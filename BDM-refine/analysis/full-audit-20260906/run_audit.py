"""Run audit reproductions using: micromamba run -n cosemu python3 -B run_audit.py.

Expected failures are evidence, not passing physics tests. Production source is
read only. The sole instrumentation records the existing binding decision.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import resource
import subprocess
import tempfile
import time

import numpy as np

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[2]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def extract(source, name, kind='subroutine'):
    pattern = rf'^\s*(?:real\s+)?{kind}\s+{name}\b.*?^\s*end\s+{kind}\s+{name}\b[^\n]*'
    match = re.search(pattern, source, re.I | re.M | re.S)
    if not match:
        raise ValueError(name)
    return match.group()


def build(destination):
    destination.mkdir(parents=True, exist_ok=True)
    source = (REPO/'PMP2linker.f90').read_text()
    modules = [re.search(rf'^module\s+{name}\b.*?^end module\s+{name}',
                         source, re.I | re.M | re.S).group()
               for name in ['BdmDuplicateRules', 'Structures']]
    modules[1] = re.sub(r'end module structures',
        'logical, allocatable :: AuditMembers(:,:)\nend module Structures', modules[1], flags=re.I)
    stub = '''module Tools
real :: Box=32.,AEXPN=.8,Om=.3,OmL=.7,ASTEP=.004,hubble=.7
integer :: NGRID=128,NROW=32,ISTEP=1,Nrealization=1
integer*8 :: Nparticles=0
character*45 :: HEADER='Controlled audit fixture'
real,allocatable :: FI(:,:,:),Xpar(:),Ypar(:),Zpar(:),VX(:),VY(:),VZ(:)
contains
real function seconds()
use omp_lib, only: omp_get_wtime
seconds=real(omp_get_wtime())
end function
real function Memory(n)
integer*8 :: n
Memory=0.
end function
end module
module LinkerList
use Structures
use Tools
use BdmDuplicateRules
contains
'''
    names = ['FindMaxima', 'GetHalo', 'FindDistinctCandidates', 'RemoveDuplicates',
             'MergeNumericalDuplicates', 'List', 'ListMaxima', 'Limits', 'EigenValues',
             'ReadParameters', 'SetParameters', 'RescaleCoords', 'AddBuffer', 'RemoveBuffer']
    routines = '\n'.join(extract(source, name) for name in names)
    tap = 'if(ee <= 0.)Ncount = Ncount +1'
    assert routines.count(tap) == 1
    routines = routines.replace(tap, 'if(ee <= 0.) AuditMembers(jp,ip)=.true.\n' + tap)
    routines += '\n' + '\n'.join(extract(source, name, 'function')
                                 for name in ['OverdenVir','OverdenAbacus','Concentration'])
    generated = '\n'.join(modules) + '\n' + stub + routines + '\nend module\n'
    path = destination/'source_cases.f90'
    path.write_text(generated + (ROOT/'audit_cases.f90').read_text())
    commands = {}
    for mode, flags in [('checked',['-O0','-g','-fcheck=all','-ffpe-trap=invalid,zero,overflow']),
                        ('optimized',['-O3'])]:
        cmd = ['gfortran',*flags,'-fopenmp','-ffree-line-length-none',
               'source_cases.f90','-o',mode]
        result = subprocess.run(cmd,cwd=destination,text=True,capture_output=True)
        commands[mode] = dict(command=cmd,returncode=result.returncode,
                              stdout=result.stdout,stderr=result.stderr)
        if result.returncode:
            raise RuntimeError(result.stderr)
    (destination/'build.json').write_text(json.dumps(commands,indent=2)+'\n')
    return commands


def run_case(build_dir, name, threads=1, mode='checked', argument=None, repeat=0):
    with tempfile.TemporaryDirectory(prefix='bdm-audit-') as scratch:
        work = Path(scratch)
        (work/'CATALOGS').mkdir()
        if name.startswith('config_'):
            config = '! header\niVirial = 2\ndLogR = 0.01\nMassMin = 1e8\n'
            if name == 'config_invalid':
                config = '! header\niVirial = 77 ! invalid\ndLogR = 0 ! invalid\nMassMin = -1 ! invalid\n'
            (work/'BDM.config').write_text(config)
        command = [str((build_dir/mode).resolve()), name]
        if argument is not None:
            command.append(str(argument))
        env = {**os.environ,'OMP_NUM_THREADS':str(threads),'OMP_DYNAMIC':'FALSE',
               'OMP_PROC_BIND':'false','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1'}
        start = time.monotonic()
        try:
            p = subprocess.run(command,cwd=work,env=env,capture_output=True,text=True,timeout=5)
            result = dict(returncode=p.returncode,stdout=p.stdout,stderr=p.stderr,
                          reached_end='AUDIT REACHED_END' in p.stdout)
        except subprocess.TimeoutExpired as error:
            result = dict(returncode=None,timeout_seconds=5,reached_end=False,
                          stdout=(error.stdout or b'').decode(),stderr=(error.stderr or b'').decode())
        result.update(case=name,mode=mode,threads=threads,argument=argument,repeat=repeat,
                      elapsed_seconds=time.monotonic()-start)
        if (work/'peaks.bin').is_file():
            result['peak_array_sha256']=sha(work/'peaks.bin')
            xyz=np.fromfile(work/'peaks.bin',dtype=np.float32).reshape(3,-1).T
            xyz=xyz[np.lexsort(xyz.T)]
            result['sorted_peak_array_sha256']=hashlib.sha256(xyz.tobytes()).hexdigest()
        if (work/'before.bin').is_file() and (work/'after.bin').is_file():
            before=np.fromfile(work/'before.bin',dtype=np.float32).reshape(6,-1)
            after=np.fromfile(work/'after.bin',dtype=np.float32).reshape(6,-1)
            result['roundtrip']=dict(changed_values_per_component=np.sum(before!=after,axis=1).tolist(),
                max_abs_change_per_component=np.max(np.abs(before.astype(float)-after),axis=1).tolist())
        if (work/'members.bin').is_file():
            raw=(work/'members.bin').read_bytes()
            phase=np.frombuffer(raw[:24000],dtype=np.float32).reshape(6,1000).T.astype(float)
            bound=np.frombuffer(raw[24000:],dtype=np.int32).astype(bool)
            retained=phase[bound]
            # Independent isolated Newtonian pair potential, excluding self.
            # G matches the finder so this tests unbinding, not constants.
            dr=retained[:,None,:3]-retained[None,:,:3]
            radius=np.linalg.norm(dr,axis=2)
            np.fill_diagonal(radius,np.inf)
            potential=4.333e-9*1.e10/.8*np.sum(1/radius,axis=1)
            hubble_a=100*np.sqrt(.3/.8**3+.7)*.8
            velocity=retained[:,3:]-retained[:,3:].mean(axis=0)+hubble_a*(retained[:,:3]-[5,5,5])
            energy=.5*np.sum(velocity**2,axis=1)-potential
            result['independent_bound_set_check']=dict(bound_count=int(bound.sum()),
                positive_energy_count=int(np.sum(energy>0)),minimum_energy_km2_s2=float(energy.min()),
                maximum_energy_km2_s2=float(energy.max()),
                interpretation='All retained particles are tested against the retained set only; self-potential excluded')
            all_dist=np.linalg.norm(phase[:,None,:3]-phase[None,:,:3],axis=2)
            np.fill_diagonal(all_dist,np.inf)
            result['independent_pair_potential_energy']=float(.5*4.333e-9*1.e20/.8*np.sum(1/all_dist))
            result['independent_vrms_all']=float(np.sqrt(np.mean(np.sum(phase[:,3:]**2,axis=1))))
            result['independent_vrms_bound']=float(np.sqrt(np.mean(np.sum(retained[:,3:]**2,axis=1))))
        return result


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--build-only',action='store_true')
    args=parser.parse_args()
    resource.setrlimit(resource.RLIMIT_CORE,(0,0))
    build_dir=ROOT/'work/cases-build'
    compilation=build(build_dir)
    if args.build_only:
        return
    results=[]
    for name in ['peak_sparse','peak_empty','peak_dynamic_range','peak_abacus',
                 'config_plain','config_invalid','shape_null_seed','shape_zero','shape_psd','shape_orthogonal',
                 'halo_central','halo_compact','halo_onepass','halo_disjoint',
                 'buffer_face','buffer_corner','buffer_search','restore_roundtrip',
                 'centering_near','centering_far']:
        result=run_case(build_dir,name)
        results.append(result)
        print(name,result['returncode'],result['reached_end'],flush=True)
    for mode in ['checked','optimized']:
        for name in ['shape_null_seed','shape_zero','halo_central','halo_compact','concentration_nan']:
            if mode=='checked' and name!='concentration_nan':
                continue
            results.append(run_case(build_dir,name,mode=mode))
    for spacing in [.04,.02,.01,.005]:
        results.append(run_case(build_dir,'halo_so',argument=spacing))
    for threads in [2,4,8,16]:
        results.append(run_case(build_dir,'peak_sparse',threads=threads))
    for threads in [1,2,4,8]:
        for repeat in range(3):
            results.append(run_case(build_dir,'peak_plateau',threads=threads,repeat=repeat))
    report=dict(checked_at_utc=datetime.now(timezone.utc).isoformat(),
        baseline_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=REPO,text=True).strip(),
        source_sha256={str(p.relative_to(REPO)):sha(p) for p in
            [REPO/'PMP2linker.f90',REPO/'PMP2bdm.f90',REPO/'PMP2mod_density.f90',
             REPO/'PMP2mod_tools.f90',ROOT/'audit_cases.f90',Path(__file__)]},
        instrumentation='GetHalo existing ee<=0 decision copied to AuditMembers; no arithmetic changed',
        compilation=compilation,results=results)
    (ROOT/'unit-results.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print('Recorded',len(results),'audit experiments; failures are preserved as findings.')


if __name__=='__main__':
    main()
