"""Assert corrected BDM peaks/configuration using extracted production routines.

Run: micromamba run -n cosemu python3 -B peaks_regression.py [--output JSON]
Builds and all run files are temporary; optional evidence is one JSON file.
The empty-BDM test uses real BDM/WriteFiles and explicit fail-on-use downstream
stubs, so it checks that no particle rescaling/buffering occurs on that path.
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

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]


def extract(source, name):
    pattern = (rf'^\s*(?:pure\s+)?(?:(?:real|logical)\s+)?(?:subroutine|function)\s+{name}\b'
               rf'.*?^\s*end\s+(?:subroutine|function)\s+{name}\b[^\n]*')
    match = re.search(pattern, source, re.I | re.M | re.S)
    if match is None:
        raise ValueError(f'production routine not found: {name}')
    return match.group()


def build(destination):
    source = (REPO / 'PMP2linker.f90').read_text()
    modules = '\n'.join(re.search(rf'^module\s+{name}\b.*?^end module\s+{name}',
                                  source, re.I | re.M | re.S).group()
                        for name in ['Structures', 'BdmDuplicateRules'])
    stub = '''module Tools
real :: Box=32.,AEXPN=1.,Om=.3,OmL=.7,ASTEP=.004,hubble=.7
integer :: NGRID=16,NROW=8,ISTEP=1,Nrealization=1
integer*8 :: Nparticles=17
character*45 :: HEADER='BDM configuration/peak regression'
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
module Density
use Tools
use Structures
contains
subroutine DENSIT
logical :: opened
if (iVirial /= 2.or.Ovdens /= 200.) error stop 'BDM configuration is late'
inquire(unit=12,opened=opened)
if (.not.opened) error stop 'BDM header is late'
end subroutine
end module
module LinkerList
use Structures
use Tools
use BdmDuplicateRules
contains
'''
    names = ['BDM','ReleaseMaxima','ReadParameters','ConfigurationError','ValidateParameters',
             'SetOverdensity','SetParameters','FindMaxima','IsDensityMaximum','WriteFiles']
    generated = modules + '\n' + stub + '\n'.join(extract(source, name) for name in names)
    for name in ['AddBuffer','SizeList','List','FindDistinctCandidates','ParametersDistinct',
                 'SizeListMaxima','ListMaxima','RemoveDuplicates']:
        generated += f"\nsubroutine {name}\nerror stop 'empty BDM called {name}'\nend subroutine {name}\n"
    generated += '''
subroutine RescaleCoords(flag)
integer :: flag
error stop 'empty BDM rescaled particles'
end subroutine RescaleCoords
subroutine RemoveBuffer(n)
integer*8 :: n
error stop 'empty BDM restored unnecessary particle buffer'
end subroutine RemoveBuffer
real function Concentration(m,r,v)
real :: m,r,v
error stop 'empty catalogue evaluated halo properties'
end function Concentration
end module
'''
    (destination / 'peaks_source.f90').write_text(generated + (HERE / 'peaks_cases.f90').read_text())
    commands = []
    for mode, flags in [('checked',['-O0','-g','-fcheck=all','-ffpe-trap=invalid,zero,overflow']),
                        ('optimized',['-O3'])]:
        command = ['gfortran', *flags, '-fopenmp', '-ffree-line-length-none',
                   'peaks_source.f90', '-o', mode]
        p = subprocess.run(command, cwd=destination, text=True, capture_output=True)
        commands.append(dict(mode=mode, command=command, returncode=p.returncode,
                             stdout=p.stdout, stderr=p.stderr))
        assert p.returncode == 0, p.stderr
    return commands


def reference_peaks(density, threshold=200/3):
    """Independent 26-neighbour reference with a global index tie order."""
    n = density.shape[0]
    index = np.arange(n**3).reshape(density.shape, order='F')
    selected = density > threshold
    for dz in [-1,0,1]:
        for dy in [-1,0,1]:
            for dx in [-1,0,1]:
                if (dx,dy,dz) == (0,0,0):
                    continue
                neighbour = np.roll(density, (dx,dy,dz), axis=(0,1,2))
                neighbour_index = np.roll(index, (dx,dy,dz), axis=(0,1,2))
                selected &= ((density > neighbour) |
                             ((density == neighbour) & (index < neighbour_index)))
    flat = np.flatnonzero(selected.ravel(order='F'))
    xyz = np.array(np.unravel_index(flat,density.shape,order='F')).T * (32/n)
    return np.column_stack([xyz,density.ravel(order='F')[flat]]).astype(np.float32)


def execute(build_dir, mode, case, threads, density=None, config=None, success=True):
    with tempfile.TemporaryDirectory(prefix='bdm-peaks-case-') as scratch:
        work = Path(scratch)
        (work / 'CATALOGS').mkdir()
        if config is not None:
            (work / 'BDM.config').write_text(config)
        if density is not None:
            density.ravel(order='F').astype(np.float32).tofile(work / 'density.bin')
        env = {**os.environ, 'OMP_NUM_THREADS':str(threads), 'OMP_DYNAMIC':'FALSE',
               'OMP_PROC_BIND':'false', 'OPENBLAS_NUM_THREADS':'1'}
        command = [str(build_dir / mode),case]
        p = subprocess.run(command,cwd=work,env=env,capture_output=True,text=True,timeout=10)
        result = dict(mode=mode,case=case,threads=threads,returncode=p.returncode,
                      stdout=p.stdout,stderr=p.stderr)
        if not success:
            assert p.returncode != 0, result
            if case == 'configuration':
                assert not list((work / 'CATALOGS').iterdir()), 'invalid config opened output files'
            return result
        assert p.returncode == 0 and 'PEAKS TEST COMPLETE' in p.stdout, result
        if case == 'peaks':
            raw = (work / 'peaks.bin').read_bytes()
            n = int(np.frombuffer(raw[:4],np.int32)[0])
            overdensity = float(np.frombuffer(raw[4:8],np.float32)[0])
            actual = np.frombuffer(raw[8:],np.float32).reshape(4,n).T
            expected = reference_peaks(density)
            assert n == len(expected) and overdensity == 200, result
            assert np.array_equal(actual,expected), (actual,expected)
            result.update(count=n, peak_sha256=hashlib.sha256(raw).hexdigest())
        else:
            catalogues = list((work / 'CATALOGS').glob('Catshort*.DAT'))
            assert len(catalogues) == 1, catalogues
            header = catalogues[0].read_text().splitlines()
            assert len(header) == 8, header
            assert header[0].endswith('[BDM finder v2]'), header
            omega_lambda = re.search(r'Omega_L=\s*([\d.]+)',header[3])
            assert omega_lambda and float(omega_lambda[1]) == .7, header[3]
            assert 'Mbound' in header[7] and 'MajorAxis' in header[7], header[7]
            result.update(catalogue=catalogues[0].name,header=header)
            if case == 'configuration':
                result['settings'] = np.loadtxt(work / 'settings.txt').tolist()
        return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    resource.setrlimit(resource.RLIMIT_CORE,(0,0))
    fixtures = {}
    for name in ['sparse','empty','dynamic_range','plateau','periodic_tie','skewed','uniform','random']:
        d = np.zeros((16,16,16),dtype=np.float32)
        if name == 'sparse': d[7,7,7] = 1000
        if name == 'dynamic_range': d[3,3,3] = 1e10; d[11,11,11] = 100
        if name == 'plateau': d[3:13,3:13,3:13] = 1000
        if name == 'periodic_tie': d[0,7,7] = 1000; d[-1,7,7] = 1000
        if name == 'skewed': d[1::2,1::2,0] = 1000
        if name == 'uniform': d[:] = 1000
        if name == 'random': d[:] = np.random.default_rng(4156).integers(0,400,size=d.shape)
        fixtures[name] = d
    invalid = ['iVirial=77','dLogR=0','dLogP=-1','NradP=0','MassMin=-1',
               'Rext=-1','SlopeR=-1','dLogR=NaN','dLogP=Inf','MassMin=garbage',
               'iVirial=2 trailing','NradP=2,3','dLogR=','MisspelledOption=1','missing equals',
               'iVirial=2'+' '*1015+'invalid']
    results = []
    with tempfile.TemporaryDirectory(prefix='bdm-peaks-build-') as scratch:
        build_dir = Path(scratch)
        compilation = build(build_dir)
        for mode in ['checked','optimized']:
            for name,density in fixtures.items():
                hashes = set()
                for threads in [1,2,4,8,16]:
                    result = execute(build_dir,mode,'peaks',threads,density=density)
                    result['fixture'] = name
                    hashes.add(result['peak_sha256'])
                    results.append(result)
                assert len(hashes) == 1, f'{mode}/{name} is thread dependent'
            nonfinite = fixtures['sparse'].copy(); nonfinite[0,0,0] = np.nan
            results.append(execute(build_dir,mode,'peaks',1,density=nonfinite,success=False))
            for invalid_config in invalid:
                result = execute(build_dir,mode,'configuration',1,config=invalid_config+'\n',success=False)
                result['invalid_config'] = invalid_config
                results.append(result)
            for ivirial,label in enumerate(['W','V','M','A']):
                for annotated in [False,True]:
                    config = f'iViRiAl={ivirial}\ndLoGr\t=\t0.01\nMinMass=1.e8\nNradP=100\n'
                    if annotated:
                        config = '! comment\n'+'\n'.join(line+' ! inline comment' for line in config.splitlines())+'\n'
                    result = execute(build_dir,mode,'configuration',1,config=config)
                    settings = result['settings']
                    assert settings[:2] == [ivirial,100]
                    assert abs(settings[2]-.01) < 1e-8 and settings[4] == 1e8
                    om = float(np.float32(.3))
                    xx = om-1
                    expected = [200/om,(178+82*xx-39*xx**2)/om,200,
                                (178+82*xx-39*xx**2)/om*200/178][ivirial]
                    assert np.isclose(settings[7],expected,rtol=1e-7)
                    assert settings[8] == (2 if ivirial == 2 else 1)
                    assert result['catalogue'] == f'Catshort{label}.0001.0001.DAT'
                    result['annotated_config'] = annotated
                    results.append(result)
            results.append(execute(build_dir,mode,'configuration',1))
            for threads in [1,8]:
                results.append(execute(build_dir,mode,'empty_bdm',threads,config='iVirial=2\n'))
    report = dict(checked_at_utc=datetime.now(timezone.utc).isoformat(),
                  source_sha256={str(p.relative_to(REPO)):hashlib.sha256(p.read_bytes()).hexdigest()
                                 for p in [REPO/'PMP2linker.f90',Path(__file__),HERE/'peaks_cases.f90']},
                  compilation=compilation,asserted_experiments=len(results),results=results)
    if args.output:
        args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(f'PASS: {len(results)} asserted peak/configuration/empty-BDM experiments (checked and optimized).')


if __name__ == '__main__':
    main()
