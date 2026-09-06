"""Independent exact-membership and LAPACK-reference numerical regressions.

Run with micromamba run -n cosemu python3 -B this_script.py.
Only actual source routines are compiled; no historical audit result is changed.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import tempfile

os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parents[2]


def routine(source, name, kind='subroutine'):
    pattern = rf'^\s*(?:real\s+)?{kind}\s+{name}\b.*?^\s*end\s+{kind}\s+{name}\b[^\n]*'
    return re.search(pattern, source, re.I | re.M | re.S).group()


DRIVER = '''program regression
use LinkerList
use, intrinsic :: ieee_arithmetic
implicit none
integer :: i,j,n,first,last,rows,which_case
real*8 :: matrix(3,3)
real :: eig(3),direction(3),mass,radius,vmaximum,concentration_value
character(30)::which
logical :: expected(10)
call get_command_argument(1,which)
open(13,status='scratch')
select case(trim(which))
case('shape')
  read(*,*)n
  do i=1,n
    read(*,*)matrix
    call EigenValues(matrix,direction,eig)
    write(*,'(6es25.16)')eig,direction
  enddo
case('concentration')
  read(*,*)n
  do i=1,n
    read(*,*)mass,radius,vmaximum
    concentration_value=Concentration(mass,radius,vmaximum)
    write(*,'(es25.16)')concentration_value
  enddo
case('invalid_concentration')
  vmaximum=ieee_value(1.,ieee_quiet_nan)
  if(Concentration(1.e13,500.,vmaximum)/=0.)error stop 'NaN concentration'
  vmaximum=ieee_value(1.,ieee_positive_inf)
  if(Concentration(1.e13,500.,vmaximum)/=0.)error stop 'Inf concentration'
  if(Concentration(-1.,500.,100.)/=0.)error stop 'negative mass'
  if(Concentration(1.e13,0.,100.)/=0.)error stop 'zero radius'
  if(Concentration(1.e13,500.,0.)/=0.)error stop 'unresolved Vmax'
  print *, 'INVALID INPUTS TERMINATE WITH DEFINED UNRESOLVED VALUE'
case('duplicates')
  Nmaxima=10;MassOne=1.e7;Cell=3.5
  Nmx=-1;Nmy=-1;Nmz=-1;Nbx=11;Nby=11;Nbz=11
  allocate(Mvir(10),Rvir(10),Mtotal(10),xMaxx(10),yMaxx(10),zMaxx(10))
  allocate(VxMaxx(10),VyMaxx(10),VzMaxx(10),Lst(10),Label(-1:11,-1:11,-1:11))
  allocate(BoundParticleIds(10))
  do i=1,10
    allocate(BoundParticleIds(i)%ids(64))
    BoundParticleIds(i)%ids=[(int(64*(i-1)+j,8),j=1,64)]
  enddo
  Mvir=64.*MassOne;Rvir=.02;Mtotal=Mvir
  xMaxx=[.01,31.99,10.00,10.15,10.30,20.,20.05,25.,25.05,27.]
  yMaxx=8.;zMaxx=8.;VxMaxx=0.;VyMaxx=0.;VzMaxx=0.
  BoundParticleIds(2)%ids=BoundParticleIds(1)%ids
  ! Exact identity is authoritative even with different apertures or velocities.
  BoundParticleIds(7)%ids=BoundParticleIds(6)%ids
  VxMaxx(7)=200.;Mtotal(7)=1.1*Mtotal(6)
  BoundParticleIds(8)%ids=2*BoundParticleIds(8)%ids
  BoundParticleIds(10)%ids=BoundParticleIds(8)%ids
  ! Equal counts, first and last IDs still do not establish equal sets.
  BoundParticleIds(9)%ids=BoundParticleIds(8)%ids
  BoundParticleIds(9)%ids(32)=BoundParticleIds(8)%ids(32)+1
  call ListMaxima
  call RemoveDuplicates
  expected=[.true.,.false.,.true.,.true.,.true.,.true.,.false.,.true.,.true.,.false.]
  if(any((Mvir>0.).neqv.expected))error stop 'exact identity or distinct-host failure'
  ! The pre-existing host-mask chain must still read immutable measurements.
  Mvir(4:10)=0.;Mvir(1:3)=[1.e12,2.e12,3.e12];Rvir(1:3)=.2
  xMaxx(1:3)=[1.,1.15,1.30]
  call ListMaxima
  call RemoveDuplicates
  if(any(Mvir(1:2)/=0.).or.Mvir(3)/=3.e12)error stop 'host-mask chain'
  print *, 'EXACT IDENTITIES, PERIODIC COPIES, DISJOINT HOSTS AND IMMUTABLE HOST MASK PASS'
case('duplicate_sort')
  read(*,*)Nmaxima
  MassOne=1.;allocate(Mvir(Nmaxima),BoundParticleIds(Nmaxima));Mvir=10.
  do i=1,Nmaxima
    read(*,*)n
    allocate(BoundParticleIds(i)%ids(n))
    read(*,*)BoundParticleIds(i)%ids
  enddo
  call MergeNumericalDuplicates
  do i=1,Nmaxima
    if(Mvir(i)>0.)write(*,'(i10)')i
  enddo
case('write_valid','write_nonfinite','write_zero_potential','write_wrong_mass')
  Nmaxima=1;MassOne=1.e10;MassMin=0.
  Xleft=0.;Yleft=0.;Zleft=0.;Xright=32.;Yright=32.;Zright=32.
  allocate(Mvir(1),Mtotal(1),Rvir(1),xMaxx(1),yMaxx(1),zMaxx(1),VxMaxx(1),VyMaxx(1),VzMaxx(1))
  allocate(EkinM(1),EpotM(1),VmaxM(1),RmaxM(1),Xoff(1),LambdaM(1),RadRms(1))
  allocate(Axba(1),Axca(1),Xax(1),Yax(1),Zax(1),BoundParticleIds(1),HaloStatus(1))
  allocate(BoundParticleIds(1)%ids(20));BoundParticleIds(1)%ids=[(int(j,8),j=1,20)]
  Mvir=20.*MassOne;Mtotal=Mvir;Rvir=.2;xMaxx=5.;yMaxx=5.;zMaxx=5.
  VxMaxx=0.;VyMaxx=0.;VzMaxx=0.;EkinM=1.e14;EpotM=1.e15
  VmaxM=0.;RmaxM=0.;Xoff=0.;LambdaM=0.;RadRms=0.;Axba=0.;Axca=0.
  Xax=1.;Yax=0.;Zax=0.;HaloStatus=4
  if(trim(which)=='write_nonfinite')Mvir=ieee_value(1.,ieee_quiet_nan)
  if(trim(which)=='write_zero_potential')EpotM=0.
  if(trim(which)=='write_wrong_mass')Mvir=2.*Mvir
  call BeginCataloguePublication('catalogue.dat')
  call WriteFiles
  if(Nhalo/=1)error stop 'writer lost valid row'
case default
  error stop 'unknown core test'
end select
end program
'''


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=ROOT/'core-results.json')
    args = parser.parse_args()
    source = (REPO/'PMP2linker.f90').read_text()
    modules = [re.search(rf'^module\s+{name}\b.*?^end module\s+{name}',
                        source, re.I | re.M | re.S).group()
               for name in ['BdmDuplicateRules', 'Structures']]
    stub = '''module Tools
real :: Box=32.,AEXPN=.8
contains
real function seconds()
call cpu_time(seconds)
end function
end module
module LinkerList
use Structures
use Tools
use BdmDuplicateRules
contains
'''
    names = ['MergeNumericalDuplicates', 'RemoveDuplicates', 'ListMaxima', 'Limits', 'EigenValues', 'WriteFiles',
             'BeginCataloguePublication', 'PublishCatalogue']
    generated = '\n'.join(modules)+'\n'+stub+'\n'.join(routine(source, name) for name in names)
    generated += routine(source, 'Concentration', 'function')+'\nend module\n'+DRIVER
    rng = np.random.default_rng(62813)
    matrices = [np.zeros((3, 3)), np.ones((3, 3)), np.eye(3),
                np.array([[2., -1., 0.], [-1., 2., 0.], [0., 0., .5]])]
    for scale in [1.e-20, 1.e-10, 1., 1.e10, 1.e20]:
        for _ in range(25):
            a = rng.normal(size=(3, 3));matrices.append(scale*(a@a.T))
    shape_input = str(len(matrices))+'\n'+'\n'.join(' '.join(map(repr, a.T.ravel().tolist())) for a in matrices)+'\n'
    c0 = 2.162581587
    expected_concentrations = np.array([3., 5., 10., 30., 100., 1000.])
    norm = (np.log1p(c0)-c0/(1+c0))/c0
    mass, radius, scale_factor = np.float32(1.e14), np.float32(1000.), float(np.float32(.8))
    virial2 = 4.333e-6*float(mass)/(float(radius)*scale_factor)
    speeds = np.sqrt(virial2*norm*expected_concentrations /
                     (np.log1p(expected_concentrations)-expected_concentrations/(1+expected_concentrations))).astype(np.float32)
    concentration_input = '6\n'+'\n'.join(f'{mass} {radius} {v}' for v in speeds)+'\n'
    sets = [tuple(sorted(rng.choice(100000, size=int(rng.integers(1, 200)), replace=False)+1)) for _ in range(70)]
    sets = [sets[int(i)] for i in rng.integers(0, len(sets), size=257)]
    expected_keep = [i+1 for i, value in enumerate(sets) if value not in sets[:i]]
    sort_input = str(len(sets))+'\n'+'\n'.join(str(len(ids))+'\n'+' '.join(map(str, ids)) for ids in sets)+'\n'
    results = []
    with tempfile.TemporaryDirectory(prefix='bdm-core-') as temp:
        work = Path(temp);(work/'cases.f90').write_text(generated)
        for mode, flags in [('checked', ['-O0', '-g', '-fcheck=all', '-ffpe-trap=invalid,zero,overflow']),
                            ('optimized', ['-O3'])]:
            command = ['gfortran', *flags, '-fopenmp', '-ffree-line-length-none', 'cases.f90', '-o', mode]
            build = subprocess.run(command, cwd=work, text=True, capture_output=True)
            assert build.returncode == 0, build.stderr
            def run(case, stdin='', threads=1):
                p = subprocess.run([str(work/mode), case], input=stdin, cwd=work, text=True,
                    capture_output=True, timeout=10, env={**os.environ, 'OMP_NUM_THREADS': str(threads)})
                assert p.returncode == 0, (case, p.stdout, p.stderr)
                return p.stdout
            values = np.fromstring(run('shape', shape_input), sep=' ').reshape(-1, 6)
            residuals = []
            for matrix, value in zip(matrices, values):
                expected = np.linalg.eigvalsh(matrix)[::-1]
                scale = max(float(np.max(np.abs(expected))), 1.e-300)
                assert np.max(abs(value[:3]-expected))/scale < 2.e-7
                assert abs(np.linalg.norm(value[3:])-1) < 1.e-7
                residual = np.linalg.norm(matrix@value[3:]-value[0]*value[3:])/scale
                assert residual < 2.e-7
                residuals.append(residual)
            concentrations = np.fromstring(run('concentration', concentration_input), sep=' ')
            assert np.allclose(concentrations, expected_concentrations, rtol=2.e-6)
            invalid = run('invalid_concentration')
            outputs = [run('duplicates', threads=t) for t in [1, 2, 4, 8]]
            keep = np.fromstring(run('duplicate_sort', sort_input), sep=' ', dtype=int).tolist()
            assert keep == expected_keep
            for case in ['write_valid', 'write_zero_potential']:
                run(case)
                row = np.loadtxt(work/'catalogue.dat', ndmin=2)
                assert row.shape == (1, 24) and np.isfinite(row).all()
                assert row[0, 12] == 0. and row[0, 13] == 20.
                assert np.isclose(row[0, 9], np.sqrt(1000.), rtol=2.e-4)
                if case == 'write_zero_potential':assert row[0, 16] == 0.
            for case, message in [('write_nonfinite', 'nonfinite'), ('write_wrong_mass', 'particle membership')]:
                p = subprocess.run([str(work/mode), case], cwd=work, text=True,
                                   capture_output=True, timeout=10)
                assert p.returncode != 0 and message in p.stderr, (case, p.stdout, p.stderr)
            results.append(dict(mode=mode,compiler_command=command,shape_matrices=len(matrices),
                maximum_relative_eigenpair_residual=max(residuals),
                concentration_expected=expected_concentrations.tolist(),concentration_measured=concentrations.tolist(),
                invalid_input_result=invalid.strip(),duplicate_controls_threads=[1, 2, 4, 8],
                duplicate_sort_candidates=len(sets),duplicate_sort_survivors=len(keep),
                writer_valid_and_invalid_controls=4,passed=True))
    report = dict(validated_at_utc=datetime.now(timezone.utc).isoformat(),
        source_sha256=hashlib.sha256(source.encode()).hexdigest(),
        driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),results=results)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
