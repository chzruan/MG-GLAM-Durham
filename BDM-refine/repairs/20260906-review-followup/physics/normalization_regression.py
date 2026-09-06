"""Unit-consistent SO and historical publication-floor controls.

Run with micromamba run -n cosemu python3 -B normalization_regression.py
--compiler gfortran|ifx --output NEW_RECEIPT.json. No historical receipt changes.
Extracts actual production routines; the sole GetHalo instrumentation prints
its double threshold before use. Fixtures include every original particle in
small complete boxes, with the unused mass well outside the local search cap.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import io
import json
import math
import os
from pathlib import Path
import re
import resource
import subprocess
import tempfile

os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
V2_COMMIT = '54f53aae2aac6700d1ddb8490bebb0719f71165d'
OLD_WRITER_COMMIT = 'a8c7715'


def routine(source, name):
    pattern = (rf'^\s*(?:pure\s+)?(?:(?:real|logical)(?:\s*\*\s*8)?\s+)?'
               rf'(?:subroutine|function)\s+{name}\b.*?^\s*end\s+'
               rf'(?:subroutine|function)\s+{name}\b[^\n]*')
    match = re.search(pattern, source, re.I | re.M | re.S)
    assert match, name
    return match.group()


def at_revision(revision):
    return subprocess.check_output(['git', 'show', f'{revision}:PMP2linker.f90'],
                                   cwd=REPO, text=True)


STUB = '''module Tools
real :: Box=32.,AEXPN=.8,Om=.3,OmL=.7,ASTEP=.004,hubble=.7
integer :: NGRID=128,NROW=8,ISTEP=1,Nrealization=1
integer*8 :: Nparticles=0
character*45 :: HEADER='BDM SO normalization regression'
real,allocatable :: Xpar(:),Ypar(:),Zpar(:),VX(:),VY(:),VZ(:)
contains
real function seconds()
seconds=0.
end function
real function Memory(n)
integer*8 :: n
Memory=0.
end function
end module
module LinkerList
use Structures
use Tools
contains
'''

FLOOR_DRIVER = '''program floor_controls
use LinkerList
implicit none
integer :: i,j
integer,parameter :: counts(8)=[9,10,19,20,21,25,26,64]
real :: requested_count
read(*,*)requested_count
open(13,status='scratch')
Nmaxima=size(counts);MassOne=1.e10;MassMin=requested_count*MassOne
Xleft=0.;Yleft=0.;Zleft=0.;Xright=32.;Yright=32.;Zright=32.
allocate(Mvir(Nmaxima),Mtotal(Nmaxima),Rvir(Nmaxima),xMaxx(Nmaxima),yMaxx(Nmaxima),zMaxx(Nmaxima))
allocate(VxMaxx(Nmaxima),VyMaxx(Nmaxima),VzMaxx(Nmaxima),EkinM(Nmaxima),EpotM(Nmaxima))
allocate(VmaxM(Nmaxima),RmaxM(Nmaxima),Xoff(Nmaxima),LambdaM(Nmaxima),RadRms(Nmaxima))
allocate(Axba(Nmaxima),Axca(Nmaxima),Xax(Nmaxima),Yax(Nmaxima),Zax(Nmaxima))
allocate(BoundParticleIds(Nmaxima),HaloStatus(Nmaxima));HaloStatus=0
do i=1,Nmaxima
  allocate(BoundParticleIds(i)%ids(counts(i)))
  BoundParticleIds(i)%ids=[(int(j,8),j=1,counts(i))]
  Mvir(i)=real(dble(counts(i))*dble(MassOne))
enddo
Mtotal=Mvir;Rvir=.2;xMaxx=5.;yMaxx=5.;zMaxx=5.
VxMaxx=0.;VyMaxx=0.;VzMaxx=0.;EkinM=1.e14;EpotM=1.e15
VmaxM=0.;RmaxM=0.;Xoff=0.;LambdaM=0.;RadRms=0.;Axba=.8;Axca=.5
Xax=1.;Yax=0.;Zax=0.
PUBLICATION_START
call WriteFiles
end program floor_controls
'''


def build(work, compiler, current, v2, old_writer, env):
    records = []
    variants = {}
    modules = re.search(r'^module\s+Structures\b.*?^end module\s+Structures',
                        current, re.I | re.M | re.S).group()
    names = ['SetParameters', 'SetOverdensity', 'ValidateParameters', 'ConfigurationError',
             'GetHalo', 'BdmHaloGather', 'BdmHaloSortRadii', 'BdmHaloSortIds',
             'BdmHaloSphericalPotential', 'BdmHaloMembershipInit', 'List', 'Limits',
             'EigenValues', 'BdmParticlePosition', 'BdmParticleCoordinate',
             'WriteFiles', 'Concentration', 'BeginCataloguePublication', 'PublishCatalogue']
    driver = (HERE / 'normalization_cases.f90').read_text()
    for version, source in [('v3', current), ('v2', v2)]:
        generated = modules + '\n' + STUB + '\n'.join(routine(source, n) for n in names)
        anchor = '        grid_size=dble(Box)/NGRID'
        assert generated.count(anchor) == 1
        generated = generated.replace(anchor, "        write(*,'(a,es25.16)')'THRESHOLD ',threshold\n"+anchor)
        variants[version] = generated + '\nend module\n' + driver
    for version, writer in [('floor_old', old_writer), ('floor_v3', current)]:
        generated = modules + '\n' + STUB + routine(writer, 'WriteFiles')
        if version == 'floor_v3':
            generated += routine(current, 'BeginCataloguePublication') + routine(current, 'PublishCatalogue')
            publication = "call BeginCataloguePublication('catalogue.dat')"
        else:
            publication = "open(12,file='catalogue.dat',status='replace')"
        # The floor control intentionally isolates writer cuts; concentration
        # is irrelevant and fixed finite for both actual writer versions.
        generated += '''\nreal function Concentration(m,r,v)
real :: m,r,v
Concentration=0.
end function Concentration
end module
'''
        variants[version] = generated + FLOOR_DRIVER.replace('PUBLICATION_START', publication)
    for variant, generated in variants.items():
        path = work / f'{variant}.f90'
        path.write_text(generated)
        for mode in ['checked', 'optimized']:
            if compiler == 'gfortran':
                flags = (['-O0', '-fcheck=all', '-ffpe-trap=invalid,zero,overflow']
                         if mode == 'checked' else ['-O3', '-fno-fast-math'])
                flags += ['-fopenmp', '-ffree-line-length-none']
            else:
                flags = ['-O0', '-check', 'bounds', '-fpe0'] if mode == 'checked' else ['-O3']
                flags += ['-fp-model', 'precise', '-qopenmp', '-free']
            command = [compiler, *flags, path.name, '-o', f'{variant}-{mode}']
            result = subprocess.run(command, cwd=work, env=env, capture_output=True, text=True)
            assert result.returncode == 0, (command, result.stderr)
            records.append(dict(variant=variant, mode=mode, command=command, stderr=result.stderr,
                                generated_sha256=hashlib.sha256(generated.encode()).hexdigest()))
    return records


def sphere(n, radius):
    q = np.arange(n, dtype=float)
    z = 1.-2.*(q+.5)/n
    phi = q*2.399963229728653
    return radius*np.c_[np.sqrt(1-z*z)*np.cos(phi), np.sqrt(1-z*z)*np.sin(phi), z]


def reference_delta(omega, a, mode):
    omega, a = float(np.float32(omega)), float(np.float32(a))
    fraction = omega/(omega+(1.-omega)*a**3)
    x = fraction-1.
    # These are the four existing BDM conventions, including the historical
    # virial constant 178 (not a replacement by 18*pi*pi).
    if mode == 0:
        critical = 200.
    elif mode == 1:
        critical = 178.+82.*x-39.*x*x
    elif mode == 2:
        return 200., fraction
    else:
        critical = (178.+82.*x-39.*x*x)*200./178.
    return float(np.float32(critical/fraction)), fraction


def reference_root(n, delta, box, nrow):
    # Independent volume equation in particle-number units. This makes no
    # use of GetHalo's mass coefficient or implementation contraction loop.
    return (3.*n/(4.*math.pi*delta)*float(box)**3/int(nrow)**3)**(1./3.)


def complete_box(local, box, nrow):
    n = nrow**3
    assert len(local) <= n
    data = np.zeros((n, 6), dtype=np.float32)
    data[:len(local), :3] = np.asarray(local)+box/4.
    remaining = n-len(local)
    q = np.arange(remaining)
    side = max(1, math.ceil(remaining**(1/3)))
    xyz = np.c_[q % side, q//side % side, q//side**2]
    data[len(local):, :3] = box*(.73+.04*(xyz+.5)/side)
    return data


def run_halo(work, env, mode, name, config, local, version='v3', invalid=False):
    box, nrow, ngrid, omega, a, virial, rext, multiplier = config
    data = complete_box(local, box, max(nrow, 8))
    with tempfile.TemporaryDirectory(prefix='bdm-so-case-') as scratch:
        wd = Path(scratch)
        data.tofile(wd / 'phase.bin')
        result = subprocess.run([str(work/f'{version}-{mode}'), 'invalid_nrow' if invalid else 'halo'],
                                input=' '.join(map(str, config))+'\n', cwd=wd, env=env,
                                capture_output=True, text=True, timeout=20)
        if invalid:
            assert result.returncode != 0 and 'Non-positive scale' in result.stderr, result.stderr
            assert not (wd/'catalogue.dat').exists()
            return dict(case=name, version=version, mode=mode, rejected=True)
        assert result.returncode == 0, (name, version, mode, result.stdout, result.stderr)
        lines = result.stdout.splitlines()
        def fields(prefix):
            return next(line.split()[1:] for line in lines if line.strip().startswith(prefix))
        threshold = float(fields('THRESHOLD')[0])
        scales = list(map(float, fields('SCALES')))
        values = list(map(float, fields('VALUES')))
        ids = list(map(int, fields('IDS')))
        header = (wd/'catalogue.dat').read_text().splitlines()
        assert header[0].endswith(f'[BDM finder {version}]'), header[0]
        assert len(header) == 9, (name, header, values)
        catalogue = np.loadtxt(io.StringIO('\n'.join(header)), skiprows=8, ndmin=2)
        assert catalogue.shape == (1, 24) and np.isfinite(catalogue).all()
        assert np.isfinite(values).all() and values[0] == 0., (name, values)
        assert ids == sorted(set(ids)) and all(1 <= q <= nrow**3 for q in ids)
        assert math.isclose(values[1], len(ids)*scales[1], rel_tol=7e-8), (name, values[1], len(ids)*scales[1])
        assert catalogue[0, 13] == len(ids)
        delta, fraction = reference_delta(omega, a, virial)
        assert scales[4] == delta, (name, scales[4], delta)
        if version == 'v3':
            expected = (4.*math.pi/3.)*delta*scales[1]*nrow**3/box**3
            assert math.isclose(threshold, expected, rel_tol=3e-15), (name, threshold, expected)
        radius = np.linalg.norm(data[:, :3].astype(float)-box/4., axis=1)
        sorted_radius = np.sort(radius)
        roots = np.cbrt(np.arange(1, len(data)+1)*scales[1]/threshold)
        # Explicit intervals between successive particles identify every
        # discrete SO crossing; select the outermost valid one inside the cap.
        upper = np.r_[sorted_radius[1:], np.inf]
        valid = ((np.arange(len(data)) >= 9) & (roots >= sorted_radius) & (roots < upper) &
                 (roots <= min(15.*scales[5], float(np.nextafter(np.float32(box/2), np.float32(-1))))))
        rso = float(max(roots[valid]))
        expected_ids = (np.flatnonzero(radius <= rso)+1).tolist()
        assert ids == expected_ids, (name, len(ids), len(expected_ids))
        grid = box/ngrid
        aperture = rso+grid*min(float(np.float32(rext))/(rso/grid)**float(np.float32(.2)), .75)
        assert math.isclose(values[3], aperture, rel_tol=7e-8), (name, values[3], aperture)
        assert math.isclose(values[2], np.count_nonzero(radius <= aperture)*scales[1], rel_tol=7e-8)
        if rext == 0.:
            density_ratio = values[2]/((4.*math.pi/3.)*values[3]**3)/(scales[1]*nrow**3/box**3)
            if version == 'v3':
                assert math.isclose(density_ratio, delta, rel_tol=2.5e-7), (name, density_ratio, delta)
        else:
            density_ratio = None
        assert 0. <= values[15] <= values[14] <= 1.
        assert math.isclose(np.linalg.norm(values[16:19]), 1., rel_tol=2e-7)
        return dict(case=name, version=version, mode=mode, config=config, full_box_particles=len(data),
                    input_phase_sha256=hashlib.sha256(data.tobytes()).hexdigest(), threshold=threshold,
                    stored_mass=scales[1], delta_mean=delta, matter_fraction=fraction,
                    expected_so_radius=rso, measured_aperture=values[3], bound_count=len(ids),
                    aperture_count=round(values[2]/scales[1]), measured_delta_mean=density_ratio,
                    ids_sha256=hashlib.sha256(np.asarray(ids, dtype='<i8').tobytes()).hexdigest(),
                    header_lines=8, columns=24, passed=True)


def run_floors(work, env, mode):
    records = []
    for cut, expected in [(0., [20,21,25,26,64]), (20., [20,21,25,26,64]),
                          (25., [25,26,64]), (25.5, [26,64])]:
        for variant in ['floor_old', 'floor_v3']:
            with tempfile.TemporaryDirectory(prefix='bdm-floor-case-') as scratch:
                result = subprocess.run([str(work/f'{variant}-{mode}')], input=f'{cut}\n',
                                        cwd=scratch, env=env, text=True, capture_output=True, timeout=10)
                assert result.returncode == 0, result.stderr
                table = np.loadtxt(Path(scratch)/'catalogue.dat', ndmin=2)
                assert table.shape == (len(expected), 24) and np.isfinite(table).all()
                actual = table[:, 13].astype(int).tolist()
                assert actual == expected, (variant, cut, actual, expected)
                records.append(dict(case='publication_floor', version=variant, mode=mode,
                                    requested_mass_in_particles=cut, published_counts=actual, passed=True))
    return records


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--compiler', choices=['gfortran', 'ifx'], default='gfortran')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    current = (REPO/'PMP2linker.f90').read_text()
    v2 = at_revision(V2_COMMIT)
    old_writer = at_revision(OLD_WRITER_COMMIT)
    env = {**os.environ, 'OMP_NUM_THREADS':'1', 'OMP_DYNAMIC':'FALSE'}
    if args.compiler == 'ifx' and env.get('BDM_AUDIT_NATIVE_LIBS'):
        env['LD_LIBRARY_PATH'] = env['BDM_AUDIT_NATIVE_LIBS']
        env.pop('LIBRARY_PATH', None)
    records = []
    with tempfile.TemporaryDirectory(prefix='bdm-so-build-') as scratch:
        work = Path(scratch)
        compilation = build(work, args.compiler, current, v2, old_writer, env)
        for mode in ['checked', 'optimized']:
            for omega in [.15, .3089, 1.]:
                for a in [.25, .5, 1.]:
                    for virial in range(4):
                        delta, _ = reference_delta(omega, a, virial)
                        local = sphere(64, .2*reference_root(64, delta, 32., 8))
                        records.append(run_halo(work, env, mode, 'overdensity_modes',
                                               [32.,8,128,omega,a,virial,0.,1.], local))
            for box, nrow, ngrid, multiplier in [(64.,8,128,1.), (32.,16,128,1.),
                                               (32.,8,64,1.), (32.,8,128,1.000001),
                                               (32.,8,128,4.)]:
                local = sphere(64, .2*reference_root(64, 200., box, nrow))
                records.append(run_halo(work, env, mode, 'box_particle_mesh_storage_scales',
                                       [box,nrow,ngrid,.3089,.8,2,0.,multiplier], local))
            r210 = reference_root(210, 200., 32., 8)
            two_roots = np.r_[sphere(10, .05*r210), sphere(200, .82*r210)]
            for ngrid in [64,128]:
                rec = run_halo(work, env, mode, 'outermost_two_crossings',
                               [32.,8,ngrid,.3089,.8,2,0.,1.], two_roots)
                assert rec['bound_count'] == 210
                records.append(rec)
            r256 = reference_root(256, 200., 32., 8)
            # Stored particle mass is obtained from the actual SetParameters
            # in an earlier run, avoiding a duplicate constructor in Python.
            mass = next(r['stored_mass'] for r in records if r.get('config') ==
                        [32.,8,128,.3089,1.,2,0.,1.])
            old_r256 = (256.*mass/(1.150e12*float(np.float32(.3089))*200.))**(1./3.)
            edge = (reference_root(257, 200., 32., 8)+old_r256)/2.
            assert edge > reference_root(257, 200., 32., 8) and edge < old_r256
            local = np.r_[sphere(256, .2*r256), [[edge,0.,0.]]]
            for version, expected in [('v2',257), ('v3',256)]:
                rec = run_halo(work, env, mode, 'catalogue_semantic_edge',
                               [32.,8,128,.3089,.8,2,0.,1.], local, version=version)
                assert rec['bound_count'] == expected, rec
                records.append(rec)
            r64 = reference_root(64, 200., 32., 8)
            # One separate cold row contributes only to the preserved Rext
            # aperture, beyond even the unextended 65-particle SO root.
            local = np.r_[sphere(64,.2*r64), [[r64+.05,0.,0.]]]
            rec = run_halo(work, env, mode, 'preserved_rext_aperture',
                           [32.,8,128,.3089,.8,2,.5,1.], local)
            assert rec['bound_count'] == 64 and rec['aperture_count'] == 65, rec
            records.append(rec)
            records.append(run_halo(work, env, mode, 'zero_nrow_guard',
                                    [32.,0,128,.3089,.8,2,0.,1.], sphere(64,.2*r64), invalid=True))
            records.extend(run_floors(work, env, mode))
    report = dict(validated_at_utc=datetime.now(timezone.utc).isoformat(), compiler=args.compiler,
                  source_sha256=hashlib.sha256(current.encode()).hexdigest(),
                  v2_commit=V2_COMMIT, v2_source_sha256=hashlib.sha256(v2.encode()).hexdigest(),
                  old_writer_commit=subprocess.check_output(['git','rev-parse',OLD_WRITER_COMMIT],
                                                            cwd=REPO,text=True).strip(),
                  old_writer_sha256=hashlib.sha256(routine(old_writer,'WriteFiles').encode()).hexdigest(),
                  test_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                               for p in [Path(__file__),HERE/'normalization_cases.f90']},
                  instrumentation='Print GetHalo double threshold immediately before grid_size assignment; no arithmetic replacement.',
                  compilation=compilation, cases=len(records), results=records, passed=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(compiler=args.compiler, cases=len(records), passed=True, output=str(args.output))))


if __name__ == '__main__':
    main()
