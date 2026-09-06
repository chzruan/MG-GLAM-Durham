"""Run preserved halo algorithm tests against the v3 density convention.

Historical Python/Fortran fixtures and receipts are never edited. This adapter
changes only fixture parameters, fixture diagnostic output and Python oracles;
every compiled production routine is extracted verbatim from current source.
Use micromamba run -n cosemu python3 -B this_script.py --output NEW.json.
"""
from __future__ import annotations
import argparse
import ast
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import resource
import subprocess
import sys
import tempfile

os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
LEGACY = REPO/'BDM-refine/repairs/20260906/tests'
sys.path.insert(0, str(LEGACY))
import halo_review as review
import halo_host_ties as ties
import normalization_regression as normalization

BOX = 32.
DELTA = 200.
OMEGA = float(np.float32(.3))
AEXPN = float(np.float32(.8))


def stored_mass(nrow):
    """Same existing real32 constructor used by SetParameters, explicit steps."""
    scale = np.float32(BOX/nrow)
    density_scale = np.float32(np.float32(2.774e11)*np.float32(scale**3))
    return float(np.float32(np.float32(OMEGA)*density_scale))


MASS = stored_mass(16)
THRESHOLD = (4.*math.pi/3.)*DELTA*MASS*(16./BOX)**3
TRANSFORMS = {}
COMPILER = 'gfortran'
SOURCE_PATH = REPO/'PMP2linker.f90'


def sha(value):
    if isinstance(value, Path):
        value = value.read_bytes()
    if isinstance(value, str):
        value = value.encode()
    return hashlib.sha256(value).hexdigest()


class FixtureConvention(ast.NodeTransformer):
    """Adapt fixture masses and SO references, never production Fortran."""
    def __init__(self, core=False):
        self.core = core
        self.mass_changes = 0
        self.threshold_changes = 0
        self.radius_changes = 0

    def visit_Assign(self, node):
        if any(isinstance(t, ast.Name) and t.id in ('threshold', 'THRESHOLD') for t in node.targets):
            if any(isinstance(n, ast.Constant) and n.value == 1.150e12 for n in ast.walk(node.value)):
                node.value = (ast.Call(ast.Name('v3_threshold', ast.Load()),
                                       [ast.Name('name', ast.Load()), ast.Name('mass', ast.Load())], [])
                              if self.core else ast.Name('V3_THRESHOLD', ast.Load()))
                self.threshold_changes += 1
                return node
        return self.generic_visit(node)

    def visit_Constant(self, node):
        if not self.core and isinstance(node.value, float) and node.value == 1.e12:
            self.mass_changes += 1
            return ast.copy_location(ast.Name('V3_MASS', ast.Load()), node)
        return node

    def visit_Assert(self, node):
        # The analytic n(<r) proportional to r profile remains a nontrivial
        # interior crossing. Its v3 radius follows from its actual box mean.
        if self.core and ast.unparse(node.test) == 'abs(so - 1.0) < 3e-05':
            node.test = ast.parse('abs(so-v3_profile_root(len(phase))) < 3.e-5', mode='eval').body
            self.radius_changes += 1
            return node
        return self.generic_visit(node)


def load_adapted(filename, core=False):
    path = LEGACY/filename
    source = path.read_text()
    tree = ast.parse(source, filename=str(path))
    transformer = FixtureConvention(core)
    tree = ast.fix_missing_locations(transformer.visit(tree))
    assert transformer.threshold_changes > 0, filename
    if core:
        assert transformer.radius_changes == 1
    else:
        assert transformer.mass_changes > 0
    name = 'v3_'+path.stem
    spec = importlib.util.spec_from_loader(name, loader=None, origin=str(path))
    module = importlib.util.module_from_spec(spec)
    module.__file__ = str(path)
    module.V3_MASS = MASS
    module.V3_THRESHOLD = THRESHOLD
    module.v3_threshold = lambda case, mass: ((4.*math.pi/3.)*DELTA*mass*
                                             ((128. if case.startswith('so') else 16.)/BOX)**3)
    module.v3_profile_root = lambda count: math.sqrt(count/(2.*(4.*math.pi/3.)*DELTA*(128./BOX)**3))
    exec(compile(tree, str(path)+' [v3 fixture adapter]', 'exec'), module.__dict__)
    TRANSFORMS[filename] = dict(original_sha256=sha(source), adapted_ast_sha256=sha(ast.unparse(tree)),
                                mass_literals=transformer.mass_changes,
                                threshold_oracles=transformer.threshold_changes,
                                analytic_radius_assertions=transformer.radius_changes)
    return module


def production_module():
    source = SOURCE_PATH.read_text()
    names = ['GetHalo', 'BdmHaloGather', 'BdmHaloSortRadii', 'BdmHaloSortIds',
             'BdmHaloSphericalPotential', 'BdmHaloMembershipInit', 'ParametersDistinct',
             'List', 'Limits', 'EigenValues', 'BdmParticlePosition', 'BdmParticleCoordinate']
    structures = re.search(r'^module\s+Structures\b.*?^end module\s+Structures',
                           source, re.I|re.M|re.S).group()
    # The box-mean metadata describes the full box; these direct-routine
    # fixtures provide only the local rows which can reach the query aperture.
    stub = normalization.STUB.replace('NROW=8', 'NROW=16')
    code = structures+'\n'+stub+'\n'.join(normalization.routine(source, n) for n in names)+'\nend module\n'
    return code


def compile_fixture(work, program):
    generated = production_module()+program
    (work/'source.f90').write_text(generated)
    records = {}
    for mode in ['checked', 'optimized']:
        if COMPILER == 'gfortran':
            flags = (['-O0', '-g', '-fcheck=all', '-ffpe-trap=invalid,zero,overflow']
                     if mode == 'checked' else ['-O3', '-fno-fast-math', '-ffp-contract=off'])
            flags += ['-fopenmp', '-ffree-line-length-none']
        else:
            flags = ['-O0', '-g', '-check', 'bounds', '-fpe0'] if mode == 'checked' else ['-O3']
            flags += ['-fp-model', 'precise', '-qopenmp', '-extend-source']
        command = [COMPILER, *flags, 'source.f90', '-o', mode]
        p = subprocess.run(command, cwd=work, capture_output=True, text=True)
        assert p.returncode == 0, (command, p.stderr)
        records[mode] = dict(command=command, returncode=p.returncode, stdout=p.stdout, stderr=p.stderr,
                             generated_sha256=sha(generated), binary_sha256=sha(work/mode))
    return records


def replace_once(text, before, after):
    assert text.count(before) == 1, before
    return text.replace(before, after)


def build_review(work, transform_program=None):
    program = (LEGACY/'halo_review_cases.f90').read_text()
    if transform_program is not None:
        program = transform_program(program)
    program = replace_once(program, "  write(*,*) 'VALUES'", "  write(*,'(a,i6,18es25.16)') 'VALUES '")
    program = replace_once(program, "    write(*,*) 'IDS',", "    write(*,'(a,*(i0,1x))') 'IDS ',")
    program = replace_once(program, "    write(*,*)'POTENTIAL',singular,energy",
                           "    write(*,'(a,l2,es25.16)')'POTENTIAL ',singular,energy")
    program = program.replace("    write(*,*)'SORTED_IDS',", "    write(*,'(a,*(i0,1x))')'SORTED_IDS ',")
    program = replace_once(program, '  call BdmHaloMembershipInit',
                           "  if(MassOne/=Om0*(2.774e11*(Box/NROW)**3))error stop 'incoherent fixture particle mass'\n"
                           '  call BdmHaloMembershipInit')
    program = replace_once(program, '  ! Call again on the same candidate',
                           "  write(*,'(a,2i20)')'UNBIND ',HaloUnbindingPasses(1),HaloUnbindingWork(1)\n"
                           '  ! Call again on the same candidate')
    program = replace_once(program, "  write(*,*) 'EMPTY_RECALL_PASS'",
                           "  if(HaloUnbindingPasses(1)/=0.or.HaloUnbindingWork(1)/=0_8)error stop 'stale unbinding counters'\n"
                           "  write(*,*) 'EMPTY_RECALL_PASS'")
    return compile_fixture(work, program)


OLD_RUN_HALO = review.run_halo


def run_review(work, name, phase, mode, cell=1., mass=MASS, identity=True):
    assert float(np.float32(mass)) == MASS, (name, mass, MASS)
    assert identity, 'real GetHalo fixtures use actual row identities'
    result = OLD_RUN_HALO(work, name, phase, mode, cell, mass, identity=True)
    parts = next(line.split()[1:] for line in result['stdout'].splitlines()
                 if line.strip().startswith('UNBIND'))
    result['unbinding_passes'], result['unbinding_work'] = map(int, parts)
    result['particle_mass'] = MASS
    result['global_nrow'] = 16
    return result


def build_core(work):
    program = (LEGACY/'halo_cases.f90').read_text()
    program = replace_once(program, '  MassOne=1.e10; Om0=.3;', '  NROW=16; Om0=.3;')
    program = replace_once(program, '  Box=32.; NGRID=128; Cell=.5',
                           '  Box=32.; NGRID=128; Cell=.5\n  MassOne=Om0*(2.774e11*(Box/NROW)**3)')
    program = replace_once(program, '    MassOne=2.*(1.150e12*Om0)*Ovdens/n',
                           '    NROW=128\n    MassOne=Om0*(2.774e11*(Box/NROW)**3)')
    # Keep the cold/hot/mixed classifications meaningful under the changed
    # mass unit. Independent references use the actual float32 phase data.
    program = replace_once(program, '  allocate(Lst(Np),Label',
                           "  if(index(which,'shell')>0)then\n"
                           '    VX=VX*sqrt(MassOne/1.e10);VY=VY*sqrt(MassOne/1.e10);VZ=VZ*sqrt(MassOne/1.e10)\n'
                           '  endif\n  allocate(Lst(Np),Label')
    program = replace_once(program, '  do ip=1,nh\n',
                           "  write(*,'(a,2i20)')'UNBIND ',HaloUnbindingPasses(1),HaloUnbindingWork(1)\n"
                           '  do ip=1,nh\n')
    return compile_fixture(work, program), []


class FollowupPaths:
    """Redirect only the historical main's hard-coded JSON output."""
    def __init__(self, destination):
        self.destination = destination

    def __truediv__(self, name):
        return self.destination if name == 'halo_followup_results.json' else LEGACY/name


def remaining_review_controls(followup):
    """Retain the review's unique central-survivor and shell-domain probes."""
    records = []
    with tempfile.TemporaryDirectory(prefix='bdm-v3-review-') as scratch:
        work = Path(scratch)
        compilation = build_review(work, followup.program_transform)
        for mode in ['checked', 'optimized']:
            phase = np.c_[review.sphere(10,.1)+5., review.sphere(10,1.e7)]
            phase[0, :3] = 5.
            phase[0, 3:] = 0.
            phase[1:, 3:] = np.tile([[1.e7,0.,0.],[-5.e6,1.e7,0.],[-5.e6,-1.e7,0.]], (3,1))
            record = run_review(work, 'one_central_survivor', phase, mode)
            assert record['ids'] == [1]
            assert record['values'][0] == 4. and record['values'][4:8] == [0.,0.,0.,0.]
            assert record['unbinding_passes'] == 2 and record['unbinding_work'] == 11
            records.append(record)
            for radii in [[], [0.], [.2], [0.,.2], [0.,0.], [.1,.1,.2,.2,.5]]:
                stdin = str(len(radii))+'\n'+' '.join(map(str,radii))+'\n'
                p = subprocess.run([str(work/mode), 'potential'], input=stdin, text=True, capture_output=True)
                assert p.returncode == 0, (radii, p.stderr)
                line = next(line for line in p.stdout.splitlines() if line.strip().startswith('POTENTIAL'))
                _, flag, energy = line.split()
                singular = len(radii)>1 and radii[1] == 0.
                assert (flag == 'T') == singular
                expected = None
                if not singular:
                    expected = 3.*4.*math.fsum(1/max(a,b) for i,a in enumerate(radii) for b in radii[i+1:])
                    assert math.isclose(float(energy), expected, rel_tol=2e-14, abs_tol=1e-14)
                records.append(dict(case='shell_potential', mode=mode, radii=radii,
                                    singular=singular, independent_pair_energy=expected,
                                    stdout=p.stdout, stderr=p.stderr))
    assert len(records) == 14
    return dict(compilation=compilation, results=records, passed=14)


def main():
    global COMPILER, SOURCE_PATH
    parser = argparse.ArgumentParser()
    parser.add_argument('--compiler', choices=['gfortran', 'ifx'], default='gfortran')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--source', type=Path, default=SOURCE_PATH,
                        help='PMP2linker.f90 in a Git worktree, read only; defaults to this worktree')
    args = parser.parse_args()
    COMPILER = args.compiler
    SOURCE_PATH = args.source.resolve()
    if SOURCE_PATH.name != 'PMP2linker.f90':
        parser.error('--source must name PMP2linker.f90 in a Git worktree')
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    if COMPILER == 'ifx':
        os.environ['LD_LIBRARY_PATH'] = os.environ['BDM_AUDIT_NATIVE_LIBS']
        os.environ.pop('LIBRARY_PATH', None)
    review.build = build_review
    review.run_halo = run_review
    review.SOURCE = SOURCE_PATH.read_text()
    ties.REPO = SOURCE_PATH.parent
    followup = load_adapted('halo_followup.py')
    # The preserved gaps module imports halo_followup by its original name.
    sys.modules['halo_followup'] = followup
    gaps = load_adapted('halo_gaps.py')
    core = load_adapted('halo_regression.py', core=True)
    core.build = build_core
    core.SOURCE = SOURCE_PATH
    historical = sorted(p for p in LEGACY.glob('halo*') if p.is_file())
    historical += [REPO/'BDM-refine/analysis/review-20260906-merged/REVIEW.md']
    preserved = {str(p.relative_to(REPO)):sha(p) for p in historical}
    source_sha = sha(SOURCE_PATH)
    argv = sys.argv[:]
    with tempfile.TemporaryDirectory(prefix='bdm-v3-suite-') as scratch:
        work = Path(scratch)
        sys.argv = ['halo_gaps.py', '--compiler', COMPILER, '--output', str(work/'gaps.json')]
        gaps.main()
        gap_report = json.loads((work/'gaps.json').read_text())
        assert gap_report['passed'] == 166
        assert gap_report['source_sha256'] == source_sha
        for record in gap_report['results']:
            if record['case'].startswith('unbinding_ladder'):
                population = record['population']
                passes = record['reference_passes']
                expected_work = passes*(population+20)//2
                assert record['unbinding_passes'] == passes
                assert record['unbinding_work'] == expected_work, record
                pairs = population//2
                expected_ids = list(range(pairs-9,pairs+1))+list(range(population-9,population+1))
                assert record['ids'] == expected_ids
                record['expected_unbinding_work'] = expected_work
        followup.HERE = FollowupPaths(work/'followup.json')
        followup.main()
        followup_report = json.loads((work/'followup.json').read_text())
        assert followup_report['passed'] == 26
        sys.argv = ['halo_regression.py', '--output', str(work/'core.json')]
        core.main()
        core_report = json.loads((work/'core.json').read_text())
        assert core_report['passed'] == 36
        expected_counts = dict(cold=64, central=64, compact=64, cold_shell=1000,
                               hot_shell=0, mixed_shell=900, bulk_shell=900)
        for record in core_report['results']:
            if record['case'] in expected_counts:
                assert record['bound_count'] == expected_counts[record['case']], record
        remaining_review = remaining_review_controls(followup)
    sys.argv = argv
    assert sha(SOURCE_PATH) == source_sha, 'production source changed during validation'
    assert all(sha(REPO/path) == digest for path, digest in preserved.items())
    report = dict(validated_at_utc=datetime.now(timezone.utc).isoformat(), compiler=COMPILER,
                  compiler_version=subprocess.check_output([COMPILER,'--version'],text=True).splitlines()[0],
                  source_sha256=source_sha, driver_sha256=sha(Path(__file__)), preserved_sha256=preserved,
                  source_path=str(SOURCE_PATH),
                  source_worktree_head=subprocess.check_output(['git','rev-parse','HEAD'],
                                                              cwd=SOURCE_PATH.parent,text=True).strip(),
                  build_helper_sha256=sha(HERE/'normalization_regression.py'),
                  fixture_convention=dict(box=BOX, delta=DELTA, omega=OMEGA, nrow_local_cases=16,
                                          stored_mass_local_cases=MASS, threshold_local_cases=THRESHOLD,
                                          nrow_radial_profile=128, stored_mass_radial_profile=stored_mass(128),
                                          local_excerpt=True),
                  fixture_adaptations=TRANSFORMS, gaps=gap_report, followup=followup_report,
                  original_algorithms=core_report, remaining_review=remaining_review, passed=242)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(f'{COMPILER}: all 242 current-source halo cases passed (166 gaps, 26 follow-up, 36 original algorithms, 14 remaining review).')


if __name__ == '__main__':
    main()
