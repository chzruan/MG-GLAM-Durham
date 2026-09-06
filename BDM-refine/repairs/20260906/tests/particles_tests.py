"""Run with micromamba run -n cosemu python3 -B particles_tests.py.

Compile actual repaired particle routines with independent numerical oracles.
All temporary modules, binaries and process files are removed on completion.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import tempfile
import time

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]


def extract(source, name, kind='subroutine'):
    pattern = rf'^\s*(?:pure\s+)?(?:real\*8\s+)?{kind}\s+{name}\b.*?^\s*end\s+{kind}\s+{name}\b[^\n]*'
    match = re.search(pattern, source, re.I | re.M | re.S)
    if not match:
        raise ValueError(name)
    return match.group()


def build(destination, source, mode, compiler="gfortran"):
    structures = re.search(r'^module\s+Structures\b.*?^end module\s+Structures',
                           source, re.I | re.M | re.S).group()
    stub = '''module Tools
real :: Box=32.,AEXPN=.8
integer :: NGRID=128,NROW=32
integer*8 :: Nparticles=0,memoryWords=0
real,allocatable :: Xpar(:),Ypar(:),Zpar(:),VX(:),VY(:),VZ(:)
contains
real function seconds()
use omp_lib, only: omp_get_wtime
seconds=real(omp_get_wtime())
end function
real function Memory(n)
integer*8 :: n
memoryWords=memoryWords+n
Memory=real(dble(memoryWords)*4.d0/1024.d0**3)
end function
end module
module LinkerList
use Structures
use Tools
contains
'''
    names = ['FindDistinctCandidates', 'List', 'Limits', 'RescaleCoords',
             'AddBuffer', 'RemoveBuffer', 'PrepareParticleSearch', 'SizeList', 'BdmParticlePosition']
    generated = structures + '\n' + stub + '\n'.join(extract(source, name) for name in names)
    generated += '\n' + extract(source, 'BdmParticleCoordinate', 'function')
    generated += '\nend module\n' + (HERE/'particles_cases.f90').read_text()
    (destination/'particle_cases.f90').write_text(generated)
    flags = ['-O0', '-g', '-fcheck=all', '-ffpe-trap=invalid,zero,overflow'] if mode == 'checked' else ['-O3']
    if compiler == 'ifx':
        flags = ['-O0', '-g', '-check', 'bounds', '-fpe0', '-fp-model', 'precise'] if mode == 'checked' else [
            '-O3', '-fp-model', 'fast=1', '-march=core-avx2', '-mfma', '-ftz', '-unroll']
        command = ['ifx', *flags, '-qopenmp', '-extend-source', 'particle_cases.f90', '-o', mode]
    else:
        command = ['gfortran', *flags, '-fopenmp', '-ffree-line-length-none', 'particle_cases.f90', '-o', mode]
    compile_env = dict(os.environ)
    if compiler == 'ifx':
        compile_env.pop('LIBRARY_PATH', None)
        compile_env['LD_LIBRARY_PATH'] = os.environ['BDM_AUDIT_NATIVE_LIBS']
    completed = subprocess.run(command, cwd=destination, env=compile_env, text=True, capture_output=True)
    if completed.returncode:
        raise RuntimeError(completed.stderr)
    return dict(command=command, stdout=completed.stdout, stderr=completed.stderr)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--benchmark', action='store_true')
    parser.add_argument('--compiler', choices=['gfortran', 'ifx'], default='gfortran')
    parser.add_argument('--output', type=Path, default=HERE/'particles_results.json')
    args = parser.parse_args()
    source = (REPO/'PMP2linker.f90').read_text()
    cases = ['centering_near', 'centering_far', 'centering_one', 'centering_empty',
             'buffer_face', 'buffer_corner', 'buffer_search', 'buffer_query', 'buffer_halfbox', 'buffer_smallbox',
             'restore_roundtrip', 'restore_buffered', 'list_check']
    records = []
    compilation = {}
    with tempfile.TemporaryDirectory(prefix='bdm-particle-repairs-') as scratch:
        work = Path(scratch)
        for mode in ['checked', 'optimized']:
            compilation[mode] = build(work, source, mode, args.compiler)
            for threads in [1, 2, 4, 8]:
                for case in (['list_benchmark'] if args.benchmark else cases):
                    command = [str(work/mode), case]
                    if args.benchmark:
                        command.append('128')
                    env = {**os.environ, 'OMP_NUM_THREADS': str(threads), 'OMP_DYNAMIC': 'FALSE',
                           'OMP_PROC_BIND': 'false', 'OPENBLAS_NUM_THREADS': '1'}
                    if args.compiler == 'ifx':
                        env['LD_LIBRARY_PATH'] = os.environ['BDM_AUDIT_NATIVE_LIBS']
                    started = time.monotonic()
                    result = subprocess.run(command, cwd=work, env=env, capture_output=True, text=True, timeout=30)
                    record = dict(case=case, mode=mode, threads=threads, returncode=result.returncode,
                                  elapsed_seconds=time.monotonic()-started, stdout=result.stdout, stderr=result.stderr)
                    records.append(record)
                    if result.returncode or 'PARTICLES_TEST_PASS' not in result.stdout:
                        print(json.dumps(record, indent=2), flush=True)
                        raise RuntimeError(f'{case} failed at {threads} threads in {mode}')
    report = dict(compiler=args.compiler, source_sha256=hashlib.sha256(source.encode()).hexdigest(),
                  test_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                               for p in [Path(__file__), HERE/'particles_cases.f90']},
                  compilation=compilation, experiments=records, passed=len(records))
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(f'Passed {len(records)} independent particle repair experiments.')


if __name__ == '__main__':
    main()
