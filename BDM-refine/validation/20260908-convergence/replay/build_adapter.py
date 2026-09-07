"""Build only campaign replay adapters from already frozen native objects.

Execute with micromamba run -n cosemu python3 -B, after loading the compiler.
The caller owns simulation builds, scheduling, input/output hashing and cleanup.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess

HERE = Path(__file__).resolve().parent
V3_SHA = '39f7e514d1db9fb3b9c7545fc7f3269b4e20b5ee2d1366945f8916b997e4cd5d'
LEGACY_SHA = '05c8d604977175568ca9e3b85d90e323c69a313781d4b2696e1d75bfc232fb5c'
OBJECTS = ['PMP2mod_tools', 'PMP2mod_fft5', 'PMP2mod_random', 'PMP2mod_density',
           'PMP2mod_power', 'PMP2mod_analyze', 'PMP2MG_subroutines', 'PMP2mod_MGbackground',
           'PMP2MGsolver_fR', 'PMP2extradof', 'PMP2MGsolver_DGP', 'PMP2MGsolver_sym',
           'PMP2MGsolver_kmf', 'PMP2MGsolver_csf']
FLAGS = ['-O3', '-g', '-traceback', '-ftz', '-unroll', '-qopenmp', '-march=core-avx2',
         '-mfma', '-shared-intel', '-mcmodel=medium', '-convert', 'big_endian', '-fp-model', 'precise']


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        while chunk := stream.read(8 * 1024**2):
            h.update(chunk)
    return h.hexdigest()


def instrument(source):
    # Two calls: valid empty catalogue and the ordinary successful path.
    assert source.count('Call WriteFiles\n') == 2
    source = source.replace('Call WriteFiles\n', 'Call WriteFiles\n      call ConvergenceDumpMembers\n')
    matches = list(re.finditer(r'^end\s+Module\s+LinkerList\s*$', source, re.I | re.M))
    assert len(matches) == 1
    start = matches[0].start()
    return source[:start] + (HERE/'membership.f90').read_text() + '\n' + source[start:]


def build(common_dir, v3_source, legacy_source, output_dir, compiler='ifx', env=None):
    common_dir, v3_source, legacy_source, output_dir = map(
        lambda p: Path(p).resolve(), (common_dir, v3_source, legacy_source, output_dir))
    assert sha(v3_source) == V3_SHA, 'Unexpected v3 production source'
    assert sha(legacy_source) == LEGACY_SHA, 'Legacy must be exact a8c7715 source'
    objects = [common_dir/(name+'.o') for name in OBJECTS]
    assert all(p.is_file() for p in objects), 'Build frozen common objects first'
    modules = [p for p in sorted(common_dir.glob('*.mod'))
               if p.name.lower() not in ['structures.mod', 'linkerlist.mod', 'bdmduplicaterules.mod']]
    assert modules
    output_dir.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ if env is None else env)
    env.update(OMP_NUM_THREADS='1', OMP_DYNAMIC='FALSE', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    if env.get('BDM_AUDIT_NATIVE_LIBS'):
        env['LD_LIBRARY_PATH'] = env['BDM_AUDIT_NATIVE_LIBS']
    env.pop('LIBRARY_PATH', None)
    report = dict(completed=False, common_dir=str(common_dir), compiler=compiler,
                  common_objects_sha256={p.name: sha(p) for p in objects},
                  common_modules_sha256={p.name: sha(p) for p in modules},
                  adapter_sha256={name: sha(HERE/name) for name in
                                  ['common_grid.f90', 'replay_entry.f90', 'membership.f90', 'build_adapter.py']},
                  variants={}, commands=[])

    def receipt():
        temporary = output_dir/'build.json.part'
        temporary.write_text(json.dumps(report, indent=2) + '\n')
        temporary.replace(output_dir/'build.json')

    def run(command, cwd):
        result = subprocess.run(list(map(str, command)), cwd=cwd, env=env, text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=300)
        report['commands'].append(dict(command=list(map(str, command)), cwd=str(cwd),
                                       returncode=result.returncode, output=result.stdout))
        receipt()
        if result.returncode:
            raise RuntimeError(result.stdout)

    run([compiler, '--version'], output_dir)
    for variant, source in [('legacy', legacy_source), ('v3', v3_source)]:
        folder = output_dir/variant
        folder.mkdir()
        for path in modules:
            shutil.copy2(path, folder/path.name)
        source_text = source.read_text()
        finder = source_text if variant == 'legacy' else instrument(source_text)
        (folder/'finder.f90').write_text(finder)
        for name in ['common_grid.f90', 'replay_entry.f90']:
            shutil.copy2(HERE/name, folder/name)
        finder_flags = list(FLAGS)
        if variant == 'legacy':
            finder_flags[-1] = 'fast=1'  # original a8c7715 compilation arithmetic
        run([compiler, *finder_flags, '-c', 'finder.f90', '-o', 'finder.o'], folder)
        run([compiler, *FLAGS, '-c', 'common_grid.f90', 'replay_entry.f90'], folder)
        binary = folder/'replay.exe'
        run([compiler, *FLAGS, '-o', binary, *objects, 'finder.o', 'common_grid.o', 'replay_entry.o'], folder)
        report['variants'][variant] = dict(binary_path=str(binary), binary_sha256=sha(binary),
            original_source_sha256=sha(source), instrumented_source_sha256=sha(folder/'finder.f90'),
            finder_flags=finder_flags, adapter_flags=FLAGS,
            membership_hook_count=0 if variant == 'legacy' else 2)
        receipt()
    assert report['common_objects_sha256'] == {p.name: sha(p) for p in objects}
    assert report['common_modules_sha256'] == {p.name: sha(p) for p in modules}
    report['completed'] = True
    receipt()
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--common-dir', type=Path, required=True)
    parser.add_argument('--v3-source', type=Path, required=True)
    parser.add_argument('--legacy-source', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--compiler', default='ifx')
    args = parser.parse_args()
    report = build(**vars(args))
    print(json.dumps({name: item['binary_path'] for name,item in report['variants'].items()}))


if __name__ == '__main__':
    main()
