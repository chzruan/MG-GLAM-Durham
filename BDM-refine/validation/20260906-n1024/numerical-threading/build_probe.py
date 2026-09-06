"""Link a diagnostic entry to frozen native objects; never rebuild the finder.

Use micromamba run -n cosemu python3 -B, with Intel modules loaded and
BDM_AUDIT_NATIVE_LIBS captured before entering the Python environment.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tarfile

ROOT = Path(__file__).resolve().parent
OBJECTS = ['PMP2mod_tools', 'PMP2mod_fft5', 'PMP2mod_random', 'PMP2mod_density',
           'PMP2mod_power', 'PMP2mod_analyze', 'PMP2MG_subroutines', 'PMP2mod_MGbackground',
           'PMP2MGsolver_fR', 'PMP2extradof', 'PMP2MGsolver_DGP', 'PMP2MGsolver_sym',
           'PMP2MGsolver_kmf', 'PMP2MGsolver_csf', 'PMP2linker']
ARCHIVE_SHA = '62a4e2cd1b7363c4a029d570801252c324fb871ec1739625a64242bc8fed134c'


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        while chunk := f.read(8 * 1024**2):
            h.update(chunk)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--repo', type=Path, default=ROOT.parents[3])
    parser.add_argument('--native-dir', type=Path, help='Alternative frozen build with matching modules and objects')
    parser.add_argument('--build-receipt', type=Path, help='Receipt describing --native-dir')
    args = parser.parse_args()
    build = ROOT / 'work/build'
    build.mkdir(parents=True, exist_ok=False)
    receipt_path = args.build_receipt or args.repo / 'BDM-refine/repairs/20260906/native-build.json'
    native_receipt = json.loads(receipt_path.read_text())
    assert native_receipt['completed']
    report = dict(started_at_utc=datetime.now(timezone.utc).isoformat(), completed=False,
                  source_commit=native_receipt['source_commit'],
                  source_sha256=native_receipt['source_sha256'],
                  native_build_receipt_sha256=sha(receipt_path), builder_sha256=sha(__file__), commands=[])
    try:
        if args.native_dir:
            source_dir = args.native_dir.resolve()
            for name, expected in native_receipt['source_sha256'].items():
                assert sha(source_dir / name) == expected, name
            for path in source_dir.iterdir():
                if path.suffix == '.mod' or path.name in [name + '.o' for name in OBJECTS]:
                    shutil.copy2(path, build / path.name)
            report['native_input_directory'] = str(source_dir)
        else:
            archive = args.repo / 'BDM-refine/repairs/20260906/work-artifacts.tar.gz'
            assert sha(archive) == ARCHIVE_SHA
            report['archive_sha256'] = ARCHIVE_SHA
            with tarfile.open(archive, 'r:gz') as tar:
                for member in tar:
                    path = Path(member.name)
                    if path.parent != Path('work/native-build'):
                        continue
                    if path.suffix == '.mod' or path.name in [name + '.o' for name in OBJECTS]:
                        with tar.extractfile(member) as src, (build / path.name).open('xb') as dst:
                            shutil.copyfileobj(src, dst)
                    if path.name in native_receipt['source_sha256']:
                        assert hashlib.sha256(tar.extractfile(member).read()).hexdigest() == native_receipt['source_sha256'][path.name]
        report['linked_input_sha256'] = {p.name: sha(p) for p in sorted(build.iterdir())}
        for name in ['thread_probe.f90', 'thread_probe_entry.f90']:
            shutil.copy2(ROOT / name, build / name)
        report['probe_source_sha256'] = {p.name: sha(p) for p in build.glob('*.f90')}
        flags = ['-g', '-traceback', '-qopenmp', '-march=core-avx2', '-shared-intel',
                 '-mcmodel=medium', '-convert', 'big_endian', '-fp-model', 'precise', '-O2']
        env = dict(os.environ, LD_LIBRARY_PATH=os.environ['BDM_AUDIT_NATIVE_LIBS'])
        env.pop('LIBRARY_PATH', None)
        for command in [['ifx', '--version'],
                        ['ifx', *flags, '-c', 'thread_probe.f90'],
                        ['ifx', *flags, '-c', 'thread_probe_entry.f90'],
                        ['ifx', *flags, '-o', 'BDM-thread-probe.exe',
                         *[name + '.o' for name in OBJECTS], 'thread_probe.o', 'thread_probe_entry.o']]:
            process = subprocess.run(command, cwd=build, env=env, capture_output=True, text=True)
            report['commands'].append(dict(command=command, returncode=process.returncode,
                                           stdout=process.stdout, stderr=process.stderr))
            if process.returncode:
                raise RuntimeError(process.stdout + process.stderr)
        report.update(completed=True, binary_sha256=sha(build / 'BDM-thread-probe.exe'),
                      finished_at_utc=datetime.now(timezone.utc).isoformat())
    finally:
        (ROOT / 'build.json').write_text(json.dumps(report, indent=2) + '\n')
    print('Linked native thread diagnostic', report['binary_sha256'])


if __name__ == '__main__':
    main()
