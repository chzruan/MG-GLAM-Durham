"""Build the campaign-only matched native IC executable in an isolated directory."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
from generate import generate


def build(repo: Path, work: Path, compiler: str, checked: bool = False) -> dict:
    work.mkdir(parents=True, exist_ok=True)
    if compiler == 'gfortran':
        flags = ['-O2', '-fopenmp', '-fno-fast-math', '-ffp-contract=off', '-std=legacy',
                 '-fallow-argument-mismatch', '-fconvert=big-endian', '-ffree-line-length-none']
        if checked:
            flags += ['-g', '-fcheck=all', '-fbacktrace']
    elif compiler == 'ifx':
        flags = ['-O2', '-qopenmp', '-fp-model', 'precise', '-convert', 'big_endian', '-assume', 'byterecl']
        if checked:
            flags += ['-g', '-check', 'bounds', '-traceback']
    else:
        raise ValueError('Supported compilers: ifx, gfortran')
    generated = generate(repo)
    (work / 'PMP2start.matched.f90').write_text(generated)
    sources = ['PMP2mod_tools.f90', 'PMP2mod_random.f90', 'PMP2mod_fft5.f90']
    for name in sources + ['luxuryp.h']:
        shutil.copy2(repo / name, work / name)
    commands = []
    with (work / 'build.log').open('w') as log:
        version = subprocess.run([compiler, '--version'], check=True, capture_output=True, text=True).stdout
        log.write(version)
        for name in sources + ['PMP2start.matched.f90']:
            command = [compiler, *flags, '-c', name]
            commands.append(command)
            subprocess.run(command, cwd=work, stdout=log, stderr=subprocess.STDOUT, check=True)
        command = [compiler, *flags, '-o', 'PMP2start.matched.exe',
                   *[Path(name).with_suffix('.o').name for name in sources + ['PMP2start.matched.f90']]]
        commands.append(command)
        subprocess.run(command, cwd=work, stdout=log, stderr=subprocess.STDOUT, check=True)
    receipt = dict(compiler=compiler, compiler_version=version, flags=flags, checked=checked, commands=commands,
                   binary_sha256=hashlib.sha256((work/'PMP2start.matched.exe').read_bytes()).hexdigest(),
                   generated_sha256=hashlib.sha256(generated.encode()).hexdigest(),
                   source_sha256={name: hashlib.sha256((repo/name).read_bytes()).hexdigest()
                                  for name in sources + ['PMP2start.f90', 'luxuryp.h']})
    (work / 'build.json').write_text(json.dumps(receipt, indent=2) + '\n')
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path, required=True)
    parser.add_argument('--work', type=Path, required=True)
    parser.add_argument('--compiler', choices=['ifx', 'gfortran'], default='ifx')
    parser.add_argument('--checked', action='store_true')
    args = parser.parse_args()
    print(json.dumps(build(args.repo.resolve(), args.work.resolve(), args.compiler, args.checked), indent=2))


if __name__ == '__main__':
    main()
