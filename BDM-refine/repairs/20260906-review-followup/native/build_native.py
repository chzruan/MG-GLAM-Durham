"""Freeze native replay variants; run only small one-core build preflights.

Invoke with micromamba run -n cosemu python3 -B after loading Intel modules.
BDM_AUDIT_NATIVE_LIBS must contain the module runtime path from outside cosemu.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import struct
import subprocess

HERE = Path(__file__).resolve().parent
REFERENCE = '800eaac76ef7fe78e6a392c85c3a88822b1c9135'
REFERENCE_SHA = 'ef8a0bdf0c3406c83404ee35f569bd89705bcbb5670b220ffa63af156c4f7598'
OBJECTS = ['PMP2mod_tools', 'PMP2mod_fft5', 'PMP2mod_random', 'PMP2mod_density',
           'PMP2mod_power', 'PMP2mod_analyze', 'PMP2MG_subroutines', 'PMP2mod_MGbackground',
           'PMP2MGsolver_fR', 'PMP2extradof', 'PMP2MGsolver_DGP', 'PMP2MGsolver_sym',
           'PMP2MGsolver_kmf', 'PMP2MGsolver_csf']
FLAGS = ['-O3', '-g', '-traceback', '-ftz', '-unroll', '-qopenmp', '-march=core-avx2',
         '-mfma', '-shared-intel', '-mcmodel=medium', '-convert', 'big_endian', '-fp-model', 'precise']


def now():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        while data := stream.read(8 * 1024**2):
            h.update(data)
    return h.hexdigest()


def write_json(path, value):
    temporary = path.with_suffix(path.suffix + '.part')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    os.replace(temporary, path)


def replace_once(source, old, new):
    assert source.count(old) == 1, old
    return source.replace(old, new)


def instrument(source, diagnostic):
    source = replace_once(source, 'integer*8, allocatable :: HaloUnbindingWork(:)',
                          'integer*8, allocatable :: HaloUnbindingWork(:)\n'
                          'logical :: RepairDumpFinalPass=.false. ! diagnostic build only')
    assert source.count('Call WriteFiles\n') == 2
    source = source.replace('Call WriteFiles\n',
                            'Call WriteFiles\n      if(RepairDumpFinalPass)call RepairDiagnosticsDump\n')
    matches = list(re.finditer(r'^end\s+Module\s+LinkerList\s*$', source, re.I | re.M))
    assert len(matches) == 1
    match = matches[0]
    return source[:match.start()] + diagnostic + '\n' + source[match.start():]


def fixture(path):
    """One small PM page with 64 real rows and sparse unused record padding."""
    path.mkdir(parents=True)
    (path/'CATALOGS').mkdir()
    (path/'BDM.config').write_text('iVirial=2\nMassMin=0\nRext=0.\n')
    header = b'Native follow-up preflight'.ljust(45) + struct.pack(
        '>4fif6f4i3fq100f', 1., .02, 1., .004, 1, 1., *([0.] * 6),
        4, 16, 1, 1, .3, .7, .7, 64, *([0.] * 99 + [32.]))
    record = struct.pack('>i', len(header))
    (path/'PMcrd.0001.DAT').write_bytes(record + header + record)
    coordinates = [[4. + .04 * (index - 1.5) for k in range(4) for j in range(4) for i in range(4)
                    for index in [[i, j, k][axis]]] for axis in range(3)]
    velocities = [[-0. if j % 2 else 0. for j in range(64)] for _ in range(3)]
    with (path/'PMcrs0.0001.DAT').open('xb') as stream:
        stream.truncate(6 * 1024**2 * 4)
        for axis, values in enumerate(coordinates + velocities):
            stream.seek(axis * 1024**2 * 4)
            stream.write(struct.pack('>64f', *values))


def link_inputs(path, target):
    target.mkdir(parents=True)
    (target/'CATALOGS').mkdir()
    shutil.copy2(path/'BDM.config', target/'BDM.config')
    for item in path.glob('PMcr*.DAT'):
        (target/item.name).symlink_to(item)


def check_tapes(path, expected_selected):
    raw = (path/'repair-members.bin').read_bytes()
    selected, candidates, mass_one = struct.unpack_from('>qqf', raw)
    assert selected == expected_selected and mass_one > 0
    pos = 20
    counts = {}
    for _ in range(selected):
        candidate, count = struct.unpack_from('>qq', raw, pos)
        properties = struct.unpack_from('>21f', raw, pos + 16)
        ids = struct.unpack_from(f'>{count}q', raw, pos + 100)
        assert count == 64 and ids == tuple(range(1, 65))
        assert abs(properties[6] - count * mass_one) <= 2.e-7 * properties[6]
        counts[candidate] = count
        pos += 100 + 8 * count
    assert pos == len(raw)
    telemetry = (path/'unbinding.bin').read_bytes()
    assert struct.unpack_from('>q', telemetry)[0] == candidates
    assert len(telemetry) == 8 + 28 * candidates
    for i in range(candidates):
        passes, status, work, count, mass = struct.unpack_from('>iiqqf', telemetry, 8 + 28 * i)
        assert passes >= 0 and status >= 0 and work >= 0 and count >= 0
        if i + 1 in counts:
            assert (passes, work, count) == (1, 64, 64)
            assert mass > 0
    return dict(selected=selected, candidates=candidates, membership_bytes=len(raw),
                unbinding_bytes=len(telemetry), exact_original_ids=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--repo', type=Path, default=HERE.parents[3])
    parser.add_argument('--root', type=Path, default=HERE.parent,
                        help='Follow-up R2 directory; outputs go to native/build.json and work/native-build')
    args = parser.parse_args()
    repo, root = args.repo.resolve(), args.root.resolve()
    receipt = root/'native/build.json'
    assert not receipt.exists(), 'Preserve the previous build receipt; use a new --root'
    build = root/'work/native-build'
    build.mkdir(parents=True, exist_ok=False)
    receipt.parent.mkdir(parents=True, exist_ok=True)
    production = build/'production'
    production.mkdir()
    sources = sorted(repo.glob('*.f90')) + sorted(repo.glob('*.h')) + [repo/'makefile']
    assert len(sources) == 43, 'Review source manifest changes explicitly'
    report = dict(completed=False, started_at_utc=now(), builder_sha256=sha(__file__),
                  source_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=repo, text=True).strip(),
                  source_worktree_status=subprocess.check_output(['git', 'status', '--porcelain', '--untracked-files=no'],
                                                                 cwd=repo, text=True),
                  source_repository=str(repo), frozen_build_path=str(build), frozen_source_path=str(production),
                  reference_commit=REFERENCE, commands=[], variants={}, preflight=[])
    env = dict(os.environ, LD_LIBRARY_PATH=os.environ['BDM_AUDIT_NATIVE_LIBS'], OMP_NUM_THREADS='1',
               OMP_DYNAMIC='FALSE', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    env.pop('LIBRARY_PATH', None)
    for variable in ['MAKEFLAGS', 'MFLAGS', 'MAKEOVERRIDES']:
        env.pop(variable, None)

    def run(command, cwd, stdin=None, success=True, timeout=300):
        process = subprocess.run([str(x) for x in command], cwd=cwd, env=env, input=stdin,
                                 text=True, capture_output=True, timeout=timeout)
        entry = dict(command=[str(x) for x in command], cwd=str(cwd), stdin=stdin,
                     returncode=process.returncode, stdout=process.stdout, stderr=process.stderr)
        report['commands'].append(entry)
        write_json(receipt, report)
        if success and process.returncode:
            raise RuntimeError(process.stdout + '\n' + process.stderr)
        if not success and process.returncode == 0:
            raise AssertionError('Expected explicit failure: ' + repr(command))
        return process

    try:
        for path in sources:
            shutil.copy2(path, production/path.name)
        report['source_sha256'] = {p.name: sha(production/p.name) for p in sources}
        assert all(sha(p) == report['source_sha256'][p.name] for p in sources)
        diagnostic_sources = ['replay_probe.f90', 'replay_entry.f90', 'publication_preflight.f90', 'diagnostics.f90']
        for name in diagnostic_sources:
            shutil.copy2(HERE/name, build/name)
        report['diagnostic_sha256'] = {name: sha(build/name) for name in diagnostic_sources}
        reference = subprocess.check_output(['git', 'show', REFERENCE + ':PMP2linker.f90'], cwd=repo)
        assert hashlib.sha256(reference).hexdigest() == REFERENCE_SHA
        # All variants deliberately share unchanged simulation/runtime objects.
        for name in OBJECTS:
            original = subprocess.check_output(['git', 'show', REFERENCE + ':' + name + '.f90'], cwd=repo)
            assert hashlib.sha256(original).hexdigest() == report['source_sha256'][name + '.f90'], name
        run(['ifx', '--version'], production)
        run(['make', '-f', 'makefile', '-j1', 'PMP2mod_tools.o', 'PMP2mod_MGbackground.o'], production)
        run(['make', '-f', 'makefile', '-j1', 'PMP2BDM'], production)
        report['production_binary'] = dict(binary_path=str(production/'PMP2BDM.exe'),
                                           binary_sha256=sha(production/'PMP2BDM.exe'))
        common = [production/(name + '.o') for name in OBJECTS]
        report['common_objects_sha256'] = {p.name: sha(p) for p in common}
        report['common_modules_sha256'] = {p.name: sha(p) for p in sorted(production.glob('*.mod'))}
        final = (production/'PMP2linker.f90').read_text()
        optimized_v2 = replace_once(final,
            'threshold=sphere_volume*dble(Ovdens)*mass*(dble(NROW)/dble(Box))**3',
            'threshold=1.150d12*dble(Om0)*dble(Ovdens)')
        optimized_v2 = replace_once(optimized_v2, '[BDM finder v3]', '[BDM finder v2]')
        diagnostic = (build/'diagnostics.f90').read_text()
        for variant, source in [('reference', reference.decode()), ('optimized-v2', optimized_v2), ('v3', final)]:
            directory = build/variant
            directory.mkdir()
            for module in production.glob('*.mod'):
                if module.name.lower() not in ['structures.mod', 'linkerlist.mod', 'bdmduplicaterules.mod']:
                    shutil.copy2(module, directory/module.name)
            (directory/'finder-original.f90').write_text(source)
            (directory/'finder.f90').write_text(instrument(source, diagnostic))
            for name in diagnostic_sources[:3]:
                shutil.copy2(build/name, directory/name)
            run(['ifx', *FLAGS, '-c', 'finder.f90', '-o', 'finder.o'], directory)
            run(['ifx', *FLAGS, '-c', 'replay_probe.f90', 'replay_entry.f90'], directory)
            binary = directory/'replay.exe'
            run(['ifx', *FLAGS, '-o', binary, *common, 'finder.o', 'replay_probe.o', 'replay_entry.o'], directory)
            run(['ifx', *FLAGS, '-o', 'publication-preflight.exe', *common,
                 'finder.o', 'replay_probe.o', 'publication_preflight.f90'], directory)
            report['variants'][variant] = dict(binary_path=str(binary), binary_sha256=sha(binary),
                source_sha256=sha(directory/'finder.f90'), finder_source_path=str(directory/'finder.f90'),
                uninstrumented_source_sha256=sha(directory/'finder-original.f90'),
                uninstrumented_source_path=str(directory/'finder-original.f90'), flags=FLAGS,
                publication_binary_sha256=sha(directory/'publication-preflight.exe'),
                diagnostic_calls_after_successful_writefiles=2, diagnostics_only_on_final_pass=True)

        preflight = build/'preflight'
        fixture(preflight/'snapshot')
        snapshots = preflight/'snapshot'
        snapshot_hashes = {p.name: sha(p) for p in snapshots.glob('PMcr*.DAT')}
        # Use a >120-character absolute density path to exercise real production paths.
        tape_dir = preflight/('density-path-' + 'x' * 140)
        tape_dir.mkdir()
        density = tape_dir/'fixed.bin'
        empty = tape_dir/'empty.bin'
        empty.write_bytes(struct.pack('>qqi', 16, 64, 1) + bytes(4 * 16**3))
        assert len(str(density)) > 120
        for variant, details in report['variants'].items():
            for case in ['bound', 'empty']:
                directory = preflight/f'{variant}-{case}'
                link_inputs(snapshots, directory)
                mode = 'write' if variant == 'reference' and case == 'bound' else 'read'
                tape = empty if case == 'empty' else density
                before = sha(tape) if tape.exists() else None
                result = run([details['binary_path'], 1, 1, 2, tape, mode], directory)
                for marker in ['REPLAY FIXED FI pass=', 'REPLAY FINDER pass=', 'REPLAY RESTORED pass=']:
                    assert result.stdout.count(marker) == 2, marker
                assert result.stdout.count('REPLAY DIAGNOSTICS candidates=') == 1
                assert result.stdout.count('REPLAY COMPLETE') == 1
                assert (directory/'p1.DAT').read_bytes() == (directory/'p2.DAT').read_bytes()
                assert len((directory/'p1.DAT').read_text().splitlines()) == (9 if case == 'bound' else 8)
                if before is not None:
                    assert sha(tape) == before, 'Read replay changed its density tape'
                assert not list(directory.glob('*.part'))
                checked = check_tapes(directory, int(case == 'bound'))
                report['preflight'].append(dict(variant=variant, case=case, completed=True,
                    repeated_catalogue_bytes_equal=True, restoration_markers=2, **checked))
            pub = build/variant/'publication-preflight.exe'
            for case in ['valid', 'empty', 'invalid']:
                directory = preflight/f'{variant}-publication-{case}'
                link_inputs(snapshots, directory)
                prior = directory/'CATALOGS/CatshortM.0001.0001.DAT'
                prior.write_bytes(b'prior complete catalogue\n')
                if case == 'invalid':
                    run([pub, 'valid'], directory)
                    (directory/'preserved.DAT').replace(prior)
                    original = prior.read_bytes()
                    assert len(original.splitlines()) == 9
                result = run([pub, case], directory, success=case != 'invalid')
                if case == 'invalid':
                    assert prior.read_bytes() == original
                    assert not (directory/'preserved.DAT').exists()
                else:
                    assert 'NATIVE PUBLICATION PREFLIGHT COMPLETE' in result.stdout
                    assert len((directory/'preserved.DAT').read_text().splitlines()) == (8 if case == 'empty' else 9)
                report['preflight'].append(dict(variant=variant, case='publication-' + case,
                    prior_preserved_on_failure=case == 'invalid', completed=True))
        pub = build/'v3/publication-preflight.exe'
        for case in ['leak-passes', 'leak-work']:
            directory = preflight/case
            directory.mkdir()
            result = run([pub, case], directory, success=False)
            assert 'Replay found retained membership or telemetry workspace' in result.stderr
            report['preflight'].append(dict(variant='v3', case=case, injected_workspace_leak_rejected=True))
        for case in ['density-truncated', 'density-metadata', 'particle-limit', 'header-step']:
            directory = preflight/case
            link_inputs(snapshots, directory)
            tape = directory/'invalid-density.bin'
            if case.startswith('density-'):
                contents = empty.read_bytes()
                if case == 'density-truncated':
                    contents = contents[:-4]
                else:
                    contents = struct.pack('>q', 15) + contents[8:]
                tape.write_bytes(contents)
                expected = 'Replay density tape byte length mismatch' if case.endswith('truncated') else 'Replay density metadata mismatch'
            else:
                header = bytearray((directory/'PMcrd.0001.DAT').read_bytes())
                if case == 'particle-limit':
                    struct.pack_into('>i', header, 4+45+12*4, 1200)
                    struct.pack_into('>q', header, 4+45+19*4, 1200**3)
                    expected = 'Replay requires NROW**3 < 1200**3'
                else:
                    struct.pack_into('>i', header, 4+45+4*4, 2)
                    expected = 'Replay snapshot/header step mismatch'
                (directory/'PMcrd.0001.DAT').unlink()
                (directory/'PMcrd.0001.DAT').write_bytes(header)
            result = run([report['variants']['v3']['binary_path'], 1, 1, 2, tape, 'read'], directory, success=False)
            assert expected in result.stderr
            assert 'REPLAY COMPLETE' not in result.stdout and not list(directory.glob('p*.DAT'))
            report['preflight'].append(dict(variant='v3', case=case, invalid_input_rejected=True))
        assert (preflight/'reference-bound/p2.DAT').read_bytes() == (preflight/'optimized-v2-bound/p2.DAT').read_bytes()
        for name in ['repair-members.bin', 'unbinding.bin']:
            assert (preflight/'reference-bound'/name).read_bytes() == (preflight/'optimized-v2-bound'/name).read_bytes()
        report['optimized_v2_preflight_catalogue_membership_telemetry_equal_reference'] = True
        directory = preflight/'production'
        link_inputs(snapshots, directory)
        result = run([production/'PMP2BDM.exe'], directory, stdin='1\n')
        assert (directory/'CATALOGS/CatshortM.0001.0001.DAT').read_bytes() == (preflight/'v3-bound/p2.DAT').read_bytes()
        report['production_make_binary_preflight_matches_instrumented_v3'] = True
        assert all(sha(snapshots/name) == expected for name, expected in snapshot_hashes.items())
        assert all(sha(p) == report['source_sha256'][p.name] for p in sources), 'Live source changed during build'
        assert all(sha(p) == report['common_objects_sha256'][p.name] for p in common)
        report.update(completed=True, finished_at_utc=now(), snapshot_preflight_sha256=snapshot_hashes,
                      preflight_density_sha256=sha(density), preflight_long_density_path=str(density))
    except BaseException as error:
        report['failure'] = repr(error)
        raise
    finally:
        write_json(receipt, report)
    print('Frozen native follow-up build and one-core preflights passed:', receipt)


if __name__ == '__main__':
    main()
