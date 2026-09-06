"""Copy the proven diagnostic, link final objects, and seal a new experiment.

With Intel modules loaded and BDM_AUDIT_NATIVE_LIBS captured, execute using
micromamba run -n cosemu python3 -B prepare_control.py --original-repo PATH.
Never writes numerical-threading or the original repository.
"""
import argparse
import json
from pathlib import Path
import shutil
import subprocess

import run_control as c


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--original-repo', type=Path, required=True)
    args = parser.parse_args()
    root = c.ROOT
    validation = args.original_repo.resolve() / 'BDM-refine/validation/20260906-n1024'
    native = validation / 'work/final-native-build'
    receipt = validation / 'native-build.json'
    preparation = validation / 'main-preparation.json'
    original_build = json.loads(receipt.read_text())
    production = json.loads(preparation.read_text())
    assert original_build['completed']
    assert original_build['source_sha256'] == production['sources_sha256']
    assert c.sha(receipt) == production['native_build_sha256']
    assert production['sources_sha256']['PMP2linker.f90'] == '1353532b32f6f3f80457208f5d023f8f359097a35a24edfc41856453a21b154c'
    historical = root.parent / 'numerical-threading'
    pilot = json.loads((historical / 'build.json').read_text())
    frozen = root / 'work/frozen'
    frozen.mkdir(parents=True, exist_ok=False)
    for name in ['thread_probe.f90', 'thread_probe_entry.f90', 'publication_preflight.f90']:
        assert c.sha(historical / name) == pilot['probe_source_sha256'][name]
        shutil.copy2(historical / name, frozen / name)
    shutil.copy2(historical / 'build_probe.py', frozen / 'build_probe.py')
    for name in ['native-build.json', 'main-preparation.json']:
        shutil.copy2(validation / name, frozen / name)
    source_dir = frozen / 'production-sources'
    source_dir.mkdir()
    for name, expected in original_build['source_sha256'].items():
        assert c.sha(native / name) == expected
        shutil.copy2(native / name, source_dir / name)
        assert c.sha(source_dir / name) == expected
    for name, expected in original_build['binaries_sha256'].items():
        assert c.sha(native / name) == expected
    command = ['micromamba', 'run', '-n', 'cosemu', 'python3', '-B', str(frozen / 'build_probe.py'),
               '--native-dir', str(native), '--build-receipt', str(frozen / 'native-build.json')]
    subprocess.run(command, check=True)
    shutil.copy2(frozen / 'build.json', root / 'build.json')
    build = json.loads((root / 'build.json').read_text())
    assert build['completed'] and build['probe_source_sha256'] == pilot['probe_source_sha256']
    for name, expected in build['linked_input_sha256'].items():
        assert c.sha(native / name) == expected, 'Native objects changed during copying'
    config = validation / 'work/main/Run1/BDM.config'
    launch = root / 'work/launch'
    launch.mkdir()
    for name in ['run_control.py', 'control.sbatch']:
        shutil.copy2(root / name, launch / name)
    files = [root / 'build.json', *[p for p in frozen.rglob('*') if p.is_file()], *launch.iterdir()]
    resources = dict(partition='cosma8-serial', account='dp004', cpus=64, mem='288G',
        time_limit_seconds=2700, expected_seconds_range=[600, 1200],
        expected_core_hours_range=[64*600/3600, 64*1200/3600], time_limit_core_hours=48,
        exclusive=False, pilot_job='11948371', pilot_batch_maxrss_kib=31325680,
        pilot_native_maxrss_kib=20359012, volume_multiplier=8, headroom_factor=1.2,
        projected_batch_gib=31325680 / 1024**2 * 8, projected_with_headroom_gib=31325680 / 1024**2 * 8 * 1.2,
        memory_accounting='N1024 adds one immutable 32 GiB FI copy and 24 GiB saved particle bits. '
            'The active 32 GiB FI is separate and BDM frees/reallocates it. '
            'Two 32 GiB output tapes are written sequentially to disk, not held as two saved FI arrays. '
            'Measured N512 batch RSS includes transient allocations and page cache; volume scaling '
            'plus 20 percent headroom gives 286.80 GiB, rounded to 288 GiB.',
        packing='Six sequential native BDM calls use at most 64 threads. The shared request leaves '
            '64 cores schedulable, allowing the independent 64-core 192 GiB normal validation concurrently.')
    plan = dict(prepared_at_utc=c.now(), validation_root=str(validation), parent_job_id='11948372',
        dependency='afterok:11948372', source_commit=build['source_commit'],
        integration_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=root, text=True).strip(),
        prepare_sha256=c.sha(__file__), build_command=command, historical_diagnostic_build_sha256=c.sha(historical / 'build.json'),
        source_sha256=build['source_sha256'], production_binaries_sha256=production['binaries_sha256'],
        config=config.read_text(), config_sha256=c.sha(config), resources=resources,
        executable=str((frozen / 'work/build/BDM-thread-probe.exe').relative_to(root)),
        runner_sha256=c.sha(launch / 'run_control.py'),
        protocol='DENSIT once at64; immutable FI -> BDM(0) at32,64,32,64; '
            'DENSIT once at32; immutable second FI -> BDM(0) at32,64. '
            'Require six exact density and particle restorations and byte equality within each field. '
            'Between-field differences and ordinary DENSIT replays remain explicit sensitivity findings.',
        frozen_files_sha256={str(p.relative_to(root)): c.sha(p) for p in sorted(files)})
    c.write_json(root / 'plan.json', plan)
    print('Prepared main control', c.sha(root / 'plan.json'), 'binary', build['binary_sha256'])


if __name__ == '__main__':
    main()
