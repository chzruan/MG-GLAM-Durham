"""Atomically consolidate the two finished probe work trees; never copy tapes.

Use micromamba run -n cosemu python3 -B. Existing historical receipts retain
their original paths. The new receipt maps these paths to preserved inodes.
"""
import csv
import io
import json
import os
from pathlib import Path
import stat
import subprocess
import traceback

import run_control as c


def metadata(path):
    value = path.lstat()
    kind = 'directory' if stat.S_ISDIR(value.st_mode) else 'symlink' if stat.S_ISLNK(value.st_mode) else 'file'
    return [value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns, kind,
            os.readlink(path) if kind == 'symlink' else None]


def tree(path):
    entries = {'.': metadata(path)}
    for directory, dirs, files in os.walk(path, followlinks=False):
        for name in sorted(dirs + files):
            item = Path(directory) / name
            entries[str(item.relative_to(path))] = metadata(item)
    return dict(sorted(entries.items()))


def main():
    validation = c.ROOT.parent
    destination = Path('/cosma8/data/dp203/dc-ruan1/mgglam_claude/MG-GLAM/BDM-refine/validation/20260906-n1024')
    receipt = c.ROOT / 'relocation.json'
    assert not receipt.exists(), 'Inspect any existing relocation receipt before resuming'
    command = ['sacct', '-j', '11948363,11948371,11948491', '--parsable2', '--format=JobID,State,ExitCode']
    accounting = subprocess.check_output(command, text=True)
    jobs = {row['JobID']:row for row in csv.DictReader(io.StringIO(accounting), delimiter='|')}
    for job, state in [('11948363', 'FAILED'), ('11948371', 'COMPLETED'), ('11948491', 'COMPLETED')]:
        assert jobs[job]['State'] == state
    references = {}
    known = {}

    def manifest(name):
        path = validation / name
        references[name] = dict(original=str(path), destination=str(destination / name), sha256=c.sha(path))
        return json.loads(path.read_text())

    def link(path, digest, reference, field):
        key = str(path)
        known.setdefault(key, []).append(dict(sha256=digest, receipt=reference, field=field))

    for part in ['numerical-threading', 'main-thread-control']:
        name = part+'/results.json'
        result = manifest(name)
        assert result['completed'] and result['fixed_density_controls_passed']
        for filename, digest in result['outputs_sha256'].items():
            link(validation / part / 'work/run' / filename, digest, name, 'outputs_sha256.'+filename)
        build_name = part+'/build.json'
        build = manifest(build_name)
        prefix = 'work/build' if part == 'numerical-threading' else 'work/frozen/work/build'
        for field in ['linked_input_sha256', 'probe_source_sha256']:
            for filename, digest in build[field].items():
                link(validation / part / prefix / filename, digest, build_name, field+'.'+filename)
        link(validation / part / prefix / 'BDM-thread-probe.exe', build['binary_sha256'], build_name, 'binary_sha256')
    failed_name = 'numerical-threading/attempt-11948363.json'
    failed = manifest(failed_name)
    for field in ['linked_input_sha256', 'probe_source_sha256']:
        for filename, digest in failed['build'][field].items():
            link(validation / 'numerical-threading/work/build-11948363' / filename,
                 digest, failed_name, 'build.'+field+'.'+filename)
    link(validation / 'numerical-threading/work/build-11948363/BDM-thread-probe.exe',
         failed['build']['binary_sha256'], failed_name, 'build.binary_sha256')
    link(validation / 'numerical-threading/work/run-11948363/probe.log',
         failed['results']['log_sha256'], failed_name, 'results.log_sha256')
    plan_name = 'main-thread-control/plan.json'
    plan = manifest(plan_name)
    for filename, digest in plan['frozen_files_sha256'].items():
        link(c.ROOT / filename, digest, plan_name, 'frozen_files_sha256.'+filename)
    report = dict(started_at_utc=c.now(), completed=False, implementation_sha256=c.sha(__file__),
        command=['micromamba', 'run', '-n', 'cosemu', 'python3', '-B', str(Path(__file__).resolve())],
        method='os.rename on the same filesystem; no file contents or historical receipts are rewritten',
        metadata_columns=['device', 'inode', 'size_bytes', 'mtime_ns', 'type', 'symlink_target'],
        sha_scope='SHA references reuse validated build/execution manifests; large tapes are not rehashed on the login node',
        sacct_command=command, sacct=accounting, existing_validated_manifests=references, moves=[])
    for part in ['numerical-threading', 'main-thread-control']:
        original, target = validation / part / 'work', destination / part / 'work'
        assert original.is_dir() and not original.is_symlink()
        assert target.parent.is_dir() and not os.path.lexists(target)
        assert original.stat().st_dev == target.parent.stat().st_dev, 'Different filesystem; no copy attempted'
        before = tree(original)
        for relative, value in before.items():
            if value[4] == 'symlink':
                assert (original / relative).exists(), 'Broken input symlink before relocation'
                assert not str((original / relative).resolve()).startswith(str(validation)), 'Moving tree contains a local absolute symlink'
        report['moves'].append(dict(original=str(original), destination=str(target), renamed=False,
                                   before=before))
    c.write_json(receipt, report)
    try:
        for move in report['moves']:
            original, target = Path(move['original']), Path(move['destination'])
            os.rename(original, target)
            move.update(renamed=True, renamed_at_utc=c.now())
            after = tree(target)
            assert not original.exists() and move['before'] == after, 'Relocated metadata changed'
            move['after'] = after
            move['all_entry_metadata_equal'] = True
            move['existing_sha_links'] = {relative:known[str(original / relative)] for relative in after
                                         if str(original / relative) in known}
            for relative, value in after.items():
                if value[4] == 'symlink':
                    assert (target / relative).exists(), 'Relocated input symlink is broken'
            move['all_symlink_targets_still_exist'] = True
            c.write_json(receipt, report)
        report.update(completed=True, finished_at_utc=c.now())
    except BaseException:
        report['failure_traceback'] = traceback.format_exc()
        raise
    finally:
        c.write_json(receipt, report)
    print('Atomic relocation complete:', [(m['destination'], len(m['after'])) for m in report['moves']])


if __name__ == '__main__':
    main()
