"""Freeze and submit the same-snapshot physics/property replays after evolution.

Run with micromamba run -n cosemu python3 -B. The separate controlled-density
job supplies the finder-thread conclusion; ordinary density replay differences
remain explicit in main-validation.json.
"""
import json
from pathlib import Path
import shutil
import subprocess

import run_validation as v


def main():
    jobs = json.loads((v.ROOT / 'jobs.json').read_text())
    assert not any(j.get('stage') == 'main-validation' for j in jobs), 'Already submitted'
    parent = next(j for j in jobs if j['job_id'] == '11948372')
    controls = json.loads((v.ROOT / 'main-driver-controls.json').read_text())
    assert controls['completed'] and controls['driver_sha256'] == v.sha(v.ROOT / 'run_validation.py')
    preparation = v.ROOT / 'main-preparation.json'
    assert v.sha(preparation) == parent['preparation_sha256']
    native = json.loads(preparation.read_text())
    for name, value in native['sources_sha256'].items():
        assert v.sha(v.REPO / name) == value
    plan = dict(partition='cosma8-serial', account='dp004', cpus=64, mem='192G',
                time_limit_seconds=2700, expected_seconds_range=[600, 1500],
                expected_core_hours_range=[64*600/3600, 64*1500/3600],
                time_limit_core_hours=48, full_node_exclusivity_required=False,
                dependencies=['afterok:11948372'], pilot_job='11948326',
                pilot_batch_maxrss_kib=19105956,
                memory_reason='N512 normal membership replay uses 12389104 KiB native RSS. '
                    'Eightfold volume gives 94.5 GiB; the repeated-state probe adds 24 GiB '
                    'of saved particle bits at N1024. The 192 GiB request also covers '
                    'transient arrays and batch page cache. Fixed-density diagnostics '
                    'need a separate measured 288 GiB allocation.',
                packing='Nine native replays run sequentially with at most 64 threads; '
                    'streaming membership checks use one CPU. Shared allocation leaves '
                    'the other 64 cores and remaining node memory schedulable.',
                scope='Membership identity and host geometry at z2,1,0; old/new catalogues '
                    'on identical snapshots; ordinary 64/32-thread comparison and two-call '
                    'bitwise particle-state restoration at z0. Normal density differences '
                    'are recorded and interpreted alongside the separate fixed-field control.')
    v.write_json(v.ROOT / 'main-validation-resource-plan.json', plan)
    frozen = v.WORK / 'launch-main-validation'
    frozen.mkdir(exist_ok=False)
    for name in ['run_validation.py', 'validation.sbatch']:
        shutil.copy2(v.ROOT / name, frozen / name)
    command = ['sbatch', '--parsable', '--job-name=bdm-n1024-validate',
               '--cpus-per-task=64', '--mem=192G', '--time=00:45:00',
               '--dependency=afterok:11948372',
               f'--export=ALL,BDM_VALIDATION_ROOT={v.ROOT},BDM_VALIDATION_PREPARATION={preparation},'
               f'BDM_VALIDATION_BIN={native["binary_directory"]}',
               str(frozen / 'validation.sbatch'), str(frozen / 'run_validation.py'),
               'validate', '--case', 'main']
    record = dict(stage='main-validation', purpose=plan['scope'], command=command,
                  prepared_at_utc=v.now(), commit=subprocess.check_output(
                      ['git', 'rev-parse', 'HEAD'], cwd=v.REPO, text=True).strip(),
                  script_sha256={p.name: v.sha(p) for p in frozen.iterdir()},
                  submitter_sha256=v.sha(__file__), dependencies=plan['dependencies'],
                  preparation_sha256=v.sha(preparation),
                  resource_plan_sha256=v.sha(v.ROOT / 'main-validation-resource-plan.json'))
    # Preserve the exact intended request even if submission fails.
    v.write_json(frozen / 'submission.json', record)
    process = subprocess.run(command, cwd=v.REPO, text=True, capture_output=True)
    record.update(submission_returncode=process.returncode, stdout=process.stdout,
                  stderr=process.stderr, submitted_at_utc=v.now())
    if process.returncode == 0:
        record['job_id'] = process.stdout.strip().split(';')[0]
        assert record['job_id'].isdigit()
        jobs.append(record)
        v.write_json(v.ROOT / 'jobs.json', jobs)
        record['scontrol'] = subprocess.check_output(
            ['scontrol', 'show', 'job', record['job_id']], text=True)
    v.write_json(frozen / 'submission.json', record)
    process.check_returncode()
    print(record['job_id'])


if __name__ == '__main__':
    main()
