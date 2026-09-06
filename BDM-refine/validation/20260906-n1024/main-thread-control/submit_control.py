"""Submit the sealed 64 CPU/288 GiB experiment after main simulation 11948372."""
import json
from pathlib import Path
import subprocess

import run_control as c


def main():
    root = c.ROOT
    assert not (root / 'submission.json').exists(), 'Do not duplicate a submitted experiment'
    plan_path = root / 'plan.json'
    plan = json.loads(plan_path.read_text())
    controls = json.loads((root / 'driver-controls.json').read_text())
    assert controls['completed'] and controls['driver_sha256'] == plan['runner_sha256'] == c.sha(root / 'run_control.py')
    for name, expected in plan['frozen_files_sha256'].items():
        assert c.sha(root / name) == expected
    command = ['sbatch', '--parsable', f'--dependency={plan["dependency"]}',
               f'--chdir={root}', f'--output={root}/work/slurm-%j.log', f'--error={root}/work/slurm-%j.log',
               f'--export=ALL,BDM_MAIN_THREAD_ROOT={root}', str(root / 'work/launch/control.sbatch'),
               str(root / 'work/launch/run_control.py'), '--plan-sha256', c.sha(plan_path)]
    record = dict(prepared_at_utc=c.now(), command=command, dependencies=[plan['dependency']],
        resources=plan['resources'], plan_sha256=c.sha(plan_path), build_sha256=c.sha(root / 'build.json'),
        script_sha256={p.name:c.sha(p) for p in (root / 'work/launch').iterdir()},
        driver_controls_sha256=c.sha(root / 'driver-controls.json'), submitter_sha256=c.sha(__file__),
        commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=root, text=True).strip())
    c.write_json(root / 'submission.json', record)
    process = subprocess.run(command, capture_output=True, text=True)
    record.update(submitted_at_utc=c.now(), submission_returncode=process.returncode,
                  stdout=process.stdout, stderr=process.stderr)
    if process.returncode == 0:
        record['job_id'] = process.stdout.strip().split(';')[0]
        assert record['job_id'].isdigit()
    c.write_json(root / 'submission.json', record)
    process.check_returncode()
    record['scontrol'] = subprocess.check_output(['scontrol', 'show', 'job', record['job_id']], text=True)
    c.write_json(root / 'submission.json', record)
    print(record['job_id'])


if __name__ == '__main__':
    main()
