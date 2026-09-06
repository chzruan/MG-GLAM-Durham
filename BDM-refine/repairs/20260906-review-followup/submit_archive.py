"""Archive completed scratch on a measured shared allocation.

Run with micromamba run -n cosemu python3 -B. Existing receipts are preserved.
"""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import subprocess


ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[2]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    receipt = ROOT / 'cleanup-submission.json'
    assert not receipt.exists()
    comparison = ROOT / 'native/comparison.json'
    accounting = ROOT / 'native/accounting.json'
    assert json.loads(comparison.read_text())['completed']
    accounts = json.loads(accounting.read_text())
    assert accounts['all_completed']
    jobs = list(accounts['jobs'])
    pilot_command = ['sacct', '-j', '11948648', '--parsable2',
                     '--format=JobID,State,ElapsedRaw,TotalCPU,ReqTRES%100,AllocTRES%100,MaxRSS,ExitCode']
    pilot = subprocess.check_output(pilot_command, text=True)
    assert '11948648|COMPLETED|' in pilot
    launch = ROOT / 'work/launch-archive'
    launch.mkdir(exist_ok=False)
    shutil.copy2(ROOT / 'archive_work.py', launch / 'archive_work.py')
    batch = launch / 'archive.sbatch'
    batch.write_text('''#!/bin/bash
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
ulimit -c 0
micromamba run -n cosemu python3 -B "$1"
''')
    command = ['sbatch', '--parsable', '--partition=cosma8-serial', '--account=dp004',
               '--nodes=1', '--ntasks=1', '--cpus-per-task=1', '--mem=4G',
               '--time=00:10:00', '--job-name=bdm-v3-archive',
               '--dependency=afterok:' + ':'.join(jobs), '--chdir=' + str(REPO),
               f'--output={ROOT}/work/archive-%j.log', f'--error={ROOT}/work/archive-%j.log',
               '--export=ALL,BDM_REVIEW_ROOT=' + str(ROOT),
               str(batch), str(launch / 'archive_work.py')]
    record = dict(prepared_at_utc=datetime.now(timezone.utc).isoformat(), command=command,
                  dependencies=jobs, submitter_sha256=sha(__file__),
                  script_sha256={p.name: sha(p) for p in launch.iterdir()},
                  comparison_sha256=sha(comparison), accounting_sha256=sha(accounting),
                  pilot_command=pilot_command, pilot_accounting=pilot,
                  resources=dict(partition='cosma8-serial', cpus=1, mem='4G', exclusive=False,
                      expected_seconds_range=[60, 240], expected_core_hours_range=[1/60, 4/60],
                      time_limit_seconds=600, time_limit_core_hours=1/6,
                      reason='Earlier verified archive job 11948648 used 87 s, 66.866 CPU s, '
                             'and 1386244 KiB batch MaxRSS on one shared core. Four GiB gives '
                             'headroom for streaming compression/readback of the follow-up '
                             'build and logs. Large science tapes remain unpacked.'))
    result = subprocess.run(command, capture_output=True, text=True)
    record.update(submitted_at_utc=datetime.now(timezone.utc).isoformat(),
                  returncode=result.returncode, stdout=result.stdout, stderr=result.stderr)
    if result.returncode == 0:
        record['job_id'] = result.stdout.strip().split(';')[0]
        record['scontrol'] = subprocess.check_output(
            ['scontrol', 'show', 'job', record['job_id'], '-o'], text=True)
    receipt.write_text(json.dumps(record, indent=2) + '\n')
    assert result.returncode == 0, result.stderr
    print('Submitted verified cleanup', record['job_id'])


if __name__ == '__main__':
    main()
