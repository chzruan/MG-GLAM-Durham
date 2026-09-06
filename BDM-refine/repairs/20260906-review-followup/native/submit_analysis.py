"""Queue one-core analysis after both replays; use micromamba run -n cosemu python3 -B."""
import json
from pathlib import Path
import shutil
import subprocess

import run_replays as r


def main():
    receipt=r.HERE/'analysis-submission.json'
    assert not receipt.exists()
    parents=json.loads((r.HERE/'submission.json').read_text())
    jobs=[s['job_id'] for s in parents]
    assert len(jobs)==2 and all(s['returncode']==0 for s in parents)
    pilot=json.loads((r.HERE/'analysis-pilot.json').read_text())
    assert pilot['completed'] and pilot['driver_sha256']==r.sha(r.HERE/'analyze_replays.py')
    launch=r.ROOT/'work/launch-analysis'
    launch.mkdir(exist_ok=False)
    for name in ['analyze_replays.py','run_replays.py']:
        shutil.copy2(r.HERE/name,launch/name)
    batch=launch/'analysis.sbatch'
    batch.write_text('''#!/bin/bash
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
ulimit -c 0
micromamba run -n cosemu python3 -B "$1"
''')
    command=['sbatch','--parsable','--partition=cosma8-serial','--account=dp004','--nodes=1','--ntasks=1',
             '--cpus-per-task=1','--mem=8G','--time=00:15:00','--job-name=bdm-v3-comparison',
             '--dependency=afterok:'+':'.join(jobs),'--chdir='+str(r.REPO),
             f'--output={r.ROOT}/work/analysis-%j.log',f'--error={r.ROOT}/work/analysis-%j.log',
             '--export=ALL,BDM_REVIEW_ROOT='+str(r.ROOT),str(batch),str(launch/'analyze_replays.py')]
    record=dict(command=command,dependencies=jobs,prepared_at_utc=r.now(),
                submitter_sha256=r.sha(__file__),script_sha256={p.name:r.sha(p) for p in launch.iterdir()},
                plan_sha256=r.sha(r.HERE/'plan.json'),pilot_sha256=r.sha(r.HERE/'analysis-pilot.json'),
                resources=dict(partition='cosma8-serial',account='dp004',cpus=1,mem='8G',time_limit_seconds=900,
                    expected_seconds_range=[120,360],expected_core_hours_range=[120/3600,360/3600],time_limit_core_hours=.25,
                    pilot='Full-size host-control job11949318 used 455556 KiB process RSS and 4181824 KiB Slurm MaxRSS. '
                      'The paired-index/NPZ analysis uses bounded per-halo membership arrays; 8 GiB doubles the measured '
                      'batch memory for extra catalogue arrays and I/O. Only one CPU is requested for its serial work.'))
    result=subprocess.run(command,capture_output=True,text=True)
    record.update(submitted_at_utc=r.now(),returncode=result.returncode,stdout=result.stdout,stderr=result.stderr)
    if result.returncode==0:
        record['job_id']=result.stdout.strip().split(';')[0]
        record['scontrol']=subprocess.check_output(['scontrol','show','job',record['job_id'],'-o'],text=True)
    r.write_json(receipt,record)
    assert result.returncode==0,result.stderr
    print('Submitted comparison',record['job_id'],'after',','.join(jobs))


if __name__=='__main__':
    main()
