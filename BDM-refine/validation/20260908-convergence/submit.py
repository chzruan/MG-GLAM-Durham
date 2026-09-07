"""Submit one explicitly sized stage; freeze Slurm script and billing estimate."""
import argparse
import json
import os
import shlex
import subprocess
import zipfile
from common import REPO, ROOT, WORK, git, now, sha, write_json


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('phase',choices=['ic','evolve'])
    parser.add_argument('case')
    parser.add_argument('--cores',type=int,required=True)
    parser.add_argument('--memory-gib',type=int,required=True)
    parser.add_argument('--minutes',type=int,required=True)
    parser.add_argument('--expected-minutes',type=float,required=True)
    parser.add_argument('--pilot-steps',type=int,default=0)
    parser.add_argument('--dependency')
    parser.add_argument('--reason',required=True)
    args=parser.parse_args()
    if not 1<=args.cores<=128:raise ValueError('Invalid shared-node core count')
    if args.minutes<args.expected_minutes:raise ValueError('Wall limit below expected runtime')
    jobs=ROOT/'jobs.json';records=json.loads(jobs.read_text()) if jobs.exists() else []
    tag=f'{args.phase}-{args.case}-t{args.cores}'+(f'-s{args.pilot_steps}' if args.pilot_steps else '')
    if any(r['tag']==tag for r in records):raise ValueError('Submission already recorded; inspect before retry')
    path=WORK/'slurm';path.mkdir(parents=True,exist_ok=True)
    script=path/(tag+'.sh')
    # A single immutable zipapp prevents queued jobs from running later edits,
    # and uses fewer inodes than copying a package tree for every submission.
    bundle=path/(tag+'.pyz')
    with zipfile.ZipFile(bundle,'x',compression=zipfile.ZIP_DEFLATED) as archive:
        archive.write(ROOT/'campaign.py','__main__.py')
        archive.write(ROOT/'common.py','common.py')
        for p in sorted((ROOT/'ic').glob('*.py')):archive.write(p,'ic/'+p.name)
        archive.writestr('ic/__init__.py','')
    command=['micromamba','run','-n','cosemu','python3','-B',str(bundle),args.phase,args.case,'--threads',str(args.cores)]
    if args.pilot_steps:command+=['--pilot-steps',str(args.pilot_steps)]
    text=f'''#!/bin/bash
#SBATCH --job-name=bdmconv-{tag}
#SBATCH --account=dp004
#SBATCH --partition=cosma8-serial
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={args.cores}
#SBATCH --mem={args.memory_gib}G
#SBATCH --time={args.minutes}
#SBATCH --output={path}/{tag}-%j.log
set -euo pipefail
export LINES=40 COLUMNS=120
module purge
module load intel_comp/2024.2.0
module load compiler-rt tbb compiler
export BDM_AUDIT_NATIVE_LIBS="$LD_LIBRARY_PATH"
export BDM_CONVERGENCE_ROOT={shlex.quote(str(ROOT))}
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OMP_DYNAMIC=FALSE OMP_PROC_BIND=close OMP_PLACES=cores OMP_STACKSIZE=128M
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
ulimit -c 0
cd {shlex.quote(str(REPO))}
scontrol show job "$SLURM_JOB_ID"
{shlex.join(command)}
'''
    script.write_text(text)
    command=['sbatch','--parsable']
    if args.dependency:command.append('--dependency=afterok:'+args.dependency)
    command.append(str(script))
    result=subprocess.run(command,capture_output=True,text=True,check=True)
    job=result.stdout.strip().split(';')[0]
    if not job.isdigit():raise RuntimeError(result.stdout)
    records.append(dict(tag=tag,job_id=job,submitted_at_utc=now(),partition='cosma8-serial',
                        cores=args.cores,memory_gib=args.memory_gib,wall_minutes=args.minutes,
                        expected_minutes=args.expected_minutes,expected_core_hours=args.cores*args.expected_minutes/60,
                        time_limit_core_hours=args.cores*args.minutes/60,reason=args.reason,
                        dependency=args.dependency,script=str(script),script_sha256=sha(script),
                        driver_sha256=sha(ROOT/'campaign.py'),bundle_sha256=sha(bundle),git_commit=git('rev-parse','HEAD'),
                        submit_command=command))
    write_json(jobs,records)
    print(job,tag,flush=True)


if __name__=='__main__':main()
