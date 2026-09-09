"""Submit small shared-node analysis jobs and record non-duplicated accounting."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import shlex
import subprocess

from measure import HERE, REPO, sha, source_hashes, write_json


def now():
    return datetime.now(timezone.utc).isoformat()


def submit(z, cores, memory, expected, reason):
    jobsfile = HERE/'jobs.json'
    jobs = json.loads(jobsfile.read_text()) if jobsfile.exists() else []
    if any(j['redshift'] == z for j in jobs):
        raise ValueError('Already submitted: inspect existing job before any rerun')
    work = HERE/'work'
    work.mkdir(exist_ok=True)
    script = work/f'F-z{z}.sh'
    hashes = source_hashes()
    command = ['micromamba', 'run', '-n', 'cosemu', 'python3', '-B', str(HERE/'measure.py'),
               '--redshift', str(z), '--threads', str(cores)]
    script.write_text(f'''#!/bin/bash
#SBATCH --job-name=bdmprop-F-z{z}
#SBATCH --account=dp004
#SBATCH --partition=cosma8-serial
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={cores}
#SBATCH --mem={memory}G
#SBATCH --time=00:15:00
#SBATCH --output={work}/F-z{z}-%j.log
set -euo pipefail
export LINES=40 COLUMNS=120
module purge
module load intel_comp/2024.2.0
module load compiler-rt tbb compiler
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OMP_DYNAMIC=FALSE OMP_PROC_BIND=close OMP_PLACES=cores
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
ulimit -c 0
cd {shlex.quote(str(REPO))}
scontrol show job "$SLURM_JOB_ID"
'''+''.join(f'test "$(sha256sum {shlex.quote(path)} | cut -d " " -f 1)" = {value}\n'
             for path, value in hashes.items()) + shlex.join(command)+'\n')
    cmd = ['sbatch', '--parsable', str(script)]
    job = subprocess.check_output(cmd, text=True).strip().split(';')[0]
    if not job.isdigit():
        raise ValueError(job)
    jobs.append(dict(job_id=job, redshift=z, submitted_at=now(), source_hashes=hashes,
                     script=str(script.relative_to(REPO)), script_sha256=sha(script),
                     partition='cosma8-serial', cores=cores, memory_gib=memory,
                     expected_minutes=expected, expected_core_hours=cores*expected/60,
                     time_limit_minutes=15, time_limit_core_hours=cores*15/60,
                     reason=reason, dependency=None, submit_command=cmd))
    write_json(jobsfile, jobs)
    print(job, flush=True)


def seconds(value):
    days = 0
    if '-' in value:
        days, value = value.split('-', 1)
    parts = list(map(float, value.split(':')))
    return int(days)*86400 + sum(x*60**i for i, x in enumerate(reversed(parts)))


def accounting():
    jobs = json.loads((HERE/'jobs.json').read_text())
    fields = ['JobIDRaw', 'State', 'ExitCode', 'ElapsedRaw', 'TotalCPU', 'MaxRSS',
              'ReqCPUS', 'AllocCPUS', 'ReqMem', 'ReqTRES', 'AllocTRES', 'TimelimitRaw', 'NodeList']
    cmd = ['sacct', '-j', ','.join(j['job_id'] for j in jobs), '-P', '-n', '--units=K',
           '--format='+','.join(fields)]
    raw = subprocess.check_output(cmd, text=True)
    rows = [dict(zip(fields, r.split('|'))) for r in raw.splitlines() if r.strip()]
    result = []
    for j in jobs:
        own = [r for r in rows if r['JobIDRaw'] == j['job_id'] or
               r['JobIDRaw'].startswith(j['job_id']+'.')]
        main = next((r for r in own if r['JobIDRaw'] == j['job_id']), None)
        if main is None:
            continue
        tres = dict(v.split('=', 1) for v in main['AllocTRES'].split(',') if '=' in v)
        billing = int(tres.get('billing', main['AllocCPUS'] or 0))
        elapsed = int(main['ElapsedRaw'] or 0)
        cpu = seconds(main['TotalCPU'])
        rss = [float(r['MaxRSS'].removesuffix('K')) for r in own if r['MaxRSS']]
        result.append(dict(job_id=j['job_id'], redshift=j['redshift'], state=main['State'],
                           exit_code=main['ExitCode'], elapsed_seconds=elapsed,
                           billing_cores=billing, billed_core_hours=billing*elapsed/3600,
                           cpu_hours=cpu/3600, cpu_efficiency=cpu/(billing*elapsed) if billing*elapsed else None,
                           maxrss_gib=max(rss)/1024**2 if rss else None,
                           request=main['ReqTRES'], allocation=main['AllocTRES']))
    report = dict(recorded_at=now(), command=cmd, raw_rows=rows, jobs=result,
                  billed_core_hours=sum(j['billed_core_hours'] for j in result),
                  cpu_hours=sum(j['cpu_hours'] for j in result),
                  expected_core_hours=sum(j['expected_core_hours'] for j in jobs),
                  time_limit_core_hours=sum(j['time_limit_core_hours'] for j in jobs))
    write_json(HERE/'accounting.json', report)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--submit', type=int, choices=[0, 1, 2])
    p.add_argument('--cores', type=int, default=2)
    p.add_argument('--memory-gib', type=int, default=8)
    p.add_argument('--expected-minutes', type=float, default=3)
    p.add_argument('--reason')
    p.add_argument('--accounting', action='store_true')
    a = p.parse_args()
    if a.submit is not None:
        if not a.reason:
            p.error('--reason is required')
        submit(a.submit, a.cores, a.memory_gib, a.expected_minutes, a.reason)
    if a.accounting:
        accounting()
