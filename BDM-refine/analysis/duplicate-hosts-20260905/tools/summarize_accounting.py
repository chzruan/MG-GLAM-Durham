"""Persist Slurm allocation, billing, peak memory, and CPU efficiency evidence."""
import argparse
import csv
import io
import json
from pathlib import Path
import subprocess


def seconds(value):
    days=0
    if '-' in value:
        d,value=value.split('-',1);days=int(d)
    parts=[float(x) for x in value.split(':')]
    total=days*86400
    for i,v in enumerate(reversed(parts)):
        total+=v*60**i
    return total


p=argparse.ArgumentParser()
p.add_argument('job_ids')
p.add_argument('output_root',type=Path)
a=p.parse_args()
cmd=['sacct','-j',a.job_ids,'--units=M','--parsable2',
     '--format=JobID,State,ElapsedRaw,TotalCPU,AllocCPUS,ReqTRES,AllocTRES,MaxRSS,ExitCode']
raw=subprocess.check_output(cmd,text=True)
(a.output_root/'all-jobs.sacct.txt').write_text(raw)
records=list(csv.DictReader(io.StringIO(raw),delimiter='|'))
steps={r['JobID']:r for r in records}
jobs=[]
for r in records:
    if '.' in r['JobID']:
        continue
    elapsed=int(r['ElapsedRaw']);cpus=int(r['AllocCPUS'])
    cpu=seconds(r['TotalCPU'])
    tres=dict(x.split('=',1) for x in r['AllocTRES'].split(',') if '=' in x)
    billing=int(tres.get('billing',cpus))
    batch=steps.get(r['JobID']+'.batch',{})
    jobs.append(dict(job_id=r['JobID'],state=r['State'],exit_code=r['ExitCode'],
        elapsed_seconds=elapsed,total_cpu_seconds=cpu,allocated_cpus=cpus,
        billed_cpus=billing,req_tres=r['ReqTRES'],alloc_tres=r['AllocTRES'],
        maxrss=batch.get('MaxRSS'),
        cpu_efficiency=cpu/(elapsed*cpus) if elapsed*cpus else None,
        allocated_core_hours=elapsed*cpus/3600,
        billed_core_hours=elapsed*billing/3600))
total_time=sum(j['elapsed_seconds']*j['allocated_cpus'] for j in jobs)
result=dict(jobs=jobs,allocated_core_hours=sum(j['allocated_core_hours'] for j in jobs),
            billed_core_hours=sum(j['billed_core_hours'] for j in jobs),
            total_cpu_hours=sum(j['total_cpu_seconds'] for j in jobs)/3600,
            aggregate_cpu_efficiency=sum(j['total_cpu_seconds'] for j in jobs)/total_time if total_time else None,
            all_completed=all(j['state']=='COMPLETED' and j['exit_code']=='0:0' for j in jobs),
            note='MaxRSS is the Slurm batch-step accounting value; per-tool ru_maxrss is also in receipts.')
(a.output_root/'accounting_summary.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
