"""Collect immutable Slurm accounting; micromamba run -n cosemu python3 -B."""
import argparse
import json
from pathlib import Path
import subprocess

import run_replays as r


def seconds(value):
    days=0
    if '-' in value:
        day,value=value.split('-',1);days=int(day)
    parts=list(map(float,value.split(':')))
    total=0.
    for number in parts:total=60*total+number
    return 86400*days+total


def tres(value):
    return dict(x.split('=',1) for x in value.split(',') if '=' in x)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--jobs',nargs='+',required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    assert not args.output.exists(),'Retain earlier accounting evidence'
    command=['sacct','-j',','.join(args.jobs),
             '--format=JobID,JobName%35,State,Partition,ElapsedRaw,TimelimitRaw,TotalCPU,ReqTRES%100,AllocTRES%100,MaxRSS,ExitCode','--parsable2']
    output=subprocess.check_output(command,text=True)
    lines=output.strip().splitlines();keys=lines[0].split('|')
    rows=[dict(zip(keys,line.split('|'))) for line in lines[1:]]
    jobs={}
    for job in args.jobs:
        row=next(x for x in rows if x['JobID']==job)
        batch=next(x for x in rows if x['JobID']==job+'.batch')
        requested,allocated=tres(row['ReqTRES']),tres(row['AllocTRES'])
        cores=int(allocated['cpu']);billing=int(allocated['billing'])
        elapsed=int(row['ElapsedRaw']);cpu=seconds(row['TotalCPU'])
        query=subprocess.run(['scontrol','show','job',job,'-o'],capture_output=True,text=True)
        jobs[job]=dict(state=row['State'],exit_code=row['ExitCode'],partition=row['Partition'],
                       requested=requested,allocated=allocated,batch_maxrss=batch['MaxRSS'],
                       elapsed_seconds=elapsed,total_cpu_seconds=cpu,
                       cpu_efficiency=cpu/(elapsed*cores) if elapsed else 0.,
                       allocated_core_hours=elapsed*cores/3600,billed_core_hours=elapsed*billing/3600,
                       time_limit_core_hours=int(row['TimelimitRaw'])*cores/60,
                       scontrol=query.stdout,scontrol_returncode=query.returncode)
    result=dict(collected_at_utc=r.now(),collector_sha256=r.sha(__file__),command=command,raw_rows=rows,jobs=jobs,
                all_completed=all(x['state']=='COMPLETED' and x['exit_code']=='0:0' for x in jobs.values()),
                allocated_core_hours=sum(x['allocated_core_hours'] for x in jobs.values()),
                billed_core_hours=sum(x['billed_core_hours'] for x in jobs.values()),
                time_limit_core_hours=sum(x['time_limit_core_hours'] for x in jobs.values()))
    r.write_json(args.output,result)
    print(json.dumps({key:result[key] for key in ['all_completed','allocated_core_hours','billed_core_hours','time_limit_core_hours']}))


if __name__=='__main__':
    main()
