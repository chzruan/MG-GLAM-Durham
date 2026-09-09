"""Record Slurm allocations and actual billing, without double-counting job steps."""
import json
import subprocess
from common import ROOT, now, write_json


def seconds(text):
    if not text or text in ['Unknown','INVALID']:return None
    days=0
    if '-' in text:
        day,text=text.split('-',1);days=int(day)
    parts=list(map(float,text.split(':')))
    if len(parts)==3:return days*86400+parts[0]*3600+parts[1]*60+parts[2]
    if len(parts)==2:return days*86400+parts[0]*60+parts[1]
    return float(text)


def main():
    jobs=json.loads((ROOT/'jobs.json').read_text())
    ids=','.join(r['job_id'] for r in jobs)
    fields=['JobIDRaw','State','ExitCode','ElapsedRaw','TotalCPU','MaxRSS','ReqCPUS','AllocCPUS',
            'ReqMem','ReqTRES','AllocTRES','TimelimitRaw','NodeList']
    command=['sacct','-j',ids,'-P','-n','--units=K','--format='+','.join(fields)]
    result=subprocess.run(command,capture_output=True,text=True,check=True)
    rows=[dict(zip(fields,line.split('|'))) for line in result.stdout.splitlines() if line.strip()]
    report=dict(recorded_at_utc=now(),command=command,jobs=[],raw_rows=rows)
    for job in jobs:
        own=[r for r in rows if r['JobIDRaw']==job['job_id'] or r['JobIDRaw'].startswith(job['job_id']+'.')]
        main=next((r for r in own if r['JobIDRaw']==job['job_id']),None)
        if main is None:continue
        tres=dict(v.split('=',1) for v in main['AllocTRES'].split(',') if '=' in v)
        billing=int(tres.get('billing',main['AllocCPUS'] or 0));elapsed=int(main['ElapsedRaw'] or 0)
        cpu=seconds(main['TotalCPU'])
        rss=[float(r['MaxRSS'].removesuffix('K')) for r in own if r['MaxRSS']]
        report['jobs'].append(dict(job_id=job['job_id'],tag=job['tag'],state=main['State'],exit_code=main['ExitCode'],
                                  elapsed_seconds=elapsed,billing_cores=billing,billed_core_hours=billing*elapsed/3600,
                                  cpu_hours=None if cpu is None else cpu/3600,
                                  cpu_efficiency=None if not cpu or not elapsed or not billing else cpu/(billing*elapsed),
                                  maxrss_gib=max(rss)/1024**2 if rss else None,
                                  request=main['ReqTRES'],allocation=main['AllocTRES'],
                                  time_limit_core_hours=job['time_limit_core_hours']))
    report['billed_core_hours_so_far']=sum(j['billed_core_hours'] for j in report['jobs'])
    report['submitted_expected_core_hours']=sum(j['expected_core_hours'] for j in jobs)
    report['submitted_time_limit_core_hours']=sum(j['time_limit_core_hours'] for j in jobs)
    write_json(ROOT/'accounting.json',report)
    for row in report['jobs']:print(json.dumps(row))


if __name__=='__main__':main()
