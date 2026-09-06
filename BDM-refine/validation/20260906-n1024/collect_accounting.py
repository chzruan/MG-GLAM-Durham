"""Record Slurm allocation, billing, CPU and resident-memory evidence."""
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import subprocess

ROOT = Path(__file__).resolve().parent


def seconds(value):
    if not value or value in ('Unknown', 'UNLIMITED', 'Partition_Limit'):
        return None
    days, _, rest = value.partition('-')
    if not rest:
        rest, days = days, '0'
    parts = [float(v) for v in rest.split(':')]
    return float(days)*86400 + sum(n * 60**i for i, n in enumerate(reversed(parts)))


def tres(value):
    return dict(part.split('=', 1) for part in value.split(',') if '=' in part)


def rss_kib(value):
    if not value:
        return 0.
    match = re.fullmatch(r'([0-9.]+)([KMGT]?)', value)
    assert match, value
    number, unit = match.groups()
    return float(number) * {'': 1/1024, 'K': 1, 'M': 1024, 'G': 1024**2, 'T': 1024**3}[unit]


def main():
    jobs = {j['job_id'] for j in json.loads((ROOT/'jobs.json').read_text())}
    for name in ['numerical-threading/results.json', 'numerical-threading/attempt-11948363.json',
                 'main-thread-control/results.json']:
        path = ROOT / name
        if path.exists():
            receipt = json.loads(path.read_text())
            jobs.add(str(receipt.get('results', receipt)['job_id']))
    fields = ['JobIDRaw', 'JobName', 'Partition', 'Account', 'State', 'ExitCode',
              'ReqTRES', 'AllocTRES', 'Elapsed', 'ElapsedRaw', 'TotalCPU', 'MaxRSS', 'Timelimit']
    command = ['sacct', '-P', '-j', ','.join(sorted(jobs)), '--units=K', '--format='+','.join(fields)]
    output = subprocess.check_output(command, text=True)
    lines = output.strip().splitlines()
    columns = lines[0].split('|')
    records = [dict(zip(columns, line.split('|'))) for line in lines[1:]]
    summaries = {}
    for job in sorted(jobs):
        parent = next(r for r in records if r['JobIDRaw'] == job)
        allocated = tres(parent['AllocTRES'])
        requested = tres(parent['ReqTRES'])
        elapsed = int(parent['ElapsedRaw'])
        cpus, billing = int(allocated.get('cpu', 0)), float(allocated.get('billing', 0))
        cpu = seconds(parent['TotalCPU']) or 0
        memory = max(rss_kib(r['MaxRSS']) for r in records if r['JobIDRaw'].split('.')[0] == job)
        limit = seconds(parent['Timelimit'])
        summaries[job] = dict(state=parent['State'], exit_code=parent['ExitCode'],
            elapsed_seconds=elapsed, total_cpu_seconds=cpu, allocated_cpus=cpus,
            requested_cpus=int(requested.get('cpu', 0)),
            billing_units=billing, maxrss_kib=memory, maxrss_gib=memory/1024**2,
            allocated_core_hours=cpus*elapsed/3600, billed_core_hours=billing*elapsed/3600,
            utilisation=cpu/(cpus*elapsed) if cpus*elapsed else None,
            time_limit_core_hours=int(requested.get('cpu', cpus))*limit/3600 if limit is not None else None,
            req_tres=parent['ReqTRES'], alloc_tres=parent['AllocTRES'])
    report = dict(collected_at_utc=datetime.now(timezone.utc).isoformat(), command=command,
                  raw_sacct=output, records=records, jobs=summaries,
                  total_billed_core_hours=sum(j['billed_core_hours'] for j in summaries.values()),
                  all_finished=all(j['state'].split()[0] not in ['RUNNING', 'PENDING', 'COMPLETING']
                                   for j in summaries.values()))
    path = ROOT/'accounting.json'
    path.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps({j: {'state': r['state'], 'elapsed': r['elapsed_seconds'],
                         'MaxRSS_GiB': r['maxrss_gib'], 'billed_core_hours': r['billed_core_hours']}
                      for j, r in summaries.items()}, indent=2))


if __name__ == '__main__':
    main()
