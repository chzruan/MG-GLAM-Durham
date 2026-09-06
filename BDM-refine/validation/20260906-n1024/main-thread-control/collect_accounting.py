"""Collect final Slurm evidence without modifying the sealed execution receipt."""
import csv
import io
import json
import re
import subprocess

import run_control as c


def cpu_seconds(value):
    days, clock = (value.split('-', 1) if '-' in value else ('0', value))
    fields = [float(piece) for piece in clock.split(':')]
    assert 1 <= len(fields) <= 3
    return 86400*int(days) + sum(part*60**index for index, part in enumerate(reversed(fields)))


def tres(value):
    return dict(item.split('=', 1) for item in value.split(',') if '=' in item)


def main():
    root = c.ROOT
    submission = json.loads((root / 'submission.json').read_text())
    job = submission['job_id']
    command = ['sacct', '-j', job, '--parsable2',
               '--format=JobID,State,ReqTRES,AllocTRES,ElapsedRaw,TotalCPU,CPUTimeRAW,MaxRSS,ExitCode']
    process = subprocess.run(command, capture_output=True, text=True, check=True)
    rows = list(csv.DictReader(io.StringIO(process.stdout), delimiter='|'))
    main = next(row for row in rows if row['JobID'] == job)
    batch = next(row for row in rows if row['JobID'] == job+'.batch')
    assert main['State'] in ['COMPLETED', 'FAILED', 'TIMEOUT', 'OUT_OF_MEMORY', 'CANCELLED'], 'Accounting is not final'
    assert batch['MaxRSS'], 'Wait for batch RSS accounting'
    requested, allocated = tres(main['ReqTRES']), tres(main['AllocTRES'])
    assert requested['cpu'] == allocated['cpu'] == allocated['billing'] == '64'
    assert requested['mem'] == allocated['mem'] == '288G'
    elapsed = int(main['ElapsedRaw'])
    cpu = cpu_seconds(main['TotalCPU'])
    matched = re.fullmatch(r'([\d.]+)([KMGTP]?)', batch['MaxRSS'])
    assert matched
    multiplier = {'':1/1024, 'K':1, 'M':1024, 'G':1024**2, 'T':1024**3, 'P':1024**4}[matched[2]]
    rss = float(matched[1])*multiplier
    control_command = ['scontrol', 'show', 'job', job]
    control = subprocess.run(control_command, capture_output=True, text=True)
    result_path = root / 'results.json'
    results = json.loads(result_path.read_text())
    report = dict(collected_at_utc=c.now(), job_id=job, state=main['State'], exit_code=main['ExitCode'],
        results_sha256=c.sha(result_path), submission_sha256=c.sha(root / 'submission.json'),
        collector_sha256=c.sha(__file__), sacct_command=command, sacct=process.stdout,
        scontrol_command=control_command, scontrol=control.stdout, scontrol_stderr=control.stderr,
        scontrol_returncode=control.returncode,
        allocation=dict(requested=requested, allocated=allocated, elapsed_seconds=elapsed,
            cpu_seconds=cpu, cpu_core_hours=cpu/3600, allocated_core_hours=64*elapsed/3600,
            cpu_utilization=cpu/(64*elapsed), batch_maxrss_kib=rss, batch_maxrss_gib=rss/1024**2,
            process_maxrss_kib=results.get('maxrss_kib'), time_limit_core_hours=48,
            headroom_fraction_of_request=1-rss/(288*1024**2)),
        scientific_completed=results['completed'], fixed_density_controls_passed=results['fixed_density_controls_passed'])
    c.write_json(root / 'accounting.json', report)
    print(json.dumps(report['allocation'], indent=2))


if __name__ == '__main__':
    main()
