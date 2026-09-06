"""Consolidate evidence and Slurm accounting; does not rerun simulations."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import subprocess

import numpy as np
from scipy.spatial import cKDTree

ROOT=Path(__file__).resolve().parent
REPO=ROOT.parents[2]


def clock_seconds(value):
    days=0
    if '-' in value:
        day,value=value.split('-');days=int(day)
    parts=[float(v) for v in value.split(':')]
    return days*86400+sum(v*60**i for i,v in enumerate(reversed(parts)))


def compare(left,right,box):
    a=cKDTree(left[:,:3]%box,boxsize=box)
    b=cKDTree(right[:,:3]%box,boxsize=box)
    distance,forward=b.query(left[:,:3]%box)
    _,backward=a.query(right[:,:3]%box)
    keep=(distance<.01)&(backward[forward]==np.arange(len(left)))
    matched=right[forward[keep]];reference=left[keep]
    fields={}
    for name,col in [('Mbound',6),('Mtot',7),('Rvir',8),('Vrms',9),('Vmax',10),('concentration',12)]:
        fields[name]=dict(changed_rows=int(np.sum(reference[:,col]!=matched[:,col])),
            max_relative_difference=float(np.max(np.abs(matched[:,col]/reference[:,col]-1))))
    return dict(reciprocal_matches=int(keep.sum()),left_unmatched=int((~keep).sum()),
        right_unmatched=len(right)-int(keep.sum()),
        matching_radius_mpc_h=.01,max_matched_displacement_mpc_h=float(distance[keep].max()),fields=fields)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--refresh-accounting',action='store_true',
                        help='Query sacct; otherwise use the preserved final accounting')
    args=parser.parse_args()
    unit=json.loads((ROOT/'unit-results.json').read_text())
    replay=json.loads((ROOT/'replay-results.json').read_text())
    followup=json.loads((ROOT/'bounds-results.json').read_text())
    arrays=np.load(ROOT/'replay-catalogues.npz')
    bounds=np.load(ROOT/'bounds-catalogues.npz')
    jobs=json.loads((ROOT/'jobs.json').read_text())
    if args.refresh_accounting:
        fields='JobID,State,ElapsedRaw,Elapsed,TotalCPU,AllocCPUS,ReqTRES,AllocTRES,MaxRSS,ExitCode'
        accounting=subprocess.check_output(['sacct','-j',','.join(j['job_id'] for j in jobs),
            '--format='+fields,'-P'],text=True)
        lines=accounting.strip().splitlines();keys=lines[0].split('|')
        rows=[dict(zip(keys,line.split('|'))) for line in lines[1:]]
        by_job={r['JobID']:r for r in rows}
        for job in jobs:
            job['final_accounting']=by_job[job['job_id']]
            job['steps']=[r for r in rows if r['JobID'].startswith(job['job_id']+'.')]
    totals=dict(billed_core_hours=0.,cpu_core_hours=0.,time_limit_core_hours=0.)
    for job in jobs:
        row=job['final_accounting']
        assert row['State'] in {'COMPLETED','FAILED','CANCELLED','TIMEOUT','OUT_OF_MEMORY'}
        job['billed_core_hours']=int(row['ElapsedRaw'])*int(row['AllocCPUS'])/3600
        job['cpu_core_hours']=clock_seconds(row['TotalCPU'])/3600
        job['cpu_utilization']=job['cpu_core_hours']/job['billed_core_hours'] if job['billed_core_hours'] else None
        for name in totals: totals[name]+=job[name]
    if args.refresh_accounting:
        (ROOT/'jobs.json').write_text(json.dumps(jobs,indent=2)+'\n')
    scale=[]
    for threads in [1,2,4,8]:
        samples=[r for r in replay['records'] if r.get('purpose') and r['threads']==threads]
        lists=[r['seconds'] for r in replay['list_benchmarks'] if r['threads']==threads]
        scale.append(dict(threads=threads,finder_median_seconds=float(np.median([r['elapsed_seconds'] for r in samples])),
            finder_median_cpu_seconds=float(np.median([r['cpu_user_seconds']+r['cpu_system_seconds'] for r in samples])),
            list_median_seconds=float(np.median(lists)),repeats=len(samples),
            distinct_canonical_catalogues=len(set(r['canonical_physical_rows_sha256'] for r in samples))))
    for row in scale:
        row['finder_speedup']=scale[0]['finder_median_seconds']/row['finder_median_seconds']
        row['finder_parallel_efficiency']=row['finder_speedup']/row['threads']
        row['finder_reserved_core_seconds']=row['finder_median_seconds']*row['threads']
        row['list_speedup']=scale[0]['list_median_seconds']/row['list_median_seconds']
    prefilter={}
    for variant in ['baseline','prefiltered']:
        samples=[r for r in followup['records'] if r['variant']==variant]
        prefilter[variant]=dict(median_seconds=float(np.median([r['elapsed_seconds'] for r in samples])),
            median_cpu_seconds=float(np.median([r['cpu_user_seconds']+r['cpu_system_seconds'] for r in samples])),
            canonical_sha256=sorted(set(r['canonical_physical_rows_sha256'] for r in samples)))
    prefilter['speedup']=prefilter['baseline']['median_seconds']/prefilter['prefiltered']['median_seconds']
    findings=[dict(id=key,priority=priority,priority_scope=scope.strip(),title=title,
                   stage=1 if int(key[1:])<=10 else 2)
              for key,priority,scope,title in re.findall(
                  r'^### (F\d+) — (P\d)([^:]*): (.+)$',(ROOT/'AUDIT.md').read_text(),re.M)]
    assert len(findings)==16
    summary=dict(completed_at_utc=datetime.now(timezone.utc).isoformat(),
        audited_commit=unit['baseline_commit'],audited_branch='audit/bdm-physics-numerics-20260905',
        verdict='changes_required',finding_groups=findings,
        source_changes_to_finder=False,unit_experiments=len(unit['results']),
        maximum_simulated_particles=128**3,particle_limit_exclusive=1200**3,
        accounting=totals,scaling=scale,prefilter_trial=prefilter,
        plateau_distinct_peak_sets=len(set(r['sorted_peak_array_sha256'] for r in unit['results'] if r['case']=='peak_plateau')),
        simulations={},native_bounds_comparisons={})
    for nrow in [64,128]:
        baseline=next(r for r in replay['records'] if r['nrow']==nrow and r['variant']=='baseline')
        member=next(r for r in replay['records'] if r['nrow']==nrow and r['variant']=='members')
        configured=next(r for r in replay['records'] if r['nrow']==nrow and r['variant']=='configured')
        summary['simulations'][str(nrow)]=dict(baseline_rows=baseline['rows'],
            configured_rows=configured['rows'],configuration_order_trial_increase_fraction=configured['rows']/baseline['rows']-1,
            instrumented_catalogue_byte_identical=baseline['catalogue_sha256']==member['catalogue_sha256'],
            membership=member['membership'])
        summary['native_bounds_comparisons'][str(nrow)]=compare(arrays[f'n{nrow}_baseline'],bounds[f'n{nrow}_bounds'],nrow)
    (ROOT/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    main()
