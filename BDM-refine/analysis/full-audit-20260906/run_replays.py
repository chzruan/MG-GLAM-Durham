"""Particle-membership validation followed by measured thread/CPU experiments."""
import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import json
import os
from pathlib import Path
import resource
import subprocess
import tempfile

import numpy as np
from scipy.spatial import cKDTree

from run_simulation import ROOT, REPO, BIN, run, sha


def catalogue_summary(data,box):
    import sys
    sys.path.insert(0,str(REPO/'BDM-refine/analysis/duplicate-hosts-20260905/tools'))
    from catalogue_core import COLUMNS,FIELDS,audit
    result=dict(rows=len(data),nonfinite_values=int((~np.isfinite(data)).sum()))
    physical=np.delete(data,11,axis=1)
    physical=physical[np.lexsort(physical[:,:3].T[::-1])]
    result['canonical_physical_rows_sha256']=hashlib.sha256(physical.tobytes()).hexdigest()
    if len(data) and np.isfinite(data).all():
        fields={key:data[:,COLUMNS.index(key)] for key in FIELDS}
        validation,_,drop,*_=audit(fields,box=box)
        result['strict_duplicate_rows']=int(drop.sum())
        result['remaining_same_bound_pairs_below_0p2']=validation['remaining_same_bound_pairs_below_0p2']
    return result


def membership_summary(work,box):
    raw=(work/'audit-survivors.bin').read_bytes()
    nmax=int(np.frombuffer(raw[:4],dtype='>i4')[0])
    mass_one=float(np.frombuffer(raw[4:8],dtype='>f4')[0])
    halos=np.frombuffer(raw[8:],dtype='>f4').reshape(5,nmax).T.astype(float)
    selected=np.flatnonzero((halos[:,0]>=max(2.5e12,20*mass_one)) &
        np.all((halos[:,2:]>=0)&(halos[:,2:]<box),axis=1))
    all_ids={}
    with (work/'audit-members.bin').open('rb') as f:
        while header:=f.read(16):
            assert len(header)==16
            candidate,n=np.frombuffer(header,dtype='>i8')
            values=np.frombuffer(f.read(int(n)*8),dtype='>i8')
            assert len(values)==n
            all_ids[int(candidate)-1]=np.sort(values)
    sets={int(i):all_ids[int(i)] for i in selected}
    grouped={};duplicates=[];repeated_ids=[];count_mismatch=[]
    for candidate,ids in sets.items():
        key=hashlib.sha256(ids.tobytes()).hexdigest()
        if key in grouped:
            assert np.array_equal(ids,sets[grouped[key]])
            duplicates.append([grouped[key],candidate])
        grouped[key]=candidate
        if len(np.unique(ids))!=len(ids): repeated_ids.append(candidate)
        if abs(halos[candidate,0]/mass_one-len(ids))>.02: count_mismatch.append(candidate)
    nearest=cKDTree(halos[selected,2:]%box,boxsize=box)
    pairs=nearest.query_pairs(max(.2,2*float(halos[selected,1].max())),output_type='ndarray')
    overlaps=[]
    for a,b in pairs:
        i,j=int(selected[a]),int(selected[b])
        overlap=len(np.intersect1d(sets[i],sets[j],assume_unique=True))
        if overlap:
            overlaps.append(dict(candidates=[i,j],shared=overlap,
                fraction_of_smaller=overlap/min(len(sets[i]),len(sets[j]))))
    return dict(candidates=nmax,selected=len(selected),mass_one=mass_one,
        exact_duplicate_member_sets=duplicates,repeated_original_ids_in_halos=repeated_ids,
        mass_count_mismatches=count_mismatch,nearby_overlapping_sets=overlaps,
        original_id_range=[min(int(ids.min()) for ids in sets.values()),max(int(ids.max()) for ids in sets.values())],
        selected_membership_sha256=hashlib.sha256(b''.join(sets[i].tobytes() for i in sorted(sets))).hexdigest())


def one_replay(nrow,threads,variant,log,repeat=0,mode=1):
    original=json.loads((ROOT/f'simulation-n{nrow}.json').read_text())
    step=original['snapshot_step'];snapshot=ROOT/'work'/f'gr-n{nrow}'/'Run1'
    assert all(sha(snapshot/name)==value for name,value in original['snapshot_sha256'].items())
    with tempfile.TemporaryDirectory(prefix='audit-replay-',dir=ROOT/'work') as tmp:
        work=Path(tmp);(work/'CATALOGS').mkdir()
        for source in snapshot.glob(f'PMcr*.{step:04d}.DAT'):
            (work/source.name).symlink_to(source.resolve())
        (work/'BDM.config').write_text(f'! audit configuration\niVirial = {mode} ! convention\n'
            'MassMin = 2.5e12 ! selection\nRext = 0.15 ! legacy correction\n')
        binary=BIN/('PMP2BDM.exe' if variant=='baseline' else f'PMP2BDM.{variant}.exe')
        record=run(binary,work,f'{step}\n',log,threads)
        record.update(nrow=nrow,variant=variant,repeat=repeat,iVirial_requested=mode)
        data=None;files=list((work/'CATALOGS').glob('Catshort*.DAT'))
        if files:
            record['catalogue_sha256']=sha(files[0])
            record['header']=files[0].read_text().splitlines()[:8]
            try:
                data=np.loadtxt(files[0],skiprows=8,ndmin=2)
                record.update(catalogue_summary(data,nrow))
            except (ValueError,IndexError) as error:
                record['catalogue_validation_error']=str(error)
        if variant=='members' and record['returncode']==0:
            record['membership']=membership_summary(work,nrow)
        return record,data


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--bounds-only',action='store_true')
    args=parser.parse_args()
    assert int(os.environ['SLURM_CPUS_PER_TASK'])>=(1 if args.bounds_only else 8)
    resource.setrlimit(resource.RLIMIT_CORE,(0,0))
    records=[];arrays={};benchmarks=[]
    report=dict(started_at_utc=datetime.now(timezone.utc).isoformat(),
        job_id=os.environ['SLURM_JOB_ID'],records=records,list_benchmarks=benchmarks)
    try:
        with gzip.open(ROOT/('bounds.log.gz' if args.bounds_only else 'replays.log.gz'),'wt') as log:
            # Finish membership/physics checks before the efficiency experiments.
            for nrow in [64,128]:
                for variant in (['bounds'] if args.bounds_only else ['baseline','members','configured','checked']):
                    record,data=one_replay(nrow,1,variant,log)
                    records.append(record)
                    if data is not None: arrays[f'n{nrow}_{variant}']=data
                    print(nrow,variant,record['returncode'],record.get('rows'),flush=True)
            if args.bounds_only:
                for repeat in range(3):
                    # Alternate order within one node/allocation to reduce cache
                    # and load biases in the isolated prefilter comparison.
                    variants=['baseline','prefiltered'] if repeat%2==0 else ['prefiltered','baseline']
                    for variant in variants:
                        record,data=one_replay(128,1,variant,log,repeat=repeat)
                        record['purpose']='paired prefilter timing'
                        records.append(record)
                        if data is not None: arrays[f'n128_{variant}_r{repeat}']=data
                report['completed']=True
                return
            for threads in [1,2,4,8]:
                for repeat in range(3):
                    record,data=one_replay(128,threads,'baseline',log,repeat=repeat)
                    record['purpose']='scaling and repeatability'
                    records.append(record)
                    if data is not None: arrays[f'n128_t{threads}_r{repeat}']=data
            for threads in [1,2,4,8]:
                for repeat in range(3):
                    env={**os.environ,'OMP_NUM_THREADS':str(threads),'OMP_DYNAMIC':'FALSE',
                         'OMP_PROC_BIND':'close','OMP_PLACES':'cores'}
                    p=subprocess.run([str(ROOT/'work/cases-build/optimized'),'list_benchmark','128'],
                        cwd=ROOT/'work',env=env,capture_output=True,text=True,check=True,timeout=60)
                    line=next(x for x in p.stdout.splitlines() if 'AUDIT LIST_SECONDS_COUNT' in x).split()
                    assert int(line[-1])==int(line[-2])==128**3
                    benchmarks.append(dict(threads=threads,repeat=repeat,seconds=float(line[-3]),
                        linked_particles=int(line[-2])))
            report['completed']=True
    finally:
        prefix='bounds' if args.bounds_only else 'replay'
        np.savez_compressed(ROOT/f'{prefix}-catalogues.npz',**arrays)
        (ROOT/f'{prefix}-results.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
    main()
