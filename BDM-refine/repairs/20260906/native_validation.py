"""Measured native replays of verified 64^3/128^3 GR snapshots.

Use one core for an authorized light local run; use the shared Slurm script
for parallel scaling. Snapshots and historical audit results are read-only.
"""
import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import json
import os
from pathlib import Path
import resource
import signal
import socket
import subprocess
import tarfile
import tempfile
import time
import traceback

os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['MKL_NUM_THREADS']='1'
import numpy as np
from scipy.spatial import cKDTree

ROOT=Path(__file__).resolve().parent
REPO=ROOT.parents[2]
AUDIT=REPO/'BDM-refine/analysis/full-audit-20260906'
BUILD=ROOT/'work/native-build'


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        while chunk:=f.read(1024*1024):h.update(chunk)
    return h.hexdigest()


def snapshots():
    archive=AUDIT/'work-artifacts.tar.gz'
    expected=json.loads((AUDIT/'validation.json').read_text())['cleanup']['archive_sha256']
    assert sha(archive)==expected
    report={}
    with tarfile.open(archive,'r:gz') as tar:
        for nrow in [64,128]:
            original=json.loads((AUDIT/f'simulation-n{nrow}.json').read_text())
            assert original['completed'] and original['particles']==nrow**3<1200**3
            destination=ROOT/f'work/snapshots/n{nrow}'
            destination.mkdir(parents=True,exist_ok=True)
            for name,value in original['snapshot_sha256'].items():
                path=destination/name
                if not path.exists():
                    with tar.extractfile(f'work/gr-n{nrow}/Run1/{name}') as src,path.open('xb') as out:
                        while chunk:=src.read(1024*1024):out.write(chunk)
                assert sha(path)==value,name
            report[nrow]=dict(directory=destination,step=original['snapshot_step'],
                              snapshot_sha256=original['snapshot_sha256'])
    return report


def memberships(path,nrow,output):
    raw=path.read_bytes()
    selected,candidates=np.frombuffer(raw[:16],dtype='>i8').astype(np.int64)
    mass_one=float(np.frombuffer(raw[16:20],dtype='>f4')[0])
    cursor=20;ids=[];offsets=[0];indices=[];properties=[];seen=set()
    for _ in range(int(selected)):
        candidate,n=np.frombuffer(raw[cursor:cursor+16],dtype='>i8');cursor+=16
        values=np.frombuffer(raw[cursor:cursor+84],dtype='>f4').astype(np.float64);cursor+=84
        members=np.frombuffer(raw[cursor:cursor+int(n)*8],dtype='>i8').astype(np.int64);cursor+=int(n)*8
        assert len(members)==n and n>=20 and np.isfinite(values).all()
        assert members.min()>=1 and members.max()<=nrow**3 and np.all(np.diff(members)>0)
        assert abs(values[6]-n*mass_one)<=2*abs(float(np.spacing(np.float32(values[6]))))
        key=members.tobytes();assert key not in seen,'identical surviving particle sets'
        seen.add(key);ids.append(members);offsets.append(offsets[-1]+len(members))
        indices.append(int(candidate));properties.append(values)
    assert cursor==len(raw)
    properties=np.asarray(properties).reshape(-1,21)
    np.savez_compressed(output,ids=np.concatenate(ids) if ids else np.empty(0,dtype=np.int64),
        offsets=np.asarray(offsets,dtype=np.int64),candidates=np.asarray(indices),properties=properties)
    overlaps=[];host_violations=[]
    if len(ids)>1:
        tree=cKDTree(properties[:,:3]%nrow,boxsize=nrow)
        for a,b in tree.query_pairs(max(.2,2*float(properties[:,8].max())),output_type='ndarray'):
            delta=properties[a,:3]-properties[b,:3]
            delta-=nrow*np.rint(delta/nrow)
            priority=lambda i:(properties[i,6],-indices[i])
            host=a if priority(a)>priority(b) else b
            if np.linalg.norm(delta)<properties[host,8]:
                host_violations.append(dict(candidates=[indices[a],indices[b]],
                    separation=float(np.linalg.norm(delta)),host_radius=float(properties[host,8])))
            shared=len(np.intersect1d(ids[a],ids[b],assume_unique=True))
            if shared:overlaps.append(dict(candidates=[indices[a],indices[b]],shared=shared,
                fraction_of_smaller=shared/min(len(ids[a]),len(ids[b]))))
    assert not host_violations,host_violations
    return dict(candidates=int(candidates),selected=int(selected),mass_one=mass_one,
        exact_duplicate_member_sets=0,repeated_original_ids=0,mass_count_mismatches=0,
        host_exclusion_violations=0,nearby_overlaps=overlaps,raw_sha256=hashlib.sha256(raw).hexdigest(),
        retained_membership_archive=output.name)


def catalogue_agreement(reference,comparison):
    """Require the same rows to within one printed unit, across all 24 fields."""
    assert reference.shape==comparison.shape,(reference.shape,comparison.shape)
    # Peak compaction and the final selection order are deterministic. A row
    # permutation is a failure here, not an opportunity for ambiguous matching.
    quantum=np.maximum(np.abs(reference),np.abs(comparison))
    quantum=10.**(np.floor(np.log10(np.maximum(quantum,1.e-30)))-3.)
    quantum[:,:3]=1.e-4;quantum[:,3:6]=.01
    quantum[:,8]/=10.  # Radius is printed with five significant figures.
    quantum[:,[11,13,14]]=0.  # Halo ID, particle count, distinct/sub flag.
    difference=np.abs(reference-comparison)
    permitted=1.01*quantum+1.e-12
    result=dict(rows=len(reference),identical=bool(np.array_equal(reference,comparison)),
        max_absolute_difference_by_column=difference.max(axis=0).tolist(),
        max_printed_units_by_column=np.max(difference/np.maximum(quantum,1.e-12),axis=0).tolist(),
        tolerance='one printed unit per continuous field; exact ID, particle count and distinct flag')
    assert np.all(difference<=permitted),result
    return result


def one_replay(snapshot,nrow,variant,threads,repeat,log,records):
    binary=BUILD/('PMP2BDM.exe' if variant=='baseline' else f'PMP2BDM.{variant}.exe')
    with tempfile.TemporaryDirectory(prefix='native-case-',dir=ROOT/'work') as temp:
        work=Path(temp);(work/'CATALOGS').mkdir()
        for source in snapshot['directory'].iterdir():(work/source.name).symlink_to(source.resolve())
        (work/'BDM.config').write_text('iVirial=1\nMassMin=2.5e12\nRext=0.15\n')
        env={**os.environ,'OMP_NUM_THREADS':str(threads),'OMP_DYNAMIC':'FALSE',
             'OMP_PROC_BIND':'close','OMP_PLACES':'cores','LD_LIBRARY_PATH':os.environ['BDM_AUDIT_NATIVE_LIBS']}
        start=time.monotonic()
        p=subprocess.Popen(['/usr/bin/time','-f','%e %U %S %M','-o',str(work/'time.txt'),str(binary)],
            stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,
            cwd=work,env=env,text=True,start_new_session=True)
        timed_out=False
        try:
            stdout,stderr=p.communicate(f"{snapshot['step']}\n",timeout=180)
        except subprocess.TimeoutExpired:
            timed_out=True
            os.killpg(p.pid,signal.SIGKILL)
            stdout,stderr=p.communicate()
        elapsed=time.monotonic()-start
        timing=(work/'time.txt').read_text().splitlines() if (work/'time.txt').exists() else []
        timing=timing[-1].split() if timing else []
        timed=not timed_out and len(timing)==4
        record=dict(nrow=nrow,variant=variant,threads=threads,repeat=repeat,returncode=p.returncode,
            elapsed_seconds=elapsed,cpu_user_seconds=float(timing[1]) if timed else None,
            cpu_system_seconds=float(timing[2]) if timed else None,
            maxrss_kib=int(timing[3]) if timed else None,timed_out=timed_out,
            binary_sha256=sha(binary),stdout_tail=stdout.splitlines()[-30:],stderr=stderr,verification_completed=False)
        records.append(record)
        log.write(f'\n=== n{nrow} {variant} threads={threads} repeat={repeat} ===\n{stdout}\n{stderr}\n')
        log.flush()
        if p.returncode!=0 or timed_out:return record,None
        files=list((work/'CATALOGS').glob('Catshort*.DAT'))
        data=None
        if files:
            assert len(files)==1
            record['catalogue_sha256']=sha(files[0]);record['header']=files[0].read_text().splitlines()[:8]
            data=np.loadtxt(files[0],skiprows=8,ndmin=2)
            if data.size==0:data=np.empty((0,24))
            assert data.shape[1]==24 and np.isfinite(data).all()
            assert np.all(data[:,6]<=data[:,7]) and np.all(data[:,8]>0) and np.all(data[:,10]>=0)
            assert np.all(data[:,:3]>=0) and np.all(data[:,:3]<nrow)
            physical=np.delete(data,11,axis=1);physical=physical[np.lexsort(physical[:,:3].T[::-1])]
            record.update(rows=len(data),nonfinite_values=0,canonical_sha256=hashlib.sha256(physical.tobytes()).hexdigest())
        if variant=='members' and p.returncode==0:
            record['membership']=memberships(work/'repair-members.bin',nrow,ROOT/f'membership-n{nrow}.npz')
            assert record['membership']['selected']==record['rows']
        if variant=='state':record['two_calls_bitwise_identical']='TWO CALLS BITWISE IDENTICAL' in stdout
        record['verification_completed']=True
        return record,data


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--phase',choices=['pilot','n128','scaling'],default='pilot')
    parser.add_argument('--threads',type=int,default=1)
    args=parser.parse_args()
    if 'SLURM_JOB_ID' not in os.environ:
        assert args.threads==1 and args.phase!='scaling','Parallel timing belongs on Slurm'
    else:assert 1<=args.threads<=int(os.environ['SLURM_CPUS_PER_TASK'])
    if args.phase=='scaling':assert args.threads in [1,2,4,8]
    resource.setrlimit(resource.RLIMIT_CORE,(0,0))
    original=snapshots()
    build=json.loads((ROOT/'native-build.json').read_text());assert build['completed']
    for name,value in build['binaries_sha256'].items():assert sha(BUILD/name)==value
    name=f'native-{args.phase}'+(f'-t{args.threads}' if args.phase=='scaling' else '')
    records=[];arrays={}
    report=dict(started_at_utc=datetime.now(timezone.utc).isoformat(),host=socket.gethostname(),
        job_id=os.environ.get('SLURM_JOB_ID'),phase=args.phase,threads=args.threads,
        source_sha256=build['source_sha256']['PMP2linker.f90'],build_commit=build['source_commit'],
        driver_sha256=sha(__file__),records=records,completed=False)
    try:
        with gzip.open(ROOT/(name+'.log.gz'),'wt') as log:
            nrow=64 if args.phase=='pilot' else 128
            tasks=[('baseline',i,t) for i in range(3) for t in [1,2,4,8] if t<=args.threads] if args.phase=='scaling' else [
                (v,0,args.threads) for v in ['baseline','members','bounds','state']]
            for variant,repeat,threads in tasks:
                row,data=one_replay(original[nrow],nrow,variant,threads,repeat,log,records)
                if data is not None:arrays[f'{variant}_{repeat}_t{threads}']=data
                print(nrow,variant,threads,row['returncode'],row.get('rows'),round(row['elapsed_seconds'],3),flush=True)
                assert row['returncode']==0 and row.get('rows',0)>0,row
                if variant=='state':assert row['two_calls_bitwise_identical']
            if args.phase!='scaling':
                assert records[0]['catalogue_sha256']==records[1]['catalogue_sha256'],'membership instrumentation changed output'
                assert records[0]['catalogue_sha256']==records[3]['catalogue_sha256'],'repeat/state probe changed output'
                report['checked_optimized_agreement']=catalogue_agreement(arrays[f'baseline_0_t{args.threads}'],arrays[f'bounds_0_t{args.threads}'])
            else:
                assert len({r['canonical_sha256'] for r in records})==1,'cross-thread/repeatability failure'
                report['cross_thread_physical_columns_identical']=True
                report['catalogue_sha256_count']=len({r['catalogue_sha256'] for r in records})
                report['full_catalogues_byte_identical']=report['catalogue_sha256_count']==1
            report['completed']=True
    except BaseException:
        report['failure_traceback']=traceback.format_exc()
        raise
    finally:
        np.savez_compressed(ROOT/(name+'-catalogues.npz'),**arrays)
        (ROOT/(name+'.json')).write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':main()
