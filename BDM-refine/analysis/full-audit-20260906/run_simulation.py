"""Small GR integration pilot. Run only in the supplied shared-queue allocation."""
import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import resource
import subprocess
import tempfile
import time

import numpy as np

ROOT=Path(__file__).resolve().parent
REPO=ROOT.parents[2]
BIN=ROOT/'work/full-build'


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(8*1024**2),b''):
            h.update(block)
    return h.hexdigest()


def run(executable,cwd,stdin,log,threads):
    env={**os.environ,'OMP_NUM_THREADS':str(threads),'OMP_DYNAMIC':'FALSE',
         'OMP_PROC_BIND':'close','OMP_PLACES':'cores','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1'}
    # micromamba remains the Python environment; native ifx subprocesses need
    # the matching compiler runtime, rather than cosemu's older libiomp5.
    env['LD_LIBRARY_PATH']=os.environ['BDM_AUDIT_NATIVE_LIBS']
    start=time.monotonic()
    timer=Path(cwd)/'timing.txt'
    p=subprocess.run(['/usr/bin/time','-f','%e %U %S %M','-o',str(timer),str(executable)],
        cwd=cwd,input=stdin,text=True,capture_output=True,env=env,timeout=700)
    elapsed=time.monotonic()-start
    log.write(f'\n=== {executable.name} {cwd} threads={threads} ===\n{p.stdout}\n{p.stderr}\n')
    values=timer.read_text().splitlines()[-1].split()
    timer.unlink()
    return dict(executable=str(executable),sha256=sha(executable),threads=threads,
        returncode=p.returncode,elapsed_seconds=elapsed,cpu_user_seconds=float(values[1]),
        cpu_system_seconds=float(values[2]),maxrss_kib=int(values[3]),
        stdout_sha256=hashlib.sha256(p.stdout.encode()).hexdigest(),
        stdout_tail=p.stdout.splitlines()[-25:],stderr=p.stderr)


def replay(snapshot_dir,step,threads,log,mode=1):
    with tempfile.TemporaryDirectory(prefix='replay-',dir=ROOT/'work') as temporary:
        work=Path(temporary);(work/'CATALOGS').mkdir()
        for source in snapshot_dir.glob(f'PMcr*.{step:04d}.DAT'):
            (work/source.name).symlink_to(source.resolve())
        (work/'BDM.config').write_text(f'! audit configuration\niVirial = {mode} ! convention\n'
            'MassMin = 2.5e12 ! selection\nRext = 0.15 ! legacy correction\n')
        record=run(BIN/'PMP2BDM.exe',work,f'{step}\n',log,threads)
        record['iVirial_requested']=mode
        files=list((work/'CATALOGS').glob('Catshort*.DAT'))
        record['catalogue_count']=len(files)
        data=None
        if files:
            catalogue=files[0]
            record['catalogue_sha256']=sha(catalogue)
            record['header']=catalogue.read_text().splitlines()[:8]
            try:
                data=np.loadtxt(catalogue,skiprows=8,ndmin=2)
                record.update(rows=len(data),columns=data.shape[1],nonfinite_values=int((~np.isfinite(data)).sum()))
                if data.size:
                    record['nonfinite_by_column']=(~np.isfinite(data)).sum(axis=0).tolist()
                    record['bound_exceeds_total_rows']=int(np.sum(data[:,6]>data[:,7]))
                    record['nonpositive_radius_rows']=int(np.sum(data[:,8]<=0))
            except ValueError as error:
                record['parse_error']=str(error)
        return record,data


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--nrow',type=int,default=64)
    p.add_argument('--threads',type=int,default=1)
    p.add_argument('--resume',action='store_true')
    args=p.parse_args()
    assert 0<args.nrow**3<1200**3
    assert 'SLURM_JOB_ID' in os.environ, 'Submit this workload with Slurm'
    assert args.threads<=int(os.environ['SLURM_CPUS_PER_TASK'])
    resource.setrlimit(resource.RLIMIT_CORE,(0,0))
    case=ROOT/'work'/f'gr-n{args.nrow}'
    case.mkdir(parents=True,exist_ok=args.resume)
    changes={'Box':float(args.nrow),'Nrow':args.nrow,'Ngrid':2*args.nrow,
        '#outputs':-1,'Steps between checkpoints':10000,'Save snapshots':1,
        'DM power spectrum':0,'Find BDM halos':0,'MG_flag':0}
    lines=[]
    for line in (REPO/'test/Init.dat').read_text().splitlines():
        key=line.split('=')[0].strip()
        if key in changes:
            line=f'{key} = {changes[key]}'
        lines.append(line)
    intended_init='\n'.join(lines)+'\n'
    if args.resume:
        assert (case/'Init.dat').read_text()==intended_init
    else:
        (case/'Init.dat').write_text(intended_init)
        for name in ['PkTable.dat','TableSeeds.dat']:
            (case/name).symlink_to((REPO/'test'/name).resolve())
    snapshot=case/'Run1';snapshot.mkdir(exist_ok=args.resume)
    (snapshot/'CATALOGS').mkdir(exist_ok=args.resume)
    records=[]
    report=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),job_id=os.environ['SLURM_JOB_ID'],
        branch=subprocess.check_output(['git','branch','--show-current'],cwd=REPO,text=True).strip(),
        baseline_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=REPO,text=True).strip(),
        nrow=args.nrow,ngrid=2*args.nrow,particles=args.nrow**3,box_mpc_h=args.nrow,
        gravity='GR',initial_conditions='repository GLAM IC generator; fixed realization 1, z_init=100',
        intended_use='integration and performance audit, not convergence calibration',
        inputs_sha256={name:sha(case/name) for name in ['Init.dat','PkTable.dat','TableSeeds.dat']},
        records=records)
    result_path=ROOT/f'simulation-n{args.nrow}.json'
    if args.resume:
        previous=json.loads(result_path.read_text())
        assert previous['inputs_sha256']==report['inputs_sha256']
        assert all(x['returncode']==0 for x in previous['records'][:2])
        assert not previous.get('completed',False)
        assert all(x['sha256']==sha(x['executable']) for x in previous['records'])
        report['previous_attempt']=previous
        report['reused_initial_conditions_sha256']={f.name:sha(f) for f in snapshot.glob('PMcr*.DAT')}
    try:
        with gzip.open(ROOT/f'simulation-n{args.nrow}.log.gz','at' if args.resume else 'wt') as log:
            for name,cwd,stdin,required in [('PMP2init.exe',case,'','Setup.dat'),
                ('PMP2start.exe',snapshot,'1\n','PMcrd.DAT'),
                ('PMP2main.exe',snapshot,'2000\n',None)]:
                if args.resume and name!='PMP2main.exe':
                    continue
                record=run(BIN/name,cwd,stdin,log,args.threads);records.append(record)
                print(name,record['returncode'],record['elapsed_seconds'],flush=True)
                if record['returncode'] or (required and not (cwd/required).is_file()):
                    raise RuntimeError(f'{name} failed or did not produce {required}')
            report['setup']= (case/'Setup.dat').read_text()
            headers=list(snapshot.glob('PMcrd.*.DAT'))
            assert len(headers)==1, f'Expected exactly one z=0 snapshot: {headers}'
            step=int(headers[0].name.split('.')[1]);report['snapshot_step']=step
            report['snapshot_sha256']={f.name:sha(f) for f in snapshot.glob(f'PMcr*.{step:04d}.DAT')}
            record,data=replay(snapshot,step,args.threads,log)
            records.append(record)
            if data is not None:
                np.savez_compressed(ROOT/f'catalogues-n{args.nrow}.npz',baseline=data)
            report['completed']=True
    finally:
        result_path.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
    main()
