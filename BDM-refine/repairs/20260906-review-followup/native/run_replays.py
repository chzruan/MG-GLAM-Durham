"""Sealed same-particle, same-density BDM v2/control/v3 replays on Slurm.

Run with micromamba run -n cosemu python3 -B. Never overwrites earlier evidence.
The original production membership checker is imported unchanged and its hash
is part of the plan. A completed stage is reusable only after input/output
hashes and its scientific checks have all been verified.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import inspect
import json
import os
from pathlib import Path
import re
import signal
import socket
import struct
import subprocess
import time
import traceback

os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
import numpy as np

ROOT = Path(os.environ.get('BDM_REVIEW_ROOT', Path(__file__).resolve().parent.parent)).resolve()
HERE = ROOT/'native'
REPO = ROOT.parents[2]
VALIDATION = REPO/'BDM-refine/validation/20260906-n1024'


def now():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        while chunk := stream.read(8*1024**2):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix+'.tmp')
    with temporary.open('w') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())
    temporary.replace(path)


def checker():
    source = VALIDATION/'run_validation.py'
    definition = importlib.util.spec_from_file_location('original_validation_checker', source)
    module = importlib.util.module_from_spec(definition)
    definition.loader.exec_module(module)
    module.ROOT = REPO  # Returned relative-path metadata only; check function unchanged.
    return module


def diagnostics(path):
    dtype = np.dtype([('passes','>i4'), ('status','>i4'), ('work','>i8'),
                     ('count','>i8'), ('mass','>f4')])
    with path.open('rb') as stream:
        count, = struct.unpack('>q',stream.read(8))
        rows = np.fromfile(stream, dtype=dtype)
    assert len(rows) == count and path.stat().st_size == 8+28*count
    assert np.all(rows['passes']>=0) and np.all(rows['work']>=0)
    assert np.all(rows['work'][rows['passes']==0]==0)
    assert np.all(rows['work']>=rows['passes'])
    active = rows['passes']>0
    bins = np.bincount(rows['passes'].astype(np.int64))
    selected = rows['mass']>=2.5e12
    long = rows['passes']>32
    top = np.argsort(rows['work'],kind='stable')[-10:][::-1]
    return dict(candidates=count, evaluated_candidates=int(active.sum()),
                max_passes=int(rows['passes'].max(initial=0)),
                active_particle_rows=int(rows['work'].sum()),
                candidates_over_32=int(long.sum()), work_over_32=int(rows['work'][long].sum()),
                selected_over_32=int((long&selected).sum()),
                selected_max_passes=int(rows['passes'][selected].max(initial=0)),
                max_retained_bound_count=int(rows['count'].max(initial=0)),
                pass_histogram={str(i):int(n) for i,n in enumerate(bins) if n},
                active_pass_percentiles={str(p):float(np.percentile(rows['passes'][active],p))
                                         for p in [50,90,95,99,99.9,100]},
                largest_work=[dict(candidate=int(i+1), passes=int(rows['passes'][i]),
                                  active_particle_rows=int(rows['work'][i]),
                                  retained_bound_count=int(rows['count'][i]),
                                  post_selection_mass=float(rows['mass'][i]),status=int(rows['status'][i]))
                              for i in top],
                sha256=sha(path))


def verify_density(path, spec, expected=None, cache=None):
    assert path.stat().st_size == 20+4*spec['ngrid']**3
    with path.open('rb') as stream:
        grid, particles, threads = struct.unpack('>qqi',stream.read(20))
    assert (grid, particles) == (spec['ngrid'], spec['nrow']**3)
    stat=path.stat()
    fingerprint=(stat.st_dev,stat.st_ino,stat.st_size,stat.st_mtime_ns,stat.st_ctime_ns)
    if cache is not None and str(path) in cache:
        # An intermediate metadata guard only: filesystem timestamp resolution
        # cannot prove byte identity. main() hashes every used tape again,
        # without this cache, before declaring the group successful.
        previous_fingerprint,value=cache[str(path)]
        assert fingerprint==previous_fingerprint, f'Density tape changed between calls: {path}'
    else:
        value=sha(path)
        if cache is not None:cache[str(path)]=(fingerprint,value)
    if expected:
        assert value == expected, f'Changed fixed density field: {path}'
    return dict(path=str(path),bytes=path.stat().st_size,sha256=value,density_threads=threads)


def run_stage(case, build, plan, spec, v, threads, runroot, density_cache):
    label = case['label']
    folder = runroot/label
    receipt = folder/'receipt.json'
    binary = Path(build['variants'][case['variant']]['binary_path'])
    expected_binary = build['variants'][case['variant']]['binary_sha256']
    assert sha(binary)==expected_binary
    density = Path(case['density_path'])
    identity = dict(case=case,binary_sha256=expected_binary,threads=threads,
                    plan_sha256=sha(HERE/'plan.json'),driver_sha256=sha(__file__))
    if receipt.exists():
        previous=json.loads(receipt.read_text())
        assert previous['completed'] and all(previous[k]==val for k,val in identity.items()), \
            f'Inspect incomplete or incompatible stage before resuming: {receipt}'
        for path,expected in previous['outputs_sha256'].items():
            assert sha(path)==expected, path
        verify_density(density,spec,previous['density']['sha256'],density_cache)
        return previous
    folder.mkdir(exist_ok=False)
    (folder/'CATALOGS').mkdir()
    (folder/'BDM.config').write_text(spec['halo_config'])
    snapshot=next(s for s in plan['snapshots'] if s['z']==case['z'])
    for name in snapshot['files_sha256']:
        (folder/name).symlink_to(Path(snapshot['directory'])/name)
    command=[str(binary),str(snapshot['header']['step']),str(threads),str(case['passes']),
             str(density),case['density_mode']]
    env=dict(os.environ,LD_LIBRARY_PATH=os.environ['BDM_AUDIT_NATIVE_LIBS'],
             OMP_NUM_THREADS=str(threads),OMP_DYNAMIC='FALSE',OMP_PROC_BIND='close',OMP_PLACES='cores')
    env.pop('LIBRARY_PATH',None)
    record=dict(**identity,completed=False,started_at_utc=now(),command=command,
                job_id=os.environ['SLURM_JOB_ID'],host=socket.gethostname(),
                config_sha256=sha(folder/'BDM.config'),outputs_sha256={})
    write_json(receipt,record)
    process=None
    start=time.monotonic()
    log=folder/'replay.log'
    timing=folder/'replay.time'
    try:
        if case['density_mode']=='read':
            density_before=verify_density(density,spec,case.get('density_sha256'),density_cache)
        else:
            assert not density.exists(), f'Refuse to replace a density realization: {density}'
            density_before=None
        with log.open('xb') as stream:
            process=subprocess.Popen(['/usr/bin/time','-f','%e %U %S %M','-o',str(timing),*command],
                                     cwd=folder,env=env,stdout=stream,stderr=subprocess.STDOUT,start_new_session=True)
            process.wait(timeout=3300)
        record.update(returncode=process.returncode,elapsed_seconds=time.monotonic()-start)
        output=log.read_text()
        assert process.returncode==0 and 'REPLAY COMPLETE' in output, output[-6000:]
        restored=re.findall(r'REPLAY RESTORED pass=\s*(\d+)',output)
        assert list(map(int,restored))==list(range(1,case['passes']+1)),restored
        timings=re.findall(r'REPLAY FINDER pass=\s*(\d+) seconds=\s*([0-9.Ee+-]+)',output)
        assert len(timings)==case['passes'],timings
        record['finder_seconds']=[float(seconds) for _,seconds in timings]
        record['stage_seconds']={key:list(map(float,re.findall(r'time for '+key+r'\s*=\s*([0-9.]+)',output)))
                                 for key in ['AddBuffer','List','ParametersDistinct','WriteFiles']}
        assert all(len(vals)==case['passes'] for vals in record['stage_seconds'].values())
        catalogues=[folder/f'p{i}.DAT' for i in range(1,case['passes']+1)]
        hashes=[sha(path) for path in catalogues]
        assert len(set(hashes))==1,'Repeated calls changed the published catalogue'
        with catalogues[-1].open() as stream:
            header=[next(stream).rstrip('\n') for _ in range(8)]
        version='v3' if case['variant']=='v3' else 'v2'
        assert header[0].endswith(f'[BDM finder {version}]'),header
        catalogue=np.loadtxt(catalogues[-1],skiprows=8,ndmin=2)
        assert len(catalogue)>0 and catalogue.shape[1]==24 and np.isfinite(catalogue).all()
        record['membership']=v.check_memberships(folder/'repair-members.bin',spec,catalogue)
        record['unbinding']=diagnostics(folder/'unbinding.bin')
        assert record['unbinding']['candidates']==record['membership']['candidates']
        record['density']=verify_density(density,spec,density_before['sha256'] if density_before else None,density_cache)
        record['density_verification']='Full hash on first group use; metadata guards between calls; mandatory final group hash'
        record['catalogue']=dict(path=str(catalogues[-1]),sha256=hashes[-1],rows=len(catalogue),version=version)
        record['repeated_catalogue_byte_identical']=True
        record['particle_state_and_workspace_restored']=True
        measured=list(map(float,timing.read_text().split()))
        record.update(native_elapsed_seconds=measured[0],cpu_user_seconds=measured[1],
                      cpu_system_seconds=measured[2],maxrss_kib=int(measured[3]))
        products=[log,timing,folder/'BDM.config',*catalogues,folder/'repair-members.bin',
                  folder/'repair-members.index.npz',folder/'unbinding.bin']
        record['outputs_sha256']={str(path):sha(path) for path in products}
        record.update(completed=True,finished_at_utc=now())
    except BaseException:
        if process and process.poll() is None:
            os.killpg(process.pid,signal.SIGTERM)
            try:process.wait(timeout=15)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid,signal.SIGKILL);process.wait()
        record.update(error=traceback.format_exc(),finished_at_utc=now())
        raise
    finally:
        write_json(receipt,record)
    print('VERIFIED',label,record['catalogue']['rows'],'haloes; max unbinding passes',
          record['unbinding']['max_passes'],flush=True)
    return record


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--threads',type=int,choices=[32,64],required=True)
    parser.add_argument('--plan-sha256',required=True)
    args=parser.parse_args()
    assert os.environ.get('SLURM_JOB_ID'),'N1024 replays require the measured Slurm allocation'
    assert int(os.environ['SLURM_CPUS_PER_TASK'])==args.threads
    assert sha(HERE/'plan.json')==args.plan_sha256
    plan=json.loads((HERE/'plan.json').read_text())
    assert sha(__file__)==plan['driver_sha256']
    assert sha(HERE/'build.json')==plan['build_sha256']
    build=json.loads((HERE/'build.json').read_text());assert build['completed']
    for path,expected in plan['evidence_sha256'].items():assert sha(path)==expected,path
    spec=plan['spec'];assert 0<spec['nrow']**3<1200**3
    v=checker()
    assert hashlib.sha256(inspect.getsource(v.check_memberships).encode()).hexdigest()==plan['checker_function_sha256']
    group=plan['groups'][str(args.threads)]
    runroot=ROOT/f'work/replays-t{args.threads}'
    runroot.mkdir(exist_ok=True)
    result_path=HERE/f'results-t{args.threads}.json'
    assert not result_path.exists(),'Preserve completed group evidence; stage recovery requires a new reviewed launch'
    record=dict(completed=False,started_at_utc=now(),job_id=os.environ['SLURM_JOB_ID'],
                host=socket.gethostname(),threads=args.threads,plan_sha256=args.plan_sha256,
                driver_sha256=sha(__file__),build_sha256=plan['build_sha256'],stages=[])
    snapshots=[s for s in plan['snapshots'] if s['z'] in {c['z'] for c in group}]
    try:
        for snapshot in snapshots:
            for name,expected in snapshot['files_sha256'].items():
                assert sha(Path(snapshot['directory'])/name)==expected,name
        density_cache={}
        for case in group:
            record['stages'].append(run_stage(case,build,plan,spec,v,args.threads,runroot,density_cache))
            write_json(result_path,record)
        by_label={s['case']['label']:s for s in record['stages']}
        before,optimized=by_label['z0-reference'],by_label['z0-optimized-v2']
        assert before['catalogue']['sha256']==optimized['catalogue']['sha256']
        assert before['membership']['raw_sha256']==optimized['membership']['raw_sha256']
        assert before['unbinding']['sha256']==optimized['unbinding']['sha256']
        record['optimization_byte_identical']=True
        for path,(_,expected) in density_cache.items():
            verify_density(Path(path),spec,expected)
        record['input_density_tapes_unchanged']=True
        for snapshot in snapshots:
            for name,expected in snapshot['files_sha256'].items():
                assert sha(Path(snapshot['directory'])/name)==expected,name
        record.update(completed=True,input_snapshots_unchanged=True,finished_at_utc=now())
    except BaseException:
        record.update(error=traceback.format_exc(),finished_at_utc=now())
        raise
    finally:
        write_json(result_path,record)
    print('REPLAY GROUP COMPLETE',args.threads,flush=True)


if __name__=='__main__':
    main()
