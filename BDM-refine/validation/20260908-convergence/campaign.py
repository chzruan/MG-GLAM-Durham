"""Resumable stage runner. Every large native invocation requires Slurm.

Use micromamba run -n cosemu python3 -B campaign.py --help.
No incomplete process is promoted to a scientifically complete simulation.
"""
import argparse
import json
import os
from pathlib import Path
import re
import shutil
import signal
import socket
import subprocess
import sys
import time
import traceback

from common import BIN, CONFIG, MATRIX, REPO, ROOT, WORK, checker, file_manifest, native_env, now, sha, verify_manifest, write_json


def spec_for(name):
    if name in MATRIX:
        nrow,ngrid=MATRIX[name]
        return dict(case=name,nrow=nrow,ngrid=ngrid,box_mpc_h=256.,epochs=[2,1,0],
                    half_steps=name=='T',cosmology=dict(Omega_m=.3089,Omega_Lambda=.6911,h=.6774,sigma8=.8159),
                    shape_repair_commit='c68c22d',analysis_ngrid=2048)
    raise ValueError(name)


def init_text(spec,save=True):
    schedule=json.loads((ROOT/'timetable.json').read_text())
    da=schedule['half_initial_step' if spec.get('half_steps') else 'normal_initial_step']
    changes={'Box':spec['box_mpc_h'],'Nrow':spec['nrow'],'Ngrid':spec['ngrid'],
             'step da':format(da,'.16e'),'#outputs':-3,'Steps between checkpoints':10000,
             'Save snapshots':int(save),'DM power spectrum':0,'Find BDM halos':0,
             'MG_flag':0,'MG_test':0,'MG_model':3}
    lines=(WORK/'reference-inputs/Init.dat').read_text().splitlines()
    result=[];i=0
    while i<len(lines):
        key=lines[i].split('=')[0].strip()
        if key=='#outputs':
            old=abs(int(lines[i].split('=')[1].split()[0]))
            result.extend(['#outputs = -3',' 2.0',' 1.0',' 0.0']);i+=old+1;continue
        result.append(f'{key} = {changes[key]}' if key in changes else lines[i]);i+=1
    return '\n'.join(result)+'\n'


def fixed_text(path,text):
    if path.exists():
        if path.read_text()!=text:raise ValueError(f'Attempt to change frozen input: {path}')
    else:path.write_text(text)


def link(path,target):
    target=target.resolve()
    if path.is_symlink():
        if path.resolve()!=target:raise ValueError(f'Changed link: {path}')
    elif path.exists():raise FileExistsError(path)
    else:path.symlink_to(target)


def initial_particles(run):
    return sorted(p for p in run.glob('PMcrs*.DAT') if re.fullmatch(r'PMcrs\d+\.DAT',p.name))


def setup_directory(case,spec,save=True):
    case.mkdir(parents=True,exist_ok=True)
    run=case/'Run1';run.mkdir(exist_ok=True);(run/'CATALOGS').mkdir(exist_ok=True)
    fixed_text(case/'Init.dat',init_text(spec,save))
    fixed_text(run/'BDM.config',CONFIG)
    table=WORK/'timetable'/('half-schedule.dat' if spec.get('half_steps') else 'native-schedule.dat')
    fixed_text(run/'campaign_schedule.dat',table.read_text())
    for name in ['PkTable.dat','TableSeeds.dat']:link(case/name,WORK/'reference-inputs'/name)
    return run


def run_native(binary,cwd,tag,threads,inputs,outputs,stdin='',arguments=(),allow_login=False,finalize=None):
    if not allow_login and 'SLURM_JOB_ID' not in os.environ:
        raise RuntimeError('Large campaign stages must run inside Slurm')
    binary=Path(binary).resolve()
    receipt=cwd/(tag+'.json');log=cwd/(tag+'.log');timing=cwd/(tag+'.time')
    identity=dict(binary_sha256=sha(binary),threads=threads,stdin=stdin,arguments=list(map(str,arguments)),
                  inputs=file_manifest(inputs))
    if receipt.exists():
        previous=json.loads(receipt.read_text())
        if previous.get('completed') and previous['identity']==identity:
            verify_manifest(previous['outputs']);verify_manifest(previous['evidence'])
            print(f'Verified completed {cwd.name}/{tag}',flush=True);return previous
        raise RuntimeError(f'Incomplete/incompatible stage requires inspection: {receipt}')
    if log.exists() or timing.exists():raise FileExistsError(f'Unreceipted artifacts: {cwd}/{tag}')
    report=dict(completed=False,identity=identity,binary=str(binary),cwd=str(cwd),started_at_utc=now(),
                job_id=os.environ.get('SLURM_JOB_ID'),host=socket.gethostname(),driver_sha256=sha(sys.argv[0]))
    write_json(receipt,report)
    start=time.monotonic();process=None
    try:
        with log.open('xb') as stream:
            process=subprocess.Popen(['/usr/bin/time','-f','%e %U %S %M','-o',str(timing),str(binary),*map(str,arguments)],
                                     cwd=cwd,env=native_env(threads),stdin=subprocess.PIPE,
                                     stdout=stream,stderr=subprocess.STDOUT,start_new_session=True)
            process.communicate(stdin.encode())
        report.update(returncode=process.returncode,elapsed_seconds=time.monotonic()-start,finished_at_utc=now())
        values=timing.read_text().splitlines()[-1].split()
        report.update(cpu_user_seconds=float(values[1]),cpu_system_seconds=float(values[2]),maxrss_kib=int(values[3]))
        if process.returncode:raise RuntimeError(f'Native failure: {log}')
        if finalize is not None:finalize()
        product=list(outputs())
        if not product:raise RuntimeError(f'No expected output: {log}')
        report.update(outputs=file_manifest(product),evidence=file_manifest([log,timing]),completed=True)
    except BaseException:
        if process is not None and process.poll() is None:
            os.killpg(process.pid,signal.SIGKILL);process.wait()
        report['failure_traceback']=traceback.format_exc()
        raise
    finally:write_json(receipt,report)
    print(f'{tag}: {report["elapsed_seconds"]:.1f} s; RSS {report["maxrss_kib"]/1024**2:.2f} GiB',flush=True)
    return report


def run_init(case,spec,allow_login=False):
    def precise_step():
        # Native es12.5 output truncates T's exact first half-step. Preserve it
        # and publish the round-trip value actually used by the IC velocities.
        original=(case/'Setup.dat').read_text()
        (case/'Setup.native.dat').write_text(original)
        lines=original.splitlines()
        table=json.loads((ROOT/'timetable.json').read_text())
        da=table['half_initial_step' if spec.get('half_steps') else 'normal_initial_step']
        lines[3]=f'{da:.16e}  Step in dAEXPN (campaign round-trip value)'
        (case/'Setup.dat').write_text('\n'.join(lines)+'\n')
    return run_native(BIN/'PMP2init.native.exe',case,'init',1,[case/'Init.dat',case/'PkTable.dat'],
                      lambda:[case/'Setup.dat',case/'Setup.native.dat'],allow_login=allow_login,finalize=precise_step)


def initialize(case,spec,threads,allow_login=False):
    run=setup_directory(case,spec)
    run_init(case,spec,allow_login)
    return run


def initial_conditions(name,threads):
    spec=spec_for(name);case=WORK/'cases'/name
    run=initialize(case,spec,threads)
    master=WORK/'cases/E/Run1'
    from ic.configure import configure
    configure(run,None if name=='E' else master,master_nrow=1024,origin_ngrid=2048)
    inputs=[case/'Setup.dat',case/'PkTable.dat',case/'TableSeeds.dat',run/'matched_ic.nml',run/'matched_ic_inputs.json']
    if name!='E':inputs.extend([master/'matched_ic_receipt.txt',master/'matched_ic_inputs.json'])
    receipt=run_native(BIN/'PMP2start.matched.exe',run,'ic',threads,inputs,
                       lambda:[run/'PMcrd.DAT',*initial_particles(run),
                               run/'matched_ic_receipt.txt',run/'matched_modes.bin'],stdin='1\n')
    header=checker().read_header(run/'PMcrd.DAT')
    assert header['nrow']==spec['nrow'] and header['ngrid']==spec['ngrid']
    assert header['particles']==spec['nrow']**3<1200**3
    assert abs(header['scale_factor']-1/101)<1.e-8
    assert sum(p.stat().st_size for p in initial_particles(run))==24*spec['nrow']**3
    write_json(ROOT/f'{name}-ic.json',dict(completed=True,spec=spec,header=header,stage_receipt=str(run/'ic.json'),
                                        outputs=receipt['outputs']))


def evolve(name,threads,pilot_steps=0):
    spec=spec_for(name);source=WORK/'cases'/name
    ic_receipt=json.loads((ROOT/f'{name}-ic.json').read_text())
    verify_manifest(ic_receipt['outputs'])
    if pilot_steps:
        case=WORK/'pilots'/f'{name}-t{threads}-s{pilot_steps}'
        run=setup_directory(case,spec,save=False)
        run_init(case,spec)
        for p in (source/'Run1').glob('PMcr*.DAT'):
            if p.name=='PMcrd.DAT' or re.fullmatch(r'PMcrs\d+\.DAT',p.name):link(run/p.name,p)
    else:case=source;run=case/'Run1'
    inputs=[case/'Setup.dat',case/'TableSeeds.dat',run/'campaign_schedule.dat',run/'PMcrd.DAT',
            *initial_particles(run)]
    products=(lambda:[run/'timing.log',run/'Run.log']) if pilot_steps else (
        lambda:[*sorted(run.glob('PMcr*.[0-9][0-9][0-9][0-9].DAT')),run/'timing.log',run/'Run.log'])
    stage=run_native(BIN/'PMP2main.schedule.exe',run,'evolve',threads,inputs,products,
                     stdin=f'{pilot_steps or 1000}\n')
    log=(run/'evolve.log').read_text()
    reached=[int(v) for v in re.findall(r'(?:Step\s*=|STEP=)\s*(\d+)',log)]
    expected=pilot_steps or (316 if spec['half_steps'] else 158)
    assert reached and max(reached)==expected, (reached[-5:],expected)
    snapshots=[]
    if not pilot_steps:
        headers=sorted(run.glob('PMcrd.*.DAT'));assert len(headers)==3
        for z,path in zip(spec['epochs'],headers):
            header=checker().read_header(path)
            assert header['particles']==spec['nrow']**3 and header['ngrid']==spec['ngrid']
            assert abs(header['scale_factor']-1/(1+z))<2.e-7
            data=sorted(run.glob(f'PMcrs*.{header["step"]:04d}.DAT'))
            assert sum(p.stat().st_size for p in data)==24*spec['nrow']**3
            snapshots.append(dict(redshift=z,header=header,header_path=str(path),data_paths=list(map(str,data))))
    report=dict(completed=True,spec=spec,pilot_steps=pilot_steps,threads=threads,snapshots=snapshots,
                stage_receipt=str(run/'evolve.json'),outputs=stage['outputs'],completed_at_utc=now())
    write_json(ROOT/(f'pilot-{name}-t{threads}-s{pilot_steps}.json' if pilot_steps else f'{name}-simulation.json'),report)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('phase',choices=['ic','evolve'])
    parser.add_argument('case',choices=list(MATRIX))
    parser.add_argument('--threads',type=int,default=int(os.environ.get('SLURM_CPUS_PER_TASK','1')))
    parser.add_argument('--pilot-steps',type=int,default=0)
    args=parser.parse_args()
    if args.phase=='ic':initial_conditions(args.case,args.threads)
    else:evolve(args.case,args.threads,args.pilot_steps)


if __name__=='__main__':main()
