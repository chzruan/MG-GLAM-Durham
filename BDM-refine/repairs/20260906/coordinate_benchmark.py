"""Pair native N128 replays to assess the original-row coordinate shortcut.

Inputs must be copied frozen native-build objects/modules, baseline/member
binaries and the verified snapshot, as recorded in copied-manifest.json.
Initialize Intel modules and BDM_AUDIT_NATIVE_LIBS, then run with cosemu Python.
All timed programs are pinned to one CPU. No simulation or Slurm job is launched.
"""
import argparse
from datetime import datetime,timezone
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import re
import signal
import statistics
import subprocess
import tempfile
import time

ROOT=Path(__file__).resolve().parent
REPO=ROOT.parents[2]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def extract(source,name,kind):
    pattern=rf'^\s*(?:pure\s+)?(?:real\*8\s+)?{kind}\s+{name}\b.*?^\s*end\s+{kind}\s+{name}\b[^\n]*'
    found=re.search(pattern,source,re.I|re.M|re.S)
    assert found,name
    return pattern,found.group()


def cpu_snapshot():
    values={}
    for line in Path('/proc/stat').read_text().splitlines():
        fields=line.split()
        if fields[0].startswith('cpu') and fields[0][3:].isdigit():
            ticks=list(map(int,fields[1:]));values[int(fields[0][3:])]=(sum(ticks[:8]),ticks[3]+ticks[4])
    return values


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--workspace',type=Path,required=True)
    parser.add_argument('--pairs',type=int,default=4)
    parser.add_argument('--output',type=Path,default=ROOT/'coordinate-results.json')
    args=parser.parse_args();assert args.pairs>=3
    workspace=args.workspace.resolve();build=workspace/'build'
    manifest=json.loads((workspace/'copied-manifest.json').read_text())
    for name,digest in manifest.items():assert sha(workspace/name)==digest,name
    source=(REPO/'PMP2linker.f90').read_text()
    baseline=(build/'PMP2linker.f90').read_text()
    env={**os.environ,'OMP_NUM_THREADS':'1','OMP_DYNAMIC':'FALSE','OMP_PROC_BIND':'close','OMP_PLACES':'cores',
         'OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1','LD_LIBRARY_PATH':os.environ['BDM_AUDIT_NATIVE_LIBS']}
    env.pop('LIBRARY_PATH',None)
    report=dict(started_at_utc=datetime.now(timezone.utc).isoformat(),
        source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=REPO,text=True).strip(),
        production_source_sha256=hashlib.sha256(source.encode()).hexdigest(),
        frozen_source_sha256=sha(build/'PMP2linker.f90'),driver_sha256=sha(__file__),
        copied_manifest=manifest,commands=[],runs=[],completed=False)
    def command(arguments):
        p=subprocess.run(arguments,cwd=build,env=env,text=True,capture_output=True)
        report['commands'].append(dict(command=arguments,returncode=p.returncode,stdout=p.stdout,stderr=p.stderr))
        assert p.returncode==0,(p.stdout,p.stderr)
    command(['ifx','--version'])
    objects=['PMP2mod_tools','PMP2mod_fft5','PMP2mod_random','PMP2mod_density','PMP2mod_power',
        'PMP2mod_analyze','PMP2MG_subroutines','PMP2mod_MGbackground','PMP2MGsolver_fR',
        'PMP2extradof','PMP2MGsolver_DGP','PMP2MGsolver_sym','PMP2MGsolver_kmf','PMP2MGsolver_csf']
    common=['-g','-traceback','-qopenmp','-march=core-avx2','-mfma','-shared-intel','-mcmodel=medium','-convert','big_endian']
    optimized=['-O3','-fp-model','precise','-ftz','-unroll']
    for variant,filename,entry,link_flags in [
        ('optimized','PMP2linker.f90','PMP2bdm.o',['-O3','-fp-model','fast=1','-ftz','-unroll']),
        ('members-optimized','finder-members.f90','entry-members.o',optimized)]:
        updated=(build/filename).read_text()
        for name,kind in [('BdmParticleCoordinate','function'),('BdmParticlePosition','subroutine')]:
            pattern,replacement=extract(source,name,kind)
            updated,count=re.subn(pattern,lambda _:replacement,updated,flags=re.I|re.M|re.S)
            assert count==1
        path=build/f'coordinate-{variant}.f90';path.write_text(updated)
        command(['ifx',*common,*optimized,'-c',path.name,'-o',f'coordinate-{variant}.o'])
        command(['ifx',*common,*link_flags,'-o',f'PMP2BDM.{variant}.exe',
                 *[name+'.o' for name in objects],f'coordinate-{variant}.o',entry])
    first=cpu_snapshot();time.sleep(.3);second=cpu_snapshot()
    utilization={cpu:1-(second[cpu][1]-first[cpu][1])/max(1,second[cpu][0]-first[cpu][0])
                 for cpu in os.sched_getaffinity(0)}
    cpu=min(utilization,key=lambda c:(utilization[c],-c))
    report['cpu_affinity']=cpu;report['cpu_busy_fraction_before']=utilization[cpu]
    golden=None;golden_members=None
    log_path=workspace/'coordinate-replays.log.gz'
    with gzip.open(log_path,'wt') as log:
        def replay(variant,pair=None):
            nonlocal golden,golden_members
            name={'baseline':'PMP2BDM.exe','members-baseline':'PMP2BDM.members.exe'}.get(variant,f'PMP2BDM.{variant}.exe')
            binary=build/name
            with tempfile.TemporaryDirectory(prefix='coordinate-run-',dir=workspace) as tmp:
                work=Path(tmp);(work/'CATALOGS').mkdir()
                for p in (workspace/'snapshot').iterdir():(work/p.name).symlink_to(p)
                (work/'BDM.config').write_text('iVirial=1\nMassMin=2.5e12\nRext=0.15\n')
                cmd=['taskset','-c',str(cpu),'/usr/bin/time','-f','%e %U %S %M','-o',str(work/'time.txt'),str(binary)]
                start=time.monotonic()
                p=subprocess.Popen(cmd,cwd=work,env=env,text=True,stdin=subprocess.PIPE,
                    stdout=subprocess.PIPE,stderr=subprocess.PIPE,start_new_session=True)
                try:stdout,stderr=p.communicate('157\n',timeout=180)
                except subprocess.TimeoutExpired:
                    os.killpg(p.pid,signal.SIGKILL);stdout,stderr=p.communicate()
                    log.write(stdout+'\n'+stderr);raise
                elapsed=time.monotonic()-start
                log.write(f'\n=== {variant} pair={pair} ===\n'+stdout+'\n'+stderr)
                assert p.returncode==0,(stdout,stderr)
                wall,user,system,rss=map(float,(work/'time.txt').read_text().splitlines()[-1].split())
                files=list((work/'CATALOGS').glob('Catshort*.DAT'));assert len(files)==1
                raw=files[0].read_bytes();rows=raw.decode().splitlines()[8:]
                assert rows and all(len(row.split())==24 and all(math.isfinite(float(x)) for x in row.split()) for row in rows)
                if golden is None:golden=raw
                assert raw==golden,f'{variant} changed full catalogue bytes'
                record=dict(variant=variant,pair=pair,wall_seconds=wall,elapsed_seconds=elapsed,
                    user_seconds=user,system_seconds=system,cpu_seconds=user+system,maxrss_kib=int(rss),
                    binary_sha256=sha(binary),catalogue_sha256=hashlib.sha256(raw).hexdigest(),rows=len(rows))
                if variant.startswith('members-'):
                    members=(work/'repair-members.bin').read_bytes()
                    if golden_members is None:golden_members=members
                    assert members==golden_members,'original-ID memberships or raw properties changed'
                    record['raw_membership_sha256']=hashlib.sha256(members).hexdigest()
                report['runs'].append(record)
                print(variant,pair,round(elapsed,3),record['cpu_seconds'],len(rows),flush=True)
        # Warm the executable and page cache before balanced AB/BA timing.
        replay('baseline',-1);replay('optimized',-1)
        for pair in range(args.pairs):
            order=['baseline','optimized'] if pair%2==0 else ['optimized','baseline']
            for variant in order:replay(variant,pair)
        replay('members-baseline');replay('members-optimized')
    timings={variant:[r for r in report['runs'] if r['variant']==variant and r['pair'] is not None and r['pair']>=0]
             for variant in ['baseline','optimized']}
    report['medians']={variant:{key:statistics.median(r[key] for r in rows)
                              for key in ['wall_seconds','elapsed_seconds','cpu_seconds']}
                       for variant,rows in timings.items()}
    report['paired_speedups']=[next(r['elapsed_seconds'] for r in timings['baseline'] if r['pair']==pair)/
                              next(r['elapsed_seconds'] for r in timings['optimized'] if r['pair']==pair)
                              for pair in range(args.pairs)]
    report['median_elapsed_speedup']=report['medians']['baseline']['elapsed_seconds']/report['medians']['optimized']['elapsed_seconds']
    report['replay_log_sha256']=sha(log_path);report['completed']=True
    args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report['medians'],indent=2),flush=True)
    print('Median elapsed speedup',report['median_elapsed_speedup'],flush=True)


if __name__=='__main__':main()
