"""Freeze once; compile serially on a quiet login node. No simulations here."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import tarfile

from common import BASE, BIN, FINDER_SHA, REPO, ROOT, WORK, file_manifest, git, native_env, now, sha, verify_manifest, write_json
import timetable


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--phase',choices=['common','adapters'],default='common')
    args=parser.parse_args()
    build=WORK/'native-build'
    production=build/'production'
    receipt=build/'common-ready.json'
    BIN.mkdir(parents=True,exist_ok=True)
    assert sha(REPO/'PMP2linker.f90')==FINDER_SHA
    if args.phase=='common':
        if receipt.exists():
            previous=json.loads(receipt.read_text());verify_manifest(previous['frozen_files'])
            print('Verified existing frozen common build');return
        production.mkdir(parents=True,exist_ok=False)
        report=dict(started_at_utc=now(),base=BASE,build_commit=git('rev-parse','HEAD'),commands=[])
        sources=sorted(REPO.glob('*.f90'))+sorted(REPO.glob('*.h'))+[REPO/'makefile']
        for path in sources:
            shutil.copy2(path,production/path.name)
        report['sources']={p.name:sha(p) for p in sources}
        audit=REPO/'BDM-refine/analysis/full-audit-20260906/work-artifacts.tar.gz'
        assert sha(audit)=='b39e39ea17aaf65a53f42a2058986717b89d8ff1513cb02fe65e14294ff9a9c9'
        inputs=WORK/'reference-inputs';inputs.mkdir()
        with tarfile.open(audit,'r:gz') as archive:
            for name in ['Init.dat','PkTable.dat','TableSeeds.dat']:
                with archive.extractfile('reference-inputs/'+name) as src,(inputs/name).open('xb') as dst:
                    shutil.copyfileobj(src,dst)
        assert sha(inputs/'PkTable.dat')=='b90b23d14a403544045ec93533402a70db8b5bc3381a3694547dfd742e0476f2'
        assert sha(inputs/'TableSeeds.dat')=='31e0ae21219dadc93da86c6cfe8cd1972e406f7bbf80c2a1046bafe5ab506e4d'
        legacy=subprocess.check_output(['git','show','a8c7715:PMP2linker.f90'],cwd=REPO)
        (build/'PMP2linker.legacy.f90').write_bytes(legacy)
        for command in [['make','-j1','PMP2mod_tools.o','PMP2mod_MGbackground.o'],
                        ['make','-j1','PMP2init','PMP2start','PMP2main','PMP2BDM']]:
            process=subprocess.run(command,cwd=production,env=native_env(1),capture_output=True,text=True,timeout=1200)
            report['commands'].append(dict(command=command,returncode=process.returncode,stdout=process.stdout,stderr=process.stderr))
            write_json(ROOT/'build.json',report)
            if process.returncode:raise RuntimeError(process.stdout+'\n'+process.stderr)
        for name in ['PMP2init','PMP2start','PMP2main','PMP2BDM']:
            shutil.copy2(production/(name+'.exe'),BIN/(name+'.native.exe'))
        timetable.generate(WORK/'timetable')
        common=[str(p) for p in production.glob('*.o') if p.name not in ['PMP2main.o','PMP2start.o','PMP2bdm.o','PMP2init.o']]
        command=['ifx','-O3','-g','-traceback','-ftz','-unroll','-qopenmp','-march=core-avx2','-mfma',
                 '-fp-model','fast=1','-shared-intel','-mcmodel=medium','-convert','big_endian',
                 str(WORK/'timetable/PMP2main.schedule.f90'),*common,'-o',str(BIN/'PMP2main.schedule.exe')]
        process=subprocess.run(command,cwd=production,env=native_env(1),capture_output=True,text=True,timeout=600)
        report['commands'].append(dict(command=command,returncode=process.returncode,stdout=process.stdout,stderr=process.stderr))
        write_json(ROOT/'build.json',report)
        if process.returncode:raise RuntimeError(process.stdout+'\n'+process.stderr)
        report.update(completed_at_utc=now(),completed=True,
                      frozen_files=file_manifest([*BIN.glob('*.exe'),*inputs.iterdir(),
                                                  build/'PMP2linker.legacy.f90',*production.glob('*.o'),*production.glob('*.mod')]))
        write_json(ROOT/'build.json',report);write_json(receipt,report)
        print('Frozen common objects and controlled main ready',flush=True)
    else:
        raise NotImplementedError('Adapters are integrated after independent preflight')


if __name__=='__main__':main()
