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
        assert receipt.exists(), 'Build the frozen common objects first'
        verify_manifest(json.loads(receipt.read_text())['frozen_files'])
        ic=WORK/'ic-build';replay=WORK/'replay-build'
        commands=[]
        if not (ic/'build.json').exists():
            commands.append(['micromamba','run','-n','cosemu','python3','-B',str(ROOT/'ic/build.py'),
                             '--repo',str(production),'--work',str(ic),'--compiler','ifx'])
        if not (replay/'build.json').exists():
            commands.append(['micromamba','run','-n','cosemu','python3','-B',str(ROOT/'replay/build_adapter.py'),
                             '--common-dir',str(production),'--v3-source',str(production/'PMP2linker.f90'),
                             '--legacy-source',str(build/'PMP2linker.legacy.f90'),'--output-dir',str(replay)])
        for command in commands:
            subprocess.run(command,cwd=REPO,env=native_env(1),check=True)
        ic_report=json.loads((ic/'build.json').read_text())
        replay_report=json.loads((replay/'build.json').read_text())
        assert replay_report['completed']
        sources={ic/'PMP2start.matched.exe':(BIN/'PMP2start.matched.exe',ic_report['binary_sha256'])}
        for variant,item in replay_report['variants'].items():
            sources[Path(item['binary_path'])]=(BIN/f'PMP2replay.{variant}.exe',item['binary_sha256'])
        for source,(destination,expected) in sources.items():
            assert sha(source)==expected
            if destination.exists():assert sha(destination)==expected
            else:shutil.copy2(source,destination)
        shutil.copy2(ic/'build.json',ROOT/'ic-build.json')
        shutil.copy2(replay/'build.json',ROOT/'replay-build.json')
        write_json(ROOT/'executables.json',dict(completed=True,prepared_at_utc=now(),
                    binaries=file_manifest(sorted(BIN.glob('*.exe'))),
                    production_finder_sha256=FINDER_SHA,
                    build_receipts=file_manifest([ROOT/'build.json',ROOT/'ic-build.json',ROOT/'replay-build.json'])))
        print('Verified frozen campaign executables and adapter receipts')


if __name__=='__main__':main()
