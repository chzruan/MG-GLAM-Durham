"""Verify catalogue replacement after successful write/close, including failures.

Run with micromamba run -n cosemu python3 -B publication_regression.py.
All staged files and builds are temporary. Tests use extracted finder routines.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import stat
import subprocess
import tempfile
import time

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]

DRIVER = '''program publication_cases
use LinkerList
use, intrinsic :: ieee_arithmetic
implicit none
integer :: j,tag,proceed
character(32) :: which,argument,ready_file
call get_command_argument(1,which)
call get_command_argument(2,argument)
tag=1
if(len_trim(argument)>0)read(argument,*)tag
call ReadParameters(1)
call SetParameters
if(trim(which)=='abort_header')error stop 'injected failure after header'
if(trim(which)=='closed_stage')then
  close(12)
  call PublishCatalogue
  error stop 'closed stage was published'
endif
Nmaxima=2;MassMin=0.
allocate(Mvir(2),Mtotal(2),Rvir(2),xMaxx(2),yMaxx(2),zMaxx(2),VxMaxx(2),VyMaxx(2),VzMaxx(2))
allocate(EkinM(2),EpotM(2),VmaxM(2),RmaxM(2),Xoff(2),LambdaM(2),RadRms(2))
allocate(Axba(2),Axca(2),Xax(2),Yax(2),Zax(2),BoundParticleIds(2))
allocate(BoundParticleIds(1)%ids(20),BoundParticleIds(2)%ids(20))
BoundParticleIds(1)%ids=[(int(j,8),j=1,20)];BoundParticleIds(2)%ids=BoundParticleIds(1)%ids+20
Mvir=20.*MassOne;Mtotal=Mvir;Rvir=.2;xMaxx=real(tag);yMaxx=5.;zMaxx=5.
VxMaxx=0.;VyMaxx=0.;VzMaxx=0.;EkinM=1.e14;EpotM=1.e15
VmaxM=0.;RmaxM=0.;Xoff=0.;LambdaM=0.;RadRms=0.;Axba=0.;Axca=0.
Xax=1.;Yax=0.;Zax=0.
if(trim(which)=='invalid_late')EkinM(2)=ieee_value(1.,ieee_quiet_nan)
if(trim(which)=='empty')Nmaxima=0
if(trim(which)=='waiting')then
  flush(12)
  write(ready_file,'(a,i0)')'ready.',tag
  open(44,file=trim(ready_file),status='new')
  write(44,'(a)')trim(CatalogueStagedPath)
  close(44)
  read(*,*)proceed
endif
call WriteFiles
if(CataloguePublicationPending)error stop 'successful publication left pending state'
print *, 'PUBLICATION COMPLETE'
end program
'''


def source_program(source):
    def routine(name):
        pattern = (rf'^\s*(?:real\s+)?(?:subroutine|function)\s+{name}\b.*?'
                   rf'^\s*end\s+(?:subroutine|function)\s+{name}\b[^\n]*')
        return re.search(pattern, source, re.I|re.M|re.S).group()
    module = re.search(r'^module\s+Structures\b.*?^end module Structures',source,re.I|re.M|re.S).group()
    tools = '''
module Tools
real :: Box=32.,AEXPN=1.,Om=.3,OmL=.7,ASTEP=.004,hubble=.7
integer :: NGRID=16,NROW=8,ISTEP=1,Nrealization=1
integer*8 :: Nparticles=40
character*45 :: HEADER='Publication regression'
end module
module LinkerList
use Structures
use Tools
contains
'''
    names = ['ReadParameters','ConfigurationError','ValidateParameters','SetOverdensity','SetParameters',
             'WriteFiles','Concentration','BeginCataloguePublication','PublishCatalogue']
    return module+tools+'\n'.join(routine(name) for name in names)+'\nend module\n'+DRIVER


def prepare(work, prior, executable):
    (work/'CATALOGS').mkdir()
    (work/'BDM.config').write_text('iVirial=2\nMassMin=0\n')
    catalogue = work/'CATALOGS/CatshortM.0001.0001.DAT'
    old = None
    if prior:
        p = subprocess.run([str(executable),'valid'],cwd=work,capture_output=True,text=True,timeout=10)
        assert p.returncode == 0,(p.stdout,p.stderr)
        old = catalogue.read_bytes()
        assert len(old.splitlines()) == 10
    return catalogue, old


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,default=HERE.parent/'publication-results.json')
    args = parser.parse_args()
    source = (REPO/'PMP2linker.f90').read_text()
    results = []
    compilation = []
    with tempfile.TemporaryDirectory(prefix='bdm-publication-build-') as scratch:
        build = Path(scratch)
        (build/'cases.f90').write_text(source_program(source))
        for mode,flags in [('checked',['-O0','-g','-fcheck=all','-ffpe-trap=invalid,zero,overflow']),
                           ('optimized',['-O3'])]:
            command = ['gfortran',*flags,'-ffree-line-length-none','cases.f90','-o',mode]
            p = subprocess.run(command,cwd=build,capture_output=True,text=True)
            assert p.returncode == 0,p.stderr
            compilation.append(dict(mode=mode,command=command,stderr=p.stderr))
            for case in ['valid','empty','invalid_late','abort_header','closed_stage','rename_failure']:
                for prior in [False,True]:
                    with tempfile.TemporaryDirectory(prefix='bdm-publication-case-') as temp:
                        work = Path(temp);catalogue,old = prepare(work,prior,build/mode)
                        protected = None
                        if case == 'rename_failure':
                            if catalogue.exists(): catalogue.unlink()
                            catalogue.mkdir();protected=catalogue/'preserved';protected.write_text('protected')
                        p = subprocess.run([str(build/mode),case],cwd=work,capture_output=True,
                                           text=True,timeout=10,umask=0o027)
                        success = case in ['valid','empty']
                        assert (p.returncode == 0) == success,(case,p.stdout,p.stderr)
                        staged = list((work/'CATALOGS').glob('.*.tmp.*'))
                        if success:
                            rows = catalogue.read_text().splitlines()
                            assert len(rows) == (8 if case == 'empty' else 10),rows
                            assert all(len(row.split()) == 24 for row in rows[8:]),rows
                            assert not staged,staged
                            assert stat.S_IMODE(catalogue.stat().st_mode) == 0o640
                        else:
                            assert len(staged) == 1,staged
                            if protected is not None:
                                assert catalogue.is_dir() and protected.read_text() == 'protected'
                            elif prior:
                                assert catalogue.read_bytes() == old,'prior catalogue changed on failure'
                            else:
                                assert not catalogue.exists(),'failure exposed a new final catalogue'
                        results.append(dict(mode=mode,case=case,prior_catalogue=prior,
                            returncode=p.returncode,staged_count=len(staged),stderr=p.stderr,passed=True))
            # Two active writers must stage distinct paths. Readers keep seeing
            # the prior final file until each complete replacement is published.
            with tempfile.TemporaryDirectory(prefix='bdm-publication-concurrent-') as temp:
                work = Path(temp);catalogue,old = prepare(work,True,build/mode)
                children = []
                try:
                    for tag in [1,2]:
                        children.append(subprocess.Popen([str(build/mode),'waiting',str(tag)],cwd=work,
                            stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True))
                    deadline=time.monotonic()+10
                    while not all((work/f'ready.{tag}').exists() and
                                  (work/f'ready.{tag}').stat().st_size>0 for tag in [1,2]):
                        assert all(p.poll() is None for p in children),'waiting writer failed'
                        assert time.monotonic()<deadline,'waiting writer timed out'
                        time.sleep(.01)
                    names=[(work/f'ready.{tag}').read_text().strip() for tag in [1,2]]
                    assert names[0] != names[1],names
                    assert catalogue.read_bytes() == old,'staging changed visible catalogue'
                    for tag,p in enumerate(children,1):
                        stdout,stderr=p.communicate('1\n',timeout=10)
                        assert p.returncode == 0,(stdout,stderr)
                        rows=catalogue.read_text().splitlines()
                        assert len(rows) == 10 and float(rows[8].split()[0]) == tag,rows
                    assert not list((work/'CATALOGS').glob('.*.tmp.*'))
                    results.append(dict(mode=mode,case='concurrent_writers',unique_staged_paths=True,passed=True))
                finally:
                    for p in children:
                        if p.poll() is None:p.kill()
                        p.communicate()
    report=dict(validated_at_utc=datetime.now(timezone.utc).isoformat(),
        source_sha256=hashlib.sha256(source.encode()).hexdigest(),
        driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        asserted_experiments=len(results),compilation=compilation,results=results)
    args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS: {len(results)} asserted catalogue-publication experiments (checked and optimized).')


if __name__ == '__main__':
    main()
