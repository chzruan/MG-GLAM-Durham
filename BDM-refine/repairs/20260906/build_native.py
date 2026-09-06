"""Freeze source and build the full native finder plus read-only audit variants."""
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess

ROOT=Path(__file__).resolve().parent
REPO=ROOT.parents[2]
BUILD=ROOT/'work/native-build'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


MEMBERSHIP='''
subroutine RepairMembershipDump
implicit none
integer :: unit,ip,nselected
logical :: selected
nselected=0
do ip=1,Nmaxima
  if(Mvir(ip)<MassMin.or.Mvir(ip)<=0.)cycle
  if(xMaxx(ip)<0..or.xMaxx(ip)>=Box.or.yMaxx(ip)<0..or.yMaxx(ip)>=Box.or.zMaxx(ip)<0..or.zMaxx(ip)>=Box)cycle
  nselected=nselected+1
enddo
open(newunit=unit,file='repair-members.bin',access='stream',form='unformatted',status='replace')
write(unit)int(nselected,8),int(Nmaxima,8),MassOne
do ip=1,Nmaxima
  if(Mvir(ip)<MassMin.or.Mvir(ip)<=0.)cycle
  if(xMaxx(ip)<0..or.xMaxx(ip)>=Box.or.yMaxx(ip)<0..or.yMaxx(ip)>=Box.or.zMaxx(ip)<0..or.zMaxx(ip)>=Box)cycle
  if(.not.allocated(BoundParticleIds(ip)%ids))error stop 'repair membership missing'
  write(unit)int(ip,8),size(BoundParticleIds(ip)%ids,kind=8)
  write(unit)xMaxx(ip),yMaxx(ip),zMaxx(ip),VxMaxx(ip),VyMaxx(ip),VzMaxx(ip), &
    Mvir(ip),Mtotal(ip),Rvir(ip),EkinM(ip),EpotM(ip),VmaxM(ip),RmaxM(ip), &
    Xoff(ip),LambdaM(ip),RadRms(ip),Axba(ip),Axca(ip),Xax(ip),Yax(ip),Zax(ip)
  write(unit)BoundParticleIds(ip)%ids
enddo
close(unit)
end subroutine RepairMembershipDump
'''

PROBE='''program BdmStateProbe
use LinkerList
use Density
implicit none
integer :: step,pass
integer*8 :: original_count
integer*4,allocatable :: saved(:,:)
read(*,*)step
call ReadDataPM(step,'')
original_count=Nparticles
allocate(saved(Nparticles,6))
saved(:,1)=transfer(Xpar,saved(:,1));saved(:,2)=transfer(Ypar,saved(:,2))
saved(:,3)=transfer(Zpar,saved(:,3));saved(:,4)=transfer(VX,saved(:,4))
saved(:,5)=transfer(VY,saved(:,5));saved(:,6)=transfer(VZ,saved(:,6))
allocate(FI(NGRID,NGRID,NGRID))
do pass=1,2
  call BDM(1)
  if(Nparticles/=original_count.or.Np/=original_count)error stop 'BDM changed particle count'
  if(size(Xpar,kind=8)/=original_count.or.size(Ypar,kind=8)/=original_count.or. &
     size(Zpar,kind=8)/=original_count.or.size(VX,kind=8)/=original_count.or. &
     size(VY,kind=8)/=original_count.or.size(VZ,kind=8)/=original_count) &
     error stop 'BDM changed particle array sizes'
  if(any(saved(:,1)/=transfer(Xpar,saved(:,1))))error stop 'BDM changed X bits'
  if(any(saved(:,2)/=transfer(Ypar,saved(:,2))))error stop 'BDM changed Y bits'
  if(any(saved(:,3)/=transfer(Zpar,saved(:,3))))error stop 'BDM changed Z bits'
  if(any(saved(:,4)/=transfer(VX,saved(:,4))))error stop 'BDM changed VX bits'
  if(any(saved(:,5)/=transfer(VY,saved(:,5))))error stop 'BDM changed VY bits'
  if(any(saved(:,6)/=transfer(VZ,saved(:,6))))error stop 'BDM changed VZ bits'
  if(.not.allocated(FI))error stop 'BDM lost density allocation'
  if(allocated(OriginalParticleId).or.allocated(BoundParticleIds).or.allocated(HaloStatus)) &
    error stop 'BDM retained per-call membership buffers'
enddo
print *, 'BDM NATIVE STATE RESTORATION: TWO CALLS BITWISE IDENTICAL'
end program BdmStateProbe
'''


def main():
    BUILD.mkdir(parents=True,exist_ok=False)
    sources=sorted(REPO.glob('*.f90'))+sorted(REPO.glob('*.h'))+[REPO/'makefile']
    for path in sources:shutil.copy2(path,BUILD/path.name)
    report=dict(started_at_utc=datetime.now(timezone.utc).isoformat(),
        source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=REPO,text=True).strip(),
        source_sha256={p.name:sha(BUILD/p.name) for p in sources},
        source_worktree_status=subprocess.check_output(['git','status','--porcelain','--untracked-files=no'],cwd=REPO,text=True),
        builder_sha256=sha(__file__),state_probe_sha256=hashlib.sha256(PROBE.encode()).hexdigest(),
        membership_probe_sha256=hashlib.sha256(MEMBERSHIP.encode()).hexdigest(),commands=[],variants={},completed=False)
    env=dict(os.environ);env['LD_LIBRARY_PATH']=os.environ['BDM_AUDIT_NATIVE_LIBS'];env.pop('LIBRARY_PATH',None)
    def run(command):
        p=subprocess.run(command,cwd=BUILD,env=env,text=True,capture_output=True)
        report['commands'].append(dict(command=command,returncode=p.returncode,stdout=p.stdout,stderr=p.stderr))
        if p.returncode:raise RuntimeError(p.stdout+'\n'+p.stderr)
    try:
        run(['ifx','--version'])
        run(['make','-f','makefile','-j1','PMP2mod_tools.o','PMP2mod_MGbackground.o'])
        run(['make','-f','makefile','-j1','PMP2init','PMP2start','PMP2main','PMP2BDM'])
        original=(BUILD/'PMP2linker.f90').read_text()
        assert original.count('Call WriteFiles\n')==2
        instrumented=original.replace('Call WriteFiles\n','Call WriteFiles\n      Call RepairMembershipDump\n')
        instrumented=re.sub(r'end module linkerlist',MEMBERSHIP+'\nend Module LinkerList',instrumented,flags=re.I)
        objects=['PMP2mod_tools','PMP2mod_fft5','PMP2mod_random','PMP2mod_density','PMP2mod_power',
                 'PMP2mod_analyze','PMP2MG_subroutines','PMP2mod_MGbackground','PMP2MGsolver_fR',
                 'PMP2extradof','PMP2MGsolver_DGP','PMP2MGsolver_sym','PMP2MGsolver_kmf','PMP2MGsolver_csf']
        common=['-g','-traceback','-qopenmp','-march=core-avx2','-shared-intel','-mcmodel=medium','-convert','big_endian']
        optimized=['-O3','-fp-model','precise','-ftz','-unroll','-mfma']
        for variant,source,flags in [('members',instrumented,optimized),
                                     ('bounds',original,['-O0','-fp-model','precise','-check','bounds','-fpe0'])]:
            path=BUILD/f'finder-{variant}.f90';path.write_text(source)
            run(['ifx',*common,*flags,'-c',path.name,'-o',f'finder-{variant}.o'])
            run(['ifx',*common,*flags,'-c','PMP2bdm.f90','-o',f'entry-{variant}.o'])
            run(['ifx',*common,*flags,'-o',f'PMP2BDM.{variant}.exe',
                 *[o+'.o' for o in objects],f'finder-{variant}.o',f'entry-{variant}.o'])
            report['variants'][variant]=dict(source_sha256=sha(path),binary_sha256=sha(BUILD/f'PMP2BDM.{variant}.exe'))
        (BUILD/'native-state-probe.f90').write_text(PROBE)
        run(['ifx',*common,*optimized,'-c','native-state-probe.f90','-o','native-state-probe.o'])
        run(['ifx',*common,*optimized,'-o','PMP2BDM.state.exe',
             *[o+'.o' for o in objects],'PMP2linker.o','native-state-probe.o'])
        report['binaries_sha256']={p.name:sha(p) for p in BUILD.glob('*.exe')}
        report['completed']=True
    finally:
        (ROOT/'native-build.json').write_text(json.dumps(report,indent=2)+'\n')
    print('Native build complete:',len(report['binaries_sha256']),'executables.')


if __name__=='__main__':main()
