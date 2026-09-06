"""Verify corrected semiaxis labels using actual GetHalo and its shape block.

Run with micromamba run -n cosemu python3 -B axis_order_regression.py --compiler gfortran.
All compilers/processes use small controlled inputs and temporary build products.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import tempfile

os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['MKL_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import numpy as np

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[3]
BASE='b6e96669b6f7a8cba2f3dcaa9853278aac7edf11'


def extract(source,name,kind='subroutine'):
    pattern=rf'^\s*(?:pure\s+)?(?:real\*8\s+)?{kind}\s+{name}\b.*?^\s*end\s+{kind}\s+{name}\b[^\n]*'
    result=re.search(pattern,source,re.I|re.M|re.S)
    assert result is not None,name
    return result.group()


def driver(source):
    get_halo=extract(source,'GetHalo')
    start=get_halo.index('              axis_ratio=sqrt(')
    finish=get_halo.index('Zax(ip)=direction(3)',start)+len('Zax(ip)=direction(3)')
    correction=get_halo[start:finish]
    return '''program axis_order_cases
use LinkerList
implicit none
character(20) :: which
integer :: n,j
integer*8 :: ip
real :: axis(3),direction(3),concentration,centre(3)
real*8 :: aperture,axis_ratio(3),concentration_proxy,slope_b,slope_c
call get_command_argument(1,which)
open(13,status='scratch')
if(trim(which)=='correction')then
  allocate(RadRms(1),Axba(1),Axca(1),Xax(1),Yax(1),Zax(1))
  read(*,*)n
  ip=1_8;aperture=1.d0
  open(50,file='corrections.bin',access='stream',form='unformatted',status='replace')
  do j=1,n
    read(*,*)axis,direction,concentration
    RadRms(1)=concentration
''' + correction + '''
    write(50)Axba(1),Axca(1),Xax(1),Yax(1),Zax(1)
  enddo
  close(50)
  stop
endif
read(*,*)n,MassOne
Np=n;Nparticles=n;Nmaxima=1;Box=32.;NGRID=128;Cell=1.
Om0=.3;Ovdens=200.;AEXPN=.8;Rext=0.;SlopeR=.2;dLogR=.02
centre=5.
Nmx=-2;Nmy=-2;Nmz=-2;Nbx=34;Nby=34;Nbz=34
allocate(Xpar(n),Ypar(n),Zpar(n),VX(n),VY(n),VZ(n),OriginalParticleId(n))
do j=1,n
  read(*,*)Xpar(j),Ypar(j),Zpar(j),VX(j),VY(j),VZ(j)
  OriginalParticleId(j)=j
enddo
allocate(Mvir(1),Rvir(1),Mtotal(1),Xoff(1),xMaxx(1),yMaxx(1),zMaxx(1), &
 VxMaxx(1),VyMaxx(1),VzMaxx(1),EpotM(1),EkinM(1),LambdaM(1),RadRms(1), &
 VmaxM(1),RmaxM(1),Xax(1),Yax(1),Zax(1),Axba(1),Axca(1))
allocate(Lst(Np),Label(Nmx:Nbx,Nmy:Nby,Nmz:Nbz))
call List
xMaxx=centre(1);yMaxx=centre(2);zMaxx=centre(3)
call BdmHaloMembershipInit
call GetHalo(centre(1),centre(2),centre(3),0.,0.,0.,1_8)
open(50,file='halo.bin',access='stream',form='unformatted',status='replace')
write(50)HaloStatus(1)
write(50)Mvir(1),Mtotal(1),Rvir(1),EkinM(1),EpotM(1),VmaxM(1),RmaxM(1), &
 VxMaxx(1),VyMaxx(1),VzMaxx(1),RadRms(1),Xoff(1),LambdaM(1),Axba(1),Axca(1), &
 Xax(1),Yax(1),Zax(1),xMaxx(1),yMaxx(1),zMaxx(1)
if(allocated(BoundParticleIds(1)%ids))then
  write(50)size(BoundParticleIds(1)%ids,kind=8),BoundParticleIds(1)%ids
else
  write(50)0_8
endif
close(50)
end program axis_order_cases
'''


def build(work,compiler):
    current=(REPO/'PMP2linker.f90').read_text()
    before=subprocess.check_output(['git','show',f'{BASE}:PMP2linker.f90'],cwd=REPO,text=True)
    compilation={};hashes={}
    for variant,source in [('before',before),('repaired',current)]:
        hashes[variant]=hashlib.sha256(source.encode()).hexdigest()
        structures=re.search(r'^module\s+Structures\b.*?^end module\s+Structures',source,re.I|re.M|re.S).group()
        stub='''module Tools
real :: Box=32.,AEXPN=.8
integer :: NGRID=128,NROW=32
integer*8 :: Nparticles=0
real,allocatable :: Xpar(:),Ypar(:),Zpar(:),VX(:),VY(:),VZ(:)
contains
real function seconds()
seconds=0.
end function
real function Memory(n)
integer*8 :: n
Memory=0.
end function
end module
module LinkerList
use Structures
use Tools
contains
'''
        names=['GetHalo','BdmHaloGather','BdmHaloSortRadii','BdmHaloSortIds','BdmHaloSphericalPotential',
               'BdmHaloMembershipInit','List','Limits','EigenValues','BdmParticlePosition']
        generated=structures+'\n'+stub+'\n'.join(extract(source,n) for n in names)
        generated+='\n'+extract(source,'BdmParticleCoordinate','function')+'\nend module\n'+driver(source)
        filename=f'{variant}.f90';(work/filename).write_text(generated)
        for mode in ['checked','optimized']:
            if compiler=='gfortran':
                flags=['-O0','-g','-fcheck=all','-ffpe-trap=invalid,zero,overflow'] if mode=='checked' else ['-O3','-fno-fast-math']
                common=['-fopenmp','-ffree-line-length-none']
            else:
                flags=['-O0','-check','bounds','-fpe0','-fp-model','precise'] if mode=='checked' else ['-O3','-fp-model','precise']
                common=['-qopenmp','-extend-source']
            binary=f'{variant}-{mode}'
            command=[compiler,*flags,*common,filename,'-o',binary]
            env=os.environ.copy()
            if compiler=='ifx':
                env['LD_LIBRARY_PATH']=env['BDM_AUDIT_NATIVE_LIBS'];env.pop('LIBRARY_PATH',None)
            result=subprocess.run(command,cwd=work,env=env,capture_output=True,text=True)
            assert result.returncode==0,result.stderr
            compilation[binary]=dict(command=command,returncode=result.returncode,stdout=result.stdout,stderr=result.stderr)
    return hashes,compilation


def run(work,compiler,variant,mode,which,stdin):
    with tempfile.TemporaryDirectory(prefix='bdm-axis-case-') as directory:
        directory=Path(directory);env={**os.environ,'OMP_NUM_THREADS':'1','OMP_DYNAMIC':'FALSE'}
        if compiler=='ifx':env['LD_LIBRARY_PATH']=env['BDM_AUDIT_NATIVE_LIBS']
        result=subprocess.run([str(work/f'{variant}-{mode}'),which],input=stdin,cwd=directory,
                              env=env,capture_output=True,text=True,timeout=20)
        assert result.returncode==0,(which,result.stdout,result.stderr)
        filename='corrections.bin' if which=='correction' else 'halo.bin'
        return (directory/filename).read_bytes()


def correction_input():
    direction=np.asarray([1.,2.,-3.]);direction/=np.linalg.norm(direction)
    entries=[(1.,.75**2,.749**2,*direction,.65)]
    for b in [0.,.1,.3,.5,.7,.75,.9,1.]:
        for fraction in [0.,.5,.99,.999,1.]:
            for concentration in [0.,.4,.400001,.6,.65,.8,1.]:
                entries.append((1.,b*b,(b*fraction)**2,*direction,concentration))
    values=np.asarray(entries,dtype=np.float32)
    text=str(len(values))+'\n'+'\n'.join(' '.join(f'{float(v):.17g}' for v in row) for row in values)+'\n'
    return values,text


def halo_inputs():
    n=512;q=np.arange(n,dtype=float);z=1.-2.*(q+.5)/n;phi=q*2.399963229728653
    sphere=np.c_[np.sqrt(1-z*z)*np.cos(phi),np.sqrt(1-z*z)*np.sin(phi),z]
    angle=.47;c=math.cos(angle);s=math.sin(angle)
    rotation=np.asarray([[c,-s,0.],[s,c,0.],[0.,0.,1.]])
    cases=[]
    for name,axes,rso,rotate in [
        ('near_degenerate_transverse',[1.,.75,.749],1.05,False),
        ('near_degenerate_transverse_rotated',[1.,.75,.749],1.05,True),
        ('nearly_spherical',[1.,.999,.998],1.05,False),
        ('separated_transverse_axes',[1.,.75,.3],1.05,True),
        ('below_empirical_threshold',[1.,.75,.749],4.,False),
        ('rank_one',[1.,0.,0.],1.05,True)]:
        position=sphere*np.asarray(axes)
        if rotate:position=position@rotation.T
        phase=np.c_[position+5.,np.tile([10000.,-5000.,2000.],(n,1))].astype(np.float32)
        mass=np.float32(1.150e12*float(np.float32(.3))*200.*rso**3/n)
        text=f'{n} {float(mass):.17g}\n'+'\n'.join(' '.join(f'{float(v):.17g}' for v in row) for row in phase)+'\n'
        cases.append((name,phase,text))
    return cases


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--compiler',choices=['gfortran','ifx'],default='gfortran')
    parser.add_argument('--output',type=Path);args=parser.parse_args();records=[];properties=[]
    with tempfile.TemporaryDirectory(prefix='bdm-axis-build-') as directory:
        work=Path(directory);hashes,compilation=build(work,args.compiler)
        entries,stdin=correction_input()
        for mode in ['checked','optimized']:
            before=np.frombuffer(run(work,args.compiler,'before',mode,'correction',stdin),dtype=np.float32).reshape(-1,5)
            after=np.frombuffer(run(work,args.compiler,'repaired',mode,'correction',stdin),dtype=np.float32).reshape(-1,5)
            assert before[0,1]>before[0,0],'explicit regression fixture must reproduce the inversion'
            assert np.array_equal(after[:,:2],np.sort(before[:,:2],axis=1)[:,::-1])
            assert np.array_equal(after[:,2:],before[:,2:])
            assert np.array_equal(after[:,2:],entries[:,3:6])
            assert np.all(np.isfinite(after)) and np.all((0<=after[:,1])&(after[:,1]<=after[:,0])&(after[:,0]<=1))
            properties.append(dict(mode=mode,cases=len(entries),before_inversions=int(np.sum(before[:,1]>before[:,0])),
                                   explicit_before=before[0].astype(float).tolist(),explicit_after=after[0].astype(float).tolist(),
                                   checks=['corrected lengths preserved bitwise as an unordered pair','0 <= c/a <= b/a <= 1',
                                           'major-axis direction preserved bitwise']))
            for name,phase,text in halo_inputs():
                old=run(work,args.compiler,'before',mode,'halo',text);new=run(work,args.compiler,'repaired',mode,'halo',text)
                old_values=np.frombuffer(old[4:88],dtype=np.float32);new_values=np.frombuffer(new[4:88],dtype=np.float32)
                assert old[:4]==new[:4] and old[88:]==new[88:],'status or membership changed'
                unchanged=[i for i in range(21) if i not in [13,14]]
                assert np.array_equal(old_values[unchanged],new_values[unchanged]),name
                assert np.array_equal(new_values[13:15],np.sort(old_values[13:15])[::-1]),name
                assert 0<=new_values[14]<=new_values[13]<=1,name
                nbound=int(np.frombuffer(new[88:96],dtype=np.int64)[0]);assert nbound==len(phase)
                if name=='near_degenerate_transverse':assert old_values[14]>old_values[13]
                offset=phase[:,:3].astype(float)-5.;r2=np.sum(offset*offset,axis=1)
                tensor=sum(np.outer(row,row)/distance2 for row,distance2 in zip(offset,r2))/len(offset)
                eigen=np.linalg.eigvalsh(tensor);direction=new_values[15:18].astype(float)
                assert abs(np.linalg.norm(direction)-1.)<2.e-7
                assert np.linalg.norm(tensor@direction-eigen[-1]*direction)<2.e-7
                records.append(dict(mode=mode,case=name,particles=len(phase),bound_count=nbound,
                                    before_ratios=old_values[13:15].astype(float).tolist(),after_ratios=new_values[13:15].astype(float).tolist(),
                                    major_direction=direction.tolist(),member_sha256=hashlib.sha256(new[96:]).hexdigest(),
                                    checks=['only the two corrected semiaxis labels can change','identical masses, dynamics, positions and memberships',
                                            'unit major direction remains a principal tensor eigenvector']))
            print(args.compiler,mode,len(entries),'shape-block cases and',len(halo_inputs()),'actual GetHalo cases passed',flush=True)
    report=dict(checked_at_utc=datetime.now(timezone.utc).isoformat(),compiler=args.compiler,baseline_commit=BASE,
                source_sha256=hashes['repaired'],baseline_source_sha256=hashes['before'],compilation=compilation,
                test_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                property_cases=properties,halo_cases=records,passed=sum(p['cases'] for p in properties)+len(records),
                interpretation='Both empirical corrected lengths and all non-shape results are preserved; only intermediate/minor labels are ordered')
    if args.output:args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(report['passed'],'axis-order regressions passed')


if __name__=='__main__':main()
