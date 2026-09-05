"""Compile actual fixed routines: periodicity, transitive groups, host-mask race."""
from pathlib import Path
import json
import os
import re
import subprocess
import sys

from source_reproduction import extract

root=Path(__file__).resolve().parents[4]
source=(root/'PMP2linker.f90').read_text()
build=Path(__file__).resolve().parent/'cases'
build.mkdir(exist_ok=True)
parts=[]
for name in ['BdmDuplicateRules','Structures']:
    parts.append(re.search(rf'^module\s+{name}\b.*?^end module\s+{name}',source,re.I|re.M|re.S).group())
parts.append('''module Tools
real :: Box=32.
contains
real function seconds()
call cpu_time(seconds)
end function
end module
module LinkerList
use Structures
use Tools
use BdmDuplicateRules
contains
''')
parts.extend(extract(source,name) for name in ['RemoveDuplicates','MergeNumericalDuplicates','Limits','ListMaxima'])
parts.append('''end module
program test
use LinkerList
implicit none
integer :: i,ip
logical :: expected(9)
Nmaxima=9;MassOne=1.e10;Cell=3.5
Nmx=-1;Nmy=-1;Nmz=-1;Nbx=11;Nby=11;Nbz=11
allocate(Mvir(9),Rvir(9),Mtotal(9),xMaxx(9),yMaxx(9),zMaxx(9))
allocate(VxMaxx(9),VyMaxx(9),VzMaxx(9),Lst(9),Label(-1:11,-1:11,-1:11))
Mvir=1.e13;Rvir=.5;Mtotal=1.e13
xMaxx=[.01,31.99,10.30,10.,10.15,20.,20.05,25.,25.05]
yMaxx=8.;zMaxx=8.;VxMaxx=0.;VyMaxx=0.;VzMaxx=0.
VxMaxx(7)=200.;Mtotal(9)=1.1e13
call ListMaxima
call RemoveDuplicates
expected=[.true.,.false.,.true.,.false.,.false.,.true.,.true.,.true.,.true.]
if(any((Mvir>0.).neqv.expected))error stop 'periodic/component/merger failure'
! Unequal-mass chain: B removes A and C removes B. Decisions must see original B.
Mvir(4:9)=0.
Mvir(1:3)=[1.e12,2.e12,3.e12];Rvir(1:3)=.2
xMaxx(1:3)=[1.,1.15,1.30];Mtotal(1:3)=Mvir(1:3)
call ListMaxima
call RemoveDuplicates
if(Mvir(1)/=0..or.Mvir(2)/=0..or.Mvir(3)/=3.e12)error stop 'host-mask race'
print *, 'PERIODIC, CHAIN, MERGER CONTROLS, READ-ONLY HOST MASK: PASS'
end program
''')
(build/'finder_cases.f90').write_text('\n'.join(parts))
cmd=['gfortran','-O0','-g','-fopenmp','-fcheck=all','-ffree-line-length-none','finder_cases.f90','-o','finder_cases']
subprocess.run(cmd,cwd=build,check=True)
results=[]
for n in [1,2,4]:
    p=subprocess.run(['./finder_cases'],cwd=build,env={**os.environ,'OMP_NUM_THREADS':str(n)},text=True,capture_output=True,check=True)
    print(n,p.stdout)
    results.append(dict(threads=n,stdout=p.stdout,stderr=p.stderr))
(build/'finder_cases_result.json').write_text(json.dumps(results,indent=2)+'\n')
