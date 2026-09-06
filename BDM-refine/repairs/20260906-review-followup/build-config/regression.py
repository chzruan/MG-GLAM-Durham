"""Bounded N1/N5/N6 compiler/runtime regressions; no simulation allocations.

Use micromamba run -n cosemu python3 -B regression.py. Load Intel 2024.2.0
modules and capture BDM_AUDIT_NATIVE_LIBS before entering the Python environment.
All compilation/run scratch is temporary; existing test receipts are untouched.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import resource
import subprocess
import tempfile

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
BASE = '54f53aa'
INPUT = '1.00000011920928955078125 0.99999988079071044921875 -1\n'
FP_MODULE = '''module fp_probe
implicit none
contains
real function muladd(a,b,c)
real, intent(in) :: a,b,c
muladd=a*b+c
end function
end module
'''
FP_MAIN = '''program probe
use fp_probe
implicit none
real :: a,b,c,value
read(*,*) a,b,c
value=muladd(a,b,c)
write(*,'(z8.8,es24.15)')transfer(value,0),value
end program
'''
CONFIG_MAIN = '''program config_cases
use LinkerList
implicit none
character(32) :: name
integer :: flag,unit
integer :: saved(2,6)
call get_command_argument(1,name)
select case(trim(name))
case('config')
  MaxMemory=17.
  call ReadParameters(1)
  print *, 'MAXMEMORY=',MaxMemory
  call PublishCatalogue
case('reset')
  call ReadParameters(1)
  if(MaxMemory/=1.25)error stop 'configured memory was not read'
  call PublishCatalogue
  open(newunit=unit,file='BDM.config',status='replace')
  write(unit,'(a)')'iVirial=1'
  close(unit)
  call ReadParameters(2)
  if(MaxMemory/=500.)error stop 'memory setting leaked between calls'
  call PublishCatalogue
case('size_list')
  call ReadParameters(1)
  ! Only the count is large. SizeList computes geometry/bytes; no particle,
  ! mesh, list, or halo arrays are allocated in this admission-gate test.
  Np=100000000_8
  call SizeList
  if(Cell/=0.25)error stop 'memory policy changed physical Cell'
  print *, 'MAXMEMORY=',MaxMemory,' CELL=',Cell
  call PublishCatalogue
case('GetProfiles')
  call GetProfiles
case('WriteProfiles')
  call WriteProfiles
case('HaloProfile')
  NradP=huge(NradP)
  call HaloProfile(0.,0.,0.,0.,0_8)
case('RemoveDuplicatesSimple')
  call RemoveDuplicatesSimple
case('rescale_zero','rescale_negative','rescale_two')
  flag=0
  if(name=='rescale_negative')flag=-1
  if(name=='rescale_two')flag=2
  call RescaleCoords(flag)
case('rescale_active')
  Nparticles=2_8;Np=2_8
  allocate(Xpar(2),Ypar(2),Zpar(2),VX(2),VY(2),VZ(2))
  Xpar=[1.,9.];Ypar=[3.,7.];Zpar=[9.,13.]
  VX=[-0.,17.3];VY=[1.,-3.];VZ=[10.,-4.]
  saved(:,1)=transfer(Xpar,saved(:,1));saved(:,2)=transfer(Ypar,saved(:,2))
  saved(:,3)=transfer(Zpar,saved(:,3));saved(:,4)=transfer(VX,saved(:,4))
  saved(:,5)=transfer(VY,saved(:,5));saved(:,6)=transfer(VZ,saved(:,6))
  call RescaleCoords(1)
  if(any(Xpar/=[0.,1.]).or.any(Ypar/=[0.25,0.75]))error stop 'active rescaling changed'
  call RemoveBuffer(2_8)
  if(any(saved(:,1)/=transfer(Xpar,saved(:,1))).or.any(saved(:,2)/=transfer(Ypar,saved(:,2))).or. &
     any(saved(:,3)/=transfer(Zpar,saved(:,3))).or.any(saved(:,4)/=transfer(VX,saved(:,4))).or. &
     any(saved(:,5)/=transfer(VY,saved(:,5))).or.any(saved(:,6)/=transfer(VZ,saved(:,6)))) &
     error stop 'active roundtrip did not preserve particle bits'
  if(allocated(BdmPMX).or.BdmPMCount/=0_8.or.Np/=2_8)error stop 'retained workspace'
case default
  error stop 'unknown regression case'
end select
print *, 'BUILD CONFIG TEST PASSED'
end program
'''


def sha(value):
    if isinstance(value, Path):
        value = value.read_bytes()
    if isinstance(value, str):
        value = value.encode()
    return hashlib.sha256(value).hexdigest()


def run(command, cwd, stdin='', success=True):
    env = {**os.environ, 'LD_LIBRARY_PATH':os.environ['BDM_AUDIT_NATIVE_LIBS'],
           'OMP_NUM_THREADS':'1', 'OMP_DYNAMIC':'FALSE', 'OMP_PROC_BIND':'FALSE',
           'OPENBLAS_NUM_THREADS':'1', 'MKL_NUM_THREADS':'1'}
    env.pop('LIBRARY_PATH', None)
    process = subprocess.run(command, cwd=cwd, input=stdin, env=env, text=True,
                             capture_output=True, timeout=90)
    record = dict(command=command, stdin=stdin, returncode=process.returncode,
                  stdout=process.stdout, stderr=process.stderr)
    if success:
        assert process.returncode == 0, record
    else:
        assert process.returncode != 0, record
    return record


def extract(source, name):
    pattern = rf'^\s*subroutine\s+{name}\b.*?^\s*end\s+subroutine\s+{name}\b[^\n]*'
    match = re.search(pattern, source, re.I|re.M|re.S)
    assert match, name
    return match.group()


def config_source(source):
    structures = re.search(r'^module\s+Structures\b.*?^end module\s+Structures', source, re.I|re.M|re.S).group()
    tools = '''module Tools
implicit none
real :: Box=2.,AEXPN=1.,Om=.3,OmL=.7,ASTEP=.004,hubble=.7
integer :: NGRID=16,NROW=8,ISTEP=1,Nrealization=1
integer*8 :: Nparticles=0,memoryWords=0
real,allocatable :: Xpar(:),Ypar(:),Zpar(:),VX(:),VY(:),VZ(:)
contains
real function Memory(n)
integer*8, intent(in) :: n
memoryWords=memoryWords+n
Memory=0.25+real(dble(memoryWords)*4.d0/1024.d0**3)
end function
end module
module LinkerList
use Structures
use Tools
contains
'''
    names = ['ReadParameters','ValidateParameters','ConfigurationError','BeginCataloguePublication',
             'PublishCatalogue','PrepareParticleSearch','SizeList','RescaleCoords','RemoveBuffer',
             'GetProfiles','WriteProfiles','HaloProfile','RemoveDuplicatesSimple']
    return structures+'\n'+tools+'\n'.join(extract(source,name) for name in names)+'\nend module\n'+CONFIG_MAIN


def config_tests(source, scratch):
    generated = config_source(source)
    (scratch/'config.f90').write_text(generated)
    cases = [(None,500.), ('mAxMeMoRy = 1.25 ! GiB\n',1.25), ('MAXMEMORY=512.5\n',512.5),
             ('MaxMemory=2.88D2\n',288.), ('MaxMemory=0.5\n',0.5)]
    invalid = ['0','-1','NaN','Inf','-Inf','1.e100','garbage','1 2','1,2','']
    records = []
    for compiler, mode, flags in [('gfortran','checked',['-O0','-g','-fcheck=all','-fopenmp']),
                                  ('gfortran','optimized',['-O3','-fopenmp','-fno-fast-math','-ffp-contract=off']),
                                  ('ifx','checked',['-O0','-g','-check','bounds','-check','pointers','-qopenmp','-fp-model','precise']),
                                  ('ifx','optimized',['-O3','-qopenmp','-march=core-avx2','-mfma','-fp-model','precise'])]:
        binary = scratch/f'{compiler}-{mode}'
        compilation = run([compiler,*flags,str(scratch/'config.f90'),'-o',str(binary)],scratch)
        results = []
        def execute(name, config=None, success=True, expected=None, message=None):
            with tempfile.TemporaryDirectory(prefix='bdm-config-case-') as tmp:
                cwd = Path(tmp);(cwd/'CATALOGS').mkdir()
                if config is not None:
                    (cwd/'BDM.config').write_text(config)
                sentinels = [cwd/'CATALOGS/CatshortV.0001.0001.DAT',cwd/'CATALOGS/outputB.0001.dat']
                if not success and name == 'config':
                    for path in sentinels:path.write_text('prior valid output\n')
                result = run([str(binary),name],cwd,success=success)
                result.update(case=name,config=config)
                if success:
                    assert 'BUILD CONFIG TEST PASSED' in result['stdout'],result
                    if expected is not None:
                        actual=float(re.search(r'MAXMEMORY=\s*(\S+)',result['stdout'])[1])
                        assert actual==expected,result
                    if name=='config' and config is None:
                        generated_config=(cwd/'BDM.config').read_text()
                        assert 'MaxMemory' in generated_config and 'GiB' in generated_config
                        result['generated_config']=generated_config
                else:
                    assert message in result['stdout']+result['stderr'],result
                    if name=='config':
                        assert sorted(p.name for p in (cwd/'CATALOGS').iterdir())==sorted(p.name for p in sentinels)
                        assert all(path.read_text()=='prior valid output\n' for path in sentinels)
                        result['existing_outputs_untouched']=True
                    if name.startswith('rescale_') or name in ['GetProfiles','WriteProfiles','HaloProfile','RemoveDuplicatesSimple']:
                        assert list((cwd/'CATALOGS').iterdir())==[]
                results.append(result)
        for config,expected in cases:execute('config',config,expected=expected)
        execute('reset','MaxMemory=1.25\n')
        for value in invalid:execute('config','MaxMemory='+value+'\n',success=False,message='BDM configuration error')
        execute('size_list','MaxMemory=1\n',expected=1.)
        execute('size_list','MaxMemory=0.5\n',success=False,message='configured memory limit')
        for name in ['GetProfiles','WriteProfiles','HaloProfile','RemoveDuplicatesSimple']:
            execute(name,success=False,message=f'BDM {name} is unsupported')
        for name in ['rescale_zero','rescale_negative','rescale_two']:
            execute(name,success=False,message='BDM RescaleCoords only supports iFlag=1')
        execute('rescale_active')
        records.append(dict(compiler=compiler,mode=mode,compilation=compilation,binary_sha256=sha(binary),results=results))
    return dict(extracted_source_sha256=sha(generated),variants=records,passed=sum(len(v['results']) for v in records))


def make_fixture(path, makefile):
    (path/'makefile').write_text(makefile)
    objects = re.search(r'^OBJ\s*=\s*(.*)$',makefile,re.M)[1].split()
    for obj in objects:
        name=Path(obj).stem
        (path/(name+'.f90')).write_text(f'subroutine stub_{name}\nend subroutine\n')
    (path/'PMP2linker.f90').write_text(FP_MODULE)
    (path/'PMP2main.f90').write_text(FP_MAIN)


def fp_tests(scratch):
    records=[]
    for label,flags,hex_value in [
            ('fast',['-fp-model','fast=1'],'A8800000'),
            ('precise',['-fp-model','precise'],'00000000'),
            ('fast_then_precise',['-fp-model','fast=1','-fp-model','precise'],'00000000'),
            ('precise_then_fast',['-fp-model','precise','-fp-model','fast=1'],'A8800000'),
            ('ofast_then_precise',['-Ofast','-fp-model','precise'],'00000000'),
            ('fastmath_then_precise',['-ffast-math','-fp-model','precise'],'00000000')]:
        path=scratch/label;path.mkdir()
        (path/'probe.f90').write_text(FP_MODULE+FP_MAIN)
        compilation=run(['ifx','-O3','-march=core-avx2','-mfma',*flags,'probe.f90','-o','probe'],path)
        result=run([str(path/'probe')],path,INPUT)
        assert result['stdout'].split()[0]==hex_value,result
        records.append(dict(label=label,compilation=compilation,runtime=result,binary_sha256=sha(path/'probe')))
    original=subprocess.check_output(['git','show',BASE+':makefile'],cwd=REPO,text=True)
    current=(REPO/'makefile').read_text()
    builds=[]
    for label,makefile,target,overrides,expected in [
            ('old_override_reproducer',original,'PMP2main',['FFLAGS=-O3 -march=core-avx2 -mfma'],'A8800000'),
            ('default',current,'PMP2main',[],'00000000'),
            ('override_without_model',current,'PMP2main',['FFLAGS=-O2 -g -qopenmp -march=core-avx2 -mfma'],'00000000'),
            ('override_conflicting_model',current,'PMP2main',['FFLAGS=-O3 -march=core-avx2 -mfma -fp-model fast=2'],'00000000'),
            ('override_bdm_flags',current,'PMP2main',['BDM_FFLAGS=-O3 -march=core-avx2 -mfma -fp-model fast=1'],'00000000'),
            ('bitmatch',current,'PMP2main-bitmatch',[],'00000000'),
            ('bitmatch_fma_sensitive',current,'PMP2main-bitmatch',
             ['FFLAGS_BITMATCH=-O2 -g -qopenmp -march=core-avx2 -mfma'],'00000000'),
            ('gnu_override',current,'PMP2main',
             ['FC=gfortran','FFLAGS=-O3 -march=core-avx2 -mfma -ffast-math','LDFLAGS=-fopenmp'],'00000000'),
            ('gnu_compiler_option',current,'PMP2main',
             ['FC=gfortran -m64','FFLAGS=-O3 -march=core-avx2 -mfma -ffast-math','LDFLAGS=-fopenmp'],'00000000')]:
        path=scratch/('make-'+label);path.mkdir();make_fixture(path,makefile)
        compilation=run(['make','-j1',*overrides,target],path)
        result=run([str(path/(target+'.exe'))],path,INPUT)
        assert result['stdout'].split()[0]==expected,result
        finder_lines=[line for line in compilation['stdout'].splitlines() if line.endswith('-c PMP2linker.f90')]
        assert len(finder_lines)==1,compilation
        if label!='old_override_reproducer':
            suffix='-fno-fast-math -ffp-contract=off -w -c PMP2linker.f90' if label.startswith('gnu_') else '-fp-model precise -w -c PMP2linker.f90'
            assert finder_lines[0].endswith(suffix),finder_lines
        builds.append(dict(label=label,makefile_sha256=sha(makefile),compilation=compilation,runtime=result,
                           binary_sha256=sha(path/(target+'.exe')),finder_command=finder_lines[0]))
    return dict(probe_source_sha256=sha(FP_MODULE+FP_MAIN),direct_compiler_controls=records,actual_make_builds=builds)


def call_graph():
    names='RescaleCoords|GetProfiles|HaloProfile|WriteProfiles|RemoveDuplicatesSimple'
    answer=dict(before=[],after=[],source_sha256={})
    for path in sorted(REPO.glob('*.f90')):
        old=subprocess.check_output(['git','show',BASE+':'+path.name],cwd=REPO,text=True)
        current=path.read_text()
        answer['source_sha256'][path.name]=sha(current)
        for stage,text in [('before',old),('after',current)]:
            for number,line in enumerate(text.splitlines(),1):
                if re.search(rf'\bcall\s+({names})\b',line.split('!')[0],re.I):
                    answer[stage].append(dict(file=path.name,line=number,code=line.strip()))
    assert len(answer['after'])==1 and re.search(r'call\s+RescaleCoords\(1\)',answer['after'][0]['code'],re.I)
    return answer


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,default=HERE/'results.json')
    args=parser.parse_args()
    resource.setrlimit(resource.RLIMIT_CORE,(0,0))
    source=(REPO/'PMP2linker.f90').read_text()
    report=dict(started_at_utc=datetime.now(timezone.utc).isoformat(),completed=False,base_commit=BASE,
                source_sha256={name:sha(REPO/name) for name in ['PMP2linker.f90','makefile']},
                driver_sha256=sha(Path(__file__)),native_runtime_libraries=os.environ['BDM_AUDIT_NATIVE_LIBS'])
    report['production_call_graph']=call_graph()
    with tempfile.TemporaryDirectory(prefix='bdm-build-config-') as tmp:
        scratch=Path(tmp)
        report['compiler_versions']=[run([compiler,'--version'],scratch) for compiler in ['ifx','gfortran']]
        report['configuration_and_entrypoints']=config_tests(source,scratch)
        report['floating_point_flags']=fp_tests(scratch)
    report.update(completed=True,finished_at_utc=datetime.now(timezone.utc).isoformat())
    args.output.write_text(json.dumps(report,indent=2)+'\n')
    print('Passed',report['configuration_and_entrypoints']['passed'],'configuration/entrypoint cases, '
          '6 compiler precedence controls and 9 actual make/runtime builds.')


if __name__=='__main__':main()
