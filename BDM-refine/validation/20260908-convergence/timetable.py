"""Generate native and subdivided schedules without altering production source."""
import argparse
import json
from pathlib import Path
import subprocess

import numpy as np

from common import REPO, ROOT, WORK, native_env, sha, write_json


def native_generator(source):
    local = source[source.index('Module LocalData'):source.index('end Module LocalData') + len('end Module LocalData')]
    schedule = source[source.index('    Alist(:) = 0.'):source.index("    write (*, *) 'Step  a_expansion")]
    exact = source[source.index('Subroutine SetExactSteps'):source.index('end Subroutine SetExactSteps') + len('end Subroutine SetExactSteps')]
    return '''module Tools
  integer :: Nexact=3,Nout=0,ISTEP=0
  real :: ASTEP0=0.0004,AEXPN0,StepFactor,AEXPN,ASTEP
  real :: zexact(1000),zout(1000)
end module
''' + local + '''
program ExportNativeSchedule
  use Tools
  use LocalData
  AEXPN0=1./101.
  AEXPN=AEXPN0;ASTEP=ASTEP0;StepFactor=ASTEP0/AEXPN0
  zexact(1:3)=[2.,1.,0.]
''' + schedule + '''
  open(71,file='native-schedule.dat',status='replace')
  write(71,*)NlastX
  do i=1,NlastX
    write(71,'(i6,2es25.16,2i3)')i,Alist(i),dAlist(i),Nlist(i),MarkX(i)
  enddo
  close(71)
end program
''' + exact + '\n'


READER = '''
subroutine ReadCampaignSchedule
  use Tools
  use LocalData
  implicit none
  integer :: io,i,index
  real :: previous
  logical :: exists
  inquire(file='campaign_schedule.dat',exist=exists)
  if(.not.exists)error stop 'Missing campaign_schedule.dat'
  if(Nexact/=3)error stop 'Campaign requires exactly three output epochs'
  open(71,file='campaign_schedule.dat',status='old',action='read')
  read(71,*,iostat=io)Ntotal
  if(io/=0)error stop 'Unreadable campaign schedule length'
  if(Ntotal<1.or.Ntotal>NstepM)error stop 'Invalid campaign schedule length'
  Alist=0.;dAlist=0.;Nlist=0;MarkX=0
  previous=AEXPN0
  do i=1,Ntotal
    read(71,*,iostat=io)index,Alist(i),dAlist(i),Nlist(i),MarkX(i)
    if(io/=0)error stop 'Incomplete campaign schedule'
    if(index/=i.or.Alist(i)<=previous.or.dAlist(i)<=0.)error stop 'Invalid campaign step'
    if(abs((Alist(i)-previous)-dAlist(i))>4.*epsilon(previous)*Alist(i)) &
      error stop 'Campaign time increment disagrees with endpoint'
    if(Nlist(i)/=0.and.Nlist(i)/=1)error stop 'Invalid campaign output flag'
    if(Nlist(i)/=MarkX(i))error stop 'Campaign output must be exact'
    previous=Alist(i)
  enddo
  read(71,*,iostat=io)index
  if(io>=0)error stop 'Unexpected extra campaign schedule rows'
  close(71)
  if(count(Nlist(1:Ntotal)==1)/=3.or.Alist(Ntotal)/=1.)error stop 'Wrong campaign outputs'
  NlastX=Ntotal
  if(ISTEP>0)then
    if(ISTEP>Ntotal)error stop 'Checkpoint beyond campaign schedule'
    if(AEXPN/=Alist(ISTEP).or.ASTEP/=dAlist(ISTEP))error stop 'Incompatible campaign checkpoint'
  else
    if(abs(ASTEP-dAlist(1))>epsilon(ASTEP)*ASTEP)error stop 'IC velocity staggering disagrees with first step'
  endif
end subroutine ReadCampaignSchedule
'''


def controlled_main(source):
    start = source.index('    !---- make table for steps')
    end = source.index("    write (*, *) 'Step  a_expansion", start)
    return source[:start] + '    Call ReadCampaignSchedule\n\n' + source[end:] + READER


def generate(destination):
    destination.mkdir(parents=True, exist_ok=True)
    source = (REPO / 'PMP2main.f90').read_text()
    generator = destination / 'export_schedule.f90'
    generator.write_text(native_generator(source))
    exe = destination / 'export_schedule.exe'
    command = ['ifx', '-O3', '-fp-model', 'fast=1', str(generator), '-o', str(exe)]
    build = subprocess.run(command, cwd=destination, env=native_env(1), capture_output=True, text=True, check=True)
    run = subprocess.run([str(exe)], cwd=destination, env=native_env(1), capture_output=True, text=True, check=True)
    raw = destination / 'native-schedule.dat'
    normal = np.loadtxt(raw, skiprows=1)
    assert len(normal) == int(raw.read_text().splitlines()[0]) == 158
    assert normal[normal[:,3] == 1, 0].tolist() == [92,109,158]
    assert np.array_equal(normal[normal[:,3] == 1,1].astype('f4'), np.array([1/3,1/2,1], dtype='f4'))
    previous = np.float32(1.) / np.float32(101.)
    half = []
    for i, endpoint, width, output, exact in normal:
        endpoint = np.float32(endpoint)
        midpoint = np.float32((float(previous)+float(endpoint))/2)
        half.append([2*i-1, midpoint, np.float32(midpoint-previous), 0, 0])
        half.append([2*i, endpoint, np.float32(endpoint-midpoint), output, exact])
        previous = endpoint
    half = np.asarray(half)
    assert np.array_equal(half[1::2,1], normal[:,1])
    assert np.array_equal(half[1::2,3:], normal[:,3:])
    assert np.all(half[:,2] > 0)
    half_path = destination / 'half-schedule.dat'
    with half_path.open('w') as stream:
        stream.write(f'{len(half)}\n')
        np.savetxt(stream, half, fmt=['%6d','%.16e','%.16e','%d','%d'])
    main = destination / 'PMP2main.schedule.f90'
    main.write_text(controlled_main(source))
    report = dict(production_main_sha256=sha(REPO/'PMP2main.f90'), generated_main_sha256=sha(main),
                  generator_sha256=sha(generator), command=command, compiler_stdout=build.stdout,
                  generator_stdout=run.stdout, normal_steps=len(normal), half_steps=len(half),
                  normal_outputs=[92,109,158], half_outputs=[184,218,316],
                  normal_initial_step=float(normal[0,2]), half_initial_step=float(half[0,2]),
                  every_normal_endpoint_preserved=True, every_output_flag_preserved=True,
                  hashes={p.name:sha(p) for p in [raw,half_path]})
    write_json(ROOT/'timetable.json',report)
    return report


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,default=WORK/'timetable')
    generate(parser.parse_args().output.resolve())
