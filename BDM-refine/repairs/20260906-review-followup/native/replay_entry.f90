program BdmReplayEntry
  use BdmReplayControl
  implicit none
  character(4096) :: arguments(5)
  integer :: step,threads,passes,j,status,io
  real :: memory_used
  step=-1;threads=0;passes=0
  if(command_argument_count()/=5)error stop 'Usage: replay step threads passes density-file read|write'
  do j=1,5
    call get_command_argument(j,arguments(j),status=status)
    if(status/=0.or.len_trim(arguments(j))==0)error stop 'Invalid or oversized replay argument'
  enddo
  read(arguments(1),*,iostat=io)step
  if(io/=0)error stop 'Invalid replay snapshot step'
  read(arguments(2),*,iostat=io)threads
  if(io/=0.or.threads<1)error stop 'Invalid replay thread count'
  read(arguments(3),*,iostat=io)passes
  if(io/=0.or.passes<1.or.passes>99)error stop 'Invalid replay pass count'
  if(trim(arguments(5))/='read'.and.trim(arguments(5))/='write')error stop 'Invalid replay density mode'
  call omp_set_dynamic(.false.)
  call omp_set_num_threads(threads)
  call ReplaySnapshotGate(step)
  call ReadDataPM(step,'')
  memory_used=Memory(1_8*NGRID*NGRID*NGRID)
  allocate(FI(NGRID,NGRID,NGRID))
  call RunBdmReplay(threads,passes,trim(arguments(4)),trim(arguments(5)))
end program BdmReplayEntry
