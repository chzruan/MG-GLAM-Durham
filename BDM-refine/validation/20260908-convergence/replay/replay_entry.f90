! Both variants use this identical entry and common DENSIT implementation.
program ConvergenceReplay
  use LinkerList, only: BDM,outputName
  use Density, only: DENSIT
  use Tools
  use ConvergenceGrid
  use iso_fortran_env, only: int32,int64
  use omp_lib, only: omp_set_num_threads,omp_set_dynamic
  implicit none
  integer :: step,target,threads,j,status,io,unit,l,nargs,original_grid
  integer(int64) :: file_bytes,grid,particles
  integer(int32) :: saved_step
  real :: memory_used,saved_a,saved_box,local_aexp0,local_partw,local_au0
  character(4096) :: arguments(5)
  character(64) :: header_path
  logical :: exists
  nargs=command_argument_count()
  step=-1;target=0;threads=0
  if(nargs/=3.and.nargs/=5)error stop 'Usage: replay STEP TARGET_GRID THREADS [write|read DENSITY_PATH]'
  arguments=''
  do j=1,nargs
    call get_command_argument(j,arguments(j),status=status)
    if(status/=0.or.len_trim(arguments(j))==0)error stop 'Invalid or oversized replay argument'
  enddo
  read(arguments(1),*,iostat=io)step
  if(io/=0.or.step<0.or.step>9999)error stop 'Invalid snapshot step'
  read(arguments(2),*,iostat=io)target
  if(io/=0.or..not.GridPowerOfTwo(target).or.target>16384)error stop 'Invalid analysis grid'
  read(arguments(3),*,iostat=io)threads
  if(io/=0.or.threads<1)error stop 'Invalid thread count'
  if(nargs==5)then
    if(arguments(4)/='write'.and.arguments(4)/='read')error stop 'Invalid density mode'
  endif
  call omp_set_dynamic(.false.)
  call omp_set_num_threads(threads)
  write(header_path,'(a,i4.4,a)')'PMcrd.',step,'.DAT'
  open(newunit=unit,file=trim(header_path),form='unformatted',status='old',action='read',iostat=io)
  if(io/=0)error stop 'Missing replay snapshot header'
  read(unit,iostat=io)HEADER,AEXPN,local_aexp0,AMPLT,ASTEP,ISTEP,local_partw, &
    TINTG,EKIN,EKIN1,EKIN2,local_au0,AEU0,NROW,NGRID,Nrealization,Nseed,Om,OmL,hubble,Nparticles,extras
  if(io/=0)error stop 'Invalid replay snapshot header'
  close(unit)
  if(ISTEP/=step)error stop 'Replay step/header mismatch'
  if(NROW<1.or.NROW>=1200.or.Nparticles/=int(NROW,int64)**3)error stop 'Replay particle count is unsupported'
  if(.not.GridPowerOfTwo(NGRID).or.NGRID>16384)error stop 'Replay source grid is unsupported'
  original_grid=NGRID
  call ReadDataPM(step,'')
  if(NGRID/=original_grid.or.ISTEP/=step)error stop 'Replay header changed during read'
  ! Main simulations use GR. Standalone old/current entries do not read Setup.dat.
  MG_flag=0;MG_test=0;MG_model=3
  call ChangeAnalysisGrid(target)
  memory_used=Memory(1_int64*NGRID*NGRID*NGRID)
  allocate(FI(NGRID,NGRID,NGRID))
  if(nargs==5.and.arguments(4)=='read')then
    open(newunit=unit,file=trim(arguments(5)),access='stream',form='unformatted',status='old', &
      action='read',iostat=io)
    if(io/=0)error stop 'Cannot read fixed density'
    inquire(unit=unit,size=file_bytes)
    if(file_bytes/=28_int64+4_int64*NGRID*NGRID*NGRID)error stop 'Fixed density byte count mismatch'
    read(unit,iostat=io)grid,particles,saved_step,saved_a,saved_box
    if(io/=0)error stop 'Invalid fixed density header'
    if(grid/=NGRID.or.particles/=Nparticles.or.saved_step/=ISTEP.or.saved_a/=AEXPN.or.saved_box/=Box) &
      error stop 'Fixed density metadata mismatch'
    do l=1,NGRID
      read(unit,iostat=io)FI(:,:,l)
      if(io/=0)error stop 'Fixed density plane read failed'
    enddo
    close(unit)
  else
    call DENSIT
    if(nargs==5)then
      inquire(file=trim(arguments(5)),exist=exists)
      if(exists)error stop 'Refuse to replace fixed density'
      open(newunit=unit,file=trim(arguments(5)),access='stream',form='unformatted',status='new', &
        action='write',iostat=io)
      if(io/=0)error stop 'Cannot create fixed density'
      write(unit,iostat=io)int(NGRID,int64),Nparticles,int(ISTEP,int32),AEXPN,Box
      if(io/=0)error stop 'Fixed density header write failed'
      do l=1,NGRID
        write(unit,iostat=io)FI(:,:,l)
        if(io/=0)error stop 'Fixed density plane write failed'
      enddo
      close(unit,iostat=io)
      if(io/=0)error stop 'Fixed density close failed'
    endif
  endif
  if(nargs==5)write(*,'(a,a)')'REPLAY DENSITY mode=',trim(arguments(4))
  ! Do not preload BDM.config: preserve the exact legacy first-call behavior.
  call BDM(0)
  write(*,'(a,a)')'REPLAY CATALOGUE ',trim(outputName)
  print '(a)','REPLAY COMPLETE'
end program
