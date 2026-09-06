program particle_followup
  use LinkerList
  use omp_lib, only: omp_get_wtime
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  implicit none
  character(80) :: which,arg
  integer :: n,repetitions,repeat,q,axis,sx,sy,sz
  integer*8 :: row,id,at,ix,iy,iz,peak,beforeWords,originalCount
  integer*8, allocatable :: expectedLst(:),expectedLabel(:,:,:)
  real, allocatable :: original(:,:)
  real*8 :: started,elapsed,point(3),stored(3),bounds(3),box64,width64
  logical :: bufferCase,benchmark
  call get_command_argument(1,which)
  n=32769; repetitions=1
  call get_command_argument(2,arg);if(len_trim(arg)>0)read(arg,*)n
  call get_command_argument(3,arg);if(len_trim(arg)>0)read(arg,*)repetitions
  Box=32.;NGRID=128;dBuffer=0.;Cell=1.
  if(trim(which)=='buffer_fractional')Box=10.3
  if(trim(which)=='buffer_halfbox')then
    Box=16.;NGRID=16
  endif
  Np=int(n,8);Nparticles=Np;originalCount=Np
  allocate(Xpar(n),Ypar(n),Zpar(n),VX(n),VY(n),VZ(n),original(6,n))
  call PrepareParticleSearch
  do q=1,n
    Xpar(q)=real(mod(37_8*(q-1),65536_8))*Box/65536.
    Ypar(q)=real(mod(101_8*(q-1),65536_8))*Box/65536.
    Zpar(q)=real(mod(271_8*(q-1),65536_8))*Box/65536.
    VX(q)=real(mod(q,29))*.125;VY(q)=-real(mod(q,17))*.25;VZ(q)=real(mod(q,31))*.5
  enddo
  if(n>=10)then
    Xpar(:10)=[0.,nearest(0.,1.),dBuffer,nearest(dBuffer,-1.),nearest(dBuffer,1.), &
      Box-dBuffer,nearest(Box-dBuffer,-1.),nearest(Box-dBuffer,1.),Box/2.,nearest(Box,-1.)]
    Ypar(:10)=Xpar(:10);Zpar(:10)=Xpar(:10)
    VX(1)=-0.;VY(2)=transfer(1,1.)
  endif
  original(1,:)=Xpar;original(2,:)=Ypar;original(3,:)=Zpar
  original(4,:)=VX;original(5,:)=VY;original(6,:)=VZ
  memoryWords=6_8*Np;peakWords=memoryWords
  open(13,status='scratch');open(17,status='scratch')
  select case(trim(which))
  case('invalid_negative')
    Xpar(1)=-.1;call AddBuffer;error stop 'missing negative-coordinate rejection'
  case('invalid_upper')
    Xpar(1)=Box;call AddBuffer;error stop 'missing upper-coordinate rejection'
  case('invalid_nan')
    Xpar(1)=ieee_value(0.,ieee_quiet_nan);call AddBuffer;error stop 'missing NaN rejection'
  case('invalid_count')
    Nparticles=Np+1;call AddBuffer;error stop 'missing count mismatch rejection'
  case('invalid_active')
    allocate(OriginalParticleId(0));call AddBuffer;error stop 'missing active-buffer rejection'
  case('invalid_capacity')
    Np=huge(0_8);Nparticles=Np;call AddBuffer;error stop 'missing safe-int64 rejection'
  case('invalid_negative_count')
    Np=-1_8;Nparticles=Np;call AddBuffer;error stop 'missing negative-count rejection'
  end select
  bufferCase=index(trim(which),'buffer')==1.or.index(trim(which),'periodic')>0.or.trim(which)=='list_benchmark'
  benchmark=index(trim(which),'benchmark')>0
  if(bufferCase)then
    elapsed=0.d0;peak=0_8
    do repeat=1,merge(repetitions,1,trim(which)=='buffer_benchmark')
      if(repeat>1)then
        deallocate(Xpar,Ypar,Zpar,VX,VY,VZ,OriginalParticleId)
        Np=originalCount;Nparticles=Np
        allocate(Xpar(n),Ypar(n),Zpar(n),VX(n),VY(n),VZ(n))
        Xpar=original(1,:);Ypar=original(2,:);Zpar=original(3,:)
        VX=original(4,:);VY=original(5,:);VZ=original(6,:)
        memoryWords=6_8*Np;peakWords=memoryWords
      endif
      started=omp_get_wtime();call AddBuffer;elapsed=elapsed+omp_get_wtime()-started
      peak=max(peak,peakWords)
    enddo
    call check_buffer
    if(trim(which)=='buffer_benchmark') &
      print *, 'FOLLOWUP_METRICS buffer',elapsed/repetitions,repetitions,originalCount,Np,peak
    call SizeList
  else
    Cell=1.;Nmx=-1;Nbx=1;Nmy=-1;Nby=1;Nmz=-3;Nbz=3
    select case(trim(which))
    case('list_i16_range');Nmz=-32768;Nbz=32767
    case('list_i16_low');Nmz=-32768;Nbz=-32760
    case('list_i16_high');Nmz=32760;Nbz=32767
    case('list_i32_low');Nmz=-32769;Nbz=-32760
    case('list_i32_high');Nmz=32760;Nbz=32768
    case('list_i32_min');Nmz=-huge(0)-1;Nbz=Nmz+2
    case('list_i32_max');Nbz=huge(0);Nmz=Nbz-2
    case('list_empty_slabs');Nmz=-1;Nbz=1
    case('invalid_cell');Cell=0.
    case('invalid_bounds');Nmz=1;Nbz=0
    case('invalid_list_nan');Zpar(1)=ieee_value(0.,ieee_quiet_nan)
    case('invalid_cache_memory');MaxMemory=0.
    end select
    if(trim(which)/='invalid_list_nan')then
      do q=1,n
        Xpar(q)=real(mod(q,9))*.5-2.;Ypar(q)=real(mod(q,7))*.5-1.5
        select case(mod(q,9))
        case(0);Zpar(q)=real(Nmz)
        case(1);Zpar(q)=real(Nmz)+1.
        case(2);Zpar(q)=nearest(real(Nmz)+1.,-1.)
        case(3);Zpar(q)=nearest(real(Nmz)+1.,1.)
        case(4);Zpar(q)=real(Nbz)
        case(5);Zpar(q)=real(Nbz)+1.
        case(6);Zpar(q)=nearest(real(Nbz)+1.,-1.)
        case(7);Zpar(q)=nearest(real(Nbz)+1.,1.)
        case(8);Zpar(q)=.5*(real(Nmz)+real(Nbz))
        end select
      enddo
    endif
  endif
  allocate(Lst(Np),Label(Nmx:Nbx,Nmy:Nby,Nmz:Nbz))
  beforeWords=memoryWords;peakWords=memoryWords
  started=omp_get_wtime()
  do repeat=1,merge(repetitions,1,trim(which)=='list_benchmark')
    call List
  enddo
  elapsed=omp_get_wtime()-started
  if(memoryWords/=beforeWords)error stop 'list cache memory-accounting leak'
  if(trim(which)=='list_benchmark') &
    print *, 'FOLLOWUP_METRICS list',elapsed/repetitions,repetitions,originalCount,Np,peakWords-beforeWords
  allocate(expectedLst(Np),expectedLabel(Nmx:Nbx,Nmy:Nby,Nmz:Nbz));expectedLabel=0_8
  do row=1,Np
    point=[dble(Xpar(row)),dble(Ypar(row)),dble(Zpar(row))]
    if(allocated(OriginalParticleId))then
      id=OriginalParticleId(row)
      do axis=1,3
        stored(axis)=point(axis)
        point(axis)=dble(original(axis,id))
        point(axis)=point(axis)+dble(nint((stored(axis)-point(axis))/dble(Box)))*dble(Box)
      enddo
    endif
    ix=reference_cell(point(1)/dble(Cell),Nmx,Nbx)
    iy=reference_cell(point(2)/dble(Cell),Nmy,Nby)
    iz=reference_cell(point(3)/dble(Cell),Nmz,Nbz)
    expectedLst(row)=expectedLabel(ix,iy,iz);expectedLabel(ix,iy,iz)=row
  enddo
  if(any(expectedLst/=Lst).or.any(expectedLabel/=Label))error stop 'independent list order/cell mismatch'
  open(77,file='result.bin',access='stream',form='unformatted',status='replace')
  write(77)Np,Nparticles,Nmx,Nbx,Nmy,Nby,Nmz,Nbz,Xpar,Ypar,Zpar,VX,VY,VZ,Lst,Label
  if(allocated(OriginalParticleId))write(77)OriginalParticleId
  close(77)
  print *, 'FOLLOWUP_PASS',trim(which),originalCount,Np
contains
  integer*8 function reference_cell(value,low,high)result(cellId)
    real*8,intent(in)::value
    integer,intent(in)::low,high
    cellId=min(max(int(low,8),ceiling(value,kind=8)-1_8),int(high,8))
  end function
  subroutine check_buffer
    integer::c
    if(size(Xpar,kind=8)/=Np.or.size(OriginalParticleId,kind=8)/=Np)error stop 'buffer capacity'
    do row=1,originalCount
      if(OriginalParticleId(row)/=row)error stop 'original identity order'
      stored=[dble(Xpar(row)),dble(Ypar(row)),dble(Zpar(row))]
      if(any(stored/=dble(original(:3,row))))error stop 'primary position copy'
    enddo
    at=originalCount;box64=dble(Box);width64=dble(dBuffer)
    do row=1,originalCount
      do sz=-1,1
      do sy=-1,1
      do sx=-1,1
        if(sx==0.and.sy==0.and.sz==0)cycle
        point=dble(original(:3,row))+dble([sx,sy,sz])*box64
        if(any(point< -width64).or.any(point>box64+width64))cycle
        at=at+1_8
        if(at>Np)error stop 'oracle images exceed allocation'
        if(OriginalParticleId(at)/=row)error stop 'ghost identity or image order'
        if(any(transfer([Xpar(at),Ypar(at),Zpar(at)],[0],3)/=transfer(real(point),[0],3))) &
          error stop 'ghost coordinate bits'
      enddo
      enddo
      enddo
    enddo
    if(at/=Np)error stop 'missing or repeated oracle image'
    do row=1,Np
      id=OriginalParticleId(row)
      if(any(transfer([VX(row),VY(row),VZ(row)],[0],3)/=transfer(original(4:6,id),[0],3)))error stop 'velocity bits'
    enddo
    if(memoryWords/=8_8*Np)error stop 'buffer scratch memory-accounting leak'
  end subroutine
end program particle_followup
