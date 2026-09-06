! Independent assertions for particle preparation, centring and restoration.
program particles_cases
  use LinkerList
  use omp_lib, only: omp_get_wtime
  implicit none
  character(80) :: which,arg
  integer :: q,h,n,ix,iy,iz,repetition,i1,i2,j1,j2,k1,k2,countFound,expectedCount
  integer*8 :: jp,id,originalCount,seenCount,lastId
  integer*8, allocatable :: seen(:),referenceLst(:),referenceLabel(:,:,:)
  real :: origin,truth(3),center(3),radius,xscale,vscale
  real*8 :: displacement,velocityError,dx,dy,dz,started,elapsed
  real, allocatable :: saved(:,:),original(:,:)
  logical :: emptyCase
  call get_command_argument(1,which)
  open(13,status='scratch');open(17,status='scratch')
  select case(trim(which))
  case('centering_near','centering_far','centering_one','centering_empty')
    n=262144; Nmaxima=1; Box=1024.; Cell=.5; MassOne=1.e10
    origin=5.125
    if(trim(which)=='centering_far')origin=1000.125
    if(trim(which)=='centering_one')n=1
    emptyCase=trim(which)=='centering_empty'
    call particles(n);call maxima()
    do q=1,n
      Xpar(q)=origin+(mod(q-1,8)-3.5)/256.
      Ypar(q)=origin+(mod((q-1)/8,8)-3.5)/256.
      Zpar(q)=origin+(mod((q-1)/64,8)-3.5)/256.
      VX(q)=100000.125+(mod(q-1,8)-3.5)/64.
      VY(q)=-300000.25; VZ(q)=7.25
    enddo
    if(emptyCase)Xpar=origin+2.
    truth=[origin,origin,origin]
    if(n==1)truth=[Xpar(1),Ypar(1),Zpar(1)]
    Nmx=int(origin/Cell)-4;Nbx=Nmx+10;Nmy=Nmx;Nby=Nbx;Nmz=Nmx;Nbz=Nbx
    allocate(Lst(Np),Label(Nmx:Nbx,Nmy:Nby,Nmz:Nbz));call List
    xMaxx=origin;yMaxx=origin;zMaxx=origin;Xoff=100.
    call FindDistinctCandidates
    displacement=sqrt(sum((dble([xMaxx(1),yMaxx(1),zMaxx(1)])-dble(truth))**2))
    if(emptyCase)then
      if(any([Mvir(1),Rvir(1),VxMaxx(1),VyMaxx(1),VzMaxx(1)]/=0.))error stop 'empty aperture values'
      if(any([xMaxx(1),yMaxx(1),zMaxx(1)]/=origin))error stop 'empty aperture centre'
    else
      velocityError=max(abs(dble(VxMaxx(1))-sum(dble(VX))/n), &
                        abs(dble(VyMaxx(1))-sum(dble(VY))/n), &
                        abs(dble(VzMaxx(1))-sum(dble(VZ))/n))
      if(displacement>0.d0)error stop 'symmetric cloud centre changed under translation'
      if(velocityError>1.d-8)error stop 'bulk velocity accumulation is inaccurate'
      if(Mvir(1)/MassOne/=real(n))error stop 'centering particle count'
      if(xMaxx(1)<minval(Xpar).or.xMaxx(1)>maxval(Xpar))error stop 'centre outside cloud'
      print *, 'CENTRE_DISPLACEMENT_VELOCITY_ERROR',displacement,velocityError
    endif
  case('buffer_face','buffer_corner','buffer_search','buffer_query','buffer_halfbox','buffer_smallbox')
    n=128; Box=32.;NGRID=128;dBuffer=5.
    if(trim(which)=='buffer_query'.or.trim(which)=='buffer_halfbox'.or.trim(which)=='buffer_smallbox')n=4096
    if(trim(which)=='buffer_halfbox')then
      Box=16.;NGRID=16
    endif
    if(trim(which)=='buffer_smallbox')then
      Box=4.;NGRID=16
    endif
    call particles(n)
    Xpar=0.;Ypar=Box/2.;Zpar=Box/2.;VX=7.25;VY=-.25;VZ=101.
    if(trim(which)=='buffer_corner')then
      Xpar=.01;Ypar=.01;Zpar=.01
    endif
    if(trim(which)=='buffer_search')Xpar=26.2
    if(trim(which)=='buffer_query'.or.trim(which)=='buffer_halfbox'.or.trim(which)=='buffer_smallbox')then
      do q=1,n
        Xpar(q)=real(mod(37_8*q,4096_8))*Box/4096.
        Ypar(q)=real(mod(101_8*q,4096_8))*Box/4096.
        Zpar(q)=real(mod(271_8*q,4096_8))*Box/4096.
      enddo
      Xpar(1)=0.;Ypar(1)=0.;Zpar(1)=0.
    endif
    originalCount=Np
    allocate(original(3,n));original(1,:)=Xpar;original(2,:)=Ypar;original(3,:)=Zpar
    call AddBuffer
    if(size(Xpar,kind=8)/=Np.or.size(OriginalParticleId,kind=8)/=Np)error stop 'exact buffer capacity'
    if(any(OriginalParticleId(:n)/=[(int(q,8),q=1,n)]))error stop 'original row order'
    if(any(VX/=7.25).or.any(VY/=-.25).or.any(VZ/=101.))error stop 'ghost velocity copying'
    do jp=1,Np
      id=OriginalParticleId(jp)
      if(id<1.or.id>n)error stop 'invalid ghost row identity'
      dx=(dble(Xpar(jp))-original(1,id))/dble(Box)
      dy=(dble(Ypar(jp))-original(2,id))/dble(Box)
      dz=(dble(Zpar(jp))-original(3,id))/dble(Box)
      if(max(abs(dx-nint(dx)),abs(dy-nint(dy)),abs(dz-nint(dz)))>1.d-7)error stop 'wrong periodic image coordinates'
    enddo
    if(trim(which)=='buffer_face')then
      if(Np/=2*n.or.count(Xpar==Box)/=n)error stop 'missing exact upper face images'
    endif
    if(trim(which)=='buffer_corner')then
      if(Np/=8*n)error stop 'corner buffer insufficient capacity'
    endif
    call SizeList
    if(ParticleSearchRadius>dBuffer)error stop 'insufficient configured buffer reach'
    if(ParticleSearchRadius>=Box/2.)error stop 'periodic aperture can count two images'
    allocate(Lst(Np),Label(Nmx:Nbx,Nmy:Nby,Nmz:Nbz));call List
    allocate(seen(n))
    do repetition=1,8
      center=[.1,.125,.25]
      if(mod(repetition,2)==0)center(1)=Box-.1
      if(mod(repetition/2,2)==0)center(2)=Box-.125
      if(mod(repetition/4,2)==0)center(3)=Box-.25
      radius=ParticleSearchRadius
      if(trim(which)=='buffer_search')then
        center=[.1,16.,16.];radius=6.
      endif
      seen=0
      call Limits(center(1),center(2),center(3),radius,i1,i2,j1,j2,k1,k2)
      do iz=k1,k2
      do iy=j1,j2
      do ix=i1,i2
        jp=Label(ix,iy,iz)
        do while(jp/=0_8)
          if(sum((dble([Xpar(jp),Ypar(jp),Zpar(jp)])-dble(center))**2)<dble(radius)**2) &
            seen(OriginalParticleId(jp))=seen(OriginalParticleId(jp))+1
          jp=Lst(jp)
        enddo
      enddo
      enddo
      enddo
      expectedCount=0
      do q=1,n
        dx=abs(dble(original(1,q))-center(1));dx=min(dx,dble(Box)-dx)
        dy=abs(dble(original(2,q))-center(2));dy=min(dy,dble(Box)-dy)
        dz=abs(dble(original(3,q))-center(3));dz=min(dz,dble(Box)-dz)
        if(dx*dx+dy*dy+dz*dz<dble(radius)**2)then
          expectedCount=expectedCount+1
          if(seen(q)/=1)error stop 'missing or repeated periodic member'
        else
          if(seen(q)/=0)error stop 'spurious periodic member'
        endif
      enddo
      if(sum(seen)/=expectedCount)error stop 'periodic membership differs from minimum image oracle'
      if(trim(which)=='buffer_search'.and.expectedCount/=n)error stop 'radius-six search oracle'
    enddo
    print *, 'EXACT_BUFFER_QUERY',n,Np,HaloSearchRadius,ParticleSearchRadius
  case('restore_roundtrip','restore_buffered')
    n=1000;Box=100.;NGRID=768;AEXPN=.7;dBuffer=5.
    call particles(n)
    do q=1,n
      Xpar(q)=1.+mod(q*37,767)+.12345
      Ypar(q)=Xpar(q);Zpar(q)=Xpar(q)
      VX(q)=q*.1234567;VY(q)=VX(q);VZ(q)=VX(q)
    enddo
    ! Include exactly representable domain and velocity edge values.
    Xpar(1)=1.;Xpar(2)=real(NGRID)+.5;VX(1)=-0.;VY(2)=transfer(1,1.)
    allocate(saved(n,6))
    saved(:,1)=Xpar;saved(:,2)=Ypar;saved(:,3)=Zpar
    saved(:,4)=VX;saved(:,5)=VY;saved(:,6)=VZ
    do repetition=1,2
      call RescaleCoords(1)
      if(trim(which)=='restore_buffered')call AddBuffer
      ! Analysis mutations must affect only its private workspace.
      Xpar(1)=17.;VX(1)=777.
      call RemoveBuffer(int(n,8))
      if(any(transfer(Xpar,[0],n)/=transfer(saved(:,1),[0],n)))error stop 'X simulation bits changed'
      if(any(transfer(Ypar,[0],n)/=transfer(saved(:,2),[0],n)))error stop 'Y simulation bits changed'
      if(any(transfer(Zpar,[0],n)/=transfer(saved(:,3),[0],n)))error stop 'Z simulation bits changed'
      if(any(transfer(VX,[0],n)/=transfer(saved(:,4),[0],n)))error stop 'VX simulation bits changed'
      if(any(transfer(VY,[0],n)/=transfer(saved(:,5),[0],n)))error stop 'VY simulation bits changed'
      if(any(transfer(VZ,[0],n)/=transfer(saved(:,6),[0],n)))error stop 'VZ simulation bits changed'
      if(allocated(OriginalParticleId).or.allocated(BdmPMX).or.allocated(BdmPMVX))error stop 'workspace lifetime leak'
      if(Nparticles/=n.or.Np/=n.or.BdmPMCount/=0_8)error stop 'restoration count mismatch'
      if(memoryWords/=0_8)error stop 'particle workspace memory accounting leak'
    enddo
    print *, 'RESTORED_EXACTLY',n,repetition-1
  case('list_check','list_benchmark')
    n=32
    call get_command_argument(2,arg)
    if(len_trim(arg)>0)read(arg,*)n
    Box=64.;Cell=1.;NGRID=128
    call particles(n**3)
    do jp=1,Np
      Xpar(jp)=mod(jp-1,int(n,8))*Box/n+.001
      Ypar(jp)=mod((jp-1)/n,int(n,8))*Box/n+.001
      Zpar(jp)=((jp-1)/(int(n,8)*n))*Box/n+.001
    enddo
    Nmx=-1;Nbx=65;Nmy=-2;Nby=66;Nmz=-3;Nbz=67
    if(trim(which)=='list_check')then
      Nmx=-20;Nbx=40;Xpar=Xpar*.5
    endif
    allocate(Lst(Np),Label(Nmx:Nbx,Nmy:Nby,Nmz:Nbz))
    started=omp_get_wtime()
    do repetition=1,10
      call List
    enddo
    elapsed=(omp_get_wtime()-started)/10
    allocate(referenceLst(Np),referenceLabel(Nmx:Nbx,Nmy:Nby,Nmz:Nbz));referenceLabel=0_8
    do jp=1,Np
      ix=min(max(Nmx,ceiling(Xpar(jp)/Cell)-1),Nbx)
      iy=min(max(Nmy,ceiling(Ypar(jp)/Cell)-1),Nby)
      iz=min(max(Nmz,ceiling(Zpar(jp)/Cell)-1),Nbz)
      referenceLst(jp)=referenceLabel(ix,iy,iz);referenceLabel(ix,iy,iz)=jp
    enddo
    if(any(referenceLst/=Lst).or.any(referenceLabel/=Label))error stop 'list differs from deterministic row-order oracle'
    allocate(seen(Np));seen=0_8
    do iz=Nmz,Nbz
    do iy=Nmy,Nby
    do ix=Nmx,Nbx
      jp=Label(ix,iy,iz);lastId=Np+1_8
      do while(jp/=0_8)
        if(jp<1_8.or.jp>Np.or.jp>=lastId)error stop 'invalid particle list order or cycle'
        seen(jp)=seen(jp)+1_8;lastId=jp;jp=Lst(jp)
      enddo
    enddo
    enddo
    enddo
    if(any(seen/=1_8))error stop 'list missed or repeated particles'
    print *, 'LIST_SECONDS_COUNT',elapsed,sum(seen),Np
  case default
    error stop 'unknown particle repair test'
  end select
  print *, 'PARTICLES_TEST_PASS',trim(which)
contains
  subroutine particles(n)
    integer,intent(in)::n
    Np=int(n,8);Nparticles=Np
    allocate(Xpar(n),Ypar(n),Zpar(n),VX(n),VY(n),VZ(n))
  end subroutine
  subroutine maxima()
    allocate(Mvir(1),Rvir(1),VxMaxx(1),VyMaxx(1),VzMaxx(1),xMaxx(1),yMaxx(1),zMaxx(1),Xoff(1))
  end subroutine
end program particles_cases
