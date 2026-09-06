! Controlled regression fixtures; linked against extracted production routines.
program halo_cases
  use LinkerList
  implicit none
  character(80) :: which,arg
  integer :: n,q,a,b,c,nh
  integer*8 :: ip
  real :: step,zz,theta,r
  call get_command_argument(1,which)
  n=64; nh=1
  if(index(which,'shell')>0)n=1000
  if(trim(which)=='so'.or.trim(which)=='so_extended')n=100000
  if(trim(which)=='parallel')nh=4
  Np=n; Nparticles=n; Nmaxima=nh
  MassOne=1.e10; Om0=.3; Ovdens=200.; dLogR=.02; Rext=0.; SlopeR=.2
  Box=32.; NGRID=128; Cell=.5
  Nmx=0;Nmy=0;Nmz=0;Nbx=32;Nby=32;Nbz=32
  allocate(Xpar(n),Ypar(n),Zpar(n),VX(n),VY(n),VZ(n),OriginalParticleId(n))
  allocate(Mvir(nh),Rvir(nh),Mtotal(nh),Xoff(nh),xMaxx(nh),yMaxx(nh),zMaxx(nh))
  allocate(VxMaxx(nh),VyMaxx(nh),VzMaxx(nh),EpotM(nh),EkinM(nh),LambdaM(nh))
  allocate(RadRms(nh),VmaxM(nh),RmaxM(nh),Xax(nh),Yax(nh),Zax(nh),Axba(nh),Axca(nh))
  xMaxx=5.;yMaxx=5.;zMaxx=5.;VxMaxx=0.;VyMaxx=0.;VzMaxx=0.
  step=.0125
  if(trim(which)=='compact')step=.0025
  q=0
  do a=1,4
  do b=1,4
  do c=1,4
    q=q+1
    Xpar(q)=5.+(a-2.5)*step
    Ypar(q)=5.+(b-2.5)*step
    Zpar(q)=5.+(c-2.5)*step
  enddo
  enddo
  enddo
  if(trim(which)=='central')then
    Xpar(1)=5.;Ypar(1)=5.;Zpar(1)=5.
  endif
  if(trim(which)=='singular')then
    Xpar(:2)=5.;Ypar(:2)=5.;Zpar(:2)=5.
  endif
  if(trim(which)=='truncated')then
    HaloSearchRadius=.02; ParticleSearchRadius=.02
  endif
  if(trim(which)=='aperture_truncated')then
    Rext=.5; ParticleSearchRadius=.22
  endif
  if(trim(which)=='so'.or.trim(which)=='so_extended')then
    call get_command_argument(2,arg); read(arg,*)dLogR
    Cell=4.;Nmx=0;Nmy=0;Nmz=0;Nbx=8;Nby=8;Nbz=8
    MassOne=2.*(1.150e12*Om0)*Ovdens/n
    if(trim(which)=='so_extended')Rext=.25
  endif
  if(index(which,'shell')>0.or.index(which,'so')==1)then
    do q=1,n
      zz=1.-2.*(q-.5)/n
      theta=2.399963229728653*(q-1)
      r=.1
      if(index(which,'so')==1)r=2.*(q-.5)/n
      Xpar(q)=5.+r*sqrt(1.-zz*zz)*cos(theta)
      Ypar(q)=5.+r*sqrt(1.-zz*zz)*sin(theta)
      Zpar(q)=5.+r*zz
    enddo
  endif
  VX=0.;VY=0.;VZ=0.
  do q=1,n
    OriginalParticleId(q)=n-q+1
    if(trim(which)=='hot_shell')then
      VX(q)=(-1.)**q*2000.
      if(q<=100)VX(q)=(-1.)**q*500.
    endif
    if(trim(which)=='mixed_shell'.or.trim(which)=='bulk_shell')then
      VX(q)=(-1.)**q*20.
      if(q<=100)VX(q)=(-1.)**q*2000.
    endif
  enddo
  if(trim(which)=='bulk_shell')then
    VX=VX+100000.;VY=100000.;VZ=100000.
  endif
  allocate(Lst(Np),Label(Nmx:Nbx,Nmy:Nby,Nmz:Nbz))
  open(13,status='scratch')
  call List
  call ParametersDistinct
  open(50,file='phase.bin',access='stream',form='unformatted',status='replace')
  write(50) Xpar,Ypar,Zpar,VX,VY,VZ
  close(50)
  open(50,file='members.bin',access='stream',form='unformatted',status='replace')
  if(allocated(BoundParticleIds(1)%ids))write(50)BoundParticleIds(1)%ids
  close(50)
  do ip=1,nh
    write(*,'(a,2i8,21es25.16)') 'HALO ',ip,HaloStatus(ip),Mvir(ip),Mtotal(ip),Rvir(ip), &
      EkinM(ip),EpotM(ip),VmaxM(ip),RmaxM(ip),VxMaxx(ip),VyMaxx(ip),VzMaxx(ip), &
      RadRms(ip),Xoff(ip),LambdaM(ip),Axba(ip),Axca(ip),Xax(ip),Yax(ip),Zax(ip), &
      MassOne,Rext,dLogR
    if(nh>1)then
      if(any(BoundParticleIds(ip)%ids/=BoundParticleIds(1)%ids))error stop 'parallel membership differs'
    endif
  enddo
end program halo_cases
