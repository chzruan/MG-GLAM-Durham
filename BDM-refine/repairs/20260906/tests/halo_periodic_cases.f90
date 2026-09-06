program periodic_precision_cases
  use LinkerList
  implicit none
  character(40) :: which
  integer :: n,nq,q,component
  integer*8 :: row,ip
  integer*8,allocatable :: rows(:)
  real,allocatable :: query(:,:)
  real*8,allocatable :: radius(:),radii(:)
  real*8 :: position(3)
  call get_command_argument(1,which)
  open(10,file='input.bin',form='unformatted',access='stream',status='old')
  read(10)n,nq
  Np=n;Nparticles=n;Nmaxima=nq
  allocate(Xpar(n),Ypar(n),Zpar(n),VX(n),VY(n),VZ(n),query(3,nq),radius(nq))
  read(10)Xpar,Ypar,Zpar,VX,VY,VZ,query,radius
  close(10)
  Box=32.;NGRID=128;dBuffer=5.;MassOne=1.e10;Om0=.3;Ovdens=200.;Rext=0.;SlopeR=.2;dLogR=.02
  if(trim(which)=='half_box')NGRID=16
  open(13,status='scratch')
  call AddBuffer
  call SizeList
  allocate(Lst(Np),Label(Nmx:Nbx,Nmy:Nby,Nmz:Nbz))
  call List
  open(10,file='images.bin',form='unformatted',access='stream',status='replace')
  write(10)Np,Xpar,Ypar,Zpar,OriginalParticleId
  do row=1,Np
    call BdmParticlePosition(row,position)
    write(10)position
  enddo
  close(10)
  open(10,file='gathers.bin',form='unformatted',access='stream',status='replace')
  do q=1,nq
    call BdmHaloGather(query(1,q),query(2,q),query(3,q),radius(q),rows,radii)
    write(10)size(rows,kind=8),OriginalParticleId(rows),radii
  enddo
  close(10)
  if(trim(which)=='centering'.or.trim(which)=='halo')then
    allocate(Mvir(nq),Rvir(nq),Mtotal(nq),Xoff(nq),xMaxx(nq),yMaxx(nq),zMaxx(nq))
    allocate(VxMaxx(nq),VyMaxx(nq),VzMaxx(nq),EpotM(nq),EkinM(nq),LambdaM(nq))
    allocate(RadRms(nq),VmaxM(nq),RmaxM(nq),Xax(nq),Yax(nq),Zax(nq),Axba(nq),Axca(nq))
    xMaxx=query(1,:);yMaxx=query(2,:);zMaxx=query(3,:);Xoff=100.
    call FindDistinctCandidates
    open(10,file='centres.bin',form='unformatted',access='stream',status='replace')
    write(10)xMaxx,yMaxx,zMaxx,VxMaxx,VyMaxx,VzMaxx,Mvir
    close(10)
    if(trim(which)=='halo')then
      call ParametersDistinct
      open(10,file='haloes.bin',form='unformatted',access='stream',status='replace')
      do ip=1,Nmaxima
        write(10)Mvir(ip),Mtotal(ip),Rvir(ip),EkinM(ip),EpotM(ip),RadRms(ip)
        write(10)HaloStatus(ip)
        if(allocated(BoundParticleIds(ip)%ids))then
          write(10)size(BoundParticleIds(ip)%ids,kind=8),BoundParticleIds(ip)%ids
        else
          write(10)0_8
        endif
      enddo
      close(10)
    endif
  endif
end program periodic_precision_cases
