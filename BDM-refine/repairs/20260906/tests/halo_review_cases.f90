program halo_review_cases
  use LinkerList
  implicit none
  character(40) :: which
  integer :: n,j,i,stat
  integer*8 :: k
  real :: cntr(3)
  real*8,allocatable :: rr(:),pp(:)
  real*8 :: energy
  logical :: singular
  call get_command_argument(1,which)
  if(trim(which)=='heap_index')then
    ! Test the actual default-integer index arithmetic at a valid Np<1200^3.
    read(*,*)n,j
    i=2*j
    write(*,*)'HEAP_CHILD',n,j,i,i<=n
    stop
  endif
  if(trim(which)=='potential')then
    read(*,*)n
    allocate(rr(n),pp(n))
    read(*,*)rr
    call BdmHaloSphericalPotential(rr,2.d0,3.d0,pp,energy,singular)
    write(*,*)'POTENTIAL',singular,energy
    write(*,*)pp
    stop
  endif
  open(31,file='input.dat',status='old')
  read(31,*)n,Cell,MassOne,Om0,Ovdens,AEXPN,Rext,SlopeR
  read(31,*)cntr
  Box=32.;NGRID=128;Np=n;Nparticles=n;Nmaxima=1
  Nmx=-2;Nbx=65;Nmy=-2;Nby=65;Nmz=-2;Nbz=65
  allocate(Xpar(n),Ypar(n),Zpar(n),VX(n),VY(n),VZ(n),OriginalParticleId(n))
  do j=1,n
    read(31,*)Xpar(j),Ypar(j),Zpar(j),VX(j),VY(j),VZ(j),OriginalParticleId(j)
  enddo
  close(31)
  allocate(Mvir(1),Rvir(1),Mtotal(1),Xoff(1),xMaxx(1),yMaxx(1),zMaxx(1), &
           VxMaxx(1),VyMaxx(1),VzMaxx(1),EpotM(1),EkinM(1),LambdaM(1),RadRms(1), &
           VmaxM(1),RmaxM(1),Xax(1),Yax(1),Zax(1),Axba(1),Axca(1))
  allocate(Lst(Np),Label(Nmx:Nbx,Nmy:Nby,Nmz:Nbz))
  open(13,status='scratch')
  call List
  xMaxx=cntr(1);yMaxx=cntr(2);zMaxx=cntr(3)
  call BdmHaloMembershipInit
  call GetHalo(cntr(1),cntr(2),cntr(3),0.,0.,0.,1_8)
  write(*,*) 'VALUES',HaloStatus(1),Mvir(1),Mtotal(1),Rvir(1),EkinM(1),EpotM(1), &
       VmaxM(1),RmaxM(1),VxMaxx(1),VyMaxx(1),VzMaxx(1),RadRms(1),Xoff(1),LambdaM(1), &
       Axba(1),Axca(1),Xax(1),Yax(1),Zax(1)
  if(allocated(BoundParticleIds(1)%ids))then
    write(*,*) 'IDS',BoundParticleIds(1)%ids
  else
    write(*,*) 'IDS'
  endif
  ! Call again on the same candidate with an empty spatial neighbourhood:
  ! membership from the previous call must be cleared before an early return.
  call GetHalo(25.,25.,25.,0.,0.,0.,1_8)
  if(allocated(BoundParticleIds(1)%ids))error stop 'stale membership after empty call'
  if(Mvir(1)/=0..or.Mtotal(1)/=0.)error stop 'stale mass after empty call'
  write(*,*) 'EMPTY_RECALL_PASS'
end program halo_review_cases
