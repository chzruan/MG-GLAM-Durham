! The Python driver supplies an entire equal-mass NROW^3 box, including the
! remote particles. Only the candidate-centred subset is a compact cold halo.
program normalization_cases
  use LinkerList
  implicit none
  integer :: j,n,input_nrow
  real :: centre(3),mass_multiplier
  character(32) :: which
  call get_command_argument(1,which)
  open(13,status='scratch')
  read(*,*)Box,NROW,NGRID,Om,AEXPN,iVirial,Rext,mass_multiplier
  OmL=1.-Om; dLogR=.02;dLogP=.02;NradP=100;MassMin=0.;SlopeR=.2;dBuffer=5.
  input_nrow=NROW
  if(trim(which)=='invalid_nrow')NROW=8
  call BeginCataloguePublication('catalogue.dat')
  call SetParameters
  if(trim(which)=='invalid_nrow')NROW=input_nrow
  ! An explicit storage perturbation case checks that the actual stored mass
  ! remains the normalization, independent of its cosmological constructor.
  MassOne=MassOne*mass_multiplier
  Cell=2.*Box/NGRID; Nparticles=int(input_nrow,8)**3
  if(trim(which)=='invalid_nrow')Nparticles=512
  Np=Nparticles;n=int(Np);Nmaxima=1;centre=Box/4.
  Nmx=-2;Nmy=-2;Nmz=-2;Nbx=NGRID/2+2;Nby=Nbx;Nbz=Nbx
  allocate(Xpar(n),Ypar(n),Zpar(n),VX(n),VY(n),VZ(n),OriginalParticleId(n))
  open(31,file='phase.bin',access='stream',form='unformatted',status='old')
  do j=1,n
    read(31)Xpar(j),Ypar(j),Zpar(j),VX(j),VY(j),VZ(j)
    OriginalParticleId(j)=j
  enddo
  close(31)
  allocate(Mvir(1),Rvir(1),Mtotal(1),Xoff(1),xMaxx(1),yMaxx(1),zMaxx(1), &
           VxMaxx(1),VyMaxx(1),VzMaxx(1),EpotM(1),EkinM(1),LambdaM(1),RadRms(1), &
           VmaxM(1),RmaxM(1),Xax(1),Yax(1),Zax(1),Axba(1),Axca(1))
  allocate(Lst(Np),Label(Nmx:Nbx,Nmy:Nby,Nmz:Nbz))
  call List
  xMaxx=centre(1);yMaxx=centre(2);zMaxx=centre(3)
  call BdmHaloMembershipInit
  call GetHalo(centre(1),centre(2),centre(3),0.,0.,0.,1_8)
  write(*,'(a,7es25.16)')'SCALES ',dble(Box),dble(MassOne),dble(Om0), &
      dble(AEXPN),dble(Ovdens),dble(Cell),dble(Rext)
  write(*,'(a,i6,18es25.16)')'VALUES ',HaloStatus(1),Mvir(1),Mtotal(1),Rvir(1),EkinM(1),EpotM(1), &
      VmaxM(1),RmaxM(1),VxMaxx(1),VyMaxx(1),VzMaxx(1),RadRms(1),Xoff(1),LambdaM(1), &
      Axba(1),Axca(1),Xax(1),Yax(1),Zax(1)
  if(allocated(BoundParticleIds(1)%ids))then
    write(*,'(a,*(i0,1x))')'IDS ',BoundParticleIds(1)%ids
  else
    write(*,*)'IDS '
  endif
  call WriteFiles
end program normalization_cases
