! Controlled tests of extracted production configuration, peak and empty BDM paths.
program peaks_cases
  use LinkerList
  implicit none
  character(40) :: test_case
  real, allocatable :: original_density(:,:,:),original_particles(:,:)
  integer :: q
  call get_command_argument(1,test_case)
  NGRID=16; NROW=8; Box=32.; Om=.3; OmL=.7; AEXPN=1.
  Nparticles=17; Np=Nparticles
  open(13,status='scratch')
  select case(trim(test_case))
  case('peaks')
    iVirial=2
    allocate(FI(NGRID,NGRID,NGRID),original_density(NGRID,NGRID,NGRID))
    open(41,file='density.bin',access='stream',form='unformatted',status='old')
    read(41) FI
    close(41)
    original_density=FI
    call FindMaxima
    if (any(FI /= original_density)) error stop 'FindMaxima modified its density input'
    open(41,file='peaks.bin',access='stream',form='unformatted',status='replace')
    write(41) Nmaxima,Ovdens,xMaxx,yMaxx,zMaxx,Xoff
    close(41)
    if (size(xMaxx) /= Nmaxima) error stop 'inexact maximum allocation'
    call ReleaseMaxima
  case('configuration')
    call ReadParameters(1)
    call SetParameters
    close(12)
    allocate(FI(NGRID,NGRID,NGRID)); FI=0.
    FI(4,4,4)=1000.; FI(12,12,12)=100.
    call FindMaxima
    open(41,file='settings.txt',status='replace')
    write(41,*) iVirial,NradP,dLogR,dLogP,MassMin,Rext,SlopeR,Ovdens,Nmaxima
    close(41)
    call ReleaseMaxima
  case('empty_bdm')
    iVirial=-1; Ovdens=-999.
    allocate(FI(NGRID,NGRID,NGRID)); FI=0.
    allocate(Xpar(Np),Ypar(Np),Zpar(Np),VX(Np),VY(Np),VZ(Np))
    allocate(original_particles(Np,6))
    do q=1,int(Np)
      Xpar(q)=1.+mod(q*37,15)+.12345; Ypar(q)=Xpar(q); Zpar(q)=Xpar(q)
      VX(q)=q*.1234567; VY(q)=VX(q); VZ(q)=VX(q)
    end do
    original_particles(:,1)=Xpar; original_particles(:,2)=Ypar; original_particles(:,3)=Zpar
    original_particles(:,4)=VX; original_particles(:,5)=VY; original_particles(:,6)=VZ
    ! DENSIT asserts that configuration and the header precede the density pass.
    call BDM(1)
    call BDM(0)
    if (.not.allocated(FI)) error stop 'empty BDM removed integrator density'
    if (any(FI /= 0.)) error stop 'empty BDM changed integrator density'
    if (allocated(Mvir).or.allocated(xMaxx).or.allocated(Axba)) error stop 'empty BDM leaked maxima'
    if (Np /= 17.or.Nparticles /= 17) error stop 'empty BDM changed particle count'
    if (any(Xpar /= original_particles(:,1)).or.any(Ypar /= original_particles(:,2)).or. &
        any(Zpar /= original_particles(:,3)).or.any(VX /= original_particles(:,4)).or. &
        any(VY /= original_particles(:,5)).or.any(VZ /= original_particles(:,6))) &
        error stop 'empty BDM modified integrator particles'
  case default
    error stop 'unknown peak repair test'
  end select
  print *, 'PEAKS TEST COMPLETE'
end program peaks_cases
