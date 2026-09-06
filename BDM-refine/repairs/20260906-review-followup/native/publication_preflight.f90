! Native production writer and teardown preflight; at most 64 synthetic rows.
program NativePublicationPreflight
  use BdmReplayControl
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  implicit none
  integer :: j,ip,unit,io
  character(32) :: mode
  call get_command_argument(1,mode)
  NROW=4;NGRID=16;Nparticles=64;Np=Nparticles
  Box=32.;Om=.3;OmL=.7;hubble=.7;AEXPN=1.;ASTEP=.004;ISTEP=1;Nrealization=1
  HEADER='Native follow-up preflight'
  call omp_set_num_threads(1)
  if(trim(mode)=='leak-passes')then
    allocate(HaloUnbindingPasses(1))
    call ReplayAssertReleased
    error stop 'Preflight failed to detect retained pass counters'
  endif
  if(trim(mode)=='leak-work')then
    allocate(HaloUnbindingWork(1))
    call ReplayAssertReleased
    error stop 'Preflight failed to detect retained work counters'
  endif
  call ReadParameters(1)
  call PrepareParticleSearch
  call SetParameters
  allocate(FI(NGRID,NGRID,NGRID));FI=0.
  if(trim(mode)/='empty')FI(4,4,4)=1000.
  if(trim(mode)=='invalid')FI(12,12,12)=1000.
  call FindMaxima
  if(trim(mode)=='empty')then
    if(Nmaxima/=0)error stop 'Preflight expected no candidate'
    call WriteFiles
    call ReleaseMaxima
  else
    if(Nmaxima/=merge(2,1,trim(mode)=='invalid'))error stop 'Preflight candidate count mismatch'
    call InitMaxima
    call BdmHaloMembershipInit
    do ip=1,Nmaxima
      allocate(BoundParticleIds(ip)%ids(20))
      BoundParticleIds(ip)%ids=[(int(j+20*(ip-1),8),j=1,20)]
    enddo
    Mvir=20.*MassOne;Mtotal=Mvir;Rvir=.2
    xMaxx=6.;yMaxx=6.;zMaxx=6.;VxMaxx=0.;VyMaxx=0.;VzMaxx=0.
    EkinM=1.e14;EpotM=1.e15;VmaxM=0.;RmaxM=0.;Xoff=0.;LambdaM=0.
    RadRms=0.;Axba=0.;Axca=0.;Xax=1.;Yax=0.;Zax=0.
    HaloUnbindingPasses=7;HaloUnbindingWork=1234_8
    ! The first row has been written when validation rejects the second.
    if(trim(mode)=='invalid')EkinM(2)=ieee_value(1.,ieee_quiet_nan)
    call WriteFiles
    call ReleaseMaxima
  endif
  call ReplayAssertReleased
  call ReplayMove(trim(outputName),'preserved.DAT')
  open(newunit=unit,file='preserved.DAT',status='old',action='read',iostat=io)
  if(io/=0)error stop 'Preflight publication preservation failed'
  close(unit)
  print '(a)','NATIVE PUBLICATION PREFLIGHT COMPLETE'
end program NativePublicationPreflight
