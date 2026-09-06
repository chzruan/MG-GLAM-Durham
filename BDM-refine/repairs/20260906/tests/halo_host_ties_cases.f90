program halo_host_ties_cases
  use LinkerList
  implicit none
  integer :: nactive,q,candidate,nmembers,i
  read(*,*)Nmaxima,nactive,Box,Cell,MassOne
  allocate(Mvir(Nmaxima),Rvir(Nmaxima),Mtotal(Nmaxima),xMaxx(Nmaxima),yMaxx(Nmaxima),zMaxx(Nmaxima))
  allocate(BoundParticleIds(Nmaxima),Lst(Nmaxima))
  Mvir=0.;Rvir=0.;Mtotal=0.;xMaxx=0.;yMaxx=0.;zMaxx=0.
  do q=1,nactive
    read(*,*)candidate,Mvir(candidate),Rvir(candidate),nmembers
    Mtotal(candidate)=Mvir(candidate)
    read(*,*)xMaxx(candidate),yMaxx(candidate),zMaxx(candidate)
    allocate(BoundParticleIds(candidate)%ids(nmembers))
    read(*,*)BoundParticleIds(candidate)%ids
  enddo
  Nmx=-2;Nmy=-2;Nmz=-2
  Nbx=ceiling(dble(Box)/dble(Cell))+2;Nby=Nbx;Nbz=Nbx
  allocate(Label(Nmx:Nbx,Nmy:Nby,Nmz:Nbz))
  open(13,status='scratch')
  call ListMaxima
  call RemoveDuplicates
  do i=1,Nmaxima
    if(Mvir(i)>MassOne)write(*,'(a,i12)')'KEEP ',i
  enddo
end program halo_host_ties_cases
