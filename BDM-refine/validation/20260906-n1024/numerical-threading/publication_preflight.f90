! Tiny preflight; no particle or mesh allocations and no simulation.
program PublicationPreflight
  use BdmThreadControl
  implicit none
  integer :: unit,io
  character(40) :: line
  Nrealization=1
  call ReadParameters(157)
  write(12,'(a)')'complete fixture'
  call PublishCatalogue
  if(len_trim(CatalogueFinalPath)/=0)error stop 'Unexpected publication lifecycle'
  call MoveThreadOutput(trim(outputName),'preserved.DAT')
  open(newunit=unit,file='preserved.DAT',status='old',iostat=io)
  if(io/=0)error stop 'Catalogue preservation failed'
  read(unit,'(a)')line
  close(unit)
  if(trim(line)/='complete fixture')error stop 'Catalogue bytes changed'
  print *, 'PUBLICATION PREFLIGHT PASSED'
end program PublicationPreflight
