program BdmThreadProbe
  use BdmThreadControl
  implicit none
  integer :: snapshot_step,high_threads,low_threads
  read(*,*)snapshot_step,high_threads,low_threads
  call omp_set_num_threads(high_threads)
  call ReadDataPM(snapshot_step,'')
  allocate(FI(NGRID,NGRID,NGRID))
  call RunBdmThreadControls(high_threads,low_threads)
end program BdmThreadProbe
