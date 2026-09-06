! Read-only diagnostic around the unchanged production BDM(0) interface.
! Link this module and the separate entry against a frozen native build.
module BdmThreadControl
  use LinkerList
  use Density, only: DENSIT
  use omp_lib, only: omp_set_num_threads,omp_set_dynamic,omp_get_wtime
  use iso_fortran_env, only: int32,int64
  use iso_c_binding, only: c_int,c_char,c_null_char
  implicit none
contains
  subroutine RunBdmThreadControls(high_threads,low_threads)
    integer,intent(in) :: high_threads,low_threads
    integer(int32),allocatable :: saved_particles(:,:)
    real,allocatable :: saved_density(:,:,:)
    integer(int64) :: original_count,ip,changed
    integer :: field,pass,passes,threads,density_threads,j,k,l,unit,io
    logical :: restored
    real(8) :: max_difference,started
    character(32) :: density_tag,tag
    if(high_threads<low_threads.or.low_threads<1)error stop 'Invalid diagnostic thread counts'
    if(Nparticles<=0_int64.or.Nparticles>=1200_int64**3)error stop 'Particle limit exceeded'
    if(.not.allocated(FI))error stop 'Allocate FI before the threading control'
    original_count=Nparticles
    allocate(saved_particles(original_count,6),saved_density(NGRID,NGRID,NGRID))
    call omp_set_dynamic(.false.)
    call omp_set_num_threads(high_threads)
!$omp parallel do private(ip)
    do ip=1,original_count
      saved_particles(ip,1)=transfer(Xpar(ip),0_int32)
      saved_particles(ip,2)=transfer(Ypar(ip),0_int32)
      saved_particles(ip,3)=transfer(Zpar(ip),0_int32)
      saved_particles(ip,4)=transfer(VX(ip),0_int32)
      saved_particles(ip,5)=transfer(VY(ip),0_int32)
      saved_particles(ip,6)=transfer(VZ(ip),0_int32)
    enddo
    do field=1,2
      density_threads=high_threads
      if(field==2)density_threads=low_threads
      call omp_set_num_threads(density_threads)
      write(density_tag,'(a,i0)')'d',density_threads
      started=omp_get_wtime()
      call DENSIT
      print *, 'PROBE DENSITY COMPLETE ',trim(density_tag),omp_get_wtime()-started
      if(field==2)then
        changed=0_int64;max_difference=0.d0
!$omp parallel do private(j,k,l) reduction(+:changed) reduction(max:max_difference)
        do l=1,NGRID
          do k=1,NGRID
            do j=1,NGRID
              if(transfer(FI(j,k,l),0_int32)/=transfer(saved_density(j,k,l),0_int32))changed=changed+1
              max_difference=max(max_difference,abs(dble(FI(j,k,l))-dble(saved_density(j,k,l))))
            enddo
          enddo
        enddo
        write(*,'(a,i0,a,es24.16)')'PROBE DENSITY DIFFERENCES cells=',changed,' maximum_absolute=',max_difference
      endif
      ! The file records the exact realization, independent of the later
      ! fingerprint calculation. The in-memory copy is the immutable control.
      open(newunit=unit,file=trim(density_tag)//'.density.bin',access='stream', &
           form='unformatted',status='new',action='write',iostat=io)
      if(io/=0)error stop 'Cannot create density tape'
      write(unit)int(NGRID,int64),original_count,int(density_threads,int32)
      do l=1,NGRID
        write(unit)FI(:,:,l)
      enddo
      close(unit,iostat=io)
      if(io/=0)error stop 'Cannot close density tape'
!$omp parallel do private(l)
      do l=1,NGRID
        saved_density(:,:,l)=FI(:,:,l)
      enddo
      call DumpThreadPeaks(trim(density_tag))
      passes=4
      if(field==2)passes=2
      do pass=1,passes
        threads=low_threads
        if(mod(pass,2)==0)threads=high_threads
        call omp_set_num_threads(threads)
        write(tag,'(a,a,i0,a,i0)')trim(density_tag),'-t',threads,'-p',pass
!$omp parallel do private(l)
        do l=1,NGRID
          FI(:,:,l)=saved_density(:,:,l)
        enddo
        changed=0_int64
!$omp parallel do private(j,k,l) reduction(+:changed)
        do l=1,NGRID
          do k=1,NGRID
            do j=1,NGRID
              if(transfer(FI(j,k,l),0_int32)/=transfer(saved_density(j,k,l),0_int32))changed=changed+1
            enddo
          enddo
        enddo
        if(changed/=0)error stop 'Restored density differs from immutable reference'
        print *, 'PROBE FIXED DENSITY VERIFIED ',trim(tag),' bit_differences=',changed
        started=omp_get_wtime()
        call BDM(0)
        write(*,'(3a,f16.6)')'PROBE FINDER COMPLETE ',trim(tag),' seconds=',omp_get_wtime()-started
        if(Nparticles/=original_count.or.Np/=original_count)error stop 'BDM changed particle counts'
        if(size(Xpar,kind=int64)/=original_count.or.size(Ypar,kind=int64)/=original_count.or. &
           size(Zpar,kind=int64)/=original_count.or.size(VX,kind=int64)/=original_count.or. &
           size(VY,kind=int64)/=original_count.or.size(VZ,kind=int64)/=original_count) &
           error stop 'BDM changed particle array sizes'
        restored=.true.
!$omp parallel do private(ip) reduction(.and.:restored)
        do ip=1,original_count
          restored=restored.and.saved_particles(ip,1)==transfer(Xpar(ip),0_int32)
          restored=restored.and.saved_particles(ip,2)==transfer(Ypar(ip),0_int32)
          restored=restored.and.saved_particles(ip,3)==transfer(Zpar(ip),0_int32)
          restored=restored.and.saved_particles(ip,4)==transfer(VX(ip),0_int32)
          restored=restored.and.saved_particles(ip,5)==transfer(VY(ip),0_int32)
          restored=restored.and.saved_particles(ip,6)==transfer(VZ(ip),0_int32)
        enddo
        if(.not.restored)error stop 'BDM changed original particle bits'
        if(allocated(OriginalParticleId).or.allocated(BoundParticleIds).or.allocated(HaloStatus)) &
          error stop 'BDM retained its particle or halo workspace'
        call MoveThreadOutput(trim(CatalogueFinalPath),trim(tag)//'.DAT')
        print *, 'PROBE PARTICLE RESTORATION VERIFIED ',trim(tag)
      enddo
    enddo
    deallocate(saved_particles,saved_density)
    print *, 'BDM FIXED DENSITY THREAD CONTROL COMPLETE'
  end subroutine RunBdmThreadControls

  subroutine DumpThreadPeaks(tag)
    character(*),intent(in) :: tag
    integer :: unit,ip,io
    real :: seed_cell,seed_radius
    ! FindMaxima is read-only on FI. Its allocated arrays are released before
    ! BDM recomputes the same list, with the configuration used by this pilot.
    iVirial=1
    close(13)
    open(13,file=tag//'.peaks.log',status='new')
    call FindMaxima
    seed_cell=2.*(Box/NGRID)
    open(newunit=unit,file=tag//'.peaks.bin',access='stream',form='unformatted',status='new',iostat=io)
    if(io/=0)error stop 'Cannot create peak tape'
    write(unit)int(Nmaxima,int64),int(NGRID,int32),Box
    do ip=1,Nmaxima
      seed_radius=min(seed_cell*max(0.5,min(2.,log10(max(Xoff(ip),0.)+10.)/2.)),nearest(0.5*Box,-1.))
      write(unit)int(ip,int64),xMaxx(ip),yMaxx(ip),zMaxx(ip),Xoff(ip),seed_radius
    enddo
    close(unit,iostat=io)
    if(io/=0)error stop 'Cannot close peak tape'
    call ReleaseMaxima
  end subroutine DumpThreadPeaks

  subroutine MoveThreadOutput(source,destination)
    character(*),intent(in) :: source,destination
    character(kind=c_char,len=:),allocatable :: old_name,new_name
    integer(c_int) :: status
    interface
      integer(c_int) function c_rename(old_name,new_name) bind(C,name='rename')
        import c_int,c_char
        character(c_char),intent(in) :: old_name(*),new_name(*)
      end function c_rename
    end interface
    old_name=source//c_null_char;new_name=destination//c_null_char
    status=c_rename(old_name,new_name)
    if(status/=0_c_int)error stop 'Cannot preserve diagnostic catalogue'
  end subroutine MoveThreadOutput
end module BdmThreadControl
