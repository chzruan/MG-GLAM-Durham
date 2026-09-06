! Diagnostic wrapper only; the compiled finder arithmetic is unchanged.
module BdmReplayControl
  use LinkerList
  use Density, only: DENSIT
  use omp_lib, only: omp_set_num_threads,omp_set_dynamic,omp_get_wtime
  use iso_fortran_env, only: int32,int64
  use iso_c_binding, only: c_int,c_char,c_null_char
  implicit none
contains
  subroutine ReplaySnapshotGate(step)
    integer,intent(in) :: step
    integer :: unit,io
    real :: local_aexp0,local_partw,local_au0
    character(64) :: path
    if(step<0.or.step>9999)error stop 'Replay step must fit the PM filename'
    write(path,'(a,i4.4,a)')'PMcrd.',step,'.DAT'
    open(newunit=unit,file=trim(path),form='unformatted',status='old',action='read',iostat=io)
    if(io/=0)error stop 'Replay snapshot header is missing'
    read(unit,iostat=io)HEADER,AEXPN,local_aexp0,AMPLT,ASTEP,ISTEP,local_partw, &
      TINTG,EKIN,EKIN1,EKIN2,local_au0,AEU0,NROW,NGRID,Nrealization,Nseed,Om,OmL,hubble,Nparticles,extras
    if(io/=0)error stop 'Replay snapshot header is invalid'
    close(unit,iostat=io)
    if(io/=0)error stop 'Replay snapshot header close failed'
    if(ISTEP/=step)error stop 'Replay snapshot/header step mismatch'
    ! Gate before ReadDataPM can allocate the particle arrays.
    if(NROW<1.or.NROW>=1200)error stop 'Replay requires NROW**3 < 1200**3'
    if(Nparticles/=int(NROW,int64)**3)error stop 'Replay requires the full equal-mass PM particle set'
    if(NGRID<1)error stop 'Replay mesh dimension is invalid'
  end subroutine ReplaySnapshotGate

  subroutine ReplayMove(source,destination)
    character(*),intent(in) :: source,destination
    character(kind=c_char,len=:),allocatable :: old_name,new_name
    integer(c_int) :: result
    logical :: exists
    interface
      integer(c_int) function c_rename(old_name,new_name) bind(C,name='rename')
        import c_int,c_char
        character(c_char),intent(in) :: old_name(*),new_name(*)
      end function c_rename
    end interface
    inquire(file=destination,exist=exists)
    if(exists)error stop 'Replay refuses to replace a preserved output'
    old_name=source//c_null_char;new_name=destination//c_null_char
    result=c_rename(old_name,new_name)
    if(result/=0_c_int)error stop 'Replay could not preserve completed output'
  end subroutine ReplayMove

  subroutine ReplaySaveDensity(path,threads)
    character(*),intent(in) :: path
    integer,intent(in) :: threads
    integer :: unit,io,l
    logical :: exists
    inquire(file=path,exist=exists)
    if(exists)error stop 'Replay refuses to replace an existing density tape'
    open(newunit=unit,file=path//'.part',access='stream',form='unformatted', &
      status='new',action='write',iostat=io)
    if(io/=0)error stop 'Replay cannot stage the density tape'
    write(unit,iostat=io)int(NGRID,int64),Nparticles,int(threads,int32)
    if(io/=0)error stop 'Replay density header write failed'
    do l=1,NGRID
      write(unit,iostat=io)FI(:,:,l)
      if(io/=0)error stop 'Replay density layer write failed'
    enddo
    close(unit,iostat=io)
    if(io/=0)error stop 'Replay density close failed'
    call ReplayMove(path//'.part',path)
  end subroutine ReplaySaveDensity

  subroutine ReplayReadDensity(path)
    character(*),intent(in) :: path
    integer :: unit,io,l
    integer(int32) :: density_threads
    integer(int64) :: grid,particles,bytes
    open(newunit=unit,file=path,access='stream',form='unformatted',status='old',action='read',iostat=io)
    if(io/=0)error stop 'Replay cannot read the density tape'
    inquire(unit=unit,size=bytes)
    if(bytes/=20_int64+4_int64*NGRID*NGRID*NGRID)error stop 'Replay density tape byte length mismatch'
    read(unit,iostat=io)grid,particles,density_threads
    if(io/=0)error stop 'Replay density header read failed'
    if(grid/=NGRID.or.particles/=Nparticles.or.density_threads<1)error stop 'Replay density metadata mismatch'
    do l=1,NGRID
      read(unit,iostat=io)FI(:,:,l)
      if(io/=0)error stop 'Replay density layer read failed'
    enddo
    close(unit,iostat=io)
    if(io/=0)error stop 'Replay density close failed'
    write(*,'(a,i0)')'REPLAY DENSITY read threads=',density_threads
  end subroutine ReplayReadDensity

  subroutine ReplayAssertReleased
    if(allocated(BdmPMX).or.allocated(BdmPMY).or.allocated(BdmPMZ).or. &
       allocated(BdmPMVX).or.allocated(BdmPMVY).or.allocated(BdmPMVZ).or.BdmPMCount/=0_int64) &
      error stop 'Replay found retained original particle workspace'
    if(allocated(OriginalParticleId).or.allocated(BoundParticleIds).or.allocated(HaloStatus).or. &
       allocated(HaloUnbindingPasses).or.allocated(HaloUnbindingWork)) &
      error stop 'Replay found retained membership or telemetry workspace'
    if(allocated(Mvir).or.allocated(Rvir).or.allocated(Xoff).or.allocated(xMaxx).or. &
       allocated(yMaxx).or.allocated(zMaxx).or.allocated(VxMaxx).or.allocated(VyMaxx).or. &
       allocated(VzMaxx).or.allocated(LstMax).or.allocated(EpotM).or.allocated(EkinM).or. &
       allocated(LambdaM).or.allocated(VmaxM).or.allocated(RmaxM).or.allocated(Mtotal).or. &
       allocated(RadRms).or.allocated(Xax).or.allocated(Yax).or.allocated(Zax).or. &
       allocated(Axba).or.allocated(Axca).or.allocated(Lst).or.allocated(Label)) &
      error stop 'Replay found retained peak or linked-list workspace'
    if(allocated(MaxIndex).or.allocated(IndexDist).or.allocated(IndexLoc).or. &
       allocated(MassProf).or.allocated(DensMax).or.allocated(DistSub).or. &
       allocated(RadH1).or.allocated(MassH1).or.allocated(VrmsH1).or.allocated(VradH1).or. &
       allocated(VrmsrH1).or.allocated(RadH2).or.allocated(MassH2).or.allocated(VrmsH2).or. &
       allocated(VradH2).or.allocated(VrmsrH2).or.allocated(NbinH1).or.allocated(NbinH2)) &
      error stop 'Replay found unexpected dormant halo workspace'
    if(CataloguePublicationPending.or.len_trim(CatalogueFinalPath)/=0.or. &
       len_trim(CatalogueStagedPath)/=0)error stop 'Replay catalogue publication is incomplete'
  end subroutine ReplayAssertReleased

  subroutine RunBdmReplay(threads,passes,density_path,mode)
    integer,intent(in) :: threads,passes
    character(*),intent(in) :: density_path,mode
    integer(int32),allocatable :: saved_particles(:,:)
    real,allocatable :: saved_density(:,:,:)
    integer(int64) :: original_count,ip,changed
    integer :: pass,j,k,l
    real(8) :: started
    logical :: restored,exists
    character(32) :: tag
    if(threads<1.or.passes<1.or.passes>99)error stop 'Replay thread/pass count is invalid'
    if(mode/='read'.and.mode/='write')error stop 'Replay density mode must be read or write'
    if(NROW<1.or.NROW>=1200.or.Nparticles/=int(NROW,int64)**3)error stop 'Replay particle limit exceeded'
    if(.not.allocated(FI))error stop 'Replay requires an allocated density mesh'
    if(any(shape(FI)/=NGRID))error stop 'Replay density mesh shape mismatch'
    call ReplayAssertReleased
    do pass=1,passes
      write(tag,'(a,i0,a)')'p',pass,'.DAT'
      inquire(file=trim(tag),exist=exists)
      if(exists)error stop 'Replay pass catalogue already exists'
    enddo
    inquire(file='repair-members.bin',exist=exists)
    if(exists)error stop 'Replay membership output already exists'
    inquire(file='unbinding.bin',exist=exists)
    if(exists)error stop 'Replay telemetry output already exists'
    original_count=Nparticles
    call omp_set_dynamic(.false.)
    call omp_set_num_threads(threads)
    allocate(saved_particles(original_count,6),saved_density(NGRID,NGRID,NGRID))
!$omp parallel do private(ip)
    do ip=1,original_count
      saved_particles(ip,1)=transfer(Xpar(ip),0_int32)
      saved_particles(ip,2)=transfer(Ypar(ip),0_int32)
      saved_particles(ip,3)=transfer(Zpar(ip),0_int32)
      saved_particles(ip,4)=transfer(VX(ip),0_int32)
      saved_particles(ip,5)=transfer(VY(ip),0_int32)
      saved_particles(ip,6)=transfer(VZ(ip),0_int32)
    enddo
    if(mode=='write')then
      started=omp_get_wtime()
      call DENSIT
      write(*,'(a,i0,a,f16.6)')'REPLAY DENSITY write threads=',threads,' seconds=',omp_get_wtime()-started
      call ReplaySaveDensity(density_path,threads)
    else
      call ReplayReadDensity(density_path)
    endif
!$omp parallel do private(l)
    do l=1,NGRID
      saved_density(:,:,l)=FI(:,:,l)
    enddo
    do pass=1,passes
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
      if(changed/=0_int64)error stop 'Replay failed to restore immutable density bits'
      write(*,'(a,i0)')'REPLAY FIXED FI pass=',pass
      RepairDumpFinalPass=pass==passes
      started=omp_get_wtime()
      call BDM(0)
      write(*,'(a,i0,a,f16.6)')'REPLAY FINDER pass=',pass,' seconds=',omp_get_wtime()-started
      RepairDumpFinalPass=.false.
      if(Nparticles/=original_count.or.Np/=original_count)error stop 'Replay particle counts changed'
      if(size(Xpar,kind=int64)/=original_count.or.size(Ypar,kind=int64)/=original_count.or. &
         size(Zpar,kind=int64)/=original_count.or.size(VX,kind=int64)/=original_count.or. &
         size(VY,kind=int64)/=original_count.or.size(VZ,kind=int64)/=original_count) &
        error stop 'Replay particle array sizes changed'
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
      if(.not.restored)error stop 'Replay original particle bits changed'
      if(.not.allocated(FI))error stop 'Replay lost density allocation'
      if(any(shape(FI)/=NGRID))error stop 'Replay density allocation shape changed'
      call ReplayAssertReleased
      write(tag,'(a,i0,a)')'p',pass,'.DAT'
      call ReplayMove(trim(outputName),trim(tag))
      write(*,'(a,i0)')'REPLAY RESTORED pass=',pass
    enddo
    deallocate(saved_particles,saved_density)
    print '(a)','REPLAY COMPLETE'
  end subroutine RunBdmReplay
end module BdmReplayControl
