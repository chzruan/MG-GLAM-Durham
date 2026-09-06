!----------------------------------------------------------------------------
!
!          BDM halo finder   A.Klypin 2010
!
! Contains:     

! Legacy numerical-field diagnostic retained for historical catalogue replay.
! Production duplicate removal compares exact original bound-particle IDs.
! Equal particle counts alone do not establish equal particle membership.
module BdmDuplicateRules
  implicit none
  real*8, parameter :: DuplicateRadius=0.2d0, DuplicateSpeed=5.d0
  real*8, parameter :: DuplicateLogMass=0.005d0
contains
  pure logical function StrictDuplicate(distance2, mass1, mass2, count1, count2, &
                                         total1, total2, speed2)
    real*8, intent(in) :: distance2,mass1,mass2,count1,count2,total1,total2,speed2
    StrictDuplicate = .false.
    if (mass1 <= 0.d0.or.total1 <= 0.d0.or.total2 <= 0.d0) return
    if (mass1 /= mass2.or.count1 /= count2) return
    if (distance2 >= DuplicateRadius**2.or.speed2 >= DuplicateSpeed**2) return
    StrictDuplicate = abs(log10(total1)-log10(total2)) < DuplicateLogMass
  end function StrictDuplicate

  integer function DuplicateRoot(parent, candidate) result(root)
    integer, intent(inout) :: parent(:)
    integer, intent(in) :: candidate
    root=candidate
    do while(parent(root) /= root)
      parent(root)=parent(parent(root))
      root=parent(root)
    end do
  end function DuplicateRoot
end module BdmDuplicateRules

Module  Structures 

Integer*4,  PARAMETER  ::                         & 
                     Nrad     =   200          ! Number of shells for halo potential
Real*4        ::                    & 
                     dLogR,         &         ! size of log binning for profiles
                     dLogP,         &         ! size of log binning for potential
                     MaxMemory=500, &         ! limit on memory for the run
                     MassMin      , &         ! minimum halo mass
                     SlopeR       , &         ! slope for extra radius for resolution correction
                     Rext                     ! extra radius shift at M=1e15
Integer*4  ::                      &
                     NradP,        &          ! Number of shells for halo profiles
                     iVirial=1                ! 0=200 critical, 1=virial, 2=200 matter, 3=Abacus

Real*4 ::            TotalMemory=0.2, t0                           ! current memory
Real*4 ::            Xleft,Xright,Yleft,Yright,Zleft,Zright,dBuffer    ! boundaries of domain
Character*120 ::     CatshortName,CatalogName,    &                ! names of halo catalogs 
                     outputName                                    ! dump file
! Publish a complete catalogue atomically; an interrupted finder may only leave
! a hidden staged file, never a partially replaced final catalogue.
character(256) :: CatalogueFinalPath='',CatalogueStagedPath=''
logical :: CataloguePublicationPending=.false.

Real*4     ::        Om0,Ovdens      ! cosmology
Real*4     ::        MassOne                                 !  current simulation
Integer*8 ::         Np                                          ! Number of particles in domain
! Analysis owns a separate particle workspace. These arrays retain the exact
! simulation state until RemoveBuffer moves it back without inverse arithmetic.
Real*4, allocatable :: BdmPMX(:),BdmPMY(:),BdmPMZ(:),BdmPMVX(:),BdmPMVY(:),BdmPMVZ(:)
Integer*8 :: BdmPMCount=0_8
Integer*8, allocatable :: OriginalParticleId(:) ! snapshot row; shared by periodic images
Real*4 :: HaloSearchRadius=0.,ParticleSearchRadius=0. ! maximum SO / aperture radii
                     ! -------------------------  Maxima ----------------------
Integer*4 ::                     Nmaxima
! Exact final bound membership is retained until the duplicate pass. IDs refer
! to original particles, so a periodic image cannot create a different identity.
type BdmHaloMembership
  integer*8, allocatable :: ids(:)
end type BdmHaloMembership
type(BdmHaloMembership), allocatable :: BoundParticleIds(:)
integer, parameter :: HaloTooFewParticles=1, HaloSearchTruncated=2, &
                      HaloUnresolvedVmax=4, HaloNoBoundParticles=8, &
                      HaloSingularCentre=16
integer, allocatable :: HaloStatus(:)
Real*4    ::                     RadMax
Real*4,        ALLOCATABLE,   DIMENSION(:) ::                        &    ! maxima of density
                                 Mvir,Rvir,Mtotal,VmaxM,RmaxM,       &
                                 xMaxx,yMaxx,zMaxx,                     &
                                 VxMaxx,VyMaxx,VzMaxx,                  &
                                 Xoff,EpotM,EkinM,LambdaM,           & 
                                 RadRms,Axba,Axca,Xax,Yax,Zax
Integer*4,   ALLOCATABLE,   DIMENSION(:) ::           & 
                                 LstMax,                             &   ! linker-list of maxima
                                 MaxIndex,                           &   !  =0 for distinct, =-2 for sub
                                 IndexDist, IndexLoc                
Real*4,        ALLOCATABLE,   DIMENSION(:,:) ::                      &    
                                 MassProf,                           &   ! mass profile
                                 DensMax,DistSub
                     ! ------------------------- Halos ----------------------
Integer*4 ::                     Nhalo
Real*4,        ALLOCATABLE,   DIMENSION(:,:) ::                      &   
                                 RadH1,MassH1,VrmsH1,VradH1,VrmsrH1, &
                                 RadH2,MassH2,VrmsH2,VradH2,VrmsrH2
Integer*4,     ALLOCATABLE,   DIMENSION(:,:) ::                      &   
                                 NbinH1,NbinH2
                     !-------------------------  Main linker-list --------------
Integer*8,  ALLOCATABLE :: Lst(:)                                   ! Linker list
Integer*8,  ALLOCATABLE :: Label(:,:,:)                         ! Head of zero-level LL 
Integer*4 ::                    Nmx,Nmy,Nmz,Nbx,Nby,Nbz         ! limits for linker-list
Real*4    ::                    Cell,Roptimal                   ! Size in Mpch of a linker-list cell

end Module Structures
!----------------------------------------------------------------------------
!
Module  LinkerList
  use Structures
  use Tools
  use BdmDuplicateRules

Contains
!----------------------------------------------------------
!                      
!                    
      SUBROUTINE BDM(mDENSIT)
!----------------------------------------------------------
   Use Density

      Integer*8  ::  idummy,ip,NpPM
      Integer*4  ::  jdummy
      integer*4  ::  OMP_GET_MAX_THREADS,OMP_GET_THREAD_NUM
      character*80 :: Path
      logical    ::  op

      NpPM = Nparticles
      Np   = Nparticles
      t0 = seconds()

      iThreads = OMP_GET_MAX_THREADS()
      write (*,'(a,i4)')        ' Number of threads      =',iThreads
      write (*,'(a,i4)')        ' mDENSIT                =',mDENSIT
      Path =''
      tstart = seconds()
      Call ReadParameters(ISTEP)
      Call PrepareParticleSearch
      Call SetParameters
      If(mDENSIT==1)Call DENSIT                 ! density on original Ng mesh
      Call FindMaxima

      ! Empty analysis is a valid catalogue, including during an inline call.
      ! No particle scaling/buffering has occurred and the PM density stays allocated.
      if (Nmaxima == 0) then
         Call WriteFiles
         Call ReleaseMaxima
         return
      end if

      myMemory =Memory(-1_8*NGRID*NGRID*NGRID)
      DeAllocate (FI)
      write(*,'(a,i11)') ' Go to RescaleCoords    : ',Nparticles
      Call RescaleCoords(1)  
      tfinish = seconds()

      write(*,'(a,i11)') ' Go to AddBuffer        : ',Nparticles
      write(13,'(10x,a,T50,2f10.2)') ' time for Dens+Maxima  =',tfinish-tstart,tfinish-t0
      write(*,'(10x,a,T50,2f10.2)')  ' time for Dens+Maxima  =',tfinish-tstart,tfinish-t0
      tstart = seconds()
      Call AddBuffer
      tfinish = seconds()
      write(13,'(10x,a,T50,2f10.2)') ' time for AddBuffer  =',tfinish-tstart,tfinish-t0
      write(*,'(10x,a,T50,2f10.2)')  ' time for AddBuffer  =',tfinish-tstart,tfinish-t0
      
      Call  SizeList
            ALLOCATE (Lst(Np),Label(Nmx:Nbx,Nmy:Nby,Nmz:Nbz))
            myMemory= Memory(2_8*(Np+(Nbx-Nmx+1_8)*(Nby-Nmy+1_8)*(Nbz-Nmz+1_8)))
      tstart = seconds()
      Call List
      tfinish = seconds()
      write(13,'(10x,a,T50,2f10.2)') ' time for List  =',tfinish-tstart,tfinish-t0      
      write(*,'(10x,a,T50,2f10.2)') ' time for List  =',tfinish-tstart,tfinish-t0      
  
      write(*,*) ' Go to FindDistinctCandidates: ',Nparticles
      Call FindDistinctCandidates    !  get estimates of Mvir and Rvir for all maxima


      write(*,'(a,i11)') ' Go to ParametersDistinct : ',Nparticles
      tstart = seconds()
        Call ParametersDistinct
      tfinish = seconds()
      write(13,'(10x,a,T50,2f10.2)') ' time for ParametersDistinct  =',tfinish-tstart,tfinish-t0      
      write(*,'(10x,a,T50,2f10.2)')  ' time for ParametersDistinct  =',tfinish-tstart,tfinish-t0      
!        Call FindSubs
!        Call ParametersSubs
              write(*,*) ' goto RemoveCloseMaxima'
              tfinish = seconds()
            DEALLOCATE (Lst,Label)
            myMemory= Memory(-2_8*(Np+(Nbx-Nmx+1_8)*(Nby-Nmy+1_8)*(Nbz-Nmz+1_8)))
              
      Call  SizeListMaxima
            ALLOCATE (Lst(Nmaxima),Label(Nmx:Nbx,Nmy:Nby,Nmz:Nbz))
            myMemory= Memory(2_8*(Nmaxima+(Nbx-Nmx+1_8)*(Nby-Nmy+1_8)*(Nbz-Nmz+1_8)))
      Call ListMaxima     
!      Call RemoveDuplicatesSimple
      Call RemoveDuplicates

      write(*,*) ' Go to Write Catshort: ',Nparticles

      tstart = seconds()
      Call WriteFiles
      tfinish = seconds()
      write(13,'(10x,a,T50,2f10.2)') ' time for WriteFiles  =',tfinish-tstart,tfinish-t0      
      write(*,'(10x,a,T50,2f10.2)')  ' time for WriteFiles  =',tfinish-tstart,tfinish-t0      
      close(12)
       !-------- deallocate temporary arrays
            DEALLOCATE (Lst,Label)
            myMemory= Memory(-2_8*(Nmaxima+(Nbx-Nmx+1_8)*(Nby-Nmy+1_8)*(Nbz-Nmz+1_8)))
      
       Call ReleaseMaxima
      !-------- restore PM structure
       Call RemoveBuffer(NpPM)
       myMemory =Memory(1_8*NGRID*NGRID*NGRID)
       Allocate (FI(NGRID,NGRID,NGRID))
             
           end SUBROUTINE BDM

      SUBROUTINE ReleaseMaxima
      implicit none
      real :: released
      DEALLOCATE(Mvir,Rvir,Xoff,xMaxx,yMaxx,zMaxx,VxMaxx,VyMaxx,VzMaxx)
      DEALLOCATE(LstMax,EpotM,EkinM,LambdaM,VmaxM,RmaxM,Mtotal,RadRms)
      DEALLOCATE(Xax,Yax,Zax,Axba,Axca)
      if(allocated(BoundParticleIds))deallocate(BoundParticleIds)
      if(allocated(HaloStatus))deallocate(HaloStatus)
      released=Memory(-22_8*Nmaxima)
      end SUBROUTINE ReleaseMaxima

!----------------------------------------------------------
!               Read/create  configuration file BDM.config       
!                    
      SUBROUTINE ReadParameters(jStep)
      use, intrinsic :: iso_fortran_env, only: iostat_end,iostat_eor
      implicit none
      integer, intent(in) :: jStep
      integer :: unit, io, i, equals, comment, line_number, nne_legacy
      character(1024) :: line, key, value
      character(10) :: CatLabel
      logical :: FileExists

      NradP=30; iVirial=1; dLogR=0.02; dLogP=0.02
      dBuffer=5.; SlopeR=0.20; Rext=0.15; MassMin=2.5e12
      inquire(file='BDM.config',exist=FileExists)
      if (FileExists) then
         open(newunit=unit,file='BDM.config',status='old',action='read',iostat=io)
         if (io /= 0) call ConfigurationError(0,'cannot open BDM.config')
         line_number=0
         do
            read(unit,'(a)',advance='no',iostat=io) line
            if (io == iostat_end) exit
            line_number=line_number+1
            if (io == 0) call ConfigurationError(line_number,'configuration line is too long')
            if (io /= iostat_eor) call ConfigurationError(line_number,'cannot read configuration')
            do i=1,len_trim(line)
               if (line(i:i) == achar(9)) line(i:i)=' '
            end do
            comment=index(line,'!')
            if (comment > 0) line(comment:)=' '
            if (len_trim(line) == 0) cycle
            equals=index(line,'=')
            if (equals <= 1) call ConfigurationError(line_number,'expected name = value')
            key=adjustl(line(:equals-1))
            value=adjustl(line(equals+1:))
            if (len_trim(value) == 0) call ConfigurationError(line_number,'missing value')
            ! A parameter has one scalar value. Do not silently accept trailing tokens.
            if (scan(trim(value),' '//achar(9)//',/') /= 0) &
               call ConfigurationError(line_number,'expected one scalar value')
            do i=1,len_trim(key)
               if (key(i:i) >= 'A'.and.key(i:i) <= 'Z') key(i:i)=achar(iachar(key(i:i))+32)
            end do
            select case(trim(key))
            case('nradp')
               read(value,*,iostat=io) NradP
            case('ivirial')
               read(value,*,iostat=io) iVirial
            case('dlogr')
               read(value,*,iostat=io) dLogR
            case('dlogp')
               read(value,*,iostat=io) dLogP
            case('rext','rextr','rextern')
               read(value,*,iostat=io) Rext
            case('sloper','slope')
               read(value,*,iostat=io) SlopeR
            case('massmin','minmass')
               read(value,*,iostat=io) MassMin
            case('nne')
               ! This historical option never affected the active finder.
               read(value,*,iostat=io) nne_legacy
               if (io == 0) then
                  if (nne_legacy <= 0) call ConfigurationError(line_number,'Nne must be positive')
                  write(*,'(a)') ' BDM.config: Nne is deprecated and has no effect'
               end if
            case default
               call ConfigurationError(line_number,'unrecognized parameter: '//trim(key))
            end select
            if (io /= 0) call ConfigurationError(line_number,'invalid value for '//trim(key))
         end do
         close(unit)
      end if
      call ValidateParameters
      if (.not.FileExists) then
         open(newunit=unit,file='BDM.config',status='new',action='write',iostat=io)
         if (io /= 0) call ConfigurationError(0,'cannot create default BDM.config')
         write(unit,'(a)') '! BDM configuration; whitespace and trailing ! comments are optional'
         write(unit,10) 'iVirial',iVirial,'! 0=200 critical, 1=virial, 2=200 matter, 3=Abacus'
         write(unit,10) 'NradP',NradP,'! Number of shells for halo profiles'
         write(unit,20) 'Rext',Rext,'! Extra radius shift at m=1e15'
         write(unit,20) 'SlopeR',SlopeR,'! Slope for extra radius shift'
         write(unit,20) 'MassMin',MassMin,'! Minimum halo mass'
         write(unit,20) 'dLogR',dLogR,'! Log bin size for potential'
         write(unit,20) 'dLogP',dLogP,'! Log bin size for profiles'
         close(unit)
      end if
      write(*,'(/a)') ' ------ current set of parameters:'
      write(*,10) 'iVirial',iVirial,'! 0=200 critical, 1=virial, 2=200 matter, 3=Abacus'
      write(*,10) 'NradP',NradP,'! Number of shells for halo profiles'
      write(*,20) 'Rext',Rext,'! Extra radius shift at m=1e15'
      write(*,20) 'SlopeR',SlopeR,'! Slope for extra radius shift'
      write(*,20) 'MassMin',MassMin,'! Minimum halo mass'
      write(*,20) 'dLogR',dLogR,'! Log bin size for potential'
      write(*,20) 'dLogP',dLogP,'! Log bin size for profiles'
10    format(10x,a,T20,' = ',i6,T40,a)
20    format(10x,a,T20,' = ',1p,g10.3,T40,a)

      ! Only validated configuration may create or replace catalogue outputs.
      write(outputName,'(a,i4.4,a)') 'CATALOGS/outputB.',jStep,'.dat'
      close(13)
      open(13,file=trim(outputName),status='replace',iostat=io)
      if (io /= 0) call ConfigurationError(0,'cannot open analysis log '//trim(outputName))
      select case(iVirial)
      case(0)
         CatLabel='W.'
      case(1)
         CatLabel='V.'
      case(2)
         CatLabel='M.'
      case(3)
         CatLabel='A.'
      end select
      write(outputName,'(2a,2(i4.4,a))') 'CATALOGS/Catshort',trim(CatLabel),jStep,'.',Nrealization,'.DAT'
      call BeginCataloguePublication(trim(outputName))
      end SUBROUTINE ReadParameters

      SUBROUTINE BeginCataloguePublication(final_path)
      use, intrinsic :: iso_c_binding, only: c_int
      implicit none
      character(*), intent(in) :: final_path
      character(96) :: suffix
      character(256) :: staged_base
      integer*8 :: stamp
      integer :: attempt,io,slash
      integer(c_int) :: process_id
      logical :: exists
      interface
         function c_getpid() bind(C,name='getpid') result(pid)
            import c_int
            integer(c_int) :: pid
         end function c_getpid
      end interface
      if (CataloguePublicationPending) error stop 'BDM catalogue publication is already pending'
      if (len_trim(final_path) == 0.or.len_trim(final_path) > len(CatalogueFinalPath)) &
         error stop 'BDM catalogue path is empty or too long'
      slash=index(trim(final_path),'/',back=.true.)
      if (slash == len_trim(final_path)) error stop 'BDM catalogue path names a directory'
      if (slash > 0) then
         staged_base=final_path(:slash)//'.'//trim(final_path(slash+1:))
      else
         staged_base='.'//trim(final_path)
      end if
      process_id=c_getpid()
      call system_clock(count=stamp)
      close(12,iostat=io)
      if (io /= 0) error stop 'BDM cannot close the previous catalogue unit'
      do attempt=0,999
         write(suffix,'(a,i0,a,i0,a,i0)') '.tmp.',process_id,'.',stamp,'.',attempt
         if (len_trim(staged_base)+len_trim(suffix) > len(CatalogueStagedPath)) &
            error stop 'BDM staged catalogue path is too long'
         CatalogueStagedPath=trim(staged_base)//trim(suffix)
         ! STATUS=NEW exclusively creates the path and keeps normal umask
         ! permissions. A collision must never truncate another staged file.
         open(12,file=trim(CatalogueStagedPath),status='new',action='write',iostat=io)
         if (io == 0) then
            CatalogueFinalPath=trim(final_path)
            CataloguePublicationPending=.true.
            return
         end if
         inquire(file=trim(CatalogueStagedPath),exist=exists)
         if (.not.exists) error stop 'BDM cannot create staged catalogue'
      end do
      error stop 'BDM cannot reserve a unique staged catalogue'
      end SUBROUTINE BeginCataloguePublication

      SUBROUTINE PublishCatalogue
      use, intrinsic :: iso_c_binding, only: c_int,c_char,c_null_char
      implicit none
      integer :: io,unit
      integer(c_int) :: rename_status
      logical :: opened
      interface
         function c_rename(old_path,new_path) bind(C,name='rename') result(status)
            import c_int,c_char
            character(kind=c_char), intent(in) :: old_path(*),new_path(*)
            integer(c_int) :: status
         end function c_rename
      end interface
      if (.not.CataloguePublicationPending) error stop 'BDM catalogue publication was not started'
      unit=-1;opened=.false.
      inquire(file=trim(CatalogueStagedPath),opened=opened,number=unit,iostat=io)
      if (io /= 0.or..not.opened.or.unit /= 12) error stop 'BDM staged catalogue is not open on unit 12'
      close(12,iostat=io)
      if (io /= 0) error stop 'BDM staged catalogue could not be completely written'
      ! C rename is atomic within this directory on the supported POSIX
      ! filesystems. Only a fully written, successfully closed file is exposed.
      rename_status=c_rename(trim(CatalogueStagedPath)//c_null_char,trim(CatalogueFinalPath)//c_null_char)
      if (rename_status /= 0_c_int) error stop 'BDM cannot publish staged catalogue'
      CataloguePublicationPending=.false.
      CatalogueFinalPath='';CatalogueStagedPath=''
      end SUBROUTINE PublishCatalogue

      SUBROUTINE ConfigurationError(line_number,message)
      use, intrinsic :: iso_fortran_env, only: error_unit
      implicit none
      integer, intent(in) :: line_number
      character(*), intent(in) :: message
      write(error_unit,'(a,i0,2a)') ' BDM configuration error, line ',line_number,': ',message
      error stop 1
      end SUBROUTINE ConfigurationError

      SUBROUTINE ValidateParameters
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      implicit none
      if (iVirial < 0.or.iVirial > 3) call ConfigurationError(0,'iVirial must be 0, 1, 2 or 3')
      if (NradP < 1) call ConfigurationError(0,'NradP must be positive')
      if (.not.all(ieee_is_finite([dLogR,dLogP,MassMin,Rext,SlopeR,dBuffer]))) &
         call ConfigurationError(0,'all real parameters must be finite')
      if (dLogR <= 0..or.dLogP <= 0.) call ConfigurationError(0,'logarithmic bin widths must be positive')
      if (MassMin < 0.) call ConfigurationError(0,'MassMin must be nonnegative')
      if (Rext < 0..or.SlopeR < 0.) call ConfigurationError(0,'radius corrections must be nonnegative')
      if (dBuffer <= 0.) call ConfigurationError(0,'buffer width must be positive')
      end SUBROUTINE ValidateParameters

      SUBROUTINE SetOverdensity
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      implicit none
      real*8 :: matter_fraction, background, density, xx
      ! Retain the historical flat-Lambda convention and normalization 178.
      ! Both peak selection and halo properties use this same calculation.
      Om0=Om
      if (.not.all(ieee_is_finite([Om0,AEXPN]))) &
         call ConfigurationError(0,'Omega_m and expansion factor must be finite')
      if (Om0 <= 0..or.AEXPN <= 0.) &
         call ConfigurationError(0,'Omega_m and expansion factor must be positive')
      background=dble(Om0)+(1.d0-dble(Om0))*dble(AEXPN)**3
      if (background <= 0.d0) call ConfigurationError(0,'flat-Lambda background must be positive')
      matter_fraction=dble(Om0)/background
      xx=matter_fraction-1.d0
      select case(iVirial)
      case(0)
         density=200.d0/matter_fraction
      case(1)
         density=(178.d0+82.d0*xx-39.d0*xx**2)/matter_fraction
      case(2)
         density=200.d0
      case(3)
         density=(178.d0+82.d0*xx-39.d0*xx**2)/matter_fraction*(200.d0/178.d0)
      case default
         call ConfigurationError(0,'iVirial must be 0, 1, 2 or 3')
      end select
      if (.not.ieee_is_finite(density)) call ConfigurationError(0,'overdensity must be finite')
      if (density <= 0.d0.or.density > dble(huge(Ovdens))) &
         call ConfigurationError(0,'overdensity is outside the representable positive range')
      Ovdens=real(density)
      end SUBROUTINE SetOverdensity

!--------------------------------------------------------------
!                        virial overdensity for cosmological model
!                        at different expansion parameter AEXPN
      Function OverdenVir()
!--------------------------------------------------------------
      xx =-(1.-Om0)*AEXPN**3/(Om0+(1.-Om0)*AEXPN**3)
      OverdenVir =(178.+82.*xx-39.*xx**2)/(1.+xx)
      !write (*,*)  '      Overdensity Delta =',OverdenVir

    END Function OverdenVir
!--------------------------------------------------------------
!                   Abacus  overdensity for cosmological model
!                        at different expansion parameter AEXPN
      Function OverdenAbacus()
!--------------------------------------------------------------
      xx =-(1.-Om0)*AEXPN**3/(Om0+(1.-Om0)*AEXPN**3)
      OverdenAbacus =(178.+82.*xx-39.*xx**2)/(1.+xx)*(200./178.)
      !write (*,*)  '      Overdensity Delta =',OverdenVir

    END Function OverdenAbacus
!----------------------------------------------------------
!                      
!                    

SUBROUTINE RescaleCoords(iFlag)
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer, intent(in) :: iFlag
  integer*8 :: ip,restoreCount
  real*8 :: Xscale,Vscale
  real :: memoryUsed
  logical :: invalid

  if(iFlag/=1)then
    restoreCount=BdmPMCount
    call RemoveBuffer(restoreCount)
    return
  endif
  if(allocated(BdmPMX))error stop 'BDM particle workspace already active'
  if(Np/=Nparticles.or.Np<0_8)error stop 'BDM original particle count is inconsistent'
  if(NGRID<=0.or.Box<=0..or.AEXPN<=0.)error stop 'Invalid BDM coordinate scales'
  Xscale=dble(Box)/dble(NGRID)
  Vscale=100.d0*Xscale/dble(AEXPN)
  BdmPMCount=Nparticles
  call move_alloc(Xpar,BdmPMX); call move_alloc(Ypar,BdmPMY); call move_alloc(Zpar,BdmPMZ)
  call move_alloc(VX,BdmPMVX); call move_alloc(VY,BdmPMVY); call move_alloc(VZ,BdmPMVZ)
  if(size(BdmPMX,kind=8)/=Np.or.size(BdmPMY,kind=8)/=Np.or.size(BdmPMZ,kind=8)/=Np.or. &
     size(BdmPMVX,kind=8)/=Np.or.size(BdmPMVY,kind=8)/=Np.or.size(BdmPMVZ,kind=8)/=Np) &
       error stop 'BDM original particle array size mismatch'
  allocate(Xpar(Np),Ypar(Np),Zpar(Np),VX(Np),VY(Np),VZ(Np))
  memoryUsed=Memory(6_8*Np)
  invalid=.false.
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ip) REDUCTION(.or.:invalid)
  do ip=1,Np
    invalid=invalid.or..not.all(ieee_is_finite([BdmPMX(ip),BdmPMY(ip),BdmPMZ(ip), &
                                             BdmPMVX(ip),BdmPMVY(ip),BdmPMVZ(ip)]))
    Xpar(ip)=real(modulo((dble(BdmPMX(ip))-1.d0)*Xscale,dble(Box)))
    Ypar(ip)=real(modulo((dble(BdmPMY(ip))-1.d0)*Xscale,dble(Box)))
    Zpar(ip)=real(modulo((dble(BdmPMZ(ip))-1.d0)*Xscale,dble(Box)))
    ! Conversion to the storage precision can round a value up to Box.
    if(Xpar(ip)>=Box)Xpar(ip)=0.
    if(Ypar(ip)>=Box)Ypar(ip)=0.
    if(Zpar(ip)>=Box)Zpar(ip)=0.
    VX(ip)=real(dble(BdmPMVX(ip))*Vscale)
    VY(ip)=real(dble(BdmPMVY(ip))*Vscale)
    VZ(ip)=real(dble(BdmPMVZ(ip))*Vscale)
  enddo
  if(invalid)error stop 'Non-finite simulation particle supplied to BDM'
end SUBROUTINE RescaleCoords
!---------------------------------------------------------------------------
!                  
!                  

SUBROUTINE WriteFiles
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   integer :: ip,iHalo
   real*4 :: x,y,z,Vrms,rr,aM,Cvir,aNpart,VirRat
   real*8 :: value
   if (.not.CataloguePublicationPending) error stop 'BDM writer requires staged catalogue publication'
   iHalo=0
   if(.not.ieee_is_finite(MassOne))error stop 'BDM particle mass is nonfinite'
   if(MassOne<=0.)error stop 'BDM particle mass must be positive'
   MassMin=max(MassMin,20.*MassOne)
   if(allocated(HaloStatus))then
     write(*,'(a,5i10)') ' BDM status [few, SO cap, Vmax, no bound, singular centre]:', &
       count(btest(HaloStatus,0)),count(btest(HaloStatus,1)), &
       count(btest(HaloStatus,2)),count(btest(HaloStatus,3)),count(btest(HaloStatus,4))
   endif
   do ip=1,Nmaxima
     if(.not.all(ieee_is_finite([xMaxx(ip),yMaxx(ip),zMaxx(ip), &
         VxMaxx(ip),VyMaxx(ip),VzMaxx(ip),Mvir(ip),Mtotal(ip),Rvir(ip), &
         EkinM(ip),EpotM(ip),VmaxM(ip),RmaxM(ip),Xoff(ip),LambdaM(ip), &
         RadRms(ip),Axba(ip),Axca(ip),Xax(ip),Yax(ip),Zax(ip)]))) &
       error stop 'BDM refuses a catalogue with nonfinite halo properties'
     x=xMaxx(ip);y=yMaxx(ip);z=zMaxx(ip)
     if(x<Xleft.or.x>=Xright.or.y<Yleft.or.y>=Yright.or.z<Zleft.or.z>=Zright)cycle
     if(Mvir(ip)<MassMin.or.Mvir(ip)<=0.)cycle
     if(EkinM(ip)<0..or.EpotM(ip)<0..or.Rvir(ip)<=0.) &
       error stop 'BDM halo has invalid mass, radius or energy'
     value=sqrt(2.d0*dble(EkinM(ip))/dble(Mvir(ip)))
     if(value>dble(huge(Vrms)))error stop 'BDM velocity dispersion is not representable'
     Vrms=real(value,4);rr=1.e3*Rvir(ip);aM=Mvir(ip)
     aNpart=Mvir(ip)/MassOne
     if(allocated(BoundParticleIds))then
       if(allocated(BoundParticleIds(ip)%ids))then
         aNpart=real(size(BoundParticleIds(ip)%ids,kind=8),4)
         value=dble(size(BoundParticleIds(ip)%ids,kind=8))*dble(MassOne)
         if(abs(value-dble(Mvir(ip)))>2.d0*dble(spacing(Mvir(ip)))) &
           error stop 'BDM bound mass disagrees with particle membership'
       endif
     endif
     if(aNpart<20.)cycle
     Cvir=Concentration(aM,rr,VmaxM(ip))
     if(Cvir<0.)then
       Cvir=0.
       if(RmaxM(ip)>0.)Cvir=2.1625816*Rvir(ip)/RmaxM(ip)
     endif
     VirRat=0.
     if(EpotM(ip)>0.)then
       value=2.d0*dble(EkinM(ip))/dble(EpotM(ip))-1.d0
       if(abs(value)>dble(huge(VirRat)))error stop 'BDM virial ratio is not representable'
       VirRat=real(value,4)
     endif
     if(.not.all(ieee_is_finite([Vrms,rr,Cvir,aNpart,VirRat]))) &
       error stop 'BDM derived catalogue property is nonfinite'
     iHalo=iHalo+1
     write(12,'(3f11.4,3x,3f10.2,1p,2g12.4,g12.5,26g12.4)') &
       x,y,z,VxMaxx(ip),VyMaxx(ip),VzMaxx(ip),Mvir(ip),Mtotal(ip), &
       rr,Vrms,VmaxM(ip),iHalo,Cvir,aNpart,0,Xoff(ip),VirRat,LambdaM(ip), &
       1.e3*RadRms(ip),Axba(ip),Axca(ip),Xax(ip),Yax(ip),Zax(ip)
   enddo
   Nhalo=iHalo
   call PublishCatalogue
end SUBROUTINE WriteFiles

!---------------------------------------------------------------------------
!                  
!                  
      SUBROUTINE WriteProfiles
integer*8 :: ic   
      iHalo = 0
      Do ip=1,Nmaxima
         If(Mvir(ip)>10*MassOne)Then
            iHalo = iHalo +1
            x    = xMaxx(ip);       y =  yMaxx(ip);    z = zMaxx(ip)
         If(Mvir(ip)>100.*MassOne)Then
         Vrms = sqrt(EkinM(ip)/Mvir(ip)*2.)
         Cvir = Concentration(Mvir(ip),1.e3*Rvir(ip),VmaxM(ip))
         If(Cvir < 0.)Cvir = Rvir(ip)/RmaxM(ip)*2.15
         iHalo= ih

         iStart  = 0
         Do i=-NradP+1,0
            If(NbinH1(i,ih)>0.and.MassH1(i,ih)>5.*MassOne)Then
               iStart = i ; exit
            EndIf
         EndDo
         Nlines  = 0 ! total lines of profile
         Do i=-NradP+1,0
            If(NbinH1(i,ih)>0.and.MassH1(i,ih)>5.*MassOne)Then
               Nlines = Nlines +1
            EndIf
         EndDo
           write(20) &
                    x,y,z,VxMaxx(ip),VyMaxx(ip),VzMaxx(ip), &
                    Mvir(ip),Mtotal(ip),1.e3*Rvir(ip),Vrms, VmaxM(ip),   & 
                    iHalo,Cvir,Mvir(ip)/MassOne,MaxIndex(ip),Xoff(ip), &
                    2.*EkinM(ip)/EpotM(ip)-1.,LambdaM(ip),1.e3*RadRms(ip),                 &
                    Axba(ip),Axca(ip),Xax(ip),Yax(ip),Zax(ip),Nlines
          Radius = 2.*Rvir(ip)
          Do i=iStart,0
            R     = Radius*10.**(i*dLogP)
            Rin   = Radius*10.**((i-1)*dLogP)
            Vcirc1 = 6.582e-5*sqrt(MassH1(i,ih)/R)/sqrt(AEXPN)
            Vcirc2 = 6.582e-5*sqrt(MassH2(i,ih)/R)/sqrt(AEXPN)
            Volume =  4.1888*(R**3-Rin**3)
            DensH1 = (MassH1(i,ih)-MassH1(i-1,ih))/Volume*1.e-9  !density Msunh/kpch**3 comoving
            DensH2 = (MassH2(i,ih)-MassH2(i-1,ih))/Volume*1.e-9
            If(NbinH1(i,ih)/= 0.and.MassH1(i,ih)>5.*MassOne)Then
 
              write(20) R*1.e3,            &
                 NbinH1(i,ih),RadH1(i,ih),MassH1(i,ih),Vcirc1, &
                 DensH1,VrmsH1(i,ih),VradH1(i,ih),VrmsrH1(i,ih), &
                 NbinH2(i,ih),RadH2(i,ih),MassH2(i,ih),Vcirc2, &
                 DensH2,VrmsH2(i,ih),VradH2(i,ih),VrmsrH2(i,ih)

            End If
         EndDo
      End If
    end If
    end do

      close (20)

    end SUBROUTINE WriteProfiles

!---------------------------------------------------------------------------
!                  Find profile of each halo and subhalos
!                  
      SUBROUTINE GetProfiles
integer*8 :: ic,ip   
      Nhalo = Nmaxima
        write(13,*) ' GetProfiles. Nhalo=',Nhalo
        iHalo = 0
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (i)
        Do i=1,Nhalo
           RadH1(:,i)  =0. ; MassH1(:,i)  =0. ; VrmsH1(:,i) =0.
           VradH1(:,i) =0. ; VrmsrH1(:,i) =0. ; NbinH1(:,i) =0
           RadH2(:,i)  =0. ; MassH2(:,i)  =0. ; VrmsH2(:,i) =0.
           VradH2(:,i) =0. ; VrmsrH2(:,i) =0. ; NbinH2(:,i) =0
        EndDo
       
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (ih,ip,ic,x,y,z,aR) 
      Do ip=1,Nmaxima
             If(Mvir(ip)>100.*MassOne)Then
                x    = xMaxx(ip);  y = yMaxx(ip);  z = zMaxx(ip)
                aR   = Rvir(ip)
                if(mod(ih,1000)==0) &
                write(13,'(i8,5f9.4)') ip,x,y,z,aR
                Call HaloProfile(x,y,z,aR,ip)
             EndIf
      EndDo         ! i
    end SUBROUTINE GetProfiles

!---------------------------------------------------------------------------
!                   Get profile of a halo
      SUBROUTINE HaloProfile(x,y,z,aR,ip)
!---------------------------------------------------------------------------
      Real*4, PARAMETER ::      fiScale  =  4.333e-9
      Real*4      :: Fi(-NradP:0)
      Real*8      :: wx,wy,wz
      integer*8   :: ic,ip,jp   
      
      Radius = 2.*aR
      d0     = Radius**2           ! get final statistics of  particles
      factorZ    = 100.*sqrt(Om0/AEXPN**3+(1.-Om0)) *AEXPN 
      

      wx = VxMaxx(ip) ; wy = VyMaxx(ip) ; wz = VzMaxx(ip) 
      Call Limits(x,y,z,Radius,i1,i2,j1,j2,k1,k2)
                                             ! Get mass profile
      Do k3 =k1, k2
      Do j3 =j1, j2
      Do i3 =i1, i2
         jp =Label(i3,j3,k3)
        Do while (jp.ne.0)
           dd =(x-Xpar(jp))**2+(y-Ypar(jp))**2+(z-Zpar(jp))**2
           If(dd< d0) Then
              r = sqrt(max(dd,1.e-20))
              ii    = max(min(INT(log10(r/Radius)/dLogP),0),-NradP)
              dx   = Xpar(jp) -x
              dy   = Ypar(jp) -y
              dz   = Zpar(jp) -z
              dvx = VX(jp) - wx +factorZ*dx    ! true velocity
              dvy = VY(jp) - wy +factorZ*dy
              dvz = VZ(jp) - wz +factorZ*dz
              vv =  dvx**2 + dvy**2 + dvz**2   ! kinetic energy
              vr =  (dvx*dx+dvy*dy+dvz*dz)/r   ! radial velocity
              MassH1(ii,ih)   = MassH1(ii,ih)  + MassOne
              RadH1(ii,ih)    = RadH1(ii,ih)   + r/aR ! radius in virial units
              VrmsH1(ii,ih)   = VrmsH1(ii,ih)  + vv
              VradH1(ii,ih)   = VradH1(ii,ih)  + vr
              VrmsrH1(ii,ih)  = VrmsrH1(ii,ih) + vr**2
              NbinH1(ii,ih)   = NbinH1(ii,ih)  + 1
           EndIf                        ! dd<d0                                 
            jp =Lst(jp)
        End Do                          !  jp/= 0
      EndDo   ! i3
      EndDo   ! j3
      EndDo   ! k3

      Do ii =-NradP+1,0
             MassH1(ii,ih) = MassH1(ii,ih) + MassH1(ii-1,ih)
      EndDo                         
      Fi           = 0.            ! get potential
      iR           = min(-INT(0.301/dLogP),0)   ! potential at R =aR
      Rin          = Radius*10.**(iR*dLogP)
      Fi(iR)       = fiScale*MassH1(iR,ih)/Rin
      Do i =iR+1,0                         ! outer part of profile: only mass
         Rin   = Radius*10.**(i*dLogP)     !  inside aR fi =GM(aR)/R
         Fi(i) = fiScale*MassH1(iR,ih)/Rin
      EndDo
      Do i =iR-1,-NradP,-1                 ! integrate inner part 
         Rin   = Radius*10.**(i*dLogP)
         Rout  = Radius*10.**((i+1)*dLogP)
         Fi(i) = Fi(i+1) + fiScale*(MassH1(i,ih)+MassH1(i+1,ih))*0.5 &
                                                 *(Rout-Rin)/(Rout*Rin)
      EndDo
      Fi = Fi/AEXPN
       !write(13,'(3g12.4)') (Fi(i),MassP(i),Radius*10.**(i*dLogR),i=-10,0)

        Do k3 =k1, k2   ! ----------- get final statistics of bound particles
        Do j3 =j1, j2
        Do i3 =i1, i2
          jp =Label(i3,j3,k3)
          Do while (jp.ne.0)
             dd =(x-Xpar(jp))**2+(y-Ypar(jp))**2+(z-Zpar(jp))**2
             If(dd< d0) Then
                r = sqrt(max(dd,1.e-20))
                dx   = Xpar(jp) -x
                dy   = Ypar(jp) -y
                dz   = Zpar(jp) -z
                dvx = VX(jp) - wx +factorZ*dx    ! true velocity
                dvy = VY(jp) - wy +factorZ*dy
                dvz = VZ(jp) - wz +factorZ*dz
                ii    = max(min(INT(log10(r/Radius)/dLogP),0),-NradP)
                vv =  dvx**2 + dvy**2 + dvz**2   ! kinetic energy
                 ee = -Fi(ii) + 0.5*vv
              if(ee <= 0.)Then                          
                 vr =  (dvx*dx+dvy*dy+dvz*dz)/r   ! radial velocity
                 MassH2(ii,ih)   = MassH2(ii,ih)  + MassOne
                 RadH2(ii,ih)    = RadH2(ii,ih)   + r/aR ! radius in virial units
                 VrmsH2(ii,ih)   = VrmsH2(ii,ih)  + vv
                 VradH2(ii,ih)   = VradH2(ii,ih)  + vr
                 VrmsrH2(ii,ih)  = VrmsrH2(ii,ih) + vr**2
                 NbinH2(ii,ih)   = NbinH2(ii,ih)  + 1 
              end if                      ! ee<0
             EndIf                        ! dd<d0                                 
            jp =Lst(jp)
          End Do                          !  jp/= 0
        EndDo   ! i3
        EndDo   ! j3
        EndDo   ! k3

        !Do ii =-NradP+1,0
        ! write(13,'(2i8,g12.4)') ii,NbinH2(ii,ih),MassH2(ii,ih)
        !EndDo
        Do ii =-NradP+1,0
                 MassH2(ii,ih) = MassH2(ii,ih) + MassH2(ii-1,ih)
        EndDo
        
        Do ii = -NradP+1,0
          If(NbinH1(ii,ih) /= 0)Then
            RadH1(ii,ih) = RadH1(ii,ih)/NbinH1(ii,ih)
            VrmsH1(ii,ih) = sqrt(VrmsH1(ii,ih)/NbinH1(ii,ih))
            VradH1(ii,ih) = VradH1(ii,ih)/NbinH1(ii,ih)
            VrmsrH1(ii,ih) = sqrt(VrmsrH1(ii,ih)/NbinH1(ii,ih))
          end If
          If(NbinH2(ii,ih) /= 0)Then
            RadH2(ii,ih)   = RadH2(ii,ih)/NbinH2(ii,ih)
            VrmsH2(ii,ih)  = sqrt(VrmsH2(ii,ih)/NbinH2(ii,ih))
            VradH2(ii,ih)  = VradH2(ii,ih)/NbinH2(ii,ih)
            VrmsrH2(ii,ih) = sqrt(VrmsrH2(ii,ih)/NbinH2(ii,ih))
          end If
       End Do
       !Do ii =-NradP+1,0
       ! write(13,'(2i8,8g12.4)') ii,NbinH1(ii,ih),MassH1(ii,ih), &
       !          RadH1(ii,ih),VrmsH1(ii,ih),VrmsrH1(ii,ih)
       !EndDo

             end SUBROUTINE HaloProfile

!---------------------------------------------------------------------------
!                 Initialize arrays for halo structure
!  
      SUBROUTINE InitMaxima
!--------------------------------------------------------------------------- 

!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (ip) 
        Do ip =1, Nmaxima
           !MassProf(:,ip)= 0.
           VxMaxx(ip)      = 0.
           VyMaxx(ip)      = 0.
           VzMaxx(ip)      = 0.
           Xoff(ip)           = 0.
           LambdaM(ip) = 0.
           EpotM(ip)       = 0.
           EkinM(ip)        = 0.
           RadRms(ip)   = 0.
           Axba(ip)          = 0.
           Axca(ip)          = 0.
           Xax(ip)            = 0.
           Yax(ip)            = 0.
           Zax(ip)            = 0.
        EndDo
      end SUBROUTINE InitMaxima

!---------------------------------------------------------------------------
!                   
!
      SUBROUTINE RemoveDuplicates
!---------------------------------------------------------------------------
! A distinct host centre lies outside every higher-priority host's aperture.
! Priority is larger bound mass, then lower stable candidate index on ties.
! Read immutable measurements in parallel, then apply the mask after the barrier:
! a suppressing host can itself be suppressed, preserving legacy host chains.
! Exact particle-set duplicates are handled separately after this host test.
        implicit none
        integer*8 :: ip,jp
        integer :: i1,i2,j1,j2,k1,k2,i3,j3,k3,sx,sy,sz
        integer :: sxlo,sxhi,sylo,syhi,szlo,szhi
        real*4 :: tstart,tfinish
        real*8 :: x,y,z,xx,yy,zz,dx,dy,dz,dd,radius,box64,cell64
        logical, allocatable :: removeHost(:)
        tstart=seconds()
        if(Nmaxima==0)return
        radius=max(3.5d0,dble(maxval(Rvir)))
        box64=dble(Box);cell64=dble(Cell)
        allocate(removeHost(Nmaxima))
        removeHost=.false.
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ip,jp,x,y,z,xx,yy,zz,dx,dy,dz,dd,i1,i2,j1,j2,k1,k2,i3,j3,k3) &
!$OMP PRIVATE(sx,sy,sz,sxlo,sxhi,sylo,syhi,szlo,szhi)
        do ip=1,Nmaxima
          if(Mvir(ip)<=MassOne)cycle
          x=dble(xMaxx(ip));y=dble(yMaxx(ip));z=dble(zMaxx(ip))
          sxlo=0;sxhi=0;sylo=0;syhi=0;szlo=0;szhi=0
          if(x-radius<0.d0)sxhi=1
          if(x+radius>=box64)sxlo=-1
          if(y-radius<0.d0)syhi=1
          if(y+radius>=box64)sylo=-1
          if(z-radius<0.d0)szhi=1
          if(z+radius>=box64)szlo=-1
          do sz=szlo,szhi
          do sy=sylo,syhi
          do sx=sxlo,sxhi
            xx=x+dble(sx)*box64;yy=y+dble(sy)*box64;zz=z+dble(sz)*box64
            ! Match ListMaxima's double-precision ceiling bins. In particular,
            ! never round a shifted x+Box query through float32 before finding
            ! its bounds: a true neighbouring host could move outside that bin.
            i1=min(max(Nmx,ceiling((xx-radius)/cell64)-1),Nbx)
            i2=min(max(Nmx,ceiling((xx+radius)/cell64)-1),Nbx)
            j1=min(max(Nmy,ceiling((yy-radius)/cell64)-1),Nby)
            j2=min(max(Nmy,ceiling((yy+radius)/cell64)-1),Nby)
            k1=min(max(Nmz,ceiling((zz-radius)/cell64)-1),Nbz)
            k2=min(max(Nmz,ceiling((zz+radius)/cell64)-1),Nbz)
            do k3=k1,k2
            do j3=j1,j2
            do i3=i1,i2
              jp=Label(i3,j3,k3)
              do while(jp/=0)
                if(Mvir(ip)<Mvir(jp).or.(Mvir(ip)==Mvir(jp).and.jp<ip))then
                  dx=x-dble(xMaxx(jp));dx=dx-box64*anint(dx/box64)
                  dy=y-dble(yMaxx(jp));dy=dy-box64*anint(dy/box64)
                  dz=z-dble(zMaxx(jp));dz=dz-box64*anint(dz/box64)
                  dd=dx*dx+dy*dy+dz*dz
                  if(dd<dble(Rvir(jp))**2)removeHost(ip)=.true.
                endif
                jp=Lst(jp)
              enddo
            enddo
            enddo
            enddo
          enddo
          enddo
          enddo
        enddo
        where(removeHost)Mvir=0.
        deallocate(removeHost)
        call MergeNumericalDuplicates
        tfinish=seconds()
        write(*,'(10x,a,T50,2f10.2)') ' time for RemoveDuplicates =',tfinish-tstart,tfinish-t0
      end SUBROUTINE RemoveDuplicates

!---------------------------------------------------------------------------
! Merge identical bound-particle sets among host-test survivors.
! The lowest original candidate index is the deterministic representative.

      SUBROUTINE MergeNumericalDuplicates
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        implicit none
        integer :: ip,jp,n,i,left,middle,right,a,b,k,width,representative
        integer*8 :: member_index
        integer, allocatable :: order(:),scratch(:)
        ! Exact survivor identity is the numerical duplicate definition. Sorting
        ! avoids an all-pairs scan and makes the lowest candidate index win.
        if(.not.all(ieee_is_finite(Mvir)))error stop 'BDM duplicate input mass is nonfinite'
        n=count(Mvir>MassOne)
        if(n==0)return
        if(.not.allocated(BoundParticleIds)) &
          error stop 'BDM duplicate removal requires bound particle identities'
        if(size(BoundParticleIds)/=Nmaxima)error stop 'BDM membership size mismatch'
        allocate(order(n),scratch(n))
        i=0
        do ip=1,Nmaxima
          if(Mvir(ip)<=MassOne)cycle
          if(.not.allocated(BoundParticleIds(ip)%ids)) &
            error stop 'BDM candidate has no bound particle identities'
          if(size(BoundParticleIds(ip)%ids)==0) &
            error stop 'BDM positive bound mass has an empty particle set'
          do member_index=1,size(BoundParticleIds(ip)%ids,kind=8)
            if(BoundParticleIds(ip)%ids(member_index)<=0_8)error stop 'BDM original particle ID must be positive'
            if(member_index==1)cycle
            if(BoundParticleIds(ip)%ids(member_index)<=BoundParticleIds(ip)%ids(member_index-1)) &
              error stop 'BDM bound identities must be sorted and unique'
          enddo
          i=i+1;order(i)=ip
        enddo
        width=1
        do while(width<n)
          left=1
          do while(left<=n)
            middle=left+min(width,n-left+1)-1
            right=middle+min(width,n-middle)
            a=left;b=middle+1
            do k=left,right
              if(a>middle)then
                scratch(k)=order(b);b=b+1
              elseif(b>right)then
                scratch(k)=order(a);a=a+1
              elseif(compare_sets(order(a),order(b))<=0)then
                scratch(k)=order(a);a=a+1
              else
                scratch(k)=order(b);b=b+1
              endif
            enddo
            left=right+1
          enddo
          order=scratch
          if(width>n/2)exit
          width=2*width
        enddo
        representative=order(1)
        do i=2,n
          jp=order(i)
          if(same_set(representative,jp))then
            Mvir(jp)=0.
          else
            representative=jp
          endif
        enddo
        deallocate(order,scratch)
      contains
        integer function compare_sets(first,second) result(comparison)
          integer,intent(in)::first,second
          integer*8 :: count_first,count_second,j
          comparison=0
          count_first=size(BoundParticleIds(first)%ids,kind=8)
          count_second=size(BoundParticleIds(second)%ids,kind=8)
          if(count_first<count_second)then
            comparison=-1;return
          elseif(count_first>count_second)then
            comparison=1;return
          endif
          do j=1,count_first
            if(BoundParticleIds(first)%ids(j)<BoundParticleIds(second)%ids(j))then
              comparison=-1;return
            elseif(BoundParticleIds(first)%ids(j)>BoundParticleIds(second)%ids(j))then
              comparison=1;return
            endif
          enddo
          ! Resolve identical sets by stable original candidate index.
          if(first<second)comparison=-1
          if(first>second)comparison=1
        end function compare_sets
        logical function same_set(first,second)
          integer,intent(in)::first,second
          same_set=.false.
          if(size(BoundParticleIds(first)%ids,kind=8)/= &
             size(BoundParticleIds(second)%ids,kind=8))return
          same_set=all(BoundParticleIds(first)%ids==BoundParticleIds(second)%ids)
        end function same_set
      end SUBROUTINE MergeNumericalDuplicates
!---------------------------------------------------------------------------
!                   
!
      SUBROUTINE RemoveDuplicatesSimple
!---------------------------------------------------------------------------
        integer*8 :: ic,ip,i
        real*4 :: m
      tstart = seconds()
!  --------------------------- 
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (ip,x,y,z,m,ic,D2) 
      Do ip=1,Nmaxima
         If(Mvir(ip)>MassOne)Then
            x   = xMaxx(ip);   y = yMaxx(ip);    z = zMaxx(ip)
            m   = Mvir(ip)
            do ic =1,Nmaxima
               If(ic/=ip)Then
                  D2 = (xMaxx(ic) -x)**2 +(yMaxx(ic) -y)**2 +(zMaxx(ic) -z)**2
                  If(D2.lt.Rvir(ic)**2.and.m.lt.Mvir(ic))Then
                     Mvir(ip) = 0.
                  End If
               end If
            end do
         end If
       End Do         ! ip
       tfinish = seconds()
      write(*,'(10x,a,T50,2f10.2)') ' time for RemoveDuplicates =',tfinish-tstart,tfinish-t0

    end SUBROUTINE RemoveDuplicatesSimple
!---------------------------------------------------------------------------
!                  Find parameters of distinct halos 
!
      SUBROUTINE ParametersDistinct
!---------------------------------------------------------------------------
integer*8 :: ic,ip,i    
      tstart = seconds()
      call BdmHaloMembershipInit

!  --------------------------- 
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (ip,x,y,z,xv,yv,zv) ! SCHEDULE(DYNAMIC,10000)
      Do ip=1,Nmaxima
            x   = xMaxx(ip);   y = yMaxx(ip);    z = zMaxx(ip)
            xv  = VxMaxx(ip); yv = VyMaxx(ip);  zv = VzMaxx(ip)
            Call GetHalo(x,y,z,xv,yv,zv,ip)
       End Do         ! ip
       tfinish = seconds()
      write(*,'(10x,a,T50,2f10.2)') ' time for ParametersDistinct =',tfinish-tstart,tfinish-t0

    end SUBROUTINE ParametersDistinct
!---------------------------------------------------------------------------
!                   Get parameters of distinct halos
      SUBROUTINE GetHalo(x,y,z,xv,yv,zv,ip)
!---------------------------------------------------------------------------
! SO uses the outermost crossing of the discrete enclosed-particle profile.
! The legacy Rext correction still defines the reported aperture and Mtotal.
! Mvir, drift, kinetic energy, shape, spin, and Vmax use only the converged bound
! population inside the unextended SO sphere. Binding uses the isolated,
! spherical Newtonian potential, with each particle's self term excluded.
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        implicit none
        real*4, intent(in) :: x,y,z,xv,yv,zv
        integer*8, intent(in) :: ip
        real*8, parameter :: gravity=4.333d-9
        integer*8, allocatable :: rows(:)
        real*8, allocatable :: radii(:),potential(:)
        real*8 :: search_cap,aperture_cap,rso,aperture,grid_size
        real*8 :: mass,threshold,hubble_a,bulk(3),offset(3),velocity(3)
        real*8 :: kinetic,rms2,centre(3),angular(3),tensor(3,3),energy
        real*8 :: shell_energy,circular2,maximum2,maximum_radius,correction
        real*8 :: axis_ratio(3),concentration_proxy,slope_b,slope_c
        real*4 :: axis(3),direction(3)
        integer :: n,nkeep,naperture,q,j,iteration
        logical :: singular

        ! ParametersDistinct allocates before entering its parallel loop. The
        ! guarded path also supports direct, sequential calls used by tests.
        if(.not.allocated(BoundParticleIds)) call BdmHaloMembershipInit
        if(size(BoundParticleIds)/=Nmaxima) error stop 'BDM membership size mismatch'
        if(allocated(BoundParticleIds(ip)%ids)) deallocate(BoundParticleIds(ip)%ids)
        HaloStatus(ip)=0
        Mvir(ip)=0.; Mtotal(ip)=0.; Rvir(ip)=0.
        VmaxM(ip)=0.; RmaxM(ip)=0.; EkinM(ip)=0.; EpotM(ip)=0.
        VxMaxx(ip)=0.; VyMaxx(ip)=0.; VzMaxx(ip)=0.
        Xoff(ip)=0.; LambdaM(ip)=0.; RadRms(ip)=0.
        Axba(ip)=0.; Axca(ip)=0.; Xax(ip)=0.; Yax(ip)=0.; Zax(ip)=0.
        if(.not.all(ieee_is_finite([x,y,z,Cell,Box,MassOne,Om0,Ovdens,AEXPN]))) &
             error stop 'Non-finite input to BDM GetHalo'
        if(min(Cell,Box,MassOne,Om0,Ovdens,AEXPN)<=0..or.NGRID<=0) &
             error stop 'Non-positive scale in BDM GetHalo'
        mass=dble(MassOne)
        threshold=1.150d12*dble(Om0)*dble(Ovdens)
        grid_size=dble(Box)/NGRID
        hubble_a=100.d0*sqrt(dble(Om0)/dble(AEXPN)**3+1.d0-dble(Om0))*dble(AEXPN)
        search_cap=min(15.d0*dble(Cell),dble(nearest(.5*Box,-1.)))
        aperture_cap=min(search_cap+.75d0*grid_size,dble(nearest(.5*Box,-1.)))
        if(HaloSearchRadius>0.) search_cap=min(search_cap,dble(HaloSearchRadius))
        if(ParticleSearchRadius>0.) aperture_cap=min(aperture_cap,dble(ParticleSearchRadius))
        ! The physical cap, not the first sampled underdense Cell, defines
        ! the SO domain. Nonmonotonic profiles can cross the threshold again
        ! outside an earlier underdense shell.
        call BdmHaloGather(x,y,z,search_cap,rows,radii)
        n=size(rows)
        if(dble(n)*mass>threshold*search_cap**3)then
          HaloStatus(ip)=HaloSearchTruncated
          return
        endif
        if(n<10)then
          HaloStatus(ip)=HaloTooFewParticles
          return
        endif
        ! With n available particles, no SO root can exceed (n*m/rho)^(1/3).
        ! Discarding rows outside that bound preserves every possible root.
        ! Repeating gives the greatest self-consistent enclosed population;
        ! most diffuse cap-neighbourhood particles need never be sorted.
        rso=min((dble(n)*mass/threshold)**(1.d0/3.d0),search_cap)
        do iteration=1,16
          nkeep=0
          do q=1,n
            if(radii(q)<=rso)then
              nkeep=nkeep+1
              rows(nkeep)=rows(q); radii(nkeep)=radii(q)
            endif
          enddo
          if(nkeep==n) exit
          n=nkeep
          if(n<10)then
            HaloStatus(ip)=HaloTooFewParticles
            return
          endif
          rso=min((dble(n)*mass/threshold)**(1.d0/3.d0),search_cap)
        enddo
        if(iteration>16)then
          ! Bound contraction can remove one row per pass in an adversarial
          ! profile. Bound that cost and finish with an O(n log n) exact scan.
          call BdmHaloSortRadii(rows(:n),radii(:n))
          rso=0.d0
          do q=n,10,-1
            energy=(dble(q)*mass/threshold)**(1.d0/3.d0)
            if(energy<radii(q)) cycle
            if(q<n)then
              if(energy>=radii(q+1)) cycle
            endif
            if(energy>search_cap) cycle
            rso=energy
            exit
          enddo
          if(rso<=0.d0)then
            HaloStatus(ip)=HaloTooFewParticles
            return
          endif
        endif
        aperture=rso+grid_size*min(dble(Rext)/(rso/grid_size)**dble(SlopeR),.75d0)
        if(aperture>aperture_cap)then
          ! A larger aperture would require an image outside the verified
          ! search domain. Reject explicitly rather than count a clipped sphere.
          HaloStatus(ip)=HaloSearchTruncated
          return
        endif
        call BdmHaloGather(x,y,z,aperture,rows,radii)
        call BdmHaloSortRadii(rows,radii)
        naperture=size(rows)
        Mtotal(ip)=real(dble(naperture)*mass)
        Rvir(ip)=real(aperture)
        n=count(radii<=rso)
        if(n==0)then
          HaloStatus(ip)=HaloNoBoundParticles
          return
        endif
        allocate(potential(n))

        ! Membership only shrinks. Every non-final pass removes at least one
        ! particle; therefore convergence needs at most the initial n+1 passes.
        ! No removed particle contributes to the next potential or drift.
        do iteration=1,n+1
          if(n==0) exit
          bulk=0.d0
          do q=1,n
            bulk=bulk+[dble(VX(rows(q))),dble(VY(rows(q))),dble(VZ(rows(q)))]
          enddo
          bulk=bulk/dble(n)
          call BdmHaloSphericalPotential(radii(:n),mass,gravity/dble(AEXPN), &
                                         potential(:n),shell_energy,singular)
          if(singular)then
            ! Multiple particles exactly at the centre have an undefined
            ! unsoftened shell potential. Do not invent a softening length.
            HaloStatus(ip)=ior(HaloStatus(ip),HaloSingularCentre)
            Mvir(ip)=0.
            return
          endif
          nkeep=0
          do q=1,n
            call BdmParticlePosition(rows(q),offset)
            offset=offset-[dble(x),dble(y),dble(z)]
            velocity=[dble(VX(rows(q))),dble(VY(rows(q))),dble(VZ(rows(q)))] &
                      -bulk+hubble_a*offset
            energy=.5d0*sum(velocity**2)-potential(q)
            if(energy<=0.d0)then
              nkeep=nkeep+1
              rows(nkeep)=rows(q); radii(nkeep)=radii(q)
            endif
          enddo
          if(nkeep==n) exit
          n=nkeep
        enddo
        if(n==0)then
          HaloStatus(ip)=ior(HaloStatus(ip),HaloNoBoundParticles)
          return
        endif

        ! The successful final pass evaluated every survivor against this same
        ! survivor set and bulk velocity. All published bound statistics below
        ! use precisely these rows and double-precision local offsets.
        kinetic=0.d0; rms2=0.d0; centre=0.d0; angular=0.d0; tensor=0.d0
        maximum2=0.d0; maximum_radius=0.d0
        do q=1,n
          call BdmParticlePosition(rows(q),offset)
          offset=offset-[dble(x),dble(y),dble(z)]
          velocity=[dble(VX(rows(q))),dble(VY(rows(q))),dble(VZ(rows(q)))] &
                    -bulk+hubble_a*offset
          kinetic=kinetic+sum(velocity**2)
          rms2=rms2+sum(offset**2)
          centre=centre+offset
          angular=angular+[offset(2)*velocity(3)-offset(3)*velocity(2), &
                           offset(3)*velocity(1)-offset(1)*velocity(3), &
                           offset(1)*velocity(2)-offset(2)*velocity(1)]
          if(radii(q)>0.d0)then
            do j=1,3
              tensor(:,j)=tensor(:,j)+offset*offset(j)/radii(q)**2
            enddo
            circular2=dble(q)*mass/radii(q)
            if(circular2>maximum2)then
              maximum2=circular2; maximum_radius=radii(q)
            endif
          endif
        enddo
        Mvir(ip)=real(dble(n)*mass)
        EkinM(ip)=real(.5d0*kinetic*mass)
        EpotM(ip)=real(shell_energy)
        VxMaxx(ip)=real(bulk(1)); VyMaxx(ip)=real(bulk(2)); VzMaxx(ip)=real(bulk(3))
        RadRms(ip)=real(sqrt(rms2/dble(n)))
        centre=centre/dble(n); angular=angular/dble(n)
        Xoff(ip)=real(sqrt(sum(centre**2))/aperture)
        ! Retain the legacy spin proxy and empirical axis corrections, but use
        ! one population for its angular momentum, RMS speed, and mass.
        LambdaM(ip)=real(sqrt(sum(angular**2))*dble(AEXPN)*sqrt(kinetic/dble(n)) &
                        /(dble(n)*mass)*1.632d8)
        tensor=tensor/dble(n)
        if(maxval(abs(tensor))>0.d0)then
          call EigenValues(tensor,direction,axis)
          if(all(ieee_is_finite(axis)).and.all(ieee_is_finite(direction)))then
            if(axis(1)>0.)then
              axis_ratio=sqrt(max(0.d0,min(1.d0,dble(axis)/dble(axis(1)))))
              concentration_proxy=dble(RadRms(ip))/aperture
              slope_b=1.d0+2.d0*max(concentration_proxy-.4d0,0.d0) &
                      +(5.7d0*max(concentration_proxy-.4d0,0.d0))**3
              slope_c=1.d0+2.d0*max(concentration_proxy-.4d0,0.d0) &
                      +(5.5d0*max(concentration_proxy-.4d0,0.d0))**3
              ! Different empirical powers can reverse the two transverse
              ! lengths. Report intermediate/minor order after correction;
              ! the principal-axis direction and both corrected lengths stay.
              axis_ratio(2)=axis_ratio(2)**slope_b
              axis_ratio(3)=axis_ratio(3)**slope_c
              Axba(ip)=real(max(axis_ratio(2),axis_ratio(3)))
              Axca(ip)=real(min(axis_ratio(2),axis_ratio(3)))
              Xax(ip)=direction(1); Yax(ip)=direction(2); Zax(ip)=direction(3)
            endif
          endif
        endif
        correction=.1d0*grid_size
        if(maximum_radius>correction.and.maximum2>0.d0)then
          VmaxM(ip)=real(sqrt(gravity*maximum2/dble(AEXPN)) &
                         /sqrt(1.d0-correction/maximum_radius))
          RmaxM(ip)=real(maximum_radius)
        else
          HaloStatus(ip)=ior(HaloStatus(ip),HaloUnresolvedVmax)
          ! Zero is the documented finite unresolved sentinel; Concentration
          ! and catalogue writing must preserve it without a 0/0 fallback.
          VmaxM(ip)=0.; RmaxM(ip)=0.
        endif
        allocate(BoundParticleIds(ip)%ids(n))
        if(allocated(OriginalParticleId))then
          BoundParticleIds(ip)%ids=OriginalParticleId(rows(:n))
        else
          ! Direct fixtures without periodic buffering use original row IDs.
          BoundParticleIds(ip)%ids=rows(:n)
        endif
        call BdmHaloSortIds(BoundParticleIds(ip)%ids)
        if(n>1)then
          if(any(BoundParticleIds(ip)%ids(2:)==BoundParticleIds(ip)%ids(:n-1))) &
               error stop 'Repeated original particle in BDM bound sphere'
        endif
      end SUBROUTINE GetHalo

! Recover the exact double-precision location of this explicit periodic image.
! Float32 storage of x+Box loses low bits on the positive face; the original
! row remains exact. Infer only the integer image shift from the stored ghost.
! Keeping that shift (rather than minimum-imaging each row independently) also
! ensures two explicit images cannot both enter one sub-half-box sphere.
      pure real*8 function BdmParticleCoordinate(row,component) result(coordinate)
        implicit none
        integer*8, intent(in) :: row
        integer, intent(in) :: component
        integer*8 :: original
        real*8 :: stored,base
        original=row
        if(allocated(OriginalParticleId)) original=OriginalParticleId(row)
        select case(component)
        case(1)
          stored=dble(Xpar(row)); base=dble(Xpar(original))
        case(2)
          stored=dble(Ypar(row)); base=dble(Ypar(original))
        case(3)
          stored=dble(Zpar(row)); base=dble(Zpar(original))
        case default
          coordinate=0.d0
          return
        end select
        coordinate=stored
        ! Original rows already store the exact analysis coordinate; only a
        ! periodic image needs an integer shift reconstructed in float64.
        if(original==row)return
        if(allocated(OriginalParticleId)) &
          coordinate=base+dble(nint((stored-base)/dble(Box)))*dble(Box)
      end function BdmParticleCoordinate

      SUBROUTINE BdmParticlePosition(row,position)
        implicit none
        integer*8, intent(in) :: row
        real*8, intent(out) :: position(3)
        integer :: component
        integer*8 :: original
        original=row
        if(allocated(OriginalParticleId))original=OriginalParticleId(row)
        if(original==row)then
          position=[dble(Xpar(row)),dble(Ypar(row)),dble(Zpar(row))]
          return
        endif
        do component=1,3
          position(component)=BdmParticleCoordinate(row,component)
        enddo
      end SUBROUTINE BdmParticlePosition

! Candidate-local workspace: two exact-count list scans, no Np-sized temporary.
      SUBROUTINE BdmHaloGather(x,y,z,radius,rows,radii)
        implicit none
        real*4, intent(in) :: x,y,z
        real*8, intent(in) :: radius
        integer*8, allocatable, intent(out) :: rows(:)
        real*8, allocatable, intent(out) :: radii(:)
        integer :: i1,i2,j1,j2,k1,k2,i,j,k,n,pass
        integer*8 :: jp
        real*8 :: distance2,offset(3)
        call Limits(x,y,z,nearest(real(radius),1.),i1,i2,j1,j2,k1,k2)
        do pass=1,2
          n=0
          do k=k1,k2
          do j=j1,j2
          do i=i1,i2
            jp=Label(i,j,k)
            do while(jp/=0)
              call BdmParticlePosition(jp,offset)
              offset=offset-[dble(x),dble(y),dble(z)]
              distance2=sum(offset**2)
              if(distance2<=radius**2)then
                n=n+1
                if(pass==2)then
                  rows(n)=jp; radii(n)=sqrt(distance2)
                endif
              endif
              jp=Lst(jp)
            enddo
          enddo
          enddo
          enddo
          if(pass==1) allocate(rows(n),radii(n))
        enddo
      end SUBROUTINE BdmHaloGather

! In-place heapsort: deterministic radius order and row-ID tie breaking.
      SUBROUTINE BdmHaloSortRadii(rows,radii)
        implicit none
        integer*8, intent(inout) :: rows(:)
        real*8, intent(inout) :: radii(:)
        integer*8 :: first,last,parent,child,n
        integer*8 :: saved_row
        real*8 :: saved_radius
        n=size(rows,kind=8)
        if(n<2) return
        first=n/2+1; last=n
        do
          if(first>1)then
            first=first-1
            saved_row=rows(first); saved_radius=radii(first)
          else
            saved_row=rows(last); saved_radius=radii(last)
            rows(last)=rows(1); radii(last)=radii(1)
            last=last-1
            if(last==1)then
              rows(1)=saved_row; radii(1)=saved_radius
              exit
            endif
          endif
          parent=first; child=2*first
          do while(child<=last)
            if(child<last)then
              if(radii(child)<radii(child+1))then
                child=child+1
              else if(radii(child)==radii(child+1))then
                if(rows(child)<rows(child+1)) child=child+1
              endif
            endif
            if(saved_radius>radii(child)) exit
            if(saved_radius==radii(child))then
              if(saved_row>=rows(child)) exit
            endif
            rows(parent)=rows(child); radii(parent)=radii(child)
            parent=child; child=2*child
          enddo
          rows(parent)=saved_row; radii(parent)=saved_radius
        enddo
      end SUBROUTINE BdmHaloSortRadii

      SUBROUTINE BdmHaloSortIds(ids)
        implicit none
        integer*8, intent(inout) :: ids(:)
        integer*8 :: saved
        integer*8 :: first,last,parent,child,n
        n=size(ids,kind=8)
        if(n<2) return
        first=n/2+1; last=n
        do
          if(first>1)then
            first=first-1; saved=ids(first)
          else
            saved=ids(last); ids(last)=ids(1); last=last-1
            if(last==1)then
              ids(1)=saved
              exit
            endif
          endif
          parent=first; child=2*first
          do while(child<=last)
            if(child<last)then
              if(ids(child)<ids(child+1)) child=child+1
            endif
            if(saved>=ids(child)) exit
            ids(parent)=ids(child); parent=child; child=2*child
          enddo
          ids(parent)=saved
        enddo
      end SUBROUTINE BdmHaloSortIds

! Exact discrete spherical-shell potential and distinct-pair energy. A pair at
! radii ri,rj contributes G*m^2/max(ri,rj); no shell includes its own particle.
! This is a monopole approximation, not the exact aspherical N-body potential.
      SUBROUTINE BdmHaloSphericalPotential(radii,mass,g_over_a,potential,energy,singular)
        implicit none
        real*8, intent(in) :: radii(:),mass,g_over_a
        real*8, intent(out) :: potential(:),energy
        logical, intent(out) :: singular
        real*8 :: exterior,interior
        integer :: q,n
        n=size(radii); exterior=0.d0; energy=0.d0; potential=0.d0
        singular=.false.
        if(n>1)then
          if(radii(2)==0.d0)then
            singular=.true.
            return
          endif
        endif
        do q=n,1,-1
          interior=0.d0
          if(q>1) interior=dble(q-1)/radii(q)
          potential(q)=g_over_a*mass*(interior+exterior)
          energy=energy+g_over_a*mass**2*interior
          if(radii(q)>0.d0) exterior=exterior+1.d0/radii(q)
        enddo
      end SUBROUTINE BdmHaloSphericalPotential

      SUBROUTINE BdmHaloMembershipInit
        implicit none
        ! Called before the production parallel loop, or from a direct fixture.
!$OMP CRITICAL (bdm_halo_membership_init)
        if(allocated(BoundParticleIds)) deallocate(BoundParticleIds)
        if(allocated(HaloStatus)) deallocate(HaloStatus)
        allocate(BoundParticleIds(Nmaxima),HaloStatus(Nmaxima))
        HaloStatus=0
!$OMP END CRITICAL (bdm_halo_membership_init)
      end SUBROUTINE BdmHaloMembershipInit


!---------------------------------------------------------------------------
!                   get centers and velocities of maxima
!
!

SUBROUTINE FindDistinctCandidates
  implicit none
  integer*8 :: jp,nn
  integer :: im,iter,i1,i2,j1,j2,k1,k2,i3,j3,k3
  real :: x,y,z,Radius,timeStart,timeFinish
  real*8 :: xc,yc,zc,xv,yv,zv,dx,dy,dz,d0,position(3)

  timeStart=seconds()
  Mvir=0.; Rvir=0.; VxMaxx=0.; VyMaxx=0.; VzMaxx=0.
  do iter=1,4
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(im,x,y,z,xc,yc,zc,xv,yv,zv,nn,i3,j3,k3,i1,i2,j1,j2,k1,k2,jp,dx,dy,dz,Radius,d0,position)
    do im=1,Nmaxima
      x=xMaxx(im); y=yMaxx(im); z=zMaxx(im)
      xc=0.d0; yc=0.d0; zc=0.d0
      xv=0.d0; yv=0.d0; zv=0.d0; nn=0_8
      Radius=min(Cell*max(0.5,min(2.,log10(max(Xoff(im),0.)+10.)/2.)),nearest(0.5*Box,-1.))
      d0=dble(Radius)**2
      call Limits(x,y,z,Radius,i1,i2,j1,j2,k1,k2)
      do k3=k1,k2
      do j3=j1,j2
      do i3=i1,i2
        jp=Label(i3,j3,k3)
        do while(jp/=0_8)
          call BdmParticlePosition(jp,position)
          dx=position(1)-dble(x)
          dy=position(2)-dble(y)
          dz=position(3)-dble(z)
          if(dx*dx+dy*dy+dz*dz<d0)then
            ! Sum local displacements: the centre must remain in its cloud,
            ! independently of the absolute box location or particle count.
            xc=xc+dx; yc=yc+dy; zc=zc+dz
            xv=xv+dble(VX(jp)); yv=yv+dble(VY(jp)); zv=zv+dble(VZ(jp))
            nn=nn+1_8
          endif
          jp=Lst(jp)
        enddo
      enddo
      enddo
      enddo
      if(nn>0_8)then
        xMaxx(im)=real(modulo(dble(x)+xc/dble(nn),dble(Box)))
        yMaxx(im)=real(modulo(dble(y)+yc/dble(nn),dble(Box)))
        zMaxx(im)=real(modulo(dble(z)+zc/dble(nn),dble(Box)))
        if(xMaxx(im)>=Box)xMaxx(im)=0.
        if(yMaxx(im)>=Box)yMaxx(im)=0.
        if(zMaxx(im)>=Box)zMaxx(im)=0.
        VxMaxx(im)=real(xv/dble(nn)); VyMaxx(im)=real(yv/dble(nn)); VzMaxx(im)=real(zv/dble(nn))
        Mvir(im)=real(dble(nn)*dble(MassOne)); Rvir(im)=Radius
      else
        ! An empty recentering aperture is not a valid mass estimate. Retain
        ! the last position, but clear derived values instead of stale mass.
        Mvir(im)=0.; Rvir(im)=0.
        VxMaxx(im)=0.; VyMaxx(im)=0.; VzMaxx(im)=0.
      endif
    enddo
  enddo
  timeFinish=seconds()
  write(13,'(10x,a,T50,2f10.2)') ' time for FindDistinctCandidates (secs) =',timeFinish-timeStart,timeFinish-t0
end SUBROUTINE FindDistinctCandidates

!---------------------------------------------------------------------------
!                  
!                  
!                  
!                 
      SUBROUTINE FindMaxima
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      use omp_lib, only: omp_get_wtime
      implicit none
      integer*8, allocatable :: PlaneOffsets(:)
      integer*8 :: plane_count
      integer :: m1,m2,m3,slot
      real :: xs,maximum_density,memory_usage
      real*8 :: started,counted,finished
      logical :: bad_density

      started=omp_get_wtime()
      call SetOverdensity
      if (NGRID < 1) error stop 'BDM peak mesh size must be positive'
      if (.not.allocated(FI)) error stop 'BDM peak density is not allocated'
      if (any(shape(FI) /= NGRID)) error stop 'BDM peak density shape does not match NGRID'
      write(*,'(a,i0,a,f12.4)') ' FindMaxima: particles=',Nparticles,', overdensity=',Ovdens

      ! Two scans avoid a full-mesh flag array. Counts and prefix offsets are
      ! per plane, independent of OpenMP scheduling and of the peak density.
      allocate(PlaneOffsets(0:NGRID)); PlaneOffsets=0_8
      maximum_density=0.; bad_density=.false.
!$OMP PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC) &
!$OMP PRIVATE(m1,m2,m3,plane_count) REDUCTION(MAX:maximum_density) REDUCTION(.OR.:bad_density)
      do m3=1,NGRID
         plane_count=0_8
         do m2=1,NGRID
            do m1=1,NGRID
               if (.not.ieee_is_finite(FI(m1,m2,m3))) then
                  bad_density=.true.
                  cycle
               end if
               maximum_density=max(maximum_density,FI(m1,m2,m3))
               if (FI(m1,m2,m3) <= Ovdens/3.) cycle
               if (IsDensityMaximum(m1,m2,m3)) plane_count=plane_count+1_8
            end do
         end do
         PlaneOffsets(m3)=plane_count
      end do
      if (bad_density) error stop 'BDM density contains a nonfinite value'
      do m3=1,NGRID
         PlaneOffsets(m3)=PlaneOffsets(m3)+PlaneOffsets(m3-1)
      end do
      if (PlaneOffsets(NGRID) > int(huge(Nmaxima),8)) &
         error stop 'BDM peak count exceeds supported integer range'
      Nmaxima=int(PlaneOffsets(NGRID))
      counted=omp_get_wtime()
      write(*,'(a,g14.6,a,i0)') ' Maximum density=',maximum_density,', number of maxima=',Nmaxima

      ! Zero-length allocations are intentional: BDM can publish an empty
      ! catalogue and release these arrays without ending an inline simulation.
      allocate(Mvir(Nmaxima),Rvir(Nmaxima),Xoff(Nmaxima))
      allocate(xMaxx(Nmaxima),yMaxx(Nmaxima),zMaxx(Nmaxima))
      allocate(VxMaxx(Nmaxima),VyMaxx(Nmaxima),VzMaxx(Nmaxima))
      allocate(LstMax(Nmaxima),EpotM(Nmaxima),EkinM(Nmaxima),LambdaM(Nmaxima))
      allocate(VmaxM(Nmaxima),RmaxM(Nmaxima),Mtotal(Nmaxima),RadRms(Nmaxima))
      allocate(Xax(Nmaxima),Yax(Nmaxima),Zax(Nmaxima),Axba(Nmaxima),Axca(Nmaxima))
      memory_usage=Memory(22_8*Nmaxima)

      xs=Box/NGRID
!$OMP PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC) PRIVATE(m1,m2,m3,slot)
      do m3=1,NGRID
         slot=int(PlaneOffsets(m3-1))
         do m2=1,NGRID
            do m1=1,NGRID
               if (FI(m1,m2,m3) <= Ovdens/3.) cycle
               if (.not.IsDensityMaximum(m1,m2,m3)) cycle
               slot=slot+1
               Xoff(slot)=FI(m1,m2,m3)
               xMaxx(slot)=(m1-1)*xs
               yMaxx(slot)=(m2-1)*xs
               zMaxx(slot)=(m3-1)*xs
            end do
         end do
         if (int(slot,8) /= PlaneOffsets(m3)) error stop 'BDM peak count changed between scans'
      end do
      deallocate(PlaneOffsets)
      finished=omp_get_wtime()
      write(*,'(a,i0,a,i0)') ' Maxima above overdensity=',count(Xoff >= Ovdens), &
                           ', maxima above 300=',count(Xoff > 300.)
      write(13,'(10x,a,2f10.3)') 'time for FindMaxima: count/list =',counted-started,finished-counted
      write(*,'(10x,a,2f10.3)') 'time for FindMaxima: count/list =',counted-started,finished-counted
      end SUBROUTINE FindMaxima

      pure logical function IsDensityMaximum(m1,m2,m3) result(selected)
      implicit none
      integer, intent(in) :: m1,m2,m3
      integer :: a,b,c,i,j,k,ix(3),iy(3),iz(3)
      integer*8 :: here,neighbour
      real :: density
      ! Immutable density and a total (z,y,x) index order break equal-density
      ! neighbour ties deterministically, including across periodic boundaries.
      ! This is a local 26-neighbour rule, not a watershed/connected-plateau fit.
      selected=.false.
      density=FI(m1,m2,m3)
      here=int(m1-1,8)+int(NGRID,8)*(int(m2-1,8)+int(NGRID,8)*int(m3-1,8))
      ix=[modulo(m1-2,NGRID)+1,m1,modulo(m1,NGRID)+1]
      iy=[modulo(m2-2,NGRID)+1,m2,modulo(m2,NGRID)+1]
      iz=[modulo(m3-2,NGRID)+1,m3,modulo(m3,NGRID)+1]
      do c=1,3
         k=iz(c)
         do b=1,3
            j=iy(b)
            do a=1,3
               i=ix(a)
               neighbour=int(i-1,8)+int(NGRID,8)*(int(j-1,8)+int(NGRID,8)*int(k-1,8))
               if (neighbour == here) cycle
               if (FI(i,j,k) > density) return
               if (FI(i,j,k) == density.and.neighbour < here) return
            end do
         end do
      end do
      selected=.true.
      end function IsDensityMaximum
!---------------------------------------------------------------------------
!                     Define size and boundaries of linked-list
      SUBROUTINE SizeListMaxima
!---------------------------------------------------------------------------
      Real*4 ::  MemList,MemMaxList,  &
                       fractionMemory  = 0.90  ! fraction of memory allocated to Lists

      Cell                  = 3.5  !Roptimal
      k    = 0
      write(13,'(2(a,g12.4))') '    SizeListMaxima: ',Cell
         Nmx         = (Xleft  -dBuffer)/Cell - 1
         Nmy         = (Yleft  -dBuffer)/Cell - 1
         Nmz         = (Zleft  -dBuffer)/Cell - 1
         Nbx         = (Xright +dBuffer)/Cell + 1
         Nby         = (Yright +dBuffer)/Cell + 1
         Nbz         = (Zright +dBuffer)/Cell + 1
      write(13,'(a,g12.3,a,6i6)') '     Size List:  Cell=',Cell, &
                          ' Limits=',Nmx,Nbx,Nmy,Nby,Nmz,Nbz

    end SUBROUTINE SizeListMaxima
!---------------------------------------------------------------------------
!                     Define size and boundaries of linked-list

SUBROUTINE PrepareParticleSearch
  implicit none
  real :: halfBox
  if(Box<=0..or.NGRID<=0)error stop 'Invalid BDM periodic box or mesh'
  ! Cell controls physical searches as well as the list geometry. It must not
  ! change after periodic images have been prepared to cover those searches.
  Cell=2.*(Box/NGRID)
  halfBox=nearest(0.5*Box,-1.)
  HaloSearchRadius=min(15.*Cell,halfBox)
  ParticleSearchRadius=min(HaloSearchRadius+0.75*(Box/NGRID),halfBox)
  ! No supported spherical query needs a wider ghost layer than half a box.
  ! This also makes the legacy five-Mpc default valid in very small boxes.
  dBuffer=min(max(dBuffer,ParticleSearchRadius),halfBox)
  Xleft=0.; Xright=Box; Yleft=0.; Yright=Box; Zleft=0.; Zright=Box
end SUBROUTINE PrepareParticleSearch

SUBROUTINE SizeList
  implicit none
  integer*8 :: cellCount,listBytes
  real*8 :: requiredGiB
  call PrepareParticleSearch
  Nmx=floor((Xleft-dBuffer)/Cell)-1; Nbx=ceiling((Xright+dBuffer)/Cell)+1
  Nmy=floor((Yleft-dBuffer)/Cell)-1; Nby=ceiling((Yright+dBuffer)/Cell)+1
  Nmz=floor((Zleft-dBuffer)/Cell)-1; Nbz=ceiling((Zright+dBuffer)/Cell)+1
  cellCount=(Nbx-Nmx+1_8)*(Nby-Nmy+1_8)*(Nbz-Nmz+1_8)
  listBytes=8_8*(Np+cellCount)
  requiredGiB=dble(listBytes)/1024.d0**3
  ! Exact int64 storage estimate. Do not coarsen physical Cell to satisfy an
  ! unrelated allocation estimate; reject insufficient memory explicitly.
  if(requiredGiB+dble(Memory(0_8))>dble(MaxMemory)) &
    error stop 'BDM linked-list allocation exceeds configured memory limit'
  write(13,'(a,g12.3,a,6i6,a,f10.4)') ' Size List: Cell=',Cell, &
       ' Limits=',Nmx,Nbx,Nmy,Nby,Nmz,Nbz,' int64 storage GiB=',requiredGiB
end SUBROUTINE SizeList
!--------------------------------------------------------------
!                          Make linker lists of particles in each cell

SUBROUTINE List
  use omp_lib, only: omp_get_max_threads,omp_get_num_threads,omp_get_thread_num
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer*2, allocatable :: zCell16(:)
  integer*4, allocatable :: zCell32(:)
  integer*8 :: jp,planeCount,i,j,k,lo,hi,cacheWords
  integer :: thread,threads
  real :: timeStart,timeFinish,memoryUsed
  logical :: compactZ,bufferedCoordinates

  timeStart=seconds()
  if(.not.ieee_is_finite(Cell).or.Cell<=0..or.Np<0_8)error stop 'Invalid BDM particle list geometry'
  if(Nmx>Nbx.or.Nmy>Nby.or.Nmz>Nbz)error stop 'Invalid BDM particle list bounds'
  ! AddBuffer validates originals and restricts image shifts to -1..1. This
  ! box/cell bound therefore proves the original fast integer conversion safe.
  ! Direct unbuffered tests/callers use the guarded clamp below instead.
  bufferedCoordinates=allocated(OriginalParticleId).and.ieee_is_finite(Box)
  if(bufferedCoordinates)bufferedCoordinates=Box>0..and. &
    2.d0*dble(Box)<dble(huge(0.)).and.2.d0*dble(Box)/dble(Cell)<dble(huge(0))
  if(omp_get_max_threads()==1.or.Np<10000_8)then
    Label=0_8
    do jp=1,Np
      if(bufferedCoordinates)then
        i=min(max(Nmx,ceiling(BdmParticleCoordinate(jp,1)/dble(Cell))-1),Nbx)
        j=min(max(Nmy,ceiling(BdmParticleCoordinate(jp,2)/dble(Cell))-1),Nby)
        k=min(max(Nmz,ceiling(BdmParticleCoordinate(jp,3)/dble(Cell))-1),Nbz)
      else
        i=ClampedCell(BdmParticleCoordinate(jp,1)/dble(Cell),Nmx,Nbx)
        j=ClampedCell(BdmParticleCoordinate(jp,2)/dble(Cell),Nmy,Nby)
        k=ClampedCell(BdmParticleCoordinate(jp,3)/dble(Cell),Nmz,Nbz)
      endif
      Lst(jp)=Label(i,j,k)
      Label(i,j,k)=jp
    enddo
  else
    planeCount=int(Nbz,8)-int(Nmz,8)+1_8
    compactZ=Nmz>=-int(huge(0_2),4)-1.and.Nbz<=int(huge(0_2),4)
    cacheWords=Np
    if(compactZ)cacheWords=Np/2_8+mod(Np,2_8)
    if(4.d0*dble(cacheWords)/1024.d0**3+dble(Memory(0_8))>dble(MaxMemory)) &
      error stop 'BDM z-cell cache exceeds configured memory limit'
    if(compactZ)then
      allocate(zCell16(Np))
    else
      allocate(zCell32(Np))
    endif
    memoryUsed=Memory(cacheWords)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(thread,threads,lo,hi,jp,i,j,k)
    thread=omp_get_thread_num()
    threads=omp_get_num_threads()
    ! Convert each exact periodic z coordinate once. Signed int16 is sufficient
    ! for ordinary grids; the int32 fallback covers every declared list bound.
!$OMP DO SCHEDULE(STATIC)
    do jp=1,Np
      if(bufferedCoordinates)then
        k=min(max(Nmz,ceiling(BdmParticleCoordinate(jp,3)/dble(Cell))-1),Nbz)
      else
        k=ClampedCell(BdmParticleCoordinate(jp,3)/dble(Cell),Nmz,Nbz)
      endif
      if(compactZ)then
        zCell16(jp)=int(k,2)
      else
        zCell32(jp)=int(k,4)
      endif
    enddo
!$OMP END DO
    ! Each cell has one writer. The slab scans keep ascending particle order,
    ! hence identical descending row-ID links, while reading only cached z cells.
    ! Int64 bounds also represent empty slabs below the minimum int32 cell.
    lo=int(Nmz,8)+planeCount*int(thread,8)/int(threads,8)
    hi=int(Nmz,8)+planeCount*(int(thread,8)+1_8)/int(threads,8)-1_8
!$OMP DO COLLAPSE(3)
    do k=int(Nmz,8),int(Nbz,8)
    do j=int(Nmy,8),int(Nby,8)
    do i=int(Nmx,8),int(Nbx,8)
      Label(i,j,k)=0_8
    enddo
    enddo
    enddo
!$OMP END DO
    if(lo<=hi)then
      do jp=1,Np
        if(compactZ)then
          k=zCell16(jp)
        else
          k=zCell32(jp)
        endif
        if(k<lo.or.k>hi)cycle
        if(bufferedCoordinates)then
          i=min(max(Nmx,ceiling(BdmParticleCoordinate(jp,1)/dble(Cell))-1),Nbx)
          j=min(max(Nmy,ceiling(BdmParticleCoordinate(jp,2)/dble(Cell))-1),Nby)
        else
          i=ClampedCell(BdmParticleCoordinate(jp,1)/dble(Cell),Nmx,Nbx)
          j=ClampedCell(BdmParticleCoordinate(jp,2)/dble(Cell),Nmy,Nby)
        endif
        Lst(jp)=Label(i,j,k)
        Label(i,j,k)=jp
      enddo
    endif
!$OMP END PARALLEL
    if(allocated(zCell16))deallocate(zCell16)
    if(allocated(zCell32))deallocate(zCell32)
    memoryUsed=Memory(-cacheWords)
  endif
  timeFinish=seconds()
  write(*,*) ' time to make list =',timeFinish-timeStart
contains
  integer*8 function ClampedCell(scaled,lower,upper) result(index)
    real*8, intent(in) :: scaled
    integer*4, intent(in) :: lower,upper
    if(.not.ieee_is_finite(scaled))error stop 'Nonfinite BDM particle cell coordinate'
    ! In an interior cell lower+1 < scaled <= upper proves that the ordinary
    ! ceiling and subtraction fit int32, even at the declared bound extrema.
    if(scaled<=dble(lower)+1.d0)then
      index=int(lower,8)
    else if(scaled>dble(upper))then
      index=int(upper,8)
    else
      index=int(ceiling(scaled)-1,8)
    endif
  end function ClampedCell
end SUBROUTINE List
!--------------------------------------------------------------
!                          Make linker lists of maxima in each cell
      SUBROUTINE ListMaxima
!--------------------------------------------------------------
        Integer*8 :: Nm,N0,N1,N2,N3,iMax,jp
        t0 =seconds()
        
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (i)
             Do jp=1,Nmaxima
                Lst(jp)=-1
             EndDo
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (i,j,k)
             Do k=Nmz,Nbz
             Do j=Nmy,Nby
             Do i=Nmx,Nbx
                Label(i,j,k)=0
             EndDo
             EndDo
          EndDo       
      Do jp=1,Nmaxima
         i=Ceiling(dble(Xmaxx(jp))/dble(Cell))-1
         j=Ceiling(dble(Ymaxx(jp))/dble(Cell))-1
         k=Ceiling(dble(Zmaxx(jp))/dble(Cell))-1
         i=MIN(MAX(Nmx,i),Nbx)
         j=MIN(MAX(Nmy,j),Nby)
         k=MIN(MAX(Nmz,k),Nbz)
         Lst(jp)      =Label(i,j,k)
         Label(i,j,k) =jp
      EndDo
      t1 = seconds()
      write(*,*) ' time to make list Max  =',t1-t0
      
    end SUBROUTINE ListMaxima
!--------------------------------------------------------------
!                eigenvalues and eigenvectors of 
!                 3x3 symmetric matrix
!                 Use Reyleigh power iterations to find
!                        maximum eigenvalue and eigenvector
!                        then use Trace(A) = sum(eigenvalues) and
!                                          product of eigenvalues = det(A)
!                        to get other two eigenvalues

  SUBROUTINE EigenValues(A,x,EigVal)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    real*8,intent(in) :: A(3,3)
    real*4,intent(out) :: EigVal(3),x(3)
    real*8 :: matrix(3,3),vectors(3,3),column(3),scale,off,angle,c,s,app,aqq,apq,u,v,value
    integer :: iteration,p,q,i,j,major,sign_component
    integer :: permutation(3),temporary
    if(.not.all(ieee_is_finite(A)))error stop 'BDM shape tensor is nonfinite'
    scale=maxval(abs(A))
    EigVal=0.;x=[1.,0.,0.]
    if(scale==0.d0)return
    matrix=0.5d0*(A/scale+transpose(A)/scale)
    vectors=0.d0
    do i=1,3
      vectors(i,i)=1.d0
    enddo
    ! Max-pivot Jacobi rotations preserve the eigenvector/eigenvalue pairing.
    ! Scaling and atan2 avoid both a fixed-seed failure and large tau overflow.
    do iteration=1,64
      p=1;q=2;off=abs(matrix(1,2))
      if(abs(matrix(1,3))>off)then
        p=1;q=3;off=abs(matrix(1,3))
      endif
      if(abs(matrix(2,3))>off)then
        p=2;q=3;off=abs(matrix(2,3))
      endif
      if(off<=16.d0*epsilon(1.d0))exit
      app=matrix(p,p);aqq=matrix(q,q);apq=matrix(p,q)
      angle=0.5d0*atan2(2.d0*apq,aqq-app)
      c=cos(angle);s=sin(angle)
      do i=1,3
        if(i==p.or.i==q)cycle
        u=matrix(i,p);v=matrix(i,q)
        matrix(i,p)=c*u-s*v;matrix(p,i)=matrix(i,p)
        matrix(i,q)=s*u+c*v;matrix(q,i)=matrix(i,q)
      enddo
      matrix(p,p)=c*c*app-2.d0*c*s*apq+s*s*aqq
      matrix(q,q)=s*s*app+2.d0*c*s*apq+c*c*aqq
      matrix(p,q)=0.d0;matrix(q,p)=0.d0
      column=vectors(:,p)
      vectors(:,p)=c*column-s*vectors(:,q)
      vectors(:,q)=s*column+c*vectors(:,q)
    enddo
    if(iteration>64)error stop 'BDM shape eigensolver failed to converge'
    permutation=[1,2,3]
    do i=1,2
      do j=i+1,3
        if(matrix(permutation(j),permutation(j))>matrix(permutation(i),permutation(i)))then
          temporary=permutation(i);permutation(i)=permutation(j);permutation(j)=temporary
        endif
      enddo
    enddo
    do i=1,3
      value=matrix(permutation(i),permutation(i))*scale
      if(abs(value)>dble(huge(EigVal)))error stop 'BDM shape eigenvalue is not representable'
      EigVal(i)=real(value,4)
    enddo
    major=permutation(1)
    sign_component=maxloc(abs(vectors(:,major)),dim=1)
    if(vectors(sign_component,major)<0.d0)vectors(:,major)=-vectors(:,major)
    x=real(vectors(:,major),4)
  end SUBROUTINE EigenValues

!--------------------------------------------------------------
!
!           Find halo concentration using M,R, and Vmax
!                 M - in Msunh, R - comoving kpch
!                 Vmax = in km/s

      real function Concentration(aM,aR,Vmax)
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        implicit none
        real*4,intent(in) :: aM,aR,Vmax
        real*8,parameter :: C0=2.162581587d0,Gkpc=4.333d-6
        real*8 :: normalization,target,virial2,lo,hi,mid,c,predicted
        integer :: iteration
        Concentration=0.
        if(.not.all(ieee_is_finite([aM,aR,Vmax,AEXPN])))return
        if(aM<=0..or.aR<=0..or.Vmax<=0..or.AEXPN<=0.)return
        virial2=Gkpc*dble(aM)/(dble(aR)*dble(AEXPN))
        target=dble(Vmax)**2/virial2
        if(target<1.d0)then
          Concentration=-1.;return
        endif
        normalization=(log(1.d0+C0)-C0/(1.d0+C0))/C0
        lo=log(C0);hi=log(dble(huge(Concentration)))
        c=exp(hi)
        predicted=normalization*c/(log(1.d0+c)-c/(1.d0+c))
        if(target>predicted)return
        ! Invert only the monotonic c >= 2.16258 NFW branch, with a finite
        ! logarithmic bracket and bounded iterations even for invalid inputs.
        do iteration=1,80
          mid=0.5d0*(lo+hi);c=exp(mid)
          predicted=normalization*c/(log(1.d0+c)-c/(1.d0+c))
          if(predicted>target)then
            hi=mid
          else
            lo=mid
          endif
          if(hi-lo<=1.d-12)exit
        enddo
        Concentration=real(exp(0.5d0*(lo+hi)),4)
      end function Concentration

!---------------------------------------------------------------------------
!              find limits for the linker-list search
      SUBROUTINE Limits(x,y,z,Radius,i1,i2,j1,j2,k1,k2)
!----------------------------------------------------------------------------
           i2=Ceiling((dble(x)+dble(Radius))/dble(Cell))-1
           j2=Ceiling((dble(y)+dble(Radius))/dble(Cell))-1
           k2=Ceiling((dble(z)+dble(Radius))/dble(Cell))-1
           i1=Ceiling((dble(x)-dble(Radius))/dble(Cell))-1
           j1=Ceiling((dble(y)-dble(Radius))/dble(Cell))-1
           k1=Ceiling((dble(z)-dble(Radius))/dble(Cell))-1
 
            i1=MIN(MAX(Nmx,i1),Nbx) 
            j1=MIN(MAX(Nmy,j1),Nby)
            k1=MIN(MAX(Nmz,k1),Nbz)
           i2=MIN(MAX(Nmx,i2),Nbx)
           j2=MIN(MAX(Nmy,j2),Nby)
           k2=MIN(MAX(Nmz,k2),Nbz)

         end SUBROUTINE Limits
!---------------------------------------------------------------------------- 
!                          Define configuration of files
!                          Read control info from the first file
!                          allocate arrays
!
      Subroutine SetParameters
     use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
     implicit none
     real :: Xscale,Vscale,Dscale
     integer :: kf,kfile

     Integer*4          :: jstep
     Integer            :: FileList(3)=(/12,20,30/)
     logical   :: exst
     character*80       :: Name,CatLabel
     Character*40       :: txt1,txt2,txt3,txt4,txt5,txt6,txt7,txt8
     Character*52       :: txt2b
     Character*10       :: sxt1,sxt2,sxt3,sxt4

           
           call ValidateParameters
           call SetOverdensity
           if (.not.ieee_is_finite(Box)) call ConfigurationError(0,'box size must be finite')
           if (Box <= 0..or.NGRID <= 0.or.NROW <= 0) &
              call ConfigurationError(0,'box, mesh size and particle-grid size must be positive')
           Xscale = Box/NGRID                 ! Scale for comoving coordinates
           Vscale = 100.*Xscale/AEXPN          ! Scale for velocities
           Dscale = 2.774e+11*(Box/NROW)**3    ! mass scale
           MassOne= Om0*Dscale                ! mass of the smallest particle


           Xleft = 0. ; Xright = Box
           Yleft = 0. ; Yright = Box
           Zleft = 0. ; Zright = Box

      write(13,*) '    Boundaries(Mpch) =',Xleft,Xright
      write(13,*) '    Boundaries       =',Yleft,Yright
      write(13,*) '    Boundaries       =',Zleft,Zright
      write(13,*) '    Buffer zone      =',dBuffer

      write(*,*) ' Scales=', Xscale, Vscale
      Do kf =1,1
                      kfile = FileList(kf)
               If(kf==2.or.kf==3)Then
                 WRITE (kfile) HEADER
                      sxt1 =' A    ='
                      sxt2 =' Step ='
                 WRITE (kfile) sxt1,AEXPN,sxt2,ASTEP 
                      sxt1 =' I    ='
                      sxt2 =' Nrow ='
                      sxt3 =' Ngrid='
                 WRITE (kfile) sxt1,ISTEP,sxt2,NROW,sxt3,NGRID 
                      sxt1 =' Omega_0='
                      sxt2 =' Omega_L='
                      sxt3 =' hubble ='
                      txt4 =' buffer width (Mpch) ='
                 WRITE (kfile) sxt1,Om0,sxt2,OmL,sxt3,hubble,txt4,dBuffer
                     txt1= ' Number of radial bins                ='
                 WRITE (kfile) txt1,NradP
                     txt1= ' Mass of smallest particle (Msunh)    ='
                 WRITE (kfile) txt1,MassOne
                     txt1= ' Overdensity limit                    ='
                 WRITE (kfile) txt1,Ovdens
                 Txt1 ='     XYZ(Mpch) '
                 Txt2 ='Vxyz(km/s)                  Mvir/Msunh '
                 Txt2b='Vxyz(km/s)                  Mbound     Mtot/Msunh '
                 Txt3 =' Rvir(kpch) Vrms(km/s) Vcirc(km/s)     '
                 Txt4 =' Nhalo  Cvir    Nparticles  Distinct/Sub'              
                 Txt5 ='   Xoff  2K/Ep-1   Lambda   RadRMS/kpch'
                 Txt6 ='  b/a  c/a MajorAxis:  x      y      z'   
                 WRITE (kfile) txt1,txt2b,txt3,txt4,txt5,txt6
               Else
                 WRITE (kfile,'(a)') trim(HEADER)//' [BDM finder v2]'
                      sxt1 =' A    ='
                      sxt2 =' Step ='
                 WRITE (kfile,'(2(a,f8.5))') sxt1,AEXPN,sxt2,ASTEP 
                      sxt1 =' I    ='
                      sxt2 =' Nrow ='
                      sxt3 =' Ngrid='
                 WRITE (kfile,'(3(a,i5))') sxt1,ISTEP,sxt2,NROW,sxt3,NGRID 
                      sxt1 =' Omega_0='
                      sxt2 =' Omega_L='
                      sxt3 =' hubble ='
                      txt4 =' buffer width (Mpch) ='
                 WRITE (kfile,'(6(a,f8.4))') sxt1,Om0,sxt2,OmL,sxt3,hubble,txt4,dBuffer
                     txt1= ' Number of radial bins                ='
                 WRITE (kfile,'(a,i4)') txt1,NradP
                     txt1= ' Mass of smallest particle (Msunh)    ='
                 WRITE (kfile,'(a,1p,g12.4)') txt1,MassOne
                     txt1= ' Overdensity limit                    ='
                 WRITE (kfile,'(a,f8.4)') txt1,Ovdens
                 Txt1 ='     XYZ(Mpch) '
                 Txt2 ='Vxyz(km/s)                  Mvir/Msunh '
                 Txt2b='Vxyz(km/s)                  Mbound     Mtot/Msunh '
                 Txt3 =' Rvir(kpch) Vrms(km/s) Vcirc(km/s)     '
                 Txt4 =' Nhalo  Cvir    Nparticles  Distinct/Sub'              
                 Txt5 ='   Xoff      2K/Ep-1   Lambda   RadRMS/kpch'
                 Txt6 ='  b/a  c/a MajorAxis:  x      y      z'   
                 WRITE (kfile,'(6a)') txt1,txt2b,txt3,txt4,txt5,txt6
              end If
             end Do

             Txt1= 'R/kpch   Npart    R/Rvir    Mass  Vcirc'
             Txt2= ' Dens/Msun/comvKpch Vrms   Vrad Vradrms'
             Txt3= ' Bound:R/kpch    Npart    R/Rvir   Mass'
             Txt4= '  Vcirc Dens/Msun/comvKpch  Vrms   Vrad'
             Txt5= '  Vradrms'
         ! WRITE(20,'(3a)') 'R/kpch   Npart    R/Rvir    Mass    Vcirc Dens/Msun/comvKpch', &

             !                 '  Vrms   Vrad Vradrms Bound:R/kpch   Npart    R/Rvir    Mass', &
         !                 '   Vcirc Dens/Msun/comvKpch  Vrms   Vrad  Vradrms'
!100      FORMAT(1X,'Header=>',A45,/     &
!               ' A=',F8.3,' Step=',F8.3,/               &
!               ' I =',I4,' Nrow=',I4,' Ngrid=',I4,/' Omega_0=',f6.3,   &
!               ' Omega_L=',f6.3,' h=',f6.3,' buffer width (Mpch) =',f9.4)
!110       FORMAT(2(a,T50,f8.3/),4(a,T50,i5/),a,T50,f8.3/,a,T50,g12.4,/4(a,T50,f8.3/))
!120       FORMAT(5x,a,T40,a,T68,a,T81,a,T92,a,T103,a,T121,a,T128,a,T137,a,       &
!                    T152,a,T166,a,T178,a,T187,a,T200,a,T213,a,T222,a,T232,a)

 
      write (13,'(6x,"Mass of smallest particle",/9x,  &
                   "in units M_sun/h is   =",3x,g10.3)')  MassOne

      write(*,*) ' Np     = ',Np,Nparticles
     

    end Subroutine SetParameters

!---------------------------------------------------------------------------- 
!         Add buffer around the computational box
!           -- 
!           -- 
!           -- move coordinates to new arrays

SUBROUTINE AddBuffer
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  real, allocatable :: Xbb(:),Ybb(:),Zbb(:),VXbb(:),VYbb(:),VZbb(:)
  integer*1, allocatable :: imageCounts(:)
  integer*8, allocatable :: blockEnd(:)
  integer*8, parameter :: blockRows=16384_8
  integer*8 :: ic,ip,originalCount,iPartMax,nx,ny,nz,block,blockCount,first,last,localImages,expectedEnd,scratchWords
  integer :: ilo,ihi,jlo,jhi,klo,khi,i,j,k
  real :: memoryUsed
  real*8 :: xx,yy,zz,box64,width64
  logical :: invalid

  if(allocated(OriginalParticleId))error stop 'BDM periodic buffer already active'
  call PrepareParticleSearch
  if(.not.ieee_is_finite(dBuffer).or.dBuffer<0..or.dBuffer>Box) &
    error stop 'BDM buffer width must be finite and between zero and Box'
  originalCount=Np
  if(Np/=Nparticles)error stop 'BDM periodic buffer original count mismatch'
  if(originalCount<0_8.or.originalCount>huge(0_8)/8_8) &
    error stop 'BDM periodic particle count exceeds safe int64 range'
  box64=dble(Box); width64=dble(dBuffer)
  ! Below half a box each axis contributes at most two images, hence 0..7
  ! ghosts per original row. Retain byte counts and only one int64 end per block.
  blockCount=originalCount/blockRows
  if(mod(originalCount,blockRows)/=0_8)blockCount=blockCount+1_8
  scratchWords=originalCount/4_8+2_8*(blockCount+1_8)
  if(mod(originalCount,4_8)/=0_8)scratchWords=scratchWords+1_8
  allocate(imageCounts(originalCount),blockEnd(0:blockCount))
  memoryUsed=Memory(scratchWords)
  invalid=.false.
!$OMP PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC) &
!$OMP PRIVATE(block,ic,first,last,localImages,xx,yy,zz,ilo,ihi,jlo,jhi,klo,khi,nx,ny,nz) REDUCTION(.or.:invalid)
  do block=1,blockCount
    first=(block-1_8)*blockRows+1_8; last=min(first+blockRows-1_8,originalCount)
    localImages=0_8
    do ic=first,last
      xx=dble(Xpar(ic)); yy=dble(Ypar(ic)); zz=dble(Zpar(ic))
      if(.not.all(ieee_is_finite([xx,yy,zz])).or.min(xx,yy,zz)<0.d0.or.max(xx,yy,zz)>=box64)then
        invalid=.true.
        imageCounts(ic)=0_1
        cycle
      endif
      ilo=ceiling((-width64-xx)/box64); ihi=floor((box64+width64-xx)/box64)
      jlo=ceiling((-width64-yy)/box64); jhi=floor((box64+width64-yy)/box64)
      klo=ceiling((-width64-zz)/box64); khi=floor((box64+width64-zz)/box64)
      nx=ihi-ilo+1_8; ny=jhi-jlo+1_8; nz=khi-klo+1_8
      if(min(nx,ny,nz)<1_8.or.max(nx,ny,nz)>2_8) &
        error stop 'BDM periodic image count outside supported half-box range'
      imageCounts(ic)=int(nx*ny*nz-1_8,1)
      localImages=localImages+int(imageCounts(ic),8)
    enddo
    blockEnd(block)=localImages
  enddo
  if(invalid)error stop 'BDM original analysis coordinates must lie in [0,Box)'
  blockEnd(0)=originalCount
  ! Prefix only the block totals. Independent blocks retain original row order
  ! during filling, including the existing k/j/i periodic-image order per row.
  do block=1,blockCount
    if(blockEnd(block)>huge(0_8)-blockEnd(block-1_8))error stop 'BDM periodic prefix exceeds int64 range'
    blockEnd(block)=blockEnd(block-1_8)+blockEnd(block)
  enddo
  iPartMax=blockEnd(blockCount)
  if(iPartMax>huge(0_8)/8_8)error stop 'BDM periodic allocation accounting exceeds int64 range'
  write(*,*) ' Allocate exact periodic particle buffer: ',originalCount,iPartMax
  allocate(Xbb(iPartMax),Ybb(iPartMax),Zbb(iPartMax),VXbb(iPartMax),VYbb(iPartMax),VZbb(iPartMax))
  allocate(OriginalParticleId(iPartMax))
  memoryUsed=Memory(8_8*iPartMax) ! six float32 fields and one int64 identity
!$OMP PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC) &
!$OMP PRIVATE(block,ic,first,last,ip,expectedEnd,xx,yy,zz,ilo,ihi,jlo,jhi,klo,khi,i,j,k)
  do block=1,blockCount
    first=(block-1_8)*blockRows+1_8; last=min(first+blockRows-1_8,originalCount)
    ip=blockEnd(block-1_8)
    do ic=first,last
      xx=dble(Xpar(ic)); yy=dble(Ypar(ic)); zz=dble(Zpar(ic))
      Xbb(ic)=Xpar(ic); Ybb(ic)=Ypar(ic); Zbb(ic)=Zpar(ic)
      VXbb(ic)=VX(ic); VYbb(ic)=VY(ic); VZbb(ic)=VZ(ic)
      OriginalParticleId(ic)=ic
      expectedEnd=ip+int(imageCounts(ic),8)
      ilo=ceiling((-width64-xx)/box64); ihi=floor((box64+width64-xx)/box64)
      jlo=ceiling((-width64-yy)/box64); jhi=floor((box64+width64-yy)/box64)
      klo=ceiling((-width64-zz)/box64); khi=floor((box64+width64-zz)/box64)
      do k=klo,khi
      do j=jlo,jhi
      do i=ilo,ihi
        if(i==0.and.j==0.and.k==0)cycle
        ip=ip+1_8
        Xbb(ip)=real(xx+dble(i)*box64); Ybb(ip)=real(yy+dble(j)*box64); Zbb(ip)=real(zz+dble(k)*box64)
        VXbb(ip)=VX(ic); VYbb(ip)=VY(ic); VZbb(ip)=VZ(ic)
        OriginalParticleId(ip)=ic
      enddo
      enddo
      enddo
      if(ip/=expectedEnd)error stop 'BDM periodic count/fill mismatch'
    enddo
    if(ip/=blockEnd(block))error stop 'BDM periodic block count/fill mismatch'
  enddo
  deallocate(imageCounts,blockEnd)
  memoryUsed=Memory(-scratchWords)
  deallocate(Xpar,Ypar,Zpar,VX,VY,VZ)
  memoryUsed=Memory(-6_8*originalCount)
  call move_alloc(Xbb,Xpar); call move_alloc(Ybb,Ypar); call move_alloc(Zbb,Zpar)
  call move_alloc(VXbb,VX); call move_alloc(VYbb,VY); call move_alloc(VZbb,VZ)
  Np=iPartMax; Nparticles=Np
end SUBROUTINE AddBuffer
!---------------------------------------------------------------------------- 
!         Remove buffer around the computational box
!
!

SUBROUTINE RemoveBuffer(NpPM)
  implicit none
  integer*8, intent(in) :: NpPM
  real :: memoryUsed
  if(.not.allocated(BdmPMX))error stop 'BDM has no retained simulation state to restore'
  if(NpPM/=BdmPMCount)error stop 'BDM restoration particle count mismatch'
  memoryUsed=Memory(-6_8*size(Xpar,kind=8))
  deallocate(Xpar,Ypar,Zpar,VX,VY,VZ)
  if(allocated(OriginalParticleId))then
    memoryUsed=Memory(-2_8*size(OriginalParticleId,kind=8))
    deallocate(OriginalParticleId)
  endif
  call move_alloc(BdmPMX,Xpar); call move_alloc(BdmPMY,Ypar); call move_alloc(BdmPMZ,Zpar)
  call move_alloc(BdmPMVX,VX); call move_alloc(BdmPMVY,VY); call move_alloc(BdmPMVZ,VZ)
  Nparticles=BdmPMCount; Np=Nparticles
  BdmPMCount=0_8
  HaloSearchRadius=0.; ParticleSearchRadius=0.
end SUBROUTINE RemoveBuffer
end Module LinkerList
