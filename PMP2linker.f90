!----------------------------------------------------------------------------
!
!          BDM halo finder   A.Klypin 2010
!
! Contains:     

! Numerical duplicate policy shared by the finder and the catalogue replay.
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
                     iVirial                  ! switch: 1 =Virial,  0=200 overdensity crit, 2=200rho_matter 

Real*4 ::            TotalMemory=0.2, t0                           ! current memory
Real*4 ::            Xleft,Xright,Yleft,Yright,Zleft,Zright,dBuffer    ! boundaries of domain
Character*120 ::     CatshortName,CatalogName,    &                ! names of halo catalogs 
                     outputName                                    ! dump file

Real*4     ::        Om0,Ovdens      ! cosmology
Real*4     ::        MassOne                                 !  current simulation
Integer*8 ::         Np                                          ! Number of particles in domain
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
      If(mDENSIT==1)Call DENSIT                 ! density on original Ng mesh
      Call FindMaxima
      

      myMemory =Memory(-1_8*NGRID*NGRID*NGRID)
      DeAllocate (FI)
      write(*,'(a,i11)') ' Go to RescaleCoords    : ',Nparticles
      Call ReadParameters(ISTEP)
      Call SetParameters
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
            myMemory= Memory(1_8*(Np+(Nbx-Nmx+1_8)*(Nby-Nmy+1_8)*(Nbz-Nmz+1_8)))
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
            myMemory= Memory(-1_8*(Np+(Nbx-Nmx+1_8)*(Nby-Nmy+1_8)*(Nbz-Nmz+1_8)))
              
      Call  SizeListMaxima
            ALLOCATE (Lst(Nmaxima),Label(Nmx:Nbx,Nmy:Nby,Nmz:Nbz))
            myMemory= Memory(1_8*(Nmaxima+(Nbx-Nmx+1_8)*(Nby-Nmy+1_8)*(Nbz-Nmz+1_8)))
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
            myMemory= Memory(-1_8*(Nmaxima+(Nbx-Nmx+1_8)*(Nby-Nmy+1_8)*(Nbz-Nmz+1_8)))
      
             DEALLOCATE (Mvir,Rvir,Xoff)
             myMemory =Memory(-3_8*Nmaxima )
             DEALLOCATE (xMaxx,yMaxx,zMaxx)
             myMemory =Memory(-3_8*Nmaxima )
             DEALLOCATE (VxMaxx,VyMaxx,VzMaxx)
             myMemory =Memory(-3_8*Nmaxima )
             DEALLOCATE(LstMax,EpotM,EkinM,LambdaM)
             myMemory =Memory(-4_8*Nmaxima )
             DEALLOCATE(VmaxM,RmaxM,Mtotal,RadRms)
             myMemory =Memory(-4_8*Nmaxima )
             DEALLOCATE(Xax,Yax,Zax)
             myMemory =Memory(-3_8*Nmaxima )
             DEALLOCATE(Axba,Axca)
             myMemory =Memory(-2_8*Nmaxima )
      !-------- restore PM structure
       Call RemoveBuffer(NpPM)
       myMemory =Memory(1_8*NGRID*NGRID*NGRID)
       Allocate (FI(NGRID,NGRID,NGRID))
             
           end SUBROUTINE BDM

!----------------------------------------------------------
!               Read/create  configuration file BDM.config       
!                    
      SUBROUTINE ReadParameters(jStep)
!----------------------------------------------------------
  Character :: txt*80,str*10,eqsign*1,Line*120,Catlabel*10
  logical   :: FileExists
                ! default set of parameters
     NradP    =   30
     iVirial  =    1
     dLogR    =  0.02
     dLogP    =  0.02
     dBuffer  = 5.
     SlopeR   = 0.20
     Rext     = 0.15
     MassMin  = 2.5e12
     
  Inquire(file='BDM.config', Exist = FileExists)
  If(FileExists)Then
     open(1,file='BDM.config')
     rewind(1)
     Read(1,'(a)')txt
     !write(*,'(a)')txt
     Do 
       Read(1,'(a)',iostat=io) Line
       If(io /= 0)exit
       !write(*,'(a)')Line
       i0 = INDEX(Line,'=')
       i1 = INDEX(Line,'!')
       !write(*,*) ' position of = and "!" in current line =',i0,i1
       If(i1.gt.i0.and.i0>0)Then  
          backspace(1)
          Read(1,*)str
          !write(*,'(a,3x,a)')str
          backspace(1)
          If(TRIM(str)=='NradP'.or.TRIM(str)=='NRADP'.or.TRIM(str)=='nradp')Then
                   Read(1,*)str,eqsign,NradP
          Else If(TRIM(str)=='Nne'  .or.TRIM(str)=='NNE'  .or.TRIM(str)=='nne')Then
                   Read(1,*)str,eqsign,Nne
          Else If(TRIM(str)=='iVirial'.or.TRIM(str)=='IVIRIAL'.or.TRIM(str)=='ivirial')Then
                   Read(1,*)str,eqsign,iVirial
          Else If(TRIM(str)=='dLogR'.or.TRIM(str)=='DLOGR'.or.TRIM(str)=='dlogr')Then
                   Read(1,*)str,eqsign,dLogR
          Else If(TRIM(str)=='Rext'.or.TRIM(str)=='Rextr'.or.TRIM(str)=='Rextern')Then
                   Read(1,*)str,eqsign,Rext
          Else If(TRIM(str)=='SlopeR'.or.TRIM(str)=='Sloper'.or.TRIM(str)=='Slope')Then
                   Read(1,*)str,eqsign,SlopeR
          Else If(TRIM(str)=='MassMin'.or.TRIM(str)=='MinMass'.or.TRIM(str)=='massmin')Then
                   Read(1,*)str,eqsign,MassMin
          Else If(TRIM(str)=='dLogP'.or.TRIM(str)=='DLOGP'.or.TRIM(str)=='dlogp')Then
                   Read(1,*)str,eqsign,dLogP
          Else
             !write(*,'(3a)') 'Unrecognized parameter: ',str,' Check spelling. I igonore it'
             Read(1,*)str
          End If
       EndIf
     EndDo
  Else
     open(1,file='BDM.config')   ! ------- create default config file

     write(1,'(a)')'! ------------- BDM configuration file. Spaces between entries are significant'
     write(1,'(a)')'!                   These are main parameters:'
     write(1,10)'iVirial', iVirial,     '! switch: 1 =Virial,  0=200 overdensity' 
     write(1,10)'NradP',   NradP,       '! Number of shells for halo profiles'
     write(1,20)'Rext',   Rext,         '! extr radius shift at m=1e15'
     write(1,20)'SlopeR',   SlopeR,     '! slope for extr radius shift'
     write(1,20)'MassMin',   MassMin,     '! Minimum halo mass'

     write(1,20)'dLogR',   dLogR,       '! size of log binning for potential'
     write(1,20)'dLogP',   dLogP,       '! size of log binning for profiles'


10   format(10x,a,T20,' = ',i6,T40,a)
20   format(10x,a,T20,' = ',1p,g10.3,T40,a)
  EndIf
     close(1)

   write(*,'(/a)') '  ------ current set of parameters:'
     write(*,'(a)')'!                   These are main parameters:'
     write(*,10)'iVirial', iVirial,     '! switch: 1 =Virial,  0=200 overdensity' 
     write(*,10)'NradP',   NradP,       '! Number of shells for halo profiles'
     write(*,20)'Rext',   Rext,         '! extr radius shift at m=1e15'
     write(*,20)'SlopeR',   SlopeR,     '! slope for extr radius shift'
     write(*,20)'MassMin',   MassMin,     '! Minimum halo mass'

     write(*,20)'dLogR',   dLogR,       '! size of log binning for potential'
     write(*,20)'dLogP',   dLogP,       '! size of log binning for profiles'

          moment = jstep
      write(outputName,'(a,i4.4,a)')'CATALOGS/outputB.',jStep,'.dat'
      !open(13,file=TRIM(outputName),buffered='NO')
      open(13,file=TRIM(outputName))
      !----------------------------  Open files ---------------------
        SELECT CASE (iVirial)
        CASE (0)
           CatLabel = 'W.' ! 200\rho_critical
        CASE  (1)
           CatLabel = 'V.' ! virial overdensity
        CASE  (2)
           CatLabel = 'M.' ! 200 matter overdensity
        CASE  (3)
           CatLabel = 'A.' ! Abacus overdensity: corrected virial
        end SELECT

      write(outputName,'(2a,2(i4.4,a))')'CATALOGS/Catshort',TRIM(CatLabel),jStep,'.',Nrealization,'.DAT'
      open(12,file=TRIM(outputName))
!      write(outputName,'(2a,i4.4,a)')'CATALOGS/Catalog',TRIM(CatLabel),jStep,'.DAT'
!         open(20,file=TRIM(outputName),form='unformatted',status ='UNKNOWN')

      end SUBROUTINE ReadParameters

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
!----------------------------------------------------------
Integer*8 :: ip

           Xscale = Box/NGRID                 ! Scale for comoving coordinates
           Vscale = 100.*Xscale/AEXPN          ! Scale for velocities
           !Dscale = 2.774e+11*(Box/NROW)**3    ! mass scale
           !MassOne= Om0*Dscale                ! mass of the smallest particle
           write(*,*) ' Inside rescale coords:', Xscale,Vscale
If(iFlag==1)Then                   ! scale to Mpc and Msun
!$OMP PARALLEL DO DEFAULT(SHARED)  PRIVATE (ip)
           Do ip =1,Np             
                Xpar(ip) = (Xpar(ip)-1.)*Xscale 
                Ypar(ip) = (Ypar(ip)-1.)*Xscale 
                Zpar(ip) = (Zpar(ip)-1.)*Xscale 
               VX(ip) =  VX(ip)*Vscale 
               VY(ip) =  VY(ip)*Vscale 
               VZ(ip) =  VZ(ip)*Vscale
            end Do
         else
!$OMP PARALLEL DO DEFAULT(SHARED)  PRIVATE (ip)
           Do ip =1,Np             
                Xpar(ip) = Xpar(ip)/Xscale+1. 
                Ypar(ip) = Ypar(ip)/Xscale+1. 
                Zpar(ip) = Zpar(ip)/Xscale+1. 
                If(Xpar(ip).ge.NGRID)Xpar(ip) =Xpar(ip)-NGRID
                If(Xpar(ip).lt.1)    Xpar(ip) =Xpar(ip)+NGRID
                If(Ypar(ip).ge.NGRID)Ypar(ip) =Ypar(ip)-NGRID
                If(Ypar(ip).lt.1)    Ypar(ip) =Ypar(ip)+NGRID
                If(Zpar(ip).ge.NGRID)Zpar(ip) =Zpar(ip)-NGRID
                If(Zpar(ip).lt.1)    Zpar(ip) =Zpar(ip)+NGRID
               VX(ip) =  VX(ip)/Vscale 
               VY(ip) =  VY(ip)/Vscale 
               VZ(ip) =  VZ(ip)/Vscale
            end Do
end If
      end SUBROUTINE RescaleCoords
!---------------------------------------------------------------------------
!                  
!                  
SUBROUTINE WriteFiles
   integer*8 :: ic    
   iHalo = 0
   MassMin = max(MassMin,20.*MassOne)
      Do ip=1,Nmaxima
         x    = xMaxx(ip);       y =   yMaxx(ip);    z = zMaxx(ip)
         If(x>=Xleft.and.x<Xright.and. &
            y>=Yleft.and.y<Yright.and. &
            z>=Zleft.and.z<Zright)Then
            if(Mvir(ip)<MassMin)cycle
            Vrms = sqrt(EkinM(ip)/Mvir(ip)*2.)
            rr   = 1.e3*Rvir(ip)
            aM   = Mvir(ip)
            vvx  = VxMaxx(ip) ;vvy  = VyMaxx(ip) ; vvz  = VzMaxx(ip)
            vvmax = VmaxM(ip)
            aNpart = Mvir(ip)/MassOne
            if(aNpart<10)cycle          !-- do not take too small halos
               iHalo = iHalo + 1
            Cvir = Concentration(aM,rr,vvmax)
            If(Cvir < 0.)Cvir = Rvir(ip)/RmaxM(ip)*2.15     ! simple estimate
            Rin   = 1.e3*RadRms(ip)
            VirRat = 2.*EkinM(ip)/EpotM(ip)-1.

           write(12,'(3f11.4,3x,3f10.2,1p,2g12.4,g12.5,26g12.4)') &
                    x,y,z,VxMaxx(ip),VyMaxx(ip),VzMaxx(ip), &
                    Mvir(ip),Mtotal(ip),1.e3*Rvir(ip),Vrms, VmaxM(ip),   & 
                    iHalo,Cvir,Mvir(ip)/MassOne,0,Xoff(ip), &
                    2.*EkinM(ip)/EpotM(ip)-1.,LambdaM(ip),1.e3*RadRms(ip),   &
                    Axba(ip),Axca(ip),Xax(ip),Yax(ip),Zax(ip)
       end If
      EndDo         ! i

      close (12)

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
! The unequal-mass host test is separate from equal-mass numerical duplicates.
! Read immutable measurements in parallel, then apply the mask after the barrier.
        implicit none
        integer*8 :: ip,jp
        integer :: i1,i2,j1,j2,k1,k2,i3,j3,k3,sx,sy,sz
        integer :: sxlo,sxhi,sylo,syhi,szlo,szhi
        real*4 :: x,y,z,xx,yy,zz,dx,dy,dz,dd,radius,tstart,tfinish
        logical, allocatable :: removeHost(:)
        tstart = seconds()
        if(Nmaxima == 0) return
        Radius = max(3.5,maxval(Rvir))
        allocate(removeHost(Nmaxima))
        removeHost=.false.
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ip,jp,x,y,z,xx,yy,zz,dx,dy,dz,dd,i1,i2,j1,j2,k1,k2,i3,j3,k3) &
!$OMP PRIVATE(sx,sy,sz,sxlo,sxhi,sylo,syhi,szlo,szhi)
      Do ip=1,Nmaxima
         If(Mvir(ip)>MassOne)Then
            x   = xMaxx(ip);   y = yMaxx(ip);    z = zMaxx(ip)
            sxlo=0;sxhi=0;sylo=0;syhi=0;szlo=0;szhi=0
            if(x-Radius<0.)sxhi=1
            if(x+Radius>=Box)sxlo=-1
            if(y-Radius<0.)syhi=1
            if(y+Radius>=Box)sylo=-1
            if(z-Radius<0.)szhi=1
            if(z+Radius>=Box)szlo=-1
            do sz=szlo,szhi
            do sy=sylo,syhi
            do sx=sxlo,sxhi
            xx=x+sx*Box; yy=y+sy*Box; zz=z+sz*Box
            Call Limits(xx,yy,zz,Radius,i1,i2,j1,j2,k1,k2)
            Do k3 =k1, k2
            Do j3 =j1, j2
            Do i3 =i1, i2
              jp =Label(i3,j3,k3)
              Do while (jp.ne.0)
                If(Mvir(ip)<Mvir(jp))Then
                  dx=x-Xmaxx(jp); dx=dx-Box*anint(dx/Box)
                  dy=y-Ymaxx(jp); dy=dy-Box*anint(dy/Box)
                  dz=z-Zmaxx(jp); dz=dz-Box*anint(dz/Box)
                  dd=dx*dx+dy*dy+dz*dz
                  If(dd.lt.Rvir(jp)**2)Then
                     removeHost(ip) = .true.
                  End If
               end If
               jp =Lst(jp)
              end do
           end do
           end do
           end Do
           end do
           end do
           end do
         end If
       End Do         ! ip
       where(removeHost) Mvir=0.
       deallocate(removeHost)
       call MergeNumericalDuplicates
       tfinish = seconds()
      write(*,'(10x,a,T50,2f10.2)') ' time for RemoveDuplicates =',tfinish-tstart,tfinish-t0

    end SUBROUTINE RemoveDuplicates

!---------------------------------------------------------------------------
! Merge connected components of the conservative strict duplicate graph.
! Measurements stay read-only until all edges have been inspected. The lowest
! original candidate index among host-test survivors wins, independent of list
! order. This catalogue-level policy is not a particle-membership assertion.
      SUBROUTINE MergeNumericalDuplicates
        implicit none
        integer :: ip,jp,i1,i2,j1,j2,k1,k2,i3,j3,k3,sx,sy,sz
        integer :: sxlo,sxhi,sylo,syhi,szlo,szhi,ri,rj
        integer, allocatable :: parent(:)
        real*4 :: x,y,z,xx,yy,zz,radius
        real*8 :: dx,dy,dz,dd,dv2
        allocate(parent(Nmaxima))
        do ip=1,Nmaxima
          parent(ip)=ip
        end do
        radius=real(DuplicateRadius)
        do ip=1,Nmaxima
          if(Mvir(ip)<=MassOne)cycle
          x=xMaxx(ip); y=yMaxx(ip); z=zMaxx(ip)
          sxlo=0;sxhi=0;sylo=0;syhi=0;szlo=0;szhi=0
          if(x-radius<0.)sxhi=1
          if(x+radius>=Box)sxlo=-1
          if(y-radius<0.)syhi=1
          if(y+radius>=Box)sylo=-1
          if(z-radius<0.)szhi=1
          if(z+radius>=Box)szlo=-1
          do sz=szlo,szhi
          do sy=sylo,syhi
          do sx=sxlo,sxhi
            xx=x+sx*Box; yy=y+sy*Box; zz=z+sz*Box
            call Limits(xx,yy,zz,radius,i1,i2,j1,j2,k1,k2)
            do k3=k1,k2
            do j3=j1,j2
            do i3=i1,i2
              jp=int(Label(i3,j3,k3))
              do while(jp/=0)
                if(jp>ip.and.Mvir(jp)>MassOne)then
                  dx=dble(x)-dble(xMaxx(jp));dx=dx-Box*anint(dx/Box)
                  dy=dble(y)-dble(yMaxx(jp));dy=dy-Box*anint(dy/Box)
                  dz=dble(z)-dble(zMaxx(jp));dz=dz-Box*anint(dz/Box)
                  dd=dx*dx+dy*dy+dz*dz
                  dv2=(dble(VxMaxx(ip))-dble(VxMaxx(jp)))**2 &
                     +(dble(VyMaxx(ip))-dble(VyMaxx(jp)))**2 &
                     +(dble(VzMaxx(ip))-dble(VzMaxx(jp)))**2
                  if(StrictDuplicate(dd,dble(Mvir(ip)),dble(Mvir(jp)), &
                         dble(Mvir(ip)/MassOne),dble(Mvir(jp)/MassOne), &
                         dble(Mtotal(ip)),dble(Mtotal(jp)),dv2))then
                    ri=DuplicateRoot(parent,ip);rj=DuplicateRoot(parent,jp)
                    parent(max(ri,rj))=min(ri,rj)
                  end if
                end if
                jp=int(Lst(jp))
              end do
            end do
            end do
            end do
          end do
          end do
          end do
        end do
        do ip=1,Nmaxima
          ri=DuplicateRoot(parent,ip)
          if(ri/=ip)Mvir(ip)=0.
        end do
        deallocate(parent)
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
        real*8 :: search,search_cap,aperture_cap,rso,aperture,grid_size
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
        search=min(dble(Cell),search_cap)

        do
          call BdmHaloGather(x,y,z,search,rows,radii)
          n=size(rows)
          if(dble(n)*mass<=threshold*search**3) exit
          if(search>=search_cap)then
            HaloStatus(ip)=HaloSearchTruncated
            return
          endif
          search=min(search+dble(Cell),search_cap)
        enddo
        call BdmHaloSortRadii(rows,radii)
        rso=0.d0
        do q=n,10,-1
          ! Between adjacent particles the enclosed mass is constant. Choose
          ! the outermost exact SO root within such an interval, including a
          ! jump onto an exactly equal-density outer particle.
          energy=(dble(q)*mass/threshold)**(1.d0/3.d0)
          if(energy<radii(q)) cycle
          if(q<n)then
            if(energy>=radii(q+1)) cycle
          endif
          if(energy>search) cycle
          rso=energy
          exit
        enddo
        if(rso<=0.d0)then
          HaloStatus(ip)=HaloTooFewParticles
          return
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
            offset=[dble(Xpar(rows(q)))-dble(x),dble(Ypar(rows(q)))-dble(y), &
                    dble(Zpar(rows(q)))-dble(z)]
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
          offset=[dble(Xpar(rows(q)))-dble(x),dble(Ypar(rows(q)))-dble(y), &
                  dble(Zpar(rows(q)))-dble(z)]
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
              Axba(ip)=real(axis_ratio(2)**slope_b); Axca(ip)=real(axis_ratio(3)**slope_c)
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

! Candidate-local workspace: two exact-count list scans, no Np-sized temporary.
      SUBROUTINE BdmHaloGather(x,y,z,radius,rows,radii)
        implicit none
        real*4, intent(in) :: x,y,z
        real*8, intent(in) :: radius
        integer*8, allocatable, intent(out) :: rows(:)
        real*8, allocatable, intent(out) :: radii(:)
        integer :: i1,i2,j1,j2,k1,k2,i,j,k,n,pass
        integer*8 :: jp
        real*8 :: distance2
        call Limits(x,y,z,nearest(real(radius),1.),i1,i2,j1,j2,k1,k2)
        do pass=1,2
          n=0
          do k=k1,k2
          do j=j1,j2
          do i=i1,i2
            jp=Label(i,j,k)
            do while(jp/=0)
              distance2=(dble(Xpar(jp))-dble(x))**2+(dble(Ypar(jp))-dble(y))**2 &
                        +(dble(Zpar(jp))-dble(z))**2
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
        integer :: first,last,parent,child,n
        integer*8 :: saved_row
        real*8 :: saved_radius
        n=size(rows)
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
        integer :: first,last,parent,child,n
        n=size(ids)
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
!---------------------------------------------------------------------------
      integer*8 :: ic,ip,i,jp 
      tstart = seconds()

!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (i)
     Do i =1, Nmaxima
        Mvir(i)   = 0.
        Rvir(i)   = 0.
        VxMaxx(i) = 0.
        VyMaxx(i) = 0.
        VzMaxx(i) = 0.
     EndDo
                              ! --- improve halo positions and halo velocities   
      Do iter =1,4
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (im,x,y,z,xc,yc,zc,xv,yv,zv,nn,i3,j3,k3,i1,i2,j1,j2,k1,k2,jp,dd,Radius,d0) 
         Do im =1, Nmaxima
            x = xMaxx(im)
            y = yMaxx(im)
            z = zMaxx(im)
            xc = 0.
            yc = 0.
            zc = 0.
            xv = 0.
            yv = 0.
            zv = 0.
            nn = 0
            Radius  = Cell*max(0.5,min(2.,log10(Xoff(im)+10.)/2.))
            d0      = Radius**2
            Call Limits(x,y,z,Radius,i1,i2,j1,j2,k1,k2)
            Do k3 =k1, k2
            Do j3 =j1, j2
            Do i3 =i1, i2
               jp =Label(i3,j3,k3)
              Do while (jp.ne.0)
                 dd = (x-Xpar(jp))**2+(y-Ypar(jp))**2+(z-Zpar(jp))**2
                     If(dd< d0) Then
                        xc = xc + Xpar(jp)
                        yc = yc + Ypar(jp)
                        zc = zc + Zpar(jp)
                        xv = xv + VX(jp)
                        yv = yv + VY(jp)
                        zv = zv + VZ(jp)
                        nn = nn + 1
                     end If
                jp =Lst(jp)
              End Do                          !  jp/= 0
            EndDo   ! i3
            EndDo   ! j3
         EndDo   ! k3
         if(nn.le.1)write(18,'(a,i9,3f9.4,3i9,g12.4)') '  no particles: ',im,x,y,z,nn,k1,k2,Xoff(im)
         !write(18,'(i9,3f9.4,6i9)') im,x,y,z,nn,k1,k2
         If(nn>1)Then
            xMaxx(im) = xc/nn
            yMaxx(im) = yc/nn
            zMaxx(im) = zc/nn
            If(xMaxx(im).lt.0.) xMaxx(im) = xMaxx(im)+Box
            If(xMaxx(im).ge.Box)xMaxx(im) = xMaxx(im)-Box
            If(yMaxx(im).lt.0.) yMaxx(im) = yMaxx(im)+Box
            If(yMaxx(im).ge.Box)yMaxx(im) = yMaxx(im)-Box
            If(zMaxx(im).lt.0.) zMaxx(im) = zMaxx(im)+Box
            If(zMaxx(im).ge.Box)zMaxx(im) = zMaxx(im)-Box
            VxMaxx(im) = xv/nn
            VyMaxx(im) = yv/nn
            VzMaxx(im) = zv/nn
            Mvir(im)   = nn*MassOne
            Rvir(im)   = Radius
         end If
      End Do
!      write(18,*) ' Iter = ',iter
!     Do im=1,Nmaxima
!        write(18,'(i9,16f10.4)') im,xMaxx(im),yMaxx(im),zMaxx(im),VxMaxx(im),VyMaxx(im),VzMaxx(im), Xoff(im)
!     end Do
     End Do         ! iter
10      tfinish = seconds()
      write(13,'(10x,a,T50,2f10.2)') ' time for FindDistinctCandidates (secs) =',tfinish-tstart,tfinish-t0

    end SUBROUTINE FindDistinctCandidates

!---------------------------------------------------------------------------
!                  
!                  
!                  
!                 
      SUBROUTINE FindMaxima
!---------------------------------------------------------------------------
  Integer*4, allocatable, dimension(:)  :: Mth
  Real*4, allocatable, dimension(:,:)  :: Xxoff,xxMax,yyMax,zzMax
  Integer*4 :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM,Nthreads

        tstart = seconds()
        Om0 = Om                                    ! Set scales
        SELECT CASE (iVirial)
        CASE (0)
           Ovdens   = 200./Om0*(Om0+(1.-Om0)*AEXPN**3) ! 200\rho_critical
        CASE  (1)
           Ovdens   = OverdenVir()                     ! virial overdensity
        CASE  (2)
           Ovdens   = 200.                             ! 200 matter overdens
        CASE  (23)
           Ovdens   = OverdenAbacus()                  ! Abacus overdens
        end SELECT

!$OMP PARALLEL
  Nthreads =OMP_GET_NUM_THREADS()
!$OMP end parallel
  write(*,*) ' Number of threads =',Nthreads

           
write(*,'(a,T40,i11,2f9.4)') ' Inside FindMaxima. Nparticles= ',Nparticles,Ovdens,Om0
Dmaxim =0.
!$OMP PARALLEL DO DEFAULT(SHARED) &         ! ------- find Maximum density
!$OMP PRIVATE (M3,M2,M1) REDUCTION(MAX:Dmaxim)
     DO M3=1,NGRID
       DO M2=1,NGRID
          DO M1=1,NGRID
             Dmaxim =max(FI(M1,M2,M3),Dmaxim)
          end DO
       end DO
    end DO
    write(*,*) ' Maximum density =',Dmaxim
    t0 = seconds()
    write(*,*) ' time for Dmax  =',t0-tstart
Nmaxima =0
!$OMP PARALLEL DO DEFAULT(SHARED) &         ! ------- count how many maxima we have
!$OMP PRIVATE (M3,M2,M1,M3p,M3m,M2p,M2m,M1p,M1m,iMax) REDUCTION(+:Nmaxima)
DO M3=1,NGRID
   M3p =M3+1 ; If(M3p>NGRID)M3p=M3p-NGRID
   M3m =M3-1 ; If(M3m<1)M3m=M3m+NGRID
       DO M2=1,NGRID
          M2p =M2+1 ; If(M2p>NGRID)M2p=M2p-NGRID
          M2m =M2-1 ; If(M2m<1)M2m=M2m+NGRID
          DO M1=1,NGRID
             M1p =M1+1 ; If(M1p>NGRID)M1p=M1p-NGRID
             M1m =M1-1 ; If(M1m<1)M1m=M1m+NGRID
             iMax = 1                   ! look for all 26 neighbors
	        If(FI(M1,M2,M3).lt.FI(M1m,M2m,M3m))iMax=0 
	        If(FI(M1,M2,M3).lt.FI(M1 ,M2m,M3m))iMax=0 
	        If(FI(M1,M2,M3).lt.FI(M1p,M2m,M3m))iMax=0 
	        If(FI(M1,M2,M3).lt.FI(M1m,M2,M3m))iMax=0 
	        If(FI(M1,M2,M3).lt.FI(M1, M2,M3m))iMax=0 
	        If(FI(M1,M2,M3).lt.FI(M1p,M2,M3m))iMax=0 
	        If(FI(M1,M2,M3).lt.FI(M1m,M2p,M3m))iMax=0 
	        If(FI(M1,M2,M3).lt.FI(M1, M2p,M3m))iMax=0 
	        If(FI(M1,M2,M3).lt.FI(M1p,M2p,M3m))iMax=0 

                If(FI(M1,M2,M3).lt.FI(M1m,M2m,M3))iMax=0 
	        If(FI(M1,M2,M3).lt.FI(M1 ,M2m,M3))iMax=0 
	        If(FI(M1,M2,M3).lt.FI(M1p,M2m,M3))iMax=0 
	        If(FI(M1,M2,M3).lt.FI(M1m,M2,M3))iMax=0 
	        If(FI(M1,M2,M3).lt.FI(M1p,M2,M3))iMax=0 
	        If(FI(M1,M2,M3).lt.FI(M1m,M2p,M3))iMax=0 
	        If(FI(M1,M2,M3).lt.FI(M1, M2p,M3))iMax=0 
	        If(FI(M1,M2,M3).lt.FI(M1p,M2p,M3))iMax=0 

                If(FI(M1,M2,M3).lt.FI(M1m,M2m,M3p))iMax=0 
	        If(FI(M1,M2,M3).lt.FI(M1 ,M2m,M3p))iMax=0 
	        If(FI(M1,M2,M3).lt.FI(M1p,M2m,M3p))iMax=0 
	        If(FI(M1,M2,M3).lt.FI(M1m,M2,M3p))iMax=0 
	        If(FI(M1,M2,M3).lt.FI(M1, M2,M3p))iMax=0 
	        If(FI(M1,M2,M3).lt.FI(M1p,M2,M3p))iMax=0 
	        If(FI(M1,M2,M3).lt.FI(M1m,M2p,M3p))iMax=0 
	        If(FI(M1,M2,M3).lt.FI(M1, M2p,M3p))iMax=0 
	        If(FI(M1,M2,M3).lt.FI(M1p,M2p,M3p))iMax=0 
         If(iMax==1.and.FI(M1,M2,M3)>Ovdens/3.)Then
            Nmaxima = Nmaxima +1
            FI(M1,M2,M3) = FI(M1,M2,M3) +Dmaxim+1.  ! assign large number to the maximum
         end If
	  END DO
       END DO
    END DO
    t1 = seconds()
    write(*,*) ' time for Find Max    =',t1-t0
    write(*,*) ' FindMaxima: 1 finished: Nmaxima= ',Nmaxima
    
    Nbuff = Nmaxima/Nthreads *5  ! length of the buffer
    Allocate(Mth(Nbuff),xxMax(Nbuff,Nthreads),yyMax(Nbuff,Nthreads),zzMax(Nbuff,Nthreads))
    Allocate(Xxoff(Nbuff,Nthreads))
      if(Nmaxima == 0)Stop ' No density maxima found'
             ALLOCATE (Mvir(Nmaxima),Rvir(Nmaxima),Xoff(Nmaxima))
             myMemory =Memory(3_8*Nmaxima )
             ALLOCATE (xMaxx(Nmaxima),yMaxx(Nmaxima),zMaxx(Nmaxima))
             myMemory =Memory(3_8*Nmaxima )
             ALLOCATE (VxMaxx(Nmaxima),VyMaxx(Nmaxima),VzMaxx(Nmaxima))
             myMemory =Memory(3_8*Nmaxima )
             ALLOCATE(LstMax(Nmaxima),EpotM(Nmaxima),EkinM(Nmaxima),LambdaM(Nmaxima))
             myMemory =Memory(4_8*Nmaxima )
             ALLOCATE(VmaxM(Nmaxima),RmaxM(Nmaxima),Mtotal(Nmaxima),RadRms(Nmaxima))
             myMemory =Memory(4_8*Nmaxima )
             ALLOCATE(Xax(Nmaxima),Yax(Nmaxima),Zax(Nmaxima))
             myMemory =Memory(3_8*Nmaxima )
             ALLOCATE(Axba(Nmaxima),Axca(Nmaxima))
             myMemory =Memory(2_8*Nmaxima )
            
    Nmaxx =0
    xs     = Box/NGRID
    Mth(:) = 0
    Mm     = 1
!$OMP PARALLEL DO PRIVATE(M1,M2,M3,Mthread) Firstprivate(Mm) 
    DO M3=1,NGRID          !--- make list of maxima
       DO M2=1,NGRID
          DO M1=1,NGRID
             If(FI(M1,M2,M3)>Dmaxim+1.)Then
                Mthread       = OMP_GET_THREAD_NUM()+1
                Mth(Mthread)  = Mm
                If(Mm>Nbuff)Then
                  write(*,'(3i6,i11,i4)') M1,M2,M3,Mm,Mthread
                  Stop ' Too many elements for buffer XX'
                end If
                xxMax(Mm,Mthread) = (M1-1)*xs 
                yyMax(Mm,Mthread) = (M2-1)*xs 
                zzMax(Mm,Mthread) = (M3-1)*xs
                Xxoff(Mm,Mthread)  = FI(M1,M2,M3) -Dmaxim-1.
                Mm = Mm +1
             end If
          end DO
       end DO
    end DO
  M = 0
  Do i=1,Nthreads
     M = M + Mth(i)
  end Do
  write(*,*) ' Elements on all buffers =',M,Nmaxima
  i = 0
  Do Mthread=1,Nthreads
     Mnow = Mth(Mthread)
     Do j = 1,Mnow
        Xoff(j+i) = Xxoff(j,Mthread)
        xMaxx(j+i) = xxMax(j,Mthread)
        yMaxx(j+i) = yyMax(j,Mthread)
        zMaxx(j+i) = zzMax(j,Mthread)
     End Do
     i = i+ Mnow
  end Do
  Deallocate(Xxoff,Mth,xxMax,yyMax,zzMax)
  

    t2 = seconds()
    iOverdens  = 0
    i300 = 0; i400 =0; i500 = 0; i600 =0; i700 =0
    write(*,*) ' time for List Max    =',t2-t1
                       !--- statistics of maxima
       Do i=1,Nmaxima
          If(Xoff(i) .ge. Ovdens)iOverdens = iOverdens +1
          If(Xoff(i)>300.)Then
                i300 = i300 +1
             if(Xoff(i)>400.)Then
                i400 = i400 +1
             if(Xoff(i)>500.)Then
                i500 = i500 +1
             if(Xoff(i)>600.)Then
                i600 = i600 +1
                if(Xoff(i)>700.)Then
                  i700 = i700 +1
                End if
             end if
             end if
             end if
          end If
          enddo
          write(*,'(a,i8,/a,i9,5(a,i8))') ' Number of maxima:',Nmaxima, &
               ' Above Overdensity =',iOverdens, &
               ' Above: 300 =',i300,' 400 =',i400,' 500 =',i500,' 600 =',i600,' 700 =',i700
       tfinish = seconds()
      write(13,'(10x,a,T50,2f10.2)') ' time for FindMaxima(secs) =',tfinish-tstart,tfinish-t0
      write(*,'(10x,a,T50,2f10.2)') ' time for FindMaxima(secs) =',tfinish-tstart,tfinish-t0
    end SUBROUTINE FindMaxima
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
      SUBROUTINE SizeList
!---------------------------------------------------------------------------
      Real*4 ::  MemList,MemMaxList,  &
                       fractionMemory  = 0.90  ! fraction of memory allocated to Lists

      !Roptimal         = (Nne*MassOne/(Om0*Ovdens))**(1./3.)*9.50e-5  ! Mpch inits
      !Roptimal         = max(Roptimal,Box/NROW) * 1.2

      MemMaxList  = (MaxMemory - TotalMemory-12.*Np/1024.**3)*fractionMemory
      If(MemMaxList<0.1)Then

              STOP ' Stop: Not enough memory !!!'
      End If
      Cell                  = (Box/NGRID*2.)  !Roptimal
      k    = 0
      write(13,'(2(a,g12.4))') '    SizeList: Roptimal=',Cell,' Allowed Memory=',MemMaxList
      Do  
         Nmx         = (Xleft  -dBuffer)/Cell - 1
         Nmy         = (Yleft  -dBuffer)/Cell - 1
         Nmz         = (Zleft  -dBuffer)/Cell - 1
         Nbx         = (Xright +dBuffer)/Cell + 1
         Nby         = (Yright +dBuffer)/Cell + 1
         Nbz         = (Zright +dBuffer)/Cell + 1
         MemList  =  (3.*4.*float(Nbx-Nmx+1)*float(Nby-Nmy+1)*float(Nbz-Nmz+1))/1024.**3
         If(MemList < MemMaxList)Exit
         k = k+1
         Cell  = Cell *1.1
         If(k > 100)STOP ' Stop: too many iterations in SizeList'
      EndDo
      write(13,'(a,g12.3,a,6i6)') '     Size List:  Cell=',Cell, &
                          ' Limits=',Nmx,Nbx,Nmy,Nby,Nmz,Nbz

    end SUBROUTINE SizeList
!--------------------------------------------------------------
!                          Make linker lists of particles in each cell
      SUBROUTINE List
!--------------------------------------------------------------
        Integer*8 :: Nm,N0,N1,N2,N3,iMax,jp
        Integer*4  :: Ndiv(0:1000)
        Integer*4 :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM,Nthreads

!$OMP PARALLEL
        Nthreads =OMP_GET_NUM_THREADS()
!$OMP end parallel
        Nsplit = Nthreads  
        d = (Nbx-Nmx)/float(Nsplit)    !--- parallelization setup
        Do i=0,Nsplit
           Ndiv(i) =  Floor(d*i+Nmx+0.5)
        End Do
        Ndiv(Nsplit) = Nbx+1
        t0 =seconds()
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (i)
             Do jp=1,Np
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
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (jp,i,j,k,Isplit)
    Do Isplit =1,Nsplit
      Do jp=1,Np
         i=Ceiling(Xpar(jp)/Cell)-1
         j=Ceiling(Ypar(jp)/Cell)-1
         k=Ceiling(Zpar(jp)/Cell)-1
         i=MIN(MAX(Nmx,i),Nbx)
         j=MIN(MAX(Nmy,j),Nby)
         k=MIN(MAX(Nmz,k),Nbz)
         If(k.ge.Ndiv(Isplit-1).and.k.lt.Ndiv(Isplit))Then
            Lst(jp)      = Label(i,j,k)
            Label(i,j,k) = jp
         End If
      EndDo
   end Do
      t1 = seconds()
      write(*,*) ' time to make list  =',t1-t0
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
         i=Ceiling(Xmaxx(jp)/Cell)-1
         j=Ceiling(Ymaxx(jp)/Cell)-1
         k=Ceiling(Zmaxx(jp)/Cell)-1
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
    Integer*4, Parameter :: Nmax =200
    Real*8, Parameter :: error =2.d-7
    Real*4  :: EigVal(3),x(3)
    Real*8 :: A(3,3),xx(3),y(3),Eig,EigNew,a1,b1,Trace,detA,x1,x2
    Integer*4 :: ind(3)
    xx       = 1.     ! initial guess
    Eig    = 0.
    Do kstep=1,Nmax
       y            = MATMUL(A,xx)
       EigNew = max(y(1),y(2),y(3))
       xx     = y/EigNew
       !write(13,'(10x,i4,1p,2g13.5,3x,4g13.5)')kstep,EigNew,Eig,xx
       If(abs(EigNew-Eig)< error*abs(EigNew))exit
       Eig = EigNew
    End Do
    Eig = EigNew
    Trace  = A(1,1) +A(2,2) +A(3,3)
    detA    = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))        &
                   -A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1))        &
                   +A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
                                        ! solve quadratic equation for other eigenvalues
    a1 = (Trace-Eig)/2.
    b1 =  detA/Eig
    d = sqrt(max(a1**2-b1,1.d-20))
    x1 = a1 + d
    x2 = a1 - d
    !EigVal(1) = max(Eig,x1,x2)
    !EigVal(3) = min(Eig,x1,x2)
    EigVal = -1.e30
    ind = 0
    If(Eig.ge.x1.and.Eig.ge.x2)Then ! Eig is the max
       EigVal(1) = Eig
       If(x1.ge.x2)Then
          EigVal(2) = x1
          EigVal(3) = x2
       Else
          EigVal(2) = x2
          EigVal(3) = x1
       EndIf
    Else  If(x1.ge.x2)Then          ! x1 is max
       EigVal(1) = x1
       If(Eig.ge.x2)Then
          EigVal(2) = Eig
          EigVal(3) = x2
       Else
          EigVal(2) = x2
          EigVal(3) = Eig
       EndIf
    Else                           ! x2 is max
       EigVal(1) = x2
       If(Eig.ge.x1)Then
          EigVal(2) = Eig
          EigVal(3) = x1
       Else
          EigVal(2) = x1
          EigVal(3) = Eig
       EndIf
    End If


!    If(EigVal(2)>EigVal(1).or.EigVal(3)>EigVal(2))Then
!          write(13,'(/5x,1p,10g13.5)')Trace,detA,a1,b1,x1,x2,EigVal
!          write(13,'(1p,3g12.4)') A
!          write(13,'(a,1p,g13.5)') '   Sum  =',EigVal(1)+EigVal(2)+EigVal(3)
!          write(13,'(a,1p,g13.5)') '   Prod =',EigVal(1)*EigVal(2)*EigVal(3)
!          write(13,'(a,1p,2g13.5)') '   Iter =',EigNew-Eig,Eig
!          write(13,'(10x,3i3)') ind
!    EndIf
    x               = xx
    d               = sqrt(SUM(x**2))
    x      =  x/d
  end SUBROUTINE EigenValues

!--------------------------------------------------------------
!
!           Find halo concentration using M,R, and Vmax
!                 M - in Msunh, R - comoving kpch
!                 Vmax = in km/s
      real function Concentration(aM,aR,Vmax)
!     ------------------------
!
        Real*8, parameter :: d0 =1.d0, dCmin =3.d-5, C0= 2.162d0
        Real*8 :: M,V,R,A,C,Fc,Vvir,Vratio,dC 

        If(aM<1.e-6.or.aR<1.e-6.or.Vmax <1.e-6)Then
           Concentration = 0. ; return
        End If
        Vvir = 2.076d-3*sqrt(aM/(aR*AEXPN))
        Vratio = Vmax/Vvir
        If(Vratio<1.)Then
           Concentration = -1.
        Else
           dC = 2.
           C= 2.*C0
           A  = Vratio**2
           !N = 0
           Do while (abs(dC)>dCmin)
              Fc = 0.2162166_8/(log(1.d0+C)/C -1.d0/(1.d0+C))
              If(Fc.gt.A.and.C.ge.C0)Then
                 dC = dC/2.d0
                 C = C - dC
              Else
                 C = C + dC
              EndIf
              !N = N +1
              !write(*,'(i6,1p,4G13.6)')N,C,dC,Fc,A
           End Do
           Concentration = C
        End If
      End function Concentration

!---------------------------------------------------------------------------
!              find limits for the linker-list search
      SUBROUTINE Limits(x,y,z,Radius,i1,i2,j1,j2,k1,k2)
!----------------------------------------------------------------------------
           i2=Ceiling((x+Radius)/Cell)-1  
           j2=Ceiling((y+Radius)/Cell)-1  
           k2=Ceiling((z+Radius)/Cell)-1 
           i1=Ceiling((x-Radius)/Cell)-1  
           j1=Ceiling((y-Radius)/Cell)-1  
           k1=Ceiling((z-Radius)/Cell)-1 
 
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

     Integer*4          :: jstep
     Integer            :: FileList(3)=(/12,20,30/)
     logical   :: exst
     character*80       :: Name,CatLabel
     Character*40       :: txt1,txt2,txt3,txt4,txt5,txt6,txt7,txt8
     Character*52       :: txt2b
     Character*10       :: sxt1,sxt2,sxt3,sxt4

           
           Om0 = Om                                    ! Set scales
           Xscale = Box/NGRID                 ! Scale for comoving coordinates
           Vscale = 100.*Xscale/AEXPN          ! Scale for velocities
           Dscale = 2.774e+11*(Box/NROW)**3    ! mass scale
           MassOne= Om0*Dscale                ! mass of the smallest particle

        SELECT CASE (iVirial)
        CASE (0)
           Ovdens   = 200./Om0*(Om0+(1.-Om0)*AEXPN**3) ! 200\rho_critical
        CASE  (1)
           Ovdens   = OverdenVir()                     ! virial overdensity
        CASE  (2)
           Ovdens   = 200.                             ! 200 matter overdens
        CASE  (3)
           Ovdens   = OverdenAbacus()                 ! Abacus overdens
        end SELECT
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
                 WRITE (kfile) sxt1,Om0,sxt2,1.-Om0,sxt3,hubble,txt4,dBuffer 
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
                 WRITE (kfile,'(a)') HEADER
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
                 WRITE (kfile,'(6(a,f8.4))') sxt1,Om0,sxt2,Oml0,sxt3,hubble,txt4,dBuffer 
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
    
    Subroutine AddBuffer

     Integer*4          :: jstep
     Integer            :: FileList(3)=(/12,20,30/)
     Real*4,       ALLOCATABLE,   DIMENSION(:) ::               &
                     Xbb,   Ybb,  Zbb,       &   !   coords
                     VXbb, VYbb, Vzbb            !    velocities
     Integer*8 :: idummy,ic,ip,iPartMax,i

     Factor   = (1.+4.*dBuffer/Box)**3
     iPartMax = Factor*Np               ! reserve extra space for total N particles
     Nparticles  = Np                   ! old number of particles
     write(*,*) ' Allocate buffers for particles. N=',iPartMax
     myMemory= Memory(6_8*iPartMax)
       ALLOCATE(Xbb(iPartMax),Ybb(iPartMax),Zbb(iPartMax))
       ALLOCATE(VXbb(iPartMax),VYbb(iPartMax),VZbb(iPartMax))
           Xleft = 0. ; Xright = Box
           Yleft = 0. ; Yright = Box
           Zleft = 0. ; Zright = Box

!$OMP PARALLEL DO DEFAULT(SHARED)  PRIVATE (ic)
       Do ic=1,iPartMax
          VXbb(ic)= 0.
          VYbb(ic)= 0.
          VZbb(ic)= 0.
          Xbb(ic)= 0.
          Ybb(ic)= 0.
          Zbb(ic)= 0.          
       End Do
       
       ip = 0
       Xl = Xleft  -dBuffer ; Xr = Xright +dBuffer
       Yl = Yleft  -dBuffer ; Yr = Yright +dBuffer
       Zl = Zleft  -dBuffer ; Zr = Zright +dBuffer
       write(*,*) '    Replicate coordinates'
!$OMP PARALLEL DO DEFAULT(SHARED)  PRIVATE (ic,xx,yy,zz)
       Do ic=1,Np
          !xx = Xp(ic) ; yy =Yp(ic) ; zz =Zp(ic)
          !Xb(ic)  =xx ; Yb(ic) = yy ; Zb(ic) = zz
          Xbb(ic)= Xpar(ic)
          Ybb(ic)= Ypar(ic)
          Zbb(ic)= Zpar(ic)
       End Do
       write(*,*) '    Replicate velocities'
!$OMP PARALLEL DO DEFAULT(SHARED)  PRIVATE (ic,xx,yy,zz)
       Do ic=1,Np
          VXbb(ic)= VX(ic)
          VYbb(ic)= VY(ic)
          VZbb(ic)= VZ(ic)
       End Do
              write(*,*) '    add buffer'

       ip = Np
!!$OMP PARALLEL DO DEFAULT(SHARED)  PRIVATE (ic,xx,yy,zz,x,y,z,k,j,i)       
       Do ic=1,Np
          xx = Xpar(ic) ; yy =Ypar(ic) ; zz =Zpar(ic)
          do k = -1,1
             z = zz+k*Box
             if(z<Zl.or.z>Zr)cycle
               do j = -1,1
                 y = yy+j*Box
                 if(y<Yl.or.y>Yr)cycle
                   do i = -1,1
                      x = xx+i*Box
                      if(x<Xl.or.x>Xr)cycle
                      if((x<Xleft.or.x>Xright) .or.    &
                           (y<Yleft.or.y>Yright) .or.    &
                           (z<Zleft.or.z>Zright))Then
                      ip = ip +1
                      If(ip>iPartMax)Stop ' Too many buffer particles. Increase Factor in AddBuffer'
                      Xbb(ip) = x
                      Ybb(ip) = y
                      Zbb(ip) = z 
                      VXbb(ip)= VX(ic)
                      VYbb(ip)= VY(ic)
                      VZbb(ip)= VZ(ic)
                   end if
                  end do
              end do
           end do
        end Do
      write(*,*) ' New number of particles = ',ip
      myMemory=Memory(-6_8*Np)  
      DEALLOCATE(Xpar,Ypar,Zpar)  
      DEALLOCATE(VX,VY,VZ)
      Np = ip
      Nparticles = Np
      myMemory=Memory(6_8*Np)  
      ALLOCATE(Xpar(Np),Ypar(Np),Zpar(Np))  
      ALLOCATE(VX(Np),VY(Np),VZ(Np))  

!$OMP PARALLEL DO DEFAULT(SHARED)  PRIVATE (ip)
      Do ip =1,Np
         Xpar(ip) = Xbb(ip) ; Ypar(ip) = Ybb(ip) ; Zpar(ip) = Zbb(ip)
         VX(ip)= VXbb(ip); VY(ip)= VYbb(ip); VZ(ip)= VZbb(ip)
      EndDo
      myMemory= Memory(-6_8*iPartMax)  
      DEALLOCATE(Xbb,Ybb,Zbb)  
      DEALLOCATE(VXbb,VYbb,VZbb)


      xmin = 1.e12; xmax = -1.e12
      ymin = 1.e12; ymax = -1.e12
      zmin = 1.e12; zmax = -1.e12

!$OMP PARALLEL DO DEFAULT(SHARED)  PRIVATE (ip)  &
!$OMP REDUCTION(MAX:xmax,ymax,zmax) REDUCTION(MIN:xmin,ymin,zmin) 
           Do ip =1,Np
              xmin =MIN(xmin,Xpar(ip))
              ymin =MIN(ymin,Ypar(ip))
              zmin =MIN(zmin,Zpar(ip))
              xmax =MAX(xmax,Xpar(ip))
              ymax =MAX(ymax,Ypar(ip))
              zmax =MAX(zmax,Zpar(ip))
              
              If(Xpar(ip).lt.Xl.or.Ypar(ip).lt.Yl.or.Zpar(ip).lt.Zl)&
                write(*,'(a,i10,1p,3g14.5)')' Error coord: ',i,Xpar(i),Ypar(i),Zpar(i)
           EndDo
!     write(*,'(a,20es13.5)') ' X= ',(Xpar(i),i=1,10)    
!     write(*,'(a,20es13.5)') ' Y= ',(Ypar(i),i=1,10)    
!     write(*,'(a,20es13.5)') ' Z= ',(Zpar(i),i=1,10)
!     write(*,'(a,20es13.5)') ' X= ',(Xpar(i),i=Np-9,Np)    
!     write(*,'(a,20es13.5)') ' Y= ',(Ypar(i),i=Np-9,Np)    
!     write(*,'(a,20es13.5)') ' Z= ',(Zpar(i),i=Np-9,Np)
     
     write(*,'(a,2es13.5)') ' X range= ',xmin,xmax
     write(*,'(a,2es13.5)') ' Y range= ',ymin,ymax
     write(*,'(a,2es13.5)') ' Z range= ',zmin,zmax
    end Subroutine AddBuffer
!---------------------------------------------------------------------------- 
!         Remove buffer around the computational box
!
!
    Subroutine RemoveBuffer(NpPM)
     Real*4,       ALLOCATABLE,   DIMENSION(:) ::               &
                     Xbb,   Ybb,  Zbb,       &   !   coords
                     VXbb, VYbb, Vzbb            !    velocities
     Integer*8 :: NpPM,ic

           Xscale = Box/NGRID                 ! Scale for comoving coordinates
           Vscale = 100.*Xscale/AEXPN          ! Scale for velocities
     
     myMemory= Memory(6_8*Nparticles)
       ALLOCATE(Xbb(Nparticles),Ybb(Nparticles),Zbb(Nparticles))
       ALLOCATE(VXbb(Nparticles),VYbb(Nparticles),VZbb(Nparticles))
       !--- copy data to buffers
       write(*,*) ' Copy to new buffer ',NpPM,Xscale
!$OMP PARALLEL DO DEFAULT(SHARED)  PRIVATE (ic)
       Do ic=1,NpPM
          Xbb(ic)= Xpar(ic)
          Ybb(ic)= Ypar(ic)
          Zbb(ic)= Zpar(ic)
          VXbb(ic)= VX(ic)
          VYbb(ic)= VY(ic)
          VZbb(ic)= VZ(ic)          
       End Do
     myMemory= Memory(-6_8*Nparticles)
     DEALLOCATE(Xpar,Ypar,Zpar,VX,VY,VZ)
     
                               !--- restore data structur, rescale Coords and Veloc
     myMemory= Memory(6_8*NpPM)
     ALLOCATE(Xpar(NpPM),Ypar(NpPM),Zpar(NpPM),VX(NpPM),VY(NpPM),VZ(NpPM))
       write(*,*) ' Copy to  old arrays ',Xscale,Vscale
!$OMP PARALLEL DO DEFAULT(SHARED)  PRIVATE (ic)
       Do ic=1,NpPM
          Xpar(ic) =Xbb(ic) 
          Ypar(ic) =Ybb(ic) 
          Zpar(ic) =Zbb(ic) 
           VX(ic)  =VXbb(ic)
           VY(ic)  =VYbb(ic)
           VZ(ic)  =VZbb(ic)         
        End Do
       DEALLOCATE(Xbb,Ybb,Zbb)
       DEALLOCATE(VXbb,VYbb,VZbb)        
       myMemory= Memory(-6_8*Nparticles)
       Nparticles = NpPM
       Np         = Nparticles
        write(*,*) ' Rescale ',Nparticles
!$OMP PARALLEL DO DEFAULT(SHARED)  PRIVATE (ic)
       Do ic=1,Nparticles      
           Xpar(ic) = Xpar(ic)/Xscale+1. 
           Ypar(ic) = Ypar(ic)/Xscale+1. 
           Zpar(ic) = Zpar(ic)/Xscale+1. 
           If(Xpar(ic).ge.NGRID)Xpar(ic) =Xpar(ic)-NGRID
           If(Xpar(ic).lt.1)    Xpar(ic) =Xpar(ic)+NGRID
           If(Ypar(ic).ge.NGRID)Ypar(ic) =Ypar(ic)-NGRID
           If(Ypar(ic).lt.1)    Ypar(ic) =Ypar(ic)+NGRID
           If(Zpar(ic).ge.NGRID)Zpar(ic) =Zpar(ic)-NGRID
           If(Zpar(ic).lt.1)    Zpar(ic) =Zpar(ic)+NGRID
           VX(ic) =  VX(ic)/Vscale 
           VY(ic) =  VY(ic)/Vscale 
           VZ(ic) =  VZ(ic)/Vscale
        End Do
       
end Subroutine RemoveBuffer
end Module LinkerList
