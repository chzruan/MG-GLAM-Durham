!-----------------------------------------------------
!     Module Tools for PMP2 code
!	 - Main shared variable and arrays
!        - Routines: ReadSetup, ReadDataPM,
!                  seconds, Timing, TimingMain, Memory
!-------------------------------------------------
!
Module Tools
      Real*4    :: Om,        &      ! cosmology
                   Omb,       &
                   OmL,       &
                   AEXPN,     &
                   ASTEP,     &
                   sigma8,    &
                   Box,       &
                   hubble
      Integer*4 :: NROW,    &      ! number of particles in 1D
                   NGRID,   &      ! number of grid points in 1D
                   Nout,    &      ! number of outputs/mocks
                   Nbiaspars,  &   ! number of bias parameters
                   Ncheckpoint, &  !'Steps between checkpoints'
                   iSave,      &   !'Save snapshots           '
                   iPower,     &   !'DM power spectrum        '
                   iPowerRSD,  &   !'Redshift distortions     '
                   iDensityDistr, &!'DM PDF                   '
                   iBias ,     &   !'Biasing model            '
                   iWriteMock,  &  !'Write Mocks              '
                   iDumpPart,  &   !'Dump random fraction     '
                   iBDM            !'Find BDM halos           '
         
      Integer*8 :: Nparticles, Ngalaxies
      Real*4    :: zinit,da,zfinal
      Real*4    :: densThr,sigV
      Real*4    :: zout(1000),BiasPars(1000)

      Real*4    :: AEXPN0,ASTEP0,AMPLT,EKIN,EKIN1,EKIN2,AEU0,  &
                   TINTG,                 &
                   extras(100),ENKIN,ENPOT
      Integer*4 :: iStep,Nrealization,Nseed
      Real*4, Allocatable, Dimension(:,:,:), target :: FI       ! Baojiu: single precision
      Real*4, Allocatable, Dimension(:,:,:), target :: FI2      ! MG array FI2: single precision
      Real*4, Allocatable, Dimension(:,:,:), target :: FI3      ! MG array FI3: single precision
      integer :: levelmax,levelmin                              ! global vars
      Complex*8, Allocatable, Dimension(:,:,:), target :: FIC
      Real*4, allocatable,  Dimension(:) ::  Xb,Yb,Zb,VXb,Vyb,Vzb
      Real*4, allocatable,  Dimension(:) ::  XPAR,YPAR,ZPAR,VX,VY,VZ
      Real*4, allocatable,  Dimension(:) ::  dens,RandP,dens2
      Real*4       :: StepFactor = 3.e-2
      Character*45 :: HEADER
      Real*8       :: CPU(0:10) = 0.      

! MG - MG flag and models
      ! Logical   :: MG_flag = .FALSE. ! Set to true to simulate modified gravity
      ! Logical   :: MG_test = .FALSE. ! Set to true for tests
      Integer*4   :: MG_flag = 0 ! Set to 1 to simulate modified gravity
      Integer*4   :: MG_test = 0 ! Set to 1 for tests

      Integer*4 :: MG_model          ! 1 - f(R) gravity
                                     ! 2 - DGP model
                                     ! 3 - symmetron model
                                     ! 4 - k-mouflage model
                                     ! 5 - coupled scalar field model

! MG - Multigrid methods
      ! Logical   :: Vcycle = .FALSE.
      ! Logical   :: Fcycle = .FALSE.
      ! Logical   :: Wcycle = .FALSE.
    Integer*4  :: Vcycle = 0
    Integer*4  :: Fcycle = 0
    Integer*4  :: Wcycle = 0

      Integer*4 :: iter_count_max = 4
      Real*8    :: res_conv = 1.0D-8

! f(R) model parameters
      Integer*4:: fr_n
      Real*8   :: fR0
! DGP model parameters
      Real*8   :: H0rc!,phi_mean,S_mean
      ! Logical  :: N_branch = .FALSE.
      ! Logical  :: S_branch = .FALSE.
      Integer*4  :: N_branch = .FALSE.
      Integer*4  :: S_branch = .FALSE.

! symmetron model parameters
      Real*8   :: sym_astar
      Real*8   :: sym_xi
      Real*8   :: sym_beta
! kmouflage model parameters
      ! Logical  :: power_law = .FALSE.
      ! Logical  :: born_infeld = .FALSE.
      Integer*4  :: power_law = .FALSE.
      Integer*4  :: born_infeld = .FALSE.

      Integer*4:: kmf_n
      Real*8   :: kmf_K0
      Real*8   :: kmf_beta
      Real*8   :: kmf_lambda
! coupled scalar field model parameters
      Integer*4:: csf_coupling
      Integer*4:: csf_potential
      Real*8   :: csf_alpha
      Real*8   :: csf_beta
      Real*8   :: csf_lambda

  Contains
!--------------------------------------------
subroutine ReadSetup
!--------------------------------------------
  character*80 :: Line
      read(11,'(a)') Line
          write(*,*)TRIM(Line)
      read(11,*) HEADER
          write(*,*) TRIM(HEADER)
      read(11,*) AEXPN0     !'Initial Expansion Parameter'
      write(*,*) ' Ainit=',AEXPN0
      
      read(11,*) ASTEP0   !'Initial Step in dAEXPN    da/a = ',da/AEXPN
      read(11,*) sigma     !'DRho/rho in box    '
      read(11,*) Box       !'Box in  Mpc/h   '
      read(11,*) sigma8    !'sigma8    '
      read(11,*) hubble    !'Hubble    '
      read(11,*) Om        !'Omega Matter'
      read(11,*) OmL       !'Omega Lambda'
      read(11,*) Omb       !'Omega Baryons or Neutrinos'
      read(11,*) NROW      !'NROW  Number Particles   '
      read(11,*) NGRID     !'NGRID Number grid points '
      read(11,*) Nseed     !'Random seed'
      read(11,*) Cell      !'Cell Size   '
      read(11,*) aMass     ! 'Particle Mass'   
      read(11,*) zinit     ! 'Initial redshift       '
      read(11,*) zfinal    ! 'Final redshift       '
      read(11,*) DensThr   ! 'Density Threshold for V correction '
      read(11,*) sigV      ! 'rms V correction factor'
      read(11,*) Nout      ! 'Number of redshifts for analysis'
      read(11,*)(zout(i),i=1,Nout)
      read(11,*) Nbiaspars !  'Number of bias parameters'
      Do i=1,Nbiaspars
         read(11,*) BiasPars(i)   ! 'Bias parameter'
      End Do
      read(11,*,iostat=ierr) Ncheckpoint  !'Steps between checkpoints'
      if(ierr/=0)Then           !--- add default values
         write (11,60) 20,   'Steps between checkpoints'
         write (11,60) 0,    'Save snapshots           '
         write (11,60) 1,    'DM power spectrum        '
         write (11,60) 0,    'Redshift distortions     '
         write (11,60) 1,    'DM PDF                   '
         write (11,60) 0,    'Biasing model            '
         write (11,60) 0,    'Write Mocks              '
         write (11,60) 0,    'Dump random fraction     '
         write (11,60) 0,    'Find BDM halos           '
60       format(i5,T20,a)
         write (*,*) ' Setup.dat was uppended. Restart the code.'
         stop
      else
         read(11,*) iSave        !'Save snapshots           '
         read(11,*) iPower       !'DM power spectrum        '
         read(11,*) iPowerRSD    !'Redshift distortions     '
         read(11,*) iDensityDistr !'DM PDF                  '
         read(11,*) iBias        !'Biasing model            '
         read(11,*) iWriteMock   !'Write Mocks              '
         read(11,*) iDumpPart    !'Dump random fraction     '
         read(11,*) iBDM         !'Find BDM halos           '
      end if

      read(11,'(a)') Line
          write(*,*)TRIM(Line)
      read(11,*) MG_flag         !'MG flag                          '
      read(11,*) MG_test         !'MG test flag                     '
      read(11,*) MG_model        !'MG model                         '
      read(11,*) Vcycle          !'V-cycle                          '
      read(11,*) Fcycle          !'F-cycle                          '
      read(11,*) Wcycle          !'W-cycle                          '
      read(11,*) iter_count_max  !'Maximum number of cycles         '
      read(11,*) res_conv        !'Relaxation convergence certerion '

      read(11,'(a)') Line
          write(*,*)TRIM(Line)
      read(11,*) fr_n            !'Hu-Sawicki power-law parameter n '
      read(11,*) fR0             !'Present value of f_R field: fR_0 '

      read(11,'(a)') Line
          write(*,*)TRIM(Line)
      read(11,*) N_branch        !'normal-branch DGP model          '
      read(11,*) S_branch        !'self-accelerate branch-DGP model '
      read(11,*) H0rc            !'DGP parameter: H0rc              '

      read(11,'(a)') Line
          write(*,*)TRIM(Line)
      read(11,*) sym_astar       !'symmetron model a_star           '
      read(11,*) sym_xi          !'symmetron model xi               '
      read(11,*) sym_beta        !'symmetron model beta_star        '

      read(11,'(a)') Line
          write(*,*)TRIM(Line)
      read(11,*) power_law       !'power-law type kmouflage model   '
      read(11,*) born_infeld     !'Born-Infeld type kmouflage model '
      read(11,*) kmf_n           !'power-law kmouflage param n      '
      read(11,*) kmf_K0          !'power-law kmouflage param K0     '
      read(11,*) kmf_beta        !'kmouflage param beta             '

      read(11,'(a)') Line
          write(*,*)TRIM(Line)
      read(11,*) csf_coupling    !'type of coupling function: 1 - exponential; 2 - quadratic '
      read(11,*) csf_potential   !'type of potential: 1 - inverse power law; 2 - SUGRA '
      read(11,*) csf_alpha       !'potential parameter alpha        '
      read(11,*) csf_beta        !'coupling function parameter beta '

      write (*,*) ' Results were read from Setup.dat'
      CLOSE (11)
      AMPLT = sigma
    end subroutine ReadSetup

!---------------------------------------
!        Read    PMfiles
!             moment <0    use PMcrd.DAT, PMcrs0.DAT ...    
!             moment >= 0  use PMcrd.xxxx.DAT, PMcrs0,XXXX.DAT ..
    SUBROUTINE ReadDataPM(moment,Path)
!      
!---------------------------------------
       Character*120 :: Name
        Logical      :: exst
        Integer*8    :: iCount,ii,ioff,ip
        Integer*8    :: Ngal,Nrecord,Jpage,NinPage
        Integer*4    :: moment
        character(len=*) :: Path 
!
        NoVelocities  = 0
        !			Read data and open files
        If(moment<0)Then
           write(Name,'(2a)')Trim(Path),'PMcrd.DAT'
           write(*,*)Trim(Name)
         Open (4,file =Trim(Name),form ='UNFORMATTED',status ='UNKNOWN')
      Else
         write(Name,'(a,i4.4,a)')'PMcrd.',moment,'.DAT'
         Open (4,file =trim(Path)//TRIM(Name),form ='UNFORMATTED',status ='UNKNOWN')
      end If
         
      READ  (4) HEADER,                        &
                       AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW, &
                       TINTG,EKIN,EKIN1,EKIN2,AU0,AEU0,     &
                       NROW,NGRID,Nrealization,Nseed,Om,OmL,hubble, &
                       Nparticles,extras
      WRITE (*,'(a,/10x,a,f8.4,4(a,i7))') HEADER,                 & 
                       ' a=',AEXPN, ' step= ',ISTEP,         &
                       ' Nrow= ', NROW, ' Ngrid=',NGRID
      WRITE (17,'(a,/10x,a,f8.4,4(a,i7))') HEADER,                 & 
                       ' a=',AEXPN, ' step= ',ISTEP,         &
                       ' Nrow= ', NROW, ' Ngrid=',NGRID
      close(4)

      Nrecord = 1024**2
      Naccess = Nrecord*6 !*4
      xR      = NGRID +1
      boxsize = extras(100)
      Box     = boxsize
      Npages   = (Nparticles-1)/Nrecord+1 ! number of records
      Nlast   = Nparticles - (Npages-1)*Nrecord ! number of particles in last record
      Jpage   = Nrecord

      write(17,'(a,i10)') ' NROW   =',NROW
      write(17,'(a,i10)') ' Ngal   =',Nparticles
      write(17,'(a,i10)') ' Ngrid  =',Ngrid
      write(17,'(a,f10.1)') ' Box    =',Box
      Allocate (Xb(Nrecord),Yb(Nrecord),Zb(Nrecord))
      Allocate (VXb(Nrecord),VYb(Nrecord),VZb(Nrecord))

      myMemory =Memory(3_8*Nparticles)
      Allocate (XPAR(Nparticles),YPAR(Nparticles),ZPAR(Nparticles))
      If(NoVelocities == 0)Then   ! --- take velocities
        Allocate (VX(Nparticles),VY(Nparticles),VZ(Nparticles))
        myMemory =Memory(3_8*Nparticles)
     End If

      iCount = 0

      ifile = 0
      jj    = 0
      If(moment<0)Then
         write(Name,'(2a,i1.1,a)')Trim(Path),'PMcrs',ifile,'.DAT'
      Else
         write(Name,'(2a,i1.1,a,i4.4,a)')Trim(Path),'PMcrs',ifile,'.',moment,'.DAT'
      End If
        INQUIRE(file=TRIM(Name),EXIST=exst)
        If(.not.exst)Stop ' File PMcrs0.DAT does not exist. Error'
        OPEN(UNIT=20,FILE=TRIM(Name),ACCESS='DIRECT', &
              FORM='unformatted',STATUS='UNKNOWN',RECL=NACCESS)
   Do ii =1,Npages
      If(ii==Npages)Then
         NinPage = Nparticles -(ii-1)*JPAGE  ! # particles in the current page
      Else
         NinPage = JPAGE
      EndIf
      jj = jj +1
      If(ii<10.or.ii==Npages)write(17,'(3(a,i9))') ' Reading page= ',ii,' record =',jj,' NinPage= ',NinPage
10      Read(20,REC=jj,iostat=ierr) Xb,Yb,Zb,VXb,VYb,VZb
      If(ierr /=0)Then
         close(20)
         ifile = ifile +1
         If(Moment<0)Then
           If(ifile<10)Then
             write(Name,'(2a,i1.1,a,i4.4,a)')Trim(Path),'PMcrs',ifile,'.DAT'
           Else
             write(Name,'(2a,i2.2,a,i4.4,a)')Trim(Path),'PMcrs',ifile,'.DAT'
          EndIf
       Else
          If(ifile<10)Then
             write(Name,'(2a,i1.1,a,i4.4,a)')Trim(Path),'PMcrs',ifile,'.',moment,'.DAT'
          Else
             write(Name,'(2a,i2.2,a,i4.4,a)')Trim(Path),'PMcrs',ifile,'.',moment,'.DAT'
          EndIf
       end If

           jj = 1
           INQUIRE(file=TRIM(Name),EXIST=exst)
           If(.not.exst)Then
             write(*,'(a,3i5)')' Attempting to read file number, record = ',ifile,ii,Npages
             write(*,'(2a)')' Attempting to read non-existing file: ',TRIM(Name)
             Stop ' Error reading PMcrs files: Did not get all files'
           End If
           Open(20,file=TRIM(Name),ACCESS='DIRECT', &
              FORM='unformatted',STATUS='UNKNOWN',RECL=NACCESS)
          write(17,'(2i7,2a,3x,i9)') ii,ifile,' Open file = ',TRIM(Name),Ninpage
          go to 10
       end If

       ioff = (ii-1)*JPAGE
!$OMP PARALLEL DO DEFAULT(SHARED)  PRIVATE (ip)       
           Do ip =1,NinPage
              If(ip+ioff > Nparticles)STOP 'Attempt to read too many particles '
              if(INT(Xb(ip))==Ngrid+1)Xb(ip)=Xb(ip)-1.e-3
              if(INT(Yb(ip))==Ngrid+1)Yb(ip)=Yb(ip)-1.e-3
              if(INT(Zb(ip))==Ngrid+1)Zb(ip)=Zb(ip)-1.e-3
              if(INT(Xb(ip))==Ngrid+1)write(*,*)'Error in boundary: ',INT(Xb(ip)),Xb(ip)
              if(INT(Yb(ip))==Ngrid+1)write(*,*)'Error in boundary: ',INT(Yb(ip)),Yb(ip)
              if(INT(Zb(ip))==Ngrid+1)write(*,*)'Error in boundary: ',INT(Zb(ip)),Zb(ip)
              
                Xpar(ip+ioff) = Xb(ip) 
                Ypar(ip+ioff) = Yb(ip) 
                Zpar(ip+ioff) = Zb(ip)
                If(NoVelocities == 0)Then                 
                  VX(ip+ioff) = VXb(ip) 
                  VY(ip+ioff) = VYb(ip) 
                  VZ(ip+ioff) = VZb(ip)
               end If
           end Do
        end DO
        close (20)

        Np = Nparticles
        xm = 1.e8 ; ym = 1.e8; zm = 1.e8
!$OMP PARALLEL DO DEFAULT(SHARED)  PRIVATE (ip) Reduction(max:xxm,yym,zzm) &
!$OMP  Reduction(min:xm,ym,zm)
        do ip =1,Nparticles
           xxm =max(xxm,Xpar(ip))
           yym =max(yym,Ypar(ip))
           zzm =max(zzm,Zpar(ip))
           xm =min(xm,Xpar(ip))
           ym =min(ym,Ypar(ip))
           zm =min(zm,Zpar(ip))
           
              If(Xpar(ip).lt.1.0.or.Ypar(ip).lt.1.0.or.Zpar(ip).lt.1.0)&
                write(17,'(a,i10,1p,3g14.5)')' Error coord: ',ip,Xpar(ip),Ypar(ip),Zpar(ip)
        EndDo
     write(17,'(a,3(2es14.5,3x))') ' Coordinates Min/Max =',xm,xxm,ym,yym,zm,zzm      
     DEALLOCATE (Xb,Yb,Zb,VXb,VYb,VZb)
     
   end SUBROUTINE ReadDataPM
!---------------------------------------
!        Write     PMfiles: PMcrd/crs   
 SUBROUTINE WriteDataPM(iFlag,Path)
!
!---------------------------------------
       Character*120 :: Name
        Logical      :: exst
        Integer*8    :: iCount,ii,Jpage,Nrecord
        Integer*8    :: Ngal,moment,i,jfirst,jlast,j0
        character(len=*) :: Path 
     Call TimingMain(4,-1)

      Npage   = 1024**2       ! number of particles per record
      Nrecord = Npage
      Jpage   = Npage
      Naccess = Nrecord*6 !!*4 !!!
      Nrecpage = 512 !256            ! max number of records per file
      moment = ISTEP
      Npages  = (Nparticles-1)/JPAGE+1
      Nfiles  = (Npages-1)/Nrecpage +1    ! number of files

      write(*,*) ' Npages      = ',Npages
      write(*,*) ' Write files = ',iFlag,moment
        !			open files
     If(iFlag == 0)Then
        Open (4,file =TRIM(Path)//'PMcrd.DAT',form ='UNFORMATTED',status ='UNKNOWN')
          write(Name,'(a,i1,a)')'PMcrs',0,'.DAT'
          OPEN(UNIT=20,FILE=TRIM(Path)//TRIM(Name),ACCESS='DIRECT', &
                 FORM='unformatted',STATUS='UNKNOWN',RECL=NACCESS)
     Else
        write(Name,'(a,i4.4,a)')'PMcrd.',moment,'.DAT'
        Open (4,file =TRIM(Path)//TRIM(Name),form ='UNFORMATTED',status ='UNKNOWN')
          write(Name,'(a,i1,a,i4.4,a)')'PMcrs',0,'.',moment,'.DAT'
          write(*,'(a)') TRIM(Name)
        OPEN(UNIT=20,FILE=TRIM(Path)//TRIM(Name),ACCESS='DIRECT', &
                 FORM='unformatted',STATUS='UNKNOWN',RECL=NACCESS)
     end If
     PARTW = float(NGRID)**3/FLOAT(Nparticles)
     AU0   = 0.
      write  (4) HEADER,                        &
                       AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW, &
                       TINTG,EKIN,EKIN1,EKIN2,AU0,AEU0,     &
                       NROW,NGRID,Nrealization,Nseed,Om,OmL,hubble, &
                       Nparticles,extras
      close(4)
      myMemory =Memory(6_8*Nrecord)
      Allocate (Xb(Nrecord),Yb(Nrecord),Zb(Nrecord))
      Allocate (VXb(Nrecord),VYb(Nrecord),VZb(Nrecord))
      
        jj = 1
Do i=1,Npages         !-------- dump into  files
   If(i==Npages)Then
      NinPage = Nparticles -(i-1)*JPAGE  ! # particles in the current page
   Else
      NinPage = JPAGE
   EndIf
    If(mod(i-1,Nrecpage)==0.and.i>1)Then   ! close old and open new file
       close(20)
       jj = 1
       ifile = (i-1)/Nrecpage       ! construct file name
       If(iFlag==0)Then
          If(ifile<10)Then
             write(Name,'(a,i1.1,a)')'PMcrs',ifile,'.DAT'
          Else
             write(Name,'(a,i2.2,a)')'PMcrs',ifile,'.DAT'
          EndIf
       Else
          If(ifile<10)Then
             write(Name,'(a,i1.1,a,i4.4,a)')'PMcrs',ifile,'.',moment,'.DAT'
          Else
             write(Name,'(a,i2.2,a,i4.4,a)')'PMcrs',ifile,'.',moment,'.DAT'
          EndIf
       end If
       Open(20,file=TRIM(Path)//TRIM(Name),ACCESS='DIRECT', &
                 FORM='unformatted',STATUS='UNKNOWN',RECL=NACCESS)
       write(*,'(2i7,2a,3x,i9)') i,ifile,' Open file = ',TRIM(Name),Ninpage
    end If

   jfirst  = (i-1)*JPAGE +1        ! first and last particles in current record
   jlast   = jfirst + NinPage-1
   If(mod(i,10)==0.or.i==Npages)    &
        write(*,'(10x,a,i5,a,4i11)')'Write page=',i,' Particles=',NinPage,jfirst,jlast

   Do j0 = jfirst,jlast
      Xb(j0-jfirst+1)   = Xpar(j0)
      Yb(j0-jfirst+1)   = Ypar(j0)
      Zb(j0-jfirst+1)   = Zpar(j0)
      Vxb(j0-jfirst+1)  = VX(j0)
      VYb(j0-jfirst+1)  = VY(j0)
      VZb(j0-jfirst+1)  = VZ(j0)
   EndDo

   !write(*,'(i10,1p,6g13.5)') (k,XPAR(k),YPAR(k),ZPAR(k),VX(k),VY(k),VZ(k),k=1,1024)
   WRITE (20,REC=jj) Xb,Yb,Zb,VXb,VYb,VZb
   jj = jj +1
EndDo                            ! end lspecies loop   
        
     myMemory =Memory(-6_8*Nrecord)
     DEALLOCATE (Xb,Yb,Zb,VXb,VYb,VZb)
     Call TimingMain(4,1)
     
   end SUBROUTINE WriteDataPM   
!
!----------------------------------------------------
function seconds ()
!----------------------------------------------------
!
!     purpose: returns elapsed time in seconds
      Integer*8, SAVE :: first=0,rate=0,i0=0
      Integer*8       :: i

      If(first==0)Then
         CALL SYSTEM_CLOCK(i,rate)
         first =1
         i0    = i
         seconds = 0.
      Else
         CALL SYSTEM_CLOCK(i)
         seconds = float(i-i0)/float(rate)
      EndIf

    end function seconds
!--------------------------------------------
subroutine Timing ( ielement , isign )
!--------------------------------------------
        ! use Timing(i,-1) to start clock
        !     Timing(i,1)  to suspend clock
        !     Timing(0,0)  to print results
    Character*80 :: FName='timing.log'
    Integer OMP_GET_THREAD_NUM, OMP_GET_MAX_THREADS
    If(isign == 0)Then
       Open(30,file=TRIM(FName),position='append')
       write(30,'(i5,1p,10G13.4)')CPU(0:10)/60.
       close(30)
       CPU(:) = 0
    EndIf
     CPU(ielement) = CPU(ielement) + float(isign) * seconds()
   end subroutine Timing
!--------------------------------------------
subroutine TimingMain ( ielement , isign )
!--------------------------------------------
        ! 0 - total
        ! 1 - force, 2- move
        ! 3 - density, 4- IO, 5- biasing, 6- BDM, 7-Analysis
    Character*80 :: FName='timing.log'
    Integer OMP_GET_THREAD_NUM, OMP_GET_MAX_THREADS
    If(isign == 0)Then
       Open(30,file=TRIM(FName),position='append')
       If(iStep==1)Then
          write(30,'(3(a,i11))') ' Ngrid    = ',Ngrid, &
                                ' Npart    = ',Nparticles, &
                                ' Nthreads = ',OMP_GET_MAX_THREADS()
          write(30,'(T3,a,T10,a,T23,a,T35,a,T49,a,T62,a,T73,a,T84,a,T95,a)')&
               'Step','Tot/min','Force','Move','Density','IO','Bias','BDM','Analysis'
       End If
       write(30,'(i5,1p,12G13.4)') iStep,CPU(0:7)/60.
       close(30)
       CPU(:) = 0
       return
    EndIf
     CPU(ielement) = CPU(ielement) + float(isign) * seconds()
   end subroutine TimingMain
!--------------------------------------------
Real*4 Function Memory(i_add)
!--------------------------------------------

    Real*8, SAVE :: mem =0.001   ! initial memory in Gb
    Integer*8 :: i_add           ! number of 4byte words
    Real*8    :: add
    mem    = mem + i_add*4./1024.**3
    Memory = mem
    write(17,'(a,f8.3,a)') ' Current Memory = ',Memory,'Gb'
  end Function Memory

end Module Tools
