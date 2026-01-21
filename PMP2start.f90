!
! _______________________ START 3-D SIMULATIONS                         
!                                                                       
!                         Klypin, August 2015                         
!
!   Uses FFT5 :  Nrow = 4*3**n*5**m
!

!-------------------------------------------------------------
Module setInitialConditions
  use Tools
  
    Integer*4, parameter  :: NtabM = 100000
    REAL*8,    PARAMETER  :: PI=3.1415926535
      Real*4   :: xkt(0:NtabM),Pkt(0:NtabM)  ! Power spectrum
      Real*4   :: StepK,alog0
      Integer*4:: Ntab,iFlip

      INTEGER*4,parameter ::   NPAGE = 1024**2  ! # particles in a record
      INTEGER*4 ::   NMAX,  &  !  = NGRID/2,     & 
                     NRECL, &  !  = NPAGE*6,   & ! # particles in a record
                     NSPEC     ! = NROW/2             ! No waves shorter than Ny
      REAL*4                ::   alpha  =0., Qscale
      Real*8 :: SKINE
      Integer*4, PARAMETER :: nbyteword = 4      ! defines length of direct-access:1or4
      Integer*4, PARAMETER :: LevelParticles = 0 


      REAL,   ALLOCATABLE, DIMENSION(:,:,:) ::  GRX,GRY,GRZ
      Real*4, dimension(NPAGE) ::  XPp,YPp,ZPp, &
     		                   VPx,VPy,VPz

!$OMP THREADPRIVATE(XPp,YPp,ZPp,VPx,VPy,VPz)

    Contains
!--------------------------------------------------                     
!                              sqrt(Power spectrum)
!                                       k = (2pi/L) 
      FUNCTION TRUNF (WK) 
!-------------------------------------------------
        real*4 :: k
      IF (WK.GE.FLOAT (NSPEC) ) THEN 
         TRUNF = 0. 
         RETURN 
      ENDIF 
          k = QSCALE * wk 
          TRUNF = sqrt (Ppk (k) )
        END FUNCTION TRUNF

!---------------------------------------                                
!            interpolate table with p(k)                       
!                                             
        FUNCTION Ppk (x) 
! 
!---------------------------------------                                                        
                              ! slope is ns =-3                         
      If (x.ge.xkt (Ntab) ) Then 
         Ppk = Pkt (Ntab) / (x / xkt (Ntab) ) **3 
         Return 
      EndIf 
                                         
      If (x.lt.xkt (1) ) Then    ! slope is ns=1                
         Ppk = Pkt (1) * (x / xkt (1) ) 
         Return 
      EndIf 
      ind = INT ( (log10 (x) - alog0) / StepK) + 1 
      dk = xkt (ind+1) - xkt (ind) 
      Ppk = (Pkt (ind) * (xkt (ind+1) - x) + Pkt (ind+1) * (x - xkt (   &
                   ind) ) ) / dk                                                     

      END FUNCTION Ppk    

!--------------------------------------------------                     
!               Store  seeds for parallelization 
!                      of random number generator
!
      SUBROUTINE SetRandomN
!                   
!-------------------------------------------------                                                     
Use Random
INCLUDE 'luxuryp.h' 

      gSet  = 0.
      Ns    = Nseed   
      lux   = 2
      Do k=1,NROW                               ! luxury
          CALL rluxgo (lux, Ns, 0, 0)      ! initialize luxury 
            !if(k/10*10==k)write(*,*) ' page number=',k
          i24A(k)    = i24
          j24A(k)    = j24
          in24A(k)  = in24
          kountA(k) = kount
          carryA(k) = carry
          gsetA(k)  = gSet
          iFlagA(k) = iFlag
             Do m=1,24
                 SeedsPage(m,k) = seeds(m)
             EndDo    
              dummy = RANDd(Ns) 
        end Do
End SUBROUTINE SetRandomN
                                                   
!-------------------------------------------------
!        
SUBROUTINE ReadPkTable
    !-- Read PkTable. Assign Omegas
      OmbT = ParseLine(10)
      OmcT = ParseLine(10)
      OmLT = ParseLine(10)
      OmT  = ParseLine(10)
      sigma8T = ParseLine(10)
           !-- check parameters
      If(abs(OmbT-Omb)>1.e-4*OmbT)Stop 'Error: OmegaBaryons is different from PkTable value'
      If(abs(OmT-Om)>1.e-4*OmT)Stop 'Error: OmegaMatter is different from PkTable value'
      write(*,'(3(a,ES12.3))')'  Ombar =', &
                       OmbT,' Om_matter =',OmT
      Ntab = 0
12   READ (10, *, end = 32, err = 32) xx, pp    !---- read P(k)
        Ntab = Ntab + 1 
        xkt (Ntab) = xx  !*hubble 
        Pkt (Ntab) = pp                !  Pk 
      GOTO 12 
32    close(10)
      Write (*,*) ' Read ', Ntab, ' lines from P(k) table '
      
      If (Ntab.le.1) stop 'wrong table for p(k)' 
      StepK = log10 (xkt (Ntab) / xkt (1) )/(Ntab-1) 
      alog0 = log10 (xkt (1) ) 
                                                                        
                            ! test that spacing is the same             
      Do k = 2, Ntab - 1 
        ss = log10 (xkt (k + 1) / xkt (k) ) 
        If (abs (ss / StepK - 1.) .gt.2.e-2) Then 
        Write (*,*) ' error in K spacing. k=', k, xkt (k + 1),xkt (k)
        STOP 
        EndIf 
     EndDo      
   END SUBROUTINE ReadPkTable
  
!-------------------------------------------------------------
!                                                                       
SUBROUTINE Initialize
!     
!-------------------------------------------------------------
                                                
   write(*,'(a,$)') ' Enter realization number for this run = '
   Read(*,*) Nrealization
   write(*,*)
   open(1,file='../TableSeeds.dat')
   read(1,*)     ! skip the first line == header
   Do ij=1,Nrealization
      read(1,*) Nseed0,Ncount
   endDo
   Nseed = Nseed0
   close(1)

   extras (:)  = 0.
   extras(100) = Box
   NMAX        = NGRID/2
   NRECL       = NPAGE*6
   NSPEC       = NROW/2
   NBYTE = NPAGE * 6 * 4 
   Nparticles  = INT(NROW,8)**3 
   ISTEP = 0                                                        
   TINTG = 0.                                                       
   AU0   = 0.                                                       
   AEU0  = 0.                                                       
   EKIN  = 0.                                                       
   EKIN1 = 0.                                                       
   EKIN2 = 0.                                                       
   QSCALE = 2.*PI/Box
      Wtotal = (float(NROW)**3*4*9       &
                + 31.*NROW*4             &
                + 6. *NPAGE*4)/1024.**3
      write (*,'(//10x,a,g13.4,a)') ' Memory required to run the code= ',Wtotal,'Gb'
      Wtotal = (float(NROW)**3*4*6) /1024.**3
      write (*,'(10x,a,g13.4,a)')   ' Disk space required            = ',Wtotal,'Gb'
      
   myMemory= Memory(9_8*Nparticles)
   Allocate(XPAR(Nparticles),VX(Nparticles))
   Allocate(YPAR(Nparticles),VY(Nparticles))
   Allocate(ZPAR(Nparticles),VZ(Nparticles))
   ALLOCATE(GRX(NROW,NROW,NROW))     ! allocate memory for spectrum
   ALLOCATE(GRY(NROW,NROW,NROW))     ! 
   ALLOCATE(GRZ(NROW,NROW,NROW))     ! 
   
   write(*,*) ' NROW   =',NROW
   write(*,*) ' NGRID  =',NGRID
   write(*,*) ' Npart  =',Nparticles
   ASTEP = ASTEP0
   AEXPN = AEXPN0
END SUBROUTINE Initialize

!---------------------------------------------------
!
!             check if all init files are present
!
   Subroutine CheckInit
!---------------------------------------------------
     logical :: exst
      Inquire(file='../PkTable.dat',exist = exst)
      if(.not.exst)Stop ' File PkTable with the power spectrum not found'
      open(10,file='../PkTable.dat')
      
      Inquire(file='../Setup.dat',exist = exst)
      if(.not.exst)Then
         write(*,*)' Error: File ../Setup.dat not found. Run PMP2init.exe'
         stop
      end if
         open(11,file='../Setup.dat')

      Inquire(file='../TableSeeds.dat',exist = exst)
      if(.not.exst)Call SetSeeds
 
end Subroutine CheckInit
    
!------------------------------
!
!             Generate a table with random seeds
!
!------------------------------
Subroutine SetSeeds
  use Random
  integer*8 :: is,ij,iM
  integer*4 :: Nseed0
  is = 1232_8**3
  
    Nseed0  = 1298302
    nslip = 137
    noff  = 2357
    Ntable =5000
    NCount =0
    open(1,file='TableSeeds.dat')
    write(1,*) 'Seeds:',Nseed0,Nslip
   Do ij=1,is
      x  =RANDd(Nseed0)
      Nn =INT(x*nslip)+1
      Do jj =1,noff+Nn
         x  =RANDd(Nseed0)
      End Do
      Ncount =Ncount +1
      write(1,*) Nseed0,Ncount
      If(Ncount>Ntable)exit
   endDo
   close(1)
 end Subroutine SetSeeds
!--------------------------------------------------
!        read line from  input file iFile
!                real format
      Function ParseLine(iFile)
      Character :: Line*120,Line2*120,Line3(120)

      Read(iFile,'(a)')Line
      Ieq =INDEX(Line,'=',BACK=.TRUE.)
               !write(*,*) '  Ieq =',Ieq
      backspace (iFile)                  !--- go to line start
      write(Line2,'(a1,i2,a)') '(',Ieq,'a1,g12.5)' ! make format
           !write(*,'(a)') Line2
      Read(iFile,Line2)(Line3(i),i=1,Ieq),dummy    ! read
      ParseLine = dummy
           !write(*,'(a,ES12.3)') ' Result =',ParseLine
      end Function ParseLine            
!------------------------------------------------                       
!                             Make a realization of spectrum of perturbations
! 
! 
SUBROUTINE SPECTR  
!------------------------------------------------                       
use Random
INCLUDE 'luxuryp.h' 
      REAL(8)   :: SUMM=0., Wi3, Wj3, Wk3, WD, TS, TRX,ss 
      REAL      ::    t0=0,t1=0.        ! counters for timing
      INTEGER*4, allocatable, DIMENSION(:) :: mapz, map3

      CALL Timing(3,-1)                  ! initialize time counter
      allocate(mapz(NROW), map3(NROW))

      SUMM = 0.
      map3(1) = 0
      mapz(1) = 1
      map3(NROW) = 0
      mapz(NROW) = NROW
      DO j=2,NROW-1,2
          j3 = j/2
          map3(j)   = -j3   ! -k   for cos
          map3(j+1) =  j3   !  k       sin
          mapz(j)   =  j+1  ! flip cos <--> sin
          mapz(j+1) =  j
      end DO
                                                           
!$OMP PARALLEL DO DEFAULT(SHARED)  &                                     
!$OMP PRIVATE(Mk3,Mk2,Mk1)                                            
      DO Mk3 = 1, NROW 
      DO Mk2 = 1, NROW 
      DO Mk1 = 1, NROW 
         GRX(MK1,MK2,MK3) =0.
	 GRY(MK1,MK2,MK3) =0.
	 GRZ(MK1,MK2,MK3) =0.
      ENDDO 
      ENDDO 
      ENDDO 
!                                      Set Spectrum                     
!$OMP PARALLEL DO DEFAULT(SHARED)  &                                     
!$OMP PRIVATE(k,ksign,kz,k3,Wk3,j,jsign,jz,j3,Wj3,i,isign,iz,i3,Wi3)  &
!$OMP PRIVATE(WD,Wk,TS,TRX,m,NRAND,gSet,iFlag)           &
!$OMP REDUCTION(+:SUMM)
      DO k = 1, NROW 
           i24      = i24A(k) 
           j24      = j24A(k) 
           in24    = in24A(k) 
           kount   = kountA(k)
           mkount  = 0 
           carry  = carryA(k)
           iflag  = iFlagA(k)
           gset   = gSetA(k)
              Do m=1,24
                  seeds(m) =SeedsPage(m,k) 
              EndDo    
         Wk3 = map3(k)**2 

         DO j = 1, NROW 
            Wj3 = map3(j)**2 
            DO i = 1, NROW 
               Wi3 = map3(i)**2 
                TS = GAUSS3 (gSet,iFlag) 
               WD  = Wi3 + Wj3 + Wk3 

               IF (WD< 1.e-3) THEN 
		  GRX(1,1,1) = 0.
		  GRY(1,1,1) = 0.
		  GRZ(1,1,1) = 0.
               ELSE 
                  Wk = SQRT (WD) 
                  TS = TRUNF (Wk) * TS
                  TRX = TS / WD 

		  GRX(mapz(i), j, k) =  TRX *map3(i)
		  GRY( i,mapz(j), k) =  TRX *map3(j)
		  GRZ( i, j,mapz(k)) =  TRX *map3(k)                  
                  SUMM = SUMM + TS**2 
               ENDIF        ! end kx=ky=kz=0 switch
            ENDDO         ! i
         ENDDO            ! j
      ENDDO               ! k
      
      IF (SUMM.LE.0.) Write (*,  * ) ' Error!!! Summ over spectrum = 0'
      ALPHA = AMPLT / SQRT (SUMM) * sqrt (8.) 
      Write (*,'(10x,a20,3g13.6)') 'SPECTR ==>  SUMM=', SUMM, ALPHA 
      
      Do i =1,Nrow             ! adjust zero-k values of FI
         Do j=1,NROW
            GRX(i,j,1) = GRX(i,j,1)/sqrt(2.)
            GRX(i,1,j) = GRX(i,1,j)/sqrt(2.)
            GRX(1,i,j) = GRX(1,i,j)/sqrt(2.)
            GRY(i,j,1) = GRY(i,j,1)/sqrt(2.)
            GRY(i,1,j) = GRY(i,1,j)/sqrt(2.)
            GRY(1,i,j) = GRY(1,i,j)/sqrt(2.)
            GRZ(i,j,1) = GRZ(i,j,1)/sqrt(2.)
            GRZ(i,1,j) = GRZ(i,1,j)/sqrt(2.)
            GRZ(1,i,j) = GRZ(1,i,j)/sqrt(2.)
         End Do
      End Do

      deallocate(mapz, map3)
      CALL Timing(3,1)                  

      END SUBROUTINE SPECTR                         
!------------------------------------------------                       
!		                        Define coordinates and velocities for        
!                                  particles with resolution = Indx     
!                                  Indx =1 = high resolution            
!                                  Indx =2 = medium resolution          
!                                  Icurrent = number of particles       
      SUBROUTINE BLOCKS (XCONS, VCONS) 
!------------------------------------------------                       
      Real(8)       :: sDispl =0.
      REAL          :: t0=0,t1=0.        ! counters for timing
      Integer*8     :: Icurrent,JROW
                                 
      CALL Timing(4,-1)                    ! initialize time counter   

      sDispl = 0.
      SKINE  = 0. 
      xShift =  0.5 
      JROW   = NROW  ! make it long integer
      XMAX = FLOAT (NGRID) + 1.
      XSHF = FLOAT (NGRID) 

      write(*,*) ' XCONS, VCONS =',XCONS, VCONS
      
      QFACT = FLOAT (NGRID) / FLOAT (NROW)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ip,jp,kp,Icurrent,Q1,Q2,Q3) &
!$OMP REDUCTION(+:sDispl,SKINE)      
      Do kp = 1, NROW
         Q3 = QFACT * (kp-1) +1.
      Do jp = 1, NROW
         Q2 = QFACT * (jp-1) +1.
      Do ip = 1, NROW
         Q1 = QFACT * (ip-1) +1.
            
         Icurrent = ip+ (jp-1_8)*JROW +(kp-1_8)*JROW*JROW 
         
         XPAR(Icurrent) = Q1- XCONS * GRX(ip,jp,kp) + xShift
         VX(Icurrent)   =     VCONS * GRX(ip,jp,kp)
         YPAR(Icurrent) = Q2- XCONS * GRY(ip,jp,kp) + xShift
         VY(Icurrent)   =     VCONS * GRY(ip,jp,kp)
         ZPAR(Icurrent) = Q3- XCONS * GRZ(ip,jp,kp) + xShift
         VZ(Icurrent)   =     VCONS * GRZ(ip,jp,kp)

         IF (XPAR(Icurrent) .GT.XMAX) XPAR(Icurrent) = XPAR(Icurrent) - XSHF 
         IF (XPAR(Icurrent) .LE.1.)   XPAR(Icurrent) = XPAR(Icurrent) + XSHF
         IF (XPAR(Icurrent) .GE.XMAX) XPAR(Icurrent) = XPAR(Icurrent) -1.e-3
         IF (XPAR(Icurrent) .GE.XMAX) XPAR(Icurrent) = XPAR(Icurrent) -1.e-3
         
         IF (YPAR(Icurrent) .GT.XMAX) YPAR(Icurrent) = YPAR(Icurrent) - XSHF 
         IF (YPAR(Icurrent) .LE.1.)   YPAR(Icurrent) = YPAR(Icurrent) + XSHF 
         IF (YPAR(Icurrent) .GE.XMAX) YPAR(Icurrent) = YPAR(Icurrent) -1.e-3
         IF (YPAR(Icurrent) .GE.XMAX) YPAR(Icurrent) = YPAR(Icurrent) -1.e-3
         
         IF (ZPAR(Icurrent) .GT.XMAX) ZPAR(Icurrent) = ZPAR(Icurrent) - XSHF 
         IF (ZPAR(Icurrent) .LE.1.)   ZPAR(Icurrent) = ZPAR(Icurrent) + XSHF 
         IF (ZPAR(Icurrent) .GE.XMAX) ZPAR(Icurrent) = ZPAR(Icurrent) -1.e-3
         IF (ZPAR(Icurrent) .GE.XMAX) ZPAR(Icurrent) = ZPAR(Icurrent) -1.e-3

         sDispl = sDispl + (XPAR(Icurrent)-Q1-xShift)**2+ &
                           (YPAR(Icurrent)-Q2-xShift)**2+ &
                           (ZPAR(Icurrent)-Q3-xShift)**2
         SKINE = SKINE + VX(Icurrent)**2 +VY(Icurrent)**2 +VZ(Icurrent)**2 
      EndDo        ! end ip
      EndDo        ! end jp
      EndDo        ! end kp
                                                                                           
      Write (16, '(10x,a20,g12.4,a,i3)') ' RMS 3d diplacement=',  &
           sqrt (sDispl/max(Icurrent,1))  

      CALL Timing(4,1)                                                                           
      END  SUBROUTINE BLOCKS                                         
!------------------------------------------------                       
!                                            FFT of the spectrum
      SUBROUTINE VECTOR  
!------------------------------------------------                       
use fft5
      integer*4, parameter :: Nlensav = 8192
      integer*4, parameter :: Nlenwrk = 8192
      real*8,    parameter :: P16 = 6.28318530718
      real*4,    parameter :: sq2 = 1. !1.41421356237
      real*8,  save        :: wsave(1:Nlensav)
      real*8,  save        :: work(1:Nlenwrk)
      REAL*8               :: XX,D1,D2,A1,A2,A3,wi,wj,wk
      Integer*4            :: OMP_GET_MAX_THREADS,OMP_GET_THREAD_NUM
      Integer*4            :: Ng,ier,lensav,lenwrk,lenr,inc
      real*8               :: r(Nlenwrk)
      REAL                 :: t0=0,t1=0.        ! counters for timing
      Real*8               :: ss
!$OMP THREADPRIVATE(work,wsave)

      CALL Timing(5,-1)

      If(NROW>Nlenwrk)Stop ' Incresase Nlenwrk in POTENTfft5'
      Ng = NROW
      lensav = NROW+int(log(real(NROW,kind = 4))/log(2.0E+00))+4
      lenwrk = NROW

      call rfft1i ( Ng, wsave, lensav, ier ) !   Initialize FFT
      inc  = 1
      lenr = Ngrid

             write(*,*)  ' XY fft'
!$OMP PARALLEL DO DEFAULT(SHARED)  copyin(wsave,work) & 
!$OMP PRIVATE ( k,j,i ,r,ier)
    Do k=1,NROW             ! fft for xy planes in x-dir
       Do j=1,NROW     
          Do i=1,NROW
             r(i) = GRX(i,j,k)
          EndDo
          call rfft1b ( Ng, inc, r, lenr, wsave, lensav, work, lenwrk, ier )
          Do i=1,NROW
            GRX(i,j,k) = r(i)
         EndDo
       EndDo
       
       Do i=1,NROW        ! fft xy planes in y-dir
          Do j=1,NROW
             r(j) = GRX(i,j,k)
          EndDo
          call rfft1b ( Ng, inc, r, lenr, wsave, lensav, work, lenwrk, ier )
          Do j=1,NROW
            GRX(i,j,k) = r(j)
          EndDo
       EndDo
    EndDo

      write (*,'(10x,a)') 'Swap k<->i'
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE ( k,j,i,aa)
      DO J=1,NROW
      DO K=1,NROW-1
            DO I=K+1,NROW
               aa = GRX(I,J,K)
               GRX(I,J,K) =GRX(K,J,I)
               GRX(K,J,I) =aa
            ENDDO
         ENDDO
      ENDDO

!$OMP PARALLEL DO DEFAULT(SHARED)  copyin(wsave,work) & 
!$OMP PRIVATE ( k,j,i ,r, ier)
     Do j=1,NROW     ! ------ z-direction
       Do i=1,NROW     
           Do k=1,NROW
             r(k) = GRX(k,j,i)
          EndDo
          call rfft1b ( Ng, inc, r, lenr, wsave, lensav, work, lenwrk, ier )
           Do k=1,NROW
             GRX(k,j,i) = r(k) 
          EndDo
       EndDo   
    end Do
      
      write (*,'(10x,a)') 'Swap i<->k'
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE ( k,j,i,aa)
      DO J=1,NROW
      DO K=1,NROW-1
            DO I=K+1,NROW
               aa = GRX(I,J,K)
               GRX(I,J,K) =GRX(K,J,I)
               GRX(K,J,I) =aa
            ENDDO
         ENDDO
      ENDDO


!--------------------- Y ---------------

             write(*,*)  ' XY fft'
!$OMP PARALLEL DO DEFAULT(SHARED)  copyin(wsave,work) & 
!$OMP PRIVATE ( k,j,i ,r,ier)
    Do k=1,NROW             ! fft for xy planes in x-dir
       Do j=1,NROW     
          Do i=1,NROW
             r(i) = GRY(i,j,k)
          EndDo
          call rfft1b ( Ng, inc, r, lenr, wsave, lensav, work, lenwrk, ier )
          Do i=1,NROW
            GRY(i,j,k) = r(i)
         EndDo
       EndDo
       
       Do i=1,NROW        ! fft xy planes in y-dir
          Do j=1,NROW
             r(j) = GRY(i,j,k)
          EndDo
          call rfft1b ( Ng, inc, r, lenr, wsave, lensav, work, lenwrk, ier )
          Do j=1,NROW
            GRY(i,j,k) = r(j)
          EndDo
       EndDo
    EndDo

      write (*,'(10x,a)') 'Swap k<->i'
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE ( k,j,i,aa)
      DO J=1,NROW
      DO K=1,NROW-1
            DO I=K+1,NROW
               aa = GRY(I,J,K)
               GRY(I,J,K) =GRY(K,J,I)
               GRY(K,J,I) =aa
            ENDDO
         ENDDO
      ENDDO

!$OMP PARALLEL DO DEFAULT(SHARED)  copyin(wsave,work) & 
!$OMP PRIVATE ( k,j,i ,r, ier)
     Do j=1,NROW     ! ------ z-direction
       Do i=1,NROW     
           Do k=1,NROW
             r(k) = GRY(k,j,i)
          EndDo
          call rfft1b ( Ng, inc, r, lenr, wsave, lensav, work, lenwrk, ier )
           Do k=1,NROW
             GRY(k,j,i) = r(k) 
          EndDo
       EndDo   
    end Do
      
      write (*,'(10x,a)') 'Swap i<->k'
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE ( k,j,i,aa)
      DO J=1,NROW
      DO K=1,NROW-1
            DO I=K+1,NROW
               aa = GRY(I,J,K)
               GRY(I,J,K) =GRY(K,J,I)
               GRY(K,J,I) =aa
            ENDDO
         ENDDO
      ENDDO
      
!--------------------- Z ---------------

             write(*,*)  ' XY fft'
!$OMP PARALLEL DO DEFAULT(SHARED)  copyin(wsave,work) & 
!$OMP PRIVATE ( k,j,i ,r,ier)
    Do k=1,NROW             ! fft for xy planes in x-dir
       Do j=1,NROW     
          Do i=1,NROW
             r(i) = GRZ(i,j,k)
          EndDo
          call rfft1b ( Ng, inc, r, lenr, wsave, lensav, work, lenwrk, ier )
          Do i=1,NROW
            GRZ(i,j,k) = r(i)
         EndDo
       EndDo
       
       Do i=1,NROW        ! fft xy planes in y-dir
          Do j=1,NROW
             r(j) = GRZ(i,j,k)
          EndDo
          call rfft1b ( Ng, inc, r, lenr, wsave, lensav, work, lenwrk, ier )
          Do j=1,NROW
            GRZ(i,j,k) = r(j)
          EndDo
       EndDo
    EndDo

      write (*,'(10x,a)') 'Swap k<->i'
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE ( k,j,i,aa)
      DO J=1,NROW
      DO K=1,NROW-1
            DO I=K+1,NROW
               aa = GRZ(I,J,K)
               GRZ(I,J,K) =GRZ(K,J,I)
               GRZ(K,J,I) =aa
            ENDDO
         ENDDO
      ENDDO

!$OMP PARALLEL DO DEFAULT(SHARED)  copyin(wsave,work) & 
!$OMP PRIVATE ( k,j,i ,r, ier)
     Do j=1,NROW     ! ------ z-direction
       Do i=1,NROW     
           Do k=1,NROW
             r(k) = GRZ(k,j,i)
          EndDo
          call rfft1b ( Ng, inc, r, lenr, wsave, lensav, work, lenwrk, ier )
           Do k=1,NROW
             GRZ(k,j,i) = r(k) 
          EndDo
       EndDo   
    end Do
      
      write (*,'(10x,a)') 'Swap i<->k'
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE ( k,j,i,aa)
      DO J=1,NROW
      DO K=1,NROW-1
            DO I=K+1,NROW
               aa = GRZ(I,J,K)
               GRZ(I,J,K) =GRZ(K,J,I)
               GRZ(K,J,I) =aa
            ENDDO
         ENDDO
      ENDDO
      
      
      Write (*,'(10x,a)') ' FFT is done' 
      CALL Timing(5,1)

      RETURN 
      END SUBROUTINE VECTOR                         

                                                                        
!---------------------------------------------------                    
!                                  Read  current data from disk/tape,   
!                                  Open files                           
!                                  Nrecl is the number of values in a re
!                                  Npage is the number of particles in a
      SUBROUTINE FilesOpen 

Character (LEN=50) :: FileName
Logical            :: FileExists 


!                                     Open control file               
      OPEN(9, FILE = 'PMcrd.DAT', form = 'unformatted') 
                            ! this clears old header file               
      Write(9) 
      CLOSE (9) 
      OPEN(9, FILE = 'PMcrd.DAT', FORM = 'UNFORMATTED', STATUS =       &
         'UNKNOWN')                                                        
      OPEN(16, FILE = 'Results.log') 
      OPEN(25, file = 'pt.dat', form = 'unformatted') 
                              

      RETURN 
      END SUBROUTINE FilesOpen                      
!---------------------------------------------                          
!                       Write current data to  disk/tape                
!                                                                       
      SUBROUTINE FilesWrite 
!----------------------------------------------                         
!                                       write header and control data   
REAL               ::    t0=0,t1=0.        ! counters for timing
Integer*8          ::    j0,jfirst,jlast,NinPage,JPAGE,Npages
Character*80       ::    fname
      CALL Timing(7,-1)                   ! initialize time counter

      Write (9) HEADER, AEXPN, AEXP0, AMPLT, ASTEP, ISTEP, 0., TINTG,&
      EKIN, EKIN1, EKIN2, AU0, AEU0, NROW, NGRID, Nrealization, Nseed,    &
      Om, OmL, hubble, Nparticles, extras                             
      REWIND 9 

      NBYTE = NPAGE*4
      NACCES= NBYTE / nbyteword
      JPAGE = NPAGE              ! int8
      Nrecpage = 512 ! 256            ! max number of records per file
      
      Npages  = (Nparticles-1)/JPAGE+1
      Nfiles  = (Npages-1)/Nrecpage +1    ! number of files
      write(*,*) ' Npages=',Npages
      write(*,*) ' Nparts=',Nparticles
      If(Npages.eq.0)stop ' Zero number of pages requested. Stop'      
      

      XMAX = FLOAT (NGRID) + 1. 
      XSHF = FLOAT (NGRID) 
      NBYTE = NRECL*4
      NACCES= NBYTE / nbyteword

      Nfiles = (Npages-1)/Nrecpage +1

      Do i=1,Nfiles                       !-- open files
         ifile = (i-1)
         If(ifile<10)Then
            write(Fname,'(a,i1.1,a)')'PMcrs',ifile,'.DAT'
         Else
            write(Fname,'(a,i2.2,a)')'PMcrs',ifile,'.DAT'
         EndIf
         Open(30+ifile,file=TRIM(Fname),ACCESS='DIRECT', &
             FORM='unformatted',STATUS='UNKNOWN',RECL=NACCES)
         write(*,'(2i7,2a,3x,i9)') i,ifile,' Open file = ',TRIM(Fname)
      End Do
         
 Do i=1,Npages        
   If(i==Npages)Then
      NinPage = Nparticles -(i-1)*JPAGE  ! # particles in the current page
   Else
      NinPage = JPAGE
   EndIf
   
   jfirst  = (i-1)*JPAGE +1
   jlast   = jfirst + NinPage-1
   If(mod(i,100)==0.or.i==Npages)    &
        write(*,'(10x,a,i5,a,4i11)')'Write page=',i,' Particles=',NinPage,jfirst,jlast

   Do j0 = jfirst,jlast                !-- move data to the buffers
      XPp(j0-jfirst+1)   = XPAR(j0)
      VPx (j0-jfirst+1)  = VX(j0)
      YPp(j0-jfirst+1)   = YPAR(j0)
      VPy (j0-jfirst+1)  = VY(j0)
      ZPp(j0-jfirst+1)   = ZPAR(j0)
      VPz (j0-jfirst+1)  = VZ(j0)
   EndDo

     WRITE ((i-1)/Nrecpage+30,REC=mod(i-1,Nrecpage)+1) XPp,YPp,ZPp,VPx,VPy,VPz

 EndDo                            ! end Npages                             
      CALL Timing(7,1)

      
      END SUBROUTINE FilesWrite                          
    end MODULE setInitialConditions
                                                                          
!-------------------------------------------------------------------------
PROGRAM  PMstartMp
  use      setInitialConditions
  use      Random
  use      FFT5
  use      Density
  use      Power
   CALL Timing(0,-1)
   
   CALL Timing(1,-1)
   iFlip = 1   ! -1
   iPower= 1   ! make power spectrum of ICs
   iWrite= 1   ! write files on disk
   Call CheckInit     ! test setup, open input files
   Call ReadSetup
   Call Initialize
   Call ReadPkTable
   Call SetRandomN
   CALL Timing(1,1)
!   write(*,'(a,$)')  &
!      ' Enter three integers: (1,-1) for Flip sims, (1,0) for making initial Pk, (1,0) for writing ICs: '
!   Read(*,*) iFlip,iPower,iWrite
   if(iFlip==0)iFlip=1
   iFlip=min(max(iFlip,-1),1)
   iWrite = max(min(iWrite,1),0)
   CALL FilesOpen                   ! this opens files on disk
      
         AEXP0 = AEXPN 
         AEXPV = AEXPN - ASTEP / 2.
         Fact   = sqrt (Om + OmL * AEXPV**3) 
         QFACT  = FLOAT (NGRID) / FLOAT (NROW) 
         Vscale = Box * 100. / NGRID 

              write(*,*) ' Go to Spectrum '
    CALL SPECTR 
          !   get the displacement vector by FFT 
          VCONS = - iFlip* ALPHA/(2.*PI/NGRID)*(AEXPV/AEXP0)*SQRT(AEXPV)*Fact
          XCONS =   iFlip* ALPHA/(2.*PI/NGRID)*(AEXPN/AEXP0) 
            Write(*,'(3x,a12,g12.4,a,g12.4)' ) 'Scaling factors:(x)=', XCONS, ' (v)=', VCONS 

    CALL VECTOR
             Write (*,'(3x,a,ES12.4)' ) 'SPECTR done: Alpha=', ALPHA 
    CALL BLOCKS (XCONS, VCONS) 
      
      DEALLOCATE (GRX,GRY,GRZ)
      PARTW = (FLOAT(NGRID)/FLOAT(NROW))**3
      EKIN = 0.5 * SKINE / AEXPV**2* PARTW
      Write (*, '('' Ekin='',E12.4,'' Particle Weight per cell='',g12.5)') EKIN,PARTW            
      Write (16,'('' Ekin='',E12.4,'' Particle Weight per cell='',g12.5)') EKIN,PARTW


      If(iPower ==1)Then
         Allocate (FI(NGRID,NGRID,NGRID))
         CALL DENSIT
         Call GetPower(0)
         Deallocate(FI)
      end If


     If(iWrite==1)Then
        CALL FilesWrite           ! write header and control data
        Write (25) astep        ! write pt.dat file: time-step for particles             
        close(25)
     End If
      CALL Timing(0,1)
            write(*,'(10x,a/,(10x,a15,f10.2))') &
                     'Wall cpu_time (seconds) ', &
                     'Init        =',CPU(1), &
                     'Spectr      =',CPU(3), &
                     'Coords      =',CPU(4), &
                     'FFT         =',CPU(5), &
                     'Write final =',CPU(7), &
                     'Total       =',CPU(0)

      STOP 
      END  PROGRAM  PMstartMp                                         
