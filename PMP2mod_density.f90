!---------------------------------------------------------
!
Module Density
!
!---------------------------------------------------------
!
  use Tools
  use ExtradofBackgroundData                            ! MG background module, only used for tests
  
Contains

!-------------------------------------------------------------
!
!         Make density contrast FI(NGRID,NGRID,NGRID)
!         for Nparticles XPAR,YPAR,ZPAR(Nparticles) in 1-Ngrid
!
SUBROUTINE DENSIT
    
!-------------------------------------------------------------
  use ExtradofBackgroundData                            ! MG background module, only used for tests

  real*8 :: XN,YN,X,Y,Z,D1,D2,D3,T1,T2,T3,T2W,D2W
  integer*8 :: IN
  integer*4 :: iSplit = 8

  real*8 :: AAA                                         ! amplitude of sine function
  real*8 :: ctilde,R_bg,R0_bg,fR_bg                     ! quantities for test of f(R) gravity
  real*8 :: beta,Orc,KK                                 ! quantities for test of DGP gravity
  real*8 :: XX1,XX2,XX0                                 ! quantities for test of kmouflage model
  real*8 :: phibar,phi,phitot,coeff                     ! quantities for test of coupled scalar field model
  integer:: ist                                         ! quantities for test of coupled scalar field model
  real*8 :: meff
  real*8 :: fct,fct1,fct2
  real*8 :: JJ = 0.02D0, aa = 0.9999D0, ww = 0.15D0     ! quantities for test of DGP gravity
  real*8 :: radius, RR                                  ! quantities for test of DGP gravity

  Call TimingMain(3,-1)
  
  XN   = FLOAT(NGRID)+1.-1.E-7
  YN   = FLOAT(NGRID)
  Wpar = YN**3/FLOAT(Nparticles)

  ! Subtract mean density
!   WRITE(*,*) ' Init density',seconds()
  
! FI(:,:,:) = -1.
  !$OMP PARALLEL DO DEFAULT(SHARED) &
  !$OMP PRIVATE (M1,M2,M3)
  DO M3=1,NGRID
     DO M2=1,NGRID
        DO M1=1,NGRID
           FI(M1,M2,M3) = -1.
        END DO
     END DO
  END DO
  
  WRITE(*,*) ' init finished: ',NGRID,seconds()

  IF(MG_flag.AND.MG_test) THEN
     !
     IF(MG_model.EQ.1) THEN
        ctilde = 2.99792458E3*DBLE(NGRID)/Box
        R_bg   = 3.0d0 * (Om/AEXPN**3 + 4.0d0*OmL)
        R0_bg  = 3.0d0 * (Om          + 4.0d0*OmL)
        fR_bg  = fR0 * (R0_bg / R_bg)**(fr_n + 1)
     ENDIF
     IF(MG_model.EQ.2) THEN
        Orc    = 1.0D0/(4.0D0*H0rc**2)
        IF(N_branch) beta = 1.0D0 + (0.5D0*Om/AEXPN**3+OmL)/(DSQRT(Orc*(Om/AEXPN**3+OmL))) ! normal branch
        IF(S_branch) beta =       - (0.5D0*Om/AEXPN**3+Orc)/(DSQRT(Orc*(Om/AEXPN**3+Orc))) ! self-accelerated branch
     ENDIF
     IF(MG_model.EQ.3.OR.MG_model.EQ.4) THEN
        ctilde = 2.99792458E3*DBLE(NGRID)/Box
     ENDIF
     IF(MG_model.EQ.5) THEN
        ctilde = 2.99792458E3*DBLE(NGRID)/Box
        ist = 1
        DO WHILE(BackgroundEvolution(ist,1).LT.AEXPN)
           ist=ist+1
        ENDDO
        !
        phibar =  BackgroundEvolution(ist-1,2)+ &
               & (BackgroundEvolution(ist,  2)-BackgroundEvolution(ist-1,2))/ &
               & (BackgroundEvolution(ist,  1)-BackgroundEvolution(ist-1,1))* &
               & (AEXPN                       -BackgroundEvolution(ist-1,1))
     ENDIF
     !
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3,radius,AAA,KK,coeff,phi,phitot)
     DO M3=1,NGRID
        DO M2=1,NGRID
           DO M1=1,NGRID
              
            ! homogeneous density
            ! FI(M3,M2,M1) = 0.0D0
            ! homogeneous density end

            ! ----------------
            ! sine field tests
            ! ----------------
            ! IF(MG_model.EQ.1) THEN
            !    AAA = 0.1D0
            !    FI(M3,M2,M1) = AEXPN/Om*ctilde**2*fR_bg*(4.0d0*DACOS(-1.0d0)**2/DBLE(NGRID)**2*AAA)* &
            !                & DSIN(2.0D0*DACOS(-1.0D0)*(DBLE(M1)-0.5D0)/DBLE(NGRID)) + &
            !                & R_bg/3.0D0*AEXPN**3/Om*(1.0/DSQRT(1.0D0+AAA*DSIN(2.0D0*DACOS(-1.0D0)*(DBLE(M1)-0.5D0)/DBLE(NGRID)))-1.0D0)
            ! ENDIF
            ! IF(MG_model.EQ.2) THEN
            !    ! the DGP sine field test seems to require a large number of iterations (~30)
            !    AAA = 0.1D0
            !    KK = 8.0D0
            !    FI(M3,M2,M1) = -4.0D0*KK*2*(DACOS(-1.0D0))**2*beta*AEXPN/Om/DBLE(NGRID)**2*AAA* &
            !                & DSIN(2.0D0*KK*DACOS(-1.0D0)*(DBLE(M1)-0.5D0)/DBLE(NGRID))
            !    fct1 = -AEXPN*beta/(Om*DBLE(NGRID)**2)
            !    FI(M3,M2,M1) = fct1*DSIN(2.0D0*DACOS(-1.0D0)*DBLE(M1)/DBLE(NGRID))
            ! ENDIF
            ! IF(MG_model.EQ.3) THEN
            !    AAA = 0.1D0
            !    KK  = 4.0D0
            !    FI(M3,M2,M1) = (AEXPN/sym_astar)**3*(1.0D0-(1.0D0+AAA*DSIN(2.0D0*DACOS(-1.0D0)*KK*(DBLE(M1)-0.5D0)/DBLE(NGRID)))**2  &
            !                 &                              -(2.0D0*DACOS(-1.0D0)*KK/DBLE(NGRID))**2*2.0D0*(ctilde*sym_xi/AEXPN)**2* &
            !                 &                                   (AAA*DSIN(2.0D0*DACOS(-1.0D0)*KK*(DBLE(M1)-0.5D0)/DBLE(NGRID)))/      &
            !                 &                             (1.0D0+AAA*DSIN(2.0D0*DACOS(-1.0D0)*KK*(DBLE(M1)-0.5D0)/DBLE(NGRID)))       &
            !                 & ) - 1.0D0
            ! ENDIF
            ! IF(MG_model.EQ.4) THEN
            !    IF(power_law) THEN
            !       AAA = 0.001D0
            !       KK  = 32.0D0
            !       XX0 = 2.0D0    *DACOS(-1.0D0)/DBLE(NGRID)*KK
            !       XX1 = 2.0D0*AAA*DACOS(-1.0D0)/DBLE(NGRID)*KK
            !       XX2 = DBLE(kmf_n)*kmf_K0*(0.5D0/AEXPN**2/kmf_lambda**2)**(kmf_n-1)
            !       FI(M3,M2,M1) = 0.0D0 - (AEXPN/(3.0D0*kmf_beta*Om))*ctilde * &
            !                    & (1.0D0+(2.0D0*DBLE(kmf_n)-1.0D0)*XX2*XX1**(2*(kmf_n-1)) * &
            !                    &        (DCOS(XX0*(DBLE(M1)-0.5D0)))**(2*(kmf_n-1)))  * &
            !                    & XX0**2*AAA*DSIN(XX0*(DBLE(M1)-0.5D0))
            !    ENDIF
            !    IF(born_infeld) THEN
            !
            !    ENDIF
            ! ENDIF
            ! IF(MG_model.EQ.5) THEN
            !    AAA = 0.1D0
            !    KK  = 4.0D0
            !    coeff = 1.0D0/(3.0D0*csf_beta*Om/AEXPN)
            !    phi = AAA*DSIN(2.0D0*DACOS(-1.0D0)/DBLE(NGRID)*KK*(DBLE(M1)-0.5D0))
            !    phitot = phibar+phi/ctilde**2
            !    FI(M3,M2,M1) = DEXP(-csf_beta*phi/ctilde**2) + &
            !                 & coeff*DEXP(-csf_beta*phitot)*csf_alpha*csf_lambda**2*AEXPN**2*(1.0D0/phitot**(1.0D0+csf_alpha)-1.0D0/phibar**(1.0D0+csf_alpha)) - &
            !                 & coeff*DEXP(-csf_beta*phitot)*(2.0D0*DACOS(-1.0D0)*KK/DBLE(NGRID))**2*phi - &
            !                 & 1.0D0
            ! ENDIF
            ! --------------------
            ! sine field tests end
            ! --------------------

            ! ----------------------------
            ! Gaussian density field tests
            ! ----------------------------
            ! IF(MG_model.EQ.1) THEN
            !    AAA = 1.0D0
            !    FI(M3,M2,M1) = 1.0D0/3.0D0 * R_bg * AEXPN**2 * (AAA**(-0.5D0) * DEXP((DBLE(M1)/DBLE(NGRID))**2 / 4.0D0) - 1.0D0) - ctilde**2 * fR_bg * AAA * DEXP(-(DBLE(M1)/DBLE(NGRID))**2 / 2.0D0) * ((DBLE(M1)/DBLE(NGRID))**2 - 1.0D0)
            !    FI(M3,M2,M1) = FI(M3,M2,M1) * AEXPN / Om
            ! ENDIF
            ! IF(MG_model.EQ.2) THEN
            !    fct1 = AEXPN*beta/(Om*DBLE(NGRID)**2)
            !    fct2 = 2.0D0 *JJ*aa/ww**2
            !    FI(M3,M2,M1) = fct1*fct2*(1.0D0 - 2.0D0*(DBLE(M1)/DBLE(NGRID)-0.5D0)**2/ww**2)*DEXP(-(DBLE(M1)/DBLE(NGRID)-0.5D0)**2/ww**2)
            ! ENDIF
            ! --------------------------------
            ! Gaussian density field tests end
            ! --------------------------------

            ! --------------------------
            ! Yukawa density field tests
            ! --------------------------
            ! AAA = 1.0D0 
            ! meff = 1.0D0 
            ! FI(M3,M2,M1) = 1.0D0/3.0D0 * R_bg * AEXPN**2 * (AAA**(-0.5) * DEXP(meff * DBLE(M1)/DBLE(NGRID) / 2.0D0) * DSQRT(DBLE(M1)/DBLE(NGRID)) - 1.0D0) - ctilde**2 * fR_bg * meff**2 * AAA * DEXP(-meff * DBLE(M1)/DBLE(NGRID)) / (DBLE(M1)/DBLE(NGRID))
            ! FI(M3,M2,M1) = AEXPN / Om * FI(M3,M2,M1)
            ! ------------------------------
            ! Yukawa density field tests end
            ! ------------------------------

            ! ---------------------------
            ! Spherical overdensity tests
            ! ---------------------------
             radius = DSQRT((DBLE(M1)-0.5D0*DBLE(NGRID))**2 + (DBLE(M2)-0.5D0*DBLE(NGRID))**2 + (DBLE(M3)-0.5D0*DBLE(NGRID))**2)
            ! RR     = 0.075D0*DBLE(NGRID)
             RR     = 0.1D0*DBLE(NGRID)
             IF(radius**2 < RR**2) THEN
            !    FI(M1,M2,M3) = 1.0D-3
                FI(M1,M2,M3) = 0.5D0
             ELSE
                FI(M1,M2,M3) = -0.5D0 * (70320.0D0 / 16706896.0D0) ! r = 0.1
            !    ! FI(M3,M2,M1) = -1.0D-3 * (29464.0D0 / 16747752.0D0) ! r = 0.075
             END IF
            ! IF(MG_model.EQ.3) THEN
            !    RR     = 0.1D0*DBLE(NGRID)
            !    radius = DSQRT((DBLE(M1)-0.5D0*DBLE(NGRID))**2+(DBLE(M2)-0.5D0*DBLE(NGRID))**2+(DBLE(M3)-0.5D0*DBLE(NGRID))**2)
            !    IF(radius**2.LT.RR**2) THEN
            !       FI(M3,M2,M1) = 5000.0D0
            !    ELSE
            !       FI(M3,M2,M1) = 1.0D0
            !    ENDIF
            ! ENDIF
            ! IF(MG_model.EQ.4) THEN
            !    RR     = 0.1D0*DBLE(NGRID)
            !    radius = DSQRT((DBLE(M1)-0.5D0*DBLE(NGRID))**2+(DBLE(M2)-0.5D0*DBLE(NGRID))**2+(DBLE(M3)-0.5D0*DBLE(NGRID))**2)
            !    IF(radius.LT.RR) THEN
            !       FI(M3,M2,M1) = 5000.0D0
            !    ELSE
            !       FI(M3,M2,M1) = 0.0D0
            !    ENDIF
            ! ENDIF
            !  IF(MG_model.EQ.5) THEN
            !     RR     = 0.1D0*DBLE(NGRID)
            !     radius = DSQRT((DBLE(M1)-0.5D0*DBLE(NGRID))**2+(DBLE(M2)-0.5D0*DBLE(NGRID))**2+(DBLE(M3)-0.5D0*DBLE(NGRID))**2)
            !     IF(radius.LT.RR) THEN
            !        FI(M3,M2,M1) = 5000.0D0
            !     ELSE
            !        FI(M3,M2,M1) = 0.0D0
            !     ENDIF
            !  ENDIF
            ! -------------------------------
            ! Spherical overdensity tests end
            ! -------------------------------

           END DO
        END DO
     END DO
     !
     RETURN
  END IF 

!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (IN,X,Y,Z,D1,D2,D3,T1,T2,T3,T2W,D2W) &
!$OMP PRIVATE (I,J,K,I1,J1,K1)    
  DO IN=1,Nparticles         ! loop over particles 
	  X=XPAR(IN)
	  Y=YPAR(IN)
	  Z=ZPAR(IN)
	  I=INT(X)
	  J=INT(Y)

	  K=INT(Z)
	  if(I<1.or.J<1.or.K<1)write(*,'(a,i12,3x,3es13.5,3x,3i8)') ' Error: out off bounds: ',IN,X,Y,Z,I,J,K
	  if(I>NGRID.or.J>NGRID.or.K>NGRID)write(*,'(a,i12,3x,3f12.4,3x,3i8)') ' Error: out off bounds: ',IN,X,Y,Z,I,J,K

	  D1=X-FLOAT(I)
	  D2=Y-FLOAT(J)
	  D3=Z-FLOAT(K)
	  T1=1.-D1
	  T2=1.-D2
	  T3=1.-D3
	  T2W =T2*WPAR
	  D2W =D2*WPAR
	  I1=I+1
	     IF(I1.GT.NGRID)I1=1
	  J1=J+1
	     IF(J1.GT.NGRID)J1=1
	  K1=K+1
          IF(K1.GT.NGRID)K1=1
!$OMP ATOMIC
          FI(I ,J ,K ) =FI(I ,J ,K ) +T3*T1*T2W
!$OMP ATOMIC
	  FI(I1,J ,K ) =FI(I1,J ,K ) +T3*D1*T2W
!$OMP ATOMIC
	  FI(I ,J1,K ) =FI(I ,J1,K ) +T3*T1*D2W
!$OMP ATOMIC
	  FI(I1,J1,K ) =FI(I1,J1,K ) +T3*D1*D2W
   
!$OMP ATOMIC
	  FI(I ,J ,K1) =FI(I ,J ,K1) +D3*T1*T2W
!$OMP ATOMIC
	  FI(I1,J ,K1) =FI(I1,J ,K1) +D3*D1*T2W
!$OMP ATOMIC
	  FI(I ,J1,K1) =FI(I ,J1,K1) +D3*T1*D2W
!$OMP ATOMIC
	  FI(I1,J1,K1) =FI(I1,J1,K1) +D3*D1*D2W
	   
  END DO 
  
  Call TimingMain(3,1)
  
   D1 = 0. ; D2 =0.
  WRITE(*,*) ' finishing density ',seconds()
!!$OMP PARALLEL DO DEFAULT(SHARED) &
!!$OMP PRIVATE (M1,M2,M3) REDUCTION(+:D1,D2)
!    DO M3=1,NGRID
!       DO M2=1,NGRID
!          DO M1=1,NGRID
!	       D1 = D1 + FI(M1,M2,M3)
!	       D2 = D2 + FI(M1,M2,M3)**2
!          END DO
!       END DO
!    END DO
!    D1 = D1/(float(NGRID))**3
!    D2 = sqrt(D2/(float(NGRID))**3)
!    write(*,*) ' Finished density: aver/rms=',D1,D2
    
END SUBROUTINE DENSIT

!-------------------------------------------------------------
!
!         Make three projections of density field
!         
!
SUBROUTINE ProjectDENSIT
!    
!-------------------------------------------------------------
Real*4, Allocatable, dimension(:,:) :: DenXY,DenXZ,DenYZ
Real*8 :: D1,D2,D3
character*120 :: Name

   Call Timing(5,-1)
ALLOCATE(DenXY(NGRID,NGRID),DenXZ(NGRID,NGRID),DenYZ(NGRID,NGRID))
moment = ISTEP
    t0 = seconds()
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2)
    DO M2=1,NGRID
       DO M1=1,NGRID
	  DenXY(M1,M2) = 0.
	  DenXZ(M1,M2) = 0.
	  DenYZ(M1,M2) = 0.
      END DO
    END DO
    t1 = seconds()
    write(*,*) ' init finished time = ',t1-t0
    write(*,*) ' Ngrid = ',NGRID
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3,D1,D2,D3)
    DO M3=1,NGRID
       DO M2=1,NGRID
          D1 = 0.
          D2 = 0.
          D3 = 0.
       DO M1=1,NGRID
          D1 = D1 + FI(M1,M2,M3)+1.
          D2 = D2 + FI(M2,M1,M3)+1.
          D3 = D3 + FI(M2,M3,M1)+1.
       end DO
	  DenYZ(M2,M3) = DenYZ(M2,M3) +D1
	  DenXZ(M2,M3) = DenXZ(M2,M3) +D2
	  DenXY(M2,M3) = DenXY(M2,M3) +D3
      END DO
    END DO
    t2 = seconds()
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3)
    DO M3=1,NGRID
       DO M2=1,NGRID
	  DenYZ(M2,M3) = DenYZ(M2,M3)/NGRID
	  DenXZ(M2,M3) = DenXZ(M2,M3)/NGRID
	  DenXY(M2,M3) = DenXY(M2,M3)/NGRID
       end DO
    end DO
    write(*,*) '      finished time = ',t2-t1
    write(*,'(a,12es12.4)') ' DenXY(1..2,1..2) : ',DenXY(1,1),DenXY(2,1),DenXY(1,2),DenXY(2,2)

      write(Name,'(2(a,i4.4),3(a,i3.3))')'ProjectionXY.',moment,'.',Nrealization,   &
           '.dat'
      OPEN(19,FILE=TRIM(Name),STATUS='UNKNOWN',form='unformatted')
      write(Name,'(2(a,i4.4),3(a,i3.3))')'ProjectionXZ.',moment,'.',Nrealization,   &
           '.dat'
      OPEN(20,FILE=TRIM(Name),STATUS='UNKNOWN',form='unformatted')
      write(Name,'(2(a,i4.4),3(a,i3.3))')'ProjectionYZ.',moment,'.',Nrealization,   &
           '.dat'
      OPEN(21,FILE=TRIM(Name),STATUS='UNKNOWN',form='unformatted')

      write(*,*) AEXPN,ISTEP,NGRID
      write(19) AEXPN,ISTEP,NGRID
      write(19) DenXY
      write(20) AEXPN,ISTEP,NGRID
      write(20) DenXZ
      write(21) AEXPN,ISTEP,NGRID
      write(21) DenYZ


      close(19)
      close(20)
      close(21)
      Call Timing(5,1)
  end SUBROUTINE ProjectDENSIT

!---------------------------------------
!                   generate a vector of gaussian numbers    
SUBROUTINE getGauss(Gg,jp,kp,N)
!
!---------------------------------------
  use Random
  use LUXURY
  Integer*8 :: kp,jp
  Real*4 :: Gg(N)

    Ns =SeedsPage(jp,kp) 
    lux = 2
    Call rluxgo(lux,Ns,0,0)
    Do i=1,N
      Gg(i) = GAUSS3(gSet,iFlag)
    end Do
END SUBROUTINE getGauss
!

!-----------------------------------------------------
!                              Density Field    
SUBROUTINE DENSITrsd(iSwitch,NGRID_old)
!
!           iSwitch = 0 - real space
!                   = 1,2,3 - redshift space 
!----------------------------------------------------
use Random
Real*8 :: ss,X,Y,Z,D1,D2,D3,T1,T2,T3
Integer*8, parameter :: Nslip= 160481183, Nf=4538127
Integer*8 :: IN, ip,jp,kp,jpp,kpp
Integer*4, save :: Nseed =198239321
Real*4, allocatable :: Gg(:)
Call Timing(7,-1)      ! start reading time

      	XN   =FLOAT(NGRID)+1.-1.E-7
        YN   =FLOAT(NGRID)
        Nseed = 121071 + Nseed
        Xscale = Float(NGRID)/NGRID_old
        write(*,*) ' Inside Densit: Ngrid =',NGRID,iSwitch
        Nthreads = OMP_GET_MAX_THREADS()
!				       Subtract mean density
!$OMP PARALLEL DO DEFAULT(SHARED) & 
!$OMP PRIVATE (M1,M2,M3)
      DO M3=1,NGRID
       DO M2=1,NGRID
	      DO M1=1,NGRID
	        FI(M1,M2,M3) = -1.
	      END DO
       END DO
      END DO

      W     = FLOAT(NGRID)**3/float(Nparticles)
      ! factor = sqrt(AEXPN/(Om0+AEXPN**3*Oml0))/100. ! V(km/s)->space
      ! Vscale = (Box/Ngrid)100/Aexpn  ! km/s
      factor = sqrt(AEXPN/(Om+AEXPN**3*OmL))/AEXPN
      sigFactor = sigv/100.*AEXPN*Ngrid/Box
      PARTW = W  
      write(*,'(a,f8.3,a,i5,a,8es12.4)') ' DENSITrsd:Particle weight = ',PARTW, &
           ' Switch=',iSwitch, &
           ' Factors= ',factor,sigFactor
      xmin =1.e+5
      xmax =-1.e+5
      ymin =1.e+5
      ymax =-1.e+5
      zmin =1.e+5
      zmax =-1.e+5
      If(iSwitch/=0)Then
         !write(*,*) ' Allocate Gg: ',Nrow
         Allocate(Gg(Nrow))
         Gg = 0.
      end If
       
!$OMP PARALLEL DO DEFAULT(SHARED) & 
!$OMP PRIVATE (IN,X,Y,Z,I,J,K,D1,D2,D3,T1,T2,T3,T2W,D2W,I1,J1,K1, &
!$OMP       vv,ip,jpp,kpp,jp,kp,Gg)
      Do kp =1,NROW
         if(mod(kp,250)==0.and.iSwitch==1)write(*,*) ' k=',kp
         kpp = (kp-1)*NROW**2
      Do jp =1,Nrow
         jpp = (jp-1)*Nrow
         If(iSwitch/=0)Then
!$OMP critical            
            Call getGauss(Gg,jp,kp,Nrow)
!$OMP end critical
        end If
      Do ip =1,Nrow
         IN = ip + jpp +kpp
         
        If(iSwitch == 0)Then
            X = (XPAR(IN)-1.)*Xscale+1.
            Y = (YPAR(IN)-1.)*Xscale+1.
            Z = (ZPAR(IN)-1.)*Xscale+1.
         Else
            If(dens(IN)>DensThr)Then
               vv = sigFactor*Gg(ip)*(dens(IN)-DensThr)**0.3333
            Else
               vv  = 0.
            end If
               SELECT CASE (iSwitch)
                CASE (1)
                   Z = (XPAR(IN)-1.)*Xscale+1. + (Vx(IN)*Xscale+vv)*factor
                   Y = (YPAR(IN)-1.)*Xscale+1.
                   X = (ZPAR(IN)-1.)*Xscale+1.
                CASE (2)
                   X = (XPAR(IN)-1.)*Xscale+1.
                   Z = (YPAR(IN)-1.)*Xscale+1. + (Vy(IN)*Xscale+vv)*factor
                   Y = (ZPAR(IN)-1.)*Xscale+1.
                CASE (3)
                   X = (XPAR(IN)-1.)*Xscale+1.
                   Y = (YPAR(IN)-1.)*Xscale+1.
                   Z = (ZPAR(IN)-1.)*Xscale+1. + (Vz(IN)*Xscale+vv)*factor
                End SELECT
      end If

           IF(X.ge.NGRID+1.)X=X-NGRID
           IF(Y.ge.NGRID+1.)Y=Y-NGRID
           IF(Z.ge.NGRID+1.)Z=Z-NGRID

           IF(X.lt.1.)X=X+NGRID
           IF(Y.lt.1.)Y=Y+NGRID
           IF(Z.lt.1.)Z=Z+NGRID

	   I=INT(X)
	   J=INT(Y)
	   K=INT(Z)
           If(I.le.0)write (*,*) ' X:',X,Y,Z,' Irow=',IROW,IN
           If(J.le.0)write (*,*) ' Y:',X,Y,Z,' Irow=',IROW,IN
           If(K.le.0)write (*,*) ' Z:',X,Y,Z,' Irow=',IROW,IN
           If(I.gt.NGRID+1)write (*,*) ' X:',X,Y,Z,' Irow=',IROW,IN
           If(J.gt.NGRID+1)write (*,*) ' Y:',X,Y,Z,' Irow=',IROW,IN
           If(K.gt.NGRID+1)write (*,*) ' Z:',X,Y,Z,' Irow=',IROW,IN
                   !---------------------------------------- CIC
	        D1=X-FLOAT(I)
	        D2=Y-FLOAT(J)
	        D3=Z-FLOAT(K)
	        T1=1.-D1
	        T2=1.-D2
	        T3=1.-D3
	        T2W =T2*W
	        D2W =D2*W
	        I1=I+1
	           IF(I1.GT.NGRID)I1=1
	        J1=J+1
	           IF(J1.GT.NGRID)J1=1
	        K1=K+1
         IF(K1.GT.NGRID)K1=1
!$OMP Atomic         
	             FI(I ,J ,K ) =FI(I ,J ,K ) +T3*T1*T2W
!$OMP Atomic         
                     FI(I1,J ,K ) =FI(I1,J ,K ) +T3*D1*T2W
!$OMP Atomic         
	             FI(I ,J1,K ) =FI(I ,J1,K ) +T3*T1*D2W
!$OMP Atomic         
	             FI(I1,J1,K ) =FI(I1,J1,K ) +T3*D1*D2W
!$OMP Atomic         
	             FI(I ,J ,K1) =FI(I ,J ,K1) +D3*T1*T2W
!$OMP Atomic         
	             FI(I1,J ,K1) =FI(I1,J ,K1) +D3*D1*T2W
!$OMP Atomic         
	             FI(I ,J1,K1) =FI(I ,J1,K1) +D3*T1*D2W
!$OMP Atomic         
	             FI(I1,J1,K1) =FI(I1,J1,K1) +D3*D1*D2W 
                !------------------------------------------ NGP
		!      FI(I ,J ,K ) =FI(I ,J ,K ) + W
           ENDDO
           ENDDO
           ENDDO

      If(iSwitch/=0)DEALLOCATE(Gg)

         Call Timing(7,1)      ! start reading time
END SUBROUTINE DENSITrsd
!--------------------------------------------------                     
!                    : Store and retrieve seeds for parallelization 
!                      of random number generator
!--------------------------------------------------                     
SUBROUTINE SetRandom
   use Tools
   use LUXURY
   use Random
   write(*,*) ' Inside SetRandom'
   write(*,*) NROW,NGRID
   
   ALLOCATE(SeedsPage(NROW,NROW))

      Ns     = Nseed 
   Do j=1,NROW 
      Do i=1,NROW 
         SeedsPage(i,j) = Ns
         dummy = RANDd(Ns) 
      End Do
   End Do
   Nseed = Ns    ! new seed
   !write(*,'(a,6i12)') ' Initialized random seeds: ',(SeedsPage(i,1),i=1,6)
End SUBROUTINE SetRandom        


SUBROUTINE DumpParticles
!
use Tools
use LUXURY
use Random
Real*8 :: ss,X,Y,Z,D1,D2,D3,T1,T2,T3
Integer*8 :: IN,im,ish,iadd,ifr,Ng3,Ng2
Integer*4,parameter :: Npd = 100
Integer*8 :: iPDF(-30:Npd)
Real*4    :: dLogR= 0.05
!Integer*4, parameter :: NpeakM = 300000
Integer*4 :: Npeak
Real*4    :: MassOne,MassOver
!Real*4, dimension(NpeakM)    :: Xpp,Ypp,Zpp,Dpp
character*120 :: Name

        write(*,'(a,i5,3x,3es12.4)') ' Inside DumpParticles: Ngrid =',NGRID,Box,Om
        Xscale = Box/NGRID
        Dscale = 2.774e+11*(Box/NROW)**3    ! mass scale
        MassOne= Om*Dscale                ! mass of the smallest particle
        MassOver = (float(NGRID)/NROW)**3
        Vscale = Xscale*100./AEXPN    ! scale V to km/s
        Nseed = 1238229

        ish=15485867; iadd = 1349973; ifr = 72421*109

        write(*,*) ' Bias         =',BiasPars(1:2)
        write(*,*) ' Mass particle=',MassOne
        write(*,*) ' Density one particle/cell  =',MassOver,seconds()
        !moment = max(INT(100.*(1./AEXPN-1.)+0.5),0)
        moment = ISTEP
        write(Name,'(2(a,i4.4),3(a,i3.3))')'ParticlesDens.',moment,   &
                                           '.',Nrealization,'.dat'
        Dth2 = BiasPars(2)
        fc   = BiasPars(3)      ! fraction for density
        OPEN(18,FILE=TRIM(Name),STATUS='UNKNOWN')
        write(Name,'(2(a,i4.4),3(a,i3.3))')'Particles.',moment,   &
                                           '.',Nrealization,'.dat'
        OPEN(19,FILE=TRIM(Name),STATUS='UNKNOWN',FORM='UNFORMATTED')

        write(18,'(a)')HEADER
        write(18,'(a,es13.4,2(a,es12.3),a,i4)') &
            ' PartThresh = ',BiasPars(2),   &
            ' Fraction = ',BiasPars(3), &
            ' Aexpn = ',AEXPN,' Step = ',ISTEP               
      ip = 0
!!$OMP PARALLEL DO DEFAULT(SHARED) & 
!!$OMP PRIVATE (IN,X,Y,Z,I,J,K,I1,J1,K1,I0,J0,K0,dd,D1,D2,D3,T1,T2,T3,Im) REDUCTION(+:ip)
      Do IN =1,Nparticles
                         !--- take random fraction above Biaspars(2) = Dth2
         X = XPAR(IN)     !--- find 8 nodes
         Y = YPAR(IN)
         Z = ZPAR(IN)
         I = INT(X)
         J = INT(Y)
         K = INT(Z)
         I1 = I + 1
         J1 = J + 1
         K1 = K +1
         If(I1>NGRID)I1=1
         If(J1>NGRID)J1=1
         If(K1>NGRID)K1=1
           D1=X-FLOAT(I)
	   D2=Y-FLOAT(J)
	   D3=Z-FLOAT(K)
	   T1=1.-D1
	   T2=1.-D2
	   T3=1.-D3
	   I1=I+1
	      IF(I1.GT.NGRID)I1=1
	   J1=J+1
	      IF(J1.GT.NGRID)J1=1
	   K1=K+1
              IF(K1.GT.NGRID)K1=1
           dd=  FI(I ,J ,K )*T3*T1*T2 + &      ! density at particle position
                FI(I1,J ,K )*T3*D1*T2 + &
                FI(I ,J1,K )*T3*T1*D2 + &
                FI(I1,J1,K )*T3*D1*D2 + &
                FI(I ,J ,K1)*D3*T1*T2 + &
                FI(I1,J ,K1)*D3*D1*T2 + &
                FI(I ,J1,K1)*D3*T1*D2 + &
                FI(I1,J1,K1)*D3*D1*D2 + 1.
           If(dd>Dth2)Then
              Im  = mod(ifr*IN+iadd,ish)           ! proxi for random number
                              !If(FI(I,J,K)>Dth2)Then
                !   ff = fc*(1-Dth2/FI(I,J,K))**0.3333*(Dth2/FI(I,J,K))
                !Else
                !   ff = 0.
                !EndIf
               If(Im/float(ish)<fc)Then
                 ip = ip +1
                 !if(ip>Ngalaxies)Stop ' Too many galaxies: change Ngalaxies parameter in mod_density.f90'
                 Xbd = (X-1.)*Xscale        ! coords in Mpch units
                 Ybd = (Y-1.)*Xscale
                 Zbd = (Z-1.)*Xscale
                 VXd = VX(IN)*Vscale  !+ dd*GAUSS(Nseed) !velocity in km/s
                 VYd = VY(IN)*Vscale  !+ dd*GAUSS(Nseed)
                 VZd = VZ(IN)*Vscale  !+ dd*GAUSS(Nseed)               
                 !! dd =  (dd+1.) 
!!$OMP critical
               if(ip<100)write(18,'(3f9.3,4es13.4)') Xbd,Ybd,Zbd,VXd,VYd,VZd,dd             
                 write(19) Xbd,Ybd,Zbd,VXd,VYd,VZd,dd             
!!$OMP end critical 
               end If
            end If
         end Do
         write(*,'(3x,a,i12,3x,a,es12.4)') 'Number of dumped particles =',ip, &
                    'Fraction of particles=',float(ip)/Nparticles
       close (18)
END SUBROUTINE DumpParticles
! 
!---------------------------------------
!               DM density at particle position 
!
    SUBROUTINE DENSPARTall
!
! 
!---------------------------------------
use Tools
use LUXURY
use Random
Real*8 :: ss,X,Y,Z,D1,D2,D3,T1,T2,T3
Integer*8 :: IN
Integer*4,parameter :: Npd = 100
Integer*8 :: iPDF(-30:Npd)
Real*4    :: dLogR= 0.05
!Integer*4, parameter :: NpeakM = 300000
Integer*4 :: Npeak
Real*4    :: MassOne,MassOver
character*120 :: Name

Call Timing(5,-1)      ! start timing
        write(*,'(a,i5,3x,3es12.4)') ' Inside DensPart: Ngrid =',NGRID,Box,Om

!$OMP PARALLEL DO DEFAULT(SHARED) & 
!$OMP PRIVATE (IN,X,Y,Z,I,J,K,I1,J1,K1,D1,D2,D3,T1,T2,T3)
      Do IN =1,Nparticles
         X = XPAR(IN)     !--- find 8 nodes
         Y = YPAR(IN)
         Z = ZPAR(IN)
         I = INT(X)
         J = INT(Y)
         K = INT(Z)
         I1 = I + 1
         J1 = J + 1
         K1 = K +1
         If(I1>NGRID)I1=1
         If(J1>NGRID)J1=1
         If(K1>NGRID)K1=1
	   D1=X-FLOAT(I)
	   D2=Y-FLOAT(J)
	   D3=Z-FLOAT(K)
	   T1=1.-D1
	   T2=1.-D2
	   T3=1.-D3
	   I1=I+1
	      IF(I1.GT.NGRID)I1=1
	   J1=J+1
	      IF(J1.GT.NGRID)J1=1
	   K1=K+1
              IF(K1.GT.NGRID)K1=1
           dens(IN)=  FI(I ,J ,K )*T3*T1*T2 + &      ! density at particle position
                FI(I1,J ,K )*T3*D1*T2 + &
                FI(I ,J1,K )*T3*T1*D2 + &
                FI(I1,J1,K )*T3*D1*D2 + &
                FI(I ,J ,K1)*D3*T1*T2 + &
                FI(I1,J ,K1)*D3*D1*T2 + &
                FI(I ,J1,K1)*D3*T1*D2 + &
                FI(I1,J1,K1)*D3*D1*D2 + 1.

       end Do
         Call Timing(5,1)    
       END SUBROUTINE DENSPARTall
!---------------------------------------
!               Biasing scheme that uses DM density at particle position 
!
SUBROUTINE DENSPART
!
! 
!---------------------------------------
use Tools
use LUXURY
use Random
Real*8 :: ss,X,Y,Z,D1,D2,D3,T1,T2,T3
Integer*8 :: IN,im,ish,iadd,ifr,Ng3,Ng2
Integer*4,parameter :: Npd = 100
Integer*8 :: iPDF(-30:Npd)
Real*4    :: dLogR= 0.05
!Integer*4, parameter :: NpeakM = 300000
Integer*4 :: Npeak
Real*4    :: MassOne,MassOver
!Real*4, dimension(NpeakM)    :: Xpp,Ypp,Zpp,Dpp
character*120 :: Name

Call Timing(5,-1)      ! start timing
        write(*,'(a,i5,3x,3es12.4)') ' Inside DensPart: Ngrid =',NGRID,Box,Om
        Xscale = Box/NGRID
        Dscale = 2.774e+11*(Box/NROW)**3    ! mass scale
        MassOne= Om*Dscale                ! mass of the smallest particle
        MassOver = (float(NGRID)/NROW)**3
        Vscale = Xscale*100./AEXPN    ! scale V to km/s
        Nseed = 1238229
        !ff    = 1.5
        !trans = 20000.

        ish=15485867; iadd = 1349973; ifr = 72421*109
        Ng3 = (Ngrid-1)**2
        Ng2 = (Ngrid-1)
        
        write(*,*) ' Bias         =',BiasPars(1:2)
        write(*,*) ' Mass particle=',MassOne
        write(*,*) ' Density one particle/cell  =',MassOver,seconds()
        !moment = max(INT(100.*(1./AEXPN-1.)+0.5),0)
        moment = ISTEP
        write(Name,'(2(a,i4.4),3(a,i3.3))')'DensDistrGal.',moment,   &
                                           '.',Nrealization,'.dat'

        OPEN(18,FILE=TRIM(Name),STATUS='UNKNOWN')
        write(18,'(a)')HEADER
        Npeak = 0
        Dth   = BiasPars(1)      ! threshold for maxima
        Dth2  = BiasPars(2)      ! threshold for density
        fc    = BiasPars(3)      ! fraction for density
!$OMP PARALLEL DO DEFAULT(SHARED) &                          ! Mark maxima density above BiasPar(1) = Dth
!$OMP PRIVATE (k,K1,Km,J,J1,Jm,I,I1,Im,Fmax,dmn,ff) REDUCTION(+:Npeak)
        Do k =1,NGRID
	   K1=K+1
              IF(K1.GT.NGRID)K1=1                 
	   Km=K-1
              IF(Km.lt.1)Km=NGRID                 
              Do j=1,NGRID
                J1=J+1
                IF(J1.GT.NGRID)J1=1                 
	        Jm=J-1
                IF(Jm.lt.1)Jm=NGRID  
              Do i=1,NGRID
                I1=I+1
                IF(I1.GT.NGRID)I1=1                 
	        Im=I-1
                IF(Im.lt.1)Im=NGRID  
                Fmax = max(FI(I1,J,K),FI(Im,J,K),FI(I,J1,K),FI(I,Jm,K),FI(I,J,K1), FI(I,J,Km) )
                Fmax = max(Fmax,FI(Im,Jm,Km),FI(I1,Jm,Km),FI(I1,Jm,Km),FI(I1,J1,Km))
                Fmax = max(Fmax,FI(Im,Jm,K1),FI(I1,Jm,K1),FI(I1,Jm,K1),FI(I1,J1,K1))
                !Im = mod(7121_8*(i+j+k)+18411_8,13447_8)     ! proxi for random number
                !Im  = mod(ifr*(i+j*Ng2+k*Ng3)+iadd,ish)       ! proxi for random number
                !If(FI(I,J,K)>Dth2)Then
                !   ff = fc*(1-Dth2/FI(I,J,K))**0.3333*(Dth2/FI(I,J,K))
                !Else
                !   ff = 0.
                !EndIf
                !If((Fmax<FI(I,J,K).and.FI(I,J,K)>Dth).or.(Im/float(ish)<fc .and.FI(I,J,K)>Dth2))Then
                If(Fmax<FI(I,J,K).and.FI(I,J,K)>Dth)Then
                   Npeak = Npeak +1
                   FI(I,J,K) = FI(I,J,K) + 2.e6          !-- mark as maximum
                End If
              End Do
            End Do
         End Do
         write(*,*) ' Number of density peaks = ',Npeak,seconds()

        
      Ngalaxies = max(2*Npeak,10000000)
      myMemory =Memory(7_8*Ngalaxies)
      ALLOCATE(dens(Ngalaxies))
      Allocate (Xb(Ngalaxies),Yb(Ngalaxies),Zb(Ngalaxies))
      Allocate (VXb(Ngalaxies),VYb(Ngalaxies),VZb(Ngalaxies))         
      ip = 0
!$OMP PARALLEL DO DEFAULT(SHARED) & 
!$OMP PRIVATE (IN,X,Y,Z,I,J,K,I1,J1,K1,I0,J0,K0,dd,D1,D2,D3,T1,T2,T3,Im)
      Do IN =1,Nparticles
         !X = XPAR(IN)+0.5     !--- find nearest node
         !Y = YPAR(IN)+0.5
         !Z = ZPAR(IN)+0.5
         !I = INT(X)
         !J = INT(Y)
         !!K = INT(Z)
         !If(I>NGRID)I=1
         !If(J>NGRID)J=1
         !If(K>NGRID)K=1
         X = XPAR(IN)     !--- find 8 nodes
         Y = YPAR(IN)
         Z = ZPAR(IN)
         I = INT(X)
         J = INT(Y)
         K = INT(Z)
         I1 = I + 1
         J1 = J + 1
         K1 = K +1
         If(I1>NGRID)I1=1
         If(J1>NGRID)J1=1
         If(K1>NGRID)K1=1
         I0 = -1 ; J0 = -1 ;  K0 = -1   ! find if any of nodes is a maximum
         If(FI(I,J,K)>1.e6)Then
            I0 = I ; J0 = J ; K0 = K
         End If
         If(FI(I1,J,K)>1.e6)Then
            I0 = I1 ; J0 = J ; K0 = K
         End If
         If(FI(I,J1,K)>1.e6)Then
            I0 = I ; J0 = J1 ; K0 = K
         End If
         If(FI(I1,J1,K)>1.e6)Then
            I0 = I1 ; J0 = J1 ; K0 = K
         End If
         
         If(FI(I,J,K1)>1.e6)Then
            I0 = I ; J0 = J ; K0 = K1
         End If
         If(FI(I1,J,K1)>1.e6)Then
            I0 = I1 ; J0 = J ; K0 = K1
         End If
         If(FI(I,J1,K1)>1.e6)Then
            I0 = I ; J0 = J1 ; K0 = K1
         End If
         If(FI(I1,J1,K1)>1.e6)Then
            I0 = I1 ; J0 = J1 ; K0 = K1
         End If
         
         If(I0>0)Then   !--- maximum
	    !dd = FI(I0,J0,K0)-2.e6
            !D1 = sqrt(abs(1.-Dth/dd))   ! statistical weight
            !Im  = mod(ifr*IN+iadd,ish)           ! proxi for random number
            !If(Im/float(ish)<D1)Then
!$OMP critical            
            FI(I0,J0,K0)= FI(I0,J0,K0)-2.e6            !--- remove maximum flag, restore density
            dd =  FI(I0,J0,K0)
            ip = ip +1
              if(ip>Ngalaxies)write(*,'(a,2i12,3x,3es13.4)') ' Too many galaxies:', ip,IN,X,Y,Z
            Xb(ip) = (X-1.)*Xscale        ! coords in Mpch units
            Yb(ip) = (Y-1.)*Xscale
            Zb(ip) = (Z-1.)*Xscale
            VXb(ip) = VX(IN)*Vscale  !+ dd*GAUSS(Nseed) !velocity in km/s
            VYb(ip) = VY(IN)*Vscale  !+ dd*GAUSS(Nseed)
            VZb(ip) = VZ(IN)*Vscale  !+ dd*GAUSS(Nseed)               
            dens(ip) = dd+1.
!$OMP end critical
	    !EndIf
         Else                   !--- take random fraction above Biaspars(2) = Dth2
	   D1=X-FLOAT(I)
	   D2=Y-FLOAT(J)
	   D3=Z-FLOAT(K)
	   T1=1.-D1
	   T2=1.-D2
	   T3=1.-D3
	   I1=I+1
	      IF(I1.GT.NGRID)I1=1
	   J1=J+1
	      IF(J1.GT.NGRID)J1=1
	   K1=K+1
              IF(K1.GT.NGRID)K1=1
           dd=  FI(I ,J ,K )*T3*T1*T2 + &      ! density at particle position
                FI(I1,J ,K )*T3*D1*T2 + &
                FI(I ,J1,K )*T3*T1*D2 + &
                FI(I1,J1,K )*T3*D1*D2 + &
                FI(I ,J ,K1)*D3*T1*T2 + &
                FI(I1,J ,K1)*D3*D1*T2 + &
                FI(I ,J1,K1)*D3*T1*D2 + &
                FI(I1,J1,K1)*D3*D1*D2 + 1.
           If(dd>Dth2)Then
              D1 = sqrt(abs(1.-Dth2/dd))*(Dth2/dd)   ! statistical weight
              Im  = mod(ifr*IN+iadd,ish)           ! proxi for random number
                              !If(FI(I,J,K)>Dth2)Then
                !   ff = fc*(1-Dth2/FI(I,J,K))**0.3333*(Dth2/FI(I,J,K))
                !Else
                !   ff = 0.
                !EndIf
               If(Im/float(ish)<fc*D1)Then
!$OMP critical
                 ip = ip +1
                 if(ip>Ngalaxies)Stop ' Too many galaxies: change Ngalaxies parameter in mod_density.f90'
                 Xb(ip) = (X-1.)*Xscale        ! coords in Mpch units
                 Yb(ip) = (Y-1.)*Xscale
                 Zb(ip) = (Z-1.)*Xscale
                 VXb(ip) = VX(IN)*Vscale  !+ dd*GAUSS(Nseed) !velocity in km/s
                 VYb(ip) = VY(IN)*Vscale  !+ dd*GAUSS(Nseed)
                 VZb(ip) = VZ(IN)*Vscale  !+ dd*GAUSS(Nseed)               
                 dens(ip) = -(dd+1.) 
           !       write(*,'(a,i12,f10.2,f9.4)') ' rand: ',ip,dens(ip),Im/float(ish)             
!$OMP end critical 
               end If
            end If
         End If
       end Do

       Ngalaxies = ip
        iPDF(:) = 0

        write(*,*) ' finished denspart =',Ngalaxies,seconds()
       write(18,'(a,i10,a,f8.3,a,es13.4,2(a,es12.3),a,I5)') &
            ' Ngalaxies=',Ngalaxies, &
            ' MaxDensThresh = ',BiasPars(1),    &
            ' PartThresh = ',BiasPars(2),   &
            ' Normalization = ',BiasPars(3), &
            ' Aexpn = ',AEXPN,' Step = ',ISTEP

!$OMP PARALLEL DO DEFAULT(SHARED) & 
!$OMP PRIVATE (IN,dd,ind)
        Do IN =1,Ngalaxies
           dd = abs(dens(IN))
           ind = INT((log10(dd)/dLogR+1000.))-1000
           ind = max(min(Npd,ind),-30)
!$OMP atomic
           iPDF(ind) =iPDF(ind) +1
        end Do
        !do IN =1,Nparticles
        !   if(dens(IN)>30..and.Xpar(IN)<300. .and. Ypar(IN)<300..and. Zpar(IN)<300.)write(*,'(i11,3x,2f9.3,3x,3f11.3)') &
        !          IN,dens(IN),dens2(IN),Xpar(IN),Ypar(IN),Zpar(IN)
        !end do
        write(18,'(5x,a)') ' density      dN       dN/(Ndrho)    rho*PDF  Particles PDF'
        Do i=-30,Npd
           If(iPDF(i)>0)Then
              dd = 10.**(i*dLogR)
              d2 = 10.**((i+1)*dLogR)
              write(18,'(3x,es12.3,i11,2es13.4)') (d2+dd)/2.,iPDF(i), &
                   iPDF(i)/(d2-dd)/FLOAT(Ngalaxies), &
                   iPDF(i)/(d2-dd)/FLOAT(Ngalaxies)*(d2+dd)/2.
          end If
       end Do
       close (18)
         Call Timing(5,1)    
       END SUBROUTINE DENSPART

!-----------------------------------------------------
!		                  Compute mean density and rms
      SUBROUTINE DENTES(DELR)
!-----------------------------------------------------
     real*8 :: SUM1,SUM2
     real*8, parameter :: dlog =0.05   
         Call Timing(3,-1)      ! start reading time
      SUM1 = 0.
      SUM2= 0.
      Nn  = 0
      Am  = 0.
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE ( k,j,i) &
!$OMP REDUCTION(+:SUM1,SUM2)
      DO K=1,NGRID
      DO J=1,NGRID
      DO I=1,NGRID
       SUM1 = SUM1 + FI(I,J,K)
       SUM2=  SUM2 + (FI(I,J,K))**2
      ENDDO
      ENDDO
      ENDDO
      Total =(FLOAT(NGRID))**3
      DENM = SUM1/Total
      DELR = DSQRT(SUM2/Total-DENM**2)
      WRITE (*,150)  DELR,DENM
      WRITE (17,150) DELR,DENM
150     format(20x,'Density is in units of average density', &
                    ' in the Box',/20x,   &
                    ' RMS Delta Rho/Rho   =',G11.4,/20x,  &
                    ' Mean Delta Rho/Rho  =',G11.4)
         Call Timing(3,1)  
     END SUBROUTINE DENTES

!-----------------------------------------------------
!		          Compute statistics of Density
      SUBROUTINE DensDistr
!-----------------------------------------------------
      real*8, parameter :: dlog =0.001   
      real*8 :: SUM1,SUM2
      Integer*8,   Allocatable,DIMENSION(:)   :: nCells
      Real*8,      Allocatable,DIMENSION(:)   :: den
      Integer*8,   Allocatable,DIMENSION(:,:) :: nCellsT
      Real*8,      Allocatable,DIMENSION(:,:) :: denT
      Integer*4  :: OMP_GET_MAX_THREADS,OMP_GET_THREAD_NUM

      Call Timing(7,-1)      ! start reading time
         
      DensMax = 1.e2    !-- assumed density maximum
      DensMin = 0.01    !-- assumed density minimum
      iMin = log10(DensMin)/dlog-1
      iMax = log10(DensMax)/dlog+1
      iThreads = OMP_GET_MAX_THREADS() 
      allocate(nCells(iMin:iMax),den(iMin:iMax))
      allocate(nCellsT(iMin:iMax,iThreads),denT(iMin:iMax,iThreads))

      nCells(:)    =0
      den(:)       =0.
      nCellsT(:,:) =0
      denT(:,:)    =0.
      
      SUM1 = 0.
      SUM2= 0.
      Nn  = 0
      Am  = 0.
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE ( k,j,i) &
!$OMP REDUCTION(+:SUM1,SUM2)
      DO K=1,NGRID
      DO J=1,NGRID
      DO I=1,NGRID
       SUM1 = SUM1 + FI(I,J,K)+1.
       SUM2=  SUM2 + (FI(I,J,K)+1.)**2
      ENDDO
      ENDDO
      ENDDO
      Total =(FLOAT(NGRID))**3
      DENM = SUM1/Total
      DELR = DSQRT(SUM2/Total-DENM**2)
      WRITE (*,150)  DELR,DENM,NGRID
      WRITE (18,150) DELR,DENM,NGRID
150     format(20x,'Density is in units of average density', &
                    ' in the Box',/20x,   &
                    ' RMS Delta Rho/Rho   =',G11.4,/20x,  &
                    ' Mean Delta Rho      =',G11.4,/20x,  &
                    ' Number grid points  =',i5)

!$OMP PARALLEL DO DEFAULT(SHARED)  &
!$OMP PRIVATE ( ind,k,j,i,iOMP)
      Do K=1,NGRID
         iOMP = OMP_GET_THREAD_NUM()+1
      DO J=1,NGRID
         DO I=1,NGRID
           ind = Floor(log10(FI(I,J,K)+1.)/dlog)
           ind = MIN(MAX(ind,iMin),iMax)
           nCellsT(ind,iOMP) = nCellsT(ind,iOMP) +1
           denT(ind,iOMP)    = denT(ind,iOMP) +FI(I,J,K)+1.
        END DO
     END DO
  END Do
  DO i=iMin,iMax            ! sum results of all threads
     DO iP=1,iThreads
        nCells(i) = nCells(i) + nCellsT(i,iP)
        den(i)    = den(i)    + denT(i,iP)
     end DO
  END DO
  dNorm = float(Ngrid)**3
  write(18,*)' dens_left   dens_right  density   dens*dN/d(dens)/Ncells  cells'
  DO i=iMin,iMax           ! print results
     d1 = 10.**(i*dlog)
     if(i==iMin)d1 =0.
    d2 = 10.**((i+1)*dlog) 
     den(i) = den(i)/max(nCells(i),1)
     write(18,'(3es12.4,3x,es13.5,i12)') d1,d2,den(i),   &
          nCells(i)/(d2-d1)/dnorm*den(i),nCells(i)     
  END DO
  
      deallocate(nCells,den)
      deallocate(nCellsT,denT)
      
         Call Timing(7,1)  
END SUBROUTINE DensDistr

!--------------------------------------------
!
!          Make density field of galaxies    
!
subroutine DensGal(iSwitch)
!
!--------------------------------------------
use Tools
character*120 :: Name

Real*8 :: ss,X,Y,Z,D1,D2,D3,T1,T2,T3
Integer*8 ::  ip
Call Timing(5,-1)      ! start reading time

      	XN   =FLOAT(NGRID)+1.-1.E-7
        YN   =FLOAT(NGRID)
!				       Subtract mean density
!$OMP PARALLEL DO DEFAULT(SHARED) & 
!$OMP PRIVATE (M1,M2,M3)
      DO M3=1,NGRID
       DO M2=1,NGRID
	      DO M1=1,NGRID
	        FI(M1,M2,M3) = -1.
	      END DO
       END DO
      END DO

      W     = FLOAT(NGRID)**3/float(Ngalaxies)

      Xscale = NGRID/Box
      Vscale = Xscale/100.*sqrt(AEXPN/(Om+AEXPN**3*OmL))
      PARTW = W  
      write(*,'(/a,es12.4,a,i5,a,8es12.4)') ' DensGal: Particle weight = ',PARTW, &
           ' Switch=',iSwitch, &
           ' Factors= ',Xscale,Vscale
      xmin =1.e+5
      xmax =-1.e+5
      ymin =1.e+5
      ymax =-1.e+5
      zmin =1.e+5
      zmax =-1.e+5

       
!$OMP PARALLEL DO DEFAULT(SHARED) & 
!$OMP PRIVATE (ip,X,Y,Z,I,J,K,D1,D2,D3,T1,T2,T3,T2W,D2W,I1,J1,K1)
      Do ip = 1,Ngalaxies
         !if(mod(ip,50000)==0)write(*,*) ' galaxy=',ip
         If(iSwitch == 0)Then
            X = Xb(ip)*Xscale+1.
            Y = Yb(ip)*Xscale+1.
            Z = Zb(ip)*Xscale+1.
         Else
               SELECT CASE (iSwitch)
               CASE (1)
                    Z = Xb(ip)*Xscale+1. +VXb(ip)*Vscale
                    Y = Yb(ip)*Xscale+1.
                    X = Zb(ip)*Xscale+1.
                 CASE (2)
                    X = Xb(ip)*Xscale+1. 
                    Z = Yb(ip)*Xscale+1. +VYb(ip)*Vscale
                    Y = Zb(ip)*Xscale+1.
                CASE (3)
                    X = Xb(ip)*Xscale+1. 
                    Y = Yb(ip)*Xscale+1.
                    Z = Zb(ip)*Xscale+1. +VZb(ip)*Vscale
                End SELECT
      end If

           IF(X.ge.NGRID+1.)X=X-NGRID
           IF(Y.ge.NGRID+1.)Y=Y-NGRID
           IF(Z.ge.NGRID+1.)Z=Z-NGRID

           IF(X.lt.1.)X=X+NGRID
           IF(Y.lt.1.)Y=Y+NGRID
           IF(Z.lt.1.)Z=Z+NGRID

	   I=INT(X)
	   J=INT(Y)
	   K=INT(Z)
           If(I.le.0)write (*,'(a,3es15.5,a,2i12)') ' X:',X,Y,Z,' Irow=',ip
           If(J.le.0)write (*,'(a,3es15.5,a,2i12)') ' Y:',X,Y,Z,' Irow=',ip
           If(K.le.0)write (*,'(a,5es15.5,a,2i12)') ' Z:',Z,Xb(ip),VXb(ip),Xb(ip)*Xscale+1.,VXb(ip)*Vscale,' Irow=',ip
           If(I.gt.NGRID+1)write (*,'(a,3es15.5,a,2i12)') ' X:',X,Y,Z,' Irow=',ip
           If(J.gt.NGRID+1)write (*,'(a,3es15.5,a,2i12)') ' Y:',X,Y,Z,' Irow=',ip
           If(K.gt.NGRID+1)write (*,'(a,5es15.5,a,2i12)') ' Z:',Z,Xb(ip),VXb(ip),Xb(ip)*Xscale+1.,VXb(ip)*Vscale,' Irow=',ip
                   !---------------------------------------- CIC
	        D1=X-FLOAT(I)
	        D2=Y-FLOAT(J)
	        D3=Z-FLOAT(K)
	        T1=1.-D1
	        T2=1.-D2
	        T3=1.-D3
	        T2W =T2*W
	        D2W =D2*W
	        I1=I+1
	           IF(I1.GT.NGRID)I1=1
	        J1=J+1
	           IF(J1.GT.NGRID)J1=1
	        K1=K+1
         IF(K1.GT.NGRID)K1=1
!$OMP Atomic         
	             FI(I ,J ,K ) =FI(I ,J ,K ) +T3*T1*T2W
!$OMP Atomic         
                     FI(I1,J ,K ) =FI(I1,J ,K ) +T3*D1*T2W
!$OMP Atomic         
	             FI(I ,J1,K ) =FI(I ,J1,K ) +T3*T1*D2W
!$OMP Atomic         
	             FI(I1,J1,K ) =FI(I1,J1,K ) +T3*D1*D2W
!$OMP Atomic         
	             FI(I ,J ,K1) =FI(I ,J ,K1) +D3*T1*T2W
!$OMP Atomic         
	             FI(I1,J ,K1) =FI(I1,J ,K1) +D3*D1*T2W
!$OMP Atomic         
	             FI(I ,J1,K1) =FI(I ,J1,K1) +D3*T1*D2W
!$OMP Atomic         
	             FI(I1,J1,K1) =FI(I1,J1,K1) +D3*D1*D2W 
                !------------------------------------------ NGP
		!      FI(I ,J ,K ) =FI(I ,J ,K ) + W
           ENDDO

         Call Timing(5,1)      ! stop reading time

end subroutine DensGal



!------------------------------------------------------------
   end Module Density
   
