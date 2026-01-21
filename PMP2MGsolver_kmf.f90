!================================================= 
!================================================= 
!
!  Scalar field
!
!================================================= 
!================================================= 

!-------------------------------------------------------------
!
! red-black Gauss-Seidel relaxation iterations on level ilevel
!
! Note: ilevel runs from 0, which means relaxation on PM grid.
! This subroutine distinguishes between PM and multgrid grids,
! since the scalar field and density fields etc. are stored in
! different places in these two cases.
!
!-------------------------------------------------------------
SUBROUTINE relaxation_iterations_kmf(ilevel,redstep)
!-------------------------------------------------------------
  use Tools
  use ExtradofBackgroundData

  implicit none

  integer :: ilevel
  logical :: redstep
  real*8  :: ctilde,dx,dx2,dx4
  real*8  :: fct1,fct2,fct3,fct4,fct5,tmp1,tmp2,Sigma1,Sigma2
  integer :: nstart,nfinal,ngrid_level
  integer :: ioffset,joffset,koffset,koffset2
  integer :: i

  real*8  :: P,Q
  integer :: M1,M2,M3,M1u,M1l,M2u,M2l,M3u,M3l

  integer :: nid,myid

  integer :: ist
  real*8  :: phibar,dphidN,HoH0,phidot,phidot2

  IF(MG_test) WRITE(*,'(A,I5,F7.4)') 'Relaxation iterations on level',levelmax-ilevel,AEXPN

  ! number of grid points on the coarse level
  ngrid_level = NGRID/(2**ilevel)
  !
  ! simulation parameters
  ctilde = 2.99792458D3*DBLE(NGRID)/Box
  dx     = DBLE(2**ilevel) 
  dx2    = dx*dx
  dx4    = dx2*dx2
  !
  ! physical and numerical quantities
  fct1   = 0.5D0*(1.0D0/(AEXPN*kmf_lambda))**2
  fct2   = 3.0D0*kmf_beta*Om/AEXPN
  fct3   = kmf_n*kmf_K0*fct1**(kmf_n-1)                                 ! fct3 is gamma in notes
  fct4   = -fct3*(    DBLE(kmf_n)-1.0D0)*2.0D0*ctilde
  fct5   =  fct3*(2.0*DBLE(kmf_n)+1.0D0)/3.0D0                          ! fct5 is (2*n+1)/3 *gamma in notes; Baojiu-01-07-2021
  !
  ist = 1
  DO WHILE(BackgroundEvolution(ist,1).LT.AEXPN)
     ist=ist+1
  ENDDO
  !
  ! background scalar field at AEXPN
  phibar =  BackgroundEvolution(ist-1,2)+ &
         & (BackgroundEvolution(ist,  2)-BackgroundEvolution(ist-1,2))/ &
         & (BackgroundEvolution(ist,  1)-BackgroundEvolution(ist-1,1))* &
         & (AEXPN                       -BackgroundEvolution(ist-1,1))
  ! d\varphi/dN at AEXPN, N=ln(a)
  dphidN =  BackgroundEvolution(ist-1,3)+ &
         & (BackgroundEvolution(ist,  3)-BackgroundEvolution(ist-1,3))/ &
         & (BackgroundEvolution(ist,  1)-BackgroundEvolution(ist-1,1))* &
         & (AEXPN                       -BackgroundEvolution(ist-1,1))
  ! H/H0 at AEXPN, here H=a'/a=dN/d\tau with '=d/d\tau
  HoH0   =  BackgroundEvolution(ist-1,5)+ &
         & (BackgroundEvolution(ist,  5)-BackgroundEvolution(ist-1,5))/ &
         & (BackgroundEvolution(ist,  1)-BackgroundEvolution(ist-1,1))* &
         & (AEXPN                       -BackgroundEvolution(ist-1,1))
  !
  fct2 = fct2*DEXP(kmf_beta*phibar)
  !
  phidot = dphidN*HoH0                                                  ! (d\varphi/d\tau)/H0
  phidot2 = phidot**2                                                   ! [(d\varphi/d\tau)/H0]^2
  !
  ! offsets for cell access to array FI3
  ioffset  = 2**(levelmax-ilevel)
  joffset  = 2**(levelmax-ilevel)
  koffset  = 0
  koffset2 = 0 
  IF(ilevel.GT.0) THEN 
     koffset  = 2**(levelmax-ilevel)*(2**ilevel-2)
     koffset2 = 2**(levelmax-ilevel)*(2**ilevel-1)
  END IF
  !
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3,M1u,M1l,M2u,M2l,M3u,M3l,i) &
!$OMP PRIVATE (tmp1,tmp2,Sigma1,Sigma2) &
!$OMP PRIVATE (P,Q)
  DO M3=1,ngrid_level
     DO M2=1,ngrid_level
        DO M1=1,ngrid_level
           ! prepare the indices of the 6 neighbours on the 3-point stencil
           ! and apply periodic boundary condition
           M1u = M1+1; IF(M1u>ngrid_level) M1u = 1
           M1l = M1-1; IF(M1l<1          ) M1l = ngrid_level
           M2u = M2+1; IF(M2u>ngrid_level) M2u = 1
           M2l = M2-1; IF(M2l<1          ) M2l = ngrid_level
           M3u = M3+1; IF(M3u>ngrid_level) M3u = 1
           M3l = M3-1; IF(M3l<1          ) M3l = ngrid_level
           ! decide on red-black update ordering
           IF((     redstep).AND.(MOD(M1+M2+M3,2).EQ.0)) CYCLE
           IF((.NOT.redstep).AND.(MOD(M1+M2+M3,2).NE.0)) CYCLE
           !
           IF(ilevel.EQ.0) THEN
              Q = (FI2(M1u,M2 ,M3 )      + &
                 & FI2(M1l,M2 ,M3 )      + &
                 & FI2(M1 ,M2u,M3 )      + &
                 & FI2(M1 ,M2l,M3 )      + &
                 & FI2(M1 ,M2 ,M3u)      + &
                 & FI2(M1 ,M2 ,M3l)) 
              ! tmp1 is \nabla^i\varphi\nabla_i\varphi
              tmp1 = (FI2(M1u,M2 ,M3 ) - FI2(M1l,M2 ,M3 ))**2 + &
                   & (FI2(M1 ,M2u,M3 ) - FI2(M1 ,M2l,M3 ))**2 + &
                   & (FI2(M1 ,M2 ,M3u) - FI2(M1 ,M2 ,M3l))**2
              tmp1 = tmp1/(4.0D0*dx2)
              !
              IF(kmf_n.EQ.2) THEN
!                Sigma1 = tmp1+phidot2                 ! Baojiu-01-07-2021
                 Sigma1 = fct5*tmp1+fct3*phidot2       ! Baojiu-01-07-2021
                 Sigma2 = 1.0D0
              ENDIF
              IF(kmf_n.EQ.3) THEN
!                Sigma1 = tmp1**2+phidot2**2           ! Baojiu-01-07-2021
                 Sigma1 = fct5*tmp1**2+fct3*phidot2**2 ! Baojiu-01-07-2021
                 Sigma2 = tmp1
              ENDIF
!             Sigma1 = Sigma1*fct3+1.0D0
              Sigma1 = Sigma1+1.0D0                    ! Baojiu-01-07-2021
              Sigma2 = Sigma2*fct4
              !
              tmp2 = (      FI2(M1u,M2 ,M3 )-      FI2(M1l,M2 ,M3 ))**2*      &
                   & (2.0D0*FI2(M1u,M2 ,M3 )+2.0D0*FI2(M1l,M2 ,M3 )    -      &
                   &        FI2(M1 ,M2u,M3 )-      FI2(M1 ,M2l,M3 )    -      &
                   &        FI2(M1 ,M2 ,M3u)-      FI2(M1 ,M2 ,M3l))/12.0D0 + &
                   & (      FI2(M1 ,M2u,M3 )-      FI2(M1 ,M2l,M3 ))**2*      &
                   & (2.0D0*FI2(M1 ,M2u,M3 )+2.0D0*FI2(M1 ,M2l,M3 )    -      &
                   &        FI2(M1u,M2 ,M3 )-      FI2(M1l,M2 ,M3 )    -      &
                   &        FI2(M1 ,M2 ,M3u)-      FI2(M1 ,M2 ,M3l))/12.0D0 + &
                   & (      FI2(M1 ,M2 ,M3u)-      FI2(M1 ,M2 ,M3l))**2*      &
                   & (2.0D0*FI2(M1 ,M2 ,M3u)+2.0D0*FI2(M1 ,M2 ,M3l)    -      &
                   &        FI2(M1u,M2 ,M3 )-      FI2(M1l,M2 ,M3 )    -      &
                   &        FI2(M1 ,M2u,M3 )-      FI2(M1 ,M2l,M3 ))/12.0D0 + &
                   & (      FI2(M1u,M2u,M3 )+      FI2(M1l,M2l,M3 )    -      &
                   &        FI2(M1u,M2l,M3 )-      FI2(M1l,M2u,M3 ))   *      &
                   & (      FI2(M1u,M2 ,M3 )-      FI2(M1l,M2 ,M3 ))   *      &
                   & (      FI2(M1 ,M2u,M3 )-      FI2(M1 ,M2l,M3 ))/8.0D0  + &
                   & (      FI2(M1u,M2 ,M3u)+      FI2(M1l,M2 ,M3l)    -      &
                   &        FI2(M1u,M2 ,M3l)-      FI2(M1l,M2 ,M3u))   *      &
                   & (      FI2(M1u,M2 ,M3 )-      FI2(M1l,M2 ,M3 ))   *      &
                   & (      FI2(M1 ,M2 ,M3u)-      FI2(M1 ,M2 ,M3l))/8.0D0  + &
                   & (      FI2(M1 ,M2u,M3u)+      FI2(M1 ,M2l,M3l)    -      &
                   &        FI2(M1 ,M2u,M3l)-      FI2(M1 ,M2l,M3u))   *      &
                   & (      FI2(M1 ,M2u,M3 )-      FI2(M1 ,M2l,M3 ))   *      &
                   & (      FI2(M1 ,M2 ,M3u)-      FI2(M1 ,M2 ,M3l))/8.0D0
              tmp2 = tmp2/dx4
              Sigma2 = Sigma2*tmp2
              !
              P = (fct2*FI(M1,M2,M3)+Sigma2)/Sigma1*dx2/ctilde
              !
              FI2(M1,M2,M3) = (Q-P)/6.0D0
           ELSE 
              Q = -FI3(M1         ,M2 +joffset,M3 +koffset2)*dx2/ctilde  + &   ! physical right-hand side  
                & (FI3(M1u+ioffset,M2 +joffset,M3 +koffset )      + &
                &  FI3(M1l+ioffset,M2 +joffset,M3 +koffset )      + &
                &  FI3(M1 +ioffset,M2u+joffset,M3 +koffset )      + &
                &  FI3(M1 +ioffset,M2l+joffset,M3 +koffset )      + &
                &  FI3(M1 +ioffset,M2 +joffset,M3u+koffset )      + &
                &  FI3(M1 +ioffset,M2 +joffset,M3l+koffset )) 
              !
              tmp1 = (FI3(M1u+ioffset,M2 +joffset,M3 +koffset) - FI3(M1l+ioffset,M2 +joffset,M3 +koffset))**2 + &
                   & (FI3(M1 +ioffset,M2u+joffset,M3 +koffset) - FI3(M1 +ioffset,M2l+joffset,M3 +koffset))**2 + &
                   & (FI3(M1 +ioffset,M2 +joffset,M3u+koffset) - FI3(M1 +ioffset,M2 +joffset,M3l+koffset))**2
              tmp1 = tmp1/(4.0D0*dx2)
              ! 
              IF(kmf_n.EQ.1) THEN
                 Sigma1 = 1.0D0 
                 Sigma2 = 0.0D0
              ENDIF
              IF(kmf_n.EQ.2) THEN
!                Sigma1 = tmp1                         ! Baojiu-01-07-2021
                 Sigma1 = fct5*tmp1+fct3*phidot2       ! Baojiu-01-07-2021
                 Sigma2 = 1.0D0
              ENDIF
              IF(kmf_n.EQ.3) THEN
!                Sigma1 = tmp1**2                      ! Baojiu-01-07-2021
                 Sigma1 = fct5*tmp1**2+fct3*phidot2**2 ! Baojiu-01-07-2021
                 Sigma2 = tmp1
              ENDIF
!             Sigma1 = Sigma1*fct3+1.0D0               ! Baojiu-01-07-2021
              Sigma1 = Sigma1+1.0D0                    ! Baojiu-01-07-2021
              Sigma2 = Sigma2*fct4
              !
              tmp2 = (      FI3(M1u+ioffset,joffset+M2 ,koffset+M3 )-      FI3(M1l+ioffset,joffset+M2 ,koffset+M3 ))**2*      &
                   & (2.0D0*FI3(M1u+ioffset,joffset+M2 ,koffset+M3 )+2.0D0*FI3(M1l+ioffset,joffset+M2 ,koffset+M3 )    -      &
                   &        FI3(M1 +ioffset,joffset+M2u,koffset+M3 )-      FI3(M1 +ioffset,joffset+M2l,koffset+M3 )    -      &
                   &        FI3(M1 +ioffset,joffset+M2 ,koffset+M3u)-      FI3(M1 +ioffset,joffset+M2 ,koffset+M3l))/12.0D0 + &
                   & (      FI3(M1 +ioffset,joffset+M2u,koffset+M3 )-      FI3(M1 +ioffset,joffset+M2l,koffset+M3 ))**2*      &
                   & (2.0D0*FI3(M1 +ioffset,joffset+M2u,koffset+M3 )+2.0D0*FI3(M1 +ioffset,joffset+M2l,koffset+M3 )    -      &
                   &        FI3(M1u+ioffset,joffset+M2 ,koffset+M3 )-      FI3(M1l+ioffset,joffset+M2 ,koffset+M3 )    -      &
                   &        FI3(M1 +ioffset,joffset+M2 ,koffset+M3u)-      FI3(M1 +ioffset,joffset+M2 ,koffset+M3l))/12.0D0 + &
                   & (      FI3(M1 +ioffset,joffset+M2 ,koffset+M3u)-      FI3(M1 +ioffset,joffset+M2 ,koffset+M3l))**2*      &
                   & (2.0D0*FI3(M1 +ioffset,joffset+M2 ,koffset+M3u)+2.0D0*FI3(M1 +ioffset,joffset+M2 ,koffset+M3l)    -      &
                   &        FI3(M1u+ioffset,joffset+M2 ,koffset+M3 )-      FI3(M1l+ioffset,joffset+M2 ,koffset+M3 )    -      &
                   &        FI3(M1 +ioffset,joffset+M2u,koffset+M3 )-      FI3(M1 +ioffset,joffset+M2l,koffset+M3 ))/12.0D0 + &
                   & (      FI3(M1u+ioffset,joffset+M2u,koffset+M3 )+      FI3(M1l+ioffset,joffset+M2l,koffset+M3 )    -      &
                   &        FI3(M1u+ioffset,joffset+M2l,koffset+M3 )-      FI3(M1l+ioffset,joffset+M2u,koffset+M3 ))   *      &
                   & (      FI3(M1u+ioffset,joffset+M2 ,koffset+M3 )-      FI3(M1l+ioffset,joffset+M2 ,koffset+M3 ))   *      &
                   & (      FI3(M1 +ioffset,joffset+M2u,koffset+M3 )-      FI3(M1 +ioffset,joffset+M2l,koffset+M3 ))/8.0D0  + &
                   & (      FI3(M1u+ioffset,joffset+M2 ,koffset+M3u)+      FI3(M1l+ioffset,joffset+M2 ,koffset+M3l)    -      &
                   &        FI3(M1u+ioffset,joffset+M2 ,koffset+M3l)-      FI3(M1l+ioffset,joffset+M2 ,koffset+M3u))   *      &
                   & (      FI3(M1u+ioffset,joffset+M2 ,koffset+M3 )-      FI3(M1l+ioffset,joffset+M2 ,koffset+M3 ))   *      &
                   & (      FI3(M1 +ioffset,joffset+M2 ,koffset+M3u)-      FI3(M1 +ioffset,joffset+M2 ,koffset+M3l))/8.0D0  + &
                   & (      FI3(M1 +ioffset,joffset+M2u,koffset+M3u)+      FI3(M1 +ioffset,joffset+M2l,koffset+M3l)    -      &
                   &        FI3(M1 +ioffset,joffset+M2u,koffset+M3l)-      FI3(M1 +ioffset,joffset+M2l,koffset+M3u))   *      &
                   & (      FI3(M1 +ioffset,joffset+M2u,koffset+M3 )-      FI3(M1 +ioffset,joffset+M2l,koffset+M3 ))   *      &
                   & (      FI3(M1 +ioffset,joffset+M2 ,koffset+M3u)-      FI3(M1 +ioffset,joffset+M2 ,koffset+M3l))/8.0D0
              tmp2 = tmp2/dx4
              Sigma2 = Sigma2*tmp2
              !
              P = (fct2*FI3(M1+ioffset,M2,M3+koffset)+Sigma2)/Sigma1*dx2/ctilde
              !
              FI3(M1+ioffset,M2+joffset,M3+koffset) = (Q-P)/6.0D0
!             write(*,*) 'aaa',P,Q,FI3(M1+ioffset,M2+joffset,M3+koffset)
!             IF(ilevel.EQ.3) write(*,*) 'bbbb',Sigma1,Sigma2
           ENDIF
        END DO
     END DO
  END DO

! IF(ilevel.EQ.3) THEN
!    OPEN(UNIT=27, FILE='test.txt', FORM='FORMATTED', STATUS='REPLACE')
!       DO M3=1,ngrid_level
!          WRITE(27,'(F20.9,F20.15,F20.15)') (DBLE(M3)-0.5d0)/DBLE(NGRID)*2.0d0*DACOS(-1.0d0),-FI3(ngrid_level/2,ngrid_level/2+joffset,M3+koffset2)*dx2/ctilde, &
!                  & FI3(ngrid_level/2+ioffset,ngrid_level/2+joffset,M3+koffset) 
!       ENDDO
!    CLOSE(27)
!    STOP
! ENDIF
 
  CALL TimingMain(3,1)
    
END SUBROUTINE relaxation_iterations_kmf

!-------------------------------------------------------------
!
! Calculate resudual and its RMS on the entire grid on ilevel.
!
!-------------------------------------------------------------
SUBROUTINE calculate_residual_kmf(ilevel,res_PM_grid)
!-------------------------------------------------------------
  use Tools
  use ExtradofBackgroundData

  integer,intent(in) :: ilevel
  real*8,intent(out) :: res_PM_grid
  
  real*8  :: ctilde,ctilde2,dx,dx2,dx4
  integer :: ngrid_level
  integer :: ioffset,joffset,koffset,koffset2
  real*8  :: fct1,fct2,fct3,fct4,fct5,tmp1,tmp2,Sigma1,Sigma2

  integer :: M1,M2,M3,M1l,M1u,M2l,M2u,M3l,M3u
  integer :: N1,N2,N3
  real*8  :: RES,OP,RES2
  integer :: i

  integer :: ist
  real*8  :: phibar,dphidN,HoH0,phidot,phidot2

  IF(MG_test) WRITE(*,'(A,I5)') 'Calculate residual on level',levelmax-ilevel

  ! number of grid points on the coarse level
  ngrid_level = NGRID/2**ilevel
  
  ctilde  = 2.99792458D3*DBLE(NGRID)/Box
  ctilde2 = ctilde**2
  dx      = DBLE(2**ilevel)
  dx2     = dx*dx
  dx4     = dx2*dx2
  !
  fct1   = 0.5D0*(1.0D0/(AEXPN*kmf_lambda))**2
  fct2   = 3.0D0*kmf_beta*Om/AEXPN
  fct3   = kmf_n*kmf_K0*fct1**(kmf_n-1)                                 ! fct3 is gamma in notes
  fct4   = -fct3*(    DBLE(kmf_n)-1.0D0)*2.0D0*ctilde
  fct5   =  fct3*(2.0*DBLE(kmf_n)+1.0D0)/3.0D0                          ! fct5 is (2*n+1)/3 *gamma in notes; Baojiu-01-07-2021
  !
  ist = 1
  DO WHILE(BackgroundEvolution(ist,1).LT.AEXPN)
     ist=ist+1
  ENDDO
  !
  ! background scalar field at AEXPN
  phibar =  BackgroundEvolution(ist-1,2)+ &
         & (BackgroundEvolution(ist,  2)-BackgroundEvolution(ist-1,2))/ &
         & (BackgroundEvolution(ist,  1)-BackgroundEvolution(ist-1,1))* &
         & (AEXPN                       -BackgroundEvolution(ist-1,1))
  ! d\varphi/dN at AEXPN, N=ln(a)
  dphidN =  BackgroundEvolution(ist-1,3)+ &
         & (BackgroundEvolution(ist,  3)-BackgroundEvolution(ist-1,3))/ &
         & (BackgroundEvolution(ist,  1)-BackgroundEvolution(ist-1,1))* &
         & (AEXPN                       -BackgroundEvolution(ist-1,1))
  ! H/H0 at AEXPN, here H=a'/a=dN/d\tau with '=d/d\tau
  HoH0   =  BackgroundEvolution(ist-1,5)+ &
         & (BackgroundEvolution(ist,  5)-BackgroundEvolution(ist-1,5))/ &
         & (BackgroundEvolution(ist,  1)-BackgroundEvolution(ist-1,1))* &
         & (AEXPN                       -BackgroundEvolution(ist-1,1))
  !
  fct2 = fct2*DEXP(kmf_beta*phibar)
  !
  phidot = dphidN*HoH0                                                  ! (d\varphi/d\tau)/H0
  phidot2 = phidot**2                                                   ! [(d\varphi/d\tau)/H0]^2
  !
  ioffset  = 2**(levelmax-ilevel)
  joffset  = 2**(levelmax-ilevel)
  koffset  = 0
  koffset2 = 0 
  IF(ilevel.GT.0) THEN 
     koffset  = 2**(levelmax-ilevel)*(2**ilevel-2)
     koffset2 = 2**(levelmax-ilevel)*(2**ilevel-1)
  END IF

!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3,M1u,M1l,M2u,M2l,M3u,M3l,i) &
!$OMP PRIVATE (tmp1,tmp2,Sigma1,Sigma2,OP)
  DO M3=1,ngrid_level  
     DO M2=1,ngrid_level  
        DO M1=1,ngrid_level 
           ! prepare the indices of the 6 neighbours on the 3-point stencil 
           M1u = M1+1; IF(M1u>ngrid_level) M1u=1              
           M1l = M1-1; IF(M1l<1          ) M1l=ngrid_level
           M2u = M2+1; IF(M2u>ngrid_level) M2u=1
           M2l = M2-1; IF(M2l<1          ) M2l=ngrid_level
           M3u = M3+1; IF(M3u>ngrid_level) M3u=1
           M3l = M3-1; IF(M3l<1          ) M3l=ngrid_level
           ! summation of u^(n+1) of all 6 neighbors on the 3-point stencil
           IF(ilevel.EQ.0) THEN
              ! calculate Laplacian of PDE
              OP = FI2(M1u,M2 ,M3 )      + &
                 & FI2(M1l,M2 ,M3 )      + &
                 & FI2(M1 ,M2u,M3 )      + &
                 & FI2(M1 ,M2l,M3 )      + &
                 & FI2(M1 ,M2 ,M3u)      + &
                 & FI2(M1 ,M2 ,M3l)      - &
                 & FI2(M1 ,M2 ,M3 )*6.0D0
              !
              tmp1 = (FI2(M1u,M2 ,M3 ) - FI2(M1l,M2 ,M3 ))**2 + &
                   & (FI2(M1 ,M2u,M3 ) - FI2(M1 ,M2l,M3 ))**2 + &
                   & (FI2(M1 ,M2 ,M3u) - FI2(M1 ,M2 ,M3l))**2
              tmp1 = tmp1/(4.0D0*dx2)
              !
              IF(kmf_n.EQ.2) THEN
!                Sigma1 = tmp1+phidot2                   ! Baojiu-01-07-2021
                 Sigma1 = fct5*tmp1+fct3*phidot2         ! Baojiu-01-07-2021
                 Sigma2 = 1.0D0
              ENDIF
              IF(kmf_n.EQ.3) THEN
!                Sigma1 = tmp1**2+phidot2**2             ! Baojiu-01-07-2021
                 Sigma1 = fct5*tmp1**2+fct3*phidot2**2   ! Baojiu-01-07-2021 
                 Sigma2 = tmp1
              ENDIF
!             Sigma1 = Sigma1*fct3+1.0D0                 ! Baojiu-01-07-2021
              Sigma1 = Sigma1+1.0D0                      ! Baojiu-01-07-2021
              Sigma2 = Sigma2*fct4
              !
              tmp2 = (      FI2(M1u,M2 ,M3 )-      FI2(M1l,M2 ,M3 ))**2*      &
                   & (2.0D0*FI2(M1u,M2 ,M3 )+2.0D0*FI2(M1l,M2 ,M3 )    -      &
                   &        FI2(M1 ,M2u,M3 )-      FI2(M1 ,M2l,M3 )    -      &
                   &        FI2(M1 ,M2 ,M3u)-      FI2(M1 ,M2 ,M3l))/12.0D0 + &
                   & (      FI2(M1 ,M2u,M3 )-      FI2(M1 ,M2l,M3 ))**2*      &
                   & (2.0D0*FI2(M1 ,M2u,M3 )+2.0D0*FI2(M1 ,M2l,M3 )    -      &
                   &        FI2(M1u,M2 ,M3 )-      FI2(M1l,M2 ,M3 )    -      &
                   &        FI2(M1 ,M2 ,M3u)-      FI2(M1 ,M2 ,M3l))/12.0D0 + &
                   & (      FI2(M1 ,M2 ,M3u)-      FI2(M1 ,M2 ,M3l))**2*      &
                   & (2.0D0*FI2(M1 ,M2 ,M3u)+2.0D0*FI2(M1 ,M2 ,M3l)    -      &
                   &        FI2(M1u,M2 ,M3 )-      FI2(M1l,M2 ,M3 )    -      &
                   &        FI2(M1 ,M2u,M3 )-      FI2(M1 ,M2l,M3 ))/12.0D0 + &
                   & (      FI2(M1u,M2u,M3 )+      FI2(M1l,M2l,M3 )    -      &
                   &        FI2(M1u,M2l,M3 )-      FI2(M1l,M2u,M3 ))   *      &
                   & (      FI2(M1u,M2 ,M3 )-      FI2(M1l,M2 ,M3 ))   *      &
                   & (      FI2(M1 ,M2u,M3 )-      FI2(M1 ,M2l,M3 ))/8.0D0  + &
                   & (      FI2(M1u,M2 ,M3u)+      FI2(M1l,M2 ,M3l)    -      &
                   &        FI2(M1u,M2 ,M3l)-      FI2(M1l,M2 ,M3u))   *      &
                   & (      FI2(M1u,M2 ,M3 )-      FI2(M1l,M2 ,M3 ))   *      &
                   & (      FI2(M1 ,M2 ,M3u)-      FI2(M1 ,M2 ,M3l))/8.0D0  + &
                   & (      FI2(M1 ,M2u,M3u)+      FI2(M1 ,M2l,M3l)    -      &
                   &        FI2(M1 ,M2u,M3l)-      FI2(M1 ,M2l,M3u))   *      &
                   & (      FI2(M1 ,M2u,M3 )-      FI2(M1 ,M2l,M3 ))   *      &
                   & (      FI2(M1 ,M2 ,M3u)-      FI2(M1 ,M2 ,M3l))/8.0D0
              tmp2 = tmp2/dx4
              Sigma2 = Sigma2*tmp2
              !
              OP = OP/dx2*ctilde
              OP = OP-(fct2*FI(M1,M2,M3)+Sigma2)/Sigma1
              FI3(M1,M2,M3) = OP
           ELSE   
              ! calculate Laplacian of PDE
              OP = FI3(M1u+ioffset,M2 +joffset,M3 +koffset)      + &
                 & FI3(M1l+ioffset,M2 +joffset,M3 +koffset)      + &
                 & FI3(M1 +ioffset,M2u+joffset,M3 +koffset)      + &
                 & FI3(M1 +ioffset,M2l+joffset,M3 +koffset)      + &
                 & FI3(M1 +ioffset,M2 +joffset,M3u+koffset)      + &
                 & FI3(M1 +ioffset,M2 +joffset,M3l+koffset)      - &
                 & FI3(M1 +ioffset,M2 +joffset,M3 +koffset)*6.0D0
              !
              tmp1 = (FI3(M1u+ioffset,M2 +joffset,M3 +koffset) - FI3(M1l+ioffset,M2 +joffset,M3 +koffset))**2 + &
                   & (FI3(M1 +ioffset,M2u+joffset,M3 +koffset) - FI3(M1 +ioffset,M2l+joffset,M3 +koffset))**2 + &
                   & (FI3(M1 +ioffset,M2 +joffset,M3u+koffset) - FI3(M1 +ioffset,M2 +joffset,M3l+koffset))**2
              tmp1 = tmp1/(4.0D0*dx2)
              !
              IF(kmf_n.EQ.2) THEN
!                Sigma1 = tmp1+phidot2                   ! Baojiu-01-07-2021
                 Sigma1 = fct5*tmp1+fct3*phidot2         ! Baojiu-01-07-2021
                 Sigma2 = 1.0D0
              ENDIF
              IF(kmf_n.EQ.3) THEN
!                Sigma1 = tmp1**2+phidot2**2             ! Baojiu-01-07-2021
                 Sigma1 = fct5*tmp1**2+fct3*phidot2**2   ! Baojiu-01-07-2021
                 Sigma2 = tmp1
              ENDIF
!             Sigma1 = Sigma1*fct3+1.0D0                 ! Baojiu-01-07-2021
              Sigma1 = Sigma1+1.0D0                      ! Baojiu-01-07-2021
              Sigma2 = Sigma2*fct4
              !
              tmp2 = (      FI3(M1u+ioffset,joffset+M2 ,koffset+M3 )-      FI3(M1l+ioffset,joffset+M2 ,koffset+M3 ))**2*      &
                   & (2.0D0*FI3(M1u+ioffset,joffset+M2 ,koffset+M3 )+2.0D0*FI3(M1l+ioffset,joffset+M2 ,koffset+M3 )    -      &
                   &        FI3(M1 +ioffset,joffset+M2u,koffset+M3 )-      FI3(M1 +ioffset,joffset+M2l,koffset+M3 )    -      &
                   &        FI3(M1 +ioffset,joffset+M2 ,koffset+M3u)-      FI3(M1 +ioffset,joffset+M2 ,koffset+M3l))/12.0D0 + &
                   & (      FI3(M1 +ioffset,joffset+M2u,koffset+M3 )-      FI3(M1 +ioffset,joffset+M2l,koffset+M3 ))**2*      &
                   & (2.0D0*FI3(M1 +ioffset,joffset+M2u,koffset+M3 )+2.0D0*FI3(M1 +ioffset,joffset+M2l,koffset+M3 )    -      &
                   &        FI3(M1u+ioffset,joffset+M2 ,koffset+M3 )-      FI3(M1l+ioffset,joffset+M2 ,koffset+M3 )    -      &
                   &        FI3(M1 +ioffset,joffset+M2 ,koffset+M3u)-      FI3(M1 +ioffset,joffset+M2 ,koffset+M3l))/12.0D0 + &
                   & (      FI3(M1 +ioffset,joffset+M2 ,koffset+M3u)-      FI3(M1 +ioffset,joffset+M2 ,koffset+M3l))**2*      &
                   & (2.0D0*FI3(M1 +ioffset,joffset+M2 ,koffset+M3u)+2.0D0*FI3(M1 +ioffset,joffset+M2 ,koffset+M3l)    -      &
                   &        FI3(M1u+ioffset,joffset+M2 ,koffset+M3 )-      FI3(M1l+ioffset,joffset+M2 ,koffset+M3 )    -      &
                   &        FI3(M1 +ioffset,joffset+M2u,koffset+M3 )-      FI3(M1 +ioffset,joffset+M2l,koffset+M3 ))/12.0D0 + &
                   & (      FI3(M1u+ioffset,joffset+M2u,koffset+M3 )+      FI3(M1l+ioffset,joffset+M2l,koffset+M3 )    -      &
                   &        FI3(M1u+ioffset,joffset+M2l,koffset+M3 )-      FI3(M1l+ioffset,joffset+M2u,koffset+M3 ))   *      &
                   & (      FI3(M1u+ioffset,joffset+M2 ,koffset+M3 )-      FI3(M1l+ioffset,joffset+M2 ,koffset+M3 ))   *      &
                   & (      FI3(M1 +ioffset,joffset+M2u,koffset+M3 )-      FI3(M1 +ioffset,joffset+M2l,koffset+M3 ))/8.0D0  + &
                   & (      FI3(M1u+ioffset,joffset+M2 ,koffset+M3u)+      FI3(M1l+ioffset,joffset+M2 ,koffset+M3l)    -      &
                   &        FI3(M1u+ioffset,joffset+M2 ,koffset+M3l)-      FI3(M1l+ioffset,joffset+M2 ,koffset+M3u))   *      &
                   & (      FI3(M1u+ioffset,joffset+M2 ,koffset+M3 )-      FI3(M1l+ioffset,joffset+M2 ,koffset+M3 ))   *      &
                   & (      FI3(M1 +ioffset,joffset+M2 ,koffset+M3u)-      FI3(M1 +ioffset,joffset+M2 ,koffset+M3l))/8.0D0  + &
                   & (      FI3(M1 +ioffset,joffset+M2u,koffset+M3u)+      FI3(M1 +ioffset,joffset+M2l,koffset+M3l)    -      &
                   &        FI3(M1 +ioffset,joffset+M2u,koffset+M3l)-      FI3(M1 +ioffset,joffset+M2l,koffset+M3u))   *      &
                   & (      FI3(M1 +ioffset,joffset+M2u,koffset+M3 )-      FI3(M1 +ioffset,joffset+M2l,koffset+M3 ))   *      &
                   & (      FI3(M1 +ioffset,joffset+M2 ,koffset+M3u)-      FI3(M1 +ioffset,joffset+M2 ,koffset+M3l))/8.0D0
              tmp2 = tmp2/dx4
              Sigma2 = Sigma2*tmp2
              !
              OP = OP/dx2*ctilde
              OP = OP-(fct2*FI3(M1+ioffset,M2,M3+koffset)+Sigma2)/Sigma1
              OP = OP-FI3(M1,M2+joffset,M3+koffset2)
              !
              FI3(M1+ioffset,M2,M3+koffset2) = OP 
              !
           END IF
        END DO
     END DO
  END DO
    
  CALL TimingMain(3,1)
  RES = 0.0D0

  IF(ilevel.EQ.0) THEN
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (N1,N2,N3) REDUCTION(+:RES)
     DO N3=1,ngrid_level
        DO N2=1,ngrid_level
           DO N1=1,ngrid_level
              RES = RES+FI3(N1,N2,N3)**2
           END DO
        END DO
     END DO
  ENDIF

  RES = DSQRT(RES/(DBLE(ngrid_level))**3)
  res_PM_grid = RES
 
! WRITE(*,'(A,F20.16,A,I5,A)') 'Residual RMS = ',res_PM_grid,' on level ',levelmax-ilevel  
! WRITE(*,'(A,F20.16,A,I5,A)') 'Residual RES2 = ',RES2,' on level ',levelmax-ilevel  
      
END SUBROUTINE calculate_residual_kmf


!-------------------------------------------------------------
!
! Calculate restricted residual field on coarse level "ilevel"
! Note: "ilevel" means the grid ilevel coarer than the PM grid
!
!-------------------------------------------------------------
SUBROUTINE restrict_residual_kmf(ilevel)
!-------------------------------------------------------------
  use Tools
  use ExtradofBackgroundData

  integer :: ilevel
  real*8  :: ctilde,ctilde2,dx,dx2,dx4
  real*8  :: fct1,fct2,fct3,fct4,fct5
  integer :: ngrid_level 
  integer :: ioffset,joffset,koffset,ioffset2,joffset2,koffset2
  integer :: i
 
  integer :: M1,M2,M3
  integer :: M1u,M1l,M2u,M2l,M3u,M3l 
  real*8  :: tmp1,tmp2,Sigma1,Sigma2,S1,S2,S3,S4,S5,S6,S7,S8,P,Q

  integer :: ist
  real*8  :: phibar,dphidN,HoH0,phidot,phidot2

  IF(MG_test) WRITE(*,'(A,I5)') 'Restrict residual to level',levelmax-ilevel
  !
  ! number of grid points on the coarse level
  ngrid_level = NGRID/2**ilevel
  !
  ctilde  = 2.99792458D3*DBLE(NGRID)/Box
  ctilde2 = ctilde**2
  dx      = 1.0d0  
  dx2     = dx*dx
  dx4     = dx2*dx2
  !
  fct1   = 0.5D0*(1.0D0/(AEXPN*kmf_lambda))**2
  fct2   = 3.0D0*kmf_beta*Om/AEXPN
  fct3   = kmf_n*kmf_K0*fct1**(kmf_n-1)                                 ! fct3 is gamma in notes
  fct4   = -fct3*(    DBLE(kmf_n)-1.0D0)*2.0D0*ctilde
  fct5   =  fct3*(2.0*DBLE(kmf_n)+1.0D0)/3.0D0                          ! fct5 is (2*n+1)/3 *gamma in notes; Baojiu-01-07-2021
  !
  ist = 1
  DO WHILE(BackgroundEvolution(ist,1).LT.AEXPN)
     ist=ist+1
  ENDDO
  !
  ! background scalar field at AEXPN
  phibar =  BackgroundEvolution(ist-1,2)+ &
         & (BackgroundEvolution(ist,  2)-BackgroundEvolution(ist-1,2))/ &
         & (BackgroundEvolution(ist,  1)-BackgroundEvolution(ist-1,1))* &
         & (AEXPN                       -BackgroundEvolution(ist-1,1))
  ! d\varphi/dN at AEXPN, N=ln(a)
  dphidN =  BackgroundEvolution(ist-1,3)+ &
         & (BackgroundEvolution(ist,  3)-BackgroundEvolution(ist-1,3))/ &
         & (BackgroundEvolution(ist,  1)-BackgroundEvolution(ist-1,1))* &
         & (AEXPN                       -BackgroundEvolution(ist-1,1))
  ! H/H0 at AEXPN, here H=a'/a=dN/d\tau with '=d/d\tau
  HoH0   =  BackgroundEvolution(ist-1,5)+ &
         & (BackgroundEvolution(ist,  5)-BackgroundEvolution(ist-1,5))/ &
         & (BackgroundEvolution(ist,  1)-BackgroundEvolution(ist-1,1))* &
         & (AEXPN                       -BackgroundEvolution(ist-1,1))
  !
  fct2 = fct2*DEXP(kmf_beta*phibar)
  !
  phidot = dphidN*HoH0                                                  ! (d\varphi/d\tau)/H0
  phidot2 = phidot**2                                                   ! [(d\varphi/d\tau)/H0]^2
  !
  IF(ilevel.EQ.1) THEN
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3,M1u,M1l,M2u,M2l,M3u,M3l,i) &
!$OMP PRIVATE (S1,S2,S3,S4,S5,S6,S7,S8) &
!$OMP PRIVATE (tmp1,tmp2,Sigma1,Sigma2,P,Q) 
     DO M3=1,ngrid_level
        DO M2=1,ngrid_level
           DO M1=1,ngrid_level
              P = 0.0d0                                         ! Laplacian terms
              Q = 0.0d0                                         ! non-Laplacian terms
              ! prepare indices of neighbour cells on the 3-point
              ! stencil and apply the periodic boundary condition
              M1u=2*M1+1; IF(M1u>NGRID) M1u=1
              M1l=2*M1-2; IF(M1l<1    ) M1l=NGRID
              M2u=2*M2+1; IF(M2u>NGRID) M2u=1
              M2l=2*M2-2; IF(M2l<1    ) M2l=NGRID
              M3u=2*M3+1; IF(M3u>NGRID) M3u=1
              M3l=2*M3-2; IF(M3l<1    ) M3l=NGRID
              ! accumulate contributions to the non-Laplacian term 
              ! of PDE from all 8 son cells
              !
              ! Source term of cell #1 (2*M1-1,2*M2-1,2*M3-1)
              tmp1 = (FI2(2*M1  ,2*M2-1,2*M3-1) - FI2(  M1l ,2*M2-1,2*M3-1))**2 + &
                   & (FI2(2*M1-1,2*M2  ,2*M3-1) - FI2(2*M1-1,  M2l ,2*M3-1))**2 + &
                   & (FI2(2*M1-1,2*M2-1,2*M3  ) - FI2(2*M1-1,2*M2-1,  M3l ))**2
              tmp1 = tmp1/(4.0D0*dx2)
              !
              IF(kmf_n.EQ.2) THEN
!                Sigma1 = tmp1+phidot2                   ! Baojiu-01-07-2021
                 Sigma1 = fct5*tmp1+fct3*phidot2         ! Baojiu-01-07-2021
                 Sigma2 = 1.0D0
              ENDIF
              IF(kmf_n.EQ.3) THEN
!                Sigma1 = tmp1**2+phidot2**2             ! Baojiu-01-07-2021
                 Sigma1 = fct5*tmp1**2+fct3*phidot2**2   ! Baojiu-01-07-2021
                 Sigma2 = tmp1
              ENDIF
!             Sigma1 = Sigma1*fct3+1.0D0                 ! Baojiu-01-07-2021
              Sigma1 = Sigma1+1.0D0                      ! Baojiu-01-07-2021
              Sigma2 = Sigma2*fct4
              !
              tmp2 = (      FI2(2*M1  ,2*M2-1,2*M3-1)-      FI2(  M1l ,2*M2-1,2*M3-1))**2*      &
                   & (2.0D0*FI2(2*M1  ,2*M2-1,2*M3-1)+2.0D0*FI2(  M1l ,2*M2-1,2*M3-1)    -      &
                   &        FI2(2*M1-1,2*M2  ,2*M3-1)-      FI2(2*M1-1,  M2l ,2*M3-1)    -      &
                   &        FI2(2*M1-1,2*M2-1,2*M3  )-      FI2(2*M1-1,2*M2-1,  M3l ))/12.0D0 + &
                   & (      FI2(2*M1-1,2*M2  ,2*M3-1)-      FI2(2*M1-1,  M2l ,2*M3-1))**2*      &
                   & (2.0D0*FI2(2*M1-1,2*M2  ,2*M3-1)+2.0D0*FI2(2*M1-1,  M2l ,2*M3-1)    -      &
                   &        FI2(2*M1  ,2*M2-1,2*M3-1)-      FI2(  M1l ,2*M2-1,2*M3-1)    -      &
                   &        FI2(2*M1-1,2*M2-1,2*M3  )-      FI2(2*M1-1,2*M2-1,  M3l ))/12.0D0 + &
                   & (      FI2(2*M1-1,2*M2-1,2*M3  )-      FI2(2*M1-1,2*M2-1,  M3l ))**2*      &
                   & (2.0D0*FI2(2*M1-1,2*M2-1,2*M3  )+2.0D0*FI2(2*M1-1,2*M2-1,  M3l )    -      &
                   &        FI2(2*M1  ,2*M2-1,2*M3-1)-      FI2(  M1l ,2*M2-1,2*M3-1)    -      &
                   &        FI2(2*M1-1,2*M2  ,2*M3-1)-      FI2(2*M1-1,  M2l ,2*M3-1))/12.0D0 + &
                   & (      FI2(2*M1  ,2*M2  ,2*M3-1)+      FI2(  M1l ,  M2l ,2*M3-1)    -      &
                   &        FI2(2*M1  ,  M2l ,2*M3-1)-      FI2(  M1l ,2*M2  ,2*M3-1))   *      &
                   & (      FI2(2*M1  ,2*M2-1,2*M3-1)-      FI2(  M1l ,2*M2-1,2*M3-1))   *      &
                   & (      FI2(2*M1-1,2*M2  ,2*M3-1)-      FI2(2*M1-1,  M2l ,2*M3-1))/8.0D0  + &
                   & (      FI2(2*M1  ,2*M2-1,2*M3  )+      FI2(  M1l ,2*M2-1,  M3l )    -      &
                   &        FI2(2*M1  ,2*M2-1,  M3l )-      FI2(  M1l ,2*M2-1,2*M3  ))   *      &
                   & (      FI2(2*M1  ,2*M2-1,2*M3-1)-      FI2(  M1l ,2*M2-1,2*M3-1))   *      &
                   & (      FI2(2*M1-1,2*M2-1,2*M3  )-      FI2(2*M1-1,2*M2-1,  M3l ))/8.0D0  + &
                   & (      FI2(2*M1-1,2*M2  ,2*M3  )+      FI2(2*M1-1,  M2l ,  M3l )    -      &
                   &        FI2(2*M1-1,2*M2  ,  M3l )-      FI2(2*M1-1,  M2l ,2*M3  ))   *      &
                   & (      FI2(2*M1-1,2*M2  ,2*M3-1)-      FI2(2*M1-1,  M2l ,2*M3-1))   *      &
                   & (      FI2(2*M1-1,2*M2-1,2*M3  )-      FI2(2*M1-1,2*M2-1,  M3l ))/8.0D0
              tmp2 = tmp2/dx4
              Sigma2 = Sigma2*tmp2
              S1   = (fct2*FI(2*M1-1,2*M2-1,2*M3-1)+Sigma2)/Sigma1
              !
              ! Source term of cell #2 (2*M1  ,2*M2-1,2*M3-1)
              tmp1 = (FI2(  M1u ,2*M2-1,2*M3-1) - FI2(2*M1-1,2*M2-1,2*M3-1))**2 + &
                   & (FI2(2*M1  ,2*M2  ,2*M3-1) - FI2(2*M1  ,  M2l ,2*M3-1))**2 + &
                   & (FI2(2*M1  ,2*M2-1,2*M3  ) - FI2(2*M1  ,2*M2-1,  M3l ))**2
              tmp1 = tmp1/(4.0D0*dx2)
              !
              IF(kmf_n.EQ.2) THEN
!                Sigma1 = tmp1+phidot2                   ! Baojiu-01-07-2021
                 Sigma1 = fct5*tmp1+fct3*phidot2         ! Baojiu-01-07-2021
                 Sigma2 = 1.0D0
              ENDIF
              IF(kmf_n.EQ.3) THEN
!                Sigma1 = tmp1**2+phidot2**2             ! Baojiu-01-07-2021
                 Sigma1 = fct5*tmp1**2+fct3*phidot2**2   ! Baojiu-01-07-2021
                 Sigma2 = tmp1
              ENDIF
!             Sigma1 = Sigma1*fct3+1.0D0                 ! Baojiu-01-07-2021
              Sigma1 = Sigma1+1.0D0                      ! Baojiu-01-07-2021
              Sigma2 = Sigma2*fct4
              !
              tmp2 = (      FI2(  M1u ,2*M2-1,2*M3-1)-      FI2(2*M1-1,2*M2-1,2*M3-1))**2*      &
                   & (2.0D0*FI2(  M1u ,2*M2-1,2*M3-1)+2.0D0*FI2(2*M1-1,2*M2-1,2*M3-1)    -      &
                   &        FI2(2*M1  ,2*M2  ,2*M3-1)-      FI2(2*M1  ,  M2l ,2*M3-1)    -      &
                   &        FI2(2*M1  ,2*M2-1,2*M3  )-      FI2(2*M1  ,2*M2-1,  M3l ))/12.0D0 + &
                   & (      FI2(2*M1  ,2*M2  ,2*M3-1)-      FI2(2*M1  ,  M2l ,2*M3-1))**2*      &
                   & (2.0D0*FI2(2*M1  ,2*M2  ,2*M3-1)+2.0D0*FI2(2*M1  ,  M2l ,2*M3-1)    -      &
                   &        FI2(  M1u ,2*M2-1,2*M3-1)-      FI2(2*M1-1,2*M2-1,2*M3-1)    -      &
                   &        FI2(2*M1  ,2*M2-1,2*M3  )-      FI2(2*M1  ,2*M2-1,  M3l ))/12.0D0 + &
                   & (      FI2(2*M1  ,2*M2-1,2*M3  )-      FI2(2*M1  ,2*M2-1,  M3l ))**2*      &
                   & (2.0D0*FI2(2*M1  ,2*M2-1,2*M3  )+2.0D0*FI2(2*M1  ,2*M2-1,  M3l )    -      &
                   &        FI2(  M1u ,2*M2-1,2*M3-1)-      FI2(2*M1-1,2*M2-1,2*M3-1)    -      &
                   &        FI2(2*M1  ,2*M2  ,2*M3-1)-      FI2(2*M1  ,  M2l ,2*M3-1))/12.0D0 + &
                   & (      FI2(  M1u ,2*M2  ,2*M3-1)+      FI2(2*M1-1,  M2l ,2*M3-1)    -      &
                   &        FI2(  M1u ,  M2l ,2*M3-1)-      FI2(2*M1-1,2*M2  ,2*M3-1))   *      &
                   & (      FI2(  M1u ,2*M2-1,2*M3-1)-      FI2(2*M1-1,2*M2-1,2*M3-1))   *      &
                   & (      FI2(2*M1  ,2*M2  ,2*M3-1)-      FI2(2*M1  ,  M2l ,2*M3-1))/8.0D0  + &
                   & (      FI2(  M1u ,2*M2-1,2*M3  )+      FI2(2*M1-1,2*M2-1,  M3l )    -      &
                   &        FI2(  M1u ,2*M2-1,  M3l )-      FI2(2*M1-1,2*M2-1,2*M3  ))   *      &
                   & (      FI2(  M1u ,2*M2-1,2*M3-1)-      FI2(2*M1-1,2*M2-1,2*M3-1))   *      &
                   & (      FI2(2*M1  ,2*M2-1,2*M3  )-      FI2(2*M1  ,2*M2-1,  M3l ))/8.0D0  + &
                   & (      FI2(2*M1  ,2*M2  ,2*M3  )+      FI2(2*M1  ,  M2l ,  M3l )    -      &
                   &        FI2(2*M1  ,2*M2  ,  M3l )-      FI2(2*M1  ,  M2l ,2*M3  ))   *      &
                   & (      FI2(2*M1  ,2*M2  ,2*M3-1)-      FI2(2*M1  ,  M2l ,2*M3-1))   *      &
                   & (      FI2(2*M1  ,2*M2-1,2*M3  )-      FI2(2*M1  ,2*M2-1,  M3l ))/8.0D0
              tmp2 = tmp2/dx4
              Sigma2 = Sigma2*tmp2
              S2   = (fct2*FI(2*M1  ,2*M2-1,2*M3-1)+Sigma2)/Sigma1
              !
              ! Source term of cell #3 (2*M1-1,2*M2  ,2*M3-1)
              tmp1 = (FI2(2*M1  ,2*M2  ,2*M3-1) - FI2(  M1l ,2*M2  ,2*M3-1))**2 + &
                   & (FI2(2*M1-1,  M2u ,2*M3-1) - FI2(2*M1-1,2*M2-1,2*M3-1))**2 + &
                   & (FI2(2*M1-1,2*M2  ,2*M3  ) - FI2(2*M1-1,2*M2  ,  M3l ))**2
              tmp1 = tmp1/(4.0D0*dx2)
              !
              IF(kmf_n.EQ.2) THEN
!                Sigma1 = tmp1+phidot2                   ! Baojiu-01-07-2021
                 Sigma1 = fct5*tmp1+fct3*phidot2         ! Baojiu-01-07-2021
                 Sigma2 = 1.0D0
              ENDIF
              IF(kmf_n.EQ.3) THEN
!                Sigma1 = tmp1**2+phidot2**2             ! Baojiu-01-07-2021
                 Sigma1 = fct5*tmp1**2+fct3*phidot2**2   ! Baojiu-01-07-2021
                 Sigma2 = tmp1
              ENDIF
!             Sigma1 = Sigma1*fct3+1.0D0                 ! Baojiu-01-07-2021
              Sigma1 = Sigma1+1.0D0                      ! Baojiu-01-07-2021
              Sigma2 = Sigma2*fct4
              !
              tmp2 = (      FI2(2*M1  ,2*M2  ,2*M3-1)-      FI2(  M1l ,2*M2  ,2*M3-1))**2*      &
                   & (2.0D0*FI2(2*M1  ,2*M2  ,2*M3-1)+2.0D0*FI2(  M1l ,2*M2  ,2*M3-1)    -      &
                   &        FI2(2*M1-1,  M2u ,2*M3-1)-      FI2(2*M1-1,2*M2-1,2*M3-1)    -      &
                   &        FI2(2*M1-1,2*M2  ,2*M3  )-      FI2(2*M1-1,2*M2  ,  M3l ))/12.0D0 + &
                   & (      FI2(2*M1-1,  M2u ,2*M3-1)-      FI2(2*M1-1,2*M2-1,2*M3-1))**2*      &
                   & (2.0D0*FI2(2*M1-1,  M2u ,2*M3-1)+2.0D0*FI2(2*M1-1,2*M2-1,2*M3-1)    -      &
                   &        FI2(2*M1  ,2*M2  ,2*M3-1)-      FI2(  M1l ,2*M2  ,2*M3-1)    -      &
                   &        FI2(2*M1-1,2*M2  ,2*M3  )-      FI2(2*M1-1,2*M2  ,  M3l ))/12.0D0 + &
                   & (      FI2(2*M1-1,2*M2  ,2*M3  )-      FI2(2*M1-1,2*M2  ,  M3l ))**2*      &
                   & (2.0D0*FI2(2*M1-1,2*M2  ,2*M3  )+2.0D0*FI2(2*M1-1,2*M2  ,  M3l )    -      &
                   &        FI2(2*M1  ,2*M2  ,2*M3-1)-      FI2(  M1l ,2*M2  ,2*M3-1)    -      &
                   &        FI2(2*M1-1,  M2u ,2*M3-1)-      FI2(2*M1-1,2*M2-1,2*M3-1))/12.0D0 + &
                   & (      FI2(2*M1  ,  M2u ,2*M3-1)+      FI2(  M1l ,2*M2-1,2*M3-1)    -      &
                   &        FI2(2*M1  ,2*M2-1,2*M3-1)-      FI2(  M1l ,  M2u ,2*M3-1))   *      &
                   & (      FI2(2*M1  ,2*M2  ,2*M3-1)-      FI2(  M1l ,2*M2  ,2*M3-1))   *      &
                   & (      FI2(2*M1-1,  M2u ,2*M3-1)-      FI2(2*M1-1,2*M2-1,2*M3-1))/8.0D0  + &
                   & (      FI2(2*M1  ,2*M2  ,2*M3  )+      FI2(  M1l ,2*M2  ,  M3l )    -      &
                   &        FI2(2*M1  ,2*M2  ,  M3l )-      FI2(  M1l ,2*M2  ,2*M3  ))   *      &
                   & (      FI2(2*M1  ,2*M2  ,2*M3-1)-      FI2(  M1l ,2*M2  ,2*M3-1))   *      &
                   & (      FI2(2*M1-1,2*M2  ,2*M3  )-      FI2(2*M1-1,2*M2  ,  M3l ))/8.0D0  + &
                   & (      FI2(2*M1-1,  M2u ,2*M3  )+      FI2(2*M1-1,2*M2-1,  M3l )    -      &
                   &        FI2(2*M1-1,  M2u ,  M3l )-      FI2(2*M1-1,2*M2-1,2*M3  ))   *      &
                   & (      FI2(2*M1-1,  M2u ,2*M3-1)-      FI2(2*M1-1,2*M2-1,2*M3-1))   *      &
                   & (      FI2(2*M1-1,2*M2  ,2*M3  )-      FI2(2*M1-1,2*M2  ,  M3l ))/8.0D0
              tmp2 = tmp2/dx4
              Sigma2 = Sigma2*tmp2
              S3   = (fct2*FI(2*M1-1,2*M2  ,2*M3-1)+Sigma2)/Sigma1
              !
              ! Source term of cell #4 (2*M1  ,2*M2  ,2*M3-1)
              tmp1 = (FI2(  M1u ,2*M2  ,2*M3-1) - FI2(2*M1-1,2*M2  ,2*M3-1))**2 + &
                   & (FI2(2*M1  ,  M2u ,2*M3-1) - FI2(2*M1  ,2*M2-1,2*M3-1))**2 + &
                   & (FI2(2*M1  ,2*M2  ,2*M3  ) - FI2(2*M1  ,2*M2  ,  M3l ))**2
              tmp1 = tmp1/(4.0D0*dx2)
              !
              IF(kmf_n.EQ.2) THEN
!                Sigma1 = tmp1+phidot2                   ! Baojiu-01-07-2021
                 Sigma1 = fct5*tmp1+fct3*phidot2         ! Baojiu-01-07-2021
                 Sigma2 = 1.0D0
              ENDIF
              IF(kmf_n.EQ.3) THEN
!                Sigma1 = tmp1**2+phidot2**2             ! Baojiu-01-07-2021
                 Sigma1 = fct5*tmp1**2+fct3*phidot2**2   ! Baojiu-01-07-2021
                 Sigma2 = tmp1
              ENDIF
!             Sigma1 = Sigma1*fct3+1.0D0                 ! Baojiu-01-07-2021
              Sigma1 = Sigma1+1.0D0                      ! Baojiu-01-07-2021
              Sigma2 = Sigma2*fct4
              !
              tmp2 = (      FI2(  M1u ,2*M2  ,2*M3-1)-      FI2(2*M1-1,2*M2  ,2*M3-1))**2*      &
                   & (2.0D0*FI2(  M1u ,2*M2  ,2*M3-1)+2.0D0*FI2(2*M1-1,2*M2  ,2*M3-1)    -      &
                   &        FI2(2*M1  ,  M2u ,2*M3-1)-      FI2(2*M1  ,2*M2-1,2*M3-1)    -      &
                   &        FI2(2*M1  ,2*M2  ,2*M3  )-      FI2(2*M1  ,2*M2  ,  M3l ))/12.0D0 + &
                   & (      FI2(2*M1  ,  M2u ,2*M3-1)-      FI2(2*M1  ,2*M2-1,2*M3-1))**2*      &
                   & (2.0D0*FI2(2*M1  ,  M2u ,2*M3-1)+2.0D0*FI2(2*M1  ,2*M2-1,2*M3-1)    -      &
                   &        FI2(  M1u ,2*M2  ,2*M3-1)-      FI2(2*M1-1,2*M2  ,2*M3-1)    -      &
                   &        FI2(2*M1  ,2*M2  ,2*M3  )-      FI2(2*M1  ,2*M2  ,  M3l ))/12.0D0 + &
                   & (      FI2(2*M1  ,2*M2  ,2*M3  )-      FI2(2*M1  ,2*M2  ,  M3l ))**2*      &
                   & (2.0D0*FI2(2*M1  ,2*M2  ,2*M3  )+2.0D0*FI2(2*M1  ,2*M2  ,  M3l )    -      &
                   &        FI2(  M1u ,2*M2  ,2*M3-1)-      FI2(2*M1-1,2*M2  ,2*M3-1)    -      &
                   &        FI2(2*M1  ,  M2u ,2*M3-1)-      FI2(2*M1  ,2*M2-1,2*M3-1))/12.0D0 + &
                   & (      FI2(  M1u ,  M2u ,2*M3-1)+      FI2(2*M1-1,2*M2-1,2*M3-1)    -      &
                   &        FI2(  M1u ,2*M2-1,2*M3-1)-      FI2(2*M1-1,  M2u ,2*M3-1))   *      &
                   & (      FI2(  M1u ,2*M2  ,2*M3-1)-      FI2(2*M1-1,2*M2  ,2*M3-1))   *      &
                   & (      FI2(2*M1  ,  M2u ,2*M3-1)-      FI2(2*M1  ,2*M2-1,2*M3-1))/8.0D0  + &
                   & (      FI2(  M1u ,2*M2  ,2*M3  )+      FI2(2*M1-1,2*M2  ,  M3l )    -      &
                   &        FI2(  M1u ,2*M2  ,  M3l )-      FI2(2*M1-1,2*M2  ,2*M3  ))   *      &
                   & (      FI2(  M1u ,2*M2  ,2*M3-1)-      FI2(2*M1-1,2*M2  ,2*M3-1))   *      &
                   & (      FI2(2*M1  ,2*M2  ,2*M3  )-      FI2(2*M1  ,2*M2  ,  M3l ))/8.0D0  + &
                   & (      FI2(2*M1  ,  M2u ,2*M3  )+      FI2(2*M1  ,2*M2-1,  M3l )    -      &
                   &        FI2(2*M1  ,  M2u ,  M3l )-      FI2(2*M1  ,2*M2-1,2*M3  ))   *      &
                   & (      FI2(2*M1  ,  M2u ,2*M3-1)-      FI2(2*M1  ,2*M2-1,2*M3-1))   *      &
                   & (      FI2(2*M1  ,2*M2  ,2*M3  )-      FI2(2*M1  ,2*M2  ,  M3l ))/8.0D0
              tmp2 = tmp2/dx4
              Sigma2 = Sigma2*tmp2
              S4   = (fct2*FI(2*M1  ,2*M2  ,2*M3-1)+Sigma2)/Sigma1
              !
              ! Source term of cell #5 (2*M1-1,2*M2-1,2*M3  )
              tmp1 = (FI2(2*M1  ,2*M2-1,2*M3  ) - FI2(  M1l ,2*M2-1,2*M3  ))**2 + &
                   & (FI2(2*M1-1,2*M2  ,2*M3  ) - FI2(2*M1-1,  M2l ,2*M3  ))**2 + &
                   & (FI2(2*M1-1,2*M2-1,  M3u ) - FI2(2*M1-1,2*M2-1,2*M3-1))**2
              tmp1 = tmp1/(4.0D0*dx2)
              !
              IF(kmf_n.EQ.2) THEN
!                Sigma1 = tmp1+phidot2                   ! Baojiu-01-07-2021
                 Sigma1 = fct5*tmp1+fct3*phidot2         ! Baojiu-01-07-2021
                 Sigma2 = 1.0D0
              ENDIF
              IF(kmf_n.EQ.3) THEN
!                Sigma1 = tmp1**2+phidot2**2             ! Baojiu-01-07-2021
                 Sigma1 = fct5*tmp1**2+fct3*phidot2**2   ! Baojiu-01-07-2021
                 Sigma2 = tmp1
              ENDIF
!             Sigma1 = Sigma1*fct3+1.0D0                 ! Baojiu-01-07-2021
              Sigma1 = Sigma1+1.0D0                      ! Baojiu-01-07-2021
              Sigma2 = Sigma2*fct4
              !
              tmp2 = (      FI2(2*M1  ,2*M2-1,2*M3  )-      FI2(  M1l ,2*M2-1,2*M3  ))**2*      &
                   & (2.0D0*FI2(2*M1  ,2*M2-1,2*M3  )+2.0D0*FI2(  M1l ,2*M2-1,2*M3  )    -      &
                   &        FI2(2*M1-1,2*M2  ,2*M3  )-      FI2(2*M1-1,  M2l ,2*M3  )    -      &
                   &        FI2(2*M1-1,2*M2-1,  M3u )-      FI2(2*M1-1,2*M2-1,2*M3-1))/12.0D0 + &
                   & (      FI2(2*M1-1,2*M2  ,2*M3  )-      FI2(2*M1-1,  M2l ,2*M3  ))**2*      &
                   & (2.0D0*FI2(2*M1-1,2*M2  ,2*M3  )+2.0D0*FI2(2*M1-1,  M2l ,2*M3  )    -      &
                   &        FI2(2*M1  ,2*M2-1,2*M3  )-      FI2(  M1l ,2*M2-1,2*M3  )    -      &
                   &        FI2(2*M1-1,2*M2-1,  M3u )-      FI2(2*M1-1,2*M2-1,2*M3-1))/12.0D0 + &
                   & (      FI2(2*M1-1,2*M2-1,  M3u )-      FI2(2*M1-1,2*M2-1,2*M3-1))**2*      &
                   & (2.0D0*FI2(2*M1-1,2*M2-1,  M3u )+2.0D0*FI2(2*M1-1,2*M2-1,2*M3-1)    -      &
                   &        FI2(2*M1  ,2*M2-1,2*M3  )-      FI2(  M1l ,2*M2-1,2*M3  )    -      &
                   &        FI2(2*M1-1,2*M2  ,2*M3  )-      FI2(2*M1-1,  M2l ,2*M3  ))/12.0D0 + &
                   & (      FI2(2*M1  ,2*M2  ,2*M3  )+      FI2(  M1l ,  M2l ,2*M3  )    -      &
                   &        FI2(2*M1  ,  M2l ,2*M3  )-      FI2(  M1l ,2*M2  ,2*M3  ))   *      &
                   & (      FI2(2*M1  ,2*M2-1,2*M3  )-      FI2(  M1l ,2*M2-1,2*M3  ))   *      &
                   & (      FI2(2*M1-1,2*M2  ,2*M3  )-      FI2(2*M1-1,  M2l ,2*M3  ))/8.0D0  + &
                   & (      FI2(2*M1  ,2*M2-1,  M3u )+      FI2(  M1l ,2*M2-1,2*M3-1)    -      &
                   &        FI2(2*M1  ,2*M2-1,2*M3-1)-      FI2(  M1l ,2*M2-1,  M3u ))   *      &
                   & (      FI2(2*M1  ,2*M2-1,2*M3  )-      FI2(  M1l ,2*M2-1,2*M3  ))   *      &
                   & (      FI2(2*M1-1,2*M2-1,  M3u )-      FI2(2*M1-1,2*M2-1,2*M3-1))/8.0D0  + &
                   & (      FI2(2*M1-1,2*M2  ,  M3u )+      FI2(2*M1-1,  M2l ,2*M3-1)    -      &
                   &        FI2(2*M1-1,2*M2  ,2*M3-1)-      FI2(2*M1-1,  M2l ,  M3u ))   *      &
                   & (      FI2(2*M1-1,2*M2  ,2*M3  )-      FI2(2*M1-1,  M2l ,2*M3  ))   *      &
                   & (      FI2(2*M1-1,2*M2-1,  M3u )-      FI2(2*M1-1,2*M2-1,2*M3-1))/8.0D0
              tmp2 = tmp2/dx4
              Sigma2 = Sigma2*tmp2
              S5   = (fct2*FI(2*M1-1,2*M2-1,2*M3  )+Sigma2)/Sigma1
              !
              ! Source term of cell #6 (2*M1  ,2*M2-1,2*M3  )
              tmp1 = (FI2(  M1u ,2*M2-1,2*M3  ) - FI2(2*M1-1,2*M2-1,2*M3  ))**2 + &
                   & (FI2(2*M1  ,2*M2  ,2*M3  ) - FI2(2*M1  ,  M2l ,2*M3  ))**2 + &
                   & (FI2(2*M1  ,2*M2-1,  M3u ) - FI2(2*M1  ,2*M2-1,2*M3-1))**2
              tmp1 = tmp1/(4.0D0*dx2)
              !
              IF(kmf_n.EQ.2) THEN
!                Sigma1 = tmp1+phidot2                   ! Baojiu-01-07-2021
                 Sigma1 = fct5*tmp1+fct3*phidot2         ! Baojiu-01-07-2021
                 Sigma2 = 1.0D0
              ENDIF
              IF(kmf_n.EQ.3) THEN
!                Sigma1 = tmp1**2+phidot2**2             ! Baojiu-01-07-2021
                 Sigma1 = fct5*tmp1**2+fct3*phidot2**2   ! Baojiu-01-07-2021
                 Sigma2 = tmp1
              ENDIF
!             Sigma1 = Sigma1*fct3+1.0D0                 ! Baojiu-01-07-2021
              Sigma1 = Sigma1+1.0D0                      ! Baojiu-01-07-2021
              Sigma2 = Sigma2*fct4
              !
              tmp2 = (      FI2(  M1u ,2*M2-1,2*M3  )-      FI2(2*M1-1,2*M2-1,2*M3  ))**2*      &
                   & (2.0D0*FI2(  M1u ,2*M2-1,2*M3  )+2.0D0*FI2(2*M1-1,2*M2-1,2*M3  )    -      &
                   &        FI2(2*M1  ,2*M2  ,2*M3  )-      FI2(2*M1  ,  M2l ,2*M3  )    -      &
                   &        FI2(2*M1  ,2*M2-1,  M3u )-      FI2(2*M1  ,2*M2-1,2*M3-1))/12.0D0 + &
                   & (      FI2(2*M1  ,2*M2  ,2*M3  )-      FI2(2*M1  ,  M2l ,2*M3  ))**2*      &
                   & (2.0D0*FI2(2*M1  ,2*M2  ,2*M3  )+2.0D0*FI2(2*M1  ,  M2l ,2*M3  )    -      &
                   &        FI2(  M1u ,2*M2-1,2*M3  )-      FI2(2*M1-1,2*M2-1,2*M3  )    -      &
                   &        FI2(2*M1  ,2*M2-1,  M3u )-      FI2(2*M1  ,2*M2-1,2*M3-1))/12.0D0 + &
                   & (      FI2(2*M1  ,2*M2-1,  M3u )-      FI2(2*M1  ,2*M2-1,2*M3-1))**2*      &
                   & (2.0D0*FI2(2*M1  ,2*M2-1,  M3u )+2.0D0*FI2(2*M1  ,2*M2-1,2*M3-1)    -      &
                   &        FI2(  M1u ,2*M2-1,2*M3  )-      FI2(2*M1-1,2*M2-1,2*M3  )    -      &
                   &        FI2(2*M1  ,2*M2  ,2*M3  )-      FI2(2*M1  ,  M2l ,2*M3  ))/12.0D0 + &
                   & (      FI2(  M1u ,2*M2  ,2*M3  )+      FI2(2*M1-1,  M2l ,2*M3  )    -      &
                   &        FI2(  M1u ,  M2l ,2*M3  )-      FI2(2*M1-1,2*M2  ,2*M3  ))   *      &
                   & (      FI2(  M1u ,2*M2-1,2*M3  )-      FI2(2*M1-1,2*M2-1,2*M3  ))   *      &
                   & (      FI2(2*M1  ,2*M2  ,2*M3  )-      FI2(2*M1  ,  M2l ,2*M3  ))/8.0D0  + &
                   & (      FI2(  M1u ,2*M2-1,  M3u )+      FI2(2*M1-1,2*M2-1,2*M3-1)    -      &
                   &        FI2(  M1u ,2*M2-1,2*M3-1)-      FI2(2*M1-1,2*M2-1,  M3u ))   *      &
                   & (      FI2(  M1u ,2*M2-1,2*M3  )-      FI2(2*M1-1,2*M2-1,2*M3  ))   *      &
                   & (      FI2(2*M1  ,2*M2-1,  M3u )-      FI2(2*M1  ,2*M2-1,2*M3-1))/8.0D0  + &
                   & (      FI2(2*M1  ,2*M2  ,  M3u )+      FI2(2*M1  ,  M2l ,2*M3-1)    -      &
                   &        FI2(2*M1  ,2*M2  ,2*M3-1)-      FI2(2*M1  ,  M2l ,  M3u ))   *      &
                   & (      FI2(2*M1  ,2*M2  ,2*M3  )-      FI2(2*M1  ,  M2l ,2*M3  ))   *      &
                   & (      FI2(2*M1  ,2*M2-1,  M3u )-      FI2(2*M1  ,2*M2-1,2*M3-1))/8.0D0
              tmp2 = tmp2/dx4
              Sigma2 = Sigma2*tmp2
              S6   = (fct2*FI(2*M1  ,2*M2-1,2*M3  )+Sigma2)/Sigma1
              !
              ! Source term of cell #7 (2*M1-1,2*M2  ,2*M3  )
              tmp1 = (FI2(2*M1  ,2*M2  ,2*M3  ) - FI2(  M1l ,2*M2  ,2*M3  ))**2 + &
                   & (FI2(2*M1-1,  M2u ,2*M3  ) - FI2(2*M1-1,2*M2-1,2*M3  ))**2 + &
                   & (FI2(2*M1-1,2*M2  ,  M3u ) - FI2(2*M1-1,2*M2  ,2*M3-1))**2
              tmp1 = tmp1/(4.0D0*dx2)
              !
              IF(kmf_n.EQ.2) THEN
!                Sigma1 = tmp1+phidot2                   ! Baojiu-01-07-2021
                 Sigma1 = fct5*tmp1+fct3*phidot2         ! Baojiu-01-07-2021
                 Sigma2 = 1.0D0
              ENDIF
              IF(kmf_n.EQ.3) THEN
!                Sigma1 = tmp1**2+phidot2**2             ! Baojiu-01-07-2021
                 Sigma1 = fct5*tmp1**2+fct3*phidot2**2   ! Baojiu-01-07-2021
                 Sigma2 = tmp1
              ENDIF
!             Sigma1 = Sigma1*fct3+1.0D0                 ! Baojiu-01-07-2021
              Sigma1 = Sigma1+1.0D0                      ! Baojiu-01-07-2021
              Sigma2 = Sigma2*fct4
              !
              tmp2 = (      FI2(2*M1  ,2*M2  ,2*M3  )-      FI2(  M1l ,2*M2  ,2*M3  ))**2*      &
                   & (2.0D0*FI2(2*M1  ,2*M2  ,2*M3  )+2.0D0*FI2(  M1l ,2*M2  ,2*M3  )    -      &
                   &        FI2(2*M1-1,  M2u ,2*M3  )-      FI2(2*M1-1,2*M2-1,2*M3  )    -      &
                   &        FI2(2*M1-1,2*M2  ,  M3u )-      FI2(2*M1-1,2*M2  ,2*M3-1))/12.0D0 + &
                   & (      FI2(2*M1-1,  M2u ,2*M3  )-      FI2(2*M1-1,2*M2-1,2*M3  ))**2*      &
                   & (2.0D0*FI2(2*M1-1,  M2u ,2*M3  )+2.0D0*FI2(2*M1-1,2*M2-1,2*M3  )    -      &
                   &        FI2(2*M1  ,2*M2  ,2*M3  )-      FI2(  M1l ,2*M2  ,2*M3  )    -      &
                   &        FI2(2*M1-1,2*M2  ,  M3u )-      FI2(2*M1-1,2*M2  ,2*M3-1))/12.0D0 + &
                   & (      FI2(2*M1-1,2*M2  ,  M3u )-      FI2(2*M1-1,2*M2  ,2*M3-1))**2*      &
                   & (2.0D0*FI2(2*M1-1,2*M2  ,  M3u )+2.0D0*FI2(2*M1-1,2*M2  ,2*M3-1)    -      &
                   &        FI2(2*M1  ,2*M2  ,2*M3  )-      FI2(  M1l ,2*M2  ,2*M3  )    -      &
                   &        FI2(2*M1-1,  M2u ,2*M3  )-      FI2(2*M1-1,2*M2-1,2*M3  ))/12.0D0 + &
                   & (      FI2(2*M1  ,  M2u ,2*M3  )+      FI2(  M1l ,2*M2-1,2*M3  )    -      &
                   &        FI2(2*M1  ,2*M2-1,2*M3  )-      FI2(  M1l ,  M2u ,2*M3  ))   *      &
                   & (      FI2(2*M1  ,2*M2  ,2*M3  )-      FI2(  M1l ,2*M2  ,2*M3  ))   *      &
                   & (      FI2(2*M1-1,  M2u ,2*M3  )-      FI2(2*M1-1,2*M2-1,2*M3  ))/8.0D0  + &
                   & (      FI2(2*M1  ,2*M2  ,  M3u )+      FI2(  M1l ,2*M2  ,2*M3-1)    -      &
                   &        FI2(2*M1  ,2*M2  ,2*M3-1)-      FI2(  M1l ,2*M2  ,  M3u ))   *      &
                   & (      FI2(2*M1  ,2*M2  ,2*M3  )-      FI2(  M1l ,2*M2  ,2*M3  ))   *      &
                   & (      FI2(2*M1-1,2*M2  ,  M3u )-      FI2(2*M1-1,2*M2  ,2*M3-1))/8.0D0  + &
                   & (      FI2(2*M1-1,  M2u ,  M3u )+      FI2(2*M1-1,2*M2-1,2*M3-1)    -      &
                   &        FI2(2*M1-1,  M2u ,2*M3-1)-      FI2(2*M1-1,2*M2-1,  M3u ))   *      &
                   & (      FI2(2*M1-1,  M2u ,2*M3  )-      FI2(2*M1-1,2*M2-1,2*M3  ))   *      &
                   & (      FI2(2*M1-1,2*M2  ,  M3u )-      FI2(2*M1-1,2*M2  ,2*M3-1))/8.0D0
              tmp2 = tmp2/dx4
              Sigma2 = Sigma2*tmp2
              S7   = (fct2*FI(2*M1-1,2*M2  ,2*M3  )+Sigma2)/Sigma1
              !
              ! Source term of cell #8 (2*M1  ,2*M2  ,2*M3  )
              tmp1 = (FI2(  M1u ,2*M2  ,2*M3  ) - FI2(2*M1-1,2*M2  ,2*M3  ))**2 + &
                   & (FI2(2*M1  ,  M2u ,2*M3  ) - FI2(2*M1  ,2*M2-1,2*M3  ))**2 + &
                   & (FI2(2*M1  ,2*M2  ,  M3u ) - FI2(2*M1  ,2*M2  ,2*M3-1))**2
              tmp1 = tmp1/(4.0D0*dx2)
              !
              IF(kmf_n.EQ.2) THEN
!                Sigma1 = tmp1+phidot2                   ! Baojiu-01-07-2021
                 Sigma1 = fct5*tmp1+fct3*phidot2         ! Baojiu-01-07-2021
                 Sigma2 = 1.0D0
              ENDIF
              IF(kmf_n.EQ.3) THEN
!                Sigma1 = tmp1**2+phidot2**2             ! Baojiu-01-07-2021
                 Sigma1 = fct5*tmp1**2+fct3*phidot2**2   ! Baojiu-01-07-2021
                 Sigma2 = tmp1
              ENDIF
!             Sigma1 = Sigma1*fct3+1.0D0                 ! Baojiu-01-07-2021
              Sigma1 = Sigma1+1.0D0                      ! Baojiu-01-07-2021
              Sigma2 = Sigma2*fct4
              !
              tmp2 = (      FI2(  M1u ,2*M2  ,2*M3  )-      FI2(2*M1-1,2*M2  ,2*M3  ))**2*      &
                   & (2.0D0*FI2(  M1u ,2*M2  ,2*M3  )+2.0D0*FI2(2*M1-1,2*M2  ,2*M3  )    -      &
                   &        FI2(2*M1  ,  M2u ,2*M3  )-      FI2(2*M1  ,2*M2-1,2*M3  )    -      &
                   &        FI2(2*M1  ,2*M2  ,  M3u )-      FI2(2*M1  ,2*M2  ,2*M3-1))/12.0D0 + &
                   & (      FI2(2*M1  ,  M2u ,2*M3  )-      FI2(2*M1  ,2*M2-1,2*M3  ))**2*      &
                   & (2.0D0*FI2(2*M1  ,  M2u ,2*M3  )+2.0D0*FI2(2*M1  ,2*M2-1,2*M3  )    -      &
                   &        FI2(  M1u ,2*M2  ,2*M3  )-      FI2(2*M1-1,2*M2  ,2*M3  )    -      &
                   &        FI2(2*M1  ,2*M2  ,  M3u )-      FI2(2*M1  ,2*M2  ,2*M3-1))/12.0D0 + &
                   & (      FI2(2*M1  ,2*M2  ,  M3u )-      FI2(2*M1  ,2*M2  ,2*M3-1))**2*      &
                   & (2.0D0*FI2(2*M1  ,2*M2  ,  M3u )+2.0D0*FI2(2*M1  ,2*M2  ,2*M3-1)    -      &
                   &        FI2(  M1u ,2*M2  ,2*M3  )-      FI2(2*M1-1,2*M2  ,2*M3  )    -      &
                   &        FI2(2*M1  ,  M2u ,2*M3  )-      FI2(2*M1  ,2*M2-1,2*M3  ))/12.0D0 + &
                   & (      FI2(  M1u ,  M2u ,2*M3  )+      FI2(2*M1-1,2*M2-1,2*M3  )    -      &
                   &        FI2(  M1u ,2*M2-1,2*M3  )-      FI2(2*M1-1,  M2u ,2*M3  ))   *      &
                   & (      FI2(  M1u ,2*M2  ,2*M3  )-      FI2(2*M1-1,2*M2  ,2*M3  ))   *      &
                   & (      FI2(2*M1  ,  M2u ,2*M3  )-      FI2(2*M1  ,2*M2-1,2*M3  ))/8.0D0  + &
                   & (      FI2(  M1u ,2*M2  ,  M3u )+      FI2(2*M1-1,2*M2  ,2*M3-1)    -      &
                   &        FI2(  M1u ,2*M2  ,2*M3-1)-      FI2(2*M1-1,2*M2  ,  M3u ))   *      &
                   & (      FI2(  M1u ,2*M2  ,2*M3  )-      FI2(2*M1-1,2*M2  ,2*M3  ))   *      &
                   & (      FI2(2*M1  ,2*M2  ,  M3u )-      FI2(2*M1  ,2*M2  ,2*M3-1))/8.0D0  + &
                   & (      FI2(2*M1  ,  M2u ,  M3u )+      FI2(2*M1  ,2*M2-1,2*M3-1)    -      &
                   &        FI2(2*M1  ,  M2u ,2*M3-1)-      FI2(2*M1  ,2*M2-1,  M3u ))   *      &
                   & (      FI2(2*M1  ,  M2u ,2*M3  )-      FI2(2*M1  ,2*M2-1,2*M3  ))   *      &
                   & (      FI2(2*M1  ,2*M2  ,  M3u )-      FI2(2*M1  ,2*M2  ,2*M3-1))/8.0D0
              tmp2 = tmp2/dx4
              Sigma2 = Sigma2*tmp2
              S8   = (fct2*FI(2*M1  ,2*M2  ,2*M3  )+Sigma2)/Sigma1
              !
              ! accumulate contributions to the non-Laplacian term 
              ! of PDE from all 8 son cells
              Q = S1+S2+S3+S4+S5+S6+S7+S8
              ! accumulate contributions to the Laplacian of PDE 
              ! from all 8 son cells
              P = P-3.0D0*FI2(2*M1-1,2*M2-1,2*M3-1)      
              P = P-3.0D0*FI2(2*M1  ,2*M2-1,2*M3-1)      
              P = P-3.0D0*FI2(2*M1-1,2*M2  ,2*M3-1)      
              P = P-3.0D0*FI2(2*M1  ,2*M2  ,2*M3-1)      
              P = P-3.0D0*FI2(2*M1-1,2*M2-1,2*M3  )      
              P = P-3.0D0*FI2(2*M1  ,2*M2-1,2*M3  )      
              P = P-3.0D0*FI2(2*M1-1,2*M2  ,2*M3  )      
              P = P-3.0D0*FI2(2*M1  ,2*M2  ,2*M3  )      
              ! accumulate further contributions to Laplacian PDE
              ! (1) contributions by 4 cells to the left
              P = P+FI2(M1l,2*M2-1,2*M3-1)
              P = P+FI2(M1l,2*M2  ,2*M3-1)      
              P = P+FI2(M1l,2*M2-1,2*M3  )      
              P = P+FI2(M1l,2*M2  ,2*M3  )      
              ! (2) contributions by 4 cells to the right
              P = P+FI2(M1u,2*M2-1,2*M3-1)      
              P = P+FI2(M1u,2*M2  ,2*M3-1)      
              P = P+FI2(M1u,2*M2-1,2*M3  )      
              P = P+FI2(M1u,2*M2  ,2*M3  )      
              ! (3) contributions by 4 cells behind
              P = P+FI2(2*M1-1,M2l,2*M3-1)      
              P = P+FI2(2*M1  ,M2l,2*M3-1)      
              P = P+FI2(2*M1-1,M2l,2*M3  )      
              P = P+FI2(2*M1  ,M2l,2*M3  )      
              ! (4) contributions by 4 cells in front
              P = P+FI2(2*M1-1,M2u,2*M3-1)      
              P = P+FI2(2*M1  ,M2u,2*M3-1)      
              P = P+FI2(2*M1-1,M2u,2*M3  )      
              P = P+FI2(2*M1  ,M2u,2*M3  )      
              ! (5) contributions by 4 cells below
              P = P+FI2(2*M1-1,2*M2-1,M3l)      
              P = P+FI2(2*M1  ,2*M2-1,M3l)      
              P = P+FI2(2*M1-1,2*M2  ,M3l)      
              P = P+FI2(2*M1  ,2*M2  ,M3l)      
              ! (6) contributions by 4 cells above
              P = P+FI2(2*M1-1,2*M2-1,M3u)      
              P = P+FI2(2*M1  ,2*M2-1,M3u)      
              P = P+FI2(2*M1-1,2*M2  ,M3u)      
              P = P+FI2(2*M1  ,2*M2  ,M3u)      
              !
              ! finalise residuals from all 8 son cells
              P = P/dx2*ctilde-Q
              !
              FI3(M1,M2,M3) = P/8.0D0  
           END DO
        END DO
     END DO
!    OPEN(UNIT=27, FILE='test.txt', FORM='FORMATTED', STATUS='REPLACE')
!       DO i=1,NGRID/2
!          WRITE(27,'(F20.9,F20.15)') (DBLE(i)-0.5d0)/DBLE(NGRID)*2.0d0*DACOS(-1.0d0),FI3(NGRID/4,NGRID/4,i)
!       ENDDO
!    CLOSE(27)
!    STOP
  ELSE
     ioffset  = 2**(levelmax-ilevel+1)
     koffset  = 2**(levelmax-ilevel)*(2**ilevel-2)
     koffset2 = 2**(levelmax-ilevel+1)*(2**(ilevel-1)-1)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3,M1u,M1l,M2u,M2l,M3u,M3l) &
!$OMP PRIVATE (P,Q) 
     DO M3=1,ngrid_level
        DO M2=1,ngrid_level
           DO M1=1,ngrid_level
              P =   FI3(ioffset+2*M1-1,2*M2-1,koffset2+2*M3-1) 
              P = P+FI3(ioffset+2*M1  ,2*M2-1,koffset2+2*M3-1) 
              P = P+FI3(ioffset+2*M1-1,2*M2  ,koffset2+2*M3-1)
              P = P+FI3(ioffset+2*M1  ,2*M2  ,koffset2+2*M3-1)
              P = P+FI3(ioffset+2*M1-1,2*M2-1,koffset2+2*M3  ) 
              P = P+FI3(ioffset+2*M1  ,2*M2-1,koffset2+2*M3  ) 
              P = P+FI3(ioffset+2*M1-1,2*M2  ,koffset2+2*M3  ) 
              P = P+FI3(ioffset+2*M1  ,2*M2  ,koffset2+2*M3  ) 
              !
              FI3(M1,M2,koffset+M3) = P/8.0D0  
              !
           END DO
        END DO
     END DO
  END IF
    
  CALL TimingMain(3,1)
END SUBROUTINE restrict_residual_kmf


!-------------------------------------------------------------
!
! Calculate the physical RHS of discrete PDE on coarser ilevel
!
!-------------------------------------------------------------
SUBROUTINE calculate_physical_right_hand_side_kmf(ilevel)
!-------------------------------------------------------------
  use Tools
  use ExtradofBackgroundData

  integer :: ilevel
  real*8  :: ctilde,dx,dx2,dx4
  real*8  :: fct1,fct2,fct3,fct4,fct5
  real*8  :: tmp1,tmp2
  real*8  :: Sigma1,Sigma2
  integer :: ngrid_level
  integer :: ioffset,joffset,koffest,koffset2
  integer :: i

  integer :: M1,M2,M3,M1l,M2l,M3l,M1u,M2u,M3u
  real*8  :: OP

  integer :: ist
  real*8  :: phibar,dphidN,HoH0,phidot,phidot2

  IF(MG_test) WRITE(*,'(A,I5)') 'Calculate physical right-hand side on level',levelmax-ilevel
  !
  ! number of grid points on the coarse level
  ngrid_level = NGRID/2**ilevel
  !
  ctilde  = 2.99792458D3*DBLE(NGRID)/Box
  dx      = DBLE(2**ilevel) 
  dx2     = dx*dx
  dx4     = dx2*dx2
  !
  fct1   = 0.5D0*(1.0D0/(AEXPN*kmf_lambda))**2
  fct2   = 3.0D0*kmf_beta*Om/AEXPN
  fct3   = kmf_n*kmf_K0*fct1**(kmf_n-1)                                 ! fct3 is gamma in notes
  fct4   = -fct3*(    DBLE(kmf_n)-1.0D0)*2.0D0*ctilde
  fct5   =  fct3*(2.0*DBLE(kmf_n)+1.0D0)/3.0D0                          ! fct5 is (2*n+1)/3 *gamma in notes; Baojiu-01-07-2021
  !
  ist = 1
  DO WHILE(BackgroundEvolution(ist,1).LT.AEXPN)
     ist=ist+1
  ENDDO
  !
  ! background scalar field at AEXPN
  phibar =  BackgroundEvolution(ist-1,2)+ &
         & (BackgroundEvolution(ist,  2)-BackgroundEvolution(ist-1,2))/ &
         & (BackgroundEvolution(ist,  1)-BackgroundEvolution(ist-1,1))* &
         & (AEXPN                       -BackgroundEvolution(ist-1,1))
  ! d\varphi/dN at AEXPN, N=ln(a)
  dphidN =  BackgroundEvolution(ist-1,3)+ &
         & (BackgroundEvolution(ist,  3)-BackgroundEvolution(ist-1,3))/ &
         & (BackgroundEvolution(ist,  1)-BackgroundEvolution(ist-1,1))* &
         & (AEXPN                       -BackgroundEvolution(ist-1,1))
  ! H/H0 at AEXPN, here H=a'/a=dN/d\tau with '=d/d\tau
  HoH0   =  BackgroundEvolution(ist-1,5)+ &
         & (BackgroundEvolution(ist,  5)-BackgroundEvolution(ist-1,5))/ &
         & (BackgroundEvolution(ist,  1)-BackgroundEvolution(ist-1,1))* &
         & (AEXPN                       -BackgroundEvolution(ist-1,1))
  !
  fct2 = fct2*DEXP(kmf_beta*phibar)
  !
  phidot = dphidN*HoH0                                                  ! (d\varphi/d\tau)/H0
  phidot2 = phidot**2                                                   ! [(d\varphi/d\tau)/H0]^2
  !
  ioffset  = 2**(levelmax-ilevel)
  joffset  = 2**(levelmax-ilevel)
  koffset  = 2**(levelmax-ilevel)*(2**ilevel-2)
  koffset2 = 2**(levelmax-ilevel)*(2**ilevel-1)

!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3,M1u,M1l,M2u,M2l,M3u,M3l,i) &
!$OMP PRIVATE (tmp1,tmp2,Sigma1,Sigma2,OP) 
  DO M3=1,ngrid_level  
     DO M2=1,ngrid_level  
        DO M1=1,ngrid_level  
           ! prepare indices of neighbour cells on the 3-point
           ! stencil and apply the periodic boundary condition
           M1u = M1+1; IF(M1u>ngrid_level) M1u=1 
           M1l = M1-1; IF(M1l<1          ) M1l=ngrid_level
           M2u = M2+1; IF(M2u>ngrid_level) M2u=1
           M2l = M2-1; IF(M2l<1          ) M2l=ngrid_level
           M3u = M3+1; IF(M3u>ngrid_level) M3u=1
           M3l = M3-1; IF(M3l<1          ) M3l=ngrid_level
           !
           ! calculate Laplacian of PDE using RESTRICTED scalar field
           OP = FI3(M1u,joffset+M2 ,koffset+M3 )      + &
              & FI3(M1l,joffset+M2 ,koffset+M3 )      + &
              & FI3(M1 ,joffset+M2u,koffset+M3 )      + &
              & FI3(M1 ,joffset+M2l,koffset+M3 )      + &
              & FI3(M1 ,joffset+M2 ,koffset+M3u)      + &
              & FI3(M1 ,joffset+M2 ,koffset+M3l)      - &
              & FI3(M1 ,joffset+M2 ,koffset+M3 )*6.0D0 
           !
           tmp1 = (FI3(M1u,joffset+M2 ,koffset+M3 ) - FI3(M1l,joffset+M2 ,koffset+M3 ))**2 + &
                & (FI3(M1 ,joffset+M2u,koffset+M3 ) - FI3(M1 ,joffset+M2l,koffset+M3 ))**2 + &
                & (FI3(M1 ,joffset+M2 ,koffset+M3u) - FI3(M1 ,joffset+M2 ,koffset+M3l))**2
           tmp1 = tmp1/(4.0D0*dx2)
           !
           IF(kmf_n.EQ.2) THEN
!             Sigma1 = tmp1+phidot2                  ! Baojiu-01-07-2021
              Sigma1 = fct5*tmp1+fct3*phidot2        ! Baojiu-01-07-2021
              Sigma2 = 1.0D0
           ENDIF
           IF(kmf_n.EQ.3) THEN
!             Sigma1 = tmp1**2+phidot2**2            ! Baojiu-01-07-2021
              Sigma1 = fct5*tmp1**2+fct3*phidot2**2  ! Baojiu-01-07-2021
              Sigma2 = tmp1
           ENDIF
!          Sigma1 = Sigma1*fct3+1.0D0                ! Baojiu-01-07-2021
           Sigma1 = Sigma1+1.0D0                     ! Baojiu-01-07-2021
           Sigma2 = Sigma2*fct4
           !
           tmp2 = (      FI3(M1u,joffset+M2 ,koffset+M3 )-      FI3(M1l,joffset+M2 ,koffset+M3 ))**2*      &
                & (2.0D0*FI3(M1u,joffset+M2 ,koffset+M3 )+2.0D0*FI3(M1l,joffset+M2 ,koffset+M3 )    -      &
                &        FI3(M1 ,joffset+M2u,koffset+M3 )-      FI3(M1 ,joffset+M2l,koffset+M3 )    -      &
                &        FI3(M1 ,joffset+M2 ,koffset+M3u)-      FI3(M1 ,joffset+M2 ,koffset+M3l))/12.0D0 + &
                & (      FI3(M1 ,joffset+M2u,koffset+M3 )-      FI3(M1 ,joffset+M2l,koffset+M3 ))**2*      &
                & (2.0D0*FI3(M1 ,joffset+M2u,koffset+M3 )+2.0D0*FI3(M1 ,joffset+M2l,koffset+M3 )    -      &
                &        FI3(M1u,joffset+M2 ,koffset+M3 )-      FI3(M1l,joffset+M2 ,koffset+M3 )    -      &
                &        FI3(M1 ,joffset+M2 ,koffset+M3u)-      FI3(M1 ,joffset+M2 ,koffset+M3l))/12.0D0 + &
                & (      FI3(M1 ,joffset+M2 ,koffset+M3u)-      FI3(M1 ,joffset+M2 ,koffset+M3l))**2*      &
                & (2.0D0*FI3(M1 ,joffset+M2 ,koffset+M3u)+2.0D0*FI3(M1 ,joffset+M2 ,koffset+M3l)    -      &
                &        FI3(M1u,joffset+M2 ,koffset+M3 )-      FI3(M1l,joffset+M2 ,koffset+M3 )    -      &
                &        FI3(M1 ,joffset+M2u,koffset+M3 )-      FI3(M1 ,joffset+M2l,koffset+M3 ))/12.0D0 + &
                & (      FI3(M1u,joffset+M2u,koffset+M3 )+      FI3(M1l,joffset+M2l,koffset+M3 )    -      &
                &        FI3(M1u,joffset+M2l,koffset+M3 )-      FI3(M1l,joffset+M2u,koffset+M3 ))   *      &
                & (      FI3(M1u,joffset+M2 ,koffset+M3 )-      FI3(M1l,joffset+M2 ,koffset+M3 ))   *      &
                & (      FI3(M1 ,joffset+M2u,koffset+M3 )-      FI3(M1 ,joffset+M2l,koffset+M3 ))/8.0D0  + &
                & (      FI3(M1u,joffset+M2 ,koffset+M3u)+      FI3(M1l,joffset+M2 ,koffset+M3l)    -      &
                &        FI3(M1u,joffset+M2 ,koffset+M3l)-      FI3(M1l,joffset+M2 ,koffset+M3u))   *      &
                & (      FI3(M1u,joffset+M2 ,koffset+M3 )-      FI3(M1l,joffset+M2 ,koffset+M3 ))   *      &
                & (      FI3(M1 ,joffset+M2 ,koffset+M3u)-      FI3(M1 ,joffset+M2 ,koffset+M3l))/8.0D0  + &
                & (      FI3(M1 ,joffset+M2u,koffset+M3u)+      FI3(M1 ,joffset+M2l,koffset+M3l)    -      &
                &        FI3(M1 ,joffset+M2u,koffset+M3l)-      FI3(M1 ,joffset+M2l,koffset+M3u))   *      &
                & (      FI3(M1 ,joffset+M2u,koffset+M3 )-      FI3(M1 ,joffset+M2l,koffset+M3 ))   *      &
                & (      FI3(M1 ,joffset+M2 ,koffset+M3u)-      FI3(M1 ,joffset+M2 ,koffset+M3l))/8.0D0
           tmp2 = tmp2/dx4
           Sigma2 = Sigma2*tmp2 
           !
           OP = OP/dx2*ctilde
           OP = OP-(fct2*FI3(M1+ioffset,M2,M3+koffset)+Sigma2)/Sigma1
           OP = OP-FI3(M1,M2,M3+koffset)                        ! add restricted residual (note the minus sign!!)
           ! 
           FI3(M1,joffset+M2,koffset2+M3) = OP  
           ! 
        END DO
     END DO
  END DO

! OPEN(UNIT=27, FILE='test.txt', FORM='FORMATTED', STATUS='REPLACE')
!    DO M3=1,NGRID
!       WRITE(27,'(F20.9,F20.15)') (DBLE(M3)-0.5d0)/DBLE(NGRID)*2.0d0*DACOS(-1.0d0),FI3(NGRID/2,joffset+NGRID/2,koffset2+M3)
!    ENDDO
! CLOSE(27)
! STOP

  CALL TimingMain(3,1)

END SUBROUTINE calculate_physical_right_hand_side_kmf
