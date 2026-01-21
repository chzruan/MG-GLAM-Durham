!================================================= 
!================================================= 
!
!  Scalar field
!
!================================================= 
!================================================= 

! Baojiu: changed to double precision throughout the file

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
SUBROUTINE relaxation_iterations_fR(ilevel,redstep)
!-------------------------------------------------------------
  use Tools

  implicit none

  integer :: ilevel
  logical :: redstep
  real*8  :: Q,Q2onethird,Q2,Q4,CC,R0_bg,R_bg,fR_bg,ctilde,dx,dx2
  real*8  :: s_fR_bg,a_fR_bg
  real*8  :: fct1,fct2                            
  real*8  :: ZERO,ONEOVERTHREE,TWOOVERTHREE,TWOPIOVERTHREE
  integer :: nstart,nfinal,ngrid_level
  integer :: ioffset,joffset,koffset,koffset2

  real*8  :: Delta10,Delta11,Delta20,Delta21,tmp1,tmp2,P,S,L,sqrt_Delta0
  integer :: M1,M2,M3,M1u,M1l,M2u,M2l,M3u,M3l

  integer :: nid,myid
! real*8  :: Pmax,Pmin

  IF(MG_test .EQ. 1) WRITE(*,'(A,I5,F7.4)') 'Relaxation iterations on level',levelmax-ilevel,AEXPN

! Pmax = 0.0D0
! Pmin = 1.0D10

  ! constants
  ZERO           = 1.0D-20
  ONEOVERTHREE   = 1.0D0/3.0D0
  TWOOVERTHREE   = 2.0D0/3.0D0
  TWOPIOVERTHREE = 2.0D0/3.0D0*DACOS(-1.0D0)

  ! number of grid points on the coarse level
  ngrid_level = NGRID/(2**ilevel)

  ! simulation parameters
  ctilde = 2.99792458D3*DBLE(NGRID)/Box
  dx     = DBLE(2**ilevel) 
  dx2    = dx*dx

  ! physical and numerical quantities
  R_bg   = 3.0d0*(Om/AEXPN**3+4.0d0*OmL)
  R0_bg  = 3.0d0*(Om         +4.0d0*OmL)
  fR_bg  = fR0*(R0_bg/R_bg)**(fr_n+1)
  fct2   = dx2/(ctilde**2*DSIGN(1.0D0,fR_bg))/6.0d0
  fct1   = -Om/AEXPN*fct2
  Q2     = R_bg*AEXPN**2*fct2/3.0d0
  Q      = Q2*(DABS(fR_bg))**(1.0D0/(DBLE(fr_n)+1.0D0))

  s_fR_bg = DSIGN(1.0D0,fR_bg)

  ! intermediate quantities
  Q2onethird = (-Q)**ONEOVERTHREE
  Q4      = 4.0D0 *Q
  Delta11 = 27.0d0*Q
  Delta20 = 12.0d0*Q

  ! offsets for cell access to array FI3
! ioffset  = 2**(levelmax-ilevel)
! joffset  = 2**(levelmax-ilevel)
! koffset  = 0
! koffset2 = 0 
! IF(ilevel.GT.0) THEN 
!    koffset  = 2**(levelmax-ilevel)*(2**ilevel-2)
!    koffset2 = 2**(levelmax-ilevel)*(2**ilevel-1)
! END IF
  ioffset  = NGRID/2**(ilevel)                                                  ! Baojiu 2025 12 27
  joffset  = NGRID/2**(ilevel)                                                  ! Baojiu 2025 12 27
  koffset  = 0                                                                  ! Baojiu 2025 12 27
  koffset2 = 0                                                                  ! Baojiu 2025 12 27
  IF(ilevel.GT.0) THEN                                                          ! Baojiu 2025 12 27
     koffset  = NGRID/2**(ilevel)*(2**ilevel-2)                                 ! Baojiu 2025 12 27
     koffset2 = NGRID/2**(ilevel)*(2**ilevel-1)                                 ! Baojiu 2025 12 27
  END IF                                                                        ! Baojiu 2025 12 27

!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3,M1u,M1l,M2u,M2l,M3u,M3l) &
!$OMP PRIVATE (Delta21,Delta10,tmp1,tmp2,sqrt_Delta0) &
!$OMP PRIVATE (P,L,S,CC)
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
           ! summation of u^(n+1) of all 6 neighbors on the 3-point stencil
           ! Case n=0
           IF(fr_n.EQ.0) THEN
              IF(ilevel.EQ.0) THEN
                 P = fct1*FI (M1 ,M2 ,M3 )-Q2-     &
                       & (FI2(M1u,M2 ,M3 )       + &
                       &  FI2(M1l,M2 ,M3 )       + &
                       &  FI2(M1 ,M2u,M3 )       + &
                       &  FI2(M1 ,M2l,M3 )       + &
                       &  FI2(M1 ,M2 ,M3u)       + &
                       &  FI2(M1 ,M2 ,M3l)       )/6.0D0
              ELSE 
                 P = fct1*FI3(M1 +ioffset,M2         ,M3 +koffset )-Q2+     &     ! add (restricted) densities from FI3 Section 2
                   & fct2*FI3(M1         ,M2 +joffset,M3 +koffset2)-        &     ! add (restricted) residuals from FI3 Section 6       
                       & (FI3(M1u+ioffset,M2 +joffset,M3 +koffset )       + &     ! add coarse l scalar fields from FI3 Section 4
                       &  FI3(M1l+ioffset,M2 +joffset,M3 +koffset )       + &
                       &  FI3(M1 +ioffset,M2u+joffset,M3 +koffset )       + &
                       &  FI3(M1 +ioffset,M2l+joffset,M3 +koffset )       + &
                       &  FI3(M1 +ioffset,M2 +joffset,M3u+koffset )       + &
                       &  FI3(M1 +ioffset,M2 +joffset,M3l+koffset )       )/6.0D0
              ENDIF
           ! Case n=1
           ELSE IF(fr_n.EQ.1) THEN
              IF(ilevel.EQ.0) THEN
                 P = fct1*FI (M1 ,M2 ,M3 )-Q2-     &
                       & (FI2(M1u,M2 ,M3 )**2    + &
                       &  FI2(M1l,M2 ,M3 )**2    + &
                       &  FI2(M1 ,M2u,M3 )**2    + &
                       &  FI2(M1 ,M2l,M3 )**2    + &
                       &  FI2(M1 ,M2 ,M3u)**2    + &
                       &  FI2(M1 ,M2 ,M3l)**2    )/6.0D0
              ELSE 
                 P = fct1*FI3(M1 +ioffset,M2         ,M3 +koffset )-Q2+     &     ! add (restricted) densities from FI3 Section 2
                   & fct2*FI3(M1         ,M2 +joffset,M3 +koffset2)-        &     ! add (restricted) residuals from FI3 Section 6       
                       & (FI3(M1u+ioffset,M2 +joffset,M3 +koffset )**2    + &     ! add coarse l scalar fields from FI3 Section 4
                       &  FI3(M1l+ioffset,M2 +joffset,M3 +koffset )**2    + &
                       &  FI3(M1 +ioffset,M2u+joffset,M3 +koffset )**2    + &
                       &  FI3(M1 +ioffset,M2l+joffset,M3 +koffset )**2    + &
                       &  FI3(M1 +ioffset,M2 +joffset,M3u+koffset )**2    + &
                       &  FI3(M1 +ioffset,M2 +joffset,M3l+koffset )**2    )/6.0D0
              ENDIF
           ! Case n=2
           ELSE IF(fr_n.EQ.2) THEN
              IF(ilevel.EQ.0) THEN
                 P = fct1*FI (M1 ,M2 ,M3 )-Q2-     &
                       & (FI2(M1u,M2 ,M3 )**3    + &
                       &  FI2(M1l,M2 ,M3 )**3    + &
                       &  FI2(M1 ,M2u,M3 )**3    + &
                       &  FI2(M1 ,M2l,M3 )**3    + &
                       &  FI2(M1 ,M2 ,M3u)**3    + &
                       &  FI2(M1 ,M2 ,M3l)**3    )/6.0D0
              ELSE 
                 P = fct1*FI3(M1 +ioffset,M2         ,M3 +koffset )-Q2+     &    ! add (restricted) densities from FI3 Section 2
                   & fct2*FI3(M1         ,M2 +joffset,M3 +koffset2)-        &    ! add (restricted) residuals from FI3 Section 6        
                       & (FI3(M1u+ioffset,M2 +joffset,M3 +koffset )**3    + &    ! add coarse l scalar fields from FI3 Section 4
                       &  FI3(M1l+ioffset,M2 +joffset,M3 +koffset )**3    + &
                       &  FI3(M1 +ioffset,M2u+joffset,M3 +koffset )**3    + &
                       &  FI3(M1 +ioffset,M2l+joffset,M3 +koffset )**3    + &
                       &  FI3(M1 +ioffset,M2 +joffset,M3u+koffset )**3    + &
                       &  FI3(M1 +ioffset,M2 +joffset,M3l+koffset )**3    )/6.0D0
              END IF
           ENDIF
           !
           ! update solution (the case of n=0)
           IF(fr_n.EQ.0) THEN
              IF(ilevel.EQ.0) THEN      ! solution on PM grid
                 FI2(M1        ,M2        ,M3        ) = 0.5D0*(-P+DSQRT(P**2-Q4))
              ELSE                      ! solution on multigrid (FI3 Section 4)
                 FI3(M1+ioffset,M2+joffset,M3+koffset) = 0.5D0*(-P+DSQRT(P**2-Q4)) 
              END IF
           ENDIF
           ! update solution (the case of n=1)
           IF(fr_n.EQ.1) THEN
              Delta10 = -3.0D0*P
              !
              CC = Delta11**2-4.0d0*Delta10**3
              !
              IF(ilevel.EQ.0) THEN
              ! solution on PM grid
                 IF(CC.GE.0.0D0) THEN
                    tmp1 = 0.5d0*(DSQRT(CC)+Delta11)
                    tmp2 = tmp1-Delta11
                    IF(DABS(P).LT.ZERO) THEN
                       FI2(M1,M2,M3) = Q2onethird 
                    ELSE
                       FI2(M1,M2,M3) = -(DSIGN(DABS(tmp1)**ONEOVERTHREE,tmp1)-DSIGN(DABS(tmp2)**ONEOVERTHREE,tmp2))*ONEOVERTHREE
                    ENDIF
                 ELSE
                    sqrt_Delta0 = DSQRT(Delta10)
                    S = DACOS(0.5d0*Delta11/(sqrt_Delta0**3))
                    FI2(M1,M2,M3) = -TWOOVERTHREE*sqrt_Delta0*DCOS(S*ONEOVERTHREE+TWOPIOVERTHREE)
                 ENDIF
              ELSE
              ! solution on multigrid (FI3 Section 4)
                 IF(CC.GE.0.0D0) THEN
                    tmp1 = 0.5d0*(DSQRT(CC)+Delta11)
                    tmp2 = tmp1-Delta11
                    IF(DABS(P).LT.ZERO) THEN
                       FI3(M1+ioffset,M2+joffset,M3+koffset) = Q2onethird
                    ELSE
                       FI3(M1+ioffset,M2+joffset,M3+koffset) = -(DSIGN(DABS(tmp1)**ONEOVERTHREE,tmp1)-DSIGN(DABS(tmp2)**ONEOVERTHREE,tmp2))*ONEOVERTHREE
                    ENDIF
                 ELSE
                    sqrt_Delta0 = DSQRT(Delta10)
                    S = DACOS(0.5d0*Delta11/sqrt_Delta0**3)
                    FI3(M1+ioffset,M2+joffset,M3+koffset) = -TWOOVERTHREE*sqrt_Delta0*DCOS(S*ONEOVERTHREE+TWOPIOVERTHREE)
                 ENDIF
              ENDIF
           END IF
           ! update solution (the case of n=2)
           IF(fr_n.EQ.2) THEN
              Delta21 = 27.0d0*P**2
              tmp1    = 0.5d0*(DSQRT(Delta21**2-4.0d0*Delta20**3)+Delta21)
              tmp2    = tmp1-Delta21
!             S       = DSQRT(tmp1**ONEOVERTHREE-tmp2**ONEOVERTHREE)/DSQRT(12.0d0)
!             S       = DSQRT(DSIGN(DABS(tmp1)**ONEOVERTHREE,tmp1)-DSIGN(DABS(tmp2)**ONEOVERTHREE,tmp2))/DSQRT(12.0d0)
              S       = DSQRT(DSIGN(DABS(tmp1)**ONEOVERTHREE, tmp1) + Delta20 / DSIGN(DABS(tmp1)**ONEOVERTHREE, tmp1)) / DSQRT(12.0d0)

              IF(ilevel.EQ.0) THEN
              ! solution on PM grid
                 IF(P.GT.ZERO) THEN
                    FI2(M1,M2,M3) = - S+0.5d0*DSQRT(-4.0d0*S**2+P/S)
                 ELSE IF(P.LT.-ZERO) THEN
                    FI2(M1,M2,M3) =   S+0.5d0*DSQRT(-4.0d0*S**2-P/S)
                 ELSE
                    FI2(M1,M2,M3) = DSQRT(DSQRT(-Q))
                 END IF
              ELSE
              ! solution on multigrid (FI3 Section 4)
                 IF(P.GT.ZERO) THEN
                    FI3(M1+ioffset,M2+joffset,M3+koffset) = - S + 0.5d0*DSQRT(-4.0d0*S**2+P/S)
                 ELSE IF(P.LT.-ZERO) THEN
                    FI3(M1+ioffset,M2+joffset,M3+koffset) =   S + 0.5d0*DSQRT(-4.0d0*S**2-P/S)
                 ELSE
                    FI3(M1+ioffset,M2+joffset,M3+koffset) = DSQRT(DSQRT(-Q))
                 END IF
              END IF
           END IF
        END DO
     END DO
  END DO

  CALL TimingMain(3,1)
    
END SUBROUTINE relaxation_iterations_fR

!-------------------------------------------------------------
!
! Calculate resudual and its RMS on the entire grid on ilevel.
!
!-------------------------------------------------------------
SUBROUTINE calculate_residual_fR(ilevel,res_PM_grid)
!-------------------------------------------------------------
use Tools

  integer,intent(in) :: ilevel
  real*8,intent(out) :: res_PM_grid
  
  real*8  :: R0_bg,R_bg,fR_bg,ctilde,ctilde2,dx,dx2
  integer :: ngrid_level
  integer :: ioffset,joffset,koffset,koffset2
  real*8  :: fct1,fct2,fct3
  real*8  :: s_fR_bg,a_fR_bg,a_fR_rt                            ! sign and amplitude of f_R(a)
  integer :: fr_n1

  integer :: M1,M2,M3,M1l,M1u,M2l,M2u,M3l,M3u
  integer :: N1,N2,N3
  real*8  :: RES,OP,RES2

  IF(MG_test .EQ. 1) WRITE(*,'(A,I5)') 'Calculate residual on level',levelmax-ilevel

  ! number of grid points on the coarse level
  ngrid_level = NGRID/2**ilevel
  
  ctilde  = 2.99792458D3*DBLE(NGRID)/Box
  ctilde2 = ctilde**2
  dx      = DBLE(2**ilevel)
  dx2     = dx*dx
  R_bg    = 3.0d0*(Om/AEXPN**3+4.0d0*OmL)
  R0_bg   = 3.0d0*(Om         +4.0d0*OmL)
  fR_bg   = fR0*(R0_bg/R_bg)**(fr_n+1)

  s_fR_bg = DSIGN(1.0D0,fR_bg)
  a_fR_bg = DABS (      fR_bg)
  fr_n1   = fr_n+1
  a_fR_rt = a_fR_bg**(1.0D0/(DBLE(fr_n1)))

  fct1    = Om/AEXPN
  fct2    = -R_bg*AEXPN**2/3.0d0
  fct3    = ctilde2*s_fR_bg/dx2
 
! ioffset  = 2**(levelmax-ilevel)
! joffset  = 2**(levelmax-ilevel)
! koffset  = 0
! koffset2 = 0 
! IF(ilevel.GT.0) THEN 
!    koffset  = 2**(levelmax-ilevel)*(2**ilevel-2)
!    koffset2 = 2**(levelmax-ilevel)*(2**ilevel-1)
! END IF
  ioffset  = NGRID/2**(ilevel)                                                  ! Baojiu 2025 12 17
  joffset  = NGRID/2**(ilevel)                                                  ! Baojiu 2025 12 17
  koffset  = 0                                                                  ! Baojiu 2025 12 17
  koffset2 = 0                                                                  ! Baojiu 2025 12 17
  IF(ilevel.GT.0) THEN                                                          ! Baojiu 2025 12 17
     koffset  = NGRID/2**(ilevel)*(2**ilevel-2)                                 ! Baojiu 2025 12 17
     koffset2 = NGRID/2**(ilevel)*(2**ilevel-1)                                 ! Baojiu 2025 12 17
  END IF                                                                        ! Baojiu 2025 12 17

!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3,M1u,M1l,M2u,M2l,M3u,M3l) &
!$OMP PRIVATE (OP)
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
           IF(fr_n.EQ.0) THEN
              IF(ilevel.EQ.0) THEN
                 ! calculate Laplacian of PDE
                 OP = FI2(M1u,M2 ,M3 )         + &
                    & FI2(M1l,M2 ,M3 )         + &
                    & FI2(M1 ,M2u,M3 )         + &
                    & FI2(M1 ,M2l,M3 )         + &
                    & FI2(M1 ,M2 ,M3u)         + &
                    & FI2(M1 ,M2 ,M3l)         - &
                    & FI2(M1 ,M2 ,M3 )*6.0D0
                 OP = OP*fct3
                 OP = OP+fct1*FI(M1,M2,M3)
                 OP = OP+fct2*(a_fR_rt/FI2(M1,M2,M3)-1.0D0)
                 FI3(M1,M2,M3) = OP
              ELSE   
                 ! calculate Laplacian of PDE - all scalar fields are from Section 4 of FI3
                 OP = FI3(M1u+ioffset,M2 +joffset,M3 +koffset)         + &
                    & FI3(M1l+ioffset,M2 +joffset,M3 +koffset)         + &
                    & FI3(M1 +ioffset,M2u+joffset,M3 +koffset)         + &
                    & FI3(M1 +ioffset,M2l+joffset,M3 +koffset)         + &
                    & FI3(M1 +ioffset,M2 +joffset,M3u+koffset)         + &
                    & FI3(M1 +ioffset,M2 +joffset,M3l+koffset)         - &
                    & FI3(M1 +ioffset,M2 +joffset,M3 +koffset)*6.0D0
                 OP = OP*fct3
                 OP = OP+fct1*FI3(M1+ioffset,M2,M3+koffset)                     ! add (restricted) density field (FI3 Section 2)
                 OP = OP+fct2*(a_fR_rt/FI3(M1+ioffset,M2+joffset,M3+koffset)-1.0D0) ! scalar field from FI3 Section 4
                 OP = OP-FI3(M1,M2+joffset,M3+koffset2)                         ! physical RHS from FI3 Section 7
                 !
                 FI3(M1+ioffset,M2,M3+koffset2) = OP                            ! store coarse level residual in FI3 Section 6
                 !
              END IF
           ELSE IF(fr_n.EQ.1) THEN
              IF(ilevel.EQ.0) THEN
                 ! calculate Laplacian of PDE - all scalar fields are from Section 4 of FI3
                 OP = FI2(M1u,M2 ,M3 )**2      + &
                    & FI2(M1l,M2 ,M3 )**2      + &
                    & FI2(M1 ,M2u,M3 )**2      + &
                    & FI2(M1 ,M2l,M3 )**2      + &
                    & FI2(M1 ,M2 ,M3u)**2      + &
                    & FI2(M1 ,M2 ,M3l)**2      - &
                    & FI2(M1 ,M2 ,M3 )**2*6.0D0
                 OP = OP*fct3
                 OP = OP+fct1*FI(M1,M2,M3)
                 OP = OP+fct2*(a_fR_rt/FI2(M1,M2,M3)-1.0D0)
                 FI3(M1,M2,M3) = OP
              ELSE   
                 ! calculate Laplacian of PDE
                 OP = FI3(M1u+ioffset,M2 +joffset,M3 +koffset)**2      + &
                    & FI3(M1l+ioffset,M2 +joffset,M3 +koffset)**2      + &
                    & FI3(M1 +ioffset,M2u+joffset,M3 +koffset)**2      + &
                    & FI3(M1 +ioffset,M2l+joffset,M3 +koffset)**2      + &
                    & FI3(M1 +ioffset,M2 +joffset,M3u+koffset)**2      + &
                    & FI3(M1 +ioffset,M2 +joffset,M3l+koffset)**2      - &
                    & FI3(M1 +ioffset,M2 +joffset,M3 +koffset)**2*6.0D0
                 OP = OP*fct3
                 OP = OP+fct1*FI3(M1+ioffset,M2,M3+koffset)                     ! add (restricted) density field (FI3 Section 2)
                 OP = OP+fct2*(a_fR_rt/FI3(M1+ioffset,M2+joffset,M3+koffset)-1.0d0) ! scalar field from FI3 Section 4
                 OP = OP-FI3(M1,M2+joffset,M3+koffset2)                         ! physical RHS from FI3 Section 7
                 !
                 FI3(M1+ioffset,M2,M3+koffset2) = OP                            ! store coarse level residual in FI3 Section 6
                 !
              END IF
           ELSE IF(fr_n.EQ.2) THEN
              IF(ilevel.EQ.0) THEN
                 ! calculate Laplacian of PDE - all scalar fields are from Section 4 of FI3
                 OP = FI2(M1u,M2 ,M3 )**3      + &
                    & FI2(M1l,M2 ,M3 )**3      + &
                    & FI2(M1 ,M2u,M3 )**3      + &
                    & FI2(M1 ,M2l,M3 )**3      + &
                    & FI2(M1 ,M2 ,M3u)**3      + &
                    & FI2(M1 ,M2 ,M3l)**3      - &
                    & FI2(M1 ,M2 ,M3 )**3*6.0D0
                 OP = OP*fct3
                 OP = OP+fct1*FI(M1,M2,M3)
                 OP = OP+fct2*(a_fR_rt/FI2(M1,M2,M3)-1.0D0)
                 FI3(M1,M2,M3) = OP
              ELSE   
                 ! calculate Laplacian of PDE
                 OP = FI3(M1u+ioffset,M2 +joffset,M3 +koffset)**3      + &
                    & FI3(M1l+ioffset,M2 +joffset,M3 +koffset)**3      + &
                    & FI3(M1 +ioffset,M2u+joffset,M3 +koffset)**3      + &
                    & FI3(M1 +ioffset,M2l+joffset,M3 +koffset)**3      + &
                    & FI3(M1 +ioffset,M2 +joffset,M3u+koffset)**3      + &
                    & FI3(M1 +ioffset,M2 +joffset,M3l+koffset)**3      - &
                    & FI3(M1 +ioffset,M2 +joffset,M3 +koffset)**3*6.0D0
                 OP = OP*fct3
                 OP = OP+fct1*FI3(M1+ioffset,M2,M3+koffset)                     ! add (restricted) density field (FI3 Section 2)
                 OP = OP+fct2*(a_fR_rt/FI3(M1+ioffset,M2+joffset,M3+koffset)-1.0D0) ! scalar field from FI3 Section 4
                 OP = OP-FI3(M1,M2+joffset,M3+koffset2)                         ! physical RHS from FI3 Section 7
                 !
                 FI3(M1+ioffset,M2,M3+koffset2) = OP                            ! store coarse level residual in FI3 Section 6
                 !
              END IF
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
! WRITE(*,'(A,F20.16,F20.16,A,I5,A)') 'Residual RES2 = ',RES2,DSQRT(DABS(fR_bg)),' on level ',levelmax-ilevel  
      
END SUBROUTINE calculate_residual_fR


!-------------------------------------------------------------
!
! Calculate restricted residual field on coarse level "ilevel"
! Note: "ilevel" means the grid ilevel coarer than the PM grid
!
!-------------------------------------------------------------
SUBROUTINE restrict_residual_fR(ilevel)
!-------------------------------------------------------------
use Tools

  integer :: ilevel
  real*8  :: R0_bg,R_bg,fR_bg,ctilde,ctilde2,dx,dx2
  real*8  :: fct1,fct2,fct3,fct4 
  integer :: ngrid_level 
  integer :: ioffset,joffset,koffset,ioffset2,joffset2,koffset2
  real*8  :: a_fR_bg,s_fR_bg                                   ! amplitude and sign of f_R(a)
 
  integer :: M1,M2,M3
  integer :: M1u,M1l,M2u,M2l,M3u,M3l 
  real*8  :: P,Q

  IF(MG_test .EQ. 1) WRITE(*,'(A,I5)') 'Restrict residual to level',levelmax-ilevel

  ! number of grid points on the coarse level
  ngrid_level = NGRID/2**ilevel

  ctilde  = 2.99792458D3*DBLE(NGRID)/Box
  ctilde2 = ctilde**2
  dx      = 1.0d0  
  dx2     = dx*dx
  R_bg    = 3.0d0*(Om/AEXPN**3+4.0d0*OmL)
  R0_bg   = 3.0d0*(Om         +4.0d0*OmL)
  fR_bg   = fR0*(R0_bg/R_bg)**(fr_n+1)

  s_fR_bg = DSIGN(1.0D0,fR_bg)
  a_fR_bg = DABS (      fR_bg)

  fct1    = Om/AEXPN
  fct2    = -R_bg*AEXPN**2/3.0d0
  fct3    = ctilde2*s_fR_bg/dx2
  fct4    = fct2*a_fR_bg**(1.0D0/(1.0D0+DBLE(fr_n)))

  IF(ilevel.EQ.1) THEN
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3,M1u,M1l,M2u,M2l,M3u,M3l) &
!$OMP PRIVATE (P,Q) 
     DO M3=1,ngrid_level
        DO M2=1,ngrid_level
           DO M1=1,ngrid_level
              P = 0.0d0 
              Q = -8.0d0*fct2
              ! accumulate contributions to the nonLaplacian term 
              ! of PDE from all 8 son cells
              Q = Q+fct1*FI (2*M1-1,2*M2-1,2*M3-1)+fct4/FI2(2*M1-1,2*M2-1,2*M3-1)
              Q = Q+fct1*FI (2*M1  ,2*M2-1,2*M3-1)+fct4/FI2(2*M1  ,2*M2-1,2*M3-1)
              Q = Q+fct1*FI (2*M1-1,2*M2  ,2*M3-1)+fct4/FI2(2*M1-1,2*M2  ,2*M3-1)
              Q = Q+fct1*FI (2*M1  ,2*M2  ,2*M3-1)+fct4/FI2(2*M1  ,2*M2  ,2*M3-1)
              Q = Q+fct1*FI (2*M1-1,2*M2-1,2*M3  )+fct4/FI2(2*M1-1,2*M2-1,2*M3  )
              Q = Q+fct1*FI (2*M1  ,2*M2-1,2*M3  )+fct4/FI2(2*M1  ,2*M2-1,2*M3  )
              Q = Q+fct1*FI (2*M1-1,2*M2  ,2*M3  )+fct4/FI2(2*M1-1,2*M2  ,2*M3  )
              Q = Q+fct1*FI (2*M1  ,2*M2  ,2*M3  )+fct4/FI2(2*M1  ,2*M2  ,2*M3  )
              ! prepare indices of neighbour cells on the 3-point
              ! stencil and apply the periodic boundary condition
              M1u=2*M1+1; IF(M1u>NGRID) M1u=1
              M1l=2*M1-2; IF(M1l<1    ) M1l=NGRID
              M2u=2*M2+1; IF(M2u>NGRID) M2u=1
              M2l=2*M2-2; IF(M2l<1    ) M2l=NGRID
              M3u=2*M3+1; IF(M3u>NGRID) M3u=1
              M3l=2*M3-2; IF(M3l<1    ) M3l=NGRID
              IF(fr_n.EQ.0) THEN
                 ! accumulate contributions to the Laplancian of PDE 
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
              ELSE IF(fr_n.EQ.1) THEN
                 ! accumulate contributions to the Laplancian of PDE 
                 ! from all 8 son cells
                 P = P-3.0D0*FI2(2*M1-1,2*M2-1,2*M3-1)**2      
                 P = P-3.0D0*FI2(2*M1  ,2*M2-1,2*M3-1)**2      
                 P = P-3.0D0*FI2(2*M1-1,2*M2  ,2*M3-1)**2      
                 P = P-3.0D0*FI2(2*M1  ,2*M2  ,2*M3-1)**2      
                 P = P-3.0D0*FI2(2*M1-1,2*M2-1,2*M3  )**2      
                 P = P-3.0D0*FI2(2*M1  ,2*M2-1,2*M3  )**2      
                 P = P-3.0D0*FI2(2*M1-1,2*M2  ,2*M3  )**2      
                 P = P-3.0D0*FI2(2*M1  ,2*M2  ,2*M3  )**2      
                 ! accumulate further contributions to Laplacian PDE
                 ! (1) contributions by 4 cells to the left
                 P = P+FI2(M1l,2*M2-1,2*M3-1)**2      
                 P = P+FI2(M1l,2*M2  ,2*M3-1)**2      
                 P = P+FI2(M1l,2*M2-1,2*M3  )**2      
                 P = P+FI2(M1l,2*M2  ,2*M3  )**2      
                 ! (2) contributions by 4 cells to the right
                 P = P+FI2(M1u,2*M2-1,2*M3-1)**2      
                 P = P+FI2(M1u,2*M2  ,2*M3-1)**2      
                 P = P+FI2(M1u,2*M2-1,2*M3  )**2      
                 P = P+FI2(M1u,2*M2  ,2*M3  )**2      
                 ! (3) contributions by 4 cells behind
                 P = P+FI2(2*M1-1,M2l,2*M3-1)**2      
                 P = P+FI2(2*M1  ,M2l,2*M3-1)**2      
                 P = P+FI2(2*M1-1,M2l,2*M3  )**2      
                 P = P+FI2(2*M1  ,M2l,2*M3  )**2      
                 ! (4) contributions by 4 cells in front
                 P = P+FI2(2*M1-1,M2u,2*M3-1)**2      
                 P = P+FI2(2*M1  ,M2u,2*M3-1)**2      
                 P = P+FI2(2*M1-1,M2u,2*M3  )**2      
                 P = P+FI2(2*M1  ,M2u,2*M3  )**2      
                 ! (5) contributions by 4 cells below
                 P = P+FI2(2*M1-1,2*M2-1,M3l)**2      
                 P = P+FI2(2*M1  ,2*M2-1,M3l)**2      
                 P = P+FI2(2*M1-1,2*M2  ,M3l)**2      
                 P = P+FI2(2*M1  ,2*M2  ,M3l)**2      
                 ! (6) contributions by 4 cells above
                 P = P+FI2(2*M1-1,2*M2-1,M3u)**2      
                 P = P+FI2(2*M1  ,2*M2-1,M3u)**2      
                 P = P+FI2(2*M1-1,2*M2  ,M3u)**2      
                 P = P+FI2(2*M1  ,2*M2  ,M3u)**2      
              ELSE IF(fr_n.EQ.2) THEN
                 ! accumulate contributions to the Laplancian of PDE 
                 ! from all 8 son cells
                 P = P-3.0D0*FI2(2*M1-1,2*M2-1,2*M3-1)**3      
                 P = P-3.0D0*FI2(2*M1  ,2*M2-1,2*M3-1)**3      
                 P = P-3.0D0*FI2(2*M1-1,2*M2  ,2*M3-1)**3      
                 P = P-3.0D0*FI2(2*M1  ,2*M2  ,2*M3-1)**3      
                 P = P-3.0D0*FI2(2*M1-1,2*M2-1,2*M3  )**3      
                 P = P-3.0D0*FI2(2*M1  ,2*M2-1,2*M3  )**3      
                 P = P-3.0D0*FI2(2*M1-1,2*M2  ,2*M3  )**3      
                 P = P-3.0D0*FI2(2*M1  ,2*M2  ,2*M3  )**3      
                 ! accumulate further contributions to Laplacian PDE
                 ! (1) contributions by 4 cells to the left
                 P = P+FI2(M1l,2*M2-1,2*M3-1)**3      
                 P = P+FI2(M1l,2*M2  ,2*M3-1)**3      
                 P = P+FI2(M1l,2*M2-1,2*M3  )**3      
                 P = P+FI2(M1l,2*M2  ,2*M3  )**3      
                 ! (2) contributions by 4 cells to the right
                 P = P+FI2(M1u,2*M2-1,2*M3-1)**3      
                 P = P+FI2(M1u,2*M2  ,2*M3-1)**3      
                 P = P+FI2(M1u,2*M2-1,2*M3  )**3      
                 P = P+FI2(M1u,2*M2  ,2*M3  )**3      
                 ! (3) contributions by 4 cells behind
                 P = P+FI2(2*M1-1,M2l,2*M3-1)**3      
                 P = P+FI2(2*M1  ,M2l,2*M3-1)**3      
                 P = P+FI2(2*M1-1,M2l,2*M3  )**3      
                 P = P+FI2(2*M1  ,M2l,2*M3  )**3      
                 ! (4) contributions by 4 cells in front
                 P = P+FI2(2*M1-1,M2u,2*M3-1)**3      
                 P = P+FI2(2*M1  ,M2u,2*M3-1)**3      
                 P = P+FI2(2*M1-1,M2u,2*M3  )**3      
                 P = P+FI2(2*M1  ,M2u,2*M3  )**3      
                 ! (5) contributions by 4 cells below
                 P = P+FI2(2*M1-1,2*M2-1,M3l)**3      
                 P = P+FI2(2*M1  ,2*M2-1,M3l)**3      
                 P = P+FI2(2*M1-1,2*M2  ,M3l)**3      
                 P = P+FI2(2*M1  ,2*M2  ,M3l)**3      
                 ! (6) contributions by 4 cells above
                 P = P+FI2(2*M1-1,2*M2-1,M3u)**3      
                 P = P+FI2(2*M1  ,2*M2-1,M3u)**3      
                 P = P+FI2(2*M1-1,2*M2  ,M3u)**3      
                 P = P+FI2(2*M1  ,2*M2  ,M3u)**3      
              ENDIF
              !
              ! finalise residuals from all 8 son cells
              P = P*fct3+Q 
              !
              FI3(M1,M2,M3) = P/8.0D0  
           END DO
        END DO
     END DO
  ELSE
!!   ioffset  = 2**(levelmax-ilevel)                                            ! Baojiu: changed ioffset 11-03-2021
!    ioffset  = 2**(levelmax-ilevel+1)                                          ! Baojiu: changed ioffset 11-03-2021
!    koffset  = 2**(levelmax-ilevel)*(2**ilevel-2)
!!   koffset2 = 2**(levelmax-ilevel+1)*(2**(ilevel-1)-2)                        ! Baojiu: changed koffset2 11-03-2021
!    koffset2 = 2**(levelmax-ilevel+1)*(2**(ilevel-1)-1)                        ! Baojiu: changed koffset2 11-03-2021
     ioffset  = NGRID/2**(ilevel-1)                                             ! Baojiu 2025 12 17
     koffset  = NGRID/2**(ilevel  )*(2**(ilevel  )-2)                           ! Baojiu 2025 12 17
     koffset2 = NGRID/2**(ilevel-1)*(2**(ilevel-1)-1)                           ! Baojiu 2025 12 17
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3,M1u,M1l,M2u,M2l,M3u,M3l) &
!$OMP PRIVATE (P,Q) 
     DO M3=1,ngrid_level
        DO M2=1,ngrid_level
           DO M1=1,ngrid_level
              P =   FI3(ioffset+2*M1-1,2*M2-1,koffset2+2*M3-1)                  ! FI3 Section 6: fine level residuals
              P = P+FI3(ioffset+2*M1  ,2*M2-1,koffset2+2*M3-1)                  ! FI3 Section 6: fine level residuals
              P = P+FI3(ioffset+2*M1-1,2*M2  ,koffset2+2*M3-1)                  ! FI3 Section 6: fine level residuals
              P = P+FI3(ioffset+2*M1  ,2*M2  ,koffset2+2*M3-1)                  ! FI3 Section 6: fine level residuals
              P = P+FI3(ioffset+2*M1-1,2*M2-1,koffset2+2*M3  )                  ! FI3 Section 6: fine level residuals
              P = P+FI3(ioffset+2*M1  ,2*M2-1,koffset2+2*M3  )                  ! FI3 Section 6: fine level residuals
              P = P+FI3(ioffset+2*M1-1,2*M2  ,koffset2+2*M3  )                  ! FI3 Section 6: fine level residuals
              P = P+FI3(ioffset+2*M1  ,2*M2  ,koffset2+2*M3  )                  ! FI3 Section 6: fine level residuals
              !
              FI3(M1,M2,koffset+M3) = P/8.0D0                                   ! FI3 Section 1: restricted residuals
              !
           END DO
        END DO
     END DO
  END IF
    
  CALL TimingMain(3,1)
END SUBROUTINE restrict_residual_fR


!-------------------------------------------------------------
!
! Calculate the physical RHS of discrete PDE on coarser ilevel
!
!-------------------------------------------------------------
SUBROUTINE calculate_physical_right_hand_side_fR(ilevel)
!-------------------------------------------------------------
use Tools

  integer :: ilevel
  real*8  :: R0_bg,R_bg,fR_bg,ctilde,dx,dx2
  real*8  :: fct1,fct2,fct3   
  integer :: ngrid_level
  integer :: ioffset,joffset,koffest,koffset2

  integer :: M1,M2,M3,M1l,M2l,M3l,M1u,M2u,M3u
  real*8  :: OP
  real*8  :: s_fR_bg,a_fR_bg,a_fR_rt                            ! sign and amplitude of f_R(a)
  integer :: fr_n1

  IF(MG_test .EQ. 1) WRITE(*,'(A,I5)') 'Calculate physical right-hand side on level',levelmax-ilevel

  ! number of grid points on the coarse level
  ngrid_level = NGRID/2**ilevel

  ctilde  = 2.99792458D3*DBLE(NGRID)/Box
  dx      = DBLE(2**ilevel) 
  dx2     = dx*dx
  R_bg    = 3.0d0*(Om/AEXPN**3+4.0d0*OmL)
  R0_bg   = 3.0d0*(Om         +4.0d0*OmL)
  fR_bg   = fR0*(R0_bg/R_bg)**(fr_n+1)

  s_fR_bg = DSIGN(1.0D0,fR_bg)
  a_fR_bg = DABS (      fR_bg)
  fr_n1   = fr_n+1
  a_fR_rt = a_fR_bg**(1.0D0/(DBLE(fr_n1)))

  fct1   = Om/AEXPN
  fct2   = -R_bg*AEXPN**2/3.0d0
  fct3   = ctilde**2*s_fR_bg/dx2

! ioffset  = 2**(levelmax-ilevel)
! joffset  = 2**(levelmax-ilevel)
! koffset  = 2**(levelmax-ilevel)*(2**ilevel-2)
! koffset2 = 2**(levelmax-ilevel)*(2**ilevel-1)
  ioffset  = NGRID/2**(ilevel)                                                  ! Baojiu 2025 12 17
  joffset  = NGRID/2**(ilevel)                                                  ! Baojiu 2025 12 17
  koffset  = NGRID/2**(ilevel)*(2**ilevel-2)                                    ! Baojiu 2025 12 17
  koffset2 = NGRID/2**(ilevel)*(2**ilevel-1)                                    ! Baojiu 2025 12 17

!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3,M1u,M1l,M2u,M2l,M3u,M3l) &
!$OMP PRIVATE (OP) 
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
           ! calculate Laplancian of PDE using RESTRICTED scalar field (Section 3 of Array FI3)
           IF(fr_n.EQ.0) THEN
              OP = FI3(M1u,joffset+M2 ,koffset+M3 )         + &
                 & FI3(M1l,joffset+M2 ,koffset+M3 )         + &
                 & FI3(M1 ,joffset+M2u,koffset+M3 )         + &
                 & FI3(M1 ,joffset+M2l,koffset+M3 )         + &
                 & FI3(M1 ,joffset+M2 ,koffset+M3u)         + &
                 & FI3(M1 ,joffset+M2 ,koffset+M3l)         - &
                 & FI3(M1 ,joffset+M2 ,koffset+M3 )*6.0D0 
           ELSE IF(fr_n.EQ.1) THEN
              OP = FI3(M1u,joffset+M2 ,koffset+M3 )**2      + &
                 & FI3(M1l,joffset+M2 ,koffset+M3 )**2      + &
                 & FI3(M1 ,joffset+M2u,koffset+M3 )**2      + &
                 & FI3(M1 ,joffset+M2l,koffset+M3 )**2      + &
                 & FI3(M1 ,joffset+M2 ,koffset+M3u)**2      + &
                 & FI3(M1 ,joffset+M2 ,koffset+M3l)**2      - &
                 & FI3(M1 ,joffset+M2 ,koffset+M3 )**2*6.0D0 
           ELSE IF(fr_n.EQ.2) THEN
              OP = FI3(M1u,joffset+M2 ,koffset+M3 )**3      + &
                 & FI3(M1l,joffset+M2 ,koffset+M3 )**3      + &
                 & FI3(M1 ,joffset+M2u,koffset+M3 )**3      + &
                 & FI3(M1 ,joffset+M2l,koffset+M3 )**3      + &
                 & FI3(M1 ,joffset+M2 ,koffset+M3u)**3      + &
                 & FI3(M1 ,joffset+M2 ,koffset+M3l)**3      - &
                 & FI3(M1 ,joffset+M2 ,koffset+M3 )**3*6.0D0 
           ENDIF
           !
           OP = OP*fct3
           OP = OP+fct1*FI3(M1+ioffset,M2,M3+koffset)                           ! add density contribution (Section 2)
           OP = OP+fct2*(a_fR_rt/FI3(M1,M2+joffset,M3+koffset)-1.0D0)           ! add term dependent on restricted scalar field (Section 3)
           OP = OP-FI3(M1,M2,M3+koffset)                                        ! add restricted residual (Section 1); note the minus sign!!!
           ! 
           FI3(M1,joffset+M2,koffset2+M3) = OP                                  ! coarse-level physical RHS in Section 7
           ! 
        END DO
     END DO
  END DO

  CALL TimingMain(3,1)

END SUBROUTINE calculate_physical_right_hand_side_fR
