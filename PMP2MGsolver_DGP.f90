!=================================================
!
!  DGP model case
!
!=================================================

!-------------------------------------------------------------
SUBROUTINE relaxation_iterations_DGP(ilevel,redstep)
!-------------------------------------------------------------
use Tools

  integer :: ilevel
  logical :: redstep
  real*8  :: alpha,beta,Orc,Rc_sq,dx,dx2,dx4
  real*8  :: fct1
  real*8  :: TWOOVERTHREE,ONEOVEREIGHT,THREEOVERFOUR,EIGHTOVERTHREE
  integer :: nstart,nfinal,ngrid_level
  integer :: ioffset,joffset,koffset,koffset2

  real*8  :: L,dLphi,Sigma!,tmp_sum
  integer :: M1,M2,M3,M1u,M1l,M2u,M2l,M3u,M3l
 
  ! constants
  ONEOVEREIGHT   = 0.125D0
  TWOOVERTHREE   = 2.0D0/3.0D0
  THREEOVERFOUR  = 0.75D0
  EIGHTOVERTHREE = 8.0D0/3.0D0

  ! number of grid points on the coarse level
  ngrid_level = NGRID/(2**ilevel)

  ! simulation parameters
  dx     = DBLE(2**ilevel) 
  dx2    = dx*dx
  dx4    = dx2*dx2

  ! physical and numerical quantities
  Orc    = 1.0D0/(4.0D0*H0rc**2)
  !
  IF(N_branch) beta = 1.0D0 + (0.5D0*Om/AEXPN**3+OmL) / (DSQRT(Orc*(Om/AEXPN**3+OmL))) ! normal branch
  IF(S_branch) beta =       - (0.5D0*Om/AEXPN**3+Orc) / (DSQRT(Orc*(Om/AEXPN**3+Orc))) ! self-accelerated branch
  !
  Rc_sq  = 1.0D0/(4.0D0*Orc)
  alpha  = 3.0D0*beta*AEXPN**2/Rc_sq
  fct1   = alpha*Om/beta/AEXPN
  dLphi  = -6.0D0/dx2

!  IF(MG_test) THEN
!     WRITE(*,'(A, I5)'    ) 'ngrid_level = ', ngrid_level
!     WRITE(*,'(A7,F20.14)') 'dx = ', dx
!     WRITE(*,'(A7,F20.14)') 'Orc = ', Orc
!     WRITE(*,'(A7,F20.14)') 'Rc_sq = ', Rc_sq
!     WRITE(*,'(A7,F20.14)') 'beta = ', beta
!     WRITE(*,'(A7,F20.14)') 'alpha = ', alpha
!     WRITE(*,'(A7,F20.14)') 'fct1 = ', fct1
!     WRITE(*,'(A7,F20.14)') 'dLphi = ', dLphi
!  ENDIF


  ! offsets for cell access to array FI3
  ioffset  = 2**(levelmax-ilevel)
  joffset  = 2**(levelmax-ilevel)
  koffset  = 0
  koffset2 = 0
  IF(ilevel.GT.0) THEN
     koffset  = 2**(levelmax-ilevel)*(2**ilevel-2)
     koffset2 = 2**(levelmax-ilevel)*(2**ilevel-1)
  END IF

!-------------------------------------------------------------
! Baojiu-13-06-2021: modification to calculate S_mean -- start
!-------------------------------------------------------------

!  tmp_sum = 0.0D0

!!$OMP PARALLEL DO DEFAULT(SHARED) &
!!$OMP PRIVATE (M1,M2,M3,M1u,M1l,M2u,M2l,M3u,M3l) &
!!$OMP PRIVATE (Sigma) REDUCTION(+:tmp_sum)
!  DO M3=1,ngrid_level
!     DO M2=1,ngrid_level
!        DO M1=1,ngrid_level
           ! prepare the indices of the 6 neighbours on the 3-point stencil
           ! and apply periodic boundary condition
!           M1u = M1+1; IF(M1u>ngrid_level) M1u = 1
!           M1l = M1-1; IF(M1l<1          ) M1l = ngrid_level
!           M2u = M2+1; IF(M2u>ngrid_level) M2u = 1
!           M2l = M2-1; IF(M2l<1          ) M2l = ngrid_level
!           M3u = M3+1; IF(M3u>ngrid_level) M3u = 1
!           M3l = M3-1; IF(M3l<1          ) M3l = ngrid_level
           ! decide on red-black update ordering
!           IF((     redstep).AND.(MOD(M1+M2+M3,2).EQ.0)) CYCLE
!           IF((.NOT.redstep).AND.(MOD(M1+M2+M3,2).NE.0)) CYCLE
           !
!           IF(ilevel.EQ.0) THEN
!              Sigma = fct1*FI (M1        ,M2,M3        )
!           ELSE
!              Sigma = fct1*FI3(M1+ioffset,M2,M3+koffset)
!           END IF
           ! summation of \phi of all 6 neighbors on the 3-point stencil
!           IF(ilevel.EQ.0) THEN
!              Sigma = Sigma                                                                               + &
!                    & ( (FI2(M1u,M2,M3) + FI2(M1l,M2,M3) - 2.0D0*FI2(M1,M2,M3))**2                        + &
!                    &   (FI2(M1,M2u,M3) + FI2(M1,M2l,M3) - 2.0D0*FI2(M1,M2,M3))**2                        + &
!                    &   (FI2(M1,M2,M3u) + FI2(M1,M2,M3l) - 2.0D0*FI2(M1,M2,M3))**2 ) * TWOOVERTHREE/dx4   - &
!                    & ( (FI2(M1u,M2,M3) + FI2(M1l,M2,M3) - 2.0D0*FI2(M1,M2,M3))                           * &
!                    &   (FI2(M1,M2u,M3) + FI2(M1,M2l,M3) - 2.0D0*FI2(M1,M2,M3))                           + &
!                    &   (FI2(M1u,M2,M3) + FI2(M1l,M2,M3) - 2.0D0*FI2(M1,M2,M3))                           * &
!                    &   (FI2(M1,M2,M3u) + FI2(M1,M2,M3l) - 2.0D0*FI2(M1,M2,M3))                           + &
!                    &   (FI2(M1,M2u,M3) + FI2(M1,M2l,M3) - 2.0D0*FI2(M1,M2,M3))                           * &
!                    &   (FI2(M1,M2,M3u) + FI2(M1,M2,M3l) - 2.0D0*FI2(M1,M2,M3))    ) * TWOOVERTHREE/dx4   + &
!                    & ( (FI2(M1u,M2u,M3) + FI2(M1l,M2l,M3) - FI2(M1u,M2l,M3) - FI2(M1l,M2u,M3))**2        + &
!                    &   (FI2(M1u,M2,M3u) + FI2(M1l,M2,M3l) - FI2(M1u,M2,M3l) - FI2(M1l,M2,M3u))**2        + &
!                    &   (FI2(M1,M2u,M3u) + FI2(M1,M2l,M3l) - FI2(M1,M2u,M3l) - FI2(M1,M2l,M3u))**2 ) * ONEOVEREIGHT/dx4
              !
!              tmp_sum = tmp_sum + (DSIGN(1.0D0,alpha)*DSQRT(alpha**2+EIGHTOVERTHREE*Sigma)-alpha)*THREEOVERFOUR
!           ELSE
!              Sigma = Sigma                                                                                                                                                         + &
!                    & ( (FI3(M1u+ioffset,M2 +joffset,M3 +koffset) + FI3(M1l+ioffset,M2 +joffset,M3 +koffset) - 2.0D0*FI3(M1+ioffset,M2+joffset,M3+koffset))**2                      + &
!                    &   (FI3(M1 +ioffset,M2u+joffset,M3 +koffset) + FI3(M1 +ioffset,M2l+joffset,M3 +koffset) - 2.0D0*FI3(M1+ioffset,M2+joffset,M3+koffset))**2                      + &
!                    &   (FI3(M1 +ioffset,M2 +joffset,M3u+koffset) + FI3(M1 +ioffset,M2 +joffset,M3l+koffset) - 2.0D0*FI3(M1+ioffset,M2+joffset,M3+koffset))**2 ) * TWOOVERTHREE/dx4 - &
!                    & ( (FI3(M1u+ioffset,M2 +joffset,M3 +koffset) + FI3(M1l+ioffset,M2 +joffset,M3 +koffset) - 2.0D0*FI3(M1+ioffset,M2+joffset,M3+koffset))                         * &
!                    &   (FI3(M1 +ioffset,M2u+joffset,M3 +koffset) + FI3(M1 +ioffset,M2l+joffset,M3 +koffset) - 2.0D0*FI3(M1+ioffset,M2+joffset,M3+koffset))                         + &
!                    &   (FI3(M1u+ioffset,M2 +joffset,M3 +koffset) + FI3(M1l+ioffset,M2 +joffset,M3 +koffset) - 2.0D0*FI3(M1+ioffset,M2+joffset,M3+koffset))                         * &
!                    &   (FI3(M1 +ioffset,M2 +joffset,M3u+koffset) + FI3(M1 +ioffset,M2 +joffset,M3l+koffset) - 2.0D0*FI3(M1+ioffset,M2+joffset,M3+koffset))                         + &
!                    &   (FI3(M1 +ioffset,M2u+joffset,M3 +koffset) + FI3(M1 +ioffset,M2l+joffset,M3 +koffset) - 2.0D0*FI3(M1+ioffset,M2+joffset,M3+koffset))                         * &
!                    &   (FI3(M1 +ioffset,M2 +joffset,M3u+koffset) + FI3(M1 +ioffset,M2 +joffset,M3l+koffset) - 2.0D0*FI3(M1+ioffset,M2+joffset,M3+koffset)) )    * TWOOVERTHREE/dx4 + &
!                    & ( (FI3(M1u+ioffset,M2u+joffset,M3 +koffset) + FI3(M1l+ioffset,M2l+joffset,M3 +koffset)                                                                        - &
!                    &    FI3(M1u+ioffset,M2l+joffset,M3 +koffset) - FI3(M1l+ioffset,M2u+joffset,M3 +koffset) )**2                                                                   + &
!                    &   (FI3(M1u+ioffset,M2 +joffset,M3u+koffset) + FI3(M1l+ioffset,M2 +joffset,M3l+koffset)                                                                        - &   
!                    &    FI3(M1u+ioffset,M2 +joffset,M3l+koffset) - FI3(M1l+ioffset,M2 +joffset,M3u+koffset) )**2                                                                   + &
!                    &   (FI3(M1 +ioffset,M2u+joffset,M3u+koffset) + FI3(M1 +ioffset,M2l+joffset,M3l+koffset)                                                                        - &
!                    &    FI3(M1 +ioffset,M2u+joffset,M3l+koffset) - FI3(M1 +ioffset,M2l+joffset,M3u+koffset) )**2 ) * ONEOVEREIGHT/dx4
              !
!              tmp_sum = tmp_sum + (DSIGN(1.0D0,alpha)*DSQRT(alpha**2+EIGHTOVERTHREE*Sigma)-alpha)*THREEOVERFOUR
!           END IF
!        END DO
!     END DO
!  END DO

!  S_mean = tmp_sum/(DBLE(ngrid_level))**3

!-----------------------------------------------------------
! Baojiu-13-06-2021: modification to calculate S_mean -- end
!-----------------------------------------------------------

!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3,M1u,M1l,M2u,M2l,M3u,M3l) &
!$OMP PRIVATE (L,Sigma)
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
              Sigma = fct1*FI (M1        ,M2,M3        )
           ELSE
              Sigma = fct1*FI3(M1+ioffset,M2,M3+koffset)
           END IF
           ! summation of \phi of all 6 neighbors on the 3-point stencil
           IF(ilevel.EQ.0) THEN
              Sigma = Sigma                                                                               + &
                    & ( (FI2(M1u,M2,M3) + FI2(M1l,M2,M3) - 2.0D0*FI2(M1,M2,M3))**2                        + &
                    &   (FI2(M1,M2u,M3) + FI2(M1,M2l,M3) - 2.0D0*FI2(M1,M2,M3))**2                        + &
                    &   (FI2(M1,M2,M3u) + FI2(M1,M2,M3l) - 2.0D0*FI2(M1,M2,M3))**2 ) * TWOOVERTHREE/dx4   - &
                    & ( (FI2(M1u,M2,M3) + FI2(M1l,M2,M3) - 2.0D0*FI2(M1,M2,M3))                           * &
                    &   (FI2(M1,M2u,M3) + FI2(M1,M2l,M3) - 2.0D0*FI2(M1,M2,M3))                           + &
                    &   (FI2(M1u,M2,M3) + FI2(M1l,M2,M3) - 2.0D0*FI2(M1,M2,M3))                           * &
                    &   (FI2(M1,M2,M3u) + FI2(M1,M2,M3l) - 2.0D0*FI2(M1,M2,M3))                           + &
                    &   (FI2(M1,M2u,M3) + FI2(M1,M2l,M3) - 2.0D0*FI2(M1,M2,M3))                           * &
                    &   (FI2(M1,M2,M3u) + FI2(M1,M2,M3l) - 2.0D0*FI2(M1,M2,M3))    ) * TWOOVERTHREE/dx4   + &
                    & ( (FI2(M1u,M2u,M3) + FI2(M1l,M2l,M3) - FI2(M1u,M2l,M3) - FI2(M1l,M2u,M3))**2        + &
                    &   (FI2(M1u,M2,M3u) + FI2(M1l,M2,M3l) - FI2(M1u,M2,M3l) - FI2(M1l,M2,M3u))**2        + &
                    &   (FI2(M1,M2u,M3u) + FI2(M1,M2l,M3l) - FI2(M1,M2u,M3l) - FI2(M1,M2l,M3u))**2 ) * ONEOVEREIGHT/dx4
              !
              L     = (  FI2(M1u,M2,M3)              + &
                    &    FI2(M1l,M2,M3)              + &
                    &    FI2(M1,M2u,M3)              + &
                    &    FI2(M1,M2l,M3)              + &
                    &    FI2(M1,M2,M3u)              + &
                    &    FI2(M1,M2,M3l)              - &
                    &    FI2(M1,M2,M3 )*6.0D0 )/dx2  - &
                    & ( DSIGN(1.0D0,alpha)*DSQRT(alpha**2+EIGHTOVERTHREE*Sigma) - alpha) * THREEOVERFOUR 
!              L     = L + S_mean ! Baojiu-13-06-2021 added S_mean contribution
           ELSE
              Sigma = Sigma                                                                                                                                                         + &
                    & ( (FI3(M1u+ioffset,M2 +joffset,M3 +koffset) + FI3(M1l+ioffset,M2 +joffset,M3 +koffset) - 2.0D0*FI3(M1+ioffset,M2+joffset,M3+koffset))**2                      + &
                    &   (FI3(M1 +ioffset,M2u+joffset,M3 +koffset) + FI3(M1 +ioffset,M2l+joffset,M3 +koffset) - 2.0D0*FI3(M1+ioffset,M2+joffset,M3+koffset))**2                      + &
                    &   (FI3(M1 +ioffset,M2 +joffset,M3u+koffset) + FI3(M1 +ioffset,M2 +joffset,M3l+koffset) - 2.0D0*FI3(M1+ioffset,M2+joffset,M3+koffset))**2 ) * TWOOVERTHREE/dx4 - &
                    & ( (FI3(M1u+ioffset,M2 +joffset,M3 +koffset) + FI3(M1l+ioffset,M2 +joffset,M3 +koffset) - 2.0D0*FI3(M1+ioffset,M2+joffset,M3+koffset))                         * &
                    &   (FI3(M1 +ioffset,M2u+joffset,M3 +koffset) + FI3(M1 +ioffset,M2l+joffset,M3 +koffset) - 2.0D0*FI3(M1+ioffset,M2+joffset,M3+koffset))                         + &
                    &   (FI3(M1u+ioffset,M2 +joffset,M3 +koffset) + FI3(M1l+ioffset,M2 +joffset,M3 +koffset) - 2.0D0*FI3(M1+ioffset,M2+joffset,M3+koffset))                         * &
                    &   (FI3(M1 +ioffset,M2 +joffset,M3u+koffset) + FI3(M1 +ioffset,M2 +joffset,M3l+koffset) - 2.0D0*FI3(M1+ioffset,M2+joffset,M3+koffset))                         + &
                    &   (FI3(M1 +ioffset,M2u+joffset,M3 +koffset) + FI3(M1 +ioffset,M2l+joffset,M3 +koffset) - 2.0D0*FI3(M1+ioffset,M2+joffset,M3+koffset))                         * &
                    &   (FI3(M1 +ioffset,M2 +joffset,M3u+koffset) + FI3(M1 +ioffset,M2 +joffset,M3l+koffset) - 2.0D0*FI3(M1+ioffset,M2+joffset,M3+koffset)) )    * TWOOVERTHREE/dx4 + &
                    & ( (FI3(M1u+ioffset,M2u+joffset,M3 +koffset) + FI3(M1l+ioffset,M2l+joffset,M3 +koffset)                                                                        - &
                    &    FI3(M1u+ioffset,M2l+joffset,M3 +koffset) - FI3(M1l+ioffset,M2u+joffset,M3 +koffset) )**2                                                                   + &
                    &   (FI3(M1u+ioffset,M2 +joffset,M3u+koffset) + FI3(M1l+ioffset,M2 +joffset,M3l+koffset)                                                                        - &   
                    &    FI3(M1u+ioffset,M2 +joffset,M3l+koffset) - FI3(M1l+ioffset,M2 +joffset,M3u+koffset) )**2                                                                   + &
                    &   (FI3(M1 +ioffset,M2u+joffset,M3u+koffset) + FI3(M1 +ioffset,M2l+joffset,M3l+koffset)                                                                        - &
                    &    FI3(M1 +ioffset,M2u+joffset,M3l+koffset) - FI3(M1 +ioffset,M2l+joffset,M3u+koffset) )**2 ) * ONEOVEREIGHT/dx4
              !
              L     = (  FI3(M1u+ioffset,M2 +joffset,M3 +koffset)               + &
                    &    FI3(M1l+ioffset,M2 +joffset,M3 +koffset)               + &
                    &    FI3(M1 +ioffset,M2u+joffset,M3 +koffset)               + &
                    &    FI3(M1 +ioffset,M2l+joffset,M3 +koffset)               + &
                    &    FI3(M1 +ioffset,M2 +joffset,M3u+koffset)               + &
                    &    FI3(M1 +ioffset,M2 +joffset,M3l+koffset)               - &
                    &    FI3(M1 +ioffset,M2 +joffset,M3 +koffset)*6.0D0 )/dx2   - &
                    & ( DSIGN(1.0D0,alpha)*DSQRT(alpha**2+EIGHTOVERTHREE*Sigma) - alpha) * THREEOVERFOUR
              L     = L - FI3(M1,M2+joffset,M3+koffset2)
!              L     = L + S_mean ! Baojiu-13-06-2021 added S_mean contribution
           END IF
           ! update solution
           IF(ilevel.EQ.0) THEN      ! solution on PM grid
!$OMP ATOMIC
               FI2(M1        ,M2        ,M3        ) = FI2(M1        ,M2        ,M3        ) - L/dLphi ! phi_old - L/(dL/dphi)
           ELSE                      ! solution on multigrid
!$OMP ATOMIC
               FI3(M1+ioffset,M2+joffset,M3+koffset) = FI3(M1+ioffset,M2+joffset,M3+koffset) - L/dLphi ! phi_old - L/(dL/dphi)
           END IF
        END DO
     END DO
  END DO
    
  CALL TimingMain(3,1)
    
END SUBROUTINE relaxation_iterations_DGP

!-------------------------------------------------------------
SUBROUTINE calculate_residual_DGP(ilevel,res_PM_grid)
!-------------------------------------------------------------
use Tools

  integer,intent(in) :: ilevel
  real*8,intent(out) :: res_PM_grid
  
  real*8  :: alpha,beta,Orc,Rc_sq,dx,dx2,dx4
  integer :: ngrid_level
  integer :: ioffset,joffset,koffset,koffset2
  real*8  :: fct1
  real*8  :: TWOOVERTHREE,ONEOVEREIGHT,THREEOVERFOUR,EIGHTOVERTHREE

  integer :: M1,M2,M3,M1l,M1u,M2l,M2u,M3l,M3u
  integer :: N1,N2,N3
  real*8  :: RES,OP,Sigma

  ! number of grid points on the coarse level
  ngrid_level = NGRID/2**ilevel

  ONEOVEREIGHT   = 1.0D0/8.0D0
  TWOOVERTHREE   = 2.0D0/3.0D0
  THREEOVERFOUR  = 3.0D0/4.0D0
  EIGHTOVERTHREE = 8.0D0/3.0D0

  dx      = DBLE(2**ilevel) 
  dx2     = dx*dx
  dx4     = dx2*dx2
  Orc     = 1.0D0/(4.0D0*H0rc**2)
  !
  IF(N_branch) beta = 1.0D0 + (0.5D0*Om/AEXPN**3+OmL) / (DSQRT(Orc*(Om/AEXPN**3+OmL))) ! normal branch
  IF(S_branch) beta =       - (0.5D0*Om/AEXPN**3+Orc) / (DSQRT(Orc*(Om/AEXPN**3+Orc))) ! self-accelerated branch
  !
  Rc_sq   = 1.0 / (4.0 * Orc)
  alpha   = 3.0 * beta * AEXPN**2 / Rc_sq
  fct1    = alpha * Om / beta / AEXPN
 
  ioffset  = 2**(levelmax-ilevel)
  joffset  = 2**(levelmax-ilevel)
  koffset  = 0
  koffset2 = 0
  IF(ilevel.GT.0) THEN
     koffset  = 2**(levelmax-ilevel)*(2**ilevel-2)
     koffset2 = 2**(levelmax-ilevel)*(2**ilevel-1)
  END IF

!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3,M1u,M1l,M2u,M2l,M3u,M3l) &
!$OMP PRIVATE (OP,Sigma)
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
              Sigma = fct1*FI(M1,M2,M3)
              Sigma = Sigma                                                                             + &
                    & ( (FI2(M1u,M2,M3) + FI2(M1l,M2,M3) - 2.0D0*FI2(M1,M2,M3))**2                      + &
                    &   (FI2(M1,M2u,M3) + FI2(M1,M2l,M3) - 2.0D0*FI2(M1,M2,M3))**2                      + &
                    &   (FI2(M1,M2,M3u) + FI2(M1,M2,M3l) - 2.0D0*FI2(M1,M2,M3))**2 ) * TWOOVERTHREE/dx4 - &
                    & ( (FI2(M1u,M2,M3) + FI2(M1l,M2,M3) - 2.0D0*FI2(M1,M2,M3))                         * &
                    &   (FI2(M1,M2u,M3) + FI2(M1,M2l,M3) - 2.0D0*FI2(M1,M2,M3))                         + &
                    &   (FI2(M1u,M2,M3) + FI2(M1l,M2,M3) - 2.0D0*FI2(M1,M2,M3))                         * &
                    &   (FI2(M1,M2,M3u) + FI2(M1,M2,M3l) - 2.0D0*FI2(M1,M2,M3))                         + &
                    &   (FI2(M1,M2u,M3) + FI2(M1,M2l,M3) - 2.0D0*FI2(M1,M2,M3))                         * &
                    &   (FI2(M1,M2,M3u) + FI2(M1,M2,M3l) - 2.0D0*FI2(M1,M2,M3))    ) * TWOOVERTHREE/dx4 + &
                    & ( (FI2(M1u,M2u,M3) + FI2(M1l,M2l,M3) - FI2(M1u,M2l,M3) - FI2(M1l,M2u,M3))**2      + &
                    &   (FI2(M1u,M2,M3u) + FI2(M1l,M2,M3l) - FI2(M1u,M2,M3l) - FI2(M1l,M2,M3u))**2      + &
                    &   (FI2(M1,M2u,M3u) + FI2(M1,M2l,M3l) - FI2(M1,M2u,M3l) - FI2(M1,M2l,M3u))**2 ) * ONEOVEREIGHT/dx4
              OP    = (  FI2(M1u,M2,M3)              + &
                    &    FI2(M1l,M2,M3)              + &
                    &    FI2(M1,M2u,M3)              + &
                    &    FI2(M1,M2l,M3)              + &
                    &    FI2(M1,M2,M3u)              + &
                    &    FI2(M1,M2,M3l)              - &
                    &    FI2(M1,M2,M3 )*6.0D0 )/dx2
              OP    = OP - (DSIGN(1.0D0,alpha)*DSQRT(alpha**2+EIGHTOVERTHREE*Sigma) - alpha) * THREEOVERFOUR
!              OP    = OP + S_mean ! Baojiu-13-06-2021 added S_mean contribution
              !
              FI3(M1,M2,M3) = OP
           ELSE
              ! calculate Sigma
              Sigma = fct1*FI3(M1+ioffset,M2,M3+koffset)
              Sigma = Sigma                                                                                                                                                         + &
                    & ( (FI3(M1u+ioffset,M2 +joffset,M3 +koffset) + FI3(M1l+ioffset,M2 +joffset,M3 +koffset) - 2.0D0*FI3(M1+ioffset,M2+joffset,M3+koffset))**2                      + &
                    &   (FI3(M1 +ioffset,M2u+joffset,M3 +koffset) + FI3(M1 +ioffset,M2l+joffset,M3 +koffset) - 2.0D0*FI3(M1+ioffset,M2+joffset,M3+koffset))**2                      + &
                    &   (FI3(M1 +ioffset,M2 +joffset,M3u+koffset) + FI3(M1 +ioffset,M2 +joffset,M3l+koffset) - 2.0D0*FI3(M1+ioffset,M2+joffset,M3+koffset))**2 ) * TWOOVERTHREE/dx4 - &
                    & ( (FI3(M1u+ioffset,M2 +joffset,M3 +koffset) + FI3(M1l+ioffset,M2 +joffset,M3 +koffset) - 2.0D0*FI3(M1+ioffset,M2+joffset,M3+koffset))                         * &
                    &   (FI3(M1 +ioffset,M2u+joffset,M3 +koffset) + FI3(M1 +ioffset,M2l+joffset,M3 +koffset) - 2.0D0*FI3(M1+ioffset,M2+joffset,M3+koffset))                         + &
                    &   (FI3(M1u+ioffset,M2 +joffset,M3 +koffset) + FI3(M1l+ioffset,M2 +joffset,M3 +koffset) - 2.0D0*FI3(M1+ioffset,M2+joffset,M3+koffset))                         * &
                    &   (FI3(M1 +ioffset,M2 +joffset,M3u+koffset) + FI3(M1 +ioffset,M2 +joffset,M3l+koffset) - 2.0D0*FI3(M1+ioffset,M2+joffset,M3+koffset))                         + &
                    &   (FI3(M1 +ioffset,M2u+joffset,M3 +koffset) + FI3(M1 +ioffset,M2l+joffset,M3 +koffset) - 2.0D0*FI3(M1+ioffset,M2+joffset,M3+koffset))                         * &
                    &   (FI3(M1 +ioffset,M2 +joffset,M3u+koffset) + FI3(M1 +ioffset,M2 +joffset,M3l+koffset) - 2.0D0*FI3(M1+ioffset,M2+joffset,M3+koffset)) )    * TWOOVERTHREE/dx4 + &
                    & ( (FI3(M1u+ioffset,M2u+joffset,M3 +koffset) + FI3(M1l+ioffset,M2l+joffset,M3 +koffset)                                                                        - &
                    &    FI3(M1u+ioffset,M2l+joffset,M3 +koffset) - FI3(M1l+ioffset,M2u+joffset,M3 +koffset) )**2                                                                   + &
                    &   (FI3(M1u+ioffset,M2 +joffset,M3u+koffset) + FI3(M1l+ioffset,M2 +joffset,M3l+koffset)                                                                        - &
                    &    FI3(M1u+ioffset,M2 +joffset,M3l+koffset) - FI3(M1l+ioffset,M2 +joffset,M3u+koffset) )**2                                                                   + &
                    &   (FI3(M1 +ioffset,M2u+joffset,M3u+koffset) + FI3(M1 +ioffset,M2l+joffset,M3l+koffset)                                                                        - &
                    &    FI3(M1 +ioffset,M2u+joffset,M3l+koffset) - FI3(M1 +ioffset,M2l+joffset,M3u+koffset) )**2 ) * ONEOVEREIGHT/dx4
              ! calculate Laplacian of PDE
              OP    = ( FI3(M1u+ioffset,M2 +joffset,M3 +koffset)               + &
                    &   FI3(M1l+ioffset,M2 +joffset,M3 +koffset)               + &
                    &   FI3(M1 +ioffset,M2u+joffset,M3 +koffset)               + &
                    &   FI3(M1 +ioffset,M2l+joffset,M3 +koffset)               + &
                    &   FI3(M1 +ioffset,M2 +joffset,M3u+koffset)               + &
                    &   FI3(M1 +ioffset,M2 +joffset,M3l+koffset)               - &
                    &   FI3(M1 +ioffset,M2 +joffset,M3 +koffset)*6.0D0 ) / dx2
              OP    = OP - (DSIGN(1.0D0,alpha)*DSQRT(alpha**2+EIGHTOVERTHREE*Sigma) - alpha) * THREEOVERFOUR
              OP    = OP - FI3(M1,M2+joffset,M3+koffset2) ! physical right-hand side
!              OP    = OP + S_mean ! Baojiu-13-06-2021 added S_mean contribution
              !
              FI3(M1+ioffset,M2,M3+koffset2) =  OP
           END IF
        END DO
     END DO
  END DO
    
  CALL TimingMain(3,1)

  RES = 0.0D0
  !
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

! WRITE(*,'(A,F20.12,A,I5,A)') 'Residual RMS = ',res_PM_grid,'on level ',levelmax-ilevel  

      
END SUBROUTINE calculate_residual_DGP

!-------------------------------------------------------------
!SUBROUTINE phi_mean_DGP
!-------------------------------------------------------------
!use Tools

!  integer :: M1,M2,M3
!  integer :: Ncells
!  real*8  :: sum

!  Ncells = NGRID**3
!  sum = 0.0D0
  !
!!$OMP PARALLEL DO DEFAULT(SHARED) &
!!$OMP PRIVATE (M1,M2,M3) REDUCTION(+:sum)
!  DO M3=1,NGRID
!     DO M2=1,NGRID
!        DO M1=1,NGRID
!           sum = sum + FI2(M1,M2,M3)
!        END DO
!     END DO
!  END DO

!  phi_mean = sum/DBLE(Ncells)
!  WRITE(*,'(A12,F20.14)') 'phi_mean = ', phi_mean

!END SUBROUTINE phi_mean_DGP

!-------------------------------------------------------------
SUBROUTINE restrict_residual_DGP(ilevel)
!-------------------------------------------------------------
use Tools

  integer :: ilevel
  real*8  :: alpha,beta,Orc,Rc_sq,dx,dx2,dx4
  integer :: ngrid_level
  integer :: ioffset,joffset,koffset,ioffset2,joffset2,koffset2
  real*8  :: fct1,fct2,fct3 
  real*8  :: TWOOVERTHREE,ONEOVEREIGHT,THREEOVERFOUR,EIGHTOVERTHREE,ONEOVERTHREE

  integer :: M1,M2,M3,M1l,M1u,M2l,M2u,M3l,M3u
  real*8  :: P,Q
  real*8  :: Sigma_1,Sigma_2,Sigma_3,Sigma_4,Sigma_5,Sigma_6,Sigma_7,Sigma_8

  ! number of grid points on the coarse level
  ngrid_level = NGRID/2**ilevel

  ONEOVEREIGHT   = 1.0D0/8.0D0
  ONEOVERTHREE   = 1.0D0/3.0D0
  TWOOVERTHREE   = 2.0D0/3.0D0
  THREEOVERFOUR  = 3.0D0/4.0D0
  EIGHTOVERTHREE = 8.0D0/3.0D0

  dx    = 1.0D0 
  dx2   = dx*dx
  dx4   = dx2*dx2
  Orc   = 1.0D0/(4.0D0*H0rc**2)
  !
  IF(N_branch) beta = 1.0D0 + (0.5D0*Om/AEXPN**3+OmL) / (DSQRT(Orc*(Om/AEXPN**3+OmL))) ! normal branch
  IF(S_branch) beta =       - (0.5D0*Om/AEXPN**3+Orc) / (DSQRT(Orc*(Om/AEXPN**3+Orc))) ! self-accelerated branch
  !
  Rc_sq = 1.0D0/(4.0D0*Orc)
  alpha = 3.0D0*beta*AEXPN**2/Rc_sq
  fct1  = alpha*Om/beta/AEXPN
  fct2  = -0.75D0*alpha
! fct3  = 8.0D0*THREEOVERFOUR*DSIGN(1.0D0,alpha)
  fct3  = DSIGN(0.75D0,alpha)

  IF(ilevel.EQ.1) THEN
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3,M1u,M1l,M2u,M2l,M3u,M3l) &
!$OMP PRIVATE (Sigma_1,Sigma_2,Sigma_3,Sigma_4,Sigma_5,Sigma_6,Sigma_7,Sigma_8) &
!$OMP PRIVATE (P,Q)
     DO M3=1,ngrid_level
        DO M2=1,ngrid_level
           DO M1=1,ngrid_level
              P      = 0.0D0
              Q      = -8.0D0*fct2
              !
              Sigma_1 = fct1*FI(2*M1-1,2*M2-1,2*M3-1) 
              Sigma_2 = fct1*FI(2*M1  ,2*M2-1,2*M3-1)
              Sigma_3 = fct1*FI(2*M1-1,2*M2  ,2*M3-1) 
              Sigma_4 = fct1*FI(2*M1  ,2*M2  ,2*M3-1)
              Sigma_5 = fct1*FI(2*M1-1,2*M2-1,2*M3  )
              Sigma_6 = fct1*FI(2*M1  ,2*M2-1,2*M3  )
              Sigma_7 = fct1*FI(2*M1-1,2*M2  ,2*M3  )
              Sigma_8 = fct1*FI(2*M1  ,2*M2  ,2*M3  ) 
              !
              ! prepare indices of neighbour cells on the 3-point
              ! stencil and apply the periodic boundary condition
              M1u=2*M1+1; IF(M1u>NGRID) M1u=1
              M1l=2*M1-2; IF(M1l<1    ) M1l=NGRID
              M2u=2*M2+1; IF(M2u>NGRID) M2u=1
              M2l=2*M2-2; IF(M2l<1    ) M2l=NGRID
              M3u=2*M3+1; IF(M3u>NGRID) M3u=1
              M3l=2*M3-2; IF(M3l<1    ) M3l=NGRID
              ! accumulate contributions to the Sigma term from all 8 son cells
              ! the following is based on Eq.(B1) of 1308.3491
              !
              ! Source term of cell (2*M1-1,2*M2-1,2*M3-1)
              Sigma_1 = Sigma_1                                                                               + &
                     ! direct neighbours in x direction
                     &        (FI2(2*M1  ,2*M2-1,2*M3-1) + FI2(  M1l ,2*M2-1,2*M3-1))                         * &
                     & (2.0D0*(FI2(2*M1  ,2*M2-1,2*M3-1) + FI2(  M1l ,2*M2-1,2*M3-1))                         - & 
                     &        (FI2(2*M1-1,2*M2  ,2*M3-1) + FI2(2*M1-1,  M2l ,2*M3-1))                         - &
                     &        (FI2(2*M1-1,2*M2-1,2*M3  ) + FI2(2*M1-1,2*M2-1,  M3l )) )    * ONEOVERTHREE/dx4 + &
                     ! direct neighbours in y direction
                     &        (FI2(2*M1-1,2*M2  ,2*M3-1) + FI2(2*M1-1,  M2l ,2*M3-1))                         * &
                     & (2.0D0*(FI2(2*M1-1,2*M2  ,2*M3-1) + FI2(2*M1-1,  M2l ,2*M3-1))                         - &
                     &        (FI2(2*M1  ,2*M2-1,2*M3-1) + FI2(  M1l ,2*M2-1,2*M3-1))                         - &
                     &        (FI2(2*M1-1,2*M2-1,2*M3  ) + FI2(2*M1-1,2*M2-1,  M3l )) )    * ONEOVERTHREE/dx4 + &
                     ! direct neighbours in z direction
                     &        (FI2(2*M1-1,2*M2-1,2*M3  ) + FI2(2*M1-1,2*M2-1,  M3l ))                         * &
                     & (2.0D0*(FI2(2*M1-1,2*M2-1,2*M3  ) + FI2(2*M1-1,2*M2-1,  M3l ))                         - &
                     &        (FI2(2*M1  ,2*M2-1,2*M3-1) + FI2(  M1l ,2*M2-1,2*M3-1))                         - &
                     &        (FI2(2*M1-1,2*M2  ,2*M3-1) + FI2(2*M1-1,  M2l ,2*M3-1)) )    * ONEOVERTHREE/dx4 + &
                     ! diagonal neighbours in x-y plane
                     & (      (FI2(2*M1  ,2*M2  ,2*M3-1) + FI2(  M1l ,  M2l ,2*M3-1))                         - &
                     &        (FI2(2*M1  ,  M2l ,2*M3-1) + FI2(  M1l ,2*M2  ,2*M3-1)) )**2 * ONEOVEREIGHT/dx4 + &
                     ! diagonal neighbours in x-z plane
                     & (      (FI2(2*M1  ,2*M2-1,2*M3  ) + FI2(  M1l ,2*M2-1,  M3l ))                         - &
                     &        (FI2(2*M1  ,2*M2-1,  M3l ) + FI2(  M1l ,2*M2-1,2*M3  )) )**2 * ONEOVEREIGHT/dx4 + &
                     ! diagonal neighbours in y-z plane
                     & (      (FI2(2*M1-1,2*M2  ,2*M3  ) + FI2(2*M1-1,  M2l ,  M3l ))                         - &
                     &        (FI2(2*M1-1,2*M2  ,  M3l ) + FI2(2*M1-1,  M2l ,2*M3  )) )**2 * ONEOVEREIGHT/dx4
              !
              ! Source term of cell (2*M1  ,2*M2-1,2*M3-1)     
              Sigma_2 = Sigma_2                                                                               + &
                     ! direct neighbours in x direction
                     &        (FI2(  M1u ,2*M2-1,2*M3-1) + FI2(2*M1-1,2*M2-1,2*M3-1))                         * &
                     & (2.0D0*(FI2(  M1u ,2*M2-1,2*M3-1) + FI2(2*M1-1,2*M2-1,2*M3-1))                         - & 
                     &        (FI2(2*M1  ,2*M2  ,2*M3-1) + FI2(2*M1  ,  M2l ,2*M3-1))                         - &
                     &        (FI2(2*M1  ,2*M2-1,2*M3  ) + FI2(2*M1  ,2*M2-1,  M3l )) )    * ONEOVERTHREE/dx4 + &
                     ! direct neighbours in y direction
                     &        (FI2(2*M1  ,2*M2  ,2*M3-1) + FI2(2*M1  ,  M2l, 2*M3-1))                         * &
                     & (2.0D0*(FI2(2*M1  ,2*M2  ,2*M3-1) + FI2(2*M1  ,  M2l ,2*M3-1))                         - &
                     &        (FI2(  M1u ,2*M2-1,2*M3-1) + FI2(2*M1-1,2*M2-1,2*M3-1))                         - &
                     &        (FI2(2*M1  ,2*M2-1,2*M3  ) + FI2(2*M1  ,2*M2-1,  M3l )) )    * ONEOVERTHREE/dx4 + &
                     ! direct neighbours in z direction
                     &        (FI2(2*M1  ,2*M2-1,2*M3  ) + FI2(2*M1  ,2*M2-1,  M3l ))                         * &
                     & (2.0D0*(FI2(2*M1  ,2*M2-1,2*M3  ) + FI2(2*M1  ,2*M2-1,  M3l ))                         - &
                     &        (FI2(  M1u ,2*M2-1,2*M3-1) + FI2(2*M1-1,2*M2-1,2*M3-1))                         - &
                     &        (FI2(2*M1  ,2*M2  ,2*M3-1) + FI2(2*M1  ,  M2l ,2*M3-1)) )    * ONEOVERTHREE/dx4 + &
                     ! diagonal neighbours in x-y plane
                     & (      (FI2(  M1u ,2*M2  ,2*M3-1) + FI2(2*M1-1,  M2l ,2*M3-1))                         - &
                     &        (FI2(  M1u ,  M2l ,2*M3-1) + FI2(2*M1-1,2*M2  ,2*M3-1)) )**2 * ONEOVEREIGHT/dx4 + &
                     ! diagonal neighbours in x-z plane
                     & (      (FI2(  M1u ,2*M2-1,2*M3  ) + FI2(2*M1-1,2*M2-1,  M3l ))                         - &
                     &        (FI2(  M1u ,2*M2-1,  M3l ) + FI2(2*M1-1,2*M2-1,2*M3  )) )**2 * ONEOVEREIGHT/dx4 + &
                     ! diagonal neighbours in y-z plane
                     & (      (FI2(2*M1  ,2*M2  ,2*M3  ) + FI2(2*M1  ,  M2l ,  M3l ))                         - &
                     &        (FI2(2*M1  ,2*M2  ,  M3l ) + FI2(2*M1  ,  M2l ,2*M3  )) )**2 * ONEOVEREIGHT/dx4
              !
              ! Source term of cell (2*M1-1,2*M2  ,2*M3-1)     
              Sigma_3 = Sigma_3                                                                               + &
                     ! direct neighbours in x direction
                     &        (FI2(2*M1  ,2*M2  ,2*M3-1) + FI2(  M1l ,2*M2  ,2*M3-1))                         * &
                     & (2.0D0*(FI2(2*M1  ,2*M2  ,2*M3-1) + FI2(  M1l ,2*M2  ,2*M3-1))                         - & 
                     &        (FI2(2*M1-1,  M2u ,2*M3-1) + FI2(2*M1-1,2*M2-1,2*M3-1))                         - &
                     &        (FI2(2*M1-1,2*M2  ,2*M3  ) + FI2(2*M1-1,2*M2  ,  M3l )) )    * ONEOVERTHREE/dx4 + &
                     ! direct neighbours in y direction
                     &        (FI2(2*M1-1,  M2u ,2*M3-1) + FI2(2*M1-1,2*M2-1,2*M3-1))                         * &
                     & (2.0D0*(FI2(2*M1-1,  M2u ,2*M3-1) + FI2(2*M1-1,2*M2-1,2*M3-1))                         - &
                     &        (FI2(2*M1  ,2*M2  ,2*M3-1) + FI2(  M1l ,2*M2  ,2*M3-1))                         - &
                     &        (FI2(2*M1-1,2*M2  ,2*M3  ) + FI2(2*M1-1,2*M2  ,  M3l )) )    * ONEOVERTHREE/dx4 + &
                     ! direct neighbours in z direction
                     &        (FI2(2*M1-1,2*M2  ,2*M3  ) + FI2(2*M1-1,2*M2  ,  M3l ))                         * &
                     & (2.0D0*(FI2(2*M1-1,2*M2  ,2*M3  ) + FI2(2*M1-1,2*M2  ,  M3l ))                         - &
                     &        (FI2(2*M1  ,2*M2  ,2*M3-1) + FI2(  M1l ,2*M2  ,2*M3-1))                         - &
                     &        (FI2(2*M1-1,  M2u ,2*M3-1) + FI2(2*M1-1,2*M2-1,2*M3-1)) )    * ONEOVERTHREE/dx4 + &
                     ! diagonal neighbours in x-y plane
                     & (      (FI2(2*M1  ,  M2u ,2*M3-1) + FI2(  M1l ,2*M2-1,2*M3-1))                         - &
                     &        (FI2(2*M1  ,2*M2-1,2*M3-1) + FI2(  M1l ,  M2u ,2*M3-1)) )**2 * ONEOVEREIGHT/dx4 + &
                     ! diagonal neighbours in x-z plane
                     & (      (FI2(2*M1  ,2*M2  ,2*M3  ) + FI2(  M1l ,2*M2  ,  M3l ))                         - &
                     &        (FI2(2*M1  ,2*M2  ,  M3l ) + FI2(  M1l ,2*M2  ,2*M3  )) )**2 * ONEOVEREIGHT/dx4 + &
                     ! diagonal neighbours in y-z plane
                     & (      (FI2(2*M1-1,  M2u ,2*M3  ) + FI2(2*M1-1,2*M2-1,  M3l ))                         - &
                     &        (FI2(2*M1-1,  M2u ,  M3l ) + FI2(2*M1-1,2*M2-1,2*M3  )) )**2 * ONEOVEREIGHT/dx4
              !
              ! Source term of cell (2*M1  ,2*M2  ,2*M3-1)     
              Sigma_4 = Sigma_4                                                                               + &
                     ! direct neighbours in x direction
                     &        (FI2(  M1u ,2*M2  ,2*M3-1) + FI2(2*M1-1,2*M2  ,2*M3-1))                         * &
                     & (2.0D0*(FI2(  M1u ,2*M2  ,2*M3-1) + FI2(2*M1-1,2*M2  ,2*M3-1))                         - & 
                     &        (FI2(2*M1  ,  M2u ,2*M3-1) + FI2(2*M1  ,2*M2-1,2*M3-1))                         - &
                     &        (FI2(2*M1  ,2*M2  ,2*M3  ) + FI2(2*M1  ,2*M2  ,  M3l )) )    * ONEOVERTHREE/dx4 + &
                     ! direct neighbours in y direction
                     &        (FI2(2*M1  ,  M2u ,2*M3-1) + FI2(2*M1  ,2*M2-1,2*M3-1))                         * &
                     & (2.0D0*(FI2(2*M1  ,  M2u ,2*M3-1) + FI2(2*M1  ,2*M2-1,2*M3-1))                         - &
                     &        (FI2(  M1u ,2*M2  ,2*M3-1) + FI2(2*M1-1,2*M2  ,2*M3-1))                         - &
                     &        (FI2(2*M1  ,2*M2  ,2*M3  ) + FI2(2*M1  ,2*M2  ,  M3l )) )    * ONEOVERTHREE/dx4 + &
                     ! direct neighbours in z direction
                     &        (FI2(2*M1  ,2*M2  ,2*M3  ) + FI2(2*M1  ,2*M2  ,  M3l ))                         * &
                     & (2.0D0*(FI2(2*M1  ,2*M2  ,2*M3  ) + FI2(2*M1  ,2*M2  ,  M3l ))                         - &
                     &        (FI2(  M1u ,2*M2  ,2*M3-1) + FI2(2*M1-1,2*M2  ,2*M3-1))                         - &
                     &        (FI2(2*M1  ,  M2u ,2*M3-1) + FI2(2*M1  ,2*M2-1,2*M3-1)) )    * ONEOVERTHREE/dx4 + &
                     ! diagonal neighbours in x-y plane
                     & (      (FI2(  M1u ,  M2u ,2*M3-1) + FI2(2*M1-1,2*M2-1,2*M3-1))                         - &
                     &        (FI2(  M1u ,2*M2-1,2*M3-1) + FI2(2*M1-1,  M2u ,2*M3-1)) )**2 * ONEOVEREIGHT/dx4 + &
                     ! diagonal neighbours in x-z plane
                     & (      (FI2(  M1u ,2*M2  ,2*M3  ) + FI2(2*M1-1,2*M2  ,  M3l ))                         - &
                     &        (FI2(  M1u ,2*M2  ,  M3l ) + FI2(2*M1-1,2*M2  ,2*M3  )) )**2 * ONEOVEREIGHT/dx4 + &
                     ! diagonal neighbours in y-z plane
                     & (      (FI2(2*M1  ,  M2u ,2*M3  ) + FI2(2*M1  ,2*M2-1,  M3l ))                         - &
                     &        (FI2(2*M1  ,  M2u ,  M3l ) + FI2(2*M1  ,2*M2-1,2*M3  )) )**2 * ONEOVEREIGHT/dx4
              !
              ! Source term of cell (2*M1-1,2*M2-1,2*M3  )     
              Sigma_5 = Sigma_5                                                                               + &
                     ! direct neighbours in x direction
                     &        (FI2(2*M1  ,2*M2-1,2*M3  ) + FI2(  M1l ,2*M2-1,2*M3  ))                         * &
                     & (2.0D0*(FI2(2*M1  ,2*M2-1,2*M3  ) + FI2(  M1l ,2*M2-1,2*M3  ))                         - & 
                     &        (FI2(2*M1-1,2*M2  ,2*M3  ) + FI2(2*M1-1,  M2l ,2*M3  ))                         - &
                     &        (FI2(2*M1-1,2*M2-1,  M3u ) + FI2(2*M1-1,2*M2-1,2*M3-1)) )    * ONEOVERTHREE/dx4 + &
                     ! direct neighbours in y direction
                     &        (FI2(2*M1-1,2*M2  ,2*M3  ) + FI2(2*M1-1,  M2l ,2*M3  ))                         * &
                     & (2.0D0*(FI2(2*M1-1,2*M2  ,2*M3  ) + FI2(2*M1-1,  M2l ,2*M3  ))                         - &
                     &        (FI2(2*M1  ,2*M2-1,2*M3  ) + FI2(  M1l ,2*M2-1,2*M3  ))                         - &
                     &        (FI2(2*M1-1,2*M2-1,  M3u ) + FI2(2*M1-1,2*M2-1,2*M3-1)) )    * ONEOVERTHREE/dx4 + &
                     ! direct neighbours in z direction
                     &        (FI2(2*M1-1,2*M2-1,  M3u ) + FI2(2*M1-1,2*M2-1,2*M3-1))                         * &
                     & (2.0D0*(FI2(2*M1-1,2*M2-1,  M3u ) + FI2(2*M1-1,2*M2-1,2*M3-1))                         - &
                     &        (FI2(2*M1  ,2*M2-1,2*M3  ) + FI2(  M1l ,2*M2-1,2*M3  ))                         - &
                     &        (FI2(2*M1-1,2*M2  ,2*M3  ) + FI2(2*M1-1,  M2l ,2*M3  )) )    * ONEOVERTHREE/dx4 + &
                     ! diagonal neighbours in x-y plane
                     & (      (FI2(2*M1  ,2*M2  ,2*M3  ) + FI2(  M1l ,  M2l ,2*M3  ))                         - &
                     &        (FI2(2*M1  ,  M2l ,2*M3  ) + FI2(  M1l ,2*M2  ,2*M3  )) )**2 * ONEOVEREIGHT/dx4 + &
                     ! diagonal neighbours in x-z plane
                     & (      (FI2(2*M1  ,2*M2-1,  M3u ) + FI2(  M1l ,2*M2-1,2*M3-1))                         - &
                     &        (FI2(2*M1  ,2*M2-1,2*M3-1) + FI2(  M1l ,2*M2-1,  M3u )) )**2 * ONEOVEREIGHT/dx4 + &
                     ! diagonal neighbours in y-z plane
                     & (      (FI2(2*M1-1,2*M2  ,  M3u ) + FI2(2*M1-1,  M2l ,2*M3-1))                         - &
                     &        (FI2(2*M1-1,2*M2  ,2*M3-1) + FI2(2*M1-1,  M2l ,  M3u )) )**2 * ONEOVEREIGHT/dx4
              !
              ! Source term of cell (2*M1  ,2*M2-1,2*M3  )     
              Sigma_6 = Sigma_6                                                                               + &
                     ! direct neighbours in x direction
                     &        (FI2(  M1u ,2*M2-1,2*M3  ) + FI2(2*M1-1,2*M2-1,2*M3  ))                         * &
                     & (2.0D0*(FI2(  M1u ,2*M2-1,2*M3  ) + FI2(2*M1-1,2*M2-1,2*M3  ))                         - & 
                     &        (FI2(2*M1  ,2*M2  ,2*M3  ) + FI2(2*M1  ,  M2l ,2*M3  ))                         - &
                     &        (FI2(2*M1  ,2*M2-1,  M3u ) + FI2(2*M1  ,2*M2-1,2*M3-1)) )    * ONEOVERTHREE/dx4 + &
                     ! direct neighbours in y direction
                     &        (FI2(2*M1  ,2*M2  ,2*M3  ) + FI2(2*M1  ,  M2l ,2*M3  ))                         * &
                     & (2.0D0*(FI2(2*M1  ,2*M2  ,2*M3  ) + FI2(2*M1  ,  M2l ,2*M3  ))                         - &
                     &        (FI2(  M1u ,2*M2-1,2*M3  ) + FI2(2*M1-1,2*M2-1,2*M3  ))                         - &
                     &        (FI2(2*M1  ,2*M2-1,  M3u ) + FI2(2*M1  ,2*M2-1,2*M3-1)) )    * ONEOVERTHREE/dx4 + &
                     ! direct neighbours in z direction
                     &        (FI2(2*M1  ,2*M2-1,  M3u ) + FI2(2*M1  ,2*M2-1,2*M3-1))                         * &
                     & (2.0D0*(FI2(2*M1  ,2*M2-1,  M3u ) + FI2(2*M1  ,2*M2-1,2*M3-1))                         - &
                     &        (FI2(  M1u ,2*M2-1,2*M3  ) + FI2(2*M1-1,2*M2-1,2*M3  ))                         - &
                     &        (FI2(2*M1  ,2*M2  ,2*M3  ) + FI2(2*M1  ,  M2l ,2*M3  )) )    * ONEOVERTHREE/dx4 + &
                     ! diagonal neighbours in x-y plane
                     & (      (FI2(  M1u ,2*M2  ,2*M3  ) + FI2(2*M1-1,  M2l ,2*M3  ))                         - &
                     &        (FI2(  M1u ,  M2l ,2*M3  ) + FI2(2*M1-1,2*M2  ,2*M3  )) )**2 * ONEOVEREIGHT/dx4 + &
                     ! diagonal neighbours in x-z plane
                     & (      (FI2(  M1u ,2*M2-1,  M3u ) + FI2(2*M1-1,2*M2-1,2*M3-1))                         - &
                     &        (FI2(  M1u ,2*M2-1,2*M3-1) + FI2(2*M1-1,2*M2-1,  M3u )) )**2 * ONEOVEREIGHT/dx4 + &
                     ! diagonal neighbours in y-z plane
                     & (      (FI2(2*M1  ,2*M2  ,  M3u ) + FI2(2*M1  ,  M2l ,2*M3-1))                         - &
                     &        (FI2(2*M1  ,2*M2  ,2*M3-1) + FI2(2*M1  ,  M2l ,  M3u )) )**2 * ONEOVEREIGHT/dx4
              !
              ! Source term of cell (2*M1-1,2*M2  ,2*M3  )     
              Sigma_7 = Sigma_7                                                                               + &
                     ! direct neighbours in x direction
                     &        (FI2(2*M1  ,2*M2  ,2*M3  ) + FI2(  M1l ,2*M2  ,2*M3  ))                         * &
                     & (2.0D0*(FI2(2*M1  ,2*M2  ,2*M3  ) + FI2(  M1l ,2*M2  ,2*M3  ))                         - & 
                     &        (FI2(2*M1-1,  M2u ,2*M3  ) + FI2(2*M1-1,2*M2-1,2*M3  ))                         - &
                     &        (FI2(2*M1-1,2*M2  ,  M3u ) + FI2(2*M1-1,2*M2  ,2*M3-1)) )    * ONEOVERTHREE/dx4 + &
                     ! direct neighbours in y direction
                     &        (FI2(2*M1-1,  M2u ,2*M3  ) + FI2(2*M1-1,2*M2-1,2*M3  ))                         * &
                     & (2.0D0*(FI2(2*M1-1,  M2u ,2*M3  ) + FI2(2*M1-1,2*M2-1,2*M3  ))                         - &
                     &        (FI2(2*M1  ,2*M2  ,2*M3  ) + FI2(  M1l ,2*M2  ,2*M3  ))                         - &
                     &        (FI2(2*M1-1,2*M2  ,  M3u ) + FI2(2*M1-1,2*M2  ,2*M3-1)) )    * ONEOVERTHREE/dx4 + &
                     ! direct neighbours in z direction
                     &        (FI2(2*M1-1,2*M2  ,  M3u ) + FI2(2*M1-1,2*M2  ,2*M3-1))                         * &
                     & (2.0D0*(FI2(2*M1-1,2*M2  ,  M3u ) + FI2(2*M1-1,2*M2  ,2*M3-1))                         - &
                     &        (FI2(2*M1  ,2*M2  ,2*M3  ) + FI2(  M1l ,2*M2  ,2*M3  ))                         - &
                     &        (FI2(2*M1-1,  M2u ,2*M3  ) + FI2(2*M1-1,2*M2-1,2*M3  )) )    * ONEOVERTHREE/dx4 + &
                     ! diagonal neighbours in x-y plane
                     & (      (FI2(2*M1  ,  M2u ,2*M3  ) + FI2(  M1l ,2*M2-1,2*M3  ))                         - &
                     &        (FI2(2*M1  ,2*M2-1,2*M3  ) + FI2(  M1l ,  M2u ,2*M3  )) )**2 * ONEOVEREIGHT/dx4 + &
                     ! diagonal neighbours in x-z plane
                     & (      (FI2(2*M1  ,2*M2  ,  M3u ) + FI2(  M1l ,2*M2  ,2*M3-1))                         - &
                     &        (FI2(2*M1  ,2*M2  ,2*M3-1) + FI2(  M1l ,2*M2  ,  M3u )) )**2 * ONEOVEREIGHT/dx4 + &
                     ! diagonal neighbours in y-z plane
                     & (      (FI2(2*M1-1,  M2u ,  M3u ) + FI2(2*M1-1,2*M2-1,2*M3-1))                         - &
                     &        (FI2(2*M1-1,  M2u ,2*M3-1) + FI2(2*M1-1,2*M2-1,  M3u )) )**2 * ONEOVEREIGHT/dx4
              !
              ! Source term of cell (2*M1  ,2*M2  ,2*M3  )     
              Sigma_8 = Sigma_8                                                                               + &
                     ! direct neighbours in x direction
                     &        (FI2(  M1u ,2*M2  ,2*M3  ) + FI2(2*M1-1,2*M2  ,2*M3  ))                         * &
                     & (2.0D0*(FI2(  M1u ,2*M2  ,2*M3  ) + FI2(2*M1-1,2*M2  ,2*M3  ))                         - & 
                     &        (FI2(2*M1  ,  M2u ,2*M3  ) + FI2(2*M1  ,2*M2-1,2*M3  ))                         - &
                     &        (FI2(2*M1  ,2*M2  ,  M3u ) + FI2(2*M1  ,2*M2  ,2*M3-1)) )    * ONEOVERTHREE/dx4 + &
                     ! direct neighbours in y direction
                     &        (FI2(2*M1  ,  M2u ,2*M3  ) + FI2(2*M1  ,2*M2-1,2*M3  ))                         * &
                     & (2.0D0*(FI2(2*M1  ,  M2u ,2*M3  ) + FI2(2*M1  ,2*M2-1,2*M3  ))                         - &
                     &        (FI2(  M1u ,2*M2  ,2*M3  ) + FI2(2*M1-1,2*M2  ,2*M3  ))                         - &
                     &        (FI2(2*M1  ,2*M2  ,  M3u ) + FI2(2*M1  ,2*M2  ,2*M3-1)) )    * ONEOVERTHREE/dx4 + &
                     ! direct neighbours in z direction
                     &        (FI2(2*M1  ,2*M2  ,  M3u ) + FI2(2*M1  ,2*M2  ,2*M3-1))                         * &
                     & (2.0D0*(FI2(2*M1  ,2*M2  ,  M3u ) + FI2(2*M1  ,2*M2  ,2*M3-1))                         - &
                     &        (FI2(  M1u ,2*M2  ,2*M3  ) + FI2(2*M1-1,2*M2  ,2*M3  ))                         - &
                     &        (FI2(2*M1  ,  M2u ,2*M3  ) + FI2(2*M1  ,2*M2-1,2*M3  )) )    * ONEOVERTHREE/dx4 + &
                     ! diagonal neighbours in x-y plane
                     & (      (FI2(  M1u ,  M2u ,2*M3  ) + FI2(2*M1-1,2*M2-1,2*M3  ))                         - &
                     &        (FI2(  M1u ,2*M2-1,2*M3  ) + FI2(2*M1-1,  M2u ,2*M3  )) )**2 * ONEOVEREIGHT/dx4 + &
                     ! diagonal neighbours in x-z plane
                     & (      (FI2(  M1u ,2*M2  ,  M3u ) + FI2(2*M1-1,2*M2  ,2*M3-1))                         - &
                     &        (FI2(  M1u ,2*M2  ,2*M3-1) + FI2(2*M1-1,2*M2  ,  M3u )) )**2 * ONEOVEREIGHT/dx4 + &
                     ! diagonal neighbours in y-z plane
                     & (      (FI2(2*M1  ,  M2u ,  M3u ) + FI2(2*M1  ,2*M2-1,2*M3-1))                         - &
                     &        (FI2(2*M1  ,  M2u ,2*M3-1) + FI2(2*M1  ,2*M2-1,  M3u )) )**2 * ONEOVEREIGHT/dx4
              !
              Q = Q - fct3*DSQRT(alpha**2+EIGHTOVERTHREE*Sigma_1)
              Q = Q - fct3*DSQRT(alpha**2+EIGHTOVERTHREE*Sigma_2)
              Q = Q - fct3*DSQRT(alpha**2+EIGHTOVERTHREE*Sigma_3)
              Q = Q - fct3*DSQRT(alpha**2+EIGHTOVERTHREE*Sigma_4)
              Q = Q - fct3*DSQRT(alpha**2+EIGHTOVERTHREE*Sigma_5)
              Q = Q - fct3*DSQRT(alpha**2+EIGHTOVERTHREE*Sigma_6)
              Q = Q - fct3*DSQRT(alpha**2+EIGHTOVERTHREE*Sigma_7)
              Q = Q - fct3*DSQRT(alpha**2+EIGHTOVERTHREE*Sigma_8)
              !
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
              !
              ! finalise residuals from all 8 son cells
              P = P/dx2 + Q 
!              P = P + 8.0D0*S_mean ! Baojiu-13-06-2021 added S_mean contribution
              FI3(M1,M2,M3) = P/8.0D0
           END DO
        END DO
     END DO
  ELSE
     koffset  = 2**(levelmax-ilevel  )*(2**(ilevel  )-2)
     ioffset  = 2**(levelmax-ilevel+1)                        
     koffset2 = 2**(levelmax-ilevel+1)*(2**(ilevel-1)-1) 
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3) &
!$OMP PRIVATE (P)
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
           END DO
        END DO
     END DO
  END IF
    
  CALL TimingMain(3,1)
END SUBROUTINE restrict_residual_DGP

!-------------------------------------------------------------
SUBROUTINE calculate_physical_right_hand_side_DGP(ilevel)
!-------------------------------------------------------------
use Tools

  integer :: ilevel
  real*8  :: alpha,beta,Orc,Rc_sq,dx,dx2,dx4
  integer :: ngrid_level
  integer :: ioffset,joffset,koffest,koffset2
  real*8  :: fct1
  real*8  :: TWOOVERTHREE,ONEOVEREIGHT,THREEOVERFOUR,EIGHTOVERTHREE

  integer :: M1,M2,M3,M1l,M2l,M3l,M1u,M2u,M3u
  real*8  :: OP,Sigma

  ! number of grid points on the coarse level
  ngrid_level = NGRID/2**ilevel

  ONEOVEREIGHT   = 1.0D0/8.0D0
  TWOOVERTHREE   = 2.0D0/3.0D0
  THREEOVERFOUR  = 3.0D0/4.0D0
  EIGHTOVERTHREE = 8.0D0/3.0D0

  dx      = DBLE(2**ilevel) 
  dx2     = dx*dx
  dx4     = dx2*dx2
  Orc     = 1.0D0/(4.0D0*H0rc**2)
  !
  IF(N_branch) beta = 1.0D0 + (0.5D0*Om/AEXPN**3+OmL) / (DSQRT(Orc*(Om/AEXPN**3+OmL))) ! normal branch
  IF(S_branch) beta =       - (0.5D0*Om/AEXPN**3+Orc) / (DSQRT(Orc*(Om/AEXPN**3+Orc))) ! self-accelerated branch
  !
  Rc_sq  = 1.0D0/(4.0D0*Orc)
  alpha  = 3.0D0*beta*AEXPN**2/Rc_sq
  fct1   = alpha*Om/beta/AEXPN

  ioffset  = 2**(levelmax-ilevel)
  joffset  = 2**(levelmax-ilevel)
  koffset  = 2**(levelmax-ilevel)*(2**ilevel-2)
  koffset2 = 2**(levelmax-ilevel)*(2**ilevel-1)
 
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3,M1u,M1l,M2u,M2l,M3u,M3l) &
!$OMP PRIVATE (OP,Sigma)
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
           ! calculate Sigma
           Sigma = fct1*FI3(M1+ioffset,M2,M3+koffset)
           Sigma = Sigma                                                                                                                               + &
                & ((FI3(M1u,M2 +joffset,M3 +koffset) + FI3(M1l,M2 +joffset,M3 +koffset) - 2.0D0*FI3(M1,M2+joffset,M3+koffset))**2                      + &
                &  (FI3(M1 ,M2u+joffset,M3 +koffset) + FI3(M1 ,M2l+joffset,M3 +koffset) - 2.0D0*FI3(M1,M2+joffset,M3+koffset))**2                      + &
                &  (FI3(M1 ,M2 +joffset,M3u+koffset) + FI3(M1 ,M2 +joffset,M3l+koffset) - 2.0D0*FI3(M1,M2+joffset,M3+koffset))**2 ) * TWOOVERTHREE/dx4 - &
                & ((FI3(M1u,M2 +joffset,M3 +koffset) + FI3(M1l,M2 +joffset,M3 +koffset) - 2.0D0*FI3(M1,M2+joffset,M3+koffset))                         * &
                &  (FI3(M1 ,M2u+joffset,M3 +koffset) + FI3(M1 ,M2l+joffset,M3 +koffset) - 2.0D0*FI3(M1,M2+joffset,M3+koffset))                         + &
                &  (FI3(M1u,M2 +joffset,M3 +koffset) + FI3(M1l,M2 +joffset,M3 +koffset) - 2.0D0*FI3(M1,M2+joffset,M3+koffset))                         * &
                &  (FI3(M1 ,M2 +joffset,M3u+koffset) + FI3(M1 ,M2 +joffset,M3l+koffset) - 2.0D0*FI3(M1,M2+joffset,M3+koffset))                         + &
                &  (FI3(M1 ,M2u+joffset,M3 +koffset) + FI3(M1 ,M2l+joffset,M3 +koffset) - 2.0D0*FI3(M1,M2+joffset,M3+koffset))                         * &
                &  (FI3(M1 ,M2 +joffset,M3u+koffset) + FI3(M1 ,M2 +joffset,M3l+koffset) - 2.0D0*FI3(M1,M2+joffset,M3+koffset)) )    * TWOOVERTHREE/dx4 + &
                & ((FI3(M1u,M2u+joffset,M3 +koffset) + FI3(M1l,M2l+joffset,M3 +koffset)                                                                - &
                &   FI3(M1u,M2l+joffset,M3 +koffset) - FI3(M1l,M2u+joffset,M3 +koffset))**2                                                            + &
                &  (FI3(M1u,M2 +joffset,M3u+koffset) + FI3(M1l,M2 +joffset,M3l+koffset)                                                                - &
                &   FI3(M1u,M2 +joffset,M3l+koffset) - FI3(M1l,M2 +joffset,M3u+koffset))**2                                                            + &
                &  (FI3(M1 ,M2u+joffset,M3u+koffset) + FI3(M1 ,M2l+joffset,M3l+koffset)                                                                - &
                &   FI3(M1 ,M2u+joffset,M3l+koffset) - FI3(M1 ,M2l+joffset,M3u+koffset))**2) * ONEOVEREIGHT/dx4
           ! calculate Laplancian of PDE
           OP    = (FI3(M1u,joffset+M2 ,koffset+M3 ) + &
                 &  FI3(M1l,joffset+M2 ,koffset+M3 ) + &
                 &  FI3(M1 ,joffset+M2u,koffset+M3 ) + &
                 &  FI3(M1 ,joffset+M2l,koffset+M3 ) + &
                 &  FI3(M1 ,joffset+M2 ,koffset+M3u) + &
                 &  FI3(M1 ,joffset+M2 ,koffset+M3l) - &
                 &  FI3(M1 ,joffset+M2 ,koffset+M3 ) * 6.0D0 ) / dx2
           OP   = OP - (DSIGN(1.0D0,alpha)*DSQRT(alpha**2+EIGHTOVERTHREE*Sigma)-alpha) * THREEOVERFOUR
           OP   = OP - FI3(M1,M2,M3+koffset)
!           OP   = OP + S_mean ! Baojiu-13-06-2021 added S_mean contribution
           !
           FI3(M1,joffset+M2,koffset2+M3) =  OP
        END DO
     END DO
  END DO
    
  CALL TimingMain(3,1)

END SUBROUTINE calculate_physical_right_hand_side_DGP
