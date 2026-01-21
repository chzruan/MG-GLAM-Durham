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
SUBROUTINE relaxation_iterations_sym(ilevel,redstep)
!-------------------------------------------------------------
  use Tools

  implicit none

  integer :: ilevel
  logical :: redstep
  real*8  :: ctilde,dx,dx2
  real*8  :: fct1,fct2,fct3,fct4 
  real*8  :: ZERO,ONEOVERTHREE,TWOOVERTHREE,TWOPIOVERTHREE
  integer :: nstart,nfinal,ngrid_level
  integer :: ioffset,joffset,koffset,koffset2

  real*8  :: Delta0,Delta1,sqrt_Delta0,tmp1,tmp2,P,Q,S,CC
  integer :: M1,M2,M3,M1u,M1l,M2u,M2l,M3u,M3l

  integer :: nid,myid

  IF(MG_test) WRITE(*,'(A,I5,F7.4)') 'Relaxation iterations on level',levelmax-ilevel,AEXPN

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
  fct1    = (sym_astar/AEXPN)**3
  fct4    = 2.0D0*(ctilde*sym_xi/AEXPN)**2
  fct2    = -fct4/dx2
  fct3    =  fct4/dx2*6.0D0

  ! offsets for cell access to array FI3
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
!$OMP PRIVATE (Delta1,Delta0,sqrt_Delta0,tmp1,tmp2) &
!$OMP PRIVATE (P,Q,S,CC)
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
              P = fct1*(FI (M1 ,M2 ,M3 )+1.0D0)-1.0D0+fct3
              Q = fct2*(FI2(M1u,M2 ,M3 )    + &
                     &  FI2(M1l,M2 ,M3 )    + &
                     &  FI2(M1 ,M2u,M3 )    + &
                     &  FI2(M1 ,M2l,M3 )    + &
                     &  FI2(M1 ,M2 ,M3u)    + &
                     &  FI2(M1 ,M2 ,M3l))   
           ELSE 
              P = fct1*(FI3(M1 +ioffset,M2         ,M3 +koffset )+1.0D0)-1.0D0+fct3
              Q = fct4* FI3(M1         ,M2 +joffset,M3 +koffset2)    + &   ! physical right-hand side  
                & fct2*(FI3(M1u+ioffset,M2 +joffset,M3 +koffset )    + &
                     &  FI3(M1l+ioffset,M2 +joffset,M3 +koffset )    + &
                     &  FI3(M1 +ioffset,M2u+joffset,M3 +koffset )    + &
                     &  FI3(M1 +ioffset,M2l+joffset,M3 +koffset )    + &
                     &  FI3(M1 +ioffset,M2 +joffset,M3u+koffset )    + &
                     &  FI3(M1 +ioffset,M2 +joffset,M3l+koffset ))
           ENDIF
           !
           Delta0 = -3.0D0*P
           Delta1 = 27.0D0*Q
           !
           CC = Delta1**2-4.0D0*Delta0**3
           !
           IF(ilevel.EQ.0) THEN
           ! solution on PM grid
              IF(CC.GE.0.0D0) THEN
                 tmp1 = 0.5D0*(DSQRT(CC)+Delta1)
                 tmp2 = tmp1-Delta1
                 IF(DABS(P).LT.ZERO) THEN
                    FI2(M1,M2,M3) = - DSIGN(DABS(Q   )**ONEOVERTHREE,Q   )
                 ELSE
                    FI2(M1,M2,M3) = -(DSIGN(DABS(tmp1)**ONEOVERTHREE,tmp1)-DSIGN(DABS(tmp2)**ONEOVERTHREE,tmp2))*ONEOVERTHREE
                 ENDIF
              ELSE
                 sqrt_Delta0 = DSQRT(Delta0)
                 S = DACOS(0.5D0*Delta1/(sqrt_Delta0**3))
                 FI2(M1,M2,M3) = -TWOOVERTHREE*sqrt_Delta0*DCOS(S*ONEOVERTHREE+TWOPIOVERTHREE)
              ENDIF
           ELSE
           ! solution on multigrid
              IF(CC.GE.0.0D0) THEN
                 tmp1 = 0.5D0*(DSQRT(CC)+Delta1)
                 tmp2 = tmp1-Delta1
                 IF(DABS(P).LT.ZERO) THEN
                    FI3(M1+ioffset,M2+joffset,M3+koffset) = - DSIGN(DABS(Q   )**ONEOVERTHREE,Q   )
                 ELSE
                    FI3(M1+ioffset,M2+joffset,M3+koffset) = -(DSIGN(DABS(tmp1)**ONEOVERTHREE,tmp1)-DSIGN(DABS(tmp2)**ONEOVERTHREE,tmp2))*ONEOVERTHREE
                 ENDIF
              ELSE
                 sqrt_Delta0 = DSQRT(Delta0)
                 S = DACOS(0.5D0*Delta1/sqrt_Delta0**3)
                 FI3(M1+ioffset,M2+joffset,M3+koffset) = -TWOOVERTHREE*sqrt_Delta0*DCOS(S*ONEOVERTHREE+TWOPIOVERTHREE)
              ENDIF
           ENDIF
        END DO
     END DO
  END DO

  CALL TimingMain(3,1)
    
END SUBROUTINE relaxation_iterations_sym

!-------------------------------------------------------------
!
! Calculate resudual and its RMS on the entire grid on ilevel.
!
!-------------------------------------------------------------
SUBROUTINE calculate_residual_sym(ilevel,res_PM_grid)
!-------------------------------------------------------------
use Tools

  integer,intent(in) :: ilevel
  real*8,intent(out) :: res_PM_grid
  
  real*8  :: ctilde,ctilde2,dx,dx2
  integer :: ngrid_level
  integer :: ioffset,joffset,koffset,koffset2
  real*8  :: fct1,fct2

  integer :: M1,M2,M3,M1l,M1u,M2l,M2u,M3l,M3u
  integer :: N1,N2,N3
  real*8  :: RES,OP,RES2

  IF(MG_test) WRITE(*,'(A,I5)') 'Calculate residual on level',levelmax-ilevel

  ! number of grid points on the coarse level
  ngrid_level = NGRID/2**ilevel
  
  ctilde  = 2.99792458D3*DBLE(NGRID)/Box
  ctilde2 = ctilde**2
  dx      = DBLE(2**ilevel)
  dx2     = dx*dx

  fct1    = (sym_astar/AEXPN)**3
  fct2    = 0.5D0*(AEXPN/(ctilde*sym_xi))**2

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
           IF(ilevel.EQ.0) THEN
              ! calculate Laplacian of PDE
              OP = FI2(M1u,M2 ,M3 )      + &
                 & FI2(M1l,M2 ,M3 )      + &
                 & FI2(M1 ,M2u,M3 )      + &
                 & FI2(M1 ,M2l,M3 )      + &
                 & FI2(M1 ,M2 ,M3u)      + &
                 & FI2(M1 ,M2 ,M3l)      - &
                 & FI2(M1 ,M2 ,M3 )*6.0D0
              OP = OP/dx2
              OP = OP-fct2*(fct1*(FI(M1,M2,M3)+1.0D0)-1.0D0)*FI2(M1,M2,M3)
              OP = OP-fct2*                                  FI2(M1,M2,M3)**3
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
              OP = OP/dx2
              OP = OP-fct2*(fct1*(FI3(M1+ioffset,M2,M3+koffset)+1.0D0)-1.0D0)*FI3(M1+ioffset,M2+joffset,M3+koffset)
              OP = OP-fct2*                                                   FI3(M1+ioffset,M2+joffset,M3+koffset)**3
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
      
END SUBROUTINE calculate_residual_sym


!-------------------------------------------------------------
!
! Calculate restricted residual field on coarse level "ilevel"
! Note: "ilevel" means the grid ilevel coarer than the PM grid
!
!-------------------------------------------------------------
SUBROUTINE restrict_residual_sym(ilevel)
!-------------------------------------------------------------
use Tools

  integer :: ilevel
  real*8  :: ctilde,ctilde2,dx,dx2
  real*8  :: fct1,fct2
  integer :: ngrid_level 
  integer :: ioffset,joffset,koffset,ioffset2,joffset2,koffset2
 
  integer :: M1,M2,M3
  integer :: M1u,M1l,M2u,M2l,M3u,M3l 
  real*8  :: P,Q

  IF(MG_test) WRITE(*,'(A,I5)') 'Restrict residual to level',levelmax-ilevel

  ! number of grid points on the coarse level
  ngrid_level = NGRID/2**ilevel

  ctilde  = 2.99792458D3*DBLE(NGRID)/Box
  ctilde2 = ctilde**2
  dx      = 1.0d0  
  dx2     = dx*dx

  fct1    = (sym_astar/AEXPN)**3
  fct2    = 0.5D0*(AEXPN/(ctilde*sym_xi))**2

  IF(ilevel.EQ.1) THEN
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3,M1u,M1l,M2u,M2l,M3u,M3l) &
!$OMP PRIVATE (P,Q) 
     DO M3=1,ngrid_level
        DO M2=1,ngrid_level
           DO M1=1,ngrid_level
              P = 0.0d0                                         ! Laplancian terms
              Q = 0.0d0                                         ! non-Laplacian terms
              ! accumulate contributions to the non-Laplacian term 
              ! of PDE from all 8 son cells
              Q = Q+fct2*(fct1*(FI (2*M1-1,2*M2-1,2*M3-1)+1.0D0)-1.0D0)* &
                &               FI2(2*M1-1,2*M2-1,2*M3-1)         &
                &  +fct2*       FI2(2*M1-1,2*M2-1,2*M3-1)**3
              Q = Q+fct2*(fct1*(FI (2*M1  ,2*M2-1,2*M3-1)+1.0D0)-1.0D0)* &
                &               FI2(2*M1  ,2*M2-1,2*M3-1)         &
                &  +fct2*       FI2(2*M1  ,2*M2-1,2*M3-1)**3
              Q = Q+fct2*(fct1*(FI (2*M1-1,2*M2  ,2*M3-1)+1.0D0)-1.0D0)* &
                &               FI2(2*M1-1,2*M2  ,2*M3-1)         &
                &  +fct2*       FI2(2*M1-1,2*M2  ,2*M3-1)**3
              Q = Q+fct2*(fct1*(FI (2*M1  ,2*M2  ,2*M3-1)+1.0D0)-1.0D0)* &
                &               FI2(2*M1  ,2*M2  ,2*M3-1)         &
                &  +fct2*       FI2(2*M1  ,2*M2  ,2*M3-1)**3
              Q = Q+fct2*(fct1*(FI (2*M1-1,2*M2-1,2*M3  )+1.0D0)-1.0D0)* &
                &               FI2(2*M1-1,2*M2-1,2*M3  )         &
                &  +fct2*       FI2(2*M1-1,2*M2-1,2*M3  )**3
              Q = Q+fct2*(fct1*(FI (2*M1  ,2*M2-1,2*M3  )+1.0D0)-1.0D0)* &
                &               FI2(2*M1  ,2*M2-1,2*M3  )         &
                &  +fct2*       FI2(2*M1  ,2*M2-1,2*M3  )**3
              Q = Q+fct2*(fct1*(FI (2*M1-1,2*M2  ,2*M3  )+1.0D0)-1.0D0)* &
                &               FI2(2*M1-1,2*M2  ,2*M3  )         &
                &  +fct2*       FI2(2*M1-1,2*M2  ,2*M3  )**3
              Q = Q+fct2*(fct1*(FI (2*M1  ,2*M2  ,2*M3  )+1.0D0)-1.0D0)* &
                &               FI2(2*M1  ,2*M2  ,2*M3  )         &
                &  +fct2*       FI2(2*M1  ,2*M2  ,2*M3  )**3
              ! prepare indices of neighbour cells on the 3-point
              ! stencil and apply the periodic boundary condition
              M1u=2*M1+1; IF(M1u>NGRID) M1u=1
              M1l=2*M1-2; IF(M1l<1    ) M1l=NGRID
              M2u=2*M2+1; IF(M2u>NGRID) M2u=1
              M2l=2*M2-2; IF(M2l<1    ) M2l=NGRID
              M3u=2*M3+1; IF(M3u>NGRID) M3u=1
              M3l=2*M3-2; IF(M3l<1    ) M3l=NGRID
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
              P = P/dx2-Q
              !
              FI3(M1,M2,M3) = P/8.0D0  
           END DO
        END DO
     END DO
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
END SUBROUTINE restrict_residual_sym


!-------------------------------------------------------------
!
! Calculate the physical RHS of discrete PDE on coarser ilevel
!
!-------------------------------------------------------------
SUBROUTINE calculate_physical_right_hand_side_sym(ilevel)
!-------------------------------------------------------------
use Tools

  integer :: ilevel
  real*8  :: ctilde,dx,dx2
  real*8  :: fct1,fct2
  integer :: ngrid_level
  integer :: ioffset,joffset,koffest,koffset2

  integer :: M1,M2,M3,M1l,M2l,M3l,M1u,M2u,M3u
  real*8  :: OP

  IF(MG_test) WRITE(*,'(A,I5)') 'Calculate physical right-hand side on level',levelmax-ilevel

  ! number of grid points on the coarse level
  ngrid_level = NGRID/2**ilevel

  ctilde  = 2.99792458D3*DBLE(NGRID)/Box
  dx      = DBLE(2**ilevel) 
  dx2     = dx*dx

  fct1   = (sym_astar/AEXPN)**3
  fct2   = 0.5D0*(AEXPN/(ctilde*sym_xi))**2

  ioffset  = 2**(levelmax-ilevel)
  joffset  = 2**(levelmax-ilevel)
  koffset  = 2**(levelmax-ilevel)*(2**ilevel-2)
  koffset2 = 2**(levelmax-ilevel)*(2**ilevel-1)

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
           ! calculate Laplancian of PDE using RESTRICTED scalar field
           OP = FI3(M1u,joffset+M2 ,koffset+M3 )      + &
              & FI3(M1l,joffset+M2 ,koffset+M3 )      + &
              & FI3(M1 ,joffset+M2u,koffset+M3 )      + &
              & FI3(M1 ,joffset+M2l,koffset+M3 )      + &
              & FI3(M1 ,joffset+M2 ,koffset+M3u)      + &
              & FI3(M1 ,joffset+M2 ,koffset+M3l)      - &
              & FI3(M1 ,joffset+M2 ,koffset+M3 )*6.0D0 
           !
           OP = OP/dx2
           OP = OP-fct2*(fct1*(FI3(M1+ioffset,M2,M3+koffset)+1.0D0)-1.0D0)*FI3(M1,joffset+M2,koffset+M3)
           OP = OP-fct2*                                                   FI3(M1,joffset+M2,koffset+M3)**3
           OP = OP-FI3(M1,M2,M3+koffset)                        ! add restricted residual (note the minus sign!!)
           ! 
           FI3(M1,joffset+M2,koffset2+M3) = OP  
           ! 
        END DO
     END DO
  END DO

  CALL TimingMain(3,1)

END SUBROUTINE calculate_physical_right_hand_side_sym
