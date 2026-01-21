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
SUBROUTINE relaxation_iterations_csf(ilevel,redstep)
!-------------------------------------------------------------
  use Tools
  use ExtradofBackgroundData

  implicit none

  integer :: ilevel
  logical :: redstep
  real*8  :: ctilde,ctilde2,dx,dx2,sixodx2
  real*8  :: fct1
  integer :: nstart,nfinal,ngrid_level
  integer :: ioffset,joffset,koffset,koffset2

  integer :: M1,M2,M3,M1u,M1l,M2u,M2l,M3u,M3l

  integer :: nid,myid

  integer :: ist
  real*8  :: phibar,A_phibar,V_phibar
  real*8  :: phi,rho,A_phi,V_phi,A_pp,V_pp,L,dL

  IF(MG_test) WRITE(*,'(A,I5,F7.4)') 'Relaxation iterations on level',levelmax-ilevel,AEXPN

  ! number of grid points on the coarse level
  ngrid_level = NGRID/(2**ilevel)
  !
  ! simulation parameters
  ctilde  = 2.99792458D3*DBLE(NGRID)/Box
  ctilde2 = ctilde**2
  dx      = DBLE(2**ilevel) 
  dx2     = dx*dx
  sixodx2 = 6.0D0/dx2
  !
  ! physical and numerical quantities
  fct1   = 3.0D0*Om/AEXPN
  !
  ist = 1
  DO WHILE(BackgroundEvolution(ist,1).LT.AEXPN)
     ist=ist+1
  ENDDO
  !
  phibar =  BackgroundEvolution(ist-1,2)+ &
         & (BackgroundEvolution(ist,  2)-BackgroundEvolution(ist-1,2))/ &
         & (BackgroundEvolution(ist,  1)-BackgroundEvolution(ist-1,1))* &
         & (AEXPN                       -BackgroundEvolution(ist-1,1))
  !
  ! background coupling function 
  IF(csf_coupling.EQ.1) THEN
     A_phibar = csf_beta*DEXP(csf_beta*phibar)
  ENDIF
  !
  ! backgrund potential derivative: (dV/d\varphi)*a^2/M_Pl^2
  IF(csf_potential.EQ.1) THEN
     V_phibar = -csf_alpha*csf_lambda**2*AEXPN**2/(phibar**(1.0D0+csf_alpha))
  ENDIF
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
!$OMP PRIVATE (M1,M2,M3,M1u,M1l,M2u,M2l,M3u,M3l) &
!$OMP PRIVATE (phi,A_phi,V_phi,A_pp,V_pp,rho,L,dL)
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
              !
              L = FI2(M1u,M2 ,M3 )    + &
                & FI2(M1l,M2 ,M3 )    + &
                & FI2(M1 ,M2u,M3 )    + &
                & FI2(M1 ,M2l,M3 )    + &
                & FI2(M1 ,M2 ,M3u)    + &
                & FI2(M1 ,M2 ,M3l)    - &
                & FI2(M1 ,M2 ,M3 )*6.0D0
              !
              phi = phibar+FI2(M1,M2,M3)/ctilde2                        ! scalar field (background+perturbation)
              rho = FI(M1,M2,M3)+1.0D0                                  ! density      (background+perturbation)
              !
           ELSE 
              !
              L = -FI3(M1         ,M2 +joffset,M3 +koffset2)*dx2 + &    ! physical right-hand side  
                &  FI3(M1u+ioffset,M2 +joffset,M3 +koffset )     + &
                &  FI3(M1l+ioffset,M2 +joffset,M3 +koffset )     + &
                &  FI3(M1 +ioffset,M2u+joffset,M3 +koffset )     + &
                &  FI3(M1 +ioffset,M2l+joffset,M3 +koffset )     + &
                &  FI3(M1 +ioffset,M2 +joffset,M3u+koffset )     + &
                &  FI3(M1 +ioffset,M2 +joffset,M3l+koffset )     - &
                &  FI3(M1 +ioffset,M2 +joffset,M3 +koffset )*6.0D0
              !
              phi = phibar+FI3(M1+ioffset,M2+joffset,M3+koffset)/ctilde2! scalar field (background+perturbation)
              rho = FI3(M1+ioffset,M2,M3+koffset)+1.0D0                 ! density      (background+perturbation)
              !
           ENDIF
           !
           IF(csf_coupling .EQ.1) THEN
              A_phi = csf_beta*DEXP(csf_beta*phi)                       ! d  A(\varphi)/d\varphi
              A_pp  = csf_beta*A_phi/ctilde2                            ! d^2A(\varphi)/d\varphi^2
           ENDIF
           IF(csf_potential.EQ.1) THEN
              V_phi = -csf_alpha*csf_lambda**2*AEXPN**2
              V_phi = V_phi/(phi**(1.0D0+csf_alpha))                    ! (d  V(\varphi)/d\varphi  )/M^2_Pl
              V_pp  = -(1.0D0+csf_alpha)*V_phi/ctilde2/phi              ! (d^2V(\varphi)/d\varphi^2)/M^2_Pl
           ENDIF
           !
           L = L/dx2
           L = L-fct1*(A_phi*rho-A_phibar)
           L = L-     (V_phi    -V_phibar)
           dL = -sixodx2-fct1*rho*A_pp-V_pp
           !
           IF(ilevel.EQ.0) THEN
              ! solution on PM grid
!$OMP ATOMIC
              FI2(M1,M2,M3)                         = FI2(M1,M2,M3)                        -L/dL
           ELSE
              ! solution on multigrid
!$OMP ATOMIC
              FI3(M1+ioffset,M2+joffset,M3+koffset) = FI3(M1+ioffset,M2+joffset,M3+koffset)-L/dL
           ENDIF
        END DO
     END DO
  END DO

! IF(ilevel.NE.0.AND..NOT.redstep) THEN
!    OPEN(UNIT=27, FILE='test.txt', FORM='FORMATTED', STATUS='REPLACE')
!       DO M3=1,ngrid_level
!          WRITE(27,'(F20.9,F20.15,F20.15)') (DBLE(M3)-0.5d0)/DBLE(ngrid_level)*2.0D0*DACOS(-1.0D0),FI3(ngrid_level/2,ngrid_level/2,M3),FI3(ngrid_level/2+ioffset,ngrid_level/2,M3+koffset)
!       ENDDO
!    CLOSE(27)
!    STOP
! ENDIF

  CALL TimingMain(3,1)
    
END SUBROUTINE relaxation_iterations_csf

!-------------------------------------------------------------
!
! Calculate resudual and its RMS on the entire grid on ilevel.
!
!-------------------------------------------------------------
SUBROUTINE calculate_residual_csf(ilevel,res_PM_grid)
!-------------------------------------------------------------
  use Tools
  use ExtradofBackgroundData 

  integer,intent(in) :: ilevel
  real*8,intent(out) :: res_PM_grid
  
  real*8  :: ctilde,ctilde2,dx,dx2
  integer :: ngrid_level
  integer :: ioffset,joffset,koffset,koffset2
  real*8  :: fct1

  integer :: M1,M2,M3,M1l,M1u,M2l,M2u,M3l,M3u
  integer :: N1,N2,N3
  real*8  :: RES,OP,RES2

  integer :: ist
  real*8  :: phibar,A_phibar,V_phibar
  real*8  :: phi,rho,A_phi,V_phi

  IF(MG_test) WRITE(*,'(A,I5)') 'Calculate residual on level',levelmax-ilevel

  ! number of grid points on the coarse level
  ngrid_level = NGRID/2**ilevel
  
  ctilde  = 2.99792458D3*DBLE(NGRID)/Box
  ctilde2 = ctilde**2
  dx      = DBLE(2**ilevel)
  dx2     = dx*dx
  !
  fct1   = 3.0D0*Om/AEXPN
  !
  ist  = 1
  DO WHILE (BackgroundEvolution(ist,1).LT.AEXPN)
     ist=ist+1
  ENDDO
  !
  phibar =  BackgroundEvolution(ist-1,2)+ &
         & (BackgroundEvolution(ist,  2)-BackgroundEvolution(ist-1,2))/ &
         & (BackgroundEvolution(ist,  1)-BackgroundEvolution(ist-1,1))* &
         & (AEXPN                       -BackgroundEvolution(ist-1,1))
  !
  ! background coupling function 
  IF(csf_coupling.EQ.1) THEN
     A_phibar = csf_beta*DEXP(csf_beta*phibar)
  ENDIF
  !
  ! backgrund potential derivative: (dV/d\varphi)*a^2/M_Pl^2
  IF(csf_potential.EQ.1) THEN
     V_phibar = -csf_alpha*csf_lambda**2*AEXPN**2/(phibar**(1.0D0+csf_alpha))
  ENDIF
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
!$OMP PRIVATE (M1,M2,M3,M1u,M1l,M2u,M2l,M3u,M3l) &
!$OMP PRIVATE (phi,rho,A_phi,V_phi,OP)
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
              !
              phi = phibar+FI2(M1,M2,M3)/ctilde2
              rho = FI(M1,M2,M3)+1.0D0
              !
              IF(csf_coupling .EQ.1) THEN 
                 A_phi = csf_beta*DEXP(csf_beta*phi)
              ENDIF
              IF(csf_potential.EQ.1) THEN
                 V_phi = -csf_alpha*csf_lambda**2*AEXPN**2/(phi**(1.0D0+csf_alpha))
              ENDIF
              !
              OP = OP-fct1*(A_phi*rho-A_phibar)
              OP = OP-     (V_phi    -V_phibar)
              !
              FI3(M1,M2,M3) = OP
!             IF(DABS(OP).LT.1.0D-8) THEN
!                WRITE(*,'(I7,I7,I7,F20.15,F20.15,F20.15)') M1,M2,M3,OP,fct1*(A_phi*rho-A_phibar)
!             ENDIF
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
              OP = OP-FI3(M1,M2+joffset,M3+koffset2)
              !
              phi = phibar+FI3(M1+ioffset,M2+joffset,M3+koffset)/ctilde2
              rho = FI3(M1+ioffset,M2,M3+koffset)+1.0D0
              !
              IF(csf_coupling .EQ.1) THEN 
                 A_phi = csf_beta*DEXP(csf_beta*phi)
              ENDIF
              IF(csf_potential.EQ.1) THEN
                 V_phi = -csf_alpha*csf_lambda**2*AEXPN**2/(phi**(1.0D0+csf_alpha))
              ENDIF
              !
              OP = OP-fct1*(A_phi*rho-A_phibar)
              OP = OP-     (V_phi    -V_phibar)
              !
              FI3(M1+ioffset,M2,M3+koffset2) = OP 
              !
           END IF

        END DO
     END DO
  END DO

! OPEN(UNIT=27, FILE='test.txt', FORM='FORMATTED', STATUS='REPLACE')
!    DO M3=1,ngrid_level
!       WRITE(27,'(F20.9,F20.15)') (DBLE(M3)-0.5d0)/DBLE(ngrid_level)*2.0D0*DACOS(-1.0D0),FI3(ngrid_level/2,ngrid_level/2,M3)
!    ENDDO
! CLOSE(27)
! STOP
    
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
      
END SUBROUTINE calculate_residual_csf


!-------------------------------------------------------------
!
! Calculate restricted residual field on coarse level "ilevel"
! Note: "ilevel" means the grid ilevel coarer than the PM grid
!
!-------------------------------------------------------------
SUBROUTINE restrict_residual_csf(ilevel)
!-------------------------------------------------------------
  use Tools
  use ExtradofBackgroundData 

  integer :: ilevel
  real*8  :: ctilde,ctilde2,dx,dx2
  real*8  :: fct1
  integer :: ngrid_level 
  integer :: ioffset,joffset,koffset,ioffset2,joffset2,koffset2
 
  integer :: M1,M2,M3
  integer :: M1u,M1l,M2u,M2l,M3u,M3l 
  real*8  :: P,Q

  integer :: ist
  real*8  :: phibar,A_phibar,V_phibar
  real*8  :: phi,rho,A_phi,V_phi

  IF(MG_test) WRITE(*,'(A,I5)') 'Restrict residual to level',levelmax-ilevel

  ! number of grid points on the coarse level
  ngrid_level = NGRID/2**ilevel

  ctilde  = 2.99792458D3*DBLE(NGRID)/Box
  ctilde2 = ctilde**2
  dx      = 1.0d0  
  dx2     = dx*dx
  !
  fct1   = 3.0D0*Om/AEXPN
  !
  ist  = 1
  DO WHILE (BackgroundEvolution(ist,1).LT.AEXPN)
     ist=ist+1
  ENDDO
  !
  phibar =  BackgroundEvolution(ist-1,2)+ &
         & (BackgroundEvolution(ist,  2)-BackgroundEvolution(ist-1,2))/ &
         & (BackgroundEvolution(ist,  1)-BackgroundEvolution(ist-1,1))* &
         & (AEXPN                       -BackgroundEvolution(ist-1,1))
  !
  ! background coupling function 
  IF(csf_coupling.EQ.1) THEN
     A_phibar = csf_beta*DEXP(csf_beta*phibar)
  ENDIF
  !
  ! backgrund potential derivative: (dV/d\varphi)*a^2/M_Pl^2
  IF(csf_potential.EQ.1) THEN
     V_phibar = -csf_alpha*csf_lambda**2*AEXPN**2/(phibar**(1.0D0+csf_alpha))
  ENDIF
  !
  IF(ilevel.EQ.1) THEN
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3,M1u,M1l,M2u,M2l,M3u,M3l) &
!$OMP PRIVATE (phi,rho,A_phi,V_phi,P,Q) 
     DO M3=1,ngrid_level
        DO M2=1,ngrid_level
           DO M1=1,ngrid_level
              P = 0.0d0                                         ! Laplancian terms
              Q = 0.0d0                                         ! non-Laplacian terms
              ! accumulate contributions to the non-Laplacian term 
              ! of PDE from all 8 son cells
              !
              phi = phibar+FI2(2*M1-1,2*M2-1,2*M3-1)/ctilde2
              rho = 1.0D0 +FI (2*M1-1,2*M2-1,2*M3-1)
              A_phi = csf_beta*DEXP(csf_beta*phi)
              V_phi = -csf_alpha*csf_lambda**2*AEXPN**2/(phi**(1.0D0+csf_alpha))
              Q = Q+fct1*(A_phi*rho-A_phibar) &
                &  +     (V_phi    -V_phibar)
              !
              phi = phibar+FI2(2*M1  ,2*M2-1,2*M3-1)/ctilde2
              rho = 1.0D0 +FI (2*M1  ,2*M2-1,2*M3-1)
              A_phi = csf_beta*DEXP(csf_beta*phi)
              V_phi = -csf_alpha*csf_lambda**2*AEXPN**2/(phi**(1.0D0+csf_alpha))
              Q = Q+fct1*(A_phi*rho-A_phibar) &
                &  +     (V_phi    -V_phibar)
              !
              phi = phibar+FI2(2*M1-1,2*M2  ,2*M3-1)/ctilde2
              rho = 1.0D0 +FI (2*M1-1,2*M2  ,2*M3-1)
              A_phi = csf_beta*DEXP(csf_beta*phi)
              V_phi = -csf_alpha*csf_lambda**2*AEXPN**2/(phi**(1.0D0+csf_alpha))
              Q = Q+fct1*(A_phi*rho-A_phibar) &
                &  +     (V_phi    -V_phibar)
              !
              phi = phibar+FI2(2*M1  ,2*M2  ,2*M3-1)/ctilde2
              rho = 1.0D0 +FI (2*M1  ,2*M2  ,2*M3-1)
              A_phi = csf_beta*DEXP(csf_beta*phi)
              V_phi = -csf_alpha*csf_lambda**2*AEXPN**2/(phi**(1.0D0+csf_alpha))
              Q = Q+fct1*(A_phi*rho-A_phibar) &
                &  +     (V_phi    -V_phibar)
              !
              phi = phibar+FI2(2*M1-1,2*M2-1,2*M3  )/ctilde2
              rho = 1.0D0 +FI (2*M1-1,2*M2-1,2*M3  )
              A_phi = csf_beta*DEXP(csf_beta*phi)
              V_phi = -csf_alpha*csf_lambda**2*AEXPN**2/(phi**(1.0D0+csf_alpha))
              Q = Q+fct1*(A_phi*rho-A_phibar) &
                &  +     (V_phi    -V_phibar)
              !
              phi = phibar+FI2(2*M1  ,2*M2-1,2*M3  )/ctilde2
              rho = 1.0D0 +FI (2*M1  ,2*M2-1,2*M3  )
              A_phi = csf_beta*DEXP(csf_beta*phi)
              V_phi = -csf_alpha*csf_lambda**2*AEXPN**2/(phi**(1.0D0+csf_alpha))
              Q = Q+fct1*(A_phi*rho-A_phibar) &
                &  +     (V_phi    -V_phibar)
              !
              phi = phibar+FI2(2*M1-1,2*M2  ,2*M3  )/ctilde2
              rho = 1.0D0 +FI (2*M1-1,2*M2  ,2*M3  )
              A_phi = csf_beta*DEXP(csf_beta*phi)
              V_phi = -csf_alpha*csf_lambda**2*AEXPN**2/(phi**(1.0D0+csf_alpha))
              Q = Q+fct1*(A_phi*rho-A_phibar) &
                &  +     (V_phi    -V_phibar)
              !
              phi = phibar+FI2(2*M1  ,2*M2  ,2*M3  )/ctilde2
              rho = 1.0D0 +FI (2*M1  ,2*M2  ,2*M3  )
              A_phi = csf_beta*DEXP(csf_beta*phi)
              V_phi = -csf_alpha*csf_lambda**2*AEXPN**2/(phi**(1.0D0+csf_alpha))
              Q = Q+fct1*(A_phi*rho-A_phibar) &
                &  +     (V_phi    -V_phibar)
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

!    OPEN(UNIT=27, FILE='test.txt', FORM='FORMATTED', STATUS='REPLACE')
!       DO M3=1,ngrid_level
!          WRITE(27,'(F20.9,F20.15)') (DBLE(M3)-0.5d0)/DBLE(ngrid_level)*2.0D0*DACOS(-1.0D0),FI3(ngrid_level/2,ngrid_level/2,M3)
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
END SUBROUTINE restrict_residual_csf


!-------------------------------------------------------------
!
! Calculate the physical RHS of discrete PDE on coarser ilevel
!
!-------------------------------------------------------------
SUBROUTINE calculate_physical_right_hand_side_csf(ilevel)
!-------------------------------------------------------------
  use Tools
  use ExtradofBackgroundData

  integer :: ilevel
  real*8  :: ctilde,ctilde2,dx,dx2
  real*8  :: fct1
  integer :: ngrid_level
  integer :: ioffset,joffset,koffest,koffset2
  integer :: ist
  real*8  :: phibar,A_phibar,V_phibar
  real*8  :: phi,rho,A_phi,V_phi

  integer :: M1,M2,M3,M1l,M2l,M3l,M1u,M2u,M3u
  real*8  :: OP
  !
  IF(MG_test) WRITE(*,'(A,I5)') 'Calculate physical right-hand side on level',levelmax-ilevel
  !
  ! number of grid points on the coarse level
  ngrid_level = NGRID/2**ilevel
  !
  ctilde  = 2.99792458D3*DBLE(NGRID)/Box
  ctilde2 = ctilde**2
  dx      = DBLE(2**ilevel) 
  dx2     = dx*dx
  !
  fct1   = 3.0D0*Om/AEXPN
  !
  ist  = 1
  DO WHILE (BackgroundEvolution(ist,1).LT.AEXPN)
     ist=ist+1
  ENDDO
  !
  phibar =  BackgroundEvolution(ist-1,2)+ &
         & (BackgroundEvolution(ist,  2)-BackgroundEvolution(ist-1,2))/ &
         & (BackgroundEvolution(ist,  1)-BackgroundEvolution(ist-1,1))* &
         & (AEXPN                       -BackgroundEvolution(ist-1,1))
  !
  ! background coupling function 
  IF(csf_coupling.EQ.1) THEN
     A_phibar = csf_beta*DEXP(csf_beta*phibar)
  ENDIF
  !
  ! backgrund potential derivative: (dV/d\varphi)*a^2/M_Pl^2
  IF(csf_potential.EQ.1) THEN
     V_phibar = -csf_alpha*csf_lambda**2*AEXPN**2/(phibar**(1.0D0+csf_alpha))
  ENDIF
  !
  ioffset  = 2**(levelmax-ilevel)
  joffset  = 2**(levelmax-ilevel)
  koffset  = 2**(levelmax-ilevel)*(2**ilevel-2)
  koffset2 = 2**(levelmax-ilevel)*(2**ilevel-1)
  !
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3,M1u,M1l,M2u,M2l,M3u,M3l) &
!$OMP PRIVATE (phi,rho,A_phi,V_phi,OP) 
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
           !
           phi = phibar+FI3(        M1,joffset+M2,koffset+M3)/ctilde2
           rho = 1.0D0 +FI3(ioffset+M1,        M2,koffset+M3)
           !
           IF(csf_coupling .EQ.1) THEN
              A_phi = csf_beta*DEXP(csf_beta*phi)
           ENDIF
           !
           IF(csf_potential.EQ.1) THEN
              V_phi = -csf_alpha*csf_lambda**2*AEXPN**2/(phi**(1.0D0+csf_alpha))
           ENDIF
           !
           OP = OP-fct1*(A_phi*rho-A_phibar)
           OP = OP-     (V_phi    -V_phibar)
           OP = OP-FI3(M1,M2,M3+koffset)                        ! add restricted residual (note the minus sign!!)
           ! 
           FI3(M1,joffset+M2,koffset2+M3) = OP  
           ! 
        END DO
     END DO
  END DO

! OPEN(UNIT=27, FILE='test.txt', FORM='FORMATTED', STATUS='REPLACE')
!    DO M3=1,ngrid_level
!       WRITE(27,'(F20.9,F20.15)') (DBLE(M3)-0.5d0)/DBLE(ngrid_level)*2.0D0*DACOS(-1.0D0),FI3(ngrid_level/2,joffset+ngrid_level/2,koffset2+M3)
!    ENDDO
! CLOSE(27)
! STOP

  CALL TimingMain(3,1)

END SUBROUTINE calculate_physical_right_hand_side_csf
