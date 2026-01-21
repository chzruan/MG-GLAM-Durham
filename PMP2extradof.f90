!=================================================
!=================================================
!
!  Main scalar field program
!
!=================================================
!=================================================

! Baojiu: changed to double precision throughout, but we also need
! a single-precision version as double precision may be unnecessary

SUBROUTINE scalar_field(i_step)
use Tools

  integer :: i_step,iter_count
  integer :: M1,M2,M3  

  real*8  :: res_PM_grid
  logical :: can_exit,exit_n_warn

  real*8  :: R_bg,R0_bg,fR_bg,u_bg
  real*8  :: res_old_1,res_old_2 

  real*8 :: FI2_max

  res_PM_grid = 0.0D0
  res_old_1   = 1.0D5
  res_old_2   = 1.0D10
  
  iter_count  = 1
  can_exit    = .FALSE.
  exit_n_warn = .FALSE.

  ! initialise scalar field values on PM grid
  IF(i_step.EQ.1) CALL initialise_extradof

   IF(MG_model.EQ.3 .AND. AEXPN.LT.sym_astar) THEN
      ! symmetron model
      ! do not need to solve the scalar field equation before the symmetry breaking time sym_astar
      RETURN 
   ENDIF 

  main_iteration_loop: DO
     !
     ! prepare initial guess for scalar field 
     ! for f(R) gravity part only
     IF(i_step.GT.1 .AND. iter_count.EQ.1) THEN
        IF(MG_model.EQ.1) THEN
           R_bg  = 3.0D0*(Om/AEXPN**3+4.0D0*OmL)
           R0_bg = 3.0D0*(Om         +4.0D0*OmL)
           fR_bg = fR0*(R0_bg/R_bg)**(fr_n+1)
           u_bg  = (DABS(fR_bg))**(1.0D0/(1.0D0+DBLE(fr_n)))
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3)
           DO M3=1,NGRID
              DO M2=1,NGRID
                 DO M1=1,NGRID
                    FI2(M1,M2,M3) = FI2(M1,M2,M3)*u_bg
                 END DO
              END DO
           END DO
        ENDIF ! IF(MG_model.EQ.1)

!        IF(MG_model.EQ.2) THEN
!!$OMP PARALLEL DO DEFAULT(SHARED) &
!!$OMP PRIVATE (M1,M2,M3)
!           DO M3=1,NGRID
!              DO M2=1,NGRID
!                 DO M1=1,NGRID
!                    FI2(M1,M2,M3) = 0.0D0
!                 END DO
!              END DO
!           END DO
!        ENDIF

         IF(MG_model.EQ.3 .AND. iter_count.EQ.1) THEN
            ! entering this part means AEXPN > sys_astar 
            u_bg = DSQRT(1.0D0-(sym_astar/AEXPN)**3)
            FI2(:,:,:) = u_bg
         ENDIF !IF(MG_model.EQ.3 .AND. iter_count.EQ.1)
     
     ENDIF ! IF(i_step.GT.1.AND.iter_count.EQ.1)
  

     ! multigrid solver
     IF(Vcycle .EQ. 1) Call V_cycle(0)
     IF(Fcycle .EQ. 1) Call F_cycle(0)
     IF(Wcycle .EQ. 1) Call W_cycle(0)


     ! update residual on PM grid after post smoothing
     SELECT CASE(MG_model)
        CASE(1)
           CALL calculate_residual_fR (0,res_PM_grid)
        CASE(2)
           CALL calculate_residual_DGP(0,res_PM_grid) 
        CASE(3)
           CALL calculate_residual_sym(0,res_PM_grid) 
        CASE(4)
           CALL calculate_residual_kmf(0,res_PM_grid) 
        CASE(5)
           CALL calculate_residual_csf(0,res_PM_grid) 
        CASE DEFAULT
     END SELECT
     WRITE(*,'(A,F20.12,A,I5,A)') 'The residual is equal to ',res_PM_grid, ' after ',iter_count, ' V-cycles.'
     !
     ! convergence criterion 1: residual small enough
     IF(res_PM_grid.LE.res_conv) THEN
        can_exit = .TRUE.
     ENDIF
     ! convergence criterion 2: residual no longer reduces
     IF((iter_count.GT.2)) THEN
        IF((DABS(res_PM_grid/res_old_1-1.0D0).LE.5.0D-3)) can_exit = .TRUE.
        IF((DABS(res_PM_grid/res_old_1-1.0D0).LE.1.0D-2) .AND. &
         & (DABS(res_old_1  /res_old_2-1.0D0).LE.1.0D-2)) can_exit = .TRUE.
     ENDIF
     ! convergence criterion 3: maximum iteration number exceeded
     IF(iter_count.GE.iter_count_max) THEN
        can_exit    = .TRUE.
        exit_n_warn = .TRUE.
     ENDIF
     !
     ! Output residual information at exit
     IF(can_exit) THEN
        WRITE(*,'(A,F20.12)') 'Residual on the PM grid = ',res_PM_grid
      !   WRITE(*,'(A,F20.12)') 'Convergence criterion = ',res_conv
        IF(exit_n_warn) WRITE(*,'(A,I5,A)') 'ATTENTION: not converged after a maximum of ',iter_count,' iterations'
     ENDIF
     !
     ! Update scalar field solution
     IF(can_exit) THEN
        IF(MG_model.EQ.1) THEN
           R_bg  = 3.0D0*(Om/AEXPN**3+4.0D0*OmL)
           R0_bg = 3.0D0*(Om         +4.0D0*OmL)
           fR_bg = fR0*(R0_bg/R_bg)**(fr_n+1)
           u_bg  = (DABS(fR_bg))**(-1.0D0/(1.0D0+DBLE(fr_n)))
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3)
           DO M3=1,NGRID
              DO M2=1,NGRID
                 DO M1=1,NGRID
                    FI2(M1,M2,M3) = FI2(M1,M2,M3)*u_bg
                 END DO
              END DO
           END DO
        ENDIF
!        IF(MG_model.EQ.2) THEN
!           CALL phi_mean_DGP
!!$OMP PARALLEL DO DEFAULT(SHARED) &
!!$OMP PRIVATE (M1,M2,M3)
!           DO M3=1,NGRID
!              DO M2=1,NGRID
!                 DO M1=1,NGRID
!                    FI2(M1,M2,M3) = FI2(M1,M2,M3) - phi_mean
!                 END DO
!              END DO
!           END DO
!        ENDIF
        !
        EXIT
     ENDIF

     ! increment cycle count
     iter_count = iter_count+1
     !
     ! update old residual values
     res_old_2 = res_old_1    ! residual of the previous previous step
     res_old_1 = res_PM_grid  ! residual of the previous step
     
  END DO main_iteration_loop

  IF(MG_test .EQ. 1) THEN
     WRITE(*, '(A,I2,A)') 'scalar field step',istep,' done.'
     CALL WriteScalarFieldToFile
     WRITE(*, '(A,I2,A)') 'save FI2 to disk', istep,' done.'
  ENDIF

end SUBROUTINE scalar_field

!--------------------------------------------
SUBROUTINE initialise_extradof
!--------------------------------------------
USE Tools

  integer :: i_step
  integer :: M1,M2,M3
  integer :: ngridmax                                                           ! Baojiu 2025-12-17
  real*8  :: r, max_n = 0.05D0, min_n = -0.05D0
  real*8  :: R_bg,R0_bg,fR_bg,u_bg,u_bg2
  
  ! Output general information
  SELECT CASE(MG_model)  
     CASE(1)
        WRITE(*,'(A)'      ) 'We are simulating f(R) gravity'
        WRITE(*,'(A,E12.6)') 'Parameter f_R0 = ',fR0
        WRITE(*,'(A,I5)'   ) 'Parameter n = ',fr_n
     CASE(2)
        IF(N_branch) WRITE(*,'(A)') 'We are simulating DGP gravity (normal branch)'
        IF(S_branch) WRITE(*,'(A)') 'We are simulating DGP gravity (self-accelerating branch)'
        WRITE(*,'(A,F12.6)') 'Parameter H0*r_c = ',H0rc
     CASE(3)
        WRITE(*,'(A)'      ) 'We are simulating symmetron model'
        WRITE(*,'(A,E12.6)') 'Parameter a_* = ',sym_astar
        WRITE(*,'(A,E12.6)') 'Parameter xi = ',sym_xi
        WRITE(*,'(A,E12.6)') 'Parameter beta_* = ',sym_beta
     CASE(4)
        IF(power_law) THEN
           WRITE(*,'(A)') 'We are simulating k-mouflage gravity (power-law type)'
           WRITE(*,'(A,I5)') 'Parameter n = ',kmf_n
           WRITE(*,'(A,E12.6)') 'Parameter K0 = ',kmf_K0
        ENDIF
        IF(born_infeld) THEN
           WRITE(*,'(A)') 'We are simulating k-mouflage gravity (Born-Infeld type)'
           WRITE(*,'(A)') 'This is currently not supported. Please try again'
           STOP
        ENDIF
        WRITE(*,'(A,E12.6)') 'Parameter beta = ',kmf_beta
     CASE(5)
        WRITE(*,'(A)'      ) 'We are simulating coupled scalar field model'
        IF(csf_coupling .EQ.1) THEN
           WRITE(*,*) 'The coupling function is exponential'
        ELSE
           WRITE(*,*) 'The coupling function is currently not supported. Please try again'
           STOP
        ENDIF
        IF(csf_potential.EQ.1) THEN
           WRITE(*,*) 'The scalar field potential is inverse power law'
        ELSE
           WRITE(*,*) 'The scalar field potential is currently not supported. Please try again'
           STOP
        ENDIF
        WRITE(*,'(A,E12.6)') 'Parameter alpha = ',csf_alpha
        WRITE(*,'(A,E12.6)') 'Parameter beta = ',csf_beta
     CASE DEFAULT  
        WRITE(*,'(A)') 'Model currently unsupported. Please specify again.'
        STOP
  END SELECT

  ! Baojiu: initialise these global variables here once and for all
  levelmax = 0   
  DO WHILE(2**levelmax.LT.NGRID) 
     levelmax = levelmax + 1
  END DO  
  ngridmax = NGRID                                                              ! Baojiu 2025 12 17
  DO WHILE(MOD(ngridmax,2).EQ.0)                                                ! Baojiu 2025 12 17
     ngridmax = ngridmax/2                                                      ! Baojiu 2025 12 17
  END DO                                                                        ! Baojiu 2025 12 17
  levelmin = 0                                                                  ! Baojiu 2025 12 17
  DO WHILE(2**levelmin.LT.ngridmax)                                             ! Baojiu 2025 12 17
     levelmin = levelmin + 1                                                    ! Baojiu 2025 12 17
  END DO                                                                        ! Baojiu 2025 12 17
  levelmin = levelmin+1                                                         ! Baojiu 2025-12-17

  WRITE(*,'(A,I5)') 'Maximum multigrid level:',levelmax
  WRITE(*,'(A,I5)') 'Minimum multigrid level:',levelmin
  ! Baojiu: end of modifications

  SELECT CASE(MG_model)
     ! 1 - f(R) gravity
     CASE(1)
        R_bg  = 3.0D0*(Om/AEXPN**3+4.0D0*OmL)
        R0_bg = 3.0D0*(Om         +4.0D0*OmL)
        fR_bg = fR0*(R0_bg/R_bg)**(fr_n+1)
        u_bg  = (DABS(fR_bg))**(1.0D0/(1.0D0+DBLE(fr_n)))
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3)
        DO M3=1,NGRID
           DO M2=1,NGRID
              DO M1=1,NGRID
                 FI2(M1,M2,M3) = u_bg
              END DO
           END DO
        END DO
     ! 2 - DGP model
     CASE(2)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3)
        DO M3=1,NGRID
           DO M2=1,NGRID
              DO M1=1,NGRID
!                IF(MG_test .EQ. 1) THEN
!                   CALL RANDOM_NUMBER(r)
!                   FI2(M1,M2,M3) = r * (max_n - min_n) + min_n
!                ELSE
                    FI2(M1,M2,M3) = 0.0D0
!                END IF
              END DO
           END DO
        END DO
     ! 3 - symmetron model
     CASE(3)
        IF(AEXPN.GT.sym_astar) THEN
           u_bg = DSQRT(1.0D0-(sym_astar/AEXPN)**3)
        ELSE
           u_bg = 0.0D0
        ENDIF
        FI2(:,:,:) = u_bg
!! $OMP PARALLEL DO DEFAULT(SHARED) &
!! $OMP PRIVATE (M1,M2,M3)
      !   DO M3=1,NGRID
      !      DO M2=1,NGRID
      !         DO M1=1,NGRID
      !            FI2(M1,M2,M3) = u_bg
      !         END DO
      !      END DO
      !   END DO
     ! 4 - k-mouflage model
     CASE(4)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3)
        DO M3=1,NGRID
           DO M2=1,NGRID
              DO M1=1,NGRID
                 FI2(M1,M2,M3) = 0.0D0
              END DO
           END DO
        END DO
     ! 5 - coupled scalar field model
     CASE(5)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3)
        DO M3=1,NGRID
           DO M2=1,NGRID
              DO M1=1,NGRID
                 FI2(M1,M2,M3) = 0.0D0
              END DO
           END DO
        END DO
     ! default case   
     CASE DEFAULT
     !
  END SELECT

!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3)
  DO M3=1,NGRID
     DO M2=1,NGRID
        DO M1=1,NGRID
           FI3(M1,M2,M3) = 0.0D0
        END DO
     END DO
  END DO

END SUBROUTINE initialise_extradof

!-------------------------------------------------------------
!
! arrangement of V-cycle multigrid relaxation (recursive call)
!
! Note: ilevel runs from 0, which means relaxation on PM grid.
!-------------------------------------------------------------
RECURSIVE SUBROUTINE V_cycle(ilevel)
!-------------------------------------------------------------
USE Tools

  integer :: ilevel
  integer :: Nsmth
  integer :: ismth
  real*8  :: res_PM_grid

  Nsmth    = 2                                                 ! max number of pre/post smoothing
  
  IF(MG_test .EQ. 1) THEN
    WRITE(*,'(A,I5)') 'V-Cycle. Current multigrid level:', levelmax-ilevel
  ENDIF

  ! exit if already on the coarsest multigrid level
  IF(levelmax-ilevel.LE.levelmin) THEN
     IF(MG_test .EQ. 1) WRITE(*,'(A)') 'Reached minimum multigrid level. Time to jump out!'
     ! pre-smoothing
     SELECT CASE(MG_model)
        CASE(1)
           ! Gauss-Seidel iterations
           DO i_smth=1,Nsmth
              CALL relaxation_iterations_fR (ilevel,.TRUE. )
              CALL relaxation_iterations_fR (ilevel,.FALSE.)
           END DO
           !
        CASE(2)
           ! Gauss-Seidel iterations
           DO i_smth=1,Nsmth
              CALL relaxation_iterations_DGP(ilevel,.TRUE. )
              CALL relaxation_iterations_DGP(ilevel,.FALSE.)
           END DO
           !
        CASE(3)
           ! Gauss-Seidel iterations
           DO i_smth=1,Nsmth
              CALL relaxation_iterations_sym(ilevel,.TRUE. )
              CALL relaxation_iterations_sym(ilevel,.FALSE.)
           END DO
           !
        CASE(4)
           ! Gauss-Seidel iterations
           DO i_smth=1,Nsmth
              CALL relaxation_iterations_kmf(ilevel,.TRUE. )
              CALL relaxation_iterations_kmf(ilevel,.FALSE.)
           END DO
           !
           !
        CASE(5)
           ! Gauss-Seidel iterations
           DO i_smth=1,Nsmth
              CALL relaxation_iterations_csf(ilevel,.TRUE. )
              CALL relaxation_iterations_csf(ilevel,.FALSE.)
           END DO
           !
        CASE DEFAULT
           ! error message
           WRITE(*,'(A)') 'Model currently unsupported. Please specify again.'
           STOP
           !
     END SELECT 
     !     
     RETURN
  ENDIF

  ! pre-smooting
  SELECT CASE(MG_model)
     CASE(1)
        ! Gauss-Seidel iterations
        DO i_smth=1,Nsmth
           CALL relaxation_iterations_fR(ilevel,.TRUE. )
           CALL relaxation_iterations_fR(ilevel,.FALSE.)
        END DO
        ! calculate residual on fine level
        CALL calculate_residual_fR(ilevel,res_PM_grid)
        ! restrict residual to coarse level
        CALL restrict_residual_fR(ilevel+1)
        !
     CASE(2)
        ! Gauss-Seidel iterations
        DO i_smth=1,Nsmth
           CALL relaxation_iterations_DGP(ilevel,.TRUE. )
           CALL relaxation_iterations_DGP(ilevel,.FALSE.)
        END DO
        ! calculate residual on fine level
        CALL calculate_residual_DGP(ilevel,res_PM_grid)
        ! restrict residual to coarse level
        CALL restrict_residual_DGP(ilevel+1)
        !
     CASE(3)
        ! Gauss-Seidel iterations
        DO i_smth=1,Nsmth
           CALL relaxation_iterations_sym(ilevel,.TRUE. )
           CALL relaxation_iterations_sym(ilevel,.FALSE.)
        END DO
        ! calculate residual on fine level
        CALL calculate_residual_sym(ilevel,res_PM_grid)
        ! restrict residual to coarse level
        CALL restrict_residual_sym(ilevel+1)
        !
     CASE(4)
        ! Gauss-Seidel iterations
        DO i_smth=1,Nsmth
           CALL relaxation_iterations_kmf(ilevel,.TRUE. )
           CALL relaxation_iterations_kmf(ilevel,.FALSE.)
        END DO
        ! calculate residual on fine level
        CALL calculate_residual_kmf(ilevel,res_PM_grid)
        ! restrict residual to coarse level
        CALL restrict_residual_kmf(ilevel+1)
        !
     CASE(5)
        ! Gauss-Seidel iterations
        DO i_smth=1,Nsmth
           CALL relaxation_iterations_csf(ilevel,.TRUE. )
           CALL relaxation_iterations_csf(ilevel,.FALSE.)
        END DO
        ! calculate residual on fine level
        CALL calculate_residual_csf(ilevel,res_PM_grid)
        ! restrict residual to coarse level
        CALL restrict_residual_csf(ilevel+1)
        !
     CASE DEFAULT
        ! error message
        WRITE(*,'(A)') 'Model currently unsupported. Please specify again.'
        STOP
        !   
  END SELECT

  ! restrict scalar field to coarse level
  CALL restrict_scalar_field(ilevel+1)
  
  ! restrict density field to coarse level
  CALL restrict_density(ilevel+1)
  
  ! calculate physical RHS on coarse level
  SELECT CASE(MG_model)
     CASE(1)
        CALL calculate_physical_right_hand_side_fR (ilevel+1)
     CASE(2)
        CALL calculate_physical_right_hand_side_DGP(ilevel+1)
     CASE(3)
        CALL calculate_physical_right_hand_side_sym(ilevel+1)
     CASE(4)
        CALL calculate_physical_right_hand_side_kmf(ilevel+1)
     CASE(5)
        CALL calculate_physical_right_hand_side_csf(ilevel+1)
     CASE DEFAULT
        WRITE(*,'(A)') 'Model currently unsupported. Please specify again.'
        STOP
  END SELECT

  ! recursive call to next coarser multigrid level
  CALL V_cycle(ilevel+1)

  ! prolongation and correction to fine level solution
  CALL correct_solution(ilevel+1)

  ! post-smoothing
  SELECT CASE(MG_model)
     CASE(1)
        ! Gauss-Seidel iterations
        DO i_smth=1,Nsmth
           CALL relaxation_iterations_fR (ilevel,.TRUE. )
           CALL relaxation_iterations_fR (ilevel,.FALSE.)
        END DO
        !
     CASE(2)
        ! Gauss-Seidel iterations
        DO i_smth=1,Nsmth
           CALL relaxation_iterations_DGP(ilevel,.TRUE. )
           CALL relaxation_iterations_DGP(ilevel,.FALSE.)
        END DO
        !
     CASE(3)
        ! Gauss-Seidel iterations
        DO i_smth=1,Nsmth
           CALL relaxation_iterations_sym(ilevel,.TRUE. )
           CALL relaxation_iterations_sym(ilevel,.FALSE.)
        END DO
        !
     CASE(4)
        ! Gauss-Seidel iterations
        DO i_smth=1,Nsmth
           CALL relaxation_iterations_kmf(ilevel,.TRUE. )
           CALL relaxation_iterations_kmf(ilevel,.FALSE.)
        END DO
        !
     CASE(5)
        ! Gauss-Seidel iterations
        DO i_smth=1,Nsmth
           CALL relaxation_iterations_csf(ilevel,.TRUE. )
           CALL relaxation_iterations_csf(ilevel,.FALSE.)
        END DO
        !
     CASE DEFAULT
        ! error message
        WRITE(*,'(A)') 'Model currently unsupported. Please specify again.'
        STOP
        !
  END SELECT

END SUBROUTINE V_cycle

!-------------------------------------------------------------
!
! arrangement of F-cycle multigrid relaxation (recursive call)
!
! Note: ilevel runs from 0, which means relaxation on PM grid.
!-------------------------------------------------------------
RECURSIVE SUBROUTINE F_cycle(ilevel)
!-------------------------------------------------------------
USE Tools

  integer :: ilevel
  integer :: Nsmth
  integer :: ismth
  real*8  :: res_PM_grid

  ! Baojiu: removed the iteration, and its associated variables

  Nsmth = 2                                                     ! max number of pre/post smoothing

  write(*,'(A)') ' '
  WRITE(*,'(A,I5)') 'F-cycle. Current multigrid level:', levelmax-ilevel

  ! exit if already on the coarsest multigrid level
  IF(levelmax-ilevel.LE.levelmin) THEN
     IF(MG_test .EQ. 1) WRITE(*,'(A)') 'Reached minimum multigrid level. Time to jump out!'
     ! re-smoothing
     SELECT CASE(MG_model)
        CASE(1)
           ! Gauss-Seidel iterations
           DO i_smth=1,Nsmth
              CALL relaxation_iterations_fR (ilevel,.TRUE. )
              CALL relaxation_iterations_fR (ilevel,.FALSE.)
           END DO
           !
        CASE(2)
           ! Gauss-Seidel iterations
           DO i_smth=1,Nsmth
              CALL relaxation_iterations_DGP(ilevel,.TRUE. )
              CALL relaxation_iterations_DGP(ilevel,.FALSE.)
           END DO
           !
        CASE(3)
           ! Gauss-Seidel iterations
           DO i_smth=1,Nsmth
              CALL relaxation_iterations_sym(ilevel,.TRUE. )
              CALL relaxation_iterations_sym(ilevel,.FALSE.)
           END DO
           !
        CASE(4)
           ! Gauss-Seidel iterations
           DO i_smth=1,Nsmth
              CALL relaxation_iterations_kmf(ilevel,.TRUE. )
              CALL relaxation_iterations_kmf(ilevel,.FALSE.)
           END DO
           !
        CASE(5)
           ! Gauss-Seidel iterations
           DO i_smth=1,Nsmth
              CALL relaxation_iterations_csf(ilevel,.TRUE. )
              CALL relaxation_iterations_csf(ilevel,.FALSE.)
           END DO
           !
        CASE DEFAULT
           ! error message
           WRITE(*,'(A)') 'Model currently unsupported. Please specify again.'
           STOP
           !
     END SELECT
     !
     RETURN
  END IF
  
  ! pre-smoothing
  SELECT CASE(MG_model)
     CASE(1)
        ! Gauss-Seidel iterations
        DO i_smth=1,Nsmth
           CALL relaxation_iterations_fR(ilevel,.TRUE. )
           CALL relaxation_iterations_fR(ilevel,.FALSE.)
        END DO
        ! calculate residual on fine level
        CALL calculate_residual_fR(ilevel,res_PM_grid)
        ! restrict residual to coarse level
        CALL restrict_residual_fR(ilevel+1)
        !
     CASE(2)
        ! Gauss-Seidel iterations
        DO i_smth=1,Nsmth
           CALL relaxation_iterations_DGP(ilevel,.TRUE. )
           CALL relaxation_iterations_DGP(ilevel,.FALSE.)
        END DO
        ! calculate residual on fine level
        CALL calculate_residual_DGP(ilevel,res_PM_grid)
        ! restrict residual to coarse level
        CALL restrict_residual_DGP(ilevel+1)
        !
     CASE(3)
        ! Gauss-Seidel iterations
        DO i_smth=1,Nsmth
           CALL relaxation_iterations_sym(ilevel,.TRUE. )
           CALL relaxation_iterations_sym(ilevel,.FALSE.)
        END DO
        ! calculate residual on fine level
        CALL calculate_residual_sym(ilevel,res_PM_grid)
        ! restrict residual to coarse level
        CALL restrict_residual_sym(ilevel+1)
        !
     CASE(4)
        ! Gauss-Seidel iterations
        DO i_smth=1,Nsmth
           CALL relaxation_iterations_kmf(ilevel,.TRUE. )
           CALL relaxation_iterations_kmf(ilevel,.FALSE.)
        END DO
        ! calculate residual on fine level
        CALL calculate_residual_kmf(ilevel,res_PM_grid)
        ! restrict residual to coarse level
        CALL restrict_residual_kmf(ilevel+1)
        !
     CASE(5)
        ! Gauss-Seidel iterations
        DO i_smth=1,Nsmth
           CALL relaxation_iterations_csf(ilevel,.TRUE. )
           CALL relaxation_iterations_csf(ilevel,.FALSE.)
        END DO
        ! calculate residual on fine level
        CALL calculate_residual_csf(ilevel,res_PM_grid)
        ! restrict residual to coarse level
        CALL restrict_residual_csf(ilevel+1)
        !
     CASE DEFAULT
        ! error message
        WRITE(*,'(A)') 'Model currently unsupported. Please specify again.'
        STOP
        !   
  END SELECT

  ! restrict scalar field to coarse level
  CALL restrict_scalar_field(ilevel+1)

  ! restrict density field to coarse level
  CALL restrict_density(ilevel+1)

  ! calculate physical RHS on coarse level
  SELECT CASE(MG_model)
     CASE(1)
        CALL calculate_physical_right_hand_side_fR (ilevel+1)
     CASE(2)
        CALL calculate_physical_right_hand_side_DGP(ilevel+1)
     CASE(3)
        CALL calculate_physical_right_hand_side_sym(ilevel+1)
     CASE(4)
        CALL calculate_physical_right_hand_side_kmf(ilevel+1)
     CASE(5)
        CALL calculate_physical_right_hand_side_csf(ilevel+1)
     CASE DEFAULT
        WRITE(*,'(A)') 'Model currently unsupported. Please specify again.'
        STOP
  END SELECT

  ! recursive call to next coarser multigrid level
  CALL F_cycle(ilevel+1)

  ! prolongation and correction to fine level solution
  CALL correct_solution(ilevel+1)

  ! post-smoothing
  SELECT CASE(MG_model)
     CASE(1)
        ! Gauss-Seidel iterations
        DO i_smth=1,Nsmth
           CALL relaxation_iterations_fR(ilevel,.TRUE. )
           CALL relaxation_iterations_fR(ilevel,.FALSE.)
        END DO
        ! calculate residual on fine level
        CALL calculate_residual_fR(ilevel,res_PM_grid)
        ! restrict residual to coarse level
        CALL restrict_residual_fR(ilevel+1)
        !
     CASE(2)
        ! Gauss-Seidel iterations
        DO i_smth=1,Nsmth
           CALL relaxation_iterations_DGP(ilevel,.TRUE. )
           CALL relaxation_iterations_DGP(ilevel,.FALSE.)
        END DO
        ! calculate residual on fine level
        CALL calculate_residual_DGP(ilevel,res_PM_grid)
        ! restrict residual to coarse level
        CALL restrict_residual_DGP(ilevel+1)
        !
     CASE(3)
        ! Gauss-Seidel iterations
        DO i_smth=1,Nsmth
           CALL relaxation_iterations_sym(ilevel,.TRUE. )
           CALL relaxation_iterations_sym(ilevel,.FALSE.)
        END DO
        ! calculate residual on fine level
        CALL calculate_residual_sym(ilevel,res_PM_grid)
        ! restrict residual to coarse level
        CALL restrict_residual_sym(ilevel+1)
        !
     CASE(4)
        ! Gauss-Seidel iterations
        DO i_smth=1,Nsmth
           CALL relaxation_iterations_kmf(ilevel,.TRUE. )
           CALL relaxation_iterations_kmf(ilevel,.FALSE.)
        END DO
        ! calculate residual on fine level
        CALL calculate_residual_kmf(ilevel,res_PM_grid)
        ! restrict residual to coarse level
        CALL restrict_residual_kmf(ilevel+1)
        !
     CASE(5)
        ! Gauss-Seidel iterations
        DO i_smth=1,Nsmth
           CALL relaxation_iterations_csf(ilevel,.TRUE. )
           CALL relaxation_iterations_csf(ilevel,.FALSE.)
        END DO
        ! calculate residual on fine level
        CALL calculate_residual_csf(ilevel,res_PM_grid)
        ! restrict residual to coarse level
        CALL restrict_residual_csf(ilevel+1)
        !
     CASE DEFAULT
        ! error message
        WRITE(*,'(A)') 'Model currently unsupported. Please specify again.'
        STOP
        !
  END SELECT

  ! restrict scalar field to coarse level
  CALL restrict_scalar_field(ilevel+1)

  ! restrict density field to coarse level
  CALL restrict_density(ilevel+1)

  ! calculate physical RHS on coarse level
  SELECT CASE(MG_model)
     CASE(1)
        CALL calculate_physical_right_hand_side_fR (ilevel+1)
     CASE(2)
        CALL calculate_physical_right_hand_side_DGP(ilevel+1)
     CASE(3)
        CALL calculate_physical_right_hand_side_sym(ilevel+1)
     CASE(4)
        CALL calculate_physical_right_hand_side_kmf(ilevel+1)
     CASE(5)
        CALL calculate_physical_right_hand_side_csf(ilevel+1)
     CASE DEFAULT
        WRITE(*,'(A)') 'Model currently unsupported. Please specify again.'
        STOP
  END SELECT

  ! recursive call to next coarser multigrid level
  CALL V_cycle(ilevel+1)

  ! prolongation and correction to fine level solution
  CALL correct_solution(ilevel+1)

  SELECT CASE(MG_model)
     CASE(1)
        ! post-smoothing
        DO i_smth=1,Nsmth
           CALL relaxation_iterations_fR (ilevel,.TRUE. )
           CALL relaxation_iterations_fR (ilevel,.FALSE.)
        END DO
        !
     CASE(2) 
        ! post-smoothing
        DO i_smth=1,Nsmth
           CALL relaxation_iterations_DGP(ilevel,.TRUE. )
           CALL relaxation_iterations_DGP(ilevel,.FALSE.)
        END DO
        !
     CASE(3) 
        ! post-smoothing
        DO i_smth=1,Nsmth
           CALL relaxation_iterations_sym(ilevel,.TRUE. )
           CALL relaxation_iterations_sym(ilevel,.FALSE.)
        END DO
        !
     CASE(4) 
        ! post-smoothing
        DO i_smth=1,Nsmth
           CALL relaxation_iterations_kmf(ilevel,.TRUE. )
           CALL relaxation_iterations_kmf(ilevel,.FALSE.)
        END DO
        !
     CASE(5) 
        ! post-smoothing
        DO i_smth=1,Nsmth
           CALL relaxation_iterations_csf(ilevel,.TRUE. )
           CALL relaxation_iterations_csf(ilevel,.FALSE.)
        END DO
        !
     CASE DEFAULT
        ! error message
        WRITE(*,'(A)') 'Model currently unsupported. Please specify again.'
        STOP
        !
  END SELECT

END SUBROUTINE F_cycle

!-------------------------------------------------------------
!
! arrangement of W-cycle multigrid relaxation (recursive call)
! Note: ilevel runs from 0, which means relaxation on PM grid.
!
!-------------------------------------------------------------
!-------------------------------------------------------------
RECURSIVE SUBROUTINE W_cycle(ilevel)
!-------------------------------------------------------------
USE Tools

  integer :: ilevel
  integer :: Nsmth
  integer :: ismth
  real*8  :: res_PM_grid

  ! Baojiu: removed the iteration, and its associated variables

  Nsmth = 2                                                     ! max number of pre/post smoothing

  write(*,'(A)') ' '
  WRITE(*,'(A,I5)') 'W-cycle. Current multigrid level:', levelmax-ilevel

  ! exit if already on the coarsest multigrid level
  IF(levelmax-ilevel.LE.levelmin) THEN
     IF(MG_test .EQ. 1) WRITE(*,'(A)') 'Reached minimum multigrid level. Time to jump out!'
     ! re-smoothing
     SELECT CASE(MG_model)
        CASE(1)
           ! Gauss-Seidel iterations
           DO i_smth=1,Nsmth
              CALL relaxation_iterations_fR (ilevel,.TRUE. )
              CALL relaxation_iterations_fR (ilevel,.FALSE.)
           END DO
        CASE(2)
           ! Gauss-Seidel iterations
           DO i_smth=1,Nsmth
              CALL relaxation_iterations_DGP(ilevel,.TRUE. )
              CALL relaxation_iterations_DGP(ilevel,.FALSE.)
           END DO
        CASE(3)
           ! Gauss-Seidel iterations
           DO i_smth=1,Nsmth
              CALL relaxation_iterations_sym(ilevel,.TRUE. )
              CALL relaxation_iterations_sym(ilevel,.FALSE.)
           END DO
        CASE(4)
           ! Gauss-Seidel iterations
           DO i_smth=1,Nsmth
              CALL relaxation_iterations_kmf(ilevel,.TRUE. )
              CALL relaxation_iterations_kmf(ilevel,.FALSE.)
           END DO
        CASE(5)
           ! Gauss-Seidel iterations
           DO i_smth=1,Nsmth
              CALL relaxation_iterations_csf(ilevel,.TRUE. )
              CALL relaxation_iterations_csf(ilevel,.FALSE.)
           END DO
        CASE DEFAULT
           ! error message
           WRITE(*,'(A)') 'Model currently unsupported. Please specify again.'
           STOP
     END SELECT
     !
     RETURN
  END IF

  ! pre-smooting
  SELECT CASE(MG_model)
     CASE(1)
        ! Gauss-Seidel iterations
        DO i_smth=1,Nsmth
           CALL relaxation_iterations_fR(ilevel,.TRUE. )
           CALL relaxation_iterations_fR(ilevel,.FALSE.)
        END DO
        ! calculate residual on fine level
        CALL calculate_residual_fR(ilevel,res_PM_grid)
        ! restrict residual to coarse level
        CALL restrict_residual_fR(ilevel+1)
        !
     CASE(2)
        ! Gauss-Seidel iterations
        DO i_smth=1,Nsmth
           CALL relaxation_iterations_DGP(ilevel,.TRUE. )
           CALL relaxation_iterations_DGP(ilevel,.FALSE.)
        END DO
        ! calculate residual on fine level
        CALL calculate_residual_DGP(ilevel,res_PM_grid)
        ! restrict residual to coarse level
        CALL restrict_residual_DGP(ilevel+1)
        !
     CASE(3)
        ! Gauss-Seidel iterations
        DO i_smth=1,Nsmth
           CALL relaxation_iterations_sym(ilevel,.TRUE. )
           CALL relaxation_iterations_sym(ilevel,.FALSE.)
        END DO
        ! calculate residual on fine level
        CALL calculate_residual_sym(ilevel,res_PM_grid)
        ! restrict residual to coarse level
        CALL restrict_residual_sym(ilevel+1)
        !
     CASE(4)
        ! Gauss-Seidel iterations
        DO i_smth=1,Nsmth
           CALL relaxation_iterations_kmf(ilevel,.TRUE. )
           CALL relaxation_iterations_kmf(ilevel,.FALSE.)
        END DO
        ! calculate residual on fine level
        CALL calculate_residual_kmf(ilevel,res_PM_grid)
        ! restrict residual to coarse level
        CALL restrict_residual_kmf(ilevel+1)
        !
     CASE(5)
        ! Gauss-Seidel iterations
        DO i_smth=1,Nsmth
           CALL relaxation_iterations_csf(ilevel,.TRUE. )
           CALL relaxation_iterations_csf(ilevel,.FALSE.)
        END DO
        ! calculate residual on fine level
        CALL calculate_residual_csf(ilevel,res_PM_grid)
        ! restrict residual to coarse level
        CALL restrict_residual_csf(ilevel+1)
        !
     CASE DEFAULT
        ! error message
        WRITE(*,'(A)') 'Model currently unsupported. Please specify again.'
        STOP
        !
  END SELECT

  ! restrict scalar field to coarse level
  CALL restrict_scalar_field(ilevel+1)

  ! restrict density field to coarse level
  CALL restrict_density(ilevel+1)

  ! calculate physical RHS on coarse level
  SELECT CASE(MG_model)
     CASE(1)
        CALL calculate_physical_right_hand_side_fR (ilevel+1)
     CASE(2)
        CALL calculate_physical_right_hand_side_DGP(ilevel+1)
     CASE(3)
        CALL calculate_physical_right_hand_side_sym(ilevel+1)
     CASE(4)
        CALL calculate_physical_right_hand_side_kmf(ilevel+1)
     CASE(5)
        CALL calculate_physical_right_hand_side_csf(ilevel+1)
     CASE DEFAULT
        WRITE(*,'(A)') 'Model currently unsupported. Please specify again.'
        STOP
  END SELECT

  ! recursive call to next coarser multigrid level
  CALL W_cycle(ilevel+1)

  ! prolongation and correction to fine level solution
  CALL correct_solution(ilevel+1)

  ! pre-smoothing
  SELECT CASE(MG_model)
     CASE(1)
        ! Gauss-Seidel iterations
        DO i_smth=1,Nsmth
           CALL relaxation_iterations_fR(ilevel,.TRUE.)
           CALL relaxation_iterations_fR(ilevel,.FALSE.)
        END DO
        ! calculate residual on fine level
        CALL calculate_residual_fR(ilevel,res_PM_grid)
        ! restrict residual to coarse level
        CALL restrict_residual_fR(ilevel+1)
        !
     CASE(2)
        ! Gauss-Seidel iterations
        DO i_smth=1,Nsmth
           CALL relaxation_iterations_DGP(ilevel,.TRUE.)
           CALL relaxation_iterations_DGP(ilevel,.FALSE.)
        END DO
        ! calculate residual on fine level
        CALL calculate_residual_DGP(ilevel,res_PM_grid)
        ! restrict residual to coarse level
        CALL restrict_residual_DGP(ilevel+1)
        !
     CASE(3)
        ! Gauss-Seidel iterations
        DO i_smth=1,Nsmth
           CALL relaxation_iterations_sym(ilevel,.TRUE.)
           CALL relaxation_iterations_sym(ilevel,.FALSE.)
        END DO
        ! calculate residual on fine level
        CALL calculate_residual_sym(ilevel,res_PM_grid)
        ! restrict residual to coarse level
        CALL restrict_residual_sym(ilevel+1)
        !
     CASE(4)
        ! Gauss-Seidel iterations
        DO i_smth=1,Nsmth
           CALL relaxation_iterations_kmf(ilevel,.TRUE.)
           CALL relaxation_iterations_kmf(ilevel,.FALSE.)
        END DO
        ! calculate residual on fine level
        CALL calculate_residual_kmf(ilevel,res_PM_grid)
        ! restrict residual to coarse level
        CALL restrict_residual_kmf(ilevel+1)
        !
     CASE(5)
        ! Gauss-Seidel iterations
        DO i_smth=1,Nsmth
           CALL relaxation_iterations_csf(ilevel,.TRUE.)
           CALL relaxation_iterations_csf(ilevel,.FALSE.)
        END DO
        ! calculate residual on fine level
        CALL calculate_residual_csf(ilevel,res_PM_grid)
        ! restrict residual to coarse level
        CALL restrict_residual_csf(ilevel+1)
        !
     CASE DEFAULT
        ! message
        WRITE(*,'(A)') 'Model currently unsupported. Please specify again.'
        STOP
        ! 
  END SELECT

  ! restrict scalar field to coarse level
  CALL restrict_scalar_field(ilevel+1)

  ! restrict density field to coarse level
  CALL restrict_density(ilevel+1)

  ! calculate physical RHS on coarse level
  SELECT CASE(MG_model)
     CASE(1)
        CALL calculate_physical_right_hand_side_fR (ilevel+1)
     CASE(2)
        CALL calculate_physical_right_hand_side_DGP(ilevel+1)
     CASE(3)
        CALL calculate_physical_right_hand_side_sym(ilevel+1)
     CASE(4)
        CALL calculate_physical_right_hand_side_kmf(ilevel+1)
     CASE(5)
        CALL calculate_physical_right_hand_side_csf(ilevel+1)
     CASE DEFAULT
        WRITE(*,'(A)') 'Model currently unsupported. Please specify again.'
        STOP
  END SELECT

  ! recursive call to next coarser multigrid level
  CALL W_cycle(ilevel+1)

  ! prolongation and correction to fine level solution
  CALL correct_solution(ilevel+1)

  ! post-smoothing
  SELECT CASE(MG_model)
     CASE(1)
        ! Gauss-Seidel iterations
        DO i_smth=1,Nsmth
           CALL relaxation_iterations_fR (ilevel,.TRUE. )
           CALL relaxation_iterations_fR (ilevel,.FALSE.)
        END DO
        !
     CASE(2)
        ! Gauss-Seidel iterations
        DO i_smth=1,Nsmth
           CALL relaxation_iterations_DGP(ilevel,.TRUE. )
           CALL relaxation_iterations_DGP(ilevel,.FALSE.)
        END DO
        !
     CASE(3)
        ! Gauss-Seidel iterations
        DO i_smth=1,Nsmth
           CALL relaxation_iterations_sym(ilevel,.TRUE. )
           CALL relaxation_iterations_sym(ilevel,.FALSE.)
        END DO
        !
     CASE(4)
        ! Gauss-Seidel iterations
        DO i_smth=1,Nsmth
           CALL relaxation_iterations_kmf(ilevel,.TRUE. )
           CALL relaxation_iterations_kmf(ilevel,.FALSE.)
        END DO
        !
     CASE(5)
        ! Gauss-Seidel iterations
        DO i_smth=1,Nsmth
           CALL relaxation_iterations_csf(ilevel,.TRUE. )
           CALL relaxation_iterations_csf(ilevel,.FALSE.)
        END DO
        !
     CASE DEFAULT
        ! error message
        WRITE(*,'(A)') 'Model currently unsupported. Please specify again.'
        STOP
        !
  END SELECT

END SUBROUTINE W_cycle

!-------------------------------------------------------------
!
! Calculate the restricted density field on coarse level ilevel
! 
! Note: ilevel is the grid that is ilevel coarser than PM grid
! and the subroutine distinguishes between the cases where the 
! fine level is the PM grid (ilevel=1) or some multigrid level
!
!-------------------------------------------------------------
SUBROUTINE restrict_density(ilevel)
!-------------------------------------------------------------
use Tools

  integer :: ilevel
  integer :: ngrid_level 
  integer :: ioffset,joffset,koffset                     
  integer :: ioffset2,joffset2,koffset2 
  real*8  :: P
 
  IF(MG_test .EQ. 1) WRITE(*,'(A,I5)') 'Call restrict_density on level',levelmax-ilevel

  ! number of grid points on the coarse level
  ngrid_level = NGRID/2**ilevel
  
! ioffset  = 2**(levelmax-ilevel)
! joffset  = 0
! koffset  = 2**(levelmax-ilevel)*(2**ilevel-2)
! ioffset2 = 2**(levelmax-ilevel+1)
! joffset2 = 0
! koffset2 = 2**(levelmax-ilevel+1)*(2**(ilevel-1)-2)
  ioffset  = NGRID/2**(ilevel  )                                                ! Baojiu 2025 12 17
  joffset  = 0                                                                  ! Baojiu 2025 12 17
  koffset  = NGRID/2**(ilevel  )*(2**(ilevel  )-2)                              ! Baojiu 2025 12 17
  ioffset2 = NGRID/2**(ilevel-1)                                                ! Baojiu 2025 12 17
  joffset2 = 0                                                                  ! Baojiu 2025 12 17
  koffset2 = NGRID/2**(ilevel-1)*(2**(ilevel-1)-2)                              ! Baojiu 2025 12 17

  IF(ilevel.EQ.1) THEN
  ! case where fine level is the PM grid, where the density field
  ! on fine level is stored in array FI
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3) &
!$OMP PRIVATE (P) 
     DO M3=1,ngrid_level
        DO M2=1,ngrid_level
           DO M1=1,ngrid_level
              ! accumulate density field values from 8 son cells
              P =   FI(2*M1-1,2*M2-1,2*M3-1)
              P = P+FI(2*M1  ,2*M2-1,2*M3-1)
              P = P+FI(2*M1-1,2*M2  ,2*M3-1)
              P = P+FI(2*M1  ,2*M2  ,2*M3-1)
              P = P+FI(2*M1-1,2*M2-1,2*M3  )
              P = P+FI(2*M1  ,2*M2-1,2*M3  )
              P = P+FI(2*M1-1,2*M2  ,2*M3  )
              P = P+FI(2*M1  ,2*M2  ,2*M3  )
              !
              FI3(ioffset+M1,M2,M3) = P/8.0d0  
              !
           END DO
        END DO
     END DO
  ELSE
  ! case where fine level is a multigrid level, where the density
  ! field on fine level is stored in a particular section of FI3
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3) &
!$OMP PRIVATE (P) 
     DO M3=1,ngrid_level
        DO M2=1,ngrid_level
           DO M1=1,ngrid_level
              ! accumulate density field values from 8 son cells
              P =   FI3(ioffset2+2*M1-1,2*M2-1,koffset2+2*M3-1)
              P = P+FI3(ioffset2+2*M1  ,2*M2-1,koffset2+2*M3-1)
              P = P+FI3(ioffset2+2*M1-1,2*M2  ,koffset2+2*M3-1)
              P = P+FI3(ioffset2+2*M1  ,2*M2  ,koffset2+2*M3-1)
              P = P+FI3(ioffset2+2*M1-1,2*M2-1,koffset2+2*M3  )
              P = P+FI3(ioffset2+2*M1  ,2*M2-1,koffset2+2*M3  )
              P = P+FI3(ioffset2+2*M1-1,2*M2  ,koffset2+2*M3  )
              P = P+FI3(ioffset2+2*M1  ,2*M2  ,koffset2+2*M3  )
              !
              FI3(ioffset+M1,M2,koffset+M3) = P/8.0d0   
              !
           END DO
        END DO
     END DO

  ENDIF
 
  CALL TimingMain(3,1)

END SUBROUTINE restrict_density

!-------------------------------------------------------------
!
! Calculate the restricted scalar field on coarse level ilevel
! 
! Note: ilevel is the grid that is ilevel coarser than PM grid
! and the subroutine distinguishes between the cases where the 
! fine level is the PM grid (ilevel=1) or some multigrid level
!
!-------------------------------------------------------------
SUBROUTINE restrict_scalar_field(ilevel)
!-------------------------------------------------------------
use Tools

  integer :: ilevel
  integer :: ngrid_level 
  integer :: ioffset,joffset,koffset       
  integer :: ioffset2,joffset2,koffset2 
  integer :: M1,M2,M3
  real*8  :: P             

  IF(MG_test .EQ. 1) WRITE(*,'(A,I5)') 'Call restrict_scalar_field on level',levelmax-ilevel

 ! number of grid points on the coarse level
  ngrid_level = NGRID/2**ilevel
  
! ioffset  = 2**(levelmax-ilevel)
! joffset  = 2**(levelmax-ilevel)
! koffset  = 2**(levelmax-ilevel)*(2**ilevel-2)
! ioffset2 = 2**(levelmax-ilevel+1)
! joffset2 = 2**(levelmax-ilevel+1)
! koffset2 = 2**(levelmax-ilevel+1)*(2**(ilevel-1)-2)
  ioffset  = NGRID/2**(ilevel  )                                                ! Baojiu 2025 12 17
  joffset  = NGRID/2**(ilevel  )                                                ! Baojiu 2025 12 17
  koffset  = NGRID/2**(ilevel  )*(2**(ilevel  )-2)                              ! Baojiu 2025 12 17
  ioffset2 = NGRID/2**(ilevel-1)                                                ! Baojiu 2025 12 17
  joffset2 = NGRID/2**(ilevel-1)                                                ! Baojiu 2025 12 17
  koffset2 = NGRID/2**(ilevel-1)*(2**(ilevel-1)-2)                              ! Baojiu 2025 12 17

  IF(ilevel.EQ.1) THEN
  ! case where fine level is the PM grid, where the scalar field
  ! on fine level is stored in array FI2
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3) &
!$OMP PRIVATE (P) 
     DO M3=1,ngrid_level
        DO M2=1,ngrid_level
           DO M1=1,ngrid_level
              ! accumulate scalar field values from 8 son cells
              P =   FI2(2*M1-1,2*M2-1,2*M3-1)
              P = P+FI2(2*M1  ,2*M2-1,2*M3-1)
              P = P+FI2(2*M1-1,2*M2  ,2*M3-1)
              P = P+FI2(2*M1  ,2*M2  ,2*M3-1)
              P = P+FI2(2*M1-1,2*M2-1,2*M3  )
              P = P+FI2(2*M1  ,2*M2-1,2*M3  )
              P = P+FI2(2*M1-1,2*M2  ,2*M3  )
              P = P+FI2(2*M1  ,2*M2  ,2*M3  )
              !
              FI3(        M1,joffset+M2,M3) = P/8.0d0
              !
              FI3(ioffset+M1,joffset+M2,M3) = P/8.0d0
              !
           END DO
        END DO
     END DO
  ELSE
  ! case where fine level is a multigrid level, where the scalar
  ! field on fine level is stored in a particular section of FI3
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3) &
!$OMP PRIVATE (P) 
     DO M3=1,ngrid_level
        DO M2=1,ngrid_level
           DO M1=1,ngrid_level
              ! accumulate scalar field values from 8 son cells
              P =   FI3(ioffset2+2*M1-1,joffset2+2*M2-1,koffset2+2*M3-1)
              P = P+FI3(ioffset2+2*M1  ,joffset2+2*M2-1,koffset2+2*M3-1)
              P = P+FI3(ioffset2+2*M1-1,joffset2+2*M2  ,koffset2+2*M3-1)
              P = P+FI3(ioffset2+2*M1  ,joffset2+2*M2  ,koffset2+2*M3-1)
              P = P+FI3(ioffset2+2*M1-1,joffset2+2*M2-1,koffset2+2*M3  )
              P = P+FI3(ioffset2+2*M1  ,joffset2+2*M2-1,koffset2+2*M3  )
              P = P+FI3(ioffset2+2*M1-1,joffset2+2*M2  ,koffset2+2*M3  )
              P = P+FI3(ioffset2+2*M1  ,joffset2+2*M2  ,koffset2+2*M3  )
              !
              FI3(        M1,joffset+M2,koffset+M3) = P/8.0d0      
              !
              FI3(ioffset+M1,joffset+M2,koffset+M3) = P/8.0d0
              !
           END DO
        END DO
     END DO
  ENDIF
    
  CALL TimingMain(3,1)

END SUBROUTINE restrict_scalar_field

!-------------------------------------------------------------
!
! Correct the fine-level solution using coarse-level solution.    
!
! Note: "coarse level" here means the multigrid level that has
! NGRID/2^ilevel grid points, which is ilevel levels below (or
! coarser than) the PM grid, while "fine level" here means the 
! multigrid or PM grid one level finer than the coarse level.
!
!-------------------------------------------------------------
SUBROUTINE correct_solution(ilevel)
!-------------------------------------------------------------
! This subroutine loops over all coarse level cells, and makes
! correction to the 8 son cells of each of these coarse cells.
! For a given fine cell, this prolongation operation will make 
! use of 8 of the 27 neighbouring coarse cells of its "parent"
! cell (including the parent cell itself). We shall call these 
! 27 coarse cells a "block" with the parent cell at the centre 
! of the block. The 8 coarse cells that make contribution to a
! fine cell have either a common face, a common edge or common
! vertex with the fine cell, or contains the fine cell itself,
! and their contributions have a weight of 9/64, 3/64, 1/64 or
! 27/64, respectively. 
!-------------------------------------------------------------
  use Tools
  use ExtradofBackgroundData

  integer :: ilevel
  integer :: ngrid_level 
  integer :: ioffset,joffset,koffset,ioffset2,joffset2,koffset2
  integer :: i1,i2,i3       
  integer :: ind,sind,pind 
  real*8  :: a,b,c,d,P
  real*8, dimension(1:8)     :: bbb
  integer,dimension(1:8,1:8) :: ccc
  integer,dimension(1:27)    :: pid1,pid2,pid3

  real*8  :: tmpt,ctilde,ctilde2
  integer :: ist
  real*8  :: phibar

  IF(MG_test .EQ. 1) WRITE(*,'(A,I5)') 'Call correct_solution for level',levelmax-ilevel
  !
  ctilde  = 2.99792458D3*DBLE(NGRID)/Box
  ctilde2 = ctilde**2
  !
  IF(MG_model.EQ.5) THEN
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
  a = 1.0d0/64.0d0
  b = 3.0d0 *a
  c = 9.0d0 *a
  d = 27.0d0*a
  !
  ! array containing weights of the 8 coarse cells
  bbb(:)   = (/a ,b ,b ,c ,b ,c ,c ,d/)
  !
  ! each array contains the indices of the 8 coarse cells that make
  ! contributions to one of the 8 fine (or son) cells of the parent 
  ! cell. ccc(:,1) for the first son cell, and so on. 1, ..., 27 is
  ! the position or index of the coarse cell in the block, with the 
  ! index of the parent cell being 14 (denotes centre of the block)
  ccc(:,1) = (/1 ,2 ,4 ,5 ,10,11,13,14/)
  ccc(:,2) = (/3 ,2 ,6 ,5 ,12,11,15,14/)
  ccc(:,3) = (/7 ,8 ,4 ,5 ,16,17,13,14/)
  ccc(:,4) = (/9 ,8 ,6 ,5 ,18,17,15,14/)
  ccc(:,5) = (/19,20,22,23,10,11,13,14/)
  ccc(:,6) = (/21,20,24,23,12,11,15,14/)
  ccc(:,7) = (/25,26,22,23,16,17,13,14/)
  ccc(:,8) = (/27,26,24,23,18,17,15,14/)

  ! number of grid points at the coarse level
  ngrid_level = NGRID/2**ilevel
  
! ioffset  = 2**(levelmax-ilevel)
! joffset  = 2**(levelmax-ilevel)
! koffset  = 2**(levelmax-ilevel)*(2**ilevel-2)
! ioffset2 = 2**(levelmax-ilevel+1)
! joffset2 = 2**(levelmax-ilevel+1)
! koffset2 = 2**(levelmax-ilevel+1)*(2**(ilevel-1)-2)
  ioffset  = NGRID/2**(ilevel  )                                                ! Baojiu 2025 12 17
  joffset  = NGRID/2**(ilevel  )                                                ! Baojiu 2025 12 17
  koffset  = NGRID/2**(ilevel  )*(2**(ilevel  )-2)                              ! Baojiu 2025 12 17
  ioffset2 = NGRID/2**(ilevel-1)                                                ! Baojiu 2025 12 17
  joffset2 = NGRID/2**(ilevel-1)                                                ! Baojiu 2025 12 17
  koffset2 = NGRID/2**(ilevel-1)*(2**(ilevel-1)-2)                              ! Baojiu 2025 12 17

! IF(ilevel.EQ.1) THEN
!    OPEN(UNIT=27, FILE='test1.txt', FORM='FORMATTED', STATUS='REPLACE')
!       DO M3=1,ngrid_level*2
!       !  WRITE(27,'(F20.9,F20.15)') (DBLE(M3)-0.5d0)/DBLE(NGRID)*2.0d0*DACOS(-1.0d0),FI3(ioffset2+ngrid_level,joffset2+ngrid_level,koffset2+M3)
!          WRITE(27,'(F20.9,F20.15)') (DBLE(M3)-0.5d0)/DBLE(NGRID)*2.0d0*DACOS(-1.0d0),FI2(ngrid_level,ngrid_level,M3)
!       ENDDO
!    CLOSE(27)
! ENDIF

  IF(ilevel.EQ.1) THEN
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3,pind,sind,ind,i1,i2,i3) &
!$OMP PRIVATE (pid1,pid2,pid3) &
!$OMP PRIVATE (P) 
     DO M3=1,ngrid_level
        DO M2=1,ngrid_level
           DO M1=1,ngrid_level
              ! get 3^3 neighbours of a parent cell (the parent cell is at centre);
              ! each coarse cell is given a 3D local index from (0,0,0) to (2,2,2),
              ! and a 3D global index (pid1,pid2,pid3) too; periodic BC is applied.
              DO pind=1,27
                 ! local index of parent cells
                 i3 = (pind-1)/9
                 i2 = (pind-1-i3*9)/3
                 i1 =  pind-1-i3*9-i2*3
                 ! global index of parent cells
                 pid1(pind) = M1+i1-1; IF(pid1(pind)<1) pid1(pind)=ngrid_level; IF(pid1(pind)>ngrid_level) pid1(pind)=1
                 pid2(pind) = M2+i2-1; IF(pid2(pind)<1) pid2(pind)=ngrid_level; IF(pid2(pind)>ngrid_level) pid2(pind)=1
                 pid3(pind) = M3+i3-1; IF(pid3(pind)<1) pid3(pind)=ngrid_level; IF(pid3(pind)>ngrid_level) pid3(pind)=1
              END DO
              ! loop over 2^3 son cells of the parent cell at the centre of block;
              ! each son cell is given a 3D local index from (0,0,0) to (1,1,1) and
              ! the corresponding 3D global index is (2*M1+i1, 2*M2+i2, 2*M3+i3) 
              DO sind=1,8
                 ! local index of son cells of parent cell
                 i3 = (sind-1)/4
                 i2 = (sind-1-i3*4)/2
                 i1 = (sind-1-i3*4-i2*2)
                 P = 0.0d0
                 ! loop over 8 (of the 27) neighbouring coarse cells
                 DO ind=1,8
                    P = P+bbb(ind)*FI3(ioffset+pid1(ccc(ind,sind)),joffset+pid2(ccc(ind,sind)),pid3(ccc(ind,sind)))
                    P = P-bbb(ind)*FI3(        pid1(ccc(ind,sind)),joffset+pid2(ccc(ind,sind)),pid3(ccc(ind,sind)))
                 END DO
                 !
                 tmpt = FI2(2*M1+i1-1,2*M2+i2-1,2*M3+i3-1)+P
                 ! correct solutions
                 ! MG_model=1: f(R) gravity, scalar field u must be positive
                 ! MG_model=3: symmetron model, scalar field u must be posistive
                 IF(MG_model.EQ.1.OR.MG_model.EQ.3) THEN
                    IF(tmpt.GT.0.0D0) FI2(2*M1+i1-1,2*M2+i2-1,2*M3+i3-1) = tmpt
                 ! MG_model=5: coupled scalar field, total scalar field (background+perturbation) must be positive
                 ELSE IF(MG_model.EQ.5) THEN
                    IF(phibar+tmpt/ctilde2.GT.0.0D0) FI2(2*M1+i1-1,2*M2+i2-1,2*M3+i3-1) = tmpt
                 ! other models
                 ELSE
                    FI2(2*M1+i1-1,2*M2+i2-1,2*M3+i3-1) = tmpt
                 ENDIF
              END DO
           END DO
        END DO
     END DO
  ELSE
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3,pind,sind,ind,i1,i2,i3) &
!$OMP PRIVATE (pid1,pid2,pid3) &
!$OMP PRIVATE (P) 
     DO M3=1,ngrid_level
        DO M2=1,ngrid_level
           DO M1=1,ngrid_level
              ! get 3^3 neighbours of a parent cell (the parent cell is at centre);
              ! each coarse cell is given a 3D local index from (0,0,0) to (2,2,2),
              ! and a 3D global index (pid1,pid2,pid3) too; periodic BC is applied.
              DO pind=1,27
                 ! local index of parent cells
                 i3 = (pind-1)/9
                 i2 = (pind-1-i3*9)/3
                 i1 =  pind-1-i3*9-i2*3
                 ! global index of parent cells
                 pid1(pind) = M1+i1-1; IF(pid1(pind)<1) pid1(pind)=ngrid_level; IF(pid1(pind)>ngrid_level) pid1(pind)=1
                 pid2(pind) = M2+i2-1; IF(pid2(pind)<1) pid2(pind)=ngrid_level; IF(pid2(pind)>ngrid_level) pid2(pind)=1
                 pid3(pind) = M3+i3-1; IF(pid3(pind)<1) pid3(pind)=ngrid_level; IF(pid3(pind)>ngrid_level) pid3(pind)=1
              END DO
              ! loop over 2^3 son cells of the parent cell at the centre of block;
              ! each son cell is given a 3D local index from (0,0,0) to (1,1,1) and
              ! the corresponding 3D global index is (2*M1+i1, 2*M2+i2, 2*M3+i3) 
              DO sind=1,8
                 ! local index of son cells of parent cell
                 i3 = (sind-1)/4
                 i2 = (sind-1-i3*4)/2
                 i1 = (sind-1-i3*4-i2*2)
                 P = 0.0d0
                 ! loop over 8 (of the 27) neighbouring coarse cells
                 DO ind=1,8                                                     ! u^H - R(u^h)
                    P = P+bbb(ind)*FI3(ioffset+pid1(ccc(ind,sind)),joffset+pid2(ccc(ind,sind)),koffset+pid3(ccc(ind,sind))) ! FI3 Section 4
                    P = P-bbb(ind)*FI3(        pid1(ccc(ind,sind)),joffset+pid2(ccc(ind,sind)),koffset+pid3(ccc(ind,sind))) ! FI3 Section 3
                 END DO
                 !
                 tmpt = FI3(ioffset2+2*M1+i1-1,joffset2+2*M2+i2-1,koffset2+2*M3+i3-1)+P
                 ! correct solutions
                 ! MG_model=1: f(R) gravity, scalar field u must be positive
                 ! MG_model=3: symmetron model, scalar field u must be posistive
                 IF(MG_model.EQ.1.OR.MG_model.EQ.3) THEN
                    IF(tmpt.GT.0.0D0) FI3(ioffset2+2*M1+i1-1,joffset2+2*M2+i2-1,koffset2+2*M3+i3-1) = tmpt
                 ! MG_model=5: coupled scalar field, total scalar field (background+perturbation) must be positive
                 ELSE IF(MG_model.EQ.5) THEN
                    IF(phibar+tmpt/ctilde2.GT.0.0D0) FI3(ioffset2+2*M1+i1-1,joffset2+2*M2+i2-1,koffset2+2*M3+i3-1) = tmpt
                 ! other models
                 ELSE
                    FI3(ioffset2+2*M1+i1-1,joffset2+2*M2+i2-1,koffset2+2*M3+i3-1) = tmpt
                 ENDIF
              END DO
           END DO
        END DO
     END DO
  ENDIF

! IF(ilevel.EQ.1) THEN
!    OPEN(UNIT=27, FILE='test2.txt', FORM='FORMATTED', STATUS='REPLACE')
!       DO M3=1,ngrid_level*2
!       !  WRITE(27,'(F20.9,F20.15)') (DBLE(M3)-0.5d0)/DBLE(NGRID)*2.0d0*DACOS(-1.0d0),FI3(ioffset2+ngrid_level,joffset2+ngrid_level,koffset2+M3)
!          WRITE(27,'(F20.9,F20.15)') (DBLE(M3)-0.5d0)/DBLE(NGRID)*2.0d0*DACOS(-1.0d0),FI2(ngrid_level,ngrid_level,M3)
!       ENDDO
!    CLOSE(27)
!    STOP
! ENDIF
    
  CALL TimingMain(3,1)

END SUBROUTINE correct_solution


SUBROUTINE WriteScalarFieldToFile
  USE Tools

  integer :: M1
  real*8  :: KK,AAA
  real*8  :: R_bg,R0_bg,fR_bg
  real*8  :: u0,uinf,ctilde,rho_in,rho_ou,m_in,m_ou,RR

  OPEN(UNIT=27, FILE='N5_L64Ng256_Vcycle_spherical_d0.5_r0.1.txt', FORM='FORMATTED', STATUS='REPLACE')

  ! sine test
  AAA = 0.1D0
  KK  = 4.0D0
  IF(MG_model.EQ.1) THEN
     R_bg  = 3.0D0*(Om/AEXPN**3+4.0D0*OmL)
     R0_bg = 3.0D0*(Om         +4.0D0*OmL)
     fR_bg = fR0*(R0_bg/R_bg)**(fr_n+1)
     DO M1=1,NGRID
!       WRITE(27,'(F20.9,F20.15,F20.15,F20.15)') (DBLE(M1)-0.5d0)/DBLE(NGRID)*2.0d0*DACOS(-1.0d0),fR0*FI2(64,64,M1)**2,FI(64,64,M1) &
!          &  ,fR0*(1.0D0+0.1D0*DSIN((DBLE(M1)-0.5d0)/DBLE(NGRID)*2.0d0*DACOS(-1.0d0)))
        WRITE(27,'(F20.9,F20.15,F20.15,F20.15)') (DBLE(M1)-0.5d0)/DBLE(NGRID)*2.0d0*DACOS(-1.0d0),FI2(64,64,M1)**(1+fr_n)/DABS(fR_bg),FI(64,64,M1) &
           &  ,(1.0D0+AAA*DSIN((DBLE(M1)-0.5d0)/DBLE(NGRID)*2.0d0*KK*DACOS(-1.0d0)))
     END DO
  ENDIF
  IF(MG_model.EQ.2) THEN
     DO M1=1,NGRID
!        WRITE(27,'(F20.9,F20.15,F20.15,F20.15)') (DBLE(M1)-0.5d0)/DBLE(NGRID)*2.0d0*DACOS(-1.0d0),FI2(64,64,M1),FI(64,64,M1) &
!           &  ,(0.1D0*DSIN((DBLE(M1)-0.5d0)/DBLE(NGRID)*2.0d0*DACOS(-1.0d0)))
        WRITE(27,'(F20.9,F20.15,F20.15)') DBLE(M1)/DBLE(NGRID),FI2(NGRID/2,NGRID/2,M1),FI(NGRID/2,NGRID/2,M1)
     END DO
  ENDIF
  IF(MG_model.EQ.3) THEN
     DO M1=1,NGRID
        WRITE(27,'(F20.9,F20.15,F20.15,F20.15)') (DBLE(M1)-0.5d0)/DBLE(NGRID)*2.0d0*DACOS(-1.0d0),FI2(NGRID/2,NGRID/2,M1),FI(NGRID/2,NGRID/2,M1), &
           &  1.0D0+AAA*DSIN((DBLE(M1)-0.5d0)/DBLE(NGRID)*2.0d0*KK*DACOS(-1.0d0))
     END DO
  ENDIF
  IF(MG_model.EQ.4) THEN
     DO M1=1,NGRID
!       WRITE(27,'(F20.9,F20.15,F20.15,F20.15)') (DBLE(M1)-0.5d0)/DBLE(NGRID)*2.0d0*DACOS(-1.0d0),FI2(NGRID/2,NGRID/2,M1),FI(NGRID/2,NGRID/2,M1), &
!          &  AAA*DSIN((DBLE(M1)-0.5d0)/DBLE(NGRID)*2.0d0*KK*DACOS(-1.0d0))
        WRITE(27,'(F23.9,F23.15,F23.15)') DBLE(M1)-DBLE(NGRID/2),FI2(NGRID/2,NGRID/2,M1),FI(NGRID/2,NGRID/2,M1)
     END DO
  ENDIF
  IF(MG_model.EQ.5) THEN
     ctilde = 2.99792458D3*DBLE(NGRID)/Box
     DO M1=1,NGRID
!       WRITE(27,'(F23.9,F23.15,F23.15)') 2.0D0*DACOS(-1.0D0)*(DBLE(M1)-0.5D0)/DBLE(NGRID),FI2(NGRID/2,NGRID/2,M1),AAA*DSIN(KK*2.0D0*DACOS(-1.0D0)*(DBLE(M1)-0.5D0)/DBLE(NGRID))
        WRITE(27,'(F23.9,F23.15,F23.15)') DBLE(M1)-DBLE(NGRID/2),FI2(NGRID/2,NGRID/2,M1)/ctilde**2,FI(NGRID/2,NGRID/2,M1)
     END DO
  ENDIF
  ! Gaussian test
! IF(MG_model.EQ.2) THEN
!    DO M1=1,NGRID
!       WRITE(27,'(F20.9,F20.15,F20.15,F20.15)') (DBLE(M1)-0.5d0),FI2(64,64,M1),FI(64,64,M1) &
!          & ,1.D0*DEXP(-((DBLE(M1)-0.5d0)-DBLE(NGRID/2))**2/  16.0D0)
!    END DO
! ENDIF
! IF(MG_model.EQ.3) THEN
!    ctilde = 2.99792458D3*DBLE(NGRID)/Box
!    rho_in = 5000.0D0
!    rho_ou = 1.0D0
!    RR   = DBLE(NGRID)*0.1D0
!    uinf = DSQRT(1.0D0-sym_astar**3)
!    m_in = DSQRT(0.5D0*(AEXPN/(sym_xi*ctilde))**2*( (sym_astar/AEXPN)**3*rho_in-1.0D0))
!    m_ou = DSQRT(1.0D0*(AEXPN/(sym_xi*ctilde))**2*(-(sym_astar/AEXPN)**3*rho_ou+1.0D0))
!    u0   = uinf*(1.0D0+m_ou*RR)/(0.5D0*(DEXP(m_in*RR)+DEXP(-m_ou*RR))+0.5D0*(m_ou/m_in)*(DEXP(m_in*RR)-DEXP(-m_in*RR)))
!    DO M1=1,NGRID
!       WRITE(27,'(F22.9,F22.15,F22.15,F22.15)') DSQRT(1.0D0)*(DBLE(M1)-DBLE(NGRID/2)),FI2(M1,NGRID/2,NGRID/2),FI(M1,NGRID/2,NGRID/2), &
!          &  m_in/m_ou
!    END DO
! ENDIF

  CLOSE(27)

END SUBROUTINE WriteScalarFieldToFile



SUBROUTINE SHOWFI2
   USE Tools 

   WRITE(*, *) 'Here are the values of some elements of FI2...'
   WRITE(*, '(F10.6)') FI2(64,64,1  ) 
   WRITE(*, '(F10.6)') FI2(64,64,10 )
   WRITE(*, '(F10.6)') FI2(64,64,20 )
   WRITE(*, '(F10.6)') FI2(64,64,30 )
   WRITE(*, '(F10.6)') FI2(64,64,40 )
   WRITE(*, '(F10.6)') FI2(64,64,50 )
   WRITE(*, '(F10.6)') FI2(64,64,60 )
   WRITE(*, '(F10.6)') FI2(64,64,70 )
   WRITE(*, '(F10.6)') FI2(64,64,80 )
   WRITE(*, '(F10.6)') FI2(64,64,90 )
   WRITE(*, '(F10.6)') FI2(64,64,100)
   WRITE(*, '(F10.6)') FI2(64,64,110)
   WRITE(*, '(F10.6)') FI2(64,64,120)
   WRITE(*, '(F10.6)') FI2(64,64,128)



END SUBROUTINE SHOWFI2
