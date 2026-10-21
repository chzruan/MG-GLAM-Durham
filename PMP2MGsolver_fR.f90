!=================================================
!=================================================
!
!  Scalar field
!
!=================================================
!=================================================

! Baojiu: changed to double precision throughout the file
!
! Storage (2026-10-21, docs/memory-audit-20261021.md, cut 2): on the PM grid
! the density is FI and the scalar field u is FI2; every coarse multigrid
! level ilevel >= 1 has its own arrays MGL(ilevel)%{res,rho,ucopy,u,rhs,lres}
! (module Tools) in place of the six sections of the former full-grid FI3.
! The arithmetic of every loop is unchanged, only the addresses differ; the
! coarse-level loops are internal procedures with explicit-shape dummies so
! that the compiler sees plain contiguous arrays. The fine-level residual is
! no longer stored (nothing read it).

!-------------------------------------------------------------
!
! red-black Gauss-Seidel relaxation iterations on level ilevel
!
! Note: ilevel runs from 0, which means relaxation on PM grid.
! This subroutine distinguishes between PM and multgrid grids,
! since the scalar field and density fields etc. are stored in
! different places in these two cases.

SUBROUTINE relaxation_iterations_fR(ilevel, redstep)
    use Tools

    implicit none

    integer :: ilevel
    logical :: redstep
    real*8  :: Q, Q2onethird, Q2, Q4, CC, R0_bg, R_bg, fR_bg, ctilde, dx, dx2
    real*8  :: s_fR_bg, a_fR_bg
    real*8  :: fct1, fct2
    real*8  :: ZERO, ONEOVERTHREE, TWOOVERTHREE, TWOPIOVERTHREE
    integer :: nstart, nfinal, ngrid_level

    real*8  :: Delta10, Delta11, Delta20, Delta21, tmp1, tmp2, P, S, L, sqrt_Delta0
    integer :: M1, M2, M3, M1u, M1l, M2u, M2l, M3u, M3l
    integer :: M1start, parity_off

    integer :: nid, myid

    CALL TimingMain(8, -1)

    IF (MG_test .EQ. 1) WRITE (*, '(A,I5,F7.4)') 'Relaxation iterations on level', levelmax - ilevel, AEXPN

    ! constants
    ZERO = 1.0D-20
    ONEOVERTHREE = 1.0D0/3.0D0
    TWOOVERTHREE = 2.0D0/3.0D0
    TWOPIOVERTHREE = 2.0D0/3.0D0*DACOS(-1.0D0)

    ! number of grid points on the coarse level
    ngrid_level = NGRID/(2**ilevel)

    ! simulation parameters
    ctilde = 2.99792458D3*DBLE(NGRID)/Box
    dx = DBLE(2**ilevel)
    dx2 = dx*dx

    ! physical and numerical quantities
    R_bg = 3.0d0*(Om/AEXPN**3 + 4.0d0*OmL)
    R0_bg = 3.0d0*(Om + 4.0d0*OmL)
    fR_bg = fR0*(R0_bg/R_bg)**(fr_n + 1)
    fct2 = dx2/(ctilde**2*DSIGN(1.0D0, fR_bg))/6.0d0
    fct1 = -Om/AEXPN*fct2
    Q2 = R_bg*AEXPN**2*fct2/3.0d0
    Q = Q2*(DABS(fR_bg))**(1.0D0/(DBLE(fr_n) + 1.0D0))

    s_fR_bg = DSIGN(1.0D0, fR_bg)

    ! intermediate quantities
    Q2onethird = (-Q)**ONEOVERTHREE
    Q4 = 4.0D0*Q
    Delta11 = 27.0d0*Q
    Delta20 = 12.0d0*Q

    ! Red-black parity offset: 0 for red pass, 1 for black pass.
    ! Combined with MOD(M2+M3+parity_off, 2) below, picks the M1 start index
    ! so the stride-2 inner loop hits the same cells the CYCLE version did.
    parity_off = MERGE(0, 1, redstep)

    IF (ilevel .GT. 0) THEN
        ! coarse level: working field u, restricted density rho, physical RHS rhs
        CALL relax_coarse(MGL(ilevel)%u, MGL(ilevel)%rho, MGL(ilevel)%rhs, ngrid_level)
    ELSE IF (fr_n .EQ. 0) THEN
!$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(STATIC) DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3,M1start,M1u,M1l,M2u,M2l,M3u,M3l) &
!$OMP PRIVATE (P)
        DO M3 = 1, ngrid_level
            DO M2 = 1, ngrid_level
                M1start = 1 + MOD(M2 + M3 + parity_off, 2)
            !DIR$ IVDEP
                DO M1 = M1start, ngrid_level, 2
                    M1u = M1 + 1; IF (M1u > ngrid_level) M1u = 1
                    M1l = M1 - 1; IF (M1l < 1) M1l = ngrid_level
                    M2u = M2 + 1; IF (M2u > ngrid_level) M2u = 1
                    M2l = M2 - 1; IF (M2l < 1) M2l = ngrid_level
                    M3u = M3 + 1; IF (M3u > ngrid_level) M3u = 1
                    M3l = M3 - 1; IF (M3l < 1) M3l = ngrid_level
                    P = fct1*FI(M1, M2, M3) - Q2 - &
                        & (FI2(M1u, M2, M3) + &
                        &  FI2(M1l, M2, M3) + &
                        &  FI2(M1, M2u, M3) + &
                        &  FI2(M1, M2l, M3) + &
                        &  FI2(M1, M2, M3u) + &
                        &  FI2(M1, M2, M3l))/6.0D0
                    FI2(M1, M2, M3) = 0.5D0*(-P + DSQRT(P*P - Q4))
                END DO
            END DO
        END DO
    ELSE IF (fr_n .EQ. 1) THEN
!$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(STATIC) DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3,M1start,M1u,M1l,M2u,M2l,M3u,M3l) &
!$OMP PRIVATE (Delta10,tmp1,tmp2,sqrt_Delta0,P,S,CC)
        DO M3 = 1, ngrid_level
            DO M2 = 1, ngrid_level
                M1start = 1 + MOD(M2 + M3 + parity_off, 2)
            !DIR$ IVDEP
                DO M1 = M1start, ngrid_level, 2
                    M1u = M1 + 1; IF (M1u > ngrid_level) M1u = 1
                    M1l = M1 - 1; IF (M1l < 1) M1l = ngrid_level
                    M2u = M2 + 1; IF (M2u > ngrid_level) M2u = 1
                    M2l = M2 - 1; IF (M2l < 1) M2l = ngrid_level
                    M3u = M3 + 1; IF (M3u > ngrid_level) M3u = 1
                    M3l = M3 - 1; IF (M3l < 1) M3l = ngrid_level
                    P = fct1*FI(M1, M2, M3) - Q2 - &
                        & (FI2(M1u, M2, M3)*FI2(M1u, M2, M3) + &
                        &  FI2(M1l, M2, M3)*FI2(M1l, M2, M3) + &
                        &  FI2(M1, M2u, M3)*FI2(M1, M2u, M3) + &
                        &  FI2(M1, M2l, M3)*FI2(M1, M2l, M3) + &
                        &  FI2(M1, M2, M3u)*FI2(M1, M2, M3u) + &
                        &  FI2(M1, M2, M3l)*FI2(M1, M2, M3l))/6.0D0
                    Delta10 = -3.0D0*P
                    CC = Delta11*Delta11 - 4.0d0*Delta10*Delta10*Delta10
                    IF (CC .GE. 0.0D0) THEN
                        tmp1 = 0.5d0*(DSQRT(CC) + Delta11)
                        tmp2 = tmp1 - Delta11
                        IF (DABS(P) .LT. ZERO) THEN
                            FI2(M1, M2, M3) = Q2onethird
                        ELSE
                            FI2(M1, M2, M3) = -(DSIGN(DABS(tmp1)**ONEOVERTHREE, tmp1) - DSIGN(DABS(tmp2)**ONEOVERTHREE, tmp2))*ONEOVERTHREE
                        END IF
                    ELSE
                        sqrt_Delta0 = DSQRT(Delta10)
                        S = DACOS(0.5d0*Delta11/(sqrt_Delta0*sqrt_Delta0*sqrt_Delta0))
                        FI2(M1, M2, M3) = -TWOOVERTHREE*sqrt_Delta0*DCOS(S*ONEOVERTHREE + TWOPIOVERTHREE)
                    END IF
                END DO
            END DO
        END DO
    ELSE IF (fr_n .EQ. 2) THEN
!$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(STATIC) DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3,M1start,M1u,M1l,M2u,M2l,M3u,M3l) &
!$OMP PRIVATE (Delta21,tmp1,tmp2,P,S)
        DO M3 = 1, ngrid_level
            DO M2 = 1, ngrid_level
                M1start = 1 + MOD(M2 + M3 + parity_off, 2)
            !DIR$ IVDEP
                DO M1 = M1start, ngrid_level, 2
                    M1u = M1 + 1; IF (M1u > ngrid_level) M1u = 1
                    M1l = M1 - 1; IF (M1l < 1) M1l = ngrid_level
                    M2u = M2 + 1; IF (M2u > ngrid_level) M2u = 1
                    M2l = M2 - 1; IF (M2l < 1) M2l = ngrid_level
                    M3u = M3 + 1; IF (M3u > ngrid_level) M3u = 1
                    M3l = M3 - 1; IF (M3l < 1) M3l = ngrid_level
                    P = fct1*FI(M1, M2, M3) - Q2 - &
                        & (FI2(M1u, M2, M3)*FI2(M1u, M2, M3)*FI2(M1u, M2, M3) + &
                        &  FI2(M1l, M2, M3)*FI2(M1l, M2, M3)*FI2(M1l, M2, M3) + &
                        &  FI2(M1, M2u, M3)*FI2(M1, M2u, M3)*FI2(M1, M2u, M3) + &
                        &  FI2(M1, M2l, M3)*FI2(M1, M2l, M3)*FI2(M1, M2l, M3) + &
                        &  FI2(M1, M2, M3u)*FI2(M1, M2, M3u)*FI2(M1, M2, M3u) + &
                        &  FI2(M1, M2, M3l)*FI2(M1, M2, M3l)*FI2(M1, M2, M3l))/6.0D0
                    Delta21 = 27.0d0*P*P
                    tmp1 = 0.5d0*(DSQRT(Delta21*Delta21 - 4.0d0*Delta20*Delta20*Delta20) + Delta21)
                    tmp2 = tmp1 - Delta21
                    S = DSQRT(DSIGN(DABS(tmp1)**ONEOVERTHREE, tmp1) + Delta20/DSIGN(DABS(tmp1)**ONEOVERTHREE, tmp1))/DSQRT(12.0d0)
                    IF (P .GT. ZERO) THEN
                        FI2(M1, M2, M3) = -S + 0.5d0*DSQRT(-4.0d0*S*S + P/S)
                    ELSE IF (P .LT. -ZERO) THEN
                        FI2(M1, M2, M3) = S + 0.5d0*DSQRT(-4.0d0*S*S - P/S)
                    ELSE
                        FI2(M1, M2, M3) = DSQRT(DSQRT(-Q))
                    END IF
                END DO
            END DO
        END DO
    END IF

    CALL TimingMain(8, 1)

CONTAINS

    SUBROUTINE relax_coarse(u, rho, rhs, n)
    ! One red or black Gauss-Seidel sweep on a coarse level: the former FI3
    ! sections 4 (u), 2 (rho) and 5 (rhs), same expressions in the same order.
        integer, intent(in)    :: n
        real*4,  intent(inout) :: u(n, n, n)
        real*4,  intent(in)    :: rho(n, n, n), rhs(n, n, n)
        integer :: M1, M2, M3, M1u, M1l, M2u, M2l, M3u, M3l, M1start
        real*8  :: Delta10, Delta21, tmp1, tmp2, P, S, CC, sqrt_Delta0

        IF (fr_n .EQ. 0) THEN
!$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(STATIC) DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3,M1start,M1u,M1l,M2u,M2l,M3u,M3l) &
!$OMP PRIVATE (P)
            DO M3 = 1, n
                DO M2 = 1, n
                    M1start = 1 + MOD(M2 + M3 + parity_off, 2)
            !DIR$ IVDEP
                    DO M1 = M1start, n, 2
                        M1u = M1 + 1; IF (M1u > n) M1u = 1
                        M1l = M1 - 1; IF (M1l < 1) M1l = n
                        M2u = M2 + 1; IF (M2u > n) M2u = 1
                        M2l = M2 - 1; IF (M2l < 1) M2l = n
                        M3u = M3 + 1; IF (M3u > n) M3u = 1
                        M3l = M3 - 1; IF (M3l < 1) M3l = n
                        P = fct1*rho(M1, M2, M3) - Q2 + &
                          & fct2*rhs(M1, M2, M3) - &
                            & (u(M1u, M2, M3) + &
                            &  u(M1l, M2, M3) + &
                            &  u(M1, M2u, M3) + &
                            &  u(M1, M2l, M3) + &
                            &  u(M1, M2, M3u) + &
                            &  u(M1, M2, M3l))/6.0D0
                        u(M1, M2, M3) = 0.5D0*(-P + DSQRT(P*P - Q4))
                    END DO
                END DO
            END DO
        ELSE IF (fr_n .EQ. 1) THEN
!$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(STATIC) DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3,M1start,M1u,M1l,M2u,M2l,M3u,M3l) &
!$OMP PRIVATE (Delta10,tmp1,tmp2,sqrt_Delta0,P,S,CC)
            DO M3 = 1, n
                DO M2 = 1, n
                    M1start = 1 + MOD(M2 + M3 + parity_off, 2)
            !DIR$ IVDEP
                    DO M1 = M1start, n, 2
                        M1u = M1 + 1; IF (M1u > n) M1u = 1
                        M1l = M1 - 1; IF (M1l < 1) M1l = n
                        M2u = M2 + 1; IF (M2u > n) M2u = 1
                        M2l = M2 - 1; IF (M2l < 1) M2l = n
                        M3u = M3 + 1; IF (M3u > n) M3u = 1
                        M3l = M3 - 1; IF (M3l < 1) M3l = n
                        P = fct1*rho(M1, M2, M3) - Q2 + &
                          & fct2*rhs(M1, M2, M3) - &
                            & (u(M1u, M2, M3)*u(M1u, M2, M3) + &
                            &  u(M1l, M2, M3)*u(M1l, M2, M3) + &
                            &  u(M1, M2u, M3)*u(M1, M2u, M3) + &
                            &  u(M1, M2l, M3)*u(M1, M2l, M3) + &
                            &  u(M1, M2, M3u)*u(M1, M2, M3u) + &
                            &  u(M1, M2, M3l)*u(M1, M2, M3l))/6.0D0
                        Delta10 = -3.0D0*P
                        CC = Delta11*Delta11 - 4.0d0*Delta10*Delta10*Delta10
                        IF (CC .GE. 0.0D0) THEN
                            tmp1 = 0.5d0*(DSQRT(CC) + Delta11)
                            tmp2 = tmp1 - Delta11
                            IF (DABS(P) .LT. ZERO) THEN
                                u(M1, M2, M3) = Q2onethird
                            ELSE
                                u(M1, M2, M3) = -(DSIGN(DABS(tmp1)**ONEOVERTHREE, tmp1) - DSIGN(DABS(tmp2)**ONEOVERTHREE, tmp2))*ONEOVERTHREE
                            END IF
                        ELSE
                            sqrt_Delta0 = DSQRT(Delta10)
                            S = DACOS(0.5d0*Delta11/(sqrt_Delta0*sqrt_Delta0*sqrt_Delta0))
                            u(M1, M2, M3) = -TWOOVERTHREE*sqrt_Delta0*DCOS(S*ONEOVERTHREE + TWOPIOVERTHREE)
                        END IF
                    END DO
                END DO
            END DO
        ELSE IF (fr_n .EQ. 2) THEN
!$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(STATIC) DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3,M1start,M1u,M1l,M2u,M2l,M3u,M3l) &
!$OMP PRIVATE (Delta21,tmp1,tmp2,P,S)
            DO M3 = 1, n
                DO M2 = 1, n
                    M1start = 1 + MOD(M2 + M3 + parity_off, 2)
            !DIR$ IVDEP
                    DO M1 = M1start, n, 2
                        M1u = M1 + 1; IF (M1u > n) M1u = 1
                        M1l = M1 - 1; IF (M1l < 1) M1l = n
                        M2u = M2 + 1; IF (M2u > n) M2u = 1
                        M2l = M2 - 1; IF (M2l < 1) M2l = n
                        M3u = M3 + 1; IF (M3u > n) M3u = 1
                        M3l = M3 - 1; IF (M3l < 1) M3l = n
                        P = fct1*rho(M1, M2, M3) - Q2 + &
                          & fct2*rhs(M1, M2, M3) - &
                            & (u(M1u, M2, M3)*u(M1u, M2, M3)*u(M1u, M2, M3) + &
                            &  u(M1l, M2, M3)*u(M1l, M2, M3)*u(M1l, M2, M3) + &
                            &  u(M1, M2u, M3)*u(M1, M2u, M3)*u(M1, M2u, M3) + &
                            &  u(M1, M2l, M3)*u(M1, M2l, M3)*u(M1, M2l, M3) + &
                            &  u(M1, M2, M3u)*u(M1, M2, M3u)*u(M1, M2, M3u) + &
                            &  u(M1, M2, M3l)*u(M1, M2, M3l)*u(M1, M2, M3l))/6.0D0
                        Delta21 = 27.0d0*P*P
                        tmp1 = 0.5d0*(DSQRT(Delta21*Delta21 - 4.0d0*Delta20*Delta20*Delta20) + Delta21)
                        tmp2 = tmp1 - Delta21
                        S = DSQRT(DSIGN(DABS(tmp1)**ONEOVERTHREE, tmp1) + Delta20/DSIGN(DABS(tmp1)**ONEOVERTHREE, tmp1))/DSQRT(12.0d0)
                        IF (P .GT. ZERO) THEN
                            u(M1, M2, M3) = -S + 0.5d0*DSQRT(-4.0d0*S*S + P/S)
                        ELSE IF (P .LT. -ZERO) THEN
                            u(M1, M2, M3) = S + 0.5d0*DSQRT(-4.0d0*S*S - P/S)
                        ELSE
                            u(M1, M2, M3) = DSQRT(DSQRT(-Q))
                        END IF
                    END DO
                END DO
            END DO
        END IF
    END SUBROUTINE relax_coarse

END SUBROUTINE relaxation_iterations_fR

!-------------------------------------------------------------
!
! Calculate resudual and its RMS on the entire grid on ilevel.
!
!-------------------------------------------------------------
SUBROUTINE calculate_residual_fR(ilevel, res_PM_grid)
!-------------------------------------------------------------
    use Tools

    integer, intent(in) :: ilevel
    real*8, intent(out) :: res_PM_grid

    real*8  :: R0_bg, R_bg, fR_bg, ctilde, ctilde2, dx, dx2
    integer :: ngrid_level
    real*8  :: fct1, fct2, fct3
    real*8  :: s_fR_bg, a_fR_bg, a_fR_rt                            ! sign and amplitude of f_R(a)
    integer :: fr_n1

    integer :: M1, M2, M3, M1l, M1u, M2l, M2u, M3l, M3u
    integer :: N1, N2, N3
    real*8  :: RES, OP, RES2

    CALL TimingMain(8, -1)

    IF (MG_test .EQ. 1) WRITE (*, '(A,I5)') 'Calculate residual on level', levelmax - ilevel

    ! number of grid points on the coarse level
    ngrid_level = NGRID/2**ilevel

    ctilde = 2.99792458D3*DBLE(NGRID)/Box
    ctilde2 = ctilde**2
    dx = DBLE(2**ilevel)
    dx2 = dx*dx
    R_bg = 3.0d0*(Om/AEXPN**3 + 4.0d0*OmL)
    R0_bg = 3.0d0*(Om + 4.0d0*OmL)
    fR_bg = fR0*(R0_bg/R_bg)**(fr_n + 1)

    s_fR_bg = DSIGN(1.0D0, fR_bg)
    a_fR_bg = DABS(fR_bg)
    fr_n1 = fr_n + 1
    a_fR_rt = a_fR_bg**(1.0D0/(DBLE(fr_n1)))

    fct1 = Om/AEXPN
    fct2 = -R_bg*AEXPN**2/3.0d0
    fct3 = ctilde2*s_fR_bg/dx2

    RES = 0.0D0

    ! Hoisted dispatch on (fr_n, ilevel). Each specialization is a clean
    ! triple-loop with a straight-line cell body so ifx can vectorize the
    ! inner M1 loop via !DIR$ IVDEP. On ilevel==0 only the L^2 norm of the
    ! residual is needed (Phase E): the fine-level residual field itself was
    ! never read again (the restriction recomputes it from FI and FI2), so it
    ! is not stored. On a coarse level the residual field feeds the next
    ! restriction and the norm is not needed (res_PM_grid = 0 there).
    IF (ilevel .GT. 0) THEN
        CALL residual_coarse(MGL(ilevel)%u, MGL(ilevel)%rho, MGL(ilevel)%rhs, MGL(ilevel)%lres, ngrid_level)
    ELSE IF (fr_n .EQ. 0) THEN
!$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(STATIC) DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3,M1u,M1l,M2u,M2l,M3u,M3l,OP) &
!$OMP REDUCTION(+:RES)
            DO M3 = 1, ngrid_level
                DO M2 = 1, ngrid_level
!DIR$ IVDEP
                    DO M1 = 1, ngrid_level
                        M1u = M1 + 1; IF (M1u > ngrid_level) M1u = 1
                        M1l = M1 - 1; IF (M1l < 1) M1l = ngrid_level
                        M2u = M2 + 1; IF (M2u > ngrid_level) M2u = 1
                        M2l = M2 - 1; IF (M2l < 1) M2l = ngrid_level
                        M3u = M3 + 1; IF (M3u > ngrid_level) M3u = 1
                        M3l = M3 - 1; IF (M3l < 1) M3l = ngrid_level
                        OP = FI2(M1u, M2, M3) + &
                           & FI2(M1l, M2, M3) + &
                           & FI2(M1, M2u, M3) + &
                           & FI2(M1, M2l, M3) + &
                           & FI2(M1, M2, M3u) + &
                           & FI2(M1, M2, M3l) - &
                           & FI2(M1, M2, M3)*6.0D0
                        OP = OP*fct3
                        OP = OP + fct1*FI(M1, M2, M3)
                        OP = OP + fct2*(a_fR_rt/FI2(M1, M2, M3) - 1.0D0)
                        RES = RES + OP*OP
                    END DO
                END DO
            END DO
    ELSE IF (fr_n .EQ. 1) THEN
!$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(STATIC) DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3,M1u,M1l,M2u,M2l,M3u,M3l,OP) &
!$OMP REDUCTION(+:RES)
            DO M3 = 1, ngrid_level
                DO M2 = 1, ngrid_level
!DIR$ IVDEP
                    DO M1 = 1, ngrid_level
                        M1u = M1 + 1; IF (M1u > ngrid_level) M1u = 1
                        M1l = M1 - 1; IF (M1l < 1) M1l = ngrid_level
                        M2u = M2 + 1; IF (M2u > ngrid_level) M2u = 1
                        M2l = M2 - 1; IF (M2l < 1) M2l = ngrid_level
                        M3u = M3 + 1; IF (M3u > ngrid_level) M3u = 1
                        M3l = M3 - 1; IF (M3l < 1) M3l = ngrid_level
                        OP = FI2(M1u, M2, M3)*FI2(M1u, M2, M3) + &
                           & FI2(M1l, M2, M3)*FI2(M1l, M2, M3) + &
                           & FI2(M1, M2u, M3)*FI2(M1, M2u, M3) + &
                           & FI2(M1, M2l, M3)*FI2(M1, M2l, M3) + &
                           & FI2(M1, M2, M3u)*FI2(M1, M2, M3u) + &
                           & FI2(M1, M2, M3l)*FI2(M1, M2, M3l) - &
                           & FI2(M1, M2, M3)*FI2(M1, M2, M3)*6.0D0
                        OP = OP*fct3
                        OP = OP + fct1*FI(M1, M2, M3)
                        OP = OP + fct2*(a_fR_rt/FI2(M1, M2, M3) - 1.0D0)
                        RES = RES + OP*OP
                    END DO
                END DO
            END DO
    ELSE IF (fr_n .EQ. 2) THEN
!$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(STATIC) DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3,M1u,M1l,M2u,M2l,M3u,M3l,OP) &
!$OMP REDUCTION(+:RES)
            DO M3 = 1, ngrid_level
                DO M2 = 1, ngrid_level
!DIR$ IVDEP
                    DO M1 = 1, ngrid_level
                        M1u = M1 + 1; IF (M1u > ngrid_level) M1u = 1
                        M1l = M1 - 1; IF (M1l < 1) M1l = ngrid_level
                        M2u = M2 + 1; IF (M2u > ngrid_level) M2u = 1
                        M2l = M2 - 1; IF (M2l < 1) M2l = ngrid_level
                        M3u = M3 + 1; IF (M3u > ngrid_level) M3u = 1
                        M3l = M3 - 1; IF (M3l < 1) M3l = ngrid_level
                        OP = FI2(M1u, M2, M3)*FI2(M1u, M2, M3)*FI2(M1u, M2, M3) + &
                           & FI2(M1l, M2, M3)*FI2(M1l, M2, M3)*FI2(M1l, M2, M3) + &
                           & FI2(M1, M2u, M3)*FI2(M1, M2u, M3)*FI2(M1, M2u, M3) + &
                           & FI2(M1, M2l, M3)*FI2(M1, M2l, M3)*FI2(M1, M2l, M3) + &
                           & FI2(M1, M2, M3u)*FI2(M1, M2, M3u)*FI2(M1, M2, M3u) + &
                           & FI2(M1, M2, M3l)*FI2(M1, M2, M3l)*FI2(M1, M2, M3l) - &
                           & FI2(M1, M2, M3)*FI2(M1, M2, M3)*FI2(M1, M2, M3)*6.0D0
                        OP = OP*fct3
                        OP = OP + fct1*FI(M1, M2, M3)
                        OP = OP + fct2*(a_fR_rt/FI2(M1, M2, M3) - 1.0D0)
                        RES = RES + OP*OP
                    END DO
                END DO
            END DO
    END IF

    CALL TimingMain(8, 1)

    RES = DSQRT(RES/(DBLE(ngrid_level))**3)
    res_PM_grid = RES

! WRITE(*,'(A,F20.16,A,I5,A)') 'Residual RMS = ',res_PM_grid,' on level ',levelmax-ilevel
! WRITE(*,'(A,F20.16,F20.16,A,I5,A)') 'Residual RES2 = ',RES2,DSQRT(DABS(fR_bg)),' on level ',levelmax-ilevel

CONTAINS

    SUBROUTINE residual_coarse(u, rho, rhs, lres, n)
    ! Residual of the coarse-level equation (former FI3 section 6) from the
    ! working field u, the restricted density rho and the physical RHS rhs.
        integer, intent(in)  :: n
        real*4,  intent(in)  :: u(n, n, n), rho(n, n, n), rhs(n, n, n)
        real*4,  intent(out) :: lres(n, n, n)
        integer :: M1, M2, M3, M1u, M1l, M2u, M2l, M3u, M3l
        real*8  :: OP

        IF (fr_n .EQ. 0) THEN
!$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(STATIC) DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3,M1u,M1l,M2u,M2l,M3u,M3l,OP)
            DO M3 = 1, n
                DO M2 = 1, n
!DIR$ IVDEP
                    DO M1 = 1, n
                        M1u = M1 + 1; IF (M1u > n) M1u = 1
                        M1l = M1 - 1; IF (M1l < 1) M1l = n
                        M2u = M2 + 1; IF (M2u > n) M2u = 1
                        M2l = M2 - 1; IF (M2l < 1) M2l = n
                        M3u = M3 + 1; IF (M3u > n) M3u = 1
                        M3l = M3 - 1; IF (M3l < 1) M3l = n
                        OP = u(M1u, M2, M3) + &
                           & u(M1l, M2, M3) + &
                           & u(M1, M2u, M3) + &
                           & u(M1, M2l, M3) + &
                           & u(M1, M2, M3u) + &
                           & u(M1, M2, M3l) - &
                           & u(M1, M2, M3)*6.0D0
                        OP = OP*fct3
                        OP = OP + fct1*rho(M1, M2, M3)
                        OP = OP + fct2*(a_fR_rt/u(M1, M2, M3) - 1.0D0)
                        OP = OP - rhs(M1, M2, M3)
                        lres(M1, M2, M3) = OP
                    END DO
                END DO
            END DO
        ELSE IF (fr_n .EQ. 1) THEN
!$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(STATIC) DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3,M1u,M1l,M2u,M2l,M3u,M3l,OP)
            DO M3 = 1, n
                DO M2 = 1, n
!DIR$ IVDEP
                    DO M1 = 1, n
                        M1u = M1 + 1; IF (M1u > n) M1u = 1
                        M1l = M1 - 1; IF (M1l < 1) M1l = n
                        M2u = M2 + 1; IF (M2u > n) M2u = 1
                        M2l = M2 - 1; IF (M2l < 1) M2l = n
                        M3u = M3 + 1; IF (M3u > n) M3u = 1
                        M3l = M3 - 1; IF (M3l < 1) M3l = n
                        OP = u(M1u, M2, M3)*u(M1u, M2, M3) + &
                           & u(M1l, M2, M3)*u(M1l, M2, M3) + &
                           & u(M1, M2u, M3)*u(M1, M2u, M3) + &
                           & u(M1, M2l, M3)*u(M1, M2l, M3) + &
                           & u(M1, M2, M3u)*u(M1, M2, M3u) + &
                           & u(M1, M2, M3l)*u(M1, M2, M3l) - &
                           & u(M1, M2, M3)*u(M1, M2, M3)*6.0D0
                        OP = OP*fct3
                        OP = OP + fct1*rho(M1, M2, M3)
                        OP = OP + fct2*(a_fR_rt/u(M1, M2, M3) - 1.0d0)
                        OP = OP - rhs(M1, M2, M3)
                        lres(M1, M2, M3) = OP
                    END DO
                END DO
            END DO
        ELSE IF (fr_n .EQ. 2) THEN
!$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(STATIC) DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3,M1u,M1l,M2u,M2l,M3u,M3l,OP)
            DO M3 = 1, n
                DO M2 = 1, n
!DIR$ IVDEP
                    DO M1 = 1, n
                        M1u = M1 + 1; IF (M1u > n) M1u = 1
                        M1l = M1 - 1; IF (M1l < 1) M1l = n
                        M2u = M2 + 1; IF (M2u > n) M2u = 1
                        M2l = M2 - 1; IF (M2l < 1) M2l = n
                        M3u = M3 + 1; IF (M3u > n) M3u = 1
                        M3l = M3 - 1; IF (M3l < 1) M3l = n
                        OP = u(M1u, M2, M3)*u(M1u, M2, M3)*u(M1u, M2, M3) + &
                           & u(M1l, M2, M3)*u(M1l, M2, M3)*u(M1l, M2, M3) + &
                           & u(M1, M2u, M3)*u(M1, M2u, M3)*u(M1, M2u, M3) + &
                           & u(M1, M2l, M3)*u(M1, M2l, M3)*u(M1, M2l, M3) + &
                           & u(M1, M2, M3u)*u(M1, M2, M3u)*u(M1, M2, M3u) + &
                           & u(M1, M2, M3l)*u(M1, M2, M3l)*u(M1, M2, M3l) - &
                           & u(M1, M2, M3)*u(M1, M2, M3)*u(M1, M2, M3)*6.0D0
                        OP = OP*fct3
                        OP = OP + fct1*rho(M1, M2, M3)
                        OP = OP + fct2*(a_fR_rt/u(M1, M2, M3) - 1.0D0)
                        OP = OP - rhs(M1, M2, M3)
                        lres(M1, M2, M3) = OP
                    END DO
                END DO
            END DO
        END IF
    END SUBROUTINE residual_coarse

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
    real*8  :: R0_bg, R_bg, fR_bg, ctilde, ctilde2, dx, dx2
    real*8  :: fct1, fct2, fct3, fct4
    integer :: ngrid_level
    real*8  :: a_fR_bg, s_fR_bg                                   ! amplitude and sign of f_R(a)

    CALL TimingMain(8, -1)

    IF (MG_test .EQ. 1) WRITE (*, '(A,I5)') 'Restrict residual to level', levelmax - ilevel

    ! number of grid points on the coarse level
    ngrid_level = NGRID/2**ilevel

    ctilde = 2.99792458D3*DBLE(NGRID)/Box
    ctilde2 = ctilde**2
    dx = 1.0d0
    dx2 = dx*dx
    R_bg = 3.0d0*(Om/AEXPN**3 + 4.0d0*OmL)
    R0_bg = 3.0d0*(Om + 4.0d0*OmL)
    fR_bg = fR0*(R0_bg/R_bg)**(fr_n + 1)

    s_fR_bg = DSIGN(1.0D0, fR_bg)
    a_fR_bg = DABS(fR_bg)

    fct1 = Om/AEXPN
    fct2 = -R_bg*AEXPN**2/3.0d0
    fct3 = ctilde2*s_fR_bg/dx2
    fct4 = fct2*a_fR_bg**(1.0D0/(1.0D0 + DBLE(fr_n)))

    IF (ilevel .EQ. 1) THEN
        ! fine level = PM grid: recomputed from FI and FI2 (the PM-grid residual is not stored)
        CALL restrict_residual_fine(MGL(1)%res, ngrid_level)
    ELSE
        ! fine level = coarse level ilevel-1: mean of its 8 son-cell residuals
        CALL restrict_residual_coarse(MGL(ilevel - 1)%lres, 2*ngrid_level, MGL(ilevel)%res, ngrid_level)
    END IF

    CALL TimingMain(8, 1)

CONTAINS

    SUBROUTINE restrict_residual_fine(res, n)
    ! Former FI3 section 1 of level 1, from the PM-grid fields FI and FI2.
        integer, intent(in)  :: n
        real*4,  intent(out) :: res(n, n, n)
        integer :: M1, M2, M3, M1u, M1l, M2u, M2l, M3u, M3l
        real*8  :: P, Q

        ! Hoisted dispatch on fr_n; Q accumulation is model-agnostic (shared).
        ! Inner M1 loop gets !DIR$ IVDEP so ifx can vectorize the 32-ary P
        ! accumulation + 8-ary 1/FI2 division across M1.
        IF (fr_n .EQ. 0) THEN
!$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(STATIC) DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3,M1u,M1l,M2u,M2l,M3u,M3l,P,Q)
            DO M3 = 1, n
                DO M2 = 1, n
!DIR$ IVDEP
                    DO M1 = 1, n
                        Q = -8.0d0*fct2
                        Q = Q + fct1*FI(2*M1 - 1, 2*M2 - 1, 2*M3 - 1) + fct4/FI2(2*M1 - 1, 2*M2 - 1, 2*M3 - 1)
                        Q = Q + fct1*FI(2*M1, 2*M2 - 1, 2*M3 - 1) + fct4/FI2(2*M1, 2*M2 - 1, 2*M3 - 1)
                        Q = Q + fct1*FI(2*M1 - 1, 2*M2, 2*M3 - 1) + fct4/FI2(2*M1 - 1, 2*M2, 2*M3 - 1)
                        Q = Q + fct1*FI(2*M1, 2*M2, 2*M3 - 1) + fct4/FI2(2*M1, 2*M2, 2*M3 - 1)
                        Q = Q + fct1*FI(2*M1 - 1, 2*M2 - 1, 2*M3) + fct4/FI2(2*M1 - 1, 2*M2 - 1, 2*M3)
                        Q = Q + fct1*FI(2*M1, 2*M2 - 1, 2*M3) + fct4/FI2(2*M1, 2*M2 - 1, 2*M3)
                        Q = Q + fct1*FI(2*M1 - 1, 2*M2, 2*M3) + fct4/FI2(2*M1 - 1, 2*M2, 2*M3)
                        Q = Q + fct1*FI(2*M1, 2*M2, 2*M3) + fct4/FI2(2*M1, 2*M2, 2*M3)
                        M1u = 2*M1 + 1; IF (M1u > NGRID) M1u = 1
                        M1l = 2*M1 - 2; IF (M1l < 1) M1l = NGRID
                        M2u = 2*M2 + 1; IF (M2u > NGRID) M2u = 1
                        M2l = 2*M2 - 2; IF (M2l < 1) M2l = NGRID
                        M3u = 2*M3 + 1; IF (M3u > NGRID) M3u = 1
                        M3l = 2*M3 - 2; IF (M3l < 1) M3l = NGRID
                        P = -3.0D0*(FI2(2*M1 - 1, 2*M2 - 1, 2*M3 - 1) + &
                                  & FI2(2*M1, 2*M2 - 1, 2*M3 - 1) + &
                                  & FI2(2*M1 - 1, 2*M2, 2*M3 - 1) + &
                                  & FI2(2*M1, 2*M2, 2*M3 - 1) + &
                                  & FI2(2*M1 - 1, 2*M2 - 1, 2*M3) + &
                                  & FI2(2*M1, 2*M2 - 1, 2*M3) + &
                                  & FI2(2*M1 - 1, 2*M2, 2*M3) + &
                                  & FI2(2*M1, 2*M2, 2*M3))
                        P = P + FI2(M1l, 2*M2 - 1, 2*M3 - 1) + FI2(M1l, 2*M2, 2*M3 - 1) + &
                              & FI2(M1l, 2*M2 - 1, 2*M3) + FI2(M1l, 2*M2, 2*M3)
                        P = P + FI2(M1u, 2*M2 - 1, 2*M3 - 1) + FI2(M1u, 2*M2, 2*M3 - 1) + &
                              & FI2(M1u, 2*M2 - 1, 2*M3) + FI2(M1u, 2*M2, 2*M3)
                        P = P + FI2(2*M1 - 1, M2l, 2*M3 - 1) + FI2(2*M1, M2l, 2*M3 - 1) + &
                              & FI2(2*M1 - 1, M2l, 2*M3) + FI2(2*M1, M2l, 2*M3)
                        P = P + FI2(2*M1 - 1, M2u, 2*M3 - 1) + FI2(2*M1, M2u, 2*M3 - 1) + &
                              & FI2(2*M1 - 1, M2u, 2*M3) + FI2(2*M1, M2u, 2*M3)
                        P = P + FI2(2*M1 - 1, 2*M2 - 1, M3l) + FI2(2*M1, 2*M2 - 1, M3l) + &
                              & FI2(2*M1 - 1, 2*M2, M3l) + FI2(2*M1, 2*M2, M3l)
                        P = P + FI2(2*M1 - 1, 2*M2 - 1, M3u) + FI2(2*M1, 2*M2 - 1, M3u) + &
                              & FI2(2*M1 - 1, 2*M2, M3u) + FI2(2*M1, 2*M2, M3u)
                        res(M1, M2, M3) = (P*fct3 + Q)/8.0D0
                    END DO
                END DO
            END DO
        ELSE IF (fr_n .EQ. 1) THEN
!$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(STATIC) DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3,M1u,M1l,M2u,M2l,M3u,M3l,P,Q)
            DO M3 = 1, n
                DO M2 = 1, n
!DIR$ IVDEP
                    DO M1 = 1, n
                        Q = -8.0d0*fct2
                        Q = Q + fct1*FI(2*M1 - 1, 2*M2 - 1, 2*M3 - 1) + fct4/FI2(2*M1 - 1, 2*M2 - 1, 2*M3 - 1)
                        Q = Q + fct1*FI(2*M1, 2*M2 - 1, 2*M3 - 1) + fct4/FI2(2*M1, 2*M2 - 1, 2*M3 - 1)
                        Q = Q + fct1*FI(2*M1 - 1, 2*M2, 2*M3 - 1) + fct4/FI2(2*M1 - 1, 2*M2, 2*M3 - 1)
                        Q = Q + fct1*FI(2*M1, 2*M2, 2*M3 - 1) + fct4/FI2(2*M1, 2*M2, 2*M3 - 1)
                        Q = Q + fct1*FI(2*M1 - 1, 2*M2 - 1, 2*M3) + fct4/FI2(2*M1 - 1, 2*M2 - 1, 2*M3)
                        Q = Q + fct1*FI(2*M1, 2*M2 - 1, 2*M3) + fct4/FI2(2*M1, 2*M2 - 1, 2*M3)
                        Q = Q + fct1*FI(2*M1 - 1, 2*M2, 2*M3) + fct4/FI2(2*M1 - 1, 2*M2, 2*M3)
                        Q = Q + fct1*FI(2*M1, 2*M2, 2*M3) + fct4/FI2(2*M1, 2*M2, 2*M3)
                        M1u = 2*M1 + 1; IF (M1u > NGRID) M1u = 1
                        M1l = 2*M1 - 2; IF (M1l < 1) M1l = NGRID
                        M2u = 2*M2 + 1; IF (M2u > NGRID) M2u = 1
                        M2l = 2*M2 - 2; IF (M2l < 1) M2l = NGRID
                        M3u = 2*M3 + 1; IF (M3u > NGRID) M3u = 1
                        M3l = 2*M3 - 2; IF (M3l < 1) M3l = NGRID
                        P = -3.0D0*(FI2(2*M1 - 1, 2*M2 - 1, 2*M3 - 1)*FI2(2*M1 - 1, 2*M2 - 1, 2*M3 - 1) + &
                                  & FI2(2*M1, 2*M2 - 1, 2*M3 - 1)*FI2(2*M1, 2*M2 - 1, 2*M3 - 1) + &
                                  & FI2(2*M1 - 1, 2*M2, 2*M3 - 1)*FI2(2*M1 - 1, 2*M2, 2*M3 - 1) + &
                                  & FI2(2*M1, 2*M2, 2*M3 - 1)*FI2(2*M1, 2*M2, 2*M3 - 1) + &
                                  & FI2(2*M1 - 1, 2*M2 - 1, 2*M3)*FI2(2*M1 - 1, 2*M2 - 1, 2*M3) + &
                                  & FI2(2*M1, 2*M2 - 1, 2*M3)*FI2(2*M1, 2*M2 - 1, 2*M3) + &
                                  & FI2(2*M1 - 1, 2*M2, 2*M3)*FI2(2*M1 - 1, 2*M2, 2*M3) + &
                                  & FI2(2*M1, 2*M2, 2*M3)*FI2(2*M1, 2*M2, 2*M3))
                        P = P + FI2(M1l, 2*M2 - 1, 2*M3 - 1)*FI2(M1l, 2*M2 - 1, 2*M3 - 1) + FI2(M1l, 2*M2, 2*M3 - 1)*FI2(M1l, 2*M2, 2*M3 - 1) + &
                              & FI2(M1l, 2*M2 - 1, 2*M3)*FI2(M1l, 2*M2 - 1, 2*M3) + FI2(M1l, 2*M2, 2*M3)*FI2(M1l, 2*M2, 2*M3)
                        P = P + FI2(M1u, 2*M2 - 1, 2*M3 - 1)*FI2(M1u, 2*M2 - 1, 2*M3 - 1) + FI2(M1u, 2*M2, 2*M3 - 1)*FI2(M1u, 2*M2, 2*M3 - 1) + &
                              & FI2(M1u, 2*M2 - 1, 2*M3)*FI2(M1u, 2*M2 - 1, 2*M3) + FI2(M1u, 2*M2, 2*M3)*FI2(M1u, 2*M2, 2*M3)
                        P = P + FI2(2*M1 - 1, M2l, 2*M3 - 1)*FI2(2*M1 - 1, M2l, 2*M3 - 1) + FI2(2*M1, M2l, 2*M3 - 1)*FI2(2*M1, M2l, 2*M3 - 1) + &
                              & FI2(2*M1 - 1, M2l, 2*M3)*FI2(2*M1 - 1, M2l, 2*M3) + FI2(2*M1, M2l, 2*M3)*FI2(2*M1, M2l, 2*M3)
                        P = P + FI2(2*M1 - 1, M2u, 2*M3 - 1)*FI2(2*M1 - 1, M2u, 2*M3 - 1) + FI2(2*M1, M2u, 2*M3 - 1)*FI2(2*M1, M2u, 2*M3 - 1) + &
                              & FI2(2*M1 - 1, M2u, 2*M3)*FI2(2*M1 - 1, M2u, 2*M3) + FI2(2*M1, M2u, 2*M3)*FI2(2*M1, M2u, 2*M3)
                        P = P + FI2(2*M1 - 1, 2*M2 - 1, M3l)*FI2(2*M1 - 1, 2*M2 - 1, M3l) + FI2(2*M1, 2*M2 - 1, M3l)*FI2(2*M1, 2*M2 - 1, M3l) + &
                              & FI2(2*M1 - 1, 2*M2, M3l)*FI2(2*M1 - 1, 2*M2, M3l) + FI2(2*M1, 2*M2, M3l)*FI2(2*M1, 2*M2, M3l)
                        P = P + FI2(2*M1 - 1, 2*M2 - 1, M3u)*FI2(2*M1 - 1, 2*M2 - 1, M3u) + FI2(2*M1, 2*M2 - 1, M3u)*FI2(2*M1, 2*M2 - 1, M3u) + &
                              & FI2(2*M1 - 1, 2*M2, M3u)*FI2(2*M1 - 1, 2*M2, M3u) + FI2(2*M1, 2*M2, M3u)*FI2(2*M1, 2*M2, M3u)
                        res(M1, M2, M3) = (P*fct3 + Q)/8.0D0
                    END DO
                END DO
            END DO
        ELSE IF (fr_n .EQ. 2) THEN
!$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(STATIC) DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3,M1u,M1l,M2u,M2l,M3u,M3l,P,Q)
            DO M3 = 1, n
                DO M2 = 1, n
!DIR$ IVDEP
                    DO M1 = 1, n
                        Q = -8.0d0*fct2
                        Q = Q + fct1*FI(2*M1 - 1, 2*M2 - 1, 2*M3 - 1) + fct4/FI2(2*M1 - 1, 2*M2 - 1, 2*M3 - 1)
                        Q = Q + fct1*FI(2*M1, 2*M2 - 1, 2*M3 - 1) + fct4/FI2(2*M1, 2*M2 - 1, 2*M3 - 1)
                        Q = Q + fct1*FI(2*M1 - 1, 2*M2, 2*M3 - 1) + fct4/FI2(2*M1 - 1, 2*M2, 2*M3 - 1)
                        Q = Q + fct1*FI(2*M1, 2*M2, 2*M3 - 1) + fct4/FI2(2*M1, 2*M2, 2*M3 - 1)
                        Q = Q + fct1*FI(2*M1 - 1, 2*M2 - 1, 2*M3) + fct4/FI2(2*M1 - 1, 2*M2 - 1, 2*M3)
                        Q = Q + fct1*FI(2*M1, 2*M2 - 1, 2*M3) + fct4/FI2(2*M1, 2*M2 - 1, 2*M3)
                        Q = Q + fct1*FI(2*M1 - 1, 2*M2, 2*M3) + fct4/FI2(2*M1 - 1, 2*M2, 2*M3)
                        Q = Q + fct1*FI(2*M1, 2*M2, 2*M3) + fct4/FI2(2*M1, 2*M2, 2*M3)
                        M1u = 2*M1 + 1; IF (M1u > NGRID) M1u = 1
                        M1l = 2*M1 - 2; IF (M1l < 1) M1l = NGRID
                        M2u = 2*M2 + 1; IF (M2u > NGRID) M2u = 1
                        M2l = 2*M2 - 2; IF (M2l < 1) M2l = NGRID
                        M3u = 2*M3 + 1; IF (M3u > NGRID) M3u = 1
                        M3l = 2*M3 - 2; IF (M3l < 1) M3l = NGRID
                        P = -3.0D0*(FI2(2*M1 - 1, 2*M2 - 1, 2*M3 - 1)*FI2(2*M1 - 1, 2*M2 - 1, 2*M3 - 1)*FI2(2*M1 - 1, 2*M2 - 1, 2*M3 - 1) + &
                                  & FI2(2*M1, 2*M2 - 1, 2*M3 - 1)*FI2(2*M1, 2*M2 - 1, 2*M3 - 1)*FI2(2*M1, 2*M2 - 1, 2*M3 - 1) + &
                                  & FI2(2*M1 - 1, 2*M2, 2*M3 - 1)*FI2(2*M1 - 1, 2*M2, 2*M3 - 1)*FI2(2*M1 - 1, 2*M2, 2*M3 - 1) + &
                                  & FI2(2*M1, 2*M2, 2*M3 - 1)*FI2(2*M1, 2*M2, 2*M3 - 1)*FI2(2*M1, 2*M2, 2*M3 - 1) + &
                                  & FI2(2*M1 - 1, 2*M2 - 1, 2*M3)*FI2(2*M1 - 1, 2*M2 - 1, 2*M3)*FI2(2*M1 - 1, 2*M2 - 1, 2*M3) + &
                                  & FI2(2*M1, 2*M2 - 1, 2*M3)*FI2(2*M1, 2*M2 - 1, 2*M3)*FI2(2*M1, 2*M2 - 1, 2*M3) + &
                                  & FI2(2*M1 - 1, 2*M2, 2*M3)*FI2(2*M1 - 1, 2*M2, 2*M3)*FI2(2*M1 - 1, 2*M2, 2*M3) + &
                                  & FI2(2*M1, 2*M2, 2*M3)*FI2(2*M1, 2*M2, 2*M3)*FI2(2*M1, 2*M2, 2*M3))
                        P = P + FI2(M1l, 2*M2 - 1, 2*M3 - 1)**3 + FI2(M1l, 2*M2, 2*M3 - 1)**3 + &
                              & FI2(M1l, 2*M2 - 1, 2*M3)**3 + FI2(M1l, 2*M2, 2*M3)**3
                        P = P + FI2(M1u, 2*M2 - 1, 2*M3 - 1)**3 + FI2(M1u, 2*M2, 2*M3 - 1)**3 + &
                              & FI2(M1u, 2*M2 - 1, 2*M3)**3 + FI2(M1u, 2*M2, 2*M3)**3
                        P = P + FI2(2*M1 - 1, M2l, 2*M3 - 1)**3 + FI2(2*M1, M2l, 2*M3 - 1)**3 + &
                              & FI2(2*M1 - 1, M2l, 2*M3)**3 + FI2(2*M1, M2l, 2*M3)**3
                        P = P + FI2(2*M1 - 1, M2u, 2*M3 - 1)**3 + FI2(2*M1, M2u, 2*M3 - 1)**3 + &
                              & FI2(2*M1 - 1, M2u, 2*M3)**3 + FI2(2*M1, M2u, 2*M3)**3
                        P = P + FI2(2*M1 - 1, 2*M2 - 1, M3l)**3 + FI2(2*M1, 2*M2 - 1, M3l)**3 + &
                              & FI2(2*M1 - 1, 2*M2, M3l)**3 + FI2(2*M1, 2*M2, M3l)**3
                        P = P + FI2(2*M1 - 1, 2*M2 - 1, M3u)**3 + FI2(2*M1, 2*M2 - 1, M3u)**3 + &
                              & FI2(2*M1 - 1, 2*M2, M3u)**3 + FI2(2*M1, 2*M2, M3u)**3
                        res(M1, M2, M3) = (P*fct3 + Q)/8.0D0
                    END DO
                END DO
            END DO
        END IF
    END SUBROUTINE restrict_residual_fine

    SUBROUTINE restrict_residual_coarse(lres, nf, res, n)
    ! Former FI3 section 1 of level ilevel from section 6 of level ilevel-1.
        integer, intent(in)  :: nf, n
        real*4,  intent(in)  :: lres(nf, nf, nf)
        real*4,  intent(out) :: res(n, n, n)
        integer :: M1, M2, M3
        real*8  :: P
!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC) DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3) &
!$OMP PRIVATE (P)
        DO M3 = 1, n
            DO M2 = 1, n
                DO M1 = 1, n
                    P = lres(2*M1 - 1, 2*M2 - 1, 2*M3 - 1)
                    P = P + lres(2*M1, 2*M2 - 1, 2*M3 - 1)
                    P = P + lres(2*M1 - 1, 2*M2, 2*M3 - 1)
                    P = P + lres(2*M1, 2*M2, 2*M3 - 1)
                    P = P + lres(2*M1 - 1, 2*M2 - 1, 2*M3)
                    P = P + lres(2*M1, 2*M2 - 1, 2*M3)
                    P = P + lres(2*M1 - 1, 2*M2, 2*M3)
                    P = P + lres(2*M1, 2*M2, 2*M3)
                    !
                    res(M1, M2, M3) = P/8.0D0
                    !
                END DO
            END DO
        END DO
    END SUBROUTINE restrict_residual_coarse

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
    real*8  :: R0_bg, R_bg, fR_bg, ctilde, dx, dx2
    real*8  :: fct1, fct2, fct3
    integer :: ngrid_level

    real*8  :: s_fR_bg, a_fR_bg, a_fR_rt                            ! sign and amplitude of f_R(a)
    integer :: fr_n1

    CALL TimingMain(8, -1)

    IF (MG_test .EQ. 1) WRITE (*, '(A,I5)') 'Calculate physical right-hand side on level', levelmax - ilevel

    ! number of grid points on the coarse level
    ngrid_level = NGRID/2**ilevel

    ctilde = 2.99792458D3*DBLE(NGRID)/Box
    dx = DBLE(2**ilevel)
    dx2 = dx*dx
    R_bg = 3.0d0*(Om/AEXPN**3 + 4.0d0*OmL)
    R0_bg = 3.0d0*(Om + 4.0d0*OmL)
    fR_bg = fR0*(R0_bg/R_bg)**(fr_n + 1)

    s_fR_bg = DSIGN(1.0D0, fR_bg)
    a_fR_bg = DABS(fR_bg)
    fr_n1 = fr_n + 1
    a_fR_rt = a_fR_bg**(1.0D0/(DBLE(fr_n1)))

    fct1 = Om/AEXPN
    fct2 = -R_bg*AEXPN**2/3.0d0
    fct3 = ctilde**2*s_fR_bg/dx2

    ! operator applied to the restricted field (ucopy, former section 3),
    ! restricted density (rho, section 2), minus the restricted residual
    ! (res, section 1) -> physical RHS (rhs, section 5)
    CALL rhs_coarse(MGL(ilevel)%ucopy, MGL(ilevel)%rho, MGL(ilevel)%res, MGL(ilevel)%rhs, ngrid_level)

    CALL TimingMain(8, 1)

CONTAINS

    SUBROUTINE rhs_coarse(ucopy, rho, res, rhs, n)
        integer, intent(in)  :: n
        real*4,  intent(in)  :: ucopy(n, n, n), rho(n, n, n), res(n, n, n)
        real*4,  intent(out) :: rhs(n, n, n)
        integer :: M1, M2, M3, M1l, M2l, M3l, M1u, M2u, M3u
        real*8  :: OP

        IF (fr_n .EQ. 0) THEN
!$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(STATIC) DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3,M1u,M1l,M2u,M2l,M3u,M3l) &
!$OMP PRIVATE (OP)
            DO M3 = 1, n
                DO M2 = 1, n
                    M2u = M2 + 1; IF (M2u > n) M2u = 1
                    M2l = M2 - 1; IF (M2l < 1) M2l = n
                    M3u = M3 + 1; IF (M3u > n) M3u = 1
                    M3l = M3 - 1; IF (M3l < 1) M3l = n
            !DIR$ IVDEP
                    DO M1 = 1, n
                        M1u = M1 + 1; IF (M1u > n) M1u = 1
                        M1l = M1 - 1; IF (M1l < 1) M1l = n
                        OP = ucopy(M1u, M2, M3) + &
                           & ucopy(M1l, M2, M3) + &
                           & ucopy(M1, M2u, M3) + &
                           & ucopy(M1, M2l, M3) + &
                           & ucopy(M1, M2, M3u) + &
                           & ucopy(M1, M2, M3l) - &
                           & ucopy(M1, M2, M3)*6.0D0
                        OP = OP*fct3
                        OP = OP + fct1*rho(M1, M2, M3)
                        OP = OP + fct2*(a_fR_rt/ucopy(M1, M2, M3) - 1.0D0)
                        OP = OP - res(M1, M2, M3)
                        rhs(M1, M2, M3) = OP
                    END DO
                END DO
            END DO
        ELSE IF (fr_n .EQ. 1) THEN
!$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(STATIC) DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3,M1u,M1l,M2u,M2l,M3u,M3l) &
!$OMP PRIVATE (OP)
            DO M3 = 1, n
                DO M2 = 1, n
                    M2u = M2 + 1; IF (M2u > n) M2u = 1
                    M2l = M2 - 1; IF (M2l < 1) M2l = n
                    M3u = M3 + 1; IF (M3u > n) M3u = 1
                    M3l = M3 - 1; IF (M3l < 1) M3l = n
            !DIR$ IVDEP
                    DO M1 = 1, n
                        M1u = M1 + 1; IF (M1u > n) M1u = 1
                        M1l = M1 - 1; IF (M1l < 1) M1l = n
                        OP = ucopy(M1u, M2, M3)*ucopy(M1u, M2, M3) + &
                           & ucopy(M1l, M2, M3)*ucopy(M1l, M2, M3) + &
                           & ucopy(M1, M2u, M3)*ucopy(M1, M2u, M3) + &
                           & ucopy(M1, M2l, M3)*ucopy(M1, M2l, M3) + &
                           & ucopy(M1, M2, M3u)*ucopy(M1, M2, M3u) + &
                           & ucopy(M1, M2, M3l)*ucopy(M1, M2, M3l) - &
                           & ucopy(M1, M2, M3)*ucopy(M1, M2, M3)*6.0D0
                        OP = OP*fct3
                        OP = OP + fct1*rho(M1, M2, M3)
                        OP = OP + fct2*(a_fR_rt/ucopy(M1, M2, M3) - 1.0D0)
                        OP = OP - res(M1, M2, M3)
                        rhs(M1, M2, M3) = OP
                    END DO
                END DO
            END DO
        ELSE IF (fr_n .EQ. 2) THEN
!$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(STATIC) DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3,M1u,M1l,M2u,M2l,M3u,M3l) &
!$OMP PRIVATE (OP)
            DO M3 = 1, n
                DO M2 = 1, n
                    M2u = M2 + 1; IF (M2u > n) M2u = 1
                    M2l = M2 - 1; IF (M2l < 1) M2l = n
                    M3u = M3 + 1; IF (M3u > n) M3u = 1
                    M3l = M3 - 1; IF (M3l < 1) M3l = n
            !DIR$ IVDEP
                    DO M1 = 1, n
                        M1u = M1 + 1; IF (M1u > n) M1u = 1
                        M1l = M1 - 1; IF (M1l < 1) M1l = n
                        OP = ucopy(M1u, M2, M3)*ucopy(M1u, M2, M3)*ucopy(M1u, M2, M3) + &
                           & ucopy(M1l, M2, M3)*ucopy(M1l, M2, M3)*ucopy(M1l, M2, M3) + &
                           & ucopy(M1, M2u, M3)*ucopy(M1, M2u, M3)*ucopy(M1, M2u, M3) + &
                           & ucopy(M1, M2l, M3)*ucopy(M1, M2l, M3)*ucopy(M1, M2l, M3) + &
                           & ucopy(M1, M2, M3u)*ucopy(M1, M2, M3u)*ucopy(M1, M2, M3u) + &
                           & ucopy(M1, M2, M3l)*ucopy(M1, M2, M3l)*ucopy(M1, M2, M3l) - &
                           & ucopy(M1, M2, M3)*ucopy(M1, M2, M3)*ucopy(M1, M2, M3)*6.0D0
                        OP = OP*fct3
                        OP = OP + fct1*rho(M1, M2, M3)
                        OP = OP + fct2*(a_fR_rt/ucopy(M1, M2, M3) - 1.0D0)
                        OP = OP - res(M1, M2, M3)
                        rhs(M1, M2, M3) = OP
                    END DO
                END DO
            END DO
        END IF
    END SUBROUTINE rhs_coarse

END SUBROUTINE calculate_physical_right_hand_side_fR
