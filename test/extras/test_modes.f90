! V1 layout unit test for the low-k mode dump (docs/degrace_companion_implementation_plan.md, section 8).
! Fills a 16^3 mesh with cos(2 pi m1.x/N) + 0.5 sin(2 pi m2.x/N), transforms it with
! FFT3Dforward and writes the mode cube with DumpModes; check_modes.py verifies amplitudes and signs.
Program test_modes
    use Tools
    use ExtraOutputs
    implicit none
    integer*4, parameter :: N = 16
    integer*4 :: ix, iy, iz
    real*8 :: ph1, ph2, twopi
    twopi = 8.d0*atan(1.d0)
    NGRID = N; NROW = 8; Box = 100.; AEXPN = 0.5; ASTEP = 0.01; ISTEP = 1; Nrealization = 1
    Om = 0.3; OmL = 0.7; hubble = 0.7; MG_flag = 0; fR0 = 0.d0; HEADER = 'mode layout unit test'
    Nparticles = int(NROW, 8)**3
    allocate (FI(N, N, N))
    Do iz = 1, N
    Do iy = 1, N
    Do ix = 1, N
        ph1 = twopi*(2.d0*(ix - 1) - 1.d0*(iy - 1) + 3.d0*(iz - 1))/N     ! m1 = ( 2,-1, 3)
        ph2 = twopi*(1.d0*(ix - 1) + 2.d0*(iy - 1) + 0.d0*(iz - 1))/N     ! m2 = ( 1, 2, 0)
        FI(ix, iy, iz) = real(cos(ph1) + 0.5d0*sin(ph2))
    End Do
    End Do
    End Do
    Call FFT3Dforward(FI, N)
    Call DumpModes(FI, N, 'vtest', 0.3, 1, 0.495)      ! n = floor(0.3*100/2pi) = 4
    write (*, *) 'D(a=0.5) =', GrowthFactor(0.5), '  chi(a=0.5) =', ComovingDistance(0.5)
End Program test_modes
