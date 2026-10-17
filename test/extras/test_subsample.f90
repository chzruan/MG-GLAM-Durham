!-------------------------------------------------
!   Unit test for module Subsample (PMP2mod_subsample.f90, item I9).
!   Stand-alone: fills the Tools particle arrays of a 32^3 toy run
!   with a deterministic pseudo-random pattern, dumps the raw arrays
!   for the Python reference (check_subsample.py) and writes
!   Subsample files for several M and chunk sizes:
!     step  7: M = 100, SubChunk = 1000 (33 chunks, partial last one)
!     step  8: M = 1   ("full" mode), SubChunk = 1000
!     step  9: M = 100, SubChunk = default -> must be bitwise equal to step 7
!     step 10: M = 0   -> must write nothing
!     step 11: M = 7,  SubChunk = 1000
!   Driven by test/extras/run_subsample_test.sh.
!-------------------------------------------------
Program test_subsample
    use Tools
    use Subsample
    use iso_c_binding, only: c_int64_t
    implicit none
    Integer*8 :: i, state, k
    Integer*4 :: iu
    Real*4 :: u(6)

    NGRID = 64
    NROW  = 32
    Box   = 256.
    AEXPN = 0.5
    ASTEP = 0.01
    ISTEP = 7
    Nrealization = 1
    Om = 0.3
    OmL = 0.7
    hubble = 0.7
    MG_flag = 0
    fR0 = 0.d0
    HEADER = 'subsample unit test'
    Nparticles = int(NROW, 8)**3
    allocate (XPAR(Nparticles), YPAR(Nparticles), ZPAR(Nparticles))
    allocate (VX(Nparticles), VY(Nparticles), VZ(Nparticles))

    !---- MINSTD generator: positions uniform in [1, 65), velocities uniform in [-5, 5)
    !     code units (-> about +-4040 km/s, well inside the int16 range)
    state = 20261017_8
    Do i = 1, Nparticles
        Call Draw(u)
        XPAR(i) = 1. + 64.*u(1)
        YPAR(i) = 1. + 64.*u(2)
        ZPAR(i) = 1. + 64.*u(3)
        VX(i) = 10.*(u(4) - 0.5)
        VY(i) = 10.*(u(5) - 0.5)
        VZ(i) = 10.*(u(6) - 0.5)
    End Do
    !---- edge cases for the position wrap/quantisation
    XPAR(1) = 1.0                      ! lower edge        -> q = 0
    YPAR(1) = 65.0                     ! upper edge, wraps -> q = 0
    ZPAR(1) = 0.5                      ! below the box, wraps -> 63.5 grid units
    XPAR(2) = 64.999996                ! largest real*4 below 65 -> q = 65535
    YPAR(2) = 33.0                     ! half box -> q = 32768 (stored 0)
    ZPAR(2) = 1.0 + 64.*(1./65536.)    ! one quantum -> q = 1
    !---- velocities that clip: the first five IDs kept at M = 100 (hence also at M = 1)
    k = 0
    Do i = 1, Nparticles
        If (glam_select(i - 1, SubSalt, 100_c_int64_t) /= 0) Then
            k = k + 1
            VX(i) = 15.*(-1)**k            ! about +-12100 km/s -> |s| = 48500 > 32767
            VY(i) = 11.                    ! about  8900 km/s
            VZ(i) = -1000.
            write (*, '(a,i8)') ' clipping velocities assigned to ID', i - 1
            If (k == 5) exit
        End If
    End Do

    !---- raw arrays for the reference decoder (little-endian float32, no markers)
    open (newunit=iu, file='subsample_test_particles.bin', status='replace', access='stream', &
          form='unformatted', convert='little_endian')
    write (iu) XPAR, YPAR, ZPAR, VX, VY, VZ
    close (iu)

    SubChunk = 1000_8
    Call WriteSubsample(7, 100)
    Call WriteSubsample(8, 1)
    SubChunk = 1048576_8
    Call WriteSubsample(9, 100)
    Call WriteSubsample(10, 0)
    SubChunk = 1000_8
    Call WriteSubsample(11, 7)
    write (*, '(a)') ' test_subsample: done'

Contains
    Subroutine Draw(u)
        Real*4, intent(out) :: u(6)
        Integer*4 :: j
        Do j = 1, 6
            state = mod(48271_8*state, 2147483647_8)
            u(j) = real(dble(state)/2147483647.d0)
            If (u(j) >= 1.) u(j) = 0.
        End Do
    End Subroutine Draw
End Program test_subsample
