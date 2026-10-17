!-------------------------------------------------
!   Extras: products added for the DEGRACE companion runs
!   (docs/degrace_companion_implementation_plan.md).
!     WriteOutputSchedule - step -> epoch index file (item I4)
!   Everything here is inert unless called; the runtime flags that
!   enable the products live in Tools (read from Setup.dat with
!   defaults = off, so old run directories behave as before).
!-------------------------------------------------
Module ExtraOutputs
    use Tools
    use fft5
    use Density
    implicit none
    integer*4, parameter, private :: Nlenmax = 16384
    real*8, save, private :: work(Nlenmax), wsave(Nlenmax)
!$OMP THREADPRIVATE(work, wsave)
    Logical, save :: fRHaloActive = .false.        ! f_R-at-haloes file open (set by fRHaloOpen/Close)
    Integer*4, save, private :: nfRrows = 0, fRmoment = 0
    Real*8, save, private :: fRbg_old = 0.d0, a_old_fR = 0.d0
Contains
!---------------------------------------------------------
!   Comoving distance chi(a) in Mpc/h for the flat matter+Lambda
!   background (Om, OmL from Setup.dat), c/H0 = 2997.92458 Mpc/h.
!   This is the convention of the DEGRACE output stations
!   chi = 256 k Mpc/h (no radiation), Simpson rule in ln a.
!---------------------------------------------------------
    Real*8 Function ComovingDistance(a)
        Real*4, intent(in) :: a
        Integer*4, parameter :: Nint = 4000
        Real*8 :: x0, h, x, s, f
        Integer*4 :: i
        If (a >= 1.) Then
            ComovingDistance = 0.d0
            Return
        End If
        x0 = log(dble(a))
        h = -x0/Nint
        s = 0.d0
        Do i = 0, Nint
            x = x0 + h*i
            f = 1.d0/(exp(x)*sqrt(dble(Om)*exp(-3.d0*x) + dble(OmL)))   ! d(ln a) / (a E(a))
            If (i == 0 .or. i == Nint) Then
                s = s + f
            Else If (mod(i, 2) == 1) Then
                s = s + 4.d0*f
            Else
                s = s + 2.d0*f
            End If
        End Do
        ComovingDistance = 2997.92458d0*s*h/3.d0
    End Function ComovingDistance
!---------------------------------------------------------
!   OutputSchedule.dat: one line per analysis output with the step
!   number that names the snapshot files, the epoch, the comoving
!   distance, whether the step is pinned on a requested redshift and
!   the per-output class. Written at start-up (and again at restart,
!   identically) in the run directory.
!---------------------------------------------------------
    Subroutine WriteOutputSchedule(Ntot, Alist, Nlist, MarkX, Nclass)
        Integer*4, intent(in) :: Ntot
        Real*4, intent(in)    :: Alist(Ntot)
        Integer*4, intent(in) :: Nlist(Ntot), MarkX(Ntot), Nclass(Ntot)
        Integer*4 :: i, n
        open (29, file='OutputSchedule.dat', status='unknown')
        write (29, '(a)') '# MG-GLAM output schedule: one line per analysis output (snapshot number = step)'
        write (29, '(a,f8.5,a,f8.5,a)') '# chi: comoving distance [Mpc/h], flat matter+Lambda, Om=', Om, &
            ' OmL=', OmL, ', c/H0 = 2997.92458 Mpc/h'
        write (29, '(a)') '# exact: 1 = step pinned on a requested redshift (#outputs < 0), 0 = nearest scheduled step'
        write (29, '(a)') '# class: 1 = full output (halo finder, extras), 2 = station (planes and P(k) only)'
        write (29, '(a)') '#  step               a               z          chi exact class'
        n = 0
        Do i = 1, Ntot
            If (Nlist(i) == 1) Then
                n = n + 1
                write (29, '(i6,2es16.8,f13.4,2i6)') i, Alist(i), 1./Alist(i) - 1., &
                    ComovingDistance(Alist(i)), MarkX(i), max(1, Nclass(i))
            End If
        End Do
        close (29)
        write (*, '(a,i4,a)') '  OutputSchedule.dat written (', n, ' outputs)'
    End Subroutine WriteOutputSchedule
!---------------------------------------------------------
!   Lens planes (item I2). FI holds the CIC overdensity delta on the
!   NGRID^3 mesh (after DENSIT; grid node i sits at (i-1)*Box/NGRID).
!   The mesh is projected along each axis in LensSlab slabs of depth
!   Box/LensSlab; every plane is the mean of (1+delta) over the slab
!   column, block-averaged to LensGrid x LensGrid (Sigma/Sigma_bar).
!   File LensPlanes.<step>.<box>.bin: raw little-endian float32, no
!   record markers, planes in the order X1..XN, Y1..YN, Z1..ZN, each
!   (LensGrid, LensGrid) column-major: X-planes indexed (y,z), Y-planes
!   (x,z), Z-planes (x,y). A .meta sidecar records the geometry.
!   Accumulation is real*8; threads own disjoint rows of the planes.
!---------------------------------------------------------
    Subroutine LensPlanes(moment)
        Integer*4, intent(in) :: moment
        Integer*4 :: Nl, Ns, f, nc, ix, iy, iz, s, ip, jp, iax, i, j
        Integer*8 :: c0, c1, cr
        Real*8, allocatable :: acc(:,:,:)
        Real*4, allocatable :: plane(:,:)
        Real*8 :: wnorm, d, pmean, amean(3)
        Character*120 :: Name
        Character*1   :: axname(3) = (/'X', 'Y', 'Z'/)

        If (LensSlab <= 0) Return
        Nl = LensGrid
        If (Nl <= 0) Nl = NGRID
        Ns = LensSlab
        If (mod(NGRID, Ns) /= 0 .or. mod(NGRID, Nl) /= 0) Stop 'LensPlanes: slabs and plane grid must divide NGRID'
        f = NGRID/Nl                                  ! block-average factor
        nc = NGRID/Ns                                 ! cells per slab along the projection axis
        wnorm = 1.d0/(dble(nc)*dble(f)*dble(f))       ! cells per plane pixel
        Call system_clock(c0, cr)

        write (Name, '(2(a,i4.4),a)') 'LensPlanes.', moment, '.', Nrealization, '.bin'
        open (29, file=TRIM(Name), status='replace', access='stream', form='unformatted', convert='little_endian')
        allocate (acc(Nl, Nl, Ns), plane(Nl, Nl))
        Do iax = 1, 3
            acc = 0.d0
            Select Case (iax)
            Case (1)        ! project along x: plane (y,z); thread owns plane row jp (from iz)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(jp,iz,iy,ip,s,ix,d) SCHEDULE(STATIC)
                Do jp = 1, Nl
                    Do iz = (jp-1)*f + 1, jp*f
                        Do iy = 1, NGRID
                            ip = (iy-1)/f + 1
                            Do s = 1, Ns
                                d = 0.d0
                                Do ix = (s-1)*nc + 1, s*nc
                                    d = d + FI(ix, iy, iz)
                                End Do
                                acc(ip, jp, s) = acc(ip, jp, s) + d
                            End Do
                        End Do
                    End Do
                End Do
            Case (2)        ! project along y: plane (x,z); thread owns plane row jp (from iz)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(jp,iz,iy,ip,s,ix) SCHEDULE(STATIC)
                Do jp = 1, Nl
                    Do iz = (jp-1)*f + 1, jp*f
                        Do iy = 1, NGRID
                            s = (iy-1)/nc + 1
                            Do ix = 1, NGRID
                                ip = (ix-1)/f + 1
                                acc(ip, jp, s) = acc(ip, jp, s) + FI(ix, iy, iz)
                            End Do
                        End Do
                    End Do
                End Do
            Case (3)        ! project along z: plane (x,y); thread owns plane row jp (from iy)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(jp,iy,iz,ip,s,ix) SCHEDULE(STATIC)
                Do jp = 1, Nl
                    Do iy = (jp-1)*f + 1, jp*f
                        Do iz = 1, NGRID
                            s = (iz-1)/nc + 1
                            Do ix = 1, NGRID
                                ip = (ix-1)/f + 1
                                acc(ip, jp, s) = acc(ip, jp, s) + FI(ix, iy, iz)
                            End Do
                        End Do
                    End Do
                End Do
            End Select
            amean(iax) = 0.d0
            Do s = 1, Ns                               ! mean(1+delta) = 1 + mean(delta)
                Do j = 1, Nl
                    Do i = 1, Nl
                        plane(i, j) = real(1.d0 + acc(i, j, s)*wnorm)
                    End Do
                End Do
                write (29) plane
                pmean = sum(dble(plane))/(dble(Nl)*dble(Nl))
                amean(iax) = amean(iax) + pmean/Ns
                write (*, '(a,a1,i3,a,f12.8,2(a,es11.3))') '   lens plane ', axname(iax), s, &
                    ': mean(1+delta) =', pmean, '  min =', minval(plane), '  max =', maxval(plane)
            End Do
        End Do
        close (29)
        deallocate (acc, plane)
        Call system_clock(c1)
        write (*, '(a,3f12.8,a,f8.2,a)') '   lens planes: axis means (must be 1) =', amean, &
            '   time =', dble(c1 - c0)/dble(cr), ' s'

        write (Name, '(2(a,i4.4),a)') 'LensPlanes.', moment, '.', Nrealization, '.meta'
        open (30, file=TRIM(Name), status='replace')
        write (30, '(a)') 'product = LensPlanes'
        write (30, '(a)') 'format = raw little-endian float32, no record markers; planes in order X1..XN, Y1..YN, Z1..ZN (projection axis, slab)'
        write (30, '(a)') 'plane_layout = (N_lens, N_lens) Fortran column-major; X-planes indexed (y,z), Y-planes (x,z), Z-planes (x,y)'
        write (30, '(a)') 'pixel = pixel (i,j) covers [(i-1), i) * Box/N_lens along the first plane axis; slab s covers [(s-1), s) * depth along the projection axis'
        write (30, '(a)') 'normalisation = mean of (1+delta) over the slab column and the pixel block (Sigma/Sigma_bar); delta from CIC on the Ngrid mesh, not deconvolved'
        write (30, '(a)') 'assignment = CIC'
        write (30, '(a)') 'interlaced = false'
        write (30, '(a,i6)') 'step = ', moment
        write (30, '(a,i6)') 'realization = ', Nrealization
        write (30, '(a,es16.8)') 'a = ', AEXPN
        write (30, '(a,es16.8)') 'z = ', 1./AEXPN - 1.
        write (30, '(a,es16.8)') 'a_x = ', AEXPN
        write (30, '(a,es16.8,a)') 'a_v = ', AEXPN - ASTEP/2., '   (velocity epoch of this output; not used by the planes)'
        write (30, '(a,f12.4)') 'Box = ', Box
        write (30, '(a)') 'units = Mpc/h'
        write (30, '(a,i6)') 'N_grid = ', NGRID
        write (30, '(a,i6)') 'N_row = ', NROW
        write (30, '(a,i20)') 'N_p = ', Nparticles
        write (30, '(a,i6)') 'N_slab = ', Ns
        write (30, '(a,f12.4)') 'depth = ', Box/Ns
        write (30, '(a,i6)') 'N_lens = ', Nl
        write (30, '(a,i6)') 'N_planes = ', 3*Ns
        write (30, '(a)') 'axis_map = plane axes follow the particle x, y, z of the halo catalogue'
        write (30, '(a,f10.6)') 'Omega_m = ', Om
        write (30, '(a,f10.6)') 'Omega_L = ', OmL
        write (30, '(a,f10.6)') 'hubble = ', hubble
        write (30, '(a,i3)') 'MG_flag = ', MG_flag
        write (30, '(a,i3)') 'MG_model = ', MG_model
        If (MG_flag == 1) Then
            write (30, '(a,es14.6)') 'fR0 = ', fR0
        Else
            write (30, '(a)') 'fR0 = none'
        End If
        write (30, '(a,a)') 'header = ', TRIM(HEADER)
        close (30)
    End Subroutine LensPlanes
!---------------------------------------------------------
!   Linear growth factor D(a) of the flat matter+Lambda background
!   (Om, OmL from Setup.dat), normalised to D(1) = 1:
!   D ~ E(a) int_0^a da'/(a' E(a'))^3, Simpson in ln a.
!---------------------------------------------------------
    Real*8 Function GrowthFactor(a)
        Real*4, intent(in) :: a
        GrowthFactor = GrowthRaw(dble(a))/GrowthRaw(1.d0)
    End Function GrowthFactor

    Real*8 Function GrowthRaw(a)
        Real*8, intent(in) :: a
        Integer*4, parameter :: Nint = 4000
        Real*8 :: x0, h, x, s, f, e
        Integer*4 :: i
        x0 = log(1.d-8)
        h = (log(a) - x0)/Nint
        s = 0.d0
        Do i = 0, Nint
            x = x0 + h*i
            f = exp(2.5d0*x)/(dble(Om) + dble(OmL)*exp(3.d0*x))**1.5d0   ! a'^{5/2}/(Om + OmL a'^3)^{3/2}
            If (i == 0 .or. i == Nint) Then
                s = s + f
            Else If (mod(i, 2) == 1) Then
                s = s + 4.d0*f
            Else
                s = s + 2.d0*f
            End If
        End Do
        e = sqrt(dble(Om)/a**3 + dble(OmL))
        GrowthRaw = e*s*h/3.d0
    End Function GrowthRaw
!---------------------------------------------------------
!   Forward 3-D real FFT of G(N,N,N) in place, replicating POWERfft5
!   exactly (three FFTPACK5 rfft1f passes along x, y, then z after a
!   transposition of indices 1 and 3; indices 1 and N of each pass are
!   multiplied by sqrt(2)). Afterwards the packed coefficient of mode
!   (px,py,pz) sits at G(pz,py,px): per dimension p=1 is sqrt(2)/N
!   times the sum, p=2m holds (2/N) sum x cos(2 pi m (j-1)/N), p=2m+1
!   holds (2/N) sum x sin(...), p=N the Nyquist term times sqrt(2).
!---------------------------------------------------------
    Subroutine FFT3Dforward(G, N)
        Integer*4, intent(in) :: N
        Real*4, intent(inout) :: G(N, N, N)
        Real*8 :: r(Nlenmax)
        Real*4 :: aa
        Integer*4 :: i, j, k, ier, lensav, lenwrk, inc, lenr
        If (N > Nlenmax) Stop 'FFT3Dforward: grid larger than the FFT work arrays'
        lensav = N + int(log(real(N, kind=4))/log(2.0E+00)) + 4
        lenwrk = N
        inc = 1
        lenr = N
        call rfft1i(N, wsave, lensav, ier)
!$OMP PARALLEL DO DEFAULT(SHARED) COPYIN(wsave,work) PRIVATE(k,j,i,r,ier)
        Do k = 1, N                ! ------ x-direction
        Do j = 1, N
            Do i = 1, N
                r(i) = G(i, j, k)
            End Do
            call rfft1f(N, inc, r, lenr, wsave, lensav, work, lenwrk, ier)
            Do i = 1, N
                G(i, j, k) = r(i)
            End Do
        End Do
        End Do
        Do k = 1, N
        Do j = 1, N
            G(1, j, k) = G(1, j, k)*sqrt(2.)
            G(N, j, k) = G(N, j, k)*sqrt(2.)
        End Do
        End Do
!$OMP PARALLEL DO DEFAULT(SHARED) COPYIN(wsave,work) PRIVATE(k,j,i,r,ier)
        Do k = 1, N                ! ------ y-direction
        Do i = 1, N
            Do j = 1, N
                r(j) = G(i, j, k)
            End Do
            call rfft1f(N, inc, r, lenr, wsave, lensav, work, lenwrk, ier)
            Do j = 1, N
                G(i, j, k) = r(j)
            End Do
        End Do
        End Do
        Do k = 1, N
        Do i = 1, N
            G(i, 1, k) = G(i, 1, k)*sqrt(2.)
            G(i, N, k) = G(i, N, k)*sqrt(2.)
        End Do
        End Do
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k,j,i,aa)
        Do j = 1, N                ! ------ transpose indices 1 and 3
        Do k = 1, N - 1
            Do i = k + 1, N
                aa = G(i, j, k)
                G(i, j, k) = G(k, j, i)
                G(k, j, i) = aa
            End Do
        End Do
        End Do
!$OMP PARALLEL DO DEFAULT(SHARED) COPYIN(wsave,work) PRIVATE(k,j,i,r,ier)
        Do j = 1, N                ! ------ z-direction (now the first index)
        Do i = 1, N
            Do k = 1, N
                r(k) = G(k, i, j)
            End Do
            call rfft1f(N, inc, r, lenr, wsave, lensav, work, lenwrk, ier)
            Do k = 1, N
                G(k, i, j) = r(k)
            End Do
        End Do
        End Do
        Do j = 1, N
        Do i = 1, N
            G(1, i, j) = G(1, i, j)*sqrt(2.)
            G(N, i, j) = G(N, i, j)*sqrt(2.)
        End Do
        End Do
    End Subroutine FFT3Dforward
!---------------------------------------------------------
!   Low-k mode dump (item I7/I8 file writer). G holds a field in the
!   packed layout produced by FFT3Dforward / POWERfft5 (before the
!   squaring). Writes delta_k = N^-3 sum_x f(x) exp(-i k.x) for the cube
!   |ix|,|iy| <= n, 0 <= iz <= n, n = floor(kc Box/2pi), k = 2pi/Box
!   (ix,iy,iz), as complex(4) C(0:n,-n:n,-n:n) (first index iz), k=0
!   set to zero, plus a .meta sidecar. File <product>.<step>.<box>.bin.
!---------------------------------------------------------
    Subroutine DumpModes(G, N, product, kc, moment, a_v)
        Integer*4, intent(in) :: N, moment
        Real*4, intent(in) :: G(N, N, N), kc, a_v
        Character(len=*), intent(in) :: product
        Complex*8, allocatable :: C(:, :, :)
        Integer*4 :: nm, mx, my, mz, ix, iy, iz, px(2), py(2), pz(2), nx, ny, nz, i, j, k
        Complex*8 :: wx(2), wy(2), wz(2), val
        Character*120 :: Name
        nm = int(dble(kc)*dble(Box)/6.283185307179586d0)
        If (nm > N/2 - 1) Then
            write (*, '(a,a,a,i5,a,i5)') ' DumpModes ', product, ': n=', nm, ' limited to N/2-1=', N/2 - 1
            nm = N/2 - 1
        End If
        allocate (C(0:nm, -nm:nm, -nm:nm))
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(mx,my,mz,px,py,pz,wx,wy,wz,nx,ny,nz,i,j,k,val)
        Do mx = -nm, nm
            Call PackedBasis(mx, N, px, wx, nx)
            Do my = -nm, nm
                Call PackedBasis(my, N, py, wy, ny)
                Do mz = 0, nm
                    Call PackedBasis(mz, N, pz, wz, nz)
                    val = (0., 0.)
                    Do i = 1, nx
                    Do j = 1, ny
                    Do k = 1, nz
                        val = val + wx(i)*wy(j)*wz(k)*G(pz(k), py(j), px(i))
                    End Do
                    End Do
                    End Do
                    C(mz, my, mx) = val
                End Do
            End Do
        End Do
        C(0, 0, 0) = (0., 0.)
        write (Name, '(a,2(a,i4.4),a)') TRIM(product), '.', moment, '.', Nrealization, '.bin'
        open (29, file=TRIM(Name), status='replace', access='stream', form='unformatted', convert='little_endian')
        write (29) C
        close (29)
        write (Name, '(a,2(a,i4.4),a)') TRIM(product), '.', moment, '.', Nrealization, '.meta'
        open (30, file=TRIM(Name), status='replace')
        write (30, '(a,a)') 'product = ', TRIM(product)
        write (30, '(a)') 'format = raw little-endian complex64 (re, im float32), no record markers'
        write (30, '(a)') 'layout = Fortran C(0:nm,-nm:nm,-nm:nm) with the first index iz; numpy: fromfile(dtype=complex64).reshape(2n+1,2n+1,nm+1) -> [ix+nm, iy+nm, iz]'
        write (30, '(a)') 'fourier_norm = delta_k = N^-3 sum_x f(x) exp(-i k.x), k = (2 pi/Box)(ix,iy,iz), P(k) = Box^3 <|delta_k|^2>; grid node i at x = (i-1) Box/N'
        write (30, '(a)') 'modes = cube |ix|,|iy| <= nm, 0 <= iz <= nm (Hermitian-redundant iz = 0 modes kept); k = 0 set to zero'
        write (30, '(a,i6)') 'n = ', nm
        write (30, '(a,f10.5)') 'k_c = ', kc
        write (30, '(a,i6)') 'N_dep = ', N
        write (30, '(a,i6)') 'N_grid = ', NGRID
        write (30, '(a)') 'assignment = CIC'
        write (30, '(a)') 'interlaced = false'
        write (30, '(a)') 'axis_map = ix, iy, iz follow the particle x, y, z of the halo catalogue'
        write (30, '(a,i6)') 'step = ', moment
        write (30, '(a,i6)') 'realization = ', Nrealization
        write (30, '(a,es16.8)') 'a = ', AEXPN
        write (30, '(a,es16.8)') 'z = ', 1./AEXPN - 1.
        write (30, '(a,es16.8)') 'a_x = ', AEXPN
        write (30, '(a,es16.8)') 'a_v = ', a_v
        write (30, '(a,es16.8)') 'D_code(a) = ', GrowthFactor(AEXPN)
        write (30, '(a,f12.4)') 'Box = ', Box
        write (30, '(a)') 'units = Mpc/h; momentum products are (1+delta) v with v the peculiar velocity in km/s at a_v'
        write (30, '(a,i6)') 'N_row = ', NROW
        write (30, '(a,i20)') 'N_p = ', Nparticles
        write (30, '(a,f10.6)') 'Omega_m = ', Om
        write (30, '(a,f10.6)') 'Omega_L = ', OmL
        write (30, '(a,f10.6)') 'hubble = ', hubble
        write (30, '(a,i3)') 'MG_flag = ', MG_flag
        write (30, '(a,i3)') 'MG_model = ', MG_model
        If (MG_flag == 1) Then
            write (30, '(a,es14.6)') 'fR0 = ', fR0
        Else
            write (30, '(a)') 'fR0 = none'
        End If
        write (30, '(a,a)') 'header = ', TRIM(HEADER)
        close (30)
        deallocate (C)
        write (*, '(a,a,a,i4,a,i10,a)') '   modes ', TRIM(product), ': n=', nm, '  (', (2*nm + 1)**2*(nm + 1), ' complex values)'
    End Subroutine DumpModes
!---------------------------------------------------------
!   Packed-basis entries of one dimension for mode m: indices p(1:np)
!   into the packed array and the complex weights w such that the 1-D
!   DFT coefficient X_m = sum_i w(i) r(p(i)) given FFT3Dforward's
!   normalisation (r(1), r(N) carry the sqrt(2) factor).
!---------------------------------------------------------
    Subroutine PackedBasis(m, N, p, w, np)
        Integer*4, intent(in) :: m, N
        Integer*4, intent(out) :: p(2), np
        Complex*8, intent(out) :: w(2)
        If (m == 0) Then
            np = 1
            p(1) = 1
            w(1) = cmplx(1./sqrt(2.), 0.)
        Else If (abs(m) == N/2) Then
            np = 1
            p(1) = N
            w(1) = cmplx(1./sqrt(2.), 0.)
        Else
            np = 2
            p(1) = 2*abs(m)
            p(2) = 2*abs(m) + 1
            w(1) = cmplx(0.5, 0.)
            w(2) = cmplx(0., -0.5*sign(1, m))       ! X_m = (r_cos - i r_sin)/2; X_{-m} = conj
        End If
    End Subroutine PackedBasis
!---------------------------------------------------------
!   Momentum deposit (item I8): CIC of w_p = (Nd^3/N_p) v_i,p onto
!   G(Nd,Nd,Nd), so that G = (1+delta) v_i with v the peculiar velocity
!   in km/s at the velocity epoch a_v (v = VX (100 Box/NGRID)/a_v).
!---------------------------------------------------------
    Subroutine MomentumDeposit(G, Nd, icomp, a_v)
        Integer*4, intent(in) :: Nd, icomp
        Real*4, intent(inout) :: G(Nd, Nd, Nd)
        Real*4, intent(in) :: a_v
        Integer*8 :: IN
        Integer*4 :: I, J, K, I1, J1, K1
        Real*4 :: X, Y, Z, D1, D2, D3, T1, T2, T3, T2W, D2W, W, vfac, sc, XN
        vfac = 100.*Box/NGRID/a_v*(dble(Nd)**3/dble(Nparticles))
        sc = real(Nd)/real(NGRID)
        XN = float(Nd) + 1. - 1.E-7
        G = 0.
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(IN,X,Y,Z,D1,D2,D3,T1,T2,T3,T2W,D2W,W,I,J,K,I1,J1,K1)
        Do IN = 1, Nparticles
            X = (XPAR(IN) - 1.)*sc + 1.
            Y = (YPAR(IN) - 1.)*sc + 1.
            Z = (ZPAR(IN) - 1.)*sc + 1.
            If (X >= XN) X = X - Nd
            If (Y >= XN) Y = Y - Nd
            If (Z >= XN) Z = Z - Nd
            Select Case (icomp)
            Case (1)
                W = VX(IN)*vfac
            Case (2)
                W = VY(IN)*vfac
            Case Default
                W = VZ(IN)*vfac
            End Select
            I = INT(X)
            J = INT(Y)
            K = INT(Z)
            D1 = X - FLOAT(I)
            D2 = Y - FLOAT(J)
            D3 = Z - FLOAT(K)
            T1 = 1. - D1
            T2 = 1. - D2
            T3 = 1. - D3
            T2W = T2*W
            D2W = D2*W
            I1 = I + 1
            IF (I1 .GT. Nd) I1 = 1
            J1 = J + 1
            IF (J1 .GT. Nd) J1 = 1
            K1 = K + 1
            IF (K1 .GT. Nd) K1 = 1
!$OMP ATOMIC
            G(I, J, K) = G(I, J, K) + T3*T1*T2W
!$OMP ATOMIC
            G(I1, J, K) = G(I1, J, K) + T3*D1*T2W
!$OMP ATOMIC
            G(I, J1, K) = G(I, J1, K) + T3*T1*D2W
!$OMP ATOMIC
            G(I1, J1, K) = G(I1, J1, K) + T3*D1*D2W
!$OMP ATOMIC
            G(I, J, K1) = G(I, J, K1) + D3*T1*T2W
!$OMP ATOMIC
            G(I1, J, K1) = G(I1, J, K1) + D3*D1*T2W
!$OMP ATOMIC
            G(I, J1, K1) = G(I, J1, K1) + D3*T1*D2W
!$OMP ATOMIC
            G(I1, J1, K1) = G(I1, J1, K1) + D3*D1*D2W
        End Do
    End Subroutine MomentumDeposit
!---------------------------------------------------------
!   Extra mode products at a class-1 output (item I8): the density
!   modes when the power spectrum was not computed (otherwise POWERfft5
!   dumps them from its own transform), and the three momentum
!   components on a transient Nd^3 grid.
!---------------------------------------------------------
    Subroutine ExtraModes(moment, mDENSIT, iPowerDone)
        Integer*4, intent(in) :: moment, iPowerDone
        Integer*4, intent(inout) :: mDENSIT
        Real*4, allocatable :: G(:, :, :)
        Integer*4 :: Nd, ic
        Integer*8 :: c0, c1, cr
        Real*4 :: a_v
        Character*5 :: pname(3) = (/'mom_x', 'mom_y', 'mom_z'/)
        If (iExtra /= 1) Return
        a_v = AEXPN - ASTEP/2.
        Call system_clock(c0, cr)
        If (iPowerDone == 0) Then                    ! density modes from FI
            If (mDENSIT == 1) Then
                Call DENSIT
                mDENSIT = 0
            End If
            Call FFT3Dforward(FI, NGRID)
            Call DumpModes(FI, NGRID, 'deltam', KcDelta, moment, a_v)
            mDENSIT = 1
        End If
        Nd = NdepMom
        If (Nd <= 0) Nd = NGRID/2
        allocate (G(Nd, Nd, Nd))
        Do ic = 1, 3
            Call MomentumDeposit(G, Nd, ic, a_v)
            Call FFT3Dforward(G, Nd)
            Call DumpModes(G, Nd, pname(ic), KcMom, moment, a_v)
        End Do
        deallocate (G)
        Call system_clock(c1)
        write (*, '(a,i5,a,f8.2,a)') '   extra modes: momentum on N_dep=', Nd, '   time =', dble(c1 - c0)/dble(cr), ' s'
    End Subroutine ExtraModes
!---------------------------------------------------------
!   f_R at halo centres (item I10, f(R) runs). FI2 holds the converged
!   scalar field u = |f_R|^(1/(n+1)) normalised by its background value,
!   solved at the pre-step epoch a_old = AEXPN - ASTEP, so that
!   f_R(x) = fR_bg(a_old) * FI2(x)^(n+1), fR_bg(a) = fR0 (R0/R(a))^(n+1),
!   R(a) = 3 (Om/a^3 + 4 OmL). Grid node i sits at code coordinate X = i.
!   fRHaloOpen is called by BDM before WriteFiles, fRHaloRow once per
!   catalogue row (same order as CatshortV), fRHaloClose afterwards.
!   File fR_halo.<step>.<box>.bin: float32 pairs (fR_centre, fR_mean).
!---------------------------------------------------------
    Subroutine fRHaloOpen(moment)
        Integer*4, intent(in) :: moment
        Real*8 :: R0, Rold
        Character*120 :: Name
        a_old_fR = dble(AEXPN) - dble(ASTEP)
        R0 = 3.d0*(dble(Om) + 4.d0*dble(OmL))
        Rold = 3.d0*(dble(Om)/a_old_fR**3 + 4.d0*dble(OmL))
        fRbg_old = fR0*(R0/Rold)**(fr_n + 1)
        write (Name, '(2(a,i4.4),a)') 'fR_halo.', moment, '.', Nrealization, '.bin'
        open (31, file=TRIM(Name), status='replace', access='stream', form='unformatted', convert='little_endian')
        nfRrows = 0
        fRmoment = moment
        fRHaloActive = .true.
    End Subroutine fRHaloOpen

    Subroutine fRHaloRow(x, y, z, rvir)
        Real*4, intent(in) :: x, y, z, rvir          ! centre [Mpc/h] in [0,Box), aperture radius [Mpc/h]
        Real*8 :: X1, Y1, Z1, D1, D2, D3, T1, T2, T3, uc, um, r, dx, dy, dz, XN
        Integer*4 :: I, J, K, I1, J1, K1, ii, jj, kk, ip, jp, kp, nr, ncell
        If (.not. fRHaloActive) Return
        XN = dble(NGRID)
        X1 = modulo(dble(x)*XN/dble(Box), XN) + 1.d0
        Y1 = modulo(dble(y)*XN/dble(Box), XN) + 1.d0
        Z1 = modulo(dble(z)*XN/dble(Box), XN) + 1.d0
        I = int(X1); J = int(Y1); K = int(Z1)
        If (I > NGRID) I = NGRID
        If (J > NGRID) J = NGRID
        If (K > NGRID) K = NGRID
        D1 = X1 - I; D2 = Y1 - J; D3 = Z1 - K
        T1 = 1.d0 - D1; T2 = 1.d0 - D2; T3 = 1.d0 - D3
        I1 = I + 1; If (I1 > NGRID) I1 = 1
        J1 = J + 1; If (J1 > NGRID) J1 = 1
        K1 = K + 1; If (K1 > NGRID) K1 = 1
        uc = T3*(T1*(T2*dble(FI2(I,J,K))**(fr_n+1) + D2*dble(FI2(I,J1,K))**(fr_n+1)) &
                + D1*(T2*dble(FI2(I1,J,K))**(fr_n+1) + D2*dble(FI2(I1,J1,K))**(fr_n+1))) &
           + D3*(T1*(T2*dble(FI2(I,J,K1))**(fr_n+1) + D2*dble(FI2(I,J1,K1))**(fr_n+1)) &
                + D1*(T2*dble(FI2(I1,J,K1))**(fr_n+1) + D2*dble(FI2(I1,J1,K1))**(fr_n+1)))
        ! mean of f_R over the grid cells whose centres lie within rvir of the centre
        r = dble(rvir)*XN/dble(Box)
        nr = int(r) + 1
        um = 0.d0
        ncell = 0
        If (r >= 0.5d0) Then
            Do kk = K - nr, K + nr
                dz = dble(kk) - Z1
                kp = modulo(kk - 1, NGRID) + 1
                Do jj = J - nr, J + nr
                    dy = dble(jj) - Y1
                    jp = modulo(jj - 1, NGRID) + 1
                    Do ii = I - nr, I + nr
                        dx = dble(ii) - X1
                        If (dx*dx + dy*dy + dz*dz <= r*r) Then
                            ip = modulo(ii - 1, NGRID) + 1
                            um = um + dble(FI2(ip, jp, kp))**(fr_n+1)
                            ncell = ncell + 1
                        End If
                    End Do
                End Do
            End Do
        End If
        If (ncell > 0) Then
            um = um/ncell
        Else
            um = uc
        End If
        write (31) real(fRbg_old*uc, 4), real(fRbg_old*um, 4)
        nfRrows = nfRrows + 1
    End Subroutine fRHaloRow

    Subroutine fRHaloClose(nhalo)
        Integer*4, intent(in) :: nhalo
        Character*120 :: Name
        If (.not. fRHaloActive) Return
        close (31)
        fRHaloActive = .false.
        If (nhalo /= nfRrows) write (*, '(a,2i10)') ' WARNING fR_halo: rows written /= catalogue rows', nfRrows, nhalo
        write (Name, '(2(a,i4.4),a)') 'fR_halo.', fRmoment, '.', Nrealization, '.meta'
        open (30, file=TRIM(Name), status='replace')
        write (30, '(a)') 'product = fR_halo'
        write (30, '(a)') 'format = raw little-endian float32 pairs (fR_centre, fR_mean), one pair per catalogue row, same order as CatshortV'
        write (Name, '(2(a,i4.4),a)') 'CATALOGS/CatshortV.', fRmoment, '.', Nrealization, '.DAT'
        write (30, '(a,a)') 'catalogue = ', TRIM(Name)
        write (30, '(a,i10)') 'rows = ', nfRrows
        write (30, '(a)') 'fR_centre = f_R interpolated (CIC) at the BDM halo centre'
        write (30, '(a)') 'fR_mean = mean of f_R over the grid cells whose centres lie within the catalogue radius (aperture Rvir) of the centre; equals fR_centre when the radius is below half a cell'
        write (30, '(a)') 'conversion = f_R(x) = fR_bg(a_old) * u(x)^(n+1), u = the code field FI2 (scalar field normalised by its background value), solved at the pre-step epoch a_old'
        write (30, '(a)') 'fR_bg(a) = fR0 * (R0/R(a))^(n+1), R(a) = 3 (Omega_m/a^3 + 4 Omega_L), R0 = R(1)'
        write (30, '(a,es16.8)') 'fR0 = ', fR0
        write (30, '(a,i3)') 'n = ', fr_n
        write (30, '(a,es16.8)') 'a_old = ', a_old_fR
        write (30, '(a,es16.8)') 'fR_bg(a_old) = ', fRbg_old
        write (30, '(a,es16.8)') 'a = ', AEXPN
        write (30, '(a,es16.8)') 'z = ', 1./AEXPN - 1.
        write (30, '(a,i6)') 'step = ', fRmoment
        write (30, '(a,i6)') 'realization = ', Nrealization
        write (30, '(a,f12.4)') 'Box = ', Box
        write (30, '(a,i6)') 'N_grid = ', NGRID
        write (30, '(a,f10.6)') 'Omega_m = ', Om
        write (30, '(a,f10.6)') 'Omega_L = ', OmL
        write (30, '(a,f10.6)') 'hubble = ', hubble
        write (30, '(a,a)') 'header = ', TRIM(HEADER)
        close (30)
        write (*, '(a,i8,a)') '   fR_halo: ', nfRrows, ' rows written'
    End Subroutine fRHaloClose
End Module ExtraOutputs
