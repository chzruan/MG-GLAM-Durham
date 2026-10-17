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
    implicit none
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
        If (MG_flag == 1) Then
            write (30, '(a,es14.6)') 'fR0 = ', fR0
        Else
            write (30, '(a)') 'fR0 = none'
        End If
        write (30, '(a,a)') 'header = ', TRIM(HEADER)
        close (30)
    End Subroutine LensPlanes
End Module ExtraOutputs
