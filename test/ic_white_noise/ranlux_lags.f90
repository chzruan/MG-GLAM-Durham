!------------------------------------------------------------------------------
!  ranlux_lags: statistics of the white noise drawn by MG-GLAM's own random
!  number routines (PMP2mod_random.f90: rluxgo, ranlux, RANDd, GAUSS3), with
!  the seeding of PMP2start.f90 SetRandomN (one RANLUX stream per plane k,
!  rluxgo(lux, Ns, 0, 0), Ns advanced by RANDd) and the draw order of SPECTR
!  (plane k, row j, column i; one GAUSS3 per (i, j, k)).  No simulation code
!  is changed or copied: the routines are called from the module.
!
!  mode 1, the Gaussian cube.  stdin:
!      1
!      Nseed NROW Nplanes Nrows Nkeep lux   (lux < 0: the code's own choice, RanluxLevel)
!      cubefile   ('none' to write no cube)
!      outfile
!    Draws planes 1..Nplanes, rows 1..Nrows, all NROW columns of each row.
!    Writes to outfile: the luxury level the generator actually uses, the
!    first values, the lag correlation along i (all columns), along j and k
!    (the kept cube i <= Nkeep), and to outfile.vprof the variance profile
!    along x of the white field synthesised row by row with rfft1b (the basis
!    of VECTOR).  Writes the kept cube (real*4, Nkeep x Nrows x Nplanes) to
!    cubefile and the basis (real*8, NROW x NROW) to cubefile.basis if asked.
!
!  mode 2, the mechanism.  stdin:
!      2
!      Nseed Npairs lux
!      outfile
!    One stream; Npairs calls of ranlux(RanNum, 2) as in GAUSS3; records the
!    uniform pairs and the Box-Muller acceptance; writes the residual of the
!    subtract-with-borrow relation between pairs p, p+7, p+12, the
!    correlations of the uniforms and of the acceptance by pair lag, the
!    correlation of the delivered Gaussians by Gaussian lag, and the
!    correlation of pairs p and p+12 split by the acceptance of pair p+7.
!------------------------------------------------------------------------------
program ranlux_lags
  implicit none
  integer*4 :: mode
  read (*, *) mode
  if (mode == 1) then
     call mode_cube
  else if (mode == 2) then
     call mode_pairs
  else
     stop 'mode must be 1 or 2'
  end if
end program ranlux_lags

!------------------------------------------------------------------------------
subroutine mode_cube
  use Random
  use FFT5
  INCLUDE 'luxuryp.h'   ! implicit typing: the include file relies on it
  integer*4, parameter :: LMAX = 40
  integer*4 :: Nseed, NROW, Nplanes, Nrows, Nkeep, lux, Ns
  integer*4 :: i, j, k, m, L, ier, lensav, lenwrk, luxlev1, nskip1, iFlag, nfirst
  integer*8 :: n1, nl(LMAX), nc
  real*4 :: gSet, TS, dummy
  real*8 :: s1, s2, sl(LMAX), first(8), cmean, cvar, sj(LMAX), sk(LMAX), chk1
  integer*8 :: nj(LMAX), nk(LMAX)
  character(len=512) :: cubefile, outfile
  real*4, allocatable :: cube(:, :, :), row(:)
  real*8, allocatable :: wsave(:), work(:), r(:), vprof(:), basis(:, :)

  read (*, *) Nseed, NROW, Nplanes, Nrows, Nkeep, lux
  read (*, '(a)') cubefile
  read (*, '(a)') outfile
  if (Nplanes > NROWmax) stop 'Nplanes > NROWmax'
  if (Nkeep > NROW) stop 'Nkeep > NROW'
  if (lux < 0) lux = RanluxLevel()   ! the code's own choice (GLAM_RANLUX_LUXURY, default 3)

  !---- SetRandomN of PMP2start.f90, planes 1..Nplanes
  gSet = 0.
  Ns = Nseed
  do k = 1, Nplanes
     call rluxgo(lux, Ns, 0, 0)
     if (k == 1) then
        luxlev1 = luxlev
        nskip1 = nskip
     end if
     i24A(k) = i24
     j24A(k) = j24
     in24A(k) = in24
     kountA(k) = kount
     carryA(k) = carry
     gSetA(k) = gSet
     iFlagA(k) = 0     ! SetRandomN stores an undefined local iFlag; 0 reproduces the runs
     do m = 1, 24
        SeedsPage(m, k) = seeds(m)
     end do
     dummy = RANDd(Ns)
  end do

  !---- SPECTR draw order: plane k, row j, column i
  allocate (cube(Nkeep, Nrows, Nplanes), row(NROW))
  lensav = NROW + int(log(real(NROW, kind=4))/log(2.0E+00)) + 4
  lenwrk = NROW
  allocate (wsave(lensav), work(lenwrk), r(NROW), vprof(NROW))
  call rfft1i(NROW, wsave, lensav, ier)
  s1 = 0.d0; s2 = 0.d0; n1 = 0; sl = 0.d0; nl = 0; vprof = 0.d0; nfirst = 0; first = 0.d0; chk1 = 0.d0
  do k = 1, Nplanes
     i24 = i24A(k)
     j24 = j24A(k)
     in24 = in24A(k)
     kount = kountA(k)
     mkount = 0
     carry = carryA(k)
     iFlag = iFlagA(k)
     gSet = gSetA(k)
     do m = 1, 24
        seeds(m) = SeedsPage(m, k)
     end do
     do j = 1, Nrows
        do i = 1, NROW
           TS = GAUSS3(gSet, iFlag)
           row(i) = TS
        end do
        if (k == 1 .and. j == 1) then
           first(1:8) = row(1:8)
           chk1 = sum(dble(row))
        end if
        do i = 1, NROW
           s1 = s1 + row(i)
           s2 = s2 + dble(row(i))**2
        end do
        n1 = n1 + NROW
        do L = 1, LMAX
           do i = 1, NROW - L
              sl(L) = sl(L) + dble(row(i))*dble(row(i + L))
           end do
           nl(L) = nl(L) + (NROW - L)
        end do
        cube(1:Nkeep, j, k) = row(1:Nkeep)
        r(1:NROW) = row(1:NROW)
        call rfft1b(NROW, 1, r, NROW, wsave, lensav, work, lenwrk, ier)
        vprof = vprof + r*r
     end do
  end do
  vprof = vprof/dble(Nrows)/dble(Nplanes)

  !---- lag correlations along j and k from the kept cube
  cmean = sum(dble(cube))/(dble(Nkeep)*Nrows*Nplanes)
  cvar = sum((dble(cube) - cmean)**2)/(dble(Nkeep)*Nrows*Nplanes)
  sj = 0.d0; sk = 0.d0; nj = 0; nk = 0
  do L = 1, LMAX
     if (L < Nrows) then
        sj(L) = sum((dble(cube(:, 1:Nrows - L, :)) - cmean)*(dble(cube(:, 1 + L:Nrows, :)) - cmean))
        nj(L) = int(Nkeep, 8)*(Nrows - L)*Nplanes
     end if
     if (L < Nplanes) then
        sk(L) = sum((dble(cube(:, :, 1:Nplanes - L)) - cmean)*(dble(cube(:, :, 1 + L:Nplanes)) - cmean))
        nk(L) = int(Nkeep, 8)*Nrows*(Nplanes - L)
     end if
  end do

  !---- outputs
  open (21, file=trim(outfile), status='replace')
  write (21, '(a,i12,a,i6,a,i6,a,i6,a,i6)') 'seed ', Nseed, ' NROW ', NROW, ' Nplanes ', Nplanes, ' Nrows ', Nrows, ' Nkeep ', Nkeep
  write (21, '(a,i3,a,i3,a,i5)') 'lux_requested ', lux, ' luxlev_set ', luxlev1, ' nskip_used ', nskip1
  write (21, '(a,8es16.8)') 'first8 ', first
  write (21, '(a,es22.14)') 'sum_row1 ', chk1
  write (21, '(a,i14,a,es14.6,a,es14.6)') 'n ', n1, ' mean ', s1/n1, ' var ', s2/n1 - (s1/n1)**2
  write (21, '(a,es14.6,a,es14.6)') 'cube_mean ', cmean, ' cube_var ', cvar
  do L = 1, LMAX
     write (21, '(a,i3,3es14.5)') 'lag ', L, (sl(L)/nl(L) - (s1/n1)**2)/(s2/n1 - (s1/n1)**2), &
          merge(sj(L)/max(nj(L), 1_8)/cvar, 0.d0, nj(L) > 0), merge(sk(L)/max(nk(L), 1_8)/cvar, 0.d0, nk(L) > 0)
  end do
  close (21)
  open (22, file=trim(outfile)//'.vprof', status='replace')
  do i = 1, NROW
     write (22, '(i6,es18.10)') i, vprof(i)
  end do
  close (22)
  if (trim(cubefile) /= 'none') then
     open (23, file=trim(cubefile), form='unformatted', access='stream', status='replace')
     write (23) cube
     close (23)
     allocate (basis(NROW, NROW))
     do i = 1, NROW
        r(:) = 0.d0
        r(i) = 1.d0
        call rfft1b(NROW, 1, r, NROW, wsave, lensav, work, lenwrk, ier)
        basis(:, i) = r(:)
     end do
     open (24, file=trim(cubefile)//'.basis', form='unformatted', access='stream', status='replace')
     write (24) basis
     close (24)
  end if
  write (*, '(a,i12,a,i3,a,i3,a,i5,a,es12.4)') ' seed ', Nseed, '  lux requested ', lux, '  level set ', luxlev1, &
       '  nskip used ', nskip1, '  first value ', first(1)
end subroutine mode_cube

!------------------------------------------------------------------------------
subroutine mode_pairs
  use Random
  INCLUDE 'luxuryp.h'   ! implicit typing: the include file relies on it
  integer*4, parameter :: MMAX = 30, LMAX = 40
  integer*4 :: Nseed, Npairs, lux, p, m, L, nacc, ng
  integer*8 :: nres, n0, n1, nc(0:1, 6:8)
  real*4 :: RanNum(2), x1, x2, R, F
  real*4, allocatable :: u1(:), u2(:), g(:)
  integer*1, allocatable :: acc(:)
  real*8 :: res, cx(MMAX), ca(MMAX), cg(LMAX), mx, vx, ma, va, mg, vg, sc(0:1, 6:8), sg(0:1), s12, twom24r
  integer*8 :: ng12(0:1)
  character(len=512) :: outfile

  read (*, *) Nseed, Npairs, lux
  read (*, '(a)') outfile
  call rluxgo(lux, Nseed, 0, 0)
  allocate (u1(Npairs), u2(Npairs), acc(Npairs), g(2*Npairs))
  nacc = 0; ng = 0
  do p = 1, Npairs
     call ranlux(RanNum, 2)
     u1(p) = RanNum(1)
     u2(p) = RanNum(2)
     x1 = 2.*RanNum(1) - 1.
     x2 = 2.*RanNum(2) - 1.
     R = x1**2 + x2**2
     if (R .ge. 1.) then
        acc(p) = 0
     else
        acc(p) = 1
        nacc = nacc + 1
        F = sqrt(-2.*LOG(R)/R)
        g(ng + 1) = x2*F      ! GAUSS3 returns x2*F first, then x1*F (gSet)
        g(ng + 2) = x1*F
        ng = ng + 2
     end if
  end do

  !---- the subtract-with-borrow relation between pairs p, p+7, p+12 (positions n, n+14, n+24)
  twom24r = 1.d0/dble(2**24)
  nres = 0
  do p = 1, Npairs - 12
     res = dble(u1(p + 12)) - (dble(u1(p + 7)) - dble(u1(p)))
     res = res - nint(res)
     if (abs(res) <= 1.5d0*twom24r) nres = nres + 1
  end do

  !---- correlations by pair lag of x1 = 2 u1 - 1 and of the acceptance; of the Gaussians by Gaussian lag
  mx = sum(2.d0*dble(u1) - 1.d0)/Npairs
  vx = sum((2.d0*dble(u1) - 1.d0 - mx)**2)/Npairs
  ma = sum(dble(acc))/Npairs
  va = sum((dble(acc) - ma)**2)/Npairs
  do m = 1, MMAX
     cx(m) = sum((2.d0*dble(u1(1:Npairs - m)) - 1.d0 - mx)*(2.d0*dble(u1(1 + m:Npairs)) - 1.d0 - mx))/(Npairs - m)/vx
     ca(m) = sum((dble(acc(1:Npairs - m)) - ma)*(dble(acc(1 + m:Npairs)) - ma))/(Npairs - m)/va
  end do
  mg = sum(dble(g(1:ng)))/ng
  vg = sum((dble(g(1:ng)) - mg)**2)/ng
  do L = 1, LMAX
     cg(L) = sum((dble(g(1:ng - L)) - mg)*(dble(g(1 + L:ng)) - mg))/(ng - L)/vg
  end do

  !---- E[x1_p x1_{p+12} | acceptance of pair p+q], q = 6, 7, 8; and the Gaussians of accepted p, p+12 by acc(p+7)
  sc = 0.d0; nc = 0; sg = 0.d0; ng12 = 0; s12 = 0.d0
  do p = 1, Npairs - 12
     x1 = 2.*u1(p) - 1.
     x2 = 2.*u1(p + 12) - 1.
     s12 = s12 + dble(x1)*dble(x2)
     do m = 6, 8
        sc(acc(p + m), m) = sc(acc(p + m), m) + dble(x1)*dble(x2)
        nc(acc(p + m), m) = nc(acc(p + m), m) + 1
     end do
     if (acc(p) == 1 .and. acc(p + 12) == 1) then
        R = x1**2 + (2.*u2(p) - 1.)**2
        F = sqrt(-2.*LOG(R)/R)
        x1 = x1*F
        R = x2**2 + (2.*u2(p + 12) - 1.)**2
        F = sqrt(-2.*LOG(R)/R)
        x2 = x2*F
        sg(acc(p + 7)) = sg(acc(p + 7)) + dble(x1)*dble(x2)
        ng12(acc(p + 7)) = ng12(acc(p + 7)) + 1
     end if
  end do

  open (21, file=trim(outfile), status='replace')
  write (21, '(a,i12,a,i10,a,i3,a,i3,a,i5)') 'seed ', Nseed, ' Npairs ', Npairs, ' lux_requested ', lux, &
       ' luxlev_set ', luxlev, ' nskip_used ', nskip
  write (21, '(a,f10.6,a,f10.6)') 'accepted_fraction ', dble(nacc)/Npairs, ' expected ', 3.14159265358979d0/4
  write (21, '(a,es12.4)') 'relation_u(p+12)=u(p+7)-u(p)_fraction_holding ', dble(nres)/(Npairs - 12)
  write (21, '(a,es12.4)') 'x1_p_x1_p+12_mean_product_over_var ', s12/(Npairs - 12)/vx
  do m = 6, 8
     write (21, '(a,i1,a,2es12.4,a,2i10)') 'x1_p_x1_p+12_given_acc_p+', m, '_rejected_accepted ', &
          sc(0, m)/max(nc(0, m), 1_8)/vx, sc(1, m)/max(nc(1, m), 1_8)/vx, '  counts ', nc(0, m), nc(1, m)
  end do
  write (21, '(a,2es12.4,a,2i10)') 'g_p_g_p+12_both_accepted_given_acc_p+7_rejected_accepted ', &
       sg(0)/max(ng12(0), 1_8), sg(1)/max(ng12(1), 1_8), '  counts ', ng12(0), ng12(1)
  do m = 1, MMAX
     write (21, '(a,i3,2es13.4)') 'pairlag ', m, cx(m), ca(m)
  end do
  do L = 1, LMAX
     write (21, '(a,i3,es13.4)') 'gausslag ', L, cg(L)
  end do
  close (21)
  write (*, '(a,i12,a,i3,a,i5,a,es10.2)') ' seed ', Nseed, '  level set ', luxlev, '  nskip used ', nskip, &
       '  relation fraction ', dble(nres)/(Npairs - 12)
end subroutine mode_pairs
