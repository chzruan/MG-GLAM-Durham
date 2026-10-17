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
End Module ExtraOutputs
