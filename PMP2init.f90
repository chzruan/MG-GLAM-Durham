!--------------------------------------------------
!
!  Generate configuration file Setup.dat
!     - If Init.dat does not exist, write a default file and quit
!     - If Init.dat exists, read it and generate Setup.dat
!
!--------------------------------------------------
Module Param
    Integer*4, parameter :: NtabM = 100000
    Real*4  ::  Om, &       ! cosmology
               Omb, &       ! Omega_massive_neutrinos
               OmL, &
               AEXPN, &
               Sn, &        ! normalization
               sigma8, &
               hubble
    Integer*4 :: NROW, &    ! number of particles in 1D
                 NGRID, &   ! number of grid points in 1D
                 Nout, &    ! number of outputs/mocks
                 Nbiaspars, &   ! number of bias parameters
                 Ncheckpoint, & ! Steps between checkpoints
                 iSave, &       ! Save snapshots
                 iPower, &      ! DM power spectrum
                 iPowerRSD, &   ! Redshift distortions
                 iDensityDistr, & ! DM PDF
                 iBias, &       ! Biasing model
                 iWriteMock, &  ! Write Mocks
                 iDumpPart, &   ! Dump random fraction
                 iBDM           ! Find BDM halos
    Real*4   :: Box, zinit, da, zfinal
    Real*4   :: densThr, sigV
    Real*4   :: zout(1000), BiasPars(1000)
    Real*8   :: rf, OmLOm0
    Real*4   :: xkt(0:NtabM), Pkt(0:NtabM)  ! Power spectrum
    Real*4   :: StepK, alog0
    INteger*4:: Ntab

    ! MG - MG flag
    Integer*4  :: MG_flag = 0       ! 0 - inactive, 1 - active
    Integer*4  :: MG_test = 0       ! 0 - inactive, 1 - active

    Integer*4:: MG_model    ! 1 - f(R) gravity model
    ! 2 - DGP model
    ! 3 - symmetron model
    ! 4 - k-mouflage model
    ! 5 - coupled scale field model

    ! MG - Multigrid methods
    Integer*4  :: Vcycle = 0    ! 0 - inactive, 1 - active
    Integer*4  :: Fcycle = 0    ! 0 - inactive, 1 - active
    Integer*4  :: Wcycle = 0    ! 0 - inactive, 1 - active

    Integer*4:: iter_count_max = 2  ! maximum number of V/F/W cycles in MG solver
    ! usually 2 is enough for V-cycle and 1 for F/W cycles

    Real*8   :: res_conv = 1.0D-8

    ! f(R) model parameters
    Integer*4:: fr_n    ! Hu-Sawicki f(R) power-law index n
    Real*8   :: fR0
    !
    ! DGP model parameters
    Integer*4  :: N_branch = 0  ! 0 - inactive, 1 - active
    Integer*4  :: S_branch = 0  ! 0 - inactive, 1 - active
    Real*8   :: H0rc
    !
    ! symmetron model parameters
    Real*8   :: sym_astar = 0.5D0
    Real*8   :: sym_xi = 1.0D-3
    Real*8   :: sym_beta = 0.1D0
    !
    ! k-mouflage model parameters
    Logical  :: power_law = .FALSE.
    Logical  :: born_infeld = .FALSE.
    Integer*4:: kmf_n = 3
    Real*8   :: kmf_K0 = 1.0D0
    Real*8   :: kmf_lambda = 1.0D0
    Real*8   :: kmf_beta = 0.1D0
    !
    ! coupled scalar field model parameters
    Integer*4 :: csf_coupling = 1
    Integer*4 :: csf_potential = 1
    Real*8    :: csf_alpha = 0.1D0
    Real*8    :: csf_beta = -0.2D0
    Real*8    :: csf_lambda = 1.0D0

Contains

    REAL*8 FUNCTION P(x)
        ! interpolate table with p(k)
        ! x is in real 1/Mpc
        !---------------------------------------
        Real*8 :: x, dk
        If (x .ge. xkt(Ntab)) Then ! slope is ns =-3
            P = Pkt(Ntab)/(x/xkt(Ntab))**3
            Return
        End If
        If (x .lt. xkt(1)) Then ! slope is ns=1
            P = Pkt(1)*(x/xkt(1))
            Return
        End If
        ind = INT((log10(x) - alog0)/StepK) + 1
        dk = xkt(ind + 1) - xkt(ind)
        P = (Pkt(ind)*(xkt(ind + 1) - x) + Pkt(ind + 1)*(x - xkt(ind)))/dk
        Return
    End FUNCTION P

    ! P*k^2*Top-Hat Filter
    REAL*8 FUNCTION Ptophat(wk)
        !---------------------------------------
        real*8 :: wk, X, TOPHAT
        IF (wk .lt. 1.d-4) THEN
            Ptophat = wk*wk**2
        ELSE
            X = wk*rf
            TOPHAT = ((SIN(X) - x*COS(X))*3./X**3)**2
            Ptophat = P(wk)*wk**2*TOPHAT
        END IF
    END FUNCTION Ptophat

    !
    !-----------------------------------  Pcold(k)*k^2
    REAL*8 FUNCTION P2(WK)
        !
        Real*8 :: WK
        P2 = WK**2*P(WK)
    END FUNCTION P2

    !
    !---------------------------------------
    REAL*8 FUNCTION Hnorm(x)
        Real*8 :: x
        Hnorm = Sqrt(1.d+0 + OmLOm0*x**3)
    End FUNCTION Hnorm

    !
    !---------------------------------------
    REAL*8 FUNCTION Hage(x)
        Real*8 :: x
        Hage = sqrt(x)/Hnorm(x)
    End FUNCTION Hage

    !
    !---------------------------------------
    REAL*8 FUNCTION Hgrow(x)
        Real*8 :: x
        Hgrow = (sqrt(x)/Hnorm(x))**3
        !Hgrow =(sqrt(x))**3*(1./sqrt(1.d+0 +OmLOm0*x**3)**3-1. )
    End FUNCTION Hgrow

    !-------------------------------------------------
    !
    !
    !
    SUBROUTINE GrFluctuations4(GrowthDen, a)
        !
        !---------------------------------------
        real*8 ::GrowthDen, a, da, ainit, Omnu, x3, a0, Dp1, D0, Dm1, C, B

        OmLOm0 = OmL/Om
        If (Omb < 0.005) Then      !--- massive neutrino
            Omnu = (1.-Omb/Om)
        Else                   !--- standard LCDM
            Omnu = 1.0
        end If

        !write(*,*) ' Omnu =',Omnu
        da = 1.e-6
        ainit = 3.e-3
        if (ainit > a) Stop ' error in initial a_exp in linear growth'

        !--- initial conditions
        D0 = 1.
        Dm1 = 1.
        a0 = ainit

        Do while (a0 < a)     !--- integrate
            x3 = OmLOm0*a0**3
            C = 0.5*(7.+10.*x3)/(1.+x3)
            B = 1.5/(1.+x3)*(Omnu - 1.-2.*x3)

            Dp1 = (D0*(2.+B*(da/a0)**2) - Dm1*(1.-C*da/a0/2.))/(1.+C*da/a0/2.)
            !write(*,'(es12.4,3x,3es12.4,5x,3es12.4,3x,3es12.4)') a0,Dm1,D0,Dp1,     C,B
            a0 = a0 + da
            Dm1 = D0
            D0 = Dp1
        end Do
        GrowthDen = D0*a0
        !write(*,*) ' Growth =',D0,a0,a
    End SUBROUTINE GrFluctuations4

    !-------------------------------------------------
    !   Age of the Universe: t0 (z=0)
    !   Expansion parameter: a =1/(1+z)
    !   Growth_Rate_Density at a: GrowthDen
    !   GrowthDen =\delta\rho/\rho
    !   normalized to GrowthDen =a for Omega=1
    !   GrowthDen < 1 at a=1
    !
    SUBROUTINE AGE(t0, GrowthDen, a)
        !
        !---------------------------------------
        real*8 :: t0, GrowthDen, a, ww, ww0, ww1, ERR, ANS
        Real*8, PARAMETER :: zero = 1.d-8
        Real*8 INTG

        OmLOm0 = OmL/Om
        t0 = 9.766/hubble/sqrt(Om)*INTG(Hage, zero, a)
        !Hubble = hubble*sqrt(Om/a**3)*Hnorm(a)
        !ww0     = INTG(Hgrow,zero,a/10.d0)
        !ww1     = INTG(Hgrow,a/10.d0,a/2.d0)
        !ww     = INTG(Hgrow,zero,a)  !+ a**2.5/2.5

        CALL DGAUS8(Hgrow, zero, a, ERR, ANS, IERR)
        ww = ANS
        GrowthDen = 2.5*Hnorm(a)/sqrt(a**3)*ww

    End SUBROUTINE AGE

End Module Param

!
!--------------------------------------------------
Program Initialize
    use Param
    CHARACTER :: Name*120, Header*45
    REAL*8 :: INTG, wk, Uklow, Ukup, Sig8, sigma, a, t0, Growthden, Growthden_0
    Real*8, parameter :: PI = 3.1415926535d0
    logical :: exst
    EXTERNAL INTG

    Call CheckInit !---- check whether  all input files exist
    Call ReadInit  !---- read all input information

    OPEN (10, file='Setup.dat')
    OPEN (1, file='lcdm.dat')

    !Omb =0.022
    H = 100.*hubble     ! Hubble constant
    a = 1.
    CALL AGE(t0, GrowthDen_0, a)
    write (*, *) ' Om = ', Om, ' Omb =', Omb
    write (*, '(3(a,es12.4))') ' T_universe= ', t0, ' a=', a, ' GrowthRate= ', GrowthDen_0
    !write (1,'(3(a,es12.4))')' T_universe= ',t0,' a=',a,' GrowthRate= ',GrowthDen_0
    Call GrFluctuations4(GrowthDen_0, a)
    write (*, '(3(a,es12.4))') ' a=', a, ' GrowthRate= ', GrowthDen_0

    rf = 8.        ! top-hat for sigma_8
    Sn = (Sigma8)**2/ &
         INTG(Ptophat, 1.d-5, 5.d+1) ! normalization of P(k)
    !          bias parameter
    Sig8 = sqrt(Sn*(INTG(Ptophat, 1.d-5, 5.d+1)))
    write (*, 10) hubble, Sigma8  ! check if all is self-consistent
    write (1, 10) hubble, Sigma8  ! check if all is self-consistent
    If (abs(Sig8 - sigma8) > 0.005*sigma8) write (*, '(a,2es14.4)') &
        ' Difference betweeen requested and real Sigma8 is:', Sig8, sigma8
10  Format(' Hubble=', f7.3, ' Sigma_8 =', G11.3)

    Uklow = 2.*pi/Box            ! frequency range for integrals
    Ukup = 2.*pi/Box*(NROW/2.)
    AEXPN = 1./(1.+zinit)            ! expansion parameter
    a = AEXPN
    !CALL AGE(t0,GrowthDen,a)
    Call GrFluctuations4(GrowthDen, a)
    sigma = (GrowthDen/GrowthDen_0)* &     !-- use 10th biaspars to tune Pk amplitude
            sqrt(Sn*INTG(P2, Uklow/sqrt(2.), Ukup))*BiasPars(10)
    WRITE (*, 40) zinit, sigma, NROW, Box
    WRITE (1, 40) zinit, sigma, NROW, Box
40  format('  z=', f8.3, ' delta\rho/rho in box=', f9.5, / &
           5X, 'Particles=', i4, ' Box=', f8.2)
    write (1, *) ' k/h        Pk*h^3   Power Spectrum at z=', 1/AEXPN - 1.
    Do i = 1, 1000
        wk = 1.d-3*10.**((i - 1)/20.)
        if (wk > 100.d0) exit
        pp = P(wk)
        write (1, '(2es14.5)') wk, pp
    end Do

    write (10, *) '--------------- Setup file for PMP2 simulations ------------------'
    write (HEADER, '(3(a,i4.4),a,i3.3,a,es8.2,a,f5.3)') &
        'N=', NROW, 'x', NGRID, 'L=', Int(Box), 'zi', INT(zinit), 'da=', da, 'f', da/AEXPN
    write (10, '(a)') TRIM(HEADER)
    write (10, 50) AEXPN, 'Expansion Parameter'
    write (10, '(es12.5,T20,a,es12.4)') da, 'Step in dAEXPN    da/a = ', da/AEXPN
    write (10, 50) sigma, 'DRho/rho in box    '
    write (10, 50) Box, 'Box in  Mpc/h   '
    write (10, 50) sigma8, 'sigma8    '
    write (10, 50) hubble, 'Hubble    '
    write (10, 50) Om, 'Omega Matter'
    write (10, 50) OmL, 'Omega Lambda'
    write (10, 50) Omb, 'Omega Baryons or Neutrinos'
    write (10, 60) NROW, 'NROW  Number Particles   '
    write (10, 60) NGRID, 'NGRID Number grid points '
    write (10, 60) 0, 'Random seed'
    write (10, 50) Box/Ngrid, 'Cell Size   '
    write (10, 50) 2.5841e11*Om*(Box/1000./(NROW/1024.))**3, 'Particle Mass'
    write (10, 50) zinit, 'Initial redshift       '
    write (10, 50) zfinal, 'Final redshift       '
    write (10, 50) DensThr, 'Density Threshold for V correction '
    write (10, 50) sigV, 'rms V correction factor'
    write (10, 60) Nout, 'Number of redshifts for analysis'
    write (10, '(100f8.3)') (zout(i), i=1, Nout)
    write (10, 60) Nbiaspars, 'Number of bias parameters'
    Do i = 1, 9
        write (10, 50) BiasPars(i), 'Bias'
    End Do
    write (10, 50) BiasPars(10), 'Pk tune'
    if (Nbiaspars > 10) Then
        Do i = 11, Nbiaspars
            write (10, 50) BiasPars(i), 'Bias'
        End Do
    end if
    write (10, 60) Ncheckpoint, 'Steps between checkpoints'
    write (10, 60) iSave, 'Save snapshots           '
    write (10, 60) iPower, 'DM power spectrum        '
    write (10, 60) iPowerRSD, 'Redshift distortions     '
    write (10, 60) iDensityDistr, 'DM PDF                   '
    write (10, 60) iBias, 'Biasing model            '
    write (10, 60) iWriteMock, 'Write Mocks              '
    write (10, 60) iDumpPart, 'Dump random fraction     '
    write (10, 60) iBDM, 'Find BDM halos           '
    write (10, *) '!--------------- MG-flag and Multigrid solvers ------------------'
    write (10, 60) MG_flag, 'MG flag'
    write (10, 60) MG_test, 'MG test flag'
    write (10, 60) MG_model, 'MG_model'
    write (10, 60) Vcycle, 'V-cycle'
    write (10, 60) Fcycle, 'F-cycle'
    write (10, 60) Wcycle, 'W-cycle'
    write (10, 60) iter_count_max, 'Maximum number of cycles'
    write (10, 50) res_conv, 'Relaxation convergence criterion'
    write (10, *) '!-------------- f(R) parameters ------------------'
    write (10, 60) fr_n, 'Hu-Sawicki power: n'
    write (10, 50) fR0, 'Present value of scalaron field: fR_0'
    write (10, *) '!--------------- DGP parameters ------------------'
    write (10, 60) N_branch, 'nDGP model'
    write (10, 60) S_branch, 'sDGP model'
    write (10, 50) H0rc, 'DGP parameter: H0rc'
    write (10, *) '!------------ symmetron parameters ---------------'
    write (10, 50) sym_astar, 'symmetron parameter: a_star'
    write (10, 50) sym_xi, 'symmetron parameter: xi'
    write (10, 50) sym_beta, 'symmetron parameter: beta_star'
    write (10, *) '!------------ kmouflage parameters ---------------'
    write (10, 60) power_law, 'power-law type kmouflage model'
    write (10, 60) born_infeld, 'Born-Infeld type kmouflage model'
    write (10, 60) kmf_n, 'power-law kmouflage parameter: n'
    write (10, 50) kmf_K0, 'power-law kmouflage parameter: K0'
    write (10, 50) kmf_beta, 'kmouflage parameter: beta'
    write (10, *) '!------------ csf model parameters ---------------'
    write (10, 60) csf_coupling, 'coupled scalar field model coupling type: 1 - exponential, 2 - quadratic'
    write (10, 60) csf_potential, 'coupled scalar field model potential type: 1 - inverse power law, 2 - SUGRA'
    write (10, 50) csf_alpha, 'coupled scalar field model potential parameter'
    write (10, 50) csf_beta, 'coupled scalar field model coupling parameter'
50  format(es12.5, T20, a)
60  format(i5, T20, a)
70  format(L, T20, a)
    write (*, *) ' Results were written to Setup.dat'
    CLOSE (10)

    a_init = AEXPN
    da_init = da
    Call TestStepping(a_init, da_init)

end Program Initialize
!
!---------------------------------------------------
!
Subroutine ReadInit
    use Param
    !-- Read PkTable. Assign Omegas and hubble
    Omb = ParseLine(10) ! note it is omega_b0 = Omega_b0 h^2 actually...
    Omc = ParseLine(10) ! note it is omega_c0 = Omega_c0 h^2 actually...
    OmL = ParseLine(10)
    Om = ParseLine(10)
    sigma8 = ParseLine(10)
    hubble = SQRT((Omb + Omc)/Om)

    If (Omb < 0.003) Then
        write (*, '(3(a,ES12.3))') ' Model with massive nu: Om_nu =', &
            Omb, ' Om_matter =', Om, ' Sigma8 =', sigma8
    else
        write (*, '(3(a,ES12.3))') ' Model from PkTable file: Ombar =', &
            Omb, ' Om_matter =', Om, ' Sigma8 =', sigma8, ' h =', hubble
    end If

    Ntab = 0
12  READ (10, *, end=32, err=32) xx, pp
    Ntab = Ntab + 1
    xkt(Ntab) = xx  !*hubble
    Pkt(Ntab) = pp  !  Pk
    GOTO 12
32  Write (*, *) ' Read ', Ntab, ' lines from P(k) table '
    close (10)
    If (Ntab .le. 1) stop 'wrong table for p(k)'
    StepK = log10(xkt(Ntab)/xkt(1))/(Ntab - 1)
    alog0 = log10(xkt(1))

    ! test that spacing is the same
    Do k = 2, Ntab - 1
        ss = log10(xkt(k + 1)/xkt(k))
        If (abs(ss/StepK - 1.) .gt. 2.e-2) Then
            Write (*, *) ' error in K spacing. k=', k, xkt(k + 1), xkt(k)
            STOP
        End If
    End Do

    !------------- Read input parameters from Init.dat
    Box = ParseLine(11)
    Nrow = iParseLine(11)
    Ngrid = iParseLine(11)
    sigma8 = ParseLine(11)
    zinit = ParseLine(11)
    da = ParseLine(11)
    zfinal = ParseLine(11)
    Nout = iParseLine(11)
    read (11, *) (zout(i), i=1, Nout)
    densThr = ParseLine(11)
    sigV = ParseLine(11)
    Nbiaspars = iParseLine(11)
    Do i = 1, Nbiaspars
        BiasPars(i) = ParseLine(11)
    End Do
    Ncheckpoint = iParseLine(11)
    iSave = iParseLine(11)
    iPower = iParseLine(11)
    iPowerRSD = iParseLine(11)
    iDensityDistr = iParseLine(11)
    iBias = iParseLine(11)
    iWriteMock = iParseLine(11)
    iDumpPart = iParseLine(11)
    iBDM = iParseLine(11)

    ! MG_flag     = logicalParseLine(11)
    ! MG_test     = logicalParseLine(11)
    MG_flag = iParseLine(11)
    MG_test = iParseLine(11)
    MG_model = iParseLine(11)

    ! Vcycle      = logicalParseLine(11)
    ! Fcycle      = logicalParseLine(11)
    ! Wcycle      = logicalParseLine(11)
    Vcycle = iParseLine(11)
    Fcycle = iParseLine(11)
    Wcycle = iParseLine(11)

    iter_count_max = iParseLine(11)
    res_conv = ParseLine(11)
    fr_n = iParseLine(11)
    fR0 = ParseLine(11)

    N_branch = iParseLine(11)
    S_branch = iParseLine(11)
    ! N_branch    = logicalParseLine(11)
    ! S_branch    = logicalParseLine(11)
    H0rc = ParseLine(11)
    sym_astar = ParseLine(11)
    sym_xi = ParseLine(11)
    sym_beta = ParseLine(11)

    power_law = iParseLine(11)
    born_infeld = iParseLine(11)
    ! power_law   = logicalParseLine(11)
    ! born_infeld = logicalParseLine(11)
    kmf_n = iParseLine(11)
    kmf_K0 = ParseLine(11)
    kmf_beta = ParseLine(11)
    csf_coupling = iParseLine(11)
    csf_potential = iParseLine(11)
    csf_alpha = ParseLine(11)
    csf_beta = ParseLine(11)

    If (BiasPars(10) < 0.1) BiasPars(10) = 1.0
    !--- make new da
    ! fr = da*(1.+zinit)*100.     ! = da/a*100
    ! frnew = INT(fr*10.)/10./100.
    ! da  = frnew/(1.+zinit)
    write (*, *) ' Done reading input from Init.dat'
    close (11)
    write (*, *) 'Box   = ', Box
    write (*, *) 'Nrow  =', Nrow
    write (*, *) 'Ngrid =', Ngrid
    write (*, *) 'Nout  =', Nout
    write (*, *) 'Nbias =', Nbiaspars
end Subroutine ReadInit
!
!---------------------------------------------------
!
!             check if all init files are present
!
Subroutine CheckInit
!---------------------------------------------------
    use Param
    logical :: exst
    Inquire (file='PkTable.dat', exist=exst)
    if (.not. exst) Stop ' File PkTable with the power spectrum not found'
    open (10, file='PkTable.dat')

    Inquire (file='Init.dat', exist=exst)
    if (.not. exst) Then
        write (*, *) ' File Init.dat with initial parameters not found.'
        write (*, *) ' I create a new one. Edit it and restart the code'
        open (11, file='Init.dat')
        write (11, '(a,f9.3)') 'Box      = ', 1000.
        write (11, '(a,i9)') 'Nrow     = ', 1000
        write (11, '(a,i9)') 'Ngrid    = ', 2000
        write (11, '(a,f9.3)') 'sig8     = ', 0.8159
        write (11, '(a,f9.3)') 'z_init   = ', 100.
        write (11, '(a,es11.4)') 'step da  = ', 4e-4
        write (11, '(a,f9.3)') 'z_final  = ', 0.
        write (11, '(a,i9)') '#outputs = ', 10
        write (11, '(20f5.2)') 2.5, 1.5, 1., 0.8, 0.7, 0.5, 0.3, 0.2, 0.1, 0.
        write (11, '(a,f9.3)') 'dens_thr = ', 30.
        write (11, '(a,f9.3)') 'Vrms     = ', sigv
        write (11, '(a,i9)') '#Params  = ', 10
        do i = 1, 9
            write (11, '(a,f9.3)') 'Parametr = ', 0.
        end do
        write (11, '(a,f9.3)') 'Pk tune  = ', 1.
        write (11, '(a,i9,a)') 'Steps between checkpoints = ', 20                                   ! Ncheckpoint
        write (11, '(a,i9,a)') 'Save snapshots            = ', 1, '  0-no, 1-yes, 2-only last'       ! iSave
        write (11, '(a,i9,a)') 'DM power spectrum         = ', 0, '  0-no'                           ! iPower
        write (11, '(a,i9,a)') 'Redshift distortions      = ', 0, '  0-no'                           ! iPowerRSD
        write (11, '(a,i9,a)') 'DM PDF                    = ', 0, '  0-no'                           ! iDensityDistr
        write (11, '(a,i9,a)') 'Biasing model             = ', 0, '  0-no 1-particles 2-density'     ! iBias
        write (11, '(a,i9,a)') 'Write Mocks               = ', 0, '  0-no 1-yes'                     ! iWriteMock
        write (11, '(a,i9,a)') 'Dump random fraction      = ', 0, '  0-no 1-yes'                     ! iDumpPart
        write (11, '(a,i9,a)') 'Find BDM halos            = ', 0, '  0-no 1-yes'                     ! iBDM
        write (11, '(a,i9,a)') 'MG_flag                   = ', 0, '  0 - inactive, 1 - active'
        write (11, '(a,i9,a)') 'MG_test                   = ', 0, '  0 - inactive, 1 - active'
        write (11, '(a,i9,a)') 'MG_model                  = ', 1, '  1 - f(R), 2 - DGP, ...'
        write (11, '(a,i9,a)') 'Vcycle                    = ', 1, '  0 - inactive, 1 - active'
        write (11, '(a,i9,a)') 'Fcycle                    = ', 0, '  0 - inactive, 1 - active'
        write (11, '(a,i9,a)') 'Wcycle                    = ', 0, '  0 - inactive, 1 - active'
        write (11, '(a,i9,a)') 'iter_count_max            = ', 2, '  Maximum Maximum number of cycles'
        write (11, '(a,es11.4,a)') 'res_conv                  = ', 1.0D-8, '  Relaxation convergence certerion'
        write (11, '(a,i9)') 'fR_n                      = ', 1         ! Hu-Sawicki f(R) parameter n
        write (11, '(a,es11.4)') 'fR0                       = ', -1.0D-5   ! Hu-Sawicki f(R) parameter f_R0
        write (11, '(a,i9,a)') 'N_branch                  = ', 0, '  0 - inactive, 1 - active' ! DGP normal branch
        write (11, '(a,i9,a)') 'S_branch                  = ', 0, '  0 - inactive, 1 - active' ! DGP self-acc branch
        write (11, '(a,f9.3)') 'H0rc                      = ', 1.0D0     ! DGP parameter H0*r_c
        write (11, '(a,f9.5)') 'sym_astar                 = ', 0.5D0     ! symmetron parameter a_*
        write (11, '(a,f9.5)') 'sym_xi                    = ', 0.001D0   ! symmetron parameter xi
        write (11, '(a,f9.5)') 'sym_beta                  = ', 0.1D0     ! symmetron parameter beta_*
        write (11, '(a,i9,a)') 'power_law                 = ', 0, '  0 - inactive, 1 - active'   ! power-law type k-mouflage
        write (11, '(a,i9,a)') 'born_infeld               = ', 0, '  0 - inactive, 1 - active'   ! Born-Infeld type k-mouflage
        write (11, '(a,i9)') 'kmf_n                     = ', 3         ! power-law kmouflage parameter n
        write (11, '(a,f9.5)') 'kmf_K0                    = ', 1.0D0     ! power-law kmouflage parameter K0
        write (11, '(a,f9.5)') 'kmf_beta                  = ', 0.1D0     ! power-law kmouflage parameter K0
        write (11, '(a,i9)') 'csf_coupling              = ', 1         ! csf model coupling function type
        write (11, '(a,i9)') 'csf_potential             = ', 1         ! csf model potential type
        write (11, '(a,f9.5)') 'csf_alpha                 = ', 0.1D0     ! csf model potential parameter alpha
        write (11, '(a,f9.5)') 'csf_beta                  = ', -0.2D0    ! csf model coupling parameter beta
        stop
    end if

    open (11, file='Init.dat')

    Inquire (file='TableSeeds.dat', exist=exst)
    if (.not. exst) Call SetSeeds

end Subroutine CheckInit
!------------------------------
!
!             Generate a table with random seeds
!
!------------------------------
Subroutine SetSeeds
    integer*8 :: is, ij, iM
    integer*4 :: Nseed
    is = 1232_8**3

    Nseed = 1298302
    nslip = 137
    Noff = 2357
    Ntable = 5000
    NCount = 0
    open (1, file='TableSeeds.dat')
    write (1, *) 'Seeds:', Nseed, Nslip
    Do ij = 1, is
        x = RANDd(Nseed)
        Nn = INT(x*nslip) + 1
        Do jj = 1, Noff + Nn
            x = RANDd(Nseed)
        End Do
        Ncount = Ncount + 1
        write (1, *) Nseed, Ncount
        If (Ncount > Ntable) exit
    end Do
    close (1)
end Subroutine SetSeeds

!------------------------------------------------
!                                          random number generator
FUNCTION RANDd(M)
!------------------------------------------------
    DATA LC, AM, KI, K1, K2, K3, K4, L1, L2, L3, L4/453815927, &
        2147483648., 2147483647, 536870912, 131072, 256, 16777216, 4, &
        16384, 8388608, 128/
    ML = M/K1*K1
    M1 = (M - ML)*L1
    ML = M/K2*K2
    M2 = (M - ML)*L2
    ML = M/K3*K3
    M3 = (M - ML)*L3
    ML = M/K4*K4
    M4 = (M - ML)*L4
    M5 = KI - M
    IF (M1 .GE. M5) M1 = M1 - KI - 1
    ML = M + M1
    M5 = KI - ML
    IF (M2 .GE. M5) M2 = M2 - KI - 1
    ML = ML + M2
    M5 = KI - ML
    IF (M3 .GE. M5) M3 = M3 - KI - 1
    ML = ML + M3
    M5 = KI - ML
    IF (M4 .GE. M5) M4 = M4 - KI - 1
    ML = ML + M4
    M5 = KI - ML
    IF (LC .GE. M5) ML = ML - KI - 1
    M = ML + LC
    RANDd = M/AM
    RETURN
END FUNCTION RANDd

!--------------------------------------------------
!        read line from  input file iFile
!                real format
Function ParseLine(iFile)
    Character :: Line*120, Line2*120, Line3(120)

    Read (iFile, '(a)') Line
    Ieq = INDEX(Line, '=', BACK=.TRUE.)
    !write(*,*) '  Ieq =',Ieq
    backspace (iFile)                  !--- go to line start
    write (Line2, '(a1,i2,a)') '(', Ieq, 'a1,g12.5)' ! make format
    !write(*,'(a)') Line2
    Read (iFile, Line2) (Line3(i), i=1, Ieq), dummy    ! read
    ParseLine = dummy
    !write(*,'(a,ES12.3)') ' Result =',ParseLine
end Function ParseLine
!--------------------------------------------------
!        read line from  input file iFile
!                          integer format
Function iParseLine(iFile)
    Character :: Line*120, Line2*120, Line3(120)

    Read (iFile, '(a)') Line
    Ieq = INDEX(Line, '=', BACK=.TRUE.)
    backspace (iFile)
    write (Line2, '(a1,i2,a)') '(', Ieq, 'a1,i10)'
    Read (iFile, Line2) (Line3(i), i=1, Ieq), idummy
    iParseLine = idummy

end Function iParseLine
!--------------------------------------------------
!        read line from  input file iFile
!                          logical format
Logical Function logicalParseLine(iFile)
    Character :: Line*120, Line2*120, Line3(120)
    Logical  :: temp_ = .FALSE.

    Line(:) = '0'
    Line2(:) = '0'
    Line3(:) = '0'

    Read (iFile, '(a)') Line
    Ieq = INDEX(Line, '=', BACK=.TRUE.)
    backspace (iFile)
    write (Line2, '(a1,i2,a)') '(', Ieq, 'a1,L,a)'
    ! Read(iFile,Line2)(Line3(i),i=1,Ieq),ldummy
    ! logicalParseLine = ldummy
    Read (iFile, Line2) (Line3(i), i=1, Ieq), temp_, Line
    logicalParseLine = temp_

    Write (*, *) 'Ieq', Ieq
    Write (*, *) 'Line', Line
    Write (*, *) 'Line2', Line2
    Write (*, *) 'Line3', Line3
    Write (*, '(L2)') 'logicalParseLine', logicalParseLine
end Function logicalParseLine
!-------------------------------------------- Simpson integration
REAL*8 FUNCTION INTG(FUNC, A, B)
!---------------------------------------
    IMPLICIT REAL*8(A - H, O - Z)
    PARAMETER(EPS=3.0d-6, JMAX=24)
    EXTERNAL FUNC
    OST = -1.D30
    OS = -1.D30
    ST = 0.
    DO J = 1, JMAX
        CALL TRAPZD(FUNC, A, B, ST, J)
        INTG = (4.0d0*ST - OST)/3.0d0
        IF (ABS(INTG - OS) .Le. EPS*ABS(OS)) RETURN
        OS = INTG
        OST = ST
    end DO
    WRITE (*, *) 'Integration did not converge'
END FUNCTION INTG
!----------------------------------------------
SUBROUTINE TRAPZD(FUNCC, A, B, S, N)
!---------------------------------------
    IMPLICIT REAL*8(A - H, O - Z)
    SAVE IT
    EXTERNAL FUNCC
    IF (N .EQ. 1) THEN
        S = 0.5d0*(B - A)*(FUNCC(A) + FUNCC(B))
        IT = 1
    ELSE
        TNM = IT
        DEL = (B - A)/TNM
        X = A + 0.5D0*DEL
        SUM = 0.0D0
        DO J = 1, IT
            SUM = SUM + FUNCC(X)
            X = X + DEL
        end DO
        S = 0.5D0*(S + (B - A)*SUM/TNM)
        IT = 2*IT
    END IF
END SUBROUTINE TRAPZD
!----------------------------------------------
SUBROUTINE TestStepping(a_init, da_init)
!----------------------------------------------
    StepFactor = da_init/a_init
    a = a_init
    da = da_init
    i = 0
    write (*, '(2(a,es12.5))') ' a_init =', a, ' z_init=', 1./a - 1.
    write (*, '(2(a,es12.5))') ' da     =', da, ' da/a  =', da/a
    write (*, '(2(a,es12.5))') ' stFact =', stepFactor
    write (*, '(a)') ' Step      a       z           da        da/a    step change'
    iCount = 0
    rMax = da/a
    Nsteps2 = 0
    Do
        ind = 0
        If (da < StepFactor/1.25*a .and. a < 0.30) Then
            write (*, *) '                  Increase step:', da, 1.25*StepFactor*a
            da = 1.5*da        ! increase step
            ind = 1
            iCount = iCount + 1
        End If
        a = a + da
        i = i + 1
        rmax = max(rMax, da/a)
        if (a > 0.333) Nsteps2 = Nsteps2 + 1
        write (*, '(i5,6es12.4,i3)') i, a, 1./a - 1., da, da/a, ind
        if (a > 1.) exit
    end Do
    write (*, '(2(a,i5))') ' Number of steps   = ', i, ' n_steps(z<2) =', Nsteps2
    write (*, '(a,i5)') ' Number of changes = ', iCount
    write (*, '(a,f8.4)') ' Maximum da/a      = ', rmax
end SUBROUTINE TestStepping

!--------------------------------------------------------------------
SUBROUTINE DGAUS8(FUN, A, B, ERR, ANS, IERR)
!***BEGIN PROLOGUE  DGAUS8
!***PURPOSE  Integrate a real function of one variable over a finite
!            interval using an adaptive 8-point Legendre-Gauss
!            algorithm.  Intended primarily for high accuracy
!            integration or integration of smooth functions.
!***LIBRARY   SLATEC
!***CATEGORY  H2A1A1
!***TYPE      DOUBLE PRECISION (GAUS8-S, DGAUS8-D)
!***KEYWORDS  ADAPTIVE QUADRATURE, AUTOMATIC INTEGRATOR,
!             GAUSS QUADRATURE, NUMERICAL INTEGRATION
!***AUTHOR  Jones, R. E., (SNLA)
!***DESCRIPTION
!
!     Abstract  *** a DOUBLE PRECISION routine ***
!        DGAUS8 integrates real functions of one variable over finite
!        intervals using an adaptive 8-point Legendre-Gauss algorithm.
!        DGAUS8 is intended primarily for high accuracy integration
!        or integration of smooth functions.
!
!        The maximum number of significant digits obtainable in ANS
!        is the smaller of 18 and the number of digits carried in
!        double precision arithmetic.
!
!     Description of Arguments
!
!        Input--* FUN, A, B, ERR are DOUBLE PRECISION *
!        FUN - name of external function to be integrated.  This name
!              must be in an EXTERNAL statement in the calling program.
!              FUN must be a DOUBLE PRECISION function of one DOUBLE
!              PRECISION argument.  The value of the argument to FUN
!              is the variable of integration which ranges from A to B.
!        A   - lower limit of integration
!        B   - upper limit of integration (may be less than A)
!        ERR - is a requested pseudorelative error tolerance.  Normally
!              pick a value of ABS(ERR) so that DTOL .LT. ABS(ERR) .LE.
!              1.0D-3 where DTOL is the larger of 1.0D-18 and the
!              double precision unit roundoff D1MACH(4).  ANS will
!              normally have no more error than ABS(ERR) times the
!              integral of the absolute value of FUN(X).  Usually,
!              smaller values of ERR yield more accuracy and require
!              more function evaluations.
!
!              A negative value for ERR causes an estimate of the
!              absolute error in ANS to be returned in ERR.  Note that
!              ERR must be a variable (not a constant) in this case.
!              Note also that the user must reset the value of ERR
!              before making any more calls that use the variable ERR.
!
!        Output--* ERR,ANS are double precision *
!        ERR - will be an estimate of the absolute error in ANS if the
!              input value of ERR was negative.  (ERR is unchanged if
!              the input value of ERR was non-negative.)  The estimated
!              error is solely for information to the user and should
!              not be used as a correction to the computed integral.
!        ANS - computed value of integral
!        IERR- a status code
!            --Normal codes
!               1 ANS most likely meets requested error tolerance,
!                 or A=B.
!              -1 A and B are too nearly equal to allow normal
!                 integration.  ANS is set to zero.
!            --Abnormal code
!               2 ANS probably does not meet requested error tolerance.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH, I1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   810223  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890911  Removed unnecessary intrinsics.  (WRB)
!   890911  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!***END PROLOGUE  DGAUS8
    INTEGER IERR, K, KML, KMX, L, LMN, LMX, LR, MXL, NBITS, &
        NIB, NLMN, NLMX
    INTEGER I1MACH
    Real*8 :: A, AA, AE, ANIB, ANS, AREA, B, C, CE, EE, EF, &
              EPS, ERR, EST, GL, GLR, GR, HH, SQ2, TOL, VL, VR, W1, W2, W3, &
              W4, X1, X2, X3, X4, X, H, &
              D1MACH5, D1MACH4
    Real*8 :: D1MACH, G8, FUN
    DIMENSION AA(60), HH(60), LR(60), VL(60), GR(60)
    SAVE X1, X2, X3, X4, W1, W2, W3, W4, SQ2, NLMN, KMX, KML
    DATA X1, X2, X3, X4/ &
        1.83434642495649805D-01, 5.25532409916328986D-01, &
        7.96666477413626740D-01, 9.60289856497536232D-01/
    DATA W1, W2, W3, W4/ &
        3.62683783378361983D-01, 3.13706645877887287D-01, &
        2.22381034453374471D-01, 1.01228536290376259D-01/
    DATA SQ2/1.41421356D0/
    DATA NLMN/1/, KMX/5000/, KML/6/
    G8(X, H) = H*((W1*(FUN(X - X1*H) + FUN(X + X1*H)) &
                   + W2*(FUN(X - X2*H) + FUN(X + X2*H))) &
                  + (W3*(FUN(X - X3*H) + FUN(X + X3*H)) &
                     + W4*(FUN(X - X4*H) + FUN(X + X4*H))))
!***FIRST EXECUTABLE STATEMENT  DGAUS8
!
!     Initialize
!
    D1MACH4 = 1.d-17
!      K = I1MACH(14)
!      ANIB = D1MACH(5)*K/0.30102000D0
!      NBITS = ANIB
!      NLMX = MIN(60,(NBITS*5)/8)
    NBITS = 16
    NLMX = MIN(60, (NBITS*5)/8)
    ANS = 0.0D0
    IERR = 1
    CE = 0.0D0
    IF (A .EQ. B) GO TO 140
    LMX = NLMX
    LMN = NLMN
    IF (B .EQ. 0.0D0) GO TO 10
    IF (SIGN(1.0D0, B)*A .LE. 0.0D0) GO TO 10
    C = ABS(1.0D0 - A/B)
    IF (C .GT. 0.1D0) GO TO 10
    IF (C .LE. 0.0D0) GO TO 140
    ANIB = 0.5D0 - LOG(C)/0.69314718D0
    NIB = ANIB
    LMX = MIN(NLMX, NBITS - NIB - 7)
    IF (LMX .LT. 1) GO TO 130
    LMN = MIN(LMN, LMX)
10  TOL = MAX(ABS(ERR), 2.0D0**(5 - NBITS))/2.0D0

    IF (ERR .EQ. 0.0D0) TOL = SQRT(D1MACH4)
    EPS = TOL
    HH(1) = (B - A)/4.0D0
    AA(1) = A
    LR(1) = 1
    L = 1
    EST = G8(AA(L) + 2.0D0*HH(L), 2.0D0*HH(L))

    K = 8
    AREA = ABS(EST)
    EF = 0.5D0
    MXL = 0
!
!     Compute refined estimates, estimate the error, etc.
!
20  GL = G8(AA(L) + HH(L), HH(L))
    GR(L) = G8(AA(L) + 3.0D0*HH(L), HH(L))
    K = K + 16
    AREA = AREA + (ABS(GL) + ABS(GR(L)) - ABS(EST))
!     IF (L .LT .LMN) GO TO 11
    GLR = GL + GR(L)
    EE = ABS(EST - GLR)*EF
    AE = MAX(EPS*AREA, TOL*ABS(GLR))
    IF (EE - AE) 40, 40, 50
30  MXL = 1
40  CE = CE + (EST - GLR)
    IF (LR(L)) 60, 60, 80
!
!     Consider the left half of this level
!
50  IF (K .GT. KMX) LMX = KML
    IF (L .GE. LMX) GO TO 30
    L = L + 1
    EPS = EPS*0.5D0
    EF = EF/SQ2
    HH(L) = HH(L - 1)*0.5D0
    LR(L) = -1
    AA(L) = AA(L - 1)
    EST = GL
    GO TO 20
!
!     Proceed to right half at this level
!
60  VL(L) = GLR
70  EST = GR(L - 1)
    LR(L) = 1
    AA(L) = AA(L) + 4.0D0*HH(L)
    GO TO 20
!
!     Return one level
!
80  VR = GLR
90  IF (L .LE. 1) GO TO 120
    L = L - 1
    EPS = EPS*2.0D0
    EF = EF*SQ2
    IF (LR(L)) 100, 100, 110
100 VL(L) = VL(L + 1) + VR
    GO TO 70
110 VR = VL(L + 1) + VR
    GO TO 90
!
!     Exit
!
120 ANS = VR
    IF ((MXL .EQ. 0) .OR. (ABS(CE) .LE. 2.0D0*TOL*AREA)) GO TO 140
    IERR = 2
!      write (*,*)  'DGAUS8: ',
!     +   'ANS is probably insufficiently accurate.'
    GO TO 140
130 IERR = -1
    write (*, *) 'DGAUS8: ', 'A and B are too nearly equal to allow normal integration.'
140 IF (ERR .LT. 0.0D0) ERR = CE
    RETURN
END SUBROUTINE DGAUS8
