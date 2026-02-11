!-------------------------------------------------
!            Parallel (OMP) PM code
!             2015, 1997  Anatoly Klypin (aklypin@nmsu.edu)
!                                   Astronomy Department, NMSU
!
!
!-------------------------------------------------

Module LocalData
    Integer*4  :: Nsteps, &   ! # steps for this run
                  Ntotal        ! # steps for the whole simulation
    Real*4     :: Alist(2000)   ! List of expansion parameter steps
    Integer*4  :: Nlist(2000)   ! List of steps for analysis
end Module LocalData
!
!-------------------------------------------------
Program PMP
    use Tools
    use fft5
    use Density
    use LocalData
    use Analyze
    use LinkerList
    character*80 :: Path
    Ncheckpoint = 30  !  save for re-starting every  Ncheckpoint step
    Path = ''
    Call Initialize(Path)
    write (*, *) ' Ncheckpoint = ', Ncheckpoint

    IF (MG_test .EQ. 1) Then
        write (*, *) ' MG_test == 1, test mode.'
        AEXPN = 1.0D0
    End If

    !  Main loop
    DO i = 1, Nsteps
        !
        Call TimingMain(0, -1)
        CALL DENSIT                ! Define density

        IF (MG_flag .EQ. 1) THEN
            CALL scalar_field(i)    ! Define scalar field
            IF (MG_test .EQ. 1) EXIT        ! do not do more loops for test
        END IF

        CALL POTENTfft5            ! Define potential

        FactV = 1.                ! normal step
        If (ASTEP < StepFactor/1.25*AEXPN .AND. AEXPN < 0.300) Then
            ASTEP = 1.5*ASTEP       ! increase step
            FactV = 1.2
        END If

        IF (MG_flag .EQ. 1) THEN
            CALL POTENT_FOR_MG
        END IF

        CALL MOVE(FactV)           ! Move particles

        CALL ADDTIME               ! Advance time

        IF (mod(ISTEP, Ncheckpoint) == 0) CALL WriteDataPM(0, Path)  ! checkpoint snapshot

        IF (Nlist(ISTEP) == 1) THEN
            mDENSIT = 1              ! flag for DM density. =0- dm is ready to use
            IF (iSave == 1) CALL WriteDataPM(1, Path)
            CALL Analysis(Path, mDENSIT)
        END IF

        Call TimingMain(0, 1)
        Call TimingMain(0, 0)

        IF (AEXPN .GE. 1.-ASTEP/2.) exit  ! Do again if a < 1
        !
    END DO

    If (iSave .ge. 1 .and. Nlist(ISTEP) /= 1) CALL WriteDataPM(1, Path)     ! make snapshot

END Program PMP

subroutine Initialize(Path)
    use Tools
    use LocalData
    use ExtradofBackgroundData  ! used for MG_model = 4 or 5

    character(len=*) :: Path
    logical :: exst

    OPEN (17, FILE=TRIM(Path)//'Run.log', STATUS='UNKNOWN', position='append')
    WRITE (*, '(A,$)') ' Enter number of steps for this run => '
    READ (*, *) Nsteps         ! Make this number of steps

    Inquire (file='../Setup.dat', exist=exst)
    if (.not. exst) Then
        write (*, *) ' Error: File ../Setup.dat not found. Run PMP2init.exe'
        stop
    end if
    open (11, file='../Setup.dat')

    Inquire (file='../TableSeeds.dat', exist=exst)
    if (.not. exst) Call SetSeeds

    Call ReadSetup
    CALL ReadDataPM(-1, Path)
    myMemory = Memory(1_8*NGRID*NGRID*NGRID)
    Allocate (FI(NGRID, NGRID, NGRID))
    IF (MG_flag .EQ. 1) THEN
        myMemory = Memory(1_8*NGRID*NGRID*NGRID)
        Allocate (FI2(NGRID, NGRID, NGRID))                              ! Allocate scalar field array FI2
        myMemory = Memory(1_8*NGRID*NGRID*NGRID)
        Allocate (FI3(NGRID, NGRID, NGRID))                              ! Allocate scalar field array FI3
    END IF
    write (*, *) ASTEP0, AEXPN0
    StepFactor = ASTEP0/AEXPN0

    IF (AEXPN .GE. 1.) THEN         ! change this if you need to run
        WRITE (*, *) ' Cannot run over a=1' !  beyond a=1
        STOP
    END IF
    write (*, '(a,T20,a,T30,i12,T45,a,i4))') ' Start running: ', &
        'Nparticles:', Nparticles, ' Ngrid= ', Ngrid
    write (*, '(T20,a,T30,i12,T45,a,es12.4)') 'Step:', ISTEP, ' da=', ASTEP
    write (*, '(T20,a,T30,es12.4)') 'StepFactor :', StepFactor

    IF (MG_flag .EQ. 1) THEN
        IF (MG_model .EQ. 4) THEN
            CALL extradof_ini_background
            OPEN (UNIT=27, FILE='background_quantities.txt', FORM='FORMATTED', STATUS='REPLACE')
            DO ii = 1, NumPointsEx
            WRITE (27, '(F20.9,F20.15,F30.15)') (BackgroundEvolution(ii, 1)), BackgroundEvolution(ii, 2), BackgroundEvolution(ii, 5)
            END DO
            CLOSE (27)
        END IF
        IF (MG_model .EQ. 5) THEN
            CALL extradof_ini_background
            OPEN (UNIT=27, FILE='background_quantities.txt', FORM='FORMATTED', STATUS='REPLACE')
            DO ii = 1, NumPointsEx
            WRITE (27, '(F20.9,F20.15,F30.15)') (BackgroundEvolution(ii, 1)), BackgroundEvolution(ii, 2), BackgroundEvolution(ii, 6)
            END DO
            CLOSE (27)
        END IF
    END IF

    !---- make table for steps
    Alist(:) = 0.
    Nlist(:) = 0
    da = ASTEP0
    a = AEXPN0
    Alist(1) = a
    i = 0
    Do
        If (da < StepFactor/1.25*a .and. a < 0.300) Then
            da = 1.5*da        ! increase step
        End If
        a = a + da
        i = i + 1
        If (i .gt. 2000) Stop 'Too many timesteps.Increase length of Alist'
        Alist(i) = a
        If (a .ge. 1.) exit
    end Do
    Ntotal = i
    Do i = 1, Nout           !-- check every zout moment
        a = 1./(1.+zout(i))
        Do j = 2, Ntotal      !-- find closest moment in all steps
            if (a .lt. Alist(j)) exit
        end Do
        da1 = Alist(j) - a
        da0 = a - Alist(j - 1)
        If (da1 < da0) Then    !-- mark closest time step for analysis
            Nlist(j) = 1
        else
            Nlist(j - 1) = 1
        end If
    end Do
    write (*, *) 'Step  a_expansion   redshift   Analyze    List moments for analysis '
    do i = 1, Ntotal
        if (Nlist(i) == 1) &
            write (*, '(i5,2es13.4,i3)') i, Alist(i), 1./Alist(i) - 1., Nlist(i)
    end Do

end subroutine Initialize
!------------------------------
!
!             Generate a table with random seeds
!
!------------------------------
Subroutine SetSeeds
    use Random
    integer*8 :: is, ij, iM
    integer*4 :: Nseed0
    is = 1232_8**3

    Nseed0 = 1298302
    nslip = 137
    noff = 2357
    Ntable = 5000
    NCount = 0
    open (1, file='TableSeeds.dat')
    write (1, *) 'Seeds:', Nseed0, Nslip
    Do ij = 1, is
        x = RANDd(Nseed0)
        Nn = INT(x*nslip) + 1
        Do jj = 1, noff + Nn
            x = RANDd(Nseed0)
        End Do
        Ncount = Ncount + 1
        write (1, *) Nseed0, Ncount
        If (Ncount > Ntable) exit
    end Do
    close (1)
end Subroutine SetSeeds

SUBROUTINE SetTest
    use Tools
    Integer*8 :: ii

    ISTEP = 0
    AEXPN = 0.01
    ASTEP = 0.001
    Om = 1.
    Oml = 0.
    write (*, *) ' Ngrid =', Ngrid
    ii = 0
    DO M3 = 1, NGRID, 2
        DO M2 = 1, NGRID, 2
            DO M1 = 1, NGRID, 2
                ii = ii + 1
                If (ii > Nparticles) Stop ' Number of particles is to big in SetTest'
                XPAR(ii) = M1 + 0.05
                YPAR(ii) = M2 + 0.05
                ZPAR(ii) = M3 + 0.05

                VX(ii) = 0.
                VY(ii) = 0.
                VZ(ii) = 0.
            END DO
        END DO
    END DO
    If (ii /= Nparticles) Stop ' Wrong number of particles in SetTest'
end SUBROUTINE SetTest
!--------------------------------------------------
!           Advance Aexpn, Istep, tIntg,...
!           AEXPN = currnet expansion parameter
!          ISTEP = current step
!          ASTEP = step in the expansion parameter
SUBROUTINE ADDTIME
!--------------------------------------------------
    use Tools

    ISTEP = ISTEP + 1
    AEXPN = AEXPN + ASTEP
    !    Energy conservation
    IF (ISTEP .EQ. 1) THEN
        EKIN1 = EKIN
        EKIN2 = 0.
        EKIN = ENKIN

        !AU0   = AEXP0*ENPOT
        !AEU0  = AEXP0*ENPOT + AEXP0*(EKIN+EKIN1)/2.
        !TINTG = 0.
        WRITE (*, 40) ISTEP, AEXPN, EKIN
        WRITE (17, 40) ISTEP, AEXPN, EKIN
40      FORMAT('**** STEP=', I3, ' A=', F10.5, ' E KIN=', E12.4)
    ELSE
        EKIN2 = EKIN1
        EKIN1 = EKIN
        EKIN = ENKIN

        WRITE (*, 50) ISTEP, AEXPN, EKIN
        WRITE (17, 50) ISTEP, AEXPN, EKIN
    END IF
50  FORMAT('Step = ', I4, ' A=', F8.5, ' Ekin=', E12.4)
END SUBROUTINE ADDTIME
!------------------------------------------------------------

!     Advance each particle:            dA          AEXPN     dA
!           by one step             I______._______I_______.______I          ->  A
!                            i-1     .            i            .          i+1        step
!                                    .          {Fi}            .           .
!                                  { vx }  { x }     .           .
!                                  { vy }  { y }     ^           .
!                                    ._______________^           .
!                                            .                   ^
!                                            .______________^
!
!                          0.5
!                  dP = - A     * Grad(Fi) * dA ; A =AEXPN
!                          i                i
!                                   3/2
!                  dX =         P(new)/A        * dA ; A      =AEXPN+dA/2
!                                 i+1/2                i+1/2
!
!------------------------------------------------
SUBROUTINE MOVE(FactV)
!------------------------------------------------
    use Tools
    use ExtradofBackgroundData ! used for MG_model = 4 (k-mouflage) or 5 (coupled scalar field)

!             PCONST = factor to change velocities
!             XCONST = factor to change coordinates
!                 Note: 0.5 is needed in Pconst because
!                    Fi(i+1)-Fi(i-1) is used as gradient
!             FactV = 1 - normal constant step
!                   = 1.2 - increase step by 1.5
    real*8 :: SVEL, SPHI, PCONST, XCONST, XN, YN, &
              D1, D2, D3, T1, T2, T3, T2T1, T2D1, D2T1, D2D1, &
              GX, GY, GZ, FP, VVx, VVY, VVZ, &
              GX111, GX211, GX121, GX221, &
              GX112, GX212, GX122, GX222, &
              GY111, GY211, GY121, GY221, &
              GY112, GY212, GY122, GY222, &
              GZ111, GZ211, GZ121, GZ221, &
              GZ112, GZ212, GZ122, GZ222, &
              X, Y, Z
    integer*8 :: IN
    integer*4 :: ist    ! MG for background quantity interpolation
    real*8    :: fric_fac, param_beta, bs_fct   ! MG additional force coefficients
    real*8    :: fct1, fct3, dphidN, HoH0, phidot, phidot2, fif2Newt    ! linearised kmoufage fifth force Baojiu-01-07-2021
    Logical   :: Linear_Kmo = .False.

    ! Variables for linear nDGP model
    Logical :: Linear_nDGP = .False.
    Real*8  :: beta, Orc

    Call TimingMain(2, -1)
    PCONST = -SQRT(AEXPN/(Om + OmL*AEXPN**3))*ASTEP*0.5/FactV
    Ahalf = AEXPN + ASTEP/2.
    XCONST = ASTEP/SQRT(Ahalf*(Om + OmL*Ahalf**3))/Ahalf
    !
    IF (MG_model .EQ. 4) THEN
        PCONST = -SQRT(AEXPN/(Om + 9.116748693E-5/AEXPN + (OmL - 9.116748693E-5)*AEXPN**3))*ASTEP*0.5/FactV ! Baojiu-03-07-2021
        XCONST = ASTEP/SQRT(Ahalf*(Om + 9.116748693E-5/Ahalf + (OmL - 9.116748693E-5)*Ahalf**3))/Ahalf       ! Baojiu-03-07-2021
    END IF
    !
    ! background expansion history modified in MG models
    ! MG_model = 4: k-mouflage
    ! MG_model = 5: coupled scalar field
    IF (MG_flag .EQ. 1) THEN
        IF (MG_model .EQ. 4 .OR. MG_model .EQ. 5) THEN
            ist = 1
            DO WHILE (BackgroundEvolution(ist, 1) < AEXPN)
                ist = ist + 1
            END DO
            !
            PCONST = BackgroundEvolution(ist - 1, 5) + &
                   & (BackgroundEvolution(ist, 5) - BackgroundEvolution(ist - 1, 5))/ &
                   & (BackgroundEvolution(ist, 1) - BackgroundEvolution(ist - 1, 1))* &
                   & (AEXPN - BackgroundEvolution(ist - 1, 1))
            XCONST = PCONST
            !
            PCONST = 1.0/PCONST
            PCONST = -PCONST*ASTEP*0.5/FactV
            !
            ! Baojiu-01/07/2021 added block start
            ist = 1
            DO WHILE (BackgroundEvolution(ist, 1) < Ahalf)
                ist = ist + 1
            END DO
            XCONST = BackgroundEvolution(ist - 1, 5) + &
                   & (BackgroundEvolution(ist, 5) - BackgroundEvolution(ist - 1, 5))/ &
                   & (BackgroundEvolution(ist, 1) - BackgroundEvolution(ist - 1, 1))* &
                   & (Ahalf - BackgroundEvolution(ist - 1, 1))
            ! Baojiu-01/07/2021 added block end
            !
            XCONST = Ahalf**2*XCONST
            XCONST = 1.0/XCONST
            XCONST = XCONST*ASTEP
            !
            ! Newtonian force boost factor due to varying particle mass
            bs_fct = BackgroundEvolution(ist - 1, 2) + &
                   & (BackgroundEvolution(ist, 2) - BackgroundEvolution(ist - 1, 2))/ &
                   & (BackgroundEvolution(ist, 1) - BackgroundEvolution(ist - 1, 1))* &
                   & (AEXPN - BackgroundEvolution(ist - 1, 1))
            IF (MG_model .EQ. 4) THEN
                bs_fct = EXP(kmf_beta*bs_fct)
            END IF
            IF (MG_model .EQ. 5) THEN
                IF (csf_coupling .EQ. 1) THEN
                    bs_fct = EXP(csf_beta*bs_fct)
                END IF
            END IF
            !
        END IF
        !
        ! Baojiu-01-07-2021 added block linearised kmoflage fifth-force-to-Newtonian-gravity ratio; start
        IF (MG_model .EQ. 4) THEN
            fct1 = 0.5D0*(1.0D0/(AEXPN*kmf_lambda))**2
            fct3 = kmf_n*kmf_K0*fct1**(kmf_n - 1)                                 ! fct3 is gamma in notes
            ist = 1
            DO WHILE (BackgroundEvolution(ist, 1) < Ahalf)
                ist = ist + 1
            END DO
            ! \varphi itself
            varphi = BackgroundEvolution(ist - 1, 2) + &
                   & (BackgroundEvolution(ist, 2) - BackgroundEvolution(ist - 1, 2))/ &
                   & (BackgroundEvolution(ist, 1) - BackgroundEvolution(ist - 1, 1))* &
                   & (AEXPN - BackgroundEvolution(ist - 1, 1))
            ! d\varphi/dN at AEXPN, N=ln(a)
            dphidN = BackgroundEvolution(ist - 1, 3) + &
                   & (BackgroundEvolution(ist, 3) - BackgroundEvolution(ist - 1, 3))/ &
                   & (BackgroundEvolution(ist, 1) - BackgroundEvolution(ist - 1, 1))* &
                   & (AEXPN - BackgroundEvolution(ist - 1, 1))
            ! H/H0 at AEXPN, here H=a'/a=dN/d\tau with '=d/d\tau
            HoH0 = BackgroundEvolution(ist - 1, 5) + &
                   & (BackgroundEvolution(ist, 5) - BackgroundEvolution(ist - 1, 5))/ &
                   & (BackgroundEvolution(ist, 1) - BackgroundEvolution(ist - 1, 1))* &
                   & (AEXPN - BackgroundEvolution(ist - 1, 1))
            phidot = dphidN*HoH0                  ! (d\varphi/d\tau) / H0
            phidot2 = phidot**2
            fif2Newt = 2.0D0*kmf_beta**2/(1.0D0 + fct3*phidot2**(kmf_n - 1))*DEXP(kmf_beta*varphi)
        END IF
        ! Baojiu-01-07-2021 added block linearised kmoflage fifth-force-to-Newtonian-gravity ratio; end
        !
    END IF
    !
    ! frictional force
    ! MG_model = 3: symmetron
    ! MG_model = 4: k-mouflage
    ! MG_model = 5: coupled scalar field
    IF (MG_flag .EQ. 1) THEN
        !
        fric_fac = 0.0
        !
        ! symmetron case
        IF (MG_model .EQ. 3 .AND. AEXPN .GT. sym_astar) THEN
            fric_fac = -9.0*Om*sym_beta**2*sym_xi**2* &
                     &      SQRT(1.0 - (sym_astar/AEXPN)**3)
        END IF
        !
        ! k-mouflage           (MG_model = 4): only support exponential coupling
        ! function for now
        ! coupled scalar field (MG_model = 5): only support exponential coupling
        ! function for now
        IF (MG_model .EQ. 4 .OR. MG_model .EQ. 5) THEN
            IF (MG_model .EQ. 4) THEN
                param_beta = kmf_beta
            END IF
            IF (MG_model .EQ. 5) THEN
                IF (csf_coupling .EQ. 1) THEN
                    param_beta = csf_beta
                END IF
            END IF
            !
            ! gather d\varphi/dN
            ist = 1
            DO WHILE (BackgroundEvolution(ist, 1) < AEXPN)
                ist = ist + 1
            END DO
            !
            fric_fac = BackgroundEvolution(ist - 1, 3) + &
                     & (BackgroundEvolution(ist, 3) - BackgroundEvolution(ist - 1, 3))/ &
                     & (BackgroundEvolution(ist, 1) - BackgroundEvolution(ist - 1, 1))* &
                     & (AEXPN - BackgroundEvolution(ist - 1, 1))
            fric_fac = fric_fac/AEXPN                                     ! d\varphi/da
            fric_fac = -fric_fac*param_beta
        END IF
        !
    END IF
    !
    SVEL = 0.                      ! counter for \Sum(v_i**2)
    SPHI = 0.                      ! counter for \Sum(phi_i)
    XN = FLOAT(NGRID) + 1.-1.E-8   ! N+1
    YN = FLOAT(NGRID)            ! N
    Wpar = YN**3/FLOAT(Nparticles)

!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (X,Y,Z,GX,GY,GZ,FP,VVx,VVY,VVZ,I,J,K) &
!$OMP PRIVATE (I1,J1,K1,K2,K3,I2,J2,I0,J0) &
!$OMP PRIVATE (D1,D2,D3,T1,T2,T3, T2T1,T2D1,D2T1,D2D1) &
!$OMP PRIVATE (F111,F211,F121,F221,F112,F212,F122,F222) &
!$OMP PRIVATE (F113,F213,F123,F223,F110,F210,F120,F220) &
!$OMP PRIVATE (F311,F321,F131,F231,F312,F322,F132,F232) &
!$OMP PRIVATE (F011,F021,F101,F201,F012,F022,F102,F202) &
!$OMP PRIVATE (GX111,GX211,GX121,GX221,GX112,GX212,GX122,GX222) &
!$OMP PRIVATE (GY111,GY211,GY121,GY221,GY112,GY212,GY122,GY222) &
!$OMP PRIVATE (GZ111,GZ211,GZ121,GZ221,GZ112,GZ212,GZ122,GZ222) &
!$OMP REDUCTION(+:SVEL,SPHI)
    DO IN = 1, Nparticles   ! Loop over particles
        X = XPAR(IN)
        Y = YPAR(IN)
        Z = ZPAR(IN)
        VVX = VX(IN)
        VVY = VY(IN)
        VVZ = VZ(IN)
        I = INT(X)
        J = INT(Y)
        K = INT(Z)
        D1 = X - FLOAT(I)
        D2 = Y - FLOAT(J)
        D3 = Z - FLOAT(K)
        T1 = 1.-D1
        T2 = 1.-D2
        T3 = 1.-D3
        T2T1 = T2*T1
        T2D1 = T2*D1
        D2T1 = D2*T1
        D2D1 = D2*D1
        I1 = I + 1
        IF (I1 .GT. NGRID) I1 = 1
        J1 = J + 1
        IF (J1 .GT. NGRID) J1 = 1
        K1 = K + 1
        IF (K1 .GT. NGRID) K1 = 1
        K2 = K + 2
        IF (K2 .GT. NGRID) K2 = K2 - NGRID
        K3 = K - 1
        IF (K3 .LT. 1) K3 = NGRID
        F111 = FI(I, J, K)  !  Read potential to Fij vars
        F211 = FI(I1, J, K)
        F121 = FI(I, J1, K)
        F221 = FI(I1, J1, K)
        !
        F112 = FI(I, J, K1)
        F212 = FI(I1, J, K1)
        F122 = FI(I, J1, K1)
        F222 = FI(I1, J1, K1)
        !
        F113 = FI(I, J, K2)
        F213 = FI(I1, J, K2)
        F123 = FI(I, J1, K2)
        F223 = FI(I1, J1, K2)
        !
        F110 = FI(I, J, K3)
        F210 = FI(I1, J, K3)
        F120 = FI(I, J1, K3)
        F220 = FI(I1, J1, K3)
        !
        I2 = I + 2
        IF (I2 .GT. NGRID) I2 = I2 - NGRID
        J2 = J + 2
        IF (J2 .GT. NGRID) J2 = J2 - NGRID
        F311 = FI(I2, J, K)
        F321 = FI(I2, J1, K)
        F131 = FI(I, J2, K)
        F231 = FI(I1, J2, K)
        !
        F312 = FI(I2, J, K1)
        F322 = FI(I2, J1, K1)
        F132 = FI(I, J2, K1)
        F232 = FI(I1, J2, K1)
        !
        I0 = I - 1
        IF (I0 .LT. 1) I0 = NGRID
        J0 = J - 1
        IF (J0 .LT. 1) J0 = NGRID
        F011 = FI(I0, J, K)
        F021 = FI(I0, J1, K)
        F101 = FI(I, J0, K)
        F201 = FI(I1, J0, K)
        !
        F012 = FI(I0, J, K1)
        F022 = FI(I0, J1, K1)
        F102 = FI(I, J0, K1)
        F202 = FI(I1, J0, K1)
        ! Find {2*gradient} in nods
        GX111 = F211 - F011
        GX211 = F311 - F111
        GX121 = F221 - F021
        GX221 = F321 - F121
        !
        GX112 = F212 - F012
        GX212 = F312 - F112
        GX122 = F222 - F022
        GX222 = F322 - F122
        !
        GY111 = F121 - F101
        GY211 = F221 - F201
        GY121 = F131 - F111
        GY221 = F231 - F211
        !
        GY112 = F122 - F102
        GY212 = F222 - F202
        GY122 = F132 - F112
        GY222 = F232 - F212
        !
        GZ111 = F112 - F110
        GZ211 = F212 - F210
        GZ121 = F122 - F120
        GZ221 = F222 - F220
        !
        GZ112 = F113 - F111
        GZ212 = F213 - F211
        GZ122 = F123 - F121
        GZ222 = F223 - F221
        ! Interpolate to the point
        GX = PCONST*(T3*(T2T1*GX111 + T2D1*GX211 + D2T1*GX121 + D2D1*GX221) + &
                     D3*(T2T1*GX112 + T2D1*GX212 + D2T1*GX122 + D2D1*GX222))

        GY = PCONST*(T3*(T2T1*GY111 + T2D1*GY211 + D2T1*GY121 + D2D1*GY221) + &
                     D3*(T2T1*GY112 + T2D1*GY212 + D2T1*GY122 + D2D1*GY222))

        GZ = PCONST*(T3*(T2T1*GZ111 + T2D1*GZ211 + D2T1*GZ121 + D2D1*GZ221) + &
                     D3*(T2T1*GZ112 + T2D1*GZ212 + D2T1*GZ122 + D2D1*GZ222))

        !         Find potential of the point
        !      FP=         T3*(T2T1*F111+T2D1*F211 +D2T1*F121+D2D1*F221 )+  &
        !                       D3*(T2T1*F112+T2D1*F212 +D2T1*F122+D2D1*F222 )
        !         SPHI = SPHI + FP*WPAR
        !
        IF (MG_flag .EQ. 1) THEN
            !
            ! frictional force
            IF ((MG_model .EQ. 3) .AND. (AEXPN .GT. sym_astar)) THEN
                VVX = VVX + fric_fac*FI2(I, J, K)*VVX*ASTEP/FactV                ! x component
                VVY = VVY + fric_fac*FI2(I, J, K)*VVY*ASTEP/FactV                ! y component
                VVZ = VVZ + fric_fac*FI2(I, J, K)*VVZ*ASTEP/FactV                ! z component
            END IF
            !
            IF ((MG_model .EQ. 4) .OR. (MG_model .EQ. 5)) THEN
                VVX = VVX + fric_fac*VVX*ASTEP/FactV                ! x component
                VVY = VVY + fric_fac*VVY*ASTEP/FactV                ! y component
                VVZ = VVZ + fric_fac*VVZ*ASTEP/FactV                ! z component
            END IF
            !
            ! boost of Newtonian force due to varying particle mass
            IF ((MG_model .EQ. 4) .OR. (MG_model .EQ. 5)) THEN
                VVX = VVX + (bs_fct - 1.0)*GX                                    ! x component
                VVY = VVY + (bs_fct - 1.0)*GY                                    ! y component
                VVZ = VVZ + (bs_fct - 1.0)*GZ                                    ! z component
            END IF
        END IF
        !
        IF (Linear_nDGP) THEN
            Orc = 1.0D0/(4.0D0*H0rc**2)
            beta = 1.0D0 + (0.5D0*Om/AEXPN**3 + OmL)/(DSQRT(Orc*(Om/AEXPN**3 + OmL))) ! normal branch
            VVX = VVX + GX*(1.0D0 + 1.0D0/(3.0D0*beta))                          ! Move points
            VVY = VVY + GY*(1.0D0 + 1.0D0/(3.0D0*beta))
            VVZ = VVZ + GZ*(1.0D0 + 1.0D0/(3.0D0*beta))
        ELSE
            VVX = VVX + GX                                                       ! Move points
            VVY = VVY + GY
            VVZ = VVZ + GZ
        END IF
        !
        ! Baojiu-01-07-2021 added block linearised kmoflage fifth force; start
        ! Comment out these lines unless you want to do a linearised kmouflage simulation!
        IF (Linear_Kmo) THEN
            VVX = VVX + GX*fif2Newt
            VVY = VVY + GY*fif2Newt
            VVZ = VVZ + GZ*fif2Newt
        END IF
        ! Baojiu-01-07-2021 added block linearised kmoflage fifth force; end
        !
        X = X + VVX*XCONST
        Y = Y + VVY*XCONST
        Z = Z + VVZ*XCONST
        IF (X .LT. 1.d0) X = X + YN                                                ! Periodical conditions
        IF (X .GE. XN) X = X - YN
        IF (Y .LT. 1.d0) Y = Y + YN
        IF (Y .GE. XN) Y = Y - YN
        IF (Z .LT. 1.d0) Z = Z + YN
        IF (Z .GE. XN) Z = Z - YN
        !
        !
        SVEL = SVEL + (VVX**2 + VVY**2 + VVZ**2)*WPAR
        !
        XPAR(IN) = X                                                         ! Write new coordinates
        YPAR(IN) = Y
        ZPAR(IN) = Z
        VX(IN) = VVX
        VY(IN) = VVY
        VZ(IN) = VVZ
        if (INT(Xpar(IN)) == Ngrid + 1) Xpar(IN) = Xpar(IN) - 1.e-3
        if (INT(Ypar(IN)) == Ngrid + 1) Ypar(IN) = Ypar(IN) - 1.e-3
        if (INT(Zpar(IN)) == Ngrid + 1) Zpar(IN) = Zpar(IN) - 1.e-3
    END DO
    ! Set energies:
    ! Kin energy now at A(i+1/2)
    ! Pot energy     at A(i)
    ENKIN = SVEL/2./(AEXPN + ASTEP/2.)**2
    ENPOT = SPHI/2.

    Call TimingMain(2, 1)

end SUBROUTINE MOVE
!-------------------------------------------------
!            Find potential on Grid FI:        DENSITY    ->        POTENTIAL
!
!                   O 1                    ^ - Fourier component
!                   |
!             1           |-4         1        ^      ^        2Pi
!             O-----O-----O     Fi    =        Rho        / (2cos(---  (i-1))+
!                   |           i,j                i,j        Ngrid
!                   |
!                   O 1                          2Pi
!                                       2cos(---  (j-1))-4)
!                       ^                        Ngrid
!                       Fi        = 1 (?)
!                         11
!                   2
!                NABLA  Fi = 3/2  /A * (Rho - <Rho>) ;
!                   X
!                              <Rho> = 1
!
SUBROUTINE POTENTfft5
!---------------------------------------------
    use Tools
    use fft5
    integer*4, parameter :: Nlensav = 16384
    integer*4, parameter :: Nlenwrk = 16384
    real*8, parameter :: P16 = 6.28318530718
    real*4, parameter :: sq2 = 1. !1.41421356237
    real*8, save        :: wsave(1:Nlensav)
    real*8, save        :: work(1:Nlenwrk)
    REAL*8               :: XX, D1, D2, A1, A2, A3, wi, wj, wk
    Integer*4            :: OMP_GET_MAX_THREADS, OMP_GET_THREAD_NUM
    Integer*4            :: Ng, ier, lensav, lenwrk, lenr, inc
    Real*8  :: GREENf(Nlenwrk)
    real*8  :: r(Nlenwrk)
!$OMP THREADPRIVATE(work,wsave)

    Call TimingMain(1, -1)
    If (Ngrid > Nlenwrk) Stop ' Incresase Nlenwrk in POTENTfft5'
    Ng = Ngrid
    lensav = Ngrid + int(log(real(Ngrid, kind=4))/log(2.0E+00)) + 4
    lenwrk = Ngrid

    If (Omb < 0.005) Then    !--- massive neutrino
        trfi = 1.5*(Om - Omb)/aexpn
    Else                 !--- standard LCDM
        trfi = 1.5*Om/aexpn
    end If

    ! Set Green function components

    GREENf(1) = 2.
    GREENf(Ngrid) = -2.
    DO i = 1, Ngrid/2 - 1
        XX = 2.*COS(P16*i/Ngrid)
        GREENf(2*i) = XX
        GREENf(2*i + 1) = XX
    End DO

    call rfft1i(Ng, wsave, lensav, ier) !   Initialize FFT
    inc = 1
    lenr = Ngrid

    write (*, *) ' time =', seconds(), ' Start POTENT'
!$OMP PARALLEL DO DEFAULT(SHARED)  copyin(wsave,work) &
!$OMP PRIVATE ( k,j,i ,r,ier)
    Do k = 1, NGRID             ! fft for xy planes in x-dir
        Do j = 1, NGRID
            Do i = 1, NGRID
                r(i) = FI(i, j, k)
            End Do
            call rfft1f(Ng, inc, r, lenr, wsave, lensav, work, lenwrk, ier)
            Do i = 1, NGRID
                FI(i, j, k) = r(i)
            End Do
        End Do

        Do i = 1, NGRID        ! fft xy planes in y-dir
            Do j = 1, NGRID
                r(j) = FI(i, j, k)
            End Do
            call rfft1f(Ng, inc, r, lenr, wsave, lensav, work, lenwrk, ier)
            Do j = 1, NGRID
                FI(i, j, k) = r(j)
            End Do
        End Do

    End Do

    write (*, *) ' time =', seconds(), ' Start transposition'

!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE ( k,j,i,aa)
    DO J = 1, Ngrid
    DO K = 1, Ngrid - 1
        DO I = K + 1, Ngrid
            aa = FI(I, J, K)
            FI(I, J, K) = FI(K, J, I)
            FI(K, J, I) = aa
        END DO
    END DO
    END DO
    write (*, *) ' time =', seconds(), ' Start z-direction'

!$OMP PARALLEL DO DEFAULT(SHARED)  copyin(wsave,work) &
!$OMP PRIVATE ( k,j,i ,r, ier,A1,A2,A3)
    Do j = 1, NGRID     ! ------ z-direction
        Do i = 1, NGRID
            Do k = 1, NGRID
                r(k) = FI(k, j, i)
            End Do
            call rfft1f(Ng, inc, r, lenr, wsave, lensav, work, lenwrk, ier)
            Do k = 1, NGRID
                FI(k, j, i) = r(k)
            End Do
        End Do
    end Do
    ww = (P16/Ngrid)**2
!$OMP PARALLEL DO DEFAULT(SHARED)  copyin(wsave,work) &
!$OMP PRIVATE ( k,j,i ,r, ier,A1,A2,A3,wi,wj,wk)
    Do j = 1, NGRID     ! ------ z-direction
        A3 = GREENf(J) - 6.
        !wj = j/2
        Do i = 1, NGRID
            A2 = GREENf(I) + A3
            !wi = i/2
            Do k = 1, NGRID
                A1 = A2 + GREENf(K)               !--- use this for descrete Poisson solver
                ! wk = k/2
                !A1 = -ww*(wi**2+wj**2+wk**2)   !--- use this for k**2 Green functions
                IF (ABS(A1) .LT. 1.d-7) A1 = 1.
                r(k) = FI(k, j, i)*trfi/A1
            End Do
            call rfft1b(Ng, inc, r, lenr, wsave, lensav, work, lenwrk, ier)
            Do k = 1, NGRID
                FI(k, j, i) = r(k)
            End Do
        End Do
    end Do
    write (*, *) ' time =', seconds(), ' end z-direction'

!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE ( k,j,i,aa)
    DO J = 1, Ngrid
    DO K = 1, Ngrid - 1
        DO I = K + 1, Ngrid
            aa = FI(I, J, K)
            FI(I, J, K) = FI(K, J, I)
            FI(K, J, I) = aa
        END DO
    END DO
    END DO
    write (*, *) ' time =', seconds(), ' start xy'

!$OMP PARALLEL DO DEFAULT(SHARED) copyin(wsave,work) &
!$OMP PRIVATE ( k,j,i ,r,ier)
    Do k = 1, NGRID             ! fft for xy planes in x-dir
        Do j = 1, NGRID
            Do i = 1, NGRID
                r(i) = FI(i, j, k)
            End Do
            call rfft1b(Ng, inc, r, lenr, wsave, lensav, work, lenwrk, ier)
            Do i = 1, NGRID
                FI(i, j, k) = r(i)
            End Do
        End Do

        Do i = 1, NGRID        ! fft xy planes in y-dir
            Do j = 1, NGRID
                r(j) = FI(i, j, k)
            End Do
            call rfft1b(Ng, inc, r, lenr, wsave, lensav, work, lenwrk, ier)
            Do j = 1, NGRID
                FI(i, j, k) = r(j)
            End Do
        End Do

    End Do
    write (*, *) ' time =', seconds(), ' Finished Potent'

    Call TimingMain(1, 1)

end SUBROUTINE POTENTfft5

SUBROUTINE POTENT_FOR_MG
    use Tools
    real*8  :: ctilde, R_bg, R0_bg, fR_bg, fct1, sf
    integer :: M1, M2, M3

!  OPEN(UNIT=29, FILE='phi_newtonian.txt', FORM='FORMATTED', STATUS='REPLACE')
!  OPEN(UNIT=30, FILE='phi_MG_N1.txt', FORM='FORMATTED', STATUS='REPLACE')

!  DO M1=1,NGRID
!     WRITE(29,'(F23.9,F23.15)') DBLE(M1)/DBLE(NGRID),FI(M1,NGRID/2,NGRID/2)
!  END DO
!  CLOSE(29)

    SELECT CASE (MG_model)
        !
    CASE (1) ! f(R) gravity part
        ctilde = 2.99792458D3*DBLE(NGRID)/Box
        R_bg = 3.0d0*(Om/AEXPN**3 + 4.0d0*OmL)
        R0_bg = 3.0d0*(Om + 4.0d0*OmL)
        fR_bg = fR0*(R0_bg/R_bg)**(fr_n + 1)
        fct1 = ctilde**2*0.5D0*fR_bg
        !
        ! Case of f(R) n=0
        IF (fr_n .EQ. 0) THEN
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3) &
!$OMP PRIVATE (sf)
            DO M3 = 1, NGRID
                DO M2 = 1, NGRID
                    DO M1 = 1, NGRID
                        sf = -fct1*FI2(M1, M2, M3)
!$OMP ATOMIC
                        FI(M1, M2, M3) = FI(M1, M2, M3) + sf
                    END DO
                END DO
            END DO
            !
            ! Case of f(R) n=1
        ELSE IF (fr_n .EQ. 1) THEN
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3) &
!$OMP PRIVATE (sf)
            DO M3 = 1, NGRID
                DO M2 = 1, NGRID
                    DO M1 = 1, NGRID
                        sf = -fct1*FI2(M1, M2, M3)**2
!$OMP ATOMIC
                        FI(M1, M2, M3) = FI(M1, M2, M3) + sf
                    END DO
                END DO
            END DO
            !
            ! Case of f(R) n=2
        ELSE IF (fr_n .EQ. 2) THEN
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3) &
!$OMP PRIVATE (sf)
            DO M3 = 1, NGRID
                DO M2 = 1, NGRID
                    DO M1 = 1, NGRID
                        sf = -fct1*FI2(M1, M2, M3)**3
!$OMP ATOMIC
                        FI(M1, M2, M3) = FI(M1, M2, M3) + sf
                    END DO
                END DO
            END DO
        END IF
        !
        ! DGP part
    CASE (2)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3) &
!$OMP PRIVATE (sf)
        DO M3 = 1, NGRID
            DO M2 = 1, NGRID
                DO M1 = 1, NGRID
                    sf = 0.5D0*FI2(M1, M2, M3)
!$OMP ATOMIC
                    FI(M1, M2, M3) = FI(M1, M2, M3) + sf
                END DO
            END DO
        END DO
        !
        ! symmetron part
    CASE (3)
        ctilde = 2.99792458D3*DBLE(NGRID)/Box
        fct1 = 3.0D0*sym_xi**2*Om*sym_beta**2*ctilde**2/sym_astar**3
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3) &
!$OMP PRIVATE (sf)
        DO M3 = 1, NGRID
            DO M2 = 1, NGRID
                DO M1 = 1, NGRID
                    sf = fct1*FI2(M1, M2, M3)**2
!$OMP ATOMIC
                    FI(M1, M2, M3) = FI(M1, M2, M3) + sf
                END DO
            END DO
        END DO
        !
        ! k-mouflage part
    CASE (4)
        ctilde = 2.99792458D3*DBLE(NGRID)/Box
        fct1 = kmf_beta*ctilde
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3) &
!$OMP PRIVATE (sf)
        DO M3 = 1, NGRID
            DO M2 = 1, NGRID
                DO M1 = 1, NGRID
                    sf = fct1*FI2(M1, M2, M3)
!$OMP ATOMIC
                    ! Comment out the line below if you want to do a linearised kmouflage simulation!
                    FI(M1, M2, M3) = FI(M1, M2, M3) + sf
                END DO
            END DO
        END DO
        !
        ! coupled scalar field part
    CASE (5)
        IF (csf_coupling .EQ. 1) fct1 = csf_beta
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3) &
!$OMP PRIVATE (sf)
        DO M3 = 1, NGRID
            DO M2 = 1, NGRID
                DO M1 = 1, NGRID
                    sf = fct1*FI2(M1, M2, M3)
!$OMP ATOMIC
                    FI(M1, M2, M3) = FI(M1, M2, M3) + sf
                END DO
            END DO
        END DO
        !
        !
    CASE DEFAULT
        WRITE (*, '(A)') 'Model currently unsupported. Please specify again.'
        !
        !
    END SELECT

!  DO M1=1,NGRID
!     WRITE(30,'(F23.9,F23.15)') DBLE(M1)/DBLE(NGRID),FI(M1,NGRID/2,NGRID/2)
!  END DO
!  CLOSE(30)

END SUBROUTINE POTENT_FOR_MG

