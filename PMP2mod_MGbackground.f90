!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                subroutines that evolve the background
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
!
! Module related with the background evolution of the scalar fields in modified gravity
!
MODULE ExtradofBackgroundData
  use Tools
  implicit none
  !
  integer,parameter :: NumPoints  =100000
  integer,parameter :: NumPointsEx=NumPoints+20
  real*8, parameter :: zst=1000000.0D0 
  real*8, parameter :: a_start_extradof_interp=1.0D-5
  real*8, parameter :: a_start_pert_evol=1.0D-5
  !
  real*8 :: BackgroundEvolution(NumPointsEx,6)                          ! Background quantities in the cosmological evolution
  real*8 :: Param_lambda                                                ! Parameter lambda, to be fixed in the subroutine Fix_Parameter_lambda
  !
  real*8, parameter :: Gnewton=6.6738D-11                               ! Newton constant
  real*8, parameter :: sigma_boltz=5.6704D-8                            ! Stefan-Boltzmann constant
  real*8, parameter :: Mpc=3.085678D22                                  ! Mpc in m
  real*8, parameter :: kappa=8.0D0*DACOS(-1.0D0)*Gnewton                ! 8 \pi G
  real*8, parameter :: col=2.99792458D8                                 ! speed of light
  !
  real*8, parameter :: tcmb=2.7255D0                                    ! CMB temperature today
  real*8, parameter :: N_nu=3.046D0                                     ! number of massless neutrino species

  CONTAINS

  ! ----------------------------------------------------------
  ! ----------------------------------------------------------
  ! Use the consistency condition H=H0 today, to fix parameter 
  ! lambda in the scalar field model by trial-and-error method 
  ! ----------------------------------------------------------
  SUBROUTINE Fix_Parameter_lambda
     !
     use Tools
     implicit none
     !        
     integer           :: ind,i,nu_i
     integer,parameter :: NumEqs=2
     real*8 :: atol
     real*8 :: c(24),w(NumEqs,9),yy(NumEqs),Nst,Nstart,Nend
     real*8 :: dum
     real*8 :: Consistency_flag,StepLength,Length,Length2,Length3,lambda,rho_psi,rho_psi_start,Previous,lambda_previous,crossing,HH_initial
     real*8 :: wrong_counter
     real*8 :: HoH0,HpoH02 
     real*8 :: rhonu
     real*8 :: AA,BB,Sigma0,Sigma1,Gama,expx
     !
     real*8 :: grhog,grhor,grhornomass,grhobc,grhom
     !
     IF(MG_model.NE.4.AND.MG_model.NE.5) RETURN                         ! only solve the background evolution for k-mouflage and coupled scalar field models
     !
     grhog       = kappa/col**2*4.0D0*sigma_boltz/col**3*tcmb**4*Mpc**2 ! present-day density of photons
     grhor       = (7.0D0/8.0D0)*(4.0D0/11.0D0)**(4.0D0/3.0D0)*grhog    ! present-day density of massless neutrinos (per species)
     grhornomass = N_nu*grhor                                           ! present-day density of massless neutrinos (all species)
     grhom       = 3.0D0*(100.0D0*hubble)**2/col**2*1000.0D0**2         ! present-day density of CDM+baryons
     grhobc      = grhom*Om                                             ! present-day density of CDM+baryons
     !
     ! Choose an appropriate initial value for lambda to start the search
     lambda = 1.0D0 
     CALL Value_of_lambda(lambda)
     !
     StepLength       = 0.1D0 
     Consistency_flag = 1.0D0 
     ! 
1234 ind     = 1
     atol    = 2.0D-13                                                  ! tolerance of the numerical solution
     w       = 0.0D0
     Nstart  = -DLOG(1.0D0+zst)                                         ! Starting N=ln(a) 
     Nst     = Nstart
     !
     IF(MG_model.EQ.4) THEN
        ! k-mouflage model
        ! For some reason, the initial condition for K0<0 cannot be
        ! too small; otherwise the code cannot solve the background
        IF(kmf_K0.LT.0.0D0) THEN
           yy(1) = 1.0D-10                                              ! \varphi at initial time
           yy(2) = 1.0D-10                                              ! d\varphi/dN at initial time
        ELSE
           yy(1) = 1.0D-30                                              ! \varphi at initial time
           yy(2) = 1.0D-30                                              ! d\varphi/dN at initial time
        ENDIF 
     ELSE IF(MG_model.EQ.5) THEN
        ! coupled scalar field model
        yy(1) = 1.0D-6                                                  ! Initial value of  \psi   , set to be very small. Final result insensitive to it 
        yy(2) = 1.0D-6                                                  ! Initial value of d\psi/dN, set to be very small. Final result insensitive to it
     ELSE
        WRITE(*,*) 'No need to solve background evolution for this model'
        STOP
     ENDIF
     !
     DO i=1,NumPoints
        Nend = Nstart+i*DABS(Nstart)/NumPoints
        call dverk(dum,NumEqs,EvolveBackground,Nst,yy,Nend,atol,ind,c,NumEqs,w)
        IF(ind.NE.3) THEN
           STOP 'Fix_Parameter_lambda: dverk error!'
        ENDIF  
     ENDDO  
     !     
     HoH0 = (grhog+grhornomass)/grhom                                   ! radiation contribution to 8*pi*G*rho*a^2 at a=1
!    IF(CP%Num_Nu_massive.NE.0) THEN
!       DO nu_i=1,CP%nu_mass_eigenstates
!          call Nu_rho(Cofpsi(yy(1),0)*nu_masses(nu_i),rhonu)
!          HoH0=HoH0+rhonu*grhormass(nu_i)/grhom                        ! massive neutrino contributions
!       ENDDO  
!    ENDIF
     HoH0 = HoH0+Cofpsi(yy(1),0)*grhobc/grhom                           ! matter contributions
     !
     ! k-mouflage
     IF(MG_model.EQ.4) THEN
        !
        HoH0 = HoH0+lambda**2/3.0D0
        !
        AA   = (DBLE(kmf_n)-0.5D0)*kmf_K0*(0.5D0/lambda**2)**(kmf_n-1)*yy(2)**(2*kmf_n)/3.0D0
        BB   = -1.0D0+yy(2)**2/6.0D0
        !
        IF(AA.GT.0.0D0) THEN
           AA = DMAX1(AA, 1.0D-20)
        ELSE
           AA = DMIN1(AA,-1.0D-20)
        ENDIF
        ! 
        IF(kmf_n.EQ.2) THEN
           HoH0   = (-BB-DSQRT(BB*BB-4.0D0*AA*HoH0))/(2.0D0*AA)
        ELSE IF(kmf_n.EQ.3) THEN
           Sigma0 = -3.0D0*AA*BB
           Sigma1 = 27.0D0*AA**2*HoH0
           IF(Sigma1**2.GE.4.0D0*Sigma0**3) THEN
              expx = DSQRT(1.0D0-4.0D0*Sigma0**3/Sigma1**2)
              Gama = (Sigma1/2.0D0*(1.0D0+expx))**(1.0D0/3.0D0)
              HoH0 = -(Gama+Sigma0/Gama)/(3.0D0*AA)
           ELSE
              expx = DACOS(Sigma1/2.0D0/DSQRT(Sigma0**3))
              HoH0 = -2.0D0*DSQRT(Sigma0)*DCOS(expx/3.0D0-2.0*DACOS(0.5D0))/(3.0D0*AA)
           ENDIF 
        ELSE
           STOP 'This n value is currently unsupported for the k-mouflage model. Please try again.'
        ENDIF 
     ENDIF
     ! coupled scalar field model 
     IF(MG_model.EQ.5) THEN
        HoH0 = HoH0+Vofpsi(yy(1),0)
        HoH0 = HoH0/(1.0D0-yy(2)**2/6.0D0)
     ENDIF
     !
     HoH0 = DSQRT(HoH0)
     !
     WRITE(*,*) '-----------------------------------------------------------------------'
     !
     IF(DABS(HoH0-1.0D0)>=1.0D-6) THEN 
        !
        IF((HoH0-1.0D0)*Consistency_flag<0.0D0) THEN
           StepLength = StepLength
        ELSE
           StepLength = -0.5D0*StepLength
        ENDIF  
        !
        Consistency_flag = -(HoH0-1.0D0)/DABS(HoH0-1.0D0)
        !
        lambda = lambda+StepLength
        CALL Value_of_lambda(lambda)
        ! 
        WRITE(*,*) 'lambda is now ',Param_lambda,' and H/H0 is ',HoH0
        !
        GO TO 1234
        !
     ELSE
        !
        CALL Value_of_lambda(lambda)
        IF(MG_model.EQ.4) kmf_lambda = Param_lambda                     ! k-mouflage model
        IF(MG_model.EQ.5) csf_lambda = Param_lambda                     ! coupled scalar field model
        WRITE(*,*) 'Trial-and-error gives Param_lambda =',Param_lambda
        WRITE(*,*) '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'
        !
     ENDIF  

  END SUBROUTINE Fix_Parameter_lambda

  ! ----------------------------------------------------------
  ! ----------------------------------------------------------
  ! Make the value of Param_lambda available in all the module
  ! ----------------------------------------------------------
  SUBROUTINE Value_of_lambda(lambda)
     !
     implicit none
     real*8 :: lambda
     !
     Param_lambda = lambda
     !
  END SUBROUTINE Value_of_lambda

  ! ----------------------------------------------------------
  ! ----------------------------------------------------------
  ! This function returns the deriv-th order derivative of the
  ! coupling function C(phi).
  ! 
  ! The input variable is the dimensionless scalar field given 
  ! by psi=phi/M_Pl, in which M_Pl is the reduced Planck mass.
  !
  ! Coupling_type = 0: No coupling:          
  !                    C(phi) = 1
  ! Coupling_type = 1: Exponential coupling:
  !                    C(phi) = exp(param_beta*phi/M_Pl)
  ! Coupling_type = 2: Quadratic coupling:   
  !                    C(phi) = 1+param_beta*(phi/M_Pl)^2
  ! ----------------------------------------------------------
  FUNCTION Cofpsi(psi,deriv)
     !
     use Tools
     implicit none
     !
     real*8  :: psi,Cofpsi
     integer :: deriv
     integer :: coupling_type
     real*8  :: param_beta
     ! 
     IF(MG_model.EQ.4) THEN
        coupling_type = 1                                       ! K-mouflage model currently only supports exponential coupling
        param_beta    = kmf_beta
     ENDIF
     IF(MG_model.EQ.5) THEN
        coupling_type = csf_coupling 
        param_beta    = csf_beta
     ENDIF
     IF(coupling_type.EQ.0) THEN
        IF     (deriv.EQ.0) THEN
           Cofpsi = 1.0D0 
        ELSE IF(deriv.EQ.1) THEN
           Cofpsi = 0.0D0 
        ELSE IF(deriv.EQ.2) THEN
           Cofpsi = 0.0D0 
        ELSE
           STOP 'Invalid deriv in Cofpsi: do not need derivs higher than second order.'
        ENDIF            
     ELSE IF(coupling_type.EQ.1) THEN
        IF     (deriv.EQ.0) THEN
           Cofpsi =               DEXP(param_beta*psi)
        ELSE IF(deriv.EQ.1) THEN
           Cofpsi = param_beta   *DEXP(param_beta*psi) 
        ELSE IF(deriv.EQ.2) THEN
           Cofpsi = param_beta**2*DEXP(param_beta*psi)
        ELSE
           STOP 'Invalid deriv in Cofpsi: do not need derivs higher than second order.'
        ENDIF 
     ELSE IF(Coupling_type.EQ.2) THEN
        IF     (deriv.EQ.0) THEN
           Cofpsi = 1.0D0+param_beta*psi**2
        ELSE IF(deriv.EQ.1) THEN
           Cofpsi = 2.0D0*param_beta*psi
        ELSE IF(deriv.EQ.2) THEN
           Cofpsi = 2.0D0*param_beta
        ELSE
           STOP 'Invalid deriv in Cofpsi: do not need derivs higher than second order.'
        ENDIF 
     ELSE
        STOP 'Unsupported type of coupling function.'
     ENDIF   

  END FUNCTION Cofpsi

  ! ----------------------------------------------------------
  ! ----------------------------------------------------------
  ! This function computes the derivatives of the scalar field
  ! potential V(\phi) normalised by 3H0^2 to be dimensionless:
  !
  ! \kappa^(1-deriv/2)*d^{deriv}V(\phi)/d^{deriv}\phi/(3*H0^2)
  !
  ! The input variable is the dimensionless scalar field given
  ! by psi=phi/M_Pl, in which M_Pl is the reduced Planck mass.
  ! Remember that kappa = 8*pi*G 
  !
  ! Potential_type = 0: Cosmological constant
  ! Potential_type = 1: Inverse power-law potential V(psi) = M^{4+alpha}/psi^alpha 
  ! Potential_type = 2: SUGRA potential V(psi) = TBD
  ! ----------------------------------------------------------
  FUNCTION Vofpsi(psi,deriv)
     ! 
     use Tools
     !
     implicit none
     real*8  :: psi,Vofpsi
     integer :: deriv
     integer :: potential_type
     real*8  :: param_alpha
     !
     IF(MG_model.EQ.5) THEN
        potential_type = csf_potential
        param_alpha    = csf_alpha
     ENDIF
     ! 
     IF(potential_type.EQ.0) THEN
        IF     (deriv.EQ.0) THEN
           Vofpsi = 1.0D0-Om
        ELSE IF(deriv.EQ.1) THEN
           Vofpsi = 0.0D0
        ELSE IF(deriv.EQ.2) THEN
           Vofpsi = 0.0D0
        ELSE
           STOP 'Invalid deriv in Vofpsi: do not need derivs higher than second.'
        ENDIF            
     ELSE IF(potential_type.EQ.1) THEN
        IF     (deriv.EQ.0) THEN
           Vofpsi = Param_lambda/psi**param_alpha
        ELSE IF(deriv.EQ.1) THEN
           Vofpsi = -Param_lambda*param_alpha                    /psi**(1.0D0+param_alpha) 
        ELSE IF(deriv.EQ.2) THEN
           Vofpsi =  Param_lambda*param_alpha*(1.0D0+param_alpha)/psi**(2.0D0+param_alpha)
        ELSE
           STOP 'Invalid deriv in Vofpsi: do not need derivs higher than second.'
        ENDIF 
     ELSE IF(potential_type.EQ.2) THEN
        IF     (deriv.EQ.0) THEN
           Vofpsi = 1.0D0-Om
        ELSE IF(deriv.EQ.1) THEN
           Vofpsi = 0.0D0
        ELSE IF(deriv.EQ.2) THEN
           Vofpsi = 0.0D0
        ELSE
           STOP 'Invalid deriv in Vofpsi: do not need derivs higher than second.'
        ENDIF 
     ELSE
        STOP 'Invalid input: go to parameter file to specify the right type of potential.'
     ENDIF 
     !
  END FUNCTION Vofpsi

  ! ----------------------------------------------------------
  ! ----------------------------------------------------------
  ! This subroutine gives background derivatives of the scalar 
  ! field wrt N=log(a) for a particular value of Param_lambda.
  ! It is also used in dverk to calculate the derivatives used 
  ! to solve scalar field equation of motion.  Note N=0 today.
  ! ----------------------------------------------------------
  subroutine EvolveBackground(dum,num,N,yy,yyprime)
     !     
     use Tools
     implicit none
     !
     real*8 :: dum
     integer num,nu_i
     real*8 :: yy(num),yyprime(num)
     real*8 :: N,lambda,AA,BB,sigma_sf,temp1,KK,K_s
     real*8 :: rhonu,pnu
     real*8 :: HoH02                                                    ! (H/H0)^2
     real*8 :: HpoH2                                                    ! H'/H^2
     real*8 :: HpoH02                                                   ! H'/H0^2
     real*8 :: drho,rho,p,term1,term2
     real*8 :: Sigma0,Sigma1,Gama,expx
     !
     real*8 :: grhog,grhor,grhornomass,grhobc,grhom
     !
     IF(MG_model.NE.4.AND.MG_model.NE.5) RETURN                         ! only solve the background evolution for k-mouflage and coupled scalar field models
     !
     grhog       = kappa/col**2*4.0D0*sigma_boltz/col**3*tcmb**4*Mpc**2 ! present-day density of photons
     grhor       = (7.0D0/8.0D0)*(4.0D0/11.0D0)**(4.0D0/3.0D0)*grhog    ! present-day density of massless neutrinos (per species)
     grhornomass = N_nu*grhor                                           ! present-day density of massless neutrinos (all species)
     grhom       = 3.0D0*(100.0D0*hubble)**2/col**2*1000.0D0**2         ! present-day density of CDM+baryons
     grhobc      = grhom*Om                                             ! present-day density of CDM+baryons
     !
     lambda = Param_lambda
     !
     ! H^2/H0^2 with H=a'/a and '=d/d\tau, tau being the conformal time
     HoH02  = (grhog+grhornomass)*DEXP(-2.0D0*N)/grhom                  ! radiation contribution to 8*pi*G*rho*a^2 at a
!    if(CP%Num_Nu_massive.NE.0) then
!       do nu_i=1,CP%nu_mass_eigenstates
!          CALL Nu_rho(Cofpsi(yy(1),0)*exp(N)*nu_masses(nu_i),rhonu)
!          HoH02=HoH02+rhonu*grhormass(nu_i)*exp(-2.0*N)/grhom ! massive neutrino contributions
!       end do
!    end if
     HoH02  = HoH02+Cofpsi(yy(1),0)*grhobc*DEXP(-N)/grhom               ! matter contributions
     !
     ! H'/H0^2
     HpoH02 = (grhog+grhornomass)*2.0D0/grhom*DEXP(-2.0D0*N) 
     temp1  = 0.0D0
!    IF(CP%Num_Nu_massive/=0) THEN
!       DO nu_i=1,CP%nu_mass_eigenstates
!          CALL Nu_background(Cofpsi(yy(1),0)*DEXP(N)*nu_masses(nu_i),rhonu,pnu)
!          HpoH02=HpoH02+(rhonu+3.0*pnu)*grhormass(nu_i)*DEXP(-2.0*N)/grhom
!          temp1=temp1+(rhonu-3.0D0*pnu)*grhormass(nu_i)*DEXP(-2.0*N)/grhom
!       ENDDO 
!    ENDIF 
     HpoH02 = HpoH02+Cofpsi(yy(1),0)*grhobc*DEXP(-N)/grhom
     HpoH02 = -HpoH02/2.0D0
     !
     ! k-mouflage model
     IF(MG_model.EQ.4) THEN
        HoH02  = HoH02+lambda**2*DEXP(2.0D0*N)/3.0D0
        AA     = (DBLE(kmf_n)-0.5D0)*kmf_K0*(0.5D0/lambda**2)**(kmf_n-1)*DEXP(-2.0D0*(kmf_n-1.0D0)*N)*yy(2)**(2*kmf_n)/3.0D0
        BB     = -1.0D0+yy(2)**2/6.0D0
        !
        IF(AA.GT.0.0D0) THEN
           AA  = DMAX1(AA, 1.0D-20)
        ELSE
           AA  = DMIN1(AA,-1.0D-20)
        ENDIF
        !  
        IF(kmf_n.EQ.2) THEN
           HoH02   = (-BB-DSQRT(BB*BB-4.0D0*AA*HoH02))/(2.0D0*AA)
        ELSE IF(kmf_n.EQ.3) THEN
           Sigma0  = -3.0D0*AA*BB
           Sigma1  = 27.0D0*AA**2*HoH02
           IF(Sigma1**2.GE.4.0D0*Sigma0**3) THEN
              expx  = DSQRT(1.0D0-4.0D0*Sigma0**3/Sigma1**2)
              Gama  = (Sigma1/2.0D0*(1.0D0+expx))**(1.0D0/3.0D0)
              HoH02 = -(Gama+Sigma0/Gama)/(3.0D0*AA)
           ELSE
              expx  = DACOS(Sigma1/2.0D0/DSQRT(Sigma0**3))
              HoH02 = -2.0D0*DSQRT(Sigma0)*DCOS(expx/3.0D0-2.0D0*DACOS(0.5D0))/(3.0D0*AA)
           ENDIF  
        ELSE
           STOP 'This n value is currently unsupported for the k-mouflage model. Please try again.'
        ENDIF 
        !
        sigma_sf = 0.5D0*HoH02*yy(2)**2*DEXP(-2.0D0*N)/lambda**2
        KK       = -1.0D0+sigma_sf  +kmf_K0*sigma_sf**(kmf_n  )
        K_s      =  1.0D0+kmf_n     *kmf_K0*sigma_sf**(kmf_n-1)
        !
        HpoH02 = HpoH02-K_s*HoH02*yy(2)**2/6.0D0-DEXP(2.0D0*N)*lambda**2*KK/3.0D0
        !
        HpoH2  = HpoH02/HoH02
        !
        ! yy(1) = \varphi
        ! yy(2) = d\varphi/dN
        yyprime(1) = yy(2)
        yyprime(2) = 2.0D0*(1.0D0+DBLE(kmf_n)*(2.0D0-DBLE(kmf_n))*kmf_K0*sigma_sf**(kmf_n-1))*yy(2)       + &
                     (1.0D0+DBLE(kmf_n)*(2.0D0*DBLE(kmf_n)-1.0D0)*kmf_K0*sigma_sf**(kmf_n-1))*HpoH2*yy(2) + &
                     3.0D0*Cofpsi(yy(1),1)*(grhobc/grhom)*DEXP(-N)/HoH02                                  + &
                     3.0D0*Cofpsi(yy(1),1)/Cofpsi(yy(1),0)*temp1/HoH02
        yyprime(2) = -yyprime(2)/(1.0D0+kmf_n*(2.0D0*kmf_n-1.0D0)*kmf_K0*sigma_sf**(kmf_n-1))
     ENDIF
     !
     ! coupled scalar field model
     IF(MG_model.EQ.5) THEN
        HoH02  = HoH02+Vofpsi(yy(1),0)*DEXP(2.0D0*N)
        HoH02  = HoH02/(1.0D0-yy(2)**2/6.0D0)
        HpoH02 = HpoH02+Vofpsi(yy(1),0)*DEXP(2.0D0*N)
        HpoH02 = HpoH02-HoH02*yy(2)**2/3.0D0
        !
        yyprime(1)   = yy(2)
        yyprime(2)   = -(yy(2)*(HpoH02+2.0D0*HoH02)+3.0D0*Vofpsi(yy(1),1)*DEXP(2.0D0*N) + &
                         3.0*(grhobc/grhom)*Cofpsi(yy(1),1)*DEXP(-N)) / HoH02
     ENDIF
     ! 
  END SUBROUTINE EvolveBackground

  ! ----------------------------------------------------------
  ! ----------------------------------------------------------
  ! This subroutine evolves the system of the backgrund scalar
  ! field equations of motion from a given initial time Nstart
  ! to some pre-set late-time grids. The results are stored in 
  ! the array BackgroundEvolution(NumPointsEx,6).  NumpointsEx 
  ! denotes the number of points in the time grid.
  !
  ! The 5 columns of the array are respectively:
  !  
  ! 1 : Values of N on the time grid;
  ! 2 : Values of \varphi on the time grid in which \varphi is 
  !     the dimensionless scalar field phi/M_Pl;
  ! 3 : Values of d\varphi/dN on the time grid;
  ! 4 : Values of d^2\varphi/dN^2 on the time grid;
  ! 5 : Values of sigma=(\varphi')^2/(2*a^2*H0^2*lambda^2) on 
  !     the time grid, where '=d/d\tau, tau the conformal time
  ! 6 : Values of H(a)/H0 where H=\dot{a}/a
  ! ----------------------------------------------------------
  SUBROUTINE extradof_ini_background
     !
     use Tools
     implicit none
     !
     integer:: ind,i,nu_i
     integer,parameter :: NumEqs=2
     real*8 :: Nend,Nfrom,Nstart
     real*8 :: c(24),w(NumEqs,9),yy(NumEqs),atol
     real*8 :: dum
     !
     real*8 :: lambda,HoH02,HpoH02,HpoH2,AA,BB,temp1,sigma_sf,y2prime
     real*8 :: rhonu,pnu
     real*8 :: rho_test,p_test,drho_test,KK,K_s,K_ss,term1,term2,ratio,r_max
     real*8 :: H0,adotoa,CdotoC,SoL
     real*8 :: Sigma0,Sigma1,Gama,expx
     !
     real*8 :: grhog,grhor,grhornomass,grhobc,grhom
     !
     IF(MG_model.NE.4.AND.MG_model.NE.5) RETURN                         ! only solve the background evolution for k-mouflage and coupled scalar field models
     !
     grhog       = kappa/col**2*4.0D0*sigma_boltz/col**3*tcmb**4*Mpc**2 ! present-day density of photons
     grhor       = (7.0D0/8.0D0)*(4.0D0/11.0D0)**(4.0D0/3.0D0)*grhog    ! present-day density of massless neutrinos (per species)
     grhornomass = N_nu*grhor                                           ! present-day density of massless neutrinos (all species)
     grhom       = 3.0D0*(100.0D0*hubble)**2/col**2*1000.0D0**2         ! present-day density of CDM+baryons
     grhobc      = grhom*Om                                             ! present-day density of CDM+baryons
     !
     r_max = 0.0D0
     !
     CALL Fix_Parameter_lambda
     !
     atol = 2.0D-13
     ind  = 1
     w    = 0.0D0
     !
     Nstart = -DLOG(1.0D0+zst)
     Nfrom  = Nstart
     ! 
     IF(MG_model.EQ.4) THEN
        ! k-mouflage model
        IF(kmf_K0.LT.0.0D0) THEN
           yy(1) = 1.0D-10
           yy(2) = 1.0D-10
        ELSE
           yy(1) = 1.0D-30
           yy(2) = 1.0D-30
        ENDIF 
     ELSE IF(MG_model.EQ.5) THEN
        ! coupled scalar field model
        yy(1) = 1.0D-6
        yy(2) = 1.0D-6
     ENDIF
     !
     lambda = Param_lambda
     !
     ! H^2/H0^2 with H=a'/a and '=d/d\tau
     HoH02  = (grhog+grhornomass)*DEXP(-2.0D0*Nstart)/grhom             ! radiation contribution to 8*pi*G*rho*a^2 at a
!    IF(CP%Num_Nu_massive/=0) THEN
!       DO nu_i=1,CP%nu_mass_eigenstates
!          call Nu_rho(Cofpsi(yy(1),0)*DEXP(Nstart)*nu_masses(nu_i),rhonu)
!          HoH02=HoH02+rhonu*grhormass(nu_i)*DEXP(-2.0D0*Nstart)/grhom  ! massive neutrino contributions
!       ENDDO 
!    ENDIF 
     HoH02 = HoH02+Cofpsi(yy(1),0)*grhobc*DEXP(-Nstart)/grhom           ! matter contributions
     !
     ! first time step
     ! H'/H0^2
     HpoH02 = (grhog+grhornomass)*2.0/grhom*DEXP(-2.0D0*Nstart) 
     !
     temp1  = 0.0D0
!    IF(CP%Num_Nu_massive.NE.0) THEN
!       DO nu_i=1,CP%nu_mass_eigenstates
!          call Nu_background(Cofpsi(yy(1),0)*DEXP(Nstart)*nu_masses(nu_i),rhonu,pnu)
!          HpoH02=HpoH02+(rhonu+3.0D0*pnu)*grhormass(nu_i)*DEXP(-2.0D0*Nstart)/grhom
!          temp1=temp1  +(rhonu-3.0D0*pnu)*grhormass(nu_i)*DEXP(-2.0D0*Nstart)/grhom
!       ENDDO 
!    ENDIF 
     !
     HpoH02 = HpoH02+Cofpsi(yy(1),0)*grhobc*DEXP(-Nstart)/grhom
     HpoH02 = -HpoH02/2.0D0
     !
     ! k-mouflage model
     IF(MG_model.EQ.4) THEN
        HoH02 = HoH02+lambda**2*DEXP(2.0D0*Nstart)/3.0D0
        AA    = (DBLE(kmf_n)-0.5D0)*kmf_K0*(0.5D0/lambda**2)**(kmf_n-1)*DEXP(-2.0D0*(kmf_n-1.0D0)*Nstart)*yy(2)**(2*kmf_n)/3.0D0
        BB    = -1.0D0+yy(2)**2/6.0D0
        !
        IF(AA.GT.0.0D0) THEN
           AA = DMAX1(AA, 1.0D-20)
        ELSE
           AA = DMIN1(AA,-1.0D-20)
        ENDIF 
        !
        IF(kmf_n.EQ.2) THEN
           HoH02   = (-BB-DSQRT(BB*BB-4.0D0*AA*HoH02))/(2.0D0*AA)
        ELSE IF(kmf_n.EQ.3) THEN
           Sigma0  = -3.0D0*AA*BB
           Sigma1  = 27.0D0*AA**2*HoH02
           IF(Sigma1**2.GE.4.0D0*Sigma0**3) THEN
              expx  = DSQRT(1.0D0-4.0D0*Sigma0**3/Sigma1**2)
              Gama  = (Sigma1/2.0D0*(1.0D0+expx))**(1.0D0/3.0D0)
              HoH02 = -(Gama+Sigma0/Gama)/(3.0D0*AA)
           ELSE
              expx  = DACOS(Sigma1/2.0D0/DSQRT(Sigma0**3))
              HoH02 = -2.0D0*DSQRT(Sigma0)*DCOS(expx/3.0D0-2.0*DACOS(0.5D0))/(3.0D0*AA)
           ENDIF 
        ELSE
           STOP 'This n value is currently unsupported for the k-mouflage model. Please try again.'
        ENDIF 
        !
        sigma_sf = 0.5D0*HoH02*yy(2)**2*DEXP(-2.0D0*Nstart)/lambda**2
        KK       = -1.0D0+sigma_sf+kmf_K0*sigma_sf**(kmf_n  )
        K_s      =  1.0D0+kmf_n   *kmf_K0*sigma_sf**(kmf_n-1)
        !
        HpoH02 = HpoH02-K_s*HoH02*yy(2)**2/6.0D0-DEXP(2.0D0*Nstart)*lambda**2*KK/3.0D0
        !
        HpoH2  = HpoH02/HoH02
        !
        y2prime = 2.0D0*(1.0D0+DBLE(kmf_n)*(2.0D0-kmf_n)*kmf_K0*sigma_sf**(kmf_n-1))*yy(2)              + &
                  (1.0D0+DBLE(kmf_n)*(2.0D0*DBLE(kmf_n)-1.0D0)*kmf_K0*sigma_sf**(kmf_n-1))*HpoH2*yy(2)  + &
                  3.0D0*Cofpsi(yy(1),1)*grhobc*DEXP(-Nstart)/grhom/HoH02                                + &
                  3.0D0*Cofpsi(yy(1),1)/Cofpsi(yy(1),0)*temp1/HoH02
        y2prime = -y2prime/(1.0D0+kmf_n*(2.0D0*kmf_n-1.0D0)*kmf_K0*sigma_sf**(kmf_n-1))
        !
        BackgroundEvolution(1,1) = DEXP(Nstart)                         ! scale factor a
        BackgroundEvolution(1,2) = yy(1)                                ! \varphi=\phi/M_Pl
        BackgroundEvolution(1,3) = yy(2)                                ! d\varphi/dN
        BackgroundEvolution(1,4) = y2prime                              ! d^2\varphi/dN^2
        BackgroundEvolution(1,5) = DSQRT(HoH02)                         ! aH/H0, H=\dot{a}/a
        BackgroundEvolution(1,6) = sigma_sf                             ! sigma (or X)
        !
     ENDIF
     !
     ! coupled scalar field model
     IF(MG_model.EQ.5) THEN
        !
        HoH02  = HoH02+Vofpsi(yy(1),0)*DEXP(2.0D0*Nstart)
        HoH02  = HoH02/(1.0D0-yy(2)**2/6.0D0)
        HpoH02 = HpoH02+Vofpsi(yy(1),0)*DEXP(2.0D0*Nstart)
        HpoH02 = HpoH02-HoH02*yy(2)**2/3.0D0
        !
        y2prime = -(yy(2)*(HpoH02+2.0D0*HoH02)+3.0D0*Vofpsi(yy(1),1)*DEXP(2.0D0*Nstart) + &
                    3.0*(grhobc/grhom)*Cofpsi(yy(1),1)*DEXP(-Nstart)) / HoH02
        !
        BackgroundEvolution(1,1) = DEXP(Nstart)                         ! scale factor a
        BackgroundEvolution(1,2) = yy(1)                                ! \varphi=\phi/M_Pl
        BackgroundEvolution(1,3) = yy(2)                                ! d\varphi/dN
        BackgroundEvolution(1,4) = y2prime                              ! d^2\varphi/dN^2
        BackgroundEvolution(1,5) = DSQRT(HoH02)                         ! aH/H0, H=\dot{a}/a so that aH=a'/a
        BackgroundEvolution(1,6) = HpoH02                               ! H'/H0, H=a'/a
     ENDIF
     !
     ! later time steps
     DO i=2,NumPointsEx
        Nend = Nstart+DBLE(i-1)*DABS(Nstart)/DBLE(NumPoints)
        CALL dverk(dum,NumEqs,EvolveBackground,Nfrom,yy,Nend,atol,ind,c,NumEqs,w)
        !
        IF(ind.NE.3) THEN
           STOP 'Bad dverk in subroutine extradof_ini_background!'
        ENDIF 
        !
        ! H^2/H0^2 with H=a'/a and '=d/d\tau
        HoH02  = (grhog+grhornomass)*EXP(-2.0*Nend)/grhom               ! radiation contribution to 8*pi*G*rho*a^2 at a
!       IF(CP%Num_Nu_massive.NE.0) THEN
!          DO nu_i=1,CP%nu_mass_eigenstates
!             CALL Nu_rho(Cofpsi(yy(1),0)*DEXP(Nend)*nu_masses(nu_i),rhonu)
!             HoH02=HoH02+rhonu*grhormass(nu_i)*DEXP(-2.0D0*Nend)/grhom ! massive neutrino contributions
!          ENDDO 
!       ENDIF
        HoH02 = HoH02+Cofpsi(yy(1),0)*grhobc*DEXP(-Nend)/grhom          ! matter contributions
        !
        ! H'/H0^2
        HpoH02 = (grhog+grhornomass)*2.0D0/grhom*DEXP(-2.0D0*Nend)
        ! 
        temp1  = 0.0D0
!       IF(CP%Num_Nu_massive.NE.0) THEN
!          DO nu_i=1,CP%nu_mass_eigenstates
!             CALL Nu_background(Cofpsi(yy(1),0)*DEXP(Nend)*nu_masses(nu_i),rhonu,pnu)
!             HpoH02 = HpoH02 + (rhonu+3.0D0*pnu)*grhormass(nu_i)*DEXP(-2.0D0*Nend)/grhom
!             temp1  = temp1  + (rhonu-3.0D0*pnu)*grhormass(nu_i)*DEXP(-2.0D0*Nend)/grhom
!          ENDDO 
!       ENDIF
        ! 
        HpoH02 = HpoH02+Cofpsi(yy(1),0)*grhobc*DEXP(-Nend)/grhom
        HpoH02 = -HpoH02/2.0D0
        !
        ! k-mouflage model
        IF(MG_model.EQ.4) THEN
           !
           HoH02 = HoH02+lambda**2*DEXP(2.0D0*Nend)/3.0D0
           AA    = (kmf_n-0.5D0)*kmf_K0*(0.5D0/lambda**2)**(kmf_n-1)*DEXP(-2.0D0*(kmf_n-1.0D0)*Nend)*yy(2)**(2*kmf_n)/3.0D0
           BB    = -1.0D0+yy(2)**2/6.0D0
           !
           IF(AA.GT.0.0D0) THEN
              AA = DMAX1(AA, 1.0D-20)
           ELSE
              AA = DMIN1(AA,-1.0D-20)
           ENDIF  
           !
           IF(kmf_n.EQ.2) THEN
              HoH02   = (-BB-DSQRT(BB*BB-4.0D0*AA*HoH02))/(2.0D0*AA)
           ELSE IF(kmf_n.EQ.3) THEN
              Sigma0  = -3.0D0*AA*BB
              Sigma1  = 27.0D0*AA**2*HoH02
              IF(Sigma1**2.GE.4.0D0*Sigma0**3) THEN
                 expx  = DSQRT(1.0D0-4.0D0*Sigma0**3/Sigma1**2)
                 Gama  = (Sigma1/2.0D0*(1.0D0+expx))**(1.0D0/3.0D0)
                 HoH02 = -(Gama+Sigma0/Gama)/(3.0D0*AA)
              ELSE
                 expx  = DACOS(Sigma1/2.0D0/DSQRT(Sigma0**3))
                 HoH02 = -2.0D0*DSQRT(Sigma0)*DCOS(expx/3.0D0-2.0D0*DACOS(0.5D0))/(3.0D0*AA)
              ENDIF 
           ELSE
              STOP 'This n value is currently unsupported for the k-mouflage model. Please try again.'
           ENDIF 
           !
           sigma_sf = 0.5D0*HoH02*yy(2)**2*DEXP(-2.0D0*Nend)/lambda**2
           KK       = -1.0D0+sigma_sf+kmf_K0*sigma_sf**(kmf_n  )
           K_s      =  1.0D0+kmf_n   *kmf_K0*sigma_sf**(kmf_n-1)
           !
           HpoH02 = HpoH02-K_s*HoH02*yy(2)**2/6.0D0-DEXP(2.0D0*Nend)*lambda**2*KK/3.0D0
           !
           HpoH2  = HpoH02/HoH02
           !
           y2prime = 2.0D0*(1.0D0+DBLE(kmf_n)*(2.0D0-DBLE(kmf_n))*kmf_K0*sigma_sf**(kmf_n-1))*yy(2) + &
                     (1.0D0+kmf_n*(2.0D0*DBLE(kmf_n)-1.0D0)*kmf_K0*sigma_sf**(kmf_n-1))*HpoH2*yy(2) + &
                     3.0D0*Cofpsi(yy(1),1)*grhobc*DEXP(-Nend)/grhom/HoH02                           + &
                     3.0D0*Cofpsi(yy(1),1)/Cofpsi(yy(1),0)*temp1/HoH02
           y2prime = -y2prime/(1.0D0+kmf_n*(2.0D0*kmf_n-1.0D0)*kmf_K0*sigma_sf**(kmf_n-1))
           !
           BackgroundEvolution(i,1) = DEXP(Nend)
           BackgroundEvolution(i,2) = yy(1)
           BackgroundEvolution(i,3) = yy(2)
           BackgroundEvolution(i,4) = y2prime
           BackgroundEvolution(i,5) = DSQRT(HoH02)               
           BackgroundEvolution(i,6) = sigma_sf
           ! 
        ENDIF
        !
        ! coupled scalar field model
        IF(MG_model.EQ.5) THEN
        !
           HoH02  = HoH02+Vofpsi(yy(1),0)*DEXP(2.0D0*Nend)
           HoH02  = HoH02/(1.0D0-yy(2)**2/6.0D0)
           HpoH02 = HpoH02+Vofpsi(yy(1),0)*DEXP(2.0D0*Nend)
           HpoH02 = HpoH02-HoH02*yy(2)**2/3.0D0
           !
           y2prime = -(yy(2)*(HpoH02+2.0D0*HoH02)+3.0D0*Vofpsi(yy(1),1)*DEXP(2.0D0*Nend) + &
                       3.0*(grhobc/grhom)*Cofpsi(yy(1),1)*DEXP(-Nend)) / HoH02
           !
           BackgroundEvolution(i,1) = DEXP(Nend)                           ! scale factor a
           BackgroundEvolution(i,2) = yy(1)                                ! \varphi=\phi/M_Pl
           BackgroundEvolution(i,3) = yy(2)                                ! d\varphi/dN
           BackgroundEvolution(i,4) = y2prime                              ! d^2\varphi/dN^2
           BackgroundEvolution(i,5) = DSQRT(HoH02)                         ! aH/H0, H=\dot{a}/a
           BackgroundEvolution(i,6) = HpoH02                               ! H'/H0, H=a'/a
           !
        ENDIF
        !
     ENDDO  
     !
  END SUBROUTINE extradof_ini_background

END MODULE ExtradofBackgroundData
