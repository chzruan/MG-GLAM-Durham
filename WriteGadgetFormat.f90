module WriteGadgetFormat

    real(kind=4), dimension(:,:), allocatable :: pos ! positions of particles in Mpc
    real(kind=4), dimension(:,:), allocatable :: vel ! velocities of particles in km/s
    integer*8,     dimension(:),  allocatable :: pid ! particle IDs

    contains


    subroutine write_gadget(fileout, npart, nparttot, &
        verbose, megaverbose, &
        Lbox, numfiles, hubble, omega0, omegaL, aexp, LengthUnit)
    !=======================================================================
    ! This routine outputs in the GADGET format positions and velocities
    ! of particles in a file.
    !
    ! INPUTS :
    ! --------
    ! fileout (string)      : name of the output file
    ! npart (integer)       : number of particles 
    ! verbose (logical)     : verbose mode (.true. or .false.)
    ! megaverbose (logical) : detailed verbose mode
    ! Lbox (real 4)         : size of the box in Mpc/h
    ! hubble (real 4)       : H0/100 where H0 is the present time Hubble 
    !                         constant expressed in km/s/Mpc
    ! omega0                : the density parameter
    ! omegaL                : the cosmological constant
    ! aexp                  : the value of the expansion factor normalized
    !                         to unity at present time
    ! LengthUnit (optional) : length unit for the positions and the box size
    !                         WRITTEN TO THE FILE. 'Mpc/h' (default) keeps the
    !                         Mpc/h that (MG-)GLAM works in; 'kpc/h' applies the
    !                         x1000 that stock Gadget snapshots conventionally
    !                         use. This only rescales the output -- pos() is
    !                         always supplied in Mpc/h. Pick to match whatever
    !                         reads the file: a reader assuming the wrong one is
    !                         off by exactly 1000 and will not complain.
    ! pos(3,npart)          : coordinates of the particles in Mpc/h
    !                         in [0,Lbox]
    ! vel(3,npart)          : velocities of particles in km/s
    !
    ! OUTPUT :
    ! --------
    ! A binary file at the GADGET format with name fileout.
    !
    ! NOTE :
    ! ------
    ! pos and vel are global module variables.
    !=======================================================================
    implicit none
    integer(kind=8) :: npart, nparttot
    integer :: numfiles

    character(len=*) :: fileout
    logical :: verbose, megaverbose
    real(kind=4) :: Lbox, hubble, omega0, omegaL, aexp
    character(len=*), optional :: LengthUnit
    character(len=8) :: unit_used
    !=======================================================================

    integer :: i,j
    integer :: NparticlesInThisFile(0:5), NparticlesTotal(0:5), numfiles_in
    real(kind=8) :: massarr_in(0:5), a_in, redshift_in
    real(kind=8) :: omega0_in, omegaL_in, hubble_in
    integer :: unused_in(64-6-12-2-2-1-1-6-1-1-2-2-2-2)
    integer :: flag_sfr_in,flag_feedback_in,flag_cooling_in
    real(kind=8) :: xLbox_in
    real(kind=8) :: facco,mass_in_kg,mass_in_sol,ctilde

    integer, parameter :: lin=10

    ! Physical constants (units : m s kg) ->
    real(kind=8), parameter :: critical_density= 1.8788d-26 ! h^2 kg/m^3
    real(kind=8), parameter :: mega_parsec=3.0857d22
    real(kind=8), parameter :: solar_mass=1.989d30
    real(kind=8), parameter :: light_speed=2.99792458d5 ! km/s
    ! <-

    !--- Output length unit. Default Mpc/h: (MG-)GLAM works in Mpc/h throughout
    !    and so does everything downstream of this converter here, so converting
    !    would just be an unrequested factor of 1000. Stock Gadget snapshots do
    !    conventionally hold kpc/h, hence the option.
    unit_used = 'Mpc/h'
    if (present(LengthUnit)) unit_used = LengthUnit

    if (verbose) write(*,*) 'Output file '//trim(fileout)
    open(unit=lin, file=fileout, form='unformatted', status='unknown', err=1)

    NparticlesInThisFile(0:5) = 0
    NparticlesInThisFile(1) = npart
    massarr_in(0:5) = 0.0d0
    hubble_in = hubble
    ! Lbox is in Mpc/h, so the physical box side is Lbox/hubble Mpc and the
    ! h^2 of critical_density cancels against two of the three 1/h's:
    !   m [Msun]   = omega0 * rho_crit*h^2 * (Lbox/h)^3 / nparttot
    !   m [Msun/h] = m * h  =>  net h-dependence cancels entirely, as it must
    !                           for a mass expressed per-h from per-h inputs.
    ! The original code treated Lbox as Mpc and multiplied by hubble at the
    ! end, writing a mass table exactly h^3 too small (verified numerically:
    ! 15.672 vs 32.443 for the h=0.785 test run). Fixed 2026-07-29.
    mass_in_kg = (DBLE(Lbox)**3/DBLE(nparttot))*mega_parsec**3 &
    &             *omega0*critical_density*hubble_in**2
    mass_in_sol = mass_in_kg/solar_mass
    massarr_in(1) = mass_in_sol/1.0d10/hubble_in**2   ! 1e10 Msun/h; = *h with the h^3 of Lbox^3 undone
    a_in = aexp
    redshift_in = 1.0d0/a_in-1.0d0
    flag_sfr_in = 0
    flag_feedback_in = 0
    NparticlesTotal(0:5) = 0
    NparticlesTotal(1) = nparttot
    flag_cooling_in = 0
    numfiles_in = numfiles
    !--- facco converts the Mpc/h that pos() and Lbox arrive in to the requested
    !    output unit. Note massarr_in above is computed from Lbox in Mpc/h and
    !    is unit-independent (1e10 Msun/h either way), so it must stay above
    !    this rescaling.
    if (trim(unit_used) == 'kpc/h') then
        facco = 1000.d0
    else if (trim(unit_used) == 'Mpc/h') then
        facco = 1.d0
    else
        write(*,*) 'ERROR in write_gadget: unknown LengthUnit "'//trim(unit_used)//'"'
        write(*,*) '  expected "Mpc/h" or "kpc/h"'
        STOP 65
    endif
    if (verbose) write(*,*) 'Positions and box size written in comoving '//trim(unit_used)

    xLbox_in = Lbox * facco ! box size in the requested unit
    omega0_in = omega0
    omegaL_in = omegaL
    unused_in = 0
    if (megaverbose) write(*,*) 'Renormalize positions'
    pos(1:3,1:npart) = pos(1:3, 1:npart) * facco
    if (megaverbose) write(*,*) 'Renormalize velocities'
      vel(1:3,1:npart)=vel(1:3,1:npart)/sqrt(a_in)
    ! vel(1:3,1:npart) = vel(1:3,1:npart) / a_in * Lbox * 100.0 * hubble_in / sqrt(a_in)

    write(lin) (NparticlesInThisFile(i),i=0,5), (massarr_in(i),i=0,5), a_in,  &
    &             redshift_in, flag_sfr_in, flag_feedback_in,  &
    &             (NparticlesTotal(i),i=0,5), flag_cooling_in, numfiles_in, &
    &             xLbox_in,omega0_in,omegaL_in,hubble_in, &
    &             unused_in

    if (megaverbose) &
    &  write(*,*) (NparticlesInThisFile(i),i=0,5), (massarr_in(i),i=0,5), a_in,  &
    &             redshift_in, flag_sfr_in, flag_feedback_in,  &
    &             (NparticlesTotal(i),i=0,5), flag_cooling_in, numfiles_in, &
    &             xLbox_in,omega0_in,omegaL_in,hubble_in, &
    &             unused_in

    if (megaverbose) write(*,*) 'Output positions'
    if (megaverbose) then
        write(*, *) 'x,y,z min = ', minval(pos)
        write(*, *) 'x,y,z max = ', maxval(pos)
    endif
    write(lin) ((pos(i,j), i=1,3), j=1,npart)
    if (megaverbose) write(*,*) 'Output velocities'
    !      if (megaverbose) then
    write(*,*) 'vx,vy,vz min=',minval(vel)
    write(*,*) 'vx,vy,vz max=',maxval(vel)
    !      endif
    write(lin) ((vel(i,j),i=1,3),j=1,npart)

    if (megaverbose) write(*,*) 'Output particle identities'
    write(lin) (pid(j),j=1,npart)
    if (megaverbose) then
    write(*,*) 'pid min=',minval(pid)
    write(*,*) 'pid max=',maxval(pid)
    endif
    close(lin)

    if (verbose) write(*,*) 'The output file has been successfully written.'
    return

1   write(*,*) 'ERROR in write_gadget : I cannot open the output file'
    STOP

    end subroutine write_gadget

end module WriteGadgetFormat


! PROGRAM test

!     use WriteGadgetFormat

!     Allocate(pos(1:3, 1:1024), vel(1:3, 1:1024), pid(1024))
!     pos(:, :) = 1
!     vel(:, :) = 1
!     pid(pid) = 1
    
!     CALL write_gadget(fileout='./test.gadget', npart=1024, nparttot=1024, &
!         verbose=.True., megaverbose=.True., &
!         Lbox=1024.0, numfiles=1, hubble=0.7, omega0=0.3, omegaL=0.4, aexp=1.0)

! END PROGRAM test