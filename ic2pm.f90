!----------------------------------------------------------------------------
!
!   ic2pm : ingest on-disk 2LPTic (Gadget format-1) initial conditions into
!           MG-GLAM PM files (PMcrd.DAT + PMcrs0.DAT, PMcrs1.DAT, ...), so a
!           run can start from pre-made 2LPTic ICs instead of PMP2start's own
!           Zel'dovich generator. Reuses Tools::ReadSetup and Tools::WriteDataPM
!           so the output is byte-identical to what PMP2MG.exe / PMP2main.exe
!           read via ReadDataPM(-1).
!
!   Run from a Run<box>/ directory: reads ../Setup.dat, writes PM files into '.'
!   (the same 2-level directory contract as PMP2start.exe).
!
!   Usage:
!       ic2pm.exe  <IC_basename>  [S_vel]  [half|sync]
!
!   - <IC_basename> is the Gadget file prefix INCLUDING the trailing '.', so the
!     files are <IC_basename>0, <IC_basename>1, ...  e.g.
!       /cosma7/data/dp004/bl267/Runs/DEGRACE/ICs/IC_data/L1024/Node_002/ics.
!       /cosma8/data/dp203/bl267/Projects/Ongoing/HEFT/ICs/IC_highres/IC_Np1d_2048_L_1024_2LPT.
!   - [S_vel] (default 1.0) scales the Gadget VEL-block value to the standard
!     Gadget unit u = v_pec/sqrt(a) in km/s (verified values, see
!     2LPTIC_Gui/README.md and VALIDATION.md):
!       ICs from our fixed FML build (2LPTIC_Gui)        :  1.0
!       original DEGRACE 2LPTic ics.* (bl267, L1024)     :  1.0
!       Gui's old HEFT IC_Np1d_2048_L_1024_2LPT.* files  :  5.168609e6
!         (= 5.12e6/0.99059529; write the NUMBER on the command line, a '/'
!          expression is rejected)
!   - [half|sync] (default half) selects the epoch of the OUTPUT velocities:
!       half : velocities are moved back half a time step, to a_v = a_init -
!              ASTEP/2, where ASTEP is the value written to the PM header
!              (= ASTEP0 from ../Setup.dat).  This is the convention GLAM's
!              kick-then-drift leapfrog expects (PMP2main::MOVE kicks the
!              momenta from a-da/2 to a+da/2 with the force at a, then drifts
!              the positions with the new momenta) and the one PMP2start
!              itself uses (AEXPV = AEXPN - ASTEP/2, VCONS built at AEXPV,
!              XCONS at AEXPN).  2LPTic/Gadget ICs are synchronous (positions
!              and velocities both at a_init), so without this shift the first
!              kick over-boosts every momentum by ~0.75*da/a_init and P(k)
!              ends up ~1.2-4.8% high at linear k (growing-mode amplitude
!              ~0.6-2.3%; da = 4e-4 .. 1.6e-3, z_i = 49; larger at nonlinear k).
!              The shift is the growing-mode rescale used by PMP2start,
!                 V(a_v)/V(a) = (a_v/a)^1.5 * F(a_v)/F(a),  F = sqrt(Om+OmL a^3)
!              (D ~ a, f ~ 1 at z_init; the 2LPT velocity term strictly scales
!              as (a_v/a)^2.5, but its rms is ~1.1% of the 1LPT term on the
!              ICs checked, so the mis-scaling is <= 4.5e-4 for da <= 1.6e-3).
!       sync : no shift.  PMcrs*.DAT are byte-identical to the output of ic2pm
!              at commit 00ea3df (the last version without the epoch argument);
!              PMcrd.DAT is identical except the AEXP0 word (0-based bytes
!              53-56), an uninitialised local of Tools::WriteDataPM that
!              differs between any two runs, even of the same binary.
!              Only for reproducing runs converted before commit ee44100.
!              Those runs were made when sync was the only behaviour and their
!              submit scripts pass NO third argument, so re-running them with
!              this binary would silently give half: append 'sync' to
!                conv_da{4,8,16}/Run1/submit_conv.sh, ic2pm_val_L1024/Run1/submit_val.sh,
!                fid2LPTIC_L512Np2048Ng4096{,_da4,_da6}/Run1/submit.sh
!              to reproduce them.  The PM header string records the mode:
!                'ic2pm: 2LPTic ingest'               (sync, and all pre-ee44100 files)
!                'ic2pm: 2LPTic ingest, v at a-da/2'  (half)
!   All failures print a message and exit with status 1; success exits 0.
!
!   Corrections vs gadget2pm.f90 (which targets z=0, kpc/h, BDM snapshots):
!     (1) positions are already Mpc/h (BoxSize=1024) -> NO /1000; keep Box=1024.
!     (2) 64-bit particle count Nall + NallHighWord<<32 (2048^3 overflows Nall).
!     (3) velocity factor  S_vel * a^1.5 * NGRID/(100*Box)  (Gadget v_pec=sqrt(a)*u,
!         GLAM V=v_pec*a*NGRID/(100*Box)); gadget2pm's a^1.0 is only valid at a=1.
!     (4) ISTEP=0 (fresh IC); NGRID / cosmology come from ../Setup.dat.
!----------------------------------------------------------------------------
Program IC2PM
   use Tools
   implicit none

   type :: GadgetHeader
      sequence
      integer*4 :: npart(6)
      real*8    :: massarr(6)
      real*8    :: time, redshift
      integer*4 :: flag_sfr, flag_feedback
      integer*4 :: nall(6)
      integer*4 :: flag_cooling, num_files
      real*8    :: BoxSize, Omega0, OmegaLambda, HubbleParam
      integer*4 :: flag_stellarage, flag_metals
      integer*4 :: nallhw(6)                 ! NallHighWord @ header byte offset 168
      character(len=64) :: fill              ! pad to 256 bytes total
   end type GadgetHeader

   type(GadgetHeader) :: gh
   character(len=256)  :: inbase, sarg, fname
   character(len=16)   :: vepoch
   integer*4           :: nfiles_in, ifile_g, np, j, ios
   character(len=*), parameter :: usage = &
      'Usage: ic2pm.exe <IC_basename incl trailing "."> [S_vel] [half|sync]'
   integer*8           :: ip, ioff, nlow
   real*4              :: xs, vfac, Svel, xx, yy, zz
   real*4              :: xmin, xmax
   real*8              :: velsq, a8, av8, fshift
   real*4, allocatable :: pos(:), vel(:)
   integer, parameter  :: uG = 50
   logical             :: ex

!--- command line -----------------------------------------------------------
   if (command_argument_count() < 1) then
      write(*,*) usage
      stop 1
   end if
   call get_command_argument(1, inbase)
   Svel = 1.0
   if (command_argument_count() >= 2) then
      call get_command_argument(2, sarg); sarg = adjustl(sarg)
      if (len_trim(sarg) == 0 .or. verify(trim(sarg), '0123456789.+-eEdD') /= 0) then
         write(*,*) ' ic2pm: bad S_vel argument "', trim(sarg), &
                    '" (must be a plain number, e.g. 1.0 or 5.168609e6)'
         write(*,*) usage
         stop 1
      end if
      read(sarg, *, iostat=ios) Svel
      if (ios /= 0 .or. Svel <= 0.0) then
         write(*,*) ' ic2pm: bad S_vel argument "', trim(sarg), '" (unreadable or <= 0)'
         write(*,*) usage
         stop 1
      end if
   end if
   vepoch = 'half'
   if (command_argument_count() >= 3) then
      call get_command_argument(3, sarg); vepoch = adjustl(sarg)
   end if
   if (trim(vepoch) /= 'half' .and. trim(vepoch) /= 'sync') then
      write(*,*) ' ic2pm: bad velocity-epoch argument "', trim(vepoch), '"'
      write(*,*) usage
      stop 1
   end if

!--- run parameters from ../Setup.dat (NGRID, NROW, Box, cosmology, AEXPN0) ---
   inquire(file='../Setup.dat', exist=ex)
   if (.not. ex) then
      write(*,*) ' ic2pm: ../Setup.dat not found (run PMP2init first, from RunN/)'
      stop 1
   end if
   open(11, file='../Setup.dat', status='old')
   call ReadSetup             ! sets AEXPN0, ASTEP0, Box, hubble, Om, OmL, NROW, NGRID, Nseed, ...

!--- file-0 Gadget header: a, Box, cosmology, total particle count -----------
   write(fname,'(2a)') trim(inbase), '0'
   inquire(file=trim(fname), exist=ex)
   if (.not. ex) then
      write(*,*) ' ic2pm: first IC file (<IC_basename>0) not found: ', trim(fname)
      stop 1
   end if
   open(uG, file=trim(fname), form='unformatted', access='sequential', &
        status='old', convert='little_endian')
   read(uG) gh
   close(uG)

   nfiles_in = gh%num_files
   nlow = int(gh%nall(2), 8)
   if (nlow < 0) nlow = nlow + 4294967296_8              ! unsigned low 32 bits
   Nparticles = nlow + ishft(int(gh%nallhw(2),8), 32)    ! 64-bit total count

!--- consistency: Gadget header vs Setup.dat (abort before a multi-hour run) --
   if (abs(Box - real(gh%BoxSize)) > 1.0e-2*Box) &
      call die('Box mismatch (Setup vs Gadget header)', Box, real(gh%BoxSize))
   if (abs(AEXPN0 - real(gh%time)) > 1.0e-4) &
      call die('AEXPN0=1/(1+z_init) vs Gadget a mismatch', AEXPN0, real(gh%time))
   if (nint(real(Nparticles,8)**(1.d0/3.d0)) /= NROW) then
      write(*,*) ' NROW(Setup)=', NROW, '  nint(Nparticles^1/3)=', &
                 nint(real(Nparticles,8)**(1.d0/3.d0))
      write(*,*) ' ic2pm: NROW mismatch (Setup vs IC)'
      stop 1
   end if
   if (abs(Om  - real(gh%Omega0))      > 1.0e-3) call die('Om mismatch',  Om,  real(gh%Omega0))
   if (abs(OmL - real(gh%OmegaLambda)) > 1.0e-3) call die('OmL mismatch', OmL, real(gh%OmegaLambda))
   if (abs(hubble - real(gh%HubbleParam)) > 1.0e-3) &
      call die('h mismatch', hubble, real(gh%HubbleParam))

!--- PM header scalars (Box/NROW/NGRID/Om/OmL/hubble/Nseed already from Setup) -
   AEXPN        = real(gh%time)         ! = AEXPN0 = 1/(1+z_init)  (e.g. 0.02 at z=49)
   ASTEP        = ASTEP0
   ISTEP        = 0                      ! fresh IC (not nint(z))
   Nrealization = 1
   EKIN = 0.0 ; EKIN1 = 0.0 ; EKIN2 = 0.0 ; TINTG = 0.0 ; AEU0 = 0.0
   extras(:)    = 0.0
   extras(100)  = Box                    ! sole Box source for ReadDataPM
   if (trim(vepoch) == 'half') then
      write(HEADER,'(a)') 'ic2pm: 2LPTic ingest, v at a-da/2'   ! records the epoch (34 chars)
   else
      write(HEADER,'(a)') 'ic2pm: 2LPTic ingest'                ! unchanged from 00ea3df
   end if

   xs   = real(NGRID) / Box
   vfac = Svel * AEXPN**1.5 * real(NGRID) / (100.0 * Box)   ! synchronous, at a_init

!--- half-step velocity epoch (GLAM leapfrog convention, see header) ---------
   a8  = real(AEXPN, 8)
   av8 = a8
   fshift = 1.0d0
   if (trim(vepoch) == 'half') then
      if (ASTEP <= 0.0 .or. ASTEP >= AEXPN) &
         call die('ASTEP (Setup.dat) must satisfy 0 < ASTEP < a_init', ASTEP, AEXPN)
      av8    = a8 - real(ASTEP, 8)/2.0d0
      fshift = (av8/a8)**1.5d0 * sqrt((real(Om,8) + real(OmL,8)*av8**3) / &
                                      (real(Om,8) + real(OmL,8)*a8**3))
      vfac   = real(real(vfac, 8) * fshift, 4)
   end if

   write(*,'(a)')          ' === ic2pm : 2LPTic -> MG-GLAM PM ==='
   write(*,'(2a)')         '  IC basename = ', trim(inbase)
   write(*,'(a,es12.4)')   '  S_vel       = ', Svel
   write(*,'(a,f12.4)')    '  Box [Mpc/h] = ', Box
   write(*,'(3(a,f9.5))')  '  Om =', Om, '  OmL =', OmL, '  h =', hubble
   write(*,'(a,f10.6)')    '  a (=AEXPN)  = ', AEXPN
   write(*,'(a,i14)')      '  Nparticles  = ', Nparticles
   write(*,'(2(a,i7))')    '  NROW =', NROW, '  NGRID =', NGRID
   write(*,'(a,i6)')       '  num_files   = ', nfiles_in
   write(*,'(2a)')         '  velocity epoch = ', trim(vepoch)
   write(*,'(a,f12.8,a,es12.4)') '  a_v (velocities) = ', av8, '   ASTEP =', ASTEP
   write(*,'(a,f12.8)')    '  velocity shift factor V(a_v)/V(a_init) = ', fshift
   write(*,'(2(a,es12.4))')'  xs =', xs, '   vfac =', vfac

!--- allocate the Tools global particle arrays -------------------------------
   allocate(XPAR(Nparticles), YPAR(Nparticles), ZPAR(Nparticles))
   allocate(VX(Nparticles),   VY(Nparticles),   VZ(Nparticles))

!--- loop over the Gadget files: scatter type-1 particles --------------------
   ioff  = 0
   velsq = 0.0d0
   do ifile_g = 0, nfiles_in - 1
      write(fname,'(a,i0)') trim(inbase), ifile_g
      open(uG, file=trim(fname), form='unformatted', access='sequential', &
           status='old', convert='little_endian')
      read(uG) gh
      np = gh%npart(2)
      allocate(pos(3*np), vel(3*np))
      read(uG) pos                                ! positions, Mpc/h
      read(uG) vel                                ! velocities, file units (IDs skipped)
      close(uG)

      if (ioff + np > Nparticles) then
         write(*,*) ' ic2pm: too many particles vs header count'
         stop 1
      end if

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j,xx,yy,zz) REDUCTION(+:velsq)
      do j = 1, np
         xx = pos(3*j-2) * xs + 1.0               ! Mpc/h -> grid units (NO /1000)
         yy = pos(3*j-1) * xs + 1.0
         zz = pos(3*j  ) * xs + 1.0
         if (xx >= NGRID+1.0) xx = xx - NGRID      ! periodic wrap -> [1, NGRID+1)
         if (yy >= NGRID+1.0) yy = yy - NGRID
         if (zz >= NGRID+1.0) zz = zz - NGRID
         if (xx < 1.0) xx = xx + NGRID
         if (yy < 1.0) yy = yy + NGRID
         if (zz < 1.0) zz = zz + NGRID
         XPAR(ioff+j) = xx
         YPAR(ioff+j) = yy
         ZPAR(ioff+j) = zz
         VX(ioff+j)   = vel(3*j-2) * vfac
         VY(ioff+j)   = vel(3*j-1) * vfac
         VZ(ioff+j)   = vel(3*j  ) * vfac
         velsq = velsq + real(vel(3*j-2),8)**2 + real(vel(3*j-1),8)**2 + real(vel(3*j),8)**2
      end do
!$OMP END PARALLEL DO

      ioff = ioff + np
      deallocate(pos, vel)
      if (mod(ifile_g,50)==0 .or. ifile_g==nfiles_in-1) &
         write(*,'(a,i5,a,i14)') '  read file ', ifile_g, '   cumulative np = ', ioff
   end do

   if (ioff /= Nparticles) then
      write(*,'(a,2i15)') ' ic2pm: FATAL scattered /= total count :', ioff, Nparticles
      stop 1
   end if

!--- diagnostics: coordinate range + implied 1D peculiar-velocity rms --------
   xmin = 1.0e30 ; xmax = -1.0e30
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ip) REDUCTION(min:xmin) REDUCTION(max:xmax)
   do ip = 1, Nparticles
      xmin = min(xmin, XPAR(ip), YPAR(ip), ZPAR(ip))
      xmax = max(xmax, XPAR(ip), YPAR(ip), ZPAR(ip))
   end do
!$OMP END PARALLEL DO
   write(*,'(a,2f12.4,a,i8,a)') '  coord min/max =', xmin, xmax, &
        '   (expect [1,', NGRID+1, '))'
   write(*,'(a,f10.3,a)') '  1D peculiar-velocity rms at a_init (IC file) = ', &
        sqrt(AEXPN) * Svel * sqrt(velsq/(3.0d0*real(Nparticles,8))), &
        ' km/s   (expect ~50 km/s at z=49)'
   ! v_pec(a_v)/v_pec(a) = fshift * a/a_v  (GLAM momentum V = v_pec*a*NGRID/(100*Box))
   write(*,'(a,f10.3,a,f10.6,a,f10.6,a)') '  1D peculiar-velocity rms written, at a_v = ', &
        sqrt(AEXPN) * Svel * sqrt(velsq/(3.0d0*real(Nparticles,8))) * fshift * a8/av8, &
        ' km/s   (a_v =', av8, ', momentum factor =', fshift, ')'

!--- write PM files (PMcrd.DAT + PMcrs0.DAT, PMcrs1.DAT, ...) into cwd --------
   write(*,'(a)') '  writing PMcrd.DAT + PMcrs*.DAT ...'
   call WriteDataPM(0, '')
   write(*,'(a)') '  ic2pm done.'

contains
   subroutine die(msg, a, b)
      character(len=*), intent(in) :: msg
      real*4,           intent(in) :: a, b
      write(*,'(3a,2es14.6)') ' ic2pm: FATAL ', msg, ' : ', a, b
      stop 1
   end subroutine die

end Program IC2PM
