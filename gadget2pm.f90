!----------------------------------------------------------------------------
!
!   gadget2pm : convert a Gadget-2 (format-1) snapshot into MG-GLAM PM files
!               (PMcrd.DAT + PMcrs0.DAT, PMcrs1.DAT, ...) so the BDM halo
!               finder (PMP2BDM.exe) can read it.
!
!   Reuses the Tools module and its WriteDataPM routine, so the output is
!   byte-for-byte the format ReadDataPM expects (same big-endian compile).
!
!   Usage:
!       gadget2pm.exe  <gadget_input_dir>  <pm_output_dir>  <NGRID>
!
!   - <gadget_input_dir> contains gadget.0 ... gadget.<numfiles-1>
!   - PM files are written into <pm_output_dir> (must exist)
!   - <NGRID> is the BDM density mesh (e.g. 2048). Coordinates are written
!     in grid units [1, NGRID+1); velocities in BDM grid units.
!
!   Assumptions (verified for the DEGRACE z=0 snapshots):
!     * Gadget format-1, little-endian, dark-matter particles are type 1,
!       masses come from the header mass table, block order HEAD/POS/VEL/ID.
!     * a = 1 (z=0), so peculiar velocity [km/s] = stored Gadget velocity.
!----------------------------------------------------------------------------
Program Gadget2PM
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
      character(len=96) :: fill
   end type GadgetHeader

   type(GadgetHeader) :: gh
   character(len=256) :: inbase, outdir, sarg, fname
   integer*4          :: ngrid_in, nfiles_in, ifile, np, j
   integer*8          :: ip, ioff
   real*4             :: xs, vfac, xx, yy, zz
   real*4             :: xmin, xmax
   real*4, allocatable :: pos(:), vel(:)
   integer, parameter :: uG = 50

!--- command line -----------------------------------------------------------
   if (command_argument_count() < 3) then
      write(*,*) 'Usage: gadget2pm.exe <gadget_input_dir> <pm_output_dir> <NGRID>'
      stop
   end if
   call get_command_argument(1, inbase)
   call get_command_argument(2, outdir)
   call get_command_argument(3, sarg);  read(sarg,*) ngrid_in

!--- read the file-0 header for the global simulation parameters ------------
   write(fname,'(2a)') trim(inbase), '/gadget.0'
   open(uG, file=trim(fname), form='unformatted', access='sequential', &
        status='old', convert='little_endian')
   read(uG) gh
   close(uG)

   Box         = real(gh%BoxSize) / 1000.0      ! kpc/h -> Mpc/h
   Om          = real(gh%Omega0)
   OmL         = real(gh%OmegaLambda)
   hubble      = real(gh%HubbleParam)
   AEXPN       = real(gh%time)
   Nparticles  = int(gh%nall(2), 8)
   nfiles_in   = gh%num_files
   NGRID       = ngrid_in
   NROW        = nint(real(Nparticles)**(1.0/3.0))
   Nrealization = 1
   Nseed       = 0
   ISTEP       = nint(gh%redshift)              ! 0 at z=0; only used as a tag
   AMPLT       = 0.0 ; ASTEP = 0.0
   EKIN = 0.0 ; EKIN1 = 0.0 ; EKIN2 = 0.0 ; TINTG = 0.0 ; AEU0 = 0.0
   extras(:)   = 0.0
   extras(100) = Box
   write(HEADER,'(a)') 'Gadget2PM converted DEGRACE snapshot'

   write(*,'(a)')          ' === gadget2pm ==='
   write(*,'(2a)')         '  input dir   = ', trim(inbase)
   write(*,'(2a)')         '  output dir  = ', trim(outdir)
   write(*,'(a,f12.4)')    '  Box [Mpc/h] = ', Box
   write(*,'(3(a,f8.5))')  '  Om =', Om, '  OmL =', OmL, '  h =', hubble
   write(*,'(a,f8.5)')     '  a           = ', AEXPN
   write(*,'(a,i12)')      '  Nparticles  = ', Nparticles
   write(*,'(2(a,i7))')    '  NROW =', NROW, '  NGRID =', NGRID
   write(*,'(a,i6)')       '  num_files   = ', nfiles_in

!--- allocate the global particle arrays (Tools module) ---------------------
   allocate(XPAR(Nparticles), YPAR(Nparticles), ZPAR(Nparticles))
   allocate(VX(Nparticles),   VY(Nparticles),   VZ(Nparticles))

   xs   = real(NGRID) / Box                     ! Mpc/h -> grid units
   vfac = AEXPN * real(NGRID) / (100.0 * Box)   ! km/s -> BDM grid velocity

!--- loop over the Gadget files, scatter type-1 particles into the arrays ----
   ioff = 0
   do ifile = 0, nfiles_in - 1
      write(fname,'(2a,i0)') trim(inbase), '/gadget.', ifile
      open(uG, file=trim(fname), form='unformatted', access='sequential', &
           status='old', convert='little_endian')
      read(uG) gh
      np = gh%npart(2)
      allocate(pos(3*np), vel(3*np))
      read(uG) pos                              ! positions, kpc/h
      read(uG) vel                              ! velocities, km/s (a=1)
      close(uG)

      if (ioff + np > Nparticles) stop ' Too many particles vs header nall'

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j,xx,yy,zz)
      do j = 1, np
         xx = pos(3*j-2)/1000.0 * xs + 1.0       ! kpc/h -> Mpc/h -> grid
         yy = pos(3*j-1)/1000.0 * xs + 1.0
         zz = pos(3*j  )/1000.0 * xs + 1.0
         if (xx >= NGRID+1.0) xx = xx - NGRID    ! periodic wrap -> [1,NGRID+1)
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
      end do
!$OMP END PARALLEL DO

      ioff = ioff + np
      deallocate(pos, vel)
      if (mod(ifile,50)==0 .or. ifile==nfiles_in-1) &
         write(*,'(a,i5,a,i12)') '  read file ', ifile, '  cumulative np = ', ioff
   end do

   if (ioff /= Nparticles) then
      write(*,'(a,2i14)') ' WARNING: scattered /= header nall :', ioff, Nparticles
   end if

!--- coordinate sanity check ------------------------------------------------
   xmin = 1.0e30 ; xmax = -1.0e30
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ip) REDUCTION(min:xmin) REDUCTION(max:xmax)
   do ip = 1, Nparticles
      xmin = min(xmin, XPAR(ip), YPAR(ip), ZPAR(ip))
      xmax = max(xmax, XPAR(ip), YPAR(ip), ZPAR(ip))
   end do
!$OMP END PARALLEL DO
   write(*,'(a,2f12.4,a,i6)') '  coord min/max =', xmin, xmax, '  expect [1,', NGRID+1
   write(*,'(a,f12.4)')       '  MassOne [1e?] = Om*2.774e11*(Box/NROW)^3 =', &
                              Om*2.774e11*(Box/real(NROW))**3

!--- write the PM files (PMcrd.DAT + PMcrs0.DAT, PMcrs1.DAT, ...) ------------
   write(*,'(a)') '  writing PM files ...'
   call WriteDataPM(0, trim(outdir)//'/')
   write(*,'(a)') '  done.'

end Program Gadget2PM
