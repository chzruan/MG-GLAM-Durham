!--------------------------------------------------
!
!     Module Rockstar: in-line Rockstar halo catalogues at full outputs
!
!     The halo finding is the C/C++ adapter in rockstar_inline/ (GNU GPL v3 or
!     later, it links MPI-Rockstar), built into GLAM with
!     make ROCKSTAR_DIR=<MPI-Rockstar source>; a default build links its stub,
!     which has no Rockstar code, and GLAM then refuses Find Rockstar halos = 1
!     at start-up. Each call reads GLAM's particle arrays only and publishes
!     CATALOGS/Rockstar.<step>.<box>.bin + .meta (distinct haloes; see
!     rockstar_inline/README.md). Nothing in the adapter stops the run: a
!     failure or a failed memory pre-check gives a catalogue whose .meta says
!     why, and the simulation continues.
!
!--------------------------------------------------
Module Rockstar
   use, intrinsic :: iso_c_binding
   implicit none
   private
   public :: RockstarHalos, CheckRockstarConfig

   character(len=*), parameter :: ConfigFile = 'Rockstar.config'   ! in the run directory, as BDM.config
   character(len=*), parameter :: CatalogDir = 'CATALOGS'

   integer, parameter :: RS_OK = 0, RS_SKIPPED = 1, RS_FAILED = 2, RS_NOT_BUILT = 3

   type, bind(C) :: rs_inline_info                                   ! struct rs_inline_info of rs_inline.h
      integer(c_int)         :: status
      character(kind=c_char) :: reason(512)
      integer(c_int64_t)     :: ngroups, nmembers, largest, nhalos_internal, nhalos, nhosts
      real(c_double)         :: t_fof, t_gather, t_kernel, t_strict, t_parents, t_write, t_total
      real(c_double)         :: rss_entry_gib, rss_max_gib, need_fof_gib, need_kernel_gib
   end type rs_inline_info

   interface
      integer(c_int) function rs_inline_available() bind(C, name='rs_inline_available')
         import :: c_int
      end function rs_inline_available

      type(c_ptr) function rs_inline_version() bind(C, name='rs_inline_version')
         import :: c_ptr
      end function rs_inline_version

      integer(c_int) function rs_inline_check_config(path, msg, msglen) bind(C, name='rs_inline_check_config')
         import :: c_int, c_char
         character(kind=c_char), dimension(*), intent(in)  :: path
         character(kind=c_char), dimension(*), intent(out) :: msg
         integer(c_int), value :: msglen
      end function rs_inline_check_config

      integer(c_int) function rs_inline_run(n, x, y, z, vx, vy, vz, box, ngrid, aexpn, astep, om, oml, hubble, &
                                            mg_flag, mg_model, fr0, header, step, box_index, config_path,     &
                                            catalog_dir, info) bind(C, name='rs_inline_run')
         import :: c_int, c_int64_t, c_float, c_double, c_char, rs_inline_info
         integer(c_int64_t), value :: n, ngrid
         real(c_float), dimension(*), intent(in) :: x, y, z, vx, vy, vz
         real(c_double), value :: box, aexpn, astep, om, oml, hubble, fr0
         integer(c_int), value :: mg_flag, mg_model, step, box_index
         character(kind=c_char), dimension(*), intent(in) :: header, config_path, catalog_dir
         type(rs_inline_info), intent(out) :: info
      end function rs_inline_run
   end interface

CONTAINS

!--------------------------------------------
!     a NUL-terminated C string (array or pointer) as a Fortran string
!--------------------------------------------
function c_text(s) result(t)
   character(kind=c_char), dimension(:), intent(in) :: s
   character(len=:), allocatable :: t
   integer :: i, n
   n = 0
   do i = 1, size(s)
      if (s(i) == c_null_char) exit
      n = i
   end do
   allocate(character(len=n) :: t)
   do i = 1, n
      t(i:i) = s(i)
   end do
end function c_text

!--------------------------------------------
!     stop the run with exit status 1 (a plain STOP returns 0, which a job script
!     cannot tell from success); the message goes to stdout (the log) and stderr
!--------------------------------------------
subroutine Fatal(msg)
   use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
   character(len=*), intent(in) :: msg
   write (output_unit, '(2a)') ' Rockstar: ', trim(msg)
   flush (output_unit)
   write (error_unit, '(2a)') ' Rockstar: ', trim(msg)
   error stop 1
end subroutine Fatal

function version_text() result(t)
   character(len=:), allocatable :: t
   character(kind=c_char), dimension(:), pointer :: p
   call c_f_pointer(rs_inline_version(), p, [1024])
   t = c_text(p)
end function version_text

!--------------------------------------------
!     Start-up check (Find Rockstar halos = 1): the build has the adapter and
!     Rockstar.config is usable; the adapter keeps the configuration it read
!     here for every output of the run.
!--------------------------------------------
subroutine CheckRockstarConfig
   character(kind=c_char) :: msg(512)
   logical :: exists
   if (rs_inline_available() == 0) &
      call Fatal('Find Rockstar halos = 1, but this executable has no in-line Rockstar: rebuild GLAM with '// &
                 'make ROCKSTAR_DIR=<MPI-Rockstar source> (rockstar_inline/README.md)')
   inquire (file=ConfigFile, exist=exists)
   if (.not. exists) call Fatal('Find Rockstar halos = 1, but '//ConfigFile//' is not in the run directory')
   if (rs_inline_check_config(ConfigFile//c_null_char, msg, int(size(msg), c_int)) /= 0) call Fatal(c_text(msg))
   write (*, '(4a)') ' ', ConfigFile, ' validated at start-up; in-line ', version_text()
end subroutine CheckRockstarConfig

!--------------------------------------------
!     The catalogue of one full output, called by Analysis after BDM. The PM
!     density FI is not needed again before DENSIT at the next step, so its
!     memory goes to Rockstar: after BDM it is allocated but untouched (not
!     resident); without BDM it is freed here. The particle arrays are only
!     read; with GLAM_ROCKSTAR_CHECKSUM=1 in the environment their checksums
!     before and after the call are compared (a test of that promise).
!--------------------------------------------
subroutine RockstarHalos(moment)
   use Tools, only: Nparticles, XPAR, YPAR, ZPAR, VX, VY, VZ, FI, NGRID, Box, AEXPN, ASTEP, Om, OmL, hubble, &
                    MG_flag, MG_model, fR0, HEADER, Nrealization, Memory
   integer, intent(in) :: moment
   type(rs_inline_info) :: info
   integer(c_int) :: rc
   integer*8 :: chk0(6), chk1(6)
   character(len=8) :: env
   logical :: check
   real*4 :: mem

   if (.not. (allocated(XPAR) .and. allocated(VX))) then
      write (*, '(a)') ' Rockstar: particle positions or velocities not allocated; no catalogue at this step'
      return
   end if
   if (allocated(FI)) then
      mem = Memory(-1_8*NGRID*NGRID*NGRID)
      deallocate (FI)
   end if
   call get_environment_variable('GLAM_ROCKSTAR_CHECKSUM', env)
   check = trim(env) == '1'
   if (check) call ParticleChecksum(chk0)

   rc = rs_inline_run(int(Nparticles, c_int64_t), XPAR, YPAR, ZPAR, VX, VY, VZ, real(Box, c_double),        &
                      int(NGRID, c_int64_t), real(AEXPN, c_double), real(ASTEP, c_double), real(Om, c_double), &
                      real(OmL, c_double), real(hubble, c_double), int(MG_flag, c_int), int(MG_model, c_int),  &
                      real(fR0, c_double), trim(HEADER)//c_null_char, int(moment, c_int),                      &
                      int(Nrealization, c_int), ConfigFile//c_null_char, CatalogDir//c_null_char, info)

   if (check) then
      call ParticleChecksum(chk1)
      if (any(chk1 /= chk0)) then
         write (*, '(a,6z17)') ' Rockstar: particle checksums changed across the call: before ', chk0
         write (*, '(a,6z17)') '                                                        after  ', chk1
         call Fatal('the particle arrays changed across the call')
      end if
      write (*, '(a,6z17)') ' Rockstar: particle checksums unchanged across the call: ', chk0
   end if
   if (rc /= RS_OK) write (*, '(a,i5,a,i2,2a)') ' Rockstar step', moment, ': no catalogue rows (status ', rc, &
                                              '): ', c_text(info%reason)

   mem = Memory(1_8*NGRID*NGRID*NGRID)
   allocate (FI(NGRID, NGRID, NGRID))
end subroutine RockstarHalos

!--------------------------------------------
!     Order-sensitive checksums of the six particle arrays (bit patterns)
!--------------------------------------------
subroutine ParticleChecksum(c)
   use Tools, only: Nparticles, XPAR, YPAR, ZPAR, VX, VY, VZ
   integer*8, intent(out) :: c(6)
   integer*8 :: i, s1, s2, s3, s4, s5, s6, w
   s1 = 0; s2 = 0; s3 = 0; s4 = 0; s5 = 0; s6 = 0
   !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, w) REDUCTION(+:s1, s2, s3, s4, s5, s6)
   do i = 1, Nparticles
      w = 2*i + 1
      s1 = s1 + w*transfer(XPAR(i), 0)
      s2 = s2 + w*transfer(YPAR(i), 0)
      s3 = s3 + w*transfer(ZPAR(i), 0)
      s4 = s4 + w*transfer(VX(i), 0)
      s5 = s5 + w*transfer(VY(i), 0)
      s6 = s6 + w*transfer(VZ(i), 0)
   end do
   !$OMP END PARALLEL DO
   c = [s1, s2, s3, s4, s5, s6]
end subroutine ParticleChecksum

end Module Rockstar
