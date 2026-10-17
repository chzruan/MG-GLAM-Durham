!-------------------------------------------------
!   Subsample: ID-hashed, quantised particle subsample for the
!   DEGRACE companion runs (item I9 of
!   docs/degrace_companion_implementation_plan.md).
!
!   Selection. Particle ID = array index - 1 (0-based; the particle
!   arrays are never permuted). For ICs ingested by ic2pm the index is
!   the Gadget ID (ic2pm places every particle at index = ID and tags
!   the lattice order in HEADER: IDx = x fastest, ID-1 = ix + Np (iy +
!   Np iz), the FML generator; IDz = z fastest, ID-1 = iz + Np (iy +
!   Np ix), 2LPTic). For native PMP2start ICs ID-1 = i + Np (j + Np k)
!   with i along x. Particle ID is kept when
!       splitmix64(ID xor SALT) mod M == 0
!   with SALT = 0x5EED00DE6C0A2026 (never change) and M the
!   subsampling factor passed by the caller (M = 1 keeps everything).
!   The hash lives in glam_hash.c (unsigned 64-bit arithmetic) and is
!   reached through iso_c_binding. IDs are not stored: the k-th
!   stored particle is the k-th ID that passes the rule, so a reader
!   recovers the IDs by re-evaluating the rule over 0 .. N_p-1
!   (reference: test/extras/check_subsample.py).
!
!   Quantisation. Positions: q = floor(x/Box * 65536) with x wrapped
!   into [0, Box), stored as int16 holding q - 32768; decode
!   x = (q + 0.5) Box/65536 (15.6 kpc/h at Box = 1024). Velocities:
!   peculiar v [km/s] = V_code * (100 Box/N_grid) / a_v with
!   a_v = AEXPN - ASTEP/2 (the leapfrog leaves the velocities half a
!   step behind the positions at every output), s = nint(v/dv),
!   dv = 0.25 km/s, clipped to [-32767, 32767]; the number of clipped
!   components is written to the .meta.
!
!   Files (run directory, same scheme as the other products):
!     Subsample.<step>.<box>.bin   raw little-endian int16, no record
!                                  markers: qx, qy, qz, vx, vy, vz,
!                                  each n_sel values, in that order
!                                  (12 bytes per particle)
!     Subsample.<step>.<box>.meta  plain text, key = value lines
!
!   Algorithm. Deterministic two-pass compaction over fixed chunks of
!   SubChunk consecutive IDs: pass 1 counts the kept IDs per chunk
!   (OpenMP over chunks), an exclusive prefix sum gives every chunk
!   its slot range, pass 2 re-evaluates the rule per chunk and fills
!   the slots in increasing ID order. The result does not depend on
!   the thread count or on SubChunk. Memory: one int16 array of
!   n_sel (one component at a time; the hash is cheap, so it is
!   re-evaluated for each of the six components) plus one int32
!   scratch array of SubChunk per thread. No per-particle mask.
!-------------------------------------------------
Module Subsample
    use Tools
    use iso_c_binding, only: c_int, c_int32_t, c_int64_t
    implicit none
    Integer(c_int64_t), parameter :: SubSalt = int(z'5EED00DE6C0A2026', c_int64_t)
    Integer*4, parameter :: SubNq = 65536            ! position quanta per box side
    Real*8, parameter    :: SubDv = 0.25d0           ! velocity quantum [km/s]
    Integer*8            :: SubChunk = 1048576_8     ! IDs per compaction chunk (output independent of it; tests shrink it)

    Interface
        ! uint64 splitmix64(uint64 x)
        Function splitmix64(x) bind(C, name='splitmix64')
            import :: c_int64_t
            Integer(c_int64_t), value :: x
            Integer(c_int64_t) :: splitmix64
        End Function splitmix64
        ! 1 if ID is kept, 0 otherwise
        Function glam_select(id, salt, m) bind(C, name='glam_select')
            import :: c_int, c_int64_t
            Integer(c_int64_t), value :: id, salt, m
            Integer(c_int) :: glam_select
        End Function glam_select
        ! number of kept IDs in [id0, id0+n)
        Function glam_count_block(id0, n, salt, m) bind(C, name='glam_count_block')
            import :: c_int64_t
            Integer(c_int64_t), value :: id0, n, salt, m
            Integer(c_int64_t) :: glam_count_block
        End Function glam_count_block
        ! offsets (0-based, from id0, increasing) of the kept IDs in [id0, id0+n) -> idx(1:count); returns count
        Function glam_select_indices(id0, n, salt, m, idx) bind(C, name='glam_select_indices')
            import :: c_int32_t, c_int64_t
            Integer(c_int64_t), value :: id0, n, salt, m
            Integer(c_int32_t), intent(out) :: idx(*)
            Integer(c_int64_t) :: glam_select_indices
        End Function glam_select_indices
    End Interface
Contains
!---------------------------------------------------------
!   Write Subsample.<moment>.<Nrealization>.{bin,meta} with every
!   M-th particle (by hash) of the current particle arrays.
!   Msub <= 0: nothing is written.
!---------------------------------------------------------
    Subroutine WriteSubsample(moment, Msub)
        Integer*4, intent(in) :: moment, Msub
        Integer*8 :: Nchunk, nsel, ic, id0, n
        Integer*8, allocatable :: cnt(:), off(:)
        Integer*2, allocatable :: q(:)
        Integer*8 :: nclip(3), nbad, ndum
        Integer*4 :: iu, ic4
        Integer*8 :: c0, c1, cr
        Real*8 :: av, xscl, vscl
        Character*120 :: Name

        If (Msub <= 0) Return
        If (Nparticles <= 0 .or. .not. allocated(XPAR)) Stop 'WriteSubsample: no particles in memory'
        If (SubChunk <= 0 .or. SubChunk > 2147483647_8) Stop 'WriteSubsample: SubChunk must be in 1 .. 2^31-1'
        Call system_clock(c0, cr)

        av   = dble(AEXPN) - dble(ASTEP)/2.d0             ! velocity epoch of this output
        xscl = dble(SubNq)/dble(NGRID)                     ! grid units (from 0) -> position quanta
        vscl = 100.d0*dble(Box)/dble(NGRID)/av/SubDv       ! code momentum -> velocity quanta

        !---- pass 1: kept IDs per chunk, exclusive prefix sum
        Nchunk = (Nparticles + SubChunk - 1)/SubChunk
        allocate (cnt(Nchunk), off(0:Nchunk))
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ic,id0,n) SCHEDULE(DYNAMIC,4)
        Do ic = 1, Nchunk
            id0 = (ic - 1)*SubChunk
            n = min(SubChunk, Nparticles - id0)
            cnt(ic) = glam_count_block(id0, n, SubSalt, int(Msub, c_int64_t))
        End Do
        off(0) = 0
        Do ic = 1, Nchunk
            off(ic) = off(ic - 1) + cnt(ic)
        End Do
        nsel = off(Nchunk)

        !---- pass 2, one component at a time: fill in ID order, write
        write (Name, '(2(a,i4.4),a)') 'Subsample.', moment, '.', Nrealization, '.bin'
        open (newunit=iu, file=TRIM(Name), status='replace', access='stream', form='unformatted', convert='little_endian')
        allocate (q(nsel))
        nclip = 0
        Do ic4 = 1, 6
            Select Case (ic4)
            Case (1)
                Call FillQuantised(XPAR, 1, xscl, q, nbad, ndum)
            Case (2)
                Call FillQuantised(YPAR, 1, xscl, q, nbad, ndum)
            Case (3)
                Call FillQuantised(ZPAR, 1, xscl, q, nbad, ndum)
            Case (4)
                Call FillQuantised(VX, 2, vscl, q, nbad, nclip(1))
            Case (5)
                Call FillQuantised(VY, 2, vscl, q, nbad, nclip(2))
            Case (6)
                Call FillQuantised(VZ, 2, vscl, q, nbad, nclip(3))
            End Select
            If (nbad /= 0) Stop 'WriteSubsample: chunk count mismatch between the two passes'
            Call WriteInt2(iu, q, nsel)
        End Do
        close (iu)
        deallocate (q, cnt, off)
        Call system_clock(c1)
        write (*, '(a,i5,a,i14,a,f12.1,a,i10,a,f8.2,a)') '   subsample: M =', Msub, '  n_sel =', nsel, &
            '  (N_p/M =', dble(Nparticles)/dble(Msub), ')  clipped v components =', sum(nclip), &
            '   time =', dble(c1 - c0)/dble(cr), ' s'

        !---- meta
        write (Name, '(2(a,i4.4),a)') 'Subsample.', moment, '.', Nrealization, '.meta'
        open (newunit=iu, file=TRIM(Name), status='replace')
        write (iu, '(a)') 'product = Subsample'
        write (iu, '(a)') 'format = raw little-endian int16, no record markers; arrays qx, qy, qz, vx, vy, vz, each n_sel values, contiguous, in that order (12 bytes per particle)'
        write (iu, '(a,i6)') 'step = ', moment
        write (iu, '(a,i6)') 'realization = ', Nrealization
        write (iu, '(a,es16.8)') 'a = ', AEXPN
        write (iu, '(a,es16.8)') 'z = ', 1./AEXPN - 1.
        write (iu, '(a,es16.8)') 'a_x = ', AEXPN
        write (iu, '(a,es16.8,a)') 'a_v = ', av, '   (velocity epoch: AEXPN - ASTEP/2, ASTEP = the step just taken)'
        write (iu, '(a,f12.4)') 'Box = ', Box
        write (iu, '(a)') 'units = Mpc/h (positions), km/s (velocities)'
        write (iu, '(a,i6)') 'N_grid = ', NGRID
        write (iu, '(a,i6)') 'N_row = ', NROW
        write (iu, '(a,i20)') 'N_p = ', Nparticles
        write (iu, '(a,i6)') 'N_side = ', NROW
        write (iu, '(a,i10)') 'M = ', Msub
        write (iu, '(a,z16.16)') 'SALT = 0x', SubSalt
        write (iu, '(a,i20)') 'n_sel = ', nsel
        write (iu, '(a,es16.8)') 'n_expected = ', dble(Nparticles)/dble(Msub)
        write (iu, '(a,es16.8,a,i6,a)') 'dx_quant = ', dble(Box)/dble(SubNq), '   (Box/', SubNq, ' Mpc/h)'
        write (iu, '(a,f8.4,a)') 'dv_quant = ', SubDv, '   (km/s)'
        write (iu, '(a,i12,a)') 'n_clip_v = ', sum(nclip), '   (velocity components clipped to +-32767, summed over x, y, z)'
        write (iu, '(a,i12)') 'n_clip_vx = ', nclip(1)
        write (iu, '(a,i12)') 'n_clip_vy = ', nclip(2)
        write (iu, '(a,i12)') 'n_clip_vz = ', nclip(3)
        write (iu, '(a)') 'v_def = peculiar km/s'
        write (iu, '(a,es16.8,a)') 'v_scale = ', 100.d0*dble(Box)/dble(NGRID)/av, '   (km/s per code momentum unit: 100 Box/N_grid/a_v)'
        write (iu, '(a)') 'x_decode = x [Mpc/h] = (q + 32768 + 0.5) * Box/65536, q the stored int16; q + 32768 = floor(x/Box * 65536) with x wrapped into [0, Box)'
        write (iu, '(a)') 'v_decode = v [km/s] = s * dv_quant, s the stored int16 = nint(v/dv_quant) clipped to [-32767, 32767]'
        write (iu, '(a)') 'id_rule = ID = index-1 in lattice order; k-th stored particle = k-th ID with splitmix64(ID xor SALT) mod M == 0'
        write (iu, '(a)') 'lattice = ID = index-1 of the particle arrays (0-based). ICs via ic2pm: index = Gadget ID; the IDx/IDz tag in header gives the order, IDx: ID = 1 + ix + Np (iy + Np iz) (FML), IDz: ID = 1 + iz + Np (iy + Np ix) (2LPTic), ID?: undetermined; no tag = pre-fix ic2pm, IC file order. Native PMP2start: ID = 1 + i + Np (j + Np k), i along x'
        write (iu, '(a)') 'hash = splitmix64(x), uint64: z = x + 0x9E3779B97F4A7C15; z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9; z = (z ^ (z >> 27)) * 0x94D049BB133111EB; return z ^ (z >> 31)'
        write (iu, '(a)') 'axis_map = qx, qy, qz follow the particle x, y, z of the halo catalogue'
        write (iu, '(a,f10.6)') 'Omega_m = ', Om
        write (iu, '(a,f10.6)') 'Omega_L = ', OmL
        write (iu, '(a,f10.6)') 'hubble = ', hubble
        write (iu, '(a,i3)') 'MG_flag = ', MG_flag
        write (iu, '(a,i3)') 'MG_model = ', MG_model
        If (MG_flag == 1) Then
            write (iu, '(a,es14.6)') 'fR0 = ', fR0
        Else
            write (iu, '(a)') 'fR0 = none'
        End If
        write (iu, '(a,a)') 'header = ', TRIM(HEADER)
        close (iu)

    Contains
        !---- pass 2 for one component: arr(ID+1) -> q(slot), slots in ID order.
        !     mode 1: position (grid units, wrapped) -> q - 32768
        !     mode 2: velocity (code momentum * scl) -> nint, clipped; nclp = clipped count
        Subroutine FillQuantised(arr, mode, scl, q, nbad, nclp)
            Real*4, intent(in)     :: arr(*)
            Integer*4, intent(in)  :: mode
            Real*8, intent(in)     :: scl
            Integer*2, intent(out) :: q(*)
            Integer*8, intent(out) :: nbad, nclp
            Integer*8 :: ic, id0, n, k, j, i, iq, o
            Real*8 :: s, r, xN
            Integer(c_int32_t), allocatable :: idx(:)

            xN = dble(NGRID)
            nbad = 0
            nclp = 0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ic,id0,n,k,j,i,iq,o,s,r,idx) REDUCTION(+:nbad,nclp)
            allocate (idx(SubChunk))
!$OMP DO SCHEDULE(DYNAMIC,4)
            Do ic = 1, Nchunk
                id0 = (ic - 1)*SubChunk
                n = min(SubChunk, Nparticles - id0)
                k = glam_select_indices(id0, n, SubSalt, int(Msub, c_int64_t), idx)
                If (k /= cnt(ic)) Then
                    nbad = nbad + 1
                    Cycle
                End If
                o = off(ic - 1)
                If (mode == 1) Then
                    Do j = 1, k
                        i = id0 + idx(j) + 1
                        s = dble(arr(i)) - 1.d0                      ! grid units from 0, nominally in [0, NGRID)
                        If (s < 0.d0 .or. s >= xN) s = modulo(s, xN) ! wrap (periodic box)
                        iq = int(s*scl, 8)                           ! floor for s >= 0
                        iq = min(max(iq, 0_8), int(SubNq - 1, 8))    ! guard the rounding edge at s*scl -> 65536
                        q(o + j) = int(iq - 32768_8, 2)
                    End Do
                Else
                    Do j = 1, k
                        i = id0 + idx(j) + 1
                        r = anint(dble(arr(i))*scl)                  ! nint in real*8: no integer overflow
                        If (r > 32767.d0) Then
                            r = 32767.d0
                            nclp = nclp + 1
                        Else If (r < -32767.d0) Then
                            r = -32767.d0
                            nclp = nclp + 1
                        End If
                        q(o + j) = int(r, 2)
                    End Do
                End If
            End Do
!$OMP END DO
            deallocate (idx)
!$OMP END PARALLEL
        End Subroutine FillQuantised

        !---- stream write of an int16 array in blocks (256 MB)
        Subroutine WriteInt2(iu, q, n)
            Integer*4, intent(in) :: iu
            Integer*2, intent(in) :: q(*)
            Integer*8, intent(in) :: n
            Integer*8, parameter :: blk = 134217728_8
            Integer*8 :: i0, nb
            Do i0 = 1, n, blk
                nb = min(blk, n - i0 + 1)
                write (iu) q(i0:i0 + nb - 1)
            End Do
        End Subroutine WriteInt2
    End Subroutine WriteSubsample
End Module Subsample
