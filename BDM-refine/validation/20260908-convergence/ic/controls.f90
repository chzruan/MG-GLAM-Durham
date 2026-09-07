! This fragment is inserted into the native setInitialConditions module.
    subroutine ReadMatchedControl
        integer :: unit, ios
        namelist /matched_ic/ ic_master_nrow, ic_origin_ngrid, ic_alpha, ic_normalize_only
        open(newunit=unit, file='matched_ic.nml', status='old', action='read', iostat=ios)
        if (ios /= 0) error stop 'Missing required matched_ic.nml'
        read(unit, nml=matched_ic, iostat=ios)
        close(unit)
        if (ios /= 0) error stop 'Invalid matched_ic.nml'
        if (NROW < 4 .or. mod(NROW,2) /= 0) error stop 'NROW must be even and >= 4'
        if (ic_master_nrow < NROW .or. mod(ic_master_nrow,2) /= 0) error stop 'Invalid master NROW'
        if (ic_master_nrow > 8192) error stop 'Master exceeds native RNG/FFT limit'
        if (NGRID < NROW .or. ic_origin_ngrid < 1) error stop 'Invalid PM/origin mesh'
        if (.not. ieee_is_finite(ic_alpha) .or. ic_alpha == 0.) error stop 'Invalid ic_alpha'
        if (ic_alpha < 0. .and. NROW /= ic_master_nrow) error stop 'Only master can normalize'
        if (ic_normalize_only .and. NROW /= ic_master_nrow) error stop 'Only master normalize-only'
        if (.not. ieee_is_finite(Box) .or. Box <= 0.) error stop 'Invalid box'
        if (.not. ieee_is_finite(AMPLT) .or. AMPLT <= 0.) error stop 'Invalid Setup amplitude'
        if (ASTEP0 <= 0. .or. ASTEP0 >= 2.*AEXPN0) error stop 'Invalid initial velocity time'
    end subroutine ReadMatchedControl

    subroutine WriteMatchedReceipt(xcons, vcons)
        real*4, intent(in) :: xcons, vcons
        integer :: unit
        open(newunit=unit, file='matched_ic_receipt.txt', status='replace')
        write(unit,'(a)') 'format = native-luxury-master-stride-v1'
        write(unit,'(a,i0)') 'nrow = ', NROW
        write(unit,'(a,i0)') 'ngrid = ', NGRID
        write(unit,'(a,i0)') 'master_nrow = ', ic_master_nrow
        write(unit,'(a,i0)') 'origin_ngrid = ', ic_origin_ngrid
        write(unit,'(a,i0)') 'realization = ', Nrealization
        write(unit,'(a,i0)') 'seed = ', Nseed
        write(unit,'(a,es25.16)') 'alpha = ', real(ALPHA,8)
        write(unit,'(a,es25.16)') 'alpha_requested = ', real(ic_alpha,8)
        write(unit,'(a,es25.16)') 'spectrum_sum = ', ic_spectrum_sum
        write(unit,'(a,es25.16)') 'setup_amplt = ', real(AMPLT,8)
        write(unit,'(a,es25.16)') 'box_mpc_h = ', real(Box,8)
        write(unit,'(a,es25.16)') 'origin_mpc_h = ', real(Box,8)/(2.d0*ic_origin_ngrid)
        write(unit,'(a,es25.16)') 'a_position = ', real(AEXPN,8)
        write(unit,'(a,es25.16)') 'a_velocity = ', real(AEXPN-ASTEP/2.,8)
        write(unit,'(a,es25.16)') 'astep = ', real(ASTEP,8)
        write(unit,'(a,es25.16)') 'xcons = ', real(xcons,8)
        write(unit,'(a,es25.16)') 'vcons = ', real(vcons,8)
        write(unit,'(a,l1)') 'normalize_only = ', ic_normalize_only
        write(unit,'(a)') 'nyquist_planes = zero'
        write(unit,'(a)') 'initial_power_mesh_allocation = disabled'
        close(unit)
    end subroutine WriteMatchedReceipt

    subroutine WriteMatchedModes
        integer :: unit, extent
        extent = min(33, NROW-1)
        open(newunit=unit,file='matched_modes.bin',access='stream',form='unformatted',status='replace')
        write(unit) extent
        write(unit) GRX(1:extent,1:extent,1:extent)
        write(unit) GRY(1:extent,1:extent,1:extent)
        write(unit) GRZ(1:extent,1:extent,1:extent)
        close(unit)
    end subroutine WriteMatchedModes
