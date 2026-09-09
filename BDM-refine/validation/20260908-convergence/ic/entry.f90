! Campaign-only standalone entry: native first-order ICs, no PM density allocation.
program matched_native_start
    use setInitialConditions
    implicit none
    real*4 :: fact, xcons, vcons, aexpv, partw
    call Timing(0, -1)
    iFlip = 1
    call CheckInit
    call ReadSetup
    call ReadMatchedControl
    call Initialize
    call ReadPkTable
    call SetRandomN
    AEXP0 = AEXPN
    aexpv = AEXPN - ASTEP/2.
    fact = sqrt(Om + OmL*aexpv**3)
    call SPECTR
    vcons = -iFlip*ALPHA/(2.*PI/NGRID)*(aexpv/AEXP0)*sqrt(aexpv)*fact
    xcons = iFlip*ALPHA/(2.*PI/NGRID)*(AEXPN/AEXP0)
    call WriteMatchedReceipt(xcons, vcons)
    if (ic_normalize_only) stop
    call WriteMatchedModes
    call FilesOpen
    call VECTOR
    call BLOCKS(xcons, vcons)
    deallocate(GRX, GRY, GRZ)
    partw = (real(NGRID)/real(NROW))**3
    EKIN = 0.5*SKINE/aexpv**2*partw
    call FilesWrite
    write(25) ASTEP
    close(25)
    close(9)
    close(16)
    call Timing(0, 1)
    write(*,'(a,f12.3)') 'Matched IC total wall seconds = ', CPU(0)
end program matched_native_start
