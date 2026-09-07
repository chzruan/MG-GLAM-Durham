program matched_probe
    use setInitialConditions
    implicit none
    character*128 :: arg
    integer :: unit
    real*4 :: xcons, vcons, aexpv, fact
    call get_command_argument(1, arg)
    read(arg,*) NROW
    call get_command_argument(2, arg)
    read(arg,*) NGRID
    call get_command_argument(3, arg)
    read(arg,*) ASTEP0
    Box = 256.
    AMPLT = 0.03
    Om = 0.307
    OmL = 0.693
    AEXPN0 = 1./101.
    call ReadMatchedControl
    call Initialize
    call SetRandomN
    Ntab = 2
    xkt(1) = 1.e-6
    xkt(2) = 1.e6
    Pkt(1:2) = 1.
    alog0 = -6.
    StepK = 12.
    AEXP0 = AEXPN
    aexpv = AEXPN-ASTEP/2.
    fact = sqrt(Om+OmL*aexpv**3)
    call SPECTR
    xcons = ALPHA/(2.*PI/NGRID)*(AEXPN/AEXP0)
    vcons = -ALPHA/(2.*PI/NGRID)*(aexpv/AEXP0)*sqrt(aexpv)*fact
    call WriteMatchedReceipt(xcons,vcons)
    if(ic_normalize_only) stop
    call WriteMatchedModes
    open(newunit=unit,file='packed.bin',access='stream',form='unformatted',status='replace')
    write(unit) GRX, GRY, GRZ
    close(unit)
    call VECTOR
    open(newunit=unit,file='real.bin',access='stream',form='unformatted',status='replace')
    write(unit) GRX, GRY, GRZ
    close(unit)
    open(16,file='Results.log',status='replace')
    call BLOCKS(xcons,vcons)
    close(16)
    open(newunit=unit,file='particles.bin',access='stream',form='unformatted',status='replace')
    write(unit) XPAR,YPAR,ZPAR,VX,VY,VZ
    close(unit)
end program matched_probe
