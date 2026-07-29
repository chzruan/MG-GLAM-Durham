!--------------------------------------------------
!   Original Authors: Anatoly Klypin and Francisco Prada
!
!   compilation:
!   ifort  -O2  -g -traceback -ftz -unroll  -qopenmp  -shared-intel -mcmodel=medium -convert big_endian -w -c WriteGadgetFormat.f90
!   ifort -O2 -g -traceback -ftz -unroll -qopenmp -shared-intel -mcmodel=medium -convert big_endian -w -c ReadGLAMFiles.f90
!   ifort  -O2   -g -traceback -ftz -unroll  -qopenmp  -shared-intel -mcmodel=medium -convert big_endian -o read.exe ReadGLAMFiles.o  WriteGadgetFormat.o
!
!   memory consumption: 
!   1 particle - [REAL*4 x,y,z,vx,vy,vz] - [4 bytes x 6 = 24 bytes]
!   particle numbers Np_3d = NROW**3
!   total memory = NROW**3 * 24 bytes = (NROW/1024)**3 * 24 GB
!   e.g. 2048**3 particles ~ 192 GB
!
!   input: 
!   - GLAMFolderPath: path to the (MG-)GLAM working folder, e.g. /cosma7/ICC-data/dc-ruan1/proj_bias_rec/GLAM_data/L512Np1024Ng2048/Run01/
!
!   output:
!   - OutputGadgetFilePath
!
!--------------------------------------------------
!
!   Read PM data GLAM-style:
!
!--------------------------------------------------
module Structures
    implicit none 

    Character  :: Name*200, HEADER*45, GLAMFolderPath*200, OutputGadgetFilePath*200
    Real*4     :: AEXPN,ASTEP0,AEXP0,AMPLT,ASTEP,EKIN,EKIN1,EKIN2, &
                TINTG,AU0,AEU0,Om,OmL,hubble, &
                extras(100),ENKIN,ENPOT, &
                Xscale,Vscale,Mscale,Box
    Real*4     :: dum1, dum2, dum3, dum4, dum5, dum6, dum7
    Integer*4  :: ISTEP,Nrealization,Nseed,NGRID,NROW
    Integer*8  :: Nparticles
    Real*4, Allocatable, Dimension(:) :: Xb,Yb,Zb,VXb,VYb,Vzb    !-- reading buffers
    Real*4, Allocatable, Dimension(:) :: Xpar,Ypar,Zpar,VX,VY,VZ !-- particles
    Integer*4 :: FileID_PMcrs, FileID_PMcrd

end module Structures
!--------------------------------------------------
!
!
!
!--------------------------------------------------
Program Read
    use Structures
    use WriteGadgetFormat ! contains pid and the (3 x Nparticles) arrays pos and vel

    Integer*4     :: SnapNum, ii, Nargs
    Character*200 :: SnapNumArg

    Nargs = Command_Argument_Count()
    IF(Nargs < 2 .or. Nargs > 3) THEN
        WRITE(*,*) 'Usage: glam2gadget.exe <GLAMFolderPath/> <SnapNum> [OutputGadgetFilePath]'
        WRITE(*,*) '   <SnapNum> < 0  : the UNNUMBERED files PMcrd.DAT / PMcrs0.DAT'
        WRITE(*,*) '                    (written by PMP2start, and rewritten in place'
        WRITE(*,*) '                     by every checkpoint -- check the reported step)'
        WRITE(*,*) '   <SnapNum> >= 0 : the numbered files PMcrd.NNNN.DAT / PMcrs0.NNNN.DAT'
        WRITE(*,*) '   GLAMFolderPath MUST end in a slash: it is concatenated directly.'
        STOP 'Wrong number of command-line arguments'
    ENDIF

    CALL Get_Command_Argument(1, GLAMFolderPath)
    CALL Get_Command_Argument(2, SnapNumArg)
    READ(SnapNumArg, *) SnapNum

    !--- guard the concatenation: 'dir' + 'PMcrd.DAT' = 'dirPMcrd.DAT'
    ii = LEN_TRIM(GLAMFolderPath)
    IF(ii == 0) STOP 'GLAMFolderPath is empty'
    IF(GLAMFolderPath(ii:ii) /= '/') THEN
        GLAMFolderPath = TRIM(GLAMFolderPath)//'/'
    ENDIF

    IF(Nargs == 3) THEN
        CALL Get_Command_Argument(3, OutputGadgetFilePath)
    ELSE IF(SnapNum < 0) THEN
        Write(OutputGadgetFilePath, '(a,a)') TRIM(GLAMFolderPath), 'particles.ic.gadget'
    ELSE
        Write(OutputGadgetFilePath, '(a,a,i4.4,a)') TRIM(GLAMFolderPath), 'particles.', SnapNum, '.gadget'
    ENDIF

    CALL ReadHeader(SnapNum)

    !--- The unnumbered files are ALSO the checkpoint target: PMP2main.f90 calls
    !    WriteDataPM(0,Path) every Ncheckpoint steps, overwriting them in place.
    !    So "the ICs" are only the ICs if no checkpoint fired. Say so loudly
    !    rather than silently converting a mid-run state labelled as ICs.
    IF(SnapNum < 0 .and. ISTEP /= 0) THEN
        WRITE(*,*) ''
        WRITE(*,*) ' *** WARNING: asked for the ICs but this header says step =', ISTEP
        WRITE(*,'(a,f8.4,a,f9.4)') '     a =', AEXPN, '   z =', 1.0/AEXPN - 1.0
        WRITE(*,*) '     A checkpoint has overwritten PMcrd.DAT/PMcrs0.DAT.'
        WRITE(*,*) '     Converting it anyway -- but it is NOT the initial conditions.'
        WRITE(*,*) ''
    ENDIF

    Allocate(&
        pos(3, Nparticles), &
        vel(3, Nparticles), &
        pid(Nparticles))

    CALL ReadParticles(SnapNum)
    ! now the particle data are stored in pos, vel and pid

!---------------------------------------------------
!   save pos, vel and pid into a Gadget format file
!---------------------------------------------------
    CALL write_gadget(&
        fileout=Trim(OutputGadgetFilePath), &
        npart=Nparticles, &
        nparttot=Nparticles, &
        verbose=.True., &
        megaverbose=.True., &
        Lbox=Box, &
        numfiles=1, &
        hubble=hubble, &
        omega0=Om, &
        omegaL=OmL, &
        aexp=AEXPN)

    write(*,*) 'x,y,z min=',minval(pos)
    write(*,*) 'x,y,z max=',maxval(pos)

end Program Read
!--------------------------------------------------
!
!
!
!--------------------------------------------------
SUBROUTINE ReadHeader(SnapNum)
    use Structures
    Integer*4 :: SnapNum
    Logical   :: exst_hdr

    IF(SnapNum < 0) THEN            !--- unnumbered files (ICs / checkpoint)
        WRITE(Name,'(a,a)') TRIM(GLAMFolderPath), 'PMcrd.DAT'
    ELSE
        WRITE(Name,'(a,a,i4.4,a)') TRIM(GLAMFolderPath), 'PMcrd.',SnapNum,'.DAT'
    ENDIF
    WRITE(*, *) TRIM(Name)
    FileID_PMcrd = 4
    !--- STATUS='OLD', not 'UNKNOWN'. An UNKNOWN open CREATES a zero-byte file
    !    when the snapshot is absent and then dies on the READ with a bare
    !    end-of-file, leaving litter behind that looks like a corrupt snapshot.
    INQUIRE(file=TRIM(Name), EXIST=exst_hdr)
    IF(.not.exst_hdr) THEN
        WRITE(*,*) ' File ', TRIM(Name), ' does not exist'
        STOP 'Header file PMcrd... does not exist. Error'
    ENDIF
    OPEN (FileID_PMcrd, file=TRIM(Name), form ='UNFORMATTED', status ='OLD')

    READ(FileID_PMcrd) HEADER, &
              AEXPN,AEXP0,AMPLT,ASTEP,ISTEP, &
              dum1, dum2, dum3, dum4, dum5, dum6, dum7, &
              NROW,NGRID,Nrealization,Nseed, &
              Om,OmL,hubble, &
              Nparticles, &
              (extras(ii), ii=1,100)
    Box     = extras(100)
    Xscale  = Box/NGRID                    ! scale for comoving coordinates
    Vscale  = 100.*Xscale/AEXPN            ! Scale for velocities
    Mscale  = Om*2.774e+11*(Box/NROW)**3  ! mass of a particle
    
    WRITE (*,'(a,/10x,a,f8.4,3(a,i7))') HEADER, & 
        'a            =', AEXPN, ', step= ',ISTEP, &
        ', Nrow= ', NROW, ', Ngrid=',NGRID

    WRITE (*,'(10x,2(a,f8.3))')   'Omega_matter =', Om, ', hubble= ', hubble
    WRITE (*,'(10x,a,f8.3,a,i5)') 'Box          =', Box, ', Ngrid =', NGRID
  
    CLOSE(FileID_PMcrd)

end SUBROUTINE ReadHeader
!--------------------------------------------------
!
!
!
!--------------------------------------------------
SUBROUTINE ReadParticles(SnapNum)
    use Structures
    use WriteGadgetFormat

    Integer*4, parameter :: Nrecord = 1024**2   !---  setup for reading data files
    Integer*4, parameter :: Naccess = Nrecord*6     !number of real*4 numbers per record  ! *6: (x,y,z,vx,vy,vz)
    Logical     :: exst
    Integer*4   :: SnapNum, ierr
    Integer*8   :: iCount, i_page, i_offset, i_par, Npages, Nlast, NinPage
    Real*4      :: xmin, ymin, zmin, xmax, ymax, zmax, &
                   xmin0, ymin0, zmin0, xmax0, ymax0, zmax0

    xR      = NGRID + 1     ! right boundary for coordinates
    Npages  = (Nparticles - 1)/Nrecord + 1  ! number of records
    WRITE(*, '(a,i12,a,i6)') 'Np_3D = ', Nparticles, ' Npages = ', Npages

    Nlast   = Nparticles - (Npages - 1) * Nrecord   ! number of particles in last record

    Allocate (Xb(Nrecord), Yb(Nrecord), Zb(Nrecord))    !--- allocate reading buffers
    Allocate (VXb(Nrecord), VYb(Nrecord), VZb(Nrecord))

    xmin = 1.e5  ; ymin = 1.e5  ; zmin = 1.e5 
    xmax = -1.e5 ; ymax = -1.e5 ; zmax = -1.e5 
    ifile = 0
    jj    = 0
    ierr  = 0

    FileID_PMcrs = 20
    Call OpenFile(ifile, SnapNum)
 
    !-------- main loop reading particles from PMcrs files -------
    ! Npages, every page has Nrecord=1024**2 particles except the last page
    Do i_page = 1, Npages
        If(i_page==Npages) THEN
            NinPage = Nparticles - (i_page - 1) * Nrecord  ! if last page, read what is left
        Else
            NinPage = Nrecord
        EndIf
        jj = jj + 1     !--- current record to read from file FileID_PMcrs
        IF(i_page<10 .or. i_page==Npages .or. mod(i_page,100)==0) THEN
            WRITE(*, '(3(a,i9))') ' Reading page= ', i_page, ' record =', jj, ' NinPage= ', NinPage
        ENDIF

10  Read(FileID_PMcrs, REC=jj, iostat=ierr) Xb,Yb,Zb,VXb,VYb,VZb    !--- read one page of data

        IF(ierr /= 0) THEN  !--- open next file when previous is finished
            ifile = ifile + 1
            Call OpenFile(ifile, SnapNum)
            jj = 1
            go to 10    ! read the record again
        END IF  !--- end open file

        i_offset = (i_page - 1) * Nrecord
        Do i_par = 1, NinPage
            !--- check for errors
            IF(i_par + i_offset > Nparticles) STOP 'Attempt to read too many particles '
            IF(INT(Xb(i_par)) == Ngrid + 1) Xb(i_par) = Xb(i_par) - 1.e-3
            IF(INT(Yb(i_par)) == Ngrid + 1) Yb(i_par) = Yb(i_par) - 1.e-3
            IF(INT(Zb(i_par)) == Ngrid + 1) Zb(i_par) = Zb(i_par) - 1.e-3
            IF(INT(Xb(i_par)) == Ngrid + 1) WRITE(*, *) 'Error in boundary: ', INT(Xb(i_par)), Xb(i_par)
            IF(INT(Yb(i_par)) == Ngrid + 1) WRITE(*, *) 'Error in boundary: ', INT(Yb(i_par)), Yb(i_par)
            IF(INT(Zb(i_par)) == Ngrid + 1) WRITE(*, *) 'Error in boundary: ', INT(Zb(i_par)), Zb(i_par)
            xmin = MIN(xmin, Xb(i_par))
            ymin = MIN(ymin, Yb(i_par))
            zmin = MIN(zmin, Zb(i_par))
            xmax = MAX(xmax, Xb(i_par))
            ymax = MAX(ymax, Yb(i_par))
            zmax = MAX(zmax, Zb(i_par))
            
            ! rescale coordinates
            ! see Eqn. (11) and (12) of the GLAM documentation
            x = (Xb(i_par) - 1.) * Xscale  !  (x,y,z) = coordinates in comoving Mpch units
            y = (Yb(i_par) - 1.) * Xscale
            z = (Zb(i_par) - 1.) * Xscale
            Wx = VXb(i_par) * Vscale   ! (Wx,Wy,Wz) = peculiar velocities in km/s
            Wy = VYb(i_par) * Vscale
            Wz = VZb(i_par) * Vscale
            
            pos(1, i_par+i_offset) = x
            pos(2, i_par+i_offset) = y
            pos(3, i_par+i_offset) = z
            vel(1, i_par+i_offset) = Wx
            vel(2, i_par+i_offset) = Wy
            vel(3, i_par+i_offset) = Wz
            pid(i_par+i_offset) = i_par + i_offset
        END DO ! Do i_par = 1, NinPage
    END DO ! Do i_page = 1, Npages

    CLOSE(FileID_PMcrs)
    
    WRITE(*,'(3(5x,a,2f10.3))') &
        'x   min/max = ', xmin, xmax, &
        'y   min/max = ', ymin, ymax, &
        'z   min/max = ', zmin, zmax

    DEALLOCATE (Xb,Yb,Zb,VXb,VYb,VZb)

end SUBROUTINE ReadParticles
!--------------------------------------------------
!
!
!
!--------------------------------------------------
SUBROUTINE OpenFile(ifile,SnapNum)
    USE Structures

    Integer*4, parameter :: Nrecord = 1024**2   !---  setup for reading data files
    Integer*4, parameter :: Naccess = Nrecord*6   !number of real*4 numbers per record  ! *6: (x,y,z,vx,vy,vz)
    logical       :: exst
    integer*4     :: ifile, SnapNum
    Character*200 :: PMcrsName

    IF(ifile/=0) THEN 
        CLOSE(FileID_PMcrs)     !-- close previous file
    ENDIF

    IF(SnapNum < 0) THEN            !--- unnumbered files (ICs / checkpoint)
        IF(ifile<10) THEN
            WRITE(PMcrsName, '(a,a,i1.1,a)') TRIM(GLAMFolderPath), 'PMcrs', ifile, '.DAT'
        Else
            WRITE(PMcrsName, '(a,a,i2.2,a)') TRIM(GLAMFolderPath), 'PMcrs', ifile, '.DAT'
        EndIf
    ELSE IF(ifile<10) THEN
        WRITE(PMcrsName, '(a,a,i1.1,a,i4.4,a)') TRIM(GLAMFolderPath), 'PMcrs', ifile, '.', SnapNum, '.DAT'
    Else
        WRITE(PMcrsName, '(a,a,i2.2,a,i4.4,a)') TRIM(GLAMFolderPath), 'PMcrs', ifile, '.', SnapNum, '.DAT'
    EndIf

    INQUIRE(file=TRIM(PMcrsName), EXIST=exst)   ! open file PMcrs
    IF(.not.exst) THEN
        WRITE(*,*) ' File',TRIM(PMcrsName),' does not exist'
        Stop ' File PMcrs... does not exist. Error'
    End IF

    FileID_PMcrs = 20
    OPEN(&
        UNIT=FileID_PMcrs, &
        FILE=TRIM(PMcrsName), &
        ACCESS='DIRECT', &
        FORM='unformatted', &
        STATUS='UNKNOWN', &
        RECL=NACCESS)

end SUBROUTINE OpenFile

