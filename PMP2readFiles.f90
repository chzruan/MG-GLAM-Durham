!--------------------------------------------------
!
!        Read PM data GLAM-style and make simple statistics
!
!--------------------------------------------------
module Structures
  Character :: Name*80, HEADER*45
  Real*4     :: AEXPN,ASTEP0,AEXP0,AMPLT,ASTEP,EKIN,EKIN1,EKIN2,  &
                TINTG,AU0,AEU0,Om,OmL,hubble,                     &
                extras(100),ENKIN,ENPOT,                          &
                Xscale,Vscale,Mscale,Box
  Integer*4  :: ISTEP,Nrealization,Nseed,NGRID,NROW
  Integer*8  :: Nparticles
  Real*4, Allocatable, Dimension(:) :: Xb,Yb,Zb,VXb,VYb,Vzb    !-- reading buffers
  Real*4, Allocatable, Dimension(:) :: Xpar,Ypar,Zpar,VX,VY,VZ !-- particles
end module Structures
!--------------------------------------------------
Program Read
  use Structures

  write(*,'(a,$)') ' Enter snapshot number: '
  read(*,*) moment
  write(*,'(a,$)') ' Do you want (1) or not (0) to read particles into RAM?  '
  read(*,*) iStore
  

  CALL ReadHeader(moment)
  If(iStore==1)Allocate(Xpar(Nparticles),Ypar(Nparticles),Zpar(Nparticles), &
                        VX(Nparticles),VY(Nparticles),VZ(Nparticles))
  CALL ReadParticles(moment,iStore)
  If(iStore==1)CALL RescaleParticles
     
     
end Program Read
!--------------------------------------------------
SUBROUTINE ReadHeader(moment)
  use Structures
  Integer*4 :: moment

  write(Name,'(a,i4.4,a)')'PMcrd.',moment,'.DAT'         !--- construct filename
  Open (4,file =TRIM(Name),form ='UNFORMATTED',status ='UNKNOWN')
  
      READ  (4) HEADER,                        &         !--- read header file
                       AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW, &
                       TINTG,EKIN,EKIN1,EKIN2,AU0,AEU0,     &
                       NROW,NGRID,Nrealization,Nseed,Om,OmL,hubble, &
                       Nparticles,extras
      Box     = extras(100)
      Xscale  = Box/NGRID                    ! scale for comoving coordinates
      Vscale  = 100.*Xscale/AEXPN            ! Scale for velocities
      Mscale  = Om*2.774e+11*(Box/NROW)**3  ! mass of a particle
      WRITE (*,'(a,/10x,a,f8.4,4(a,i7))') HEADER,            & 
                       ' a=',AEXPN, ' step= ',ISTEP,         &
                       ' Nrow= ', NROW, ' Ngrid=',NGRID
      WRITE (*,'(10x,2(a,f8.3))')   'Omega_matter =',Om,' hubble= ',hubble
      WRITE (*,'(10x,a,f8.3,a,i5)') 'Box          =',Box,' Ngrid =',NGRID
      close(4)

    end SUBROUTINE ReadHeader
!--------------------------------------------------
SUBROUTINE ReadParticles(moment,iStore)
  use Structures
  Integer*4, parameter :: Nrecord = 1024**2, &     !---  setup for reading data files
                          Naccess = Nrecord*6!*4      !     number of real*4 numbers per record  Logical      :: exst
  Integer*4    :: moment,iStore
  Integer*8    :: iCount,ii,ioff,ip,Npages,Nlast,Jpage,NinPage
  Real*4       :: xmin,ymin,zmin,xmax,ymax,zmax, &
                  xmin0,ymin0,zmin0,xmax0,ymax0,zmax0
      xR      = NGRID +1                        ! right bondary for coordinates
      Npages   = (Nparticles-1)/Nrecord+1       ! number of records
      Nlast   = Nparticles - (Npages-1)*Nrecord ! number of particles in last record
      Jpage   = Nrecord
 
      Allocate (Xb(Nrecord),Yb(Nrecord),Zb(Nrecord))     !--- allocate reading buffers
      Allocate (VXb(Nrecord),VYb(Nrecord),VZb(Nrecord))

      xmin = 1.e5  ; ymin = 1.e5  ; zmin = 1.e5 
      xmax = -1.e5 ; ymax = -1.e5 ; zmax = -1.e5 
  ifile = 0
  jj    = 0
  Call OpenFile(ifile,moment)
  Open(50,file='particles.dat')
  write(50,'(a,i4,2(a,es12.4))') ' Particles for snapshot= ',ISTEP,' Aexpn= ',AEXPN,' Box= ',Box
  
                        !-------- main loop reading particles from Npages -------      
   Do ii =1,Npages
      If(ii==Npages)Then
         NinPage = Nparticles -(ii-1)*JPAGE  ! if last page, read what is left
      Else
         NinPage = JPAGE
      EndIf
      jj = jj +1                                     !--- current record to read from file 20
      If(ii<10.or.ii==Npages.or.mod(ii,100)==0)write(*,'(3(a,i9))') ' Reading page= ',ii,' record =',jj,' NinPage= ',NinPage
10    Read(20,REC=jj,iostat=ierr) Xb,Yb,Zb,VXb,VYb,VZb        !--- read one page of data
      
      If(ierr /=0)Then                               !--- open next file when previous is finished
         ifile = ifile +1
         Call OpenFile(ifile,moment)
         jj = 1
         go to 10                                    !     read the record again
       end If                                        !--- end open file

       ioff = (ii-1)*JPAGE
       Do ip =1,NinPage
                                                     !--- check for errors
              If(ip+ioff > Nparticles)STOP 'Attempt to read too many particles '
              if(INT(Xb(ip))==Ngrid+1)Xb(ip)=Xb(ip)-1.e-3
              if(INT(Yb(ip))==Ngrid+1)Yb(ip)=Yb(ip)-1.e-3
              if(INT(Zb(ip))==Ngrid+1)Zb(ip)=Zb(ip)-1.e-3
              if(INT(Xb(ip))==Ngrid+1)write(*,*)'Error in boundary: ',INT(Xb(ip)),Xb(ip)
              if(INT(Yb(ip))==Ngrid+1)write(*,*)'Error in boundary: ',INT(Yb(ip)),Yb(ip)
              if(INT(Zb(ip))==Ngrid+1)write(*,*)'Error in boundary: ',INT(Zb(ip)),Zb(ip)
              xmin = min(xmin,Xb(ip))
              ymin = min(ymin,Yb(ip))
              zmin = min(ymin,Zb(ip))
              xmax = max(xmax,Xb(ip))
              ymax = max(ymax,Yb(ip))
              zmax = max(zmax,Zb(ip))
                                                    !--- rescale coordinates
              x = (Xb(ip)-1.)*Xscale
              y = (Yb(ip)-1.)*Xscale
              z = (Zb(ip)-1.)*Xscale
                                                    !--- write some coordiantes to a file
              If(x<50..and.y<50..and.z<50.)write(50,'(i12,3x,3f10.4)')ip,x,y,z
              If(iStore == 1)Then        !--- store the data
                Xpar(ip+ioff) = Xb(ip) 
                Ypar(ip+ioff) = Yb(ip) 
                Zpar(ip+ioff) = Zb(ip) 
                  VX(ip+ioff) = VXb(ip) 
                  VY(ip+ioff) = VYb(ip) 
                  VZ(ip+ioff) = VZb(ip)
              EndIf
           end Do
        end DO
        close (20)
        close(50)
           write(*,'(3(5x,a,2f10.3))')'x   min/max= ',xmin,xmax,'y   min/max= ',ymin,ymax,'z   min/max= ',zmin,zmax

      DEALLOCATE (Xb,Yb,Zb,VXb,VYb,VZb)
           
   end SUBROUTINE ReadParticles
!-----------------------------------------------
SUBROUTINE OpenFile(ifile,moment)
  Integer*4, parameter :: Nrecord = 1024**2, &     !---  setup for reading data files
                          Naccess = Nrecord*6!*4      !     number of real*4 numbers per record  Logical      :: exst
  logical :: exst
  integer*4    :: ifile,moment
  Character*80 :: Name

  if(ifile/=0)close(20)     !-- close previous file
     If(ifile<10)Then
        write(Name,'(a,i1.1,a,i4.4,a)')'PMcrs',ifile,'.',moment,'.DAT'
     Else
        write(Name,'(a,i2.2,a,i4.4,a)')'PMcrs',ifile,'.',moment,'.DAT'
     EndIf

   INQUIRE(file=TRIM(Name),EXIST=exst)                         !    Open file PMcrs
   If(.not.exst)Then
      write(*,*) ' File',TRIM(Name),' does not exist'
      Stop ' File PMcrs... does not exist. Error'
   End If
   OPEN(UNIT=20,FILE=TRIM(Name),ACCESS='DIRECT', &
        FORM='unformatted',STATUS='UNKNOWN',RECL=NACCESS)
   write(*,*) ' Open file ',ifile
   end SUBROUTINE OpenFile
!--------------------------------------------------
SUBROUTINE RescaleParticles
  use Structures
  Integer*8 :: ip
!$OMP PARALLEL DO DEFAULT(SHARED)  PRIVATE (ip)
    Do ip =1,Nparticles            
         Xpar(ip) = (Xpar(ip)-1.)*Xscale 
         Ypar(ip) = (Ypar(ip)-1.)*Xscale 
         Zpar(ip) = (Zpar(ip)-1.)*Xscale 
        VX(ip) =  VX(ip)*Vscale 
        VY(ip) =  VY(ip)*Vscale 
        VZ(ip) =  VZ(ip)*Vscale
     end Do
   end SUBROUTINE RescaleParticles
