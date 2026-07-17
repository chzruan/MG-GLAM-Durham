!--------------------------------------------------
!
!        Read particle sample, make 3D overdensity grid,
!        Calculate power spectrum
!          
!          
!          
!          
!          
!          
!
!--------------------------------------------------
Module Parameters
  Integer*4 ::  Ng, NGRIDnew, NGRIDold
  character(len=80) :: Fname
end Module Parameters

!--------------------------------------------------
Program Spectrum
  use Tools
  use Density
  use Power
  use Parameters
  integer*8 :: Ntot
  Path = ''
  
  Call Initialize(moment,iFlag)
   NGRIDold  = NGRID                  ! store old value of NGRID
   NGRID     = NGRIDnew
   Ndiv      = NGRIDold/NGRIDnew
   Mem_current = Memory(1_8*NGRID*NGRID*NGRID)
   ALLOCATE(FI(NGRID,NGRID,NGRID))
   AvWeight = Float(NGRIDold)**3/Float(NROW)**3
    
   write(*,*) ' NGRIDold = ',NGRIDold
   write(*,*) ' NGRIDnew = ',NGRID
   
   Do ksub =1,Ndiv
   Do jsub =1,Ndiv
   Do isub =1,Ndiv
      mybox = isub +(jsub-1)*Ndiv + (ksub-1)*Ndiv**2
      write(Fname,'(a,i4.4,a,i3.3,a,i4.4,a,i4.4,a)')'PowerDM.',moment,'.',Ndiv,'.',Nrealization,'.',mybox,'.dat'
      Open(12,file=TRIM(Fname))                                 ! output
      write(Fname,'(a,i4.4,a,i3.3,a,i4.4,a,i4.4,a)')'DensDisrt.',moment,'.',Ndiv,'.',Nrealization,'.',mybox,'.dat'
      Open(18,file=TRIM(Fname))                                 ! output
   
      Call DensSubbox(isub,jsub,ksub,iFlag,Ntot)         !- density in real space with new Ngrid
      write(12,'(a)') HEADER
      write(12,'(a,f9.4,a,i4,a,f6.1,a,2i5,a,i4,a,3i3,a,i3,a,es12.4)')  &
                      'Aexpn= ',AEXPN,' Istep= ',ISTEP,' Box =',Box,&
                      ' Ngrids= ',NGRIDold,NGRIDnew,' Realization= ',Nrealization, &
           ' subbox =',isub,jsub,ksub,' Weight= ',iFlag,' AvDens= ',AvWeight*float(Ntot)/float(NGRID)**3
      write(18,'(a)') HEADER
      write(18,'(a,f9.4,a,i4,a,f6.1,a,2i5,a,i4,a,3i3,a,i3,a,es12.4)')  &
                      'Aexpn= ',AEXPN,' Istep= ',ISTEP,' Box =',Box,&
                      ' Ngrids= ',NGRIDold,NGRIDnew,' Realization= ',Nrealization, &
                      ' subbox =',isub,jsub,ksub,' Weight= ',iFlag,' AvDens= ',AvWeight*float(Ntot)/float(NGRID)**3
      Call DensDistr
      Call GetPowerSSC
      close(12)
      close(18)
   end Do
   end Do
   end Do

   DEALLOCATE(FI)
            
end Program Spectrum

!--------------------------------------------
subroutine Initialize(moment,iFlag)
!--------------------------------------------
  use Tools
  use Parameters
  logical :: FileExists
  character*80 :: Fout
           write(*,'(a,$)') ' Enter  snapshot    number = '
           read(*,*) moment
           write(*,'(a,$)') ' Enter  new Mesh size      = '
           read(*,*) NGRIDnew
           write(*,'(a,$)') ' Enter  weight flag        = '
           read(*,*) iFlag
           iFlag = max(min(iFlag,1),0)

           If(mod(NGRID,NGRIDnew)/=0)Stop ' New Ngrid is not multiple of Old Ngrid'
           CALL ReadDataPMcoords(moment,Path)
           
      end subroutine Initialize
!---------------------------------------
!        Read    PMfiles
!             moment <0    use PMcrd.DAT, PMcrs0.DAT ...    
!             moment >= 0  use PMcrd.xxxx.DAT, PMcrs0,XXXX.DAT ..
    SUBROUTINE ReadDataPMcoords(moment,Path)
!      
!---------------------------------------
use Tools
use Parameters
      Character*80 :: Name
        Logical      :: exst
        Integer*8    :: iCount,ii,ioff,ip
        Integer*8    :: Ngal,Nrecord,Jpage,NinPage
        character(len=*) :: Path 
!
        !			Read data and open files
        If(moment<0)Then
           write(Name,'(2a)')Trim(Path),'PMcrd.DAT'
           write(*,*)Trim(Name)
         Open (4,file =Trim(Name),form ='UNFORMATTED',status ='UNKNOWN')
      Else
         write(Name,'(a,i4.4,a)')'PMcrd.',moment,'.DAT'
         Open (4,file =Path//TRIM(Name),form ='UNFORMATTED',status ='UNKNOWN')
      end If
         
      READ  (4) HEADER,                        &
                       AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW, &
                       TINTG,EKIN,EKIN1,EKIN2,AU0,AEU0,     &
                       NROW,NGRID,Nrealization,Nseed,Om,OmL,hubble, &
                       Nparticles,extras
      WRITE (*,'(a,/10x,a,f8.4,4(a,i7))') HEADER,                 & 
                       ' a=',AEXPN, ' step= ',ISTEP,         &
                       ' Nrow= ', NROW, ' Ngrid=',NGRID
      WRITE (17,'(a,/10x,a,f8.4,4(a,i7))') HEADER,                 & 
                       ' a=',AEXPN, ' step= ',ISTEP,         &
                       ' Nrow= ', NROW, ' Ngrid=',NGRID
      close(4)

      Nrecord = 1024**2
      Naccess = Nrecord*6 !*4
      xR      = NGRID +1
      boxsize = extras(100)
      Box     = boxsize
      Npages   = (Nparticles-1)/Nrecord+1 ! number of records
      Nlast   = Nparticles - (Npages-1)*Nrecord ! number of particles in last record
      Jpage   = Nrecord

      write(17,'(a,i10)') ' NROW   =',NROW
      write(17,'(a,i10)') ' Ngal   =',Nparticles
      write(17,'(a,i10)') ' Ngrid  =',Ngrid
      write(17,'(a,f10.1)') ' Box    =',Box
      Allocate (Xb(Nrecord),Yb(Nrecord),Zb(Nrecord))
      Allocate (VXb(Nrecord),VYb(Nrecord),VZb(Nrecord))

      myMemory =Memory(6_8*Nparticles)
      Allocate (XPAR(Nparticles),YPAR(Nparticles),ZPAR(Nparticles))

      iCount = 0

      ifile = 0
      jj    = 0
      If(moment<0)Then
         write(Name,'(2a,i1.1,a)')Trim(Path),'PMcrs',ifile,'.DAT'
      Else
         write(Name,'(2a,i1.1,a,i4.4,a)')Trim(Path),'PMcrs',ifile,'.',moment,'.DAT'
      End If
        INQUIRE(file=TRIM(Name),EXIST=exst)
        If(.not.exst)Stop ' File PMcrs0.DAT does not exist. Error'
        OPEN(UNIT=20,FILE=TRIM(Name),ACCESS='DIRECT', &
              FORM='unformatted',STATUS='UNKNOWN',RECL=NACCESS)
   Do ii =1,Npages
      If(ii==Npages)Then
         NinPage = Nparticles -(ii-1)*JPAGE  ! # particles in the current page
      Else
         NinPage = JPAGE
      EndIf
      jj = jj +1
      If(ii<10.or.ii==Npages)write(17,'(3(a,i9))') ' Reading page= ',ii,' record =',jj,' NinPage= ',NinPage
10      Read(20,REC=jj,iostat=ierr) Xb,Yb,Zb,VXb,VYb,VZb
      If(ierr /=0)Then
         close(20)
         ifile = ifile +1
         If(Moment<0)Then
           If(ifile<10)Then
             write(Name,'(2a,i1.1,a,i4.4,a)')Trim(Path),'PMcrs',ifile,'.DAT'
           Else
             write(Name,'(2a,i2.2,a,i4.4,a)')Trim(Path),'PMcrs',ifile,'.DAT'
          EndIf
       Else
          If(ifile<10)Then
             write(Name,'(2a,i1.1,a,i4.4,a)')Trim(Path),'PMcrs',ifile,'.',moment,'.DAT'
          Else
             write(Name,'(2a,i2.2,a,i4.4,a)')Trim(Path),'PMcrs',ifile,'.',moment,'.DAT'
          EndIf
       end If

           jj = 1
           INQUIRE(file=TRIM(Name),EXIST=exst)
           If(.not.exst)Then
             write(*,'(a,3i5)')' Attempting to read file number, record = ',ifile,ii,Npages
             write(*,'(2a)')' Attempting to read non-existing file: ',TRIM(Name)
             Stop ' Error reading PMcrs files: Did not get all files'
           End If
           Open(20,file=TRIM(Name),ACCESS='DIRECT', &
              FORM='unformatted',STATUS='UNKNOWN',RECL=NACCESS)
          write(17,'(2i7,2a,3x,i9)') ii,ifile,' Open file = ',TRIM(Name),Ninpage
          go to 10
       end If

       ioff = (ii-1)*JPAGE
!$OMP PARALLEL DO DEFAULT(SHARED)  PRIVATE (ip)       
           Do ip =1,NinPage
              If(ip+ioff > Nparticles)STOP 'Attempt to read too many particles '
              if(INT(Xb(ip))==Ngrid+1)Xb(ip)=Xb(ip)-1.e-3
              if(INT(Yb(ip))==Ngrid+1)Yb(ip)=Yb(ip)-1.e-3
              if(INT(Zb(ip))==Ngrid+1)Zb(ip)=Zb(ip)-1.e-3
              if(INT(Xb(ip))==Ngrid+1)write(*,*)'Error in boundary: ',INT(Xb(ip)),Xb(ip)
              if(INT(Yb(ip))==Ngrid+1)write(*,*)'Error in boundary: ',INT(Yb(ip)),Yb(ip)
              if(INT(Zb(ip))==Ngrid+1)write(*,*)'Error in boundary: ',INT(Zb(ip)),Zb(ip)
              
                Xpar(ip+ioff) = Xb(ip) 
                Ypar(ip+ioff) = Yb(ip) 
                Zpar(ip+ioff) = Zb(ip) 
           end Do
        end DO
        close (20)
      
           xx = MAXVAL(Xpar)
           xm = MINVAL(Xpar)
           write(17,*)' x     min/max= ',xm,xx
           xx = MAXVAL(Ypar)
           xm = MINVAL(Ypar)
           write(17,*)' y     min/max= ',xm,xx
           xx = MAXVAL(Zpar)
           xm = MINVAL(Zpar)
           write(17,*)' z     min/max= ',xm,xx
           Np = Nparticles
           Do ii =1,Nparticles
              If(Xpar(ii).lt.1.0.or.Ypar(ii).lt.1.0.or.Zpar(ii).lt.1.0)&
                write(17,'(a,i10,1p,3g14.5)')' Error coord: ',ii,Xpar(ii),Ypar(ii),Zpar(ii)
           EndDo
     DEALLOCATE (Xb,Yb,Zb,VXb,VYb,VZb)
     
   end SUBROUTINE ReadDataPMcoords
!----------------------------------------------------
function seconds ()
!----------------------------------------------------
!
!     purpose: returns elapsed time in seconds
      Integer*8, SAVE :: first=0,rate=0,i0=0
      Integer*8       :: i

      If(first==0)Then
         CALL SYSTEM_CLOCK(i,rate)
         first =1
         i0    = i
         seconds = 0.
      Else
         CALL SYSTEM_CLOCK(i)
         seconds = float(i-i0)/float(rate)
      EndIf

    end function seconds

!--------------------------------------------
!
!          Make density field of galaxies    
!
subroutine DensSubbox(isub,jsub,ksub,iFlag,Ntot)
!
!--------------------------------------------
use Tools
use Parameters
character*120 :: Name

Real*8 :: ss,X,Y,Z,D1,D2,D3,T1,T2,T3
Integer*8 ::  ip,Ntot

        Ngalaxies  = Nparticles
        XN   =FLOAT(NGRID)+1.-1.E-7
        YN   =FLOAT(NGRID)
        il   = (isub-1)*NGRID+1
        ir   =  isub   *NGRID+1
        jl   = (jsub-1)*NGRID+1
        jr   =  jsub   *NGRID+1
        kl   = (ksub-1)*NGRID+1
        kr   =  ksub   *NGRID+1
        write(*,'(a,3i5)') ' Subbox= ',isub,jsub,ksub
!				       Subtract mean density
!$OMP PARALLEL DO DEFAULT(SHARED) & 
!$OMP PRIVATE (M1,M2,M3)
      DO M3=1,NGRID
       DO M2=1,NGRID
	      DO M1=1,NGRID
	        FI(M1,M2,M3) = 0.
	      END DO
       END DO
      END DO

!$OMP PARALLEL DO DEFAULT(SHARED) & 
!$OMP PRIVATE (ip,X,Y,Z,I,J,K,D1,D2,D3,T1,T2,T3,I1,J1,K1)
      Do ip = 1,Ngalaxies
         !if(mod(ip,50000)==0)write(*,*) ' galaxy=',ip
            X = Xpar(ip)
            Y = Ypar(ip)
            Z = Zpar(ip)
              I = INT(X)
              J = INT(Y)
              K = INT(Z)
	        I1=I+1
	        J1=J+1
                K1=K+1
                If(I1.ge.NGRIDold+1)I1=1
                If(J1.ge.NGRIDold+1)J1=1
                If(K1.ge.NGRIDold+1)K1=1
                !If(I.ge.1.and.I.le.NGRID.and.J.ge.1.and.J.le.NGRID.and.K.ge.1.and.K.le.NGRID)Then
                !   FI(I ,J ,K ) =FI(I,J ,K ) +1.
                !end If
                
                   !---------------------------------------- CIC
	        D1=X-FLOAT(I)
	        D2=Y-FLOAT(J)
	        D3=Z-FLOAT(K)
	        T1=1.-D1
	        T2=1.-D2
	        T3=1.-D3
         
         I = I -il
         J = J -jl
         K = K -kl
         I1 = I1 -il
         J1 = J1 -jl
         K1 = K1 -kl
         !If(X>=il.and.X<=ir.and.Y>=jl.and.Y<=jr.and.Z>=kl.and.Z<=kr) &
         !     write(*,'(/i11,3x,3f9.3,3x,3i6,3i6)') ip,X,Y,Z,I,J,K,I1,J1,K1
         
               If(K.ge.1.and.K.le.NGRID)Then         !-- K
                  if(J.ge.1.and.J.le.NGRID)Then
                     if(I.ge.1.and.I.le.NGRID)Then
                       ! write(*,'(i11,3x,3f9.3,3x,3i5,i3)') ip,X,Y,Z,I,J,K,1
!$OMP Atomic         
                       FI(I ,J ,K ) =FI(I,J ,K ) +T3*T1*T2
                    end if
                    if(I1.ge.1.and.I1.le.NGRID)Then                     
                       ! write(*,'(i11,3x,3f9.3,3x,3i5,i3)') ip,X,Y,Z,I1,J,K,2
!$OMP Atomic         
                       FI(I1,J ,K ) =FI(I1,J ,K ) +T3*D1*T2
                    end if
                 end if
                  if(J1.ge.1.and.J1.le.NGRID)Then                    
                    if(I.ge.1.and.I.le.NGRID)Then                                         
                       ! write(*,'(i11,3x,3f9.3,3x,3i5,i3)') ip,X,Y,Z,I,J1,K,3
!$OMP Atomic         
                       FI(I ,J1,K ) =FI(I ,J1,K ) +T3*T1*D2
                    end if
                    if(I1.ge.1.and.I1.le.NGRID)Then                                         
                       ! write(*,'(i11,3x,3f9.3,3x,3i5,i3)') ip,X,Y,Z,I1,J1,K,4
!$OMP Atomic         
                     FI(I1,J1,K ) =FI(I1,J1,K ) +T3*D1*D2
                    End If
                  end if     !J1
               end If         ! K
               If(K1.ge.1.and.K1.le.NGRID)Then         !-- K1
                  if(J.ge.1.and.J.le.NGRID)Then
                    if(I.ge.1.and.I.le.NGRID)Then                     
                       ! write(*,'(i11,3x,3f9.3,3x,3i5,i3)') ip,X,Y,Z,I,J,K1,5
!$OMP Atomic         
                       FI(I ,J ,K1 ) =FI(I ,J ,K1 ) +D3*T1*T2
                    end if
                    if(I1.ge.1.and.I1.le.NGRID)Then                     
                      !  write(*,'(i11,3x,3f9.3,3x,3i5,i3)') ip,X,Y,Z,I1,J,K1,6
!$OMP Atomic         
                       FI(I1,J ,K1 ) =FI(I1,J ,K1 ) +D3*D1*T2
                    end if
                 end if
                  if(J1.ge.1.and.J1.le.NGRID)Then                    
                    if(I.ge.1.and.I.le.NGRID)Then                                         
                       ! write(*,'(i11,3x,3f9.3,3x,3i5,i3)') ip,X,Y,Z,I,J1,K1,7
!$OMP Atomic         
                       FI(I ,J1,K1 ) =FI(I ,J1,K1 ) +D3*T1*D2
                    end if
                    if(I1.ge.1.and.I1.le.NGRID)Then                                         
                      !  write(*,'(i11,3x,3f9.3,3x,3i5,i3)') ip,X,Y,Z,I1,J1,K1,8
!$OMP Atomic         
                     FI(I1,J1,K1 ) =FI(I1,J1,K1 ) +D3*D1*D2
                    End If
                  end if   ! J1
               end If       ! K1
            END Do

        SS =0.
        D2 =0.
!$OMP PARALLEL DO DEFAULT(SHARED) &         !--- find mass inside subbox
!$OMP PRIVATE (M1,M2,M3) Reduction(+:SS,D2)
      DO M3=1,NGRID
       DO M2=1,NGRID
	      DO M1=1,NGRID
	        SS = SS+ FI(M1,M2,M3)
	        D2 = D2+ FI(M1,M2,M3)**2
	      END DO
       END DO
    END DO
      write(*,'(a,es12.4)') ' Npart = ',SS
      write(*,*) ' Density=',SS/float(NGRID)**3
      If(iFlag ==1)Then
         W   = FLOAT(NGRID)**3/SS
      Else
         W  = FLOAT(NGRIDold)**3/float(Ngalaxies)
      End If
                  
      write(18,*) ' Weight = ',W
      Ntot = INT(SS)
        SS =0.
        D2 =0.
!$OMP PARALLEL DO DEFAULT(SHARED) &           !-- rescale density
!$OMP PRIVATE (M1,M2,M3) Reduction(+:SS,D2)
      DO M3=1,NGRID
       DO M2=1,NGRID
          DO M1=1,NGRID
             FI(M1,M2,M3) = FI(M1,M2,M3) *W -1.
	        SS = SS+ FI(M1,M2,M3)
	        D2 = D2+ FI(M1,M2,M3)**2
	      END DO
       END DO
    END DO
    SS  = SS/float(NGRID)**3
    write(*,*) ' Average density =',SS/float(NGRID)**3
    write(*,*) ' RMS     density =',sqrt(D2/float(NGRID)**3 -SS**2)
    !write(18,*) ' Average density =',SS/float(NGRID)**3
    !write(18,*) ' RMS     density =',sqrt(D2/float(NGRID)**3 -SS**2)


       end subroutine DensSubbox

!----------------------------------------------------
!           power spectrum   P(k) for tracers
!
!
SUBROUTINE GetPowerSSC
!
!----------------------------------------------------
  use Tools
  use Density
  use Power
  use Parameters
  Real*4      :: kNyq
  Character*120 :: Name
    Allocate(Pk(NGRID,0:6),dNharm(NGRID,0:6),kharm(NGRID,0:6))

    Pk(:,:)     = 0.
    dNharm(:,:) = 0.
    kharm(:,:)  = 0.
    boxsize = Box
    kNyq = 2.*3.1415926/Box*NGRID/2.
    ReScalePk = (float(NGRIDnew)/NGRIDold)**3
    Box_old =Box
    Box = Box*(float(NGRIDnew)/NGRIDold)
       CALL POWERfft5(0,0)
       WRITE (12,'(3x,a,T9,a,T19,a,T30,a,T43,a)') 'bin','k/Mpch',      &
                 'N','Preal'
       DO I=1,NGRID/2
         WRITE (12,'(i6,f9.5,es11.4,7es13.5)')I,kharm(I,0),dNharm(I,0),Pk(I,0)
      EndDo
      

    DeAllocate(Pk,dNharm,kharm)
    close(18)
    Box = Box_old
  end SUBROUTINE GetPowerSSC

  
