  !--------------------------------------------------
  !
  !           Correlation function of mock  catalog
  !
  !--------------------------------------------------
Module Corr
Integer*4,  PARAMETER  ::                         & 
                     Nrad     =   500         ! Number of radial bins
Real*4        ::                    & 
                     dLogR,         &         ! size of log binning 
                     Rmax,          &         ! maximum radius of xi(R)
                     Mmin,          &         ! minimum mass of halos
                     Vmin                     ! minimum vcirc of halos
Real*4 ::            Xleft,Xright,Yleft,Yright,Zleft,Zright,dBuffer    ! boundaries of domain
Integer*8 ::         Np,iPartMax                      ! Number of particles in domain
!-------------------------  Main linker-list --------------
Integer*8,  ALLOCATABLE :: Lst(:)                                   ! Linker list
Integer*8,  ALLOCATABLE :: Label(:,:,:)                         ! Head of zero-level LL 
Real*4,     ALLOCATABLE :: MassHalo(:) 
Integer*4 ::                    Nmx,Nmy,Nmz,Nbx,Nby,Nbz         ! limits for linker-list
Real*4    ::                    Cell  , t0                      ! Size in Mpch of a linker-list cell
Real*8    ::                    Ncorr(-Nrad:Nrad)
character :: Letter*1            !--- for halo catalogs
Integer*4 :: iSnap
Contains
!---------------------------------------------------------------------------
!                     Define size and boundaries of linked-list
      SUBROUTINE SizeList
!---------------------------------------------------------------------------
       Cell       = Rmax
       Nmx        = (Xleft  -dBuffer)/Cell - 1
       Nmy        = (Yleft  -dBuffer)/Cell - 1
       Nmz        = (Zleft  -dBuffer)/Cell - 1
       Nbx        = (Xright +dBuffer)/Cell + 1
       Nby        = (Yright +dBuffer)/Cell + 1
       Nbz        = (Zright +dBuffer)/Cell + 1
      write(*,'(a,g12.3,a,6i6)') '     Size List:  Cell=',Cell, &
                                 ' Limits=',Nmx,Nbx,Nmy,Nby,Nmz,Nbz
      Allocate(Lst(Np),Label(Nmx:Nbx,Nmy:Nby,Nmz:Nbz))
    end SUBROUTINE SizeList
!--------------------------------------------------------------
!                          Make linker lists of particles in each cell
      SUBROUTINE List
!--------------------------------------------------------------
use Tools
        Integer*8 :: Nm,N0,N1,N2,N3,iMax,jp
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (i)
             Do jp=1,Np
                Lst(jp)=-1
             EndDo
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (i,j,k)
             Do k=Nmz,Nbz
             Do j=Nmy,Nby
             Do i=Nmx,Nbx
                Label(i,j,k)=0
             EndDo
             EndDo
          EndDo
          write(*,*) ' List: done init '
       
      Do jp=1,Np
         i=Ceiling(Xpar(jp)/Cell)-1
         j=Ceiling(Ypar(jp)/Cell)-1
         k=Ceiling(Zpar(jp)/Cell)-1
         i=MIN(MAX(Nmx,i),Nbx)
         j=MIN(MAX(Nmy,j),Nby)
         k=MIN(MAX(Nmz,k),Nbz)
         Lst(jp)      =Label(i,j,k)
         Label(i,j,k) =jp
      EndDo
      write(*,*) ' List: done all '
    end SUBROUTINE List
!---------------------------------------------------------------------------
!              find limits for the linker-list search
      SUBROUTINE Limits(x,y,z,Radius,i1,i2,j1,j2,k1,k2)
!----------------------------------------------------------------------------
           i2=Ceiling((x+Radius)/Cell)-1  
           j2=Ceiling((y+Radius)/Cell)-1  
           k2=Ceiling((z+Radius)/Cell)-1 
           i1=Ceiling((x-Radius)/Cell)-1  
           j1=Ceiling((y-Radius)/Cell)-1  
           k1=Ceiling((z-Radius)/Cell)-1 
 
            i1=MIN(MAX(Nmx,i1),Nbx) 
            j1=MIN(MAX(Nmy,j1),Nby)
            k1=MIN(MAX(Nmz,k1),Nbz)
           i2=MIN(MAX(Nmx,i2),Nbx)
           j2=MIN(MAX(Nmy,j2),Nby)
           k2=MIN(MAX(Nmz,k2),Nbz)

         end SUBROUTINE Limits
!---------------------------------------------------------------------------- 
!         Add buffer around the computational box
!           -- 
!           -- 
!           -- move coordinates to new arrays
    
Subroutine AddBuffer
!---------------------------------------------------------------------------- 
use Tools
     Real*4,       ALLOCATABLE,   DIMENSION(:) ::               &
                     Xbb,   Ybb,  Zbb,       &   !   coords
                     VXbb, VYbb, Vzbb,       &   !    velocities
                     MassHalob
     Integer*8 :: idummy,ic,ip,iPartMax,i

     Factor   = (1.+4.*dBuffer/Box)**3 !*2 !!! remove *2
     iPartMax = Factor*Np               ! reserve extra space for total N particles
     Nparticles  = Np                   ! old number of particles
     write(*,*) ' Allocate buffers for particles. N=',iPartMax
       ALLOCATE(Xbb(iPartMax),Ybb(iPartMax),Zbb(iPartMax))
       ALLOCATE(VXbb(iPartMax),VYbb(iPartMax),VZbb(iPartMax))
       ALLOCATE(MassHalob(iPartMax))
       
       Xleft = 0. ; Xright = Box
       
           Yleft = 0. ; Yright = Box
           Zleft = 0. ; Zright = Box

!$OMP PARALLEL DO DEFAULT(SHARED)  PRIVATE (ic)
       Do ic=1,iPartMax
                    VXbb(ic)= 0.
                    VYbb(ic)= 0.
                    VZbb(ic)= 0.
          Xbb(ic)= 0.
          Ybb(ic)= 0.
          Zbb(ic)= 0.
          MassHalob(ic) = 0.
       End Do
       
       ip = 0
       Xl = Xleft  -dBuffer ; Xr = Xright +dBuffer
       Yl = Yleft  -dBuffer ; Yr = Yright +dBuffer
       Zl = Zleft  -dBuffer ; Zr = Zright +dBuffer
       write(*,*) '    Replicate coordinates'
!$OMP PARALLEL DO DEFAULT(SHARED)  PRIVATE (ic,xx,yy,zz)
       Do ic=1,Np
          Xbb(ic)= Xpar(ic)
          Ybb(ic)= Ypar(ic)
          Zbb(ic)= Zpar(ic)
          MassHalob(ic) = MassHalo(ic)
       End Do
       write(*,*) '    Replicate velocities'
!$OMP PARALLEL DO DEFAULT(SHARED)  PRIVATE (ic,xx,yy,zz)
       Do ic=1,Np
          VXbb(ic)= VX(ic)
          VYbb(ic)= VY(ic)
          VZbb(ic)= VZ(ic)
       End Do
              write(*,*) '    add buffer'

       ip = Np
       Do ic=1,Np
          xx = Xpar(ic) ; yy =Ypar(ic) ; zz =Zpar(ic)
          do k = -1,1
             z = zz+k*Box
             if(z<Zl.or.z>Zr)cycle
               do j = -1,1
                 y = yy+j*Box
                 if(y<Yl.or.y>Yr)cycle
                   do i = -1,1
                      x = xx+i*Box
                      if(x<Xl.or.x>Xr)cycle
                      if((x<Xleft.or.x>Xright) .or.    &
                           (y<Yleft.or.y>Yright) .or.    &
                           (z<Zleft.or.z>Zright))Then
                      ip = ip +1
                      If(ip>iPartMax)Stop ' Too many buffer particles. Increase Factor in AddBuffer'
                      Xbb(ip) = x
                      Ybb(ip) = y
                      Zbb(ip) = z
                      Masshalob(ip) = MassHalo(ic)
                      VXbb(ip)= VX(ic)
                      VYbb(ip)= VY(ic)
                      VZbb(ip)= VZ(ic)
                   end if
                  end do
              end do
           end do
        end Do
      write(*,*) ' New number of particles = ',ip
      DEALLOCATE(Xpar,Ypar,Zpar)  
      DEALLOCATE(VX,VY,VZ)
      DEALLOCATE(MassHalo)
      Np = ip
      ALLOCATE(Xpar(Np),Ypar(Np),Zpar(Np))  
      ALLOCATE(VX(Np),VY(Np),VZ(Np))  
      ALLOCATE(MassHalo(Np))

!$OMP PARALLEL DO DEFAULT(SHARED)  PRIVATE (ip)
      Do ip =1,Np
         Xpar(ip) = Xbb(ip) ; Ypar(ip) = Ybb(ip) ; Zpar(ip) = Zbb(ip)
         VX(ip)= VXbb(ip); VY(ip)= VYbb(ip); VZ(ip)= VZbb(ip)
         MassHalo(ip) = MassHalob(ip)
      EndDo
      DEALLOCATE(Xbb,Ybb,Zbb)  
      DEALLOCATE(VXbb,VYbb,VZbb)
      DEALLOCATE(MassHalob)


      xmin = 1.e12; xmax = -1.e12
      ymin = 1.e12; ymax = -1.e12
      zmin = 1.e12; zmax = -1.e12

!$OMP PARALLEL DO DEFAULT(SHARED)  PRIVATE (ip)  &
!$OMP REDUCTION(MAX:xmax,ymax,zmax) REDUCTION(MIN:xmin,ymin,zmin) 
           Do ip =1,Np
              xmin =MIN(xmin,Xpar(ip))
              ymin =MIN(ymin,Ypar(ip))
              zmin =MIN(zmin,Zpar(ip))
              xmax =MAX(xmax,Xpar(ip))
              ymax =MAX(ymax,Ypar(ip))
              zmax =MAX(zmax,Zpar(ip))
              
              If(Xpar(ip).lt.Xl.or.Ypar(ip).lt.Yl.or.Zpar(ip).lt.Zl)&
                write(*,'(a,i10,1p,3g14.5)')' Error coord: ',i,Xpar(i),Ypar(i),Zpar(i)
           EndDo
!     write(*,'(a,20es13.5)') ' X= ',(Xpar(i),i=1,10)    
!     write(*,'(a,20es13.5)') ' Y= ',(Ypar(i),i=1,10)    
!     write(*,'(a,20es13.5)') ' Z= ',(Zpar(i),i=1,10)
!     write(*,'(a,20es13.5)') ' X= ',(Xpar(i),i=Np-9,Np)    
!     write(*,'(a,20es13.5)') ' Y= ',(Ypar(i),i=Np-9,Np)    
!     write(*,'(a,20es13.5)') ' Z= ',(Zpar(i),i=Np-9,Np)
     
     write(*,'(a,2es13.5)') ' X range= ',xmin,xmax
     write(*,'(a,2es13.5)') ' Y range= ',ymin,ymax
     write(*,'(a,2es13.5)') ' Z range= ',zmin,zmax
    end Subroutine AddBuffer

!---------------------------------------------------------------------------
!                   
!
SUBROUTINE CorrGet
!---------------------------------------------------------------------------
use Tools
  integer*8 :: ic,ip,i
  Integer OMP_GET_THREAD_NUM, OMP_GET_MAX_THREADS
  Real*8, Allocatable :: Ncc(:,:)
  Integer*8, Allocatable :: iNcc(:,:)
  Integer*8 :: iCorr(-1000:1000)
  Real*4 :: MassH1,MassH2
  Real*8 :: Msum

  Ncorr(:) = 0
  iCorr(:) = 0
     Nthr     = OMP_GET_MAX_THREADS()
     write(*,*) ' Number of Threads =',Nthr
     Allocate(Ncc(-Nrad:Nrad,Nthr))
     Allocate(iNcc(-Nrad:Nrad,Nthr))
     Ncc(:,:)  = 0
     iNcc(:,:) = 0
     R2max    = Rmax**2
     Radius   = Rmax
     tstart   = seconds()
     write(*,'(2(a,i12))') ' Nparticles= ',Nparticles,' Nparticles total= ',Np
! 
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (ip,x,y,z,jp,dd,k1,k2,k3,j1,j2,j3,i1,i2,i3,ind,mThr,MassH1) 
     Do ip=1, Nparticles
        mThr = OMP_GET_THREAD_NUM()+1
        x   = Xpar(ip);   y = Ypar(ip);    z = Zpar(ip)
        MassH1 =  MassHalo(ip)
            Call Limits(x,y,z,Radius,i1,i2,j1,j2,k1,k2)
            if(mod(ip,10000)==0.and.mThr==1)write(*,'(a,i11,3x,12i4)') '  ip=',ip,i1,i2,j1,j2,k1,k2
            Do k3 =k1, k2
            Do j3 =j1, j2
            Do i3 =i1, i2
              jp =Label(i3,j3,k3)
              Do while (jp.ne.0)
                If(jp/=ip)Then
                  dd =(x-Xpar(jp))**2+(y-Ypar(jp))**2+(z-Zpar(jp))**2
                  If(dd.lt.R2max)Then
                     ind = max(-Nrad,min(Nrad,INT(Log10(sqrt(dd))/dLogR+1000.)-1000))
                !write(*,'(i9,es12.4,i6,2es13.4)') ip,sqrt(dd),ind,10.**((ind)*dLogR),10.**((ind+1)*dLogR)
                     Ncc(ind,mThr)  = Ncc(ind,mThr)  +  MassH1 *MassHalo(jp)
                     iNcc(ind,mThr) = iNcc(ind,mThr) +  1
                  End If
               end If
               jp =Lst(jp)
              end do
           end do
           end do
           end Do
       End Do         ! ip
       tfinish = seconds()
      write(*,'(10x,a,T50,2f10.2)') ' time for CorrGet =',tfinish-tstart,tfinish-t0

      Do i=-Nrad,Nrad
         Do j=1,Nthr
            Ncorr(i) = Ncorr(i) +Ncc(i,j)
            iCorr(i) = iCorr(i) +iNcc(i,j)
         End Do
      EndDo
      Msum =0.
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (ip) REDUCTION(+:Msum) 
     Do ip=1, Nparticles
        Msum = Msum +  MassHalo(ip)
     EndDo
      
     !fact = Box**3/Nparticles**2
     fact = Box**3/Msum**2
     factN = Box**3/(float(Nparticles))**2
     write(*,'(a,es12.4)') ' Average mass =',Msum/Nparticles
      write(12,'(3x,a)') ' R            dN           xi_mass(R)   R**2xi_mass(R) xi_numb    rR**2xi_numb'
      Do i =-Nrad+1,Nrad
         R1 = 10.**((i)*dLogR)
         R2 = 10.**((i+1)*dLogR)
         dV = 4.188*(R2**3-R1**3)
         xi = Ncorr(i)*fact/dV -1.              ! 2Npairs/(dV/V) -1
         xN = iCorr(i)*factN/dV -1.              ! 2Npairs/(dV/V) -1
         Rr = (R1+R2)/2.
         if(R2<Rmax.and.Rr>1.00)write(12,'(f9.4,3x,i9,2(3x,2es12.4))') &
                                    Rr,iCorr(i),xi,xi*Rr**2,xN,xN*Rr**2
      EndDo
      deallocate(Ncc)
    end SUBROUTINE CorrGet

!---------------------------------------------------------
!                         Read BDM  halos catalog
!                       
      SUBROUTINE ReadHalos
!--------------------------------------------------------------
use Tools
     Character     :: file1*80,Line*120,txt*80
     Logical       :: inside
     Real*8        :: Massav,Vcmin,Vcav
     NpMax = 5e7
    ALLOCATE(Xpar(NpMax),Ypar(NpMax),Zpar(NpMax),VX(NpMax),VY(NpMax),VZ(NpMax))
    ALLOCATE(MassHalo(NpMax))
    Ncenter  =0
    write(*,*) ' Mmin = ',Mmin,' Rmax = ',Rmax
                    ! --- read header of Catshort catalog
   Call GetFile(0)
   read(3,'(a)')HEADER
   write(*,'(a)')HEADER
   !write(12,'(a)')HEADER

   If(Vmin<10..or.Vmin>1000.)Then   !-- selection by Mass_tot
      Vmin  = 0.
   else                          !-- selection by Vcirc
      Mmin  = 0.
   EndIf
   Do i=1, 7
     Read(3,'(a)')Line 
     if(i==1.or.(i>2.and.i<14))write(12,'(a)') Line   
     write(*,'(a)') Line 
   EndDo
   iFile   = 1 
   Massav =0. ;  Vcmin =10000. ; Vcav =0.
   Np = 0
      Ndistinct = 0
   Do     !i =1,18000000            ! -------- read halos
     !read (3,*,iostat=iStat) X0, Y0 , Z0, VvX,VvY,VvZ,       &
     !                        aM,rr,vrmss,vcrc,vacc,aMacc,aAcc,iBDM,Conc,    &
     !                        iDist,Offset,VirRat,aLam,Rin,       &
     !                        ba,ca,xax,yax,zax
     read (3,*,iostat=iStat) X0, Y0 , Z0, VvX,VvY,VvZ,       &
                             aM,aMtot,rr,vrmss,vcrc,iHalo,Conc,aNpart,    &
                             iDist,Offset,VirRat,aLam,Rin,       &
                             ba,ca,xax,yax,zax
     !read (3,*,iostat=iStat) ii,X0, Y0 , Z0, &
     !                                          aM0,aM,vnow,vcrc
     !rr = aM0**0.33333/4.82e1      !  virial radius in kpch for Om=0.27 delta=360.8
     Np = Np +1
     If(mod(Np,10000)==0)&
      write (*,'(a,i11,i9,3f10.4,2g12.3,3f9.2)') ' Reading: ',Np,Ncenter, X0, Y0 , Z0,aM,aMtot,vcrc,Conc !,aM,rr,vnow,vcrc
45          Format(F9.4,2F9.4,3F8.1,g11.3,f8.2,2f7.1,I9,f8.1,f8.2,i4)

        If(iStat.ne.0)Call GetFile(iFile)
        If(iFile ==0)Exit
        If(iDist/=0)vcrc =vcrc*1.15
        If(aMtot> Mmin.and.vcrc>Vmin)Then    !!! selection conditions
        Ncenter = Ncenter +1
        iHalo = Ncenter
        Massav = Massav + aMtot
        Vcmin  = min(Vcmin,vcrc)
        Vcav   = Vcav + vcrc
        !Npart = INT(aNpart)
        !If(Ncenter /= iHalo)  &
        !     Call HandleError(' --- Catshort :',' Error in sequence of halos:',&
        !     iHalo,X0, Y0 , Z0, VvX,VvY,VvZ,aM,rr,vrmss,vcrc,Ncenter,Conc)
        If(X0 > Box .or. Y0 > Box .or. Z0 > Box)  &
             Call HandleError(' --- Catshort :',' Error in upper limit of coordinates:',&
             iHalo,X0, Y0 , Z0, VvX,VvY,VvZ,aM,rr,vrmss,vcrc,Npart,Conc)
        If(X0<0. .or. Y0 < 0. .or. Z0<0.)  &
             Call HandleError(' --- Catshort :',' Error in low limit of coordinates:',&
             iHalo,X0, Y0 , Z0, VvX,VvY,VvZ,aM,rr,vrmss,vcrc,Npart,Conc)
        If(aM>1.e16.or.aM<1.e7) &
             Call HandleError(' --- Catshort :',' Error in halo mass     :',&
             iHalo,X0, Y0 , Z0, VvX,VvY,VvZ,aM,rr,vrmss,vcrc,Npart,Conc)
        If(vcrc<0.or.vcrc>10000.) &
             Call HandleError(' --- Catshort :',' Error in halo Vcirc     :',&
             iHalo,X0, Y0 , Z0, VvX,VvY,VvZ,aM,rr,vrmss,vcrc,Npart,Conc)
        !If(aNpart<0.or.aNpart>1.e10) &
        !     Call HandleError(' --- Catshort :',' Error in number of particles:',&
        !     iHalo,X0, Y0 , Z0, VvX,VvY,VvZ,aM,rr,vrmss,vcrc,Npart,Conc)
        If(iDist<0 .or.iDist > 2000000000) &
             Call HandleError(' --- HaloList Error :',' Error in distinct index:',&
             iHalo,X0, Y0 , Z0, VvX,VvY,VvZ,aM,rr,vrmss,vcrc,iDist,Conc)
        If(Ncenter .le. NpMax)Then
              Xpar(Ncenter) = X0
              Ypar(Ncenter) = Y0
              Zpar(Ncenter) = Z0
              VX(Ncenter) = VvX
              VY(Ncenter) = VvY
              VZ(Ncenter) = VvZ
              MassHalo(Ncenter) = aMtot
        EndIf
        EndIf     !  test mass
   Enddo
40          Format(g11.4,I9,g11.3,g10.3,2f7.1,g10.3,i8,i7,i5,g11.3)
   close (3)
   If(Ncenter == 0)Stop ' !---- No halos found'
   Vcav   = Vcav/Ncenter
   Massav = Massav/Ncenter
    write(*,'(2x,a,i11,3(a,es12.4))') ' Read Halos = ',Ncenter,' Average mass=',Massav, &
                                                    ' Vcirc min=',Vcmin,' Vcirc aver=',Vcav
    write(*,'(2x,a,es12.3,3(a,es12.4))') ' Selection:  Mass > ',Mmin,' Vcirc > ',Vmin
    write(12,'(2x,a,i11,3(a,es12.4))') ' Read Halos = ',Ncenter,' Average mass=',Massav, &
                                                    ' Vcirc min=',Vcmin,' Vcirc aver=',Vcav
    write(12,'(2x,a,es12.3,3(a,es12.4))') ' Selection: Mass >',Mmin,' Vcirc >',Vmin
    If(Ncenter >= NpMax)Stop ' Not enough space for halos. Increase NpMax.'
    Nparticles = Ncenter
    Np         = Nparticles
      Return
10    write(*,*) ' Could not read radial bins'
      stop
20    write(*,*) ' Could not read text R(kpc/h)'
      stop
30    write(*,*) ' Unexpected End of file'
      stop

    End SUBROUTINE ReadHalos

!---------------------------------------------------------
!                         Read RockStar  halos catalog
!                       
      SUBROUTINE ReadRockHalos
!--------------------------------------------------------------
use Tools
     Character     :: file1*80,FileName*120,Line*120,txt*80
     Logical       :: inside,FileExists
     Real*8        :: Massav,Vcmin,Vcav,dum(11),aMtot, m200b, m200c, m500c, m2500c
     Integer*8     :: i0, i1, i2, i3,i4,i5,i6,i7,i8, nextc, lastc, lastm, llm,id1
    ALLOCATE(Xpar(NpMax),Ypar(NpMax),Zpar(NpMax),VX(NpMax),VY(NpMax),VZ(NpMax))

    Ncenter  =0
    write(*,'(a,$)') ' Enter Min mass and min Vpeak = '
    read(*,*) Mmin,Vmin
    write(*,*) ' Mmin = ',Mmin,' Rmax = ',Rmax
    write(*,'(a,$)') ' Enter input file name           => '
    read(*,'(a)') FileName
    Inquire(file=TRIM(FileName), Exist = FileExists)
    If(.not.FileExists)Then
       write(*,'(2a)') ' Input file does not exit. Check file name:',TRIM(FileName)
       stop
    End If

  Open(3,file=TRIM(FileName))
  FileName = TRIM(FileName)//'.corr'
  Open(12,file=TRIM(FileName))


    ! --- read header 
   If(Vmin<10..or.Vmin>1000.)Then   !-- selection by Mass_tot
      Vmin  = 0.
   else                          !-- selection by Vcirc
      Mmin  = 0.
   EndIf
   Do i=1, 63
     Read(3,'(a)')Line 
     write(*,'(a)') Line 
   EndDo
   Massav =0. ;  Vcmin =10000. ; Vcav =0.
   Np = 0
      Ndistinct = 0
   Do    !i =1,2800000           ! -------- read halos
      !read (3,*,iostat=iStat) Aexpn, i0, a0, i1, i2, i3, i4,i5, i6,        &  ! 1-9   scale(0) id(1) desc_scale(2) desc_id(3) num_prog(4) pid(5) upid(6) desc_pid(7) phantom(8)
      !                       d0, aMvir, Rvir, Rs, Vrms, mmp, a1, Vmax,      &  ! 10-17 sam_mvir(9) mvir(10) rvir(11) rs(12) vrms(13) mmp?(14) scale_of_last_MM(15) vmax(16)
      !                       X0, Y0 , Z0, VvX,VvY,VvZ,  aJx, aJy, aJz,     &  ! 18-26 x(17) y(18) z(19) vx(20) vy(21) vz(22) Jx(23) Jy(24) Jz(25)                             
      !                       alam, id1, i7, i8, isnap, nextc, lastc, lastm, llm, &  ! 27-35 Spin(26) (27) (28) ID(29) (30) (31) NextID(32) LastID(33) Last_mainleafID(34)
      !                       Rsk, aMtot, m200b, m200c, m500c, m2500c,      &  ! 36-41  Rs_Klypin(35) Mmvir_all(36) M200b(37) M200c(38) M500c(39) M2500c(40) 
      !                       Xoff, Voff, dum(1:11),                        &  ! 42-54 Xoff, Voff
      !                       ToverU, b1, b2, b3, aMacc, aMpeak, Vacc, Vpeak
      read (3,*,iostat=iStat) aMvir,id1,Rvir,Vmax,X0, Y0 , Z0, VvX,VvY,VvZ,alam,Rs,aMtot,aM200b,aM200c,Xoff,ToverU,aMacc,aMpeak,Vacc,Vpeak
     Np = Np +1
     If(iStat.ne.0)exit
     If(mod(Np,1000000)==0)&
      write (*,'(a,i12,i11,4f10.4,2es12.3,4f9.2,i12,f9.4)') ' Reading: ',Np,Ncenter,Aexpn, X0, Y0 , Z0,aMvir,aMtot,Vmax,Vpeak,Rs,Rvir,id1,alam !,aM,rr,vnow,vcrc

        If(aMtot> Mmin.and.Vpeak>Vmin)Then    !!! selection conditions
        Ncenter = Ncenter +1
        iHalo = Ncenter
        Massav = Massav + aMtot
        Vcmin  = min(Vcmin,Vpeak)
        Vcav   = Vcav + Vpeak

        If(Ncenter .le. NpMax)Then
              Xpar(Ncenter) = X0
              Ypar(Ncenter) = Y0
              Zpar(Ncenter) = Z0
              VX(Ncenter) = VvX
              VY(Ncenter) = VvY
              VZ(Ncenter) = VvZ
           Else
              write(*,*) ' Not enough space for particles. Increase NpMax =',NpMax
          Stop  ' Not enough space for particles'    
        EndIf
        EndIf     !  test mass
   Enddo
40          Format(g11.4,I9,g11.3,g10.3,2f7.1,g10.3,i8,i7,i5,g11.3)
   close (3)
   If(Ncenter == 0)Stop ' !---- No halos found'
   Vcav   = Vcav/Ncenter
   Massav = Massav/Ncenter
    write(*,'(2x,a,i11,3(a,es12.4))') ' Read Halos = ',Ncenter,' Average mass=',Massav, &
                                                    ' Vcirc min=',Vcmin,' Vcirc aver=',Vcav
    write(*,'(2x,a,es12.3,3(a,es12.4))') ' Selection:  Mass > ',Mmin,' Vcirc > ',Vmin
    write(12,'(2x,a,i11,3(a,es12.4))') ' Read Halos = ',Ncenter,' Average mass=',Massav, &
                                                    ' Vcirc min=',Vcmin,' Vcirc aver=',Vcav
    write(12,'(2x,a,es12.3,3(a,es12.4))') ' Selection: Mass >',Mmin,' Vcirc >',Vmin
    If(Ncenter >= NpMax)Stop ' Not enough space for halos. Increase NpMax.'
    Nparticles = Ncenter
    Np         = Nparticles
    
      xmin = 1.e12; xmax = -1.e12
      ymin = 1.e12; ymax = -1.e12
      zmin = 1.e12; zmax = -1.e12
      Xl = 0.
      Yl = 0.
      Zl = 0.
!$OMP PARALLEL DO DEFAULT(SHARED)  PRIVATE (ip)  &
!$OMP REDUCTION(MAX:xmax,ymax,zmax) REDUCTION(MIN:xmin,ymin,zmin) 
           Do ip =1,Np
              xmin =MIN(xmin,Xpar(ip))
              ymin =MIN(ymin,Ypar(ip))
              zmin =MIN(zmin,Zpar(ip))
              xmax =MAX(xmax,Xpar(ip))
              ymax =MAX(ymax,Ypar(ip))
              zmax =MAX(zmax,Zpar(ip))
              
              If(Xpar(ip).lt.Xl.or.Ypar(ip).lt.Yl.or.Zpar(ip).lt.Zl)&
                write(*,'(a,i10,1p,3g14.5)')' Error coord: ',ip,Xpar(ip),Ypar(ip),Zpar(ip)
           EndDo
     write(*,'(a,20es13.5)') ' X= ',(Xpar(i),i=1,10)    
     write(*,'(a,20es13.5)') ' Y= ',(Ypar(i),i=1,10)    
     write(*,'(a,20es13.5)') ' Z= ',(Zpar(i),i=1,10)
     write(*,'(a,20es13.5)') ' X= ',(Xpar(i),i=Np-9,Np)    
     write(*,'(a,20es13.5)') ' Y= ',(Ypar(i),i=Np-9,Np)    
     write(*,'(a,20es13.5)') ' Z= ',(Zpar(i),i=Np-9,Np)
     
     write(*,'(a,2es13.5)') ' X range= ',xmin,xmax
     write(*,'(a,2es13.5)') ' Y range= ',ymin,ymax
     write(*,'(a,2es13.5)') ' Z range= ',zmin,zmax
      Return
10    write(*,*) ' Could not read radial bins'
      stop
20    write(*,*) ' Could not read text R(kpc/h)'
      stop
30    write(*,*) ' Unexpected End of file'
      stop

    End SUBROUTINE ReadRockHalos    
!--------------------------------------------------------------
!                       Reads and recognizes files
!  
      SUBROUTINE GetFile(iFile)
!--------------------------------------------------------------
use Tools
     Integer*4 , parameter :: Nbin = 1000
     Integer*4, SAVE  :: iFileout,iFlag
     Logical    :: exst
     Character*120 :: file1,file2,catname,listname,line,line2     

   Open(20,file='ErrorAnalyzeCatalogs.dat',position='append')

   If(iFile == 0) Then    ! recongize configuration

      write(*,'(a,$)')' Enter snapshot, realization, Catalog letter, Minimum halo mass , Minimum Vcirc => '
      read(*,*)iSnap,ir,Letter,Mmin,Vmin

       If(Vmin<10. .or. Vmin >1000.)Then
         write(file2,'(a,2a1,i4.4,a6,i3.3,a,i4.4,a)')'HaloCorrelation',TRIM(Letter),'.', &
              iSnap,'.mass.',INT(10.*log10(Mmin)),'.',ir,'.DAT'
       else
         write(file2,'(a,2a1,i4.4,a6,i4.4,a)')'HaloCorrelation',TRIM(Letter),'.', &
              iSnap,'.vcir.',INT(Vmin+0.5),'.DAT'
       end If
         Open(12,file=TRIM(file2))                                 ! output
         iFileout = ir
         write(file1,'(a,2a1,i4.4,a1,i4.4,a)')'Catshort',TRIM(Letter),'.',iSnap,  &
                                              '.',iFileout,'.DAT'
         write(*,'(/a)') TRIM(file1)
         Inquire(file=TRIM(file1),exist = exst)
         write(*,*) ' Status =',exst
         If(exst)Then              ! multiple files CatshortX.XXX.XX.DAT
            iFlag = 0
            open(3,file = TRIM(file1))
         Else
            iFlag = 1              !---  single file CatshortX.xxx.DAT
            write(file1,'(a,2a1,i4.4,a)')'Catshort',TRIM(Letter),'.',iSnap,  &
                                                 '.DAT'      
            Inquire(file=TRIM(file1),exist = exst)
         write(*,'(a)') TRIM(file1)
         write(*,*) ' Status =',exst
            If(exst)Then
               open(3,file = TRIM(file1))               
            Else
               Stop '----  Input Catshort file not found '
            End If
         End If
      read(3,'(a)')HEADER
      !read(3,'(2(a7,f8.5))')  line(1:7),AEXPN,Line2(1:7),Step 
      read(3,'(2(a8,f8.5))')  line(1:8),AEXPN  !,Line2(1:7),Step 
      write(*,*) AEXPN
      rewind(3)
Else

      If(iFlag == 1)Then  ! single file
         iFile = 0
         Return
      Else
         iFileout = iFileout +1
         write(*,*)  ' -------- open file =',iFileout
         write(file1,'(a,2a1,i4.4,a1,i2.2,a)')'Catshort',TRIM(Letter),'.',iSnap,  &
              '.',iFileout,'.DAT'
         write(*,*) TRIM(file1)
         Inquire(file=TRIM(file1),exist = exst)
         If(exst)Then              ! multiple files CatshortX.XXX.XX.DAT
            open(3,file = TRIM(file1))
            write(*,'(a)') TRIM(file1)
         Else                      ! no more files to read
            iFile = 0
            Return
         End If
      End If
End If
      end SUBROUTINE GetFile
!---------------------------------------------------------------
!             
!                       
      SUBROUTINE HandleError(Text1,Text2,iHalo,X0, Y0 , Z0, & 
                             VvX,VvY,VvZ,aM,rr,vrmss,vcrc,Nhalo,Conc)
!---------------------------------------------------------------
   character(*) :: Text1, Text2
   integer :: ListFiles(2) =(/0,20/)
   Do j=1,2
      write (ListFiles(j),'(2a,/3x,a,i9,a,1p,3g12.4,a,3g12.4, &
            /3x,2(a,g12.4),g12.4,a,0p,i9,a,1p,g12.4)') &
             Trim(Text1),TRIM(Text2), &
           ' Line in file=',iHalo,' Coords=', X0, Y0 , Z0, & 
           ' Veloc=',VvX,VvY,VvZ,        & 
           ' M=',aM,' R=',rr,            &
           ' Vrms=',vrmss,vcrc, ' Nhalo=',Nhalo,&
           ' C =',Conc
   EndDo
      stop
      End SUBROUTINE HandleError
    
!---------------------------------------------------------------------------
!                   
!
SUBROUTINE ReadMock
!---------------------------------------------------------------------------
use Tools
  character*150 :: Line,FileName
  logical FileExists

    ALLOCATE(Xpar(NpMax),Ypar(NpMax),Zpar(NpMax),VX(NpMax),VY(NpMax),VZ(NpMax))
    ALLOCATE(Dens(NpMax))

    write(*,'(a,$)') ' Enter input file name           => '
    read(*,'(a)') FileName
    Inquire(file=TRIM(FileName), Exist = FileExists)
    If(.not.FileExists)Then
       write(*,'(2a)') ' Input file does not exit. Check file name:',TRIM(FileName)
       stop
    End If

  Open(10,file=TRIM(FileName))
  FileName = TRIM(FileName)//'.corr'
  Open(12,file=TRIM(FileName))

    write(*,*)
    do i=1,2              !--- header
      Read(10,'(a)') Line
      write(12,'(a)') TRIM(Line)
      write(*,'(a)') TRIM(Line)
    end do
    Read(10,'(a)') Line
    Np = 0
    do
       read(10,*,iostat=io)x,y,z,vvx,vvy,vvz,den
       if(io/=0)exit
       Np = Np+1
       If(Np>Npmax) Stop ' Not enough space for particles. Increase NpMax'
       Xpar(Np) = x ; Ypar(Np) = y ; Zpar(Np) = z ;
       VX(Np) = vvx ; VY(Np) = vvy ; VZ(Np) = vvz ;
       Dens(Np) = den
    EndDo
    Nparticles = Np
    write(*,*) ' Read N particles = ',Nparticles
    close(10)
  end SUBROUTINE ReadMock
!
end Module Corr
!---------------------------------------------
Program CorrFunction
  use Corr
  use Tools
  character*120 :: FileName
  logical        FileExists  

  t0 = seconds()
  write(*,'(a,$)') ' Enter Box size and Max Radius   =  '
  read(*,*) Box,Rmax
  

  dBuffer = Rmax
  dLogR   = 0.03
  !Call ReadRockHalos
  !Call ReadMock
  Call ReadHalos
  Call AddBuffer
  !Call GetDMparticles(moment,'',fraction)
  !Call RandomDMparticles(2000000)
  Call SizeList
  Call List
  Call CorrGet
  
end Program CorrFunction

!!--------------------------------------------
!subroutine RandomDMparticles(N)
!!--------------------------------------------
!  use Tools
!  use Corr
!  use Random
!  INCLUDE 'luxuryp.h'
!  logical :: exst, in1, in2
!  character*80 :: Filename
!  Real*4, Allocatable  :: Ggx(:),Ggy(:),Ggz(:)

!      Inquire(file='../Setup.dat',exist = exst)
!      if(.not.exst)Then
!         write(*,*)' Error: File ../Setup.dat not found. Run PMP2init.exe'
!         stop
!      end if
!         open(11,file='../Setup.dat')
!      write(*,'(a,$)') ' Enter moment            =  '
!      read(*,*) moment
!      write(*,'(a,$)') ' Enter particle fraction =  '
!      read(*,*) fraction

!      write(*,*) ' enter setup '
!      Call ReadSetup
!      Ngrid = 1000
!      AEXPN = 1.0
!      ISTEP = 0
!      Box   = 1000.
!      Nrecord = Ngrid**2
!      iPartMax= N*(1.+2.5*dBuffer/Box)**3
!      write(*,*) ' done setup ',iPartMax
!       Nseed  = 1210831          !--- initialize luxury 
!       lux    = 2
!       Ns     = Nseed
!       CALL rluxgo (lux, Ns, 0, 0)
!       Npages = (N-1)/Nrecord+1
!       Allocate(Ggx(Nrecord),Ggy(Nrecord),Ggz(Nrecord))
!       ALLOCATE(XPAR(iPartMax),YPAR(iPartMax),ZPAR(iPartMax))
!       write(*,*) ' Npages = ',Npages
!       write(*,*) ' Nrecord = ',Nrecord
!       iCount = 0
!       Do j=1,Npages
!        Call ranlux(Ggx,Nrecord)
!        Call ranlux(Ggy,Nrecord)
!        Call ranlux(Ggz,Nrecord)
!        write(*,'(4(3f8.4,3x),i9)')(Ggx(i),Ggy(i),Ggz(i),i=1,4),j
!          Do ip =1,Nrecord
!             iCount = iCount +1
!             If(ip>N)exit
!                If(iCount > iPartMax)STOP 'Attempt to read too many particles '
!                Xpar(iCount) = Ggx(ip)*Box
!                Ypar(iCount) = Ggy(ip)*Box
!                Zpar(iCount) = Ggz(ip)*Box
!          end Do ! ip
!       end Do
!           Np = iCount -1
!           Nparticles = iCount -1
!     ip = Np                 !--- add particles in buffer width dBuffer

!     write(*,'(a,6i12)') ' Np  =', Np,ip,iPartMax
!     write(*,'(a,f8.3)') ' Box =', Box,' dBuffer =',dBuffer
     
!     Xleft = 0. ; Xright = Box
!     Yleft = 0. ; Yright = Box
!     Zleft = 0. ; Zright = Box
!       Xl =  -dBuffer ; Xr = Box +dBuffer
!       Yl =  -dBuffer ; Yr = Box +dBuffer
!       Zl =  -dBuffer ; Zr = Box +dBuffer
!       Do ic=1,Np
!          xx = Xpar(ic) ; yy =Ypar(ic) ; zz =Zpar(ic)
!          do k = -1,1
!             z = zz+k*Box
!               do j = -1,1
!                 y = yy+j*Box
!                   do i = -1,1
!                     x = xx+i*Box
!                     in1 = (k==0).and.(j==0).and.(i==0)
!                     If(.not.in1)Then
!                       in2 = (x>Xl.and.x<Xr).and.(y>Yl.and.y<Yr).and.(z>Zl.and.z<Zr)
!                       if(in2)Then
!                         ip = ip +1
!                         If(ip>iPartMax)Stop ' Too many buffer particles. Increase Factor in AddBuffer'
!                         Xpar(ip) = x
!                         Ypar(ip) = y
!                         Zpar(ip) = z 
!                       end if ! in2
!                    end If ! in1
!                    end do ! i
!               end do ! j
!           end do ! k
!        end Do ! ic
!      write(*,*) ' New number of particles = ',ip
!      Np = ip           
!      write(FileName,'(a,i4.4,a,i8.8,a)') 'Corr.step.',ISTEP,'.N.',Np,'.dat'
!      open(12,name=TRIM(FileName))
!   end subroutine RandomDMparticles
!--------------------------------------------
subroutine GetDMparticles(moment,Path,fraction)
!--------------------------------------------
  use Tools
  use Corr
  character(len=*) :: Path
  character*80 :: FileName
     logical :: exst, in1, in2
     integer*8 :: ic,ip

      Inquire(file='../Setup.dat',exist = exst)
      if(.not.exst)Then
         write(*,*)' Error: File ../Setup.dat not found. Run PMP2init.exe'
         stop
      end if
         open(11,file='../Setup.dat')
      write(*,'(a,$)') ' Enter moment                     =  '
      read(*,*) moment
      write(*,'(a,$)') ' Enter particle fraction          =  '
      read(*,*) fraction
      write(*,'(a,$)') ' Enter RSD flag (0-no,1,2,3-xyz)  =  '
      read(*,*) iFlag

      write(*,*) ' enter setup '
      Call ReadSetup
      write(*,*) ' done setup '
      CALL ReadDataPMcoords(moment,Path,fraction,iFlag)

      write(*,'(a,T20,a,T30,i12,T45,a,i4)')' Start running: ', &
                      'Nparticles:',Nparticles,' Ngrid= ',Ngrid
      write(*,'(T20,a,T30,i12,T45,a,es12.4)')  'Step:',ISTEP,' Aexp =',AEXPN
      Np = Nparticles
      Xscale = Box/NGRID

      xmin = 1.e12; xmax = -1.e12
      ymin = 1.e12; ymax = -1.e12
      zmin = 1.e12; zmax = -1.e12
      Xl = 0.
      Yl = 0.
      Zl = 0.
!$OMP PARALLEL DO DEFAULT(SHARED)  PRIVATE (ip)  &
!$OMP REDUCTION(MAX:xmax,ymax,zmax) REDUCTION(MIN:xmin,ymin,zmin) 
      Do ip =1,Np
        ! Xpar(ip) =(Xpar(ip)-1.)*Xscale
        ! Ypar(ip) =(Ypar(ip)-1.)*Xscale
        ! Zpar(ip) =(Zpar(ip)-1.)*Xscale
              xmin =MIN(xmin,Xpar(ip))
              ymin =MIN(ymin,Ypar(ip))
              zmin =MIN(zmin,Zpar(ip))
              xmax =MAX(xmax,Xpar(ip))
              ymax =MAX(ymax,Ypar(ip))
              zmax =MAX(zmax,Zpar(ip))
              
              If(Xpar(ip).lt.Xl.or.Ypar(ip).lt.Yl.or.Zpar(ip).lt.Zl)&
                write(*,'(a,i10,1p,3g14.5)')' Error coord: ',ip,Xpar(ip),Ypar(ip),Zpar(ip)
           EndDo
     write(*,'(a,20es13.5)') ' X= ',(Xpar(i),i=1,10)    
     write(*,'(a,20es13.5)') ' Y= ',(Ypar(i),i=1,10)    
     write(*,'(a,20es13.5)') ' Z= ',(Zpar(i),i=1,10)
     write(*,'(a,20es13.5)') ' X= ',(Xpar(i),i=Np-9,Np)    
     write(*,'(a,20es13.5)') ' Y= ',(Ypar(i),i=Np-9,Np)    
     write(*,'(a,20es13.5)') ' Z= ',(Zpar(i),i=Np-9,Np)
     
     write(*,'(a,2es13.5)') ' X range= ',xmin,xmax
     write(*,'(a,2es13.5)') ' Y range= ',ymin,ymax
     write(*,'(a,2es13.5)') ' Z range= ',zmin,zmax
  
     ip = Np                 !--- add particles in buffer width dBuffer

     write(*,'(a,6i12)') ' Np  =', Np,ip,iPartMax
     write(*,'(2(a,f8.3))') ' Box =', Box,' dBuffer = ',dBuffer
     
     Xleft = 0. ; Xright = Box
     Yleft = 0. ; Yright = Box
     Zleft = 0. ; Zright = Box
       Xl =  -dBuffer ; Xr = Box +dBuffer
       Yl =  -dBuffer ; Yr = Box +dBuffer
       Zl =  -dBuffer ; Zr = Box +dBuffer
       Do ic=1,Np
          xx = Xpar(ic) ; yy =Ypar(ic) ; zz =Zpar(ic)
          do k = -1,1
             z = zz+k*Box
               do j = -1,1
                 y = yy+j*Box
                   do i = -1,1
                     x = xx+i*Box
                     in1 = (k==0).and.(j==0).and.(i==0)
                     If(.not.in1)Then
                       in2 = (x>Xl.and.x<Xr).and.(y>Yl.and.y<Yr).and.(z>Zl.and.z<Zr)
                       if(in2)Then
                         ip = ip +1
                         If(ip>iPartMax)Stop ' Too many buffer particles. Increase Factor in AddBuffer'
                         Xpar(ip) = x
                         Ypar(ip) = y
                         Zpar(ip) = z 
                       end if
                      end if
                    end do ! i
               end do ! j
           end do ! k
        end Do ! ic
      write(*,*) ' New number of particles = ',ip
      Np = ip
      write(12,'(a,f10.4,2(a,i12),a,i2)') ' Fraction =',fraction, &
           ' Nparticles = ',Nparticles,' Np = ',Np, &
           ' RSD =',iFlag
     
   end subroutine GetDMparticles
!---------------------------------------
!        Read    PMfiles
!             moment <0    use PMcrd.DAT, PMcrs0.DAT ...    
!             moment >= 0  use PMcrd.xxxx.DAT, PMcrs0,XXXX.DAT ..
!        iFlag =  0 -- real space
!                 1 -- redshift along z
    SUBROUTINE ReadDataPMcoords(moment,Path,fraction,iFlag)
!      
!---------------------------------------
  use Tools
  use Corr
  use Random
  INCLUDE 'luxuryp.h'
      Character*80 :: Name
        Logical      :: exst
        Integer*8    :: iCount,ii,ioff,ip
        Integer*8    :: Ngal,Jpage,NinPage
        Integer*8    :: idummy,ic,i
        Integer*4    :: moment,Nrecord
        Real*4, Allocatable  :: Gg(:)
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
      close(4)


      Nrecord = 1024**2
      Naccess = Nrecord*6 !*4
      xR      = NGRID +1
      boxsize = extras(100)
      Box     = boxsize
      Npages  = (Nparticles-1)/Nrecord+1 ! number of records
      Nlast   = Nparticles - (Npages-1)*Nrecord ! number of particles in last record
      Jpage   = Nrecord
      Xscale = Box/NGRID
      factorV = sqrt(AEXPN/(Om+AEXPN**3*OmL))/AEXPN*Xscale
      If(iFlag ==0)Then
         iFx = 0 ; iFy =0 ; iFz =0
      EndIf
      If(iFlag ==1)Then
         iFx = 1 ; iFy =0 ; iFz =0
      EndIf
      If(iFlag ==2)Then
         iFx = 0 ; iFy =1 ; iFz =0
      EndIf
      If(iFlag ==3)Then
         iFx = 0 ; iFy =0 ; iFz =1
      EndIf

      write(*,'(a,4i4,7es12.4)') ' iFlag =', iFlag,iFx,iFy,iFz,factorV,Xscale,AEXPN,Om,OmL
      Allocate(Gg(Nrecord))

      write(17,'(a,i10)') ' NROW   =',NROW
      write(17,'(a,i10)') ' Ngal   =',Nparticles
      write(17,'(a,i10)') ' Ngrid  =',Ngrid
      write(17,'(a,f10.1)') ' Box    =',Box


      Allocate (Xb(Nrecord),Yb(Nrecord),Zb(Nrecord))
      Allocate (VXb(Nrecord),VYb(Nrecord),VZb(Nrecord))

     Factor   = (1.+3.*dBuffer/Box)**3*fraction 
     Np       = Nparticles
     iPartMax = Factor*Np               ! reserve extra space for particles included in buffer
     write(*,*) ' Allocate buffers for particles. N=',iPartMax
      write(Name,'(a,i4.4,a,i2.2,a,i4.4,a)') 'Corr.step.',ISTEP,&
                '.RSD.',iFlag,'.',Nrealization,'.dat'
      open(12,name=TRIM(Name))
      WRITE (12,'(a,/10x,a,f8.4,4(a,i7))') HEADER,                 & 
                       ' a=',AEXPN, ' step= ',ISTEP,         &
                       ' Nrow= ', NROW, ' Ngrid=',NGRID
       ALLOCATE(XPAR(iPartMax),YPAR(iPartMax),ZPAR(iPartMax))
       myMemory =Memory(3_8*iPartMax)
       
       Nseed  = 1210831          !--- initialize luxury 
       lux    = 2
       Ns     = Nseed
       CALL rluxgo (lux, Ns, 0, 0)   

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
        
   Do ii =1,Npages                !--- read particles
      If(ii==Npages)Then
         NinPage = Nparticles -(ii-1)*JPAGE  ! # particles in the current page
      Else
         NinPage = JPAGE
      EndIf
      jj = jj +1
      If(ii<10.or.ii==Npages)write(*,'(3(a,i9))') ' Reading page= ',ii,' record =',jj,' NinPage= ',NinPage
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
          go to 10
       end If

       Call ranlux(Gg,Nrecord)
       Do ip =1,NinPage
          If(Gg(ip).lt.fraction)Then
             iCount = iCount +1
              If(iCount > iPartMax)STOP 'Attempt to read too many particles '
              
                Xpar(iCount) = (Xb(ip)-1.)*Xscale  + iFx*Vxb(ip)*factorV
                Ypar(iCount) = (Yb(ip)-1.)*Xscale  + iFy*Vyb(ip)*factorV
                Zpar(iCount) = (Zb(ip)-1.)*Xscale  + iFz*Vzb(ip)*factorV
                
              if(Xpar(iCount)>Box)Xpar(iCount)=Xpar(iCount)-Box
              if(Ypar(iCount)>Box)Ypar(iCount)=Ypar(iCount)-Box
              if(Zpar(iCount)>Box)Zpar(iCount)=Zpar(iCount)-Box
              if(Xpar(iCount)<0.)Xpar(iCount)=Xpar(iCount)+Box
              if(Ypar(iCount)<0.)Ypar(iCount)=Ypar(iCount)+Box
              if(Zpar(iCount)<0.)Zpar(iCount)=Zpar(iCount)+Box

              
                 !VX(iCount) = VXb(ip) 
                  !VY(iCount) = VYb(ip) 
                  !VZ(iCount) = VZb(ip)
             end If  ! fraction
           end Do ! ip
        end DO
        close (20)
        Nparticles = iCount
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
     DEALLOCATE(Xb,Yb,Zb,VXb,VYb,VZb)
     DEALLOCATE(Gg)
   end SUBROUTINE ReadDataPMcoords
