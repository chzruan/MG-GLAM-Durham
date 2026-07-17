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
Integer*4,  ALLOCATABLE :: Indx(:)                               ! Index of distinct halo
Real*4,     ALLOCATABLE, dimension(:) :: Rvir,aMvir,aMtot,Rvirb,aMvirb,aMtotb 
Integer*4 ::                    Nmx,Nmy,Nmz,Nbx,Nby,Nbz         ! limits for linker-list
Real*4    ::                    Cell  , t0                      ! Size in Mpch of a linker-list cell
Integer*8 ::                    Ncorr(-Nrad:Nrad)
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
!!$OMP PARALLEL DO DEFAULT(SHARED) &
!!$OMP PRIVATE (jp,i,j,k)          
      Do jp=1,Np
         i=Ceiling(Xpar(jp)/Cell)-1
         j=Ceiling(Ypar(jp)/Cell)-1
         k=Ceiling(Zpar(jp)/Cell)-1
         i=MIN(MAX(Nmx,i),Nbx)
         j=MIN(MAX(Nmy,j),Nby)
         k=MIN(MAX(Nmz,k),Nbz)
!!$OMP atomic write
         Lst(jp)      =Label(i,j,k)
!!$OMP atomic write
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
                     VXbb, VYbb, Vzbb            !    velocities
     Integer*8 :: idummy,ic,ip,iPartMax,i

     Factor   = (1.+4.*dBuffer/Box)**3 !*2 !!! remove *2
     iPartMax = Factor*Np               ! reserve extra space for total N particles
     Nparticles  = Np                   ! old number of particles
     write(*,*) ' Allocate buffers for particles. N=',iPartMax
       ALLOCATE(Xbb(iPartMax),Ybb(iPartMax),Zbb(iPartMax))
       ALLOCATE(VXbb(iPartMax),VYbb(iPartMax),VZbb(iPartMax))
       ALLOCATE(Rvirb(iPartMax),aMvirb(iPartMax),aMtotb(iPartMax))
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
       End Do
       write(*,*) '    Replicate velocities'
!$OMP PARALLEL DO DEFAULT(SHARED)  PRIVATE (ic,xx,yy,zz)
       Do ic=1,Np
          VXbb(ic)= VX(ic)
          VYbb(ic)= VY(ic)
          VZbb(ic)= VZ(ic)
       End Do
       write(*,*) '    Replicate the rest'
!$OMP PARALLEL DO DEFAULT(SHARED)  PRIVATE (ic,xx,yy,zz)
       Do ic=1,Np
          Rvirb(ic) = Rvir(ic)
          aMvirb(ic)= aMvir(ic)
          aMtotb(ic)= aMtot(ic)
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
                      Xbb(ip)   = x
                      Ybb(ip)   = y
                      Zbb(ip)   = z 
                      VXbb(ip)  = VX(ic)
                      VYbb(ip)  = VY(ic)
                      VZbb(ip)  = VZ(ic)
                      Rvirb(ip) = Rvir(ic)
                      aMvirb(ip)= aMvir(ic)
                      aMtotb(ip)= aMtot(ic)
                   end if
                  end do
              end do
           end do
        end Do
      write(*,*) ' New number of particles = ',ip
      DEALLOCATE(Xpar,Ypar,Zpar)  
      DEALLOCATE(VX,VY,VZ)
      DEALLOCATE(Rvir,aMtot,aMvir)
      Np = ip
      ALLOCATE(Xpar(Np),Ypar(Np),Zpar(Np))  
      ALLOCATE(VX(Np),VY(Np),VZ(Np))  
      ALLOCATE(Rvir(Np),aMvir(Np),aMtot(Np))  

!$OMP PARALLEL DO DEFAULT(SHARED)  PRIVATE (ip)
      Do ip =1,Np
         Xpar(ip) = Xbb(ip) ;  Ypar(ip) = Ybb(ip) ;   Zpar(ip) = Zbb(ip)
         VX(ip)   = VXbb(ip);  VY(ip)   = VYbb(ip);   VZ(ip)   = VZbb(ip)
         Rvir(ip) = Rvirb(ip); aMtot(ip)= aMtotb(ip); aMvir(ip)= aMvirb(ip)
      EndDo
      DEALLOCATE(Xbb,Ybb,Zbb)  
      DEALLOCATE(VXbb,VYbb,VZbb)
      DEALLOCATE(Rvirb,Amvirb,aMtotb)


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
SUBROUTINE DistinctCorrect
!---------------------------------------------------------------------------
use Tools
  integer*8 :: ic,ip,i
  Integer OMP_GET_THREAD_NUM, OMP_GET_MAX_THREADS

     tstart   = seconds()

  Allocate(Indx(Np))
  Indx(:) = 1
     write(*,'(2(a,i12))') ' Nparticles= ',Nparticles,' Nparticles total= ',Np
! 
!!$OMP PARALLEL DO DEFAULT(SHARED) &
!!$OMP PRIVATE (ip,x,y,z,jp,dd,k1,k2,k3,j1,j2,j3,i1,i2,i3,ind,Radius,R2max) 
     Do ip=1, Nparticles
        x   = Xpar(ip);   y = Ypar(ip);    z = Zpar(ip)
        Radius =  Rvir(ip)/1000.
        R2max  =  Radius**2
            Call Limits(x,y,z,Radius,i1,i2,j1,j2,k1,k2)
            if(mod(ip,10000)==0)write(*,'(a,i11,3x,12i4)') '  ip=',ip,i1,i2,j1,j2,k1,k2
            Do k3 =k1, k2
            Do j3 =j1, j2
            Do i3 =i1, i2
              jp =Label(i3,j3,k3)
              Do while (jp.ne.0)
                If(jp/=ip)Then
                  dd =(x-Xpar(jp))**2+(y-Ypar(jp))**2+(z-Zpar(jp))**2
                  If(dd.lt.R2max.and.aMvir(ip).lt.aMvir(jp))Then
                     Indx(ip) =0
               write(*,'(a,2es13.4,3x,3f9.3,3x,3f10.3)') ' Not clean= ',aMtot(ip),aMtot(jp),Rvir(ip),Rvir(jp), &
                                    1000.*sqrt(dd),Xpar(ip),Ypar(ip),Zpar(ip)                    
                  End If
               end If
               jp =Lst(jp)
              end do
           end do
           end do
           end Do
        End Do         ! ip
        Nclean = 0
        Do ip=1,Nparticles
           If(Indx(ip)==1)Then
              Nclean = Nclean +1
           Else
              write(*,'(a,es12.4,f9.3,3f9.3)') ' Not clean distinct halo= ',aMtot(ip),Rvir(ip), &
                                       Xpar(ip),Ypar(ip),Zpar(ip)
           end If
        end Do
        write(*,'(2(a,i8))') ' Number of clean Distinct halos= ',Nclean,' Total in the catalog= ',Nparticles
       tfinish = seconds()
      write(*,'(10x,a,T50,2f10.2)') ' time for clean distinct =',tfinish-tstart,tfinish-t0

    end SUBROUTINE DistinctCorrect
!---------------------------------------------------------------------------
!                   
!
SUBROUTINE CorrGet
!---------------------------------------------------------------------------
use Tools
  integer*8 :: ic,ip,i
  Integer OMP_GET_THREAD_NUM, OMP_GET_MAX_THREADS
  Integer*4, Allocatable :: Ncc(:,:)

     Ncorr(:) = 0
     Nthr     = OMP_GET_MAX_THREADS()
     write(*,*) ' Number of Threads =',Nthr
     Allocate(Ncc(-Nrad:Nrad,Nthr))
     Ncc(:,:) = 0
     R2max    = Rmax**2
     Radius   = Rmax
     tstart   = seconds()
     write(*,'(2(a,i12))') ' Nparticles= ',Nparticles,' Nparticles total= ',Np
! 
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (ip,x,y,z,jp,dd,k1,k2,k3,j1,j2,j3,i1,i2,i3,ind,mThr) 
     Do ip=1, Nparticles
        mThr = OMP_GET_THREAD_NUM()+1
            x   = Xpar(ip);   y = Ypar(ip);    z = Zpar(ip)
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
                     Ncc(ind,mThr) = Ncc(ind,mThr) +1
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
         End Do
      EndDo
      fact = Box**3/Nparticles**2
      write(12,'(3x,a)') ' R            dN           xi(R)        R**2xi(R)'
      Do i =-Nrad+1,Nrad
         R1 = 10.**((i)*dLogR)
         R2 = 10.**((i+1)*dLogR)
         dV = 4.188*(R2**3-R1**3)
         xi = Ncorr(i)*fact/dV -1.              ! 2Npairs/(dV/V) -1
         Rr = (R1+R2)/2.
         if(Ncorr(i)>10.and.R2<Rmax.and.Rr>0.200)write(12,'(f9.4,3x,i11,3x,2es12.4)') &
                                    Rr,Ncorr(i),xi,xi*Rr**2
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
     NpMax = 5e6
    ALLOCATE(Xpar(NpMax),Ypar(NpMax),Zpar(NpMax),VX(NpMax),VY(NpMax),VZ(NpMax))
    ALLOCATE(Rvir(NpMax),aMvir(NpMax),aMtot(NpMax))

    write(*,'(a,$)') ' Enter Catalog File Name= '
    read(*,'(a)') Line
    open(3,file=TRIM(Line))
    
    write(*,'(a,$)') ' Enter minimum Mass= '
    read(*,*) Mmin
    
    Ncenter  =0
    write(*,*) ' Mmin = '
                    ! --- read header of Catshort catalog


    read(3,'(a)')HEADER
   write(*,'(a)')HEADER
   !write(12,'(a)')HEADER

   Do i=1, 7
     Read(3,'(a)')Line 
     if(i==1.or.(i>2.and.i<14))write(12,'(a)') Line   
     write(*,'(a)') Line 
   EndDo
   Massav =0. ;  Vcmin =10000. ; Vcav =0.
   Np = 0
      Ndistinct = 0
   Do     !i =1,18000000            ! -------- read halos
     !read (3,*,iostat=iStat) X0, Y0 , Z0, VvX,VvY,VvZ,       &
     !                        aM,rr,vrmss,vcrc,vacc,aMacc,aAcc,iBDM,Conc,    &
     !                        iDist,Offset,VirRat,aLam,Rin,       &
     !                        ba,ca,xax,yax,zax
     read (3,*,iostat=iStat) X0, Y0 , Z0, VvX,VvY,VvZ,       &
                             aM,aMt,rr,vrmss,vcrc,iHalo,Conc,aNpart,    &
                             iDist,Offset,VirRat,aLam,Rin,       &
                             ba,ca,xax,yax,zax
     If(iStat/=0)exit
     !read (3,*,iostat=iStat) ii,X0, Y0 , Z0, &
     !                                          aM0,aM,vnow,vcrc
     !rr = aM0**0.33333/4.82e1      !  virial radius in kpch for Om=0.27 delta=360.8
     Np = Np +1
     If(mod(Np,1000000)==0)&
      write (*,'(a,i11,i9,3f10.4,2g12.3,3f9.2)') ' Reading: ',Np,Ncenter, X0, Y0 , Z0,aM,aMt,vcrc,Conc !,aM,rr,vnow,vcrc
45          Format(F9.4,2F9.4,3F8.1,g11.3,f8.2,2f7.1,I9,f8.1,f8.2,i4)

        If(iDist/=0)vcrc =vcrc*1.15
        If(aMt> Mmin)Then    !!! selection conditions
        Ncenter = Ncenter +1
        iHalo = Ncenter
        Massav = Massav + aMt
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
              Rvir(Ncenter) = rr
              aMvir(Ncenter) = aM
              aMtot(Ncenter) = aMt
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
  dLogR   = 0.005

  Call ReadHalos
  Call AddBuffer

  Call SizeList
  Call List
  Call DistinctCorrect
  !Call CorrGet
  
end Program CorrFunction

