!----------------------------------------------------------------------------------------
!                    Find properties of large distinct halos
!                   
!     Input:  Halo Catalog (short list)
!                
!     Output:  Distribution functions
!     Parameters:
!                     
!                      
!
MODULE   setHalos
PUBLIC      

Integer*4, Parameter :: NhaloMax = 19500000, &  ! maximum number of halos
                        NinBin   = 600,       &  ! threshold to start increasing of bin size 
                        NtabM    = 10000         ! Max size of Pk table
Real,   PARAMETER    ::                               &
                         dlog      = 0.1,             & 
                         dVLog     = 0.05,            & 
                         dMlog     = 0.10,            &
                         factorlog = 1.5,             & ! increment for variable binning
                         slog_start= 0.02,            & ! starting bin size for Vcirc/Vvir mass selection
                         vlog_start= 0.010,           & ! starting bin size for Vcirc/Vvir  velocity selection
                         dBuffer   = 2.5
Integer*8, DIMENSION(NhaloMax)        :: HaloNumber,indR
Real,      DIMENSION (NhaloMax)     :: Xc,Yc,Zc,         &
                                       VXc,VYc,Vzc,      &
                                       Amc,Rmc,Vrms,Vcirc, &
                                       Xoff,VirRatio,CAratio,aLambda,Concen, &
                                       DD                ! buffers

Real*4 ::                              deltac22,deltac,Uklow,Ukup , &
                                       xkt(0:NtabM),Pkt(0:NtabM),StepK,alog0

Integer    :: Ncenter,iFlag,iSnap,Ntab
Real*8    :: Box, Rmax, Aexpn, Om0,Omb, Ovdens, sigma8, hsmall, &
             OmLOm0
Real*8    ::  Grf,rf 
Character :: Letter*2                                         
Contains

!----------------------------------------------------

SUBROUTINE Gr(Growth,a)
!----------------------------------------------------
      !  use SetHalos
      Real*8 :: a,ww,x,Grow0,Hgrow,Growth,err
      External Hgrow
      
      x  = (OmLOm0)**0.33333     ! Normalize to growth at a=1
      Call QAG(Hgrow, 1.d-6, x, 1.d-6,1.d-6,2,ww,err,Neval,ier)
      Grow0 =sqrt(1.+x**3)/x**1.5*ww

      x  = (OmLOm0)**0.33333*a
      Call QAG(Hgrow, 1.d-6, x, 1.d-6,1.d-6,2,ww,err,Neval,ier)

      Growth =sqrt(1.+x**3)/x**1.5*ww/Grow0
    End SUBROUTINE Gr


!-------------------------------------- rms fluctuations
!                                       for mass aMass
!                              aMass is in Msun/h units
REAL*8 FUNCTION LogSlope(aMass)
!---------------------------------------
!Use SetHalos
     Real*8, PARAMETER :: pi = 3.1415926535d0, xx =0.691d0, dM =1.e-1
     REAL*8     ::        a,aMass,dd,err,dm1,ds1,M1,M2,s1,s2,s
     REAL*8, EXTERNAL ::         Ptophat

      rf= (aMass/(1.162e12*Om0))**0.33333333  ! in Mpc/h units
      Call QAG(Ptophat, 1.d-6, 1.d3, 1.d-6,1.d-6,2,dd,err,Neval,ier)
      s =GrF*sqrt((dd)/rf**3/(2.*pi**2))

      M1 = aMass*(1.+dM)     
      rf= (M1/(1.162e12*Om0))**0.33333333  ! in Mpc/h units
      Call QAG(Ptophat, 1.d-6, 1.d3, 1.d-6,1.d-6,2,dd,err,Neval,ier)
      s1 =GrF*sqrt((dd)/rf**3/(2.*pi**2))

      M2 = aMass*(1.-dM)     
      rf= (M2/(1.162e12*Om0))**0.33333333  ! in Mpc/h units
      Call QAG(Ptophat, 1.d-6, 1.d3, 1.d-6,1.d-6,2,dd,err,Neval,ier)
      s2 =GrF*sqrt((dd)/rf**3/(2.*pi**2))

      LogSlope = -aMass/s*(s1-s2)/(M1-M2)

    END FUNCTION LogSlope
!-------------------------------------- rms fluctuations
!                                       for mass aMass
!                              aMass is in Msun/h units
REAL*8 FUNCTION SigMass2(aMass)
!---------------------------------------
!Use SetHalos
     Real*8, PARAMETER :: pi = 3.1415926535d0, xx =0.691d0
     REAL*8     ::        a,aMass,dd,err
     REAL*8, EXTERNAL ::         Ptophat

      rf= (aMass/(1.162e12*Om0))**0.33333333  ! in Mpc/h units

      Call QAG(Ptophat, 1.d-6, 1.d3, 1.d-6,1.d-6,2,dd,err,Neval,ier)
      SigMass2 =GrF*sqrt((dd)/rf**3/(2.*pi**2))

    END FUNCTION SigMass2
!--------------------------------------------------- 
!                                     
!                   Setup for the power sapectrum         
!                              
Subroutine SetPk
!----------------------------------------------------        
        Real*8, parameter :: pi=3.141592653
        Real*8            :: sig8,a,wk
        REAL*8  :: Ptophat,Ppk,FF,Growth,aM
        external Ptophat

         Call ReadTable
      !   Pkt(:) = Pkt(:)*1.023   !!!! scale up Pk for bigMD sim
         a = Aexpn

         Call Gr(Growth,a)  

         GrF = Growth
         dM =1.05
         rf = 8.
         aM = 1.162e12*Om0*rf**3
         Sig8  = SigMass2(aM)
         Write(*,'(1p,3(a,g14.4))') 'a = ',a,' GrFactor =',GrF
         Write(*,'(1p,3(a,g14.4))') 'M = ',aM, ' R = ',rf,' Sigma8 = ',Sig8

        write(*,'(a,$)') ' Do you want to change Sigma8?  [New value or "return" for "no"] '
        READ  (*,'(g12.4)',IOSTAT=iStat) Sig8new
           if(iStat.eq.0.and.Sig8new.gt.0.) Then
              Pkt(:) = (Sig8new/Sig8)**2*Pkt(:)
              sig8 = SigMass2(aM)
              write(*,*) ' New sigma_8 = ', sig8,Sig8New
           Else
              write(*,*) ' No change in sigma_8'
           End if

      end Subroutine SetPk
!---------------------------------------
SUBROUTINE  ReadTable
!                     interpolate table with p(k)
!---------------------------------------
!Use SetHalos
   Character*120 :: Name,Line
   Real*8 :: a,ww,x,Grow0,Hgrow,Growth,Omc,OmL,B
   Logical :: exst
        Name = 'PkTable.dat'
        Open(1,file=TRIM(Name))
        Call ParseHeader(Omb,B)
        Call ParseHeader(Omc,B)
        Call ParseHeader(OmL,B)
        Call ParseHeader(Om0,B)
        Call ParseHeader(A,sigma8)
        write(*,'(a,4g12.4)') 'Omb =',Omb,Omc,A,sigma8
        OmLOm0 = OmL/Om0
        j =  0
        Do i=1,NtabM
           Read(1,*,iostat=ierr) ak,Pk
           If(ierr/=0)exit
           xkt(i)  = ak
           Pkt(i)  = Pk 
           j       = j+1
        EndDo
        Ntab = j
      StepK = log10 (xkt (Ntab) / xkt (1) )/(Ntab-1) 
      alog0 = log10 (xkt (1) ) 
                                                                        
                            ! test that spacing is the same             
      Do k = 2, Ntab - 1 
        ss = log10 (xkt (k + 1) / xkt (k) ) 
        If (abs (ss / StepK - 1.) .gt.2.e-2) Then 
        Write (*,*) ' error in K spacing. k=', k, xkt (k + 1),xkt (k)
        STOP 
        EndIf 
      EndDo 


   write(*,*) ' Number of lines read =',Ntab
                                                                        
      END SUBROUTINE ReadTable                      
!---------------------------------------------
!
subroutine ParseHeader(A,B)
      CHARACTER*150 :: Str1
      Real*8 :: A,B
       read(1,'(a)') Str1

       ipos = INDEX(Str1,'=')   ! First equal position
       Str1 = Str1(ipos+1:)     ! shift left
       read(Str1,*) A
       ipos = INDEX(Str1,'=')   ! Second equal position
       if(ipos==0)Then
          B = 0.
       Else
          Str1 = Str1(ipos+1:)     ! shift left
          read(Str1,*) B
       End if
end subroutine ParseHeader

!--------------------------------------------------------------
!                       Reads and recognizes files
!  
      SUBROUTINE GetFile(iFile)
!--------------------------------------------------------------
Integer*4 , parameter :: Nbin = 1000
Integer*4, SAVE  :: iFileout
Logical    :: exst

Character*120 :: file1,file2,catname,listname,HEADER*45,line,line2     
Open(20,file='ErrorAnalyzeCatalogs.dat',position='append')

If(iFile == 0) Then    ! recongize configuration

      write(*,'(a,$)')' Enter snapshot number, Catalog letter, Realization and Box size => '
      read(*,*)iSnap,Letter,iRealiz,Box
      If(iRealiz.le.0)Then
         write(file2,'(a,2a1,i4.4,a)')'HaloStaTotM',TRIM(Letter),'.',iSnap,'.DAT'
         Open( 2,file=TRIM(file2))                                 ! output
         iFileout = 0
         write(file1,'(a,2a1,i4.4,a1,i2.2,a)')'Catshort',TRIM(Letter),'.',iSnap,  &
              '.',iFileout,'.DAT'
         write(*,'(a)') TRIM(file1)
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
      Else
         iFlag = 1
         write(file2,'(a,2a1,i4.4,a1,i4.4,a)')'HaloStaTot',TRIM(Letter),'.',iSnap,'.',iRealiz,'.DAT'
         Open( 2,file=TRIM(file2))                                 ! output
         iFileout = 0
         write(file1,'(a,2a1,i4.4,a1,i4.4,a)')'Catshort',TRIM(Letter),'.',iSnap,'.',iRealiz,  &
              '.DAT'
         
      EndIf
         Inquire(file=TRIM(file1),exist = exst)
         If(exst)Then              
            open(3,file = TRIM(file1))
         Else                     
            iFile = 0
            Return
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
         
         write(file1,'(a,2a1,i4.4,a1,i2.2,a)')'Catshort',TRIM(Letter),'.',iSnap,  &
                                              '.',iFileout,'.DAT'
         Inquire(file=TRIM(file1),exist = exst)
         If(exst)Then              ! multiple files CatshortX.XXX.XX.DAT
            open(3,file = TRIM(file1))
         Else                      ! no more files to read
            iFile = 0
            Return
         End If
      End If
End If
      end SUBROUTINE GetFile
!--------------------------------------------------------------
!                         Velocity function of distinct halos
!  
      SUBROUTINE VelocityFunction
!--------------------------------------------------------------
Integer, parameter :: Nbin = 1000
Integer, Dimension(-Nbin:Nbin) ::         &
                             mHalosD, iHalosD, vHalosD
Real, Dimension(-Nbin:Nbin) ::  massbinD,vbinD


vHalosD  = 0  ; vbinD        = 0.
        nD=0            ! count all halos
        Do i=1,Ncenter
           vv = sqrt(VXc(i)**2+VYc(i)**2+VZc(i)**2)+1.e-5
                     nD =nD +1

                     ind =INT(log10(vv)/dVlog+Nbin)-Nbin
                     If(ind<Nbin)Then
                            vHalosD(ind) = vHalosD(ind) + 1
                            vbinD(ind)      = vbinD(ind)      + vv
                     EndIf
        EndDo

        vbinD = vbinD/(max(vHalosD,1))

                             ! ------------------------------- print mass function of Distinct and Subhalos
        iHalosD  = vHalosD

        istart =-Nbin      ! find first and last non-empty bin
        Do i =-Nbin,Nbin
           If(vHalosD(i)>0.)Then
              istart =i ; exit
           EndIf
        EndDo
        ifinish =Nbin
        Do i =Nbin,-Nbin,-1
           If(vHalosD(i)>0)Then
              ifinish =i ; exit
           EndIf
        EndDo
        Do i =ifinish-1,istart,-1
           vHalosD(i) = vHalosD(i) +vHalosD(i+1)
        EndDo
        write(2,'(/2a,T90,a)') 'Vdrift:V(km/s)     N(>V)     n(>V)Mpch3    V(km/s)        dN     VdN/dV'
        Do i=istart,ifinish         ! print results
              aVm =10.**(i*dVlog)
              dV    = 10.**((i+1)*dVlog)-aVm
              If(aVm>30.)Then
              write(2,'(3x,1P,G12.5,2(2g11.4,2x,1P,4g11.4))')aVm, vHalosD(i),vHalosD(i)/(Box)**3
           End If
        EndDo

      End SUBROUTINE VelocityFunction

!--------------------------------------------------------------
!                        
!  
      SUBROUTINE VelocityDistr
!--------------------------------------------------------------
Integer, parameter :: Nbin = 1000
Real*4,  Dimension(-Nbin:Nbin) ::         &
                             mHalosDx, iHalosDx, vHalosDx, &
                             mHalosDy, iHalosDy, vHalosDy, &
                             mHalosDz, iHalosDz, vHalosDz

Real*4, Dimension(-Nbin:Nbin) ::  vbinDx,vbinDy,vbinDz


vHalosDx  = 0  ; vbinDx        = 0.
vHalosDy  = 0  ; vbinDy        = 0.
vHalosDz  = 0  ; vbinDz        = 0.
         nD=0             ! count all halos
        Do i=1,Ncenter
           vvx = abs(VXc(i))+1.e-10
           vvy = abs(VYc(i))+1.e-10
           vvz = abs(VZc(i))+1.e-10 

                     nD =nD +1

                     ind =INT(log10(vvx)/dVlog+Nbin)-Nbin
                     if(ind<-1000)write(*,'(i9,3f9.3,i12)') i,VXc(i),VYc(i),VZc(i),ind
                     If(ind<Nbin)Then
                            vHalosDx(ind) = vHalosDx(ind) + 1
                            vbinDx(ind)   = vbinDx(ind)   + vvx
                     EndIf
                     ind =INT(log10(vvy)/dVlog+Nbin)-Nbin
                     If(ind<Nbin)Then
                            vHalosDy(ind) = vHalosDy(ind) + 1
                            vbinDy(ind)   = vbinDy(ind)   + vvy
                     EndIf
                     ind =INT(log10(vvz)/dVlog+Nbin)-Nbin
                     If(ind<Nbin)Then
                            vHalosDz(ind) = vHalosDz(ind) + 1
                            vbinDz(ind)   = vbinDz(ind)   + vvz
                     EndIf

        EndDo

        Do i=-Nbin,Nbin
            vbinDx(i) = vbinDx(i)/(max(vHalosDx(i),1.))
            vbinDy(i) = vbinDy(i)/(max(vHalosDy(i),1.))
            vbinDz(i) = vbinDz(i)/(max(vHalosDz(i),1.))
         End Do

                             ! ------------------------------- print 

        write(2,'(/a,5x,2a)') '     V(km/s)    Vx       dN            Vy          dN          Vz        dN'
        Do i=0,100         ! print results
              aVm =10.**(i*dVlog)
              dV    = 10.**((i+1)*dVlog)-aVm
              
              If(aVm<5000.)Then
              write(2,'(3x,f8.2,1p,2(6g12.4,3x))')aVm, vbinDx(i),(vHalosDx(i)/Ncenter/dV), &
                       vbinDy(i),(vHalosDy(i)/Ncenter/dV),vbinDz(i),(vHalosDz(i)/Ncenter/dV)    
              End If
        EndDo

      End SUBROUTINE VelocityDistr


!--------------------------------------------------------------
!                         Mass function of distinct halos
!  
      SUBROUTINE MassFunction
!--------------------------------------------------------------
Integer, parameter :: Nbin = 1000
Integer, Dimension(-Nbin:Nbin) ::         &
                             mHalosD, iHalosD, vHalosD
Real*8  :: aMm,sigma,aMd,sigmad,Slope
Real, Dimension(-Nbin:Nbin) ::  massbinD,vbinD

mHalosD  = 0 ; massbinD = 0.
vHalosD  = 0  ; vbinD        = 0.
         nD=0            ! count all halos
        Do i=1,Ncenter
                     nD =nD +1
                     ind =INT(log10(Amc(i))/dMlog+Nbin)-Nbin
                     If(ind<Nbin)Then
                            mHalosD(ind)  = mHalosD(ind)   + 1
                            massbinD(ind) = massbinD(ind) + Amc(i)
                     EndIf
                     ind =INT(log10(Vcirc(i))/dVlog+Nbin)-Nbin
                     If(ind<Nbin)Then
                            vHalosD(ind) = vHalosD(ind) + 1
                            vbinD(ind)      = vbinD(ind)      + Vcirc(i)
                     EndIf
        EndDo
        massbinD = massbinD/(max(mHalosD,1))
        vbinD = vbinD/(max(vHalosD,1))

                 ! ------------------------------- print mass function of Distinct
        istart =-Nbin      ! find first and last non-empty bin
        Do i =-Nbin,Nbin
           If(mHalosD(i)>0.)Then
              istart =i ; exit
           EndIf
        EndDo
        ifinish =Nbin
        Do i =Nbin,-Nbin,-1
           If(mHalosD(i)>0)Then
              ifinish =i ; exit
           EndIf
        EndDo
        iHalosD  = mHalosD
        Do i =ifinish-1,istart,-1
           mHalosD(i) = mHalosD(i) +mHalosD(i+1)
        EndDo

        write(2,'(/3a)') ' Distinct Mass(Msunh)   sigma      N(>M)    n(>M)/Mpch3',&
                          '  Mass       sigma      Slope          dN       MdN/dM'
        Do i=istart,ifinish         ! print results
              aMm =10.**(i*dMlog)
              dM    = 10.**((i+1)*dMlog)-aMm
              if(massbinD(i)<1.e-7)massbinD(i)=aMm
              sigma = SigMass2(aMm)
              aMd   =   massbinD(i)
              sigmaD=  SigMass2(aMd) 
              Slope =  LogSlope(aMd) 
              write(2,'(9x,es12.4,f9.4,i11,2es12.4,2f9.4,3x,i11,3es12.4)')  &
                   aMm,         sigma, mHalosD(i),mHalosD(i)/(Box)**3, &
                   massbinD(i),sigmaD,Slope,iHalosD(i), massbinD(i)*iHalosD(i)/dM/(Box)**3
        EndDo

                             ! ------------------------------- print mass function of Distinct and Subhalos
        iHalosD  = vHalosD

        istart =-Nbin      ! find first and last non-empty bin
        Do i =-Nbin,Nbin
           If(vHalosD(i)>0.)Then
              istart =i ; exit
           EndIf
        EndDo
        ifinish =Nbin
        Do i =Nbin,-Nbin,-1
           If(vHalosD(i)>0)Then
              ifinish =i ; exit
           EndIf
        EndDo
        Do i =ifinish-1,istart,-1
           vHalosD(i) = vHalosD(i) +vHalosD(i+1)
        EndDo
        write(2,'(/2a,T90,a)') 'Dist:V(km/s)     N(>V)     n(>V)Mpch3    V(km/s)        dN     VdN/dV'
        Do i=istart,ifinish         ! print results
              aVm =10.**(i*dVlog)
              dV    = 10.**((i+1)*dVlog)-aVm
              write(2,'(3x,1P,G12.5,2(2g11.4,2x,1P,4g11.4))')aVm, vHalosD(i),vHalosD(i)/(Box)**3,   &
                          vbinD(i),iHalosD(i),vbinD(i)*iHalosD(i)/dV/(Box)**3
        EndDo
        write(*,*) ' exit Mass Function'
      End SUBROUTINE  MassFunction

!--------------------------------------------------------------
!                     Concentration and offset functions of distinct halos
!  
      SUBROUTINE ConFunction
!--------------------------------------------------------------
Integer, parameter :: Nbin = 1000
Integer, Dimension(-Nbin:Nbin) ::         &
                             mHalosD 
Real*8  :: aMm
Real, Dimension(-Nbin:Nbin) ::  massbinD,cbin,c2bin,offbin,off2bin

mHalosD  = 0  ; massbinD = 0.
cbin     = 0  ; c2bin    =0.; offbin  = 0. ; off2bin =0.
         nD=0            ! count all halos
        Do i=1,Ncenter
                     nD =nD +1
                     ind =INT(log10(Amc(i))/dMlog+Nbin)-Nbin
                     If(ind<Nbin)Then
                            mHalosD(ind)  = mHalosD(ind)   + 1
                            massbinD(ind) = massbinD(ind)  + Amc(i)
                            cbin(ind)     = cbin(ind)      + Concen(i)
                            c2bin(ind)    = c2bin(ind)     + Concen(i)**2
                            offbin(ind)   = offbin(ind)    + Xoff(i)
                            off2bin(ind)  = off2bin(ind)   + Xoff(i)**2
                     EndIf
        EndDo
        massbinD = massbinD/(max(mHalosD,1))
        cbin     = cbin/max(mHalosD,1)
        c2bin    = sqrt(c2bin/max(mHalosD,1)-cbin**2 )
        offbin   = offbin/max(mHalosD,1)
        off2bin  = sqrt(off2bin/max(mHalosD,1)-offbin**2 )

                 ! ------------------------------- print mass function of Distinct
        istart =-Nbin      ! find first and last non-empty bin
        Do i =-Nbin,Nbin
           If(mHalosD(i)>0.)Then
              istart =i ; exit
           EndIf
        EndDo
        ifinish =Nbin
        Do i =Nbin,-Nbin,-1
           If(mHalosD(i)>0)Then
              ifinish =i ; exit
           EndIf
        EndDo


        write(2,'(/3a)') ' Distinct Mass(Msunh)     Nhalos   ',&
                          '  <C>       rms C       <Xoff>      rms Xoff'
        Do i=istart,ifinish         ! print results
              aMm =10.**(i*dMlog)
              dM    = 10.**((i+1)*dMlog)-aMm
              if(massbinD(i)<1.e-7)massbinD(i)=aMm
              write(2,'(9x,es12.4,i11,2f9.4,3x,2f9.4)')  &
                   aMm, mHalosD(i), cbin(i),c2bin(i),offbin(i),off2bin(i)
        EndDo

        write(*,*) ' exit Concentration Function'
      End SUBROUTINE ConFunction



!---------------------------------------------------------
!                         Read halos catalog
!                       
      SUBROUTINE ReadHalos
!--------------------------------------------------------------
     Character     :: file1*80,Line*120,txt*80
     CHARACTER*45  :: HEADER
     Logical       :: inside


   Ncenter  =0
   Rdmax    =0.
   factor   =1.162e12*Om0*Ovdens 
                            ! --- read header of Catshort catalog
   read(3,'(a)')HEADER
   write(*,'(a)')HEADER
   write(2,'(a)')HEADER

  ! read(3,'(2(a7,f8.5))')  txt(1:7),AEXPN,Line(1:7),Step 
  ! write(*,'(2(a7,f8.5))') txt(1:7),AEXPN,Line(1:7),Step
  ! write(2,'(2(a7,f8.5))') txt(1:7),AEXPN,Line(1:7),Step
  ! read(3,'(a7,i5,a)')     txt(1:7),ISTEP,Line
  ! write(*,'(a7,i5,a)')    txt(1:7),ISTEP,Trim(Line)
  ! write(2,'(a7,i5,a)')    txt(1:7),ISTEP,Trim(Line)
  ! Do i=1, 14
  !   Read(3,'(a)')Line 
  !   if(i .le.13)write(2,'(a)') Line   
  !   write(*,'(a)') Line 
  ! EndDo

   Do i=1, 7
     Read(3,'(a)')Line 
     if(i .le.13)write(2,'(a)') Line   
     write(*,'(a)') Line 
   EndDo
   Vclimit = 0.
   iFile   = 1 


      Ndistinct = 0
   Do  !  i =1,1800000            ! -------- read halos
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
     If(mod(Ncenter,100000)==0)&
      write (*,'(a,i9,3f10.4,2g12.3,3f9.2)') ' Reading: ',Ncenter, X0, Y0 , Z0,aM,aMtot,vcrc,Conc !,aM,rr,vnow,vcrc
45          Format(F9.4,2F9.4,3F8.1,g11.3,f8.2,2f7.1,I9,f8.1,f8.2,i4)

        If(iStat.ne.0)Call GetFile(iFile)
        If(iFile ==0)Exit
        !inside = (X0>dBuffer.and.X0<Box-dBuffer)
        !inside = (inside.and.(Y0>dBuffer.and.Y0<Box-dBuffer))
        !inside = (inside.and.(Z0>dBuffer.and.Z0<Box-dBuffer))
               rvr  = (aMtot/factor)**0.3333*1000.
               !Vvir = 2.081e-3*sqrt(aMtot/rvr/Aexpn) 
!       If(vcrc> vcLimit .and. rr > 0.95*rvr)Then    !!! selection conditions
       If(vcrc> vcLimit)Then    !!! selection conditions
        Ncenter = Ncenter +1
        iHalo = Ncenter
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
        
        If(Ncenter .le. NhaloMax)Then
              Xc(Ncenter) = X0
              Yc(Ncenter) = Y0
              Zc(Ncenter) = Z0
              VXc(Ncenter) = VvX
              VYc(Ncenter) = VvY
              VZc(Ncenter) = VvZ
              Vcirc(Ncenter) = vcrc
               rvr  = (aMtot/factor)**0.3333*1000.
               !Vvir = 2.081e-3*sqrt(aMtot/rvr/Aexpn)

               !if(aMtot>3.5e13)write(*,'(1p,2g12.4,3x,8g12.4)') aM,aMtot,rr,rvr,rvr/rr,Offset,vcrc,Vvir
              Rmc(Ncenter) = rr  !rvr
              Amc(Ncenter) =  aMtot ! aM !
              HaloNumber(Ncenter) = iHalo
              Xoff(Ncenter) = Offset
              Concen(Ncenter) = Conc
                !----- Correct concentration
              r0 = rr/Conc
              d0 = min(280./8./r0,3.)
              Concen(Ncenter) = Conc *(1+d0)/exp((d0/1.5)**2)
              !CAratio(Ncenter) = ca
              VirRatio(Ncenter) = VirRat
              aLambda(Ncenter) = aLam
              !Vrms(Ncenter) = vrmss
              Rdmax  = MAX(rr/1000.,Rdmax)
              Ndistinct = Ndistinct +1
        EndIf
        EndIf     !  test mass
   Enddo
40          Format(g11.4,I9,g11.3,g10.3,2f7.1,g10.3,i8,i7,i5,g11.3)
   close (3)
    write(*,'(2(a,T39,a,i10/),a,T39,a,1p,g12.4,a)') ' Read Halos','= ',Ncenter,' Ndistinct','= ',Ndistinct, &
                                            ' Halo max Radius','= ',Rdmax,'Mpch'
    write(2,'(2(a,T39,a,i10/),a,T39,a,1p,g12.4,a)') ' Number of Halos','= ',Ncenter,' Ndistinct','= ',Ndistinct, &
                                            ' Halo max Radius','= ',Rdmax,'Mpch'
    If(Ncenter >= NhaloMax)Stop ' Not enough space for halos. Increase NhaloMax.'
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

!     ------------------------
      real function Concentration(Vratio)
!     ------------------------
!
        Real*8 :: A,Aa,C,dC,x,f

        If(Vratio<1.)Then
           Concentration = 1.
        Else
           C= 2.162
           dC =3.d-5
           A  = Vratio**2
           Aa = 1.d0
           Do while(Aa<A)
              C = C +dC
              Fc = log(1.d0+C) -C/(1.d0+C)
              Aa = 0.216216*C/Fc
           end Do
           Concentration  = C
           !write(*,'(3x,1p,5g12.4)') Vratio,C,dC,Aa,A
        End If
      End function Concentration
!     ------------------------
      real function seconds ()
!     ------------------------
!
!     purpose: returns elapsed time in seconds
      Integer, SAVE :: first=0,rate=0
      Real   , SAVE :: t0=0.
      !real*8 dummy
      !real tarray(2)
      If(first==0)Then
         CALL SYSTEM_CLOCK(i,rate)
         first =1
         t0    =float(i)/float(rate)
         seconds = 0.
      Else
!        seconds = etime(tarray)
         CALL SYSTEM_CLOCK(i)
         seconds = float(i)/float(rate)-t0
      EndIf
      return
      end function seconds
!--------------------------------------------------------------
!                        virial overdensity for cosmological model
!                        at different expansion parameter AEXPN
      Function OverdenVir()
!--------------------------------------------------------------
      xx =-(1.-Om0)*AEXPN**3/(Om0+(1.-Om0)*AEXPN**3)
      OverdenVir =(178.+82.*xx-39.*xx**2)/(1.+xx)
      write (13,*)  '      Overdensity Delta =',OverdenVir

    END Function OverdenVir

end module SetHalos
!--------------------------------------------------------------------------------------
!
!     
PROGRAM HaloSatellites
use   SetHalos 
      Character     :: file1*80,Line*120,txt*80

      CALL GetFile(0)              ! set configuration. Open first file
      !CALL ReadROCK(1)
      CALL SetPk
      iVirial = 0
      write(*,*) ' Aexpn =',AEXPN
      write(*,*) ' Om0   =',Om0
      If(TRIM(Letter) == 'V')iVirial = 1
           If(iVirial ==1)Then
               Ovdens   = OverdenVir()               ! overdensity relative to Om0\rho_crit
           Else
               Ovdens   = 200./Om0*(Om0+(1.-Om0)*AEXPN**3)
           EndIf
           write(*,'(a,i2,a,f8.3)') ' Virial =',iVirial,' Overdens =',Ovdens
      CALL ReadHalos                  ! read halos
      write(*,*) '    Go to Mass Function'     
      Call MassFunction                    
      Call ConFunction                    
      Call VelocityFunction
      !Call VcircMvir                           ! Vcirc-Mvir relation for Distinct and Subs
      !Call VmaxVvir
      !Call VmaxVvirM
      Call VelocityDistr

End PROGRAM HaloSatellites

Subroutine R_mrgref (XVALT, IRNGT,N)
!   Ranks array XVALT into index array IRNGT, using merge-sort
! __________________________________________________________
!   This version is not optimized for performance, and is thus
!   not as difficult to read as some other ones.
!   Michel Olagnon - April 2000
! __________________________________________________________
! _________________________________________________________
      Real, Dimension (N), Intent (In) :: XVALT
      Integer, Dimension (N), Intent (Out) :: IRNGT
! __________________________________________________________
!
      Integer, Dimension (:), Allocatable :: JWRKT
      Integer :: LMTNA, LMTNC
      Integer :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
!
      NVAL = Min (SIZE(XVALT), SIZE(IRNGT))
      If (NVAL <= 0) Then
         Return
      End If
!
!  Fill-in the index array, creating ordered couples
!
      Do IIND = 2, NVAL, 2
         If (XVALT(IIND-1) < XVALT(IIND)) Then
            IRNGT (IIND-1) = IIND - 1
            IRNGT (IIND) = IIND
         Else
            IRNGT (IIND-1) = IIND
            IRNGT (IIND) = IIND - 1
         End If
      End Do
      If (Modulo (NVAL, 2) /= 0) Then
         IRNGT (NVAL) = NVAL
      End If
!
!  We will now have ordered subsets A - B - A - B - ...
!  and merge A and B couples into     C   -   C   - ...
!
      Allocate (JWRKT(1:NVAL))
      LMTNC = 2
      LMTNA = 2
!
!  Iteration. Each time, the length of the ordered subsets
!  is doubled.
!
      Do
         If (LMTNA >= NVAL) Exit
         IWRKF = 0
         LMTNC = 2 * LMTNC
         IWRK = 0
!
!   Loop on merges of A and B into C
!
         Do
            IINDA = IWRKF
            IWRKD = IWRKF + 1
            IWRKF = IINDA + LMTNC
            JINDA = IINDA + LMTNA
            If (IWRKF >= NVAL) Then
               If (JINDA >= NVAL) Exit
               IWRKF = NVAL
            End If
            IINDB = JINDA
!
!   Shortcut for the case when the max of A is smaller
!   than the min of B (no need to do anything)
!
            If (XVALT(IRNGT(JINDA)) <= XVALT(IRNGT(JINDA+1))) Then
               IWRK = IWRKF
               Cycle
            End If
!
!  One steps in the C subset, that we create in the final rank array
!
            Do
               If (IWRK >= IWRKF) Then
!
!  Make a copy of the rank array for next iteration
!
                  IRNGT (IWRKD:IWRKF) = JWRKT (IWRKD:IWRKF)
                  Exit
               End If
!
               IWRK = IWRK + 1
!
!  We still have unprocessed values in both A and B
!
               If (IINDA < JINDA) Then
                  If (IINDB < IWRKF) Then
                     If (XVALT(IRNGT(IINDA+1)) > XVALT(IRNGT(IINDB+1))) &
                    & Then
                        IINDB = IINDB + 1
                        JWRKT (IWRK) = IRNGT (IINDB)
                     Else
                        IINDA = IINDA + 1
                        JWRKT (IWRK) = IRNGT (IINDA)
                     End If
                  Else
!
!  Only A still with unprocessed values
!
                     IINDA = IINDA + 1
                     JWRKT (IWRK) = IRNGT (IINDA)
                  End If
               Else
!
!  Only B still with unprocessed values
!
                  IRNGT (IWRKD:IINDB) = JWRKT (IWRKD:IINDB)
                  IWRK = IWRKF
                  Exit
               End If
!
            End Do
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 2 * LMTNA
      End Do
!
!  Clean up
!
      Deallocate (JWRKT)
      Return
!
End Subroutine R_mrgref 

!---------------------------------------
      REAL*8 FUNCTION Ppk(x)
!                     interpolate table with p(k)
!                       x is in real 1/Mpc
!---------------------------------------
        use SetHalos
      Real*8 :: x,dk

      If(x.ge.xkt(Ntab))Then  ! slope is ns =-3
         Ppk =Pkt(Ntab)/(x/xkt(Ntab))**3
         Return
      EndIf
      If(x.lt.xkt(1))Then                ! slope is ns=1
         Ppk =Pkt(1)*(x/xkt(1))
         Return
      EndIf
      ind = INT((log10(x)-alog0)/StepK) +1 
      dk  = xkt(ind+1)-xkt(ind)
        Ppk   = (Pkt(ind)*(xkt(ind+1)-x)+Pkt(ind+1)*(x-xkt(ind)))/dk
    End FUNCTION Ppk
subroutine qag ( f, a, b, epsabs, epsrel, key, result, abserr, neval, ier )

!*****************************************************************************80
!
!! QAG approximates an integral over a finite interval.
!
!  Discussion:
!
!    The routine calculates an approximation RESULT to a definite integral   
!      I = integral of F over (A,B),
!    hopefully satisfying
!      || I - RESULT || <= max ( EPSABS, EPSREL * ||I|| ).
!
!    QAG is a simple globally adaptive integrator using the strategy of 
!    Aind (Piessens, 1973).  It is possible to choose between 6 pairs of
!    Gauss-Kronrod quadrature formulae for the rule evaluation component. 
!    The pairs of high degree of precision are suitable for handling
!    integration difficulties due to a strongly oscillating integrand.
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, external real F, the name of the function routine, of the form
!      function f ( x )
!      real f
!      real x
!    which evaluates the integrand function.
!
!    Input, real A, B, the limits of integration.
!
!    Input, real EPSABS, EPSREL, the absolute and relative accuracy requested.
!
!    Input, integer KEY, chooses the order of the local integration rule:
!    1,  7 Gauss points, 15 Gauss-Kronrod points,
!    2, 10 Gauss points, 21 Gauss-Kronrod points,
!    3, 15 Gauss points, 31 Gauss-Kronrod points,
!    4, 20 Gauss points, 41 Gauss-Kronrod points,
!    5, 25 Gauss points, 51 Gauss-Kronrod points,
!    6, 30 Gauss points, 61 Gauss-Kronrod points.
!
!    Output, real RESULT, the estimated value of the integral.
!
!    Output, real ABSERR, an estimate of || I - RESULT ||.
!
!    Output, integer NEVAL, the number of times the integral was evaluated.
!
!    Output, integer IER, return code.
!    0, normal and reliable termination of the routine.  It is assumed that the 
!      requested accuracy has been achieved.
!    1, maximum number of subdivisions allowed has been achieved.  One can 
!      allow more subdivisions by increasing the value of LIMIT in QAG. 
!      However, if this yields no improvement it is advised to analyze the
!      integrand to determine the integration difficulties.  If the position
!      of a local difficulty can be determined, such as a singularity or
!      discontinuity within the interval) one will probably gain from 
!      splitting up the interval at this point and calling the integrator 
!      on the subranges.  If possible, an appropriate special-purpose 
!      integrator should be used which is designed for handling the type 
!      of difficulty involved.
!    2, the occurrence of roundoff error is detected, which prevents the
!      requested tolerance from being achieved.
!    3, extremely bad integrand behavior occurs at some points of the
!      integration interval.
!    6, the input is invalid, because EPSABS < 0 and EPSREL < 0.
!
!  Local parameters:
!
!    LIMIT is the maximum number of subintervals allowed in
!    the subdivision process of QAGE.
!
  implicit none

  integer, parameter :: limit = 500

  real*8 a
  real*8 abserr
  real*8 alist(limit)
  real*8 b
  real*8 blist(limit)
  real*8 elist(limit)
  real*8 epsabs
  real*8 epsrel
  real*8, external :: f
  integer ier
  integer iord(limit)
  integer key
  integer last
  integer neval
  real*8 result
  real*8 rlist(limit)
  call qage ( f, a, b, epsabs, epsrel, key, limit, result, abserr, neval, &
    ier, alist, blist, rlist, elist, iord, last )

end subroutine qag

!****************************************************************************
subroutine qage ( f, a, b, epsabs, epsrel, key, limit, result, abserr, neval, &
  ier, alist, blist, rlist, elist, iord, last )

!*****************************************************************************80
!
!! QAGE estimates a definite integral.
!
!  Discussion:
!
!    The routine calculates an approximation RESULT to a definite integral   
!      I = integral of F over (A,B),
!    hopefully satisfying
!      || I - RESULT || <= max ( EPSABS, EPSREL * ||I|| ).
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, external real F, the name of the function routine, of the form
!      function f ( x )
!      real f
!      real x
!    which evaluates the integrand function.
!
!    Input, real A, B, the limits of integration.
!
!    Input, real EPSABS, EPSREL, the absolute and relative accuracy requested.
!
!    Input, integer KEY, chooses the order of the local integration rule:
!    1,  7 Gauss points, 15 Gauss-Kronrod points,
!    2, 10 Gauss points, 21 Gauss-Kronrod points,
!    3, 15 Gauss points, 31 Gauss-Kronrod points,
!    4, 20 Gauss points, 41 Gauss-Kronrod points,
!    5, 25 Gauss points, 51 Gauss-Kronrod points,
!    6, 30 Gauss points, 61 Gauss-Kronrod points.
!
!    Input, integer LIMIT, the maximum number of subintervals that
!    can be used.
!
!    Output, real RESULT, the estimated value of the integral.
!
!    Output, real ABSERR, an estimate of || I - RESULT ||.
!
!    Output, integer NEVAL, the number of times the integral was evaluated.
!
!    Output, integer IER, return code.
!    0, normal and reliable termination of the routine.  It is assumed that the 
!      requested accuracy has been achieved.
!    1, maximum number of subdivisions allowed has been achieved.  One can 
!      allow more subdivisions by increasing the value of LIMIT in QAG. 
!      However, if this yields no improvement it is advised to analyze the
!      integrand to determine the integration difficulties.  If the position
!      of a local difficulty can be determined, such as a singularity or
!      discontinuity within the interval) one will probably gain from 
!      splitting up the interval at this point and calling the integrator 
!      on the subranges.  If possible, an appropriate special-purpose 
!      integrator should be used which is designed for handling the type 
!      of difficulty involved.
!    2, the occurrence of roundoff error is detected, which prevents the
!      requested tolerance from being achieved.
!    3, extremely bad integrand behavior occurs at some points of the
!      integration interval.
!    6, the input is invalid, because EPSABS < 0 and EPSREL < 0.
!
!    Workspace, real ALIST(LIMIT), BLIST(LIMIT), contains in entries 1 
!    through LAST the left and right ends of the partition subintervals.
!
!    Workspace, real RLIST(LIMIT), contains in entries 1 through LAST
!    the integral approximations on the subintervals.
!
!    Workspace, real ELIST(LIMIT), contains in entries 1 through LAST
!    the absolute error estimates on the subintervals.
!
!    Output, integer IORD(LIMIT), the first K elements of which are pointers 
!    to the error estimates over the subintervals, such that
!    elist(iord(1)), ..., elist(iord(k)) form a decreasing sequence, with
!    k = last if last <= (limit/2+2), and k = limit+1-last otherwise.
!
!    Output, integer LAST, the number of subintervals actually produced 
!    in the subdivision process.
!
!  Local parameters:
!
!    alist     - list of left end points of all subintervals
!                       considered up to now
!    blist     - list of right end points of all subintervals
!                       considered up to now
!    elist(i)  - error estimate applying to rlist(i)
!    maxerr    - pointer to the interval with largest error estimate
!    errmax    - elist(maxerr)
!    area      - sum of the integrals over the subintervals
!    errsum    - sum of the errors over the subintervals
!    errbnd    - requested accuracy max(epsabs,epsrel*abs(result))
!    *****1    - variable for the left subinterval
!    *****2    - variable for the right subinterval
!    last      - index for subdivision
!
  implicit none

  integer limit

  real*8 a
  real*8 abserr
  real*8 alist(limit)
  real*8 area
  real*8 area1
  real*8 area12
  real*8 area2
  real*8 a1
  real*8 a2
  real*8 b
  real*8 blist(limit)
  real*8 b1
  real*8 b2
  real*8 c
  real*8 defabs
  real*8 defab1
  real*8 defab2
  real*8 elist(limit)
  real*8 epsabs
  real*8 epsrel
  real*8 errbnd
  real*8 errmax
  real*8 error1
  real*8 error2
  real*8 erro12
  real*8 errsum
  real*8, external :: f
  integer ier
  integer iord(limit)
  integer iroff1
  integer iroff2
  integer key
  integer keyf
  integer last
  integer maxerr
  integer neval
  integer nrmax
  real*8 resabs
  real*8 result
  real*8 rlist(limit)
!
!  Test on validity of parameters.
!
  ier = 0
  neval = 0
  last = 0
  result = 0.0e+00
  abserr = 0.0e+00
  alist(1) = a
  blist(1) = b
  rlist(1) = 0.0e+00
  elist(1) = 0.0e+00
  iord(1) = 0

  if ( epsabs < 0.0e+00 .and. epsrel < 0.0e+00 ) then
    ier = 6
    return
  end if
!
!  First approximation to the integral.
!
  keyf = key
  keyf = max ( keyf, 1 )
  keyf = min ( keyf, 2 )

  c = keyf
  neval = 0

  if ( keyf == 1 ) then
    call qk15 ( f, a, b, result, abserr, defabs, resabs )
  else if ( keyf == 2 ) then
    call qk21 ( f, a, b, result, abserr, defabs, resabs )
  end if

  last = 1
  rlist(1) = result
  elist(1) = abserr
  iord(1) = 1
!
!  Test on accuracy.
!
  errbnd = max ( epsabs, epsrel * abs ( result ) )

  if ( abserr <= 5.0e+01 * epsilon ( defabs ) * defabs .and. &
    errbnd < abserr ) then
    ier = 2
  end if

  if ( limit == 1 ) then
    ier = 1
  end if

  if ( ier /= 0 .or. &
    ( abserr <= errbnd .and. abserr /= resabs ) .or. &
    abserr == 0.0e+00 ) then

    if ( keyf /= 1 ) then
      neval = (10*keyf+1) * (2*neval+1)
    else
      neval = 30 * neval + 15
    end if

    return

  end if
!
!  Initialization.
!
  errmax = abserr
  maxerr = 1
  area = result
  errsum = abserr
  nrmax = 1
  iroff1 = 0
  iroff2 = 0

  do last = 2, limit
!
!  Bisect the subinterval with the largest error estimate.
!
    a1 = alist(maxerr)
    b1 = 0.5E+00 * ( alist(maxerr) + blist(maxerr) )
    a2 = b1
    b2 = blist(maxerr)

    if ( keyf == 1 ) then
      call qk15 ( f, a1, b1, area1, error1, resabs, defab1 )
    else if ( keyf == 2 ) then
      call qk21 ( f, a1, b1, area1, error1, resabs, defab1 )
    end if

    if ( keyf == 1 ) then
      call qk15 ( f, a2, b2, area2, error2, resabs, defab2 )
    else if ( keyf == 2 ) then
      call qk21 ( f, a2, b2, area2, error2, resabs, defab2 )

    end if
!
!  Improve previous approximations to integral and error and
!  test for accuracy.
!
    neval = neval + 1
    area12 = area1 + area2
    erro12 = error1 + error2
    errsum = errsum + erro12 - errmax
    area = area + area12 - rlist(maxerr)

    if ( defab1 /= error1 .and. defab2 /= error2 ) then

      if ( abs ( rlist(maxerr) - area12 ) <= 1.0e-05 * abs ( area12 ) &
        .and. 9.9e-01 * errmax <= erro12 ) then
        iroff1 = iroff1 + 1
      end if

      if ( 10 < last .and. errmax < erro12 ) then
        iroff2 = iroff2 + 1
      end if

    end if

    rlist(maxerr) = area1
    rlist(last) = area2
    errbnd = max ( epsabs, epsrel * abs ( area ) )
!
!  Test for roundoff error and eventually set error flag.
!
    if ( errbnd < errsum ) then

      if ( 6 <= iroff1 .or. 20 <= iroff2 ) then
        ier = 2
      end if
!
!  Set error flag in the case that the number of subintervals
!  equals limit.
!
      if ( last == limit ) then
        ier = 1
      end if
!
!  Set error flag in the case of bad integrand behavior
!  at a point of the integration range.
!
      if ( max ( abs ( a1 ), abs ( b2 ) ) <= ( 1.0e+00 + c * 1.0e+03 * &
        epsilon ( a1 ) ) * ( abs ( a2 ) + 1.0e+04 * tiny ( a2 ) ) ) then
        ier = 3
      end if

    end if
!
!  Append the newly-created intervals to the list.
!
    if ( error2 <= error1 ) then
      alist(last) = a2
      blist(maxerr) = b1
      blist(last) = b2
      elist(maxerr) = error1
      elist(last) = error2
    else
      alist(maxerr) = a2
      alist(last) = a1
      blist(last) = b1
      rlist(maxerr) = area2
      rlist(last) = area1
      elist(maxerr) = error2
      elist(last) = error1
    end if
!
!  Call QSORT to maintain the descending ordering
!  in the list of error estimates and select the subinterval
!  with the largest error estimate (to be bisected next).
!
    call qsort ( limit, last, maxerr, errmax, elist, iord, nrmax )
 
    if ( ier /= 0 .or. errsum <= errbnd ) then
      exit
    end if

  end do
!
!  Compute final result.
!
  result = sum ( rlist(1:last) )

  abserr = errsum

  if ( keyf /= 1 ) then
    neval = ( 10 * keyf + 1 ) * ( 2 * neval + 1 )
  else
    neval = 30 * neval + 15
  end if

end subroutine qage
subroutine qk15 ( f, a, b, result, abserr, resabs, resasc )

!*****************************************************************************80
!
!! QK15 carries out a 15 point Gauss-Kronrod quadrature rule.
!
!  Discussion:
!
!    This routine approximates
!      I = integral ( A <= X <= B ) F(X) dx
!    with an error estimate, and
!      J = integral ( A <= X <= B ) | F(X) | dx
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, external real*8 F, the name of the function routine, of the form
!      function f ( x )
!      real*8 f
!      real*8 x
!    which evaluates the integrand function.
!
!    Input, real*8 A, B, the limits of integration.
!
!    Output, real*8 RESULT, the estimated value of the integral.
!    RESULT is computed by applying the 15-point Kronrod rule (RESK) 
!    obtained by optimal addition of abscissae to the 7-point Gauss rule 
!    (RESG).
!
!    Output, real*8 ABSERR, an estimate of | I - RESULT |.
!
!    Output, real*8 RESABS, approximation to the integral of the absolute
!    value of F.
!
!    Output, real*8 RESASC, approximation to the integral | F-I/(B-A) | 
!    over [A,B].
!
!  Local Parameters:
!
!           the abscissae and weights are given for the interval (-1,1).
!           because of symmetry only the positive abscissae and their
!           corresponding weights are given.
!
!           xgk    - abscissae of the 15-point Kronrod rule
!                    xgk(2), xgk(4), ...  abscissae of the 7-point
!                    Gauss rule
!                    xgk(1), xgk(3), ...  abscissae which are optimally
!                    added to the 7-point Gauss rule
!
!           wgk    - weights of the 15-point Kronrod rule
!
!           wg     - weights of the 7-point Gauss rule
!
!           centr  - mid point of the interval
!           hlgth  - half-length of the interval
!           absc   - abscissa
!           fval*  - function value
!           resg   - result of the 7-point Gauss formula
!           resk   - result of the 15-point Kronrod formula
!           reskh  - approximation to the mean value of f over (a,b),
!                    i.e. to i/(b-a)
!
  implicit none

  real*8 a
  real*8 absc
  real*8 abserr
  real*8 b
  real*8 centr
  real*8 dhlgth
  real*8, external :: f
  real*8 fc
  real*8 fsum
  real*8 fval1
  real*8 fval2
  real*8 fv1(7)
  real*8 fv2(7)
  real*8 hlgth
  integer j
  integer jtw
  integer jtwm1
  real*8 resabs
  real*8 resasc
  real*8 resg
  real*8 resk
  real*8 reskh
  real*8 result
  real*8 wg(4)
  real*8 wgk(8)
  real*8 xgk(8)

  data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8)/ &
       9.914553711208126e-01,   9.491079123427585e-01, &
       8.648644233597691e-01,   7.415311855993944e-01, &
       5.860872354676911e-01,   4.058451513773972e-01, &
       2.077849550078985e-01,   0.0e+00              /
  data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8)/ &
       2.293532201052922e-02,   6.309209262997855e-02, &
       1.047900103222502e-01,   1.406532597155259e-01, &
       1.690047266392679e-01,   1.903505780647854e-01, &
       2.044329400752989e-01,   2.094821410847278e-01/
  data wg(1),wg(2),wg(3),wg(4)/ &
       1.294849661688697e-01,   2.797053914892767e-01, &
       3.818300505051189e-01,   4.179591836734694e-01/
!
  centr = 5.0e-01*(a+b)
  hlgth = 5.0e-01*(b-a)
  dhlgth = abs(hlgth)
!
!  Compute the 15-point Kronrod approximation to the integral,
!  and estimate the absolute error.
!
write(*,'(1p,a,3g15.4)') '      qk15: ',centr,hlgth

  fc = f(centr)
write(*,'(1p,a,3g15.4)') '      qk15 2: ',centr,hlgth
  resg = fc*wg(4)
  resk = fc*wgk(8)
  resabs = abs(resk)
write(*,'(1p,a,3g15.4)') '      qk15 2: ',centr,hlgth

  do j = 1, 3
    jtw = j*2
    absc = hlgth*xgk(jtw)
    fval1 = f(centr-absc)
    fval2 = f(centr+absc)
    fv1(jtw) = fval1
    fv2(jtw) = fval2
    fsum = fval1+fval2
    resg = resg+wg(j)*fsum
    resk = resk+wgk(jtw)*fsum
    resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
  end do
write(*,'(1p,a,3g15.4)') '      qk15 3: ',centr,hlgth

  do j = 1, 4
    jtwm1 = j*2-1
    absc = hlgth*xgk(jtwm1)
    fval1 = f(centr-absc)
    fval2 = f(centr+absc)
    fv1(jtwm1) = fval1
    fv2(jtwm1) = fval2
    fsum = fval1+fval2
    resk = resk+wgk(jtwm1)*fsum
    resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
  end do

  reskh = resk * 5.0e-01
  resasc = wgk(8)*abs(fc-reskh)
write(*,'(1p,a,3g15.4)') '      qk15 4: ',centr,hlgth

  do j = 1, 7
    resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
  end do

  result = resk*hlgth
  resabs = resabs*dhlgth
  resasc = resasc*dhlgth
  abserr = abs((resk-resg)*hlgth)
write(*,'(1p,a,3g15.4)') '      qk15 5: ',centr,hlgth

  if ( resasc /= 0.0e+00.and.abserr /= 0.0e+00 ) then
    abserr = resasc*min ( 1.0e+00,(2.0e+02*abserr/resasc)**1.5e+00)
  end if

  if ( resabs > tiny ( resabs ) / (5.0e+01* epsilon ( resabs ) ) ) then
    abserr = max (( epsilon ( resabs ) *5.0e+01)*resabs,abserr)
  end if

  return
end subroutine qk15
subroutine qk21 ( f, a, b, result, abserr, resabs, resasc )

!*****************************************************************************80
!
!! QK21 carries out a 21 point Gauss-Kronrod quadrature rule.
!
!  Discussion:
!
!    This routine approximates
!      I = integral ( A <= X <= B ) F(X) dx
!    with an error estimate, and
!      J = integral ( A <= X <= B ) | F(X) | dx
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, external real*8 F, the name of the function routine, of the form
!      function f ( x )
!      real*8 f
!      real*8 x
!    which evaluates the integrand function.
!
!    Input, real*8 A, B, the limits of integration.
!
!    Output, real*8 RESULT, the estimated value of the integral.
!    RESULT is computed by applying the 21-point Kronrod rule (resk) 
!    obtained by optimal addition of abscissae to the 10-point Gauss 
!    rule (resg).
!
!    Output, real*8 ABSERR, an estimate of | I - RESULT |.
!
!    Output, real*8 RESABS, approximation to the integral of the absolute
!    value of F.
!
!    Output, real*8 RESASC, approximation to the integral | F-I/(B-A) | 
!    over [A,B].
!
  implicit none

  real*8 a
  real*8 absc
  real*8 abserr
  real*8 b
  real*8 centr
  real*8 dhlgth
  real*8, external :: f
  real*8 fc
  real*8 fsum
  real*8 fval1
  real*8 fval2
  real*8 fv1(10)
  real*8 fv2(10)
  real*8 hlgth
  integer j
  integer jtw
  integer jtwm1
  real*8 resabs
  real*8 resasc
  real*8 resg
  real*8 resk
  real*8 reskh
  real*8 result
  real*8 wg(5)
  real*8 wgk(11)
  real*8 xgk(11)
!
!           the abscissae and weights are given for the interval (-1,1).
!           because of symmetry only the positive abscissae and their
!           corresponding weights are given.
!
!           xgk    - abscissae of the 21-point Kronrod rule
!                    xgk(2), xgk(4), ...  abscissae of the 10-point
!                    Gauss rule
!                    xgk(1), xgk(3), ...  abscissae which are optimally
!                    added to the 10-point Gauss rule
!
!           wgk    - weights of the 21-point Kronrod rule
!
!           wg     - weights of the 10-point Gauss rule
!
  data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8), &
    xgk(9),xgk(10),xgk(11)/ &
       9.956571630258081e-01,     9.739065285171717e-01, &
       9.301574913557082e-01,     8.650633666889845e-01, &
       7.808177265864169e-01,     6.794095682990244e-01, &
       5.627571346686047e-01,     4.333953941292472e-01, &
       2.943928627014602e-01,     1.488743389816312e-01, &
       0.000000000000000e+00/
!
  data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8), &
    wgk(9),wgk(10),wgk(11)/ &
       1.169463886737187e-02,     3.255816230796473e-02, &
       5.475589657435200e-02,     7.503967481091995e-02, &
       9.312545458369761e-02,     1.093871588022976e-01, &
       1.234919762620659e-01,     1.347092173114733e-01, &
       1.427759385770601e-01,     1.477391049013385e-01, &
       1.494455540029169e-01/
!
  data wg(1),wg(2),wg(3),wg(4),wg(5)/ &
       6.667134430868814e-02,     1.494513491505806e-01, &
       2.190863625159820e-01,     2.692667193099964e-01, &
       2.955242247147529e-01/
!
!
!           list of major variables
!
!           centr  - mid point of the interval
!           hlgth  - half-length of the interval
!           absc   - abscissa
!           fval*  - function value
!           resg   - result of the 10-point Gauss formula
!           resk   - result of the 21-point Kronrod formula
!           reskh  - approximation to the mean value of f over (a,b),
!                    i.e. to i/(b-a)
!
  centr = 5.0e-01*(a+b)
  hlgth = 5.0e-01*(b-a)
  dhlgth = abs(hlgth)
!
!  Compute the 21-point Kronrod approximation to the
!  integral, and estimate the absolute error.
!
  resg = 0.0e+00
  fc = f(centr)
  resk = wgk(11)*fc
  resabs = abs(resk)

  do j = 1, 5
    jtw = 2*j
    absc = hlgth*xgk(jtw)
    fval1 = f(centr-absc)
    fval2 = f(centr+absc)
    fv1(jtw) = fval1
    fv2(jtw) = fval2
    fsum = fval1+fval2
    resg = resg+wg(j)*fsum
    resk = resk+wgk(jtw)*fsum
    resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
  end do

  do j = 1, 5
    jtwm1 = 2*j-1
    absc = hlgth*xgk(jtwm1)
    fval1 = f(centr-absc)
    fval2 = f(centr+absc)
    fv1(jtwm1) = fval1
    fv2(jtwm1) = fval2
    fsum = fval1+fval2
    resk = resk+wgk(jtwm1)*fsum
    resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
  end do

  reskh = resk*5.0e-01
  resasc = wgk(11)*abs(fc-reskh)

  do j = 1, 10
    resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
  end do

  result = resk*hlgth
  resabs = resabs*dhlgth
  resasc = resasc*dhlgth
  abserr = abs((resk-resg)*hlgth)

  if ( resasc /= 0.0e+00.and.abserr /= 0.0e+00) then
    abserr = resasc*min ( 1.0e+00,(2.0e+02*abserr/resasc)**1.5e+00)
  end if

  if ( resabs > tiny ( resabs ) /(5.0e+01* epsilon ( resabs ) )) then
    abserr = max (( epsilon ( resabs ) *5.0e+01)*resabs,abserr)
  end if

  return
end subroutine qk21
subroutine qsort ( limit, last, maxerr, ermax, elist, iord, nrmax )

!*****************************************************************************80
!
!! QSORT maintains the order of a list of local error estimates.
!
!  Discussion:
!
!    This routine maintains the descending ordering in the list of the 
!    local error estimates resulting from the interval subdivision process. 
!    At each call two error estimates are inserted using the sequential 
!    search top-down for the largest error estimate and bottom-up for the
!    smallest error estimate.
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, integer LIMIT, the maximum number of error estimates the list can
!    contain.
!
!    Input, integer LAST, the current number of error estimates.
!
!    Input/output, integer MAXERR, the index in the list of the NRMAX-th 
!    largest error.
!
!    Output, real*8 ERMAX, the NRMAX-th largest error = ELIST(MAXERR).
!
!    Input, real*8 ELIST(LIMIT), contains the error estimates.
!
!    Input/output, integer IORD(LAST).  The first K elements contain 
!    pointers to the error estimates such that ELIST(IORD(1)) through
!    ELIST(IORD(K)) form a decreasing sequence, with
!      K = LAST 
!    if 
!      LAST <= (LIMIT/2+2), 
!    and otherwise
!      K = LIMIT+1-LAST.
!
!    Input/output, integer NRMAX.
!
  implicit none

  integer last

  real*8 elist(last)
  real*8 ermax
  real*8 errmax
  real*8 errmin
  integer i
  integer ibeg
  integer iord(last)
  integer isucc
  integer j
  integer jbnd
  integer jupbn
  integer k
  integer limit
  integer maxerr
  integer nrmax
!
!  Check whether the list contains more than two error estimates.
!
  if ( last <= 2 ) then
    iord(1) = 1
    iord(2) = 2
    go to 90
  end if
!
!  This part of the routine is only executed if, due to a
!  difficult integrand, subdivision increased the error
!  estimate. in the normal case the insert procedure should
!  start after the nrmax-th largest error estimate.
!
  errmax = elist(maxerr)

  do i = 1, nrmax-1

    isucc = iord(nrmax-1)

    if ( errmax <= elist(isucc) ) then
      exit
    end if

    iord(nrmax) = isucc
    nrmax = nrmax-1

  end do
!
!  Compute the number of elements in the list to be maintained
!  in descending order.  This number depends on the number of
!  subdivisions still allowed.
!
  jupbn = last

  if ( (limit/2+2) < last ) then
    jupbn = limit+3-last
  end if

  errmin = elist(last)
!
!  Insert errmax by traversing the list top-down, starting
!  comparison from the element elist(iord(nrmax+1)).
!
  jbnd = jupbn-1
  ibeg = nrmax+1

  do i = ibeg, jbnd
    isucc = iord(i)
    if ( elist(isucc) <= errmax ) then
      go to 60
    end if
    iord(i-1) = isucc
  end do

  iord(jbnd) = maxerr
  iord(jupbn) = last
  go to 90
!
!  Insert errmin by traversing the list bottom-up.
!
60 continue

  iord(i-1) = maxerr
  k = jbnd

  do j = i, jbnd
    isucc = iord(k)
    if ( errmin < elist(isucc) ) then
      go to 80
    end if
    iord(k+1) = isucc
    k = k-1
  end do

  iord(i) = last
  go to 90

80 continue

  iord(k+1) = last
!
!  Set maxerr and ermax.
!
90 continue

  maxerr = iord(nrmax)
  ermax = elist(maxerr)

  return
end subroutine qsort

!----------------------------------------------------
Real*8      FUNCTION Hgrow(x)
     Real*8 :: x
         Hgrow =(sqrt(x/(1.+x**3)))**3 
     End FUNCTION Hgrow

!-------------------------------------- P*k^2*Top-Hat Filter
REAL*8 FUNCTION Ptophat(x)
!---------------------------------------
     use setHalos
     Real*8 :: wk,x,TOPHAT,Ppk
        wk = x/rf
	    TOPHAT =( (SIN(x)-x*COS(x))*3.d0/x**3 )**2
         !   write(*,*) '       ',wk,TOPHAT
        
            Ptophat=Ppk(wk)*x**2*TOPHAT
    END FUNCTION Ptophat
