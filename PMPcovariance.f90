  !-------------------------------------------
  !
  !   Read PM Power.DAT table and interpolate to
  !               k values given in another Power table
  !
!-------------------------------------------
Module Structures
  Real*8,    allocatable, dimension(:,:) :: wKd,wNd,Prd,Cov,CovF,CovE
  Integer*4, allocatable, dimension(:,:) :: Nd
  Real*8,    allocatable, dimension(:) :: Prf,wk,CovD,CovA
  Integer*4, allocatable, dimension(:) :: Nf
  
  Character*120 :: Line(10000),Name

  Integer, parameter :: Nmax = 500
  Integer*4  :: Nfiles,iFr,iSm
  Integer*4  :: iFilter = 1   ! width of filter in k
end Module Structures
!-------------------------------------------
Program Covariance
  use Structures

   write(*,'(a,$)') ' Enter name of file with the list of spectra : '
   read(*,'(a)') Name
   write(*,*) Trim(Name)
   open(40,file=TRIM(Name))
  
  read(40,*) Nfiles
  !Nfiles = 1800
  write(*,'(a,i4)') ' Nfiles to analyze= ',Nfiles
  allocate(Nd(Nmax,Nfiles),wKd(Nmax,Nfiles))
  allocate(wNd(Nmax,Nfiles),Prd(Nmax,Nfiles))
  allocate(Cov(Nmax,Nmax),CovF(Nmax,Nmax),CovE(Nmax,Nmax),CovD(Nmax),CovA(Nmax))
  allocate(Prf(Nmax),wk(Nmax),Nf(Nmax))

  
  Do i=1,Nfiles             !--- read data
     read(40,'(a)') Line(i)
     Open(20,file=TRIM(Line(i)))
     Call ReadData(i)
     close(20)
  EndDo
  close(40)
      Nf(:) = wNd(:,1)
      wk(:) = wKd(:,1)
   iFr = 96 !40              !- central bin
   iSm = 0
   write(*,'(a,$)') ' Enter k and width in k for averaging: '
   read(*,*) ValK,iSm
   iFilter  = iSm
   do i=1,Nmax
      If(wk(i)>ValK)exit
   end do
   iFr = i
   if(wk(iFr+1)-ValK < ValK-wk(iFr))iFr = iFr+1
   write(*,*) ' iFr= ',iFr,iSm
   write(Name,'(a,f5.3,a,i2.2,a)')'Covar.',ValK,'.',iSm,'.dat'
  open(30,file=TRIM(Name))
  Open(20,file=TRIM(Line(1)))
  open(10,file='PowerAver.dat')
  Do i=1,3
     Read(20,'(a)')Name
     if(i<3)write(30,'(a)')TRIM(Name)
     if(i<3)write(10,'(a)')TRIM(Name)
     if(i<3)write(*,'(a)')TRIM(Name)
  EndDo
  close(20)

!  write(1,'(7x,a,3x,a,3(10x,a),3x,6(8x,a))')'k/Mpch','Nharmon',&
!          a  'Pk','errPk'
  Call MakeCov

end Program Covariance
!------------------------------------------------
!
!           Read One File
!
Subroutine  ReadData(iFile)
  use Structures
  Real, SAVE ::  Box,Res
  Integer, save ::  Np
  character*80 :: rLine
  real*8 :: a,b
  
  Do i=1,100      !--- find the first line
     Read(20,*,iostat=ierr) n
     if(ierr ==0)Then
        !write(*,*) ' Line=',i, ' n= ',n
        exit
     End if
  end Do
  iL = i
  backspace(20)
  If(iFile == 1)Then
     Box =2500.
     Res =3300.
     Np  =1000
     write(*,'(a,$)') ' Enter: Box, Mesh, Nparticles ='
     read(*,*) Box,Res,Np
  End If
  Pk_corr = (Box/Np)**3
  fNy  = 2.*3.14159/Box*Res/2.
  !write(*,*) ' fNy =',fNy
  Do i =1,Nmax
     !write(*,*) '       read i=',i
     Read(20,*,iostat=ierr) n,w,wn,Pr  !!,P0,P2
     if(ierr ==0)Then
        Nd(i,iFile) = n       ! Line number
        wKd(i,iFile)= w      ! k 
        wNd(i,iFile)= wn      ! N harmonics
        Prd(i,iFile)= Pr/((1-0.667*(sin(3.1415*w/fNy/2.))**2))-Pk_corr      ! Pk      
        !write(*,'(12es12.4)')w,((1-0.667*(sin(3.1415*w/fNy/2.))**2))
     Else
        exit
     End if
  end Do
  !write(*,*) ' finished reading iFile= ',iFile

end Subroutine ReadData
!------------------------------------------------
!
!           Make Covariance matrix
!
Subroutine  makeCov
  use Structures
  Real*8 :: a,b,c,ea,eb,ec

  Prf(:)  = 0.
  CovD(:) = 0.
  !write(*,*) ' make covariance maxtrix'
  Do i= 1,Nmax     !--- find <P(i)>
        a =0.
        b =0.
     Do iFile=1,Nfiles        
        a = a + Prd(i,iFile)
        b = b + Prd(i,iFile)**2
        !write(*,'(3i6,es13.3,3x,6es12.4)') i,j,iFile,a,Prd(i,iFile),Prd(j,iFile)
     EndDo
     Prf(i) = a/Nfiles    !--- Prf = average Pk over all realizations
     CovD(i) = sqrt(b/Nfiles -(a/Nfiles)**2)    !--- Prf = average Pk over all realizations
  end Do
  write(10,'(6x,a)') ' k        N          Pk           errPk'
  Do i=  1,Nmax
     write(10,'(es13.5,i9,2e13.5)') wk(i), Nf(i), Prf(i), CovD(i)
  end Do
  close(10)
  CovD(:) =0.
  
  Do i= 1,Nmax    !--- find <P_i P_j>
  Do j= 1,Nmax
     a =0.
     Do iFile=1,Nfiles
        a = a + Prd(i,iFile)*Prd(j,iFile)
     EndDo
     Cov(i,j) = a/(Nfiles)-(Prf(i)*Prf(j))
  end Do
  end Do

  write(10,'(6x,a)') ' k        N          Pk           errPk'
  Do i=  1,Nmax
     write(10,'(es13.5,i9,es12.5,18f13.6)') wk(i), Nf(i), Prf(i), (wk(j),1.e6*Cov(i,j)/Prf(i)/Prf(j),j=86,100,10)
  end Do
  close(10)
  
  Do i=1,Nmax     !--- store the diagonal components
     CovD(i) = Cov(i,i)
  EnDdo

  
  Do i= 1,Nmax    !--- normalize <P_i P_j>/(<P_i> <P_j>)
     Do j= 1,Nmax
        a = Cov(i,j)
        b = sqrt((Cov(i,j)**2+CovD(i)*CovD(j))/Nfiles/(2*iFilter+1)**2)
        CovE(i,j) = sqrt((abs(a)+b)/(Prf(i)*Prf(j))) - sqrt(abs(a)/(Prf(i)*Prf(j)))
        Cov(i,j) = a/(Prf(i)*Prf(j))
     end Do
  end Do
  CovA(:) = Cov(:,iFr)
  
  If(iFilter/=0)Then  !---- Smooth Cov
     
     CovA(:) =0.      !--- smooth Cov(i,:) over iNfilter values of k
     Do i=1,iFilter              ! store first iFilter elements
        CovA(i) = Cov(i,iFr)
     EndDo
     Do i=Nmax-iFilter+1,Nmax    ! store last iFilter elements
        CovA(i) = Cov(i,iFr)
     EndDo
     
     Do i=1+iFilter,Nmax-iFilter-1          ! sum all elements along diagonal         
       Do j=-iFilter,iFilter
          CovA(i) =CovA(i) +Cov(i+j,iFr+j)
       EndDo
     End Do
     CovA(:) = CovA(:)/(2.*iFilter+1)       ! find average
     CovD(:) = CovA(:)

     Do i=1+iFilter,Nmax-iFilter-1
       If(abs(i-iFr) > 4 )Then
        CovD(i) =0.
        Do j=-iFilter,iFilter
          CovD(i) =CovD(i) +CovA(i+j)
       EndDo
       CovD(i) = CovD(i)/(2.*iFilter+1)
      End If
     End Do
     CovA(:) = CovD(:)
 End If
   
   write(30,'(a,f8.4,a,i3,a,2f8.4,a,i5)')      &
        ' k=',wk(iFr),' width= ',iSm,' k_limits=',wk(iFr-iSm),wk(iFr+iSm), &
        ' Nfiles = ',Nfiles
   write(30,*) '      k        CovarCoefficient'
  do i=1,Nmax
     write(30,*) wk(i),sqrt(abs(CovA(i))),sqrt(CovE(i,iFr))
  end do
  !write(*,'(9x,20es12.4)') (wk(i),i=iFr-iSm,iFr+iSm,1)
  !   do j=1,200
  !      write(*,'(f9.4,20es12.4)') wk(j),(Cov(i,j),i=iFr-iSm,iFr+iSm,1)
  !   End do
end Subroutine makeCov
