!------------------------------------------
!
!
!   Read tables of Power Spectra
!   find average and rms
!   If required, include linear prediction
!
!
!------------------------------------------
Module Param
      Integer*4, parameter :: NtabM = 100000
      Real*4  :: Om,        &      ! cosmology
                 Omb,       &
                 Omc,       &
                 OmL,       &
                 AEXPN,     &
                 Sn,        &      ! normalization
                 sigma8,    &
                 hubble

      Real*8   :: rf,OmLOm0
      Real*4   :: xkt(0:NtabM),Pkt(0:NtabM)  ! Power spectrum
      Real*4   :: StepK,alog0
      INteger*4:: Ntab

    end Module Param
!
!__________________________________________________________
!
Program Average
  integer, parameter :: Nmax=300000
  Real*8 :: Aver(0:Nmax), RMS(0:Nmax),k1(0:Nmax),k2(0:Nmax),k3(0:Nmax),aa,rr,Plin
  Integer*8 :: nr(0:Nmax),N
  Character*80 :: Names,Fname,Line
  
  Aver = 0.; RMS = 0.
  rad = 0. ; nr = 0

  open(1,file='list.dat')
  read(1,*) Ntab

  Do i=1,Ntab                !--- read files
     Read(1,'(a)') Names
     Open(10,file=TRIM(Names))
     If(i==1)Then
        write(Fname,'(a,a)')Trim(Names),'.Aver'
        open(20,file=TRIM(Fname))
        write(*,'(2a)') ' Writing averaged bispectrum to file:',TRIM(Fname)
     End If
     Do j=1,4                 !--- read the header
        Read(10,'(a)')Line
        If(i==1)Then
           If(j<=3)write(20,*)TRIM(Line)
           If(j==4)Then
              write(20,*) ' Realizations  = ',Ntab
              write(20,'(2a)')'      k1          k2          k3 ', &
              '         dN      Bispectrum      rms'
           endIf
        end If
     EndDo ! j
     k =0
     Do                      !--- read data from each file
        read(10,*,iostat=ierr) ak1, ak2, ak3,dN,P
        If(ierr/=0)exit
           k = k +1
           if(i==1)Then
              k1(k) = ak1
              k2(k) = ak2
              k3(k) = ak3
              nr(k)  = dN
           end if
           Aver(k) = Aver(k) + P
           RMS(k)  = RMS(k)  + P**2
     end Do
     close(10)
  EndDo  ! Ntab

  write(*,*) ' k    =',k
  write(*,*) ' Ntab =' ,Ntab
     Nn = k

     Do i=1,Nn
        aa = Aver(i)/Ntab
        rr = sqrt(max(RMS(i)/Ntab - aa**2,1.d-20))
        write(20,'(3f11.7,i13,3es12.4)') k1(i),k2(i),k3(i),nr(i),aa,rr 
     end Do
     
end Program Average
!
!
!
