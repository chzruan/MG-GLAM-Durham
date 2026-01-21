  !------------------------------------------
  !
  !
  !   Read tables of Power Spectra
  !   find average and rms
  !------------------------------------------
Program Average
  integer, parameter :: Nmax=10000
  Real*8 :: Aver(0:Nmax), RMS(0:Nmax),ak(0:Nmax),aa,rr
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
        write(Fname,'(a15,a)')Trim(Names),'Aver.dat'
        open(20,file=TRIM(Fname))
     End If
     Do j=1,3                 !--- read the header
        Read(10,'(a)')Line
        If(i==1)Then
           If(j<=2)write(20,*)TRIM(Line)
           If(j==3)Then
              write(20,*) ' Realizations  = ',Ntab
              write(20,'(a)')'   k         dN         Pk        rms  '
           endIf
        end If
     EndDo ! j
     k =0
     Do                      !--- read data from each file
        read(10,*,iostat=ierr)ii,wk,N,P
        If(ierr/=0)exit
           k = k +1
           if(i==1)Then
              ak(k) = wk
              nr(k)  = N
           end if
           Aver(k) = Aver(k) + P
           RMS(k)  = RMS(k)  + P**2
     end Do
     close(10)
  EndDo  ! Ntab

  Nn = k
  Do i=1,Nn
     aa = Aver(i)/Ntab
     rr = sqrt(max(RMS(i)/Ntab - aa**2,1.d-20))
     write(20,'(f9.5,i13,3es13.5)') ak(i),nr(i),aa,rr 
  end Do
end Program Average
  
