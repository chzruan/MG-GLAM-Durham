  !--------------------
  !
  ! OMP parallization of linker list
  !
  !---------------------
Program OOO
  Integer*4, parameter :: N =7000
  Real*4, save :: A(N,N,N)
  Real*4, allocatable, dimension(:)    :: X
  Integer*4, allocatable, dimension(:)  :: Mth
  Real*4, allocatable, dimension(:,:)  :: XX
  Real*8 :: pi =3.14159265
  Real*8 :: D
  Integer*8 :: M,L,i,j,Mnow
  Integer*4 :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM,Nthreads

  Do iter =1,4
  tstart = seconds()
  Alimit =0.99991 ; D =0.
  M  =0
!$OMP PARALLEL
  Nthreads =OMP_GET_NUM_THREADS()
!$OMP end parallel
  write(*,*) ' Number of threads =',Nthreads

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k)
  Do k=1,N
     Do j=1,N
        Do i=1,N
           A(i,j,k) =sin(2.d0*pi/N*float(i+j+k)/3.528)
        end Do
     end Do
  end Do

 t0 = seconds()
  write(*,*) ' time for init array    =',t0-tstart
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k) REDUCTION(+:M,D)
  Do k=1,N
     Do j=1,N
        Do i=1,N
           If(A(i,j,k)>Alimit)Then
              M = M + 1
              D = D + A(i,j,k)
              !write(*,'(3i4,f8.3,i4,3x,i4)') i,j,k,A(i,j,k),Mthread,Mm
           end If
        end Do
     end Do
  end Do
  t1 =  seconds()
  write(*,*) ' time for counting cells    =',t1-t0
  write(*,'(a,f8.4,a,i12,a,es13.5)') ' Limit=',Alimit,' Number of elements= ',M, &
       ' Average of elements= ',D/max(M,1)
  L = M
  Nbuff = L/Nthreads *5  ! length of the buffer
  write(*,'(a,i10,a,es12.4)') ' Length of buffer=',Nbuff, &
       ' Memory needed/MB=',4.*Nbuff*Nthreads/1024.**3
  
  Allocate(X(L),Mth(Nthreads),XX(Nbuff,Nthreads))
  M      = 0
  X(:)   = -100.
  Mth(:) = 0
  Mm     = 1
  write(*,*) Nbuff,Nthreads
!$OMP PARALLEL DO PRIVATE(i,j,k,Mthread) Firstprivate(Mm) 
  Do k=1,N
     Do j=1,N
        Do i=1,N
           If(A(i,j,k)>Alimit)Then
               Mthread       = OMP_GET_THREAD_NUM()+1
              Mth(Mthread)   =   Mm
              !write(*,'(3i4,f8.3,i4,3x,i4)') i,j,k,A(i,j,k),Mthread,Mm
              If(Mm>Nbuff)Then
                 write(*,'(3i6,i11,i4)') i,j,k,Mm,Mthread
                 Stop ' Too many elements for buffer XX'
              end If
              XX(Mm,Mthread) = A(i,j,k)
              Mm             =   Mm +1
           end If
        end Do
     end Do
  end Do
  t2 =  seconds()
  write(*,*) ' time for parallel counting cells =',t2-t1  
  write(*,'(a,/(20i11))') ' Number of elements in each thread =',Mth
  M = 0
  Do i=1,Nthreads
     M = M + Mth(i)
  end Do
  write(*,*) ' Elements =',M,L
  
  i = 0
  Do Mthread=1,Nthreads
     Mnow = Mth(Mthread)
     Do j = 1,Mnow
        X(j+i) = XX(j,Mthread)
     End Do
     i = i+ Mnow
  end Do
  write(*,*) ' Elements =',M,i,L
  D =0.
  Do i=1,M
     D= D+X(i)
  end Do
  Deallocate(XX,Mth,X)

  write(*,'(a,i12,a,es13.5//)') '  Number of elements= ',M, &
               ' Average of elements= ',D/max(M,1)
  
end Do
end Program OOO
!
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
