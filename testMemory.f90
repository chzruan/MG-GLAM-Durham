  !--------- quasi random
  Integer*4 :: Ngrid =5000

  Real*4, allocatable :: Fi(:,:,:)
  real*8 :: rms
  Allocate(Fi(Ngrid,Ngrid,Ngrid))

  write(*,*) ' Started: Ngrid=',Ngrid
  write(*,*) '     Memory = ',4.*(Ngrid/1024.)**3,' Gb'
  Do m= 1,5
    t0 = seconds()
!$OMP PARALLEL DO DEFAULT(SHARED) &         
!$OMP PRIVATE (k,j,i)
     Do k=1,Ngrid
        Do j=1,Ngrid
           Do i=1,Ngrid
              Fi(i,j,k) =-1.
           EndDo
        end Do
     end Do
    t1 = seconds()
    write(*,*) ' Finished :',t1-t0
  End Do
write(*,*) ' Finished '
  stop
end program

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
