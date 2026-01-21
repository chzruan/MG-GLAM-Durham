!--------------------------------------------------
!
!        Read PM data and make 3D overdensity grid
!
!--------------------------------------------------
Program Density
  use Tools
  
  Path = ''
  open(17,file='output')
  
   Call Initialize(Path)
   Call DENSIT
   
end Program Density

!--------------------------------------------
subroutine Initialize(Path)
!--------------------------------------------
use Tools

        logical :: exst
        character(len=*) :: Path

      WRITE (*,'(A,$)') ' Enter snapshot number to analyze => '
      READ (*,*) moment	 ! read snapshot number
      
      CALL ReadDataPM(moment,Path)
      myMemory =Memory(1_8*NGRID*NGRID*NGRID)
      Allocate (FI(NGRID,NGRID,NGRID))
         
end subroutine Initialize

!-------------------------------------------------------------
!
!         Make density contrast FI(NGRID,NGRID,NGRID)
!         for Nparticles XPAR,YPAR,ZPAR(Nparticles) in 1-Ngrid
!
SUBROUTINE DENSIT
    
!-------------------------------------------------------------
use Tools
  real*8 :: XN,YN, &
     X,Y,Z,D1,D2,D3,T1,T2,T3,T2W,D2W
integer*8 :: IN

    XN   =FLOAT(NGRID)+1.-1.E-7
    YN   =FLOAT(NGRID)
    Wpar = YN**3/FLOAT(Nparticles)
				!       Subtract mean density
write(*,'(a,i5,a,i12)') ' Init density: Ngrid= ',Ngrid,' Nparticles= ',Nparticles
write(17,'(a,i5,a,i12)') ' Init density: Ngrid= ',Ngrid,' Nparticles= ',Nparticles
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3)
    DO M3=1,NGRID
       DO M2=1,NGRID
          DO M1=1,NGRID
	        FI(M1,M2,M3) = -1.
	  END DO
       END DO
    END DO
    
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (IN,X,Y,Z,D1,D2,D3,T1,T2,T3,T2W,D2W) &
!$OMP PRIVATE (I,J,K,I1,J1,K1)    
     DO   IN=1,Nparticles         ! loop over particles 
	  X=XPAR(IN)
	  Y=YPAR(IN)
	  Z=ZPAR(IN)
	  I=INT(X)
	  J=INT(Y)
	  K=INT(Z)

	  D1=X-FLOAT(I)
	  D2=Y-FLOAT(J)
	  D3=Z-FLOAT(K)
	  T1=1.-D1
	  T2=1.-D2
	  T3=1.-D3
	  T2W =T2*WPAR
	  D2W =D2*WPAR
	  I1=I+1
	     IF(I1.GT.NGRID)I1=1
	  J1=J+1
	     IF(J1.GT.NGRID)J1=1
	  K1=K+1
          IF(K1.GT.NGRID)K1=1
!$OMP ATOMIC
          FI(I ,J ,K ) =FI(I ,J ,K ) +T3*T1*T2W
!$OMP ATOMIC
	  FI(I1,J ,K ) =FI(I1,J ,K ) +T3*D1*T2W
!$OMP ATOMIC
	  FI(I ,J1,K ) =FI(I ,J1,K ) +T3*T1*D2W
!$OMP ATOMIC
	  FI(I1,J1,K ) =FI(I1,J1,K ) +T3*D1*D2W
   
!$OMP ATOMIC
	  FI(I ,J ,K1) =FI(I ,J ,K1) +D3*T1*T2W
!$OMP ATOMIC
	  FI(I1,J ,K1) =FI(I1,J ,K1) +D3*D1*T2W
!$OMP ATOMIC
	  FI(I ,J1,K1) =FI(I ,J1,K1) +D3*T1*D2W
!$OMP ATOMIC
	  FI(I1,J1,K1) =FI(I1,J1,K1) +D3*D1*D2W
	   
      ENDDO 
    Call TimingMain(3,1)
    D1 = 0. ; D2 =0.
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3) REDUCTION(+:D1,D2)
    DO M3=1,NGRID
       DO M2=1,NGRID
          DO M1=1,NGRID
	       D1 = D1 + FI(M1,M2,M3)
	       D2 = D2 + FI(M1,M2,M3)**2
          END DO
       END DO
    END DO
    D1 = D1/(float(NGRID))**3
    D2 = sqrt(D2/(float(NGRID))**3)
    write(*,*) ' Finished density: average/rms= ',D1,D2
    write(17,*) ' Finished density: average/rms= ',D1,D2
    
  END SUBROUTINE DENSIT
