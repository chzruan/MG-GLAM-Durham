!--------------------------------------------------
!
!        Read halo/particle sample, make 3D overdensity grid,
!        Calculate power spectrum
!           
!
!--------------------------------------------------
Module Parameters
  Integer*4 :: iFlag = 0, Ng, NGRID_old, moment
  Integer*8 :: Nparticles,Ntotal
  Real*4, Allocatable,  Dimension(:,:,:) :: Rho1,Rho2,Rho3,Rho4,Rho5,Rho6
  Real*4, allocatable,  Dimension(:) ::  XPAR,YPAR,ZPAR,VX,VY,VZ
  Real*4   :: AEXPN2,Box2
end Module Parameters

!--------------------------------------------------
Program Spectrum
  use Parameters
  
     t0 = seconds()
  write(*,'(a,$)') ' Enter grid size Ng        = '
  read(*,*) Ng
  
  Call InitializeGLAM        !--- comment out this line if gadget format is used
  Call RescaleGLAM
  NGRID = Ng
  ! Call InitializeGadget    !--- uncomment these two line for Gadget format
  ! NGRID = Ng
     t1 = seconds()
     write(*,*) ' time to read data =',t1-t0,' secs'
     
                             !--- read particles and update meshes as many times as needed
   ! do Iread =1,Ineed  
   Call DensMeshes         
     t2 = seconds()
     write(*,*) ' time to make density =',t2-t1,' secs'
   ! End Do

   ! CALL MPI_REDUCE( ... MPI_SUM ...) !--- sum up all density meshes
   !                                        result is placed in root MPI task  

   !  if(mpi_node == root)Then
        Call NormalizeMeshes
        Call GetPowerMulti
   !   endif
     t3 = seconds()
     write(*,*) ' time to make Pk      =',t3-t2,' secs'

     
     write(*,*) ' time for the run     =',t3-t0,' secs'

      DEALLOCATE(Rho1,Rho2,Rho3,Rho4,Rho5,Rho6)
            
end Program Spectrum

!--------------------------------------------
subroutine InitializeGLAM
!--------------------------------------------
  use Parameters

  write(*,'(a,$)') ' Enter  moment    number = '
  read(*,*) moment
  OPEN(17,file='Out.dat')
  CALL ReadDataPM          !--- read coordinates and velocities
                           !    initialize parameters:
                           !       Ng,Box2, AEXPN2, Nparticles, Ntotal
   NGRID    = Ng
   Box      = Box2  
   AEXPN    = AEXPN2
   Ntotal   = Nparticles   !--- for GLAM all particles are read in one call

     ALLOCATE(Rho1(NGRID,NGRID,NGRID),Rho2(NGRID,NGRID,NGRID))
     ALLOCATE(Rho3(NGRID,NGRID,NGRID),Rho4(NGRID,NGRID,NGRID))
     ALLOCATE(Rho5(NGRID,NGRID,NGRID),Rho6(NGRID,NGRID,NGRID))
!				       initiate density
!$OMP PARALLEL DO DEFAULT(SHARED) & 
!$OMP PRIVATE (M1,M2,M3)
      DO M3=1,NGRID
       DO M2=1,NGRID
	 DO M1=1,NGRID
          Rho1(M1,M2,M3) = 0.
          Rho2(M1,M2,M3) = 0.
          Rho3(M1,M2,M3) = 0.
          Rho4(M1,M2,M3) = 0.
          Rho5(M1,M2,M3) = 0.
          Rho6(M1,M2,M3) = 0.
	 END DO
       END DO
      END DO

   write(*,'(a,f8.2,f8.4,i6)') ' read: Box,a,Ngrid = ',Box,AEXPN,NGRID
end subroutine InitializeGLAM

!--------------------------------------------
subroutine RescaleGadget
!        
!      Input coordinates Xpar,Ypar,Zpar are
!            in the range of 0--Box
!      Output coordinate are 1--(Ngrid+1)
!  
!--------------------------------------------
  use Parameters

  Box    = Box2
  NGRID  = Ng
  Xscale = NGRID/Box
!$OMP PARALLEL DO DEFAULT(SHARED) & 
!$OMP PRIVATE (ip,X,Y,Z)
      Do ip = 1,Nparticles
            X = Xpar(ip)*Xscale +1. !--- rescale
            Y = Ypar(ip)*Xscale +1.
            Z = Zpar(ip)*Xscale +1.
           IF(X.ge.NGRID+1.)X=X-NGRID  !--- periodical conditions
           IF(Y.ge.NGRID+1.)Y=Y-NGRID
           IF(Z.ge.NGRID+1.)Z=Z-NGRID

           IF(X.lt.1.)X=X+NGRID
           IF(Y.lt.1.)Y=Y+NGRID
           IF(Z.lt.1.)Z=Z+NGRID
           Xpar(ip) = X
           Ypar(ip) = Y
           Zpar(ip) = Z
        end Do
      end subroutine RescaleGadget
!--------------------------------------------
subroutine RescaleGLAM
!        
!      Input coordinates Xpar,Ypar,Zpar are
!            in the range of 1--(Ngid_old+1)
!      Output coordinate are 1--(Ngrid+1)
!  
!--------------------------------------------
  use Parameters
  Integer*8 :: ip
  
  NGRID  = Ng
  Xscale = NGRID/float(NGRID_old)
  write(*,*) ' NGRID    =',NGRID
  write(*,*) ' NGRIDold =',NGRID_old
!$OMP PARALLEL DO DEFAULT(SHARED) & 
!$OMP PRIVATE (ip,X,Y,Z)
      Do ip = 1,Nparticles
            X = (Xpar(ip)-1.)*Xscale +1.
            Y = (Ypar(ip)-1.)*Xscale +1.
            Z = (Zpar(ip)-1.)*Xscale +1.
           IF(X.ge.NGRID+1.)X=X-NGRID
           IF(Y.ge.NGRID+1.)Y=Y-NGRID
           IF(Z.ge.NGRID+1.)Z=Z-NGRID

           IF(X.lt.1.)X=X+NGRID
           IF(Y.lt.1.)Y=Y+NGRID
           IF(Z.lt.1.)Z=Z+NGRID
           Xpar(ip) = X
           Ypar(ip) = Y
           Zpar(ip) = Z
        end Do
                xm = 1.e8 ; ym = 1.e8; zm = 1.e8
!$OMP PARALLEL DO DEFAULT(SHARED)  PRIVATE (ip) Reduction(max:xxm,yym,zzm) &
!$OMP  Reduction(min:xm,ym,zm)
        do ip =1,Nparticles
           xxm =max(xxm,Xpar(ip))
           yym =max(yym,Ypar(ip))
           zzm =max(zzm,Zpar(ip))
           xm =min(xm,Xpar(ip))
           ym =min(ym,Ypar(ip))
           zm =min(zm,Zpar(ip))
           
          ! If(Xpar(ip).lt.1.0.or.Ypar(ip).lt.1.0.or.Zpar(ip).lt.1.0) &
          !   write(*,'(a,i10,1p,3g14.5)')' Error coord: ',ip,Xpar(ip),Ypar(ip),Zpar(ip)
        EndDo
     write(*,'(a,3(2es14.5,3x))') ' Coordinates Min/Max =',xm,xxm,ym,yym,zm,zzm      

      end subroutine RescaleGLAM
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
!          Make density field for stacked meshes
!          Nparticles = current number of particles
!
!
subroutine DensMeshes
!
!--------------------------------------------
  use Parameters
character*120 :: Name

Real*8 :: ss,X,Y,Z,D1,D2,D3,T1,T2,T3, &
             D4,D5,D6,T4,T5,T6,  &
             T1a,T2a,T3a,T4a,T5a,T6a, &
             W,Wmesh,Xscale
Integer*8 ::  ip

        NGRID = Ng
        XN   =FLOAT(NGRID)+1.-1.E-7
        YN   =FLOAT(NGRID)
        write(*,*) ' New Ngrid= ',NGRID
        write(*,*) ' Old Ngrid= ',NGRID_old

      W     = FLOAT(NGRID)**3/float(Ntotal)  !-- weight for first mesh
      Xscale = Ngrid/float(Ngrid_old)

Do Kmesh=1,6
      Imesh = 2**Kmesh/2
      Nmesh = Ngrid/Imesh  ! size of mesh at this level
      Wmesh = W *Imesh**1.5   ! weight of each particle
      write(*,'(i3,a,i4,a,es12.4)') Kmesh,' mesh size=',Nmesh,' Weight= ',Wmesh
!$OMP PARALLEL DO DEFAULT(SHARED) & 
!$OMP PRIVATE (ip,X,Y,Z,I,J,K,D1,D2,D3,T1,T2,T3,T2W,D2W,I1,J1,K1)
      Do ip = 1,Nparticles
         !if(mod(ip,50000)==0)write(*,*) ' galaxy=',ip
            X = Xpar(ip)-1.
            Y = Ypar(ip)-1.
            Z = Zpar(ip)-1.
            X = X -INT(X/Nmesh)*Nmesh   !-- 0...Nmesh
            Y = Y -INT(Y/Nmesh)*Nmesh 
            Z = Z -INT(Z/Nmesh)*Nmesh
            X = X*Imesh+1.              !-- 1...Ngrid
            Y = Y*Imesh+1. 
            Z = Z*Imesh+1.
            
            
           IF(X.ge.NGRID+1.)X=X-NGRID
           IF(Y.ge.NGRID+1.)Y=Y-NGRID
           IF(Z.ge.NGRID+1.)Z=Z-NGRID

           IF(X.lt.1.)X=X+NGRID
           IF(Y.lt.1.)Y=Y+NGRID
           IF(Z.lt.1.)Z=Z+NGRID

	   I=INT(X)
	   J=INT(Y)
	   K=INT(Z)
           If(I.le.0)write (*,'(a,3es15.5,a,2i12)') ' X:',X,Y,Z,' Irow=',ip
           If(J.le.0)write (*,'(a,3es15.5,a,2i12)') ' Y:',X,Y,Z,' Irow=',ip
           If(I.gt.NGRID+1)write (*,'(a,3es15.5,a,2i12)') ' X:',X,Y,Z,' Irow=',ip
           If(J.gt.NGRID+1)write (*,'(a,3es15.5,a,2i12)') ' Y:',X,Y,Z,' Irow=',ip
                   !---------------------------------------- CIC
	        D1=X-FLOAT(I)
	        D2=Y-FLOAT(J)
	        D3=Z-FLOAT(K)
	        T1=1.-D1
	        T2=1.-D2
	        T3=1.-D3
	        T2W =T2*Wmesh
	        D2W =D2*Wmesh
	        I1=I+1
	           IF(I1.GT.NGRID)I1=1
	        J1=J+1
	           IF(J1.GT.NGRID)J1=1
	        K1=K+1
         IF(K1.GT.NGRID)K1=1
         Select CASE (Kmesh)
            CASE (1)
!$OMP Atomic         
	             Rho1(I ,J ,K ) =Rho1(I ,J ,K ) +T3*T1*T2W
!$OMP Atomic         
                     Rho1(I1,J ,K ) =Rho1(I1,J ,K ) +T3*D1*T2W
!$OMP Atomic         
	             Rho1(I ,J1,K ) =Rho1(I ,J1,K ) +T3*T1*D2W
!$OMP Atomic         
	             Rho1(I1,J1,K ) =Rho1(I1,J1,K ) +T3*D1*D2W
!$OMP Atomic         
	             Rho1(I ,J ,K1) =Rho1(I ,J ,K1) +D3*T1*T2W
!$OMP Atomic         
	             Rho1(I1,J ,K1) =Rho1(I1,J ,K1) +D3*D1*T2W
!$OMP Atomic         
	             Rho1(I ,J1,K1) =Rho1(I ,J1,K1) +D3*T1*D2W
!$OMP Atomic         
                    Rho1(I1,J1,K1) =Rho1(I1,J1,K1) +D3*D1*D2W
            CASE (2)
!$OMP Atomic         
	             Rho2(I ,J ,K ) =Rho2(I ,J ,K ) +T3*T1*T2W
!$OMP Atomic         
                     Rho2(I1,J ,K ) =Rho2(I1,J ,K ) +T3*D1*T2W
!$OMP Atomic         
	             Rho2(I ,J1,K ) =Rho2(I ,J1,K ) +T3*T1*D2W
!$OMP Atomic         
	             Rho2(I1,J1,K ) =Rho2(I1,J1,K ) +T3*D1*D2W
!$OMP Atomic         
	             Rho2(I ,J ,K1) =Rho2(I ,J ,K1) +D3*T1*T2W
!$OMP Atomic         
	             Rho2(I1,J ,K1) =Rho2(I1,J ,K1) +D3*D1*T2W
!$OMP Atomic         
	             Rho2(I ,J1,K1) =Rho2(I ,J1,K1) +D3*T1*D2W
!$OMP Atomic         
                    Rho2(I1,J1,K1) =Rho2(I1,J1,K1) +D3*D1*D2W
            CASE (3)
!$OMP Atomic         
	             Rho3(I ,J ,K ) =Rho3(I ,J ,K ) +T3*T1*T2W
!$OMP Atomic         
                     Rho3(I1,J ,K ) =Rho3(I1,J ,K ) +T3*D1*T2W
!$OMP Atomic         
	             Rho3(I ,J1,K ) =Rho3(I ,J1,K ) +T3*T1*D2W
!$OMP Atomic         
	             Rho3(I1,J1,K ) =Rho3(I1,J1,K ) +T3*D1*D2W
!$OMP Atomic         
	             Rho3(I ,J ,K1) =Rho3(I ,J ,K1) +D3*T1*T2W
!$OMP Atomic         
	             Rho3(I1,J ,K1) =Rho3(I1,J ,K1) +D3*D1*T2W
!$OMP Atomic         
	             Rho3(I ,J1,K1) =Rho3(I ,J1,K1) +D3*T1*D2W
!$OMP Atomic         
                    Rho3(I1,J1,K1) =Rho3(I1,J1,K1) +D3*D1*D2W
            CASE (4)
!$OMP Atomic         
	             Rho4(I ,J ,K ) =Rho4(I ,J ,K ) +T3*T1*T2W
!$OMP Atomic         
                     Rho4(I1,J ,K ) =Rho4(I1,J ,K ) +T3*D1*T2W
!$OMP Atomic         
	             Rho4(I ,J1,K ) =Rho4(I ,J1,K ) +T3*T1*D2W
!$OMP Atomic         
	             Rho4(I1,J1,K ) =Rho4(I1,J1,K ) +T3*D1*D2W
!$OMP Atomic         
	             Rho4(I ,J ,K1) =Rho4(I ,J ,K1) +D3*T1*T2W
!$OMP Atomic         
	             Rho4(I1,J ,K1) =Rho4(I1,J ,K1) +D3*D1*T2W
!$OMP Atomic         
	             Rho4(I ,J1,K1) =Rho4(I ,J1,K1) +D3*T1*D2W
!$OMP Atomic         
                    Rho4(I1,J1,K1) =Rho4(I1,J1,K1) +D3*D1*D2W
            CASE (5)
!$OMP Atomic         
	             Rho5(I ,J ,K ) =Rho5(I ,J ,K ) +T3*T1*T2W
!$OMP Atomic         
                     Rho5(I1,J ,K ) =Rho5(I1,J ,K ) +T3*D1*T2W
!$OMP Atomic         
	             Rho5(I ,J1,K ) =Rho5(I ,J1,K ) +T3*T1*D2W
!$OMP Atomic         
	             Rho5(I1,J1,K ) =Rho5(I1,J1,K ) +T3*D1*D2W
!$OMP Atomic         
	             Rho5(I ,J ,K1) =Rho5(I ,J ,K1) +D3*T1*T2W
!$OMP Atomic         
	             Rho5(I1,J ,K1) =Rho5(I1,J ,K1) +D3*D1*T2W
!$OMP Atomic         
	             Rho5(I ,J1,K1) =Rho5(I ,J1,K1) +D3*T1*D2W
!$OMP Atomic         
                    Rho5(I1,J1,K1) =Rho5(I1,J1,K1) +D3*D1*D2W
            CASE (6)
!$OMP Atomic         
	             Rho6(I ,J ,K ) =Rho6(I ,J ,K ) +T3*T1*T2W
!$OMP Atomic         
                     Rho6(I1,J ,K ) =Rho6(I1,J ,K ) +T3*D1*T2W
!$OMP Atomic         
	             Rho6(I ,J1,K ) =Rho6(I ,J1,K ) +T3*T1*D2W
!$OMP Atomic         
	             Rho6(I1,J1,K ) =Rho6(I1,J1,K ) +T3*D1*D2W
!$OMP Atomic         
	             Rho6(I ,J ,K1) =Rho6(I ,J ,K1) +D3*T1*T2W
!$OMP Atomic         
	             Rho6(I1,J ,K1) =Rho6(I1,J ,K1) +D3*D1*T2W
!$OMP Atomic         
	             Rho6(I ,J1,K1) =Rho6(I ,J1,K1) +D3*T1*D2W
!$OMP Atomic         
                    Rho6(I1,J1,K1) =Rho6(I1,J1,K1) +D3*D1*D2W
              
                !------------------------------------------ NGP
              !      FI(I ,J ,K ) =FI(I ,J ,K ) + W
           end Select
        End Do
     End Do

   end subroutine DensMeshes

!---------------------------------------------
!
subroutine NormalizeMeshes
!
!--------------------------------------------
  use Parameters
character*120 :: Name

Real*8 ::    D1,D2,D3,T1,T2,T3,     &
             D4,D5,D6,T4,T5,T6,     &
             T1a,T2a,T3a,T4a,T5a,T6a
Integer*8 ::  ip

        NGRID = Ng
   
!				       Subtract mean density
           T1 =0. ; T2 =0.; T3=0.
           T4 =0. ; T5 =0.; T6=0.
           D1 =0. ; D2 =0.; D3=0.
           D4 =0. ; D5 =0.; D6=0.
           
!$OMP PARALLEL DO DEFAULT(SHARED) & 
!$OMP PRIVATE (M1,M2,M3) REDUCTION(+:T1,T2,T3,T4,T5,T6)
      DO M3=1,NGRID
       DO M2=1,NGRID
	 DO M1=1,NGRID
          T1 = T1 + Rho1(M1,M2,M3)
          T2 = T2 + Rho2(M1,M2,M3)
          T3 = T3 + Rho3(M1,M2,M3)
          T4 = T4 + Rho4(M1,M2,M3)
          T5 = T5 + Rho5(M1,M2,M3)
          T6 = T6 + Rho6(M1,M2,M3)
             
	 END DO
       END DO
    END DO
      T1a = T1/(float(Ngrid))**3
      T2a = T2/(float(Ngrid))**3
      T3a = T3/(float(Ngrid))**3
      T4a = T4/(float(Ngrid))**3
      T5a = T5/(float(Ngrid))**3
      T6a = T6/(float(Ngrid))**3
      T1 =0. ; T2 =0.; T3=0.
      T4 =0. ; T5 =0.; T6=0.

!$OMP PARALLEL DO DEFAULT(SHARED) & 
!$OMP PRIVATE (M1,M2,M3) REDUCTION(+:T1,T2,T3,T4,T5,T6,D1,D2,D3,D4,D5,D6)
      DO M3=1,NGRID
       DO M2=1,NGRID
	 DO M1=1,NGRID
          Rho1(M1,M2,M3) = Rho1(M1,M2,M3) -T1a
          Rho2(M1,M2,M3) = Rho2(M1,M2,M3) -T2a
          Rho3(M1,M2,M3) = Rho3(M1,M2,M3) -T3a
          Rho4(M1,M2,M3) = Rho4(M1,M2,M3) -T4a
          Rho5(M1,M2,M3) = Rho5(M1,M2,M3) -T5a
          Rho6(M1,M2,M3) = Rho6(M1,M2,M3) -T6a
          T1 = T1 + Rho1(M1,M2,M3)
          T2 = T2 + Rho2(M1,M2,M3)
          T3 = T3 + Rho3(M1,M2,M3)
          T4 = T4 + Rho4(M1,M2,M3)
          T5 = T5 + Rho5(M1,M2,M3)
          T6 = T6 + Rho6(M1,M2,M3)
          
          D1 = D1 + Rho1(M1,M2,M3)**2
          D2 = D2 + Rho2(M1,M2,M3)**2
          D3 = D3 + Rho3(M1,M2,M3)**2          
          D4 = D4 + Rho4(M1,M2,M3)**2
          D5 = D5 + Rho5(M1,M2,M3)**2
          D6 = D6 + Rho6(M1,M2,M3)**2          
	 END DO
       END DO
    END DO
      
      T1 = T1/(float(Ngrid))**3
      T2 = T2/(float(Ngrid))**3
      T3 = T3/(float(Ngrid))**3
      T4 = T4/(float(Ngrid))**3
      T5 = T5/(float(Ngrid))**3
      T6 = T6/(float(Ngrid))**3      
      D1 = sqrt(D1/(float(Ngrid))**3 -T1**2)          
      D2 = sqrt(D2/(float(Ngrid))**3 -T2**2)          
      D3 = sqrt(D3/(float(Ngrid))**3 -T3**2)          
      D4 = sqrt(D4/(float(Ngrid))**3 -T4**2)          
      D5 = sqrt(D5/(float(Ngrid))**3 -T5**2)          
      D6 = sqrt(D6/(float(Ngrid))**3 -T6**2)          

         write(*,'(a,i5)') ' Rho/Rms Meshes with Ngrid= ',Ngrid
         write(*,'(a,3es12.4)') ' 1: ',T1,D1,T1a
         write(*,'(a,3es12.4)') ' 2: ',T2,D2,T2a
         write(*,'(a,3es12.4)') ' 3: ',T3,D3,T3a
         write(*,'(a,3es12.4)') ' 4: ',T4,D4,T4a
         write(*,'(a,3es12.4)') ' 5: ',T5,D5,T5a
         write(*,'(a,3es12.4)') ' 6: ',T6,D6,Y6a

       end subroutine NormalizeMeshes
!----------------------------------------------------
!           power spectrum   P(k) for multi Mesh 6 lrvrls
!
!              
SUBROUTINE GetPowerMulti
!
!----------------------------------------------------
  use Parameters
  use Power
  Real*8        :: kNy,Alias
  Character*120 :: Name
  Real*8, parameter :: pi =3.1415926535

  NGRID = Ng
  Allocate(Pk(NGRID,0:6),dNharm(NGRID,0:6),kharm(NGRID,0:6))
  Allocate(FI(NGRID,NGRID,NGRID))
  Box   = Box2
  AEXPN = AEXPN2
    moment = max(INT(100.*(1./AEXPN-1.)+0.5),0)  ! make integer out of Z_moment
    Pk(:,:)     = 0.
    dNharm(:,:) = 0.
    kharm(:,:)  = 0.
    boxsize = Box
    kNy     = 2.*3.1415926/Box*NGRID/2.
       write(Name,'(2(a,i5.5),3(a,i3.3))')'PowerDM.all.',moment,'.dat'
       OPEN(18,FILE=TRIM(Name),STATUS='UNKNOWN')
       write(Name,'(2(a,i5.5),3(a,i3.3))')'PowerDM.',moment,'.dat'
       OPEN(20,FILE=TRIM(Name),STATUS='UNKNOWN')
       write(20,'(a,f7.4,a,f8.3,a,i4,a,f8.3)')' Aexpn =',AEXPN,    &
            ' Redshift= ',1./AEXPN-1.,' Ngrid= ',Ngrid,' kNyq= ',kNy
       FI = Rho1      !--- find Pk for the first mesh
       CALL POWERfft5(0,0)
       kharm(:,1)  = kharm(:,0)     !--- store results
       dNharm(:,1) = dNharm(:,0)
       Pk(:,1)     = Pk(:,0)
       
       FI = Rho2      !--- find Pk for the second mesh
       Box = Box/2
       CALL POWERfft5(0,0)
       kharm(:,2)  = kharm(:,0)     !--- store results
       dNharm(:,2) = dNharm(:,0)
       Pk(:,2)     = Pk(:,0)
       
       FI = Rho3      !--- find Pk for the third mesh
       Box = Box/2
       CALL POWERfft5(0,0)
       kharm(:,3)  = kharm(:,0)     !--- store results
       dNharm(:,3) = dNharm(:,0)
       Pk(:,3)     = Pk(:,0)
       
       FI = Rho4      !--- find Pk for the forth mesh
       Box = Box/2
       CALL POWERfft5(0,0)
       kharm(:,4)  = kharm(:,0)     !--- store results
       dNharm(:,4) = dNharm(:,0)
       Pk(:,4)     = Pk(:,0)
       
       FI = Rho5      !--- find Pk for the fifth mesh
       Box = Box/2
       CALL POWERfft5(0,0)
       kharm(:,5)  = kharm(:,0)     !--- store results
       dNharm(:,5) = dNharm(:,0)
       Pk(:,5)     = Pk(:,0)
       
       FI = Rho6      !--- find Pk for the six mesh
       Box = Box/2
       CALL POWERfft5(0,0)
       kharm(:,6)  = kharm(:,0)     !--- store results
       dNharm(:,6) = dNharm(:,0)
       Pk(:,6)     = Pk(:,0)


       
       
       WRITE (18,'(3x,a,T9,a,T19,a,T30,a,T43,a)') 'bin','k/Mpch',      &
                 'N','Preal'
       DO I=1,NGRID/2
         WRITE (18,'(i6,6(f9.5,es11.4,es13.5))')I,(kharm(I,j),dNharm(I,j),Pk(I,j),j=1,6)
      EndDo

       WRITE (20,'(3x,a,T15,a,T25,a,T35,a,T47,a)') 'k/Mpch',      &
                 'N','Pk','Alias','Mesh'

       N1   = NGRID/8
       Box  = boxsize
       exk  = 2.5
       kNy = 2.*pi/Box*NGRID/2.
       write(*,'(a,i4,2(a,es12.4))') 'N1 =',N1,' KNy= ',kNy,' Box= ',Box
       DO I=1,N1                        !--- first mesh: take first N1 harmonics
         Alias  = (1.-0.6667*sin(pi/2.*kharm(I,1)/kNy)**2) 
         WRITE (20,'(f9.5,es11.4,2es13.5,i3)')kharm(I,1),dNharm(I,1),Pk(I,1)/Alias,Alias,1
      EndDo

      ak = kharm(N1,1)
      write(*,'(a,i4,a,es12.4)') ' Number of harmonics for the first mesh =',N1,' last k= ',ak
      Do j=2,6                          !--- other meshes: take k = (ak,bk)
        bk = ak*exk
        kNy = kNy*2.
        write(*,'(a,i4,a,2es12.4)') ' Mesh number                          =',j,' first and last k= ',ak,bk
        DO I=1,NGRID/2
            Alias  = 1.-0.6667*sin(pi/2.*kharm(I,j)/kNy)**2
            If(kharm(I,j).gt.ak)Then
              WRITE (20,'(f9.5,es11.4,2es13.5,i3)')kharm(I,j),dNharm(I,j),Pk(I,j)/Alias,Alias,j
            end If
            If(kharm(I,j).gt.bk.or.I==NGRID/2)Then
               ak = kharm(I,j)
               write(*,'(a,i4,a,2es12.4)') '    I = ',I,' ak/bk= ',ak,bk
               exit
            end If
         EndDo
         
      end Do
      DEALLOCATE(FI)
      
end SUBROUTINE GetPowerMulti
!---------------------------------------
!        Read    PMfiles
!             moment <0    use PMcrd.DAT, PMcrs0.DAT ...    
!             moment >= 0  use PMcrd.xxxx.DAT, PMcrs0,XXXX.DAT ..
    SUBROUTINE ReadDataPM
!      
!---------------------------------------
use Parameters
       Character*120 :: Name
        Logical      :: exst
        Integer*8    :: iCount,ii,ioff,ip
        Integer*8    :: Nrecord,Jpage,NinPage
        Real*4       :: extras(100)
        Real*4, allocatable,  Dimension(:) ::  Xb,Yb,Zb,VXb,Vyb,Vzb
        Character*45 :: HEADER
!
        !			Read data and open files
        If(moment<0)Then
           write(Name,'(2a)')'PMcrd.DAT'
           write(*,*)Trim(Name)
         Open (4,file =Trim(Name),form ='UNFORMATTED',status ='UNKNOWN')
      Else
         write(Name,'(a,i4.4,a)')'PMcrd.',moment,'.DAT'
         Open (4,file =TRIM(Name),form ='UNFORMATTED',status ='UNKNOWN')
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
      Box2   = Box
      AEXPN2 = AEXPN
      NGRID_old = NGRID
      write(17,'(a,i10)') ' NROW   =',NROW
      write(17,'(a,i10)') ' Ngal   =',Nparticles
      write(17,'(a,i10)') ' Ngrid  =',Ngrid
      write(17,'(a,f10.1)') ' Box    =',Box
      Allocate (Xb(Nrecord),Yb(Nrecord),Zb(Nrecord))
      Allocate (VXb(Nrecord),VYb(Nrecord),VZb(Nrecord))

      Allocate (XPAR(Nparticles),YPAR(Nparticles),ZPAR(Nparticles))
      Allocate (VX(Nparticles),VY(Nparticles),VZ(Nparticles))

      iCount = 0

      ifile = 0
      jj    = 0
      If(moment<0)Then
         write(Name,'(a,i1.1,a)')'PMcrs',ifile,'.DAT'
      Else
         write(Name,'(a,i1.1,a,i4.4,a)')'PMcrs',ifile,'.',moment,'.DAT'
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
             write(Name,'(a,i1.1,a,i4.4,a)')'PMcrs',ifile,'.DAT'
           Else
             write(Name,'(a,i2.2,a,i4.4,a)')'PMcrs',ifile,'.DAT'
          EndIf
       Else
          If(ifile<10)Then
             write(Name,'(a,i1.1,a,i4.4,a)')'PMcrs',ifile,'.',moment,'.DAT'
          Else
             write(Name,'(a,i2.2,a,i4.4,a)')'PMcrs',ifile,'.',moment,'.DAT'
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
                  VX(ip+ioff) = VXb(ip) 
                  VY(ip+ioff) = VYb(ip) 
                  VZ(ip+ioff) = VZb(ip)
           end Do
        end DO
        close (20)

        Np = Nparticles
        xm = 1.e8 ; ym = 1.e8; zm = 1.e8
!$OMP PARALLEL DO DEFAULT(SHARED)  PRIVATE (ip) Reduction(max:xxm,yym,zzm) &
!$OMP  Reduction(min:xm,ym,zm)
        do ip =1,Nparticles
           xxm =max(xxm,Xpar(ip))
           yym =max(yym,Ypar(ip))
           zzm =max(zzm,Zpar(ip))
           xm =min(xm,Xpar(ip))
           ym =min(ym,Ypar(ip))
           zm =min(zm,Zpar(ip))
           
              If(Xpar(ip).lt.1.0.or.Ypar(ip).lt.1.0.or.Zpar(ip).lt.1.0)&
                write(17,'(a,i10,1p,3g14.5)')' Error coord: ',ip,Xpar(ip),Ypar(ip),Zpar(ip)
        EndDo
     write(17,'(a,3(2es14.5,3x))') ' Coordinates Min/Max =',xm,xxm,ym,yym,zm,zzm      
     DEALLOCATE (Xb,Yb,Zb,VXb,VYb,VZb)
     
   end SUBROUTINE ReadDataPM

