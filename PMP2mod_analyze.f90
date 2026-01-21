!--------------------------------------------------
!
!      Routines to analyze results on PM simulations
!
!--------------------------------------------------
Module Analyze

CONTAINS

!--------------------------------------------
!
subroutine Analysis(Path,mDENSIT)
!
!--------------------------------------------
use Tools
use Density
use LUXURY
use Random
use Power
use LinkerList
        character*120 :: NameF
        character(len=*) :: Path


           !-----  Fill the table for things to analyze
    !iPower         = 1   !-- all particles Power spectrum 
    !iPowerRSD      = 0   !       redshift distortions for DM
    !iDensityDistr  = 1   !-- PDF for DM for different cell-sizes
    !iBias          = 0   !-- Biasing model: =1 for bias particles. =2 for density
    !iSave          = 0   !-- Save snapshot
    !iWriteMock     = 0   !-- Save mock catalog
    !iDumpPart      = 0   !-- Dump random fraction
    ! iBDM  = 0
    Call Timing(7,-1)      ! start reading time

    NGRID_old = NGRID                  ! store old value of NGRID
                                       ! mDENSIT DM density index:
                                       !  = 1 if density was overwritten
                                       !      and needs to be restored
    !moment = 100.*(1./AEXPN-1.)+0.5 ! redshift*100.
    !moment = max(moment,0)
    moment = ISTEP
    Nseed = 13979821
               !------------ Density Distribution in real space
    If(iDensityDistr == 1)Then
       write(NameF,'(2(a,i4.4),3(a,i3.3))')'DensDistrDM.',moment,   &
                                        '.',Nrealization,'.dat'
       OPEN(18,FILE=TRIM(NameF),STATUS='UNKNOWN')
       write(18,'(a)')HEADER
       write(18,'(a,f7.4,a,i4,a,f8.3,a,i4,a,f8.2)')' Aexpn =',AEXPN,' Step=',ISTEP,    &
            ' Redshift= ',1./AEXPN-1.,' Ngrid= ',Ngrid
       If(mDENSIT==1)Then
          CALL DENSIT
          mDENSIT = 0
       End If
      Call DensDistr
       Mem_current = Memory(-1_8*NGRID*NGRID*NGRID)
       Deallocate (FI)
       NGRID_old = NGRID                  ! store old value of NGRID
       NGRID = NGRID/2                    ! increase cell-size twice
       Mem_current = Memory(1_8*NGRID*NGRID*NGRID)
       ALLOCATE(FI(NGRID,NGRID,NGRID))
       mDENSIT = 1                        ! flip the density flag
      Call DENSITrsd(0,NGRID_old)         !- density in real space with new Ngrid
      Call DensDistr                      !  statistics of PDF
       Mem_current = Memory(-1_8*NGRID*NGRID*NGRID)
       Deallocate (FI)
       NGRID = NGRID_old/4
       Mem_current = Memory(1_8*NGRID*NGRID*NGRID)
       ALLOCATE(FI(NGRID,NGRID,NGRID))
      Call DENSITrsd(0,NGRID_old)
      Call DensDistr
       Mem_current = Memory(-1_8*NGRID*NGRID*NGRID)
       Deallocate (FI)
       NGRID = NGRID_old/8
       Mem_current = Memory(1_8*NGRID*NGRID*NGRID)
       ALLOCATE(FI(NGRID,NGRID,NGRID))
      Call DENSITrsd(0,NGRID_old)
      Call DensDistr
       Mem_current = Memory(-1_8*NGRID*NGRID*NGRID)
       Deallocate (FI)
       NGRID = NGRID_old/16
       Mem_current = Memory(1_8*NGRID*NGRID*NGRID)
       ALLOCATE(FI(NGRID,NGRID,NGRID))
      Call DENSITrsd(0,NGRID_old)
      Call DensDistr
      close(18)
                                         ! restore NGRID
       Mem_current = Memory(-1_8*NGRID*NGRID*NGRID)
       DEALLOCATE(FI)
      NGRID = NGRID_old
        Mem_current = Memory(1_8*NGRID*NGRID*NGRID)
       ALLOCATE(FI(NGRID,NGRID,NGRID))
    End If
                      !----------------- Power spectrum
    If(iPower == 1)Then
       If(mDENSIT==1)Then
          CALL DENSIT
          mDENSIT = 0
       EndIf
       If(iPowerRSD == 1)Then
         Mem_current = Memory(1_8*Nparticles)
         ALLOCATE(dens(Nparticles))
         CALL DENSPARTall     ! get density for each particle. Need it for corrections
      end If
       

      Call GetPower(0)  ! pk only real space, Ngrid mesh
      mDENSIT = 1
       If(iPowerRSD == 1)Then
         CALL SetRandom    ! set seeds for random numbers
         Call GetPower(1)  ! pk in real and redshift space, Ngrid/2 mesh
       end If
    End If
                      !----------------- Save snapshot    
    If(iDumpPart == 1)Then 
           If(mDENSIT==1)Then
             CALL DENSIT
             mDENSIT = 0
            EndIf
            Call ParticlesRandom
            Call ProjectDENSIT
    end If
                      !----------------- Biasing model  1  
    If(iBias == 1)Then
       If(mDENSIT==1)Then              ! restore DM density
          CALL DENSIT
          mDENSIT = 0
       EndIf       
       CALL DENSPART                   ! biasing scheme
       mDENSIT = 1                     ! flip flag because denspart corrupts FI
       Call BiasParticles  ! write galaxies into a file
         myMemory =Memory(-7_8*Ngalaxies)
         DeAllocate (Xb,Yb,Zb)
         DeAllocate (VXb,VYb,VZb)
         If(ALLOCATED(dens))Then
            myMemory =Memory(-2_8*Ngalaxies)
            DeAllocate (dens)
         End If
     end If      ! end Biasing model 1
                        !----------------- Biasing model  2
                        !          BiasPars(4) == celsize for Biasing model
    If(iBias == 2)Then
       Mem_current = Memory(-1_8*NGRID*NGRID*NGRID)
       DEALLOCATE(FI)             ! release FI to make room for other arrays      
                                  !--- deallocate arrays which are not needed
       mDENSIT = 1
       If(ALLOCATED(dens))Then
            myMemory =Memory(-2_8*Nparticles)
            DeAllocate (dens)
       End If
    If(ALLOCATED(SeedsPage))Then
            DeAllocate (SeedsPage)
         End If
    If(BiasPars(4)<1.e-5.or.BiasPars(4)>Box)Then
       write(*,*) ' Cell-size for BiasModel=2 is wrong: BiasPars(4)=',BiasPars(4)
       goto 10
    End If
    CellSize =  BiasPars(4)              ! new ngrid 
    NGRID = INT(Box/CellSize+0.5)
    If(.not.ALLOCATED(FI))Then 
       Mem_current = Memory(1_8*NGRID*NGRID*NGRID)
       ALLOCATE(FI(NGRID,NGRID,NGRID))
    end If
    Call DENSITrsd(0,NGRID_old)          ! new density field

    Call BiasDensity                     ! find number of galaxies in a cell
    !Call SelectGalaxies                 ! find particles using biased density

                                             ! restore NGRID and Density FI
       Mem_current = Memory(-1_8*NGRID*NGRID*NGRID)
       DEALLOCATE(FI)
       If(ALLOCATED(dens))Then
            myMemory =Memory(-2_8*Nparticles)
            DeAllocate (dens,dens2)
       End If
      NGRID = NGRID_old
        Mem_current = Memory(1_8*NGRID*NGRID*NGRID)
       ALLOCATE(FI(NGRID,NGRID,NGRID))
       mDENSIT = 1
 End If        !--- end BiasModel=2
 
10 continue
 If(iBDM==1)Call BDM(mDENSIT)
     Call Timing(7,1)      ! stop reading time

  end subroutine Analysis
!--------------------------------------------
!
!            create 3d mesh with number of galaxies in each node
!
subroutine BiasDensity
!
!           BiasPars(5) = density threshold
!                  (6) = normalization
!                  (7) = power law slope
!--------------------------------------------
use Tools
use LUXURY
use Random
character*120 :: Name
Real*4, allocatable :: Gg(:),SeedsP(:,:)
Integer*8 :: ip,jp,kp,jpp,kpp,IN,iGal

write(*,*) ' Inside BiasDensity   Allocate Gg and SeedsP Ngrid= ',NGRID
Call Timing(5,-1)
      ! moment = 100.*(1./AEXPN-1.) ! redshift*100.
      ! moment = max(moment,0)
      moment = ISTEP
      Allocate(Gg(NGRID))
      ALLOCATE(SeedsP(NGRID,NGRID))

      Ns     = Nseed      ! store initial seed
   Do j=1,NGRID 
      Do i=1,NGRID 
         SeedsP(i,j) = Ns
         dummy = RANDd(Ns) 
      End Do
   End Do
   Nseed = Ns    ! new seed
   !write(*,'(a,6i12)') ' Initialized random seeds: ',(SeedsPage(i,1),i=1,6)
      
      Gg = 0.
      iGal =0
      Do kp =1,NGRID
         if(mod(kp,50)==0)write(*,*) ' k=',kp
         kpp = (kp-1)*NGRID**2
      Do jp =1,NGRID
         jpp = (jp-1)*NGRID
            !Call getRandomG(Gg,jp,kp,NGRID)
            Call getRandom(Gg,jp,kp,NGRID)
       Do ip =1,NGRID
                        !---- Biasing model
         dd = FI(ip,jp,kp)+1.-BiasPars(1)
         If(dd>0.)Then             ! take only above density threshold
            
            ff = BiasPars(2)*dd**BiasPars(3)
            
            If(ff>1.)write(*,'(2(a,es13.4))') ' Too big bias: probability= ',ff, &
                    ' DM density= ',dens(IN)
         end If
      end Do
      end Do
      end Do

      Ngalaxies = iGal
      write(*,*)' Selected Galaxies=',Ngalaxies
      Call Timing(5,1)
      If(Ngalaxies == 0)Then
         write(*,*) '  No selected galaxies. Skip analysis.'
         deallocate(Gg,SeedsP)
         return
      End If
    end subroutine BiasDensity
!--------------------------------------------
!
!            randomly select a fraction of particles
!
subroutine ParticlesRandom
!
!                      
!--------------------------------------------
use Tools
use Random
character*120 :: Name
Integer*8 ::   ip,IN,iGal,ib,im,ish,iadd,ifr
Integer*8, parameter :: Nbuff = 100000
Real*4, dimension(Nbuff) :: Xbb,Ybb,Zbb,Vxbb,Vybb,Vzbb,DenP

write(*,*) ' Inside Particles Random :',NROW,NGRID
write(*,*) '              BiasPar 3  :',BiasPars(3)
   Nparticles = NROW*NROW*NROW
   Fraction = BiasPars(3)
   If(Fraction.le.0)Return
   IN = Fraction*1.1*Nparticles
   write(*,*) ' Inside ParticlesRandom. Fraction= ',Fraction
    write(*,*) '                      Nparticles= ',Nparticles
    write(*,*) '                      Expected  = ',IN
    ALLOCATE(Xb(IN),Yb(IN),Zb(IN))
    ALLOCATE(VXb(IN),VYb(IN),VZb(IN))
    ALLOCATE(dens(IN))

            ish=15485867; iadd = 1349973; ifr = 72421*109
    
   Call Timing(5,-1)
   !moment = 100.*(1./AEXPN-1.) ! redshift*100.
   !moment = max(moment,0)
   moment = ISTEP
      write(Name,'(2(a,i4.4),3(a,i3.3))')'Particles.',moment,'.',Nrealization,   &
           '.dat'
      OPEN(19,FILE=TRIM(Name),STATUS='UNKNOWN',form='unformatted')

      Xscale = Box/NGRID
      Vscale = Xscale*100./AEXPN    ! scale V to km/s
      Nseed = 10234909
      Ns     = Nseed      ! store initial seed
      t0 = seconds()
      iGal =0
      Do ip = 1,Nparticles
         Im  = mod(ifr*ip+iadd,ish)    !--- random numbers
         Rand = Im/float(ish)
         !if(Randd(Ns)>Fraction)cycle
         if(Rand>Fraction)cycle
         iGal = iGal +1
         if(iGal>IN)Stop 'Too many selected particles Increase IN buffer'
         Xb(iGal) = (XPAR(ip)-1.)*Xscale  !--- coordinates
         Yb(iGal) = (YPAR(ip)-1.)*Xscale
         Zb(iGal) = (ZPAR(ip)-1.)*Xscale
         VXb(iGal) = VX(ip)*Vscale        !--- velocity
         VYb(iGal) = VY(ip)*Vscale
         VZb(iGal) = VZ(ip)*Vscale
         
         X = XPAR(ip)     !--- find density
         Y = YPAR(ip)
         Z = ZPAR(ip)
         I = INT(X)
         J = INT(Y)
         K = INT(Z)
         I1 = I + 1
         J1 = J + 1
         K1 = K +1
         If(I1>NGRID)I1=1
         If(J1>NGRID)J1=1
         If(K1>NGRID)K1=1
           D1=X-FLOAT(I)
	   D2=Y-FLOAT(J)
	   D3=Z-FLOAT(K)
	   T1=1.-D1
	   T2=1.-D2
	   T3=1.-D3
	   I1=I+1
	      IF(I1.GT.NGRID)I1=1
	   J1=J+1
	      IF(J1.GT.NGRID)J1=1
	   K1=K+1
              IF(K1.GT.NGRID)K1=1
              dens(iGal)=               &
                FI(I ,J ,K )*T3*T1*T2 + &      ! density at particle position
                FI(I1,J ,K )*T3*D1*T2 + &
                FI(I ,J1,K )*T3*T1*D2 + &
                FI(I1,J1,K )*T3*D1*D2 + &
                FI(I ,J ,K1)*D3*T1*T2 + &
                FI(I1,J ,K1)*D3*D1*T2 + &
                FI(I ,J1,K1)*D3*T1*D2 + &
                FI(I1,J1,K1)*D3*D1*D2 + 1.
      end Do
      t1 = seconds()
      Ngalaxies = iGal
      write(*,*)' Selected Particles= ',Ngalaxies
      write(*,*)' time              = ',t1-t0
      If(Ngalaxies == 0)Then
         write(*,*) '  No selected galaxies. Skip analysis.'
      End If
      write(19) AEXPN,ISTEP,Fraction,Ngalaxies,Nbuff
      
      ib = 0
      Do ip =1,Ngalaxies
         ib = ib +1
         Xbb(ib) = Xb(ip)
         Ybb(ib) = Yb(ip)
         Zbb(ib) = Zb(ip)
         VXbb(ib) = VXb(ip)
         VYbb(ib) = VYb(ip)
         VZbb(ib) = VZb(ip)
         denP(ib) = dens(ip)
         If(ib.ge.Nbuff)Then
            write(19) Xbb,Ybb,Zbb,VXbb,VYbb,VZbb,denP
            ib = 0
         end If
      end Do
      If(ib/=0)write(19) Xbb,Ybb,Zbb,VXbb,VYbb,VZbb,denP
      
      deallocate(Xb,Yb,Zb)
      deallocate(VXb,VYb,VZb,dens)

      close(19)
      Call Timing(5,1)
    end subroutine ParticlesRandom
!--------------------------------------------
!
!            create a set of biased particles
!
subroutine BiasParticles
!
!           BiasPars(1) = density threshold
!                  (2) = normalization
!                  (3) = power law slope
!--------------------------------------------
use Tools
use LUXURY
use Random
character*120 :: Name
Real*4, allocatable :: Gg(:),SelGal(:)
Integer*8 :: ip,jp,kp,jpp,kpp,IN,iGal

write(*,*) ' Inside BiasParticles   Nrow= ',Nrow
      Call Timing(5,-1)

       !moment = 100.*(1./AEXPN-1.) ! redshift*100.
       !moment = max(moment,0)
       moment = ISTEP
      write(*,*)' Selected Galaxies=',Ngalaxies
      If(Ngalaxies == 0)Then
         write(*,*) '  No selected galaxies. Skip analysis.'
         return
      End If
      
      write(*,*) ' iWriteMock =',iWriteMock
      If(iWriteMock ==1)Then
       write(Name,'(2(a,i4.4),3(a,i3.3))')'GalaxiesZ.',moment,   &
                                        '.',Nrealization,'.dat'
       OPEN(18,FILE=TRIM(Name),STATUS='UNKNOWN')
       write(18,'(a)')HEADER
       write(18,'(a,i10,a,f8.3,a,es13.4,2(a,es12.3)),a,I5') &
            ' Ngalaxies=',Ngalaxies, &
            ' MaxDensThresh = ',BiasPars(1),    &
            ' PartThresh = ',BiasPars(2),   &
            ' Normalization = ',BiasPars(3), &
            ' Aexpn = ',AEXPN,' Step = ',ISTEP
       write(18,'(T10,a,T43,a,T67,a)')'XYZ(Mpch)','Vxyz(km/s)','DensDM rho/<rho>'
       Do i=1,Ngalaxies
          write(18,'(3f11.5,3x,3f10.3,es12.3)') Xb(i),Yb(i),Zb(i), &
                        VXb(i),VYb(i),VZb(i),dens(i)
       End Do
       close(18)
      end If

      Call Timing(5,1)
      
end subroutine BiasParticles


end Module Analyze
