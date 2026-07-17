  !----------------------------------------------
  !
  !  Read and re-format CatshortV catalog
  !
!----------------------------------------------
Module Structures
  real*8     :: C,delta_vir,r_vir, aMvir, Om0, r_s,aMhalo,Aexpn
  Integer*4  :: iSnap,iRealize
End module Structures
!----------------------------------------------    
Program Reformat
  
      use Structures
      Real*8 :: MassDelta,delta
      Character*120 :: File1,File2,Line

      Call Initialize
      nn = 0
      Do 
         read(1,*,iostat = ierr) x,y,z,vx,vy,vz,aMb,aMvir,Rvir,Vrms,Vcirc,Nhalo,C,aNp,iD,Xoff,te,aLambda
         if(ierr /= 0)exit
         aM500 = MassDelta(500.d0)
         aM200 = MassDelta(200.d0)
         write(10,'(3f10.4,3f10.2,4es12.4,f8.2,2f8.1,f9.3,3f9.4)')  &
              x,y,z,vx,vy,vz,aMb,aMvir,aM500,aM200,Rvir,Vrms,Vcirc,C,Xoff,aLambda
         nn = nn +1
      End Do
      close (1)
      close(10)
      write(*,*) ' Number of halos read = ',nn
  
end Program Reformat
!-------------------------------------------------
!
   Real*8 Function MassDelta(delta)
!
   use Structures
      real*8 :: x,FF,dx,xold,derivative,gold,delta,f
      FF = delta/delta_vir*f(C)/C**3
      xold = C
      !write(*,'(3(a,es12.4))') ' delta =',delta,' aMvir = ',aMvir,' delta_vir= ',delta_vir
      !write(*,*) ' x               xold      g_old       error derivative '
         derivative = (3.+4.*xold)/xold**3/(1.+xold)**2 - 3.*log(1.+xold)/xold**4
         gold  = f(xold)/xold**3
         x = xold +         (FF-gold)/derivative/3.
         xold = x
       Do i= 1,15
         derivative = (3.+4.*xold)/xold**3/(1.+xold)**2 - 3.*log(1.+xold)/xold**4
         gold  = f(xold)/xold**3
         x = xold +         (FF-gold)/derivative
         !write(*,'(10es12.4)') x,xold,gold,FF-gold,derivative
         xold = x
         if(abs(FF-gold)<1.e-5)exit
      end Do
    MassDelta = aMvir *(x/C)**3*delta/delta_vir
    end Function MassDelta
!
!------------------------
      Real*8 Function f(x)
      real*8 :: x
	f = log(1.+x)-x/(1.+x)
      end Function f
!-------------------------------------------------
!
   Subroutine  Initialize
!
     use Structures
     Character*120 :: File1,File2,Line,Str1*11,Str2*10,Str3*10,Str4*4,Str5*4,Line2*104,Line3*60
     Character :: FF*120
     logical  exst

      write(*,'(a,$)') ' Catshort name: '
      read(*,*)FF
      write(File1,'(a)')TRIM(FF)
      
   Inquire(file=TRIM(File1),exist = exst)
   if(.not.exst)Stop ' CatshortV not found. wrong file name'
   Open(1,file=TRIM(File1))

   write(File2,*)'HaloList',FF(9:23)
   Open(10,file=TRIM(ADJUSTL(File2)))

   Read(1,'(a)')Line
   Write(10,*)TRIM(Line)
   Read(1,*)Str4,Str5,Aexpn
   write(*,*) ' A=', Aexpn
   backspace(1)
   Do i=1,2
      Read(1,'(a)')Line
      Write(10,*)TRIM(Line)
   End Do
   Read(1,*)Str1,Om0
   write(*,*) TRIM(Str1),Om0
   backspace(1)
   Do i=1,3
     Read(1,'(a)')Line
     Write(10,'(a)')TRIM(Line)
  end Do
   
   Read(1,*)Str1,Str2,Str3  ,delta_vir0
    delta_vir = delta_vir0*Om0/(Om0+(1.-Om0)*AEXPN**3)    ! overdensity relative to the critical
   
   write(10,'(3a,16x,a,f10.4,a,f10.4)')'  ',TRIM(Str1),TRIM(Str2),TRIM(Str3),delta_vir0,' overdens/crit= ',delta_vir
   Read(1,'(2a)') Line2,Line3
   write(10,'(3x,a,22x,a,18x,4a)') 'XYZ(Mpch)','Vxyz(km/s)','Mbound      Mtot/Msunh', &
                               '   M500c','       M200c       Rvir    Vrms    Vcirc', &
                               '    Cvir     Xoff    Lambda' 
 end Subroutine Initialize





  
