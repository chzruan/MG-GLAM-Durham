!--------------------------------------------------
!
!  Tune time-stepping to fit constraints on certain redshifts 
!    
!    Input three redshifts to match and limit number of steps
!
!--------------------------------------------------
Module Param
      Real*4   :: zout(100000)

    Contains


End Module Param
!
!--------------------------------------------------
Program Initialize
  use Param
      CHARACTER  ::    Name*120
      Real*8 :: ee,wN,dW
      Nda    = 1000
      Nainit = 1000
      write(*,'(a,$)') 'Enter three redshifts to match from largest to smallest : '
      read(*,*) z1,z2,z3
      ze3 =3.e-4
      ze2 =0.01
      ze1 =0.02
      Err = 1.e10
      write(*,'(a,$)') 'Enter optimum and spread in Number of steps  : '
      read(*,*) Nopt,Nwidth
      dW = Nwidth
      write(*,'(a,$)') 'Enter range of z_init from largest to smallest  : '
      read(*,*) z_init_large,z_init_small
       write(*,'(a,$)') 'Enter  da_init  : '
      read(*,*) da_init_0
      OPEN(10,file='Tune.dat')
      write(10,*) '---- Tuning time-stepping ----'
      write(10,'(a,3f9.3)') ' Redshifts to match     = ',z1,z2,z3
      write(10,'(a,3f9.3)') ' Initial Redshift range = ',z_init_large,z_init_small
      write(10,'(a,3f9.3)') ' Initial da guess       = ',da_init_0
      write(10,'(a,2i5)') ' Number of steps+range  = ',Nopt,Nwidth
      

      a_1 = 1./(1.+z_init_large)
      a_2 = 1./(1.+z_init_small)
      dda = (a_2-a_1)/Nainit
      Do i=0,Nainit
         a_init = a_1 + i*dda
         !write(*,*) ' a_init= ',a_init
         If(a_init>5.e-3.and.a_init.lt.0.02)Then
           Do j=-Nda/2,Nda/2
              da_init = da_init_0 + 2.e-6*j
              If(da_init>1.e-5.and.da_init<0.02)Then
              !write(*,*) '        da_init= ',da_init
                 Call TestStepping(a_init,da_init,Nn,0)
                 If(Nn> Nopt-Nwidth.and.Nn<Nopt+Nwidth)Then   
                      dz1 =1.e3 ; dz2 = 1.e3; dz3 =1.e3
                      Do k=2,Nn
                         zm = zout(k-1)
                         zp = zout(k)
                         If(zm>z1.and.zp.le.z1)dz1 = min(zm-z1,z1-zp)
                         If(zm>z2.and.zp.le.z2)dz2 = min(zm-z2,z2-zp)
                         If(zm>z3.and.zp.le.z3)dz3 = min(zm-z3,z3-zp)
                         !write(*,'(20x,2i4,3x,2es12.4,3x,3x,6es12.4)')Nn,k,zm,zp,z1,z2,z3,dz1,dz2,dz3
                      end Do
                      wN = 10.*exp((float(Nn-Nopt)/dW)**2)
                      ee = dz1/ze1 + dz2/ze2 + dz3/ze3 + wN   !--- goodness of fit
                      !write(*,'(f9.4,3x,es12.4,i5,3x,8es12.4 )') a_init,da_init,Nn, ee,Err,dz1,dz2,dz3
                      if(ee<Err)Then
                         Err = ee
                         a_init_1 = a_init
                         da_init_1 = da_init
                         write(*,'(a,es12.4,a,es12.4,a,i5,a,6es12.4)') &
                              'Zinit=',1./a_init_1-1,' da=',da_init_1,' N=',Nn,' Errs=',ee, dz1,dz2,dz3
                      end if
                   end If
                end If
             end Do
        end If
     end Do
     a_init  = a_init_1
     da_init = da_init_1
     Call TestStepping(a_init,da_init,Nn,1)
     
     
    end Program Initialize
!
!
!----------------------------------------------
SUBROUTINE TestStepping(a_init,da_init,Nn,iFlag)
  !----------------------------------------------
  use Param
  StepFactor = da_init/a_init
  zout(:) =1.e10
     a = a_init
         da= da_init
         i = 0
         If(iFlag ==1)Then
            write(*,'(a)') ' Step      z         a       da        da/a    step change'
         End If
         iCount = 0
         rMax   = da/a
         Nsteps2= 0
         Do 
            ind = 0
            If(da<StepFactor/1.25*a.and.a<0.5)Then
               If(iFlag ==1)write(*,*) '                  Increase step:',da,1.25*StepFactor*a
               da =1.5*da        ! increase step
               ind =1
               iCount = iCount +1
            EndIf
            a = a +da
            i = i +1
            zout(i) = 1./a-1.
            rmax = max(rMax,da/a)
            if(a>0.50)Nsteps2 = Nsteps2 +1
            If(iFlag ==1.and.a>0.15)write(*,'(i5,4es12.4,i3)') i,1./a-1.,a,da,da/a,ind
            if(a>1.)exit
         end Do
         Nn = i
         If(iFlag ==1)Then
            write(*,'(2(a,i5))')    ' Number of steps   = ',i, ' n_steps(z<2) =',Nsteps2 
            write(*,'(a,i5)')       ' Number of changes = ',iCount 
            write(*,'(a,f8.4)')     ' Maximum da/a      = ',rmax
            write(*,'(2(a,es12.5))') ' a_init    = ',a_init, ' z_init      = ',1./a_init-1.
            write(*,'(2(a,f11.5))') ' da_final  = ',da,     ' da/a_final  = ',da/a
            write(*,'(2(a,es12.5))') ' da_init   = ',da_init,' da/a_init = ',stepFactor
            write(10,'(2(a,i5))')    ' Number of steps   = ',i, ' n_steps(z<2) =',Nsteps2 
            write(10,'(a,i5)')       ' Number of changes = ',iCount 
            write(10,'(a,f8.4)')     ' Maximum da/a      = ',rmax
            write(10,'(2(a,es12.5))') ' a_init    = ',a_init, ' z_init      = ',1./a_init-1.
            write(10,'(2(a,f11.5))') ' da_final  = ',da,     ' da/a_final  = ',da/a
            write(10,'(2(a,es12.5))') ' da_init   = ',da_init,' da/a_init = ',stepFactor
         end If
end SUBROUTINE TestStepping

