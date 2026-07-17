!--------------------------------------------------
!
!        Read PM data and make 3D overdensity grid
!        Calculate spectrum
!
!--------------------------------------------------
Program Spectrum
  use Tools
  use Density
  use Power
  !!! Real*4, ALLOCATABLE, DIMENSION(:,:,:) :: Fi1,Fi2,Fi3 ....
        character(len=80) :: Path,Fname
  
  Path = ''
  open(17,file='output')
  
   Call Initialize(Path)
   NGRID_old = NGRID                  ! store old value of NGRID
      WRITE (*,'(A,$)') ' Enter Mesh size for spectrum     => '
      READ (*,*) NGRID	 ! read snapshot number
   Mem_current = Memory(1_8*NGRID*NGRID*NGRID)
   ALLOCATE(FI(NGRID,NGRID,NGRID))

   write(Fname,'(a,i4.4,a,i4.4,a)') 'DensDistr.',ISTEP,'.',NGRID,'.dat'
   open(18,file=TRIM(Fname))
    
   write(*,*) ' NGRIDold = ',NGRID_old
   write(*,*) ' NGRIDnew = ',NGRID
   
      
      Call DENSITrsd(0,NGRID_old)         !- density in real space with new Ngrid
      Call DensDistr                      !  statistics of PDF
      Call GetPower(0)

   DEALLOCATE(FI)
   Mem_current = Memory(-1_8*NGRID*NGRID*NGRID)
            
end Program Spectrum

!--------------------------------------------
subroutine Initialize(Path)
!--------------------------------------------
use Tools

        logical :: exst
        character(len=*) :: Path

      WRITE (*,'(A,$)') ' Enter snapshot number to analyze => '
      READ (*,*) moment	         ! read snapshot number
      
      CALL ReadDataPM(moment,Path)

end subroutine Initialize

