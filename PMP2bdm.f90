!----------------------------------------------------------------------------
!
!          BDM halo finder   A.Klypin 2010
!
!          PM files format
!          No subhalos
!
!----------------------------------------------------------------------------
 Program BDMrun
   use LinkerList
   Use Density

      Integer*8  ::  idummy,ip
      Integer*4  ::  jdummy
      integer*4  ::  OMP_GET_MAX_THREADS,OMP_GET_THREAD_NUM
      character*80 :: Path
      logical    ::  op
      t0 = seconds()

      iThreads = OMP_GET_MAX_THREADS()
      write (*,'(a,i4)')        ' Number of threads      =',iThreads
      Path =''
         write(*,'(a,$)')'  Enter snapshot number = '
         read(*,*)jStep

      write(*,'(a,i4)')  '  Go to ReadDataPM      =',jstep
      tstart = seconds()
      Call ReadDataPM(jstep,Path)             ! read PM data
      tfinish = seconds()
      write(13,'(10x,a,T50,2f10.2)') ' time for Reading  =',tfinish-tstart,tfinish-t0
      write(*,'(10x,a,T50,2f10.2)')  ' time for Reading  =',tfinish-tstart,tfinish-t0
      Np = Nparticles
      myMemory =Memory(1_8*NGRID*NGRID*NGRID)
      Allocate (FI(NGRID,NGRID,NGRID))

      tstart = seconds()

      Call BDM

 end Program BDMrun
