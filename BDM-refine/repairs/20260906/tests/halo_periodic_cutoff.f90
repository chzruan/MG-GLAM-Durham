
program periodic_cutoff
use LinkerList
implicit none
integer*8,allocatable :: rows(:)
real*8,allocatable :: radii(:)
real*8 :: dx
integer :: q
Np=1;Nparticles=1; dBuffer=5.
allocate(Xpar(1),Ypar(1),Zpar(1),VX(1),VY(1),VZ(1))
Xpar=.1;Ypar=16.;Zpar=16.;VX=0.;VY=0.;VZ=0.
call AddBuffer
call SizeList
allocate(Lst(Np),Label(Nmx:Nbx,Nmy:Nby,Nmz:Nbz))
call List
call BdmHaloGather(31.9,16.,16.,.2d0,rows,radii)
dx=dble(Xpar(1))-dble(31.9)
dx=dx-dble(Box)*anint(dx/dble(Box))
print *, 'OBSERVED_GATHER_COUNT',size(rows)
print *, 'EXACT_ORIGINAL_MINIMUM_IMAGE_DISTANCE',abs(dx)
do q=1,size(rows)
 print *, 'GHOST_ROW_ID_DISTANCE',rows(q),OriginalParticleId(rows(q)),radii(q)
enddo
end program
