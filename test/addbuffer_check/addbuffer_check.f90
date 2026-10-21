program addbuffer_check
  ! Cut 4 check: RetainParticles + AddBuffer(rescaleRetained=.true.) must give
  ! bit-for-bit the same buffered set as RescaleCoords(1) + AddBuffer, restore
  ! the simulation arrays exactly, and leave the memory counter where it was.
  use Tools
  use Structures
  use LinkerList
  implicit none
  integer*8 :: n,ip,npA,npB
  integer :: bad
  real, allocatable :: x0(:),y0(:),z0(:),vx0(:),vy0(:),vz0(:)
  real, allocatable :: xa(:),ya(:),za(:),vxa(:),vya(:),vza(:)
  integer*8, allocatable :: ida(:)
  real*8 :: state
  real :: m0,mA,mB

  n=300000_8; Box=256.; NGRID=1024; AEXPN=0.7; dBuffer=7.6875; Cell=3.5
  Nparticles=n; Np=n
  allocate(Xpar(n),Ypar(n),Zpar(n),VX(n),VY(n),VZ(n))
  state=0.123456789d0
  do ip=1,n
    Xpar(ip)=1.+real(NGRID)*real(lcg()); Ypar(ip)=1.+real(NGRID)*real(lcg()); Zpar(ip)=1.+real(NGRID)*real(lcg())
    VX(ip)=real((lcg()-0.5d0)*2000.d0); VY(ip)=real((lcg()-0.5d0)*2000.d0); VZ(ip)=real((lcg()-0.5d0)*2000.d0)
  enddo
  ! edge values: exact lower edge, values that round up to Box in real*4, far corners
  Xpar(1)=1.; Ypar(1)=1.; Zpar(1)=1.
  Xpar(2)=real(NGRID)+0.9999999; Ypar(2)=real(NGRID)+0.5; Zpar(2)=1.0000001
  Xpar(3)=real(NGRID)*0.999999+1.; Ypar(3)=1.+real(NGRID)*0.03; Zpar(3)=1.+real(NGRID)*0.97
  VX(1)=-0.; VY(2)=transfer(1,1.)
  allocate(x0(n),y0(n),z0(n),vx0(n),vy0(n),vz0(n))
  x0=Xpar; y0=Ypar; z0=Zpar; vx0=VX; vy0=VY; vz0=VZ
  m0=Memory(0_8)
  bad=0

  ! path A: the historical two-step path (rescaled copy, then the buffer)
  call RescaleCoords(1)
  call AddBuffer
  npA=Np
  allocate(xa(npA),ya(npA),za(npA),vxa(npA),vya(npA),vza(npA),ida(npA))
  xa=Xpar; ya=Ypar; za=Zpar; vxa=VX; vya=VY; vza=VZ; ida=OriginalParticleId
  call RemoveBuffer(n)
  mA=Memory(0_8)
  call check_restored('A')

  ! path B: retained workspace, buffer built straight from it
  call RetainParticles
  call AddBuffer(rescaleRetained=.true.)
  npB=Np
  if(npB/=npA)then; print *,'FAIL: buffered counts differ',npA,npB; bad=bad+1; endif
  if(npB==npA)then
    if(any(transfer(Xpar,[0],npB)/=transfer(xa,[0],npA)))then; print *,'FAIL: X bits differ'; bad=bad+1; endif
    if(any(transfer(Ypar,[0],npB)/=transfer(ya,[0],npA)))then; print *,'FAIL: Y bits differ'; bad=bad+1; endif
    if(any(transfer(Zpar,[0],npB)/=transfer(za,[0],npA)))then; print *,'FAIL: Z bits differ'; bad=bad+1; endif
    if(any(transfer(VX,[0],npB)/=transfer(vxa,[0],npA)))then; print *,'FAIL: VX bits differ'; bad=bad+1; endif
    if(any(transfer(VY,[0],npB)/=transfer(vya,[0],npA)))then; print *,'FAIL: VY bits differ'; bad=bad+1; endif
    if(any(transfer(VZ,[0],npB)/=transfer(vza,[0],npA)))then; print *,'FAIL: VZ bits differ'; bad=bad+1; endif
    if(any(OriginalParticleId/=ida))then; print *,'FAIL: ids differ'; bad=bad+1; endif
  endif
  call RemoveBuffer(n)
  mB=Memory(0_8)
  call check_restored('B')
  if(abs(mA-m0)>1e-6.or.abs(mB-m0)>1e-6)then; print *,'FAIL: memory counter',m0,mA,mB; bad=bad+1; endif
  print '(a,i0,a,i0,a,i0,a,3f10.6)',' originals ',n,' buffered A ',npA,' B ',npB,'  memory counter start/A/B',m0,mA,mB
  if(bad==0)then
    print *,'ADDBUFFER_CHECK PASS'
  else
    print *,'ADDBUFFER_CHECK FAIL',bad
    error stop 1
  endif
contains
  real*8 function lcg()
    state=mod(state*16807.d0,2147483647.d0)
    lcg=state/2147483647.d0
  end function lcg
  subroutine check_restored(tag)
    character(*), intent(in) :: tag
    if(any(transfer(Xpar,[0],n)/=transfer(x0,[0],n)).or.any(transfer(Ypar,[0],n)/=transfer(y0,[0],n)).or. &
       any(transfer(Zpar,[0],n)/=transfer(z0,[0],n)).or.any(transfer(VX,[0],n)/=transfer(vx0,[0],n)).or. &
       any(transfer(VY,[0],n)/=transfer(vy0,[0],n)).or.any(transfer(VZ,[0],n)/=transfer(vz0,[0],n)))then
      print *,'FAIL: simulation arrays not restored exactly, path ',tag; bad=bad+1
    endif
    if(allocated(OriginalParticleId).or.allocated(BdmPMX).or.allocated(BdmPMVX).or.BdmPMCount/=0_8.or.Np/=n.or.Nparticles/=n)then
      print *,'FAIL: workspace lifetime, path ',tag; bad=bad+1
    endif
  end subroutine check_restored
end program addbuffer_check
