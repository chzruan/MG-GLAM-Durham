! Campaign-only change of PM coordinate units. No finder source is modified.
module ConvergenceGrid
  use Tools, only: NGRID,NROW,Nparticles,Box,AEXPN,Xpar,Ypar,Zpar,VX,VY,VZ
  use iso_fortran_env, only: int64,real32,real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
contains
  logical function GridPowerOfTwo(n)
    integer,intent(in) :: n
    GridPowerOfTwo=.false.
    if(n>0)GridPowerOfTwo=iand(n,n-1)==0
  end function

  pure real(real32) function GridPosition(value,source,target) result(mapped)
    real(real32),intent(in) :: value
    integer,intent(in) :: source,target
    mapped=real(1.d0+modulo(real(value,real64)-1.d0,real(source,real64))* &
      real(target,real64)/real(source,real64),real32)
    ! Rounding can land on the upper periodic image. Canonicalize that image.
    if(mapped>=real(target,real32)+1.)mapped=1.
  end function

  subroutine GridComponent(values,source,target,position,max_error,changed)
    real(real32),intent(inout) :: values(:)
    integer,intent(in) :: source,target
    logical,intent(in) :: position
    real(real64),intent(out) :: max_error
    integer(int64),intent(out) :: changed
    real(real64) :: ratio,old_physical,new_physical,error,scale
    real(real32) :: old,new
    integer(int64) :: ip
    logical :: invalid
    ratio=real(target,real64)/real(source,real64)
    scale=real(Box,real64)
    if(.not.position)scale=100.d0*scale/real(AEXPN,real64)
    invalid=.false.;changed=0_int64;max_error=0.d0
!$omp parallel do private(ip,old,new,old_physical,new_physical,error) &
!$omp reduction(.or.:invalid) reduction(+:changed) reduction(max:max_error)
    do ip=1,size(values,kind=int64)
      old=values(ip)
      if(.not.ieee_is_finite(old))then
        invalid=.true.;cycle
      endif
      if(position)then
        if(old<1..or.real(old,real64)>real(source,real64)+1.d0)then
          invalid=.true.;cycle
        endif
        new=GridPosition(old,source,target)
        old_physical=modulo(real(old,real64)-1.d0,real(source,real64))*scale/source
        new_physical=(real(new,real64)-1.d0)*scale/target
        error=abs(new_physical-old_physical)
        error=min(error,scale-error)
      else
        if(abs(real(old,real64)*ratio)>real(huge(old),real64))then
          invalid=.true.;cycle
        endif
        new=real(real(old,real64)*ratio,real32)
        error=abs(real(new,real64)*scale/target-real(old,real64)*scale/source)
      endif
      if(.not.ieee_is_finite(new))then
        invalid=.true.;cycle
      endif
      if(transfer(new,0)/=transfer(old,0))changed=changed+1_int64
      max_error=max(max_error,error)
      values(ip)=new
    enddo
    if(invalid)error stop 'Invalid or overflowing replay particle coordinate/velocity'
  end subroutine

  subroutine ChangeAnalysisGrid(target)
    integer,intent(in) :: target
    integer :: source
    integer(int64) :: changed(6)
    real(real64) :: error(6),position_bound
    source=NGRID
    if(.not.GridPowerOfTwo(source).or..not.GridPowerOfTwo(target)) &
      error stop 'Replay supports only power-of-two source and analysis grids'
    if(max(source,target)>16384)error stop 'Replay grid exceeds supported FFT workspace'
    if(NROW<1.or.NROW>=1200.or.Nparticles/=int(NROW,int64)**3) &
      error stop 'Replay requires a full PM particle set below 1200 cubed'
    if(.not.ieee_is_finite(Box).or..not.ieee_is_finite(AEXPN))error stop 'Invalid replay scales'
    if(Box<=0..or.AEXPN<=0.)error stop 'Nonpositive replay scales'
    if(.not.allocated(Xpar).or..not.allocated(Ypar).or..not.allocated(Zpar).or. &
       .not.allocated(VX).or..not.allocated(VY).or..not.allocated(VZ))error stop 'Missing replay particles'
    if(size(Xpar,kind=int64)/=Nparticles.or.size(Ypar,kind=int64)/=Nparticles.or. &
       size(Zpar,kind=int64)/=Nparticles.or.size(VX,kind=int64)/=Nparticles.or. &
       size(VY,kind=int64)/=Nparticles.or.size(VZ,kind=int64)/=Nparticles) &
      error stop 'Replay particle-array length mismatch'
    call GridComponent(Xpar,source,target,.true.,error(1),changed(1))
    call GridComponent(Ypar,source,target,.true.,error(2),changed(2))
    call GridComponent(Zpar,source,target,.true.,error(3),changed(3))
    call GridComponent(VX,source,target,.false.,error(4),changed(4))
    call GridComponent(VY,source,target,.false.,error(5),changed(5))
    call GridComponent(VZ,source,target,.false.,error(6),changed(6))
    position_bound=real(spacing(real(target,real32)+1.),real64)*real(Box,real64)/target
    if(maxval(error(:3))>position_bound)error stop 'Replay grid conversion exceeds coordinate rounding bound'
    NGRID=target
    write(*,'(a,2i8,a,i0)')'REPLAY GRID source/analysis=',source,target,' particles=',Nparticles
    write(*,'(a,6es24.15)')'REPLAY GRID physical errors xyz Mpc/h, velocity km/s=',error
    write(*,'(a,es24.15)')'REPLAY GRID position error bound Mpc/h=',position_bound
    write(*,'(a,6i20)')'REPLAY GRID changed storage values=',changed
  end subroutine
end module
