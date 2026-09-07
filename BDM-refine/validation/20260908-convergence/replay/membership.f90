! Diagnostic-only hook after successful publication. Published haloes only.
logical function ConvergenceSelected(ip) result(selected)
  implicit none
  integer,intent(in) :: ip
  selected=.false.
  if(Mvir(ip)<MassMin.or.Mvir(ip)<=0.)return
  if(xMaxx(ip)<Xleft.or.xMaxx(ip)>=Xright.or.yMaxx(ip)<Yleft.or.yMaxx(ip)>=Yright.or. &
     zMaxx(ip)<Zleft.or.zMaxx(ip)>=Zright)return
  if(.not.allocated(BoundParticleIds(ip)%ids))error stop 'Published membership is missing'
  if(size(BoundParticleIds(ip)%ids,kind=8)<20_8)return
  selected=.true.
end function

subroutine ConvergenceDumpMembers
  use iso_c_binding, only: c_int,c_char,c_null_char
  implicit none
  integer :: unit,ip,nselected,io
  integer(c_int) :: rename_result
  logical :: exists
  interface
    integer(c_int) function conv_rename(old_name,new_name) bind(C,name='rename')
      import c_int,c_char
      character(c_char),intent(in) :: old_name(*),new_name(*)
    end function
  end interface
  inquire(file='repair-members.bin',exist=exists)
  if(exists)error stop 'Refuse to replace published membership tape'
  nselected=0
  do ip=1,Nmaxima
    if(ConvergenceSelected(ip))nselected=nselected+1
  enddo
  if(nselected/=Nhalo)error stop 'Membership selection differs from published catalogue'
  open(newunit=unit,file='repair-members.bin.part',access='stream',form='unformatted', &
    status='new',action='write',iostat=io)
  if(io/=0)error stop 'Cannot stage membership tape'
  write(unit,iostat=io)int(nselected,8),int(Nmaxima,8),MassOne
  if(io/=0)error stop 'Membership header write failed'
  do ip=1,Nmaxima
    if(.not.ConvergenceSelected(ip))cycle
    write(unit,iostat=io)int(ip,8),size(BoundParticleIds(ip)%ids,kind=8)
    if(io/=0)error stop 'Membership row header write failed'
    write(unit,iostat=io)xMaxx(ip),yMaxx(ip),zMaxx(ip),VxMaxx(ip),VyMaxx(ip),VzMaxx(ip), &
      Mvir(ip),Mtotal(ip),Rvir(ip),EkinM(ip),EpotM(ip),VmaxM(ip),RmaxM(ip), &
      Xoff(ip),LambdaM(ip),RadRms(ip),Axba(ip),Axca(ip),Xax(ip),Yax(ip),Zax(ip)
    if(io/=0)error stop 'Membership property write failed'
    write(unit,iostat=io)BoundParticleIds(ip)%ids
    if(io/=0)error stop 'Membership ID write failed'
  enddo
  close(unit,iostat=io)
  if(io/=0)error stop 'Membership tape close failed'
  rename_result=conv_rename('repair-members.bin.part'//c_null_char,'repair-members.bin'//c_null_char)
  if(rename_result/=0_c_int)error stop 'Membership publication failed'
  print '(a,i0,a,i0)','REPLAY MEMBERSHIP candidates=',Nmaxima,' selected=',nselected
end subroutine
