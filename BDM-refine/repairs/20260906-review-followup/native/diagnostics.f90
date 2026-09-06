! Appended only to frozen diagnostic finder sources, after successful WriteFiles.
logical function RepairSelected(ip) result(selected)
  implicit none
  integer,intent(in) :: ip
  selected=.false.
  if(Mvir(ip)<MassMin.or.Mvir(ip)<=0.)return
  if(xMaxx(ip)<Xleft.or.xMaxx(ip)>=Xright.or.yMaxx(ip)<Yleft.or.yMaxx(ip)>=Yright.or. &
     zMaxx(ip)<Zleft.or.zMaxx(ip)>=Zright)return
  if(.not.allocated(BoundParticleIds(ip)%ids))error stop 'Replay selected membership is missing'
  if(size(BoundParticleIds(ip)%ids,kind=8)<20_8)return
  selected=.true.
end function RepairSelected

subroutine RepairDiagnosticPublish(source,destination)
  use iso_c_binding, only: c_int,c_char,c_null_char
  implicit none
  character(*),intent(in) :: source,destination
  character(kind=c_char,len=:),allocatable :: old_name,new_name
  integer(c_int) :: result
  logical :: exists
  interface
    integer(c_int) function c_rename(old_name,new_name) bind(C,name='rename')
      import c_int,c_char
      character(c_char),intent(in) :: old_name(*),new_name(*)
    end function c_rename
  end interface
  inquire(file=destination,exist=exists)
  if(exists)error stop 'Replay refuses to replace a diagnostic tape'
  old_name=source//c_null_char;new_name=destination//c_null_char
  result=c_rename(old_name,new_name)
  if(result/=0_c_int)error stop 'Replay diagnostic tape publication failed'
end subroutine RepairDiagnosticPublish

subroutine RepairDiagnosticsDump
  implicit none
  integer :: unit,ip,nselected,io
  integer*8 :: bound_count
  if(Nmaxima>0)then
    if(.not.allocated(BoundParticleIds).or..not.allocated(HaloStatus).or. &
       .not.allocated(HaloUnbindingPasses).or..not.allocated(HaloUnbindingWork)) &
      error stop 'Replay final diagnostic workspace is missing'
    if(size(BoundParticleIds)/=Nmaxima.or.size(HaloStatus)/=Nmaxima.or. &
       size(HaloUnbindingPasses)/=Nmaxima.or.size(HaloUnbindingWork)/=Nmaxima) &
      error stop 'Replay final diagnostic workspace size mismatch'
  endif
  nselected=0
  do ip=1,Nmaxima
    if(RepairSelected(ip))nselected=nselected+1
  enddo
  if(nselected/=Nhalo)error stop 'Replay membership selection differs from the published catalogue'
  open(newunit=unit,file='repair-members.bin.part',access='stream',form='unformatted', &
    status='new',action='write',iostat=io)
  if(io/=0)error stop 'Replay membership tape staging failed'
  write(unit,iostat=io)int(nselected,8),int(Nmaxima,8),MassOne
  if(io/=0)error stop 'Replay membership header write failed'
  do ip=1,Nmaxima
    if(.not.RepairSelected(ip))cycle
    write(unit,iostat=io)int(ip,8),size(BoundParticleIds(ip)%ids,kind=8)
    if(io/=0)error stop 'Replay membership row header write failed'
    write(unit,iostat=io)xMaxx(ip),yMaxx(ip),zMaxx(ip),VxMaxx(ip),VyMaxx(ip),VzMaxx(ip), &
      Mvir(ip),Mtotal(ip),Rvir(ip),EkinM(ip),EpotM(ip),VmaxM(ip),RmaxM(ip), &
      Xoff(ip),LambdaM(ip),RadRms(ip),Axba(ip),Axca(ip),Xax(ip),Yax(ip),Zax(ip)
    if(io/=0)error stop 'Replay membership properties write failed'
    write(unit,iostat=io)BoundParticleIds(ip)%ids
    if(io/=0)error stop 'Replay membership ID write failed'
  enddo
  close(unit,iostat=io)
  if(io/=0)error stop 'Replay membership close failed'
  call RepairDiagnosticPublish('repair-members.bin.part','repair-members.bin')
  open(newunit=unit,file='unbinding.bin.part',access='stream',form='unformatted', &
    status='new',action='write',iostat=io)
  if(io/=0)error stop 'Replay unbinding tape staging failed'
  write(unit,iostat=io)int(Nmaxima,8)
  if(io/=0)error stop 'Replay unbinding header write failed'
  do ip=1,Nmaxima
    bound_count=0_8
    if(allocated(BoundParticleIds(ip)%ids))bound_count=size(BoundParticleIds(ip)%ids,kind=8)
    write(unit,iostat=io)int(HaloUnbindingPasses(ip),4),int(HaloStatus(ip),4), &
      HaloUnbindingWork(ip),bound_count,Mvir(ip)
    if(io/=0)error stop 'Replay unbinding row write failed'
  enddo
  close(unit,iostat=io)
  if(io/=0)error stop 'Replay unbinding close failed'
  call RepairDiagnosticPublish('unbinding.bin.part','unbinding.bin')
  print '(a,i0,a,i0)','REPLAY DIAGNOSTICS candidates=',Nmaxima,' selected=',nselected
end subroutine RepairDiagnosticsDump
