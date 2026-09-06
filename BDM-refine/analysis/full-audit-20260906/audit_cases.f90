! Test drivers only. LinkerList routines are extracted from the audited source.
program audit_cases
  use LinkerList
  use, intrinsic :: ieee_arithmetic
  use omp_lib, only: omp_get_wtime
  implicit none
  character(80) :: which, arg
  integer :: q,h,j,k,n,a,b,c,repeat_count
  integer*8 :: ip,jp,original_n
  real :: origin, r, theta, zz, step, matvec(3), eig(3), direction(3)
  real :: saved(6), delta_r, vout
  real*8 :: tensor(3,3), tbegin, elapsed, norm, expected, ratio
  call get_command_argument(1,which)
  open(13,status='scratch'); open(18,status='scratch')
  select case(trim(which))
  case('peak_sparse','peak_empty','peak_dynamic_range','peak_plateau','peak_abacus')
    NGRID=16; Box=32.; AEXPN=1.; Om=.3; Nparticles=4096
    iVirial=2; Ovdens=-999.
    allocate(FI(NGRID,NGRID,NGRID)); FI=0.
    if(trim(which)=='peak_sparse') FI(8,8,8)=1000.
    if(trim(which)=='peak_dynamic_range')then
      FI(4,4,4)=1.e10; FI(12,12,12)=100.
    endif
    if(trim(which)=='peak_plateau') FI(4:13,4:13,4:13)=1000.
    if(trim(which)=='peak_abacus')then
      iVirial=3; FI(8,8,8)=1000.
    endif
    call FindMaxima
    print *, 'AUDIT PEAKS',Nmaxima,Ovdens
    open(50,file='peaks.bin',access='stream',form='unformatted',status='replace')
    write(50) xMaxx,yMaxx,zMaxx
    close(50)
  case('config_plain','config_invalid')
    call ReadParameters(1)
    print *, 'AUDIT CONFIG',iVirial,dLogR,MassMin
  case('shape_null_seed','shape_zero','shape_psd','shape_orthogonal')
    tensor=0.d0
    if(trim(which)=='shape_null_seed')then
      tensor(1,1)=1.d0; tensor(2,2)=1.d0
      tensor(1,2)=-1.d0; tensor(2,1)=-1.d0
    endif
    if(trim(which)=='shape_psd')then
      tensor(1,1)=3.d0; tensor(2,2)=2.d0; tensor(3,3)=1.d0
      tensor(1,2)=.2d0; tensor(2,1)=.2d0
    endif
    if(trim(which)=='shape_orthogonal')then
      tensor(1,1)=2.d0; tensor(2,2)=2.d0; tensor(3,3)=.5d0
      tensor(1,2)=-1.d0; tensor(2,1)=-1.d0
    endif
    call EigenValues(tensor,direction,eig)
    print *, 'AUDIT SHAPE',eig,direction
    print *, 'AUDIT PRINCIPAL_RESIDUAL',sqrt(sum((matmul(tensor,dble(direction))-eig(1)*direction)**2))
  case('concentration_nan')
    vout=Concentration(1.e13,500.,ieee_value(1.,ieee_quiet_nan))
    print *, 'AUDIT CONCENTRATION',vout
  case('buffer_face','buffer_corner','buffer_search')
    n=128; Np=n; Nparticles=n; original_n=Np; dBuffer=5.; Box=32.
    call particles(n)
    Xpar=0.; Ypar=16.; Zpar=16.; VX=0.; VY=0.; VZ=0.
    if(trim(which)=='buffer_corner')then
      Xpar=.01; Ypar=.01; Zpar=.01
    endif
    if(trim(which)=='buffer_search') Xpar=26.2
    call AddBuffer
    print *, 'AUDIT BUFFER',original_n,Np,count(Xpar==Box)
    if(trim(which)=='buffer_search')then
      Cell=1.; Nmx=-6;Nmy=-6;Nmz=-6;Nbx=38;Nby=38;Nbz=38
      allocate(Lst(Np),Label(Nmx:Nbx,Nmy:Nby,Nmz:Nbz)); call List
      n=0
      do jp=1,Np
        if((Xpar(jp)-.1)**2+(Ypar(jp)-16.)**2+(Zpar(jp)-16.)**2<36.) n=n+1
      enddo
      print *, 'AUDIT PERIODIC_FOUND_EXPECTED',n,original_n
    endif
  case('restore_roundtrip')
    n=1000; Np=n; Nparticles=n; Box=100.;NGRID=768
    call particles(n)
    do q=1,n
      Xpar(q)=1.+mod(q*37,767)+.12345
      Ypar(q)=Xpar(q); Zpar(q)=Xpar(q)
      VX(q)=q*.1234567; VY(q)=VX(q); VZ(q)=VX(q)
    enddo
    open(50,file='before.bin',access='stream',form='unformatted',status='replace')
    write(50) Xpar,Ypar,Zpar,VX,VY,VZ;close(50)
    call RescaleCoords(1)
    call RemoveBuffer(int(n,8))
    open(50,file='after.bin',access='stream',form='unformatted',status='replace')
    write(50) Xpar,Ypar,Zpar,VX,VY,VZ;close(50)
    print *, 'AUDIT RESTORED',Np
  case('halo_central','halo_compact','halo_onepass','halo_so','halo_disjoint')
    n=64
    if(trim(which)=='halo_onepass')n=1000
    if(trim(which)=='halo_so')n=100000
    if(trim(which)=='halo_disjoint')n=128
    Np=n;Nparticles=n;Nmaxima=1
    if(trim(which)=='halo_disjoint')Nmaxima=2
    call particles(n);call maxima(Nmaxima)
    MassOne=1.e10;Om0=.3;Ovdens=200.;dLogR=.02;Rext=0.;SlopeR=.2
    Box=32.;NGRID=128;Cell=.5
    Nmx=0;Nmy=0;Nmz=0;Nbx=32;Nby=32;Nbz=32
    xMaxx=5.;yMaxx=5.;zMaxx=5.;Xoff=100.
    step=.0125
    if(trim(which)=='halo_compact')step=.0025
    if(trim(which)=='halo_disjoint')then
      MassOne=1.e7;MassMin=1.e8;NGRID=4096;Cell=.05;step=.001
      xMaxx=[5.,5.15]
      Nmx=95;Nmy=95;Nmz=95;Nbx=110;Nby=105;Nbz=105
    endif
    if(n==64.or.trim(which)=='halo_disjoint')then
      q=0
      do h=1,Nmaxima
      do a=1,4
      do b=1,4
      do c=1,4
        q=q+1
        Xpar(q)=xMaxx(h)+(a-2.5)*step
        Ypar(q)=5.+(b-2.5)*step
        Zpar(q)=5.+(c-2.5)*step
      enddo
      enddo
      enddo
      enddo
    endif
    if(trim(which)=='halo_central')then
      Xpar(1)=5.;Ypar(1)=5.;Zpar(1)=5.
    endif
    if(trim(which)=='halo_so')then
      call get_command_argument(2,arg);read(arg,*)dLogR
      Cell=4.;Nmx=0;Nmy=0;Nmz=0;Nbx=8;Nby=8;Nbz=8
      MassOne=2.*(1.150e12*Om0)*Ovdens/n
    endif
    if(trim(which)=='halo_onepass'.or.trim(which)=='halo_so')then
      do q=1,n
        zz=1.-2.*(q-.5)/n
        theta=2.399963229728653*(q-1)
        r=.1
        if(trim(which)=='halo_so')r=2.*(q-.5)/n
        Xpar(q)=5.+r*sqrt(1.-zz*zz)*cos(theta)
        Ypar(q)=5.+r*sqrt(1.-zz*zz)*sin(theta)
        Zpar(q)=5.+r*zz
      enddo
    endif
    VX=0.;VY=0.;VZ=0.
    if(trim(which)=='halo_onepass')then
      do q=1,n
        VX(q)=(-1.)**q*2000.
        if(q<=100)VX(q)=(-1.)**q*500.
      enddo
    endif
    allocate(Lst(Np),Label(Nmx:Nbx,Nmy:Nby,Nmz:Nbz))
    allocate(AuditMembers(Np,Nmaxima));AuditMembers=.false.
    call List
    do ip=1,Nmaxima
      call GetHalo(xMaxx(ip),yMaxx(ip),zMaxx(ip),0.,0.,0.,ip)
      print *, 'AUDIT HALO',ip,Mvir(ip),Mtotal(ip),Rvir(ip),VmaxM(ip), &
          EkinM(ip),EpotM(ip),RadRms(ip),Axba(ip),Axca(ip)
    enddo
    if(trim(which)=='halo_so')then
      ratio=dble(Mtotal(1))/(dble(1.150e12*Om0)*Ovdens*dble(Rvir(1))**3)
      print *, 'AUDIT SO_RADIUS_DENSITY',Rvir(1),ratio
    endif
    if(trim(which)=='halo_onepass')then
      open(50,file='members.bin',access='stream',form='unformatted',status='replace')
      write(50) Xpar,Ypar,Zpar,VX,VY,VZ,AuditMembers(:,1)
      close(50)
      print *, 'AUDIT BOUND',count(AuditMembers(:,1))
    endif
    if(trim(which)=='halo_disjoint')then
      print *, 'AUDIT DISJOINT_BEFORE',count(Mvir>0.),count(AuditMembers(:,1)), &
          count(AuditMembers(:,2)),count(AuditMembers(:,1).and.AuditMembers(:,2))
      call ListMaxima
      call RemoveDuplicates
      print *, 'AUDIT DISJOINT_AFTER',count(Mvir>0.)
    endif
  case('centering_near','centering_far')
    n=262144;Np=n;Nparticles=n;Nmaxima=1;Box=1024.;Cell=.5;MassOne=1.e10
    origin=5.125
    if(trim(which)=='centering_far')origin=1000.125
    call particles(n);call maxima(1)
    do q=1,n
      Xpar(q)=origin+(mod(q-1,8)-3.5)/256.
      Ypar(q)=origin+(mod((q-1)/8,8)-3.5)/256.
      Zpar(q)=origin+(mod((q-1)/64,8)-3.5)/256.
    enddo
    VX=7.25;VY=0.;VZ=0.
    Nmx=int(origin/Cell)-4;Nbx=Nmx+8
    Nmy=Nmx;Nby=Nbx;Nmz=Nmx;Nbz=Nbx
    allocate(Lst(Np),Label(Nmx:Nbx,Nmy:Nby,Nmz:Nbz));call List
    xMaxx=origin;yMaxx=origin;zMaxx=origin;Xoff=100.
    call FindDistinctCandidates
    print *, 'AUDIT CENTRE_TRUE_MEASURED',origin,xMaxx,yMaxx,zMaxx
    print *, 'AUDIT CENTRE_DISPLACEMENT',sqrt((xMaxx(1)-origin)**2+(yMaxx(1)-origin)**2+(zMaxx(1)-origin)**2)
  case('list_benchmark')
    call get_command_argument(2,arg);read(arg,*)n
    Np=int(n,8)**3;Nparticles=Np;Box=64.;Cell=1.;NGRID=128
    Nmx=-1;Nmy=-1;Nmz=-1;Nbx=65;Nby=65;Nbz=65
    call particles(int(Np))
    do jp=1,Np
      Xpar(jp)=mod(jp-1,int(n,8))*Box/n+.001
      Ypar(jp)=mod((jp-1)/n,int(n,8))*Box/n+.001
      Zpar(jp)=((jp-1)/(int(n,8)*n))*Box/n+.001
    enddo
    allocate(Lst(Np),Label(Nmx:Nbx,Nmy:Nby,Nmz:Nbz))
    tbegin=omp_get_wtime()
    do repeat_count=1,10
      call List
    enddo
    elapsed=(omp_get_wtime()-tbegin)/10
    n=0
    do k=Nmz,Nbz
    do j=Nmy,Nby
    do q=Nmx,Nbx
      jp=Label(q,j,k)
      do while(jp/=0)
        n=n+1;jp=Lst(jp)
      enddo
    enddo
    enddo
    enddo
    print *, 'AUDIT LIST_SECONDS_COUNT',elapsed,n,Np
  case default
    error stop 'unknown audit case'
  end select
  print *, 'AUDIT REACHED_END'
contains
  subroutine particles(n)
    integer,intent(in)::n
    allocate(Xpar(n),Ypar(n),Zpar(n),VX(n),VY(n),VZ(n))
  end subroutine
  subroutine maxima(n)
    integer,intent(in)::n
    allocate(Mvir(n),Rvir(n),Mtotal(n),Xoff(n),xMaxx(n),yMaxx(n),zMaxx(n))
    allocate(VxMaxx(n),VyMaxx(n),VzMaxx(n),EpotM(n),EkinM(n),LambdaM(n))
    allocate(RadRms(n),VmaxM(n),RmaxM(n),Xax(n),Yax(n),Zax(n),Axba(n),Axca(n))
    Mvir=0.;Rvir=0.;Mtotal=0.;VxMaxx=0.;VyMaxx=0.;VzMaxx=0.
    EpotM=0.;EkinM=0.;LambdaM=0.;RadRms=0.;VmaxM=0.;RmaxM=0.
    Xax=0.;Yax=0.;Zax=0.;Axba=0.;Axca=0.
  end subroutine
end program
