module input_files
  character*1024 :: eqfile, cfile, gfile,pfile,convexfile,fluxdatapath
  integer :: iunit=1738
!
  data eqfile  /'d3d.equ'/
  data cfile   /'DATA/ccoil.dat'/
!  data gfile   /'gfiles/shot115452/g115452.03525'/
!  data pfile   /'Conly/probe_g129_bfield.out'/
end module input_files
!
module theta_rz_mod
  integer :: icall=0
  integer :: nsqp,nlab,nthe,icp_pt
  integer, dimension(:,:), allocatable :: ipoint_pt
  real(kind=8) :: hsqpsi,hlabel,htheqt,psiaxis,sigma_qt,raxis,zaxis
  real(kind=8), dimension(:,:),   allocatable :: spllabel
  real(kind=8), dimension(:,:,:), allocatable :: splthet
  real(kind=8), dimension(:),     allocatable :: sqpsi,flab,theqt
end module theta_rz_mod
!
module extract_fluxcoord_mod
  integer :: load_extract_fluxcoord=1
  integer :: nphinorm
  double precision :: psif_extract,theta_extract,psifmin,hpsif
  double precision :: psifmax,phifmax,sigcos
  double precision, dimension(:), allocatable :: phinorm_arr
end module extract_fluxcoord_mod
!
module field_eq_mod
  logical :: use_fpol = .true.                                      !<=18.12.18
  integer :: icall_eq=0
  integer :: nrad,nzet,icp,nwindow_r,nwindow_z
  real(kind=8), parameter                      :: pi=3.14159265358979d0
  real(kind=8) :: psib,btf,rtf,hrad,hzet
  real(kind=8) :: psi_axis,psi_sep,hfpol                            !<=18.12.18
  real(kind=8), dimension(:,:), allocatable    :: psi, psi0
  real(kind=8), dimension(:,:), allocatable    :: splfpol           !<=18.12.18
  real(kind=8), dimension(:,:,:), allocatable  :: splpsi
  real(kind=8), dimension(:), allocatable      :: rad, zet, xi,f
  integer, dimension(:), allocatable           :: imi,ima,jmi,jma
  integer, dimension(:,:), allocatable         :: ipoint
  double precision :: psif,dpsidr,dpsidz,d2psidr2,d2psidrdz,d2psidz2
end module field_eq_mod

module utils_bdivfree_mod

  implicit none

  contains

subroutine indef_bdf(u,umin,dum1,nup,indu)
! defines interval for 1D interpolation on uniform mesh, normally 
! looks for the central interval of stencil, but
! stops moving of stencil at the boundary (works for mp=4 only!)
! Input:
!    u - coordinate of a point (to be interpolated)
!    umin - minimal value of u
!    dum1 = 1./h reciprocal to length of mesh interval
!    nup - total number of mesh points
! Output:
!    indu(mp) - relative index of stencil points
!
! the power 3 of polinomial is fixed strictly:
!
    integer :: nup, i
    integer, parameter :: mp = 4
    integer, dimension(mp) :: indu
    double precision :: u, umin, dum1
                            
    indu(1) = int((u-umin)*dum1)
    if( indu(1) .le. 0 ) indu(1) = 1
    indu(mp) = indu(1) + mp - 1
    if( indu(mp) .gt. nup ) then
        indu(mp) = nup
        indu(1) = indu(mp) - mp + 1
    endif
    do i=2,mp-1
        indu(i) = indu(i-1) + 1
    enddo

    return 
end subroutine indef_bdf
!
!
subroutine plag2d_bdf(x,y,fp,dxm1,dym1,xp,yp,polyl2d)
!
! 2D interpolation by means of Lagrange polynomial
! the power 3 is fixed strictly:
  integer, parameter :: mp = 4
! uniform mesh (increasingly ordered) in all dimensions is implied
!
! Input parameters:
!   x,y - coordinates of the point for interpolation
!   dxm1,dym1 - inverse steps in each direction
!   xp,yp - vertices of stencil
!
! Output parameters:
! polyl2d - polynomial itself
  double precision :: x,y,dxm1,dym1,polyl2d
  double precision, dimension(mp) :: cx,cy,xp,yp
  double precision, dimension(mp,mp) :: fp
  integer :: i, j
!
  call coefs_bdf(x,xp,dxm1,cx)
  call coefs_bdf(y,yp,dym1,cy)
!
  polyl2d = 0.d0
  do j=1,mp
    do i=1,mp
      polyl2d = polyl2d + fp(i,j)*cx(i)*cy(j)
    enddo
  enddo
!
  return
end subroutine plag2d_bdf
!
subroutine coefs_bdf(u,up,dum1,cu)
!
    integer, parameter :: mp = 4
    double precision :: u, dum1, one6, du3
    double precision, dimension(mp) :: up, cu
    data one6/0.16666666666667d0/
    du3 = dum1**3
    cu(1) = (u - up(2)) * (u - up(3)) * (u - up(4)) * (-one6*du3)
    cu(2) = (u - up(1)) * (u - up(3)) * (u - up(4)) * (0.5d0*du3)
    cu(3) = (u - up(1)) * (u - up(2)) * (u - up(4)) * (-0.5d0*du3)
    cu(4) = (u - up(1)) * (u - up(2)) * (u - up(3)) * (one6*du3)
    return
end subroutine coefs_bdf
!
subroutine spl_five_per(n,h,a,b,c,d,e,f)
!
! Periodic spline of the 5-th order. First and last values of function must
! be the same.
!
  integer :: n,i,ip1,ip2
  double precision :: h,rhop,rhom,fac,xplu,xmin,gammao_m,gammao_p
  double precision :: c_gammao_m,c_gammao_p
  double precision, dimension(n) :: a,b,c,d,e,f
  double precision, dimension(:), allocatable :: alp,bet,gam
!
  rhop=13.d0+sqrt(105.d0)
  rhom=13.d0-sqrt(105.d0)
!
  allocate(alp(n),bet(n),gam(n))
!
  alp(1)=0.0d0
  bet(1)=0.0d0
!
  do i=1,n-4
    ip1=i+1
    alp(ip1)=-1.d0/(rhop+alp(i))
    bet(ip1)=alp(ip1)*(bet(i)- &
             5.d0*(a(i+4)-4.d0*a(i+3)+6.d0*a(i+2)-4.d0*a(ip1)+a(i)))
  enddo
  alp(n-2)=-1.d0/(rhop+alp(n-3))
  bet(n-2)=alp(n-2)*(bet(n-3)- &
           5.d0*(a(2)-4.d0*a(1)+6.d0*a(n-1)-4.d0*a(n-2)+a(n-3)))
  alp(n-1)=-1.d0/(rhop+alp(n-2))
  bet(n-1)=alp(n-1)*(bet(n-2)- &
           5.d0*(a(3)-4.d0*a(2)+6.d0*a(1)-4.d0*a(n-1)+a(n-2)))
  alp(n)=-1.d0/(rhop+alp(n-1))
  bet(n)=alp(n)*(bet(n-1)- &
           5.d0*(a(4)-4.d0*a(3)+6.d0*a(2)-4.d0*a(1)+a(n-1)))
!
  gam(n)=bet(n)
  do i=n-1,1,-1
    gam(i)=gam(i+1)*alp(i)+bet(i)
  enddo
!
  xplu=sqrt(0.25d0*rhop**2-1.d0)-0.5d0*rhop
  xmin=-sqrt(0.25d0*rhop**2-1.d0)-0.5d0*rhop
  gammao_m=(gam(2)+xplu*gam(n))/(xmin-xplu)
  gammao_p=(gam(2)+xmin*gam(n))/(xplu-xmin)
  if(abs(xmin).lt.1) then
    c_gammao_m=gammao_m/(xmin**(n-1)-1.d0)
  else
    c_gammao_m=gammao_m*(1.d0/xmin)**(n-1)/(1.d0-(1.d0/xmin)**(n-1))
  endif
  if(abs(xplu).lt.1) then
    c_gammao_p=gammao_p/(xplu**(n-1)-1.d0)
  else
    c_gammao_p=gammao_p*(1.d0/xplu)**(n-1)/(1.d0-(1.d0/xplu)**(n-1))
  endif
  gam(1)=gam(1)+c_gammao_m+c_gammao_p
  do i=2,n
    if(abs(xmin).lt.1) then
      c_gammao_m=gammao_m*xmin**(i-1)/(xmin**(n-1)-1.d0)
    else
      c_gammao_m=gammao_m*(1.d0/xmin)**(n-i)/(1.d0-(1.d0/xmin)**(n-1))
    endif
    if(abs(xplu).lt.1) then
      c_gammao_p=gammao_p*xplu**(i-1)/(xplu**(n-1)-1.d0)
    else
      c_gammao_p=gammao_p*(1.d0/xplu)**(n-i)/(1.d0-(1.d0/xplu)**(n-1))
    endif
    gam(i)=gam(i)+c_gammao_m+c_gammao_p
  enddo
!
  alp(1)=0.0d0
  bet(1)=0.d0
!
  do i=1,n-1
    ip1=i+1
    alp(ip1)=-1.d0/(rhom+alp(i))
    bet(ip1)=alp(ip1)*(bet(i)-gam(i))
  enddo
!
  e(n)=bet(n)
  do i=n-1,1,-1
    e(i)=e(i+1)*alp(i)+bet(i)
  enddo
!
  xplu=sqrt(0.25d0*rhom**2-1.d0)-0.5d0*rhom
  xmin=-sqrt(0.25d0*rhom**2-1.d0)-0.5d0*rhom
  gammao_m=(e(2)+xplu*e(n))/(xmin-xplu)
  gammao_p=(e(2)+xmin*e(n))/(xplu-xmin)
  if(abs(xmin).lt.1) then
    c_gammao_m=gammao_m/(xmin**(n-1)-1.d0)
  else
    c_gammao_m=gammao_m*(1.d0/xmin)**(n-1)/(1.d0-(1.d0/xmin)**(n-1))
  endif
  if(abs(xplu).lt.1) then
    c_gammao_p=gammao_p/(xplu**(n-1)-1.d0)
  else
    c_gammao_p=gammao_p*(1.d0/xplu)**(n-1)/(1.d0-(1.d0/xplu)**(n-1))
  endif
  e(1)=e(1)+c_gammao_m+c_gammao_p
  do i=2,n
    if(abs(xmin).lt.1) then
      c_gammao_m=gammao_m*xmin**(i-1)/(xmin**(n-1)-1.d0)
    else
      c_gammao_m=gammao_m*(1.d0/xmin)**(n-i)/(1.d0-(1.d0/xmin)**(n-1))
    endif
    if(abs(xplu).lt.1) then
      c_gammao_p=gammao_p*xplu**(i-1)/(xplu**(n-1)-1.d0)
    else
      c_gammao_p=gammao_p*(1.d0/xplu)**(n-i)/(1.d0-(1.d0/xplu)**(n-1))
    endif
    e(i)=e(i)+c_gammao_m+c_gammao_p
  enddo
!
  do i=n-1,1,-1
    f(i)=(e(i+1)-e(i))/5.d0
  enddo
  f(n)=f(1)
!
  d(n-1)=(a(3)-3.d0*a(2)+3.d0*a(1)-a(n-1))/6.d0 &
      -(e(3)+27.d0*e(2)+93.d0*e(1)+59.d0*e(n-1))/30.d0
  d(n-2)=(a(2)-3.d0*a(1)+3.d0*a(n-1)-a(n-2))/6.d0 &
      -(e(2)+27.d0*e(1)+93.d0*e(n-1)+59.d0*e(n-2))/30.d0
  do i=n-3,1,-1
    d(i)=(a(i+3)-3.d0*a(i+2)+3.d0*a(i+1)-a(i))/6.d0 &
        -(e(i+3)+27.d0*e(i+2)+93.d0*e(i+1)+59.d0*e(i))/30.d0
  enddo
  d(n)=d(1)
  c(n-1)=0.5d0*(a(2)+a(n-1))-a(1)-0.5d0*d(1)-2.5d0*d(n-1) &
      -0.1d0*(e(2)+18.d0*e(1)+31.d0*e(n-1))
  b(n-1)=a(1)-a(n-1)-c(n-1)-d(n-1)-0.2d0*(4.d0*e(n-1)+e(1))
!
  do i=n-2,1,-1
    c(i)=0.5d0*(a(i+2)+a(i))-a(i+1)-0.5d0*d(i+1)-2.5d0*d(i) &
        -0.1d0*(e(i+2)+18.d0*e(i+1)+31.d0*e(i))
    b(i)=a(i+1)-a(i)-c(i)-d(i)-0.2d0*(4.d0*e(i)+e(i+1))
  enddo
  b(n)=b(1)
  c(n)=c(1)
!
  fac=1.d0/h
  b=b*fac
  fac=fac/h
  c=c*fac
  fac=fac/h
  d=d*fac
  fac=fac/h
  e=e*fac
  fac=fac/h
  f=f*fac
!
  deallocate(alp,bet,gam)
!
  return
end subroutine spl_five_per
!
subroutine s2dring(nx,ny,hx,hy,f,icount,spl,ipoint)
!
! Calculates coefficients of a 2D spline for a ring domain
! (periodic over y variable)
! equidistant mesh, but hx must not be = hy
!
!  Input parameters:
!                    nx           - horizontal size of the mesh (over x)
!                    ny           - vertical size of the mesh (over y)
!                    hx           - step of the mesh over x
!                    hy           - step of the mesh over y
!                    f(i,j)       - array of values to be interpolated
!                                   (i = 1, ..., nx; j = 1, ..., ny).
!
!                                   For the case of non-rectangular domain:
!                                   numbers of the mesh points which
!                                   correspond to the boundaries of the
!                                   interpolation region:
!                    icount       - maximum number of entries in spl
! Output parameters:
!                    spl(l,m,k)   - spline coefficients (i,j = 1, ... , n;
!                    ipoint(i,j)    l,m = 1, ..., 4; i,j - numbers of the
!                                   mesh point in horizontal and vertical
!                                   direction (over x and over y), l,m -
!                                   the numbers of expansion power over x
!                                   and y (~ dx**(l-1)*dy**(m-1) ))
!                                   ipoint(i,j) contains the pointer to k
!
    use spline5_RZ_mod, only: spl_five_reg
! 
  double precision, dimension(nx,ny)          :: f
  double precision, dimension(6,6,icount)     :: spl(6,6,icount)
  double precision, dimension(:), allocatable :: ai,bi,ci,di,ei,fi
  double precision                            :: hx, hy
  integer,          dimension(nx,ny)          :: ipoint(nx,ny)
  integer,          dimension(:), allocatable :: imi,ima,jmi,jma
  integer                                     :: nx,ny,icount, i, ic, ii, j, jj, l, nmax, nsi
! 
  nmax=max(nx,ny)
! 
  allocate( ai(nmax),bi(nmax),ci(nmax),di(nmax),ei(nmax),fi(nmax) )
  allocate(imi(ny),ima(ny),jmi(nx),jma(nx))
!
  imi=1
  ima=nx
  jmi=1
  jma=ny
! 
  spl=0.d0
  ipoint=-1
!
!  spline along Y-axis
!
  ic = 0
  do i=1,nx
    if(jmi(i).gt.0) then
      nsi=jma(i)-jmi(i)+1
      do j=jmi(i),jma(i)
        ai(j-jmi(i)+1)=f(i,j)
      enddo
      call spl_five_per(nsi,hy,ai,bi,ci,di,ei,fi)
      do j=jmi(i),jma(i)
        jj=j-jmi(i)+1
        ic = ic+1
        ipoint(i,j)=ic
        spl(1,1,ic)=ai(jj)
        spl(1,2,ic)=bi(jj)
        spl(1,3,ic)=ci(jj)
        spl(1,4,ic)=di(jj)
        spl(1,5,ic)=ei(jj)
        spl(1,6,ic)=fi(jj)
      enddo
    endif
  enddo
!
  if (ic .ne. icount) then
    write (6,*) 'Warning, ic, icount:  ',ic,icount
  endif
!
!  spline along X-axis
!
  do j=1,ny
    if(imi(j).gt.0) then
      nsi=ima(j)-imi(j)+1
      do l=1,6
        do i=imi(j),ima(j)
          ai(i-imi(j)+1)=spl(1,l,ipoint(i,j))
        enddo
        call spl_five_reg(nsi,hx,ai,bi,ci,di,ei,fi)
        do i=imi(j),ima(j)
          ii=i-imi(j)+1
          spl(2,l,ipoint(i,j))=bi(ii)
          spl(3,l,ipoint(i,j))=ci(ii)
          spl(4,l,ipoint(i,j))=di(ii)
          spl(5,l,ipoint(i,j))=ei(ii)
          spl(6,l,ipoint(i,j))=fi(ii)
        enddo
      enddo
    endif
  enddo
!
  deallocate( ai,bi,ci,di,ei,fi )
!
  return
end subroutine s2dring
!
subroutine load_theta
!
  use theta_rz_mod
  use input_files, only : iunit,fluxdatapath
  use spline5_RZ_mod, only: spl_five_reg
!
  double precision, parameter :: pi=3.14159265358979d0
!
  integer :: nsqpsi,nlabel,ntheta,i
  real(kind=8) :: sqpsimin,sqpsimax
  real(kind=8) :: flabel_min,flabel_max
!
  real(kind=8), dimension(:),   allocatable :: flabel
  real(kind=8), dimension(:,:), allocatable :: theta_of_theta_qt
!
  open(iunit,form='unformatted',                                 &
       file=trim(fluxdatapath)//'/theta_of_theta_qt_flabel.dat')
  read (iunit) nsqpsi,nlabel,ntheta,sqpsimin,sqpsimax,flabel_min,flabel_max &
              ,raxis,zaxis,psiaxis,sigma_qt
  allocate(theta_of_theta_qt(nlabel,0:ntheta),flabel(0:nsqpsi))
  read (iunit) theta_of_theta_qt
  read (iunit) flabel
  close(iunit)
!
  nsqp=nsqpsi
  nlab=nlabel
  nthe=ntheta+1
!
  hsqpsi=(sqpsimax-sqpsimin)/(nsqp-1)
  hlabel=(flabel_max-flabel_min)/(nlab-1)
  htheqt=2.d0*pi/ntheta
!
  allocate(sqpsi(nsqp),flab(nlab),theqt(nthe))
!
  do i=1,nsqp
    sqpsi(i)=sqpsimin+hsqpsi*(i-1)
  enddo
!
  do i=1,nlab
    flab(i)=flabel_min+hlabel*(i-1)
  enddo
!
  do i=1,nthe
    theqt(i)=htheqt*(i-1)
  enddo
!
  icp_pt=nthe*nlab
  allocate( splthet(6,6,icp_pt), ipoint_pt(nlab,nthe) )
!
  call s2dring(nlab,nthe,hlabel,htheqt,theta_of_theta_qt(:,0:ntheta), &
               icp_pt,splthet,ipoint_pt)
!
  allocate(spllabel(6,nsqpsi))
  spllabel(1,:)=flabel(1:nsqpsi)
  call spl_five_reg(nsqpsi,hsqpsi,spllabel(1,:),spllabel(2,:),spllabel(3,:) &
                   ,              spllabel(4,:),spllabel(5,:),spllabel(6,:))
!
  return
end subroutine load_theta
!
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine psithet_rz(rrr,zzz,                                          &
                        theta,theta_r,theta_z,theta_rr,theta_rz,theta_zz, &
                        flabel,s_r,s_z,s_rr,s_rz,s_zz,s0,ds0ds,dds0ds)
!
  use theta_rz_mod
  use field_eq_mod, only : nrad,nzet,rad,zet,hrad,hzet,icp,splpsi,ipoint  &
                         , psif,dpsidr,dpsidz,d2psidr2,d2psidrdz,d2psidz2
  use extract_fluxcoord_mod, only : psif_extract,theta_extract
  use spline5_RZ_mod, only: spline
! 
  real(kind=8), parameter :: pi=3.14159265358979d0
! 
  integer :: npoint,i,j,ierr,k
  real(kind=8) :: rrr,zzz,theta,theta_r,theta_z,theta_rr,theta_rz,theta_zz
  real(kind=8) :: theta_s,theta_t,theta_ss,theta_st,theta_tt
  real(kind=8) :: sqpsi_qt,s_r,s_z,s_rr,s_rz,s_zz
  real(kind=8) :: theta_qt,t_r,t_z,t_rr,t_rz,t_zz
  real(kind=8) :: rho2,rho4,dr,dz,flabel,dflabel,ddflabel,dx,dfl_dpsi,ddfl_dpsi
  real(kind=8) :: s0,ds0ds,dds0ds
! 
  if(icall.eq.0) then
    icall=1
    call load_theta
  endif
! 
  call spline(nrad,nzet,rad,zet,hrad,hzet,icp,splpsi,ipoint,rrr,zzz, &
              psif,dpsidr,dpsidz,d2psidr2,d2psidrdz,d2psidz2,ierr)
! 
  sqpsi_qt=sqrt(abs(psif-psiaxis))
!
  k=min(int(sqpsi_qt/hsqpsi),nsqp)
  dx=sqpsi_qt-hsqpsi*k
  flabel=spllabel(1,k)+dx*(spllabel(2,k)+dx*(spllabel(3,k)           &
        +dx*(spllabel(4,k)+dx*(spllabel(5,k)+dx*spllabel(6,k)))))
  dflabel=spllabel(2,k)+dx*(2.d0*spllabel(3,k)                       &
          +dx*(3.d0*spllabel(4,k)+dx*(4.d0*spllabel(5,k)             &
          +dx*5.d0*spllabel(6,k))))
  ddflabel=2.d0*spllabel(3,k)+dx*(6.d0*spllabel(4,k)                 &
           +dx*(12.d0*spllabel(5,k)+dx*20.d0*spllabel(6,k)))
!
!  dfl_dpsi=0.5d0*dflabel/sqpsi_qt
  dfl_dpsi=sign(0.5d0,psif-psiaxis)*dflabel/sqpsi_qt
  ddfl_dpsi=0.25d0*(ddflabel-dflabel/sqpsi_qt)/abs(psif-psiaxis)
  s_r=dpsidr*dfl_dpsi
  s_z=dpsidz*dfl_dpsi
  s_rr=d2psidr2*dfl_dpsi+dpsidr**2*ddfl_dpsi
  s_rz=d2psidrdz*dfl_dpsi+dpsidr*dpsidz*ddfl_dpsi
  s_zz=d2psidz2*dfl_dpsi+dpsidz**2*ddfl_dpsi
!
  s0=sqpsi_qt
  ds0ds=1.d0/dflabel
  dds0ds=-ds0ds**3*ddflabel
!
  dr=rrr-raxis
  dz=zzz-zaxis
  rho2=dr**2+dz**2
  rho4=rho2**2
  theta_qt=mod(sigma_qt*atan2(dz,dr)+2*pi,2*pi)
  t_r=-sigma_qt*dz/rho2
  t_z=sigma_qt*dr/rho2
  t_rr=2.d0*sigma_qt*dr*dz/rho4
  t_zz=-t_rr
  t_rz=sigma_qt*(dz**2-dr**2)/rho4
! 
  call spline(nlab,nthe,flab,theqt,hlabel,htheqt,icp_pt,splthet,ipoint_pt, &
              flabel,theta_qt,                                             &
              theta,theta_s,theta_t,theta_ss,theta_st,theta_tt,ierr)
!
  theta=theta+theta_qt
  theta_r=theta_s*s_r+(theta_t+1.d0)*t_r
  theta_z=theta_s*s_z+(theta_t+1.d0)*t_z
  theta_rr=theta_ss*s_r**2+2.d0*theta_st*s_r*t_r+theta_tt*t_r**2 &
          +theta_s*s_rr+(theta_t+1.d0)*t_rr
  theta_rz=theta_ss*s_r*s_z+theta_st*(s_r*t_z+s_z*t_r)+theta_tt*t_r*t_z &
          +theta_s*s_rz+(theta_t+1.d0)*t_rz
  theta_zz=theta_ss*s_z**2+2.d0*theta_st*s_z*t_z+theta_tt*t_z**2 &
          +theta_s*s_zz+(theta_t+1.d0)*t_zz
!
  psif_extract=psif
  theta_extract=theta
!
  return
end subroutine psithet_rz
!
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine cspl_five_reg(n,h,a,b,c,d,e,f)

  implicit none
!
  integer :: n,i,ip1,ip2
  double precision :: h,rhop,rhom,fac,fpl31,fpl40,fmn31,fmn40          ,x
  double precision :: a11,a12,a13,a21,a22,a23,a31,a32,a33,det
  complex(kind=8) :: abeg,bbeg,cbeg,dbeg,ebeg,fbeg
  complex(kind=8) :: aend,bend,cend,dend,eend,fend
  complex(kind=8) :: b1,b2,b3
  complex(kind=8), dimension(n) :: a,b,c,d,e,f
  complex(kind=8), dimension(:), allocatable :: alp,bet,gam
!
  rhop=13.d0+sqrt(105.d0)
  rhom=13.d0-sqrt(105.d0)
!
  a11=1.d0
  a12=1.d0/4.d0
  a13=1.d0/16.d0
  a21=3.d0
  a22=27.d0/4.d0
  a23=9.d0*27.d0/16.d0
  a31=5.d0
  a32=125.d0/4.d0
  a33=5.d0**5/16.d0
  det=a11*a22*a33+a12*a23*a31+a13*a21*a32-a12*a21*a33-a13*a22*a31-a11*a23*a32
  b1=a(4)-a(3)
  b2=a(5)-a(2)
  b3=a(6)-a(1)
  bbeg=b1*a22*a33+a12*a23*b3+a13*b2*a32-a12*b2*a33-a13*a22*b3-b1*a23*a32
  bbeg=bbeg/det
  dbeg=a11*b2*a33+b1*a23*a31+a13*a21*b3-b1*a21*a33-a13*b2*a31-a11*a23*b3
  dbeg=dbeg/det
  fbeg=a11*a22*b3+a12*b2*a31+b1*a21*a32-a12*a21*b3-b1*a22*a31-a11*b2*a32
  fbeg=fbeg/det
  b1=a(n-2)-a(n-3)
  b2=a(n-1)-a(n-4)
  b3=a(n)-a(n-5)
  bend=b1*a22*a33+a12*a23*b3+a13*b2*a32-a12*b2*a33-a13*a22*b3-b1*a23*a32
  bend=bend/det
  dend=a11*b2*a33+b1*a23*a31+a13*a21*b3-b1*a21*a33-a13*b2*a31-a11*a23*b3
  dend=dend/det
  fend=a11*a22*b3+a12*b2*a31+b1*a21*a32-a12*a21*b3-b1*a22*a31-a11*b2*a32
  fend=fend/det
  a11=2.d0
  a12=1.d0/2.d0
  a13=1.d0/8.d0
  a21=2.d0
  a22=9.d0/2.d0
  a23=81.d0/8.d0
  a31=2.d0
  a32=25.d0/2.d0
  a33=625.d0/8.d0
  det=a11*a22*a33+a12*a23*a31+a13*a21*a32-a12*a21*a33-a13*a22*a31-a11*a23*a32
  b1=a(4)+a(3)
  b2=a(5)+a(2)
  b3=a(6)+a(1)
  abeg=b1*a22*a33+a12*a23*b3+a13*b2*a32-a12*b2*a33-a13*a22*b3-b1*a23*a32
  abeg=abeg/det
  cbeg=a11*b2*a33+b1*a23*a31+a13*a21*b3-b1*a21*a33-a13*b2*a31-a11*a23*b3
  cbeg=cbeg/det
  ebeg=a11*a22*b3+a12*b2*a31+b1*a21*a32-a12*a21*b3-b1*a22*a31-a11*b2*a32
  ebeg=ebeg/det
  b1=a(n-2)+a(n-3)
  b2=a(n-1)+a(n-4)
  b3=a(n)+a(n-5)
  aend=b1*a22*a33+a12*a23*b3+a13*b2*a32-a12*b2*a33-a13*a22*b3-b1*a23*a32
  aend=aend/det
  cend=a11*b2*a33+b1*a23*a31+a13*a21*b3-b1*a21*a33-a13*b2*a31-a11*a23*b3
  cend=cend/det
  eend=a11*a22*b3+a12*b2*a31+b1*a21*a32-a12*a21*b3-b1*a22*a31-a11*b2*a32
  eend=eend/det
!
  allocate(alp(n),bet(n),gam(n))
!
  alp(1)=0.0d0
  bet(1)=ebeg*(2.d0+rhom)-5.d0*fbeg*(3.d0+1.5d0*rhom) !gamma1
!
  do i=1,n-4
    ip1=i+1
    alp(ip1)=-1.d0/(rhop+alp(i))
    bet(ip1)=alp(ip1)*(bet(i)- &
             5.d0*(a(i+4)-4.d0*a(i+3)+6.d0*a(i+2)-4.d0*a(ip1)+a(i)))
  enddo
!
  gam(n-2)=eend*(2.d0+rhom)+5.d0*fend*(3.d0+1.5d0*rhom) !gamma
  do i=n-3,1,-1
    gam(i)=gam(i+1)*alp(i)+bet(i)
  enddo
!
  alp(1)=0.0d0
  bet(1)=ebeg-2.5d0*5.d0*fbeg !e1
!
  do i=1,n-2
    ip1=i+1
    alp(ip1)=-1.d0/(rhom+alp(i))
    bet(ip1)=alp(ip1)*(bet(i)-gam(i))
  enddo
!
  e(n)=eend+2.5d0*5.d0*fend
  e(n-1)=e(n)*alp(n-1)+bet(n-1)
  f(n-1)=(e(n)-e(n-1))/5.d0
  e(n-2)=e(n-1)*alp(n-2)+bet(n-2)
  f(n-2)=(e(n-1)-e(n-2))/5.d0
  d(n-2)=dend+1.5d0*4.d0*eend+1.5d0**2*10.d0*fend
!
  do i=n-3,1,-1
    e(i)=e(i+1)*alp(i)+bet(i)
    f(i)=(e(i+1)-e(i))/5.d0
    d(i)=(a(i+3)-3.d0*a(i+2)+3.d0*a(i+1)-a(i))/6.d0 &
        -(e(i+3)+27.d0*e(i+2)+93.d0*e(i+1)+59.d0*e(i))/30.d0
    c(i)=0.5d0*(a(i+2)+a(i))-a(i+1)-0.5d0*d(i+1)-2.5d0*d(i) &
        -0.1d0*(e(i+2)+18.d0*e(i+1)+31.d0*e(i))
    b(i)=a(i+1)-a(i)-c(i)-d(i)-0.2d0*(4.d0*e(i)+e(i+1))
  enddo
!
  do i=n-3,n
    b(i)=b(i-1)+2.d0*c(i-1)+3.d0*d(i-1)+4.d0*e(i-1)+5.d0*f(i-1)
    c(i)=c(i-1)+3.d0*d(i-1)+6.d0*e(i-1)+10.d0*f(i-1)
    d(i)=d(i-1)+4.d0*e(i-1)+10.d0*f(i-1)
    if(i.ne.n) f(i)= a(i+1)-a(i)-b(i)-c(i)-d(i)-e(i)
  enddo
  f(n)=f(n-1)
!
  fac=1.d0/h
  b=b*fac
  fac=fac/h
  c=c*fac
  fac=fac/h
  d=d*fac
  fac=fac/h
  e=e*fac
  fac=fac/h
  f=f*fac
!
  deallocate(alp,bet,gam)
!
  return
end subroutine cspl_five_reg
!
!
subroutine localizer(x1,x2,x,weight,dweight,ddweight)
!
  implicit none
!
  double precision, parameter :: c1=-6.283185307179586d0,c2=-1.414213562373095d0
!
  double precision :: x1,x2,x,t,weight,dweight,ddweight,exin
!
  t=(x-x1)/(x2-x1)
!
  if(t.le.0.d0) then
    weight=1.d0
    dweight=0.d0
    ddweight=0.d0
  elseif(t.ge.1.d0) then
    weight=0.d0
    dweight=0.d0
    ddweight=0.d0
  else
    exin=exp(c2/t)
    weight=exp(c1/(1.d0-t)*exin)
    dweight=weight*c1*(1.d0/(1.d0-t)-c2/t**2)*exin/(1.d0-t)
    ddweight=dweight*c1*(1.d0/(1.d0-t)-c2/t**2)*exin/(1.d0-t) &
            +weight*c1*(1.d0/(1.d0-t)**2+2.d0*c2/t**3)*exin/(1.d0-t) &
            +weight*c1*(1.d0/(1.d0-t)-c2/t**2)**2*exin/(1.d0-t)
  endif
!
  dweight=dweight/(x2-x1)
  ddweight=ddweight/(x2-x1)**2
!
  return
end subroutine localizer
!
subroutine read_eqfile_west(nrad, nzet, psib, btf, rtf, rad, zet, psi)
  !
    use input_files, only : iunit,gfile
  !
    implicit none
  !
    integer :: nrad,nzet,ir
    real(kind=8) :: psib, btf, rtf
    real(kind=8) :: rad(nrad), zet(nzet)
    real(kind=8) :: psi(nrad,nzet)
  !
    psib=0.d0
  !
    open(unit=iunit,file=trim(gfile),status='old',action='read')
    read(iunit,*) nrad,nzet
    read(iunit,*) btf
    read(iunit,*) rad
    read(iunit,*) zet
    do ir=1,nrad
      read(iunit,*) psi(ir,:)
    enddo
    close(iunit)
    rtf=0.5d0*(rad(1)+rad(nrad))
    btf=btf/rtf
  !
end subroutine read_eqfile_west
!
subroutine read_eqfile2(nwEQD,nhEQD,psiAxis,psiSep,bt0,rzero,fpol,rad,zet,psiRZ)
  use input_files
  implicit none
  integer :: nwEQD, nhEQD
  integer :: gunit, idum
  character*10 case(6)
  integer :: i,j
  real (kind=8) :: xdim,zdim,r1,zmid,rmaxis,zmaxis,xdum
  real (kind=8) :: bt0, rzero, plas_cur, psiAxis, psiSep
  real (kind=8), dimension(nwEQD) :: fpol,pres,ffprim,pprime,qpsi
  real (kind=8), dimension(nwEQD,nhEQD) :: psiRZ
  real (kind=8) :: rad(nwEQD), zet(nhEQD)

  integer :: n_bndyxy,nlimEQD
  real (kind=8), dimension(:), allocatable :: LCFS, limEQD

      gunit=iunit

      open(unit=gunit,file=trim(gfile),status='old',action='read')

! Equilibrium Parameters
      read(gunit,2000)(case(i),i=1,6),idum,nwEQD,nhEQD
      write(*,*) 'READ_EQFILE1: ',trim(gfile),nwEQD,nhEQD
      read(gunit,2010,end=55,err=250)xdim,zdim,rzero,r1,zmid
      write(*,*) xdim, zdim, rzero, r1, zmid
      read(gunit,2010,end=55,err=250)rmaxis,zmaxis,psiAxis,psiSep,bt0
      write(*,*) rmaxis,zmaxis,psiAxis,psiSep,bt0
      read(gunit,2010,end=55,err=250)plas_cur,psiAxis,xdum,rmaxis,xdum
      write(*,*) plas_cur,psiAxis,xdum,rmaxis,xdum
      read(gunit,2010,end=55,err=250)zmaxis,xdum,psiSep,xdum,xdum
      write(*,*) zmaxis,xdum,psiSep,xdum,xdum
      read(gunit,2010,end=55,err=250)(fpol(i),i=1,nwEQD)
      read(gunit,2010,end=55,err=250)(pres(i),i=1,nwEQD)
      read(gunit,2010,end=55,err=250)(ffprim(i),i=1,nwEQD)
      read(gunit,2010,end=55,err=250)(pprime(i),i=1,nwEQD)
      read(gunit,2010,end=55,err=250)((psiRZ(i,j),i=1,nwEQD),j=1,nhEQD)
      read(gunit,2010,end=55,err=250)(qpsi(i),i=1,nwEQD)
      print *, 'Equilibrium Done.', trim(gfile)
! Boundary Data
      read(gunit,*,end=55,err=250)n_bndyxy,nlimEQD
      allocate(LCFS(2*n_bndyxy))
      allocate(limEQD(2*nlimEQD))
      read(gunit,2010,end=55,err=250)(LCFS(i),i=1,2*n_bndyxy)
      read(gunit,2010,end=55,err=250)(limEQD(i),i=1,2*nlimEQD)
!      print *, 'Boundary Done.'
      close(gunit)

      call set_eqcoords(nwEQD,nhEQD,xdim,zdim,r1,zmid,rad,zet)
  return

2000  format(6a8,3i4)
2010  format(5(e16.9))
55    print *, 'READ_EQFILE1: Early EOF in',trim(gfile); STOP
250   print *, 'READ_EQFILE1: Error reading ',trim(gfile); STOP

end subroutine read_eqfile2
!
subroutine spline_fpol
!
  use field_eq_mod, only : nrad,hfpol,splfpol
  use spline5_RZ_mod, only: spl_five_reg
!
  implicit none
!
  double precision, dimension(:), allocatable :: b,c,d,e,f
!
  allocate(b(nrad),c(nrad),d(nrad),e(nrad),f(nrad))
!
  hfpol=1.d0/dble(nrad-1)
!
  call spl_five_reg(nrad,hfpol,splfpol(0,:),b,c,d,e,f)
!
  splfpol(1,:)=b
  splfpol(2,:)=c
  splfpol(3,:)=d
  splfpol(4,:)=e
  splfpol(5,:)=f
!
  deallocate(b,c,d,e,f)
!
end subroutine spline_fpol
!
subroutine splint_fpol(x,f,fp)
  !
    use field_eq_mod, only : nrad,hfpol,splfpol
  !
    implicit none
  !
    integer :: j,k
    double precision :: x,f,fp,dx
  !
    k=max(0,int(x/hfpol))
    dx=x-k*hfpol
    k=k+1
  !
    f=splfpol(5,k)
    fp=0.d0
    do j=4,0,-1
      f=f*dx+splfpol(j,k)
      fp=fp*dx+splfpol(j+1,k)*dble(j+1)
    enddo
  !
end subroutine splint_fpol
!
subroutine read_dimeq1(nwEQD,nhEQD)
  use input_files
  implicit none
  integer :: nwEQD, nhEQD,i
  integer :: idum
  character*10 case(6)
!
     open(unit=iunit,file=trim(gfile),status='old',action='read')
     read(iunit,2000)(case(i),i=1,6),idum,nwEQD,nhEQD
     close(iunit)
  return

2000  format(6a8,3i4)
55    print *, 'READ_EQDIM1: Early EOF in',trim(gfile); STOP
250   print *, 'READ_EQDIM1: Error reading ',trim(gfile); STOP
end subroutine read_dimeq1
!
subroutine window_filter(n,nw,arr_in,arr_out)
!
  implicit none
!
  integer :: n,nw,nwa,i
  double precision, dimension(n) :: arr_in,arr_out
!
  do i=1,n
    nwa=min(nw,i-1,n-i)
    arr_out(i)=sum(arr_in(i-nwa:i+nwa))/(2*nwa+1)
  enddo
!
  return
end subroutine window_filter
!
subroutine read_eqfile1(nwEQD,nhEQD,psiSep, bt0, rzero, rad, zet, psiRZ)
  use input_files
  implicit none
  integer :: nwEQD, nhEQD
  integer :: gunit, idum
  character*10 case(6)
  integer :: i,j
  real (kind=8) :: xdim,zdim,r1,zmid,rmaxis,zmaxis,xdum
  real (kind=8) :: bt0, rzero, plas_cur, psiAxis, psiSep
  real (kind=8), dimension(nwEQD) :: fpol,pres,ffprim,pprime,qpsi
  real (kind=8), dimension(nwEQD,nhEQD) :: psiRZ
  real (kind=8) :: rad(nwEQD), zet(nhEQD)

  integer :: n_bndyxy,nlimEQD
  real (kind=8), dimension(:), allocatable :: LCFS, limEQD

      gunit=iunit

      open(unit=gunit,file=trim(gfile),status='old',action='read')

! Equilibrium Parameters
      read(gunit,2000)(case(i),i=1,6),idum,nwEQD,nhEQD
      write(*,*) 'READ_EQFILE1: ',trim(gfile),nwEQD,nhEQD
      read(gunit,2010,end=55,err=250)xdim,zdim,rzero,r1,zmid
      write(*,*) xdim, zdim, rzero, r1, zmid
      read(gunit,2010,end=55,err=250)rmaxis,zmaxis,psiAxis,psiSep,bt0
      write(*,*) rmaxis,zmaxis,psiAxis,psiSep,bt0
      read(gunit,2010,end=55,err=250)plas_cur,psiAxis,xdum,rmaxis,xdum
      write(*,*) plas_cur,psiAxis,xdum,rmaxis,xdum
      read(gunit,2010,end=55,err=250)zmaxis,xdum,psiSep,xdum,xdum
      write(*,*) zmaxis,xdum,psiSep,xdum,xdum
      read(gunit,2010,end=55,err=250)(fpol(i),i=1,nwEQD)
      read(gunit,2010,end=55,err=250)(pres(i),i=1,nwEQD)
      read(gunit,2010,end=55,err=250)(ffprim(i),i=1,nwEQD)
      read(gunit,2010,end=55,err=250)(pprime(i),i=1,nwEQD)
      read(gunit,2010,end=55,err=250)((psiRZ(i,j),i=1,nwEQD),j=1,nhEQD)
      read(gunit,2010,end=55,err=250)(qpsi(i),i=1,nwEQD)
      print *, 'Equilibrium Done.', trim(gfile)
! Boundary Data
      read(gunit,*,end=55,err=250)n_bndyxy,nlimEQD

      if (n_bndyxy > 0) then
        allocate(LCFS(2*n_bndyxy))
        read(gunit,2010,end=55,err=250)(LCFS(i),i=1,2*n_bndyxy)
      end if

      if (nlimEQD > 0) then
        allocate(limEQD(2*nlimEQD))
        read(gunit,2010,end=55,err=250)(limEQD(i),i=1,2*nlimEQD)
      end if
!      print *, 'Boundary Done.'
      close(gunit)

      call set_eqcoords(nwEQD,nhEQD,xdim,zdim,r1,zmid,rad,zet)
  return

2000  format(6a8,3i4)
2010  format(5(e16.9))
55    print *, 'READ_EQFILE1: Early EOF in',trim(gfile); STOP
250   print *, 'READ_EQFILE1: Error reading ',trim(gfile); STOP

end subroutine read_eqfile1
!
subroutine set_eqcoords(nwEQD,nhEQD,xdim,zdim,r1,zmid,rad,zet)
  implicit none
  integer :: j,k,nwEQD,nhEQD
  real (kind=8) :: xdim,zdim,r1,zmid,z1
  real (kind=8) :: rad(nwEQD), zet(nhEQD)

  do j=1,nwEQD
    rad(j) = r1 + (j-1)*(xdim/(nwEQD-1))
  end do

  z1 = zmid - zdim/2.0		! check this definition wrt zmid
  do k=1,nhEQD			! runov chooses lower, probe chooses upper
    zet(k) = z1 + (k-1)*(zdim/(nhEQD-1))
  end do

!      print *, 'set_coords done.'
!      print *, 'rad'
!      print 2010, (rad(j),j=1,nwEQD)
!      print *, ''
!      print *, 'zet'
!      print 2010, (zet(k),k=1,nhEQD)

  return
2010  format(5(e16.9))
end subroutine set_eqcoords
! Input of axisymmetric equilibrium for WEST tokamak:
!
subroutine read_dimeq_west(nrad,nzet)
!
  use input_files, only : iunit,gfile
!
  implicit none
!
  integer :: nrad,nzet
!
  open(unit=iunit,file=trim(gfile),status='old',action='read')
  read(iunit,*) nrad,nzet
  close(iunit)
!
  end subroutine read_dimeq_west


!All the remaining subroutines in this module are currently used nowhere:

!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine invert_mono_per(nx,arry_in,xmin,xmax,ny,arrx,ymin,ymax)
  !
  ! Inverts the monotonous function y(x) given on the equidistant grid
  ! of x values on the interval [xmin,xmax] by the array y_i=arry(i).
  ! Special case: y(x) is a sum of linear and periodic functions.
  ! The result, function x(y), is given on the equdistant grid of y values
  ! at the interval [ymin,ymax] by the array x_i=arrx(i).
  !
    implicit none
  !
    integer :: ny,nx,iy,ix,ixfix,ix1,ix2,ix3,ix4
  !
    double precision :: xmin,xmax,ymin,ymax,hy,y,hx,x1,x2,x3,x4,y1,y2,y3,y4
    double precision, dimension(0:nx) :: arry_in
    double precision, dimension(0:ny) :: arrx
    double precision, dimension(:), allocatable :: arry
  !
    allocate(arry(-1:nx+1))
    arry(0:nx)=arry_in
  !
    ymin=arry(0)
    ymax=arry(nx)
    arry(-1)=arry(nx-1)-ymax+ymin
    arry(nx+1)=arry(1)+ymax-ymin
  !
    hy=(ymax-ymin)/ny
    hx=(xmax-xmin)/nx
  !
    arrx(0)=xmin
    arrx(ny)=xmax
  !
    do iy=1,ny-1
      y=ymin+iy*hy
      do ix=0,nx
        if(arry(ix).gt.y) then
          ixfix=ix-3
          exit
        endif
      enddo
  !    ixfix=max(ixfix,-1)
  !    ixfix=min(ixfix,nx-4)
      ix1=ixfix+1
      ix2=ixfix+2
      ix3=ixfix+3
      ix4=ixfix+4
      x1=xmin+ix1*hx
      x2=xmin+ix2*hx
      x3=xmin+ix3*hx
      x4=xmin+ix4*hx
      y1=arry(ix1)
      y2=arry(ix2)
      y3=arry(ix3)
      y4=arry(ix4)
      arrx(iy) = x1*(y-y2)/(y1-y2)*(y-y3)/(y1-y3)*(y-y4)/(y1-y4)    &
                + x2*(y-y3)/(y2-y3)*(y-y4)/(y2-y4)*(y-y1)/(y2-y1)    &
                + x3*(y-y4)/(y3-y4)*(y-y1)/(y3-y1)*(y-y2)/(y3-y2)    &
                + x4*(y-y1)/(y4-y1)*(y-y2)/(y4-y2)*(y-y3)/(y4-y3)
    enddo
  !
    return
  end subroutine invert_mono_per
  
  !
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  subroutine extract_fluxcoord(phinorm,theta)
  !
    use extract_fluxcoord_mod
    use input_files, only : iunit,fluxdatapath
  !
    implicit none
  !
    integer :: k
    double precision :: phinorm,theta,xpsif
  !
    if(load_extract_fluxcoord.eq.1) then
      load_extract_fluxcoord=0
      open(iunit,file=trim(fluxdatapath)//'/phinorm_arr.dat')
      read (iunit,*) nphinorm,psifmin,hpsif
      allocate(phinorm_arr(nphinorm))
      do k=1,nphinorm
        read (iunit,*) phinorm_arr(k)
      enddo
      close(iunit)
    endif
  !
    xpsif=(psif_extract-psifmin)/hpsif
    k=min(nphinorm-2,max(0,int(xpsif)))
    phinorm=phinorm_arr(k+1)*(k+1-xpsif)+phinorm_arr(k+2)*(xpsif-k)
  !
    theta=theta_extract
  !
  end subroutine extract_fluxcoord
  
  !---------------------------------------------------------------------
  !
  subroutine invert_mono_reg(nx,arry,xmin,xmax,ny,arrx,ymin,ymax)
  !
  ! Inverts the monotonous function y(x) given on the equidistant grid 
  ! of x values on the interval [xmin,xmax] by the array y_i=arry(i). 
  ! The result, function x(y), is given on the equidistant grid of y values 
  ! at the interval [ymin,ymax] by the array x_i=arrx(i).
  !
    implicit none
  !
    integer :: ny,nx,iy,ix,ixfix,ix1,ix2,ix3,ix4
  !
    double precision :: xmin,xmax,ymin,ymax,hy,y,hx,x1,x2,x3,x4,y1,y2,y3,y4
    double precision, dimension(0:nx) :: arry
    double precision, dimension(0:ny) :: arrx
  !
    ymin=arry(0)
    ymax=arry(nx)
  !
    hy=(ymax-ymin)/ny
    hx=(xmax-xmin)/nx
  !
    arrx(0)=xmin
    arrx(ny)=xmax
  !
    do iy=1,ny-1
      y=ymin+iy*hy
      do ix=0,nx
        if(arry(ix).gt.y) then
          ixfix=ix-3
          exit
        endif
      enddo
      ixfix=max(ixfix,-1)
      ixfix=min(ixfix,nx-4)
      ix1=ixfix+1
      ix2=ixfix+2
      ix3=ixfix+3
      ix4=ixfix+4
      x1=xmin+ix1*hx
      x2=xmin+ix2*hx
      x3=xmin+ix3*hx
      x4=xmin+ix4*hx
      y1=arry(ix1)
      y2=arry(ix2)
      y3=arry(ix3)
      y4=arry(ix4)
      arrx(iy) = x1*(y-y2)/(y1-y2)*(y-y3)/(y1-y3)*(y-y4)/(y1-y4)    & 
                + x2*(y-y3)/(y2-y3)*(y-y4)/(y2-y4)*(y-y1)/(y2-y1)    &
                + x3*(y-y4)/(y3-y4)*(y-y1)/(y3-y1)*(y-y2)/(y3-y2)    &
                + x4*(y-y1)/(y4-y1)*(y-y2)/(y4-y2)*(y-y3)/(y4-y3)
    enddo
  !
    return
  end subroutine invert_mono_reg
  
  !---------------------------------------------------------------------
subroutine indsmp_bdf(index,nup,indu)

  ! defines interval for 1D interpolation on uniform mesh
  ! by known index.
  ! Normally looks for the central interval of stencil, but
  ! stops moving of stencil at the boundary (works for mp=4 only!)
  ! Input:
  !    index - number of a cell on the mesh
  !    nup - total number of mesh points
  ! Output:
  !    indu(mp) - relative index of stencil points
  
  ! the power 3 of polinomial is fixed strictly:

  integer, parameter :: mp = 4
  integer, dimension(mp) :: indu
  integer :: index, nup, i
                          
  indu(1) = index - 1
  if( indu(1) .le. 0 ) indu(1) = 1
  indu(mp) = indu(1) + mp - 1
  if( indu(mp) .gt. nup ) then
      indu(mp) = nup
      indu(1) = indu(mp) - mp + 1
  endif
  do i=2,mp-1
      indu(i) = indu(i-1) + 1
  enddo

  return 
end subroutine indsmp_bdf

end module utils_bdivfree_mod