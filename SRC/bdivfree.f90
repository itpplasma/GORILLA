!
  module bdivfree_mod
    integer :: nr,nz,ntor,icp
    integer, dimension(:,:), allocatable :: ipoint
    double precision :: rmin,zmin,hr,hz,pmin,pfac
    double precision, dimension(:),       allocatable :: rpoi,zpoi
    double precision, dimension(:,:,:),   allocatable :: apav,rbpav_coef
    double precision, dimension(:,:,:,:), allocatable :: aznre,aznim,arnre,arnim
  end module bdivfree_mod
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  module amn_mod
! Fourier ampitudes of the original field:
    integer :: ntor_amn=1,mpol,ntor_ff,mpol_ff,nsqpsi,icall=0
    double precision :: sqpsimin,sqpsimax,hsqpsi
    complex(kind=8), dimension(:,:,:,:), allocatable :: splapsi,splatet
    complex(kind=8), dimension(:,:), allocatable :: amnpsi,   amntet,     &
                                                   amnpsi_s, amntet_s,   &
                                                   amnpsi_ss,amntet_ss
    complex(kind=8), dimension(:),   allocatable :: expthe,expphi
! Formfactors:
    integer :: nsqpsi_ff,nmodes_ff
    double precision :: sqpsimin_ff,sqpsimax_ff,hsqpsi_ff
    integer,        dimension(:,:), allocatable :: ipoi_ff
    complex(kind=8), dimension(:,:,:), allocatable :: splffp,splfft
    complex(kind=8), dimension(:),   allocatable :: fmnpsi,   fmntet,     &
                                                   fmnpsi_s, fmntet_s,   &
                                                   fmnpsi_ss,fmntet_ss
  end module amn_mod
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  module extract_fluxcoord_mod
    integer :: load_extract_fluxcoord=1
    integer :: nphinorm
    double precision :: psif_extract,theta_extract,psifmin,hpsif
    double precision :: psifmax,phifmax,sigcos
    double precision, dimension(:), allocatable :: phinorm_arr
  end module extract_fluxcoord_mod
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  module getout_vector_potentials_mod
    double precision :: ar,az
  end module getout_vector_potentials_mod
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine vector_potentials(nr_in,np_in,nz_in,ntor_in,      &
             rmin_in,rmax_in,pmin_in,pmax_in,zmin_in,zmax_in,  &
             br,bp,bz)
!
  use bdivfree_mod
!
  implicit none
!
  double precision, parameter :: pi=3.14159265358979d0
!
  integer :: nr_in,np_in,nz_in,ntor_in,ip,np,n,ir,iz  
  integer, dimension(:), allocatable :: imi,ima,jmi,jma
!
  integer :: nashli_rukami
  integer :: irmin, irmax, i,j
  double precision, dimension(4), parameter :: weight=(/-1., 13., 13., -1./)/24.
!
  double precision :: rmin_in,rmax_in,pmin_in,pmax_in,zmin_in,zmax_in
  double precision :: hp,r,rm,zm,sumbz,hrm1,hzm1
  double precision, dimension(nr_in,np_in,nz_in)  :: br,bp,bz
  double precision, dimension(:),     allocatable :: dummy
  double precision, dimension(:,:),   allocatable :: a_re, a_im, rbpav_dummy
  double precision, dimension(:,:),   allocatable :: brm,bpm,bzm
!
  complex(kind=8) :: four_ampl
  complex(kind=8), dimension(:,:), allocatable :: expon
!
  integer, parameter :: mp=4 ! power of Lagrange's polynomial =3
  integer,          dimension(mp)    :: indx,indy
  double precision, dimension(mp)    :: xp,yp
  double precision, dimension(mp,mp) :: fp
!
  nr=nr_in
  nz=nz_in
  np=np_in-1
  ntor=ntor_in
  nashli_rukami=(nr_in+1)/2
!
  rmin=rmin_in
  zmin=zmin_in
  hr=(rmax_in-rmin_in)/(nr-1)
  hz=(zmax_in-zmin_in)/(nz-1)
  hp=2.d0*pi/np
  pmin=pmin_in
  pfac = dble(nint(2.d0*pi/(pmax_in-pmin_in)))
!
  allocate(expon(np,ntor),a_re(nr,nz),a_im(nr,nz),rbpav_dummy(nr,nz))
  allocate(imi(nz),ima(nz),jmi(nr),jma(nr), dummy(nr))
  allocate(rpoi(nr),zpoi(nz))
  allocate(brm(nr,nz),bpm(nr,nz),bzm(nr,nz))
!
  imi=1
  ima=nr
  jmi=1
  jma=nz
  do ir=1,nr
    rpoi(ir)=rmin+hr*(ir-1)
  enddo
  do iz=1,nz
    zpoi(iz)=zmin+hz*(iz-1)
  enddo
!
! Truncation of data outside the limiting convex:
!
  hrm1=1.d0/hr
  hzm1=1.d0/hz
  do ip=1,np
    do ir=1,nr
      do iz=1,nz
        call stretch_coords(rpoi(ir),zpoi(iz),rm,zm)
        call indef_bdf(rm,rmin,hrm1,nr,indx)
        call indef_bdf(zm,zmin,hzm1,nz,indy)
!
        do i=1,mp
          xp(i) = rpoi(indx(i))
          yp(i) = zpoi(indy(i))
        enddo
!
        do j=1,mp
          do i=1,mp
            fp(i,j) = Br(indx(i),ip,indy(j))
          enddo
        enddo
        call plag2d_bdf(rm,zm,fp,hrm1,hzm1,xp,yp,brm(ir,iz))
        do j=1,mp
          do i=1,mp
            fp(i,j) = Bp(indx(i),ip,indy(j))
          enddo
        enddo
        call plag2d_bdf(rm,zm,fp,hrm1,hzm1,xp,yp,bpm(ir,iz))
        bpm(ir,iz)=bpm(ir,iz)*rm/rpoi(ir)
        do j=1,mp
          do i=1,mp
            fp(i,j) = Bz(indx(i),ip,indy(j))
          enddo
        enddo
        call plag2d_bdf(rm,zm,fp,hrm1,hzm1,xp,yp,bzm(ir,iz))
!
      enddo
    enddo
    Br(:,ip,:)=brm
    Bp(:,ip,:)=bpm
    Bz(:,ip,:)=bzm
  enddo
!
! End of data truncation
!
  allocate(ipoint(nr,nz))
  icp=nr*nz
  allocate(aznre(6,6,icp,ntor),aznim(6,6,icp,ntor))
  allocate(arnre(6,6,icp,ntor),arnim(6,6,icp,ntor))
  allocate(apav(6,6,icp),rbpav_coef(6,6,icp))
!
  do n=1,ntor
    do ip=1,np
      expon(ip,n)=exp(dcmplx(0.d0,-n*(ip-1)*hp))/np
    enddo
  enddo
!
  do n=1,ntor
    do ir=1,nr
      r=rmin+hr*(ir-1)
      do iz=1,nz
        four_ampl=sum(br(ir,1:np,iz)*expon(:,n))*dcmplx(0.d0,-r/(n*pfac))
        a_re(ir,iz)=dble(four_ampl)
        a_im(ir,iz)=aimag(four_ampl)
      enddo
    enddo
    call s2dcut(nr,nz,hr,hz,a_re,imi,ima,jmi,jma,icp,aznre(:,:,:,n),ipoint)
    call s2dcut(nr,nz,hr,hz,a_im,imi,ima,jmi,jma,icp,aznim(:,:,:,n),ipoint)
    do ir=1,nr
      r=rmin+hr*(ir-1)
      do iz=1,nz
        four_ampl=sum(bz(ir,1:np,iz)*expon(:,n))*dcmplx(0.d0,r/(n*pfac))
        a_re(ir,iz)=dble(four_ampl)
        a_im(ir,iz)=aimag(four_ampl)
      enddo
    enddo
    call s2dcut(nr,nz,hr,hz,a_re,imi,ima,jmi,jma,icp,arnre(:,:,:,n),ipoint)
    call s2dcut(nr,nz,hr,hz,a_im,imi,ima,jmi,jma,icp,arnim(:,:,:,n),ipoint)
  enddo
!
  do iz=1,nz
     do ir=1,nr
        r=rmin+hr*(ir-1)
        dummy(ir) = sum(bz(ir,1:np,iz))*hr*r/np
     enddo
     a_re(nashli_rukami,iz) = 0.
     sumbz=0.d0
     do ir=nashli_rukami+1,nr
        irmax = min(ir+1,nr) 
        irmin = irmax - 3
        sumbz = sumbz + sum(dummy(irmin:irmax)*weight)
        a_re(ir,iz)=sumbz
     enddo
     sumbz=0.d0
     do ir=nashli_rukami-1,1,-1
        irmin = max(ir-1,1) 
        irmax = irmin + 3
        sumbz = sumbz - sum(dummy(irmin:irmax)*weight)
        a_re(ir,iz)=sumbz
     enddo
  enddo
!
  call s2dcut(nr,nz,hr,hz,a_re,imi,ima,jmi,jma,icp,apav,ipoint)
!
  do iz=1,nz
     do ir=1,nr
        r=rmin+hr*(ir-1)
        rbpav_dummy(ir,iz) = r*sum(bp(ir,1:np,iz))/np
     enddo
  enddo
!
  call s2dcut(nr,nz,hr,hz,rbpav_dummy,imi,ima,jmi,jma,icp,rbpav_coef,ipoint)
!
  deallocate(expon,a_re,a_im,rbpav_dummy,imi,ima,jmi,jma,dummy,brm,bpm,bzm)
!
102 format(1000e15.7)
!
  return
  end subroutine vector_potentials
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine field_divfree(r,phi,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ    &
                          ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
!
  use bdivfree_mod
  use inthecore_mod, only : incore                                            &
                          , vacf,dvacdr,dvacdz,d2vacdr2,d2vacdrdz,d2vacdz2
  use amn_mod, only : ntor_amn
  use getout_vector_potentials_mod, only : ar,az
!
  implicit none
!
  integer :: n,ierr
  double precision :: r,phi,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ    &
                     ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ
  double precision :: f,fr,fz,frr,frz,fzz
  double precision :: g,gr,gz,grr,grz,gzz
  double precision :: delbr,delbz,delbp
  double precision :: deldBrdR,deldBrdp,deldBrdZ
  double precision :: deldBpdR,deldBpdp,deldBpdZ
  double precision :: deldBzdR,deldBzdp,deldBzdZ
!  double precision :: ar,az,dar_dr,dar_dz,dar_dp,daz_dr,daz_dz,daz_dp
  double precision :: dar_dr,dar_dz,dar_dp,daz_dr,daz_dz,daz_dp
  complex(kind=8) :: expon,anr,anz,anr_r,anr_z,anz_r,anz_z
  complex(kind=8) :: anr_rr,anr_rz,anr_zz,anz_rr,anz_rz,anz_zz
!
!
  call spline(nr,nz,rpoi,zpoi,hr,hz,icp,rbpav_coef,ipoint,r,z,         &
              f,fr,fz,frr,frz,fzz,ierr)
  Bp = f/r
  dBpdR = fr/r - Bp/r
  dBpdZ = fz/r
  dBpdp = 0.d0
!
  call spline(nr,nz,rpoi,zpoi,hr,hz,icp,apav,ipoint,r,z,         &
              f,fr,fz,frr,frz,fzz,ierr)
!
  Br=-fz/r
  Bz=fr/r
  dBrdR=fz/r**2-frz/r
  dBrdZ=-fzz/r
  dBzdR=-fr/r**2+frr/r
  dBzdZ=frz/r
  dBrdp=0.d0
  dBzdp=0.d0
!
  ar=0.d0
  az=0.d0
!
  do n=1,ntor
    call spline(nr,nz,rpoi,zpoi,hr,hz,icp,arnre(:,:,:,n),ipoint,r,z,   &
                f,fr,fz,frr,frz,fzz,ierr)
    call spline(nr,nz,rpoi,zpoi,hr,hz,icp,arnim(:,:,:,n),ipoint,r,z,   &
                g,gr,gz,grr,grz,gzz,ierr)
    anr=dcmplx(f,g)
    anr_r=dcmplx(fr,gr)
    anr_z=dcmplx(fz,gz)
    anr_rr=dcmplx(frr,grr)
    anr_rz=dcmplx(frz,grz)
    anr_zz=dcmplx(fzz,gzz)
    call spline(nr,nz,rpoi,zpoi,hr,hz,icp,aznre(:,:,:,n),ipoint,r,z,   &
                f,fr,fz,frr,frz,fzz,ierr)
    call spline(nr,nz,rpoi,zpoi,hr,hz,icp,aznim(:,:,:,n),ipoint,r,z,   &
                g,gr,gz,grr,grz,gzz,ierr)
    anz=dcmplx(f,g)
    anz_r=dcmplx(fr,gr)
    anz_z=dcmplx(fz,gz)
    anz_rr=dcmplx(frr,grr)
    anz_rz=dcmplx(frz,grz)
    anz_zz=dcmplx(fzz,gzz)
!
    expon=exp(dcmplx(0.d0,n*pfac*phi))
    delbr=2.d0*dble(dcmplx(0.d0,pfac)*n*anz*expon/r)
    delbz=-2.d0*dble(dcmplx(0.d0,pfac)*n*anr*expon/r)
    delbp=2.d0*dble((anr_z-anz_r)*expon)
    deldBrdR=-delbr/r+2.d0*dble(dcmplx(0.d0,pfac)*n*anz_r*expon/r)
    deldBrdZ=2.d0*dble(dcmplx(0.d0,pfac)*n*anz_z*expon/r)
    deldBrdp=-2.d0*dble((pfac*n)**2*anz*expon/r)
    deldBzdR=-delbz/r-2.d0*dble(dcmplx(0.d0,pfac)*n*anr_r*expon/r)
    deldBzdZ=-2.d0*dble(dcmplx(0.d0,pfac)*n*anr_z*expon/r)
    deldBzdp=2.d0*dble((pfac*n)**2*anr*expon/r)
    deldBpdR=2.d0*dble((anr_rz-anz_rr)*expon)
    deldBpdZ=2.d0*dble((anr_zz-anz_rz)*expon)
    deldBpdp=2.d0*dble(dcmplx(0.d0,pfac)*n*(anr_z-anz_r)*expon)
!
!    if(incore.eq.-1.or.n.gt.ntor_amn) then
      br=br+delbr
      bz=bz+delbz
      bp=bp+delbp
      dBrdR=dBrdR+deldBrdR
      dBrdZ=dBrdZ+deldBrdZ
      dBrdp=dBrdp+deldBrdp
      dBzdR=dBzdR+deldBzdR
      dBzdZ=dBzdZ+deldBzdZ
      dBzdp=dBzdp+deldBzdp
      dBpdR=dBpdR+deldBpdR
      dBpdZ=dBpdZ+deldBpdZ
      dBpdp=dBpdp+deldBpdp
!    elseif(incore.eq.0) then
!!
      ar=ar+2.d0*dble(anr*expon)
      az=az+2.d0*dble(anz*expon)
!      dar_dr=2.d0*dble(anr_r*expon)
!      dar_dz=2.d0*dble(anr_z*expon)
!      dar_dp=-delbz*r
!      daz_dr=2.d0*dble(anz_r*expon)
!      daz_dz=2.d0*dble(anz_z*expon)
!      daz_dp=delbr*r
!!
!      br=br+delbr*vacf
!      bz=bz+delbz*vacf
!      bp=bp+delbp*vacf+ar*dvacdz-az*dvacdr
!      dBrdR=dBrdR+deldBrdR*vacf+delbr*dvacdr
!      dBrdZ=dBrdZ+deldBrdZ*vacf+delbr*dvacdz
!      dBrdp=dBrdp+deldBrdp*vacf
!      dBzdR=dBzdR+deldBzdR*vacf+delbz*dvacdr
!      dBzdZ=dBzdZ+deldBzdZ*vacf+delbz*dvacdz
!      dBzdp=dBzdp+deldBzdp*vacf
!      dBpdR=dBpdR+deldBpdR*vacf+delbp*dvacdr &
!           +dar_dr*dvacdz-daz_dr*dvacdr+ar*d2vacdrdz-az*d2vacdr2
!      dBpdZ=dBpdZ+deldBpdZ*vacf+delbp*dvacdz &
!           +dar_dz*dvacdz-daz_dz*dvacdr+ar*d2vacdz2-az*d2vacdrdz-az
!      dBpdp=dBpdp+deldBpdp*vacf              &
!           +dar_dp*dvacdz-daz_dp*dvacdr
!    endif
!
  enddo
!
  return
  end subroutine field_divfree
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
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
      implicit double precision (a-h,o-z)
!
      parameter(mp=4)
      integer indu(mp)  
                             
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
      end
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
      parameter(mp=4)
      integer indu(mp)  
                             
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
      end
!---------------------------------------------------------------------
      subroutine plag2d_bdf(x,y,fp,dxm1,dym1,xp,yp,polyl2d)
!
      implicit double precision (a-h,o-z)
!
! 2D interpolation by means of Lagrange polynomial
! the power 3 is fixed strictly:
      parameter(mp=4)
! uniform mesh (increasingly ordered) in all dimensions is implied
!
! Input parameters:
!   x,y - coordinates of the point for interpolation
!   dxm1,dym1 - inverse steps in each direction
!   xp,yp - vertices of stencil
!
! Output parameters:
! polyl2d - polynomial itself
      dimension cx(mp),cy(mp),fp(mp,mp),xp(mp),yp(mp)
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
      end
!---------------------------------------------------------------------
      subroutine coefs_bdf(u,up,dum1,cu)
!
      implicit double precision (a-h,o-z)
!
      parameter(mp=4)
      dimension up(mp),cu(mp)
      data one6/0.16666666666667d0/
      du3 = dum1**3
      cu(1) = (u - up(2)) * (u - up(3)) * (u - up(4)) * (-one6*du3)
      cu(2) = (u - up(1)) * (u - up(3)) * (u - up(4)) * (0.5d0*du3)
      cu(3) = (u - up(1)) * (u - up(2)) * (u - up(4)) * (-0.5d0*du3)
      cu(4) = (u - up(1)) * (u - up(2)) * (u - up(3)) * (one6*du3)
      return
      end
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
  end
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
  end
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine spl_five_per(n,h,a,b,c,d,e,f)
!
! Periodic spline of the 5-th order. First and last values of function must
! be the same.
!
  implicit none
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
  end
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
  implicit double precision (a-h,o-z)
! 
  dimension f(nx,ny),spl(6,6,icount),ipoint(nx,ny)
! 
  integer,          dimension(:), allocatable :: imi,ima,jmi,jma
  double precision, dimension(:), allocatable :: ai,bi,ci,di,ei,fi
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
  end
!
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine load_theta
!
  use theta_rz_mod
  use input_files, only : iunit,fluxdatapath
!
  implicit none
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
  end
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
! 
  implicit none
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
!
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
  end
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine field_fourier(r,phi,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ              &
                          ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
!
! Caution: derivatives are not computed, for derivatives call 
! a driver routine "field_fourier_derivs"
!
  use amn_mod
  use input_files,           only : iunit,fluxdatapath
  use inthecore_mod, only : incore,psi_sep                                 &
                          , plaf,dpladr,dpladz,d2pladr2,d2pladrdz,d2pladz2
  use field_eq_mod,  only : psif,dpsidr,dpsidz,d2psidr2,d2psidrdz,d2psidz2
  use theta_rz_mod,  only : psiaxis
  use bdivfree_mod,  only : pfac
!
  implicit none
!
  integer :: m,n,i,k,ierr,ntor
  double precision :: r,phi,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ                &
                     ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ
  double precision :: sqpsi,dx,g11,g12,g11_r,g11_z,g12_r,g12_z
  double precision :: theta,theta_r,theta_z,theta_rr,theta_rz,theta_zz, &
                      s_r,s_z,s_rr,s_rz,s_zz
  double precision :: fun,fr,fz,frr,frz,fzz
  double precision :: apsi,apsi_s,apsi_t,apsi_p
  double precision :: athe,athe_s,athe_t,athe_p
  double precision :: delbr,delbz,delbp,delar,delaz
  double precision :: deldBrdR,deldBrdp,deldBrdZ
  double precision :: deldBpdR,deldBpdp,deldBpdZ
  double precision :: deldBzdR,deldBzdp,deldBzdZ
  double precision :: delardR,delazdR,delardZ,delazdZ
  double precision :: fcjac,g11_t,g12_t,s0,ds0ds,dds0ds,sqpsi_sep
!
  integer, dimension(:,:), allocatable :: idummy2
!
  complex(kind=8) :: expon
  complex(kind=8), dimension(:), allocatable :: a,b,c,d,e,f
  complex(kind=8), dimension(:,:,:), allocatable :: apsimn,athetmn
!
  integer, save :: nper
!
! Initialization ------------------------------------------------------------
!
  if(icall.eq.0) then
    icall=1
!
    nper=nint(pfac)
    print *,'nper = ',nper
! Toroidally symetric part of the vacuum perturbation field - comes now
! from the cylindrical vacuum field routine
!
!
! Fourier ampitudes of the non-axisymmetric vacuum perturbation field:
!
    open(iunit,form='unformatted',file=trim(fluxdatapath)//'/amn.dat')
    read (iunit) ntor,mpol,nsqpsi,sqpsimin,sqpsimax
    allocate(apsimn(-mpol:mpol,ntor,nsqpsi))
    allocate(athetmn(-mpol:mpol,ntor,nsqpsi))
    read (iunit) apsimn,athetmn
    close(iunit)
!
    call psithet_rz(r,z,                                              &
                    theta,theta_r,theta_z,theta_rr,theta_rz,theta_zz, &
                    sqpsi,s_r,s_z,s_rr,s_rz,s_zz,s0,ds0ds,dds0ds)
!
    nsqpsi_ff=0
    open(iunit,form='unformatted',file=trim(fluxdatapath)//'/formfactors.dat')
    read (iunit,end=1) nmodes_ff,nsqpsi_ff,mpol_ff,mpol_ff,ntor_ff,ntor_ff
    read (iunit) sqpsimin_ff,sqpsimax_ff
    print *,'nsqpsi_ff = ',nsqpsi_ff
    if(ntor_ff.ne.ntor.or.mpol_ff.ne.mpol) then
      print *,'Number of harmonics in formfactors differs from original field'
      print *,'ntor = ',ntor,'ntor_ff = ',ntor_ff, &
              'mpol = ',mpol,'mpol_ff = ',mpol_ff
    endif
    allocate(ipoi_ff(-mpol_ff:mpol_ff,ntor_ff))
    allocate(splfft(nmodes_ff,6,nsqpsi_ff))
    read (iunit) ipoi_ff(-mpol_ff:mpol_ff,1:ntor_ff)
    read (iunit) splfft(1:nmodes_ff,1,1:nsqpsi_ff)
!splfft(1:nmodes_ff,1,1:nsqpsi_ff)=(0.d0,0.d0)
!ipoi_ff(1,1)=-1
!ipoi_ff(-1,1)=-1
1   close(iunit)
!
    if(nsqpsi_ff.gt.0) then
!
      call smear_formfactors(nmodes_ff,nsqpsi_ff,sqpsimin_ff,sqpsimax_ff, &
                             splfft(1:nmodes_ff,1,1:nsqpsi_ff))
!
! use those formfactors which available and necessary
      mpol_ff=min(mpol,mpol_ff)
      ntor_ff=min(ntor,ntor_ff)
      allocate(splffp(nmodes_ff,6,nsqpsi_ff))
      splffp(:,1,:)=splfft(:,1,:)
      allocate(idummy2(-mpol_ff:mpol_ff,ntor_ff))
      idummy2=ipoi_ff(-mpol_ff:mpol_ff,1:ntor_ff)
      deallocate(ipoi_ff)
      allocate(ipoi_ff(-mpol_ff:mpol_ff,ntor_ff))
      ipoi_ff=idummy2
      deallocate(idummy2)
    else
      print *,'Formfactor file formfactors.dat empty or absent,' &
             //' compute vacuum field'
      nmodes_ff=1
      mpol_ff=1
      ntor_ff=1
      nsqpsi_ff=10
      sqpsimin_ff=0.d0
      sqpsimax_ff=sqrt(abs(psi_sep-psiaxis))
      allocate(ipoi_ff(-mpol_ff:mpol_ff,ntor_ff))
      ipoi_ff=-1
      allocate(splffp(nmodes_ff,6,nsqpsi_ff))
      allocate(splfft(nmodes_ff,6,nsqpsi_ff))
      splffp(:,1,:)=1.d0
      splfft(:,1,:)=splffp(:,1,:)
    endif
!
    mpol=mpol_ff
    ntor=ntor_ff
    ntor_amn=ntor
!
    allocate(splapsi(-mpol:mpol,ntor,6,nsqpsi))
    allocate(splatet(-mpol:mpol,ntor,6,nsqpsi))
!
    allocate(amnpsi(-mpol:mpol,ntor),amntet(-mpol:mpol,ntor))
    allocate(amnpsi_s(-mpol:mpol,ntor),amntet_s(-mpol:mpol,ntor))
    allocate(amnpsi_ss(-mpol:mpol,ntor),amntet_ss(-mpol:mpol,ntor))
!
    allocate(expthe(-mpol:mpol),expphi(ntor))
!
    hsqpsi=(sqpsimax-sqpsimin)/(nsqpsi-1)
!
    allocate(a(nsqpsi),b(nsqpsi),c(nsqpsi),d(nsqpsi),e(nsqpsi),f(nsqpsi))
!
    do m=-mpol,mpol
      do n=1,ntor
        a=apsimn(m,n,:)
        call cspl_five_reg(nsqpsi,hsqpsi,a,b,c,d,e,f)
        splapsi(m,n,1,:)=a
        splapsi(m,n,2,:)=b
        splapsi(m,n,3,:)=c
        splapsi(m,n,4,:)=d
        splapsi(m,n,5,:)=e
        splapsi(m,n,6,:)=f
        a=athetmn(m,n,:)
        call cspl_five_reg(nsqpsi,hsqpsi,a,b,c,d,e,f)
        splatet(m,n,1,:)=a
        splatet(m,n,2,:)=b
        splatet(m,n,3,:)=c
        splatet(m,n,4,:)=d
        splatet(m,n,5,:)=e
        splatet(m,n,6,:)=f
      enddo
    enddo
!   
! Formfactors:
!
!
    allocate(fmnpsi(nmodes_ff))
    allocate(fmntet(nmodes_ff))
    allocate(fmnpsi_s(nmodes_ff))
    allocate(fmntet_s(nmodes_ff))
    allocate(fmnpsi_ss(nmodes_ff))
    allocate(fmntet_ss(nmodes_ff))
!
    hsqpsi_ff=(sqpsimax_ff-sqpsimin_ff)/(nsqpsi_ff-1)
!
    deallocate(a,b,c,d,e,f)
    allocate(a(nsqpsi_ff),b(nsqpsi_ff),c(nsqpsi_ff))
    allocate(d(nsqpsi_ff),e(nsqpsi_ff),f(nsqpsi_ff))
!
    do i=1,nmodes_ff
      a=splffp(i,1,:)
      call cspl_five_reg(nsqpsi_ff,hsqpsi_ff,a,b,c,d,e,f)
      splffp(i,1,:)=a
      splffp(i,2,:)=b
      splffp(i,3,:)=c
      splffp(i,4,:)=d
      splffp(i,5,:)=e
      splffp(i,6,:)=f
      a=splfft(i,1,:)
      call cspl_five_reg(nsqpsi_ff,hsqpsi_ff,a,b,c,d,e,f)
      splfft(i,1,:)=a
      splfft(i,2,:)=b
      splfft(i,3,:)=c
      splfft(i,4,:)=d
      splfft(i,5,:)=e
      splfft(i,6,:)=f
    enddo
!
! Normalize formfactors to 1 at the separatrix:
    sqpsi_sep=sqrt(abs(psi_sep-psiaxis))
    k=min(nsqpsi_ff,max(1,ceiling((sqpsi_sep-sqpsimin_ff)/hsqpsi_ff)))
    dx=sqpsi_sep-sqpsimin_ff-hsqpsi_ff*(k-1)

    fmnpsi=splffp(:,1,k)+dx*(splffp(:,2,k)+dx*(splffp(:,3,k)        &
          +dx*(splffp(:,4,k)+dx*(splffp(:,5,k)+dx*splffp(:,6,k)))))
    fmntet=splfft(:,1,k)+dx*(splfft(:,2,k)+dx*(splfft(:,3,k)        &
          +dx*(splfft(:,4,k)+dx*(splfft(:,5,k)+dx*splfft(:,6,k)))))
    do i=1,6
      do k=1,nsqpsi_ff
        splffp(:,i,k)=splffp(:,i,k)/fmnpsi
        splfft(:,i,k)=splfft(:,i,k)/fmntet
      enddo
    enddo
    splffp(:,1,:)=splffp(:,1,:)-1.d0
    splfft(:,1,:)=splfft(:,1,:)-1.d0
!
  endif
!
! End of initialization ------------------------------------------------------
!
! Toroidally symmetric part of the perturbation field - not computed, comes
! from the vacuum routine
!
  Br=0.d0
  Bp=0.d0
  Bz=0.d0
  dBrdR=0.d0
  dBrdZ=0.d0
  dBrdp=0.d0
  dBpdR=0.d0
  dBpdZ=0.d0
  dBpdp=0.d0
  dBzdR=0.d0
  dBzdZ=0.d0
  dBzdp=0.d0
!
! Asymmetric part of the perturbation field:
!
  call psithet_rz(r,z,                                              &
                  theta,theta_r,theta_z,theta_rr,theta_rz,theta_zz, &
                  sqpsi,s_r,s_z,s_rr,s_rz,s_zz,s0,ds0ds,dds0ds)
!
  g11=dpsidr**2+dpsidz**2
  g12=dpsidr*theta_r+dpsidz*theta_z
  g11_r=2.d0*dpsidr*d2psidr2+2.d0*dpsidz*d2psidrdz
  g11_z=2.d0*dpsidr*d2psidrdz+2.d0*dpsidz*d2psidz2
  g12_r=d2psidr2*theta_r+dpsidr*theta_rr+d2psidrdz*theta_z+dpsidz*theta_rz
  g12_z=d2psidrdz*theta_r+dpsidr*theta_rz+d2psidz2*theta_z+dpsidz*theta_zz
  fcjac=dpsidr*theta_z-dpsidz*theta_r
  g11_t=(dpsidr*g11_z-dpsidz*g11_r)/fcjac
  g12_t=(dpsidr*g12_z-dpsidz*g12_r)/fcjac
!
  k=min(nsqpsi,max(1,ceiling((sqpsi-sqpsimin)/hsqpsi)))
  dx=sqpsi-sqpsimin-hsqpsi*(k-1)
!
  amnpsi=splapsi(:,:,1,k)+dx*(splapsi(:,:,2,k)+dx*(splapsi(:,:,3,k)       &
        +dx*(splapsi(:,:,4,k)+dx*(splapsi(:,:,5,k)+dx*splapsi(:,:,6,k)))))
  amntet=splatet(:,:,1,k)+dx*(splatet(:,:,2,k)+dx*(splatet(:,:,3,k)       &
        +dx*(splatet(:,:,4,k)+dx*(splatet(:,:,5,k)+dx*splatet(:,:,6,k)))))
  amnpsi_s=splapsi(:,:,2,k)+dx*(2.d0*splapsi(:,:,3,k)                     &
          +dx*(3.d0*splapsi(:,:,4,k)+dx*(4.d0*splapsi(:,:,5,k)            &
          +dx*5.d0*splapsi(:,:,6,k))))
  amntet_s=splatet(:,:,2,k)+dx*(2.d0*splatet(:,:,3,k)                     &
          +dx*(3.d0*splatet(:,:,4,k)+dx*(4.d0*splatet(:,:,5,k)            &
          +dx*5.d0*splatet(:,:,6,k))))
  amnpsi_ss=2.d0*splapsi(:,:,3,k)+dx*(6.d0*splapsi(:,:,4,k)               &
           +dx*(12.d0*splapsi(:,:,5,k)+dx*20.d0*splapsi(:,:,6,k)))
  amntet_ss=2.d0*splatet(:,:,3,k)+dx*(6.d0*splatet(:,:,4,k)               &
           +dx*(12.d0*splatet(:,:,5,k)+dx*20.d0*splatet(:,:,6,k)))
!   
! Formfactors:
!
  k=min(nsqpsi_ff,max(1,ceiling((s0-sqpsimin_ff)/hsqpsi_ff)))
  dx=s0-sqpsimin_ff-hsqpsi_ff*(k-1)
!
  fmnpsi=splffp(:,1,k)+dx*(splffp(:,2,k)+dx*(splffp(:,3,k)        &
        +dx*(splffp(:,4,k)+dx*(splffp(:,5,k)+dx*splffp(:,6,k)))))
  fmntet=splfft(:,1,k)+dx*(splfft(:,2,k)+dx*(splfft(:,3,k)        &
        +dx*(splfft(:,4,k)+dx*(splfft(:,5,k)+dx*splfft(:,6,k)))))
  fmnpsi_s=splffp(:,2,k)+dx*(2.d0*splffp(:,3,k)                     &
          +dx*(3.d0*splffp(:,4,k)+dx*(4.d0*splffp(:,5,k)            &
          +dx*5.d0*splffp(:,6,k))))
  fmntet_s=splfft(:,2,k)+dx*(2.d0*splfft(:,3,k)                     &
          +dx*(3.d0*splfft(:,4,k)+dx*(4.d0*splfft(:,5,k)            &
          +dx*5.d0*splfft(:,6,k))))
  fmnpsi_ss=2.d0*splffp(:,3,k)+dx*(6.d0*splffp(:,4,k)               &
           +dx*(12.d0*splffp(:,5,k)+dx*20.d0*splffp(:,6,k)))
  fmntet_ss=2.d0*splfft(:,3,k)+dx*(6.d0*splfft(:,4,k)               &
           +dx*(12.d0*splfft(:,5,k)+dx*20.d0*splfft(:,6,k)))
!
! convert forfactor derivatives to derivatives over new label:
!
  fmnpsi_ss=fmnpsi_ss*ds0ds**2+fmnpsi_s*dds0ds
  fmnpsi_s=fmnpsi_s*ds0ds
  fmntet_ss=fmntet_ss*ds0ds**2+fmntet_s*dds0ds
  fmntet_s=fmntet_s*ds0ds
!
! Product:
!
  do m=-mpol_ff,mpol_ff
    do n=1,ntor_ff
      if(n*nper.gt.ntor_ff) cycle
      if(ipoi_ff(m,n*nper).gt.0) then
        k=ipoi_ff(m,n*nper)
        amnpsi_ss(m,n)=amnpsi_ss(m,n)*fmnpsi(k)         &
                      +2.d0*amnpsi_s(m,n)*fmnpsi_s(k)   &
                      +amnpsi(m,n)*fmnpsi_ss(k)
        amntet_ss(m,n)=amntet_ss(m,n)*fmntet(k)         &
                      +2.d0*amntet_s(m,n)*fmntet_s(k)   &
                      +amntet(m,n)*fmntet_ss(k)
        amnpsi_s(m,n) =amnpsi_s(m,n)*fmnpsi(k)          &
                      +amnpsi(m,n)*fmnpsi_s(k)
        amntet_s(m,n) =amntet_s(m,n)*fmntet(k)          &
                      +amntet(m,n)*fmntet_s(k)
        amnpsi(m,n)   =amnpsi(m,n)*fmnpsi(k)
        amntet(m,n)   =amntet(m,n)*fmntet(k)
      endif
    enddo
  enddo
!
  expthe(0)=(1.d0,0.d0)
  expthe(1)=exp(dcmplx(0.d0,theta))
  expthe(-1)=conjg(expthe(1))
  do m=2,mpol
    expthe(m)=expthe(m-1)*expthe(1)
    expthe(-m)=conjg(expthe(m))
  enddo
!  expphi(1)=exp(dcmplx(0.d0,phi))
  expphi(1)=exp(dcmplx(0.d0,pfac*phi))
  do n=2,ntor_amn
    expphi(n)=expphi(n-1)*expphi(1)
  enddo
!
  apsi=0.d0
  apsi_s=0.d0
  apsi_t=0.d0
  apsi_p=0.d0
  athe=0.d0
  athe_s=0.d0
  athe_t=0.d0
  athe_p=0.d0
  do m=-mpol,mpol
    do n=1,ntor_amn
      if(n*nper.gt.ntor_ff) cycle
      if(ipoi_ff(m,n*nper).gt.0) then
        expon=expthe(m)*expphi(n)
        apsi=apsi+2.d0*dble(expon*amnpsi(m,n))
        apsi_s=apsi_s+2.d0*dble(expon*amnpsi_s(m,n))
        apsi_t=apsi_t+2.d0*dble((0.d0,1.d0)*m*expon*amnpsi(m,n))
        apsi_p=apsi_p+2.d0*dble((0.d0,1.d0)*n*expon*amnpsi(m,n))*pfac
        athe=athe+2.d0*dble(expon*amntet(m,n))
        athe_s=athe_s+2.d0*dble(expon*amntet_s(m,n))
        athe_t=athe_t+2.d0*dble((0.d0,1.d0)*m*expon*amntet(m,n))
        athe_p=athe_p+2.d0*dble((0.d0,1.d0)*n*expon*amntet(m,n))*pfac
      endif
    enddo
  enddo
!
  delar=(apsi-g12*athe)/g11*dpsidr+athe*theta_r
  delaz=(apsi-g12*athe)/g11*dpsidz+athe*theta_z
  delbr=((apsi_p-g12*athe_p)/g11*dpsidz+athe_p*theta_z)/r
  delbz=-((apsi_p-g12*athe_p)/g11*dpsidr+athe_p*theta_r)/r
  delbp=fcjac*( (apsi_t-g12*athe_t-g12_t*athe)/g11  &
       -(apsi-g12*athe)*g11_t/g11**2 )              &
       +athe_s*(theta_r*s_z-theta_z*s_r)
!
  if(incore.eq.1) then
    Br=Br+delbr
    Bz=Bz+delbz
    Bp=Bp+delbp
  else
    Br=Br+delbr*plaf
    Bz=Bz+delbz*plaf
    Bp=Bp+delbp*plaf+delar*dpladz-delaz*dpladr
  endif
!
  end subroutine field_fourier
!
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine field_fourier_derivs(r,phi,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ    &
                                 ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
!
! Computes the field and its derivatives using central differences 
! for the field components computed by "field_fourier".
!
  implicit none
!
  double precision, parameter :: eps=1.d-7
  double precision :: rrr,ppp,zzz,r,phi,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ       &
                     ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ,del              &
                     ,rm,zm,Br0,Bp0,Bz0,dBrdR0,dBrdp0,dBrdZ0               &
                     ,dBpdR0,dBpdp0,dBpdZ0,dBzdR0,dBzdp0,dBzdZ0 
!
    del=eps*r
!
    rrr=r-del
    zzz=z
    ppp=phi
    call stretch_coords(rrr,zzz,rm,zm)
    rrr=rm
    zzz=zm
    call field_eq(rrr,ppp,zzz,Br0,Bp0,Bz0,dBrdR0,dBrdp0,dBrdZ0   &
                  ,dBpdR0,dBpdp0,dBpdZ0,dBzdR0,dBzdp0,dBzdZ0)
    call inthecore(rrr,zzz)
    call field_fourier(rrr,ppp,zzz,Br0,Bp0,Bz0,dBrdR0,dBrdp0,dBrdZ0   &
                      ,dBpdR0,dBpdp0,dBpdZ0,dBzdR0,dBzdp0,dBzdZ0)
    rrr=r+del
    zzz=z
    ppp=phi
    call stretch_coords(rrr,zzz,rm,zm)
    rrr=rm
    zzz=zm
    call field_eq(rrr,ppp,zzz,Br,Bp,Bz,dBrdR0,dBrdp0,dBrdZ0   &
                  ,dBpdR0,dBpdp0,dBpdZ0,dBzdR0,dBzdp0,dBzdZ0)
    call inthecore(rrr,zzz)
    call field_fourier(rrr,ppp,zzz,Br,Bp,Bz,dBrdR0,dBrdp0,dBrdZ0      &
                      ,dBpdR0,dBpdp0,dBpdZ0,dBzdR0,dBzdp0,dBzdZ0)
    dBrdR=0.5d0*(Br-Br0)/del
    dBpdR=0.5d0*(Bp-Bp0)/del
    dBzdR=0.5d0*(Bz-Bz0)/del
!
    rrr=r
    zzz=z-del
    ppp=phi
    call stretch_coords(rrr,zzz,rm,zm)
    rrr=rm
    zzz=zm
    call field_eq(rrr,ppp,zzz,Br0,Bp0,Bz0,dBrdR0,dBrdp0,dBrdZ0   &
                  ,dBpdR0,dBpdp0,dBpdZ0,dBzdR0,dBzdp0,dBzdZ0)
    call inthecore(rrr,zzz)
    call field_fourier(rrr,ppp,zzz,Br0,Bp0,Bz0,dBrdR0,dBrdp0,dBrdZ0   &
                      ,dBpdR0,dBpdp0,dBpdZ0,dBzdR0,dBzdp0,dBzdZ0)
    rrr=r
    zzz=z+del
    ppp=phi
    call stretch_coords(rrr,zzz,rm,zm)
    rrr=rm
    zzz=zm
    call field_eq(rrr,ppp,zzz,Br,Bp,Bz,dBrdR0,dBrdp0,dBrdZ0   &
                  ,dBpdR0,dBpdp0,dBpdZ0,dBzdR0,dBzdp0,dBzdZ0)
    call inthecore(rrr,zzz)
    call field_fourier(rrr,ppp,zzz,Br,Bp,Bz,dBrdR0,dBrdp0,dBrdZ0      &
                      ,dBpdR0,dBpdp0,dBpdZ0,dBzdR0,dBzdp0,dBzdZ0)
    dBrdZ=0.5d0*(Br-Br0)/del
    dBpdZ=0.5d0*(Bp-Bp0)/del
    dBzdZ=0.5d0*(Bz-Bz0)/del
!
    del=eps
!
    rrr=r
    zzz=z
    ppp=phi-del
    call stretch_coords(rrr,zzz,rm,zm)
    rrr=rm
    zzz=zm
    call field_eq(rrr,ppp,zzz,Br0,Bp0,Bz0,dBrdR0,dBrdp0,dBrdZ0   &
                  ,dBpdR0,dBpdp0,dBpdZ0,dBzdR0,dBzdp0,dBzdZ0)
    call inthecore(rrr,zzz)
    call field_fourier(rrr,ppp,zzz,Br0,Bp0,Bz0,dBrdR0,dBrdp0,dBrdZ0   &
                      ,dBpdR0,dBpdp0,dBpdZ0,dBzdR0,dBzdp0,dBzdZ0)
    rrr=r
    zzz=z
    ppp=phi+del
    call stretch_coords(rrr,zzz,rm,zm)
    rrr=rm
    zzz=zm
    call field_eq(rrr,ppp,zzz,Br,Bp,Bz,dBrdR0,dBrdp0,dBrdZ0   &
                  ,dBpdR0,dBpdp0,dBpdZ0,dBzdR0,dBzdp0,dBzdZ0)
    call inthecore(rrr,zzz)
    call field_fourier(rrr,ppp,zzz,Br,Bp,Bz,dBrdR0,dBrdp0,dBrdZ0      &
                      ,dBpdR0,dBpdp0,dBpdZ0,dBzdR0,dBzdp0,dBzdZ0)
    dBrdp=0.5d0*(Br-Br0)/del
    dBpdp=0.5d0*(Bp-Bp0)/del
    dBzdp=0.5d0*(Bz-Bz0)/del
!
    call field_eq(r,phi,z,Br0,Bp0,Bz0,dBrdR0,dBrdp0,dBrdZ0   &
                  ,dBpdR0,dBpdp0,dBpdZ0,dBzdR0,dBzdp0,dBzdZ0)
    call inthecore(r,z)
    call field_fourier(r,phi,z,Br,Bp,Bz,dBrdR0,dBrdp0,dBrdZ0          &
                      ,dBpdR0,dBpdp0,dBpdZ0,dBzdR0,dBzdp0,dBzdZ0)
!
  end subroutine field_fourier_derivs
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
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine smear_formfactors(nmodes_ff,nsqpsi_ff,sqpsimin_ff,sqpsimax_ff, &
                               formfactors)
!
  use inthecore_mod, only : psi_sep,psi_cut
  use theta_rz_mod,  only : psiaxis
!
  implicit none
!
  integer :: nmodes_ff,nsqpsi_ff,i
  double precision :: sqpsimin_ff,sqpsimax_ff,hsqpsi_ff,apsif
  double precision :: apsi_sep,apsi_cut,weight,dweight,ddweight
  double precision :: R=1.d0,Z=0.d0
  complex(kind=8), dimension(nmodes_ff,nsqpsi_ff) :: formfactors
!
  call inthecore(R,Z)
!
  hsqpsi_ff=(sqpsimax_ff-sqpsimin_ff)/dfloat(nsqpsi_ff-1)
  apsi_sep=abs(psi_sep-psiaxis)
  apsi_cut=abs(psi_cut-psiaxis)
!
  do i=1,nsqpsi_ff
    apsif=(sqpsimin_ff+hsqpsi_ff*dfloat(i-1))**2
    call localizer(apsi_cut,apsi_sep,apsif,weight,dweight,ddweight)
    formfactors(:,i)=weight*formfactors(:,i)+1.d0-weight
  enddo
!
  end subroutine smear_formfactors
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
