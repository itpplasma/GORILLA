!
  module exchange_get_cancoord_mod
    logical :: onlytheta
    double precision :: vartheta_c,varphi_c,sqg,aiota,Bcovar_vartheta,Bcovar_varphi,theta
  end module exchange_get_cancoord_mod
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
module get_canonical_coordinates_mod

  implicit none

  contains

  subroutine get_canonical_coordinates
!
  use odeint_mod, only: odeint_allroutines
  use canonical_coordinates_mod, only : ns_c,n_theta_c,n_phi_c,           &
                                        hs_c,h_theta_c,h_phi_c,           &
                                        ns_s_c,ns_tp_c,                   &
                                        nh_stencil,G_c,sqg_c,             &
                                        B_vartheta_c,B_varphi_c
  use vector_potential_mod, only : ns,hs
  use exchange_get_cancoord_mod, only : vartheta_c,varphi_c,sqg,aiota,Bcovar_vartheta,Bcovar_varphi, &
                                        theta,onlytheta
  use new_vmec_stuff_mod, only : n_theta,n_phi,h_theta,h_phi,ns_s,ns_tp
  use spline_vmec_data_mod, only: deallocate_vmec_spline
!
  implicit none
!
  logical :: fullset
  double precision, parameter :: relerr=1d-10
  integer :: i_theta,i_phi,i_y,i_sten,ndim,is_beg
  integer,          dimension(:),     allocatable :: ipoi_t,ipoi_p
  double precision, dimension(:),     allocatable :: y,dy
  double precision, dimension(:),     allocatable :: dstencil_theta,dstencil_phi
!
  double precision :: r,r1,r2,G_beg,dG_c_dt,dG_c_dp
  integer :: is
!
  !external rhs_cancoord
!
!
  ns_c=ns
  n_theta_c=n_theta
  n_phi_c=n_phi
  h_theta_c=h_theta
  h_phi_c=h_phi
  hs_c=hs
!
  nh_stencil=3
!
  allocate(dstencil_theta(-nh_stencil:nh_stencil),dstencil_phi(-nh_stencil:nh_stencil))
!
  if(nh_stencil.eq.1) then
    dstencil_theta(-1)=-0.5d0
    dstencil_theta(0)=0.0d0
    dstencil_theta(1)=0.5d0
  elseif(nh_stencil.eq.2) then
    dstencil_theta(-2)=1.d0/12.d0
    dstencil_theta(-1)=-2.d0/3.d0
    dstencil_theta(0)=0.0d0
    dstencil_theta(1)=2.d0/3.d0
    dstencil_theta(2)=-1.d0/12.d0
  elseif(nh_stencil.eq.3) then
    dstencil_theta(-3)=-1.d0/60.d0
    dstencil_theta(-2)=0.15d0
    dstencil_theta(-1)=-0.75d0
    dstencil_theta(0)=0.0d0
    dstencil_theta(1)=0.75d0
    dstencil_theta(2)=-0.15d0
    dstencil_theta(3)=1.d0/60.d0
  endif
!
  dstencil_phi=dstencil_theta
  dstencil_theta=dstencil_theta/h_theta_c
  dstencil_phi=dstencil_phi/h_phi_c
!
  allocate(ipoi_t(1-nh_stencil:n_theta_c+nh_stencil))
  allocate(ipoi_p(1-nh_stencil:n_phi_c+nh_stencil))
!
  do i_theta=1,n_theta_c
    ipoi_t(i_theta)=i_theta
  enddo
!
  do i_phi=1,n_phi_c
    ipoi_p(i_phi)=i_phi
  enddo
!
  do i_sten=1,nh_stencil
    ipoi_t(1-i_sten)=ipoi_t(n_theta-i_sten)
    ipoi_t(n_theta_c+i_sten)=ipoi_t(1+i_sten)
    ipoi_p(1-i_sten)=ipoi_p(n_phi_c-i_sten)
    ipoi_p(n_phi_c+i_sten)=ipoi_p(1+i_sten)
  enddo
!
  allocate(G_c(ns_c,n_theta_c,n_phi_c))
  allocate(sqg_c(ns_c,n_theta_c,n_phi_c))
  allocate(B_vartheta_c(ns_c,n_theta_c,n_phi_c))
  allocate(B_varphi_c(ns_c,n_theta_c,n_phi_c))
!
  onlytheta=.false.
  ndim=1
  allocate(y(ndim),dy(ndim))
  is_beg=ns_c/2
  G_beg=1.d-5
!
  do i_theta=1,n_theta_c
    print *,'integrate ODE: ',i_theta,' of ',n_theta_c
    vartheta_c=h_theta_c*dfloat(i_theta-1)
    do i_phi=1,n_phi_c
      varphi_c=h_phi_c*dfloat(i_phi-1)
!
      G_c(is_beg,i_theta,i_phi)=G_beg
      y(1)=G_beg
!
!      do is=is_beg-1,1,-1
      do is=is_beg-1,2,-1
        r1=hs_c*dfloat(is)
        r2=hs_c*dfloat(is-1)
!        if(is.eq.1) r2=hs_c*1d-5
!
        call odeint_allroutines(y,ndim,r1,r2,relerr,rhs_cancoord)
!
        G_c(is,i_theta,i_phi)=y(1)
      enddo
!
      y(1)=G_beg
!
      do is=is_beg+1,ns_c
        r1=hs_c*dfloat(is-2)
        r2=hs_c*dfloat(is-1)
!
        call odeint_allroutines(y,ndim,r1,r2,relerr,rhs_cancoord)
!
        G_c(is,i_theta,i_phi)=y(1)
      enddo
    enddo
  enddo
!
  do i_theta=1,n_theta_c
    print *,'compute components: ',i_theta,' of ',n_theta_c
    vartheta_c=h_theta_c*dfloat(i_theta-1)
    do i_phi=1,n_phi_c
      varphi_c=h_phi_c*dfloat(i_phi-1)
!      do is=1,ns_c
      do is=2,ns_c
        r=hs_c*dfloat(is-1)
        y(1)=G_c(is,i_theta,i_phi)
!
        call rhs_cancoord(r,y,dy)
!
        dG_c_dt=sum(dstencil_theta*G_c(is,ipoi_t(i_theta-nh_stencil:i_theta+nh_stencil),i_phi))
        dG_c_dp=sum(dstencil_phi*G_c(is,i_theta,ipoi_p(i_phi-nh_stencil:i_phi+nh_stencil)))
        sqg_c(is,i_theta,i_phi)=sqg*(1.d0+aiota*dG_c_dt+dG_c_dp)
        B_vartheta_c(is,i_theta,i_phi)=Bcovar_vartheta+(aiota*Bcovar_vartheta+Bcovar_varphi)*dG_c_dt
        B_varphi_c(is,i_theta,i_phi)=Bcovar_varphi+(aiota*Bcovar_vartheta+Bcovar_varphi)*dG_c_dp
      enddo
!First point is=1 (on axis) is bad, extrapolate with parabola:
      sqg_c(1,i_theta,i_phi) = 3.d0*(sqg_c(2,i_theta,i_phi)-sqg_c(3,i_theta,i_phi))                      &
                             + sqg_c(4,i_theta,i_phi)
      B_vartheta_c(1,i_theta,i_phi) = 3.d0*(B_vartheta_c(2,i_theta,i_phi)-B_vartheta_c(3,i_theta,i_phi)) &
                                    + B_vartheta_c(4,i_theta,i_phi)
      B_varphi_c(1,i_theta,i_phi) = 3.d0*(B_varphi_c(2,i_theta,i_phi)-B_varphi_c(3,i_theta,i_phi))       &
                                  + B_varphi_c(4,i_theta,i_phi)
    enddo
  enddo
!
  ns_s_c=ns_s
  ns_tp_c=ns_tp
  fullset=.true.
!
  call deallocate_vmec_spline(1)
!
  onlytheta=.true.
!
  call spline_can_coord(fullset)
!
  deallocate(dstencil_theta,dstencil_phi,ipoi_t,ipoi_p,y,dy,sqg_c,B_vartheta_c,B_varphi_c,G_c)
!
  end subroutine get_canonical_coordinates
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine rhs_cancoord(r,y,dy)
!
  use exchange_get_cancoord_mod, only : vartheta_c,varphi_c,sqg,aiota,Bcovar_vartheta,Bcovar_varphi, &
                                        theta,onlytheta
  use spline_vmec_data_mod, only: splint_iota, splint_lambda
  use splint_vmec_data_mod, only: vmec_field
!
  implicit none
!
  double precision, parameter :: epserr=1.d-14
  integer :: iter
  double precision :: s,varphi,A_theta,A_phi,dA_theta_ds,dA_phi_ds,                 &
                      alam,dl_ds,dl_dt,dl_dp,Bctrvr_vartheta,Bctrvr_varphi,Bcovar_r
!
  double precision :: r,vartheta,daiota_ds,deltheta
  double precision, dimension(1) :: y,dy
!
  s=r
!
  call splint_iota(s,aiota,daiota_ds)
!
  vartheta=vartheta_c+aiota*y(1)
  varphi=varphi_c+y(1)
!
! Begin Newton iteration to find VMEC theta
!
  theta = vartheta
!
  do iter=1,100
!
    call splint_lambda(s,theta,varphi,alam,dl_dt)
!
    deltheta = (vartheta-theta-alam)/(1.d0+dl_dt)
    theta = theta + deltheta
    if(abs(deltheta).lt.epserr) exit
  enddo
!
! End Newton iteration to find VMEC theta
!
  if(onlytheta) return
!
  call vmec_field(s,theta,varphi,A_theta,A_phi,dA_theta_ds,dA_phi_ds,aiota,     &
                  sqg,alam,dl_ds,dl_dt,dl_dp,Bctrvr_vartheta,Bctrvr_varphi,     &
                  Bcovar_r,Bcovar_vartheta,Bcovar_varphi)
!
  dy(1)=-(Bcovar_r+daiota_ds*Bcovar_vartheta*y(1))/(aiota*Bcovar_vartheta+Bcovar_varphi)
!
  end subroutine rhs_cancoord
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine spline_can_coord(fullset)
!
  use spl_three_to_five_mod
  use canonical_coordinates_mod, only : ns_c,n_theta_c,n_phi_c,hs_c,h_theta_c,h_phi_c,    &
                                        ns_s_c,ns_tp_c,G_c,sqg_c,B_vartheta_c,B_varphi_c, &
                                        s_sqg_c,s_B_vartheta_c,s_B_varphi_c,s_G_c
!
  implicit none
!
  logical :: fullset
  integer :: k,is,i_theta,i_phi
  integer :: iss,ist,isp
  double precision, dimension(:,:), allocatable :: splcoe
!
  allocate(s_sqg_c(ns_s_c+1,ns_tp_c+1,ns_tp_c+1,ns_c,n_theta_c,n_phi_c))
  allocate(s_B_vartheta_c(ns_s_c+1,ns_tp_c+1,ns_tp_c+1,ns_c,n_theta_c,n_phi_c))
  allocate(s_B_varphi_c(ns_s_c+1,ns_tp_c+1,ns_tp_c+1,ns_c,n_theta_c,n_phi_c))
  if(fullset) allocate(s_G_c(ns_s_c+1,ns_tp_c+1,ns_tp_c+1,ns_c,n_theta_c,n_phi_c))
!
  s_sqg_c(1,1,1,:,:,:)=sqg_c
  s_B_vartheta_c(1,1,1,:,:,:)=B_vartheta_c
  s_B_varphi_c(1,1,1,:,:,:)=B_varphi_c
  if(fullset) s_G_c(1,1,1,:,:,:)=G_c
!
! splining over $\varphi$:
!
  allocate(splcoe(0:ns_tp_c,n_phi_c))
!
  do is=1,ns_c
    do i_theta=1,n_theta_c
!
      splcoe(0,:)=s_sqg_c(1,1,1,is,i_theta,:)
!
      call spl_per(ns_tp_c,n_phi_c,h_phi_c,splcoe)
!
      do k=1,ns_tp_c
        s_sqg_c(1,1,k+1,is,i_theta,:)=splcoe(k,:)
      enddo
!
      splcoe(0,:)=s_B_vartheta_c(1,1,1,is,i_theta,:)
!
      call spl_per(ns_tp_c,n_phi_c,h_phi_c,splcoe)
!
      do k=1,ns_tp_c
        s_B_vartheta_c(1,1,k+1,is,i_theta,:)=splcoe(k,:)
      enddo
!
      splcoe(0,:)=s_B_varphi_c(1,1,1,is,i_theta,:)
!
      call spl_per(ns_tp_c,n_phi_c,h_phi_c,splcoe)
!
      do k=1,ns_tp_c
        s_B_varphi_c(1,1,k+1,is,i_theta,:)=splcoe(k,:)
      enddo
!
      if(fullset) then
!
        splcoe(0,:)=s_G_c(1,1,1,is,i_theta,:)
!
        call spl_per(ns_tp_c,n_phi_c,h_phi_c,splcoe)
!
        do k=1,ns_tp_c
          s_G_c(1,1,k+1,is,i_theta,:)=splcoe(k,:)
        enddo
!
      endif
    enddo
  enddo
!
  deallocate(splcoe)
!
! splining over $\vartheta$:
!
  allocate(splcoe(0:ns_tp_c,n_theta_c))
!
  do is=1,ns_c
    do i_phi=1,n_phi_c
      do isp=1,ns_tp_c+1
!
        splcoe(0,:)=s_sqg_c(1,1,isp,is,:,i_phi)
!
        call spl_per(ns_tp_c,n_theta_c,h_theta_c,splcoe)
!
        do k=1,ns_tp_c
          s_sqg_c(1,k+1,isp,is,:,i_phi)=splcoe(k,:)
        enddo
!
        splcoe(0,:)=s_B_vartheta_c(1,1,isp,is,:,i_phi)
!
        call spl_per(ns_tp_c,n_theta_c,h_theta_c,splcoe)
!
        do k=1,ns_tp_c
          s_B_vartheta_c(1,k+1,isp,is,:,i_phi)=splcoe(k,:)
        enddo
!
        splcoe(0,:)=s_B_varphi_c(1,1,isp,is,:,i_phi)
!
        call spl_per(ns_tp_c,n_theta_c,h_theta_c,splcoe)
!
        do k=1,ns_tp_c
          s_B_varphi_c(1,k+1,isp,is,:,i_phi)=splcoe(k,:)
        enddo
!
        if(fullset) then
!
          splcoe(0,:)=s_G_c(1,1,isp,is,:,i_phi)
!
          call spl_per(ns_tp_c,n_theta_c,h_theta_c,splcoe)
!
          do k=1,ns_tp_c
            s_G_c(1,k+1,isp,is,:,i_phi)=splcoe(k,:)
          enddo
!
        endif
!
      enddo
    enddo
  enddo
!
  deallocate(splcoe)
!
! splining over $s$:
!
  allocate(splcoe(0:ns_s_c,ns_c))
!
  do i_theta=1,n_theta_c
    do i_phi=1,n_phi_c
      do ist=1,ns_tp_c+1
        do isp=1,ns_tp_c+1
!
          splcoe(0,:)=s_sqg_c(1,ist,isp,:,i_theta,i_phi)
!
          call spl_reg(ns_s_c,ns_c,hs_c,splcoe)
!
          do k=1,ns_s_c
            s_sqg_c(k+1,ist,isp,:,i_theta,i_phi)=splcoe(k,:)
          enddo
!
          splcoe(0,:)=s_B_vartheta_c(1,ist,isp,:,i_theta,i_phi)
!
          call spl_reg(ns_s_c,ns_c,hs_c,splcoe)
!
          do k=1,ns_s_c
            s_B_vartheta_c(k+1,ist,isp,:,i_theta,i_phi)=splcoe(k,:)
          enddo
!
          splcoe(0,:)=s_B_varphi_c(1,ist,isp,:,i_theta,i_phi)
!
          call spl_reg(ns_s_c,ns_c,hs_c,splcoe)
!
          do k=1,ns_s_c
            s_B_varphi_c(k+1,ist,isp,:,i_theta,i_phi)=splcoe(k,:)
          enddo
!
          if(fullset) then
!
            splcoe(0,:)=s_G_c(1,ist,isp,:,i_theta,i_phi)
!
            call spl_reg(ns_s_c,ns_c,hs_c,splcoe)
!
            do k=1,ns_s_c
              s_G_c(k+1,ist,isp,:,i_theta,i_phi)=splcoe(k,:)
            enddo
!
          endif
!
        enddo
      enddo
    enddo
  enddo
!
  deallocate(splcoe)
!
  end subroutine spline_can_coord
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine can_to_vmec(r,vartheta_c_in,varphi_c_in,theta_vmec,varphi_vmec)
!
  use exchange_get_cancoord_mod, only : vartheta_c,varphi_c,theta
  use splint_can_coord_mod, only: splint_can_coord
!
  implicit none
!
  logical :: fullset
  double precision :: theta_vmec,varphi_vmec
  double precision :: r,vartheta_c_in,varphi_c_in,                                     &
                      A_phi,A_theta,dA_phi_dr,dA_theta_dr,d2A_phi_dr2,                 &
                      sqg_c,dsqg_c_dr,dsqg_c_dt,dsqg_c_dp,                             &
                      B_vartheta_c,dB_vartheta_c_dr,dB_vartheta_c_dt,dB_vartheta_c_dp, &
                      B_varphi_c,dB_varphi_c_dr,dB_varphi_c_dt,dB_varphi_c_dp,G_c
  double precision, dimension(1) :: y,dy
!
  fullset=.true.
!
  call splint_can_coord(r,vartheta_c_in,varphi_c_in,                                     &
                        A_theta,A_phi,dA_theta_dr,dA_phi_dr,d2A_phi_dr2,                 &
                        sqg_c,dsqg_c_dr,dsqg_c_dt,dsqg_c_dp,                             &
                        B_vartheta_c,dB_vartheta_c_dr,dB_vartheta_c_dt,dB_vartheta_c_dp, &
                        B_varphi_c,dB_varphi_c_dr,dB_varphi_c_dt,dB_varphi_c_dp,         &
                        fullset,G_c)
!
  vartheta_c=vartheta_c_in
  varphi_c=varphi_c_in
  y(1)=G_c
!
  call rhs_cancoord(r,y,dy)
!
  theta_vmec=theta
  varphi_vmec=varphi_c_in+G_c
!
  end subroutine can_to_vmec
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine deallocate_can_coord
!
  use canonical_coordinates_mod, only : s_sqg_c,s_B_vartheta_c,s_B_varphi_c,s_G_c
!
  implicit none
!
  deallocate(s_sqg_c,s_B_vartheta_c,s_B_varphi_c)
  if(allocated(s_G_c)) deallocate(s_G_c)
!
  end subroutine deallocate_can_coord

end module get_canonical_coordinates_mod