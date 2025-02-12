module splint_can_coord_mod

  implicit none

  contains

  subroutine splint_can_coord(r,vartheta_c,varphi_c,                                           &
                              A_theta,A_phi,dA_theta_dr,dA_phi_dr,d2A_phi_dr2,                 &
                              sqg_c,dsqg_c_dr,dsqg_c_dt,dsqg_c_dp,                             &
                              B_vartheta_c,dB_vartheta_c_dr,dB_vartheta_c_dt,dB_vartheta_c_dp, &
                              B_varphi_c,dB_varphi_c_dr,dB_varphi_c_dt,dB_varphi_c_dp,         &
                              fullset,G_c)
!
  use canonical_coordinates_mod, only : ns_c,n_theta_c,n_phi_c,hs_c,h_theta_c,h_phi_c,    &
                                        ns_s_c,ns_tp_c,                                   &
                                        s_sqg_c,s_B_vartheta_c,s_B_varphi_c,s_G_c
  use vector_potential_mod, only : ns,hs,torflux,sA_phi
  use new_vmec_stuff_mod,   only : nper,ns_A
!
  implicit none
!
  integer, parameter :: ns_max=6
  double precision, parameter :: twopi=2.d0*3.14159265358979d0
!
  logical :: fullset
!
  integer :: nstp
  integer :: k,is,i_theta,i_phi
  integer :: iss,ist,isp
!
  double precision :: r,vartheta_c,varphi_c,                                           &
                      A_phi,A_theta,dA_phi_dr,dA_theta_dr,d2A_phi_dr2,                 &
                      sqg_c,dsqg_c_dr,dsqg_c_dt,dsqg_c_dp,                             &
                      B_vartheta_c,dB_vartheta_c_dr,dB_vartheta_c_dt,dB_vartheta_c_dp, &
                      B_varphi_c,dB_varphi_c_dr,dB_varphi_c_dt,dB_varphi_c_dp,G_c
  double precision :: s,ds,dtheta,dphi
!
  double precision, dimension(ns_max)        :: sp_sqg,sp_bt,sp_bp,sp_G
  double precision, dimension(ns_max)        :: dsp_sqg_ds,dsp_bt_ds,dsp_bp_ds
  double precision, dimension(ns_max)        :: dsp_sqg_dt,dsp_bt_dt,dsp_bp_dt
  double precision, dimension(ns_max,ns_max) :: stp_sqg,stp_bt,stp_bp,stp_G
  double precision, dimension(ns_max,ns_max) :: dstp_sqg_ds,dstp_bt_ds,dstp_bp_ds
!
  A_theta=torflux*r
  dA_theta_dr=torflux
!
  dtheta=modulo(vartheta_c,twopi)/h_theta_c
  i_theta=max(0,min(n_theta_c-1,int(dtheta)))
  dtheta=(dtheta-dfloat(i_theta))*h_theta_c
  i_theta=i_theta+1
!
  dphi=modulo(varphi_c,twopi/dfloat(nper))/h_phi_c
  i_phi=max(0,min(n_phi_c-1,int(dphi)))
  dphi=(dphi-dfloat(i_phi))*h_phi_c
  i_phi=i_phi+1
!
! Begin interpolation over $s$
!
  ds=r/hs
  is=max(0,min(ns-1,int(ds)))
  ds=(ds-dfloat(is))*hs
  is=is+1
!
  A_phi=sA_phi(ns_A+1,is)
  dA_phi_dr=0.d0
!
  do k=ns_A,1,-1
    A_phi=sA_phi(k,is)+ds*A_phi
    dA_phi_dr=sA_phi(k+1,is)*dfloat(k)+ds*dA_phi_dr
  enddo
!
  d2A_phi_dr2=0.d0
!
  do k=ns_A,2,-1
    d2A_phi_dr2=sA_phi(k+1,is)*dfloat(k)*dfloat(k-1)+ds*d2A_phi_dr2
  enddo
!
  if(ns_c.ne.ns) then
    ds=r/hs_c
    is=max(0,min(ns_c-1,int(ds)))
    ds=(ds-dfloat(is))*hs_c
    is=is+1
  endif
!
  nstp=ns_tp_c+1
!
  stp_sqg(1:nstp,1:nstp)=s_sqg_c(ns_s_c+1,:,:,is,i_theta,i_phi)
  dstp_sqg_ds(1:nstp,1:nstp)=0.d0
!
  stp_bt(1:nstp,1:nstp)=s_B_vartheta_c(ns_s_c+1,:,:,is,i_theta,i_phi)
  dstp_bt_ds(1:nstp,1:nstp)=0.d0
!
  stp_bp(1:nstp,1:nstp)=s_B_varphi_c(ns_s_c+1,:,:,is,i_theta,i_phi)
  dstp_bp_ds(1:nstp,1:nstp)=0.d0
!
  if(fullset) stp_G(1:nstp,1:nstp)=s_G_c(ns_s_c+1,:,:,is,i_theta,i_phi)
!
  do k=ns_s_c,1,-1
    stp_sqg(1:nstp,1:nstp)=s_sqg_c(k,:,:,is,i_theta,i_phi)+ds*stp_sqg(1:nstp,1:nstp)
    dstp_sqg_ds(1:nstp,1:nstp)=s_sqg_c(k+1,:,:,is,i_theta,i_phi)*dfloat(k)+ds*dstp_sqg_ds(1:nstp,1:nstp)
!
    stp_bt(1:nstp,1:nstp)=s_B_vartheta_c(k,:,:,is,i_theta,i_phi)+ds*stp_bt(1:nstp,1:nstp)
    dstp_bt_ds(1:nstp,1:nstp)=s_B_vartheta_c(k+1,:,:,is,i_theta,i_phi)*dfloat(k)+ds*dstp_bt_ds(1:nstp,1:nstp)
!
    stp_bp(1:nstp,1:nstp)=s_B_varphi_c(k,:,:,is,i_theta,i_phi)+ds*stp_bp(1:nstp,1:nstp)
    dstp_bp_ds(1:nstp,1:nstp)=s_B_varphi_c(k+1,:,:,is,i_theta,i_phi)*dfloat(k)+ds*dstp_bp_ds(1:nstp,1:nstp)
!
    if(fullset) stp_G(1:nstp,1:nstp)=s_G_c(k,:,:,is,i_theta,i_phi)+ds*stp_G(1:nstp,1:nstp)
  enddo
!
! End interpolation over $s$
!----------------------------
!
! Begin interpolation over $\theta$
!
  sp_sqg(1:nstp)=stp_sqg(nstp,1:nstp)
  dsp_sqg_ds(1:nstp)=dstp_sqg_ds(nstp,1:nstp)
  dsp_sqg_dt(1:nstp)=0.d0
  sp_bt(1:nstp)=stp_bt(nstp,1:nstp)
  dsp_bt_ds(1:nstp)=dstp_bt_ds(nstp,1:nstp)
  dsp_bt_dt(1:nstp)=0.d0
  sp_bp(1:nstp)=stp_bp(nstp,1:nstp)
  dsp_bp_ds(1:nstp)=dstp_bp_ds(nstp,1:nstp)
  dsp_bp_dt(1:nstp)=0.d0
  if(fullset) sp_G(1:nstp)=stp_G(nstp,1:nstp)
!
  do k=ns_tp_c,1,-1
    sp_sqg(1:nstp)=stp_sqg(k,1:nstp)+dtheta*sp_sqg(1:nstp)
    dsp_sqg_ds(1:nstp)=dstp_sqg_ds(k,1:nstp)+dtheta*dsp_sqg_ds(1:nstp)
    dsp_sqg_dt(1:nstp)=stp_sqg(k+1,1:nstp)*dfloat(k)+dtheta*dsp_sqg_dt(1:nstp)
!
    sp_bt(1:nstp)=stp_bt(k,1:nstp)+dtheta*sp_bt(1:nstp)
    dsp_bt_ds(1:nstp)=dstp_bt_ds(k,1:nstp)+dtheta*dsp_bt_ds(1:nstp)
    dsp_bt_dt(1:nstp)=stp_bt(k+1,1:nstp)*dfloat(k)+dtheta*dsp_bt_dt(1:nstp)
!
    sp_bp(1:nstp)=stp_bp(k,1:nstp)+dtheta*sp_bp(1:nstp)
    dsp_bp_ds(1:nstp)=dstp_bp_ds(k,1:nstp)+dtheta*dsp_bp_ds(1:nstp)
    dsp_bp_dt(1:nstp)=stp_bp(k+1,1:nstp)*dfloat(k)+dtheta*dsp_bp_dt(1:nstp)
!
    if(fullset) sp_G(1:nstp)=stp_G(k,1:nstp)+dtheta*sp_G(1:nstp)
  enddo
!
! End interpolation over $\theta$
!--------------------------------
!
! Begin interpolation over $\varphi$
!
  sqg_c=sp_sqg(nstp)
  dsqg_c_dr=dsp_sqg_ds(nstp)
  dsqg_c_dt=dsp_sqg_dt(nstp)
  dsqg_c_dp=0.d0
!
  B_vartheta_c=sp_bt(nstp)
  dB_vartheta_c_dr=dsp_bt_ds(nstp)
  dB_vartheta_c_dt=dsp_bt_dt(nstp)
  dB_vartheta_c_dp=0.d0
!
  B_varphi_c=sp_bp(nstp)
  dB_varphi_c_dr=dsp_bp_ds(nstp)
  dB_varphi_c_dt=dsp_bp_dt(nstp)
  dB_varphi_c_dp=0.d0
!
  if(fullset) G_c=sp_G(nstp)
!
  do k=ns_tp_c,1,-1
    sqg_c=sp_sqg(k)+dphi*sqg_c
    dsqg_c_dr=dsp_sqg_ds(k)+dphi*dsqg_c_dr
    dsqg_c_dt=dsp_sqg_dt(k)+dphi*dsqg_c_dt
    dsqg_c_dp=sp_sqg(k+1)*dfloat(k)+dphi*dsqg_c_dp
!
    B_vartheta_c=sp_bt(k)+dphi*B_vartheta_c
    dB_vartheta_c_dr=dsp_bt_ds(k)+dphi*dB_vartheta_c_dr
    dB_vartheta_c_dt=dsp_bt_dt(k)+dphi*dB_vartheta_c_dt
    dB_vartheta_c_dp=sp_bt(k+1)*dfloat(k)+dphi*dB_vartheta_c_dp
!
    B_varphi_c=sp_bp(k)+dphi*B_varphi_c
    dB_varphi_c_dr=dsp_bp_ds(k)+dphi*dB_varphi_c_dr
    dB_varphi_c_dt=dsp_bp_dt(k)+dphi*dB_varphi_c_dt
    dB_varphi_c_dp=sp_bp(k+1)*dfloat(k)+dphi*dB_varphi_c_dp
!
    if(fullset) G_c=sp_G(k)+dphi*G_c
  enddo
!
! End interpolation over $\varphi$
!
  end subroutine splint_can_coord

end module splint_can_coord_mod