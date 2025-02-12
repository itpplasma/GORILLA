module splint_vmec_data_mod

  implicit none

  contains

subroutine splint_vmec_data(s,theta,varphi,A_phi,A_theta,dA_phi_ds,dA_theta_ds,aiota,       &
  R,Z,alam,dR_ds,dR_dt,dR_dp,dZ_ds,dZ_dt,dZ_dp,dl_ds,dl_dt,dl_dp)
!
  use new_vmec_stuff_mod,   only : n_theta,n_phi,h_theta,h_phi,sR,sZ,slam,nper,ns_A,ns_s,ns_tp
  use vector_potential_mod, only : ns,hs,torflux,sA_phi
!
  implicit none
!
  double precision, parameter :: twopi=2.d0*3.14159265358979d0
!
  integer :: is,i_theta,i_phi,k
  double precision :: ds,dtheta,dphi,rho_tor
  double precision :: s,theta,varphi,A_phi,A_theta,dA_phi_ds,dA_theta_ds,aiota,       &
  R,Z,alam,dR_ds,dR_dt,dR_dp,dZ_ds,dZ_dt,dZ_dp,dl_ds,dl_dt,dl_dp
!
  integer, parameter :: ns_max=6
!
  integer :: nstp
!
  double precision, dimension(ns_max)        :: sp_R,sp_Z,sp_lam
  double precision, dimension(ns_max)        :: dsp_R_ds,dsp_Z_ds,dsp_lam_ds
  double precision, dimension(ns_max)        :: dsp_R_dt,dsp_Z_dt,dsp_lam_dt
  double precision, dimension(ns_max,ns_max) :: stp_R,stp_Z,stp_lam
  double precision, dimension(ns_max,ns_max) :: dstp_R_ds,dstp_Z_ds,dstp_lam_ds
!
  nstp=ns_tp+1
!
  A_theta=torflux*s
  dA_theta_ds=torflux
!
  ds=s/hs
  is=max(0,min(ns-1,int(ds)))
  ds=(ds-dble(is))*hs
  is=is+1
!
  dtheta=modulo(theta,twopi)/h_theta
  i_theta=max(0,min(n_theta-1,int(dtheta)))
  dtheta=(dtheta-dble(i_theta))*h_theta
  i_theta=i_theta+1
!
  dphi=modulo(varphi,twopi/dble(nper))/h_phi
  i_phi=max(0,min(n_phi-1,int(dphi)))
  dphi=(dphi-dble(i_phi))*h_phi
  i_phi=i_phi+1
!
! Begin interpolation over $s$
!
! Vector potential $A_\varphi$ and its derivative:
  A_phi=sA_phi(ns_A+1,is)
  dA_phi_ds=0.d0
!
  do k=ns_A,1,-1
  A_phi=sA_phi(k,is)+ds*A_phi
  dA_phi_ds=sA_phi(k+1,is)*dble(k)+ds*dA_phi_ds
  enddo
!
! R, Z and $\lambda$ and their derivatives over $s$:
!
  rho_tor=sqrt(s)
  ds=rho_tor/hs
  is=max(0,min(ns-2,int(ds)))
  ds=(ds-dble(is))*hs
  is=is+1
!
  stp_R(1:nstp,1:nstp)=sR(ns_s+1,:,:,is,i_theta,i_phi)
  dstp_R_ds(1:nstp,1:nstp)=0.d0
  stp_Z(1:nstp,1:nstp)=sZ(ns_s+1,:,:,is,i_theta,i_phi)
  dstp_Z_ds(1:nstp,1:nstp)=0.d0
  stp_lam(1:nstp,1:nstp)=slam(ns_s+1,:,:,is,i_theta,i_phi)
  dstp_lam_ds(1:nstp,1:nstp)=0.d0
!
  do k=ns_s,1,-1
  stp_R(1:nstp,1:nstp)=sR(k,:,:,is,i_theta,i_phi)+ds*stp_R(1:nstp,1:nstp)
  dstp_R_ds(1:nstp,1:nstp)=sR(k+1,:,:,is,i_theta,i_phi)*dble(k)+ds*dstp_R_ds(1:nstp,1:nstp)
  stp_Z(1:nstp,1:nstp)=sZ(k,:,:,is,i_theta,i_phi)+ds*stp_Z(1:nstp,1:nstp)
  dstp_Z_ds(1:nstp,1:nstp)=sZ(k+1,:,:,is,i_theta,i_phi)*dble(k)+ds*dstp_Z_ds(1:nstp,1:nstp)
  stp_lam(1:nstp,1:nstp)=slam(k,:,:,is,i_theta,i_phi)+ds*stp_lam(1:nstp,1:nstp)
  dstp_lam_ds(1:nstp,1:nstp)=slam(k+1,:,:,is,i_theta,i_phi)*dble(k)+ds*dstp_lam_ds(1:nstp,1:nstp)
  enddo
!
! End interpolation over $s$
!----------------------------
!
  aiota=-dA_phi_ds/dA_theta_ds
!
! Begin interpolation over $\theta$
!
  sp_R(1:nstp)=stp_R(nstp,1:nstp)
  dsp_R_ds(1:nstp)=dstp_R_ds(nstp,1:nstp)
  dsp_R_dt(1:nstp)=0.d0
  sp_Z(1:nstp)=stp_Z(nstp,1:nstp)
  dsp_Z_ds(1:nstp)=dstp_Z_ds(nstp,1:nstp)
  dsp_Z_dt(1:nstp)=0.d0
  sp_lam(1:nstp)=stp_lam(nstp,1:nstp)
  dsp_lam_ds(1:nstp)=dstp_lam_ds(nstp,1:nstp)
  dsp_lam_dt(1:nstp)=0.d0
!
  do k=ns_tp,1,-1
  sp_R(1:nstp)=stp_R(k,1:nstp)+dtheta*sp_R(1:nstp)
  dsp_R_ds(1:nstp)=dstp_R_ds(k,1:nstp)+dtheta*dsp_R_ds(1:nstp)
  dsp_R_dt(1:nstp)=stp_R(k+1,1:nstp)*dble(k)+dtheta*dsp_R_dt(1:nstp)
!
  sp_Z(1:nstp)=stp_Z(k,1:nstp)+dtheta*sp_Z(1:nstp)
  dsp_Z_ds(1:nstp)=dstp_Z_ds(k,1:nstp)+dtheta*dsp_Z_ds(1:nstp)
  dsp_Z_dt(1:nstp)=stp_Z(k+1,1:nstp)*dble(k)+dtheta*dsp_Z_dt(1:nstp)
!
  sp_lam(1:nstp)=stp_lam(k,1:nstp)+dtheta*sp_lam(1:nstp)
  dsp_lam_ds(1:nstp)=dstp_lam_ds(k,1:nstp)+dtheta*dsp_lam_ds(1:nstp)
  dsp_lam_dt(1:nstp)=stp_lam(k+1,1:nstp)*dble(k)+dtheta*dsp_lam_dt(1:nstp)
  enddo
!
! End interpolation over $\theta$
!--------------------------------
!
! Begin interpolation over $\varphi$
!
  R=sp_R(nstp)
  dR_ds=dsp_R_ds(nstp)
  dR_dt=dsp_R_dt(nstp)
  dR_dp=0.d0
  Z=sp_Z(nstp)
  dZ_ds=dsp_Z_ds(nstp)
  dZ_dt=dsp_Z_dt(nstp)
  dZ_dp=0.d0
  alam=sp_lam(nstp)
  dl_ds=dsp_lam_ds(nstp)
  dl_dt=dsp_lam_dt(nstp)
  dl_dp=0.d0
!
  do k=ns_tp,1,-1
  R=sp_R(k)+dphi*R
  dR_ds=dsp_R_ds(k)+dphi*dR_ds
  dR_dt=dsp_R_dt(k)+dphi*dR_dt
  dR_dp=sp_R(k+1)*dble(k)+dphi*dR_dp
!
  Z=sp_Z(k)+dphi*Z
  dZ_ds=dsp_Z_ds(k)+dphi*dZ_ds
  dZ_dt=dsp_Z_dt(k)+dphi*dZ_dt
  dZ_dp=sp_Z(k+1)*dble(k)+dphi*dZ_dp
!
  alam=sp_lam(k)+dphi*alam
  dl_ds=dsp_lam_ds(k)+dphi*dl_ds
  dl_dt=dsp_lam_dt(k)+dphi*dl_dt
  dl_dp=sp_lam(k+1)*dble(k)+dphi*dl_dp
  enddo
!
! End interpolation over $\varphi$
!
  dR_ds=0.5d0*dR_ds/rho_tor
  dZ_ds=0.5d0*dZ_ds/rho_tor
  dl_ds=0.5d0*dl_ds/rho_tor
!
end subroutine splint_vmec_data

subroutine vmec_field(s,theta,varphi,A_theta,A_phi,dA_theta_ds,dA_phi_ds,aiota,     &
  sqg,alam,dl_ds,dl_dt,dl_dp,Bctrvr_vartheta,Bctrvr_varphi,     &
  Bcovar_r,Bcovar_vartheta,Bcovar_varphi)    
!
  implicit none
!
  double precision :: s,theta,varphi,A_phi,A_theta,dA_phi_ds,dA_theta_ds,aiota,        &
    R,Z,alam,dR_ds,dR_dt,dR_dp,dZ_ds,dZ_dt,dZ_dp,dl_ds,dl_dt,dl_dp
  double precision :: Bctrvr_vartheta,Bctrvr_varphi,Bcovar_vartheta,Bcovar_varphi,sqg
  double precision :: cjac,sqgV,Bcovar_r
  double precision, dimension(3,3) :: cmat,gV,g
!
!
  call splint_vmec_data(s,theta,varphi,A_phi,A_theta,dA_phi_ds,dA_theta_ds,aiota,      &
      R,Z,alam,dR_ds,dR_dt,dR_dp,dZ_ds,dZ_dt,dZ_dp,dl_ds,dl_dt,dl_dp)
!
  gV(1,1)=dR_ds**2+dZ_ds**2
  gV(1,2)=dR_ds*dR_dt+dZ_ds*dZ_dt
  gV(1,3)=dR_ds*dR_dp+dZ_ds*dZ_dp
  gV(2,1)=gV(1,2)
  gV(2,2)=dR_dt**2+dZ_dt**2
  gV(2,3)=dR_dt*dR_dp+dZ_dt*dZ_dp
  gV(3,1)=gV(1,3)
  gV(3,2)=gV(2,3)
  gV(3,3)=R**2+dR_dp**2+dZ_dp**2
  sqgV=R*(dR_dt*dZ_ds-dR_ds*dZ_dt)
!
  cjac=1.d0/(1.d0+dl_dt)
  sqg=sqgV*cjac
  Bctrvr_vartheta=-dA_phi_ds/sqg
  Bctrvr_varphi=dA_theta_ds/sqg
!
  cmat(1,2:3)=0.d0
  cmat(3,1:2)=0.d0
  cmat(1,1)=1.d0
  cmat(3,3)=1.d0
  cmat(2,1)=-dl_ds*cjac
  cmat(2,2)=cjac
  cmat(2,3)=-dl_dp*cjac
!
  g=matmul(transpose(cmat),matmul(gV,cmat))
!
  Bcovar_r=g(1,2)*Bctrvr_vartheta+g(1,3)*Bctrvr_varphi
  Bcovar_vartheta=g(2,2)*Bctrvr_vartheta+g(2,3)*Bctrvr_varphi
  Bcovar_varphi=g(3,2)*Bctrvr_vartheta+g(3,3)*Bctrvr_varphi
!
end subroutine vmec_field

end module splint_vmec_data_mod