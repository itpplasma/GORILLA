!
  module magdata_in_symfluxcoor_mod
    integer, parameter :: nspl  = 3 !spline order in poloidal angle interpolation
    integer, parameter :: nplag = 4 !stencil for Largange polynomial interpolation (polynomial order + 1)
    integer, parameter :: nder  = 1 !number of derivatives from Largange polynomial interpolation
    logical :: load=.true.
    integer :: nlabel,ntheta
    double precision :: rmn,rmx,zmn,zmx,raxis,zaxis,h_theta,twopi,psipol_max,psitor_max
    double precision, dimension(nplag)        :: R_lag,Z_lag,sqrtg_lag,bmod_lag,dbmod_dt_lag, &
                                                 dR_dt_lag, dZ_dt_lag
    double precision, dimension(0:nder,nplag) :: coef
    double precision, dimension(:),     allocatable :: rbeg,rsmall,qsaf,psisurf,phitor
    double precision, dimension(:,:,:), allocatable :: R_st,Z_st,bmod_st,sqgnorm_st
  end module magdata_in_symfluxcoor_mod
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine load_magdata_in_symfluxcoord
!
! Reads and splines magnetic data over theta
!
  use magdata_in_symfluxcoor_mod
  use spl_three_to_five_mod, only : spl_per
!
  implicit none
!
  integer :: i,nthetap1
  double precision, dimension(:,:), allocatable :: splcoe
!
  twopi = atan(1.d0)*8.d0
!
!-----------------------------------------------------------------------
!
  open(1,form='formatted',file='box_size_axis.dat')
  read (1,*) rmn,rmx     !'<= rmn, rmx (cm)'
  read (1,*) zmn,zmx     !'<= zmn, zmx (cm)'
  read (1,*) raxis,zaxis !'<= raxis, zaxis (cm)'
  close(1)
!
!-----------------------------------------------------------------------
!
  open(1,form='formatted',file='twodim_functions.dat')
  read (1,*) nlabel, ntheta   !'<= nlabel, ntheta'
!
  allocate(rbeg(nlabel),rsmall(nlabel),qsaf(nlabel),psisurf(0:nlabel),phitor(0:nlabel))
  allocate(R_st(0:nspl,0:ntheta,nlabel))
  allocate(Z_st(0:nspl,0:ntheta,nlabel))
  allocate(bmod_st(0:nspl,0:ntheta,nlabel))
  allocate(sqgnorm_st(0:nspl,0:ntheta,nlabel))
!
  allocate(splcoe(0:nspl,0:ntheta))
  h_theta = twopi/dfloat(ntheta)
  nthetap1=ntheta+1
!
  read (1,*)
  do i=1,nlabel
    read (1,*) splcoe(0,1:ntheta)
    splcoe(0,0)=splcoe(0,ntheta)
!
    call spl_per(nspl,nthetap1,h_theta,splcoe)
!
    R_st(:,:,i)=splcoe
  enddo
!
  read (1,*)
  do i=1,nlabel
    read (1,*) splcoe(0,1:ntheta)
    splcoe(0,0)=splcoe(0,ntheta)
!
    call spl_per(nspl,nthetap1,h_theta,splcoe)
!
    Z_st(:,:,i)=splcoe
  enddo
!
  read (1,*)
  do i=1,nlabel
    read (1,*) splcoe(0,1:ntheta)
    splcoe(0,0)=splcoe(0,ntheta)
!
    call spl_per(nspl,nthetap1,h_theta,splcoe)
!
    bmod_st(:,:,i)=splcoe
  enddo
!
  read (1,*)
  do i=1,nlabel
    read (1,*) splcoe(0,1:ntheta)
    splcoe(0,0)=splcoe(0,ntheta)
! 
    call spl_per(nspl,nthetap1,h_theta,splcoe)
!
    sqgnorm_st(:,:,i)=splcoe
  enddo
  close(1)
!
  deallocate(splcoe)
!
!-----------------------------------------------------------------------
!
  open(1,form='formatted',file='flux_functions.dat')
  read (1,*)
  do i=1,nlabel
    read (1,*) rbeg(i),rsmall(i),qsaf(i),psisurf(i),phitor(i)
  enddo
  close(1)
!
  psisurf(0)=0.d0
  phitor(0)=0.d0
  psipol_max=psisurf(nlabel)
  psitor_max=phitor(nlabel)
  psisurf=psisurf/psipol_max
  phitor=phitor/psitor_max
!
  load=.false.
!
  end subroutine load_magdata_in_symfluxcoord
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine magdata_in_symfluxcoord_ext(inp_label,s,psi,theta,q,dq_ds, &
                                         sqrtg,bmod,dbmod_dtheta,R,dR_ds,dR_dtheta,       &
                                         Z,dZ_ds,dZ_dtheta)
  !
  ! Computes safety factor, sqrt(g) of symmetry flux coordinates, module of B, cylindrical 
  ! coordinates R and Z
  ! and poloidal angle of symmetry flux coordinates. 
  ! Computes also derivative of module-B over the poloidal angle and derivatives of R, Z 
  ! Uses periodic spline interpolation over theta and Lagrange polynomial interpolation over flux srface label.
  ! Two flux surface labels can be used as an input: normalized toroidal flux s and dimensional poloidal
  ! flux psi, depending on the input switch inp_label. If inp_label=1 s is an input variable and psi=psi(s) is computed. 
  ! If inp_label=2, psi is an input variable, and s=s(psi) is computed.
  ! 
  ! Input (inout) arguments:
  !                 inp_label - input switch: 1 for s and 2 for psi
  !                 s         - normalized toroidal flux (inout)
  !                 psi       - poloidal flux (inout)
  !                 theta     - poloidal angle of symmetry flux coordinates
  ! Output arguments:
  !                 q         - safety factor
  !                 dq_ds     - derivative of safety factor over normalized toroidal or poloidal flux
  !                 sqrtg     - normalized metric determinant sqrt(g)
  !                 bmod      - module of B
  !                 dbmod_dt  - derivative of module of B over theta
  !                 R         - radius R of cylindrical coordinates
  !                 dR_ds     - derivative of R over normalized toroidal or poloidal flux
  !                 dR_dtheta - derivative of R over polidal angle of symmetry flux coordinates theta
  !                 Z         - height Z of cylindrical coordinates
  !                 dZ_ds     - derivative of Z over normalized toroidal or poloidal flux
  !                 dZ_dtheta - derivative of Z over polidal angle of symmetry flux coordinates theta
  !
  use magdata_in_symfluxcoor_mod
  !
  implicit none
  !
  integer :: inp_label,ir,it,k,km1,ibeg,iend
  double precision :: s,psi,theta,q,dq_ds,sqrtg,bmod,dbmod_dtheta,dtheta, &
                      R,dR_ds,dR_dtheta,Z,dZ_ds,dZ_dtheta
!
  if(inp_label.eq.1) then
!
    call binsrc(phitor(0:nlabel),0,nlabel,s,ibeg)
  !
    ibeg=max(1,ibeg-nplag/2)
    iend=ibeg+nplag-1
    if(iend.gt.nlabel) then
      iend=nlabel
      ibeg=iend+1-nplag
    endif
  !
    call plag_coeff(nplag,nder,s,phitor(ibeg:iend),coef)
  !
    psi=sum(coef(0,:)*psisurf(ibeg:iend))*psipol_max
  elseif(inp_label.eq.2) then
    s=psi/psipol_max
  !
    call binsrc(psisurf(0:nlabel),0,nlabel,s,ibeg)
  !
    ibeg=max(1,ibeg-nplag/2)
    iend=ibeg+nplag-1
    if(iend.gt.nlabel) then
      iend=nlabel
      ibeg=iend+1-nplag
    endif
!
    call plag_coeff(nplag,nder,s,psisurf(ibeg:iend),coef)
  !
    s=sum(coef(0,:)*phitor(ibeg:iend))
  else
    print *,'unknown mode for inp_label =',inp_label
    return
  endif
!
  q=sum(coef(0,:)*qsaf(ibeg:iend))
  dq_ds=sum(coef(1,:)*qsaf(ibeg:iend))
!
  dtheta=modulo(theta,twopi)/h_theta
  it=max(0,min(ntheta-1,int(dtheta)))
  dtheta=(dtheta-dfloat(it))*h_theta
!
  sqrtg_lag=sqgnorm_st(nspl,it,ibeg:iend)
  bmod_lag=bmod_st(nspl,it,ibeg:iend)
  R_lag=R_st(nspl,it,ibeg:iend)
  Z_lag=Z_st(nspl,it,ibeg:iend)
!
  dbmod_dt_lag=0.d0
  dR_dt_lag=0.d0
  dZ_dt_lag=0.d0
!
  do k=nspl,1,-1
    km1=k-1
    sqrtg_lag=sqrtg_lag*dtheta+sqgnorm_st(km1,it,ibeg:iend)
    bmod_lag=bmod_lag*dtheta+bmod_st(km1,it,ibeg:iend)
    R_lag=R_lag*dtheta+R_st(km1,it,ibeg:iend)
    Z_lag=Z_lag*dtheta+Z_st(km1,it,ibeg:iend)
  !
    dbmod_dt_lag=dbmod_dt_lag*dtheta+bmod_st(k,it,ibeg:iend)*dfloat(k)
    dR_dt_lag=dR_dt_lag*dtheta+R_st(k,it,ibeg:iend)*dfloat(k)
    dZ_dt_lag=dZ_dt_lag*dtheta+Z_st(k,it,ibeg:iend)*dfloat(k)
  enddo
!
  sqrtg=sum(coef(0,:)*sqrtg_lag)
  bmod=sum(coef(0,:)*bmod_lag)
  R=sum(coef(0,:)*R_lag)
  Z=sum(coef(0,:)*Z_lag)
!
  dR_ds=sum(coef(1,:)*R_lag)
  dZ_ds=sum(coef(1,:)*Z_lag)
!
  dbmod_dtheta=sum(coef(0,:)*dbmod_dt_lag)
  dR_dtheta=sum(coef(0,:)*dR_dt_lag)
  dZ_dtheta=sum(coef(0,:)*dZ_dt_lag)
!
end subroutine magdata_in_symfluxcoord_ext
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
