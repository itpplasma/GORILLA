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
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
module getout_vector_potentials_mod
  double precision :: ar,az
end module getout_vector_potentials_mod
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
module inthecore_mod
  logical :: prop=.true.
  integer :: npoi,ijumpb,ibeg,iend
  double precision, parameter :: epssep=1.d-6
  double precision :: rc,zc,twopi,sig,psi_sep,psi_cut,sigpsi,cutoff
  double precision, dimension(:), allocatable :: rho2i,theti
  integer          :: incore
  double precision :: vacf,dvacdr,dvacdz,d2vacdr2,d2vacdrdz,d2vacdz2
  double precision :: plaf,dpladr,dpladz,d2pladr2,d2pladrdz,d2pladz2
end module inthecore_mod
!

module field_mod
  integer          :: icall=0
  integer          :: ipert,iequil,iaxieq=0
  double precision :: ampl
end module field_mod
!
module utils_field_divB0_mod

  implicit none

  contains

subroutine field_divfree(r,phi,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ    &
                          ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
!
  use bdivfree_mod
  use inthecore_mod, only : incore                                            &
                          , vacf,dvacdr,dvacdz,d2vacdr2,d2vacdrdz,d2vacdz2
  use amn_mod, only : ntor_amn
  use getout_vector_potentials_mod, only : ar,az
  use spline5_RZ_mod, only: spline
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
subroutine vector_potentials(nr_in,np_in,nz_in,ntor_in,      &
  rmin_in,rmax_in,pmin_in,pmax_in,zmin_in,zmax_in,  &
  br,bp,bz)
!
  use bdivfree_mod
  use spline5_RZ_mod, only: s2dcut
  use utils_bdivfree_mod, only: indef_bdf, plag2d_bdf
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
subroutine field_fourier(r,phi,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ              &
  ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
!
! Caution: derivatives are not computed, for derivatives call 
! a driver routine "field_fourier_derivs"
!
  use amn_mod
  use input_files,   only : iunit,fluxdatapath
  use inthecore_mod, only : incore,psi_sep                                 &
    , plaf,dpladr,dpladz,d2pladr2,d2pladrdz,d2pladz2
  use field_eq_mod,  only : psif,dpsidr,dpsidz,d2psidr2,d2psidrdz,d2psidz2
  use theta_rz_mod,  only : psiaxis
  use bdivfree_mod,  only : pfac
  use utils_bdivfree_mod, only:  psithet_rz, cspl_five_reg
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
subroutine inthecore(R,Z)
!
  use inthecore_mod
  use input_files,  only : iunit,fluxdatapath
  use field_eq_mod, only : psif,dpsidr,dpsidz,d2psidr2,d2psidrdz,d2psidz2
  use utils_bdivfree_mod, only: localizer
!
  implicit none
!
  integer :: i
  double precision :: R,Z,rho2,thet,scalp,xx,yy
  double precision :: weight,dweight,ddweight
  double precision, dimension(4) :: x,y
  double precision, dimension(:), allocatable :: ri,zi
!
  if(prop) then
    prop=.false.
    open(iunit,file=trim(fluxdatapath)//'/separ.dat')
    read(iunit,*) x(1:2)
    psi_sep=x(2)+(x(1)-x(2))*(1.d0-epssep)
    psi_cut=x(2)+(x(1)-x(2))*cutoff
    sigpsi=sign(1.d0,psi_sep-psi_cut)
    npoi=0
    do while(npoi.ge.0)
      npoi=npoi+1
      read(iunit,*,end=1) rc
    enddo
  1  allocate(ri(0:npoi),zi(0:npoi),rho2i(0:npoi),theti(0:npoi))
    ri=0.d0
    zi=0.d0
    close(iunit)
    open(iunit,file=trim(fluxdatapath)//'/separ.dat')
    read(iunit,*)
    do i=1,npoi-1
      read(iunit,*) ri(i),zi(i)
    enddo
    close(iunit)
    rc=sum(ri(1:npoi-1))/(npoi-1)
    zc=sum(zi(1:npoi-1))/(npoi-1)
    rho2i=(ri-rc)**2+(zi-zc)**2
    theti=atan2(zi-zc,ri-rc)
    sig=theti(2)-theti(1)
    do i=2,npoi-2
      if((theti(i+1)-theti(i))*sig.lt.0.d0) then
        ijumpb=i
        exit
      endif
    enddo
    twopi=8.d0*atan2(1.d0,1.d0)
    ri(1:npoi-1-ijumpb)=rho2i(ijumpb+1:npoi-1)
    ri(npoi-ijumpb:npoi-1)=rho2i(1:ijumpb)
    rho2i=ri
    ri(1:npoi-1-ijumpb)=theti(ijumpb+1:npoi-1)
    ri(npoi-ijumpb:npoi-1)=theti(1:ijumpb)
    theti=ri
    deallocate(ri,zi)
    sig=theti(2)-theti(1)
    rho2i(npoi)=rho2i(1)
    theti(npoi)=theti(1)+sign(twopi,sig)
    rho2i(0)=rho2i(npoi-1)
    theti(0)=theti(npoi-1)-sign(twopi,sig)
  endif
!
  rho2=(r-rc)**2+(z-zc)**2
  thet=atan2(z-zc,r-rc)
!
  ibeg=0
  iend=npoi
  do while(ibeg+1.lt.iend)
    i=(ibeg+iend)/2
    if((thet-theti(i))*sig.gt.0.d0) then
      ibeg=i
    else
      iend=i
    endif
  enddo
  iend=min(iend,npoi-1)
  ibeg=iend-1
  x=theti(ibeg-1:iend+1)
  y=rho2i(ibeg-1:iend+1)
!
  xx=thet
  yy=y(1)*(xx-x(2))/(x(1)-x(2))*(xx-x(3))/(x(1)-x(3))*(xx-x(4))/(x(1)-x(4)) &
    +y(2)*(xx-x(3))/(x(2)-x(3))*(xx-x(4))/(x(2)-x(4))*(xx-x(1))/(x(2)-x(1)) &
    +y(3)*(xx-x(4))/(x(3)-x(4))*(xx-x(1))/(x(3)-x(1))*(xx-x(2))/(x(3)-x(2)) &
    +y(4)*(xx-x(1))/(x(4)-x(1))*(xx-x(2))/(x(4)-x(2))*(xx-x(3))/(x(4)-x(3))
!
  if(rho2.gt.yy) then
    incore=-1
    return
  elseif((psif-psi_cut)*sigpsi.lt.0.d0) then
    incore=1
    return
  endif
!
  incore=0
!
  call localizer(psi_cut,psi_sep,psif,weight,dweight,ddweight)
!
  plaf=weight
  dpladr=dweight*dpsidr
  dpladz=dweight*dpsidz
  d2pladr2=ddweight*dpsidr**2+dweight*d2psidr2
  d2pladrdz=ddweight*dpsidr*dpsidz+dweight*d2psidrdz
  d2pladz2=ddweight*dpsidz**2+dweight*d2psidz2
!
  vacf=1.d0-plaf
  dvacdr=-dpladr
  dvacdz=-dpladz
  d2vacdr2=-d2pladr2
  d2vacdrdz=-d2pladrdz
  d2vacdz2=-d2pladz2
!
  return
end subroutine inthecore

subroutine stretch_coords(r,z,rm,zm)
  use input_files, only : iunit,convexfile
  implicit none
  integer icall, i, j, nrz ! number of points "convex wall" in input file
  integer, parameter :: nrzmx=100 ! possible max. of nrz
  integer, parameter :: nrhotht=360
  real(kind=8), parameter :: pi = 3.14159265358979d0
  real(kind=8) R0,Rw, Zw, htht, Rl, Zl, a, b, r, z, rm, zm, rho, tht, rho_c, delta
  real(kind=8), dimension(100):: rad_w, zet_w ! points "convex wall"
  real(kind=8), dimension(:), allocatable :: rho_w, tht_w
  real(kind=8), dimension(nrhotht) :: rho_wall, tht_wall ! polar coords of CW
  data icall /0/, delta/1./
  save
!----------- 1st call --------------------------------------------------------
  if(icall .eq. 0) then
      icall = 1
      nrz = 0
      rad_w = 0.
      zet_w = 0.
      open(iunit,file=trim(convexfile))
      do i=1,nrzmx
        read(iunit,*,END=10)rad_w(i),zet_w(i)
        nrz = nrz + 1
      enddo
10   continue
      close(iunit)

      nrz = nrz+1
      rad_w(nrz) = rad_w(1)
      zet_w(nrz) = zet_w(1)
      allocate(rho_w(nrz), tht_w(nrz))
      R0 = (maxval(rad_w(1:nrz)) +  minval(rad_w(1:nrz)))*0.5
      do i=1,nrz
        rho_w(i) = sqrt( (rad_w(i)-R0)**2 + zet_w(i)**2 )
        tht_w(i) = atan2(zet_w(i),(rad_w(i)-R0))
        if(tht_w(i) .lt. 0.) tht_w(i) = tht_w(i) + 2.*pi
      enddo
      htht = 2.*pi/(nrhotht-1)
      do i=2,nrhotht
        tht_wall(i) = htht*(i-1)
        do j=1,nrz-1
            if(tht_wall(i).ge.tht_w(j) .and. tht_wall(i).le.tht_w(j+1)) then
              if( abs((rad_w(j+1) - rad_w(j))/rad_w(j)) .gt. 1.e-3) then
                  a = (zet_w(j+1) - zet_w(j))/(rad_w(j+1) - rad_w(j))
                  b = zet_w(j) - a*(rad_w(j) - R0)
                  Rw = b/(tan(tht_wall(i)) - a) + R0
                  Zw = a*(Rw - R0) + b
              else
                  a = (rad_w(j+1) - rad_w(j))/(zet_w(j+1) - zet_w(j))
                  b = rad_w(j)-R0 - a*zet_w(j)
                  Zw = b/(1./tan(tht_wall(i)) - a)
                  Rw = a*Zw + b + R0
              endif
            endif
        enddo
        rho_wall(i) = sqrt((Rw-R0)**2 + Zw**2)
      enddo
      tht_wall(1) = 0.
      rho_wall(1) = rho_wall(nrhotht)
!!$  do i=1,nrhotht
!!$     write(19,*)tht_wall(i), rho_wall(i)
!!$  enddo
  endif
!----------- end of the 1st call --------------------------------------------
  rm = r
  zm = z
  rho = sqrt((r-R0)**2 + z**2)
  tht = atan2(z,(r-R0))
  if(tht .lt. 0.) tht = tht + 2.*pi
  i = modulo(int(tht/htht), nrhotht-1) + 1
  rho_c = (rho_wall(i+1) - rho_wall(i))/(tht_wall(i+1) - tht_wall(i))   &
        *(tht - tht_wall(i)) + rho_wall(i)
!print *,rho,rho_c,i,tht
  if(rho .ge. rho_c) then
      rho = rho_c + delta*atan2((rho-rho_c), delta)
      rm = rho*cos(tht) + R0
      zm = rho*sin(tht)
  endif

  return
end subroutine stretch_coords
!
subroutine smear_formfactors(nmodes_ff,nsqpsi_ff,sqpsimin_ff,sqpsimax_ff, &
  formfactors)
!
  use inthecore_mod, only : psi_sep,psi_cut
  use theta_rz_mod,  only : psiaxis
  use utils_bdivfree_mod, only: localizer
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
!
subroutine field_eq(r,ppp,z,Brad,Bphi,Bzet,dBrdR,dBrdp,dBrdZ  &
  ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
!
  use input_files
  use field_eq_mod
  use field_mod, only : iaxieq
  use spline5_RZ_mod, only: spline, s2dcut
  use utils_bdivfree_mod, only: read_eqfile_west, read_eqfile2, spline_fpol, splint_fpol, read_dimeq1, &
                                window_filter, read_eqfile1, read_dimeq_west
!
  implicit none
!
  integer :: ierr,i,j
!
  double precision :: rrr,ppp,zzz,Brad,Bphi,Bzet,dBrdR,dBrdp,dBrdZ  &
      ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ,r,z
  double precision :: psihat,fpol,fpol_prime                                         !<=18.12.18
!
!-------first call: read data from disk-------------------------------
  if(icall_eq .lt. 1) then
!
  select case(iaxieq)
  case(0)
  print *,'axisymmetric equilibrium: EFIT format'
!
  call read_dimeq1(nrad,nzet)
!
  allocate(rad(nrad),zet(nzet))
  allocate(psi0(nrad,nzet),psi(nrad,nzet))
  allocate(splfpol(0:5,nrad))
!
  case(1)
  print *,'axisymmetric equilibrium: WEST format'
  use_fpol=.false.
!
  call read_dimeq_west(nrad,nzet)
!
  allocate(rad(nrad),zet(nzet))
  allocate(psi0(nrad,nzet),psi(nrad,nzet))
!
  call read_eqfile_west(nrad, nzet, psib, btf, rtf, rad, zet, psi)
!
  case default
  print *,'axisymmetric equilibrium: unknown format'
  stop
  end select
!
  if(use_fpol) then                                                                !<=18.12.18
  call read_eqfile2(nrad, nzet, psi_axis, psi_sep, btf, rtf,    &                !<=18.12.18
        splfpol(0,:), rad, zet, psi)                                 !<=18.12.18
  psib=-psi_axis                                                                 !<=18.12.18
  psi_sep=(psi_sep-psi_axis)*1.d8                                                !<=18.12.18
  splfpol(0,:)=splfpol(0,:)*1.d6                                                 !<=18.12.18
  call spline_fpol                                                               !<=18.12.18
  elseif(iaxieq.eq.0) then
  call read_eqfile1(nrad, nzet, psib, btf, rtf, rad, zet, psi)
  endif                                                                            !<=18.12.18
!
  ! Filtering:
  !open(191,file='psi_orig.dat')
  !do i=1,nrad
  !write (191,*) psi(i,:)+psib
  !enddo
  !close(191)
  !stop
  do i=1,nzet
  call window_filter(nrad,nwindow_r,psi(:,i),psi0(:,i))
  enddo
!
  do i=1,nrad
  call window_filter(nzet,nwindow_z,psi0(i,:),psi(i,:))
  enddo
  !open(191,file='psi_filt.dat')
  !do i=1,nrad
  !write (191,*) psi(i,:)
  !enddo
  !close(191)
  !stop
  ! End filtering
  !     allocate(xi(nzet),f(nzet))
  !     npoint = nzet
  !     xi = zet
  !     do i=1,nrad
  !        f = psi(i,:)
  !        call leastsq(npoint,xi,f)
  !        psi0(i,:) = f
  !     enddo
  !     deallocate(xi,f)

  !     allocate(xi(nrad),f(nrad))
  !     npoint = nrad
  !     xi = rad
  !     do i=1,nzet
  !        f = psi0(:,i)
  !        call leastsq(npoint,xi,f)
  !        psi(:,i) = f
  !     enddo
!
  rad = rad*1.d2 ! cm
  zet = zet*1.d2 ! cm
  rtf = rtf*1.d2 ! cm
  psi = psi*1.d8
  psib= psib*1.d8
  btf = btf*1.d4
!
  psi=psi+psib
!
  hrad = rad(2) - rad(1)
  hzet = zet(2) - zet(1)
!
! rectangular domain:
  allocate( imi(nzet),ima(nzet),jmi(nrad),jma(nrad) )
  imi = 1
  ima = nrad
  jmi = 1
  jma = nzet
!
!  Computation of the number of data in splpsi
  icp = 0
  do i=1,nzet
  if ( imi(i) .gt. 0 .and. ima(i) .gt. 0 ) then
  icp = icp + ima(i) - imi(i) + 1
  endif
  enddo
  write(6,*) 'number of points in the table:  ',icp
!
  allocate( splpsi(6,6,icp), ipoint(nrad,nzet) )
!
  call s2dcut(nrad,nzet,hrad,hzet,psi,imi,ima,jmi,jma,icp,splpsi,ipoint)
!
  if(icall_eq.eq.-1) then
! Quit after initialization with zero field
  Brad=0.d0
  Bphi=0.d0
  Bzet=0.d0
  dBrdR=0.d0
  dBrdp=0.d0
  dBrdZ=0.d0
  dBpdR=0.d0
  dBpdp=0.d0
  dBpdZ=0.d0
  dBzdR=0.d0
  dBzdp=0.d0
  dBzdZ=0.d0
  icall_eq = 1
  return
  endif
  icall_eq = 1
  endif
!
! ------- end first call ----------------------------------------------
  rrr=max(rad(1),min(rad(nrad),r))
  zzz=max(zet(1),min(zet(nzet),z))
!
  call spline(nrad,nzet,rad,zet,hrad,hzet,icp,splpsi,ipoint,rrr,zzz, &
  psif,dpsidr,dpsidz,d2psidr2,d2psidrdz,d2psidz2,ierr)
!
  Brad = -dpsidz/rrr
  Bzet =  dpsidr/rrr
!
! axisymmetric case:
  dBrdp = 0.
  dBpdp = 0.
  dBzdp = 0.
!
  dBrdR = -d2psidrdz/rrr+dpsidz/rrr**2
  dBzdZ =  d2psidrdz/rrr
  dBrdZ = -d2psidz2/rrr
  dBzdR =  d2psidr2/rrr-dpsidr/rrr**2
!
  if(use_fpol) then                                                                  !<=18.12.18
  psihat=psif/psi_sep                                                              !<=18.12.18
  if(psihat.gt.1.d0) then                                                          !<=18.12.18
  fpol=splfpol(0,nrad)                                                           !<=18.12.18
  fpol_prime=0.d0                                                                !<=18.12.18
  else                                                                             !<=18.12.18
  call splint_fpol(psihat,fpol,fpol_prime)                                       !<=18.12.18
  endif                                                                            !<=18.12.18
  Bphi = fpol/rrr                                                                  !<=18.12.18
  dBpdR = fpol_prime*dpsidr/(psi_sep*rrr)-fpol/rrr**2                              !<=18.12.18
  dBpdZ = fpol_prime*dpsidz/(psi_sep*rrr)                                          !<=18.12.18
  else                                                                               !<=18.12.18
  Bphi = btf*rtf/rrr
  dBpdR = -btf*rtf/rrr**2
  dBpdZ = 0.
  endif                                                                              !<=18.12.18
!
  return
end subroutine field_eq

end module utils_field_divB0_mod