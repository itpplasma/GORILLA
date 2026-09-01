module vmecin_mod

  implicit none

  contains
!
! Usage:
!
!    call vmecin(rmnc,zmns,almns,rmns,zmnc,almnc,aiota,phi,sps,axm,axn,s, &
!               nsurfm,nstrm,kpar,torflux)
!
!  where scalars are:

!  nstrm - (integer) number of harmonics
!  kpar  - (integer) number of radial points

!  vectors are:
!  aiota(0:kpar)             - (double precision) iota profile
!  sps(0:kpar) = 0:kpar      - (double precision) radial index as dble number
!  phi(0:kpar)               - (double precision) toroidal flux
!  s(0:kpar)                 - (double precision) normalized toroidal flux
!  axm(nstrm)                - (double precision) poloidal mode numbers
!  axn(nstrm)                - (double precision) toroidal mode numbers

!  matrices are:
!  rmnc/rmns(nstrm,0:kpar) - cosine/sine amplitudes of R
!  zmnc/zmns(nstrm,0:kpar) - cosine/sine amplitudes of Z
!  almnc/almns(nstrm,0:kpar) - cosine/sine amplitudes of lambda
!
!

  subroutine vmecin(rmnc,zmns,almns,rmns,zmnc,almnc,aiota,phi,sps,axm,axn,s, &
                    nsurfb,nstrb,kparb,flux)
!
  use new_vmec_stuff_mod, only : netcdffile
  use nctools_module
!
  implicit none
!
  double precision, parameter :: pi=3.14159265358979d0
  double precision, parameter :: fac_b=1d4,fac_r=1d2
!
  integer :: nsurfb,nstrb,kparb,ncid,i,lasym_int
  logical :: lasym
  double precision :: flux
  double precision, dimension(nstrb)         :: axm,axn
  double precision, dimension(0:kparb)       :: sps,aiota,phi,s
  double precision, dimension(nstrb,0:kparb) :: rmnc,zmnc,almnc,lmnc
  double precision, dimension(nstrb,0:kparb) :: rmns,zmns,almns,lmns
!
  do i=0,kparb
    sps(i)=dble(i)
  enddo
!
  call nc_open(netcdffile, ncid)
!
  lasym = .false.
  if (nc_has_var(ncid, 'lasym__logical__')) then
    call nc_get(ncid, 'lasym__logical__', lasym_int)
    lasym = lasym_int.eq.1
  endif
!
  call nc_get(ncid, 'phi', phi)
  phi = phi/(2*pi)  ! added by Christopher Albert, 2019-09-16 for correct normalization
!
  flux=phi(kparb)
  flux=flux*fac_b*fac_r**2
  phi=phi*fac_b*fac_r**2
  s=phi/flux
!
  call nc_get(ncid, 'xm', axm)
  call nc_get(ncid, 'xn', axn)
  call nc_get(ncid, 'iotaf', aiota)
  call nc_get(ncid, 'rmnc', rmnc)
  call nc_get(ncid, 'zmns', zmns)
  call nc_get(ncid, 'lmns', lmns)
!
  if (lasym) then
    call nc_get(ncid, 'rmns', rmns)
    call nc_get(ncid, 'zmnc', zmnc)
    call nc_get(ncid, 'lmnc', lmnc)
  else
    rmns=0.d0
    zmnc=0.d0
    lmnc=0.d0
  endif

! Convert half-mesh to full mesh for lambda
! added by Christopher Albert, 2020-02-11
  almnc(:,0) = 0d0
  almns(:,0) = 0d0
  do i=1,kparb-1
    almnc(:,i) = 0.5d0*(lmnc(:,i+1) + lmnc(:,i))
    almns(:,i) = 0.5d0*(lmns(:,i+1) + lmns(:,i))
  enddo
  almnc(:,kparb) = lmnc(:,kparb) &
  + 0.5d0*(lmnc(:,kparb)-lmnc(:,kparb-1))
  almns(:,kparb) = lmns(:,kparb) &
  + 0.5d0*(lmns(:,kparb)-lmns(:,kparb-1))
!
  rmnc=rmnc*fac_r
  zmnc=zmnc*fac_r
  rmns=rmns*fac_r
  zmns=zmns*fac_r
!
  call nc_close(ncid)
!
  end subroutine vmecin
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine stevvo(RT0,R0i,L1i,cbfi,bz0i,bf0)
!
  use new_vmec_stuff_mod, only : nper,rmajor
!
  implicit none
!
  integer :: L1i
  double precision :: RT0,R0i,cbfi,bz0i,bf0
!
  L1i=nper
  RT0=rmajor*1.d2
!
  end subroutine stevvo

end module vmecin_mod
