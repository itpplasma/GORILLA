module period_mod
  double precision :: per_phi, per_tht
end module period_mod

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
module field_c_mod
  integer :: icall_c=0
  integer :: ntor=16
  integer :: nr,np,nz
  integer :: icftype
  double precision :: rmin,pmin,zmin,rmax,pmax,zmax
end module field_c_mod
!
module field_eq_mod
  logical :: use_fpol = .true.                                      !<=18.12.18
  logical :: skip_read = .false.
  integer :: iaxieq = 0
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
!
module field_mod
  integer          :: icall=0
  integer          :: ipert,iequil
  double precision :: ampl
end module field_mod
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
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! -----------------------------------------------------------------
subroutine read_field_input
  use input_files, only: iunit, gfile, pfile, convexfile, fluxdatapath
  use field_c_mod, only: ntor, icftype
  use field_eq_mod, only: nwindow_r, nwindow_z
  use field_mod, only: ipert, iequil, ampl
  use inthecore_mod, only: cutoff

  open(iunit, file='field_divB0.inp')
  read(iunit,*) ipert        ! 0=eq only, 1=vac, 2=vac+plas no derivatives,
                             ! 3=plas+vac with derivatives
  read(iunit,*) iequil       ! 0=perturbation alone, 1=with equilibrium
  read(iunit,*) ampl         ! amplitude of perturbation, a.u.
  read(iunit,*) ntor         ! number of toroidal harmonics
  read(iunit,*) cutoff       ! inner cutoff in psi/psi_a units
  read(iunit,*) icftype      ! type of coil file
  read(iunit,*) gfile        ! equilibrium file
  read(iunit,*) pfile        ! coil        file
  read(iunit,*) convexfile   ! convex file for stretchcoords
  read(iunit,*) fluxdatapath ! directory with data in flux coord.
  read(iunit,*) nwindow_r    ! widow size for filtering of psi array over R
  read(iunit,*) nwindow_z    ! widow size for filtering of psi array over Z
  close(iunit)
end subroutine read_field_input

subroutine field(r,p,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ   &
                ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
!
  use period_mod
  use field_c_mod, only: icall_c
  use field_mod
  use inthecore_mod, only: incore
!
  implicit none
!
  double precision :: r,p,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ         &
                     ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ
  double precision :: rm,zm,Brc,Bpc,Bzc,dBrdRc,dBrdpc,dBrdZc   &
                     ,dBpdRc,dBpdpc,dBpdZc,dBzdRc,dBzdpc,dBzdZc
!
  if(icall .eq. 0) then
     icall = 1
     call read_field_input
     print *, 'Perturbation field',ipert,'Ampl',ampl
     if(icall_c.eq.-1) ipert=1
  endif

  call stretch_coords(r,z,rm,zm)

  call field_eq(rm,p,zm,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ   &
               ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
!
  if(iequil.eq.0) then
    Br=0.d0
    Bp=0.d0
    Bz=0.d0
    dBrdR=0.d0
    dBrdp=0.d0
    dBrdZ=0.d0
    dBpdR=0.d0
    dBpdp=0.d0
    dBpdZ=0.d0
    dBzdR=0.d0
    dBzdp=0.d0
    dBzdZ=0.d0
  endif
!
  if(ipert.gt.0) then
!
    if(ipert.gt.1) then
      call inthecore(rm,zm)
    else
      incore=-1
    endif
!
! vacuum perturbation coil field:
!
    call field_c(rm,p,zm,Brc,Bpc,Bzc,dBrdRc,dBrdpc,dBrdZc   &
                ,dBpdRc,dBpdpc,dBpdZc,dBzdRc,dBzdpc,dBzdZc)
!
    Br = Br + Brc*ampl
    Bp = Bp + Bpc*ampl
    Bz = Bz + Bzc*ampl
    dBrdR = dBrdR + dBrdRc*ampl
    dBrdp = dBrdp + dBrdpc*ampl
    dBrdZ = dBrdZ + dBrdZc*ampl
    dBpdR = dBpdR + dBpdRc*ampl
    dBpdp = dBpdp + dBpdpc*ampl
    dBpdZ = dBpdZ + dBpdZc*ampl
    dBzdR = dBzdR + dBzdRc*ampl
    dBzdp = dBzdp + dBzdpc*ampl
    dBzdZ = dBzdZ + dBzdZc*ampl
!
    if(incore.gt.-1) then
! perturbation coil field with plasma shielding:
!
      if(ipert.eq.2) then
!
        call field_fourier(rm,p,zm,Brc,Bpc,Bzc,dBrdRc,dBrdpc,dBrdZc           &
                          ,dBpdRc,dBpdpc,dBpdZc,dBzdRc,dBzdpc,dBzdZc)
!
        Br = Br + Brc*ampl
        Bp = Bp + Bpc*ampl
        Bz = Bz + Bzc*ampl
!
      elseif(ipert.eq.3) then
!
        call field_fourier_derivs(rm,p,zm,Brc,Bpc,Bzc,dBrdRc,dBrdpc,dBrdZc    &
                                 ,dBpdRc,dBpdpc,dBpdZc,dBzdRc,dBzdpc,dBzdZc)
!
        Br = Br + Brc*ampl
        Bp = Bp + Bpc*ampl
        Bz = Bz + Bzc*ampl
        dBrdR = dBrdR + dBrdRc*ampl
        dBrdp = dBrdp + dBrdpc*ampl
        dBrdZ = dBrdZ + dBrdZc*ampl
        dBpdR = dBpdR + dBpdRc*ampl
        dBpdp = dBpdp + dBpdpc*ampl
        dBpdZ = dBpdZ + dBpdZc*ampl
        dBzdR = dBzdR + dBzdRc*ampl
        dBzdp = dBzdp + dBzdpc*ampl
        dBzdZ = dBzdZ + dBzdZc*ampl
!
      endif
!
    endif
!
   end if
!
   return
 end subroutine field
! ========================================================================
subroutine field_eq(r,ppp,z,Brad,Bphi,Bzet,dBrdR,dBrdp,dBrdZ  &
                   ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
!
  use input_files
  use field_eq_mod
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
    if (.not. skip_read) then
        if (iaxieq .eq. 1) then
           call read_dimeq_west(nrad, nzet)
        else
!          call read_dimeq0(nrad, nzet)
           call read_dimeq1(nrad, nzet)
        end if
        allocate(rad(nrad), zet(nzet))
        allocate(psi0(nrad, nzet), psi(nrad, nzet))
    end if
    if (use_fpol) then                                                               !<=18.12.18
      if (.not. skip_read) then
         allocate(splfpol(0:5, nrad))                                                !<=18.12.18
!        call read_eqfile0(nrad, nzet, psib, btf, rtf, rad, zet, psi)
         call read_eqfile2(nrad, nzet, psi_axis, psi_sep, btf, rtf, &                !<=18.12.18
              splfpol(0, :), rad, zet, psi)                                          !<=18.12.18
      end if
      psib=-psi_axis                                                                 !<=18.12.18
      psi_sep=(psi_sep-psi_axis)*1.d8                                                !<=18.12.18
      splfpol(0,:)=splfpol(0,:)*1.d6                                                 !<=18.12.18
      call spline_fpol                                                               !<=18.12.18
    else if (.not. skip_read) then                                                   !<=18.12.18
      if (iaxieq .eq. 1) then
         call read_eqfile_west(nrad, nzet, psib, btf, rtf, rad, zet, psi)
      else
         call read_eqfile1(nrad, nzet, psib, btf, rtf, rad, zet, psi)
      end if
    end if                                                                           !<=18.12.18
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

! ----------- Runov's Original Method --------------------------------
subroutine read_dimeq0(nrad,nzet)
  use input_files
  integer :: nrad, nzet

     open(11,file=eqfile)
     read(11,*)
     read(11,*)
     read(11,*)
     read(11,*)
     read(11,*)
     read(11,*)
     read(11,*)
     read(11,*)
     read(11,*)
     read(11,*)

     read(11,111) nrad
     read(11,111) nzet
111  format(12x,i3)

     close(11)
  return
end subroutine read_dimeq0

subroutine read_eqfile0(nrad, nzet, psib, btf, rtf, rad, zet, psi)
  use input_files
  integer :: nrad, nzet, dum
  real(kind=8) :: psib, btf, rtf
  real(kind=8) :: rad(nrad), zet(nzet)
  real(kind=8) :: psi(nrad,nzet)

     open(11,file=eqfile)
     read(11,*)
     read(11,*)
     read(11,*)
     read(11,*)
     read(11,*)
     read(11,*)
     read(11,*)
     read(11,*)
     read(11,*)
     read(11,*)

     read(11,111) dum !nrad
     read(11,111) dum !nzet
     read(11,112) psib
     read(11,112) btf
     read(11,112) rtf
     read(11,*)
     read(11,*)

     print *, nrad, nzet, psib, btf, rtf

!     nrad = nrad - 3
!     nzet = nzet - 4

!!$     read(11,113)dummy,dummy,(rad(i),i=1,nrad)
!!$     read(11,*)
!!$     read(11,*)
!!$     read(11,113)dummy,dummy,(zet(i),i=1,nzet)
!!$     read(11,*)
!!$     read(11,*)
!!$     read(11,113)(dummy,dummy,(psi(j,k),j=1,nrad),dummy,k=1,2)
!!$     read(11,113)(dummy,dummy,(psi(j,k),j=1,nrad),dummy,k=1,nzet)

     read(11,113)(rad(i),i=1,nrad)
     read(11,*)
     read(11,*)
     read(11,113)(zet(i),i=1,nzet)
     read(11,*)
     read(11,*)
     read(11,113)((psi(j,k),j=1,nrad),k=1,nzet)

!!$     do k=1,nzet
!!$        write(41,*)(psi(j,k),j=1,nrad)
!!$     enddo
    close(11)
    return

111  format(12x,i3)
112  format(12x,f21.2)
113  format(5(e17.4))
end subroutine read_eqfile0


! ----------- Read gfile directly --------------------------------
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

! ===========================================================================

! Input of axisymmetric equilibrium for WEST tokamak

  subroutine read_dimeq_west(nrad, nzet)

  use input_files, only: iunit, gfile

  implicit none

  integer :: nrad, nzet

  open(unit = iunit,file = trim(gfile), status = 'old', action = 'read')
  read(iunit, *) nrad, nzet
  close(iunit)

  end subroutine read_dimeq_west

!----------------------------------------------------------------------

  subroutine read_eqfile_west(nrad, nzet, psib, btf, rtf, rad, zet, psi)

  use input_files, only: iunit, gfile

  implicit none

  integer :: nrad, nzet, ir
  double precision :: psib, btf, rtf
  double precision :: rad(nrad), zet(nzet)
  double precision :: psi(nrad, nzet)

  psib = 0.0d0

  open(unit = iunit, file = trim(gfile), status = 'old', action = 'read')
  read (iunit, *) nrad, nzet
  read (iunit, *) btf
  read (iunit, *) rad
  read (iunit, *) zet
  do ir = 1, nrad
    read (iunit, *) psi(ir, :)
  end do
  close(iunit)
  rtf = 0.5d0 * (rad(1) + rad(nrad))
  btf = btf / rtf

  end subroutine read_eqfile_west

! ===========================================================================
!
subroutine field_c(rrr,ppp,zzz,Brad,Bphi,Bzet,dBrdR,dBrdp,dBrdZ  &
                  ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
!
  use input_files
  use field_c_mod
!
  implicit none
!
  double precision, parameter :: pi=3.14159265358979d0
!
  double precision :: rrr,ppp,zzz,Brad,Bphi,Bzet,dBrdR,dBrdp,dBrdZ  &
                     ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ
  double precision :: hrm1,hpm1,hzm1
  double precision, dimension(:,:,:), allocatable :: Br,Bp,Bz
!
!-------first call: read data from disk-------------------------------
  if(icall_c .lt. 1) then
    print *,'coils: file type = ',icftype
    if(icftype.eq.1) then
      nr=129  !64
      np=37   !37
      nz=129  !64
    elseif(icftype.eq.2) then
      nr=129
      np=33
      nz=129
    elseif(icftype.eq.3) then
      nr=129
      np=37
      nz=131
      icftype=1
    elseif(icftype.eq.4) then
      call read_sizes(nr,np,nz)
    else
      print *,'field_c: unknown coil file type'
      stop
    endif
    allocate(Br(nr,np,nz),Bp(nr,np,nz),Bz(nr,np,nz))
!
!    call read_field0(rad,phi,zet,rmin,pmin,zmin,hrm1,hpm1,hzm1,Br,Bp,Bz)
!    call read_field1(icftype,nr,np,nz,rmin,rmax,pmin,pmax,zmin,zmax,Br,Bp,Bz)
!
!
    if(icftype.lt.4) then
      call read_field2(icftype,nr,np,nz,rmin,rmax,pmin,pmax,zmin,zmax,Br,Bp,Bz)
    else
      call read_field4(nr,np,nz,rmin,rmax,pmin,pmax,zmin,zmax,Br,Bp,Bz)
    endif
!
    print *,'coils: nr,np,nz = ',nr,np,nz
    print *,'coils: rmin,rmax = ',rmin,rmax
    print *,'coils: zmin,zmax = ',zmin,zmax
    print *,'coils: pmin,pmax = ',pmin,pmax
!
    call vector_potentials(nr,np,nz,ntor,rmin,rmax,pmin,pmax,zmin,zmax,br,bp,bz)
!
    deallocate(Br,Bp,Bz)
!
    if(icall_c.eq.-1) then
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
      icall_c = 1
      return
    endif
    icall_c = 1
  endif
  !------- end first call ----------------------------------------------
!
print *,'field_divfree'
  call field_divfree(rrr,ppp,zzz,Brad,Bphi,Bzet,dBrdR,dBrdp,dBrdZ    &
                    ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
print *,'end field_divfree'
!
  return
end subroutine field_c



! ===========================================================================
subroutine read_field0(rad,phi,zet,rmin,pmin,zmin,hrm1,hpm1,hzm1,Br,Bp,Bz)
!
  use input_files
  parameter(nr=64,np=37,nz=64)
  real, parameter :: pi=3.14159265358979d0
  parameter (mp=4) ! power of Lagrange's polynomial =3
  dimension Bz(nr,np,nz)
  dimension Br(nr,np,nz),Bp(nr,np,nz)
  dimension rad(nr), phi(np), zet(nz)
  dimension xp(mp),yp(mp),zp(mp),fp(mp,mp,mp)
  integer indx(mp), indy(mp), indz(mp)
  data icall/0/
  save
!
!-------first call: read data from disk-------------------------------
     open(1,file=cfile,status='old',action='read')
     read(1,*)
     read(1,*)
     read(1,*)
     read(1,*)
     read(1,*)

!---Input B      -->T = V*s/m/m
     do j=1,np-1	 !only npmax-1 points are given
        do k=nz,1,-1  !reverse order of probe data
           do i=1,nr
              read(1,*) Br(i,j,k), Bp(i,j,k), Bz(i,j,k)

              Br(i,j,k) = Br(i,j,k)*1.d4
              Bp(i,j,k) = Bp(i,j,k)*1.d4
              Bz(i,j,k) = Bz(i,j,k)*1.d4

           enddo
           read(1,*)
        enddo
        read(1,*)
     enddo
     close(1)
     !
     rmin = 84.
     rmax = 254.
     zmin = -160.
     zmax = 160.
     pmin = 0.
     pmax = 2.*pi


     hrad = (rmax - rmin)/(nr-1)
     hphi = (pmax - pmin)/(np-1)
     hzet = (zmax - zmin)/(nz-1)

     do i=1,nr
        rad(i) = rmin + hrad*(i-1)
     enddo
     do i=1,np
        phi(i) = pmin + hphi*(i-1)
     enddo
     do i=1,nz
        zet(i) = zmin + hzet*(i-1)
     enddo

     do i=1,nr
        do k=1,nz
           Br(i,np,k) = Br(i,1,k)
           Bp(i,np,k) = Bp(i,1,k)
           Bz(i,np,k) = Bz(i,1,k)
        enddo
     enddo
end subroutine read_field0
!
subroutine read_field1(icftype,nr,np,nz,rmin,rmax,pmin,pmax,zmin,zmax,Br,Bp,Bz)
!
  use input_files
!
  implicit none
!
  double precision, parameter :: pi=3.14159265358979d0
!
  integer :: nr,np,nz,i,j,k,icftype
  double precision :: rmin,pmin,zmin,rmax,pmax,zmax,xdim,zdim,zmid,dum
  double precision, dimension(nr,np,nz) :: Br,Bp,Bz
!
  open(iunit,file=trim(pfile),status='old',action='read')
  read(iunit,*)
  read(iunit,*)
  read(iunit,*)
  read(iunit,*)  !PROBE
  read(iunit,*)
  read(iunit,*)
  if(icftype.eq.2) then
    read(iunit,*)    !New Format
    read(iunit,*)    !New Format
  endif

!---Input B      -->T = V*s/m/m
  do j=1,np-1   !only npmax-1 points are given
     do k=nz,1,-1  !reverse order of probe data
        do i=1,nr
           if(icftype.eq.1) then
!					Old Format
             read(iunit,*) Br(i,j,k),Bp(i,j,k),Bz(i,j,k)
           elseif(icftype.eq.2) then
!					New Format
             read(iunit,*) dum,dum,dum,Br(i,j,k),Bp(i,j,k),Bz(i,j,k),dum,dum
           endif
!
				  !Convert to CGS
           Br(i,j,k) = Br(i,j,k)*1.d4
           Bp(i,j,k) = Bp(i,j,k)*1.d4
           Bz(i,j,k) = Bz(i,j,k)*1.d4
        enddo
        read(iunit,*)
     enddo
     read(iunit,*)
  enddo
  close(iunit)
!
  xdim=170.d0
  rmin=84.d0
  rmax=rmin+xdim
!
  pmin = 0.
  pmax = 2.*pi
!
  zdim=320.d0
  zmid=0.d0
  zmin=zmid - zdim/2.d0
  zmax=zmid + zdim/2.d0
!
  do i=1,nr
     do k=1,nz
        Br(i,np,k) = Br(i,1,k)
        Bp(i,np,k) = Bp(i,1,k)
        Bz(i,np,k) = Bz(i,1,k)
     enddo
  enddo
end subroutine read_field1
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
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
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine inthecore(R,Z)
!
  use inthecore_mod
  use input_files,  only : iunit,fluxdatapath
  use field_eq_mod, only : psif,dpsidr,dpsidz,d2psidr2,d2psidrdz,d2psidz2
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
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine read_field2(icftype,nr,np,nz,rmin,rmax,pmin,pmax,zmin,zmax,Br,Bp,Bz)
!
  use input_files
!
  implicit none
!
  double precision, parameter :: pi=3.14159265358979d0
!
  integer :: nr,np,nz,i,j,k,icftype
  double precision :: rmin,pmin,zmin,rmax,pmax,zmax,xdim,zdim,zmid,dum
  double precision, dimension(nr,np,nz) :: Br,Bp,Bz
!
  open(iunit,file=trim(pfile),status='old',action='read')

!---Input B      -->T = V*s/m/m
  do j=1,np-1   !only npmax-1 points are given
     do k=nz,1,-1  !reverse order of probe data
        do i=1,nr
           if(icftype.eq.1) then
!                                       Old Format
             read(iunit,*) Br(i,j,k),Bp(i,j,k),Bz(i,j,k)
           elseif(icftype.eq.2) then
!                                       New Format
             read(iunit,*) dum,dum,dum,Br(i,j,k),Bp(i,j,k),Bz(i,j,k),dum,dum
           endif
!
                                  !Convert to CGS
           Br(i,j,k) = Br(i,j,k)*1.d4
           Bp(i,j,k) = Bp(i,j,k)*1.d4
           Bz(i,j,k) = Bz(i,j,k)*1.d4
        enddo
     enddo
  enddo
  close(iunit)
!
  xdim=300.d0
  rmin=100.d0
  rmax=rmin+xdim
!
  pmin = 0.
  pmax = 2.*pi
!
  zdim=400.d0
  zmid=0.d0
  zmin=zmid - zdim/2.d0
  zmax=zmid + zdim/2.d0
!
  do i=1,nr
     do k=1,nz
        Br(i,np,k) = Br(i,1,k)
        Bp(i,np,k) = Bp(i,1,k)
        Bz(i,np,k) = Bz(i,1,k)
     enddo
  enddo
end subroutine read_field2
!
subroutine read_sizes(nr,np,nz)
!
  use input_files, only : iunit,pfile
!
  implicit none
  integer :: nr,np,nz
!
  open(iunit,file=pfile)
  read(iunit,*) nr,np,nz
  close(iunit)
!
end subroutine read_sizes
!
subroutine read_field4(nr,np,nz,rmin,rmax,pmin,pmax,zmin,zmax,Br,Bp,Bz)
!
  use input_files, only : iunit,pfile
!
  implicit none
!
  integer :: nr,np,nz,i,j,k
  double precision :: rmin,rmax,pmin,pmax,zmin,zmax
  double precision, dimension(nr,np,nz)       :: Br,Bp,Bz
!
  open(iunit,file=pfile)
  read(iunit,*) nr,np,nz
  read(iunit,*) rmin,rmax
  read(iunit,*) pmin,pmax
  read(iunit,*) zmin,zmax
  do i=1,nr
     do j=1,np
        do k=1,nz
           read(iunit,*) Br(i,j,k),Bp(i,j,k),Bz(i,j,k)
        enddo
     enddo
  enddo
  close(iunit)
!
end subroutine read_field4

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine spline_fpol
!
  use field_eq_mod, only : nrad,hfpol,splfpol
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
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
