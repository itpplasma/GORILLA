module period_mod
  double precision :: per_phi, per_tht
end module period_mod
!
module field_c_mod
  integer :: icall_c=0
  integer :: ntor=16
  integer :: nr,np,nz
  integer :: icftype
  double precision :: rmin,pmin,zmin,rmax,pmax,zmax
end module field_c_mod
!
module field_divB0_mod

  implicit none

  contains

subroutine field(r,p,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ   &
                ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
!
  use period_mod
  use input_files
  use field_c_mod,   only : ntor,icftype,icall_c
  use field_mod
  use inthecore_mod, only : incore,cutoff
  use field_eq_mod, only : nwindow_r,nwindow_z
  use tetra_grid_settings_mod, only: g_file_filename, convex_wall_filename, iaxieq_in !=> Michael Eder, 05 May 2022
  use utils_field_divB0_mod, only: field_fourier, field_fourier_derivs, inthecore, stretch_coords, field_eq
!
  implicit none
!
  double precision :: r,p,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ         &
                     ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ
  double precision :: rm,zm,Brc,Bpc,Bzc,dBrdRc,dBrdpc,dBrdZc   &
                     ,dBpdRc,dBpdpc,dBpdZc,dBzdRc,dBzdpc,dBzdZc
  character*1024 :: dummy
  ! logical :: boole_reverse_field = .false.

!
  if(icall .eq. 0) then
     icall = 1
     open(iunit, file='field_divB0.inp')
     read(iunit,*) ipert        ! 0=eq only, 1=vac, 2=vac+plas no derivatives,
                                ! 3=plas+vac with derivatives
     read(iunit,*) iequil       ! 0=perturbation alone, 1=with equilibrium
     read(iunit,*) ampl         ! amplitude of perturbation, a.u.
     read(iunit,*) ntor         ! number of toroidal harmonics
     read(iunit,*) cutoff       ! inner cutoff in psi/psi_a units
     read(iunit,*) icftype      ! type of coil file
     read(iunit,*) dummy !gfile        ! equilibrium file => Michael Eder, 08 March 2021
     read(iunit,*) pfile        ! coil        file
     read(iunit,*) dummy !convexfile   ! convex file for stretchcoords => Michael Eder, 08 March 2021
     read(iunit,*) fluxdatapath ! directory with data in flux coord.
     read(iunit,*) nwindow_r    ! widow size for filtering of psi array over R
     read(iunit,*) nwindow_z    ! widow size for filtering of psi array over Z
     read(iunit,*,err=7) dummy !iaxieq ! format of axisymmetric equilibrium file, 0 for EFIT (default), 1 for WEST => Michael Eder, 05 May 2022
7    close(iunit)
     print *, 'Perturbation field',ipert,'Ampl',ampl
     if(icall_c.eq.-1) ipert=1

     !Use input data from tetra_grid_settings_mod => Michael Eder, 05 May 2022
     gfile = g_file_filename
     convexfile = convex_wall_filename
     iaxieq = iaxieq_in
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
  !  if (boole_reverse_field) then
  !   Brc = -Brc
  !   Bpc = -Bpc
  !   Bzc = -Bzc
  !   dBrdRc = -dBrdRc
  !   dBrdpc = -dBrdpc
  !   dBrdZc = -dBrdZc
  !   dBpdRc = -dBpdRc
  !   dBpdpc = -dBpdpc
  !   dBpdZc = -dBpdZc
  !   dBzdRc = -dBzdRc
  !   dBzdpc = -dBzdpc
  !   dBzdZc = -dBzdZc
  !  endif
!
   return
 end subroutine field


! ----------- Read gfile directly --------------------------------

! ===========================================================================
!
subroutine field_c(rrr,ppp,zzz,Brad,Bphi,Bzet,dBrdR,dBrdp,dBrdZ  &
                  ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
!
  use input_files
  use field_c_mod
  use utils_field_divB0_mod, only: field_divfree, vector_potentials
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
  call field_divfree(rrr,ppp,zzz,Brad,Bphi,Bzet,dBrdR,dBrdp,dBrdZ    &
                    ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
!
  return
end subroutine field_c
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

!!! These subroutines are currently unused

! ========================================================================

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
  integer :: nrad, nzet, dum, i, j, k
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

end module field_divB0_mod