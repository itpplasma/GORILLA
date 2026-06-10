! Analytic magnetic field for a large-aspect-ratio circular tokamak.
!
! Concentric circular flux surfaces.  Two q-profile modes:
!   (a) analytic   q(r) = q0 + q1*(r/a)^p   (p = 2 or 4), closed-form psi.
!   (b) tabulated  q(r) read from a 2-column file (r[cm], q), cubic-splined,
!                  with psi_pol = B0*int_0^r r'/q dr' integrated numerically
!                  and splined.  Lets GORILLA use KIM's exact q-profile.
! Field is divergence-free by construction (psi = B0 int r/q dr).
!
! Field components in cylindrical (R, phi, Z):
!   B_phi = B0*R0/R
!   B_R   = -B0*Z       / (R*q(r))
!   B_Z   = +B0*(R-R0)  / (R*q(r))
! where r = sqrt((R-R0)^2 + Z^2).
!
module field_analytic_circ_mod

  implicit none

  double precision :: R0_ac   ! major radius [cm]
  double precision :: a_ac    ! minor radius [cm]
  double precision :: B0_ac   ! toroidal field at R=R0 [Gauss]
  double precision :: q0_ac   ! safety factor on axis (analytic mode)
  double precision :: q1_ac   ! safety factor gradient coefficient (analytic mode)
  integer          :: p_ac    ! q-profile exponent: q = q0 + q1*(r/a)^p (2 or 4)

  ! --- tabulated q-profile state (mode b) ---
  logical                       :: tabulated_ac = .false.
  integer                       :: ntab_ac = 0
  double precision, allocatable :: rt_ac(:)    ! minor radius grid [cm]
  double precision, allocatable :: qt_ac(:)    ! q on the grid
  double precision, allocatable :: qy2_ac(:)   ! q spline 2nd derivs
  double precision, allocatable :: pt_ac(:)    ! psi_pol on the grid
  double precision, allocatable :: py2_ac(:)   ! psi spline 2nd derivs

contains

  subroutine set_field_analytic_circ(R0, a, B0, q0, q1, q_exp, q_file)
    double precision, intent(in) :: R0, a, B0, q0, q1
    integer,          intent(in) :: q_exp
    character(len=*), intent(in), optional :: q_file
    R0_ac = R0
    a_ac  = a
    B0_ac = B0
    q0_ac = q0
    q1_ac = q1
    p_ac  = q_exp
    tabulated_ac = .false.
    if (present(q_file)) then
      if (len_trim(q_file) > 0) call load_q_table(q_file)
    end if
  end subroutine set_field_analytic_circ


  subroutine load_q_table(fname)
    ! Read (r[cm], q) table, spline q(r), integrate psi_pol = B0*int_0^r r'/q dr'
    ! (trapezoidal on the table grid; psi(r1) seeded assuming q~const below r1),
    ! and spline psi(r).  Sets tabulated_ac = .true.
    character(len=*), intent(in) :: fname
    integer :: iu, ios, n, i
    character(len=512) :: line
    double precision :: rr, qq

    open(newunit=iu, file=fname, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      print *, 'ERROR field_analytic_circ: cannot open q_profile_file: ', trim(fname); stop
    end if
    n = 0
    do
      read(iu,'(A)',iostat=ios) line
      if (ios /= 0) exit
      line = adjustl(line)
      if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
      n = n + 1
    end do
    rewind(iu)
    ntab_ac = n
    allocate(rt_ac(n), qt_ac(n), qy2_ac(n), pt_ac(n), py2_ac(n))
    i = 0
    do
      read(iu,'(A)',iostat=ios) line
      if (ios /= 0) exit
      line = adjustl(line)
      if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
      i = i + 1
      read(line,*) rr, qq
      rt_ac(i) = rr
      qt_ac(i) = abs(qq)            ! use |q|; sign handled by mode/field convention
    end do
    close(iu)

    call spline_setup(rt_ac, qt_ac, ntab_ac, qy2_ac)

    ! psi_pol = B0 * int_0^r r'/q dr' ; seed psi(r1) ~ B0*r1^2/(2 q(r1))
    pt_ac(1) = B0_ac * rt_ac(1)**2 / (2.d0*qt_ac(1))
    do i = 2, ntab_ac
      pt_ac(i) = pt_ac(i-1) + B0_ac * 0.5d0 * &
                 (rt_ac(i)/qt_ac(i) + rt_ac(i-1)/qt_ac(i-1)) * (rt_ac(i)-rt_ac(i-1))
    end do
    call spline_setup(rt_ac, pt_ac, ntab_ac, py2_ac)

    tabulated_ac = .true.
    print *, '[field_analytic_circ] tabulated q-profile loaded: ', trim(fname), &
             '  npts=', ntab_ac
  end subroutine load_q_table


  double precision function psi_tabulated(rho)
    double precision, intent(in) :: rho
    psi_tabulated = splint_val(rt_ac, pt_ac, py2_ac, ntab_ac, rho)
  end function psi_tabulated


  ! Safety factor q at minor radius rho (tabulated spline or analytic q0+q1(r/a)^p).
  double precision function q_analytic_circ(rho) result(q)
    double precision, intent(in) :: rho
    if (tabulated_ac) then
      q = splint_val(rt_ac, qt_ac, qy2_ac, ntab_ac, rho)
    else
      q = q0_ac + q1_ac*(rho/a_ac)**p_ac
    end if
  end function q_analytic_circ


  ! Poloidal flux psi_pol(rho) for the analytic circular tokamak.  Single source
  ! of truth used by BOTH field_divB0 (-> tetra_physics%Aphi1) and the internal
  ! flux_functions writer, so the equil-mapping is consistent with Aphi1 by
  ! construction.  Tabulated -> psi_tabulated; analytic -> closed-form parabola
  ! (p=2) or arctan (p=4) of psi_pol = B0 int_0^r r'/q dr'.
  double precision function psi_pol_analytic_circ(rho) result(psi)
    double precision, intent(in) :: rho
    double precision :: q_loc
    if (tabulated_ac) then
      psi = psi_tabulated(rho)
    else if (p_ac .eq. 2) then
      q_loc = q0_ac + q1_ac*(rho/a_ac)**2
      if (q1_ac .ne. 0.d0 .and. abs(q0_ac) .gt. 1.d-10 .and. q_loc/q0_ac .gt. 0.d0) then
        psi = B0_ac * a_ac**2 / (2.d0*q1_ac) * log(q_loc/q0_ac)
      else
        psi = 0.d0
      end if
    else
      psi = B0_ac * a_ac**2 / (2.d0*sqrt(q0_ac*q1_ac)) &
            * atan( (rho/a_ac)**2 * sqrt(q1_ac/q0_ac) )
    end if
  end function psi_pol_analytic_circ


  ! Write a SELF-CONSISTENT flux_functions.dat (R_beg, r, q, psi_pol, psi_tor)
  ! for grid_kind=5 directly from the analytic field -- so grid_kind=5 never
  ! depends on an externally provided (and possibly inconsistent) equil-mapping.
  ! Mirrors preload_for_SYNCH for EFIT/grid_kind=2.  Columns:
  !   r       : minor radius grid [0, a]
  !   s       : geometric flux label (R0-sqrt(R0^2-r^2))/s_edge, encoded in psi_tor
  !             so s_grid = psi_tor/psi_tor(edge) == the mesh eval_s_local s
  !   psi_pol : = psi_pol_analytic_circ(r) == field_divB0 psif == tetra_physics%Aphi1
  !   q       : field q(r)
  subroutine write_flux_functions_analytic_circ()
    integer, parameter :: nf = 1001
    integer :: iu, i
    double precision :: r, s_edge, s_geom, psi_tor_scale
    s_edge = R0_ac - sqrt(R0_ac**2 - a_ac**2)
    psi_tor_scale = B0_ac * a_ac**2 / 2.d0
    open(newunit=iu, file='flux_functions.dat', status='replace', action='write')
    write(iu,'(A)') ' # R_beg, r,  q, psi_pol, psi_tor'
    do i = 1, nf
      r = dble(i-1)/dble(nf-1) * a_ac
      if (i .eq. 1) r = 1.d-6
      s_geom = (R0_ac - sqrt(R0_ac**2 - r**2)) / s_edge
      write(iu,'(5(1x,ES23.16))') R0_ac + r, r, q_analytic_circ(r), &
            psi_pol_analytic_circ(r), psi_tor_scale * s_geom
    end do
    close(iu)
    print *, '[field_analytic_circ] wrote self-consistent flux_functions.dat, nf=', nf
  end subroutine write_flux_functions_analytic_circ


  subroutine field_analytic_circ(r, p, z,                          &
      Br, Bp, Bz,                                                   &
      dBrdR, dBrdp, dBrdZ, dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)

    double precision, intent(in)  :: r, p, z
    double precision, intent(out) :: Br, Bp, Bz
    double precision, intent(out) :: dBrdR, dBrdp, dBrdZ
    double precision, intent(out) :: dBpdR, dBpdp, dBpdZ
    double precision, intent(out) :: dBzdR, dBzdp, dBzdZ

    double precision :: rho, q, dqdR, dqdZ, Rshift, dqdrho

    Rshift = r - R0_ac                    ! R - R0
    rho    = sqrt(Rshift**2 + z**2)       ! minor radius

    if (tabulated_ac) then
      q      = splint_val(rt_ac, qt_ac, qy2_ac, ntab_ac, rho)
      dqdrho = splint_der(rt_ac, qt_ac, qy2_ac, ntab_ac, rho)
      dqdR   = dqdrho * Rshift / rho      ! drho/dR = Rshift/rho (rho>0 in annulus)
      dqdZ   = dqdrho * z      / rho
    else
      q    = q0_ac + q1_ac*(rho/a_ac)**p_ac
      dqdR = q1_ac*dble(p_ac)*rho**(p_ac-2)*Rshift / a_ac**p_ac   ! dq/dR
      dqdZ = q1_ac*dble(p_ac)*rho**(p_ac-2)*z      / a_ac**p_ac   ! dq/dZ
    end if

    ! --- field components ---
    Bp = B0_ac*R0_ac / r                  ! toroidal
    Br = -B0_ac*z          / (r*q)        ! radial
    Bz = B0_ac*Rshift      / (r*q)        ! vertical

    ! --- phi derivatives vanish (axisymmetric) ---
    dBrdp = 0.d0
    dBpdp = 0.d0
    dBzdp = 0.d0

    ! --- toroidal field derivatives ---
    dBpdR = -B0_ac*R0_ac / r**2
    dBpdZ = 0.d0

    ! --- radial field derivatives ---
    dBrdR = B0_ac*z * (1.d0/(r**2*q) + dqdR/(r*q**2))
    dBrdZ = -B0_ac/(r*q) + B0_ac*z*dqdZ/(r*q**2)

    ! --- vertical field derivatives ---
    dBzdR = B0_ac*(1.d0/(r*q) - Rshift/(r**2*q) - Rshift*dqdR/(r*q**2))
    dBzdZ = -B0_ac*Rshift*dqdZ / (r*q**2)

  end subroutine field_analytic_circ


  ! ---------- self-contained natural cubic spline (value + 1st derivative) ----------
  subroutine spline_setup(x, y, n, y2)
    integer, intent(in) :: n
    double precision, intent(in)  :: x(n), y(n)
    double precision, intent(out) :: y2(n)
    integer :: i, k
    double precision :: sig, pp, u(n)
    y2(1) = 0.d0; u(1) = 0.d0                      ! natural BC
    do i = 2, n-1
      sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
      pp  = sig*y2(i-1) + 2.d0
      y2(i) = (sig-1.d0)/pp
      u(i) = (6.d0*((y(i+1)-y(i))/(x(i+1)-x(i)) - (y(i)-y(i-1))/(x(i)-x(i-1))) &
              /(x(i+1)-x(i-1)) - sig*u(i-1))/pp
    end do
    y2(n) = 0.d0                                    ! natural BC
    do k = n-1, 1, -1
      y2(k) = y2(k)*y2(k+1) + u(k)
    end do
  end subroutine spline_setup

  subroutine locate(x, n, xx, klo)
    integer, intent(in) :: n
    double precision, intent(in) :: x(n), xx
    integer, intent(out) :: klo
    integer :: khi, k
    if (xx <= x(1)) then; klo = 1; return; end if
    if (xx >= x(n)) then; klo = n-1; return; end if
    klo = 1; khi = n
    do while (khi-klo > 1)
      k = (khi+klo)/2
      if (x(k) > xx) then; khi = k; else; klo = k; end if
    end do
  end subroutine locate

  double precision function splint_val(x, y, y2, n, xx)
    integer, intent(in) :: n
    double precision, intent(in) :: x(n), y(n), y2(n), xx
    integer :: klo, khi
    double precision :: h, a, b
    call locate(x, n, xx, klo); khi = klo+1
    h = x(khi)-x(klo); a = (x(khi)-xx)/h; b = (xx-x(klo))/h
    splint_val = a*y(klo) + b*y(khi) + ((a**3-a)*y2(klo) + (b**3-b)*y2(khi))*h*h/6.d0
  end function splint_val

  double precision function splint_der(x, y, y2, n, xx)
    integer, intent(in) :: n
    double precision, intent(in) :: x(n), y(n), y2(n), xx
    integer :: klo, khi
    double precision :: h, a, b
    call locate(x, n, xx, klo); khi = klo+1
    h = x(khi)-x(klo); a = (x(khi)-xx)/h; b = (xx-x(klo))/h
    splint_der = (y(khi)-y(klo))/h - (3.d0*a*a-1.d0)/6.d0*h*y2(klo) &
                 + (3.d0*b*b-1.d0)/6.d0*h*y2(khi)
  end function splint_der

end module field_analytic_circ_mod
