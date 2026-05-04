! Analytic magnetic field for a large-aspect-ratio circular tokamak.
!
! Concentric circular flux surfaces, q(r) = q0 + q1*(r/a)^2.
! Derived from the poloidal flux function psi(r) = B0 * int_0^r r'/q(r') dr',
! so the field is divergence-free by construction.
!
! Units: lengths in cm, B in Gauss (consistent with the rest of GORILLA).
!
! Field components in cylindrical (R, phi, Z):
!   B_phi = B0*R0/R
!   B_R   = -B0*Z       / (R*q(r))
!   B_Z   = +B0*(R-R0)  / (R*q(r))
! where r = sqrt((R-R0)^2 + Z^2),  q(r) = q0 + q1*(r/a)^2.
!
module field_analytic_circ_mod

  implicit none

  double precision :: R0_ac   ! major radius [cm]
  double precision :: a_ac    ! minor radius [cm]
  double precision :: B0_ac   ! toroidal field at R=R0 [Gauss]
  double precision :: q0_ac   ! safety factor on axis
  double precision :: q1_ac   ! safety factor gradient coefficient

contains

  subroutine set_field_analytic_circ(R0, a, B0, q0, q1)
    use field_eq_mod, only: rtf, btf
    double precision, intent(in) :: R0, a, B0, q0, q1
    R0_ac = R0
    a_ac  = a
    B0_ac = B0
    q0_ac = q0
    q1_ac = q1
    ! Set rtf*btf = B0*R0 so that A_z = -rtf*btf*log(R) gives B_phi = B0*R0/R
    rtf = R0_ac
    btf = B0_ac
  end subroutine set_field_analytic_circ


  subroutine field_analytic_circ(r, p, z,                          &
      Br, Bp, Bz,                                                   &
      dBrdR, dBrdp, dBrdZ, dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)

    use field_eq_mod, only: psif
    double precision, intent(in)  :: r, p, z
    double precision, intent(out) :: Br, Bp, Bz
    double precision, intent(out) :: dBrdR, dBrdp, dBrdZ
    double precision, intent(out) :: dBpdR, dBpdp, dBpdZ
    double precision, intent(out) :: dBzdR, dBzdp, dBzdZ

    double precision :: rho, q, dqdR, dqdZ, Rshift

    Rshift = r - R0_ac                    ! R - R0
    rho    = sqrt(Rshift**2 + z**2)       ! minor radius

    q    = q0_ac + q1_ac*(rho/a_ac)**2
    dqdR = 2.d0*q1_ac*Rshift / a_ac**2   ! dq/dR
    dqdZ = 2.d0*q1_ac*z      / a_ac**2   ! dq/dZ

    ! Poloidal flux: psif = integral_0^rho B0*r'/q(r') dr'
    ! = B0*(a^2/(2*q1))*ln(q(rho)/q0), or B0*rho^2/(2*q0) when q1=0
    if (q1_ac .gt. 0.d0) then
      psif = B0_ac * (a_ac**2 / (2.d0*q1_ac)) * log(q/q0_ac)
    else
      psif = B0_ac * rho**2 / (2.d0*q0_ac)
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

end module field_analytic_circ_mod
