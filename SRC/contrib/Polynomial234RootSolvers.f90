!! ----------------------------------------------------------------------------------
!! MODULE Polynomial234RootSolvers
!!
!!    Minimal in-tree replacement for the ACM TOMS Algorithm 954 polynomial root
!!    solver.  Implements quadraticRoots and cubicRoots with classical
!!    (numerically-stable) formulas; quarticRoots still aborts with a clear
!!    error so that exotic code paths that rely on the original library remain
!!    visible.  For high-precision benchmarks the real Algorithm 954 file should
!!    be dropped in to replace this one.
!!
!!    Reference equations:
!!      quadraticRoots:   x**2 + q1*x + q0 = 0
!!      cubicRoots:       x**3 + c2*x**2 + c1*x + c0 = 0
!!      quarticRoots:     x**4 + q3*x**3 + q2*x**2 + q1*x + q0 = 0  (NOT IMPLEMENTED)
!!
!! ----------------------------------------------------------------------------------

module Polynomial234RootSolvers

use SetWorkingPrecision, ONLY : wp, wpformat

implicit none

real(wp), parameter, private :: pi_ = 3.141592653589793238462643383279502884_wp
real(wp), parameter, private :: third_ = 1.0_wp / 3.0_wp
real(wp), parameter, private :: two_pi_over_three_ = 2.0_wp * pi_ / 3.0_wp

contains

! ---------------------------------------------------------------------
! Quadratic: x**2 + q1*x + q0 = 0
! Uses the numerically-stable form
!     q = -0.5*(q1 + sign(q1)*sqrt(discr))
!     root1 = q;  root2 = q0 / q
! to avoid cancellation when q1 and sqrt(discr) are close in magnitude.
! ---------------------------------------------------------------------
subroutine quadraticRoots (q1, q0, nReal, root)

  implicit none

  integer            , intent (out) :: nReal
  real    (kind = wp), intent (in)  :: q1, q0
  real    (kind = wp), intent (out) :: root (1:2,1:2)

  real(wp) :: discr, sqrt_discr, qhelp

  root = 0.0_wp
  discr = q1*q1 - 4.0_wp * q0

  if (discr > 0.0_wp) then
    nReal = 2
    sqrt_discr = sqrt(discr)
    if (q1 >= 0.0_wp) then
      qhelp = -0.5_wp * (q1 + sqrt_discr)
    else
      qhelp = -0.5_wp * (q1 - sqrt_discr)
    end if
    root(1, 1) = qhelp
    root(1, 2) = 0.0_wp
    if (qhelp /= 0.0_wp) then
      root(2, 1) = q0 / qhelp
    else
      root(2, 1) = 0.0_wp
    end if
    root(2, 2) = 0.0_wp
  else if (discr == 0.0_wp) then
    nReal = 2
    root(1, 1) = -0.5_wp * q1
    root(1, 2) = 0.0_wp
    root(2, 1) = -0.5_wp * q1
    root(2, 2) = 0.0_wp
  else
    nReal = 0
    sqrt_discr = sqrt(-discr)
    root(1, 1) = -0.5_wp * q1
    root(1, 2) =  0.5_wp * sqrt_discr
    root(2, 1) = -0.5_wp * q1
    root(2, 2) = -0.5_wp * sqrt_discr
  end if

end subroutine quadraticRoots

! ---------------------------------------------------------------------
! Cubic: x**3 + c2*x**2 + c1*x + c0 = 0
!   1. depress with y = x + c2/3:  y**3 + p*y + q = 0
!   2. branch on sign of discriminant disc = q**2/4 + p**3/27
!        disc > 0: one real, two complex conjugate (Cardano)
!        disc = 0: three real (some equal)
!        disc < 0: three distinct real (trigonometric)
!   3. shift back x = y - c2/3.
! ---------------------------------------------------------------------
subroutine cubicRoots (c2, c1, c0, nReal, root, printInfo)

  implicit none

  logical, optional  , intent (in)  :: printInfo
  integer            , intent (out) :: nReal
  real    (kind = wp), intent (in)  :: c2, c1, c0
  real    (kind = wp), intent (out) :: root (1:3,1:2)

  real(wp) :: shift, p, q, disc
  real(wp) :: A_, B_, t1, t2
  real(wp) :: m_, theta, cos_arg

  root = 0.0_wp
  shift = c2 * third_
  p = c1 - c2*c2 * third_
  q = (2.0_wp * c2*c2*c2) / 27.0_wp - (c2 * c1) * third_ + c0
  disc = (q*q) * 0.25_wp + (p*p*p) / 27.0_wp

  if (disc > 0.0_wp) then
    ! one real, two complex conjugate
    nReal = 1
    A_ = -q * 0.5_wp + sqrt(disc)
    B_ = -q * 0.5_wp - sqrt(disc)
    t1 = sign(abs(A_)**third_, A_)
    t2 = sign(abs(B_)**third_, B_)
    root(1, 1) = t1 + t2 - shift
    root(1, 2) = 0.0_wp
    root(2, 1) = -0.5_wp * (t1 + t2) - shift
    root(2, 2) = 0.5_wp * sqrt(3.0_wp) * (t1 - t2)
    root(3, 1) = root(2, 1)
    root(3, 2) = -root(2, 2)
  else if (disc == 0.0_wp) then
    nReal = 3
    if (q == 0.0_wp) then
      root(1, 1) = -shift
      root(2, 1) = -shift
      root(3, 1) = -shift
    else
      t1 = sign(abs(q * 0.5_wp)**third_, q)
      root(1, 1) = -2.0_wp * t1 - shift
      root(2, 1) =  t1 - shift
      root(3, 1) =  t1 - shift
    end if
  else
    ! three distinct real roots, trigonometric form (p < 0 here)
    nReal = 3
    m_ = sqrt(-p * third_)
    cos_arg = -q / (2.0_wp * m_ * m_ * m_)
    cos_arg = max(-1.0_wp, min(1.0_wp, cos_arg))
    theta = acos(cos_arg) * third_
    root(1, 1) = 2.0_wp * m_ * cos(theta)                          - shift
    root(2, 1) = 2.0_wp * m_ * cos(theta - two_pi_over_three_)     - shift
    root(3, 1) = 2.0_wp * m_ * cos(theta + two_pi_over_three_)     - shift
  end if

  if (present(printInfo)) then
    if (printInfo) then
      write(*, '(a, 3es16.8)') 'cubicRoots: real parts = ', root(1,1), root(2,1), root(3,1)
    end if
  end if

end subroutine cubicRoots

! ---------------------------------------------------------------------
! Quartic: x**4 + q3*x**3 + q2*x**2 + q1*x + q0 = 0
! NOT implemented in this in-tree replacement.  Calls are caught at
! runtime with a clear error so the missing path is obvious.
! ---------------------------------------------------------------------
subroutine quarticRoots (q3, q2, q1, q0, nReal, root, printInfo)

  implicit none

  logical, optional  , intent (in)  :: printInfo
  integer            , intent (out) :: nReal
  real    (kind = wp), intent (in)  :: q3, q2, q1, q0
  real    (kind = wp), intent (out) :: root (1:4,1:2)

  print *, 'ERROR: quarticRoots is not implemented in the in-tree Polynomial234RootSolvers.'
  print *, '       Replace SRC/contrib/Polynomial234RootSolvers.f90 with the real ACM TOMS'
  print *, '       Algorithm 954 file from <https://doi.org/10.1145/2699468>.'
  print *, '       coefficients: q3=', q3, ' q2=', q2, ' q1=', q1, ' q0=', q0
  nReal = 0
  root  = 0.0_wp
  stop

end subroutine quarticRoots

end module Polynomial234RootSolvers
