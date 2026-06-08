!! ----------------------------------------------------------------------------------
!! MODULE Polynomial234RootSolvers
!!
!! Wrapper around the Skowron & Gould (2012) general complex polynomial root
!! solver `cmplx_roots_gen` (see SRC/contrib/cmplx_roots_sg.f90), providing the
!! quadratic/cubic/quartic interfaces previously supplied by ACM Algorithm 954
!! (Flocke 2015). External call sites are unchanged: arguments and the
!! `root(:,1)=Re, root(:,2)=Im` layout are preserved, real roots are emitted
!! with an exactly-zero imaginary part, and real roots are placed before
!! complex ones (largest real part first), matching the prior convention.
!!
!! Skowron J., Gould A., 2012, "General Complex Polynomial Root Solver and Its
!! Further Optimization for Binary Microlenses", arXiv:1203.1034.
!! License: Apache 2.0 (see SRC/contrib/LICENSE.cmplx_roots_sg).
!!
!! ----------------------------------------------------------------------------------

module Polynomial234RootSolvers

  use SetWorkingPrecision, only: wp
  use cmplx_roots_sg_mod,  only: cmplx_roots_gen

  implicit none

  private
  public :: quadraticRoots, cubicRoots, quarticRoots

  ! A root with |Im(z)| <= cmplx_tol * max(1, |Re(z)|) is taken to be real
  ! and emitted with Im = 0 exactly. Sits well above double-precision
  ! macheps (~2.2e-16) while still rejecting genuinely complex roots whose
  ! imaginary parts are merely small.
  real(kind = wp), parameter :: cmplx_tol = 1.0e-12_wp

contains

  subroutine quadraticRoots(q1, q0, nReal, root, printInfo)
    real(kind = wp),   intent(in)  :: q1, q0
    integer,           intent(out) :: nReal
    real(kind = wp),   intent(out) :: root(1:2, 1:2)
    logical, optional, intent(in)  :: printInfo

    complex(kind = wp) :: poly(3), croots(2)

    poly(1) = cmplx(q0,     0.0_wp, kind = wp)
    poly(2) = cmplx(q1,     0.0_wp, kind = wp)
    poly(3) = cmplx(1.0_wp, 0.0_wp, kind = wp)

    call cmplx_roots_gen(croots, poly, 2, .true., .false.)
    call pack_roots(croots, 2, nReal, root)
  end subroutine quadraticRoots

  subroutine cubicRoots(c2, c1, c0, nReal, root, printInfo)
    real(kind = wp),   intent(in)  :: c2, c1, c0
    integer,           intent(out) :: nReal
    real(kind = wp),   intent(out) :: root(1:3, 1:2)
    logical, optional, intent(in)  :: printInfo

    complex(kind = wp) :: poly(4), croots(3)

    poly(1) = cmplx(c0,     0.0_wp, kind = wp)
    poly(2) = cmplx(c1,     0.0_wp, kind = wp)
    poly(3) = cmplx(c2,     0.0_wp, kind = wp)
    poly(4) = cmplx(1.0_wp, 0.0_wp, kind = wp)

    call cmplx_roots_gen(croots, poly, 3, .true., .false.)
    call pack_roots(croots, 3, nReal, root)
  end subroutine cubicRoots

  subroutine quarticRoots(q3, q2, q1, q0, nReal, root, printInfo)
    real(kind = wp),   intent(in)  :: q3, q2, q1, q0
    integer,           intent(out) :: nReal
    real(kind = wp),   intent(out) :: root(1:4, 1:2)
    logical, optional, intent(in)  :: printInfo

    complex(kind = wp) :: poly(5), croots(4)

    poly(1) = cmplx(q0,     0.0_wp, kind = wp)
    poly(2) = cmplx(q1,     0.0_wp, kind = wp)
    poly(3) = cmplx(q2,     0.0_wp, kind = wp)
    poly(4) = cmplx(q3,     0.0_wp, kind = wp)
    poly(5) = cmplx(1.0_wp, 0.0_wp, kind = wp)

    call cmplx_roots_gen(croots, poly, 4, .true., .false.)
    call pack_roots(croots, 4, nReal, root)
  end subroutine quarticRoots

  subroutine pack_roots(croots, n, nReal, root)
    integer,            intent(in)  :: n
    complex(kind = wp), intent(in)  :: croots(n)
    integer,            intent(out) :: nReal
    real(kind = wp),    intent(out) :: root(n, 2)

    real(kind = wp) :: re_part(n), im_part(n), tol_i
    logical         :: is_real(n)
    integer         :: i, idx

    do i = 1, n
      re_part(i) = real(croots(i),  kind = wp)
      im_part(i) = aimag(croots(i))
      tol_i      = cmplx_tol * max(1.0_wp, abs(re_part(i)))
      is_real(i) = abs(im_part(i)) .le. tol_i
    end do

    nReal = count(is_real)
    root  = 0.0_wp

    idx = 0
    do i = 1, n
      if (is_real(i)) then
        idx = idx + 1
        root(idx, 1) = re_part(i)
        root(idx, 2) = 0.0_wp
      end if
    end do
    call sort_real_descending(root, nReal)

    do i = 1, n
      if (.not. is_real(i)) then
        idx = idx + 1
        root(idx, 1) = re_part(i)
        root(idx, 2) = im_part(i)
      end if
    end do
  end subroutine pack_roots

  subroutine sort_real_descending(arr, m)
    real(kind = wp), intent(inout) :: arr(:, :)
    integer,         intent(in)    :: m
    integer         :: i, j
    real(kind = wp) :: tmp_re, tmp_im

    do i = 2, m
      tmp_re = arr(i, 1)
      tmp_im = arr(i, 2)
      j = i - 1
      do while (j >= 1)
        if (arr(j, 1) >= tmp_re) exit
        arr(j+1, 1) = arr(j, 1)
        arr(j+1, 2) = arr(j, 2)
        j = j - 1
      end do
      arr(j+1, 1) = tmp_re
      arr(j+1, 2) = tmp_im
    end do
  end subroutine sort_real_descending

end module Polynomial234RootSolvers
