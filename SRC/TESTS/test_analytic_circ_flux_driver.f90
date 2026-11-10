! Fluxes of the analytic circular tokamak against quadratures of its field:
! psi_pol(rho) = int_R0^{R0+rho} R B_Z dR (all q branches, both signs) and the
! flux_functions.dat table (q = field-line winding, psi_tor = int B_phi dA/2pi).
program test_analytic_circ_flux_driver
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_analytic_circ_mod, only: set_field_analytic_circ, field_analytic_circ, &
        psi_pol_analytic_circ, write_flux_functions_analytic_circ
    implicit none
    real(dp), parameter :: r0 = 165.0_dp, minor = 70.0_dp
    integer, parameter :: n = 2048, nth = 512
    real(dp) :: pi, b0, q0, rho, acc, wind, th, row(5), b(3)
    integer :: ib, iq, p, grad, k, l, iu, i

    pi = acos(-1.0_dp)
    do ib = 1, 2
        b0 = (-1.0_dp)**ib*20000.0_dp
        do iq = 1, 2
            q0 = (-1.0_dp)**iq
            do p = 2, 4, 2
                do grad = 0, 1
                    call set_field_analytic_circ(r0, minor, b0, q0, 2.0_dp*grad*q0, p)
                    do i = 1, 4
                        rho = minor*i/4.0_dp
                        acc = 0.0_dp
                        do k = 1, n
                            b = field_at(r0 + k*rho/n, 0.0_dp)
                            acc = acc + merge(1, merge(4, 2, mod(k, 2) == 1), k == n) &
                                *(r0 + k*rho/n)*b(3)*rho/n/3.0_dp
                        end do
                        if (abs(psi_pol_analytic_circ(rho)/acc - 1.0_dp) > 1.0e-10_dp) &
                            error stop 'psi_pol differs from the field quadrature'
                    end do
                end do
            end do
        end do

        call set_field_analytic_circ(r0, minor, b0, 1.0_dp, 2.0_dp, 2)
        call write_flux_functions_analytic_circ()
        open(newunit=iu, file='flux_functions.dat', status='old')
        read(iu, *)
        do i = 1, 1001
            read(iu, *) row
            if (i == 1 .or. mod(i - 1, 200) /= 0) cycle
            wind = 0.0_dp
            acc = 0.0_dp
            do l = 1, nth
                th = (l - 0.5_dp)*2.0_dp*pi/nth
                b = field_at(r0 + row(2)*cos(th), row(2)*sin(th))
                wind = wind + row(2)*b(2)/((r0 + row(2)*cos(th)) &
                    *(-b(1)*sin(th) + b(3)*cos(th)))/nth
                do k = 1, 256
                    rho = (k - 0.5_dp)*row(2)/256
                    b = field_at(r0 + rho*cos(th), rho*sin(th))
                    acc = acc + b(2)*rho*row(2)/256/nth
                end do
            end do
            if (abs(row(3)/wind - 1.0_dp) > 1.0e-12_dp) error stop 'q is not the winding'
            if (abs(row(5)/acc - 1.0_dp) > 1.0e-5_dp) error stop 'psi_tor differs'
        end do
        close(iu)
    end do
    print *, 'PASS: analytic circular fluxes agree with field quadratures'

contains

    function field_at(r, z) result(b)
        real(dp), intent(in) :: r, z
        real(dp) :: b(3), d(9)
        call field_analytic_circ(r, 0.0_dp, z, b(1), b(2), b(3), d(1), d(2), d(3), &
            d(4), d(5), d(6), d(7), d(8), d(9))
    end function field_at

end program test_analytic_circ_flux_driver
