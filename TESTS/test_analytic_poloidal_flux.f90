program test_analytic_poloidal_flux
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_analytic_circ_mod, only: set_field_analytic_circ, &
        psi_pol_analytic_circ, field_analytic_circ
    implicit none
    integer, parameter :: steps = 2048
    real(dp), parameter :: r0 = 165.0_dp, minor = 70.0_dp
    real(dp) :: b0, q0, q1, rho, dr, r, br, bp, bz, derivatives(9)
    real(dp) :: measured, expected, weight, flux_derivative, delta
    integer :: ib, iq, exponent, gradient, radial, k

    do ib = 1, 2
        b0 = (-1.0_dp)**ib*20000.0_dp
        do iq = 1, 2
            q0 = (-1.0_dp)**iq
            do exponent = 2, 4, 2
                do gradient = 0, 1
                    q1 = 2.0_dp*gradient*q0
                    call set_field_analytic_circ(r0, minor, b0, q0, q1, exponent)
                    do radial = 1, 4
                        rho = minor*real(radial, dp)/4.0_dp
                        dr = rho/steps
                        expected = 0.0_dp
                        ! Independent physical poloidal-flux oracle: integrate
                        ! R*B_Z at theta=0, using the actual magnetic field.
                        do k = 0, steps
                            r = r0+k*dr
                            call field_analytic_circ(r, 0.0_dp, 0.0_dp, &
                                br, bp, bz, derivatives(1), derivatives(2), &
                                derivatives(3), derivatives(4), derivatives(5), &
                                derivatives(6), derivatives(7), derivatives(8), &
                                derivatives(9))
                            weight = 2.0_dp
                            if (mod(k, 2) /= 0) weight = 4.0_dp
                            if (k == 0 .or. k == steps) weight = 1.0_dp
                            expected = expected + weight*r*bz*dr/3.0_dp
                        end do
                        measured = psi_pol_analytic_circ(rho)
                        if (abs(measured/expected - 1.0_dp) > 1.e-10_dp) &
                            error stop 'flux integral does not match actual field'
                        delta = 1.e-4_dp*rho
                        flux_derivative = (psi_pol_analytic_circ(rho + delta) &
                            - psi_pol_analytic_circ(rho - delta))/(2.0_dp*delta)
                        call field_analytic_circ(r0 + rho, 0.0_dp, 0.0_dp, &
                            br, bp, bz, derivatives(1), derivatives(2), &
                            derivatives(3), derivatives(4), derivatives(5), &
                            derivatives(6), derivatives(7), derivatives(8), &
                            derivatives(9))
                        if (abs(flux_derivative/((r0 + rho)*bz) - 1.0_dp) &
                            > 1.e-7_dp) error stop 'flux derivative mismatch'
                    end do
                end do
            end do
        end do
    end do
    print *, 'PASS: physical poloidal flux, constant q, p=2/4, both signs'
end program test_analytic_poloidal_flux
