program test_analytic_circular_flux
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_analytic_circ_mod, only: set_field_analytic_circ, &
        write_flux_functions_analytic_circ, field_analytic_circ
    implicit none
    real(dp), parameter :: r0 = 165.0_dp, minor = 70.0_dp
    integer, parameter :: ntheta = 512, nradial = 256
    real(dp) :: row(5), field0, theta, rho, radius, z, br, bp, bz, deriv(9)
    real(dp) :: flux, winding, btheta, pi, dr, dtheta
    integer :: sign_case, iu, i, j, k
    character(len=100) :: header

    pi = acos(-1.0_dp)
    dtheta = 2.0_dp*pi/ntheta
    do sign_case = 1, 2
        field0 = (-1.0_dp)**sign_case*20000.0_dp
        call set_field_analytic_circ(r0, minor, field0, 1.0_dp, 2.0_dp, 2)
        call write_flux_functions_analytic_circ()
        open(newunit=iu, file='flux_functions.dat', status='old')
        read(iu, '(A)') header
        do i = 1, 1001
            read(iu, *) row
            if (mod(i-1, 200) /= 0) cycle
            if (i == 1) cycle
            dr = row(2)/nradial
            flux = 0.0_dp
            winding = 0.0_dp
            ! Independent surface-flux quadrature and field-line dphi/dtheta.
            do j = 1, ntheta
                theta = (real(j, dp)-0.5_dp)*dtheta
                radius = r0+row(2)*cos(theta)
                z = row(2)*sin(theta)
                call field_analytic_circ(radius, 0.0_dp, z, br, bp, bz, &
                    deriv(1), deriv(2), deriv(3), deriv(4), deriv(5), deriv(6), &
                    deriv(7), deriv(8), deriv(9))
                btheta = -br*sin(theta)+bz*cos(theta)
                winding = winding+row(2)*bp/(radius*btheta)/ntheta
                do k = 1, nradial
                    rho = (real(k, dp)-0.5_dp)*dr
                    radius = r0+rho*cos(theta)
                    z = rho*sin(theta)
                    call field_analytic_circ(radius, 0.0_dp, z, br, bp, bz, &
                        deriv(1), deriv(2), deriv(3), deriv(4), deriv(5), deriv(6), &
                        deriv(7), deriv(8), deriv(9))
                    flux = flux+bp*rho*dr/ntheta
                end do
            end do
            if (abs(row(3)/winding-1.0_dp) > 1.0e-12_dp) error stop 1
            if (abs(row(5)/flux-1.0_dp) > 1.0e-5_dp) error stop 2
        end do
        close(iu)
    end do
    print *, 'PASS: flux table agrees with magnetic surface flux and field winding'
end program test_analytic_circular_flux
