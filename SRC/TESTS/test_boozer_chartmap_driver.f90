program test_boozer_chartmap_driver
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use boozer_chartmap_mod, only: load_boozer_chartmap, evaluate_boozer_chartmap, &
                                   evaluate_chartmap_geometry, chartmap_iota, &
                                   chartmap_nfp, chartmap_rmajor, chartmap_zaxis
    use netcdf

    implicit none

    real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp)
    real(dp) :: first(8), periodic(8), r, phi, z, drds, dzds
    real(dp) :: btheta, bphi, expected_sqrtg
    integer :: failures
    character(len=1024) :: external_file

    failures = 0
    if (command_argument_count() > 0) then
        call get_command_argument(1, external_file)
        call load_boozer_chartmap(trim(external_file))
        call evaluate_boozer_chartmap(0.25_dp, 0.37_dp, 0.42_dp, first(1), &
                                      first(2), first(3), first(4), first(5), &
                                      first(6), first(7), first(8))
        if (.not. all(ieee_is_finite(first))) error stop 'Non-finite chartmap field'
        if (first(7) <= 0.0_dp) error stop 'Non-positive chartmap Bmod'
        if (.not. ieee_is_finite(chartmap_iota(0.25_dp))) then
            error stop 'Non-finite chartmap iota'
        end if
        stop
    end if
    call write_fixture('manufactured_boozer_chartmap.nc')
    call load_boozer_chartmap('manufactured_boozer_chartmap.nc')
    call evaluate_boozer_chartmap(0.25_dp, 0.37_dp, 0.42_dp, first(1), first(2), &
                                  first(3), first(4), first(5), first(6), &
                                  first(7), first(8))
    btheta = 105.0_dp
    bphi = 490.0_dp
    expected_sqrtg = (bphi*2.0e8_dp - btheta*(-0.8e8_dp))/1.0e6_dp
    call check_close('A_s', first(1), 0.0_dp, 1.0e-13_dp, failures)
    call check_close('A_theta', first(2), 0.5e8_dp, 1.0e-5_dp, failures)
    call check_close('A_phi', first(3), -0.2e8_dp, 1.0e-5_dp, failures)
    call check_close('h_s', first(4), 0.0_dp, 1.0e-13_dp, failures)
    call check_close('h_theta', first(5), btheta/1000.0_dp, 1.0e-13_dp, failures)
    call check_close('h_phi', first(6), bphi/1000.0_dp, 1.0e-13_dp, failures)
    call check_close('Bmod', first(7), 1000.0_dp, 1.0e-13_dp, failures)
    call check_close('sqrtg', first(8), expected_sqrtg, 1.0e-7_dp, failures)
    call check_close('iota', chartmap_iota(0.73_dp), 0.4_dp, 1.0e-13_dp, failures)
    if (chartmap_nfp /= 2) failures = failures + 1
    call check_close('axis R', chartmap_rmajor, 200.0_dp, 1.0e-12_dp, failures)
    call check_close('axis Z', chartmap_zaxis, 0.0_dp, 1.0e-12_dp, failures)

    call evaluate_boozer_chartmap(0.25_dp, 0.25_dp, 0.31_dp, first(1), first(2), &
                                  first(3), first(4), first(5), first(6), &
                                  first(7), first(8))
    call evaluate_boozer_chartmap(0.25_dp, 0.25_dp + 2.0_dp*pi, 0.31_dp + pi, &
                                  periodic(1), periodic(2), periodic(3), &
                                  periodic(4), periodic(5), periodic(6), &
                                  periodic(7), periodic(8))
    call check_close('periodicity', maxval(abs(first - periodic)), 0.0_dp, &
                     1.0e-10_dp, failures)
    call evaluate_chartmap_geometry(0.25_dp, 0.0_dp, 0.0_dp, r, phi, z, drds, dzds)
    call check_close('R', r, 205.0_dp, 1.0e-12_dp, failures)
    call check_close('phi', phi, 0.0_dp, 1.0e-12_dp, failures)
    call check_close('Z', z, 0.0_dp, 1.0e-12_dp, failures)
    call check_close('dR/ds', drds, 10.0_dp, 1.0e-6_dp, failures)
    call check_close('dZ/ds', dzds, 0.0_dp, 1.0e-6_dp, failures)

    if (failures /= 0) error stop 'Boozer chartmap oracle failed'

contains

    subroutine check_close(name, actual, expected, tolerance, count)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: actual, expected, tolerance
        integer, intent(inout) :: count

        if (abs(actual - expected) > tolerance) then
            write(*, '(A,2ES24.15)') trim(name)//' mismatch: ', actual, expected
            count = count + 1
        end if
    end subroutine check_close

    subroutine write_fixture(filename)
        character(len=*), intent(in) :: filename
        integer :: ncid, dims(4), vars(12), ir, it, iz
        integer :: dim1(1), dim3(3)
        real(dp) :: rho(3), s(3), theta(4), zeta(2), radius
        real(dp) :: x(3, 4, 2), y(3, 4, 2), zz(3, 4, 2), bmod(3, 4, 2)
        real(dp) :: aphi(3), btheta_profile(3), bphi_profile(3)

        rho = [0.0_dp, 0.5_dp, 1.0_dp]
        s = rho**2
        theta = [0.0_dp, 0.5_dp*pi, pi, 1.5_dp*pi]
        zeta = [0.0_dp, 0.5_dp*pi]
        bmod = 1000.0_dp
        aphi = -0.4_dp*2.0e8_dp*s
        btheta_profile = 100.0_dp + 10.0_dp*rho
        bphi_profile = 500.0_dp - 20.0_dp*rho
        do iz = 1, 2
            do it = 1, 4
                do ir = 1, 3
                    radius = 200.0_dp + 10.0_dp*rho(ir)*cos(theta(it))
                    x(ir, it, iz) = radius*cos(zeta(iz))
                    y(ir, it, iz) = radius*sin(zeta(iz))
                    zz(ir, it, iz) = 10.0_dp*rho(ir)*sin(theta(it))
                end do
            end do
        end do

        call ok(nf90_create(filename, nf90_clobber, ncid))
        call ok(nf90_def_dim(ncid, 'rho', 3, dims(1)))
        call ok(nf90_def_dim(ncid, 's', 3, dims(2)))
        call ok(nf90_def_dim(ncid, 'theta', 4, dims(3)))
        call ok(nf90_def_dim(ncid, 'zeta', 2, dims(4)))
        dim1 = dims(1)
        call ok(nf90_def_var(ncid, 'rho', nf90_double, dim1, vars(1)))
        dim1 = dims(2)
        call ok(nf90_def_var(ncid, 's', nf90_double, dim1, vars(2)))
        dim1 = dims(3)
        call ok(nf90_def_var(ncid, 'theta', nf90_double, dim1, vars(3)))
        dim1 = dims(4)
        call ok(nf90_def_var(ncid, 'zeta', nf90_double, dim1, vars(4)))
        dim1 = dims(2)
        call ok(nf90_def_var(ncid, 'A_phi', nf90_double, dim1, vars(5)))
        dim1 = dims(1)
        call ok(nf90_def_var(ncid, 'B_theta', nf90_double, dim1, vars(6)))
        call ok(nf90_def_var(ncid, 'B_phi', nf90_double, dim1, vars(7)))
        dim3 = dims([1, 3, 4])
        call ok(nf90_def_var(ncid, 'x', nf90_double, dim3, vars(8)))
        call ok(nf90_def_var(ncid, 'y', nf90_double, dim3, vars(9)))
        call ok(nf90_def_var(ncid, 'z', nf90_double, dim3, vars(10)))
        call ok(nf90_def_var(ncid, 'Bmod', nf90_double, dim3, vars(11)))
        call ok(nf90_def_var(ncid, 'num_field_periods', nf90_int, vars(12)))
        call ok(nf90_put_att(ncid, nf90_global, 'torflux', 2.0e8_dp))
        call ok(nf90_put_att(ncid, nf90_global, 'boozer_field', 1))
        call ok(nf90_enddef(ncid))
        call ok(nf90_put_var(ncid, vars(1), rho))
        call ok(nf90_put_var(ncid, vars(2), s))
        call ok(nf90_put_var(ncid, vars(3), theta))
        call ok(nf90_put_var(ncid, vars(4), zeta))
        call ok(nf90_put_var(ncid, vars(5), aphi))
        call ok(nf90_put_var(ncid, vars(6), btheta_profile))
        call ok(nf90_put_var(ncid, vars(7), bphi_profile))
        call ok(nf90_put_var(ncid, vars(8), x))
        call ok(nf90_put_var(ncid, vars(9), y))
        call ok(nf90_put_var(ncid, vars(10), zz))
        call ok(nf90_put_var(ncid, vars(11), bmod))
        call ok(nf90_put_var(ncid, vars(12), 2))
        call ok(nf90_close(ncid))
    end subroutine write_fixture

    subroutine ok(status)
        integer, intent(in) :: status
        if (status /= nf90_noerr) error stop nf90_strerror(status)
    end subroutine ok

end program test_boozer_chartmap_driver
