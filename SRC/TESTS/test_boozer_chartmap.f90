module test_boozer_chartmap
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use boozer_chartmap_mod, only: load_boozer_chartmap, evaluate_boozer_chartmap, &
                                   evaluate_chartmap_geometry, chartmap_iota, &
                                   chartmap_nfp, chartmap_rmajor, chartmap_zaxis
    use funit
    use netcdf

    implicit none
    real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp)
    character(len=*), parameter :: fixture = 'manufactured_boozer_chartmap.nc'

contains

    @before
    subroutine make_fixture()
        call write_manufactured_chartmap(fixture)
        call load_boozer_chartmap(fixture)
    end subroutine make_fixture

    @test
    subroutine test_field_and_jacobian_oracle()
        real(dp) :: as, atheta, aphi, hs, htheta, hphi, bmod, sqrtg
        real(dp) :: expected_btheta, expected_bphi, expected_sqrtg

        call evaluate_boozer_chartmap(0.25_dp, 0.37_dp, 0.42_dp, as, atheta, &
                                      aphi, hs, htheta, hphi, bmod, sqrtg)
        expected_btheta = 105.0_dp
        expected_bphi = 490.0_dp
        expected_sqrtg = (expected_bphi*2.0e8_dp &
                          - expected_btheta*(-0.8e8_dp))/1.0e6_dp

        @assertEqual(0.0_dp, as, tolerance=1.0e-13_dp)
        @assertEqual(0.5e8_dp, atheta, tolerance=1.0e-5_dp)
        @assertEqual(-0.2e8_dp, aphi, tolerance=1.0e-5_dp)
        @assertEqual(0.0_dp, hs, tolerance=1.0e-13_dp)
        @assertEqual(expected_btheta/1000.0_dp, htheta, tolerance=1.0e-13_dp)
        @assertEqual(expected_bphi/1000.0_dp, hphi, tolerance=1.0e-13_dp)
        @assertEqual(1000.0_dp, bmod, tolerance=1.0e-13_dp)
        @assertEqual(expected_sqrtg, sqrtg, tolerance=1.0e-7_dp)
        @assertEqual(0.4_dp, chartmap_iota(0.73_dp), tolerance=1.0e-13_dp)
        @assertEqual(2, chartmap_nfp)
        @assertEqual(200.0_dp, chartmap_rmajor, tolerance=1.0e-12_dp)
        @assertEqual(0.0_dp, chartmap_zaxis, tolerance=1.0e-12_dp)
    end subroutine test_field_and_jacobian_oracle

    @test
    subroutine test_periodicity_and_geometry_oracle()
        real(dp) :: a1(9), a2(9), r, phi, z, drds, dzds

        call evaluate_boozer_chartmap(0.25_dp, 0.25_dp, 0.31_dp, a1(1), a1(2), &
                                      a1(3), a1(4), a1(5), a1(6), a1(7), a1(8))
        call evaluate_boozer_chartmap(0.25_dp, 0.25_dp + 2.0_dp*pi, &
                                      0.31_dp + pi, a2(1), a2(2), a2(3), &
                                      a2(4), a2(5), a2(6), a2(7), a2(8))
        @assertEqual(a1(1:8), a2(1:8), tolerance=1.0e-10_dp)

        call evaluate_chartmap_geometry(0.25_dp, 0.0_dp, 0.0_dp, r, phi, z, &
                                        drds, dzds)
        @assertEqual(205.0_dp, r, tolerance=1.0e-12_dp)
        @assertEqual(0.0_dp, phi, tolerance=1.0e-12_dp)
        @assertEqual(0.0_dp, z, tolerance=1.0e-12_dp)
        @assertEqual(10.0_dp, drds, tolerance=1.0e-6_dp)
        @assertEqual(0.0_dp, dzds, tolerance=1.0e-6_dp)
    end subroutine test_periodicity_and_geometry_oracle

    subroutine write_manufactured_chartmap(filename)
        character(len=*), intent(in) :: filename
        integer :: ncid, dim_rho, dim_s, dim_theta, dim_zeta
        integer :: vrho, vs, vtheta, vzeta, vaphi, vbtheta, vbphi
        integer :: vx, vy, vz, vbmod, vnfp, ir, it, iz
        real(dp) :: rho(3), s(3), theta(4), zeta(2)
        real(dp) :: x(3, 4, 2), y(3, 4, 2), z(3, 4, 2), bmod(3, 4, 2)
        real(dp) :: aphi(3), btheta(3), bphi(3), radius

        rho = [0.0_dp, 0.5_dp, 1.0_dp]
        s = rho**2
        theta = [0.0_dp, 0.5_dp*pi, pi, 1.5_dp*pi]
        zeta = [0.0_dp, 0.5_dp*pi]
        aphi = -0.4_dp*2.0e8_dp*s
        btheta = 100.0_dp + 10.0_dp*rho
        bphi = 500.0_dp - 20.0_dp*rho
        bmod = 1000.0_dp
        do iz = 1, 2
            do it = 1, 4
                do ir = 1, 3
                    radius = 200.0_dp + 10.0_dp*rho(ir)*cos(theta(it))
                    x(ir, it, iz) = radius*cos(zeta(iz))
                    y(ir, it, iz) = radius*sin(zeta(iz))
                    z(ir, it, iz) = 10.0_dp*rho(ir)*sin(theta(it))
                end do
            end do
        end do

        call nc_ok(nf90_create(filename, nf90_clobber, ncid))
        call nc_ok(nf90_def_dim(ncid, 'rho', 3, dim_rho))
        call nc_ok(nf90_def_dim(ncid, 's', 3, dim_s))
        call nc_ok(nf90_def_dim(ncid, 'theta', 4, dim_theta))
        call nc_ok(nf90_def_dim(ncid, 'zeta', 2, dim_zeta))
        call nc_ok(nf90_def_var(ncid, 'rho', nf90_double, [dim_rho], vrho))
        call nc_ok(nf90_def_var(ncid, 's', nf90_double, [dim_s], vs))
        call nc_ok(nf90_def_var(ncid, 'theta', nf90_double, [dim_theta], vtheta))
        call nc_ok(nf90_def_var(ncid, 'zeta', nf90_double, [dim_zeta], vzeta))
        call nc_ok(nf90_def_var(ncid, 'A_phi', nf90_double, [dim_s], vaphi))
        call nc_ok(nf90_def_var(ncid, 'B_theta', nf90_double, [dim_rho], vbtheta))
        call nc_ok(nf90_def_var(ncid, 'B_phi', nf90_double, [dim_rho], vbphi))
        call nc_ok(nf90_def_var(ncid, 'x', nf90_double, &
                                [dim_rho, dim_theta, dim_zeta], vx))
        call nc_ok(nf90_def_var(ncid, 'y', nf90_double, &
                                [dim_rho, dim_theta, dim_zeta], vy))
        call nc_ok(nf90_def_var(ncid, 'z', nf90_double, &
                                [dim_rho, dim_theta, dim_zeta], vz))
        call nc_ok(nf90_def_var(ncid, 'Bmod', nf90_double, &
                                [dim_rho, dim_theta, dim_zeta], vbmod))
        call nc_ok(nf90_def_var(ncid, 'num_field_periods', nf90_int, vnfp))
        call nc_ok(nf90_put_att(ncid, nf90_global, 'torflux', 2.0e8_dp))
        call nc_ok(nf90_put_att(ncid, nf90_global, 'boozer_field', 1))
        call nc_ok(nf90_enddef(ncid))
        call nc_ok(nf90_put_var(ncid, vrho, rho))
        call nc_ok(nf90_put_var(ncid, vs, s))
        call nc_ok(nf90_put_var(ncid, vtheta, theta))
        call nc_ok(nf90_put_var(ncid, vzeta, zeta))
        call nc_ok(nf90_put_var(ncid, vaphi, aphi))
        call nc_ok(nf90_put_var(ncid, vbtheta, btheta))
        call nc_ok(nf90_put_var(ncid, vbphi, bphi))
        call nc_ok(nf90_put_var(ncid, vx, x))
        call nc_ok(nf90_put_var(ncid, vy, y))
        call nc_ok(nf90_put_var(ncid, vz, z))
        call nc_ok(nf90_put_var(ncid, vbmod, bmod))
        call nc_ok(nf90_put_var(ncid, vnfp, 2))
        call nc_ok(nf90_close(ncid))
    end subroutine write_manufactured_chartmap

    subroutine nc_ok(status)
        integer, intent(in) :: status
        if (status /= nf90_noerr) error stop nf90_strerror(status)
    end subroutine nc_ok

end module test_boozer_chartmap
