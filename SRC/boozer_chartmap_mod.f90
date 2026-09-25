module boozer_chartmap_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use netcdf

    implicit none
    private

    real(dp), parameter :: twopi = 8.0_dp*atan(1.0_dp)

    integer, public, protected :: chartmap_nfp = 1
    real(dp), public, protected :: chartmap_rmajor = 0.0_dp
    real(dp), public, protected :: chartmap_zaxis = 0.0_dp
    real(dp), public, protected :: chartmap_torflux = 0.0_dp

    real(dp), allocatable :: rho_grid(:), s_grid(:), aphi_grid(:)
    real(dp), allocatable :: btheta_grid(:), bphi_grid(:)
    real(dp), allocatable :: x_grid(:, :, :), y_grid(:, :, :), z_grid(:, :, :)
    real(dp), allocatable :: bmod_grid(:, :, :)
    integer :: nrho = 0, ns = 0, ntheta = 0, nzeta = 0
    real(dp) :: theta_period = twopi, zeta_period = twopi
    logical :: theta_has_endpoint = .false., zeta_has_endpoint = .false.
    logical :: initialized = .false.

    public :: load_boozer_chartmap, evaluate_boozer_chartmap
    public :: evaluate_chartmap_geometry, chartmap_iota

contains

    subroutine load_boozer_chartmap(filename)
        character(len=*), intent(in) :: filename

        integer :: ncid, status
        integer :: boozer_field
        real(dp), allocatable :: theta(:), zeta(:)

        call clear_chartmap()
        call check_nc(nf90_open(trim(filename), nf90_nowrite, ncid), 'open chartmap')

        call read_dimension(ncid, 'rho', nrho)
        call read_dimension(ncid, 's', ns)
        call read_dimension(ncid, 'theta', ntheta)
        call read_dimension(ncid, 'zeta', nzeta)

        if (nrho < 2 .or. ns < 2 .or. ntheta < 2 .or. nzeta < 1) then
            error stop 'Boozer chartmap dimensions are too small'
        end if
        allocate(rho_grid(nrho), s_grid(ns), aphi_grid(ns))
        allocate(btheta_grid(nrho), bphi_grid(nrho))
        allocate(theta(ntheta), zeta(nzeta))
        allocate(x_grid(nrho, ntheta, nzeta), y_grid(nrho, ntheta, nzeta))
        allocate(z_grid(nrho, ntheta, nzeta), bmod_grid(nrho, ntheta, nzeta))

        call read_real_1d(ncid, 'rho', rho_grid)
        call read_real_1d(ncid, 's', s_grid)
        call read_real_1d(ncid, 'theta', theta)
        call read_real_1d(ncid, 'zeta', zeta)
        call read_real_1d(ncid, 'A_phi', aphi_grid)
        call read_real_1d(ncid, 'B_theta', btheta_grid)
        call read_real_1d(ncid, 'B_phi', bphi_grid)
        call read_real_3d(ncid, 'x', x_grid)
        call read_real_3d(ncid, 'y', y_grid)
        call read_real_3d(ncid, 'z', z_grid)
        call read_real_3d(ncid, 'Bmod', bmod_grid)
        call read_int_scalar(ncid, 'num_field_periods', chartmap_nfp)
        call check_nc(nf90_get_att(ncid, nf90_global, 'torflux', chartmap_torflux), &
                      'read torflux attribute')

        boozer_field = 0
        status = nf90_get_att(ncid, nf90_global, 'boozer_field', boozer_field)
        if (status /= nf90_noerr .or. boozer_field /= 1) then
            error stop 'NetCDF file is not an extended Boozer chartmap'
        end if
        call check_nc(nf90_close(ncid), 'close chartmap')

        if (chartmap_nfp < 1) error stop 'Invalid chartmap field-period count'
        theta_period = twopi
        zeta_period = twopi/real(chartmap_nfp, dp)
        call validate_axis(rho_grid, 'rho')
        call validate_axis(s_grid, 's')
        call validate_periodic_axis(theta, theta_period, theta_has_endpoint, 'theta')
        call validate_periodic_axis(zeta, zeta_period, zeta_has_endpoint, 'zeta')
        if (minval(bmod_grid) <= 0.0_dp) error stop 'Chartmap Bmod must be positive'
        if (abs(chartmap_torflux) <= tiny(1.0_dp)) then
            error stop 'Chartmap torflux must be nonzero'
        end if

        chartmap_rmajor = sum(sqrt(x_grid(1, :, :)**2 + y_grid(1, :, :)**2)) &
                          /real(ntheta*nzeta, dp)
        chartmap_zaxis = sum(z_grid(1, :, :))/real(ntheta*nzeta, dp)
        initialized = .true.
        deallocate(theta, zeta)
    end subroutine load_boozer_chartmap

    subroutine evaluate_boozer_chartmap(s, theta, zeta, a_s, a_theta, a_phi, &
                                        h_s, h_theta, h_phi, bmod, sqrtg)
        real(dp), intent(in) :: s, theta, zeta
        real(dp), intent(out) :: a_s, a_theta, a_phi
        real(dp), intent(out) :: h_s, h_theta, h_phi, bmod, sqrtg

        real(dp) :: rho, daphids, btheta, bphi

        call assert_loaded()
        rho = sqrt(max(0.0_dp, s))
        call interpolate_profile(s_grid, aphi_grid, s, a_phi, daphids)
        call interpolate_profile(rho_grid, btheta_grid, rho, btheta)
        call interpolate_profile(rho_grid, bphi_grid, rho, bphi)
        call interpolate_periodic_3d(bmod_grid, rho, theta, zeta, bmod, .false.)

        a_s = 0.0_dp
        a_theta = chartmap_torflux*s
        h_s = 0.0_dp
        h_theta = btheta/bmod
        h_phi = bphi/bmod
        sqrtg = (bphi*chartmap_torflux - btheta*daphids)/(bmod*bmod)
    end subroutine evaluate_boozer_chartmap

    subroutine evaluate_chartmap_geometry(s, theta, zeta, r_cyl, phi_cyl, z_cyl, &
                                          dr_ds, dz_ds)
        real(dp), intent(in) :: s, theta, zeta
        real(dp), intent(out) :: r_cyl, phi_cyl, z_cyl, dr_ds, dz_ds

        real(dp) :: rho, x, y, ds_probe, r_plus, r_minus, z_plus, z_minus
        real(dp) :: x_dummy, y_dummy

        call assert_loaded()
        rho = sqrt(max(0.0_dp, s))
        call interpolate_periodic_3d(x_grid, rho, theta, zeta, x, .true., y_grid, y)
        call interpolate_periodic_3d(z_grid, rho, theta, zeta, z_cyl, .false.)
        r_cyl = sqrt(x*x + y*y)
        phi_cyl = atan2(y, x)
        if (phi_cyl < 0.0_dp) phi_cyl = phi_cyl + twopi

        ds_probe = max(1.0e-8_dp, 1.0e-6_dp*minval(s_grid(2:) - s_grid(:ns - 1)))
        call geometry_rz(min(s_grid(ns), s + ds_probe), theta, zeta, r_plus, z_plus)
        call geometry_rz(max(s_grid(1), s - ds_probe), theta, zeta, r_minus, z_minus)
        dr_ds = (r_plus - r_minus) &
                / (min(s_grid(ns), s + ds_probe) - max(s_grid(1), s - ds_probe))
        dz_ds = (z_plus - z_minus) &
                / (min(s_grid(ns), s + ds_probe) - max(s_grid(1), s - ds_probe))
    end subroutine evaluate_chartmap_geometry

    function chartmap_iota(s) result(iota)
        real(dp), intent(in) :: s
        real(dp) :: iota, aphi, daphids

        call assert_loaded()
        call interpolate_profile(s_grid, aphi_grid, s, aphi, daphids)
        iota = -daphids/chartmap_torflux
    end function chartmap_iota

    subroutine geometry_rz(s, theta, zeta, r_cyl, z_cyl)
        real(dp), intent(in) :: s, theta, zeta
        real(dp), intent(out) :: r_cyl, z_cyl
        real(dp) :: x, y, rho

        rho = sqrt(max(0.0_dp, s))
        call interpolate_periodic_3d(x_grid, rho, theta, zeta, x, .true., y_grid, y)
        call interpolate_periodic_3d(z_grid, rho, theta, zeta, z_cyl, .false.)
        r_cyl = sqrt(x*x + y*y)
    end subroutine geometry_rz

    subroutine interpolate_profile(axis, values, coordinate, value, derivative)
        real(dp), intent(in) :: axis(:), values(:), coordinate
        real(dp), intent(out) :: value
        real(dp), intent(out), optional :: derivative
        integer :: i
        real(dp) :: u, x_clamped, slope

        x_clamped = min(axis(size(axis)), max(axis(1), coordinate))
        i = locate_interval(axis, x_clamped)
        slope = (values(i + 1) - values(i))/(axis(i + 1) - axis(i))
        u = (x_clamped - axis(i))/(axis(i + 1) - axis(i))
        value = (1.0_dp - u)*values(i) + u*values(i + 1)
        if (present(derivative)) derivative = slope
    end subroutine interpolate_profile

    subroutine interpolate_periodic_3d(values, rho, theta, zeta, value, rotate_xy, &
                                       y_values, y_value)
        real(dp), intent(in) :: values(:, :, :), rho, theta, zeta
        real(dp), intent(out) :: value
        logical, intent(in) :: rotate_xy
        real(dp), intent(in), optional :: y_values(:, :, :)
        real(dp), intent(out), optional :: y_value

        integer :: ir, it0, it1, iz0, iz1, dr, dt, dz
        real(dp) :: ur, ut, uz, weight, sample_x, sample_y, angle
        real(dp) :: theta_wrapped, zeta_wrapped
        logical :: zeta_wrap

        ir = locate_interval(rho_grid, min(rho_grid(nrho), max(rho_grid(1), rho)))
        ur = (min(rho_grid(nrho), max(rho_grid(1), rho)) - rho_grid(ir)) &
             /(rho_grid(ir + 1) - rho_grid(ir))
        theta_wrapped = modulo(theta, theta_period)
        zeta_wrapped = modulo(zeta, zeta_period)
        call periodic_indices(theta_wrapped, theta_period, ntheta, theta_has_endpoint, &
                              it0, it1, ut, zeta_wrap)
        call periodic_indices(zeta_wrapped, zeta_period, nzeta, zeta_has_endpoint, &
                              iz0, iz1, uz, zeta_wrap)

        value = 0.0_dp
        if (present(y_value)) y_value = 0.0_dp
        do dz = 0, 1
            do dt = 0, 1
                do dr = 0, 1
                    weight = merge(ur, 1.0_dp - ur, dr == 1) &
                             *merge(ut, 1.0_dp - ut, dt == 1) &
                             *merge(uz, 1.0_dp - uz, dz == 1)
                    sample_x = values(ir + dr, merge(it1, it0, dt == 1), &
                                      merge(iz1, iz0, dz == 1))
                    if (rotate_xy) then
                        if (.not. present(y_values) .or. .not. present(y_value)) then
                            error stop 'Both Cartesian chartmap components are required'
                        end if
                        sample_y = y_values(ir + dr, merge(it1, it0, dt == 1), &
                                            merge(iz1, iz0, dz == 1))
                        if (dz == 1 .and. iz1 == 1 .and. .not. zeta_has_endpoint) then
                            angle = zeta_period
                            call rotate_cartesian(sample_x, sample_y, angle)
                        end if
                        y_value = y_value + weight*sample_y
                    end if
                    value = value + weight*sample_x
                end do
            end do
        end do
    end subroutine interpolate_periodic_3d

    subroutine periodic_indices(coordinate, period, count, has_endpoint, i0, i1, u, wrapped)
        real(dp), intent(in) :: coordinate, period
        integer, intent(in) :: count
        logical, intent(in) :: has_endpoint
        integer, intent(out) :: i0, i1
        real(dp), intent(out) :: u
        logical, intent(out) :: wrapped
        integer :: cells
        real(dp) :: scaled

        cells = merge(count - 1, count, has_endpoint)
        scaled = coordinate*real(cells, dp)/period
        i0 = min(cells - 1, int(scaled)) + 1
        u = scaled - real(i0 - 1, dp)
        i1 = i0 + 1
        wrapped = .false.
        if (i1 > count) then
            i1 = 1
            wrapped = .true.
        end if
    end subroutine periodic_indices

    integer function locate_interval(axis, coordinate) result(index)
        real(dp), intent(in) :: axis(:), coordinate
        integer :: lo, hi, mid

        lo = 1
        hi = size(axis)
        do while (hi - lo > 1)
            mid = (lo + hi)/2
            if (coordinate >= axis(mid)) then
                lo = mid
            else
                hi = mid
            end if
        end do
        index = min(size(axis) - 1, lo)
    end function locate_interval

    subroutine rotate_cartesian(x, y, angle)
        real(dp), intent(inout) :: x, y
        real(dp), intent(in) :: angle
        real(dp) :: old_x

        old_x = x
        x = cos(angle)*old_x - sin(angle)*y
        y = sin(angle)*old_x + cos(angle)*y
    end subroutine rotate_cartesian

    subroutine validate_axis(axis, name)
        real(dp), intent(in) :: axis(:)
        character(len=*), intent(in) :: name

        if (any(axis(2:) <= axis(:size(axis) - 1))) then
            error stop 'Chartmap axis is not strictly increasing: '//name
        end if
    end subroutine validate_axis

    subroutine validate_periodic_axis(axis, period, has_endpoint, name)
        real(dp), intent(in) :: axis(:), period
        logical, intent(out) :: has_endpoint
        character(len=*), intent(in) :: name
        real(dp) :: spacing, tolerance

        call validate_axis(axis, name)
        if (abs(axis(1)) > 1.0e-10_dp) then
            error stop 'Periodic chartmap axis must start at zero: '//name
        end if
        spacing = axis(2) - axis(1)
        tolerance = max(1.0e-10_dp, 1.0e-6_dp*spacing)
        has_endpoint = abs(axis(size(axis)) - period) <= tolerance
        if (.not. has_endpoint) then
            if (abs(spacing*real(size(axis), dp) - period) > tolerance) then
                error stop 'Periodic chartmap axis must be uniform: '//name
            end if
        else
            if (abs(spacing*real(size(axis) - 1, dp) - period) > tolerance) then
                error stop 'Periodic chartmap axis must be uniform: '//name
            end if
        end if
    end subroutine validate_periodic_axis

    subroutine read_dimension(ncid, name, length)
        integer, intent(in) :: ncid
        character(len=*), intent(in) :: name
        integer, intent(out) :: length
        integer :: dimid

        call check_nc(nf90_inq_dimid(ncid, name, dimid), 'find dimension '//name)
        call check_nc(nf90_inquire_dimension(ncid, dimid, len=length), &
                      'read dimension '//name)
    end subroutine read_dimension

    subroutine read_real_1d(ncid, name, values)
        integer, intent(in) :: ncid
        character(len=*), intent(in) :: name
        real(dp), intent(out) :: values(:)
        integer :: varid

        call check_nc(nf90_inq_varid(ncid, name, varid), 'find variable '//name)
        call check_nc(nf90_get_var(ncid, varid, values), 'read variable '//name)
    end subroutine read_real_1d

    subroutine read_real_3d(ncid, name, values)
        integer, intent(in) :: ncid
        character(len=*), intent(in) :: name
        real(dp), intent(out) :: values(:, :, :)
        integer :: varid

        call check_nc(nf90_inq_varid(ncid, name, varid), 'find variable '//name)
        call check_nc(nf90_get_var(ncid, varid, values), 'read variable '//name)
    end subroutine read_real_3d

    subroutine read_int_scalar(ncid, name, value)
        integer, intent(in) :: ncid
        character(len=*), intent(in) :: name
        integer, intent(out) :: value
        integer :: varid

        call check_nc(nf90_inq_varid(ncid, name, varid), 'find variable '//name)
        call check_nc(nf90_get_var(ncid, varid, value), 'read variable '//name)
    end subroutine read_int_scalar

    subroutine check_nc(status, operation)
        integer, intent(in) :: status
        character(len=*), intent(in) :: operation

        if (status /= nf90_noerr) then
            write(*, '(A)') trim(operation)//': '//trim(nf90_strerror(status))
            error stop 'Boozer chartmap NetCDF error'
        end if
    end subroutine check_nc

    subroutine assert_loaded()
        if (.not. initialized) error stop 'Boozer chartmap has not been loaded'
    end subroutine assert_loaded

    subroutine clear_chartmap()
        if (allocated(rho_grid)) deallocate(rho_grid, s_grid, aphi_grid)
        if (allocated(btheta_grid)) deallocate(btheta_grid, bphi_grid)
        if (allocated(x_grid)) deallocate(x_grid, y_grid, z_grid, bmod_grid)
        initialized = .false.
    end subroutine clear_chartmap

end module boozer_chartmap_mod
