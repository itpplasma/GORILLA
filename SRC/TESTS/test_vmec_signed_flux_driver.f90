program test_vmec_signed_flux
    use netcdf
    use constants, only: pi
    use tetra_grid_settings_mod, only: load_tetra_grid_inp
    use spline_vmec_data_mod, only: spline_vmec_data
    use splint_vmec_data_mod, only: vmec_field
    implicit none

    real(8), parameter :: radius_m = 2.d0, minor_m = .5d0, iota = .4d0
    real(8), parameter :: s = .25d0, theta = 0.d0, phi = 0.d0
    character(32) :: arg, mode
    integer :: orientation, field_sign, geom_sign, unit
    real(8) :: expected, b_tor_g, b_pol_g, bmod_g, r_cm
    real(8) :: a_theta, a_phi, dat, dap, iota_out, jac, alam, ls, lt, lp
    real(8) :: bt, bp, bc_s, bc_t, bc_p

    call get_command_argument(1, arg)
    read(arg,*) orientation
    call get_command_argument(2, arg)
    read(arg,*) field_sign
    call get_command_argument(3, mode)
    geom_sign = merge(-1, orientation, orientation == 0)
    call write_wout(geom_sign, field_sign, orientation /= 0)
    open(newunit=unit, file='tetra_grid.inp', status='replace')
    write(unit,'(a)') '&TETRA_GRID_NML', ' grid_kind=3, n1=8, n2=8, n3=8,', &
        " netcdf_filename='signed_flux.wout.nc', sfc_s_min=.01, theta_geom_flux=1,"
    if (mode == 'legacy') write(unit,'(a)') ' boole_vmec_legacy_flux=.true.,'
    write(unit,'(a)') '/'
    close(unit)
    call load_tetra_grid_inp()
    call spline_vmec_data
    call vmec_field(s,theta,phi,a_theta,a_phi,dat,dap,iota_out,jac,alam,ls,lt,lp, &
                    bt,bp,bc_s,bc_t,bc_p)

    r_cm = 100.d0*(radius_m + minor_m*sqrt(s))
    b_tor_g = r_cm*bp
    bmod_g = sqrt(bt*bc_t + bp*bc_p)
    expected = 2.d4*field_sign
    if (orientation == 0 .or. mode == 'legacy') expected = expected*geom_sign
    b_pol_g = expected*iota*100.d0*minor_m*sqrt(s)/r_cm
    print '(a,f12.3,a,f12.3,a)', 'VMEC B_phi = ', b_tor_g, ' G, expected ', expected, ' G'
    if (sign(1.d0,jac) /= geom_sign) error stop 'VMEC Jacobian handedness changed'
    if (abs(iota_out + geom_sign*iota) > 1.d-12) error stop 'VMEC iota changed'
    if (abs(b_tor_g/expected - 1.d0) > 2.d-4) error stop 'Wrong physical toroidal field'
    if (abs(bmod_g/hypot(expected,b_pol_g) - 1.d0) > 2.d-4) error stop 'Wrong field magnitude'

contains

    subroutine write_wout(handedness, polarity, write_signgs)
        integer, intent(in) :: handedness, polarity
        logical, intent(in) :: write_signgs
        integer, parameter :: ns = 17, nm = 3
        real(8) :: radial(ns), rc(nm,ns), zs(nm,ns), lambda(nm,ns)
        integer :: i, file, m, r, nfp_id, major_id, sign_id
        integer :: phi_id, iota_id, xm_id, xn_id, r_id, z_id, l_id

        radial = [(real(i-1,8)/real(ns-1,8), i=1,ns)]
        rc(1,:) = radius_m
        rc(2,:) = minor_m*sqrt(radial)
        rc(3,:) = 0.d0
        zs = 0.d0
        zs(2,:) = -handedness*minor_m*sqrt(radial)
        lambda = 0.d0
        call nf90_check(nf90_create('signed_flux.wout.nc', nf90_clobber, file))
        call nf90_check(nf90_def_dim(file, 'mn_mode', nm, m))
        call nf90_check(nf90_def_dim(file, 'radius', ns, r))
        call nf90_check(nf90_def_var(file, 'nfp', nf90_int, nfp_id))
        call nf90_check(nf90_def_var(file, 'Rmajor_p', nf90_double, major_id))
        if (write_signgs) call nf90_check(nf90_def_var(file, 'signgs', nf90_int, sign_id))
        call nf90_check(nf90_def_var(file, 'phi', nf90_double, [r], phi_id))
        call nf90_check(nf90_def_var(file, 'iotaf', nf90_double, [r], iota_id))
        call nf90_check(nf90_def_var(file, 'xm', nf90_double, [m], xm_id))
        call nf90_check(nf90_def_var(file, 'xn', nf90_double, [m], xn_id))
        call nf90_check(nf90_def_var(file, 'rmnc', nf90_double, [m,r], r_id))
        call nf90_check(nf90_def_var(file, 'zmns', nf90_double, [m,r], z_id))
        call nf90_check(nf90_def_var(file, 'lmns', nf90_double, [m,r], l_id))
        call nf90_check(nf90_enddef(file))
        call nf90_check(nf90_put_var(file, nfp_id, 1))
        call nf90_check(nf90_put_var(file, major_id, radius_m))
        if (write_signgs) call nf90_check(nf90_put_var(file, sign_id, handedness))
        call nf90_check(nf90_put_var(file, phi_id, polarity*2.d0*pi*minor_m**2*radial))
        call nf90_check(nf90_put_var(file, iota_id, [( -handedness*iota, i=1,ns)]))
        call nf90_check(nf90_put_var(file, xm_id, [0.d0, 1.d0, 8.d0]))
        call nf90_check(nf90_put_var(file, xn_id, [0.d0, 0.d0, 2.d0]))
        call nf90_check(nf90_put_var(file, r_id, rc))
        call nf90_check(nf90_put_var(file, z_id, zs))
        call nf90_check(nf90_put_var(file, l_id, lambda))
        call nf90_check(nf90_close(file))
    end subroutine write_wout

    subroutine nf90_check(status)
        integer, intent(in) :: status

        if (status /= nf90_noerr) error stop nf90_strerror(status)
    end subroutine nf90_check
end program test_vmec_signed_flux
