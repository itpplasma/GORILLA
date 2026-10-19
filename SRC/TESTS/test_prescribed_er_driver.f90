! Oracle for the prescribed radial electric field loader.
! Er = E0*r_eff/r_max integrates to Phi_0(r_eff) = -E0*r_eff**2/(2*r_max).
! The mapping file gives psi_pol = c*r_eff**2, once ascending and once
! descending in psi_pol, so Phi_0 is linear in psi_pol. Linear interpolation
! of Er, the trapezoidal rule and the cubic spline are then all exact, and
! Phi_0 must match at and between the nodes to roundoff.
program test_prescribed_er_driver
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use strong_electric_field_mod, only: load_phi0_from_er, eval_phi0
    implicit none
    real(dp), parameter :: e0 = 2.5_dp, c_psi = 3.0e4_dp, r_max = 50.0_dp
    integer, parameter :: n_er = 7, n_eq = 101
    real(dp) :: sign_psi

    call write_er_file('er_const.dat')
    sign_psi = 1.0_dp
    call check_mapping('mapping_ascending.dat', sign_psi)
    sign_psi = -1.0_dp
    call check_mapping('mapping_descending.dat', sign_psi)
    print *, 'prescribed Er oracle: PASS'

contains

    subroutine write_er_file(filename)
        character(len=*), intent(in) :: filename
        integer :: iu, i
        open(newunit=iu, file=filename, status='replace', action='write')
        write(iu, '(a)') '# r_eff [cm]  Er [statV/cm]'
        do i = 1, n_er
            write(iu, '(2es25.16)') r_max*real(i - 1, dp)/real(n_er - 1, dp), &
                e0*real(i - 1, dp)/real(n_er - 1, dp)
        end do
        close(iu)
    end subroutine write_er_file

    subroutine check_mapping(filename, psi_sign)
        character(len=*), intent(in) :: filename
        real(dp), intent(in) :: psi_sign
        integer :: iu, i
        real(dp) :: r, psi, expected, err_node, err_mid

        open(newunit=iu, file=filename, status='replace', action='write')
        do i = 1, n_eq
            r = r_max*real(i - 1, dp)/real(n_eq - 1, dp)
            write(iu, '(5es25.16)') 200.0_dp + r, r, 1.5_dp, psi_sign*c_psi*r**2, 0.0_dp
        end do
        close(iu)

        call load_phi0_from_er('er_const.dat', filename)

        err_node = 0.0_dp
        err_mid = 0.0_dp
        do i = 1, n_eq
            r = r_max*real(i - 1, dp)/real(n_eq - 1, dp)
            psi = psi_sign*c_psi*r**2
            expected = -e0*r**2/(2.0_dp*r_max)
            err_node = max(err_node, abs(eval_phi0(psi) - expected))
            if (i < n_eq) then
                r = r + 0.5_dp*r_max/real(n_eq - 1, dp)
                psi = psi_sign*c_psi*r**2
                err_mid = max(err_mid, abs(eval_phi0(psi) + e0*r**2/(2.0_dp*r_max)))
            end if
        end do
        print *, trim(filename), ' max node error =', err_node, &
            ' max midpoint error =', err_mid
        if (err_node > 1.0e-12_dp*e0*r_max) error stop 'Phi_0 differs at the nodes'
        if (err_mid > 1.0e-12_dp*e0*r_max) error stop 'Phi_0 differs between the nodes'
        ! Outside the table the potential is held at the end values.
        if (abs(eval_phi0(psi_sign*2.0_dp*c_psi*r_max**2) + 0.5_dp*e0*r_max) > 1.0e-12_dp*e0*r_max) &
            error stop 'Phi_0 not clamped beyond the table'
    end subroutine check_mapping

end program test_prescribed_er_driver
