program test_jorek_restart_smoke
    !> Smoke test of the JOREK restart chain against the committed circular
    !> fixture SRC/TESTS/jorek_data/circular_equilibrium_restart.h5 (see the
    !> README there for provenance). It checks the reader metadata and the
    !> restart-backed field plumbing without golden field values, so it works
    !> on any revision that provides the JOREK backend.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use jorek_field_backend_mod, only: load_jorek_field_backend, &
        evaluate_jorek_field_backend, free_jorek_field_backend
    use jorek_restart, only: jorek_restart_t, load_jorek_restart, &
        free_jorek_restart

    implicit none

    type(jorek_restart_t) :: restart
    character(len=1024) :: filename
    real(dp) :: a_covariant(3), b_covariant(3), bmod
    integer :: ierr

    call get_command_argument(1, filename)
    if (len_trim(filename) == 0) error stop 'fixture restart path required'

    call load_jorek_restart(trim(filename), restart, ierr)
    if (ierr /= 0) error stop 'fixture restart load failed'
    if (restart%jorek_model /= 303) error stop 'fixture model mismatch'
    if (restart%n_nodes /= 490) error stop 'fixture node count mismatch'
    if (restart%n_elements /= 476) error stop 'fixture element count mismatch'
    if (restart%n_var /= 7) error stop 'fixture variable count mismatch'
    if (restart%n_degrees /= 4) error stop 'fixture degree count mismatch'
    if (restart%n_tor /= 1) error stop 'fixture toroidal harmonic mismatch'
    if (restart%n_order /= 3) error stop 'fixture Bezier order mismatch'
    if (abs(restart%F0 - 10.0_dp) > 1.0e-12_dp) error stop 'fixture F0 mismatch'
    call free_jorek_restart(restart)

    call load_jorek_field_backend(trim(filename), ierr)
    if (ierr /= 0) error stop 'fixture backend initialization failed'

    ! The fixture is centred at R_geo = 10 with minor radius 1 (JOREK units).
    call evaluate_jorek_field_backend(10.2_dp, 0.3_dp, 0.1_dp, &
        a_covariant, b_covariant, bmod, ierr)
    if (ierr /= 0) error stop 'inside-domain fixture query failed'
    if (.not. all(ieee_is_finite(a_covariant))) &
        error stop 'fixture vector potential is not finite'
    if (.not. all(ieee_is_finite(b_covariant))) &
        error stop 'fixture magnetic field is not finite'
    if (.not. ieee_is_finite(bmod)) error stop 'fixture bmod is not finite'
    if (bmod <= 0.0_dp) error stop 'fixture field magnitude is not positive'

    call evaluate_jorek_field_backend(30.0_dp, 0.0_dp, 0.0_dp, &
        a_covariant, b_covariant, bmod, ierr)
    if (ierr == 0) error stop 'outside-domain fixture query was accepted'

    call free_jorek_field_backend()
    print '(A)', 'PASS: JOREK circular-fixture restart smoke'
end program test_jorek_restart_smoke
