program test_jorek_restart_backend
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use jorek_field_backend_mod, only: load_jorek_field_backend, &
        evaluate_jorek_field_backend, free_jorek_field_backend
    use jorek_restart, only: jorek_restart_t, load_jorek_restart, &
        free_jorek_restart

    implicit none

    type(jorek_restart_t) :: restart
    character(len=1024) :: filename
    real(dp) :: a(3), b(3), bmod
    integer :: ierr

    call get_command_argument(1, filename)
    if (len_trim(filename) == 0) error stop 'golden restart path required'

    call load_jorek_restart(trim(filename), restart, ierr)
    if (ierr /= 0) error stop 'golden restart load failed'
    if (restart%jorek_model /= 303 .or. restart%n_nodes /= 6281 &
            .or. restart%n_elements /= 6129 .or. restart%n_tor /= 3) &
        error stop 'golden restart metadata mismatch'
    call require_close(restart%F0, 2.930189997653415_dp, 2.0e-12_dp)
    call free_jorek_restart(restart)

    call load_jorek_field_backend(trim(filename), ierr, 100.0_dp, 1.0e4_dp)
    if (ierr /= 0) error stop 'golden backend initialization failed'

    call evaluate_jorek_field_backend(171.047_dp, 0.0_dp, 3.863_dp, &
        a, b, bmod, ierr)
    call require_success(ierr)
    call require_vector(a, [0.0_dp, 2.6376570704692677e7_dp, &
        -1.5728327716843248e6_dp])
    call require_vector(b, [-7.8091258771564851_dp, &
        2.9301899972574934e6_dp, 7.3835745814089364e1_dp])
    call require_close(bmod, 1.7131066423046243e4_dp, 2.0e-9_dp)

    call evaluate_jorek_field_backend(171.047_dp, 1.0_dp, 3.863_dp, &
        a, b, bmod, ierr)
    call require_success(ierr)
    call require_vector(a, [0.0_dp, 2.6376954320458438e7_dp, &
        -1.5728327716843248e6_dp])
    call require_vector(b, [-7.8191018919215480_dp, &
        2.9301899972574934e6_dp, 7.3868651584250500e1_dp])
    call require_close(bmod, 1.7131066569453778e4_dp, 2.0e-9_dp)

    call evaluate_jorek_field_backend(190.0_dp, 1.0_dp, 0.0_dp, &
        a, b, bmod, ierr)
    call require_success(ierr)
    call require_vector(a, [0.0_dp, 2.2037693420978863e7_dp, &
        -1.8807538394451437e6_dp])
    call require_vector(b, [-3.7080170447481675e2_dp, &
        2.9301899954255959e6_dp, -2.3842691701895715e3_dp])
    call require_close(bmod, 1.5609674564466492e4_dp, 2.0e-9_dp)

    call evaluate_jorek_field_backend(300.0_dp, 0.0_dp, 0.0_dp, &
        a, b, bmod, ierr)
    if (ierr == 0 .or. any(a /= 0.0_dp) .or. any(b /= 0.0_dp) &
            .or. bmod /= 0.0_dp) &
        error stop 'outside-domain restart query was accepted'

    call free_jorek_field_backend()
    print '(A)', 'PASS: JOREK restart-backed field golden values'

contains

    subroutine require_success(status)
        integer, intent(in) :: status

        if (status /= 0) error stop 'restart-backed field query failed'
    end subroutine require_success

    subroutine require_vector(actual, expected)
        real(dp), intent(in) :: actual(3), expected(3)
        integer :: i

        do i = 1, 3
            call require_close(actual(i), expected(i), 2.0e-9_dp)
        end do
    end subroutine require_vector

    subroutine require_close(actual, expected, tolerance)
        real(dp), intent(in) :: actual, expected, tolerance

        if (abs(actual - expected) > tolerance*max(1.0_dp, abs(expected))) &
            error stop 'restart-backed field value mismatch'
    end subroutine require_close

end program test_jorek_restart_backend
