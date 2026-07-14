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
    call require_vector(a, [0.0_dp, 2.6379834772598185e7_dp, &
        -1.5728327712887882e6_dp])
    call require_vector(b, [-3.9579526048319588e-2_dp, &
        2.9301899976530299e6_dp, -1.8625597037654090e-1_dp])
    call require_close(bmod, 1.7130905527919473e4_dp, 2.0e-9_dp)

    call evaluate_jorek_field_backend(171.047_dp, 1.0_dp, 3.863_dp, &
        a, b, bmod, ierr)
    call require_success(ierr)
    call require_vector(a, [0.0_dp, 2.6380218701328982e7_dp, &
        -1.5728327712887882e6_dp])
    call require_vector(b, [-5.2962610719154167e-2_dp, &
        2.9301899976530299e6_dp, -1.6139042786267069e-1_dp])
    call require_close(bmod, 1.7130905527703315e4_dp, 2.0e-9_dp)

    call evaluate_jorek_field_backend(190.0_dp, 1.0_dp, 0.0_dp, &
        a, b, bmod, ierr)
    call require_success(ierr)
    call require_vector(a, [0.0_dp, 2.2046295916477572e7_dp, &
        -1.8807538368041692e6_dp])
    call require_vector(b, [-3.6767339193728509e2_dp, &
        2.9301899980665701e6_dp, -2.3250116147547578e3_dp])
    call require_close(bmod, 1.5600659274184742e4_dp, 2.0e-9_dp)

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

        if (maxval(abs(actual - expected) &
                /max(1.0_dp, abs(expected))) > 2.0e-9_dp) then
            write (*, '(A, 3(ES24.16E3, 1X))') 'actual:   ', actual
            write (*, '(A, 3(ES24.16E3, 1X))') 'expected: ', expected
            error stop 'restart-backed field vector mismatch'
        end if
    end subroutine require_vector

    subroutine require_close(actual, expected, tolerance)
        real(dp), intent(in) :: actual, expected, tolerance

        if (abs(actual - expected) > tolerance*max(1.0_dp, abs(expected))) then
            write (*, '(A, 2(ES24.16E3, 1X))') 'actual, expected: ', &
                actual, expected
            error stop 'restart-backed field value mismatch'
        end if
    end subroutine require_close

end program test_jorek_restart_backend
