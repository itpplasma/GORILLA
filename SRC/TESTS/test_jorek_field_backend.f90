program test_jorek_field_backend
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use jorek_field_backend_mod, only: evaluate_jorek_field_backend, &
        evaluate_jorek_model303_gorilla
    use jorek_restart, only: jorek_restart_t

    implicit none

    type(jorek_restart_t) :: data
    real(dp) :: a_covariant(3), b_covariant(3), bmod
    integer :: ierr, node
    real(dp), parameter :: r_node(4) = [10.0_dp, 11.0_dp, 11.0_dp, 10.0_dp]
    real(dp), parameter :: z_node(4) = [0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp]

    call evaluate_jorek_field_backend(10.25_dp, 0.4_dp, 0.6_dp, &
        a_covariant, b_covariant, bmod, ierr)
    if (ierr /= 9) error stop 'uninitialized JOREK backend was accepted'

    data%jorek_model = 303
    data%n_order = 3
    data%n_degrees = 4
    data%n_coord_tor = 1
    data%n_dim = 2
    data%n_tor = 1
    data%n_period = 1
    data%n_var = 7
    data%n_nodes = 4
    data%n_elements = 1
    data%n_vertex_max = 4
    data%F0 = 20.0_dp
    allocate (data%x(4, 1, 4, 2), data%values(4, 1, 4, 7))
    allocate (data%vertex(1, 4), data%size(1, 4, 4))
    data%x = 0.0_dp
    data%values = 0.0_dp
    data%vertex(1, :) = [1, 2, 3, 4]
    data%size = 1.0_dp
    do node = 1, 4
        data%x(node, 1, 1, 1) = r_node(node)
        data%x(node, 1, 1, 2) = z_node(node)
        data%x(node, 1, 2, 1) = 1.0_dp/3.0_dp
        data%x(node, 1, 3, 2) = 1.0_dp/3.0_dp
        data%values(node, 1, 1, 1) = r_node(node)*z_node(node)
        data%values(node, 1, 2, 1) = z_node(node)/3.0_dp
        data%values(node, 1, 3, 1) = r_node(node)/3.0_dp
        data%values(node, 1, 4, 1) = 1.0_dp/9.0_dp
    end do

    call evaluate_jorek_model303_gorilla(data, 10.25_dp, 0.4_dp, 0.6_dp, &
        a_covariant, b_covariant, bmod, ierr)
    if (ierr /= 0) error stop 'JOREK backend returned an error'
    call check(a_covariant, [0.0_dp, -6.15_dp, -20.0_dp*log(10.25_dp)])
    call check(b_covariant, [1.0_dp, 20.0_dp, -0.6_dp/10.25_dp])
    if (abs(bmod - sqrt(1.0_dp + (20.0_dp/10.25_dp)**2 &
            + (0.6_dp/10.25_dp)**2)) > 3.0e-14_dp) &
        error stop 'JOREK backend returned the wrong field magnitude'

    print '(A)', 'PASS: JOREK model-303 GORILLA field conventions'

contains

    subroutine check(actual, expected)
        real(dp), intent(in) :: actual(3), expected(3)

        if (maxval(abs(actual - expected)) > 3.0e-14_dp) &
            error stop 'JOREK backend component mismatch'
    end subroutine check

end program test_jorek_field_backend
