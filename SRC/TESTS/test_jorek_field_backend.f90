program test_jorek_field_backend
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use jorek_bezier, only: jorek_locator_t, build_jorek_locator
    use jorek_field_backend_mod, only: evaluate_jorek_field_backend, &
        evaluate_jorek_model303_gorilla, evaluate_jorek_model303_gorilla_element
    use jorek_restart, only: jorek_restart_t

    implicit none

    type(jorek_restart_t) :: data, overlap
    type(jorek_locator_t) :: locator
    real(dp) :: a_covariant(3), b_covariant(3), bmod
    integer :: ierr, node
    real(dp), parameter :: r_node(4) = [10.0_dp, 11.0_dp, 11.0_dp, 10.0_dp]
    real(dp), parameter :: z_node(4) = [0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp]

    call evaluate_jorek_field_backend(10.25_dp, 0.4_dp, 0.6_dp, &
        a_covariant, b_covariant, bmod, ierr)
    if (ierr /= 9) error stop 'uninitialized JOREK backend was accepted'

    call evaluate_jorek_model303_gorilla(data, 10.25_dp, 0.4_dp, 0.6_dp, &
        a_covariant, b_covariant, bmod, ierr, length_scale_in=0.0_dp)
    if (ierr /= 10) error stop 'invalid JOREK conversion scale was accepted'

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
    data%size(:, 2:3, 2) = -1.0_dp
    data%size(:, 3:4, 3) = -1.0_dp
    data%size(:, 2, 4) = -1.0_dp
    data%size(:, 4, 4) = -1.0_dp
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

    call build_jorek_locator(data, locator, ierr)
    if (ierr /= 0) error stop 'JOREK locator construction failed'
    call evaluate_jorek_model303_gorilla(data, 1025.0_dp, 0.4_dp, 60.0_dp, &
        a_covariant, b_covariant, bmod, ierr, locator, 100.0_dp, 1.0e4_dp)
    if (ierr /= 0) error stop 'JOREK backend returned an error'
    call check(a_covariant, [0.0_dp, -6.15e8_dp, &
        -20.0e6_dp*log(10.25_dp)])
    call check(b_covariant, [1.0e4_dp, 20.0e6_dp, &
        -0.6e4_dp/10.25_dp])
    if (abs(bmod - 1.0e4_dp*sqrt(1.0_dp + (20.0_dp/10.25_dp)**2 &
            + (0.6_dp/10.25_dp)**2)) > 3.0e-10_dp) &
        error stop 'JOREK backend returned the wrong field magnitude'

    call initialize_overlap(overlap)
    call evaluate_jorek_model303_gorilla_element(overlap, 2, 0.25_dp, &
        0.6_dp, 0.4_dp, 1025.0_dp, a_covariant, b_covariant, bmod, ierr, &
        100.0_dp, 1.0e4_dp)
    if (ierr /= 0) error stop 'explicit overlapping owner returned an error'
    call check(a_covariant, [0.0_dp, -12.3e8_dp, &
        -20.0e6_dp*log(10.25_dp)])
    call check(b_covariant, [2.0e4_dp, 20.0e6_dp, &
        -1.2e4_dp/10.25_dp])
    if (abs(bmod - 1.0e4_dp*sqrt(4.0_dp + (20.0_dp/10.25_dp)**2 &
            + (1.2_dp/10.25_dp)**2)) > 3.0e-10_dp) &
        error stop 'explicit overlapping owner returned the wrong magnitude'

    call evaluate_jorek_model303_gorilla_element(data, 1, 0.25_dp, 0.6_dp, &
        0.4_dp, 1025.0_dp, a_covariant, b_covariant, bmod, ierr, &
        100.0_dp, 1.0e4_dp)
    if (ierr /= 0) error stop 'owner-aware JOREK backend returned an error'
    call check(a_covariant, [0.0_dp, -6.15e8_dp, &
        -20.0e6_dp*log(10.25_dp)])
    call check(b_covariant, [1.0e4_dp, 20.0e6_dp, &
        -0.6e4_dp/10.25_dp])

    print '(A)', 'PASS: JOREK model-303 GORILLA field conventions'

contains

    subroutine initialize_overlap(overlap)
        type(jorek_restart_t), intent(out) :: overlap

        real(dp) :: factor
        integer :: element, local_node, node

        overlap%jorek_model = 303
        overlap%n_order = 3
        overlap%n_degrees = 4
        overlap%n_coord_tor = 1
        overlap%n_dim = 2
        overlap%n_tor = 1
        overlap%n_period = 1
        overlap%n_var = 7
        overlap%n_nodes = 8
        overlap%n_elements = 2
        overlap%n_vertex_max = 4
        overlap%F0 = 20.0_dp
        allocate (overlap%x(8, 1, 4, 2), overlap%values(8, 1, 4, 7))
        allocate (overlap%vertex(2, 4), overlap%size(2, 4, 4))
        overlap%x = 0.0_dp
        overlap%values = 0.0_dp
        overlap%size = 1.0_dp
        overlap%size(:, 2:3, 2) = -1.0_dp
        overlap%size(:, 3:4, 3) = -1.0_dp
        overlap%size(:, 2, 4) = -1.0_dp
        overlap%size(:, 4, 4) = -1.0_dp
        do element = 1, 2
            factor = real(element, dp)
            do local_node = 1, 4
                node = 4*(element - 1) + local_node
                overlap%vertex(element, local_node) = node
                overlap%x(node, 1, 1, 1) = r_node(local_node)
                overlap%x(node, 1, 1, 2) = z_node(local_node)
                overlap%x(node, 1, 2, 1) = 1.0_dp/3.0_dp
                overlap%x(node, 1, 3, 2) = 1.0_dp/3.0_dp
                overlap%values(node, 1, 1, 1) = &
                    factor*r_node(local_node)*z_node(local_node)
                overlap%values(node, 1, 2, 1) = &
                    factor*z_node(local_node)/3.0_dp
                overlap%values(node, 1, 3, 1) = &
                    factor*r_node(local_node)/3.0_dp
                overlap%values(node, 1, 4, 1) = factor/9.0_dp
            end do
        end do
    end subroutine initialize_overlap

    subroutine check(actual, expected)
        real(dp), intent(in) :: actual(3), expected(3)

        if (maxval(abs(actual - expected)) > &
                1.0e-12_dp*max(1.0_dp, maxval(abs(expected)))) &
            error stop 'JOREK backend component mismatch'
    end subroutine check

end program test_jorek_field_backend
