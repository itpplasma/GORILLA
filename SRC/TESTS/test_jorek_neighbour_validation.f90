program test_jorek_neighbour_validation
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use jorek_refined_plane_mod, only: extract_refined_jorek_plane
    use jorek_restart, only: jorek_restart_t

    implicit none

    type(jorek_restart_t) :: data
    real(dp), allocatable :: plane_rz(:, :), psi(:), vertex_st(:, :)
    integer, allocatable :: triangles(:, :), vertex_element(:)
    integer :: ierr

    call build_synthetic_restart(data)

    call extract_refined_jorek_plane(data, 2, plane_rz, psi, triangles, &
        vertex_element, vertex_st, ierr)
    if (ierr /= 0) error stop 'baseline extraction without neighbours failed'
    if (size(plane_rz, 2) /= 9 .or. size(triangles, 1) /= 8) &
        error stop 'baseline refined plane has wrong dimensions'

    allocate(data%neighbours(data%n_elements, 3), source=0)
    call extract_refined_jorek_plane(data, 2, plane_rz, psi, triangles, &
        vertex_element, vertex_st, ierr)
    if (ierr /= 12) error stop 'neighbours with malformed side count was accepted'
    deallocate(data%neighbours)

    allocate(data%neighbours(data%n_elements + 1, 4), source=0)
    call extract_refined_jorek_plane(data, 2, plane_rz, psi, triangles, &
        vertex_element, vertex_st, ierr)
    if (ierr /= 12) error stop 'neighbours with malformed element count was accepted'
    deallocate(data%neighbours)

    allocate(data%neighbours(data%n_elements, 4), source=0)
    call extract_refined_jorek_plane(data, 2, plane_rz, psi, triangles, &
        vertex_element, vertex_st, ierr)
    if (ierr /= 0) error stop 'well-formed boundary-only neighbours were rejected'

    data%neighbours(1, 1) = data%n_elements + 1
    call extract_refined_jorek_plane(data, 2, plane_rz, psi, triangles, &
        vertex_element, vertex_st, ierr)
    if (ierr /= 10) error stop 'out-of-range neighbour reference was accepted'

    print '(A)', 'PASS: JOREK neighbour-array validation'

contains

    subroutine build_synthetic_restart(data)
        type(jorek_restart_t), intent(out) :: data

        real(dp), parameter :: r_node(4) = [10.0_dp, 11.0_dp, 11.0_dp, 10.0_dp]
        real(dp), parameter :: z_node(4) = [0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp]
        integer :: node

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
        allocate(data%x(4, 1, 4, 2), data%values(4, 1, 4, 7))
        allocate(data%vertex(1, 4), data%size(1, 4, 4))
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
    end subroutine build_synthetic_restart

end program test_jorek_neighbour_validation
