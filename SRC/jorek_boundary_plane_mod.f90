module jorek_boundary_plane_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use jorek_bezier, only: evaluate_jorek_geometry
    use jorek_field_values, only: evaluate_jorek_variable
    use jorek_refined_plane_mod, only: extract_refined_jorek_plane, &
        orient_triangles
    use jorek_restart, only: jorek_restart_t

    implicit none
    private

    public :: extract_boundary_refined_jorek_plane

contains

    subroutine extract_boundary_refined_jorek_plane(data, subdivisions, &
            plane_rz, psi, triangles, vertex_element, vertex_st, &
            element_triangle_range, ierr)
        type(jorek_restart_t), intent(in) :: data
        integer, intent(in) :: subdivisions
        real(dp), allocatable, intent(out) :: plane_rz(:, :), psi(:)
        integer, allocatable, intent(out) :: triangles(:, :), vertex_element(:)
        real(dp), allocatable, intent(out) :: vertex_st(:, :)
        integer, allocatable, intent(out) :: element_triangle_range(:, :)
        integer, intent(out) :: ierr

        integer, allocatable :: base_triangles(:, :), element_vertex(:, :, :)
        integer, allocatable :: cell_center(:, :, :), extra_node(:, :, :)

        ierr = 1
        if (subdivisions < 4 .or. mod(subdivisions, 2) /= 0) return
        if (.not. allocated(data%neighbours)) return
        if (any(shape(data%neighbours) /= [data%n_elements, 4])) return
        call extract_refined_jorek_plane(data, 2, plane_rz, psi, &
            base_triangles, vertex_element, vertex_st, ierr, element_vertex)
        if (ierr /= 0) return
        deallocate(base_triangles)
        call append_boundary_nodes(data, subdivisions, &
            plane_rz, psi, vertex_element, vertex_st, extra_node, ierr)
        if (ierr /= 0) return
        call append_boundary_cell_centers(data, plane_rz, psi, vertex_element, &
            vertex_st, cell_center, ierr)
        if (ierr /= 0) return
        call build_boundary_triangles(data, subdivisions, element_vertex, &
            extra_node, cell_center, triangles, element_triangle_range)
        call orient_triangles(plane_rz, triangles, ierr)
    end subroutine extract_boundary_refined_jorek_plane

    subroutine append_boundary_nodes(data, subdivisions, &
            plane_rz, psi, vertex_element, vertex_st, extra_node, ierr)
        type(jorek_restart_t), intent(in) :: data
        integer, intent(in) :: subdivisions
        real(dp), allocatable, intent(inout) :: plane_rz(:, :), psi(:)
        integer, allocatable, intent(inout) :: vertex_element(:)
        real(dp), allocatable, intent(inout) :: vertex_st(:, :)
        integer, allocatable, intent(out) :: extra_node(:, :, :)
        integer, intent(out) :: ierr

        real(dp), allocatable :: extended_rz(:, :), extended_psi(:)
        real(dp), allocatable :: extended_st(:, :)
        integer, allocatable :: extended_element(:)
        integer :: new_count, node

        new_count = count(data%neighbours == 0)*(subdivisions - 2)
        allocate(extended_rz(2, size(psi) + new_count))
        allocate(extended_psi(size(psi) + new_count))
        allocate(extended_element(size(psi) + new_count))
        allocate(extended_st(2, size(psi) + new_count))
        extended_rz(:, :size(psi)) = plane_rz
        extended_psi(:size(psi)) = psi
        extended_element(:size(psi)) = vertex_element
        extended_st(:, :size(psi)) = vertex_st
        allocate(extra_node(4, subdivisions - 1, data%n_elements), source=0)
        node = size(psi)
        call populate_boundary_nodes(data, subdivisions, extra_node, &
            extended_rz, extended_psi, extended_element, extended_st, &
            node, ierr)
        if (ierr /= 0) return
        if (node /= size(extended_psi)) then
            ierr = 2
            return
        end if
        call move_alloc(extended_rz, plane_rz)
        call move_alloc(extended_psi, psi)
        call move_alloc(extended_element, vertex_element)
        call move_alloc(extended_st, vertex_st)
        ierr = 0
    end subroutine append_boundary_nodes

    subroutine populate_boundary_nodes(data, subdivisions, extra_node, &
            plane_rz, psi, vertex_element, vertex_st, node, ierr)
        type(jorek_restart_t), intent(in) :: data
        integer, intent(in) :: subdivisions
        integer, intent(inout) :: extra_node(:, :, :), vertex_element(:), node
        real(dp), intent(inout) :: plane_rz(:, :), psi(:), vertex_st(:, :)
        integer, intent(out) :: ierr

        integer :: element, k, side

        do element = 1, data%n_elements
            do side = 1, 4
                if (data%neighbours(element, side) /= 0) cycle
                do k = 1, subdivisions - 1
                    if (k == subdivisions/2) cycle
                    node = node + 1
                    extra_node(side, k, element) = node
                    call store_boundary_node(data, element, side, k, &
                        subdivisions, plane_rz(:, node), psi(node), &
                        vertex_st(:, node), ierr)
                    if (ierr /= 0) return
                    vertex_element(node) = element
                end do
            end do
        end do
        ierr = 0
    end subroutine populate_boundary_nodes

    subroutine append_boundary_cell_centers(data, plane_rz, psi, &
            vertex_element, vertex_st, cell_center, ierr)
        type(jorek_restart_t), intent(in) :: data
        real(dp), allocatable, intent(inout) :: plane_rz(:, :), psi(:)
        integer, allocatable, intent(inout) :: vertex_element(:)
        real(dp), allocatable, intent(inout) :: vertex_st(:, :)
        integer, allocatable, intent(out) :: cell_center(:, :, :)
        integer, intent(out) :: ierr

        real(dp), allocatable :: extended_rz(:, :), extended_psi(:)
        real(dp), allocatable :: extended_st(:, :)
        integer, allocatable :: extended_element(:)
        real(dp) :: st(2)
        integer :: element, i, j, new_count, node

        new_count = count_boundary_cells(data)
        allocate(extended_rz(2, size(psi) + new_count))
        allocate(extended_psi(size(psi) + new_count))
        allocate(extended_element(size(psi) + new_count))
        allocate(extended_st(2, size(psi) + new_count))
        extended_rz(:, :size(psi)) = plane_rz
        extended_psi(:size(psi)) = psi
        extended_element(:size(psi)) = vertex_element
        extended_st(:, :size(psi)) = vertex_st
        allocate(cell_center(0:1, 0:1, data%n_elements), source=0)
        node = size(psi)
        do element = 1, data%n_elements
            do j = 0, 1
                do i = 0, 1
                    if (.not. is_boundary_cell(data, element, i, j)) cycle
                    node = node + 1
                    cell_center(i, j, element) = node
                    st = [0.5_dp*real(i, dp) + 0.25_dp, &
                        0.5_dp*real(j, dp) + 0.25_dp]
                    call store_element_node(data, element, st, &
                        extended_rz(:, node), extended_psi(node), ierr)
                    if (ierr /= 0) return
                    extended_element(node) = element
                    extended_st(:, node) = st
                end do
            end do
        end do
        call move_alloc(extended_rz, plane_rz)
        call move_alloc(extended_psi, psi)
        call move_alloc(extended_element, vertex_element)
        call move_alloc(extended_st, vertex_st)
        ierr = 0
    end subroutine append_boundary_cell_centers

    subroutine store_boundary_node(data, element, side, k, subdivisions, &
            rz, value, st, ierr)
        type(jorek_restart_t), intent(in) :: data
        integer, intent(in) :: element, side, k, subdivisions
        real(dp), intent(out) :: rz(2), value, st(2)
        integer, intent(out) :: ierr

        real(dp) :: fraction

        fraction = real(k, dp)/subdivisions
        call side_parameters(side, fraction, st)
        call store_element_node(data, element, st, rz, value, ierr)
    end subroutine store_boundary_node

    subroutine store_element_node(data, element, st, rz, value, ierr)
        type(jorek_restart_t), intent(in) :: data
        integer, intent(in) :: element
        real(dp), intent(in) :: st(2)
        real(dp), intent(out) :: rz(2), value
        integer, intent(out) :: ierr

        real(dp) :: derivatives(3), rz_st(2, 2)

        call evaluate_jorek_geometry(data, element, st(1), st(2), rz, &
            rz_st, ierr)
        if (ierr /= 0) return
        call evaluate_jorek_variable(data, element, 1, st(1), st(2), &
            0.0_dp, value, derivatives, ierr)
    end subroutine store_element_node

    subroutine build_boundary_triangles(data, subdivisions, element_vertex, &
            extra_node, cell_center, triangles, element_triangle_range)
        type(jorek_restart_t), intent(in) :: data
        integer, intent(in) :: subdivisions, element_vertex(0:, 0:, :)
        integer, intent(in) :: extra_node(:, :, :), cell_center(0:, 0:, :)
        integer, allocatable, intent(out) :: triangles(:, :)
        integer, allocatable, intent(out) :: element_triangle_range(:, :)

        integer, allocatable :: work(:, :)
        integer :: element, maximum, triangle

        maximum = 8*data%n_elements + 10*count(data%neighbours == 0)
        allocate(work(maximum, 3))
        allocate(element_triangle_range(2, data%n_elements))
        triangle = 0
        do element = 1, data%n_elements
            element_triangle_range(1, element) = triangle + 1
            call append_element_cells(data, element, subdivisions, &
                element_vertex, extra_node, cell_center, work, triangle)
            element_triangle_range(2, element) = triangle
        end do
        allocate(triangles(triangle, 3), source=work(:triangle, :))
    end subroutine build_boundary_triangles

    subroutine append_element_cells(data, element, subdivisions, element_vertex, &
            extra_node, cell_center, triangles, triangle)
        type(jorek_restart_t), intent(in) :: data
        integer, intent(in) :: element, subdivisions
        integer, intent(in) :: element_vertex(0:, 0:, :)
        integer, intent(in) :: extra_node(:, :, :), cell_center(0:, 0:, :)
        integer, intent(inout) :: triangles(:, :), triangle

        integer :: i, j

        do j = 0, 1
            do i = 0, 1
                if (is_boundary_cell(data, element, i, j)) then
                    call append_boundary_cell_fan(data, element, i, j, &
                        subdivisions, element_vertex, extra_node, &
                        cell_center(i, j, element), triangles, triangle)
                else
                    call append_factor_two_cell(element, i, j, &
                        element_vertex, triangles, triangle)
                end if
            end do
        end do
    end subroutine append_element_cells

    subroutine append_factor_two_cell(element, i, j, element_vertex, &
            triangles, triangle)
        integer, intent(in) :: element, i, j, element_vertex(0:, 0:, :)
        integer, intent(inout) :: triangles(:, :), triangle

        call append_triangle([element_vertex(i, j, element), &
            element_vertex(i + 1, j, element), &
            element_vertex(i + 1, j + 1, element)], triangles, triangle)
        call append_triangle([element_vertex(i, j, element), &
            element_vertex(i + 1, j + 1, element), &
            element_vertex(i, j + 1, element)], triangles, triangle)
    end subroutine append_factor_two_cell

    subroutine append_boundary_cell_fan(data, element, i, j, subdivisions, &
            element_vertex, extra_node, center, triangles, triangle)
        type(jorek_restart_t), intent(in) :: data
        integer, intent(in) :: element, i, j, subdivisions, center
        integer, intent(in) :: element_vertex(0:, 0:, :), extra_node(:, :, :)
        integer, intent(inout) :: triangles(:, :), triangle

        integer :: current, side

        current = element_vertex(i, j, element)
        do side = 1, 4
            call append_cell_edge(data, element, i, j, side, subdivisions, &
                element_vertex, extra_node, center, current, triangles, triangle)
        end do
    end subroutine append_boundary_cell_fan

    subroutine append_cell_edge(data, element, i, j, side, subdivisions, &
            element_vertex, extra_node, center, current, triangles, triangle)
        type(jorek_restart_t), intent(in) :: data
        integer, intent(in) :: element, i, j, side, subdivisions, center
        integer, intent(in) :: element_vertex(0:, 0:, :), extra_node(:, :, :)
        integer, intent(inout) :: current, triangles(:, :), triangle

        integer :: first, k, last, next

        if (cell_edge_is_physical(data, element, i, j, side)) then
            select case (side)
            case (1)
                first = i*subdivisions/2
            case (2)
                first = j*subdivisions/2
            case (3)
                first = (1 - i)*subdivisions/2
            case (4)
                first = (1 - j)*subdivisions/2
            end select
            last = first + subdivisions/2
            do k = first + 1, last
                next = boundary_node(element, side, k, subdivisions, &
                    element_vertex, extra_node)
                call append_triangle([center, current, next], triangles, triangle)
                current = next
            end do
        else
            next = cell_edge_endpoint(element, i, j, side, element_vertex)
            call append_triangle([center, current, next], triangles, triangle)
            current = next
        end if
    end subroutine append_cell_edge

    logical function cell_edge_is_physical(data, element, i, j, side)
        type(jorek_restart_t), intent(in) :: data
        integer, intent(in) :: element, i, j, side

        select case (side)
        case (1)
            cell_edge_is_physical = j == 0
        case (2)
            cell_edge_is_physical = i == 1
        case (3)
            cell_edge_is_physical = j == 1
        case (4)
            cell_edge_is_physical = i == 0
        end select
        cell_edge_is_physical = cell_edge_is_physical &
            .and. data%neighbours(element, side) == 0
    end function cell_edge_is_physical

    integer function cell_edge_endpoint(element, i, j, side, element_vertex)
        integer, intent(in) :: element, i, j, side
        integer, intent(in) :: element_vertex(0:, 0:, :)

        select case (side)
        case (1)
            cell_edge_endpoint = element_vertex(i + 1, j, element)
        case (2)
            cell_edge_endpoint = element_vertex(i + 1, j + 1, element)
        case (3)
            cell_edge_endpoint = element_vertex(i, j + 1, element)
        case (4)
            cell_edge_endpoint = element_vertex(i, j, element)
        end select
    end function cell_edge_endpoint

    integer function count_boundary_cells(data)
        type(jorek_restart_t), intent(in) :: data

        integer :: element, i, j

        count_boundary_cells = 0
        do element = 1, data%n_elements
            do j = 0, 1
                do i = 0, 1
                    if (is_boundary_cell(data, element, i, j)) &
                        count_boundary_cells = count_boundary_cells + 1
                end do
            end do
        end do
    end function count_boundary_cells

    logical function is_boundary_cell(data, element, i, j)
        type(jorek_restart_t), intent(in) :: data
        integer, intent(in) :: element, i, j

        is_boundary_cell = (j == 0 .and. data%neighbours(element, 1) == 0) &
            .or. (i == 1 .and. data%neighbours(element, 2) == 0) &
            .or. (j == 1 .and. data%neighbours(element, 3) == 0) &
            .or. (i == 0 .and. data%neighbours(element, 4) == 0)
    end function is_boundary_cell

    subroutine append_triangle(nodes, triangles, triangle)
        integer, intent(in) :: nodes(3)
        integer, intent(inout) :: triangles(:, :), triangle

        if (any(nodes <= 0)) return
        if (nodes(1) == nodes(2) .or. nodes(1) == nodes(3) &
                .or. nodes(2) == nodes(3)) return
        triangle = triangle + 1
        triangles(triangle, :) = nodes
    end subroutine append_triangle

    integer function boundary_node(element, side, k, subdivisions, &
            element_vertex, extra_node)
        integer, intent(in) :: element, side, k, subdivisions
        integer, intent(in) :: element_vertex(0:, 0:, :), extra_node(:, :, :)

        if (k == subdivisions/2) then
            boundary_node = side_midpoint(element, side, element_vertex)
        else if (k == subdivisions) then
            boundary_node = side_endpoint(element, side, element_vertex)
        else
            boundary_node = extra_node(side, k, element)
        end if
    end function boundary_node

    integer function side_midpoint(element, side, element_vertex)
        integer, intent(in) :: element, side, element_vertex(0:, 0:, :)

        select case (side)
        case (1)
            side_midpoint = element_vertex(1, 0, element)
        case (2)
            side_midpoint = element_vertex(2, 1, element)
        case (3)
            side_midpoint = element_vertex(1, 2, element)
        case (4)
            side_midpoint = element_vertex(0, 1, element)
        end select
    end function side_midpoint

    integer function side_endpoint(element, side, element_vertex)
        integer, intent(in) :: element, side, element_vertex(0:, 0:, :)

        select case (side)
        case (1)
            side_endpoint = element_vertex(2, 0, element)
        case (2)
            side_endpoint = element_vertex(2, 2, element)
        case (3)
            side_endpoint = element_vertex(0, 2, element)
        case (4)
            side_endpoint = element_vertex(0, 0, element)
        end select
    end function side_endpoint

    subroutine side_parameters(side, fraction, st)
        integer, intent(in) :: side
        real(dp), intent(in) :: fraction
        real(dp), intent(out) :: st(2)

        select case (side)
        case (1)
            st = [fraction, 0.0_dp]
        case (2)
            st = [1.0_dp, fraction]
        case (3)
            st = [1.0_dp - fraction, 1.0_dp]
        case (4)
            st = [0.0_dp, 1.0_dp - fraction]
        end select
    end subroutine side_parameters

end module jorek_boundary_plane_mod
