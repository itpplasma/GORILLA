module jorek_inserted_plane_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use jorek_bezier, only: evaluate_jorek_geometry
    use jorek_field_values, only: evaluate_jorek_variable
    use jorek_refined_plane_mod, only: extract_refined_jorek_plane, &
        orient_triangles
    use jorek_restart, only: jorek_restart_t

    implicit none
    private

    public :: extract_inserted_side_jorek_plane

contains

    subroutine extract_inserted_side_jorek_plane(data, subdivisions, &
            selected_side, plane_rz, psi, triangles, vertex_element, vertex_st, &
            triangle_parent, element_vertex, ierr)
        type(jorek_restart_t), intent(in) :: data
        integer, intent(in) :: subdivisions
        logical, intent(in) :: selected_side(:, :)
        real(dp), allocatable, intent(out) :: plane_rz(:, :), psi(:)
        integer, allocatable, intent(out) :: triangles(:, :), vertex_element(:)
        real(dp), allocatable, intent(out) :: vertex_st(:, :)
        integer, allocatable, intent(out) :: triangle_parent(:)
        integer, allocatable, intent(out) :: element_vertex(:, :, :)
        integer, intent(out) :: ierr

        integer, allocatable :: parent_work(:), triangle_work(:, :)
        integer :: extra_count, node, triangle_count

        ierr = 1
        if (subdivisions < 4 .or. mod(subdivisions, 2) /= 0) return
        if (.not. allocated(data%neighbours)) return
        if (any(shape(selected_side) /= [4, data%n_elements])) return
        call extract_refined_jorek_plane(data, 2, plane_rz, psi, triangles, &
            vertex_element, vertex_st, ierr, element_vertex)
        if (ierr /= 0) return
        allocate(triangle_parent(size(triangles, 1)))
        call assign_factor_two_parents(data, element_vertex, triangles, &
            triangle_parent, ierr)
        if (ierr /= 0) return
        extra_count = count_inserted_nodes(data, subdivisions, selected_side, ierr)
        if (ierr /= 0) return
        call extend_vertices(extra_count, plane_rz, psi, vertex_element, vertex_st)
        allocate(triangle_work(size(triangles, 1) + 2*extra_count, 3))
        allocate(parent_work(size(triangle_parent) + 2*extra_count))
        triangle_count = size(triangles, 1)
        triangle_work(:triangle_count, :) = triangles
        parent_work(:triangle_count) = triangle_parent
        node = size(psi) - extra_count
        call insert_selected_nodes(data, subdivisions, selected_side, plane_rz, &
            psi, vertex_element, vertex_st, triangle_work, parent_work, &
            triangle_count, node, ierr)
        if (ierr /= 0 .or. node /= size(psi)) return
        deallocate(triangles, triangle_parent)
        allocate(triangles(triangle_count, 3), &
            source=triangle_work(:triangle_count, :))
        allocate(triangle_parent(triangle_count), &
            source=parent_work(:triangle_count))
        call orient_triangles(plane_rz, triangles, ierr)
    end subroutine extract_inserted_side_jorek_plane

    integer function count_inserted_nodes(data, subdivisions, selected_side, &
            ierr) result(count_nodes)
        type(jorek_restart_t), intent(in) :: data
        integer, intent(in) :: subdivisions
        logical, intent(in) :: selected_side(:, :)
        integer, intent(out) :: ierr

        integer :: element, neighbor, neighbor_side, side

        count_nodes = 0
        ierr = 1
        do element = 1, data%n_elements
            do side = 1, 4
                if (.not. selected_side(side, element)) cycle
                neighbor = data%neighbours(element, side)
                if (neighbor > 0) then
                    neighbor_side = reciprocal_side(data, neighbor, element)
                    if (neighbor_side == 0) return
                    if (.not. selected_side(neighbor_side, neighbor)) return
                    if (neighbor < element) cycle
                end if
                count_nodes = count_nodes + subdivisions - 2
            end do
        end do
        ierr = 0
    end function count_inserted_nodes

    subroutine extend_vertices(extra_count, plane_rz, psi, vertex_element, &
            vertex_st)
        integer, intent(in) :: extra_count
        real(dp), allocatable, intent(inout) :: plane_rz(:, :), psi(:)
        integer, allocatable, intent(inout) :: vertex_element(:)
        real(dp), allocatable, intent(inout) :: vertex_st(:, :)

        real(dp), allocatable :: extended_rz(:, :), extended_psi(:)
        real(dp), allocatable :: extended_st(:, :)
        integer, allocatable :: extended_element(:)
        integer :: old_count

        old_count = size(psi)
        allocate(extended_rz(2, old_count + extra_count), &
            extended_psi(old_count + extra_count), &
            extended_element(old_count + extra_count), &
            extended_st(2, old_count + extra_count))
        extended_rz(:, :old_count) = plane_rz
        extended_psi(:old_count) = psi
        extended_element(:old_count) = vertex_element
        extended_st(:, :old_count) = vertex_st
        call move_alloc(extended_rz, plane_rz)
        call move_alloc(extended_psi, psi)
        call move_alloc(extended_element, vertex_element)
        call move_alloc(extended_st, vertex_st)
    end subroutine extend_vertices

    subroutine insert_selected_nodes(data, subdivisions, selected_side, &
            plane_rz, psi, vertex_element, vertex_st, triangles, parents, &
            triangle_count, node, ierr)
        type(jorek_restart_t), intent(in) :: data
        integer, intent(in) :: subdivisions
        logical, intent(in) :: selected_side(:, :)
        real(dp), intent(inout) :: plane_rz(:, :), psi(:), vertex_st(:, :)
        integer, intent(inout) :: vertex_element(:), triangles(:, :), parents(:)
        integer, intent(inout) :: triangle_count, node
        integer, intent(out) :: ierr

        integer :: element, k, neighbor, side

        do element = 1, data%n_elements
            do side = 1, 4
                if (.not. selected_side(side, element)) cycle
                neighbor = data%neighbours(element, side)
                if (neighbor > 0 .and. neighbor < element) cycle
                do k = 1, subdivisions - 1
                    if (k == subdivisions/2) cycle
                    node = node + 1
                    call store_side_node(data, element, side, k, subdivisions, &
                        plane_rz(:, node), psi(node), vertex_st(:, node), ierr)
                    if (ierr /= 0) return
                    vertex_element(node) = element
                    call insert_node_in_plane(plane_rz, node, triangles, parents, &
                        triangle_count, ierr)
                    if (ierr /= 0) return
                end do
            end do
        end do
        ierr = 0
    end subroutine insert_selected_nodes

    subroutine store_side_node(data, element, side, k, subdivisions, rz, &
            value, st, ierr)
        type(jorek_restart_t), intent(in) :: data
        integer, intent(in) :: element, side, k, subdivisions
        real(dp), intent(out) :: rz(2), value, st(2)
        integer, intent(out) :: ierr

        real(dp) :: derivatives(3), fraction, rz_st(2, 2)

        fraction = real(k, dp)/real(subdivisions, dp)
        call side_parameters(side, fraction, st)
        call evaluate_jorek_geometry(data, element, st(1), st(2), rz, &
            rz_st, ierr)
        if (ierr /= 0) return
        call evaluate_jorek_variable(data, element, 1, st(1), st(2), 0.0_dp, &
            value, derivatives, ierr)
    end subroutine store_side_node

    subroutine insert_node_in_plane(plane_rz, node, triangles, parents, &
            triangle_count, ierr)
        real(dp), intent(in) :: plane_rz(:, :)
        integer, intent(in) :: node
        integer, intent(inout) :: triangles(:, :), parents(:), triangle_count
        integer, intent(out) :: ierr

        real(dp), parameter :: tolerance = 1.0e-10_dp
        real(dp) :: weight(3)
        integer :: edge_matches(32), edge_zero(32), edge_count
        integer :: strict_match, trial, zero_count

        edge_count = 0
        strict_match = 0
        do trial = 1, triangle_count
            call barycentric_weights(plane_rz(:, triangles(trial, :)), &
                plane_rz(:, node), weight, ierr)
            if (ierr /= 0) return
            if (minval(weight) < -tolerance) cycle
            zero_count = count(abs(weight) <= tolerance)
            if (zero_count == 0) then
                if (strict_match /= 0) then
                    ierr = 2
                    return
                end if
                strict_match = trial
            else if (zero_count == 1) then
                edge_count = edge_count + 1
                if (edge_count > size(edge_matches)) then
                    ierr = 3
                    return
                end if
                edge_matches(edge_count) = trial
                edge_zero(edge_count) = minloc(abs(weight), dim=1)
            else
                ierr = 4
                return
            end if
        end do
        if (strict_match > 0 .and. edge_count > 0) then
            ierr = 5
            return
        else if (strict_match > 0) then
            call split_interior_triangle(strict_match, node, triangles, parents, &
                triangle_count)
        else if (edge_count > 0) then
            do trial = 1, edge_count
                call split_edge_triangle(edge_matches(trial), &
                    edge_zero(trial), node, triangles, parents, triangle_count)
            end do
        else
            ierr = 6
            return
        end if
        ierr = 0
    end subroutine insert_node_in_plane

    subroutine split_interior_triangle(triangle, node, triangles, parents, &
            triangle_count)
        integer, intent(in) :: triangle, node
        integer, intent(inout) :: triangles(:, :), parents(:), triangle_count

        integer :: old_nodes(3), parent

        old_nodes = triangles(triangle, :)
        parent = parents(triangle)
        triangles(triangle, :) = [old_nodes(1), old_nodes(2), node]
        triangle_count = triangle_count + 1
        triangles(triangle_count, :) = [old_nodes(2), old_nodes(3), node]
        parents(triangle_count) = parent
        triangle_count = triangle_count + 1
        triangles(triangle_count, :) = [old_nodes(3), old_nodes(1), node]
        parents(triangle_count) = parent
    end subroutine split_interior_triangle

    subroutine split_edge_triangle(triangle, zero, node, triangles, parents, &
            triangle_count)
        integer, intent(in) :: triangle, zero, node
        integer, intent(inout) :: triangles(:, :), parents(:), triangle_count

        integer :: edge1, edge2, old_nodes(3), parent

        old_nodes = triangles(triangle, :)
        parent = parents(triangle)
        edge1 = modulo(zero, 3) + 1
        edge2 = modulo(zero + 1, 3) + 1
        triangles(triangle, :) = [old_nodes(zero), old_nodes(edge1), node]
        triangle_count = triangle_count + 1
        triangles(triangle_count, :) = [old_nodes(zero), node, old_nodes(edge2)]
        parents(triangle_count) = parent
    end subroutine split_edge_triangle

    subroutine assign_factor_two_parents(data, element_vertex, triangles, &
            parents, ierr)
        type(jorek_restart_t), intent(in) :: data
        integer, intent(in) :: element_vertex(0:, 0:, :), triangles(:, :)
        integer, intent(out) :: parents(:), ierr

        integer :: cell_i, cell_j, element, nodes(4), split, triangle
        integer, parameter :: triangle_vertices(3, 2) = reshape( &
            [1, 2, 3, 1, 3, 4], [3, 2])

        triangle = 0
        do element = 1, data%n_elements
            do cell_j = 0, 1
                do cell_i = 0, 1
                    nodes = [element_vertex(cell_i, cell_j, element), &
                        element_vertex(cell_i + 1, cell_j, element), &
                        element_vertex(cell_i + 1, cell_j + 1, element), &
                        element_vertex(cell_i, cell_j + 1, element)]
                    do split = 1, 2
                        if (has_duplicate(nodes(triangle_vertices(:, split)))) cycle
                        triangle = triangle + 1
                        parents(triangle) = element
                    end do
                end do
            end do
        end do
        ierr = merge(0, 1, triangle == size(triangles, 1))
    end subroutine assign_factor_two_parents

    logical function has_duplicate(nodes)
        integer, intent(in) :: nodes(3)

        has_duplicate = nodes(1) == nodes(2) .or. nodes(1) == nodes(3) &
            .or. nodes(2) == nodes(3)
    end function has_duplicate

    subroutine barycentric_weights(vertices, point, weight, ierr)
        real(dp), intent(in) :: vertices(2, 3), point(2)
        real(dp), intent(out) :: weight(3)
        integer, intent(out) :: ierr

        real(dp) :: denominator

        denominator = (vertices(2, 2) - vertices(2, 3)) &
            *(vertices(1, 1) - vertices(1, 3)) &
            + (vertices(1, 3) - vertices(1, 2)) &
            *(vertices(2, 1) - vertices(2, 3))
        if (abs(denominator) <= tiny(denominator)) then
            ierr = 1
            return
        end if
        weight(1) = ((vertices(2, 2) - vertices(2, 3)) &
            *(point(1) - vertices(1, 3)) &
            + (vertices(1, 3) - vertices(1, 2)) &
            *(point(2) - vertices(2, 3)))/denominator
        weight(2) = ((vertices(2, 3) - vertices(2, 1)) &
            *(point(1) - vertices(1, 3)) &
            + (vertices(1, 1) - vertices(1, 3)) &
            *(point(2) - vertices(2, 3)))/denominator
        weight(3) = 1.0_dp - weight(1) - weight(2)
        ierr = 0
    end subroutine barycentric_weights

    integer function reciprocal_side(data, element, neighbor) result(side)
        type(jorek_restart_t), intent(in) :: data
        integer, intent(in) :: element, neighbor

        do side = 1, 4
            if (data%neighbours(element, side) == neighbor) return
        end do
        side = 0
    end function reciprocal_side

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

end module jorek_inserted_plane_mod
