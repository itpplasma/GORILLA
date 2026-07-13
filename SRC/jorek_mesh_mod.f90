module jorek_mesh_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use circular_mesh_SOLEDGE3X_EIRENE, only: calc_mesh_SOLEDGE3X_EIRENE
    use jorek_bezier, only: jorek_locator_t, build_jorek_locator, &
        evaluate_jorek_geometry, is_jorek_axis_target, &
        locate_jorek_element_indexed
    use jorek_field_values, only: evaluate_jorek_variable
    use jorek_restart, only: jorek_restart_t

    implicit none
    private

    public :: build_jorek_mesh_arrays

contains

    subroutine build_jorek_mesh_arrays(data, n_slices, length_scale, &
            points_rphiz, verts, neighbours, neighbour_faces, perbou_phi, &
            ierr, vertex_element, vertex_st)
        type(jorek_restart_t), intent(in) :: data
        integer, intent(in) :: n_slices
        real(dp), intent(in) :: length_scale
        real(dp), allocatable, intent(out) :: points_rphiz(:, :)
        integer, allocatable, intent(out) :: verts(:, :), neighbours(:, :)
        integer, allocatable, intent(out) :: neighbour_faces(:, :), perbou_phi(:, :)
        integer, intent(out) :: ierr
        integer, allocatable, intent(out), optional :: vertex_element(:)
        real(dp), allocatable, intent(out), optional :: vertex_st(:, :)

        real(dp), allocatable :: plane_rz(:, :), plane_st(:, :), psi(:)
        integer, allocatable :: plane_element(:), triangles(:, :)
        integer :: n_tetras

        call extract_jorek_plane(data, plane_rz, psi, triangles, &
            plane_element, plane_st, ierr)
        if (ierr /= 0) return
        if (n_slices < 3 .or. length_scale <= 0.0_dp .or. data%n_period < 1) then
            ierr = 5
            return
        end if
        call extrude_jorek_plane(plane_rz*length_scale, n_slices, &
            data%n_period, points_rphiz)
        if (present(vertex_element)) &
            call extrude_integer_plane(plane_element, n_slices, vertex_element)
        if (present(vertex_st)) call extrude_real_plane(plane_st, n_slices, vertex_st)
        call calc_mesh_SOLEDGE3X_EIRENE(n_slices, points_rphiz, size(plane_rz, 2), &
            n_tetras, verts, neighbours, neighbour_faces, perbou_phi, &
            triangles, psi)
        call validate_face_orientations(points_rphiz, verts, neighbours, &
            neighbour_faces, n_slices, data%n_period, ierr)
    end subroutine build_jorek_mesh_arrays

    subroutine extract_jorek_plane(data, plane_rz, psi, triangles, &
            vertex_element, vertex_st, ierr)
        type(jorek_restart_t), intent(in) :: data
        real(dp), allocatable, intent(out) :: plane_rz(:, :), psi(:)
        integer, allocatable, intent(out) :: triangles(:, :)
        integer, allocatable, intent(out) :: vertex_element(:)
        real(dp), allocatable, intent(out) :: vertex_st(:, :)
        integer, intent(out) :: ierr

        real(dp), parameter :: corner_s(4) = [0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp]
        real(dp), parameter :: corner_t(4) = [0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp]
        logical, allocatable :: seen(:)
        real(dp) :: derivatives(3), rz(2), rz_st(2, 2), value
        integer :: element, node, vertex

        ierr = 0
        if (data%n_nodes < 1 .or. data%n_elements < 1 &
                .or. data%n_vertex_max /= 4 .or. data%n_degrees /= 4 &
                .or. data%n_coord_tor < 1 .or. data%n_dim < 2) then
            ierr = 1
            return
        end if
        allocate(plane_rz(2, data%n_nodes), psi(data%n_nodes), seen(data%n_nodes))
        allocate(triangles(2*data%n_elements, 3))
        plane_rz = 0.0_dp
        psi = 0.0_dp
        seen = .false.
        do element = 1, data%n_elements
            call set_element_triangles(data, element, triangles, ierr)
            if (ierr /= 0) return
            do vertex = 1, 4
                node = data%vertex(element, vertex)
                call evaluate_jorek_geometry(data, element, corner_s(vertex), &
                    corner_t(vertex), rz, rz_st, ierr)
                if (ierr /= 0) return
                call evaluate_jorek_variable(data, element, 1, corner_s(vertex), &
                    corner_t(vertex), 0.0_dp, value, derivatives, ierr)
                if (ierr /= 0) return
                call store_node(node, rz, value, plane_rz, psi, seen, ierr)
                if (ierr /= 0) return
            end do
        end do
        if (.not. all(seen)) then
            ierr = 4
            return
        end if
        call collapse_coincident_nodes(data, plane_rz, psi, triangles, ierr)
        if (ierr /= 0) return
        call assign_vertex_owners(data, plane_rz, vertex_element, vertex_st, ierr)
        if (ierr /= 0) return
        call orient_triangles(plane_rz, triangles, ierr)
    end subroutine extract_jorek_plane

    subroutine assign_vertex_owners(data, plane_rz, vertex_element, vertex_st, ierr)
        type(jorek_restart_t), intent(in) :: data
        real(dp), intent(in) :: plane_rz(:, :)
        integer, allocatable, intent(out) :: vertex_element(:)
        real(dp), allocatable, intent(out) :: vertex_st(:, :)
        integer, intent(out) :: ierr

        type(jorek_locator_t) :: locator
        integer :: node

        call build_jorek_locator(data, locator, ierr)
        if (ierr /= 0) return
        allocate(vertex_element(size(plane_rz, 2)))
        allocate(vertex_st(2, size(plane_rz, 2)))
        do node = 1, size(vertex_element)
            if (is_jorek_axis_target(locator, plane_rz(:, node))) then
                vertex_element(node) = 0
                vertex_st(:, node) = 0.0_dp
                cycle
            end if
            call locate_jorek_element_indexed(data, locator, plane_rz(:, node), &
                vertex_element(node), vertex_st(1, node), vertex_st(2, node), ierr)
            if (ierr /= 0) return
        end do
    end subroutine assign_vertex_owners

    subroutine set_element_triangles(data, element, triangles, ierr)
        type(jorek_restart_t), intent(in) :: data
        integer, intent(in) :: element
        integer, intent(inout) :: triangles(:, :)
        integer, intent(out) :: ierr

        integer :: nodes(4)

        nodes = data%vertex(element, :)
        ierr = 0
        if (any(nodes < 1) .or. any(nodes > data%n_nodes)) then
            ierr = 2
            return
        end if
        triangles(2*element - 1, :) = nodes([1, 2, 3])
        triangles(2*element, :) = nodes([1, 3, 4])
    end subroutine set_element_triangles

    subroutine collapse_coincident_nodes(data, plane_rz, psi, triangles, ierr)
        type(jorek_restart_t), intent(in) :: data
        real(dp), allocatable, intent(inout) :: plane_rz(:, :), psi(:)
        integer, allocatable, intent(inout) :: triangles(:, :)
        integer, intent(out) :: ierr

        logical, allocatable :: valid(:)
        integer, allocatable :: compact(:), parent(:), reduced(:, :)
        real(dp), allocatable :: reduced_rz(:, :), reduced_psi(:)
        integer :: column, count_nodes, count_triangles, element, i, j, left, right

        parent = [(i, i=1, data%n_nodes)]
        ierr = 0
        do element = 1, data%n_elements
            do i = 1, 3
                do j = i + 1, 4
                    left = data%vertex(element, i)
                    right = data%vertex(element, j)
                    if (same_point(plane_rz(:, left), plane_rz(:, right))) then
                        if (.not. same_value(psi(left), psi(right))) then
                            ierr = 8
                            return
                        end if
                        call join_nodes(parent, left, right)
                    end if
                end do
            end do
        end do
        do i = 1, size(parent)
            parent(i) = node_root(parent, i)
        end do
        do i = 1, size(triangles, 1)
            do column = 1, 3
                triangles(i, column) = parent(triangles(i, column))
            end do
        end do
        allocate(valid(size(triangles, 1)))
        do i = 1, size(triangles, 1)
            valid(i) = all(triangles(i, 1) /= triangles(i, 2:3)) &
                .and. triangles(i, 2) /= triangles(i, 3)
        end do
        count_triangles = count(valid)
        allocate(compact(size(parent)))
        compact = 0
        count_nodes = 0
        do i = 1, size(parent)
            if (parent(i) /= i .or. .not. any(triangles == i)) cycle
            count_nodes = count_nodes + 1
            compact(i) = count_nodes
        end do
        allocate(reduced_rz(2, count_nodes), reduced_psi(count_nodes))
        do i = 1, size(parent)
            if (compact(i) == 0) cycle
            reduced_rz(:, compact(i)) = plane_rz(:, i)
            reduced_psi(compact(i)) = psi(i)
        end do
        allocate(reduced(count_triangles, 3))
        j = 0
        do i = 1, size(triangles, 1)
            if (.not. valid(i)) cycle
            j = j + 1
            reduced(j, :) = compact(triangles(i, :))
        end do
        call move_alloc(reduced_rz, plane_rz)
        call move_alloc(reduced_psi, psi)
        call move_alloc(reduced, triangles)
    end subroutine collapse_coincident_nodes

    subroutine join_nodes(parent, left, right)
        integer, intent(inout) :: parent(:)
        integer, intent(in) :: left, right

        integer :: left_root, right_root

        left_root = node_root(parent, left)
        right_root = node_root(parent, right)
        if (left_root == right_root) return
        parent(max(left_root, right_root)) = min(left_root, right_root)
    end subroutine join_nodes

    integer function node_root(parent, node)
        integer, intent(in) :: parent(:), node

        node_root = node
        do while (parent(node_root) /= node_root)
            node_root = parent(node_root)
        end do
    end function node_root

    logical function same_point(left, right)
        real(dp), intent(in) :: left(2), right(2)

        real(dp) :: scale

        scale = max(1.0_dp, maxval(abs(left)), maxval(abs(right)))
        same_point = maxval(abs(left - right)) <= 1.0e-10_dp*scale
    end function same_point

    logical function same_value(left, right)
        real(dp), intent(in) :: left, right

        same_value = abs(left - right) <= 1.0e-10_dp*max(1.0_dp, &
            abs(left), abs(right))
    end function same_value

    subroutine orient_triangles(plane_rz, triangles, ierr)
        real(dp), intent(in) :: plane_rz(:, :)
        integer, intent(inout) :: triangles(:, :)
        integer, intent(out) :: ierr

        real(dp) :: area, scale
        integer :: index, swap

        ierr = 0
        do index = 1, size(triangles, 1)
            area = (plane_rz(1, triangles(index, 2)) &
                    - plane_rz(1, triangles(index, 1))) &
                *(plane_rz(2, triangles(index, 3)) &
                    - plane_rz(2, triangles(index, 1))) &
                - (plane_rz(2, triangles(index, 2)) &
                    - plane_rz(2, triangles(index, 1))) &
                *(plane_rz(1, triangles(index, 3)) &
                    - plane_rz(1, triangles(index, 1)))
            scale = max(1.0_dp, maxval(abs(plane_rz(:, triangles(index, :)))))
            if (abs(area) <= 64.0_dp*epsilon(1.0_dp)*scale**2) then
                ierr = 6
                return
            end if
            if (area < 0.0_dp) then
                swap = triangles(index, 2)
                triangles(index, 2) = triangles(index, 3)
                triangles(index, 3) = swap
            end if
        end do
    end subroutine orient_triangles

    subroutine store_node(node, rz, value, plane_rz, psi, seen, ierr)
        integer, intent(in) :: node
        real(dp), intent(in) :: rz(2), value
        real(dp), intent(inout) :: plane_rz(:, :), psi(:)
        logical, intent(inout) :: seen(:)
        integer, intent(out) :: ierr

        real(dp) :: scale

        ierr = 0
        if (.not. seen(node)) then
            plane_rz(:, node) = rz
            psi(node) = value
            seen(node) = .true.
            return
        end if
        scale = max(1.0_dp, maxval(abs(rz)), maxval(abs(plane_rz(:, node))))
        if (maxval(abs(plane_rz(:, node) - rz)) > 1.0e-10_dp*scale &
                .or. abs(psi(node) - value) > 1.0e-10_dp*max(1.0_dp, &
                abs(psi(node)), abs(value))) ierr = 3
    end subroutine store_node

    subroutine extrude_jorek_plane(plane_rz, n_slices, n_period, points_rphiz)
        real(dp), intent(in) :: plane_rz(:, :)
        integer, intent(in) :: n_slices, n_period
        real(dp), allocatable, intent(out) :: points_rphiz(:, :)

        real(dp), parameter :: pi = acos(-1.0_dp)
        integer :: first, n_nodes, slice

        n_nodes = size(plane_rz, 2)
        allocate(points_rphiz(3, n_nodes*n_slices))
        do slice = 0, n_slices - 1
            first = slice*n_nodes + 1
            points_rphiz(1, first:first + n_nodes - 1) = plane_rz(1, :)
            points_rphiz(2, first:first + n_nodes - 1) = &
                2.0_dp*pi*slice/(n_period*n_slices)
            points_rphiz(3, first:first + n_nodes - 1) = plane_rz(2, :)
        end do
    end subroutine extrude_jorek_plane

    subroutine extrude_integer_plane(plane_values, n_slices, values)
        integer, intent(in) :: plane_values(:), n_slices
        integer, allocatable, intent(out) :: values(:)

        integer :: first, n_nodes, slice

        n_nodes = size(plane_values)
        allocate(values(n_nodes*n_slices))
        do slice = 0, n_slices - 1
            first = slice*n_nodes + 1
            values(first:first + n_nodes - 1) = plane_values
        end do
    end subroutine extrude_integer_plane

    subroutine extrude_real_plane(plane_values, n_slices, values)
        real(dp), intent(in) :: plane_values(:, :)
        integer, intent(in) :: n_slices
        real(dp), allocatable, intent(out) :: values(:, :)

        integer :: first, n_nodes, slice

        n_nodes = size(plane_values, 2)
        allocate(values(size(plane_values, 1), n_nodes*n_slices))
        do slice = 0, n_slices - 1
            first = slice*n_nodes + 1
            values(:, first:first + n_nodes - 1) = plane_values
        end do
    end subroutine extrude_real_plane

    subroutine validate_face_orientations(points, verts, neighbours, faces, &
            n_slices, n_period, ierr)
        real(dp), intent(in) :: points(:, :)
        integer, intent(in) :: verts(:, :), neighbours(:, :), faces(:, :)
        integer, intent(in) :: n_slices, n_period
        integer, intent(out) :: ierr

        real(dp) :: left(3), right(3)
        integer :: face, neighbour, tetra

        ierr = 0
        do tetra = 1, size(verts, 2)
            do face = 1, 4
                neighbour = neighbours(face, tetra)
                if (neighbour == -1) cycle
                call inward_face_normal(points, verts, tetra, face, n_slices, &
                    n_period, left)
                call inward_face_normal(points, verts, neighbour, &
                    faces(face, tetra), n_slices, n_period, right)
                if (dot_product(left, right) >= 0.0_dp) then
                    ierr = 7
                    return
                end if
            end do
        end do
    end subroutine validate_face_orientations

    subroutine inward_face_normal(points, verts, tetra, face, n_slices, &
            n_period, normal)
        real(dp), intent(in) :: points(:, :)
        integer, intent(in) :: verts(:, :), tetra, face, n_slices, n_period
        real(dp), intent(out) :: normal(3)

        real(dp), parameter :: pi = acos(-1.0_dp)
        real(dp) :: coordinates(3, 4), edge_a(3), edge_b(3), opposite(3)
        integer :: index, selected(3), tetras_per_slice

        coordinates = points(:, verts(:, tetra))
        tetras_per_slice = size(verts, 2)/n_slices
        if (tetra > size(verts, 2) - tetras_per_slice) then
            where (coordinates(2, :) == 0.0_dp)
                coordinates(2, :) = 2.0_dp*pi/n_period
            end where
        end if
        selected = pack([(index, index=1, 4)], [(index /= face, index=1, 4)])
        edge_a = coordinates(:, selected(1)) - coordinates(:, selected(3))
        edge_b = coordinates(:, selected(2)) - coordinates(:, selected(3))
        normal = [edge_a(2)*edge_b(3) - edge_a(3)*edge_b(2), &
            edge_a(3)*edge_b(1) - edge_a(1)*edge_b(3), &
            edge_a(1)*edge_b(2) - edge_a(2)*edge_b(1)]
        opposite = coordinates(:, face) - coordinates(:, selected(3))
        if (dot_product(normal, opposite) < 0.0_dp) normal = -normal
    end subroutine inward_face_normal

end module jorek_mesh_mod
