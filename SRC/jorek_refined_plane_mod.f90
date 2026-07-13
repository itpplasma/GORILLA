module jorek_refined_plane_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use jorek_bezier, only: jorek_locator_t, build_jorek_locator, &
        evaluate_jorek_geometry, is_jorek_axis_target, &
        locate_jorek_element_indexed
    use jorek_field_values, only: evaluate_jorek_variable
    use jorek_refined_topology_mod, only: build_refined_element_nodes
    use jorek_restart, only: jorek_restart_t

    implicit none
    private

    public :: extract_refined_jorek_plane

contains

    subroutine extract_refined_jorek_plane(data, subdivisions, plane_rz, psi, &
            triangles, vertex_element, vertex_st, ierr, element_vertex)
        type(jorek_restart_t), intent(in) :: data
        integer, intent(in) :: subdivisions
        real(dp), allocatable, intent(out) :: plane_rz(:, :), psi(:)
        integer, allocatable, intent(out) :: triangles(:, :)
        integer, allocatable, intent(out) :: vertex_element(:)
        real(dp), allocatable, intent(out) :: vertex_st(:, :)
        integer, intent(out) :: ierr
        integer, allocatable, intent(out), optional :: element_vertex(:, :, :)

        integer, allocatable :: element_nodes(:, :, :)
        logical, allocatable :: is_corner(:), seen(:)
        real(dp) :: derivatives(3), rz(2), rz_st(2, 2), s, t, value
        integer :: element, i, j, n_nodes, node, triangle

        call validate_layout(data, subdivisions, ierr)
        if (ierr /= 0) return
        allocate(element_nodes(0:subdivisions, 0:subdivisions, &
            data%n_elements))
        call build_refined_element_nodes(data, subdivisions, element_nodes)
        n_nodes = maxval(element_nodes)
        allocate(plane_rz(2, n_nodes), psi(n_nodes), seen(n_nodes))
        allocate(vertex_element(n_nodes), vertex_st(2, n_nodes))
        allocate(is_corner(n_nodes), source=.false.)
        is_corner(1:data%n_nodes) = .true.
        plane_rz = 0.0_dp
        psi = 0.0_dp
        seen = .false.
        vertex_element = huge(1)
        vertex_st = 0.0_dp
        allocate(triangles(2*data%n_elements*subdivisions**2, 3))
        triangle = 0
        do element = 1, data%n_elements
            do j = 0, subdivisions
                t = real(j, dp)/subdivisions
                do i = 0, subdivisions
                    s = real(i, dp)/subdivisions
                    node = element_nodes(i, j, element)
                    call evaluate_jorek_geometry(data, element, s, t, rz, &
                        rz_st, ierr)
                    if (ierr /= 0) return
                    call evaluate_jorek_variable(data, element, 1, s, t, &
                        0.0_dp, value, derivatives, ierr)
                    if (ierr /= 0) return
                    call store_node(node, element, s, t, rz, value, plane_rz, &
                        psi, seen, vertex_element, vertex_st, ierr)
                    if (ierr /= 0) return
                end do
            end do
            do j = 0, subdivisions - 1
                do i = 0, subdivisions - 1
                    triangle = triangle + 1
                    triangles(triangle, :) = [element_nodes(i, j, element), &
                        element_nodes(i + 1, j, element), &
                        element_nodes(i + 1, j + 1, element)]
                    triangle = triangle + 1
                    triangles(triangle, :) = [element_nodes(i, j, element), &
                        element_nodes(i + 1, j + 1, element), &
                        element_nodes(i, j + 1, element)]
                end do
            end do
        end do
        if (.not. all(seen)) then
            ierr = 4
            return
        end if
        call collapse_axis_nodes(data, subdivisions, element_nodes, plane_rz, &
            psi, triangles, vertex_element, vertex_st, is_corner, ierr)
        if (ierr /= 0) return
        call assign_special_owners(data, plane_rz, is_corner, vertex_element, &
            vertex_st, ierr)
        if (ierr /= 0) return
        call orient_triangles(plane_rz, triangles, ierr)
        if (present(element_vertex)) then
            allocate(element_vertex, source=element_nodes)
        end if
    end subroutine extract_refined_jorek_plane

    subroutine validate_layout(data, subdivisions, ierr)
        type(jorek_restart_t), intent(in) :: data
        integer, intent(in) :: subdivisions
        integer, intent(out) :: ierr

        ierr = 0
        if (subdivisions < 2) then
            ierr = 9
        else if (data%n_nodes < 1 .or. data%n_elements < 1 &
                .or. data%n_vertex_max /= 4 .or. data%n_degrees /= 4 &
                .or. data%n_coord_tor < 1 .or. data%n_dim < 2 &
                .or. .not. allocated(data%vertex) &
                .or. any(shape(data%vertex) /= &
                    [data%n_elements, data%n_vertex_max]) &
                .or. any(data%vertex < 1) &
                .or. any(data%vertex > data%n_nodes)) then
            ierr = 1
        end if
    end subroutine validate_layout

    subroutine store_node(node, element, s, t, rz, value, plane_rz, psi, &
            seen, vertex_element, vertex_st, ierr)
        integer, intent(in) :: node, element
        real(dp), intent(in) :: s, t, rz(2), value
        real(dp), intent(inout) :: plane_rz(:, :), psi(:), vertex_st(:, :)
        logical, intent(inout) :: seen(:)
        integer, intent(inout) :: vertex_element(:)
        integer, intent(out) :: ierr

        real(dp) :: scale

        ierr = 0
        if (.not. seen(node)) then
            plane_rz(:, node) = rz
            psi(node) = value
            seen(node) = .true.
        else
            scale = max(1.0_dp, maxval(abs(rz)), &
                maxval(abs(plane_rz(:, node))))
            if (maxval(abs(plane_rz(:, node) - rz)) > 1.0e-10_dp*scale &
                    .or. abs(psi(node) - value) > 1.0e-10_dp*max(1.0_dp, &
                    abs(psi(node)), abs(value))) then
                ierr = 3
                return
            end if
        end if
        if (element < vertex_element(node)) then
            vertex_element(node) = element
            vertex_st(:, node) = [s, t]
        end if
    end subroutine store_node

    subroutine collapse_axis_nodes(data, subdivisions, element_nodes, &
            plane_rz, psi, triangles, vertex_element, vertex_st, is_corner, &
            ierr)
        type(jorek_restart_t), intent(in) :: data
        integer, intent(in) :: subdivisions
        integer, intent(inout) :: element_nodes(0:, 0:, :)
        real(dp), allocatable, intent(inout) :: plane_rz(:, :), psi(:)
        integer, allocatable, intent(inout) :: triangles(:, :), vertex_element(:)
        real(dp), allocatable, intent(inout) :: vertex_st(:, :)
        logical, allocatable, intent(inout) :: is_corner(:)
        integer, intent(out) :: ierr

        integer, allocatable :: parent(:), root_element(:)
        logical, allocatable :: root_corner(:)
        real(dp), allocatable :: root_st(:, :)
        integer :: element, i, j

        allocate(parent(size(psi)))
        parent = [(i, i=1, size(parent))]
        ierr = 0
        do element = 1, data%n_elements
            do i = 1, 3
                do j = i + 1, 4
                    call join_if_coincident(data%vertex(element, i), &
                        data%vertex(element, j), plane_rz, psi, parent, ierr)
                    if (ierr /= 0) return
                end do
            end do
            do i = 0, subdivisions - 1
                call join_if_coincident(element_nodes(i, 0, element), &
                    element_nodes(i + 1, 0, element), plane_rz, psi, parent, ierr)
                if (ierr /= 0) return
                call join_if_coincident(element_nodes(subdivisions, i, element), &
                    element_nodes(subdivisions, i + 1, element), plane_rz, &
                    psi, parent, ierr)
                if (ierr /= 0) return
                call join_if_coincident(element_nodes(subdivisions - i, &
                    subdivisions, element), element_nodes(subdivisions - i - 1, &
                    subdivisions, element), plane_rz, psi, parent, ierr)
                if (ierr /= 0) return
                call join_if_coincident(element_nodes(0, subdivisions - i, &
                    element), element_nodes(0, subdivisions - i - 1, element), &
                    plane_rz, psi, parent, ierr)
                if (ierr /= 0) return
            end do
        end do
        do i = 1, size(parent)
            parent(i) = node_root(parent, i)
        end do
        allocate(root_element(size(parent)), source=huge(1))
        allocate(root_st(2, size(parent)), source=0.0_dp)
        allocate(root_corner(size(parent)), source=.false.)
        do i = 1, size(parent)
            root_corner(parent(i)) = root_corner(parent(i)) .or. is_corner(i)
            if (vertex_element(i) < root_element(parent(i))) then
                root_element(parent(i)) = vertex_element(i)
                root_st(:, parent(i)) = vertex_st(:, i)
            end if
        end do
        call compact_refined_plane(parent, subdivisions, element_nodes, &
            plane_rz, psi, triangles, vertex_element, vertex_st, is_corner, &
            root_element, root_st, root_corner)
    end subroutine collapse_axis_nodes

    subroutine compact_refined_plane(parent, subdivisions, element_nodes, &
            plane_rz, psi, triangles, vertex_element, vertex_st, is_corner, &
            root_element, root_st, root_corner)
        integer, intent(in) :: parent(:), subdivisions, root_element(:)
        integer, intent(inout) :: element_nodes(0:, 0:, :)
        integer, allocatable, intent(inout) :: triangles(:, :), vertex_element(:)
        real(dp), allocatable, intent(inout) :: plane_rz(:, :), psi(:)
        real(dp), allocatable, intent(inout) :: vertex_st(:, :)
        logical, allocatable, intent(inout) :: is_corner(:)
        real(dp), intent(in) :: root_st(:, :)
        logical, intent(in) :: root_corner(:)

        integer, allocatable :: compact(:), reduced(:, :), reduced_element(:)
        logical, allocatable :: reduced_corner(:), used(:), valid(:)
        real(dp), allocatable :: reduced_rz(:, :), reduced_psi(:), reduced_st(:, :)
        integer :: count_nodes, count_triangles, element, i, j

        do i = 1, size(triangles, 1)
            triangles(i, :) = parent(triangles(i, :))
        end do
        allocate(valid(size(triangles, 1)))
        valid = triangles(:, 1) /= triangles(:, 2) &
            .and. triangles(:, 1) /= triangles(:, 3) &
            .and. triangles(:, 2) /= triangles(:, 3)
        count_triangles = count(valid)
        allocate(used(size(parent)), source=.false.)
        do i = 1, size(triangles, 1)
            if (.not. valid(i)) cycle
            used(triangles(i, :)) = .true.
        end do
        allocate(compact(size(parent)), source=0)
        count_nodes = 0
        do i = 1, size(parent)
            if (parent(i) /= i .or. .not. used(i)) cycle
            count_nodes = count_nodes + 1
            compact(i) = count_nodes
        end do
        allocate(reduced_rz(2, count_nodes), reduced_psi(count_nodes))
        allocate(reduced_element(count_nodes), reduced_st(2, count_nodes))
        allocate(reduced_corner(count_nodes))
        do i = 1, size(parent)
            if (compact(i) == 0) cycle
            reduced_rz(:, compact(i)) = plane_rz(:, i)
            reduced_psi(compact(i)) = psi(i)
            reduced_element(compact(i)) = root_element(i)
            reduced_st(:, compact(i)) = root_st(:, i)
            reduced_corner(compact(i)) = root_corner(i)
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
        call move_alloc(reduced_element, vertex_element)
        call move_alloc(reduced_st, vertex_st)
        call move_alloc(reduced_corner, is_corner)
        call move_alloc(reduced, triangles)
        do element = 1, size(element_nodes, 3)
            do j = 0, subdivisions
                do i = 0, subdivisions
                    element_nodes(i, j, element) = &
                        compact(parent(element_nodes(i, j, element)))
                end do
            end do
        end do
    end subroutine compact_refined_plane

    subroutine join_if_coincident(left, right, plane_rz, psi, parent, ierr)
        integer, intent(in) :: left, right
        real(dp), intent(in) :: plane_rz(:, :), psi(:)
        integer, intent(inout) :: parent(:)
        integer, intent(out) :: ierr

        ierr = 0
        if (.not. same_point(plane_rz(:, left), plane_rz(:, right))) return
        if (.not. same_value(psi(left), psi(right))) then
            ierr = 8
            return
        end if
        call join_nodes(parent, left, right)
    end subroutine join_if_coincident

    subroutine assign_special_owners(data, plane_rz, is_corner, &
            vertex_element, vertex_st, ierr)
        type(jorek_restart_t), intent(in) :: data
        real(dp), intent(in) :: plane_rz(:, :)
        logical, intent(in) :: is_corner(:)
        integer, intent(inout) :: vertex_element(:)
        real(dp), intent(inout) :: vertex_st(:, :)
        integer, intent(out) :: ierr

        type(jorek_locator_t) :: locator
        integer :: node

        call build_jorek_locator(data, locator, ierr)
        if (ierr /= 0) return
        do node = 1, size(vertex_element)
            if (is_jorek_axis_target(locator, plane_rz(:, node))) then
                vertex_element(node) = 0
                vertex_st(:, node) = 0.0_dp
            else if (is_corner(node)) then
                call locate_jorek_element_indexed(data, locator, &
                    plane_rz(:, node), vertex_element(node), &
                    vertex_st(1, node), vertex_st(2, node), ierr)
                if (ierr /= 0) return
            end if
        end do
    end subroutine assign_special_owners

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

        same_value = abs(left - right) <= 1.0e-9_dp*max(1.0_dp, &
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

end module jorek_refined_plane_mod
