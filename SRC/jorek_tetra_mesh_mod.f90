module jorek_tetra_mesh_mod
    use, intrinsic :: iso_fortran_env, only: int64

    implicit none
    private

    public :: extrude_jorek_triangles

contains

    subroutine extrude_jorek_triangles(n_nodes, n_slices, triangles, verts, &
            neighbours, neighbour_faces, perbou_phi, ierr)
        integer, intent(in) :: n_nodes, n_slices
        integer, intent(in) :: triangles(:, :)
        integer, allocatable, intent(out) :: verts(:, :), neighbours(:, :)
        integer, allocatable, intent(out) :: neighbour_faces(:, :), perbou_phi(:, :)
        integer, intent(out) :: ierr

        integer, allocatable :: edge_left(:), edge_right(:), edge_triangle(:)
        integer :: edge, left, prism, right, slot, table_size

        ierr = 0
        if (n_nodes < 3 .or. n_slices < 3 .or. size(triangles, 2) /= 3 &
                .or. size(triangles, 1) < 1 &
                .or. any(triangles < 1) .or. any(triangles > n_nodes)) then
            ierr = 1
            return
        end if
        if (any(triangles(:, 1) == triangles(:, 2)) &
                .or. any(triangles(:, 1) == triangles(:, 3)) &
                .or. any(triangles(:, 2) == triangles(:, 3))) then
            ierr = 1
            return
        end if

        call allocate_mesh(n_nodes, n_slices, triangles, verts, neighbours, &
            neighbour_faces, perbou_phi)
        table_size = edge_table_size(size(triangles, 1))
        allocate(edge_left(table_size), edge_right(table_size), &
            edge_triangle(table_size))
        edge_left = 0
        edge_right = 0
        edge_triangle = 0
        do prism = 1, size(triangles, 1)
            do edge = 1, 3
                left = min(triangles(prism, edge), &
                    triangles(prism, modulo(edge, 3) + 1))
                right = max(triangles(prism, edge), &
                    triangles(prism, modulo(edge, 3) + 1))
                slot = edge_slot(left, right, table_size, edge_left, edge_right)
                if (slot == 0) then
                    ierr = 2
                    return
                end if
                if (edge_triangle(slot) == 0) then
                    edge_left(slot) = left
                    edge_right(slot) = right
                    edge_triangle(slot) = prism
                else if (edge_triangle(slot) > 0) then
                    call connect_prism_edges(edge_triangle(slot), prism, left, &
                        right, n_slices, triangles, neighbours, neighbour_faces)
                    edge_triangle(slot) = -edge_triangle(slot)
                else
                    ierr = 3
                    return
                end if
            end do
        end do
    end subroutine extrude_jorek_triangles

    subroutine allocate_mesh(n_nodes, n_slices, triangles, verts, neighbours, &
            neighbour_faces, perbou_phi)
        integer, intent(in) :: n_nodes, n_slices, triangles(:, :)
        integer, allocatable, intent(out) :: verts(:, :), neighbours(:, :)
        integer, allocatable, intent(out) :: neighbour_faces(:, :), perbou_phi(:, :)

        integer :: bottom(3), first, n_prisms, n_tetras, prism, slice, top(3)

        n_prisms = size(triangles, 1)
        n_tetras = 3*n_prisms*n_slices
        allocate(verts(4, n_tetras), neighbours(4, n_tetras), &
            neighbour_faces(4, n_tetras), perbou_phi(4, n_tetras))
        neighbours = -1
        neighbour_faces = -1
        perbou_phi = 0
        do slice = 0, n_slices - 1
            do prism = 1, n_prisms
                bottom = slice*n_nodes + sorted_nodes(triangles(prism, :))
                if (slice == n_slices - 1) then
                    top = sorted_nodes(triangles(prism, :))
                else
                    top = bottom + n_nodes
                end if
                first = 3*(slice*n_prisms + prism - 1) + 1
                verts(:, first) = [bottom, top(3)]
                verts(:, first + 1) = [bottom(1:2), top(2:3)]
                verts(:, first + 2) = [bottom(1), top]
                call connect_faces(first, 3, first + 1, 3, &
                    neighbours, neighbour_faces)
                call connect_faces(first + 1, 2, first + 2, 2, &
                    neighbours, neighbour_faces)
                call connect_toroidal_faces(prism, slice, n_prisms, n_slices, &
                    first, neighbours, neighbour_faces, perbou_phi)
            end do
        end do
    end subroutine allocate_mesh

    subroutine connect_toroidal_faces(prism, slice, n_prisms, n_slices, &
            first, neighbours, neighbour_faces, perbou_phi)
        integer, intent(in) :: prism, slice, n_prisms, n_slices, first
        integer, intent(inout) :: neighbours(:, :), neighbour_faces(:, :)
        integer, intent(inout) :: perbou_phi(:, :)

        integer :: next_first, previous_first

        previous_first = 3*(modulo(slice - 1, n_slices)*n_prisms + prism - 1) + 1
        next_first = 3*(modulo(slice + 1, n_slices)*n_prisms + prism - 1) + 1
        neighbours(4, first) = previous_first + 2
        neighbour_faces(4, first) = 1
        neighbours(1, first + 2) = next_first
        neighbour_faces(1, first + 2) = 4
        if (slice == 0) perbou_phi(4, first) = -1
        if (slice == n_slices - 1) perbou_phi(1, first + 2) = 1
    end subroutine connect_toroidal_faces

    subroutine connect_prism_edges(left_prism, right_prism, left_node, &
            right_node, n_slices, triangles, neighbours, neighbour_faces)
        integer, intent(in) :: left_prism, right_prism, left_node, right_node
        integer, intent(in) :: n_slices, triangles(:, :)
        integer, intent(inout) :: neighbours(:, :), neighbour_faces(:, :)

        integer :: left_faces(2), left_tetras(2), right_faces(2), right_tetras(2)
        integer :: slice

        do slice = 0, n_slices - 1
            call prism_edge_faces(left_prism, left_node, right_node, slice, &
                triangles, left_tetras, left_faces)
            call prism_edge_faces(right_prism, left_node, right_node, slice, &
                triangles, right_tetras, right_faces)
            call connect_faces(left_tetras(1), left_faces(1), right_tetras(1), &
                right_faces(1), neighbours, neighbour_faces)
            call connect_faces(left_tetras(2), left_faces(2), right_tetras(2), &
                right_faces(2), neighbours, neighbour_faces)
        end do
    end subroutine connect_prism_edges

    subroutine prism_edge_faces(prism, left_node, right_node, slice, triangles, &
            tetras, faces)
        integer, intent(in) :: prism, left_node, right_node, slice, triangles(:, :)
        integer, intent(out) :: tetras(2), faces(2)

        integer :: first, nodes(3)

        nodes = sorted_nodes(triangles(prism, :))
        first = 3*(slice*size(triangles, 1) + prism - 1) + 1
        if (left_node == nodes(1) .and. right_node == nodes(2)) then
            tetras = [first + 1, first + 2]
            faces = [4, 4]
        else if (left_node == nodes(2) .and. right_node == nodes(3)) then
            tetras = [first, first + 1]
            faces = [1, 1]
        else
            tetras = [first, first + 2]
            faces = [2, 3]
        end if
    end subroutine prism_edge_faces

    subroutine connect_faces(left_tetra, left_face, right_tetra, right_face, &
            neighbours, neighbour_faces)
        integer, intent(in) :: left_tetra, left_face, right_tetra, right_face
        integer, intent(inout) :: neighbours(:, :), neighbour_faces(:, :)

        neighbours(left_face, left_tetra) = right_tetra
        neighbour_faces(left_face, left_tetra) = right_face
        neighbours(right_face, right_tetra) = left_tetra
        neighbour_faces(right_face, right_tetra) = left_face
    end subroutine connect_faces

    integer function edge_table_size(n_triangles)
        integer, intent(in) :: n_triangles

        edge_table_size = 1
        do while (edge_table_size < 4*n_triangles)
            edge_table_size = 2*edge_table_size
        end do
    end function edge_table_size

    integer function edge_slot(left, right, table_size, edge_left, edge_right)
        integer, intent(in) :: left, right, table_size, edge_left(:), edge_right(:)

        integer :: probes
        integer(int64) :: key

        key = 104729_int64*int(left, int64) + 130363_int64*int(right, int64)
        edge_slot = 1 + int(modulo(key, int(table_size, int64)))
        do probes = 1, table_size
            if (edge_left(edge_slot) == 0 &
                    .or. (edge_left(edge_slot) == left &
                    .and. edge_right(edge_slot) == right)) return
            edge_slot = modulo(edge_slot, table_size) + 1
        end do
        edge_slot = 0
    end function edge_slot

    pure function sorted_nodes(nodes) result(sorted)
        integer, intent(in) :: nodes(3)
        integer :: sorted(3), swap

        sorted = nodes
        if (sorted(1) > sorted(2)) then
            swap = sorted(1)
            sorted(1) = sorted(2)
            sorted(2) = swap
        end if
        if (sorted(2) > sorted(3)) then
            swap = sorted(2)
            sorted(2) = sorted(3)
            sorted(3) = swap
        end if
        if (sorted(1) > sorted(2)) then
            swap = sorted(1)
            sorted(1) = sorted(2)
            sorted(2) = swap
        end if
    end function sorted_nodes

end module jorek_tetra_mesh_mod
