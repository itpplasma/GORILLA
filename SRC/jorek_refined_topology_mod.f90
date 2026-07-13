module jorek_refined_topology_mod
    use, intrinsic :: iso_fortran_env, only: int64
    use jorek_restart, only: jorek_restart_t

    implicit none
    private

    public :: build_refined_element_nodes

contains

    subroutine build_refined_element_nodes(data, subdivisions, nodes)
        type(jorek_restart_t), intent(in) :: data
        integer, intent(in) :: subdivisions
        integer, intent(out) :: nodes(0:, 0:, :)

        integer, allocatable :: side_edge(:, :)

        call build_side_edges(data, side_edge)
        call assign_element_nodes(data, subdivisions, side_edge, nodes)
    end subroutine build_refined_element_nodes

    subroutine build_side_edges(data, side_edge)
        type(jorek_restart_t), intent(in) :: data
        integer, allocatable, intent(out) :: side_edge(:, :)

        integer, parameter :: side_start(4) = [1, 2, 3, 4]
        integer, parameter :: side_end(4) = [2, 3, 4, 1]
        integer(int64), allocatable :: keys(:)
        integer, allocatable :: order(:)
        integer :: element, high, index, low, n_edges, previous, side

        allocate(side_edge(4, data%n_elements))
        allocate(keys(4*data%n_elements), order(4*data%n_elements))
        index = 0
        do element = 1, data%n_elements
            do side = 1, 4
                index = index + 1
                low = min(data%vertex(element, side_start(side)), &
                    data%vertex(element, side_end(side)))
                high = max(data%vertex(element, side_start(side)), &
                    data%vertex(element, side_end(side)))
                keys(index) = int(low, int64)*int(data%n_nodes + 1, int64) &
                    + int(high, int64)
                order(index) = index
            end do
        end do
        call sort_edge_order(order, keys, 1, size(order))
        n_edges = 0
        previous = 0
        do index = 1, size(order)
            if (index == 1) then
                n_edges = n_edges + 1
            else if (keys(order(index)) /= keys(previous)) then
                n_edges = n_edges + 1
            end if
            previous = order(index)
            element = (order(index) - 1)/4 + 1
            side = mod(order(index) - 1, 4) + 1
            side_edge(side, element) = n_edges
        end do
    end subroutine build_side_edges

    recursive subroutine sort_edge_order(order, keys, left, right)
        integer, intent(inout) :: order(:)
        integer(int64), intent(in) :: keys(:)
        integer, intent(in) :: left, right

        integer :: i, j, swap
        integer(int64) :: pivot

        i = left
        j = right
        pivot = keys(order((left + right)/2))
        do
            do while (keys(order(i)) < pivot)
                i = i + 1
            end do
            do while (keys(order(j)) > pivot)
                j = j - 1
            end do
            if (i <= j) then
                swap = order(i)
                order(i) = order(j)
                order(j) = swap
                i = i + 1
                j = j - 1
            end if
            if (i > j) exit
        end do
        if (left < j) call sort_edge_order(order, keys, left, j)
        if (i < right) call sort_edge_order(order, keys, i, right)
    end subroutine sort_edge_order

    subroutine assign_element_nodes(data, subdivisions, side_edge, nodes)
        type(jorek_restart_t), intent(in) :: data
        integer, intent(in) :: subdivisions, side_edge(:, :)
        integer, intent(out) :: nodes(0:, 0:, :)

        integer, parameter :: side_start(4) = [1, 2, 3, 4]
        integer, parameter :: side_end(4) = [2, 3, 4, 1]
        integer :: canonical, edge, element, i, j, node, offset, side

        offset = data%n_nodes + maxval(side_edge)*(subdivisions - 1)
        do element = 1, data%n_elements
            nodes(0, 0, element) = data%vertex(element, 1)
            nodes(subdivisions, 0, element) = data%vertex(element, 2)
            nodes(subdivisions, subdivisions, element) = data%vertex(element, 3)
            nodes(0, subdivisions, element) = data%vertex(element, 4)
            do side = 1, 4
                edge = side_edge(side, element)
                do i = 1, subdivisions - 1
                    if (data%vertex(element, side_start(side)) &
                            < data%vertex(element, side_end(side))) then
                        canonical = i
                    else
                        canonical = subdivisions - i
                    end if
                    node = data%n_nodes + (edge - 1)*(subdivisions - 1) &
                        + canonical
                    select case (side)
                    case (1)
                        nodes(i, 0, element) = node
                    case (2)
                        nodes(subdivisions, i, element) = node
                    case (3)
                        nodes(subdivisions - i, subdivisions, element) = node
                    case (4)
                        nodes(0, subdivisions - i, element) = node
                    end select
                end do
            end do
            do j = 1, subdivisions - 1
                do i = 1, subdivisions - 1
                    nodes(i, j, element) = offset &
                        + (element - 1)*(subdivisions - 1)**2 &
                        + (j - 1)*(subdivisions - 1) + i
                end do
            end do
        end do
    end subroutine assign_element_nodes

end module jorek_refined_topology_mod
