program test_jorek_edge_sagitta
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use jorek_bezier, only: evaluate_jorek_geometry
    use jorek_refined_plane_mod, only: extract_refined_jorek_plane
    use jorek_restart, only: jorek_restart_t, load_jorek_restart

    implicit none

    type(jorek_restart_t) :: data
    real(dp), allocatable :: plane_rz(:, :), psi(:), sagitta(:), vertex_st(:, :)
    integer, allocatable :: element_vertex(:, :, :), triangles(:, :)
    integer, allocatable :: vertex_element(:)
    logical, allocatable :: aliased(:)
    character(len=1024) :: filename, output_filename
    real(dp) :: exact_midpoint(2), maximum_sagitta
    real(dp) :: rz_st(2, 2), s, t
    integer :: element, ierr, neighbor, output_unit, record, segment, side

    if (command_argument_count() /= 2) &
        error stop 'restart path and output path are required'
    call get_command_argument(1, filename)
    call get_command_argument(2, output_filename)
    call load_jorek_restart(trim(filename), data, ierr)
    if (ierr /= 0) error stop 'JOREK restart load failed'
    call extract_refined_jorek_plane(data, 2, plane_rz, psi, triangles, &
        vertex_element, vertex_st, ierr, element_vertex)
    if (ierr /= 0) error stop 'JOREK factor-2 plane extraction failed'
    allocate(sagitta(4*data%n_elements), aliased(4*data%n_elements))
    open(newunit=output_unit, file=trim(output_filename), status='replace', &
        action='write', iostat=ierr)
    if (ierr /= 0) error stop 'cannot open JOREK edge-sagitta output'
    write(output_unit, '(A)') &
        'node1,node2,neighbor_node1,neighbor_node2,element,side,segment,' // &
        'neighbor,neighbor_side,sagitta_m,aliased'

    record = 0
    do element = 1, data%n_elements
        do side = 1, 4
            neighbor = data%neighbours(element, side)
            if (neighbor <= element) cycle
            do segment = 0, 1
                record = record + 1
                call side_parameters(side, 0.25_dp + 0.5_dp*segment, s, t)
                call evaluate_jorek_geometry(data, element, s, t, &
                    exact_midpoint, rz_st, ierr)
                if (ierr /= 0) error stop 'JOREK edge midpoint evaluation failed'
                call store_segment(element, side, segment, neighbor, record, &
                    exact_midpoint)
            end do
        end do
    end do
    close(output_unit)
    if (record /= 24218) error stop 'JOREK shared-edge population changed'
    maximum_sagitta = maxval(sagitta(:record))
    if (.not. all(ieee_is_finite(sagitta(:record))) &
            .or. any(sagitta(:record) < 0.0_dp) &
            .or. maximum_sagitta > 2.0e-3_dp) &
        error stop 'JOREK shared-edge sagitta is invalid'
    if (count(aliased(:record)) /= 0) &
        error stop 'JOREK shared edge retains duplicate node identities'
    print '(A, I0)', 'shared factor-2 edge segments=', record
    print '(A, ES14.6)', 'maximum shared-edge sagitta [m]=', maximum_sagitta
    print '(A, I0)', 'aliased shared segments=', count(aliased(:record))
    print '(A)', 'PASS: conforming JOREK shared-edge geometry'

contains

    subroutine store_segment(element, side, segment, neighbor, record, midpoint)
        integer, intent(in) :: element, side, segment, neighbor, record
        real(dp), intent(in) :: midpoint(2)

        real(dp) :: chord_midpoint(2)
        integer :: endpoint_nodes(2), key(2), neighbor_nodes(2), neighbor_side
        integer :: neighbor_key(2)

        endpoint_nodes = [side_node(element, side, segment), &
            side_node(element, side, segment + 1)]
        key = [minval(endpoint_nodes), maxval(endpoint_nodes)]
        neighbor_side = reciprocal_side(neighbor, element)
        call matched_neighbor_nodes(endpoint_nodes, neighbor, neighbor_side, &
            neighbor_nodes)
        neighbor_key = [minval(neighbor_nodes), maxval(neighbor_nodes)]
        chord_midpoint = 0.5_dp*(plane_rz(:, endpoint_nodes(1)) &
            + plane_rz(:, endpoint_nodes(2)))
        sagitta(record) = sqrt(sum((midpoint - chord_midpoint)**2))
        aliased(record) = any(key /= neighbor_key)
        write(output_unit, '(*(g0,:,","))') key, neighbor_key, element, side, &
            segment, neighbor, neighbor_side, sagitta(record), &
            merge(1, 0, aliased(record))
    end subroutine store_segment

    integer function reciprocal_side(element, neighbor) result(side)
        integer, intent(in) :: element, neighbor

        do side = 1, 4
            if (data%neighbours(element, side) == neighbor) return
        end do
        error stop 'JOREK shared side is not reciprocal'
    end function reciprocal_side

    subroutine matched_neighbor_nodes(nodes, neighbor, side, matched)
        integer, intent(in) :: nodes(2), neighbor, side
        integer, intent(out) :: matched(2)

        real(dp) :: best, direct, reverse
        integer :: candidate, trial(2)

        best = huge(1.0_dp)
        matched = 0
        do candidate = 0, 1
            trial = [side_node(neighbor, side, candidate), &
                side_node(neighbor, side, candidate + 1)]
            direct = sum((plane_rz(:, nodes(1)) - plane_rz(:, trial(1)))**2) &
                + sum((plane_rz(:, nodes(2)) - plane_rz(:, trial(2)))**2)
            reverse = sum((plane_rz(:, nodes(1)) - plane_rz(:, trial(2)))**2) &
                + sum((plane_rz(:, nodes(2)) - plane_rz(:, trial(1)))**2)
            if (min(direct, reverse) < best) then
                best = min(direct, reverse)
                matched = trial
            end if
        end do
        if (sqrt(best) > 1.0e-8_dp) &
            error stop 'JOREK neighbor segment geometry does not match'
    end subroutine matched_neighbor_nodes

    integer function side_node(element, side, index) result(node)
        integer, intent(in) :: element, side, index

        select case (side)
        case (1)
            node = element_vertex(index, 0, element)
        case (2)
            node = element_vertex(2, index, element)
        case (3)
            node = element_vertex(2 - index, 2, element)
        case (4)
            node = element_vertex(0, 2 - index, element)
        end select
    end function side_node

    subroutine side_parameters(side, u, s, t)
        integer, intent(in) :: side
        real(dp), intent(in) :: u
        real(dp), intent(out) :: s, t

        select case (side)
        case (1)
            s = u
            t = 0.0_dp
        case (2)
            s = 1.0_dp
            t = u
        case (3)
            s = 1.0_dp - u
            t = 1.0_dp
        case (4)
            s = 0.0_dp
            t = 1.0_dp - u
        end select
    end subroutine side_parameters

end program test_jorek_edge_sagitta
