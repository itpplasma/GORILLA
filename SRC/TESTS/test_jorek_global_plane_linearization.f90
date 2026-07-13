program test_jorek_global_plane_linearization
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use jorek_bezier, only: jorek_locator_t, build_jorek_locator, &
        evaluate_jorek_geometry
    use jorek_boundary_plane_mod, only: extract_boundary_refined_jorek_plane, &
        extract_selected_refined_jorek_plane
    use jorek_field_backend_mod, only: evaluate_jorek_model303_gorilla, &
        evaluate_jorek_model303_gorilla_element
    use jorek_model303_field, only: evaluate_jorek_model303_b
    use jorek_refined_plane_mod, only: extract_refined_jorek_plane
    use jorek_restart, only: jorek_restart_t, load_jorek_restart
    use jorek_triangle_locator_mod, only: build_jorek_triangle_locator, &
        jorek_triangle_locator_t, locate_jorek_triangle

    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: phases(4) = [0.0_dp, 0.5_dp*pi, pi, 1.5_dp*pi]
    real(dp), parameter :: weights(3, 4) = reshape([ &
        1.0_dp/3, 1.0_dp/3, 1.0_dp/3, 0.6_dp, 0.2_dp, 0.2_dp, &
        0.2_dp, 0.6_dp, 0.2_dp, 0.2_dp, 0.2_dp, 0.6_dp], [3, 4])
    integer, parameter :: triangle_vertices(3, 2) = &
        reshape([1, 2, 3, 1, 3, 4], [3, 2])
    integer, parameter :: max_tail_records = 100
    real(dp), parameter :: corner_s(4) = [0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp]
    real(dp), parameter :: corner_t(4) = [0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp]
    character(len=11), parameter :: relation_names(4) = &
        [character(len=11) :: 'same', 'face', 'vertex_only', 'remote']
    type :: field_maximum_t
        real(dp) :: error = 0.0_dp
        real(dp) :: s = 0.0_dp, t = 0.0_dp, phi = 0.0_dp
        real(dp) :: rz(2) = 0.0_dp, weight(3) = 0.0_dp
        real(dp) :: exact(3) = 0.0_dp, interpolated(3) = 0.0_dp
        real(dp) :: vertex_rz(2, 3) = 0.0_dp, vertex_st(2, 3) = 0.0_dp
        integer :: element = 0, triangle = 0, triangle_element = 0, matches = 0
        integer :: relation = 0
        integer :: element_distance = -1
        integer :: nodes(3) = 0, owners(3) = 0, neighbors(4) = 0
        integer :: owner_neighbors(4, 3) = 0
    end type field_maximum_t
    type(jorek_restart_t) :: data
    type(jorek_locator_t) :: field_locator
    type(jorek_triangle_locator_t) :: plane_locator
    type(field_maximum_t) :: maximum(2)
    real(dp), allocatable :: plane_rz(:, :), psi(:), vertex_st(:, :)
    integer, allocatable :: element_vertex(:, :, :), triangles(:, :)
    integer, allocatable :: element_triangle_range(:, :)
    integer, allocatable :: triangle_parent(:), vertex_element(:)
    logical, allocatable :: hole_element(:), selected_side(:, :)
    real(dp) :: barycentric(3), corner_rz(2, 4), metrics(2, 2), rz(2)
    real(dp) :: recovered_relation_metrics(2, 4)
    real(dp) :: rz_st(2, 2), s, t
    character(len=1024) :: argument, filename, output_filename
    character(len=1024) :: maximum_output_filename, overlap_output_filename
    character(len=1024) :: tail_output_filename, topology_output_filename
    integer :: ambiguous, boundary_outside, element, global_outside, i, ierr, j
    integer :: interior_outside, located, refinement
    integer :: matches, neighbor_recovered
    integer :: phase, sample, source_covered, source_outside, triangle
    integer :: strict_ambiguous, unique_sample
    integer :: metric_counts(2), recovered_high_counts(4)
    integer :: recovered_relation_counts(4)
    integer :: tail_edge_count, tail_edge_nodes(2, max_tail_records)
    integer :: tail_edge_phase_counts(max_tail_records)
    integer :: tail_point_count, tail_point_element(max_tail_records)
    real(dp) :: tail_point_s(max_tail_records), tail_point_t(max_tail_records)
    integer :: output_unit, overlap_output_unit, tail_output_unit
    logical :: boundary_mode, custom_mode, source_hit, tail_mode, write_output
    logical :: write_overlap_output, write_tail_output

    if (command_argument_count() < 2 .or. command_argument_count() > 7) &
        error stop 'restart path, refinement, and optional output paths are required'
    call get_command_argument(1, filename)
    call get_command_argument(2, argument)
    boundary_mode = trim(argument) == 'boundary8'
    tail_mode = trim(argument) == 'tail4'
    custom_mode = boundary_mode .or. tail_mode
    if (boundary_mode) then
        refinement = -8
    else if (tail_mode) then
        refinement = -4
    else
        read(argument, *, iostat=ierr) refinement
        if (ierr /= 0 .or. .not. any(refinement == [2, 4, 8, 16])) &
            error stop 'refinement must be 2, 4, 8, 16, boundary8, or tail4'
    end if
    call load_jorek_restart(trim(filename), data, ierr)
    if (ierr /= 0) error stop 'JOREK restart load failed'
    call build_jorek_locator(data, field_locator, ierr)
    if (ierr /= 0) error stop 'JOREK field locator build failed'
    if (boundary_mode) then
        call extract_boundary_refined_jorek_plane(data, 8, plane_rz, psi, &
            triangles, vertex_element, vertex_st, element_triangle_range, ierr)
    else if (tail_mode) then
        call build_tail_side_mask
        call extract_selected_refined_jorek_plane(data, 4, selected_side, &
            plane_rz, psi, triangles, vertex_element, vertex_st, &
            element_triangle_range, ierr)
    else
        call extract_refined_jorek_plane(data, refinement, plane_rz, psi, &
            triangles, vertex_element, vertex_st, ierr, element_vertex)
    end if
    if (ierr /= 0) error stop 'JOREK plane refinement failed'
    allocate(triangle_parent(size(triangles, 1)))
    if (custom_mode) then
        call build_boundary_triangle_parents
    else
        call build_triangle_parents
    end if
    call build_jorek_triangle_locator(plane_rz, triangles, plane_locator, ierr)
    if (ierr /= 0) error stop 'JOREK triangle locator build failed'
    allocate(hole_element(data%n_elements), source=.false.)
    output_unit = -1
    overlap_output_unit = -1
    tail_output_unit = -1
    write_output = .false.
    write_overlap_output = .false.
    write_tail_output = .false.
    if (command_argument_count() >= 3) then
        call get_command_argument(3, output_filename)
        open(newunit=output_unit, file=trim(output_filename), status='replace', &
            action='write', iostat=ierr)
        if (ierr /= 0) error stop 'cannot open global-hole output'
        write_output = .true.
        write(output_unit, '(A)') 'element,s,t,r_m,z_m,' // &
            'neighbor_1,neighbor_2,neighbor_3,neighbor_4'
    end if
    if (command_argument_count() >= 4) then
        call get_command_argument(4, overlap_output_filename)
        open(newunit=overlap_output_unit, file=trim(overlap_output_filename), &
            status='replace', action='write', iostat=ierr)
        if (ierr /= 0) error stop 'cannot open global-overlap output'
        write_overlap_output = .true.
        write(overlap_output_unit, '(A)') &
            'element,s,t,r_m,z_m,first_min_weight,matches,' // &
            'source_hit,neighbor_1,neighbor_2,neighbor_3,neighbor_4'
    end if
    if (command_argument_count() == 7) then
        call get_command_argument(7, tail_output_filename)
        open(newunit=tail_output_unit, file=trim(tail_output_filename), &
            status='replace', action='write', iostat=ierr)
        if (ierr /= 0) error stop 'cannot open recovered-tail output'
        write_tail_output = .true.
        write(tail_output_unit, '(A)') &
            'relation,error,element,s,t,phi_rad,r_m,z_m,triangle,' // &
            'triangle_element,element_distance,matches,w1,w2,w3,' // &
            'node1,node2,node3,owner1,owner2,owner3,' // &
            'cell_i,cell_j,crossed_edge,edge_node1,edge_node2,' // &
            'signed_chord_distance_m,element_side,side_neighbor,' // &
            'exact_br_g,exact_bphi_g,exact_bz_g,' // &
            'interpolated_br_g,interpolated_bphi_g,interpolated_bz_g'
    end if
    global_outside = 0
    boundary_outside = 0
    interior_outside = 0
    ambiguous = 0
    neighbor_recovered = 0
    source_covered = 0
    source_outside = 0
    strict_ambiguous = 0
    metrics = 0.0_dp
    metric_counts = 0
    recovered_high_counts = 0
    recovered_relation_counts = 0
    recovered_relation_metrics = 0.0_dp
    tail_edge_count = 0
    tail_edge_nodes = 0
    tail_edge_phase_counts = 0
    tail_point_count = 0
    tail_point_element = 0
    tail_point_s = 0.0_dp
    tail_point_t = 0.0_dp
    unique_sample = 0

    do element = 1, data%n_elements
        call element_corners(element, corner_rz)
        if (has_coincident_corners(corner_rz)) cycle
        do j = 0, 1
            do i = 0, 1
                do triangle = 1, 2
                    do sample = 1, size(weights, 2)
                        unique_sample = unique_sample + 1
                        call reference_parameters(i, j, triangle, &
                            weights(:, sample), s, t)
                        call evaluate_jorek_geometry(data, element, s, t, rz, &
                            rz_st, ierr)
                        if (ierr /= 0) error stop 'JOREK sample geometry failed'
                        call locate_jorek_triangle(plane_rz, triangles, &
                            plane_locator, rz, located, barycentric, ierr, matches)
                        if (modulo(unique_sample, 12000) == 0) &
                            call verify_indexed_location(rz, located, &
                                barycentric, matches, ierr)
                        source_hit = source_contains(element, s, t, rz)
                        if (matches > 1) then
                            ambiguous = ambiguous + size(phases)
                            if (minval(barycentric) > 1.0e-8_dp) &
                                strict_ambiguous = strict_ambiguous + size(phases)
                            if (write_overlap_output) write(overlap_output_unit, &
                                '(I0,5(",",ES24.16E3),6(",",I0))') element, &
                                s, t, rz, minval(barycentric), matches, &
                                merge(1, 0, source_hit), &
                                data%neighbours(element, :)
                        end if
                        if (source_hit) then
                            source_covered = source_covered + size(phases)
                        else
                            source_outside = source_outside + size(phases)
                            if (ierr == 0) &
                                neighbor_recovered = neighbor_recovered + size(phases)
                        end if
                        if (ierr /= 0) then
                            hole_element(element) = .true.
                            if (write_output) write(output_unit, &
                                '(I0,4(",",ES24.16E3),4(",",I0))') element, &
                                s, t, rz, data%neighbours(element, :)
                            global_outside = global_outside + size(phases)
                            if (any(data%neighbours(element, :) <= 0)) then
                                boundary_outside = boundary_outside + size(phases)
                            else
                                interior_outside = interior_outside + size(phases)
                            end if
                            cycle
                        end if
                        do phase = 1, size(phases)
                            call compare_field(element, s, t, phases(phase), rz, &
                                located, barycentric, matches, &
                                merge(1, 2, source_hit))
                        end do
                    end do
                end do
            end do
        end do
    end do
    if (.not. boundary_mode) then
        do i = 1, 2
            maximum(i)%element_distance = neighbor_distance( &
                maximum(i)%element, maximum(i)%triangle_element)
        end do
    end if
    select case (refinement)
    case (-8)
        if (source_covered /= 420684 .or. source_outside /= 349108 &
                .or. neighbor_recovered /= 344324 .or. global_outside /= 4784 &
                .or. ambiguous /= 0 .or. strict_ambiguous /= 0 &
                .or. boundary_outside /= 4784 .or. interior_outside /= 0 &
                .or. count(hole_element) /= 148) &
            error stop 'boundary-8 global-plane fixture changed'
    case (-4)
        if (source_covered /= 420320 .or. source_outside /= 349472 &
                .or. neighbor_recovered /= 344272 .or. global_outside /= 5200 &
                .or. ambiguous /= 92 .or. strict_ambiguous /= 92 &
                .or. boundary_outside /= 5200 .or. interior_outside /= 0 &
                .or. count(hole_element) /= 148) &
            error stop 'tail-4 global-plane fixture changed'
    case (2)
        if (source_covered /= 420264 .or. source_outside /= 349528 &
                .or. neighbor_recovered /= 344328 .or. global_outside /= 5200 &
                .or. ambiguous /= 0 .or. boundary_outside /= 5200 &
                .or. interior_outside /= 0 .or. count(hole_element) /= 148) &
            error stop 'factor-2 global-plane fixture changed'
        call verify_factor2_maximum
    case (4)
        if (source_covered /= 432832 .or. source_outside /= 336960 &
                .or. neighbor_recovered /= 332200 .or. global_outside /= 4760 &
                .or. ambiguous /= 7628 .or. strict_ambiguous /= 7628 &
                .or. boundary_outside /= 4760 &
                .or. interior_outside /= 0 .or. count(hole_element) /= 148) &
            error stop 'factor-4 global-plane fixture changed'
    case (8)
        if (source_covered /= 765124 .or. source_outside /= 4668 &
                .or. neighbor_recovered /= 4372 .or. global_outside /= 296 &
                .or. ambiguous /= 334792 .or. boundary_outside /= 296 &
                .or. interior_outside /= 0 .or. count(hole_element) /= 59) &
            error stop 'factor-8 global-plane fixture changed'
    case (16)
        if (source_covered /= 769032 .or. source_outside /= 760 &
                .or. neighbor_recovered /= 756 .or. global_outside /= 4 &
                .or. ambiguous /= 338408 .or. boundary_outside /= 4 &
                .or. interior_outside /= 0 .or. count(hole_element) /= 1) &
            error stop 'factor-16 global-plane fixture changed'
    end select
    if (source_covered + source_outside /= 769792 &
            .or. source_covered + source_outside - global_outside &
                /= sum(metric_counts) &
            .or. sum(recovered_relation_counts) /= neighbor_recovered &
            .or. sum(recovered_high_counts) > neighbor_recovered &
            .or. boundary_outside + interior_outside /= global_outside &
            .or. .not. all(ieee_is_finite(metrics)) &
            .or. .not. all(ieee_is_finite(recovered_relation_metrics))) &
        error stop 'global-plane sample partition failed'
    if (refinement == 2) then
        if (sum(tail_edge_phase_counts) /= sum(recovered_high_counts) &
                .or. tail_point_count /= 18 .or. tail_edge_count /= 16) &
            error stop 'factor-2 recovered chord partition changed'
    end if
    if (boundary_mode) then
        print '(A)', 'refinement=boundary8'
    else if (tail_mode) then
        print '(A)', 'refinement=tail4'
    else
        print '(A, I0)', 'refinement=', refinement
    end if
    print '(A, I0, A, I0, A, I0)', 'source outside=', source_outside, &
        ' recovered by global plane=', neighbor_recovered, &
        ' global outside=', global_outside
    print '(A, I0)', 'multiple containing global triangles=', ambiguous
    print '(A, I0)', 'strict first-triangle overlaps=', strict_ambiguous
    print '(A, I0, A, I0)', 'global outside from boundary elements=', &
        boundary_outside, ' from interior elements=', interior_outside
    print '(A, I0)', 'elements containing global holes=', count(hole_element)
    print '(A, 2(ES14.6, 1X))', 'source-contained max/rms B error: ', &
        metrics(1, 1), sqrt(metrics(2, 1)/metric_counts(1))
    print '(A, 2(ES14.6, 1X))', 'neighbor-recovered max/rms B error: ', &
        metrics(1, 2), sqrt(metrics(2, 2)/neighbor_recovered)
    if (refinement == 2) print '(A, I0, A, I0)', &
        'above-gate unique points=', tail_point_count, &
        ' crossed chords=', tail_edge_count
    if (command_argument_count() >= 6) then
        call get_command_argument(6, topology_output_filename)
        call write_topology_output(trim(topology_output_filename))
    end if
    if (command_argument_count() >= 5) then
        call get_command_argument(5, maximum_output_filename)
        call write_maximum_output(trim(maximum_output_filename))
    end if
    print '(A)', 'PASS: global JOREK plane containment and interpolation'
    if (write_output) close(output_unit)
    if (write_overlap_output) close(overlap_output_unit)
    if (write_tail_output) close(tail_output_unit)

contains

    subroutine build_tail_side_mask
        integer, parameter :: tail_element(10) = [ &
            5061, 5062, 5063, 5174, 5175, 5175, 5986, 5986, 5987, 5988]
        integer, parameter :: tail_side(10) = [2, 2, 2, 2, 3, 2, 2, 3, 2, 2]
        integer :: index, neighbor, neighbor_side, side

        allocate(selected_side(4, data%n_elements), source=.false.)
        do index = 1, size(tail_element)
            selected_side(tail_side(index), tail_element(index)) = .true.
            neighbor = data%neighbours(tail_element(index), tail_side(index))
            if (neighbor <= 0) error stop 'tail side is not an interior interface'
            neighbor_side = 0
            do side = 1, 4
                if (data%neighbours(neighbor, side) == tail_element(index)) &
                    neighbor_side = side
            end do
            if (neighbor_side == 0) error stop 'tail side is not reciprocal'
            selected_side(neighbor_side, neighbor) = .true.
        end do
        if (count(selected_side) /= 20) &
            error stop 'tail side closure changed'
    end subroutine build_tail_side_mask

    subroutine element_corners(element, corners)
        integer, intent(in) :: element
        real(dp), intent(out) :: corners(2, 4)

        real(dp) :: rz_st(2, 2)
        integer :: corner, ierr

        do corner = 1, 4
            call evaluate_jorek_geometry(data, element, &
                corner_s(corner), corner_t(corner), corners(:, corner), &
                rz_st, ierr)
            if (ierr /= 0) error stop 'JOREK corner geometry failed'
        end do
    end subroutine element_corners

    subroutine reference_parameters(cell_i, cell_j, triangle, weight, s, t)
        integer, intent(in) :: cell_i, cell_j, triangle
        real(dp), intent(in) :: weight(3)
        real(dp), intent(out) :: s, t

        real(dp) :: local_s(4), local_t(4)
        integer :: index, vertex

        local_s = 0.5_dp*real([cell_i, cell_i + 1, cell_i + 1, cell_i], dp)
        local_t = 0.5_dp*real([cell_j, cell_j, cell_j + 1, cell_j + 1], dp)
        s = 0.0_dp
        t = 0.0_dp
        do index = 1, 3
            vertex = triangle_vertices(index, triangle)
            s = s + weight(index)*local_s(vertex)
            t = t + weight(index)*local_t(vertex)
        end do
    end subroutine reference_parameters

    logical function source_contains(element, s, t, point)
        integer, intent(in) :: element
        real(dp), intent(in) :: s, t, point(2)

        real(dp) :: corners(2, 4), weight(3)
        integer :: cell_i, cell_j, grid_i(4), grid_j(4), ierr, node, trial

        if (custom_mode) then
            source_contains = .false.
            do trial = element_triangle_range(1, element), &
                    element_triangle_range(2, element)
                call barycentric_weights(plane_rz(:, triangles(trial, :)), &
                    point, weight, ierr)
                if (ierr == 0) then
                    if (all(weight >= -1.0e-10_dp) &
                            .and. all(weight <= 1.0_dp + 1.0e-10_dp)) then
                        source_contains = .true.
                        return
                    end if
                end if
            end do
            return
        end if
        cell_i = min(refinement - 1, int(s*refinement))
        cell_j = min(refinement - 1, int(t*refinement))
        grid_i = [cell_i, cell_i + 1, cell_i + 1, cell_i]
        grid_j = [cell_j, cell_j, cell_j + 1, cell_j + 1]
        do trial = 1, 4
            node = element_vertex(grid_i(trial), grid_j(trial), element)
            corners(:, trial) = plane_rz(:, node)
        end do
        source_contains = .false.
        do trial = 1, 2
            call barycentric_weights(corners(:, triangle_vertices(:, trial)), &
                point, weight, ierr)
            if (ierr == 0) then
                if (all(weight >= -1.0e-10_dp) &
                        .and. all(weight <= 1.0_dp + 1.0e-10_dp)) then
                    source_contains = .true.
                    return
                end if
            end if
        end do
    end function source_contains

    subroutine compare_field(element, s, t, phi, rz, triangle, weight, &
            matches, group)
        integer, intent(in) :: element, triangle, matches, group
        real(dp), intent(in) :: s, t, phi, rz(2), weight(3)

        real(dp) :: a(3), b_covariant(3), b_exact_jorek(3), b_exact(3)
        real(dp) :: b_interp(3), bmod, bmod_interp, error, h_interp(3), r_cm
        integer :: corner, ierr, node, owner

        h_interp = 0.0_dp
        bmod_interp = 0.0_dp
        do corner = 1, 3
            node = triangles(triangle, corner)
            owner = vertex_element(node)
            r_cm = 100.0_dp*plane_rz(1, node)
            if (owner == 0) then
                call evaluate_jorek_model303_gorilla(data, r_cm, phi, &
                    100.0_dp*plane_rz(2, node), a, b_covariant, bmod, ierr, &
                    field_locator, 100.0_dp, 1.0e4_dp)
            else
                call evaluate_jorek_model303_gorilla_element(data, owner, &
                    vertex_st(1, node), vertex_st(2, node), phi, r_cm, a, &
                    b_covariant, bmod, ierr, 100.0_dp, 1.0e4_dp, field_locator)
            end if
            if (ierr /= 0 .or. bmod <= 0.0_dp) &
                error stop 'JOREK plane vertex field failed'
            h_interp = h_interp + weight(corner)*b_covariant/bmod
            bmod_interp = bmod_interp + weight(corner)*bmod
        end do
        r_cm = 100.0_dp*rz(1)
        b_interp = [h_interp(1)*bmod_interp, &
            h_interp(2)*bmod_interp/r_cm, h_interp(3)*bmod_interp]
        call evaluate_jorek_model303_b(data, element, s, t, phi, &
            b_exact_jorek, ierr)
        if (ierr /= 0) error stop 'JOREK exact sample field failed'
        b_exact = 1.0e4_dp*[b_exact_jorek(1), b_exact_jorek(3), b_exact_jorek(2)]
        error = sqrt(sum((b_interp - b_exact)**2)) &
            /max(sqrt(sum(b_exact**2)), tiny(1.0_dp))
        if (error > metrics(1, group)) call record_maximum(element, s, t, phi, &
            rz, triangle, weight, matches, b_exact, b_interp, error, group)
        if (group == 2) then
            call record_recovered_relation(element, triangle_parent(triangle), &
                error)
            if (error > 0.02_dp) then
                if (refinement == 2) call record_tail_chord(element, s, t, rz)
                if (write_tail_output) call write_recovered_tail(element, s, &
                    t, phi, rz, triangle, weight, matches, b_exact, b_interp, &
                    error)
            end if
        end if
        metrics(1, group) = max(metrics(1, group), error)
        metrics(2, group) = metrics(2, group) + error**2
        metric_counts(group) = metric_counts(group) + 1
    end subroutine compare_field

    subroutine record_maximum(element, s, t, phi, rz, triangle, weight, &
            matches, b_exact, b_interp, error, group)
        integer, intent(in) :: element, triangle, matches, group
        real(dp), intent(in) :: s, t, phi, rz(2), weight(3)
        real(dp), intent(in) :: b_exact(3), b_interp(3), error

        integer :: corner, node

        maximum(group)%error = error
        maximum(group)%s = s
        maximum(group)%t = t
        maximum(group)%phi = phi
        maximum(group)%rz = rz
        maximum(group)%weight = weight
        maximum(group)%exact = b_exact
        maximum(group)%interpolated = b_interp
        maximum(group)%element = element
        maximum(group)%triangle = triangle
        maximum(group)%triangle_element = triangle_parent(triangle)
        maximum(group)%relation = parent_relation(element, &
            maximum(group)%triangle_element)
        maximum(group)%matches = matches
        maximum(group)%neighbors = data%neighbours(element, :)
        do corner = 1, 3
            node = triangles(triangle, corner)
            maximum(group)%nodes(corner) = node
            maximum(group)%owners(corner) = vertex_element(node)
            maximum(group)%vertex_rz(:, corner) = plane_rz(:, node)
            maximum(group)%vertex_st(:, corner) = vertex_st(:, node)
            if (maximum(group)%owners(corner) > 0) &
                maximum(group)%owner_neighbors(:, corner) = &
                    data%neighbours(maximum(group)%owners(corner), :)
        end do
    end subroutine record_maximum

    subroutine verify_factor2_maximum
        real(dp), parameter :: tolerance = 1.0e-12_dp

        if (maximum(2)%element /= 5175 .or. maximum(2)%triangle /= 47661 &
                .or. maximum(2)%triangle_element /= 5987 &
                .or. maximum(2)%element_distance /= 3 &
                .or. maximum(2)%relation /= 4 &
                .or. maximum(2)%matches /= 1 &
                .or. any(maximum(2)%nodes /= [16404, 5062, 16403]) &
                .or. any(maximum(2)%owners /= [5986, 5061, 5987]) &
                .or. any(maximum(2)%neighbors <= 0) &
                .or. abs(maximum(2)%s - 5.0_dp/6.0_dp) > tolerance &
                .or. abs(maximum(2)%t - 2.0_dp/3.0_dp) > tolerance &
                .or. abs(maximum(2)%phi) > tolerance &
                .or. abs(maximum(2)%error - 0.043064226571950817_dp) &
                    > tolerance &
                .or. any(recovered_relation_counts &
                    /= [4, 297216, 47088, 20]) &
                .or. any(recovered_high_counts /= [0, 43, 12, 8]) &
                .or. maxval(abs(recovered_relation_metrics(1, :) - [ &
                    0.0051750424482673758_dp, 0.033642681741412578_dp, &
                    0.028862752750893582_dp, 0.043064226571950817_dp])) &
                    > tolerance) &
            error stop 'factor-2 recovered maximum changed'
    end subroutine verify_factor2_maximum

    subroutine build_triangle_parents
        integer :: candidate, cell_i, cell_j, split, valid_triangle
        integer :: grid_node(4), trial_node(3)

        valid_triangle = 0
        do candidate = 1, data%n_elements
            do cell_j = 0, refinement - 1
                do cell_i = 0, refinement - 1
                    grid_node = [ &
                        element_vertex(cell_i, cell_j, candidate), &
                        element_vertex(cell_i + 1, cell_j, candidate), &
                        element_vertex(cell_i + 1, cell_j + 1, candidate), &
                        element_vertex(cell_i, cell_j + 1, candidate)]
                    do split = 1, 2
                        trial_node = grid_node(triangle_vertices(:, split))
                        if (trial_node(1) == trial_node(2) &
                                .or. trial_node(1) == trial_node(3) &
                                .or. trial_node(2) == trial_node(3)) cycle
                        valid_triangle = valid_triangle + 1
                        triangle_parent(valid_triangle) = candidate
                    end do
                end do
            end do
        end do
        if (valid_triangle /= size(triangles, 1)) &
            error stop 'global triangle parent count changed'
    end subroutine build_triangle_parents

    subroutine build_boundary_triangle_parents
        integer :: candidate

        triangle_parent = 0
        do candidate = 1, data%n_elements
            triangle_parent(element_triangle_range(1, candidate): &
                element_triangle_range(2, candidate)) = candidate
        end do
        if (any(triangle_parent == 0)) &
            error stop 'boundary triangle parent map changed'
    end subroutine build_boundary_triangle_parents

    integer function parent_relation(source, parent) result(relation)
        integer, intent(in) :: source, parent

        integer :: corner

        if (parent == source) then
            relation = 1
        else if (any(data%neighbours(source, :) == parent)) then
            relation = 2
        else
            relation = 4
            do corner = 1, 4
                if (any(data%vertex(source, :) &
                        == data%vertex(parent, corner))) relation = 3
            end do
        end if
    end function parent_relation

    subroutine record_recovered_relation(source, parent, error)
        integer, intent(in) :: source, parent
        real(dp), intent(in) :: error

        integer :: relation

        relation = parent_relation(source, parent)
        recovered_relation_counts(relation) = &
            recovered_relation_counts(relation) + 1
        recovered_relation_metrics(1, relation) = &
            max(recovered_relation_metrics(1, relation), error)
        recovered_relation_metrics(2, relation) = &
            recovered_relation_metrics(2, relation) + error**2
        if (error > 0.02_dp) recovered_high_counts(relation) = &
            recovered_high_counts(relation) + 1
    end subroutine record_recovered_relation

    integer function neighbor_distance(source, target) result(distance)
        integer, intent(in) :: source, target

        integer, allocatable :: depth(:), queue(:)
        integer :: current, head, neighbor, side, tail

        allocate(depth(data%n_elements), source=-1)
        allocate(queue(data%n_elements))
        head = 1
        tail = 1
        queue(1) = source
        depth(source) = 0
        do while (head <= tail)
            current = queue(head)
            head = head + 1
            if (current == target) exit
            do side = 1, 4
                neighbor = data%neighbours(current, side)
                if (neighbor <= 0 .or. depth(neighbor) >= 0) cycle
                tail = tail + 1
                queue(tail) = neighbor
                depth(neighbor) = depth(current) + 1
            end do
        end do
        distance = depth(target)
    end function neighbor_distance

    subroutine write_maximum_output(path)
        character(len=*), intent(in) :: path

        character(len=9), parameter :: group_name(2) = &
            [character(len=9) :: 'source', 'recovered']
        integer :: group, ierr, unit

        open(newunit=unit, file=path, status='replace', action='write', &
            iostat=ierr)
        if (ierr /= 0) error stop 'cannot open field-maximum output'
        write(unit, '(A)') 'group,error,element,s,t,phi_rad,r_m,z_m,' // &
            'triangle,triangle_element,element_distance,relation,matches,' // &
            'w1,w2,w3,node1,node2,node3,owner1,owner2,owner3,' // &
            'v1_r_m,v1_z_m,v2_r_m,v2_z_m,v3_r_m,v3_z_m,' // &
            'v1_s,v1_t,v2_s,v2_t,v3_s,v3_t,' // &
            'neighbor1,neighbor2,neighbor3,neighbor4,' // &
            'owner1_neighbor1,owner1_neighbor2,owner1_neighbor3,' // &
            'owner1_neighbor4,owner2_neighbor1,owner2_neighbor2,' // &
            'owner2_neighbor3,owner2_neighbor4,owner3_neighbor1,' // &
            'owner3_neighbor2,owner3_neighbor3,owner3_neighbor4,' // &
            'exact_br_g,exact_bphi_g,exact_bz_g,' // &
            'interpolated_br_g,interpolated_bphi_g,interpolated_bz_g'
        do group = 1, 2
            write(unit, '(*(g0,:,","))') trim(group_name(group)), &
                maximum(group)%error, maximum(group)%element, &
                maximum(group)%s, maximum(group)%t, maximum(group)%phi, &
                maximum(group)%rz, maximum(group)%triangle, &
                maximum(group)%triangle_element, maximum(group)%element_distance, &
                maximum(group)%relation, maximum(group)%matches, &
                maximum(group)%weight, &
                maximum(group)%nodes, maximum(group)%owners, &
                maximum(group)%vertex_rz, maximum(group)%vertex_st, &
                maximum(group)%neighbors, maximum(group)%owner_neighbors, &
                maximum(group)%exact, maximum(group)%interpolated
        end do
        close(unit)
    end subroutine write_maximum_output

    subroutine write_topology_output(path)
        character(len=*), intent(in) :: path

        real(dp) :: rms
        integer :: ierr, relation, unit

        open(newunit=unit, file=path, status='replace', action='write', &
            iostat=ierr)
        if (ierr /= 0) error stop 'cannot open recovered-topology output'
        write(unit, '(A)') 'relation,count,above_2pct,max_error,rms_error'
        do relation = 1, 4
            rms = 0.0_dp
            if (recovered_relation_counts(relation) > 0) &
                rms = sqrt(recovered_relation_metrics(2, relation) &
                    /recovered_relation_counts(relation))
            write(unit, '(*(g0,:,","))') trim(relation_names(relation)), &
                recovered_relation_counts(relation), &
                recovered_high_counts(relation), &
                recovered_relation_metrics(1, relation), rms
        end do
        close(unit)
    end subroutine write_topology_output

    subroutine write_recovered_tail(element, s, t, phi, rz, triangle, weight, &
            matches, b_exact, b_interp, error)
        integer, intent(in) :: element, triangle, matches
        real(dp), intent(in) :: s, t, phi, rz(2), weight(3)
        real(dp), intent(in) :: b_exact(3), b_interp(3), error

        real(dp) :: chord_distance
        integer :: cell_i, cell_j, chord_nodes(2), crossed_edge, distance
        integer :: element_side, nodes(3), owners(3), parent, relation
        integer :: side_neighbor

        parent = triangle_parent(triangle)
        relation = parent_relation(element, parent)
        distance = neighbor_distance(element, parent)
        nodes = triangles(triangle, :)
        owners = vertex_element(nodes)
        call source_crossed_chord(element, s, t, rz, cell_i, cell_j, &
            crossed_edge, chord_nodes, chord_distance, element_side, &
            side_neighbor)
        write(tail_output_unit, '(*(g0,:,","))') trim(relation_names(relation)), &
            error, element, s, t, phi, rz, triangle, parent, distance, matches, &
            weight, nodes, owners, cell_i, cell_j, crossed_edge, chord_nodes, &
            chord_distance, element_side, side_neighbor, b_exact, b_interp
    end subroutine write_recovered_tail

    subroutine source_crossed_chord(element, s, t, point, cell_i, cell_j, &
            crossed_edge, chord_nodes, distance, element_side, side_neighbor)
        integer, intent(in) :: element
        real(dp), intent(in) :: s, t, point(2)
        integer, intent(out) :: cell_i, cell_j, crossed_edge, chord_nodes(2)
        real(dp), intent(out) :: distance
        integer, intent(out) :: element_side, side_neighbor

        real(dp) :: area, chord_distance, corners(2, 4), edge_vector(2)
        real(dp) :: orientation
        integer :: corner, grid_i(4), grid_j(4), next_corner, nodes(4)

        cell_i = min(refinement - 1, int(s*refinement))
        cell_j = min(refinement - 1, int(t*refinement))
        grid_i = [cell_i, cell_i + 1, cell_i + 1, cell_i]
        grid_j = [cell_j, cell_j, cell_j + 1, cell_j + 1]
        do corner = 1, 4
            nodes(corner) = element_vertex(grid_i(corner), grid_j(corner), &
                element)
            corners(:, corner) = plane_rz(:, nodes(corner))
        end do
        area = 0.0_dp
        do corner = 1, 4
            next_corner = modulo(corner, 4) + 1
            area = area + corners(1, corner)*corners(2, next_corner) &
                - corners(1, next_corner)*corners(2, corner)
        end do
        orientation = merge(1.0_dp, -1.0_dp, area >= 0.0_dp)
        distance = huge(1.0_dp)
        crossed_edge = 0
        do corner = 1, 4
            next_corner = modulo(corner, 4) + 1
            edge_vector = corners(:, next_corner) - corners(:, corner)
            chord_distance = orientation*(edge_vector(1) &
                *(point(2) - corners(2, corner)) - edge_vector(2) &
                *(point(1) - corners(1, corner)))/sqrt(sum(edge_vector**2))
            if (chord_distance < distance) then
                distance = chord_distance
                crossed_edge = corner
            end if
        end do
        chord_nodes = [nodes(crossed_edge), &
            nodes(modulo(crossed_edge, 4) + 1)]
        call source_element_side(cell_i, cell_j, crossed_edge, element_side)
        side_neighbor = 0
        if (element_side > 0) side_neighbor = &
            data%neighbours(element, element_side)
    end subroutine source_crossed_chord

    subroutine source_element_side(cell_i, cell_j, crossed_edge, element_side)
        integer, intent(in) :: cell_i, cell_j, crossed_edge
        integer, intent(out) :: element_side

        element_side = 0
        if (crossed_edge == 1 .and. cell_j == 0) element_side = 1
        if (crossed_edge == 2 .and. cell_i == refinement - 1) element_side = 2
        if (crossed_edge == 3 .and. cell_j == refinement - 1) element_side = 3
        if (crossed_edge == 4 .and. cell_i == 0) element_side = 4
    end subroutine source_element_side

    subroutine record_tail_chord(element, s, t, point)
        integer, intent(in) :: element
        real(dp), intent(in) :: s, t, point(2)

        real(dp) :: distance
        integer :: cell_i, cell_j, chord_nodes(2), crossed_edge, edge_record
        integer :: element_side, point_record, side_neighbor, sorted_nodes(2)

        call source_crossed_chord(element, s, t, point, cell_i, cell_j, &
            crossed_edge, chord_nodes, distance, element_side, side_neighbor)
        if (distance >= 0.0_dp .or. element_side <= 0 .or. side_neighbor <= 0) &
            error stop 'recovered tail does not cross an interior element side'
        sorted_nodes = [minval(chord_nodes), maxval(chord_nodes)]
        edge_record = 0
        do point_record = 1, tail_edge_count
            if (all(tail_edge_nodes(:, point_record) == sorted_nodes)) &
                edge_record = point_record
        end do
        if (edge_record == 0) then
            tail_edge_count = tail_edge_count + 1
            if (tail_edge_count > max_tail_records) &
                error stop 'too many recovered tail chords'
            edge_record = tail_edge_count
            tail_edge_nodes(:, edge_record) = sorted_nodes
        end if
        tail_edge_phase_counts(edge_record) = &
            tail_edge_phase_counts(edge_record) + 1
        do point_record = 1, tail_point_count
            if (tail_point_element(point_record) == element &
                    .and. tail_point_s(point_record) == s &
                    .and. tail_point_t(point_record) == t) return
        end do
        tail_point_count = tail_point_count + 1
        if (tail_point_count > max_tail_records) &
            error stop 'too many recovered tail points'
        tail_point_element(tail_point_count) = element
        tail_point_s(tail_point_count) = s
        tail_point_t(tail_point_count) = t
    end subroutine record_tail_chord

    logical function has_coincident_corners(corners)
        real(dp), intent(in) :: corners(2, 4)

        integer :: i, j

        has_coincident_corners = .false.
        do i = 1, 3
            do j = i + 1, 4
                if (maxval(abs(corners(:, i) - corners(:, j))) <= 1.0e-12_dp) &
                    has_coincident_corners = .true.
            end do
        end do
    end function has_coincident_corners

    subroutine barycentric_weights(vertices, point, weight, ierr)
        real(dp), intent(in) :: vertices(2, 3), point(2)
        real(dp), intent(out) :: weight(3)
        integer, intent(out) :: ierr

        real(dp) :: determinant

        determinant = (vertices(1, 2) - vertices(1, 1)) &
            *(vertices(2, 3) - vertices(2, 1)) &
            - (vertices(1, 3) - vertices(1, 1)) &
            *(vertices(2, 2) - vertices(2, 1))
        weight = 0.0_dp
        ierr = 1
        if (abs(determinant) <= 64.0_dp*epsilon(1.0_dp) &
                *max(1.0_dp, maxval(abs(vertices)))**2) return
        weight(2) = ((point(1) - vertices(1, 1)) &
            *(vertices(2, 3) - vertices(2, 1)) &
            - (vertices(1, 3) - vertices(1, 1)) &
            *(point(2) - vertices(2, 1)))/determinant
        weight(3) = ((vertices(1, 2) - vertices(1, 1)) &
            *(point(2) - vertices(2, 1)) &
            - (point(1) - vertices(1, 1)) &
            *(vertices(2, 2) - vertices(2, 1)))/determinant
        weight(1) = 1.0_dp - weight(2) - weight(3)
        ierr = 0
    end subroutine barycentric_weights

    subroutine verify_indexed_location(point, indexed, indexed_weight, &
            indexed_matches, indexed_ierr)
        real(dp), intent(in) :: point(2), indexed_weight(3)
        integer, intent(in) :: indexed, indexed_matches, indexed_ierr

        real(dp) :: brute_weight(3), first_weight(3)
        integer :: brute, brute_ierr, brute_matches, first

        first = 0
        first_weight = 0.0_dp
        brute_matches = 0
        do brute = 1, size(triangles, 1)
            call barycentric_weights(plane_rz(:, triangles(brute, :)), point, &
                brute_weight, brute_ierr)
            if (brute_ierr == 0) then
                if (all(brute_weight >= -1.0e-10_dp) &
                        .and. all(brute_weight <= 1.0_dp + 1.0e-10_dp)) then
                    brute_matches = brute_matches + 1
                    if (first == 0) then
                        first = brute
                        first_weight = brute_weight
                    end if
                end if
            end if
        end do
        if (indexed_matches /= brute_matches .or. indexed /= first) &
            error stop 'indexed triangle match differs from brute force'
        if (brute_matches == 0) then
            if (indexed_ierr == 0) error stop 'indexed outside status changed'
        else
            if (indexed_ierr /= 0 &
                    .or. maxval(abs(indexed_weight - first_weight)) > 1.0e-12_dp) &
                error stop 'indexed barycentric weight differs from brute force'
        end if
    end subroutine verify_indexed_location

end program test_jorek_global_plane_linearization
