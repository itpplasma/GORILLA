program test_jorek_poloidal_linearization
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use jorek_bezier, only: jorek_locator_t, build_jorek_locator, &
        evaluate_jorek_geometry
    use jorek_field_backend_mod, only: jorek_chart_requires_global
    use jorek_model303_field, only: evaluate_jorek_model303_a, &
        evaluate_jorek_model303_at, evaluate_jorek_model303_b
    use jorek_refined_plane_mod, only: extract_refined_jorek_plane
    use jorek_restart, only: jorek_restart_t, load_jorek_restart

    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: corner_s(4) = [0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp]
    real(dp), parameter :: corner_t(4) = [0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp]
    real(dp), parameter :: weights(3, 4) = reshape([ &
        1.0_dp/3, 1.0_dp/3, 1.0_dp/3, 0.6_dp, 0.2_dp, 0.2_dp, &
        0.2_dp, 0.6_dp, 0.2_dp, 0.2_dp, 0.2_dp, 0.6_dp], [3, 4])
    integer, parameter :: triangle_vertices(3, 2) = reshape([1, 2, 3, 1, 3, 4], [3, 2])
    integer, parameter :: refinements(4) = [1, 2, 4, 8]
    integer, parameter :: reference_refinement = 2

    type(jorek_restart_t) :: data
    type(jorek_locator_t) :: locator
    logical, allocatable :: axis_mask(:), uncovered_elements(:)
    logical :: write_metrics
    real(dp), allocatable :: element_metrics(:, :)
    integer, allocatable :: element_covered(:), element_outside(:)
    real(dp), allocatable :: mesh_psi(:), mesh_rz(:, :), mesh_st(:, :)
    integer, allocatable :: mesh_nodes(:, :, :), mesh_owner(:), mesh_triangles(:, :)
    real(dp) :: corner_rz(2, 4), metrics(6), phi, rz_st(2, 2), target_metrics(2)
    real(dp) :: split_metrics(2, 2)
    real(dp) :: outside_distance(2)
    real(dp) :: worst_b_exact(3), worst_b_interp(3), worst_rz(2)
    character(len=1024) :: filename, metrics_filename
    integer :: axis_elements, element, i, ierr, j, level, metrics_unit, outside, phase
    integer :: refinement, sample, samples, triangle, vertex
    integer :: fallback_samples, regular_samples
    real(dp) :: target_s, target_t
    integer :: worst_case(4), worst_found
    integer :: worst_corner_owner(4)
    real(dp) :: worst_st(2)

    if (command_argument_count() < 1 .or. command_argument_count() > 2) &
        error stop 'restart path and optional metrics path are required'
    call get_command_argument(1, filename)
    call load_jorek_restart(trim(filename), data, ierr)
    if (ierr /= 0) error stop 'JOREK restart load failed'
    call build_jorek_locator(data, locator, ierr)
    if (ierr /= 0) error stop 'JOREK locator build failed'
    allocate(axis_mask(data%n_elements), uncovered_elements(data%n_elements))
    allocate(element_metrics(6, data%n_elements))
    allocate(element_covered(data%n_elements), element_outside(data%n_elements))
    metrics_unit = -1
    write_metrics = .false.
    if (command_argument_count() == 2) then
        call get_command_argument(2, metrics_filename)
        open(newunit=metrics_unit, file=trim(metrics_filename), &
            status='replace', action='write', iostat=ierr)
        if (ierr /= 0) error stop 'cannot open JOREK metrics output'
        write_metrics = .true.
        write(metrics_unit, '(A)') 'refinement,element,axis,covered,outside,' // &
            'b_max_relative,b_rms_relative,miss_max_cm,miss_rms_cm,' // &
            'target_b_max_relative,target_b_rms_relative'
    end if

    do level = 1, size(refinements)
        refinement = refinements(level)
        metrics = 0.0_dp
        target_metrics = 0.0_dp
        split_metrics = 0.0_dp
        worst_case = 0
        axis_elements = 0
        outside = 0
        outside_distance = 0.0_dp
        samples = 0
        fallback_samples = 0
        regular_samples = 0
        axis_mask = .false.
        element_covered = 0
        element_outside = 0
        element_metrics = 0.0_dp
        uncovered_elements = .false.
        if (refinement > 1) then
            call extract_refined_jorek_plane(data, refinement, mesh_rz, &
                mesh_psi, mesh_triangles, mesh_owner, mesh_st, ierr, mesh_nodes)
            if (ierr /= 0) error stop 'production JOREK refinement failed'
        end if
        do element = 1, data%n_elements
            do vertex = 1, 4
                call evaluate_jorek_geometry(data, element, corner_s(vertex), &
                    corner_t(vertex), corner_rz(:, vertex), rz_st, ierr)
                if (ierr /= 0) error stop 'JOREK corner geometry failed'
            end do
            if (has_coincident_corners(corner_rz)) then
                axis_elements = axis_elements + 1
                axis_mask(element) = .true.
                cycle
            end if
            do phase = 0, 3
                phi = 0.5_dp*pi*phase
                do j = 0, reference_refinement - 1
                    do i = 0, reference_refinement - 1
                        do triangle = 1, 2
                            do sample = 1, size(weights, 2)
                                call reference_parameters(i, j, triangle, &
                                    weights(:, sample), target_s, target_t)
                                call compare_fixed_sample(data, phi, &
                                    target_s, target_t, refinement, metrics, &
                                    samples, outside)
                            end do
                        end do
                    end do
                end do
            end do
        end do
        if (samples == 0) error stop 'no interior linearization samples'
        metrics(2) = sqrt(metrics(2)/samples)
        metrics(4) = sqrt(metrics(4)/samples)
        metrics(6) = sqrt(metrics(6)/samples)
        target_metrics(2) = sqrt(target_metrics(2)/samples)
        if (refinement == 1 &
                .and. maxval(abs(target_metrics - metrics(1:2))) > 1.0e-14_dp) &
            error stop 'factor-1 owner interpolation identity failed'
        if (samples + outside /= 769792) &
            error stop 'fixed poloidal sample count changed'
        if (sum(element_covered) /= samples &
                .or. sum(element_outside) /= outside &
                .or. .not. all(ieee_is_finite(element_metrics))) &
            error stop 'per-element metric partition failed'
        if (fallback_samples + regular_samples /= samples &
                .or. .not. all(ieee_is_finite(split_metrics))) &
            error stop 'chart-fallback metric partition failed'
        if (outside == 0) error stop 'coverage-failure fixture changed'
        print '(A, I0, A, I0, A, I0, A, I0)', 'refinement=', refinement, &
            ' samples=', samples, ' outside=', outside, &
            ' axis_elements=', axis_elements
        print '(A, 6(ES14.6, 1X))', &
            'max/rms B, A, |B| relative errors: ', metrics
        print '(A, 2(ES14.6, 1X))', &
            'target-owner max/rms B relative errors: ', target_metrics
        print '(A, I0, A, 2(ES14.6, 1X))', 'fallback samples=', &
            fallback_samples, ' max/rms target-owner B errors: ', &
            split_metrics(1, 1), &
            sqrt(split_metrics(2, 1)/max(1, fallback_samples))
        print '(A, I0, A, 2(ES14.6, 1X))', 'regular samples=', &
            regular_samples, ' max/rms target-owner B errors: ', &
            split_metrics(1, 2), &
            sqrt(split_metrics(2, 2)/max(1, regular_samples))
        print '(A, I0)', 'uncovered elements=', count(uncovered_elements)
        if (outside > 0) then
            outside_distance(2) = sqrt(outside_distance(2)/outside)
            print '(A, 2(ES14.6, 1X))', &
                'max/rms uncovered distance [cm]: ', 100.0_dp*outside_distance
        end if
        print '(A, 4(I0, 1X), A, I0, A, 2(ES14.6, 1X), A, 2(ES14.6, 1X))', &
            'worst B element/phase/i/j ', worst_case, &
            'found ', worst_found, ' st ', worst_st, ' RZ ', worst_rz
        print '(A, 3(ES14.6, 1X))', 'worst B exact ', worst_b_exact
        print '(A, 3(ES14.6, 1X))', 'worst B interpolated ', worst_b_interp
        print '(A, 4(I0, 1X))', 'worst corner owners ', worst_corner_owner
        if (write_metrics) call write_element_metrics(metrics_unit, &
            refinement, axis_mask, element_covered, element_outside, &
            element_metrics)
    end do
    if (write_metrics) close(metrics_unit)
    call print_axis_limit(data, 0.0_dp)
    call print_axis_limit(data, pi)
    call print_regularized_axis(data, locator, 0.0_dp)
    call print_regularized_axis(data, locator, pi)

contains

    subroutine reference_parameters(cell_i, cell_j, triangle, weight, s, t)
        integer, intent(in) :: cell_i, cell_j, triangle
        real(dp), intent(in) :: weight(3)
        real(dp), intent(out) :: s, t

        real(dp) :: local_s(4), local_t(4)
        integer :: index, vertex

        local_s = real([cell_i, cell_i + 1, cell_i + 1, cell_i], dp) &
            /reference_refinement
        local_t = real([cell_j, cell_j, cell_j + 1, cell_j + 1], dp) &
            /reference_refinement
        s = 0.0_dp
        t = 0.0_dp
        do index = 1, 3
            vertex = triangle_vertices(index, triangle)
            s = s + weight(index)*local_s(vertex)
            t = t + weight(index)*local_t(vertex)
        end do
    end subroutine reference_parameters

    subroutine compare_fixed_sample(data, phi, target_s, target_t, &
            refinement, metrics, samples, uncovered)
        type(jorek_restart_t), intent(in) :: data
        real(dp), intent(in) :: phi, target_s, target_t
        integer, intent(in) :: refinement
        real(dp), intent(inout) :: metrics(6)
        integer, intent(inout) :: samples, uncovered

        real(dp) :: a_corner(3, 4), a_exact(3), a_interp(3), a_jorek(3)
        real(dp) :: a_target_corner(3, 4), a_target_interp(3)
        real(dp) :: barycentric(3), b_exact(3), b_interp(3), b_jorek(3)
        real(dp) :: b_target_interp(3)
        real(dp) :: b_error, miss_distance
        real(dp) :: bmod_corner(4), bmod_exact, bmod_interp, bmod_target(4)
        real(dp) :: bmod_target_interp
        real(dp) :: corner_rz(2, 4), h_corner(3, 4), h_exact(3)
        real(dp) :: h_target_corner(3, 4), rz(2), rz_st(2, 2), st(2)
        integer :: cell_i, cell_j, corner_owner(4), found, ierr, vertices(3)
        logical :: target_fallback(4)

        cell_i = min(refinement - 1, int(target_s*refinement))
        cell_j = min(refinement - 1, int(target_t*refinement))
        call prepare_sample_corners(data, phi, refinement, cell_i, cell_j, &
            corner_rz, corner_owner, a_corner, h_corner, bmod_corner, &
            a_target_corner, h_target_corner, bmod_target, target_fallback)
        call evaluate_jorek_geometry(data, element, target_s, target_t, rz, &
            rz_st, ierr)
        if (ierr /= 0) error stop 'JOREK fixed-point geometry failed'
        call find_sample_triangle(corner_rz, rz, target_s, target_t, &
            refinement, cell_i, cell_j, vertices, barycentric, ierr)
        if (ierr /= 0) then
            uncovered = uncovered + 1
            element_outside(element) = element_outside(element) + 1
            uncovered_elements(element) = .true.
            miss_distance = distance_to_triangles(corner_rz, rz)
            call update_metrics(miss_distance, outside_distance)
            call update_metrics(miss_distance, element_metrics(3:4, element))
            return
        end if
        call evaluate_jorek_model303_a(data, element, target_s, target_t, &
            phi, a_jorek, ierr)
        if (ierr == 0) call evaluate_jorek_model303_b(data, element, target_s, &
            target_t, phi, b_jorek, ierr)
        if (ierr /= 0) error stop 'JOREK fixed-point field evaluation failed'
        found = element
        st = [target_s, target_t]
        call convert_components(rz(1)*100.0_dp, a_jorek, b_jorek, a_exact, &
            h_exact, bmod_exact)
        call interpolate_corner_field(vertices, barycentric, rz(1), a_corner, &
            h_corner, bmod_corner, a_interp, b_interp, bmod_interp)
        call interpolate_corner_field(vertices, barycentric, rz(1), &
            a_target_corner, h_target_corner, bmod_target, a_target_interp, &
            b_target_interp, bmod_target_interp)
        a_exact = [a_exact(1), a_exact(2)/(rz(1)*100.0_dp), a_exact(3)]
        b_exact = [h_exact(1)*bmod_exact, &
            h_exact(2)*bmod_exact/(rz(1)*100.0_dp), h_exact(3)*bmod_exact]
        b_error = relative_error(b_interp, b_exact)
        if (b_error > metrics(1)) then
            worst_case = [element, phase, cell_i, cell_j]
            worst_found = found
            worst_st = st
            worst_rz = rz
            worst_b_exact = b_exact
            worst_b_interp = b_interp
            worst_corner_owner = corner_owner
        end if
        call update_metrics(b_error, metrics(1:2))
        call update_metrics(b_error, element_metrics(1:2, element))
        call update_metrics(relative_error(b_target_interp, b_exact), &
            target_metrics)
        call update_metrics(relative_error(b_target_interp, b_exact), &
            element_metrics(5:6, element))
        if (any(target_fallback(vertices))) then
            call update_metrics(relative_error(b_target_interp, b_exact), &
                split_metrics(:, 1))
            fallback_samples = fallback_samples + 1
        else
            call update_metrics(relative_error(b_target_interp, b_exact), &
                split_metrics(:, 2))
            regular_samples = regular_samples + 1
        end if
        call update_metrics(relative_error(a_interp, a_exact), metrics(3:4))
        call update_metrics(abs(bmod_interp/bmod_exact - 1.0_dp), metrics(5:6))
        samples = samples + 1
        element_covered(element) = element_covered(element) + 1
    end subroutine compare_fixed_sample

    subroutine prepare_sample_corners(data, phi, refinement, cell_i, cell_j, &
            corner_rz, corner_owner, a_corner, h_corner, bmod_corner, &
            a_target, h_target, bmod_target, target_fallback)
        type(jorek_restart_t), intent(in) :: data
        real(dp), intent(in) :: phi
        integer, intent(in) :: refinement, cell_i, cell_j
        real(dp), intent(out) :: corner_rz(2, 4), a_corner(3, 4), h_corner(3, 4)
        real(dp), intent(out) :: bmod_corner(4), a_target(3, 4), h_target(3, 4)
        real(dp), intent(out) :: bmod_target(4)
        integer, intent(out) :: corner_owner(4)
        logical, intent(out) :: target_fallback(4)

        real(dp) :: rz_st(2, 2), sub_s(4), sub_t(4)
        integer :: grid_i(4), grid_j(4), ierr, index, mesh_node, owner

        sub_s = real([cell_i, cell_i + 1, cell_i + 1, cell_i], dp)/refinement
        sub_t = real([cell_j, cell_j, cell_j + 1, cell_j + 1], dp)/refinement
        grid_i = [cell_i, cell_i + 1, cell_i + 1, cell_i]
        grid_j = [cell_j, cell_j, cell_j + 1, cell_j + 1]
        do index = 1, 4
            if (refinement == 1) then
                call evaluate_jorek_geometry(data, element, sub_s(index), &
                    sub_t(index), corner_rz(:, index), rz_st, ierr)
                if (ierr /= 0) error stop 'JOREK corner geometry failed'
                owner = element
            else
                mesh_node = mesh_nodes(grid_i(index), grid_j(index), element)
                corner_rz(:, index) = mesh_rz(:, mesh_node)
                owner = mesh_owner(mesh_node)
                if (owner == 0) error stop 'non-axis sample uses axis owner'
            end if
            corner_owner(index) = owner
            call jorek_chart_requires_global(data, element, sub_s(index), &
                sub_t(index), target_fallback(index), ierr)
            if (ierr /= 0) error stop 'JOREK target chart test failed'
            if (refinement == 1) then
                call evaluate_owner_corner(data, owner, sub_s(index), &
                    sub_t(index), phi, corner_rz(:, index), a_corner(:, index), &
                    h_corner(:, index), bmod_corner(index))
            else
                call evaluate_owner_corner(data, owner, mesh_st(1, mesh_node), &
                    mesh_st(2, mesh_node), phi, corner_rz(:, index), &
                    a_corner(:, index), h_corner(:, index), bmod_corner(index))
            end if
            if (owner == element) then
                a_target(:, index) = a_corner(:, index)
                h_target(:, index) = h_corner(:, index)
                bmod_target(index) = bmod_corner(index)
            else
                call evaluate_owner_corner(data, element, sub_s(index), &
                    sub_t(index), phi, corner_rz(:, index), a_target(:, index), &
                    h_target(:, index), bmod_target(index))
            end if
        end do
    end subroutine prepare_sample_corners

    subroutine find_sample_triangle(corners, point, target_s, target_t, &
            refinement, cell_i, cell_j, vertices, barycentric, ierr)
        real(dp), intent(in) :: corners(2, 4), point(2), target_s, target_t
        integer, intent(in) :: refinement, cell_i, cell_j
        integer, intent(out) :: vertices(3), ierr
        real(dp), intent(out) :: barycentric(3)

        if ((target_t*refinement - cell_j) &
                <= (target_s*refinement - cell_i)) then
            vertices = triangle_vertices(:, 1)
        else
            vertices = triangle_vertices(:, 2)
        end if
        call physical_barycentric(corners(:, vertices), point, barycentric, ierr)
        if (valid_barycentric(barycentric, ierr)) return
        if (all(vertices == triangle_vertices(:, 1))) then
            vertices = triangle_vertices(:, 2)
        else
            vertices = triangle_vertices(:, 1)
        end if
        call physical_barycentric(corners(:, vertices), point, barycentric, ierr)
        if (.not. valid_barycentric(barycentric, ierr)) ierr = 1
    end subroutine find_sample_triangle

    logical function valid_barycentric(weight, ierr)
        real(dp), intent(in) :: weight(3)
        integer, intent(in) :: ierr

        valid_barycentric = ierr == 0 .and. all(weight >= -1.0e-10_dp) &
            .and. all(weight <= 1.0_dp + 1.0e-10_dp)
    end function valid_barycentric

    subroutine interpolate_corner_field(vertices, weight, r, a_corner, &
            h_corner, bmod_corner, a_interp, b_interp, bmod_interp)
        integer, intent(in) :: vertices(3)
        real(dp), intent(in) :: weight(3), r, a_corner(3, 4), h_corner(3, 4)
        real(dp), intent(in) :: bmod_corner(4)
        real(dp), intent(out) :: a_interp(3), b_interp(3), bmod_interp

        real(dp) :: h_interp(3), r_cm
        integer :: index

        a_interp = 0.0_dp
        h_interp = 0.0_dp
        bmod_interp = 0.0_dp
        do index = 1, 3
            a_interp = a_interp + weight(index)*a_corner(:, vertices(index))
            h_interp = h_interp + weight(index)*h_corner(:, vertices(index))
            bmod_interp = bmod_interp + weight(index)*bmod_corner(vertices(index))
        end do
        r_cm = 100.0_dp*r
        a_interp = [a_interp(1), a_interp(2)/r_cm, a_interp(3)]
        b_interp = [h_interp(1)*bmod_interp, &
            h_interp(2)*bmod_interp/r_cm, h_interp(3)*bmod_interp]
    end subroutine interpolate_corner_field

    subroutine write_element_metrics(unit, refinement, axis, covered, &
            outside, values)
        integer, intent(in) :: unit, refinement, covered(:), outside(:)
        logical, intent(in) :: axis(:)
        real(dp), intent(in) :: values(:, :)

        real(dp) :: b_rms, miss_rms, target_b_rms
        integer :: element

        do element = 1, size(axis)
            b_rms = 0.0_dp
            miss_rms = 0.0_dp
            target_b_rms = 0.0_dp
            if (covered(element) > 0) &
                b_rms = sqrt(values(2, element)/covered(element))
            if (outside(element) > 0) &
                miss_rms = 100.0_dp*sqrt(values(4, element)/outside(element))
            if (covered(element) > 0) &
                target_b_rms = sqrt(values(6, element)/covered(element))
            write(unit, '(I0,",",I0,",",I0,",",I0,",",I0,6(",",ES24.16E3))') &
                refinement, element, merge(1, 0, axis(element)), &
                covered(element), outside(element), values(1, element), &
                b_rms, 100.0_dp*values(3, element), miss_rms, &
                values(5, element), target_b_rms
        end do
    end subroutine write_element_metrics

    subroutine evaluate_owner_corner(data, element, s, t, phi, rz, a, h, bmod)
        type(jorek_restart_t), intent(in) :: data
        integer, intent(in) :: element
        real(dp), intent(in) :: s, t, phi, rz(2)
        real(dp), intent(out) :: a(3), h(3), bmod

        real(dp) :: a_jorek(3), b_jorek(3)
        real(dp) :: st_found(2)
        integer :: found, ierr
        logical :: use_global

        call jorek_chart_requires_global(data, element, s, t, use_global, ierr)
        if (ierr /= 0) error stop 'JOREK chart test failed'
        if (use_global) then
            call evaluate_jorek_model303_at(data, rz, phi, a_jorek, b_jorek, &
                found, st_found, ierr, locator)
        else
            call evaluate_jorek_model303_a(data, element, s, t, phi, &
                a_jorek, ierr)
            if (ierr == 0) call evaluate_jorek_model303_b(data, element, s, t, &
                phi, b_jorek, ierr)
        end if
        if (ierr /= 0) error stop 'JOREK owner corner field failed'
        call convert_components(rz(1)*100.0_dp, a_jorek, b_jorek, a, h, bmod)
    end subroutine evaluate_owner_corner

    subroutine physical_barycentric(vertices, point, weight, ierr)
        real(dp), intent(in) :: vertices(2, 3), point(2)
        real(dp), intent(out) :: weight(3)
        integer, intent(out) :: ierr

        real(dp) :: determinant, scale

        determinant = (vertices(1, 2) - vertices(1, 1)) &
            *(vertices(2, 3) - vertices(2, 1)) &
            - (vertices(1, 3) - vertices(1, 1)) &
            *(vertices(2, 2) - vertices(2, 1))
        scale = max(1.0_dp, maxval(abs(vertices)))
        ierr = 0
        weight = 0.0_dp
        if (abs(determinant) <= 64.0_dp*epsilon(1.0_dp)*scale**2) then
            ierr = 1
            return
        end if
        weight(2) = ((point(1) - vertices(1, 1)) &
            *(vertices(2, 3) - vertices(2, 1)) &
            - (vertices(1, 3) - vertices(1, 1)) &
            *(point(2) - vertices(2, 1)))/determinant
        weight(3) = ((vertices(1, 2) - vertices(1, 1)) &
            *(point(2) - vertices(2, 1)) &
            - (point(1) - vertices(1, 1)) &
            *(vertices(2, 2) - vertices(2, 1)))/determinant
        weight(1) = 1.0_dp - weight(2) - weight(3)
    end subroutine physical_barycentric

    real(dp) function distance_to_triangles(corners, point)
        real(dp), intent(in) :: corners(2, 4), point(2)

        integer, parameter :: edges(2, 5) = reshape([ &
            1, 2, 2, 3, 3, 4, 4, 1, 1, 3], [2, 5])
        integer :: edge

        distance_to_triangles = huge(1.0_dp)
        do edge = 1, size(edges, 2)
            distance_to_triangles = min(distance_to_triangles, &
                point_segment_distance(point, corners(:, edges(1, edge)), &
                    corners(:, edges(2, edge))))
        end do
    end function distance_to_triangles

    real(dp) function point_segment_distance(point, left, right)
        real(dp), intent(in) :: point(2), left(2), right(2)

        real(dp) :: fraction, segment(2), scale

        segment = right - left
        scale = sum(segment**2)
        if (scale <= tiny(1.0_dp)) then
            point_segment_distance = sqrt(sum((point - left)**2))
            return
        end if
        fraction = dot_product(point - left, segment)/scale
        fraction = min(1.0_dp, max(0.0_dp, fraction))
        point_segment_distance = &
            sqrt(sum((point - left - fraction*segment)**2))
    end function point_segment_distance

    subroutine convert_components(r, a_jorek, b_jorek, a, h, bmod)
        real(dp), intent(in) :: r, a_jorek(3), b_jorek(3)
        real(dp), intent(out) :: a(3), h(3), bmod

        real(dp) :: b_physical(3)

        a = a_jorek*1.0e6_dp
        a(2) = a(2)*100.0_dp
        b_physical = b_jorek*1.0e4_dp
        h = [b_physical(1), r*b_physical(3), b_physical(2)]
        bmod = sqrt(sum(b_physical**2))
        h = h/bmod
    end subroutine convert_components

    real(dp) function relative_error(actual, expected)
        real(dp), intent(in) :: actual(3), expected(3)

        relative_error = sqrt(sum((actual - expected)**2)) &
            /max(sqrt(sum(expected**2)), tiny(1.0_dp))
    end function relative_error

    subroutine update_metrics(error, metric)
        real(dp), intent(in) :: error
        real(dp), intent(inout) :: metric(2)

        metric(1) = max(metric(1), error)
        metric(2) = metric(2) + error**2
    end subroutine update_metrics

    logical function has_coincident_corners(rz)
        real(dp), intent(in) :: rz(2, 4)

        integer :: i, j

        has_coincident_corners = .false.
        do i = 1, 3
            do j = i + 1, 4
                if (maxval(abs(rz(:, i) - rz(:, j))) <= 1.0e-12_dp) &
                    has_coincident_corners = .true.
            end do
        end do
    end function has_coincident_corners

    subroutine print_axis_limit(data, phi)
        type(jorek_restart_t), intent(in) :: data
        real(dp), intent(in) :: phi

        real(dp) :: b(3), b2_sum(3), b_sum(3), corners(2, 4), mean(3)
        real(dp) :: rz(2), rz_st(2, 2), s
        integer :: count_axis, element, ierr, power, vertex

        print '(A, ES12.4)', 'axis-limit phi=', phi
        do power = 1, 5
            s = 10.0_dp**(-power)
            b_sum = 0.0_dp
            b2_sum = 0.0_dp
            count_axis = 0
            do element = 1, data%n_elements
                do vertex = 1, 4
                    call evaluate_jorek_geometry(data, element, corner_s(vertex), &
                        corner_t(vertex), corners(:, vertex), rz_st, ierr)
                end do
                if (.not. has_coincident_corners(corners)) cycle
                call evaluate_jorek_geometry(data, element, s, 0.5_dp, rz, &
                    rz_st, ierr)
                call evaluate_jorek_model303_b(data, element, s, 0.5_dp, phi, &
                    b, ierr)
                b_sum = b_sum + b
                b2_sum = b2_sum + b**2
                count_axis = count_axis + 1
            end do
            mean = b_sum/count_axis
            print '(ES10.2, 1X, 3(ES12.4, 1X), A, 3(ES12.4, 1X))', &
                s, mean, 'std ', sqrt(max(0.0_dp, b2_sum/count_axis - mean**2))
        end do
    end subroutine print_axis_limit

    subroutine print_regularized_axis(data, locator, phi)
        type(jorek_restart_t), intent(in) :: data
        type(jorek_locator_t), intent(in) :: locator
        real(dp), intent(in) :: phi

        real(dp) :: a(3), b(3), st(2)
        integer :: element, ierr

        call evaluate_jorek_model303_at(data, locator%axis_rz, phi, a, b, &
            element, st, ierr, locator)
        print '(A, ES12.4, A, I0, A, I0, A, 3(ES14.6, 1X))', &
            'regularized axis phi=', phi, ' ierr=', ierr, ' element=', element, &
            ' B=', b
    end subroutine print_regularized_axis

end program test_jorek_poloidal_linearization
