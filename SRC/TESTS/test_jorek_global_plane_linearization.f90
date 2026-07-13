program test_jorek_global_plane_linearization
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use jorek_bezier, only: jorek_locator_t, build_jorek_locator, &
        evaluate_jorek_geometry
    use jorek_field_backend_mod, only: evaluate_jorek_model303_gorilla, &
        evaluate_jorek_model303_gorilla_element
    use jorek_model303_field, only: evaluate_jorek_model303_b
    use jorek_refined_plane_mod, only: extract_refined_jorek_plane
    use jorek_restart, only: jorek_restart_t, load_jorek_restart
    use jorek_triangle_locator_mod, only: build_jorek_triangle_locator, &
        jorek_triangle_locator_t, locate_jorek_triangle

    implicit none

    integer, parameter :: refinement = 8
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: phases(4) = [0.0_dp, 0.5_dp*pi, pi, 1.5_dp*pi]
    real(dp), parameter :: weights(3, 4) = reshape([ &
        1.0_dp/3, 1.0_dp/3, 1.0_dp/3, 0.6_dp, 0.2_dp, 0.2_dp, &
        0.2_dp, 0.6_dp, 0.2_dp, 0.2_dp, 0.2_dp, 0.6_dp], [3, 4])
    integer, parameter :: triangle_vertices(3, 2) = &
        reshape([1, 2, 3, 1, 3, 4], [3, 2])
    real(dp), parameter :: corner_s(4) = [0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp]
    real(dp), parameter :: corner_t(4) = [0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp]
    type(jorek_restart_t) :: data
    type(jorek_locator_t) :: field_locator
    type(jorek_triangle_locator_t) :: plane_locator
    real(dp), allocatable :: plane_rz(:, :), psi(:), vertex_st(:, :)
    integer, allocatable :: element_vertex(:, :, :), triangles(:, :)
    integer, allocatable :: vertex_element(:)
    real(dp) :: barycentric(3), corner_rz(2, 4), metrics(2, 2), rz(2)
    real(dp) :: rz_st(2, 2), s, t
    character(len=1024) :: filename
    integer :: ambiguous, element, global_outside, i, ierr, j, located
    integer :: matches, neighbor_recovered
    integer :: phase, sample, source_covered, source_outside, triangle
    integer :: unique_sample
    integer :: metric_counts(2)
    logical :: source_hit

    if (command_argument_count() /= 1) error stop 'restart path argument is required'
    call get_command_argument(1, filename)
    call load_jorek_restart(trim(filename), data, ierr)
    if (ierr /= 0) error stop 'JOREK restart load failed'
    call build_jorek_locator(data, field_locator, ierr)
    if (ierr /= 0) error stop 'JOREK field locator build failed'
    call extract_refined_jorek_plane(data, refinement, plane_rz, psi, &
        triangles, vertex_element, vertex_st, ierr, element_vertex)
    if (ierr /= 0) error stop 'JOREK plane refinement failed'
    call build_jorek_triangle_locator(plane_rz, triangles, plane_locator, ierr)
    if (ierr /= 0) error stop 'JOREK triangle locator build failed'
    global_outside = 0
    ambiguous = 0
    neighbor_recovered = 0
    source_covered = 0
    source_outside = 0
    metrics = 0.0_dp
    metric_counts = 0
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
                        if (matches > 1) ambiguous = ambiguous + size(phases)
                        source_hit = source_contains(element, s, t, rz)
                        if (source_hit) then
                            source_covered = source_covered + size(phases)
                        else
                            source_outside = source_outside + size(phases)
                            if (ierr == 0) &
                                neighbor_recovered = neighbor_recovered + size(phases)
                        end if
                        if (ierr /= 0) then
                            global_outside = global_outside + size(phases)
                            cycle
                        end if
                        do phase = 1, size(phases)
                            call compare_field(element, s, t, phases(phase), rz, &
                                located, barycentric, merge(1, 2, source_hit))
                        end do
                    end do
                end do
            end do
        end do
    end do
    if (source_covered /= 765124 .or. source_outside /= 4668) &
        error stop 'source-element coverage fixture changed'
    if (source_covered + source_outside /= 769792 &
            .or. source_covered + source_outside - global_outside &
                /= sum(metric_counts) &
            .or. .not. all(ieee_is_finite(metrics))) &
        error stop 'global-plane sample partition failed'
    print '(A, I0, A, I0, A, I0)', 'source outside=', source_outside, &
        ' recovered by global plane=', neighbor_recovered, &
        ' global outside=', global_outside
    print '(A, I0)', 'multiple containing global triangles=', ambiguous
    print '(A, 2(ES14.6, 1X))', 'source-contained max/rms B error: ', &
        metrics(1, 1), sqrt(metrics(2, 1)/metric_counts(1))
    print '(A, 2(ES14.6, 1X))', 'neighbor-recovered max/rms B error: ', &
        metrics(1, 2), sqrt(metrics(2, 2)/neighbor_recovered)
    print '(A)', 'PASS: global JOREK plane containment and interpolation'

contains

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

    subroutine compare_field(element, s, t, phi, rz, triangle, weight, group)
        integer, intent(in) :: element, triangle, group
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
        metrics(1, group) = max(metrics(1, group), error)
        metrics(2, group) = metrics(2, group) + error**2
        metric_counts(group) = metric_counts(group) + 1
    end subroutine compare_field

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
