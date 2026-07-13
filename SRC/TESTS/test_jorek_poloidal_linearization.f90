program test_jorek_poloidal_linearization
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use jorek_bezier, only: jorek_locator_t, build_jorek_locator, &
        evaluate_jorek_geometry
    use jorek_model303_field, only: evaluate_jorek_model303_at, &
        evaluate_jorek_model303_b
    use jorek_restart, only: jorek_restart_t, load_jorek_restart

    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: corner_s(4) = [0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp]
    real(dp), parameter :: corner_t(4) = [0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp]
    real(dp), parameter :: weights(3, 4) = reshape([ &
        1.0_dp/3, 1.0_dp/3, 1.0_dp/3, 0.6_dp, 0.2_dp, 0.2_dp, &
        0.2_dp, 0.6_dp, 0.2_dp, 0.2_dp, 0.2_dp, 0.6_dp], [3, 4])
    integer, parameter :: triangle_vertices(3, 2) = reshape([1, 2, 3, 1, 3, 4], [3, 2])
    integer, parameter :: refinements(3) = [1, 2, 4]
    integer, parameter :: reference_refinement = 2

    type(jorek_restart_t) :: data
    type(jorek_locator_t) :: locator
    logical, allocatable :: owner_mismatch_elements(:), uncovered_elements(:)
    real(dp) :: corner_rz(2, 4), metrics(6), phi, rz_st(2, 2)
    real(dp) :: worst_b_exact(3), worst_b_interp(3), worst_rz(2)
    character(len=1024) :: filename
    integer :: axis_elements, element, i, ierr, j, level, outside, phase
    integer :: refinement, sample, samples, triangle, vertex
    real(dp) :: target_s, target_t
    integer :: worst_case(4), worst_found
    real(dp) :: worst_st(2)

    if (command_argument_count() /= 1) error stop 'restart path argument is required'
    call get_command_argument(1, filename)
    call load_jorek_restart(trim(filename), data, ierr)
    if (ierr /= 0) error stop 'JOREK restart load failed'
    call build_jorek_locator(data, locator, ierr)
    if (ierr /= 0) error stop 'JOREK locator build failed'
    allocate (owner_mismatch_elements(data%n_elements))
    allocate (uncovered_elements(data%n_elements))

    do level = 1, size(refinements)
        refinement = refinements(level)
        metrics = 0.0_dp
        worst_case = 0
        axis_elements = 0
        outside = 0
        samples = 0
        owner_mismatch_elements = .false.
        uncovered_elements = .false.
        do element = 1, data%n_elements
            do vertex = 1, 4
                call evaluate_jorek_geometry(data, element, corner_s(vertex), &
                    corner_t(vertex), corner_rz(:, vertex), rz_st, ierr)
                if (ierr /= 0) error stop 'JOREK corner geometry failed'
            end do
            if (has_coincident_corners(corner_rz)) then
                axis_elements = axis_elements + 1
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
                                call compare_fixed_sample(data, locator, phi, &
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
        if (samples + outside /= 769792) &
            error stop 'fixed poloidal sample count changed'
        if (outside == 0) error stop 'coverage-failure fixture changed'
        if (level == size(refinements) &
                .and. metrics(1) <= 2.0e-2_dp) &
            error stop 'maximum field-error rejection unexpectedly cleared'
        print '(A, I0, A, I0, A, I0, A, I0)', 'refinement=', refinement, &
            ' samples=', samples, ' outside=', outside, &
            ' axis_elements=', axis_elements
        print '(A, 6(ES14.6, 1X))', &
            'max/rms B, A, |B| relative errors: ', metrics
        print '(A, I0, A, I0)', 'uncovered elements=', &
            count(uncovered_elements), ' interior owner mismatches=', &
            count(owner_mismatch_elements)
        print '(A, 4(I0, 1X), A, I0, A, 2(ES14.6, 1X), A, 2(ES14.6, 1X))', &
            'worst B element/phase/i/j ', worst_case, &
            'found ', worst_found, ' st ', worst_st, ' RZ ', worst_rz
        print '(A, 3(ES14.6, 1X))', 'worst B exact ', worst_b_exact
        print '(A, 3(ES14.6, 1X))', 'worst B interpolated ', worst_b_interp
    end do
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

    subroutine compare_fixed_sample(data, locator, phi, target_s, target_t, &
            refinement, metrics, samples, uncovered)
        type(jorek_restart_t), intent(in) :: data
        type(jorek_locator_t), intent(in) :: locator
        real(dp), intent(in) :: phi, target_s, target_t
        integer, intent(in) :: refinement
        real(dp), intent(inout) :: metrics(6)
        integer, intent(inout) :: samples, uncovered

        real(dp) :: a_corner(3, 4), a_exact(3), a_interp(3), a_jorek(3)
        real(dp) :: barycentric(3), b_exact(3), b_interp(3), b_jorek(3)
        real(dp) :: bmod_corner(4), bmod_exact, bmod_interp
        real(dp) :: corner_rz(2, 4), h_corner(3, 4), h_exact(3), h_interp(3)
        real(dp) :: rz(2), rz_st(2, 2), st(2), sub_s(4), sub_t(4)
        integer :: cell_i, cell_j, corner_owner, found, ierr, index, vertices(3)

        cell_i = min(refinement - 1, int(target_s*refinement))
        cell_j = min(refinement - 1, int(target_t*refinement))
        sub_s = real([cell_i, cell_i + 1, cell_i + 1, cell_i], dp)/refinement
        sub_t = real([cell_j, cell_j, cell_j + 1, cell_j + 1], dp)/refinement
        do index = 1, 4
            call evaluate_jorek_geometry(data, element, sub_s(index), &
                sub_t(index), corner_rz(:, index), rz_st, ierr)
            if (ierr /= 0) error stop 'JOREK refined corner geometry failed'
            call evaluate_global_corner(data, locator, corner_rz(:, index), &
                phi, a_corner(:, index), h_corner(:, index), &
                bmod_corner(index), corner_owner)
            if (sub_s(index) > 0.0_dp .and. sub_s(index) < 1.0_dp &
                    .and. sub_t(index) > 0.0_dp &
                    .and. sub_t(index) < 1.0_dp &
                    .and. corner_owner /= element) &
                owner_mismatch_elements(element) = .true.
        end do
        call evaluate_jorek_geometry(data, element, target_s, target_t, rz, &
            rz_st, ierr)
        if (ierr /= 0) error stop 'JOREK fixed-point geometry failed'
        if ((target_t*refinement - cell_j) &
                <= (target_s*refinement - cell_i)) then
            vertices = triangle_vertices(:, 1)
        else
            vertices = triangle_vertices(:, 2)
        end if
        call physical_barycentric(corner_rz(:, vertices), rz, barycentric, ierr)
        if (ierr /= 0 .or. any(barycentric < -1.0e-10_dp) &
                .or. any(barycentric > 1.0_dp + 1.0e-10_dp)) then
            if (all(vertices == triangle_vertices(:, 1))) then
                vertices = triangle_vertices(:, 2)
            else
                vertices = triangle_vertices(:, 1)
            end if
            call physical_barycentric(corner_rz(:, vertices), rz, &
                barycentric, ierr)
            if (ierr /= 0 .or. any(barycentric < -1.0e-10_dp) &
                    .or. any(barycentric > 1.0_dp + 1.0e-10_dp)) then
                uncovered = uncovered + 1
                uncovered_elements(element) = .true.
                return
            end if
        end if
        call evaluate_jorek_model303_at(data, rz, phi, a_jorek, b_jorek, &
            found, st, ierr, locator)
        if (ierr /= 0) error stop 'JOREK fixed-point field evaluation failed'
        call convert_components(rz(1)*100.0_dp, a_jorek, b_jorek, a_exact, &
            h_exact, bmod_exact)
        a_interp = 0.0_dp
        h_interp = 0.0_dp
        bmod_interp = 0.0_dp
        do index = 1, 3
            a_interp = a_interp + barycentric(index)*a_corner(:, vertices(index))
            h_interp = h_interp + barycentric(index)*h_corner(:, vertices(index))
            bmod_interp = bmod_interp &
                + barycentric(index)*bmod_corner(vertices(index))
        end do
        a_interp = [a_interp(1), a_interp(2)/(rz(1)*100.0_dp), a_interp(3)]
        a_exact = [a_exact(1), a_exact(2)/(rz(1)*100.0_dp), a_exact(3)]
        b_interp = [h_interp(1)*bmod_interp, &
            h_interp(2)*bmod_interp/(rz(1)*100.0_dp), h_interp(3)*bmod_interp]
        b_exact = [h_exact(1)*bmod_exact, &
            h_exact(2)*bmod_exact/(rz(1)*100.0_dp), h_exact(3)*bmod_exact]
        if (relative_error(b_interp, b_exact) > metrics(1)) then
            worst_case = [element, phase, cell_i, cell_j]
            worst_found = found
            worst_st = st
            worst_rz = rz
            worst_b_exact = b_exact
            worst_b_interp = b_interp
        end if
        call update_metrics(relative_error(b_interp, b_exact), metrics(1:2))
        call update_metrics(relative_error(a_interp, a_exact), metrics(3:4))
        call update_metrics(abs(bmod_interp/bmod_exact - 1.0_dp), metrics(5:6))
        samples = samples + 1
    end subroutine compare_fixed_sample

    subroutine evaluate_global_corner(data, locator, rz, phi, a, h, bmod, &
            found)
        type(jorek_restart_t), intent(in) :: data
        type(jorek_locator_t), intent(in) :: locator
        real(dp), intent(in) :: rz(2), phi
        real(dp), intent(out) :: a(3), h(3), bmod
        integer, intent(out) :: found

        real(dp) :: a_jorek(3), b_jorek(3), st(2)
        integer :: ierr

        call evaluate_jorek_model303_at(data, rz, phi, a_jorek, b_jorek, &
            found, st, ierr, locator)
        if (ierr /= 0) error stop 'JOREK global corner field failed'
        call convert_components(rz(1)*100.0_dp, a_jorek, b_jorek, a, h, bmod)
    end subroutine evaluate_global_corner

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
