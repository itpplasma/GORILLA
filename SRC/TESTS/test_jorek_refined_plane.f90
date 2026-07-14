program test_jorek_refined_plane
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use jorek_bezier, only: jorek_locator_t, build_jorek_locator, &
        evaluate_jorek_geometry
    use jorek_field_backend_mod, only: evaluate_jorek_model303_gorilla, &
        evaluate_jorek_model303_gorilla_element, jorek_chart_requires_global
    use jorek_refined_plane_mod, only: extract_refined_jorek_plane
    use jorek_restart, only: jorek_restart_t, load_jorek_restart

    implicit none

    integer, parameter :: levels(3) = [2, 4, 8]
    integer, parameter :: expected_nodes(3) = [24585, 98201, 392529]
    integer, parameter :: expected_triangles(3) = [48802, 195668, 783592]
    type(jorek_restart_t) :: data
    type(jorek_locator_t) :: locator
    real(dp), allocatable :: plane_rz(:, :), psi(:), vertex_st(:, :)
    integer, allocatable :: triangles(:, :), vertex_element(:)
    real(dp) :: a(3), b(3), bmod, max_bmod, rebuilt_rz(2), rz_st(2, 2)
    character(len=1024) :: filename
    integer :: fallback_count, ierr, level, node
    logical :: use_global

    if (command_argument_count() /= 1) error stop 'restart path argument is required'
    call get_command_argument(1, filename)
    call load_jorek_restart(trim(filename), data, ierr)
    if (ierr /= 0) error stop 'JOREK restart load failed'
    call build_jorek_locator(data, locator, ierr)
    if (ierr /= 0) error stop 'JOREK locator build failed'
    do level = 1, size(levels)
        call extract_refined_jorek_plane(data, levels(level), plane_rz, psi, &
            triangles, vertex_element, vertex_st, ierr)
        if (ierr /= 0) then
            write (*, '(A, I0, A, I0)') 'refinement=', levels(level), &
                ' error=', ierr
            error stop 'JOREK plane refinement failed'
        end if
        print '(A, I0, A, I0, A, I0)', 'refinement=', levels(level), &
            ' nodes=', size(plane_rz, 2), ' triangles=', size(triangles, 1)
        if (size(plane_rz, 2) /= expected_nodes(level) &
                .or. size(triangles, 1) /= expected_triangles(level)) &
            error stop 'JOREK refined plane dimensions changed'
        if (count(vertex_element == 0) /= 1 &
                .or. any(vertex_element < 0) &
                .or. any(vertex_element > data%n_elements)) &
            error stop 'JOREK refined plane owner is invalid'
        if (.not. all(ieee_is_finite(plane_rz)) &
                .or. .not. all(ieee_is_finite(psi)) &
                .or. any(vertex_st < 0.0_dp) .or. any(vertex_st > 1.0_dp)) &
            error stop 'JOREK refined plane metadata is invalid'
        do node = 1, size(vertex_element), 997
            if (vertex_element(node) == 0) cycle
            call evaluate_jorek_geometry(data, vertex_element(node), &
                vertex_st(1, node), vertex_st(2, node), rebuilt_rz, rz_st, ierr)
            if (ierr /= 0 .or. maxval(abs(rebuilt_rz - plane_rz(:, node))) &
                    > 1.0e-10_dp) &
                error stop 'JOREK refined vertex owner does not reproduce geometry'
        end do
        if (levels(level) == 8) then
            fallback_count = 0
            max_bmod = 0.0_dp
            do node = 1, size(vertex_element)
                if (vertex_element(node) == 0) then
                    call evaluate_jorek_model303_gorilla(data, &
                        100.0_dp*plane_rz(1, node), 0.0_dp, &
                        100.0_dp*plane_rz(2, node), a, b, bmod, ierr, &
                        locator, 100.0_dp, 1.0e4_dp)
                else
                    call jorek_chart_requires_global(data, vertex_element(node), &
                        vertex_st(1, node), vertex_st(2, node), use_global, ierr)
                    if (use_global) fallback_count = fallback_count + 1
                    call evaluate_jorek_model303_gorilla_element(data, &
                        vertex_element(node), vertex_st(1, node), &
                        vertex_st(2, node), 0.0_dp, &
                        100.0_dp*plane_rz(1, node), a, b, bmod, ierr, &
                        100.0_dp, 1.0e4_dp, locator)
                end if
                if (ierr /= 0 .or. .not. ieee_is_finite(bmod)) &
                    error stop 'JOREK refined vertex field is invalid'
                max_bmod = max(max_bmod, bmod)
            end do
            print '(A, I0, A, ES12.4)', 'chart fallbacks=', fallback_count, &
                ' max |B| [G]=', max_bmod
            if (fallback_count /= 0 .or. max_bmod > 3.0e4_dp) &
                error stop 'JOREK refined chart fallback did not bound the field'
        end if
    end do
    print '(A)', 'PASS: owner-labelled JOREK poloidal refinement'
end program test_jorek_refined_plane
