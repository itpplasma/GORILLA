program test_jorek_owner_field_consistency
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use jorek_bezier, only: jorek_locator_t, build_jorek_locator
    use jorek_field_backend_mod, only: evaluate_jorek_model303_gorilla_element
    use jorek_refined_plane_mod, only: extract_refined_jorek_plane
    use jorek_restart, only: jorek_restart_t, load_jorek_restart

    implicit none

    integer, parameter :: levels(3) = [2, 4, 8]
    integer, parameter :: expected_samples(3) = [120928, 217800, 411544]
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: phases(4) = [0.0_dp, 0.5_dp*pi, pi, 1.5_dp*pi]
    type(jorek_restart_t) :: data
    type(jorek_locator_t) :: locator
    integer, allocatable :: element_vertex(:, :, :), foreign_vertices(:)
    integer, allocatable :: jumps_above_limit(:), triangles(:, :), vertex_element(:)
    real(dp), allocatable :: element_metrics(:, :), plane_rz(:, :), psi(:)
    real(dp), allocatable :: vertex_st(:, :)
    real(dp) :: a(3), b_covariant(3), b_reference(3), b_target(3)
    real(dp) :: bmod, error, metrics(2), r
    character(len=1024) :: filename, output_filename
    integer :: element, i, ierr, j, level, node, output_unit, phase, samples
    integer :: source_element
    real(dp) :: s, t

    if (command_argument_count() /= 2) &
        error stop 'restart path and output path are required'
    call get_command_argument(1, filename)
    call get_command_argument(2, output_filename)
    call load_jorek_restart(trim(filename), data, ierr)
    if (ierr /= 0) error stop 'JOREK restart load failed'
    call build_jorek_locator(data, locator, ierr)
    if (ierr /= 0) error stop 'JOREK locator build failed'
    allocate(element_metrics(2, data%n_elements), foreign_vertices(data%n_elements))
    allocate(jumps_above_limit(data%n_elements))
    open(newunit=output_unit, file=trim(output_filename), status='replace', &
        action='write', iostat=ierr)
    if (ierr /= 0) error stop 'cannot open owner-consistency output'
    write(output_unit, '(A)') 'refinement,element,foreign_vertices,' // &
        'samples,jumps_above_2pct,max_relative,rms_relative'

    do level = 1, size(levels)
        call extract_refined_jorek_plane(data, levels(level), plane_rz, psi, &
            triangles, vertex_element, vertex_st, ierr, element_vertex)
        if (ierr /= 0) error stop 'JOREK plane refinement failed'
        element_metrics = 0.0_dp
        foreign_vertices = 0
        jumps_above_limit = 0
        metrics = 0.0_dp
        samples = 0
        do element = 1, data%n_elements
            do j = 0, levels(level)
                t = real(j, dp)/levels(level)
                do i = 0, levels(level)
                    s = real(i, dp)/levels(level)
                    node = element_vertex(i, j, element)
                    source_element = vertex_element(node)
                    if (source_element == 0 .or. source_element == element) cycle
                    foreign_vertices(element) = foreign_vertices(element) + 1
                    r = 100.0_dp*plane_rz(1, node)
                    do phase = 1, size(phases)
                        call evaluate_field(element, s, t, phases(phase), r, &
                            b_target)
                        call evaluate_field(source_element, vertex_st(1, node), &
                            vertex_st(2, node), phases(phase), r, b_reference)
                        error = relative_error(b_target, b_reference)
                        call update_metrics(error, metrics)
                        call update_metrics(error, element_metrics(:, element))
                        if (error > 0.02_dp) &
                            jumps_above_limit(element) = &
                                jumps_above_limit(element) + 1
                        samples = samples + 1
                    end do
                end do
            end do
        end do
        print '(A, I0, A, I0, A, ES12.4, A, ES12.4, A, I0)', &
            'refinement=', levels(level), ' samples=', samples, &
            ' max=', metrics(1), ' rms=', sqrt(metrics(2)/samples), &
            ' above 2%=', sum(jumps_above_limit)
        if (samples /= expected_samples(level) &
                .or. samples /= size(phases)*sum(foreign_vertices) &
                .or. samples == 0 &
                .or. .not. all(ieee_is_finite(element_metrics))) &
            error stop 'owner-consistency sample partition failed'
        call write_element_metrics(output_unit, levels(level), &
            foreign_vertices, jumps_above_limit, element_metrics)
    end do
    close(output_unit)
    print '(A)', 'PASS: JOREK shared-vertex owner field comparison'

contains

    subroutine evaluate_field(element, s, t, phi, r, b)
        integer, intent(in) :: element
        real(dp), intent(in) :: s, t, phi, r
        real(dp), intent(out) :: b(3)

        call evaluate_jorek_model303_gorilla_element(data, element, s, t, &
            phi, r, a, b_covariant, bmod, ierr, 100.0_dp, 1.0e4_dp, locator)
        if (ierr /= 0 .or. .not. ieee_is_finite(bmod)) &
            error stop 'JOREK owner field evaluation failed'
        b = [b_covariant(1), b_covariant(3), b_covariant(2)/r]
    end subroutine evaluate_field

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

    subroutine write_element_metrics(unit, refinement, foreign, above, values)
        integer, intent(in) :: unit, refinement, foreign(:), above(:)
        real(dp), intent(in) :: values(:, :)

        integer :: element, samples
        real(dp) :: rms

        do element = 1, size(foreign)
            samples = size(phases)*foreign(element)
            rms = 0.0_dp
            if (samples > 0) rms = sqrt(values(2, element)/samples)
            write(unit, '(I0,",",I0,",",I0,",",I0,",",I0,2(",",ES24.16E3))') &
                refinement, element, foreign(element), samples, above(element), &
                values(1, element), rms
        end do
    end subroutine write_element_metrics

end program test_jorek_owner_field_consistency
