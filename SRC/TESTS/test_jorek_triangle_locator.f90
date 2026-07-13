program test_jorek_triangle_locator
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use jorek_triangle_locator_mod, only: build_jorek_triangle_locator, &
        jorek_triangle_locator_t, locate_jorek_triangle

    implicit none

    type(jorek_triangle_locator_t) :: locator
    real(dp), parameter :: vertices(2, 6) = reshape([ &
        0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, &
        0.0_dp, 1.0_dp, 2.0_dp, 0.0_dp, 2.0_dp, 1.0_dp], [2, 6])
    integer, parameter :: triangles(4, 3) = reshape([ &
        1, 1, 2, 2, 2, 3, 5, 6, 3, 4, 6, 3], [4, 3])
    real(dp) :: weight(3)
    integer :: ierr, matches, triangle

    call build_jorek_triangle_locator(vertices, triangles, locator, ierr)
    if (ierr /= 0) error stop 'triangle locator build failed'
    call expect_location([0.75_dp, 0.25_dp], 1)
    call expect_location([0.25_dp, 0.75_dp], 2)
    call expect_location([1.75_dp, 0.25_dp], 3)
    call expect_location([1.25_dp, 0.75_dp], 4)
    call expect_location([0.5_dp, 0.5_dp], 1)
    call locate_jorek_triangle(vertices, triangles, locator, [0.5_dp, 0.5_dp], &
        triangle, weight, ierr, matches)
    if (ierr /= 0 .or. triangle /= 1 .or. matches /= 2) &
        error stop 'shared-edge triangle multiplicity changed'
    call locate_jorek_triangle(vertices, triangles, locator, [-0.1_dp, 0.5_dp], &
        triangle, weight, ierr)
    if (ierr == 0 .or. triangle /= 0) error stop 'outside point was accepted'
    call build_jorek_triangle_locator(vertices, reshape([1, 2, 7], [1, 3]), &
        locator, ierr)
    if (ierr == 0) error stop 'invalid triangle index was accepted'
    print '(A)', 'PASS: deterministic indexed JOREK triangle location'

contains

    subroutine expect_location(point, expected)
        real(dp), intent(in) :: point(2)
        integer, intent(in) :: expected

        call locate_jorek_triangle(vertices, triangles, locator, point, &
            triangle, weight, ierr)
        if (ierr /= 0 .or. triangle /= expected &
                .or. abs(sum(weight) - 1.0_dp) > 1.0e-14_dp) &
            error stop 'triangle location mismatch'
    end subroutine expect_location

end program test_jorek_triangle_locator
