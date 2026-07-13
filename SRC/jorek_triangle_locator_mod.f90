module jorek_triangle_locator_mod
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none
    private

    type, public :: jorek_triangle_locator_t
        real(dp) :: lower(2) = 0.0_dp
        real(dp) :: scale(2) = 0.0_dp
        integer :: n_bins(2) = 0
        integer, allocatable :: offsets(:), entries(:)
    end type jorek_triangle_locator_t

    public :: build_jorek_triangle_locator, locate_jorek_triangle

contains

    subroutine build_jorek_triangle_locator(vertices, triangles, locator, ierr)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :)
        type(jorek_triangle_locator_t), intent(out) :: locator
        integer, intent(out) :: ierr

        integer, allocatable :: counts(:), cursor(:)
        real(dp) :: extent(2), upper(2)
        integer :: bin, first(2), last(2), target_bins, triangle, x, y

        locator = jorek_triangle_locator_t()
        ierr = 1
        if (size(vertices, 1) /= 2 .or. size(vertices, 2) < 3 &
                .or. size(triangles, 2) /= 3 .or. size(triangles, 1) < 1) return
        if (.not. all(ieee_is_finite(vertices))) return
        if (any(triangles < 1) .or. any(triangles > size(vertices, 2))) return
        locator%lower = minval(vertices, dim=2)
        upper = maxval(vertices, dim=2)
        extent = upper - locator%lower
        if (any(extent <= 0.0_dp)) return
        target_bins = max(1, size(triangles, 1)/4)
        locator%n_bins(1) = max(1, min(2048, &
            nint(sqrt(real(target_bins, dp)*extent(1)/extent(2)))))
        locator%n_bins(2) = max(1, min(2048, &
            (target_bins + locator%n_bins(1) - 1)/locator%n_bins(1)))
        locator%scale = real(locator%n_bins, dp)/extent
        allocate(counts(product(locator%n_bins)), source=0)
        do triangle = 1, size(triangles, 1)
            call triangle_bins(vertices(:, triangles(triangle, :)), locator, &
                first, last)
            do y = first(2), last(2)
                do x = first(1), last(1)
                    bin = linear_bin(locator, x, y)
                    counts(bin) = counts(bin) + 1
                end do
            end do
        end do
        allocate(locator%offsets(size(counts) + 1))
        locator%offsets(1) = 1
        do bin = 1, size(counts)
            locator%offsets(bin + 1) = locator%offsets(bin) + counts(bin)
        end do
        allocate(locator%entries(locator%offsets(size(counts) + 1) - 1))
        allocate(cursor, source=locator%offsets(1:size(counts)))
        do triangle = 1, size(triangles, 1)
            call triangle_bins(vertices(:, triangles(triangle, :)), locator, &
                first, last)
            do y = first(2), last(2)
                do x = first(1), last(1)
                    bin = linear_bin(locator, x, y)
                    locator%entries(cursor(bin)) = triangle
                    cursor(bin) = cursor(bin) + 1
                end do
            end do
        end do
        ierr = 0
    end subroutine build_jorek_triangle_locator

    subroutine locate_jorek_triangle(vertices, triangles, locator, point, &
            triangle, barycentric, ierr, match_count)
        real(dp), intent(in) :: vertices(:, :), point(2)
        integer, intent(in) :: triangles(:, :)
        type(jorek_triangle_locator_t), intent(in) :: locator
        integer, intent(out) :: triangle, ierr
        real(dp), intent(out) :: barycentric(3)
        integer, intent(out), optional :: match_count

        real(dp) :: trial_weight(3)
        integer :: bin, entry, matches, trial, x, y

        triangle = 0
        barycentric = 0.0_dp
        ierr = 1
        matches = 0
        if (present(match_count)) match_count = 0
        if (.not. locator_ready(locator) .or. .not. all(ieee_is_finite(point))) &
            return
        if (any(point < locator%lower) &
                .or. any(point > locator%lower &
                    + real(locator%n_bins, dp)/locator%scale)) return
        call point_bin(locator, point, x, y)
        bin = linear_bin(locator, x, y)
        do entry = locator%offsets(bin), locator%offsets(bin + 1) - 1
            trial = locator%entries(entry)
            call triangle_barycentric(vertices(:, triangles(trial, :)), &
                point, trial_weight, ierr)
            if (ierr == 0) then
                if (all(trial_weight >= -1.0e-10_dp) &
                        .and. all(trial_weight <= 1.0_dp + 1.0e-10_dp)) then
                    matches = matches + 1
                    if (triangle == 0) then
                        triangle = trial
                        barycentric = trial_weight
                    end if
                end if
            end if
        end do
        if (present(match_count)) match_count = matches
        if (matches > 0) then
            ierr = 0
        else
            barycentric = 0.0_dp
            ierr = 1
        end if
    end subroutine locate_jorek_triangle

    subroutine triangle_bins(vertices, locator, first, last)
        real(dp), intent(in) :: vertices(2, 3)
        type(jorek_triangle_locator_t), intent(in) :: locator
        integer, intent(out) :: first(2), last(2)

        call point_bin(locator, minval(vertices, dim=2), first(1), first(2))
        call point_bin(locator, maxval(vertices, dim=2), last(1), last(2))
    end subroutine triangle_bins

    subroutine point_bin(locator, point, x, y)
        type(jorek_triangle_locator_t), intent(in) :: locator
        real(dp), intent(in) :: point(2)
        integer, intent(out) :: x, y

        x = min(locator%n_bins(1), max(1, &
            1 + int((point(1) - locator%lower(1))*locator%scale(1))))
        y = min(locator%n_bins(2), max(1, &
            1 + int((point(2) - locator%lower(2))*locator%scale(2))))
    end subroutine point_bin

    integer function linear_bin(locator, x, y)
        type(jorek_triangle_locator_t), intent(in) :: locator
        integer, intent(in) :: x, y

        linear_bin = x + locator%n_bins(1)*(y - 1)
    end function linear_bin

    logical function locator_ready(locator)
        type(jorek_triangle_locator_t), intent(in) :: locator

        locator_ready = all(locator%n_bins > 0) .and. allocated(locator%offsets) &
            .and. allocated(locator%entries)
    end function locator_ready

    subroutine triangle_barycentric(vertices, point, weight, ierr)
        real(dp), intent(in) :: vertices(2, 3), point(2)
        real(dp), intent(out) :: weight(3)
        integer, intent(out) :: ierr

        real(dp) :: determinant, scale

        determinant = (vertices(1, 2) - vertices(1, 1)) &
            *(vertices(2, 3) - vertices(2, 1)) &
            - (vertices(1, 3) - vertices(1, 1)) &
            *(vertices(2, 2) - vertices(2, 1))
        scale = max(1.0_dp, maxval(abs(vertices)))
        weight = 0.0_dp
        ierr = 1
        if (abs(determinant) <= 64.0_dp*epsilon(1.0_dp)*scale**2) return
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
    end subroutine triangle_barycentric

end module jorek_triangle_locator_mod
