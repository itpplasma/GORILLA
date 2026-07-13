module jorek_field_backend_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use jorek_bezier, only: jorek_locator_t, build_jorek_locator
    use jorek_model303_field, only: evaluate_jorek_model303_a, &
        evaluate_jorek_model303_b, evaluate_jorek_model303_at
    use jorek_restart, only: jorek_restart_t, load_jorek_restart, &
        free_jorek_restart

    implicit none
    private

    type(jorek_restart_t), save :: field_data
    type(jorek_locator_t), save :: field_locator
    logical, save :: field_loaded = .false.
    real(dp), save :: length_scale = 1.0_dp
    real(dp), save :: magnetic_field_scale = 1.0_dp

    public :: load_jorek_field_backend, free_jorek_field_backend
    public :: evaluate_jorek_field_backend, evaluate_jorek_field_backend_element
    public :: evaluate_jorek_model303_gorilla, &
        evaluate_jorek_model303_gorilla_element

contains

    subroutine load_jorek_field_backend(filename, ierr, length_scale_in, &
            magnetic_field_scale_in)
        character(len=*), intent(in) :: filename
        integer, intent(out) :: ierr
        real(dp), intent(in), optional :: length_scale_in
        real(dp), intent(in), optional :: magnetic_field_scale_in

        if (field_loaded) call free_jorek_restart(field_data)
        field_locator = jorek_locator_t()
        length_scale = 1.0_dp
        magnetic_field_scale = 1.0_dp
        if (present(length_scale_in)) length_scale = length_scale_in
        if (present(magnetic_field_scale_in)) &
            magnetic_field_scale = magnetic_field_scale_in
        if (.not. ieee_is_finite(length_scale) .or. length_scale <= 0.0_dp &
                .or. .not. ieee_is_finite(magnetic_field_scale) &
                .or. magnetic_field_scale <= 0.0_dp) then
            ierr = 10
            field_loaded = .false.
            return
        end if
        call load_jorek_restart(filename, field_data, ierr)
        if (ierr == 0) call build_jorek_locator(field_data, field_locator, ierr)
        field_loaded = ierr == 0
        if (.not. field_loaded) call free_jorek_restart(field_data)
    end subroutine load_jorek_field_backend

    subroutine free_jorek_field_backend()
        if (field_loaded) call free_jorek_restart(field_data)
        field_locator = jorek_locator_t()
        field_loaded = .false.
    end subroutine free_jorek_field_backend

    pure subroutine evaluate_jorek_model303_gorilla(data, r, phi, z, &
            a_covariant, b_covariant, bmod, ierr, locator, length_scale_in, &
            magnetic_field_scale_in)
        type(jorek_restart_t), intent(in) :: data
        real(dp), intent(in) :: r, phi, z
        real(dp), intent(out) :: a_covariant(3), b_covariant(3), bmod
        integer, intent(out) :: ierr
        type(jorek_locator_t), intent(in), optional :: locator
        real(dp), intent(in), optional :: length_scale_in
        real(dp), intent(in), optional :: magnetic_field_scale_in

        real(dp) :: a_jorek(3), b_physical(3), st(2)
        real(dp) :: length_conversion, field_conversion
        integer :: element

        length_conversion = 1.0_dp
        field_conversion = 1.0_dp
        if (present(length_scale_in)) length_conversion = length_scale_in
        if (present(magnetic_field_scale_in)) &
            field_conversion = magnetic_field_scale_in
        if (.not. ieee_is_finite(length_conversion) &
                .or. length_conversion <= 0.0_dp &
                .or. .not. ieee_is_finite(field_conversion) &
                .or. field_conversion <= 0.0_dp) then
            a_covariant = 0.0_dp
            b_covariant = 0.0_dp
            bmod = 0.0_dp
            ierr = 10
            return
        end if

        if (present(locator)) then
            call evaluate_jorek_model303_at(data, &
                [r, z]/length_conversion, phi, a_jorek, &
                b_physical, element, st, ierr, locator)
        else
            call evaluate_jorek_model303_at(data, &
                [r, z]/length_conversion, phi, a_jorek, &
                b_physical, element, st, ierr)
        end if
        if (ierr /= 0) then
            a_covariant = 0.0_dp
            b_covariant = 0.0_dp
            bmod = 0.0_dp
            return
        end if
        a_covariant = a_jorek*field_conversion*length_conversion
        a_covariant(2) = a_covariant(2)*length_conversion
        b_physical = b_physical*field_conversion
        b_covariant = [b_physical(1), r*b_physical(3), b_physical(2)]
        bmod = sqrt(sum(b_physical**2))
    end subroutine evaluate_jorek_model303_gorilla

    pure subroutine evaluate_jorek_model303_gorilla_element(data, element, &
            s, t, phi, r, a_covariant, b_covariant, bmod, ierr, &
            length_scale_in, magnetic_field_scale_in)
        type(jorek_restart_t), intent(in) :: data
        integer, intent(in) :: element
        real(dp), intent(in) :: s, t, phi, r
        real(dp), intent(out) :: a_covariant(3), b_covariant(3), bmod
        integer, intent(out) :: ierr
        real(dp), intent(in), optional :: length_scale_in
        real(dp), intent(in), optional :: magnetic_field_scale_in

        real(dp) :: a_jorek(3), b_physical(3)
        real(dp) :: length_conversion, field_conversion

        length_conversion = 1.0_dp
        field_conversion = 1.0_dp
        if (present(length_scale_in)) length_conversion = length_scale_in
        if (present(magnetic_field_scale_in)) &
            field_conversion = magnetic_field_scale_in
        if (.not. ieee_is_finite(length_conversion) &
                .or. length_conversion <= 0.0_dp &
                .or. .not. ieee_is_finite(field_conversion) &
                .or. field_conversion <= 0.0_dp &
                .or. .not. ieee_is_finite(s) &
                .or. .not. ieee_is_finite(t) &
                .or. .not. ieee_is_finite(phi) &
                .or. .not. ieee_is_finite(r) .or. r <= 0.0_dp &
                .or. element < 1 .or. element > data%n_elements &
                .or. s < 0.0_dp .or. s > 1.0_dp &
                .or. t < 0.0_dp .or. t > 1.0_dp) then
            a_covariant = 0.0_dp
            b_covariant = 0.0_dp
            bmod = 0.0_dp
            ierr = 10
            return
        end if
        call evaluate_jorek_model303_a(data, element, s, t, phi, a_jorek, ierr)
        if (ierr == 0) call evaluate_jorek_model303_b(data, element, s, t, &
            phi, b_physical, ierr)
        if (ierr /= 0) then
            a_covariant = 0.0_dp
            b_covariant = 0.0_dp
            bmod = 0.0_dp
            return
        end if
        a_covariant = a_jorek*field_conversion*length_conversion
        a_covariant(2) = a_covariant(2)*length_conversion
        b_physical = b_physical*field_conversion
        b_covariant = [b_physical(1), r*b_physical(3), b_physical(2)]
        bmod = sqrt(sum(b_physical**2))
    end subroutine evaluate_jorek_model303_gorilla_element

    subroutine evaluate_jorek_field_backend(r, phi, z, a_covariant, &
            b_covariant, bmod, ierr)
        real(dp), intent(in) :: r, phi, z
        real(dp), intent(out) :: a_covariant(3), b_covariant(3), bmod
        integer, intent(out) :: ierr

        if (.not. field_loaded) then
            a_covariant = 0.0_dp
            b_covariant = 0.0_dp
            bmod = 0.0_dp
            ierr = 9
            return
        end if
        call evaluate_jorek_model303_gorilla(field_data, r, phi, z, &
            a_covariant, b_covariant, bmod, ierr, field_locator, &
            length_scale, magnetic_field_scale)
    end subroutine evaluate_jorek_field_backend

    subroutine evaluate_jorek_field_backend_element(r, phi, z, element, st, &
            a_covariant, b_covariant, bmod, ierr)
        real(dp), intent(in) :: r, phi, z, st(2)
        integer, intent(in) :: element
        real(dp), intent(out) :: a_covariant(3), b_covariant(3), bmod
        integer, intent(out) :: ierr

        if (element == 0) then
            call evaluate_jorek_field_backend(r, phi, z, a_covariant, &
                b_covariant, bmod, ierr)
            return
        end if
        if (.not. field_loaded) then
            a_covariant = 0.0_dp
            b_covariant = 0.0_dp
            bmod = 0.0_dp
            ierr = 9
            return
        end if
        call evaluate_jorek_model303_gorilla_element(field_data, element, &
            st(1), st(2), phi, r, a_covariant, b_covariant, bmod, ierr, &
            length_scale, magnetic_field_scale)
    end subroutine evaluate_jorek_field_backend_element

end module jorek_field_backend_mod
