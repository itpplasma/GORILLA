module jorek_field_backend_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use jorek_model303_field, only: evaluate_jorek_model303_at
    use jorek_restart, only: jorek_restart_t, load_jorek_restart, &
        free_jorek_restart

    implicit none
    private

    type(jorek_restart_t), save :: field_data
    logical, save :: field_loaded = .false.

    public :: load_jorek_field_backend, free_jorek_field_backend
    public :: evaluate_jorek_field_backend, evaluate_jorek_model303_gorilla

contains

    subroutine load_jorek_field_backend(filename, ierr)
        character(len=*), intent(in) :: filename
        integer, intent(out) :: ierr

        if (field_loaded) call free_jorek_restart(field_data)
        call load_jorek_restart(filename, field_data, ierr)
        field_loaded = ierr == 0
    end subroutine load_jorek_field_backend

    subroutine free_jorek_field_backend()
        if (field_loaded) call free_jorek_restart(field_data)
        field_loaded = .false.
    end subroutine free_jorek_field_backend

    pure subroutine evaluate_jorek_model303_gorilla(data, r, phi, z, &
            a_covariant, b_covariant, bmod, ierr)
        type(jorek_restart_t), intent(in) :: data
        real(dp), intent(in) :: r, phi, z
        real(dp), intent(out) :: a_covariant(3), b_covariant(3), bmod
        integer, intent(out) :: ierr

        real(dp) :: b_physical(3), st(2)
        integer :: element

        call evaluate_jorek_model303_at(data, [r, z], phi, a_covariant, &
            b_physical, element, st, ierr)
        if (ierr /= 0) then
            a_covariant = 0.0_dp
            b_covariant = 0.0_dp
            bmod = 0.0_dp
            return
        end if
        b_covariant = [b_physical(1), r*b_physical(3), b_physical(2)]
        bmod = sqrt(sum(b_physical**2))
    end subroutine evaluate_jorek_model303_gorilla

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
            a_covariant, b_covariant, bmod, ierr)
    end subroutine evaluate_jorek_field_backend

end module jorek_field_backend_mod
