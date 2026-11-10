module hamiltonian_time_coefficients_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private
    public :: hamiltonian_time_coefficients
contains
    pure subroutine hamiltonian_time_coefficients(h, gh, ca, ch, ce, cmq, &
                                                 strong_e, a, b, da, db)
        real(dp), intent(in) :: h(3), gh(3, 3), ca(3), ch(3), ce(3), cmq
        logical, intent(in) :: strong_e
        real(dp), intent(out) :: a, b, da(3), db(3)
        real(dp) :: symplectic_curl(3)

        symplectic_curl = ca
        if (strong_e) symplectic_curl = symplectic_curl + cmq*ce
        a = dot_product(h, symplectic_curl)
        b = dot_product(h, ch)
        da = matmul(gh, symplectic_curl)
        db = matmul(gh, ch)
    end subroutine hamiltonian_time_coefficients
end module hamiltonian_time_coefficients_mod
