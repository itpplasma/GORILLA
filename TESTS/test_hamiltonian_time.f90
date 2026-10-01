program test_hamiltonian_time
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use hamiltonian_time_coefficients_mod, only: hamiltonian_time_coefficients
    implicit none
    real(dp) :: h(3), gh(3, 3), ca(3), ch(3), ce(3), x(3), v, cmq
    real(dp) :: a, b, da(3), db(3), hx(3), bs(3), omega(4, 4), jac, expected
    integer :: i, j, mode, sign_case
    logical :: strong_e

    h = [0.2_dp, 1.4_dp, -0.3_dp]
    gh = reshape([0.1_dp, -0.2_dp, 0.3_dp, 0.05_dp, &
        0.01_dp, -0.03_dp, 0.2_dp, 0.1_dp, 0.04_dp], [3, 3])
    ca = [2.0_dp, 4.0_dp, -1.0_dp]
    ch = [0.3_dp, -0.2_dp, 0.1_dp]
    ce = [0.4_dp, 1.1_dp, -0.2_dp]
    do sign_case = -1, 1, 2
        cmq = sign_case*0.7_dp
        do mode = 0, 1
            strong_e = mode == 1
            call hamiltonian_time_coefficients(h, gh, ca, ch, ce, cmq, &
                                               strong_e, a, b, da, db)
            do i = 1, 41
                x(1) = 0.01_dp*i
                x(2) = -0.02_dp*i
                x(3) = 0.003_dp*i
                v = -2.0_dp + 0.1_dp*i
                do j = 1, 3
                    hx(j) = h(j) + dot_product(gh(:, j), x)
                end do
                bs = ca + cmq*v*ch
                if (strong_e) bs = bs + cmq*ce
                ! Independent Pfaffian of the local 4-D symplectic form.
                omega = 0.0_dp
                omega(1, 2) = bs(3)
                omega(1, 3) = -bs(2)
                omega(2, 3) = bs(1)
                omega(1:3, 4) = hx
                do j = 1, 4
                    omega(j, 1:j-1) = -omega(1:j-1, j)
                end do
                expected = omega(1, 2)*omega(3, 4) &
                    - omega(1, 3)*omega(2, 4) + omega(1, 4)*omega(2, 3)
                jac = a + dot_product(da, x) + cmq*v*(b + dot_product(db, x))
                if (abs(jac-expected) > 2.0e-14_dp*abs(expected)) error stop 1
            end do
        end do
    end do
    print *, 'PASS: Hamiltonian clock matches symplectic Pfaffian, both signs'
end program test_hamiltonian_time
