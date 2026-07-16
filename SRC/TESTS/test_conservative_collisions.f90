program test_conservative_collisions
    ! Deterministic tests of the conservative BKP collision kernel:
    ! analytic nu(v) values, exact ensemble conservation on a fixed seed,
    ! beam isotropization, thermalization to the background Maxwellian
    ! (H-theorem sanity), and Maxwellian stationarity.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use constants, only: ame, echarge, ev2erg
    use xorshift_rng_mod, only: xorshift_rng_t, rng_seed, rng_uniform, &
        rng_normal, rng_sign
    use conservative_collisions_mod, only: chandrasekhar_g, nu_d_hat, &
        nu_e_hat, dln_nu_e_dxsq, nu_ref_self, scatter_pitch_bkp, &
        collide_marker, apply_collisions_conservative, collision_correction_t

    implicit none

    integer, parameter :: n_markers = 4096

    call test_rng_moments()
    call test_frequencies_analytic()
    call test_conservation_fixed_seed()
    call test_isotropization()
    call test_thermalization()
    call test_maxwellian_stationarity()
    print '(A)', 'PASS: conservative BKP collision kernel'

contains

    subroutine expect_close(val, ref, tol, label)
        real(dp), intent(in) :: val, ref, tol
        character(len=*), intent(in) :: label

        if (abs(val - ref) > tol*max(abs(ref), 1.0_dp)) then
            print '(A,ES23.15,A,ES23.15)', label//': got ', val, ' want ', ref
            error stop 'analytic value mismatch'
        end if
    end subroutine expect_close

    subroutine test_rng_moments()
        type(xorshift_rng_t) :: rng
        real(dp) :: mean_u, mean_g, var_g, g
        integer :: k

        call rng_seed(rng, 12345_int64)
        mean_u = 0.0_dp
        do k = 1, 100000
            mean_u = mean_u + rng_uniform(rng)
        end do
        mean_u = mean_u/100000.0_dp
        if (abs(mean_u - 0.5_dp) > 0.01_dp) error stop 'uniform mean off'
        mean_g = 0.0_dp
        var_g = 0.0_dp
        do k = 1, 100000
            g = rng_normal(rng)
            mean_g = mean_g + g
            var_g = var_g + g**2
        end do
        mean_g = mean_g/100000.0_dp
        var_g = var_g/100000.0_dp - mean_g**2
        if (abs(mean_g) > 0.02_dp) error stop 'normal mean off'
        if (abs(var_g - 1.0_dp) > 0.03_dp) error stop 'normal variance off'
    end subroutine test_rng_moments

    subroutine test_frequencies_analytic()
        ! Reference values computed with 30-digit arithmetic (mpmath) from
        ! the same definitions; the Chandrasekhar/deflection conventions
        ! match NEO-2 COMMON/collop_compute.f90.
        real(dp) :: nu

        call expect_close(chandrasekhar_g(0.5_dp), &
                          0.16221717669064828_dp, 5.0e-13_dp, 'G(0.5)')
        call expect_close(chandrasekhar_g(1.0_dp), &
                          0.21379664776456008_dp, 5.0e-13_dp, 'G(1)')
        call expect_close(chandrasekhar_g(2.0_dp), &
                          0.11924853678884608_dp, 5.0e-13_dp, 'G(2)')
        call expect_close(chandrasekhar_g(3.0_dp), &
                          0.055531119463731176_dp, 5.0e-13_dp, 'G(3)')
        call expect_close(nu_d_hat(0.5_dp), &
                          3.810237319903349_dp, 5.0e-13_dp, 'nu_d_hat(0.5)')
        call expect_close(nu_d_hat(1.0_dp), &
                          0.83602768048790201_dp, 5.0e-13_dp, 'nu_d_hat(1)')
        call expect_close(nu_d_hat(2.0_dp), &
                          0.14557502374486922_dp, 5.0e-13_dp, 'nu_d_hat(2)')
        call expect_close(nu_d_hat(3.0_dp), &
                          0.046499676388346053_dp, 5.0e-13_dp, 'nu_d_hat(3)')
        call expect_close(nu_e_hat(0.5_dp), &
                          0.86256737852508018_dp, 5.0e-13_dp, 'nu_e_hat(0.5)')
        call expect_close(nu_e_hat(1.0_dp), &
                          0.56841703746147706_dp, 5.0e-13_dp, 'nu_e_hat(1)')
        call expect_close(nu_e_hat(2.0_dp), &
                          0.15852189618467875_dp, 5.0e-13_dp, 'nu_e_hat(2)')
        call expect_close(nu_e_hat(3.0_dp), &
                          0.049213173269292289_dp, 5.0e-13_dp, 'nu_e_hat(3)')
        call expect_close(dln_nu_e_dxsq(0.5_dp), &
                          -0.58267862341543231_dp, 5.0e-12_dp, 'dlnnuE(0.5)')
        call expect_close(dln_nu_e_dxsq(1.0_dp), &
                          -0.52920017277884365_dp, 5.0e-12_dp, 'dlnnuE(1)')
        call expect_close(dln_nu_e_dxsq(2.0_dp), &
                          -0.33167245504511477_dp, 5.0e-12_dp, 'dlnnuE(2)')
        ! Small-x limits: nu_E -> 1, drag log-derivative -> -3/5.
        call expect_close(nu_e_hat(1.0e-3_dp), 1.0_dp, 1.0e-5_dp, &
                          'nu_e_hat(0+)')
        call expect_close(dln_nu_e_dxsq(1.0e-3_dp), -0.6_dp, 1.0e-4_dp, &
                          'dlnnuE(0+)')
        ! Large-x scaling nu_D ~ 3 sqrt(pi)/4 x^-3 (1 percent at x = 8).
        call expect_close(nu_d_hat(8.0_dp)*8.0_dp**3, &
                          1.329340388179137_dp, 1.0e-2_dp, 'nu_d_hat tail')
        ! Reference frequency 1/tau_aa, electrons at n = 1e13 cm^-3,
        ! T = 1 keV, ln Lambda = 15, with GORILLA's CGS constants.
        nu = nu_ref_self(1.0e13_dp, 1.0e3_dp*ev2erg, ame, echarge, 15.0_dp)
        call expect_close(nu, 13785.503557619853_dp, 1.0e-12_dp, 'nu_ref')
        ! First-principles pairing check: nu_ref_self*nu_d_hat(x) must be
        ! the physical deflection frequency nu_D = nu_perp/2 of Trubnikov/
        ! NRL, references computed with 40-digit mpmath from the
        ! numerically integrated Maxwell integral psi(u) -- independent of
        ! the erf-G closed form inside nu_d_hat, so it pins the
        ! normalization pairing, not just the formula.
        call expect_close(nu*nu_d_hat(0.5_dp), 52526.040128903552_dp, &
                          1.0e-12_dp, 'nu_D(0.5 v_T)')
        call expect_close(nu*nu_d_hat(1.0_dp), 11525.062563634647_dp, &
                          1.0e-12_dp, 'nu_D(v_T)')
        call expect_close(nu*nu_d_hat(2.0_dp), 2006.8250077354892_dp, &
                          1.0e-12_dp, 'nu_D(2 v_T)')
    end subroutine test_frequencies_analytic

    subroutine sample_maxwellian(rng, v_thermal, v, lam)
        type(xorshift_rng_t), intent(inout) :: rng
        real(dp), intent(in) :: v_thermal
        real(dp), intent(out) :: v(:), lam(:)

        real(dp) :: sigma, vx, vy, vz
        integer :: i

        sigma = v_thermal/sqrt(2.0_dp)
        do i = 1, size(v)
            vx = sigma*rng_normal(rng)
            vy = sigma*rng_normal(rng)
            vz = sigma*rng_normal(rng)
            v(i) = sqrt(vx**2 + vy**2 + vz**2)
            if (v(i) > 0.0_dp) then
                lam(i) = vz/v(i)
            else
                lam(i) = 0.0_dp
            end if
        end do
    end subroutine sample_maxwellian

    subroutine test_conservation_fixed_seed()
        ! 200 conservative steps on a weighted Maxwellian ensemble: the
        ! momentum and energy moments must return to their initial values
        ! to accumulated roundoff, while the raw test-particle defects are
        ! nonzero (the operator actually scatters).
        type(xorshift_rng_t) :: rng
        type(collision_correction_t) :: corr
        real(dp) :: w(n_markers), v(n_markers), lam(n_markers)
        real(dp), parameter :: v_thermal = 1.3e8_dp
        real(dp), parameter :: nu_ref = 1.0e4_dp, dt = 5.0e-6_dp
        real(dp) :: p_init, k_init, wsum, max_defect
        integer :: step, i, ierr

        call rng_seed(rng, 20260716_int64)
        call sample_maxwellian(rng, v_thermal, v, lam)
        do i = 1, n_markers
            w(i) = 1.0_dp + 0.5_dp*rng_uniform(rng)
        end do
        wsum = sum(w)
        p_init = sum(w*v*lam)
        k_init = sum(w*v**2)
        max_defect = 0.0_dp
        do step = 1, 200
            call apply_collisions_conservative(w, v, lam, dt, nu_ref, &
                                               v_thermal, rng, ierr, corr)
            if (ierr /= 0) error stop 'conservative step failed'
            max_defect = max(max_defect, abs(corr%momentum_defect))
        end do
        if (abs(sum(w*v*lam) - p_init)/(wsum*v_thermal) > 1.0e-11_dp) &
            error stop 'parallel momentum not conserved'
        if (abs(sum(w*v**2) - k_init)/k_init > 1.0e-11_dp) &
            error stop 'kinetic energy not conserved'
        if (max_defect <= 0.0_dp) error stop 'operator did not scatter'
    end subroutine test_conservation_fixed_seed

    subroutine test_isotropization()
        ! H-theorem sanity for the Lorentz part: a lam = 1 beam relaxes
        ! monotonically toward the isotropic state <lam> = 0, <lam^2> = 1/3.
        type(xorshift_rng_t) :: rng
        real(dp) :: lam(n_markers)
        real(dp), parameter :: nu_d_dt = 0.05_dp
        real(dp) :: mean_lam_20, mean_lam, mean_lam2
        integer :: step, i

        call rng_seed(rng, 777_int64)
        lam = 1.0_dp
        mean_lam_20 = 0.0_dp
        do step = 1, 400
            do i = 1, n_markers
                call scatter_pitch_bkp(lam(i), nu_d_dt, rng_sign(rng))
            end do
            if (step == 20) mean_lam_20 = sum(lam)/n_markers
        end do
        mean_lam = sum(lam)/n_markers
        mean_lam2 = sum(lam**2)/n_markers
        ! After 20 steps (nu_D t = 1) the beam mean sits near exp(-1).
        if (mean_lam_20 < 0.25_dp) error stop 'beam decayed too fast'
        if (mean_lam_20 > 0.50_dp) error stop 'beam decayed too slowly'
        if (abs(mean_lam) > 0.05_dp) error stop 'beam did not isotropize'
        if (abs(mean_lam2 - 1.0_dp/3.0_dp) > 0.03_dp) &
            error stop 'pitch distribution not isotropic'
    end subroutine test_isotropization

    subroutine test_thermalization()
        ! H-theorem sanity for the energy part: a monoenergetic population
        ! at u = E/T = 4 relaxes under the raw test-particle operator to
        ! the background Maxwellian moments <u> = 3/2, <u^2> = 15/4.
        type(xorshift_rng_t) :: rng
        real(dp) :: v(n_markers), lam(n_markers)
        real(dp), parameter :: v_thermal = 1.0_dp
        real(dp), parameter :: nu_ref = 1.0_dp, dt = 0.05_dp
        real(dp) :: mean_u, mean_u2
        integer :: step, i

        call rng_seed(rng, 4242_int64)
        v = 2.0_dp*v_thermal
        do i = 1, n_markers
            lam(i) = 2.0_dp*rng_uniform(rng) - 1.0_dp
        end do
        do step = 1, 2000
            do i = 1, n_markers
                call collide_marker(v(i), lam(i), dt, nu_ref, v_thermal, rng)
            end do
        end do
        mean_u = sum((v/v_thermal)**2)/n_markers
        mean_u2 = sum((v/v_thermal)**4)/n_markers
        if (abs(mean_u - 1.5_dp) > 0.1_dp) &
            error stop 'energy did not thermalize to background T'
        if (abs(mean_u2 - 3.75_dp) > 0.5_dp) &
            error stop 'energy spread is not Maxwellian'
    end subroutine test_thermalization

    subroutine test_maxwellian_stationarity()
        ! The background Maxwellian is a fixed point of the conservative
        ! operator: moments stay at their Maxwellian values.
        type(xorshift_rng_t) :: rng
        real(dp) :: w(n_markers), v(n_markers), lam(n_markers)
        real(dp), parameter :: v_thermal = 1.0_dp
        real(dp), parameter :: nu_ref = 1.0_dp, dt = 0.05_dp
        real(dp) :: mean_u, mean_u2, mean_lam
        integer :: step, ierr

        call rng_seed(rng, 99991_int64)
        call sample_maxwellian(rng, v_thermal, v, lam)
        w = 1.0_dp
        do step = 1, 100
            call apply_collisions_conservative(w, v, lam, dt, nu_ref, &
                                               v_thermal, rng, ierr)
            if (ierr /= 0) error stop 'stationarity step failed'
        end do
        mean_u = sum((v/v_thermal)**2)/n_markers
        mean_u2 = sum((v/v_thermal)**4)/n_markers
        mean_lam = sum(lam)/n_markers
        if (abs(mean_u - 1.5_dp) > 0.1_dp) &
            error stop 'Maxwellian mean energy drifted'
        if (abs(mean_u2 - 3.75_dp) > 0.5_dp) &
            error stop 'Maxwellian energy spread drifted'
        if (abs(mean_lam) > 0.05_dp) error stop 'Maxwellian anisotropized'
    end subroutine test_maxwellian_stationarity

end program test_conservative_collisions
