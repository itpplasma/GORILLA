program test_kinetic_sources
    ! Deterministic tests of the WP4 source skeleton: injection fills
    ! marker slots with Maxwellian samples and the balance accounting
    ! matches the stored markers exactly; heating adds the requested
    ! energy to roundoff; error paths leave the state untouched.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use xorshift_rng_mod, only: xorshift_rng_t, rng_seed
    use kinetic_sources_mod, only: source_balance_t, reset_balance, &
        sample_maxwellian_marker, inject_maxwellian_markers, heat_markers

    implicit none

    integer, parameter :: capacity = 8192
    real(dp), parameter :: v_thermal = 1.3e8_dp, mass = 9.1094e-28_dp
    real(dp), parameter :: weight_each = 0.7_dp
    integer, parameter :: n_pre = 100, n_new = 5000

    type(xorshift_rng_t) :: rng
    type(source_balance_t) :: balance
    real(dp) :: w(capacity), v(capacity), lam(capacity)
    real(dp) :: kin_direct, mom_direct, kin0, kin1, delta_energy, mean_u
    real(dp) :: mean_lam, energy_before
    integer :: n_active, i, ierr

    call rng_seed(rng, 314159_int64)
    call reset_balance(balance)

    ! Pre-existing population outside the source accounting.
    n_active = n_pre
    do i = 1, n_pre
        w(i) = 1.0_dp
        call sample_maxwellian_marker(rng, v_thermal, v(i), lam(i))
    end do

    call inject_maxwellian_markers(balance, rng, v_thermal, mass, &
                                   weight_each, n_new, w, v, lam, n_active, &
                                   ierr)
    if (ierr /= 0) error stop 'injection failed'
    if (n_active /= n_pre + n_new) error stop 'active count wrong'
    if (balance%n_injections /= 1) error stop 'injection count wrong'

    ! Balance must equal the moments of the markers actually stored.
    kin_direct = 0.0_dp
    mom_direct = 0.0_dp
    do i = n_pre + 1, n_active
        kin_direct = kin_direct + 0.5_dp*w(i)*mass*v(i)**2
        mom_direct = mom_direct + w(i)*mass*v(i)*lam(i)
    end do
    if (abs(balance%particles_injected - n_new*weight_each) &
        > 1.0e-12_dp*n_new*weight_each) error stop 'particle balance wrong'
    if (abs(balance%energy_injected - kin_direct) > 1.0e-12_dp*kin_direct) &
        error stop 'energy balance wrong'
    if (abs(balance%momentum_injected - mom_direct) &
        > 1.0e-12_dp*mass*v_thermal*n_new) error stop 'momentum balance wrong'

    ! Sampled moments close to Maxwellian: <u> = 3/2, <lam> = 0.
    mean_u = sum((v(n_pre + 1:n_active)/v_thermal)**2)/n_new
    mean_lam = sum(lam(n_pre + 1:n_active))/n_new
    if (abs(mean_u - 1.5_dp) > 0.1_dp) error stop 'source spectrum wrong'
    if (abs(mean_lam) > 0.05_dp) error stop 'source not isotropic'

    ! Capacity overflow must not modify anything.
    energy_before = balance%energy_injected
    call inject_maxwellian_markers(balance, rng, v_thermal, mass, &
                                   weight_each, capacity, w, v, lam, &
                                   n_active, ierr)
    if (ierr /= 3) error stop 'capacity overflow not detected'
    if (n_active /= n_pre + n_new) error stop 'overflow changed count'
    if (balance%energy_injected /= energy_before) &
        error stop 'overflow changed balance'

    ! Invalid arguments are rejected.
    call inject_maxwellian_markers(balance, rng, v_thermal, mass, -1.0_dp, &
                                   1, w, v, lam, n_active, ierr)
    if (ierr /= 2) error stop 'negative weight accepted'

    ! Heating adds exactly the requested energy and records it.
    kin0 = 0.5_dp*mass*sum(w(1:n_active)*v(1:n_active)**2)
    delta_energy = 0.25_dp*kin0
    energy_before = balance%energy_injected
    call heat_markers(balance, mass, w, v, lam, n_active, delta_energy, ierr)
    if (ierr /= 0) error stop 'heating failed'
    kin1 = 0.5_dp*mass*sum(w(1:n_active)*v(1:n_active)**2)
    if (abs(kin1 - kin0 - delta_energy) > 1.0e-12_dp*kin0) &
        error stop 'heating energy mismatch'
    if (abs(balance%energy_injected - energy_before - (kin1 - kin0)) &
        > 1.0e-12_dp*kin0) error stop 'heating balance mismatch'
    if (balance%n_heatings /= 1) error stop 'heating count wrong'

    ! Cooling past zero energy is rejected without modification.
    call heat_markers(balance, mass, w, v, lam, n_active, -2.0_dp*kin1, ierr)
    if (ierr /= 3) error stop 'overdraining accepted'
    if (abs(0.5_dp*mass*sum(w(1:n_active)*v(1:n_active)**2) - kin1) &
        > 1.0e-12_dp*kin1) error stop 'rejected cooling changed state'

    ! Reset clears all accounting.
    call reset_balance(balance)
    if (balance%particles_injected /= 0.0_dp) error stop 'reset failed'
    if (balance%energy_injected /= 0.0_dp) error stop 'reset failed'
    if (balance%momentum_injected /= 0.0_dp) error stop 'reset failed'
    if (balance%n_injections /= 0) error stop 'reset failed'
    if (balance%n_heatings /= 0) error stop 'reset failed'

    print '(A)', 'PASS: kinetic source skeleton balance accounting'

end program test_kinetic_sources
