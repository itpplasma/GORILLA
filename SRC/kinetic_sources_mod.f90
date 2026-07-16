module kinetic_sources_mod
    ! Particle and energy source skeleton with balance accounting for the
    ! WP4 kinetic loop. Particle sources add markers drawn from a Maxwellian
    ! at the local thermal speed; the energy source heats existing markers.
    ! Every action records the injected particle, parallel-momentum, and
    ! energy moments so a transport run can audit its balance exactly.
    !
    ! CGS units: mass g, velocity cm/s, momentum g cm/s, energy erg.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use xorshift_rng_mod, only: xorshift_rng_t, rng_normal

    implicit none
    private

    type, public :: source_balance_t
        real(dp) :: particles_injected = 0.0_dp ! sum of injected weights
        real(dp) :: momentum_injected = 0.0_dp  ! sum of w m vpar, g cm/s
        real(dp) :: energy_injected = 0.0_dp    ! sum of w m v^2/2, erg
        integer :: n_injections = 0
        integer :: n_heatings = 0
    end type source_balance_t

    public :: reset_balance, sample_maxwellian_marker
    public :: inject_maxwellian_markers, heat_markers

contains

    subroutine reset_balance(balance)
        type(source_balance_t), intent(out) :: balance

        balance = source_balance_t()
    end subroutine reset_balance

    subroutine sample_maxwellian_marker(rng, v_thermal, v, lam)
        ! Draw (v, lam) from an isotropic Maxwellian with thermal speed
        ! v_thermal = sqrt(2 T/m): three Cartesian components of variance
        ! T/m each.
        type(xorshift_rng_t), intent(inout) :: rng
        real(dp), intent(in) :: v_thermal
        real(dp), intent(out) :: v, lam

        real(dp) :: sigma, vx, vy, vz

        sigma = v_thermal/sqrt(2.0_dp)
        vx = sigma*rng_normal(rng)
        vy = sigma*rng_normal(rng)
        vz = sigma*rng_normal(rng)
        v = sqrt(vx**2 + vy**2 + vz**2)
        if (v > 0.0_dp) then
            lam = vz/v
        else
            lam = 0.0_dp
        end if
    end subroutine sample_maxwellian_marker

    subroutine inject_maxwellian_markers(balance, rng, v_thermal, mass, &
                                         weight_each, n_new, w, v, lam, &
                                         n_active, ierr)
        ! Append n_new Maxwellian markers of weight weight_each to the
        ! marker slots n_active+1 .. n_active+n_new and record the injected
        ! particle, parallel-momentum, and energy moments from the markers
        ! actually stored.
        !
        ! ierr: 0 ok; 1 inconsistent array sizes; 2 invalid arguments;
        ! 3 capacity exceeded. On nonzero ierr nothing is modified.
        type(source_balance_t), intent(inout) :: balance
        type(xorshift_rng_t), intent(inout) :: rng
        real(dp), intent(in) :: v_thermal, mass, weight_each
        integer, intent(in) :: n_new
        real(dp), intent(inout) :: w(:), v(:), lam(:)
        integer, intent(inout) :: n_active
        integer, intent(out) :: ierr

        integer :: i

        ierr = 1
        if (size(v) /= size(w)) return
        if (size(lam) /= size(w)) return
        ierr = 2
        if (v_thermal <= 0.0_dp) return
        if (mass <= 0.0_dp) return
        if (weight_each <= 0.0_dp) return
        if (n_new < 1) return
        if (n_active < 0) return
        ierr = 3
        if (n_active + n_new > size(w)) return

        do i = n_active + 1, n_active + n_new
            w(i) = weight_each
            call sample_maxwellian_marker(rng, v_thermal, v(i), lam(i))
            balance%particles_injected = balance%particles_injected + w(i)
            balance%momentum_injected = balance%momentum_injected &
                                        + w(i)*mass*v(i)*lam(i)
            balance%energy_injected = balance%energy_injected &
                                      + 0.5_dp*w(i)*mass*v(i)**2
        end do
        n_active = n_active + n_new
        balance%n_injections = balance%n_injections + 1
        ierr = 0
    end subroutine inject_maxwellian_markers

    subroutine heat_markers(balance, mass, w, v, lam, n_active, delta_energy, &
                            ierr)
        ! Add delta_energy (erg, may be negative for cooling) to the active
        ! marker population by a uniform speed rescale, and record the
        ! measured energy and parallel-momentum change. The lab-frame
        ! rescale changes the momentum moment by (alpha - 1) times its
        ! value; a flow-frame-preserving heating is deferred to the orbit
        ! coupling.
        !
        ! ierr: 0 ok; 1 inconsistent sizes or no active markers; 2 invalid
        ! arguments; 3 heating undefined (nonpositive kinetic energy or
        ! cooling past zero). On nonzero ierr nothing is modified.
        type(source_balance_t), intent(inout) :: balance
        real(dp), intent(in) :: mass
        real(dp), intent(in) :: w(:)
        real(dp), intent(inout) :: v(:)
        real(dp), intent(in) :: lam(:)
        integer, intent(in) :: n_active
        real(dp), intent(in) :: delta_energy
        integer, intent(out) :: ierr

        real(dp) :: kin0, kin1, mom0, alpha

        ierr = 1
        if (size(v) /= size(w)) return
        if (size(lam) /= size(w)) return
        if (n_active < 1) return
        if (n_active > size(w)) return
        ierr = 2
        if (mass <= 0.0_dp) return
        if (any(w(1:n_active) <= 0.0_dp)) return
        ierr = 3
        kin0 = 0.5_dp*mass*sum(w(1:n_active)*v(1:n_active)**2)
        if (kin0 <= 0.0_dp) return
        if (kin0 + delta_energy <= 0.0_dp) return

        alpha = sqrt(1.0_dp + delta_energy/kin0)
        mom0 = mass*sum(w(1:n_active)*v(1:n_active)*lam(1:n_active))
        v(1:n_active) = alpha*v(1:n_active)
        kin1 = 0.5_dp*mass*sum(w(1:n_active)*v(1:n_active)**2)

        balance%energy_injected = balance%energy_injected + (kin1 - kin0)
        balance%momentum_injected = balance%momentum_injected &
                                    + (alpha - 1.0_dp)*mom0
        balance%n_heatings = balance%n_heatings + 1
        ierr = 0
    end subroutine heat_markers

end module kinetic_sources_mod
