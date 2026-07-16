module conservative_collisions_mod
    ! Conservative Monte Carlo collision kernel for guiding-center delta-f
    ! markers (WP4): Boozer-Kuo-Petravic pitch-angle (Lorentz) scattering
    ! plus energy scattering with drag against a Maxwellian field-particle
    ! background [A. H. Boozer, G. Kuo-Petravic, Phys. Fluids 24, 851
    ! (1981)], followed by an ensemble correction that returns the
    ! parallel-momentum and energy transfer to the background, so
    ! like-particle collisions conserve both moments.
    !
    ! Conventions (CGS/Gaussian): velocity cm/s, temperature erg, density
    ! 1/cm^3, mass g, charge statcoul. x = v/v_T with v_T = sqrt(2 T/m).
    ! The normalized deflection frequency matches NEO-2
    ! (COMMON/collop_compute.f90, nu_D_hat):
    !     nu_d_hat(x) = 3 sqrt(pi)/4 * (erf(x) - G(x))/x^3,
    ! with G the Chandrasekhar function. The energy-scattering frequency is
    !     nu_e_hat(x) = x^2 * nu_par_hat(x) = 3 sqrt(pi)/2 * G(x)/x,
    ! the normalized parallel velocity diffusion nu_par_hat = 3 sqrt(pi)/2
    ! * G(x)/x^3 mapped to the energy variable, so the energy diffusion
    ! coefficient is D_EE = 2 T E nu_E and the BKP drift term makes the
    ! background Maxwellian the exact stationary state.
    !
    ! Both hats are normalized to the inverse self-collision time
    !     1/tau_aa = (4/(3 sqrt(pi))) * 4 pi n q^4 ln(Lambda)/(m^2 v_T^3)
    ! (NEO-2 convention), returned by nu_ref_self, so that
    ! nu_ref*nu_d_hat(v/v_T) is the physical deflection frequency
    ! nu_D = nu_perp/2 (Trubnikov/NRL) on a Maxwellian background.
    !
    ! Self-contained: no coupling to the orbit loop in this increment. All
    ! sampling weights must be positive; delta-f perturbation weights evolve
    ! separately and ride on top of the scattered markers.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use xorshift_rng_mod, only: xorshift_rng_t, rng_sign

    implicit none
    private

    type, public :: collision_correction_t
        ! Diagnostics of one conservative step: raw test-particle drift of
        ! the conserved moments and the correction actually applied.
        real(dp) :: momentum_defect = 0.0_dp ! change of sum(w v lam), weight*cm/s
        real(dp) :: energy_defect = 0.0_dp ! change of sum(w v^2), weight*cm^2/s^2
        real(dp) :: vpar_shift = 0.0_dp ! uniform parallel shift, cm/s
        real(dp) :: speed_scale = 1.0_dp ! thermal-frame rescale factor
    end type collision_correction_t

    public :: chandrasekhar_g, nu_d_hat, nu_e_hat, dln_nu_e_dxsq, nu_ref_self
    public :: scatter_pitch_bkp, scatter_energy_bkp, collide_marker
    public :: apply_collisions_conservative

    real(dp), parameter :: pi = 3.14159265358979323846_dp
    real(dp), parameter :: sqrt_pi = 1.7724538509055160273_dp
    ! Switch to Taylor expansions below this x (numerical stability of the
    ! cancellations in G and in d ln nu_E / d x^2, as in NEO-2).
    real(dp), parameter :: x_taylor = 1.0e-2_dp
    ! Speed floor (units of v_T) protecting the diverging Lorentz rate and
    ! the 1/x factors; markers reflected at u = 0 are re-floored here.
    real(dp), parameter, public :: collision_x_floor = 1.0e-3_dp
    ! Cap on nu*dt per substep; an operator step beyond this is not a valid
    ! diffusive Monte Carlo step anyway.
    real(dp), parameter, public :: collision_nu_dt_cap = 0.5_dp

contains

    pure function chandrasekhar_g(x) result(g)
        ! Chandrasekhar function G(x) = (erf(x) - x erf'(x))/(2 x^2) with a
        ! fifth-order Taylor expansion near x = 0 (as in NEO-2).
        real(dp), intent(in) :: x
        real(dp) :: g

        if (x <= x_taylor) then
            g = 2.0_dp*x/(3.0_dp*sqrt_pi) - 2.0_dp*x**3/(5.0_dp*sqrt_pi) &
                + x**5/(7.0_dp*sqrt_pi)
        else
            g = (erf(x) - x*2.0_dp/sqrt_pi*exp(-x**2))/(2.0_dp*x**2)
        end if
    end function chandrasekhar_g

    pure function nu_d_hat(x) result(nu)
        ! Normalized deflection frequency nu_D/nu_ref for like-particle
        ! collisions; identical to NEO-2's nu_D_hat. Requires x > 0.
        real(dp), intent(in) :: x
        real(dp) :: nu

        nu = 0.75_dp*sqrt_pi*(erf(x) - chandrasekhar_g(x))/x**3
    end function nu_d_hat

    pure function nu_e_hat(x) result(nu)
        ! Normalized energy-scattering frequency nu_E/nu_ref = x^2 times the
        ! parallel velocity diffusion frequency. Limit 1 at x -> 0; requires
        ! x > 0.
        real(dp), intent(in) :: x
        real(dp) :: nu

        nu = 1.5_dp*sqrt_pi*chandrasekhar_g(x)/x
    end function nu_e_hat

    pure function dln_nu_e_dxsq(x) result(d)
        ! d ln(nu_e_hat)/d(x^2) for the BKP drag correction; smooth limit
        ! -3/5 at x -> 0.
        real(dp), intent(in) :: x
        real(dp) :: d

        real(dp) :: g, dg

        if (x <= x_taylor) then
            d = -0.6_dp + (12.0_dp/175.0_dp)*x**2
        else
            g = chandrasekhar_g(x)
            dg = d_chandrasekhar_g(x)
            d = (dg/g - 1.0_dp/x)/(2.0_dp*x)
        end if
    end function dln_nu_e_dxsq

    pure function d_chandrasekhar_g(x) result(dg)
        ! Derivative of the Chandrasekhar function for x > x_taylor.
        real(dp), intent(in) :: x
        real(dp) :: dg

        real(dp) :: derf

        derf = 2.0_dp/sqrt_pi*exp(-x**2)
        dg = ((1.0_dp + x**2)*derf - erf(x)/x)/x**2
    end function d_chandrasekhar_g

    pure function nu_ref_self(dens, temp, mass, charge, coulomb_log) result(nu)
        ! Reference like-particle collision frequency (CGS): the inverse
        ! self-collision time
        !     nu_ref = 1/tau_aa
        !            = (4/(3 sqrt(pi))) * 4 pi n q^4 ln(Lambda)/(m^2 v_T^3),
        ! v_T = sqrt(2 T/m), the normalization of nu_d_hat and nu_e_hat, so
        ! that nu_D(v) = nu_ref*nu_d_hat(v/v_T) is the physical deflection
        ! frequency on a Maxwellian background of the same species
        ! (Helander & Sigmar nu_D = nu_hat (erf - G)/x^3 with
        ! nu_hat = 4 pi n q^4 ln(Lambda)/(m^2 v_T^3), i.e. Trubnikov/NRL
        ! nu_perp/2).
        real(dp), intent(in) :: dens, temp, mass, charge, coulomb_log
        real(dp) :: nu

        real(dp) :: v_thermal

        v_thermal = sqrt(2.0_dp*temp/mass)
        nu = (4.0_dp/(3.0_dp*sqrt_pi))*4.0_dp*pi*dens*charge**4*coulomb_log &
             /(mass**2*v_thermal**3)
    end function nu_ref_self

    pure subroutine scatter_pitch_bkp(lam, nu_d_dt, rand_sign)
        ! BKP Lorentz pitch-angle step:
        !     lam' = lam (1 - nu_D dt) +/- sqrt((1 - lam^2) nu_D dt).
        real(dp), intent(inout) :: lam
        real(dp), intent(in) :: nu_d_dt, rand_sign

        lam = lam*(1.0_dp - nu_d_dt) &
              + rand_sign*sqrt(max((1.0_dp - lam**2)*nu_d_dt, 0.0_dp))
        lam = max(-1.0_dp, min(1.0_dp, lam))
    end subroutine scatter_pitch_bkp

    pure subroutine scatter_energy_bkp(u, nu_e_dt, dln_nu_e, rand_sign)
        ! BKP energy step in u = E/T = x^2:
        !     u' = u - 2 nu_E dt (u - 3/2 - u dln nu_E/dx^2)
        !            +/- 2 sqrt(u nu_E dt),
        ! reflected at u = 0. The drag term keeps the background Maxwellian
        ! stationary for the diffusion coefficient D_uu = 2 u nu_E.
        real(dp), intent(inout) :: u
        real(dp), intent(in) :: nu_e_dt, dln_nu_e, rand_sign

        u = u - 2.0_dp*nu_e_dt*(u - 1.5_dp - u*dln_nu_e) &
            + 2.0_dp*rand_sign*sqrt(max(u*nu_e_dt, 0.0_dp))
        u = abs(u)
    end subroutine scatter_energy_bkp

    subroutine collide_marker(v, lam, dt, nu_ref, v_thermal, rng)
        ! One BKP test-particle step (pitch + energy) for a single marker.
        ! No conservation correction; use apply_collisions_conservative for
        ! the conserving ensemble operator.
        real(dp), intent(inout) :: v, lam
        real(dp), intent(in) :: dt, nu_ref, v_thermal
        type(xorshift_rng_t), intent(inout) :: rng

        real(dp) :: x, u, nu_d_dt, nu_e_dt

        x = max(v/v_thermal, collision_x_floor)
        nu_d_dt = min(nu_ref*dt*nu_d_hat(x), collision_nu_dt_cap)
        call scatter_pitch_bkp(lam, nu_d_dt, rng_sign(rng))
        u = x**2
        nu_e_dt = min(nu_ref*dt*nu_e_hat(x), collision_nu_dt_cap)
        call scatter_energy_bkp(u, nu_e_dt, dln_nu_e_dxsq(x), rng_sign(rng))
        v = v_thermal*sqrt(max(u, collision_x_floor**2))
    end subroutine collide_marker

    subroutine apply_collisions_conservative(w, v, lam, dt, nu_ref, v_thermal, &
                                             rng, ierr, correction)
        ! BKP pitch+energy step on all markers, then a two-stage correction
        ! restoring the weighted parallel momentum sum(w v lam) and kinetic
        ! energy sum(w v^2) to their pre-collision values: a uniform
        ! parallel shift fixes the momentum, then an energy rescale about
        ! the mean parallel flow fixes the energy without disturbing the
        ! momentum. Mass factors are constant across markers, so these
        ! moments per unit mass are the conserved like-particle momentum
        ! and energy that the field-particle response must restore.
        !
        ! ierr: 0 ok; 1 inconsistent array sizes or no markers; 2 invalid
        ! dt, nu_ref, v_thermal, or nonpositive weight; 3 degenerate
        ! ensemble (nonpositive thermal energy, rescale undefined).
        real(dp), intent(in) :: w(:)
        real(dp), intent(inout) :: v(:), lam(:)
        real(dp), intent(in) :: dt, nu_ref, v_thermal
        type(xorshift_rng_t), intent(inout) :: rng
        integer, intent(out) :: ierr
        type(collision_correction_t), intent(out), optional :: correction

        integer :: i, n
        real(dp) :: wsum, p0, k0, p1, k1, shift, u_mean, k_thermal, k_target
        real(dp) :: alpha, vpar, vperp2

        ierr = 1
        n = size(w)
        if (n == 0) return
        if (size(v) /= n) return
        if (size(lam) /= n) return
        ierr = 2
        if (dt <= 0.0_dp) return
        if (nu_ref <= 0.0_dp) return
        if (v_thermal <= 0.0_dp) return
        if (any(w <= 0.0_dp)) return

        wsum = sum(w)
        p0 = sum(w*v*lam)
        k0 = sum(w*v**2)

        do i = 1, n
            call collide_marker(v(i), lam(i), dt, nu_ref, v_thermal, rng)
        end do
        p1 = sum(w*v*lam)
        k1 = sum(w*v**2)

        ! Stage 1: uniform parallel shift restores the momentum moment.
        shift = (p0 - p1)/wsum
        do i = 1, n
            vpar = v(i)*lam(i) + shift
            vperp2 = v(i)**2*(1.0_dp - lam(i)**2)
            v(i) = sqrt(vpar**2 + vperp2)
            if (v(i) > 0.0_dp) then
                lam(i) = vpar/v(i)
            else
                lam(i) = 0.0_dp
            end if
        end do

        ! Stage 2: rescale about the mean parallel flow; the flow-frame
        ! momentum is zero, so the momentum moment stays fixed while the
        ! energy moment returns to k0 exactly.
        u_mean = p0/wsum
        k_thermal = sum(w*v**2) - wsum*u_mean**2
        k_target = k0 - wsum*u_mean**2
        ierr = 3
        if (k_thermal <= 0.0_dp) return
        if (k_target <= 0.0_dp) return
        alpha = sqrt(k_target/k_thermal)
        do i = 1, n
            vpar = u_mean + alpha*(v(i)*lam(i) - u_mean)
            vperp2 = alpha**2*v(i)**2*(1.0_dp - lam(i)**2)
            v(i) = sqrt(vpar**2 + vperp2)
            if (v(i) > 0.0_dp) then
                lam(i) = vpar/v(i)
            else
                lam(i) = 0.0_dp
            end if
        end do
        ierr = 0

        if (present(correction)) then
            correction%momentum_defect = p1 - p0
            correction%energy_defect = k1 - k0
            correction%vpar_shift = shift
            correction%speed_scale = alpha
        end if
    end subroutine apply_collisions_conservative

end module conservative_collisions_mod
