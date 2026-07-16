module xorshift_rng_mod
    ! Deterministic xorshift64 pseudo-random stream for Monte Carlo kernels.
    ! The state update uses only XOR and bit shifts (no wrapping multiplies),
    ! so a fixed seed reproduces the stream bit-for-bit across compilers.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64

    implicit none
    private

    type, public :: xorshift_rng_t
        integer(int64) :: state = 88172645463325252_int64
        logical :: have_cached_normal = .false.
        real(dp) :: cached_normal = 0.0_dp
    end type xorshift_rng_t

    public :: rng_seed, rng_uniform, rng_normal, rng_sign

    real(dp), parameter :: two_pi = 8.0_dp*atan(1.0_dp)
    real(dp), parameter :: scale_53 = 1.0_dp/9007199254740992.0_dp ! 2^-53

contains

    subroutine rng_seed(rng, seed)
        type(xorshift_rng_t), intent(out) :: rng
        integer(int64), intent(in) :: seed

        integer :: k

        rng%state = ieor(seed, 88172645463325252_int64)
        if (rng%state == 0_int64) rng%state = 88172645463325252_int64
        do k = 1, 8
            call advance(rng%state)
        end do
    end subroutine rng_seed

    pure subroutine advance(state)
        integer(int64), intent(inout) :: state

        state = ieor(state, shiftl(state, 13))
        state = ieor(state, shiftr(state, 7))
        state = ieor(state, shiftl(state, 17))
    end subroutine advance

    function rng_uniform(rng) result(u)
        ! Uniform deviate strictly inside (0, 1): the top 53 state bits are
        ! offset by one half unit, so u is never exactly 0 or 1.
        type(xorshift_rng_t), intent(inout) :: rng
        real(dp) :: u

        call advance(rng%state)
        u = (real(shiftr(rng%state, 11), dp) + 0.5_dp)*scale_53
    end function rng_uniform

    function rng_normal(rng) result(g)
        ! Standard normal deviate via Box-Muller with one cached value.
        type(xorshift_rng_t), intent(inout) :: rng
        real(dp) :: g

        real(dp) :: amplitude, angle

        if (rng%have_cached_normal) then
            rng%have_cached_normal = .false.
            g = rng%cached_normal
        else
            amplitude = sqrt(-2.0_dp*log(rng_uniform(rng)))
            angle = two_pi*rng_uniform(rng)
            g = amplitude*cos(angle)
            rng%cached_normal = amplitude*sin(angle)
            rng%have_cached_normal = .true.
        end if
    end function rng_normal

    function rng_sign(rng) result(s)
        ! Random sign, +1 or -1 with equal probability.
        type(xorshift_rng_t), intent(inout) :: rng
        real(dp) :: s

        s = sign(1.0_dp, rng_uniform(rng) - 0.5_dp)
    end function rng_sign

end module xorshift_rng_mod
