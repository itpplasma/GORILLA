module cylindrical_characteristics_mod
    ! Exact Gaussian averages of cylindrical guiding-centre characteristics.
    ! u=v_parallel/sqrt(T/m), du=-u*d(tau)+sqrt(2)*dW, tau=nu*t.
    ! The helical phase is integral(x1*u-x2) d(tau). No magnetic drifts or FLR.
    ! I(m,n)=integral E[u(tau)**m*u(0)**n*exp(-i*phase)] d(tau).
    ! This is the number-conserving OU resolvent, not energy-restoring collisions.
    use, intrinsic :: iso_fortran_env, only: real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: cylindrical_ou_resolvent
contains
    subroutine cylindrical_ou_resolvent(x1, x2, moments, info, refinement)
        real(real64), intent(in) :: x1, x2
        complex(real64), intent(out) :: moments(0:3,0:3)
        integer, intent(out) :: info
        integer, intent(in), optional :: refinement
        real(real64), parameter :: nodes(8) = [ &
            -0.9602898564975363_real64, -0.7966664774136267_real64, &
            -0.5255324099163290_real64, -0.1834346424956498_real64, &
             0.1834346424956498_real64,  0.5255324099163290_real64, &
             0.7966664774136267_real64,  0.9602898564975363_real64 ]
        real(real64), parameter :: weights(8) = [ &
            0.1012285362903763_real64, 0.2223810344533745_real64, &
            0.3137066458778873_real64, 0.3626837833783620_real64, &
            0.3626837833783620_real64, 0.3137066458778873_real64, &
            0.2223810344533745_real64, 0.1012285362903763_real64 ]
        complex(real64), parameter :: ii=(0.0_real64,1.0_real64)
        real(real64), parameter :: tail_time=36.0_real64
        real(real64) :: h, centre, tau, decay, omdecay, variance_half
        real(real64) :: coeff(0:3,0:3)
        complex(real64) :: shifted_mean, phase, polynomials(0:3,0:3), denom
        integer :: panels, panel, k, m, n, factor
        moments=0.0_real64
        info=1
        if (.not.ieee_is_finite(x1) .or. .not.ieee_is_finite(x2)) return
        factor=1
        if (present(refinement)) factor=refinement
        if (factor<1) return
        ! Exactly parallel-resonant radius: use the analytic Hermite expansion.
        ! At x1=x2=0 the density resolvent has no stationary solution.
        if (x1==0.0_real64) then
            if (x2==0.0_real64) return
            coeff=0.0_real64
            coeff(0,0)=1.0_real64
            coeff(1,1)=1.0_real64
            coeff(2,0)=1.0_real64
            coeff(2,2)=sqrt(2.0_real64)
            coeff(3,1)=3.0_real64
            coeff(3,3)=sqrt(6.0_real64)
            do n=0,3
                do m=0,3
                    do k=0,3
                        moments(m,n)=moments(m,n)+coeff(m,k)*coeff(n,k)/ &
                            (real(k,real64)-ii*x2)
                    end do
                end do
            end do
        else
            ! Eight-point Gauss panels resolve both streaming and Doppler phase.
            ! The OU covariance has saturated at tau=36. Integrate the remaining
            ! exponential tail analytically; do not add a spurious Krook damping.
            if (max(abs(x1),abs(x2))>1.0e5_real64) return
            panels=ceiling(2.0_real64*tail_time*max(1.0_real64,abs(x1),abs(x2)))
            if (panels>huge(panels)/factor) return
            panels=panels*factor
            h=tail_time/real(panels,real64)
            do panel=1,panels
                centre=(real(panel,real64)-0.5_real64)*h
                do k=1,8
                    tau=centre+0.5_real64*h*nodes(k)
                    decay=exp(-tau)
                    if (tau<1.0e-3_real64) then
                        omdecay=tau*(1.0_real64-tau/2.0_real64+tau**2/6.0_real64 &
                            -tau**3/24.0_real64+tau**4/120.0_real64)
                        variance_half=tau**2*(0.5_real64-tau/6.0_real64 &
                            +tau**2/24.0_real64-tau**3/120.0_real64+tau**4/720.0_real64)
                    else
                        omdecay=1.0_real64-decay
                        variance_half=tau-omdecay
                    end if
                    shifted_mean=-ii*x1*omdecay
                    phase=exp(-x1*x1*variance_half+ii*x2*tau)
                    call gaussian_polynomials(shifted_mean,decay,polynomials)
                    moments=moments+(0.5_real64*h*weights(k))*phase*polynomials
                end do
            end do
            call gaussian_polynomials(-ii*x1,0.0_real64,polynomials)
            denom=cmplx(x1*x1,-x2,real64)
            phase=exp(-x1*x1*(tail_time-1.0_real64)+ii*x2*tail_time)
            moments=moments+phase*polynomials/denom
        end if
        if (all(ieee_is_finite(real(moments))) .and. &
            all(ieee_is_finite(aimag(moments)))) info=0
    end subroutine

    pure subroutine gaussian_polynomials(mu,cov,polys)
        complex(real64), intent(in) :: mu
        real(real64), intent(in) :: cov
        complex(real64), intent(out) :: polys(0:3,0:3)
        integer :: m,n
        polys=0.0_real64
        polys(0,0)=1.0_real64
        polys(0,1)=mu
        do n=2,3
            polys(0,n)=mu*polys(0,n-1)+real(n-1,real64)*polys(0,n-2)
        end do
        do m=1,3
            do n=0,3
                polys(m,n)=mu*polys(m-1,n)
                if (m>1) polys(m,n)=polys(m,n)+real(m-1,real64)*polys(m-2,n)
                if (n>0) polys(m,n)=polys(m,n)+real(n,real64)*cov*polys(m-1,n-1)
            end do
        end do
    end subroutine
end module
