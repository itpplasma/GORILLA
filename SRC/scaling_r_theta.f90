module scaling_r_theta
implicit none
! 
  contains
    function scaling_r(r) result(r_scaled) !restrictions: f(0)=0, f(1)=1
        double precision, dimension(:), intent(in) :: r !r_frac from points_2d
        double precision, dimension(size(r)) :: r_scaled
!         
        r_scaled = r !explicit scaling function is defined here, function should fulfill restrictions stated above
!         
!         
    end function
! 
    function scaling_theta(theta) result(theta_scaled) !restrictions: f(0)=0, f(1)=1
        use tetra_grid_settings_mod, only: grid_kind, theta_xpoint_clustering, theta_xpoint_clustering_width
        use constants, only: pi
        double precision, dimension(:), intent(in) :: theta
        double precision, dimension(size(theta)) :: theta_scaled
        double precision :: a, w
!         
        theta_scaled = theta !explicit scaling function is defined here, function should fulfill restrictions stated above
!
        !grid_kind = 2: cluster points toward theta = 0 (the X-point for theta0_at_xpoint = .true.).
        !Point spacing proportional to 1 - a*g(t), g = exp(-(t/w)^2) + exp(-((t-1)/w)^2) (periodic
        !images at t = 0 and 1); f(t) = (t - a*G(t)) / (1 - a*G(1)) with G the integral of g, so f(0)=0, f(1)=1,
        !and f is monotonic for a < 1, w <= 0.25.
        if ((grid_kind .eq. 2) .and. (theta_xpoint_clustering .gt. 0.d0)) then
            a = theta_xpoint_clustering
            w = theta_xpoint_clustering_width
            theta_scaled = (theta - a * big_g(theta)) / (1.d0 - a * big_g_scalar(1.d0))
        endif
!
        contains
!
        function big_g(t) result(g)
            double precision, dimension(:), intent(in) :: t
            double precision, dimension(size(t)) :: g
            g = 0.5d0 * w * sqrt(pi) * (erf(t / w) + erf((t - 1.d0) / w) - erf(-1.d0 / w))
        end function big_g
!
        function big_g_scalar(t) result(g)
            double precision, intent(in) :: t
            double precision :: g
            g = 0.5d0 * w * sqrt(pi) * (erf(t / w) + erf((t - 1.d0) / w) - erf(-1.d0 / w))
        end function big_g_scalar
!
    end function
! 
end module scaling_r_theta
