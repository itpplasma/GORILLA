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
        double precision, dimension(:), intent(in) :: theta
        double precision, dimension(size(theta)) :: theta_scaled
!         
        theta_scaled = theta !explicit scaling function is defined here, function should fulfill restrictions stated above
!         
    end function
! 
end module scaling_r_theta
