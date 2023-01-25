module test_scaling_r_theta
   use scaling_r_theta, only: scaling_r, scaling_theta
   use funit
   implicit none
!
! Background:
! The grid on which GORILLA operates is by default made out of equidistant points placed in in the corresponding 
! coordinates space (e.g. 10 equidistant values between R_min and R_max). One can choose to increase the number
! of points to refine the grid, but if one wants to resolve effects only in a destinct region of space, one opts
! to instead use scaling functions that places more of the grid points in the region of interest, while still
! keeping the total number of points the same (for performance). The desired scaling functions can be set in the
! module scaling_r_theta. By default, there is no scaling and the scaling functions just return the inputs unchanged.
!
! Georg Graßler (25.01.2023)
!       
contains

	@test
    ! Description of test_scaling_r:
    ! By default no scaling is set. The provided trial array of "radius"-values in test_scaling_r should therefore be
    ! the same as the output of the scaling routine.
    !
    ! Georg Graßler (25.01.2023)
    !  
	subroutine test_scaling_r()
        double precision, dimension(3):: r
        double precision, dimension(size(r)) :: r_scaled
        
        r = [1.d0,2.d0,3.d0]
        r_scaled = scaling_r(r)
        
        @assertEqual(r_scaled,r,message = "Testing scaling_r")
        
	end subroutine test_scaling_r
	
	@test
    ! Description of test_scaling_theta:
    ! By default no scaling is set. The provided trial array of "angle"-values in test_scaling_theta should therefore be
    ! the same as the output of the scaling routine.
    !
    ! Georg Graßler (25.01.2023)
    !  
	subroutine test_scaling_theta()
        double precision, dimension(3):: theta
        double precision, dimension(size(theta)) :: theta_scaled
        
        theta = [1.d0,2.d0,3.d0]
        theta_scaled = scaling_theta(theta)
        
        @assertEqual(theta_scaled,theta,message = "Testing scaling_theta")
        
	end subroutine test_scaling_theta
	
end module test_scaling_r_theta
