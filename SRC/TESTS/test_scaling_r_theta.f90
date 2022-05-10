module test_scaling_r_theta
   use scaling_r_theta, only: scaling_r, scaling_theta
   use funit
   implicit none
   
contains

	@test
	subroutine test_scaling_r()
        double precision, dimension(3):: r
        double precision, dimension(size(r)) :: r_scaled
        
        r = [1.d0,2.d0,3.d0]
        r_scaled = scaling_r(r)
        
        @assertEqual(r_scaled,r,message = "Testing scaling_r")
        
	end subroutine test_scaling_r
	
	@test
	subroutine test_scaling_theta()
        double precision, dimension(3):: theta
        double precision, dimension(size(theta)) :: theta_scaled
        
        theta = [1.d0,2.d0,3.d0]
        theta_scaled = scaling_theta(theta)
        
        @assertEqual(theta_scaled,theta,message = "Testing scaling_theta")
        
	end subroutine test_scaling_theta
	
end module test_scaling_r_theta
