module test_pusher_tetra_rk
	use pusher_tetra_rk_mod, only: initialize_const_motion_rk !, perpinv, perpinv2
   use funit
   implicit none
   
contains
	
	@test
	!@disable
	subroutine test_initialize_const_motion_rk()
	
		double precision  :: perpinv_in,perpinv2_in
		perpinv_in = 3.0d0
		perpinv2_in = 2.5d0
		
		call initialize_const_motion_rk(perpinv_in, perpinv2_in)
		
		!perpinv and perpinv2 are not public, due to this they are not able to inherit from pusher_tetra_rk_mod, so they are not testable in this way
		!@assertEqual(perpinv, perpinv_in, tolerance = 1e-13, message = "Testing initialize_const_motion_rk perpinv")
		!@assertEqual(perpinv2, perpinv2_in, tolerance = 1e-13, message = "Testing initialize_const_motion_rk perpinv2")
		
	end subroutine test_initialize_const_motion_rk
	
end module
