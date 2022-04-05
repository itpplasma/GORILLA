module test_gorilla_settings_mod
   use gorilla_settings_mod, only: set_eps_Phi, eps_Phi, load_gorilla_inp
   use funit
   implicit none
   
contains

	@before
	subroutine test_load_gorilla_inp()
		call load_gorilla_inp()
	end subroutine

	@test
	subroutine test_set_eps_Phi()
		
		double precision :: compare_eps_from_file
		double precision :: eps_Phi_in
		compare_eps_from_file = 0.d0
		
		@assertEqual(eps_Phi, compare_eps_from_file, tolerance = 1e-13, message = "Check if read data in OK")
		
		eps_Phi_in = 0.8
		
		call set_eps_Phi(eps_Phi_in)
		@assertEqual(eps_Phi, eps_Phi_in, tolerance = 1e-13, message = "Testing test_set_eps_Phi")
	
	end subroutine test_set_eps_Phi

end module test_gorilla_settings_mod
