module test_gorilla_settings_mod
   use gorilla_settings_mod, only: set_eps_Phi, eps_Phi, load_gorilla_inp
   use funit
   implicit none
!
! Background:
! GORILLA is controlled via *.inp files that designate states of boolean variables, integers for switches,
! values of parameters and so on. The *.inp files are read in by corresponding modules, in the case considered
! here: gorilla_settings_mod for the gorilla.inp file, which lists amongst others a parameter eps_Phi, which
! is used to parameterize the electric field. As the main focus of GORILLA are strong magnetic fields, the electric
! contribution is currently implemented rather simple: The value of the electric potential Phi_elect is calculated by
!
! Phi_elect = eps_Phi*A_x2
!
! where A_x2 denotes the second component of the vector potential. In flux coodinates this leads to the contours of
! constant electric potential Phi_elect aligning with the flux surfaces of the magnetic field. The parameter eps_Phi
! is like all the other inputs public protected in the module, so is usually only intended to be read. However, there
! is a writting routine in place, that lets one specifically change the value of eps_Phi during runtime (set_eps_Phi).
! By default, eps_Phi = 0 in the blueprint gorilla.inp file.
!
! Georg Graßler (25.01.2023)
!    
contains

	@before
	subroutine test_load_gorilla_inp()
		call load_gorilla_inp()
	end subroutine

	@test
    ! Description of test_set_eps_Phi:
    ! The test checks first whether if the default case of no electric field (eps_Phi = 0) is set,
    ! by loading the blueprint input file and checking the public protected variable eps_Phi. After 
    ! that it tries to change its value with the writer routine set_eps_Phi(ep_Phi_in).
    !
    ! Georg Graßler (25.01.2023)
    !    
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
