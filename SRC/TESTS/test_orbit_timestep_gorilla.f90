module test_orbit_timestep_gorilla
	use orbit_timestep_gorilla_mod, only: initialize_gorilla
	use funit
	implicit none

contains

	@before
	subroutine prepare_for_initialize_gorilla()
    	use tetra_grid_settings_mod, only:      load_tetra_grid_inp
    	use gorilla_settings_mod, only:         load_gorilla_inp
		
		call load_tetra_grid_inp()
		call load_gorilla_inp()
		
	end subroutine prepare_for_initialize_gorilla
	
	@test
	@disable
	subroutine test_initialize_gorilla()
		call initialize_gorilla()
	end subroutine


end module test_orbit_timestep_gorilla
