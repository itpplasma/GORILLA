module test_orbit_timestep_gorilla
	use orbit_timestep_gorilla_mod, only: initialize_gorilla
	use funit
	implicit none
!
! Background:
! When GORILLA is loaded, several aspects have to be set up before the actual orbit
! integration scheme can be performed. This includes constructing the mesh, linearizing
! the fields over the tetrahedrons, setting the species-specific values for mass, charge
! and so on. All that is taken care of by the subroutine initialize_gorilla.
!
! Georg Graßler (25.01.2023)
!    
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
    ! Description of test_initialize_gorilla:
    ! The intend of the test is to check via the integrated check_tetra_overlaps if
    ! the generated tetrahedron mesh is sane in the used coordinate system. See in
    ! DOCUMENTATION/SUPPLEMENTAL_MATERIAL/Bauer for more details about issues concerning
    ! using different coordinate systems.
    !
    ! However, in the current version the program structure does not support the reading
    ! out of the sanity check result. Until this is changed this test stays incomplete
    ! and therefore DISABLED!
    !
    ! Georg Graßler (25.01.2023)
    !    
	subroutine test_initialize_gorilla()
		call initialize_gorilla()
	end subroutine


end module test_orbit_timestep_gorilla
