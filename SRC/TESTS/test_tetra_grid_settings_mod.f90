module test_tetra_grid_settings_mod
   use tetra_grid_settings_mod
   use funit
   implicit none
!
! Background:
! Similar to test_gorilla_settings_mod, there exist a seperate input file for the generation of the
! tetrahedron mesh with a corresponding blueprint tetra_grid.inp file. As for all input files, a short
! explanation of the various parameters can be found in the blueprint input file itself. All input 
! variables are public protected in the module tetra_grid_settings_mod, but there exist writing routines
! for some of them.
!
! Georg Graßler (25.01.2023)
!       
contains

	@before
	subroutine test_load_tetra_grid_inp()
		call load_tetra_grid_inp()
	end subroutine test_load_tetra_grid_inp
	
	@test
    ! Description of test_load_tetra_grid_inp_values:
    ! The test checks whether the default values in the blueprint file
    ! are unchanged (grid_kind = 3, filename_mesh_rphiz = 'mesh_rphiz.obj', ...).
    !
    ! Georg Graßler (25.01.2023)
    !  
   subroutine test_load_tetra_grid_inp_values()
   	
   	double precision :: sfc_s_min_compare
   	sfc_s_min_compare = 0.1d0
   	
   	@assertTrue(boole_n_field_periods)
   	@assertFalse(boole_write_mesh_obj)
   	@assertEqual(filename_mesh_rphiz, 'mesh_rphiz.obj')
   	@assertEqual(filename_mesh_sthetaphi, 'mesh_sthetaphi.obj')
   	@assertEqual(g_file_filename, 'MHD_EQUILIBRIA/g_file_for_test')
   	@assertEqual(grid_kind, 3)
   	@assertEqual(n1, 100)
   	@assertEqual(n2, 40)
   	@assertEqual(n3, 40)
   	@assertEqual(n_field_periods_manual, 1)
   	@assertEqual(netcdf_filename, 'MHD_EQUILIBRIA/netcdf_file_for_test.nc')
   	@assertEqual(sfc_s_min, sfc_s_min_compare, tolerance=1e-13) 
   	@assertTrue(theta0_at_xpoint)
   	@assertEqual(theta_geom_flux,1)
   	
   end subroutine test_load_tetra_grid_inp_values
   
   @test
   ! Description of test_set_n_field_periods:
   ! The number of periods along the toroidal direction in the field configuration 
   ! is given by n_field_periods. This quantity can either be determined from the
   ! field data automatically, or been set directly by the input variable 
   ! n_field_periods_manual. The latter variant is checked in this test.
   !
   ! Georg Graßler (25.01.2023)
   !  
   subroutine test_set_n_field_periods()
   
   	call set_n_field_periods(n_field_periods_manual)
   	@assertEqual(n_field_periods, n_field_periods_manual, message = "Testing set_n_field_periods")
   
   end subroutine test_set_n_field_periods
   
end module test_tetra_grid_settings_mod
