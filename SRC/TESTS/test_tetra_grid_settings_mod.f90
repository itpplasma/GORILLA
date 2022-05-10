module test_tetra_grid_settings_mod
   use tetra_grid_settings_mod
   use funit
   implicit none
   
contains

	@before
	subroutine test_load_tetra_grid_inp()
		call load_tetra_grid_inp()
	end subroutine test_load_tetra_grid_inp
	
	@test
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
   subroutine test_set_n_field_periods()
   
   	call set_n_field_periods(n_field_periods_manual)
   	@assertEqual(n_field_periods, n_field_periods_manual, message = "Testing set_n_field_periods")
   
   end subroutine test_set_n_field_periods
   
   @test
   subroutine test_set_grid_size()
   
   	call set_grid_size([n1,n2,n3])
   	@assertEqual(grid_size, [n1,n2,n3], message = "Testing set_grid_size")
   
   end subroutine test_set_grid_size
   
end module test_tetra_grid_settings_mod
