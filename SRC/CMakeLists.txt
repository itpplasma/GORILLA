add_library(GORILLA STATIC
	SetWorkingPrecision.f90
	contrib/Polynomial234RootSolvers.f90 
	constants_mod.f90 
	tetra_grid_settings_mod.f90
	gorilla_settings_mod.f90
	various_functions_mod.f90 
	gorilla_diag_mod.f90 
	canonical_coordinates_mod.f90 
	nctools_module.f90 
	contrib/rkf45.f90 
	odeint_rkf45.f90 
	runge_kutta_mod.f90 
	magfie.f90 
	chamb_m.f90 
	vmecinm_m.f90 
	spl_three_to_five_mod.f90 
	spline_vmec_data.f90 
	new_vmec_allocation_stuff.f90 
	binsrc.f90 
	field_divB0.f90 
	scaling_r_theta.f90
	field_line_integration_for_SYNCH.f90 
	preload_for_SYNCH.f90
	plag_coeff.f90
	magdata_in_symfluxcoord.f90 
	points_2d.f90
	circular_mesh.f90
	circular_mesh_SOLEDGE3X_EIRENE.f90
	tetra_grid_mod.f90 
	make_grid_rect.f90 
	bdivfree.f90 
	tetra_physics_mod.f90
	tetra_physics_poly_precomp_mod.f90 
	differentiate.f90 
	spline5_RZ.f90 
	supporting_functions_mod.f90 
	pusher_tetra_func_mod.f90 
	pusher_tetra_poly.f90 
	pusher_tetra_rk.f90 
	get_canonical_coordinates.f90 
	orbit_timestep_gorilla.f90 
	gorilla_plot_mod.f90
	)

if ($ENV{GORILLA_COVERAGE} STREQUAL "TRUE")
	target_compile_options(GORILLA
		PRIVATE --coverage
		)
endif()

add_executable(test_gorilla_main.x
  test_gorilla_main.f90
  )

target_link_libraries(test_gorilla_main.x GORILLA netcdff netcdf lapack)

set_target_properties (GORILLA PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

target_include_directories(GORILLA PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

if ($ENV{GORILLA_COVERAGE} STREQUAL "TRUE")
	add_subdirectory(TESTS)
endif()


