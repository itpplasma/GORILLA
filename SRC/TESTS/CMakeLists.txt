add_pfunit_ctest (test_various_functions_mod
  TEST_SOURCES test_various_functions_mod.f90
  LINK_LIBRARIES GORILLA
  )

add_pfunit_ctest (test_binsrc
  TEST_SOURCES test_binsrc.f90
  LINK_LIBRARIES GORILLA
  )
  
add_pfunit_ctest (test_chamb_m
  TEST_SOURCES test_chamb_m.f90
  LINK_LIBRARIES GORILLA
  )
  
add_pfunit_ctest (test_tetra_grid_settings_mod
  TEST_SOURCES test_tetra_grid_settings_mod.f90
  LINK_LIBRARIES GORILLA
  )
  
add_pfunit_ctest (test_scaling_r_theta
  TEST_SOURCES test_scaling_r_theta.f90
  LINK_LIBRARIES GORILLA
  )
  
add_pfunit_ctest (test_supporting_functions_mod
  TEST_SOURCES test_supporting_functions_mod.f90
  LINK_LIBRARIES GORILLA lapack netcdf netcdff
  )
  
add_pfunit_ctest (test_gorilla_settings_mod
  TEST_SOURCES test_gorilla_settings_mod.f90
  LINK_LIBRARIES GORILLA
  )
  
add_pfunit_ctest (test_pusher_tetra_rk
  TEST_SOURCES test_pusher_tetra_rk.f90
  LINK_LIBRARIES GORILLA lapack netcdf netcdff
  )
  
add_pfunit_ctest (test_orbit_timestep_gorilla
  TEST_SOURCES test_orbit_timestep_gorilla.f90
  LINK_LIBRARIES GORILLA lapack netcdf netcdff
  )
  
file(COPY ../../INPUT/tetra_grid.inp ../../INPUT/gorilla.inp DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
