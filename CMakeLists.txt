cmake_minimum_required(VERSION 3.12)

project (GORILLA
  VERSION 1.0.0)

enable_language(Fortran)
enable_testing()

include_directories(/usr/local/include)
if(UNIX AND NOT APPLE)
	include_directories(/usr/include)
elseif(APPLE AND ${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "arm64")
	include_directories(/opt/homebrew/include)
endif()

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(netCDF REQUIRED)

if ("$ENV{GORILLA_COVERAGE}" STREQUAL "TRUE")
	find_package(PFUNIT REQUIRED)
	enable_testing()
endif()

add_compile_options(-g -fbacktrace -ffpe-trap=zero,overflow,invalid -fbounds-check -fopenmp)
add_link_options(-g -fbacktrace -ffpe-trap=zero,overflow,invalid -fbounds-check -fopenmp)

add_link_options(-L/usr/local/lib)
if(APPLE AND ${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "arm64")
	add_link_options(-L/opt/homebrew/lib)
endif()

if ("$ENV{GORILLA_COVERAGE}" STREQUAL "TRUE")
	add_link_options(--coverage)
endif()

add_subdirectory(SRC)

add_executable(test_gorilla_main.x
  SRC/test_gorilla_main.f90
)

target_link_libraries(test_gorilla_main.x GORILLA netcdff netcdf lapack)

add_executable(test_boozer_chartmap_driver.x
  SRC/TESTS/test_boozer_chartmap_driver.f90
)
target_link_libraries(test_boozer_chartmap_driver.x GORILLA netcdff netcdf lapack)
add_test(NAME boozer_chartmap_oracle COMMAND test_boozer_chartmap_driver.x)

add_executable(test_radial_electric_projection_driver.x
  SRC/TESTS/test_radial_electric_projection_driver.f90
)
target_link_libraries(test_radial_electric_projection_driver.x
  GORILLA netcdff netcdf lapack)
add_test(NAME radial_electric_projection_oracle
  COMMAND test_radial_electric_projection_driver.x)

add_executable(test_prescribed_er_driver.x
  SRC/TESTS/test_prescribed_er_driver.f90
)
target_link_libraries(test_prescribed_er_driver.x GORILLA netcdff netcdf lapack)
add_test(NAME prescribed_er_oracle COMMAND test_prescribed_er_driver.x)
set_tests_properties(prescribed_er_oracle PROPERTIES
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
add_executable(test_hamiltonian_time_driver.x
  SRC/TESTS/test_hamiltonian_time_driver.f90
)
target_link_libraries(test_hamiltonian_time_driver.x GORILLA netcdff netcdf lapack)
add_test(NAME hamiltonian_time_oracle COMMAND test_hamiltonian_time_driver.x)

add_executable(test_analytic_circ_flux_driver.x
  SRC/TESTS/test_analytic_circ_flux_driver.f90
)
target_link_libraries(test_analytic_circ_flux_driver.x GORILLA netcdff netcdf lapack)
add_test(NAME analytic_circ_flux_oracle COMMAND test_analytic_circ_flux_driver.x
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

add_executable(test_analytic_circ_mesh_driver.x
  SRC/TESTS/test_analytic_circ_mesh_driver.f90
)
target_link_libraries(test_analytic_circ_mesh_driver.x GORILLA netcdff netcdf lapack)
foreach(cs 1 2)
  foreach(b0sign 1 -1)
    set(mesh_dir ${CMAKE_CURRENT_BINARY_DIR}/analytic_circ_mesh_cs${cs}_b${b0sign})
    file(MAKE_DIRECTORY ${mesh_dir})
    add_test(NAME analytic_circ_mesh_cs${cs}_b${b0sign}
      COMMAND test_analytic_circ_mesh_driver.x ${cs} ${b0sign}
      WORKING_DIRECTORY ${mesh_dir})
  endforeach()
endforeach()

message(${CMAKE_HOST_SYSTEM_PROCESSOR})

# Exact species contract used by external namelist writers and native table readers.
find_package(Python3 REQUIRED COMPONENTS Interpreter)
add_executable(test_marker_species.x SRC/TESTS/test_marker_species.f90)
target_link_libraries(test_marker_species.x GORILLA netcdff netcdf lapack)
add_test(NAME marker_species
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/SRC/TESTS/test_marker_species.py
    $<TARGET_FILE:test_marker_species.x> ${CMAKE_CURRENT_SOURCE_DIR}/INPUT/gorilla.inp)
