cmake_minimum_required(VERSION 3.22)

project (GORILLA
  VERSION 1.0.0)

option(GORILLA_ENABLE_JOREK "Build the JOREK model-303 field backend" OFF)

enable_language(Fortran)

include_directories(/usr/local/include)
if(UNIX AND NOT APPLE)
	include_directories(/usr/include)
elseif(APPLE AND ${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "arm64")
	include_directories(/opt/homebrew/include)
endif()

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(netCDF REQUIRED)

if(GORILLA_ENABLE_JOREK)
	include(FetchContent)
	set(LIBNEO_BUILD_TESTING OFF CACHE BOOL "" FORCE)
	set(LIBNEO_BUILD_NUMERICS OFF CACHE BOOL "" FORCE)
	set(LIBNEO_BUILD_JOREK ON CACHE BOOL "" FORCE)
	FetchContent_Declare(libneo
		GIT_REPOSITORY https://github.com/itpplasma/libneo.git
		GIT_TAG 7b07b8121dff371bea67e657f298f900ce9df496
	)
	FetchContent_MakeAvailable(libneo)
endif()

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

if(GORILLA_ENABLE_JOREK)
	enable_testing()
	add_executable(test_jorek_field_backend.x
		SRC/TESTS/test_jorek_field_backend.f90)
	target_link_libraries(test_jorek_field_backend.x GORILLA LIBNEO::jorek)
	add_test(NAME test_jorek_field_backend COMMAND test_jorek_field_backend.x)
	add_executable(test_jorek_mesh.x SRC/TESTS/test_jorek_mesh.f90)
	target_link_libraries(test_jorek_mesh.x GORILLA LIBNEO::jorek netcdff netcdf lapack)
	add_test(NAME test_jorek_mesh COMMAND test_jorek_mesh.x)
	set(GORILLA_JOREK_GOLDEN_RESTART "" CACHE FILEPATH
		"JOREK restart with SHA256 ac152cd345cd233664b6fa659e1000200519158ba385c2acab09494929bc315c")
	if(GORILLA_JOREK_GOLDEN_RESTART)
		if(NOT EXISTS "${GORILLA_JOREK_GOLDEN_RESTART}")
			message(FATAL_ERROR "GORILLA_JOREK_GOLDEN_RESTART does not exist")
		endif()
		add_executable(test_jorek_restart_backend.x
			SRC/TESTS/test_jorek_restart_backend.f90)
		target_link_libraries(test_jorek_restart_backend.x GORILLA LIBNEO::jorek)
		add_test(NAME test_jorek_restart_backend
			COMMAND test_jorek_restart_backend.x "${GORILLA_JOREK_GOLDEN_RESTART}")
		set_tests_properties(test_jorek_restart_backend PROPERTIES
			LABELS "jorek;golden")
		add_executable(test_jorek_restart_mesh.x
			SRC/TESTS/test_jorek_restart_mesh.f90)
		target_link_libraries(test_jorek_restart_mesh.x GORILLA LIBNEO::jorek
			netcdff netcdf lapack)
		add_test(NAME test_jorek_restart_mesh
			COMMAND test_jorek_restart_mesh.x "${GORILLA_JOREK_GOLDEN_RESTART}")
		set_tests_properties(test_jorek_restart_mesh PROPERTIES
			DEPENDS test_jorek_restart_backend LABELS "jorek;golden")
		add_executable(test_jorek_poloidal_linearization.x
			SRC/TESTS/test_jorek_poloidal_linearization.f90)
		target_link_libraries(test_jorek_poloidal_linearization.x
			GORILLA LIBNEO::jorek netcdff netcdf lapack)
		add_test(NAME test_jorek_poloidal_linearization
			COMMAND test_jorek_poloidal_linearization.x
				"${GORILLA_JOREK_GOLDEN_RESTART}"
				"${CMAKE_CURRENT_BINARY_DIR}/jorek_poloidal_element_metrics.csv")
		set_tests_properties(test_jorek_poloidal_linearization PROPERTIES
			DEPENDS test_jorek_restart_mesh LABELS "jorek;golden")
		add_executable(test_jorek_refined_plane.x
			SRC/TESTS/test_jorek_refined_plane.f90)
		target_link_libraries(test_jorek_refined_plane.x
			GORILLA LIBNEO::jorek netcdff netcdf lapack)
		add_test(NAME test_jorek_refined_plane
			COMMAND test_jorek_refined_plane.x
				"${GORILLA_JOREK_GOLDEN_RESTART}")
		set_tests_properties(test_jorek_refined_plane PROPERTIES
			DEPENDS test_jorek_restart_mesh LABELS "jorek;golden")
		set(JOREK_ORBIT_DIR "${CMAKE_CURRENT_BINARY_DIR}/jorek_snapshot_orbit")
		file(MAKE_DIRECTORY "${JOREK_ORBIT_DIR}")
		configure_file(SRC/TESTS/jorek_snapshot_orbit/gorilla.inp.in
			"${JOREK_ORBIT_DIR}/gorilla.inp" @ONLY)
		foreach(INPUT_FILE tetra_grid.inp gorilla_plot.inp orbit_start_rphizlambda.dat)
			configure_file("SRC/TESTS/jorek_snapshot_orbit/${INPUT_FILE}"
				"${JOREK_ORBIT_DIR}/${INPUT_FILE}" COPYONLY)
		endforeach()
		add_executable(test_jorek_orbit_output.x
			SRC/TESTS/test_jorek_orbit_output.f90)
		add_test(NAME test_jorek_snapshot_orbit
			COMMAND ${CMAKE_COMMAND}
				-DGORILLA_MAIN=$<TARGET_FILE:test_gorilla_main.x>
				-DORBIT_VALIDATOR=$<TARGET_FILE:test_jorek_orbit_output.x>
				-DORBIT_DIR=${JOREK_ORBIT_DIR}
				-P ${CMAKE_CURRENT_SOURCE_DIR}/SRC/TESTS/run_jorek_snapshot_orbit.cmake)
		set_tests_properties(test_jorek_snapshot_orbit PROPERTIES
			DEPENDS test_jorek_restart_backend LABELS "jorek;golden")
		set(JOREK_DERIVED_ORBIT_DIR
			"${CMAKE_CURRENT_BINARY_DIR}/jorek_derived_orbit")
		file(MAKE_DIRECTORY "${JOREK_DERIVED_ORBIT_DIR}")
		configure_file(SRC/TESTS/jorek_snapshot_orbit/gorilla.inp.in
			"${JOREK_DERIVED_ORBIT_DIR}/gorilla.inp" @ONLY)
		configure_file(SRC/TESTS/jorek_derived_orbit/tetra_grid.inp
			"${JOREK_DERIVED_ORBIT_DIR}/tetra_grid.inp" COPYONLY)
		foreach(INPUT_FILE gorilla_plot.inp orbit_start_rphizlambda.dat)
			configure_file("SRC/TESTS/jorek_snapshot_orbit/${INPUT_FILE}"
				"${JOREK_DERIVED_ORBIT_DIR}/${INPUT_FILE}" COPYONLY)
		endforeach()
		add_test(NAME test_jorek_derived_mesh_orbit
			COMMAND ${CMAKE_COMMAND}
				-DGORILLA_MAIN=$<TARGET_FILE:test_gorilla_main.x>
				-DORBIT_VALIDATOR=$<TARGET_FILE:test_jorek_orbit_output.x>
				-DORBIT_DIR=${JOREK_DERIVED_ORBIT_DIR}
				-DEXPECTED_RECORDS=99
				-DMIN_PPHI_RESPONSE=1.0e-7
				-P ${CMAKE_CURRENT_SOURCE_DIR}/SRC/TESTS/run_jorek_snapshot_orbit.cmake)
		set_tests_properties(test_jorek_derived_mesh_orbit PROPERTIES
			DEPENDS test_jorek_restart_mesh LABELS "jorek;golden;slow")
		set(JOREK_REFINED_ORBIT_DIR
			"${CMAKE_CURRENT_BINARY_DIR}/jorek_refined_orbit")
		file(MAKE_DIRECTORY "${JOREK_REFINED_ORBIT_DIR}")
		configure_file(SRC/TESTS/jorek_snapshot_orbit/gorilla.inp.in
			"${JOREK_REFINED_ORBIT_DIR}/gorilla.inp" @ONLY)
		configure_file(SRC/TESTS/jorek_refined_orbit/tetra_grid.inp
			"${JOREK_REFINED_ORBIT_DIR}/tetra_grid.inp" COPYONLY)
		foreach(INPUT_FILE gorilla_plot.inp orbit_start_rphizlambda.dat)
			configure_file("SRC/TESTS/jorek_snapshot_orbit/${INPUT_FILE}"
				"${JOREK_REFINED_ORBIT_DIR}/${INPUT_FILE}" COPYONLY)
		endforeach()
		add_test(NAME test_jorek_refined_mesh_orbit
			COMMAND ${CMAKE_COMMAND}
				-DGORILLA_MAIN=$<TARGET_FILE:test_gorilla_main.x>
				-DORBIT_VALIDATOR=$<TARGET_FILE:test_jorek_orbit_output.x>
				-DORBIT_DIR=${JOREK_REFINED_ORBIT_DIR}
				-DEXPECTED_RECORDS=246
				-DMIN_PPHI_RESPONSE=1.0e-7
				-P ${CMAKE_CURRENT_SOURCE_DIR}/SRC/TESTS/run_jorek_snapshot_orbit.cmake)
		set_tests_properties(test_jorek_refined_mesh_orbit PROPERTIES
			DEPENDS test_jorek_refined_plane LABELS "jorek;golden;slow")
	endif()
endif()

message(${CMAKE_HOST_SYSTEM_PROCESSOR})
