cmake_minimum_required(VERSION 3.12)

project (GORILLA
  VERSION 1.0.0)

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

message(${CMAKE_HOST_SYSTEM_PROCESSOR})

# JOREK reader smoke test against the committed circular fixture in
# SRC/TESTS/jorek_data (see the README there for provenance). This block is a
# no-op unless the JOREK backend infrastructure (option GORILLA_ENABLE_JOREK
# plus the libneo FetchContent block providing LIBNEO::jorek) is present in
# this source tree and enabled.
if(GORILLA_ENABLE_JOREK AND TARGET LIBNEO::jorek)
	enable_testing()
	add_executable(test_jorek_restart_smoke.x
		SRC/TESTS/test_jorek_restart_smoke.f90)
	target_link_libraries(test_jorek_restart_smoke.x GORILLA LIBNEO::jorek)
	add_test(NAME test_jorek_restart_smoke
		COMMAND test_jorek_restart_smoke.x
			"${CMAKE_CURRENT_SOURCE_DIR}/SRC/TESTS/jorek_data/circular_equilibrium_restart.h5")
endif()
