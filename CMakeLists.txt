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
		GIT_TAG 3ec46e07c6c2b6ea1725af0ac5189971880ff34b
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
endif()

message(${CMAKE_HOST_SYSTEM_PROCESSOR})
