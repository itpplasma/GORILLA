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

include(CTest)
if(BUILD_TESTING)
    add_executable(test_analytic_circular_flux.x TESTS/test_analytic_circular_flux.f90)
    target_link_libraries(test_analytic_circular_flux.x GORILLA netcdff netcdf lapack)
    add_test(NAME analytic_circular_flux COMMAND test_analytic_circular_flux.x)
endif()

if(BUILD_TESTING)
    add_executable(test_analytic_poloidal_flux.x
        SRC/field_analytic_circ_mod.f90
        TESTS/test_analytic_poloidal_flux.f90)
    set_target_properties(test_analytic_poloidal_flux.x PROPERTIES
        Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/poloidal_test_modules)
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/poloidal_test_modules)
    target_include_directories(test_analytic_poloidal_flux.x PRIVATE
        ${CMAKE_CURRENT_BINARY_DIR}/poloidal_test_modules)
    add_test(NAME analytic_poloidal_flux COMMAND test_analytic_poloidal_flux.x)
endif()

if(BUILD_TESTING)
    add_executable(test_hamiltonian_time.x
        SRC/hamiltonian_time_coefficients_mod.f90
        TESTS/test_hamiltonian_time.f90)
    set_target_properties(test_hamiltonian_time.x PROPERTIES
        Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/time_test_modules)
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/time_test_modules)
    target_include_directories(test_hamiltonian_time.x PRIVATE
        ${CMAKE_CURRENT_BINARY_DIR}/time_test_modules)
    add_test(NAME hamiltonian_time COMMAND test_hamiltonian_time.x)
endif()
