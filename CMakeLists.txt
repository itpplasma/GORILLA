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
find_package(OpenMP REQUIRED COMPONENTS Fortran)

if ("$ENV{GORILLA_COVERAGE}" STREQUAL "TRUE")
	find_package(PFUNIT REQUIRED)
	enable_testing()
endif()

# Keep performance builds free of runtime checks unless explicitly requested.
if(CMAKE_BUILD_TYPE MATCHES "^(Release|MinSizeRel)$")
  set(_runtime_checks_default OFF)
else()
  set(_runtime_checks_default ON)
endif()
option(GORILLA_RUNTIME_CHECKS
  "Fortran bounds checks, FPE traps and backtraces" ${_runtime_checks_default})

if(GORILLA_RUNTIME_CHECKS)
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    add_compile_options(
      "$<$<COMPILE_LANGUAGE:Fortran>:-g;-fbacktrace;-ffpe-trap=zero,overflow,invalid;-fbounds-check>")
  elseif(CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC|PGI")
    add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-g;-C;-Ktrap=inv,divz,ovf>")
  else()
    message(FATAL_ERROR "GORILLA_RUNTIME_CHECKS is not supported for "
                        "${CMAKE_Fortran_COMPILER_ID}; use a supported compiler "
                        "or set GORILLA_RUNTIME_CHECKS=OFF.")
  endif()
endif()

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

add_executable(test_cell_generated_driver.x
  SRC/TESTS/test_cell_generated_driver.f90
)
target_link_libraries(test_cell_generated_driver.x GORILLA)
add_test(NAME cell_generated_oracle COMMAND test_cell_generated_driver.x)

message(${CMAKE_HOST_SYSTEM_PROCESSOR})
