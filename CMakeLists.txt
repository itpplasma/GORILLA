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

# On macOS, FindBLAS defaults to the Accelerate framework. Its vectorised
# kernels compute on masked-off SIMD lanes and discard the results, which
# raises spurious FP exceptions. Combined with the -ffpe-trap flags set below
# that aborts with SIGFPE inside libBLAS.dylib before any physics runs.
# Prefer OpenBLAS instead; set -DBLA_VENDOR explicitly to override.
if(APPLE AND NOT DEFINED BLA_VENDOR)
	set(BLA_VENDOR OpenBLAS)
	# Homebrew keeps openblas keg-only, so it is not on the default search
	# path. Ask brew where it lives; this resolves both the Intel
	# (/usr/local) and Apple silicon (/opt/homebrew) prefixes.
	find_program(BREW_EXECUTABLE brew)
	if(BREW_EXECUTABLE)
		execute_process(COMMAND ${BREW_EXECUTABLE} --prefix openblas
			OUTPUT_VARIABLE OPENBLAS_PREFIX
			OUTPUT_STRIP_TRAILING_WHITESPACE
			ERROR_QUIET)
		if(OPENBLAS_PREFIX AND IS_DIRECTORY "${OPENBLAS_PREFIX}")
			list(APPEND CMAKE_PREFIX_PATH "${OPENBLAS_PREFIX}")
		endif()
	endif()
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

target_link_libraries(test_gorilla_main.x GORILLA netcdff netcdf ${LAPACK_LIBRARIES})

add_executable(test_boozer_chartmap_driver.x
  SRC/TESTS/test_boozer_chartmap_driver.f90
)
target_link_libraries(test_boozer_chartmap_driver.x GORILLA netcdff netcdf ${LAPACK_LIBRARIES})
add_test(NAME boozer_chartmap_oracle COMMAND test_boozer_chartmap_driver.x)

add_executable(test_radial_electric_projection_driver.x
  SRC/TESTS/test_radial_electric_projection_driver.f90
)
target_link_libraries(test_radial_electric_projection_driver.x
  GORILLA netcdff netcdf ${LAPACK_LIBRARIES})
add_test(NAME radial_electric_projection_oracle
  COMMAND test_radial_electric_projection_driver.x)

message(${CMAKE_HOST_SYSTEM_PROCESSOR})
