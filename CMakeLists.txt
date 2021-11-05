cmake_minimum_required(VERSION 3.12)

project (GORILLA
  VERSION 1.0.0
  LANGUAGES Fortran)

set(CMAKE_VERBOSE_MAKEFILE ON)

include_directories(/usr/include)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(netCDF REQUIRED)

if ($ENV{GORILLA_COVERAGE} STREQUAL "TRUE")
	find_package(PFUNIT REQUIRED)
	enable_testing()
endif()

add_compile_options(-g -fbacktrace -ffpe-trap=zero,overflow,invalid -fbounds-check -fopenmp)
add_link_options(-g -fbacktrace -ffpe-trap=zero,overflow,invalid -fbounds-check -fopenmp)

if ($ENV{GORILLA_COVERAGE} STREQUAL "TRUE")
	add_link_options(--coverage)
endif()

add_subdirectory(SRC)

