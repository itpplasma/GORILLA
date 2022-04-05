cmake_minimum_required(VERSION 3.12)

project (GORILLA
  VERSION 1.0.0)
  
enable_language(Fortran)

if(UNIX AND NOT APPLE)
	include_directories(/usr/include)
elseif(APPLE)
	include_directories(/opt/local/include)
endif()

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(netCDF REQUIRED)

if ($ENV{GORILLA_COVERAGE} STREQUAL "TRUE")
	find_package(PFUNIT REQUIRED)
	enable_testing()
endif()

add_compile_options(-g -fbacktrace -ffpe-trap=zero,overflow,invalid -fbounds-check -fopenmp)
add_link_options(-g -fbacktrace -ffpe-trap=zero,overflow,invalid -fbounds-check -fopenmp)

if(APPLE)
	add_link_options(-L/opt/local/lib)
endif()

if ($ENV{GORILLA_COVERAGE} STREQUAL "TRUE")
	add_link_options(--coverage)
endif()

add_subdirectory(SRC)

