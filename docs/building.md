# Building GORILLA

GORILLA is built with CMake. The top-level `Makefile` is a thin wrapper that
configures a Ninja build in `BUILD/`.

## Requirements

- GNU Fortran (the only supported compiler)
- CMake (3.12 or newer) and Ninja
- NetCDF-C and NetCDF-Fortran
- BLAS and LAPACK

Ubuntu:

```bash
sudo apt install git gfortran make cmake ninja-build liblapack-dev libnetcdff-dev
```

macOS with [Homebrew](https://brew.sh):

```bash
brew install git gcc make cmake ninja netcdf netcdf-fortran libomp
```

On Apple silicon, `CMakeLists.txt` adds `/opt/homebrew/include` and
`/opt/homebrew/lib` to the search paths.

## Building with make

```bash
make                    # Release build in BUILD/
make CONFIG=Debug       # other CMake build type (on first configuration)
make reconfigure CONFIG=Debug
make test               # ctest in BUILD/
make clean              # remove BUILD/
```

The executable is `BUILD/test_gorilla_main.x`. The example folders link to
this path.

## Building with CMake directly

```bash
cmake -S . -B BUILD
cmake --build BUILD
```

Any generator works; with the default Unix Makefiles generator no Ninja is
needed.

## Tests

`make test` runs the CTest checks registered in `CMakeLists.txt`:
`boozer_chartmap_oracle` and `radial_electric_projection_oracle`.

The unit tests in `SRC/TESTS` use [pFUnit](https://github.com/Goddard-Fortran-Ecosystem/pFUnit).
They are built only when the environment variable `GORILLA_COVERAGE=TRUE` is
set, which also enables coverage instrumentation. Each test source contains a
description of what it checks. The tests focus on low-level functions with
trial or default inputs compared against expected results. Tests that cannot
yet be completed because of the current program structure are disabled but
kept in `SRC/TESTS`.

## Coverage

Install pFUnit and [lcov](https://github.com/linux-test-project/lcov)
(`sudo apt install lcov` or `brew install lcov`), point `PFUNIT_DIR` to the
pFUnit installation and run

```bash
./build_coverage.sh
```

The script rebuilds `BUILD/` with `GORILLA_COVERAGE=TRUE`, runs `ctest` and
writes an HTML report to `BUILD/COVERAGE/index.html`.

## Third-party code

- `SRC/contrib/rkf45.f90`: Runge-Kutta-Fehlberg routines from
  <https://people.sc.fsu.edu/~jburkardt/f_src/rkf45/rkf45.html> (GNU LGPL).
- `SRC/contrib/cmplx_roots_sg.f90`: the complex polynomial root solver of
  Skowron and Gould (Apache 2.0, see
  [SRC/contrib/LICENSE.cmplx_roots_sg](../SRC/contrib/LICENSE.cmplx_roots_sg)),
  used for the quadratic, cubic and quartic roots in the polynomial pusher.
  J. Skowron and A. Gould, "General Complex Polynomial Root Solver and Its
  Further Optimization for Binary Microlenses", arXiv:[1203.1034](https://arxiv.org/abs/1203.1034) (2012).

Both are part of the repository; no download is needed.

## CINECA Leonardo

The following recipe builds GORILLA on the Leonardo CPU login nodes with GCC
12.2.0, CMake 3.27.9, OpenBLAS 0.3.26, NetCDF-C 4.9.2 and NetCDF-Fortran
4.6.1. Two workarounds are needed: the Spack NetCDF modules ship no
`netCDFConfig.cmake`, and `CMakeLists.txt` links `-llapack`, which OpenBLAS
does not provide as a separate library. Both are replaced by small shims.

```bash
module purge
module load profile/base gcc/12.2.0 cmake/3.27.9 ninja/1.11.1 \
    openblas/0.3.26--gcc--12.2.0 python/3.11.7 \
    netcdf-c/4.9.2--gcc--12.2.0-spack0.22 \
    netcdf-fortran/4.6.1--gcc--12.2.0-spack0.22
export B=$CINECA_SCRATCH/builds TMPDIR=$CINECA_SCRATCH/builds/tmp
mkdir -p $TMPDIR $B/shim/netCDF $B/shim/lib
export CMAKE_BUILD_PARALLEL_LEVEL=4 OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1
echo 'set(netCDF_FOUND TRUE)' > $B/shim/netCDF/netCDFConfig.cmake
for l in lapack blas; do ln -sf $OPENBLAS_LIB/libopenblas.so $B/shim/lib/lib$l.so; done
LINK="-L$B/shim/lib -Wl,-rpath,$OPENBLAS_LIB -L$NETCDF_FORTRAN_LIB -L$NETCDF_C_LIB -Wl,-rpath,$NETCDF_FORTRAN_LIB -Wl,-rpath,$NETCDF_C_LIB"
INC="-I$NETCDF_FORTRAN_INC -I$NETCDF_C_INC"
cd $B
git clone https://github.com/itpplasma/GORILLA.git && cd GORILLA
nice cmake -S . -B build-lrd -G Ninja -DCMAKE_BUILD_TYPE=Release \
    -DnetCDF_DIR=$B/shim/netCDF "-DCMAKE_Fortran_FLAGS=$INC" \
    "-DCMAKE_EXE_LINKER_FLAGS=$LINK" "-DCMAKE_SHARED_LINKER_FLAGS=$LINK"
nice cmake --build build-lrd -j 4
(cd build-lrd && ctest)
```

Examples can be run on the serial partition, e.g. in a batch job with
`--partition=lrd_all_serial --cpus-per-task=4`:

```bash
cd EXAMPLES/example_1 && ln -sf ../../build-lrd/test_gorilla_main.x . && ./test_gorilla_main.x
```
