# Local polynomial kernels

The checked-in Fortran module and CUDA header come from the same FortSym
expressions in `generate_cell_kernels.f90`. Ordinary GORILLA builds compile the
checked-in sources and do not require FortSym.

With FortSym installed, build and check the generator separately:

```sh
cmake -S gen -B /tmp/gorilla-generator -Dfortsym_DIR=/path/to/fortsym/lib/cmake/fortsym
cmake --build /tmp/gorilla-generator
/tmp/gorilla-generator/generate_cell_kernels SRC/generated
/tmp/gorilla-generator/generate_gpu_physics SRC/generated
sh gen/check_generated.sh /tmp/gorilla-generator/generate_cell_kernels \
  /tmp/gorilla-generator/generate_gpu_physics
```

`gorilla_affine_coefficients2` produces the constant, linear and quadratic
coefficients of the affine system `dz/dtau = A*z + b`. Here `z` contains the
three local position coordinates and parallel velocity, `tau` is the evolution
parameter, `A` is the constant four-by-four system matrix, and `b` is the
constant four-component forcing vector.
`gorilla_poly_evaluate2` evaluates these coefficients. Both languages use
Fortran column-major storage: `A[row + 4*column]` and
`coef[row + 4*power]` with zero-based C++ indices.

The Fortran `gorilla_poly_step1` through `gorilla_poly_step4` accept cached
matrix-vector products. They preserve the pusher's sequential accumulation of
Taylor terms, including when the cached products differ from fresh products of
the current state. Input and output arrays must not overlap.

The native `cell_generated_oracle` CTest checks exact nilpotent trajectories,
dense affine coefficients, matrix Taylor series and independent cached products.

`generate_gpu_physics.f90` derives the affine operator and four face-distance
polynomials used by the CUDA pusher. It emits both Fortran and CUDA from the
same FortSym expressions. The regeneration check compares both output pairs;
it requires no GPU or CUDA compiler.
