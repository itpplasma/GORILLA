# Local polynomial kernels

The checked-in Fortran module and CUDA header come from the same FortSym
expressions in `generate_cell_kernels.f90`. Ordinary GORILLA builds compile the
checked-in sources and do not require FortSym.

With FortSym installed, build and check the generator separately:

```sh
cmake -S gen -B /tmp/gorilla-generator -Dfortsym_DIR=/path/to/fortsym/lib/cmake/fortsym
cmake --build /tmp/gorilla-generator
/tmp/gorilla-generator/generate_cell_kernels SRC/generated
sh gen/check_generated.sh /tmp/gorilla-generator/generate_cell_kernels
```

`gorilla_affine_coefficients2` produces the constant, linear and quadratic
coefficients of the four-component affine system `dz/dtau = A*z + b`.
`gorilla_poly_evaluate2` evaluates these coefficients. Both languages use
Fortran column-major storage: `A[row + 4*column]` and
`coef[row + 4*power]` with zero-based C++ indices.

The Fortran `gorilla_poly_step1` through `gorilla_poly_step4` accept cached
matrix-vector products. They preserve the pusher's sequential accumulation of
Taylor terms, including when the cached products differ from fresh products of
the current state. Input and output arrays must not overlap.

The native `cell_generated_oracle` CTest checks exact nilpotent trajectories,
dense affine coefficients, matrix Taylor series and independent cached products.
