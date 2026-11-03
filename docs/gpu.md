# Polynomial CUDA batches

GORILLA's optional CUDA backend advances an ensemble through the quadratic
polynomial pusher. The tetrahedral mesh stays on the device until the batch
handle is finalized. Marker state is transferred for each bounded device
advance, and uncertain cell crossings are replayed with the existing CPU
pusher. End-to-end timings therefore include transfers and CPU recovery.

The supported settings are `ipusher=2`, `poly_order=2`, `i_precomp=0`,
`i_time_tracing_option=1`, and `handover_processing_kind=1`, with adaptive time
steps, strong electric fields, and optional trajectory quantities disabled.
Both coordinate systems are supported. Other settings are rejected during
initialization. The batch uses the species selected by the ordinary GORILLA
initialization; it does not select a different particle mass or charge.

## Build

The default build remains CPU only and does not require a CUDA toolkit.
To enable CUDA, use CMake 3.18 or newer, a CUDA toolkit supporting the target
GPU, and the usual Fortran, BLAS, LAPACK, and NetCDF dependencies:

```sh
cmake -S . -B build-gpu -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_Fortran_COMPILER=gfortran -DGORILLA_ENABLE_CUDA=ON \
  -DCMAKE_CUDA_ARCHITECTURES=75
cmake --build build-gpu --parallel
ctest --test-dir build-gpu --output-on-failure
```

Choose the architecture for the GPU that will run the executable: for example,
`75` for Turing or `120` for Blackwell, with a sufficiently recent toolkit.
With CMake 3.24 or newer, `native` selects the build host's GPU. Use an explicit
architecture when compiling on a host without a GPU. NetCDF's Fortran modules
must be compatible with the chosen Fortran compiler.

`GORILLA_RUNTIME_CHECKS=ON` enables the existing host bounds and floating-point
checks for diagnostic builds. CUDA qualification requires an accessible device;
the low-level CUDA oracle reports a missing device as a skipped test.

## Batch interface

The module `orbit_timestep_gorilla_batch_mod` provides
`initialize_gorilla_batch`, `orbit_timestep_gorilla_batch`, and
`finalize_gorilla_batch`. Initialize the usual GORILLA mesh and physics first,
then initialize a `gorilla_batch_t` with backend `cpu` or `cuda`. The batch call
accepts positions shaped `(3,n)`, parallel and perpendicular velocity vectors,
a shared signed time step, and each marker's initialization, tetrahedron, and
entry-face state. Reuse the handle between calls and finalize it when finished.
Reinitialize it after rebuilding the mesh or changing the species or settings.

Each marker returns `GORILLA_COMPLETE`, `GORILLA_PHYSICAL_LOSS`, or
`GORILLA_NUMERICAL_FAILURE`, defined in `gorilla_status_mod`. A rejected device
crossing starts CPU recovery at its unchanged cell-entry state. Recovery
preserves the perpendicular invariant and remaining signed time; a failed CPU
solve remains a numerical failure. Consumers must remove physically lost
markers from later batches and must not count numerical failures as physical
losses or successful trajectories.

The optional `gorilla_batch_profile_t` output separates setup, upload, kernel,
download, CPU replay, and total host times. Both backends count completed
markers, physical losses, numerical failures, and cell crossings; CUDA also
counts replayed cells. Setup is reported once per initialized handle.
Compare complete application wall times
on the same host, with identical inputs, numerical settings, and CPU allocation;
kernel time alone does not measure an application speedup.
