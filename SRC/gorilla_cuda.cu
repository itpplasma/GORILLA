#include "gorilla_cuda.h"
#include "generated/gorilla_cell_generated.cuh"
#include "generated/gorilla_cuda_physics.cuh"

#include <cuda_runtime.h>
#include <cfloat>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <new>

namespace {
constexpr int active = NGORILLA_STATUS_ACTIVE;
constexpr int complete = NGORILLA_STATUS_COMPLETE;
constexpr int physical_loss = NGORILLA_STATUS_PHYSICAL_LOSS;
constexpr int replay_cell = NGORILLA_STATUS_REPLAY_CELL;
constexpr int numerical_failure = NGORILLA_STATUS_NUMERICAL_FAILURE;
constexpr int mesh_fields = NGORILLA_MESH_NFIELD;
constexpr int topology_fields = NGORILLA_TOPO_NFIELD;
constexpr int real_fields = NGORILLA_RSTATE_NFIELD;
constexpr int integer_fields = NGORILLA_ISTATE_NFIELD;
constexpr double pi = 3.141592653589793238462643383279502884;
using clock_type = std::chrono::steady_clock;

struct device_mesh {
  const double *values;
  const int64_t *topology;
  int64_t count;
  int coordinates, periods, sign_sqg;
  double cm_over_e, clight, eps;
};

struct mesh_handle {
  device_mesh mesh{};
  cudaStream_t stream = nullptr;
  cudaEvent_t event[4]{};
  double *mesh_values = nullptr, *state = nullptr;
  int64_t *topology = nullptr, *indices = nullptr;
  int64_t *crossings = nullptr;
  int64_t capacity = 0;
  int device = -1;
};

double milliseconds(clock_type::time_point start) {
  return std::chrono::duration<double, std::milli>(clock_type::now() - start).count();
}

bool valid_size(int64_t n, int fields, size_t element_size) {
  return n >= 0 && static_cast<size_t>(n) <=
      std::numeric_limits<size_t>::max() / element_size / fields;
}

int abi_error(cudaError_t error) {
  return error == cudaSuccess ? NGORILLA_ERR_OK : NGORILLA_ERR_LAUNCH;
}

cudaError_t release(mesh_handle *handle) {
  if (!handle) return cudaSuccess;
  cudaError_t result = cudaSuccess;
  auto save_error = [&result](cudaError_t error) {
    if (result == cudaSuccess) result = error;
  };
  save_error(cudaFree(handle->state));
  save_error(cudaFree(handle->indices));
  save_error(cudaFree(handle->crossings));
  save_error(cudaFree(handle->mesh_values));
  save_error(cudaFree(handle->topology));
  for (auto event : handle->event)
    if (event) save_error(cudaEventDestroy(event));
  if (handle->stream) save_error(cudaStreamDestroy(handle->stream));
  delete handle;
  return result;
}

cudaError_t resize_particles(mesh_handle &handle, int64_t n) {
  if (n <= handle.capacity) return cudaSuccess;
  double *state = nullptr;
  int64_t *indices = nullptr;
  int64_t *crossings = nullptr;
  cudaError_t error = cudaMalloc(&state, sizeof(double) * real_fields * size_t(n));
  if (error == cudaSuccess)
    error = cudaMalloc(&indices, sizeof(int64_t) * integer_fields * size_t(n));
  if (error == cudaSuccess)
    error = cudaMalloc(&crossings, sizeof(int64_t) * size_t(n));
  if (error != cudaSuccess) {
    cudaFree(state);
    cudaFree(indices);
    cudaFree(crossings);
    return error;
  }
  cudaFree(handle.state);
  cudaFree(handle.indices);
  cudaFree(handle.crossings);
  handle.state = state;
  handle.indices = indices;
  handle.crossings = crossings;
  handle.capacity = n;
  return cudaSuccess;
}

// This branch structure is the first-attempt solver in pusher_tetra_poly.f90.
__device__ double quadratic_root(double a, double b, double c, double eps) {
  double root = DBL_MAX, discriminant, denominator;
  if (c > 0.0) {
    if (a > 0.0) {
      if (b < 0.0) {
        discriminant = b*b - 2.0*a*c;
        if (discriminant > 0.0) {
          denominator = -b + sqrt(discriminant);
          root = fabs(denominator) > eps ? 2.0*c/denominator
                                         : (-sqrt(discriminant) - b)/a;
        } else if (discriminant == 0.0) root = -b/a;
      }
    } else if (a < 0.0) {
      discriminant = b*b - 2.0*a*c;
      denominator = -b + sqrt(discriminant);
      root = fabs(denominator) > eps ? 2.0*c/denominator
                                     : (-sqrt(discriminant) - b)/a;
    } else if (b < 0.0) root = -c/b;
  } else if (c < 0.0) {
    if (a < 0.0) {
      if (b > 0.0) {
        discriminant = b*b - 2.0*a*c;
        if (discriminant > 0.0) root = (sqrt(discriminant) - b)/a;
        else if (discriminant == 0.0) root = -b/a;
      }
    } else if (a > 0.0) {
      discriminant = b*b - 2.0*a*c;
      root = (sqrt(discriminant) - b)/a;
    } else if (b > 0.0) root = -c/b;
  } else if ((a > 0.0 && b < 0.0) || (a < 0.0 && b > 0.0)) {
    root = -2.0*b/a;
  }
  return root;
}

__device__ void cell_system(const double *m, const double *z, double perpinv,
                            const device_mesh &mesh, double sign_rhs,
                            double *matrix, double *affine) {
  gorilla_cell_system(z[0], z[1], z[2], z[3], perpinv, mesh.cm_over_e,
      mesh.clight, sign_rhs, m[17], m[18], m[19], m[20], m[21], m[22],
      m[23], m[24], m[25], m[26], m[27], m[28], m[29], m[30], m[31],
      m[32], m[33], m[34], m[35], m[36], m[37], m[38], m[39], m[40],
      m[41], m[42], m[43], m[44], m[45], m[46], m[47], m[48], m[49],
      m[50], m[51], m[52], m[53], m[54], m[55], m[56], m[57], m[58],
      &matrix[0], &matrix[1], &matrix[2], &matrix[3], &matrix[4], &matrix[5],
      &matrix[6], &matrix[7], &matrix[8], &matrix[9], &matrix[10], &matrix[11],
      &matrix[12], &matrix[13], &matrix[14], &matrix[15],
      &affine[0], &affine[1], &affine[2], &affine[3]);
}

__device__ double plane_distance(const double *cell, int face, const double *z) {
  double distance;
  const int offset = 3 + 3*face;
  gorilla_plane_distance(cell[offset], cell[offset + 1], cell[offset + 2],
      face == 0 ? cell[15] : 0.0, z[0], z[1], z[2], &distance);
  return distance;
}

__global__ void advance_cells(device_mesh mesh, int64_t n, double *state,
                              int64_t *indices, int64_t *crossings,
                              int64_t max_crossings) {
  const int64_t marker = int64_t(blockIdx.x)*blockDim.x + threadIdx.x;
  if (marker >= n || indices[3*size_t(n) + marker] != active) return;
  for (int64_t crossing = 0; crossing < max_crossings; ++crossing) {
    const int64_t tetra = indices[marker], entry_face = indices[size_t(n) + marker];
    if (tetra < 1 || tetra > mesh.count || entry_face < 0 || entry_face > 4) {
      indices[3*size_t(n) + marker] = numerical_failure;
      return;
    }
    double cell[mesh_fields], z[4], matrix[16], affine[4], coefficient[12];
    for (int field = 0; field < mesh_fields; ++field)
      cell[field] = mesh.values[field*size_t(mesh.count) + tetra - 1];
    for (int k = 0; k < 3; ++k) z[k] = state[k*size_t(n) + marker] - cell[k];
    z[3] = state[3*size_t(n) + marker];
    const double remaining = state[5*size_t(n) + marker];
    if (remaining == 0.0) {
      indices[2*size_t(n) + marker] = tetra;
      for (int k = 0; k < 3; ++k) state[(6 + k)*size_t(n) + marker] = z[k];
      indices[3*size_t(n) + marker] = complete;
      return;
    }
    const double sign_rhs = mesh.sign_sqg*(remaining < 0.0 ? -1.0 : 1.0);
    const double time_scale = sign_rhs*cell[16];
    if (!isfinite(remaining) || !isfinite(time_scale) || time_scale == 0.0) {
      indices[3*size_t(n) + marker] = replay_cell;
      return;
    }
    cell_system(cell, z, state[4*size_t(n) + marker], mesh, sign_rhs,
                matrix, affine);
    gorilla_affine_coefficients2(matrix, affine, z, coefficient);
    for (double value : coefficient) {
      if (!isfinite(value)) {
        indices[3*size_t(n) + marker] = replay_cell;
        return;
      }
    }
    double tau = DBL_MAX;
    int exit_face = -1;
    for (int face = 0; face < 4; ++face) {
      const int offset = 3 + 3*face;
      double distance, linear, quadratic, candidate = DBL_MAX;
      gorilla_face_coefficients2(cell[offset], cell[offset + 1], cell[offset + 2],
          face == 0 ? cell[15] : 0.0, z[0], z[1], z[2],
          coefficient[4], coefficient[5], coefficient[6],
          coefficient[8], coefficient[9], coefficient[10],
          &distance, &linear, &quadratic);
      if (face + 1 == entry_face || distance == 0.0) {
        if (quadratic != 0.0) candidate = -linear/(quadratic/2.0);
      } else if (quadratic == 0.0) {
        if (linear != 0.0) candidate = -distance/linear;
      } else candidate = quadratic_root(quadratic, linear, distance, mesh.eps);
      if (isfinite(candidate) && candidate > 0.0 && candidate < tau) {
        tau = candidate;
        exit_face = face;
      }
    }
    if (exit_face < 0) {
      indices[3*size_t(n) + marker] = replay_cell;
      return;
    }
    double endpoint[4];
    gorilla_poly_evaluate2(coefficient, tau, endpoint);
    bool valid = true;
    for (double value : endpoint) valid = valid && isfinite(value);
    for (int face = 0; face < 4; ++face) {
      const double distance = plane_distance(cell, face, endpoint);
      valid = valid && isfinite(distance) &&
          (face == exit_face ? fabs(distance) <= 1.0e-11 : distance >= 0.0);
    }
    double velocity;
    const int normal_offset = 3 + 3*exit_face;
    gorilla_normal_velocity2(cell[normal_offset], cell[normal_offset + 1],
        cell[normal_offset + 2], coefficient[4], coefficient[5], coefficient[6],
        coefficient[8], coefficient[9], coefficient[10], tau, &velocity);
    valid = valid && isfinite(velocity) && velocity <= 0.0;
    double passed = tau*time_scale;
    valid = valid && isfinite(passed) && passed != 0.0 &&
        ((remaining < 0.0) == (passed < 0.0));
    if (!valid) {
      indices[3*size_t(n) + marker] = replay_cell;
      return;
    }
    const bool finished = fabs(passed) >= fabs(remaining);
    if (finished) {
      gorilla_poly_evaluate2(coefficient, remaining/time_scale, endpoint);
      for (double value : endpoint) valid = valid && isfinite(value);
      for (int face = 0; face < 4; ++face) {
        const double distance = plane_distance(cell, face, endpoint);
        valid = valid && isfinite(distance) && distance >= 0.0;
      }
      if (!valid) {
        indices[3*size_t(n) + marker] = replay_cell;
        return;
      }
      passed = remaining;
    }
    int64_t neighbor = tetra, next_face = 0, phi_shift = 0, theta_shift = 0;
    if (!finished) {
      neighbor = mesh.topology[exit_face*size_t(mesh.count) + tetra - 1];
      next_face = mesh.topology[(4 + exit_face)*size_t(mesh.count) + tetra - 1];
      phi_shift = mesh.topology[(8 + exit_face)*size_t(mesh.count) + tetra - 1];
      theta_shift = mesh.topology[(12 + exit_face)*size_t(mesh.count) + tetra - 1];
      if (neighbor == 0 || neighbor > mesh.count ||
          (neighbor > 0 && (next_face < 1 || next_face > 4))) {
        indices[3*size_t(n) + marker] = replay_cell;
        return;
      }
    }
    for (int k = 0; k < 3; ++k) {
      state[k*size_t(n) + marker] = endpoint[k] + cell[k];
      state[(6 + k)*size_t(n) + marker] = endpoint[k];
    }
    state[3*size_t(n) + marker] = endpoint[3];
    state[5*size_t(n) + marker] = remaining - passed;
    indices[2*size_t(n) + marker] = tetra;
    indices[marker] = neighbor;
    indices[size_t(n) + marker] = next_face;
    if (finished) {
      indices[3*size_t(n) + marker] = complete;
      return;
    }
    if (phi_shift == 1 || phi_shift == -1) {
      const int phi_index = mesh.coordinates == 1 ? 1 : 2;
      state[phi_index*size_t(n) + marker] -= phi_shift*(2.0*pi/mesh.periods);
    }
    if (mesh.coordinates == 2 && (theta_shift == 1 || theta_shift == -1))
      state[size_t(n) + marker] -= theta_shift*(2.0*pi);
    ++crossings[marker];
    if (neighbor < 0) {
      indices[3*size_t(n) + marker] = physical_loss;
      return;
    }
  }
}
} // namespace

extern "C" int gorilla_cuda_create(int64_t ntetr, const double *mesh_values,
    const int64_t *topology, int coordinates, int periods, int sign_sqg,
    double cm_over_e, double clight, double eps, void **output) {
  if (!output) return NGORILLA_ERR_BAD_ARG;
  *output = nullptr;
  if (ntetr < 1 || !mesh_values || !topology ||
      !valid_size(ntetr, mesh_fields, sizeof(double)) ||
      !valid_size(ntetr, topology_fields, sizeof(int64_t)) ||
      (coordinates != 1 && coordinates != 2) || periods < 1 ||
      (sign_sqg != 1 && sign_sqg != -1) || !std::isfinite(cm_over_e) ||
      cm_over_e == 0.0 || !std::isfinite(clight) || clight <= 0.0 ||
      !std::isfinite(eps) || eps < 0.0) return NGORILLA_ERR_BAD_ARG;
  auto *handle = new (std::nothrow) mesh_handle;
  if (!handle) return NGORILLA_ERR_LAUNCH;
  cudaError_t error = cudaGetDevice(&handle->device);
  if (error == cudaSuccess)
    error = cudaStreamCreateWithFlags(&handle->stream, cudaStreamNonBlocking);
  for (int i = 0; i < 4 && error == cudaSuccess; ++i)
    error = cudaEventCreate(&handle->event[i]);
  const size_t real_bytes = sizeof(double)*mesh_fields*size_t(ntetr);
  const size_t integer_bytes = sizeof(int64_t)*topology_fields*size_t(ntetr);
  if (error == cudaSuccess) error = cudaMalloc(&handle->mesh_values, real_bytes);
  if (error == cudaSuccess) error = cudaMalloc(&handle->topology, integer_bytes);
  if (error == cudaSuccess)
    error = cudaMemcpyAsync(handle->mesh_values, mesh_values, real_bytes,
                            cudaMemcpyHostToDevice, handle->stream);
  if (error == cudaSuccess)
    error = cudaMemcpyAsync(handle->topology, topology, integer_bytes,
                            cudaMemcpyHostToDevice, handle->stream);
  if (error == cudaSuccess) error = cudaStreamSynchronize(handle->stream);
  if (error != cudaSuccess) {
    release(handle);
    return abi_error(error);
  }
  handle->mesh = {handle->mesh_values, handle->topology, ntetr, coordinates,
                  periods, sign_sqg, cm_over_e, clight, eps};
  *output = handle;
  return 0;
}

extern "C" int gorilla_cuda_advance(void *opaque, int64_t n, double *state,
    int64_t *indices, int64_t *crossings, int64_t max_crossings, double profile_ms[5]) {
  const auto start = clock_type::now();
  if (profile_ms) for (int i = 0; i < 5; ++i) profile_ms[i] = 0.0;
  if (!opaque || !profile_ms || n < 0 || max_crossings < 1 ||
      !valid_size(n, real_fields, sizeof(double)) ||
      !valid_size(n, integer_fields, sizeof(int64_t)) ||
      n > int64_t(std::numeric_limits<int>::max())) return NGORILLA_ERR_BAD_ARG;
  if (n == 0) return 0;
  if (!state || !indices || !crossings) return NGORILLA_ERR_BAD_ARG;
  auto &handle = *static_cast<mesh_handle *>(opaque);
  int device;
  cudaError_t error = cudaGetDevice(&device);
  if (error != cudaSuccess) return abi_error(error);
  if (device != handle.device) return NGORILLA_ERR_BAD_ARG;
  const auto allocation_start = clock_type::now();
  error = resize_particles(handle, n);
  profile_ms[0] = milliseconds(allocation_start);
  const size_t real_bytes = sizeof(double)*real_fields*size_t(n);
  const size_t integer_bytes = sizeof(int64_t)*integer_fields*size_t(n);
  const size_t count_bytes = sizeof(int64_t)*size_t(n);
  if (error == cudaSuccess) error = cudaEventRecord(handle.event[0], handle.stream);
  if (error == cudaSuccess)
    error = cudaMemcpyAsync(handle.state, state, real_bytes,
                            cudaMemcpyHostToDevice, handle.stream);
  if (error == cudaSuccess)
    error = cudaMemcpyAsync(handle.indices, indices, integer_bytes,
                            cudaMemcpyHostToDevice, handle.stream);
  if (error == cudaSuccess)
    error = cudaMemcpyAsync(handle.crossings, crossings, count_bytes,
                            cudaMemcpyHostToDevice, handle.stream);
  if (error == cudaSuccess) error = cudaEventRecord(handle.event[1], handle.stream);
  if (error == cudaSuccess) {
    advance_cells<<<(n - 1)/128 + 1, 128, 0, handle.stream>>>(handle.mesh,
        n, handle.state, handle.indices, handle.crossings, max_crossings);
    error = cudaGetLastError();
  }
  if (error == cudaSuccess) error = cudaEventRecord(handle.event[2], handle.stream);
  if (error == cudaSuccess)
    error = cudaMemcpyAsync(state, handle.state, real_bytes,
                            cudaMemcpyDeviceToHost, handle.stream);
  if (error == cudaSuccess)
    error = cudaMemcpyAsync(indices, handle.indices, integer_bytes,
                            cudaMemcpyDeviceToHost, handle.stream);
  if (error == cudaSuccess)
    error = cudaMemcpyAsync(crossings, handle.crossings, count_bytes,
                            cudaMemcpyDeviceToHost, handle.stream);
  if (error == cudaSuccess) error = cudaEventRecord(handle.event[3], handle.stream);
  const cudaError_t completion = cudaStreamSynchronize(handle.stream);
  if (error == cudaSuccess) error = completion;
  if (error == cudaSuccess) {
    for (int phase = 1; phase <= 3; ++phase) {
      float elapsed;
      error = cudaEventElapsedTime(&elapsed, handle.event[phase - 1], handle.event[phase]);
      if (error != cudaSuccess) break;
      profile_ms[phase] = elapsed;
    }
  }
  profile_ms[4] = milliseconds(start);
  return abi_error(error);
}

extern "C" int gorilla_cuda_destroy(void **opaque) {
  if (!opaque || !*opaque) return 0;
  auto *handle = static_cast<mesh_handle *>(*opaque);
  int device;
  cudaError_t error = cudaGetDevice(&device);
  if (error != cudaSuccess) return abi_error(error);
  if (device != handle->device) return NGORILLA_ERR_BAD_ARG;
  const int result = abi_error(release(handle));
  *opaque = nullptr;
  return result;
}
