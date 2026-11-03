#include "gorilla_cuda.h"
#include <cuda_runtime.h>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <vector>

namespace {
constexpr double pi = 3.141592653589793238462643383279502884;

void require(bool condition, const char *message) {
  if (!condition) {
    std::fprintf(stderr, "CUDA poly2: %s\n", message);
    std::exit(1);
  }
}

void close(double actual, double expected, const char *message) {
  require(std::isfinite(actual) && std::abs(actual - expected) <
      2.0e-12*(1.0 + std::abs(expected)), message);
}

struct mesh {
  int64_t n;
  std::vector<double> values;
  std::vector<int64_t> topology;
  explicit mesh(int64_t count) : n(count), values(59*count), topology(16*count, -1) {}
  double &at(int cell, int field) { return values[field*n + cell]; }
  int64_t &neighbor(int cell, int field) { return topology[field*n + cell]; }
  void simplex(int cell) {
    at(cell, 3) = -1.0;
    at(cell, 4) = -1.0;
    at(cell, 5) = -1.0;
    at(cell, 6) = 1.0;
    at(cell, 10) = 1.0;
    at(cell, 14) = 1.0;
    at(cell, 15) = 1.0;
    at(cell, 16) = 1.0;
    at(cell, 17) = 1.0;
    at(cell, 56) = 1.0;
    for (int field = 4; field < 16; ++field) neighbor(cell, field) = 0;
  }
  void *create(int coordinates = 1) {
    void *handle = nullptr;
    require(gorilla_cuda_create(n, values.data(), topology.data(), coordinates,
        1, 1, 1.0, 1.0, 2.2204460492503131e-16, &handle) == NGORILLA_ERR_OK,
        "create");
    require(handle != nullptr, "missing handle after successful create");
    return handle;
  }
};

struct marker {
  std::array<double, 9> real{0.1, 0.1, 0.1, 1.0, 0.0, 0.2, 0.1, 0.1, 0.1};
  std::array<int64_t, 4> integer{1, 0, 1, NGORILLA_STATUS_ACTIVE};
  int64_t crossings = 0;
  void advance(void *handle, int64_t budget = 64) {
    double profile[5];
    require(gorilla_cuda_advance(handle, 1, real.data(), integer.data(),
        &crossings, budget, profile) == NGORILLA_ERR_OK, "advance");
    for (double duration : profile)
      require(std::isfinite(duration) && duration >= 0.0, "invalid phase duration");
  }
};

void destroy(void *&handle) {
  require(gorilla_cuda_destroy(&handle) == NGORILLA_ERR_OK, "destroy");
  require(handle == nullptr, "destroy must clear the handle");
  require(gorilla_cuda_destroy(&handle) == NGORILLA_ERR_OK, "repeat destroy");
}

void partial_and_loss() {
  mesh cells(1);
  cells.simplex(0);
  void *handle = cells.create();
  marker partial;
  partial.advance(handle);
  require(partial.integer[3] == NGORILLA_STATUS_COMPLETE, "partial status");
  close(partial.real[0], 0.3, "partial position");
  close(partial.real[5], 0.0, "partial remaining time");
  require(partial.integer[1] == 0 && partial.crossings == 0, "partial face/count");
  marker loss;
  loss.real[5] = 1.0;
  loss.advance(handle);
  require(loss.integer[3] == NGORILLA_STATUS_PHYSICAL_LOSS, "boundary status");
  require(loss.integer[0] == -1 && loss.integer[2] == 1, "last cell at boundary");
  close(loss.real[0], 0.8, "boundary position");
  close(loss.real[5], 0.3, "boundary remaining time");
  require(loss.crossings == 1, "boundary crossing count");
  marker backwards;
  backwards.real[0] = 0.5;
  backwards.real[5] = -0.2;
  backwards.advance(handle);
  require(backwards.integer[3] == NGORILLA_STATUS_COMPLETE, "backward status");
  close(backwards.real[0], 0.3, "backward position");
  double profile[5];
  require(gorilla_cuda_advance(handle, 0, nullptr, nullptr, nullptr, 1, profile) ==
      NGORILLA_ERR_OK, "empty batch");
  require(gorilla_cuda_advance(handle, -1, nullptr, nullptr, nullptr, 1, profile) ==
      NGORILLA_ERR_BAD_ARG, "negative batch length");
  destroy(handle);
}

void acceleration_and_handles() {
  mesh constant(1), accelerated(1);
  constant.simplex(0);
  accelerated.simplex(0);
  accelerated.at(0, 22) = -2.0;
  accelerated.at(0, 35) = -2.0;
  void *first = constant.create(), *second = accelerated.create();
  marker a, b;
  b.real[3] = 0.0;
  b.advance(second);
  a.advance(first);
  require(b.integer[3] == NGORILLA_STATUS_COMPLETE, "acceleration status");
  close(b.real[0], 0.14, "constant acceleration position");
  close(b.real[3], 0.4, "constant acceleration velocity");
  close(a.real[0], 0.3, "independent handles");
  destroy(first);
  b = marker{};
  b.real[3] = 0.0;
  b.advance(second);
  close(b.real[0], 0.14, "surviving independent handle");
  destroy(second);
}

void two_cells_and_budget() {
  mesh cells(2);
  cells.simplex(0);
  cells.simplex(1);
  cells.neighbor(0, 0) = 2;
  cells.neighbor(0, 4) = 2;
  cells.at(1, 0) = 1.0;
  const double normals[12] = {1, -1, -1, 1, 1, 1, -1, -1, 1, -1, 1, -1};
  for (int k = 0; k < 12; ++k) cells.at(1, 3 + k) = normals[k];
  cells.at(1, 15) = 2.0;
  void *handle = cells.create();
  marker m;
  m.real[2] = 0.2;
  m.real[5] = 0.65;
  m.advance(handle, 1);
  require(m.integer[3] == NGORILLA_STATUS_ACTIVE, "budget must leave active marker");
  require(m.integer[0] == 2 && m.integer[1] == 2, "neighbor topology");
  close(m.real[0], 0.7, "first-cell exit");
  close(m.real[5], 0.05, "budget remaining time");
  m.advance(handle, 1);
  require(m.integer[3] == NGORILLA_STATUS_COMPLETE, "budget resume status");
  require(m.integer[0] == 2 && m.integer[2] == 2, "resume last cell");
  close(m.real[0], 0.75, "second-cell final position");
  close(m.real[6], -0.25, "second-cell local position");
  require(m.crossings == 1, "resumed crossing count");
  destroy(handle);
}

void periodic_handover(int coordinates, bool theta = false) {
  mesh cells(1);
  const int phi = coordinates == 1 || theta ? 1 : 2;
  const int other = phi == 1 ? 2 : 1;
  cells.at(0, 3 + phi) = -1;
  cells.at(0, 6 + phi) = 1;
  cells.at(0, 9) = 1;
  cells.at(0, 12 + other) = 1;
  cells.at(0, 15) = 2*pi;
  cells.at(0, 16) = 1;
  cells.at(0, 17) = 1;
  cells.at(0, 56 + phi) = 1;
  for (int field = 0; field < 16; ++field) cells.neighbor(0, field) = 0;
  cells.neighbor(0, 0) = 1;
  cells.neighbor(0, 4) = 2;
  cells.neighbor(0, theta ? 12 : 8) = 1;
  void *handle = cells.create(coordinates);
  marker m;
  m.real[5] = 2*pi + 0.1;
  m.advance(handle, 1);
  require(m.integer[3] == NGORILLA_STATUS_ACTIVE, "periodic budget status");
  close(m.real[phi], 0.0, "periodic position exchange");
  close(m.real[6 + phi], 2*pi, "local coordinates before periodic exchange");
  m.advance(handle, 1);
  require(m.integer[3] == NGORILLA_STATUS_COMPLETE, "periodic resume");
  close(m.real[phi], 0.2, "periodic final position");
  destroy(handle);
}

void recovery_and_invalid_state() {
  mesh cells(1);
  cells.simplex(0);
  void *handle = cells.create();
  marker stationary;
  stationary.real[3] = 0.0;
  const auto entry = stationary.real;
  stationary.advance(handle);
  require(stationary.integer[3] == NGORILLA_STATUS_REPLAY_CELL, "no root needs replay");
  require(stationary.real == entry && stationary.crossings == 0,
      "replay must preserve cell-entry state");
  marker invalid;
  invalid.integer[0] = 2;
  invalid.advance(handle);
  require(invalid.integer[3] == NGORILLA_STATUS_NUMERICAL_FAILURE,
      "invalid cell is numerical failure");
  marker outside;
  outside.real[1] = -0.1;
  const auto outside_entry = outside.real;
  outside.advance(handle);
  require(outside.integer[3] == NGORILLA_STATUS_REPLAY_CELL &&
      outside.real == outside_entry, "failed geometry must preserve entry state");
  marker inward;
  inward.real[0] = 1.2;
  inward.real[3] = -1.0;
  const auto inward_entry = inward.real;
  inward.advance(handle);
  require(inward.integer[3] == NGORILLA_STATUS_REPLAY_CELL &&
      inward.real == inward_entry, "inward crossing must preserve entry state");
  destroy(handle);
}

void batched_layout() {
  mesh cells(1);
  cells.simplex(0);
  void *handle = cells.create();
  constexpr int n = 5;
  std::array<double, 9*n> real;
  std::array<int64_t, 4*n> integer;
  std::array<int64_t, n> crossings{};
  const double times[n] = {0.1, -0.1, 1.0, 0.0, 0.4};
  for (int j = 0; j < n; ++j) {
    marker m;
    m.real[0] = 0.2;
    m.real[5] = times[j];
    if (j == 4) m.integer[3] = NGORILLA_STATUS_COMPLETE;
    for (int k = 0; k < 9; ++k) real[k*n + j] = m.real[k];
    for (int k = 0; k < 4; ++k) integer[k*n + j] = m.integer[k];
  }
  double profile[5];
  require(gorilla_cuda_advance(handle, n, real.data(), integer.data(),
      crossings.data(), 64, profile) == NGORILLA_ERR_OK, "batch advance");
  close(real[0], 0.3, "batched forward marker");
  close(real[1], 0.1, "batched backward marker");
  close(real[2], 0.8, "batched loss marker");
  close(real[3], 0.2, "batched zero-time marker");
  close(real[4], 0.2, "batched terminal marker");
  close(real[5*n + 4], 0.4, "terminal marker retains time");
  for (int j = 0; j < n; ++j) {
    require(integer[3*n + j] == (j == 2 ? NGORILLA_STATUS_PHYSICAL_LOSS :
        NGORILLA_STATUS_COMPLETE), "batched status");
    require(crossings[j] == (j == 2 ? 1 : 0), "batched crossing count");
  }
  destroy(handle);
}
} // namespace

int main() {
  int devices = 0;
  if (cudaGetDeviceCount(&devices) != cudaSuccess || devices == 0) return 77;
  partial_and_loss();
  acceleration_and_handles();
  two_cells_and_budget();
  periodic_handover(1);
  periodic_handover(2);
  periodic_handover(2, true);
  recovery_and_invalid_state();
  batched_layout();
  std::puts("CUDA poly2 analytic, recovery, periodic and lifecycle checks passed");
  return 0;
}
