#!/bin/sh
set -eu
if [ "$#" -lt 1 ] || [ "$#" -gt 2 ]; then
    echo "usage: $0 /path/to/generate_cell_kernels [/path/to/generate_gpu_physics]" >&2
    exit 2
fi
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
output=$(mktemp -d)
trap 'rm -rf "$output"' EXIT HUP INT TERM
"$1" "$output"
diff -u "$root/SRC/generated/gorilla_cell_generated.f90" "$output/gorilla_cell_generated.f90"
diff -u "$root/SRC/generated/gorilla_cell_generated.cuh" "$output/gorilla_cell_generated.cuh"
if [ "$#" -eq 2 ]; then
    "$2" "$output"
    diff -u "$root/SRC/generated/gorilla_cuda_physics.f90" "$output/gorilla_cuda_physics.f90"
    diff -u "$root/SRC/generated/gorilla_cuda_physics.cuh" "$output/gorilla_cuda_physics.cuh"
fi
