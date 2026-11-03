#!/bin/sh
set -eu
if [ "$#" -ne 1 ]; then
    echo "usage: $0 /path/to/generate_cell_kernels" >&2
    exit 2
fi
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
output=$(mktemp -d)
trap 'rm -rf "$output"' EXIT HUP INT TERM
"$1" "$output"
diff -u "$root/SRC/generated/gorilla_cell_generated.f90" "$output/gorilla_cell_generated.f90"
diff -u "$root/SRC/generated/gorilla_cell_generated.cuh" "$output/gorilla_cell_generated.cuh"
