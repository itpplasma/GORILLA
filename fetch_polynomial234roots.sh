#!/usr/bin/env bash
set -euo pipefail

repo_root="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
archive_path="${repo_root}/954.zip"
extract_dir="${repo_root}/954"
source_path="${extract_dir}/F90/Src/Polynomial234RootSolvers.f90"

if [[ ! -f "${archive_path}" ]]; then
    wget -O "${archive_path}" "https://netlib.org/toms/954.zip"
fi

rm -rf "${extract_dir}"
unzip -q -o "${archive_path}" -d "${repo_root}"

if [[ ! -f "${source_path}" ]]; then
    echo "downloaded archive does not contain ${source_path}" >&2
    exit 1
fi

echo "downloaded external Polynomial234RootSolvers.f90 to ${source_path}"
echo "the tracked placeholder stays untouched; CMake will use this external file when present"
