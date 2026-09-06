#!/bin/bash
set -euo pipefail
export LINES=40 COLUMNS=120
audit_dir=$(cd -- "$(dirname -- "$0")" && pwd)
repo_dir=$(git -C "$audit_dir" rev-parse --show-toplevel)
build_dir="$audit_dir/work/full-build"
mkdir -p "$build_dir"
module purge
module load intel_comp/2024.2.0
module load compiler-rt tbb compiler
cd "$build_dir"
# Only sources are linked: never reuse objects from the main working directory.
for source in "$repo_dir"/*.f90 "$repo_dir"/*.h; do
    ln -sfn "$source" "${source##*/}"
done
ifx --version
make -f "$repo_dir/makefile" -j1 PMP2mod_tools.o PMP2mod_MGbackground.o
make -f "$repo_dir/makefile" -j1 PMP2init PMP2start PMP2main PMP2BDM
