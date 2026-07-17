#!/bin/bash
# Build MG-GLAM on Cosma8 with the Intel oneAPI toolchain.
# Usage: source Cosma_make.sh [make targets...]
#   Must be *sourced* (not executed) so `module load` affects your shell.
#   Defaults to: PMP2start PMP2MG
#
# Deliberately no `set -e`: this runs in your interactive shell via `source`,
# so an `exit` on error would close your whole SSH session instead of just
# stopping the script. Each risky step is checked explicitly with `return`
# (safe under `source`) instead.

_cosma_make_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)" || return 1
cd "$_cosma_make_dir" || return 1

module purge
module load intel_comp/2024.2.0 || return 1
export I_MPI_F90=ifx
module load compiler-rt tbb compiler mpi || return 1

# Build order matters here: PMP2mod_tools.o produces Tools.mod, which almost
# every other source file USEs (including PMP2mod_MGbackground.f90). Build it
# first, then PMP2mod_MGbackground.o (-> ExtradofBackgroundData.mod, needed by
# density/extradof/csf/kmf), before letting make handle the rest in whatever
# order it likes.
make PMP2mod_tools.o || return 1
make PMP2mod_MGbackground.o || return 1

_cosma_make_targets=("$@")
if [ ${#_cosma_make_targets[@]} -eq 0 ]; then
    _cosma_make_targets=(PMP2start PMP2MG)
fi

make "${_cosma_make_targets[@]}"
