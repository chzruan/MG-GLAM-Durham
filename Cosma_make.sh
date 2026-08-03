#!/bin/bash
# Build MG-GLAM on Cosma8 with the Intel oneAPI toolchain.
# Usage: source Cosma_make.sh [make targets...]
#   Must be *sourced* (not executed) so `module load` affects your shell.
#   Defaults to: PMP2init PMP2start PMP2MG glam2gadget
#   The first three targets are the full simulation pipeline (see e.g.
#   symmetron/prepare/sanity/*/submit.sh) and between them cover every object
#   file in the makefile's $(OBJ) list plus the two extra entry points
#   (PMP2init.o, PMP2start.o). glam2gadget (the PM->Gadget converter used by
#   ic2gadget/final2gadget) is a separate entry point not covered by $(OBJ),
#   so it's listed explicitly -- otherwise a no-args `source Cosma_make.sh`
#   silently stops rebuilding it.
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
    _cosma_make_targets=(PMP2init PMP2start PMP2MG glam2gadget)
fi

make "${_cosma_make_targets[@]}"
