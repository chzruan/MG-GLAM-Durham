#!/bin/bash
# V1 layout unit test for the low-k mode dump. Builds in a scratch dir against the
# repository's module objects (build the code first: make PMP2main), runs in seconds.
set -u
HERE=$(cd "$(dirname "$0")" && pwd); ROOT=$(cd "$HERE/../.." && pwd)
W=${TMPDIR:-/tmp}/glam_modes_test_$$; mkdir -p "$W"; trap 'rm -rf "$W"' EXIT
if [ -z "${SKIP_MODULES:-}" ] && command -v module >/dev/null 2>&1; then
    module purge 2>/dev/null; module load intel_comp/2024.2.0 2>/dev/null; module load compiler-rt tbb compiler mpi 2>/dev/null
fi
FFLAGS="-O2 -qopenmp -march=core-avx2 -mcmodel=medium -convert big_endian -I$ROOT"
ifx $FFLAGS -o "$W/test_modes.exe" "$HERE/test_modes.f90" "$ROOT"/PMP2mod_tools.o "$ROOT"/PMP2mod_fft5.o \
    "$ROOT"/PMP2mod_random.o "$ROOT"/PMP2MG_subroutines.o "$ROOT"/PMP2mod_MGbackground.o "$ROOT"/PMP2mod_density.o "$ROOT"/PMP2mod_extras.o -module "$W" || exit 2
( cd "$W" && OMP_NUM_THREADS=2 ./test_modes.exe ) || exit 2
python3 "$HERE/check_modes.py" "$W/vtest.0001.0001.bin"
