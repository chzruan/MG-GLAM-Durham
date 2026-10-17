#!/bin/bash
# Unit test for module Subsample (PMP2mod_subsample.f90 + glam_hash.c, item I9
# of docs/degrace_companion_implementation_plan.md).
#
# Compiles glam_hash.c, PMP2mod_tools.f90, PMP2mod_subsample.f90 and
# test_subsample.f90 in a work directory (nothing is written to the
# repository root), runs the 32^3 toy test with 1 and 4 OpenMP threads,
# then checks:
#   - bitwise identity of every Subsample file between 1 and 4 threads
#     (deterministic two-pass compaction)
#   - bitwise identity of step 7 (SubChunk = 1000) and step 9 (default
#     chunk): the output does not depend on the chunk size
#   - hash_ref.c: block routines agree with glam_select per ID
#   - check_subsample.py: numpy splitmix64 == C on 10^7 IDs, decoded
#     positions/velocities vs the raw arrays, n_sel vs N_p/M, M = 1 complete
# Light: a few seconds, < 200 MB, 4 threads, under nice.
#
# Usage:  bash run_subsample_test.sh [workdir]      (default: test/extras/work)
#         SKIP_MODULES=1 bash run_subsample_test.sh  (Intel env already loaded)
#         PYTHON="python3" bash run_subsample_test.sh (default: micromamba run -n cosemu python3)
set -u

HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../.." && pwd)
WORK=${1:-$HERE/work}
PYTHON=${PYTHON:-"micromamba run -n cosemu python3"}
NHASH=${NHASH:-10000000}
# makefile FFLAGS without the opt-report switches (they only produce .optrpt files)
FFLAGS="-O3 -g -traceback -ftz -unroll -qopenmp -march=core-avx2 -mfma -fp-model fast=1 -shared-intel -mcmodel=medium -convert big_endian"
CFLAGS="-O2 -g -march=core-avx2 -mcmodel=medium"
NICE="nice -n 19"

# ---------------------------------------------------------------- environment
if [ -z "${SKIP_MODULES:-}" ] && command -v module >/dev/null 2>&1; then
    module purge 2>/dev/null
    module load intel_comp/2024.2.0 2>/dev/null
    export I_MPI_F90=ifx
    module load compiler-rt 2>/dev/null
    module load tbb 2>/dev/null
    module load compiler 2>/dev/null
    module load mpi 2>/dev/null
fi
command -v ifx >/dev/null 2>&1 || { echo "FATAL: ifx not found (load intel_comp/2024.2.0 + compiler)"; exit 2; }
command -v icx >/dev/null 2>&1 || { echo "FATAL: icx not found"; exit 2; }

nfail=0
report() {   # report <0|1> <label>
    if [ "$1" -eq 0 ]; then echo "PASS  $2"; else echo "FAIL  $2"; nfail=$((nfail + 1)); fi
}

# ---------------------------------------------------------------- build
mkdir -p "$WORK" || exit 2
cd "$WORK" || exit 2
rm -rf t1 t4 hash_ref.bin ./*.o ./*.mod ./*.exe
echo "== building in $WORK"
$NICE icx $CFLAGS -c "$ROOT/glam_hash.c"                         || { echo "FATAL: glam_hash.c"; exit 2; }
$NICE ifx $FFLAGS -w -c "$ROOT/PMP2mod_tools.f90"                || { echo "FATAL: PMP2mod_tools.f90"; exit 2; }
$NICE ifx $FFLAGS -c "$ROOT/PMP2mod_subsample.f90"               || { echo "FATAL: PMP2mod_subsample.f90"; exit 2; }
$NICE ifx $FFLAGS -c "$HERE/test_subsample.f90"                  || { echo "FATAL: test_subsample.f90"; exit 2; }
$NICE ifx $FFLAGS -o test_subsample.exe PMP2mod_tools.o glam_hash.o PMP2mod_subsample.o test_subsample.o \
                                                                 || { echo "FATAL: link test_subsample"; exit 2; }
$NICE icx $CFLAGS -o hash_ref.exe "$HERE/hash_ref.c" glam_hash.o || { echo "FATAL: hash_ref.c"; exit 2; }
report 0 "build: glam_hash.c, PMP2mod_subsample.f90, test_subsample.f90, hash_ref.c"

# ---------------------------------------------------------------- run
echo "== hash reference ($NHASH IDs)"
$NICE ./hash_ref.exe hash_ref.bin "$NHASH" > hash_ref.log 2>&1
report $? "hash_ref: block routines agree with glam_select ($(tail -1 hash_ref.log))"

for nt in 1 4; do
    mkdir -p "t$nt"
    ( cd "t$nt" && OMP_NUM_THREADS=$nt $NICE ../test_subsample.exe > run.log 2>&1 )
    report $? "test_subsample.exe with $nt thread(s)"
done
grep '^   subsample:' t4/run.log

# ---------------------------------------------------------------- determinism
ok=0
for f in t4/Subsample.*; do
    cmp -s "$f" "t1/$(basename "$f")" || { ok=1; echo "   differs: $(basename "$f")"; }
done
[ "$(ls t4/Subsample.* 2>/dev/null | wc -l)" -ge 8 ] || ok=1
report $ok "1 vs 4 threads: all Subsample files bitwise identical"

ok=0
cmp -s t4/Subsample.0007.0001.bin t4/Subsample.0009.0001.bin || ok=1
diff -q <(grep -v '^step' t4/Subsample.0007.0001.meta) <(grep -v '^step' t4/Subsample.0009.0001.meta) >/dev/null || ok=1
report $ok "SubChunk 1000 vs default: step 7 and step 9 files identical"

# ---------------------------------------------------------------- reference decoder
echo "== check_subsample.py"
$PYTHON -I "$HERE/check_subsample.py" "$WORK/t4" "$WORK/hash_ref.bin"
report $? "check_subsample.py (all checks above)"

# ---------------------------------------------------------------- cleanup (keep logs and the small outputs)
rm -f hash_ref.bin ./*.o ./*.mod
echo
if [ "$nfail" -eq 0 ]; then
    echo "ALL PASS  (outputs and logs in $WORK)"
    exit 0
else
    echo "$nfail FAILURE(S)  (see $WORK)"
    exit 1
fi
