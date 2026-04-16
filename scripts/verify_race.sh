#!/bin/bash
# Task A driver: race-detection regression on MG relaxation.
#
# Builds PMP2main-bitmatch.exe (pre-Pass-4 flags: -O2, no AVX2/FMA/fast=1),
# then runs a short (Nsteps=3) simulation from the fid_F5/Run11 checkpoint
# at OMP_NUM_THREADS ∈ {1, 16, 64, 128}, capturing stdout to
# verify/race_T{N}.log. Finally invokes scripts/compare_residuals.py which
# compares per-V-cycle residuals across thread counts (T=1 reference).
#
# Pass bar: max|Δres| ≤ 1e-14 (exact bit-match expected at bitmatch flags).
#
# Run from the repo root on a Cosma8 compute node with the Intel toolchain
# already loaded:
#   module purge && module load intel_comp/2024.2.0
#   export I_MPI_F90=ifx
#   module load compiler-rt tbb compiler mpi
#   bash scripts/verify_race.sh

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"

VERIFY_DIR="$REPO_ROOT/verify"
RUN_DIR="$REPO_ROOT/fid_F5/Run11"
IC_SNAPSHOT="$RUN_DIR/.race_ic_snapshot"
NSTEPS=3
THREAD_COUNTS=(1 16 64 128)

# Files that define the starting state — restored before each T run.
IC_FILES=(PMcrd.DAT pt.dat fort.17)

mkdir -p "$VERIFY_DIR"

echo "=== Task A: MG race-detection regression ==="
echo "Repo:     $REPO_ROOT"
echo "Run dir:  $RUN_DIR"
echo "Logs:     $VERIFY_DIR/race_T{N}.log"
echo "Nsteps:   $NSTEPS"
echo "Threads:  ${THREAD_COUNTS[*]}"
echo

# ---- 1. Build the bitmatch binary ------------------------------------------
echo "--- Building PMP2main-bitmatch.exe ---"
make PMP2main-bitmatch
if [[ ! -x "$REPO_ROOT/PMP2main-bitmatch.exe" ]]; then
    echo "ERROR: PMP2main-bitmatch.exe not produced" >&2
    exit 2
fi
echo

# ---- 2. Snapshot IC state from Run11 (once) --------------------------------
if [[ ! -d "$IC_SNAPSHOT" ]]; then
    echo "--- Snapshotting IC state from $RUN_DIR ---"
    mkdir -p "$IC_SNAPSHOT"
    for f in "${IC_FILES[@]}"; do
        if [[ ! -f "$RUN_DIR/$f" ]]; then
            echo "ERROR: missing IC file $RUN_DIR/$f" >&2
            exit 2
        fi
        cp -p "$RUN_DIR/$f" "$IC_SNAPSHOT/$f"
    done
    echo "Snapshot: $IC_SNAPSHOT"
    echo
fi

restore_ic() {
    for f in "${IC_FILES[@]}"; do
        cp -p "$IC_SNAPSHOT/$f" "$RUN_DIR/$f"
    done
}

# ---- 3. Run at each OMP_NUM_THREADS ----------------------------------------
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export OMP_STACKSIZE=512M

for T in "${THREAD_COUNTS[@]}"; do
    LOG="$VERIFY_DIR/race_T${T}.log"
    echo "--- Run T=$T -> $LOG ---"
    restore_ic
    pushd "$RUN_DIR" >/dev/null
    OMP_NUM_THREADS="$T" \
        "$REPO_ROOT/PMP2main-bitmatch.exe" <<<"$NSTEPS" \
        2>&1 | tee "$LOG" >/dev/null
    popd >/dev/null
    nres=$(grep -c 'The residual is equal to' "$LOG" || true)
    echo "    V-cycles logged: $nres"
done
echo

# ---- 4. Compare --------------------------------------------------------------
echo "--- Comparing residuals ---"
LOG_ARGS=()
for T in "${THREAD_COUNTS[@]}"; do
    LOG_ARGS+=("$VERIFY_DIR/race_T${T}.log")
done
python3 "$REPO_ROOT/scripts/compare_residuals.py" \
    "${LOG_ARGS[@]}" 2>&1 | tee "$VERIFY_DIR/race_compare.md"
rc=${PIPESTATUS[0]}

# Restore IC one last time so Run11 is left in its original state.
restore_ic

if [[ $rc -eq 0 ]]; then
    echo
    echo "=== Task A: PASS ==="
else
    echo
    echo "=== Task A: FAIL (rc=$rc) ===" >&2
fi
exit $rc
