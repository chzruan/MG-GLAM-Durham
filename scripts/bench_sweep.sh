#!/bin/bash
# Task C driver: benchmark sweep for MG relaxation time vs OMP thread count.
#
# Uses the production PMP2MG.exe (must exist — built via `make PMP2MG`)
# so we measure what the user actually runs for MG (fid_F5). PMP2main.exe
# is LCDM-only and skips the MG code path. At each OMP_NUM_THREADS ∈
# {16,32,64,128} runs 3 repeats of Nsteps=5, preserving each run's
# timing.log as verify/bench_T{N}_r{R}.log. Before each run the IC state
# is restored from the same Run11 snapshot verify_race.sh uses.
#
# Run from repo root on a Cosma8 compute node with Intel toolchain loaded.
#
#   bash scripts/bench_sweep.sh
#
# Final step: bench_report.py is invoked to emit a markdown summary.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"

VERIFY_DIR="$REPO_ROOT/verify"
RUN_DIR="$REPO_ROOT/fid_F5/Run11"
IC_SNAPSHOT="$RUN_DIR/.race_ic_snapshot"
NSTEPS=5
REPEATS=3
THREAD_COUNTS=(16 32 64 128)

IC_FILES=(PMcrd.DAT pt.dat fort.17)
# PMcrs*.DAT are the big particle-data files; globbed at snapshot time.

mkdir -p "$VERIFY_DIR"

echo "=== Task C: MG benchmark sweep ==="
echo "Binary:   $REPO_ROOT/PMP2MG.exe"
echo "Run dir:  $RUN_DIR"
echo "Nsteps:   $NSTEPS   Repeats: $REPEATS"
echo "Threads:  ${THREAD_COUNTS[*]}"
echo

if [[ ! -x "$REPO_ROOT/PMP2MG.exe" ]]; then
    echo "ERROR: $REPO_ROOT/PMP2MG.exe missing — run 'make PMP2MG' first" >&2
    exit 2
fi

# Snapshot Run11 IC if the race driver hasn't already. Invalidates any
# snapshot that predates the current PMcrs0.DAT (IC was regenerated).
snapshot_stale=0
if [[ -d "$IC_SNAPSHOT" && -f "$RUN_DIR/PMcrs0.DAT" ]]; then
    if [[ ! -f "$IC_SNAPSHOT/PMcrs0.DAT" \
          || "$RUN_DIR/PMcrs0.DAT" -nt "$IC_SNAPSHOT/PMcrs0.DAT" ]]; then
        snapshot_stale=1
    fi
fi
if [[ ! -d "$IC_SNAPSHOT" || "$snapshot_stale" -eq 1 ]]; then
    echo "--- Snapshotting IC state from $RUN_DIR ---"
    rm -rf "$IC_SNAPSHOT"
    mkdir -p "$IC_SNAPSHOT"
    for f in "${IC_FILES[@]}"; do
        if [[ ! -f "$RUN_DIR/$f" ]]; then
            echo "ERROR: missing IC file $RUN_DIR/$f" >&2
            exit 2
        fi
        cp -p "$RUN_DIR/$f" "$IC_SNAPSHOT/$f"
    done
    # Snapshot every PMcrs*.DAT present (at least PMcrs0.DAT; some
    # configurations also write PMcrs1.DAT etc.).
    shopt -s nullglob
    pmcrs_files=("$RUN_DIR"/PMcrs*.DAT)
    shopt -u nullglob
    if [[ ${#pmcrs_files[@]} -eq 0 ]]; then
        echo "ERROR: no PMcrs*.DAT in $RUN_DIR — run IC generator first" >&2
        exit 2
    fi
    for src in "${pmcrs_files[@]}"; do
        cp -p "$src" "$IC_SNAPSHOT/$(basename "$src")"
    done
    echo "    snapshotted: ${IC_FILES[*]} $(basename -a "${pmcrs_files[@]}" | tr '\n' ' ')"
    echo
fi

restore_ic() {
    for f in "${IC_FILES[@]}"; do
        cp -p "$IC_SNAPSHOT/$f" "$RUN_DIR/$f"
    done
    shopt -s nullglob
    local src
    for src in "$IC_SNAPSHOT"/PMcrs*.DAT; do
        cp -p "$src" "$RUN_DIR/$(basename "$src")"
    done
    shopt -u nullglob
}

export OMP_PROC_BIND=close
export OMP_PLACES=cores
export OMP_STACKSIZE=512M

for T in "${THREAD_COUNTS[@]}"; do
    for R in $(seq 1 "$REPEATS"); do
        LOG_COPY="$VERIFY_DIR/bench_T${T}_r${R}.log"
        STDOUT_COPY="$VERIFY_DIR/bench_T${T}_r${R}.stdout"
        echo "--- T=$T  r=$R ---"
        restore_ic
        pushd "$RUN_DIR" >/dev/null
        # Delete any prior timing.log so a stale copy can't masquerade as
        # the output of a run that silently failed at IC load.
        rm -f timing.log
        OMP_NUM_THREADS="$T" \
            "$REPO_ROOT/PMP2MG.exe" <<<"$NSTEPS" \
            >"$STDOUT_COPY" 2>&1
        if [[ ! -f timing.log ]]; then
            echo "ERROR: no timing.log written for T=$T r=$R" >&2
            echo "       tail of stdout:" >&2
            tail -5 "$STDOUT_COPY" >&2
            popd >/dev/null
            exit 3
        fi
        # PMP2MG exits 0 even on IC-read failure; detect that case here.
        if grep -qE "Error reading|Did not get all files" "$STDOUT_COPY"; then
            echo "ERROR: PMP2MG aborted at startup (IC read failure) for T=$T r=$R" >&2
            tail -5 "$STDOUT_COPY" >&2
            popd >/dev/null
            exit 3
        fi
        cp -p timing.log "$LOG_COPY"
        popd >/dev/null
        nsteps=$(awk 'NR>2 && $1 ~ /^[0-9]+$/ {n++} END {print n+0}' "$LOG_COPY")
        echo "    steps logged: $nsteps -> $LOG_COPY"
    done
done
echo

# Restore IC so Run11 is left in pristine state.
restore_ic

# ---- Summarize -------------------------------------------------------------
echo "--- Writing summary report ---"
LOG_ARGS=()
for T in "${THREAD_COUNTS[@]}"; do
    for R in $(seq 1 "$REPEATS"); do
        LOG_ARGS+=("$VERIFY_DIR/bench_T${T}_r${R}.log")
    done
done
python3 "$REPO_ROOT/scripts/bench_report.py" \
    --baseline 16 \
    "${LOG_ARGS[@]}" | tee "$VERIFY_DIR/bench_report.md"

echo
echo "=== Task C: done ==="
echo "Summary: $VERIFY_DIR/bench_report.md"
