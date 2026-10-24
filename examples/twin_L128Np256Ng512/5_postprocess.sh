#!/bin/bash
# Step 5: post-process both runs with DEGRACE-companion/pack_run.sh (one small Slurm job each): pack Run<N>/ into one
# HDF5 file per output step, rebuild every original file and compare it byte for byte, then delete the originals.
# Afterwards Run<N>/ holds Output.<step>.<tag>.h5 per output, RunInfo.<tag>.h5 and packed.md5; glamio reads them as
# it read the native files. KEEP_ORIGINALS=1 packs and checks without deleting. pack_run.sh reads the state of each
# run, so a rerun resumes (or only re-verifies a finished one).
set -euo pipefail
source "$(dirname "$0")/config.sh"
trap stop_waiting INT TERM
prune=yes
[ "${KEEP_ORIGINALS:-0}" = 1 ] && prune=no
pids=()
for twin in $TWINS; do
    run=$WORK/$twin
    submit_wait "ex_pack_$twin" "$run" -p "$PARTITION" -A "$ACCOUNT" -o "$run/pack_%j.out" -e "$run/pack_%j.err" \
                --export=ALL,PRUNE=$prune "$REPO/DEGRACE-companion/pack_run.sh" "$run/Run$REALISATION" &
    pids+=($!)
done
rc=0
for p in "${pids[@]}"; do wait "$p" || rc=1; done
for twin in $TWINS; do
    echo "=== $twin: $(ls "$WORK/$twin/Run$REALISATION" | wc -l) files in Run$REALISATION/"
    log=$(ls -t "$WORK/$twin"/pack_*.out 2> /dev/null | head -1)
    [ -n "$log" ] && grep -E "^status|byte check passed|deleted|originals kept|verify only|error" "$log" || true
done
exit $rc
