#!/bin/bash
# Step 3: run both twins at once (one Slurm job each, 16 threads, ~10-20 min) and wait for them.
# Each run writes into $WORK/<twin>/Run<N>/: P(k) at every output, lens planes, BDM halo catalogues, the low-k density
# and momentum cubes, the particle subsample (and f_R at haloes for f(R)), the schedule and the logs; job_done.txt
# marks a run that finished. A rerun waits for jobs still queued or running and skips finished runs.
set -euo pipefail
source "$(dirname "$0")/config.sh"
trap stop_waiting INT TERM
IC=$WORK/ICs/IC_L${BOX}Np${NP}_seed${SEED}.
[ -f "$WORK/ICs/md5sum.txt" ] || { echo "no ICs in $WORK/ICs: run 1_ics.sh first"; exit 1; }
pids=()
rc=0
for twin in $TWINS; do
    run=$WORK/$twin
    out=$run/Run$REALISATION
    [ -f "$run/Init.dat" ] || { echo "no inputs in $run: run 2_inputs.sh first"; exit 1; }
    if [ -z "$(jobs_in "ex_$twin" "$run")" ]; then
        if [ -f "$out/job_done.txt" ] || compgen -G "$out/RunInfo.*.h5" > /dev/null; then
            echo "$twin: finished earlier: skipped"; continue
        fi
        if [ -e "$out/main.log" ]; then
            echo "$twin: an earlier run did not finish (see $run/slurm_*.out); remove $run, then rerun 2_inputs.sh and this step"
            rc=1; continue
        fi
    fi
    submit_wait "ex_$twin" "$run" -p "$PARTITION" -A "$ACCOUNT" -o "$run/slurm_%j.out" -e "$run/slurm_%j.err" \
                --export=ALL,IC="$IC",REALISATION="$REALISATION" "$EXAMPLE/run_job.sh" &
    pids+=($!)
done
for p in "${pids[@]}"; do wait "$p" || rc=1; done
[ $rc = 0 ] && echo "runs finished; job logs: $WORK/<twin>/slurm_<job>.out" || echo "a run failed or is incomplete: see the messages above"
exit $rc
