#!/bin/bash
# The whole example in one go: ICs, inputs, both runs, checks, post-processing and the figure (~30 min, mostly
# waiting for Slurm). Each step skips work that is already done and waits for jobs still queued or running, so
# the script can be rerun after a failure or a Ctrl-C.
set -euo pipefail
here=$(dirname "$0")
source "$here/config.sh"
for step in 1_ics 2_inputs 3_run 4_check 5_postprocess; do
    echo "################ $step $(date +%H:%M:%S)"
    bash "$here/$step.sh"
done
echo "################ analyse $(date +%H:%M:%S)"
$PYTHON "$here/analyse.py"
