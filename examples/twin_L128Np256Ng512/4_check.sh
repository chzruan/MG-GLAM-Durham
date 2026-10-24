#!/bin/bash
# Step 4: verify both runs with DEGRACE-companion/check_run.py: the schedule against Init.dat (exact redshifts,
# output classes), every P(k) and lens-plane file, the catalogues at the full outputs, the ic2pm velocity epoch,
# the timing and the Slurm accounting. Run it before step 5 deletes the native files (it reads both layouts).
set -euo pipefail
source "$(dirname "$0")/config.sh"
rc=0
for twin in $TWINS; do
    echo "=== $twin"
    $PYTHON "$REPO/DEGRACE-companion/check_run.py" "$WORK/$twin" --box "$REALISATION" --planes all || rc=1
done
exit $rc
