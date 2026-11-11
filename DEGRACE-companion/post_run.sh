#!/bin/bash
# Post-processing of a finished GLAM box in one cosma8-serial job: check_run.py, then pack_run.sh (pack, byte check,
# prune). Packing starts only when check_run.py reports 0 FAIL and 0 WARN; otherwise the originals stay untouched and
# the job exits with an error (a WARN, such as a skipped Rockstar catalogue, needs a decision first: AGENTS.md).
#
# Chain it to the simulation job, so that it starts when the run has finished successfully (one core instead of the
# simulation's whole node):
#   j=$(sbatch --parsable submit.sh)
#   sbatch --dependency=afterok:$j --kill-on-invalid-dep=yes -J post_<name> \
#          /path/to/MG-GLAM/DEGRACE-companion/post_run.sh "$PWD/Run<box>" $j
# The second argument (the simulation's job ID) lets check_run.py report its sacct record. PRUNE=no is passed on to
# pack_run.sh. A rerun is safe: pack_run.sh continues from the state of the box.
#SBATCH -p cosma8-serial
#SBATCH -A dp203
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH -t 12:00:00
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
set -euo pipefail
export MAMBA_USE_LOCKFILES=false
RUN=$(readlink -f "${1:?usage: post_run.sh /abs/path/Run<box> [simulation job ID]}")
SIMJOB=${2:-}
PY=${PYTHON:-micromamba run -n cosemu python3}
fail() { echo "error: $*" >&2; exit 1; }
[[ $(basename "$RUN") =~ ^Run([0-9]+)$ ]] || fail "$RUN is not a Run<box> directory"
BOX=$((10#${BASH_REMATCH[1]}))
RUNDIR=$(dirname "$RUN")
# the repository: MGGLAM_REPO, or the directory of this script (under sbatch, the path it was submitted with)
if [[ -z ${MGGLAM_REPO:-} ]]; then
    self=$0
    if [[ -n ${SLURM_JOB_ID:-} ]]; then
        self=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/^ *Command=/{print $2; exit}' | awk '{print $1}')
    fi
    MGGLAM_REPO=$(cd "$(dirname "$self")/.." && pwd)
fi
[[ -f $MGGLAM_REPO/DEGRACE-companion/check_run.py && -f $MGGLAM_REPO/DEGRACE-companion/pack_run.sh ]] ||
    fail "MG-GLAM repository not found at $MGGLAM_REPO (set MGGLAM_REPO)"
echo "== $(date -Is) $(hostname) post-processing $RUN (box $BOX) with $MGGLAM_REPO ($(git -C "$MGGLAM_REPO" rev-parse --short HEAD 2>/dev/null || echo 'no git'))"

REPORT=$RUNDIR/check_run_${SIMJOB:-${SLURM_JOB_ID:-$$}}.txt
echo "== $(date -Is) check_run.py -> $REPORT"
set +e
$PY "$MGGLAM_REPO/DEGRACE-companion/check_run.py" "$RUNDIR" --box "$BOX" --planes all ${SIMJOB:+--job "$SIMJOB"} > "$REPORT" 2>&1
rc=$?
set -e
tail -1 "$REPORT"
[[ $rc -eq 0 ]] || fail "check_run.py failed (exit $rc); nothing packed, see $REPORT"
tail -1 "$REPORT" | grep -q '^0 FAIL, 0 WARN$' || fail "check_run.py reports warnings; nothing packed, see $REPORT"

echo "== $(date -Is) pack_run.sh"
bash "$MGGLAM_REPO/DEGRACE-companion/pack_run.sh" "$RUN"
echo "== $(date -Is) post-processing done"
