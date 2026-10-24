#!/bin/bash
# Post-processing of a finished GLAM box: pack Run<box>/ into one HDF5 file per output step (glamio.pack),
# prove that the packed files hold every original byte for byte, then delete the originals.
#
# usage (from the run directory, after the simulation job has finished, check_run.py has passed and every script that
# needs the native files has run):
#   sbatch -J pack_<name> /path/to/MG-GLAM/DEGRACE-companion/pack_run.sh "$PWD/Run<box>"                        # pack, check, prune
#   sbatch -J pack_<name> --export=ALL,PRUNE=no /path/to/MG-GLAM/DEGRACE-companion/pack_run.sh "$PWD/Run<box>"  # keep the originals
# The originals are deleted only after the byte check below has passed; unpack rebuilds them byte for byte.
#
# The job reads `python -m glamio.pack status` and continues from where the box is:
#   native / incomplete / packed:
#     1. pack in place: Output.<step>.<box>.h5 per output step, RunInfo.<box>.h5 and packed.md5 beside the originals
#        (each file is verified against the md5 of its originals, flushed, then renamed; an existing packed file is
#        kept only if it still matches the current originals);
#     2. byte check: unpack EVERY packed file into $CHECK_DIR (default pack-check-<job> beside Run<box>/, i.e. in the run directory),
#        require the same file list and `cmp` equality of every file with Run<box>, then remove the rebuilt copy
#        (needs free space equal to the box: 75 GB for L1024 Np2048; removed on failure too);
#     3. unless PRUNE=no: `glamio.pack prune --yes` (refuses while any original is in no packed file; re-verifies every
#        packed file and the md5 of every original; deletes nothing if anything differs), then `md5sum -c packed.md5`
#        and a read of the packed run through glamio;
#   partially-pruned (an interrupted prune): unless PRUNE=no, byte-checks the remaining originals (step 2 on the files still
#     there) and finishes the prune (step 3); otherwise the job stops;
#   pruned: verify only; busy (a .partial exists) or inconsistent: the job stops and changes nothing.
# Any failure stops the job before an original is touched. Cost for an L1024 Np2048 box (75 GB, one core): pack 20 min,
# check 15-30 min, prune 30 min, measured on 2026-10-24 with the equivalent scratch jobs 12217626, 12217719, 12217742
# (peak RSS 0.8 GB; sacct MaxRSS shows the 8 GB limit, which is page cache).
# Restore the native files later (75 GB: under Slurm) with: python -m glamio.pack unpack Run<box> --out DIR
#SBATCH -p cosma8-serial
#SBATCH -A dp203
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH -t 03:00:00
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
set -euo pipefail
export MAMBA_USE_LOCKFILES=false
RUN=$(readlink -f "${1:?usage: pack_run.sh /abs/path/Run<box>}")
PRUNE=${PRUNE:-yes}
JOB=${SLURM_JOB_ID:-$$}
CHECK_DIR=${CHECK_DIR:-$(dirname "$RUN")/pack-check-$JOB}
PY=${PYTHON:-micromamba run -n cosemu python3}     # needs numpy and h5py; glamio from PYTHONPATH or the environment
CREATED_CHECK_DIR=""
cleanup() { if [[ -n $CREATED_CHECK_DIR && -d $CREATED_CHECK_DIR ]]; then rm -r "${CREATED_CHECK_DIR:?}"; echo "removed $CREATED_CHECK_DIR"; fi; }
trap cleanup EXIT
fail() { echo "error: $*" >&2; exit 1; }

[[ $(basename "$RUN") =~ ^Run[0-9]+$ ]] || fail "$RUN is not a Run<box> directory"
case "$(readlink -m "$CHECK_DIR")/" in "$RUN"/*) fail "CHECK_DIR $CHECK_DIR lies inside $RUN" ;; esac
code_sha() { $PY -c 'from glamio.pack import _code_sha256; print(_code_sha256())'; }   # = glamio_code_sha256 in the files
CODE_SHA=$(code_sha)
echo "== $(date -Is) $(hostname) $RUN, PRUNE=$PRUNE"
echo "glamio $($PY -c 'import glamio; print(glamio.__version__, glamio.__file__)'), sources sha256 $CODE_SHA"
STATUS_OUT=$($PY -m glamio.pack status "$RUN") || true
echo "$STATUS_OUT"
STATE=$(awk '/^status:/{print $2}' <<< "$STATUS_OUT")
[[ -n $STATE ]] || fail "glamio.pack status failed"

do_pack() {
    echo "== $(date -Is) 1. pack: $(find "$RUN" -type f | wc -l) files, $(du -s --apparent-size -B1 "$RUN" | cut -f1) bytes"
    $PY -u -m glamio.pack pack "$RUN"
}
do_check() {    # $1 = full: the rebuilt file list must equal the originals'; remaining: every original still present is checked
    echo "== $(date -Is) 2. byte check ($1): unpack everything into $CHECK_DIR and compare"
    [[ ! -e $CHECK_DIR ]] || fail "$CHECK_DIR exists"
    mkdir -p "$CHECK_DIR"
    CREATED_CHECK_DIR=$CHECK_DIR
    $PY -m glamio.pack unpack "$RUN" --out "$CHECK_DIR/rebuilt" | tail -1
    (cd "$RUN" && find . -type f ! -name 'Output.*.h5' ! -name 'RunInfo.*.h5' ! -name packed.md5 | sort) > "$CHECK_DIR/files_original.txt"
    (cd "$CHECK_DIR/rebuilt" && find . -type f | sort) > "$CHECK_DIR/files_rebuilt.txt"
    if [[ $1 == full ]]; then
        cmp "$CHECK_DIR/files_original.txt" "$CHECK_DIR/files_rebuilt.txt" || fail "file lists differ"
    else
        [[ -z $(comm -23 "$CHECK_DIR/files_original.txt" "$CHECK_DIR/files_rebuilt.txt") ]] || fail "originals in no packed file"
    fi
    while IFS= read -r f; do
        cmp -- "$RUN/$f" "$CHECK_DIR/rebuilt/$f" || fail "$f differs"
    done < "$CHECK_DIR/files_original.txt"
    echo "byte check passed ($1): $(wc -l < "$CHECK_DIR/files_original.txt") original files, every one cmp-identical to its rebuilt copy"
    rm -r "${CHECK_DIR:?}"
    CREATED_CHECK_DIR=""
}
do_prune() {
    [[ $(code_sha) == "$CODE_SHA" ]] || fail "the glamio sources changed during the job; nothing deleted"
    echo "== $(date -Is) 3. prune"
    $PY -m glamio.pack prune "$RUN" --yes
    echo "after: $(find "$RUN" -type f | wc -l) files, $(find "$RUN" -mindepth 1 -type d | wc -l) subdirectories"
    (cd "$RUN" && md5sum -c --quiet packed.md5 && echo "packed.md5 OK: $(wc -l < packed.md5) files")
    $PY - "$RUN" <<'EOF'
import sys
import glamio
run = sys.argv[1]
s = glamio.read_schedule(glamio.product_path(run, "schedule", 0))
box = int(glamio.split_packed(s.path)[0].rsplit(".", 2)[-2])
n = 0
for step in s.steps:
    found = glamio.products_at(run, int(step), box)
    assert found and all(".h5/" in p for p in found.values()), (int(step), found)
    n += len(found)
msg = f"glamio reads the packed run: {len(s)} outputs, {n} products"
last = glamio.products_at(run, int(s.steps[-1]), box)
if "catalogue" in last:
    msg += f", {len(glamio.read_catalogue(last['catalogue'], mmap=True))} haloes at step {int(s.steps[-1])}"
print(msg)
EOF
}

case $STATE in
    native|incomplete|packed)
        do_pack
        do_check full
        if [[ $PRUNE != no ]]; then do_prune; else echo "originals kept (PRUNE=no); resubmit without it to delete them"; fi ;;
    partially-pruned)
        [[ $PRUNE != no ]] || fail "$RUN is partially pruned (an interrupted prune?); resubmit without PRUNE=no to finish it"
        do_check remaining
        do_prune ;;
    pruned)
        echo "== already packed and pruned: verify only"
        $PY -m glamio.pack verify "$RUN"
        (cd "$RUN" && md5sum -c --quiet packed.md5 && echo "packed.md5 OK") ;;
    *)
        fail "state '$STATE': nothing changed; inspect $RUN (python -m glamio.pack status)" ;;
esac
echo "== $(date -Is) done"
