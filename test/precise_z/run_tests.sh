#!/bin/bash
# Automated tests for the exact-redshift snapshot feature (#outputs < 0).
#
# Runs a tiny LCDM box (36^3 particles, 72^3 mesh, 150 Mpc/h) through
# PMP2init -> PMP2start -> PMP2main for several configurations, then checks
# every produced PMcrd.NNNN.DAT header against an independent Real*4
# replication of the schedule arithmetic (asserts.py).
#
# Config convention under test: the existing '#outputs' line + redshift list.
#   #outputs =  N   -> legacy: analyze at the step CLOSEST to each redshift
#   #outputs = -N   -> exact:  hit each listed redshift exactly (new steps are
#                      inserted, or a step within 0.5% in 1+z is moved onto it)
#
# Cases (see asserts.py for the exact expectations):
#   case5a_legacy    - positive #outputs: legacy schedule (baseline)
#   case5b_trailing  - Setup.dat with stale trailing lines appended (e.g. the
#                      interim Nexact block): must behave identically to 5a
#   case1_merge      - exact z equal to an existing output moment: merged,
#                      single snapshot, no duplicate step
#   case2_insert     - exact z between two scheduled steps: new step inserted,
#                      snapshot lands bit-exactly on 1/(1+z)
#   case3_multi      - three exact z given unsorted: all produced, in order
#   case4a/4b        - invalid z (>= z_init / negative) rejected by PMP2init
#   case4c/4d        - close pair / out-of-range z in a hand-edited Setup.dat
#                      rejected by PMP2main
#   case6_tailinsert - z=0 inserted after the legacy stop step (NlastX fix)
#
# Usage:   bash run_tests.sh          (interactive; ~1-2 min on a login node)
#          SKIP_MODULES=1 bash run_tests.sh   (if the Intel env is already loaded)
set -u

HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../.." && pwd)
WORK=$HERE/work
NSTEPS=100
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}

# ---------------------------------------------------------------- environment
if [ -z "${SKIP_MODULES:-}" ] && command -v module >/dev/null 2>&1; then
    module purge 2>/dev/null
    module load intel_comp/2024.2.0 2>/dev/null
    module load tbb compiler-rt 2>/dev/null
    module load compiler mpi 2>/dev/null
fi
for exe in PMP2init.exe PMP2start.exe PMP2main.exe; do
    if [ ! -x "$ROOT/$exe" ]; then
        echo "== building $exe"
        ( cd "$ROOT" && make "${exe%.exe}" ) || { echo "FATAL: cannot build $exe"; exit 2; }
    fi
done

rm -rf "$WORK"; mkdir -p "$WORK"

# ------------------------------------------------------------ input data files
# P(k) table: any consistent table works (cosmology is read from its header).
for src in "$ROOT/test/PkTable.dat" "$ROOT/fid2LPTIC_L512Np2048Ng4096/PkTable.dat"; do
    [ -f "$src" ] && { cp "$src" "$WORK/PkTable.dat"; break; }
done
[ -f "$WORK/PkTable.dat" ] || { echo "FATAL: no PkTable.dat found"; exit 2; }
# Seed table: reuse an existing one, else write a minimal valid table.
for src in "$ROOT/test/TableSeeds.dat" "$ROOT/LCDM_L1024Np2048Ng4096/TableSeeds.dat"; do
    [ -f "$src" ] && { cp "$src" "$WORK/TableSeeds.dat"; break; }
done
if [ ! -f "$WORK/TableSeeds.dat" ]; then
    { echo ' Seeds:     1298302         137'
      for i in $(seq 1 20); do echo "     $((1298302 + 7919*i))           $i"; done
    } > "$WORK/TableSeeds.dat"
fi

# ------------------------------------------------------------------ helpers
prep_case () {   # $1 = case dir
    local d="$WORK/$1"
    mkdir -p "$d/Run1"
    cp "$HERE/Init.base.dat" "$d/Init.dat"
    cp "$WORK/PkTable.dat" "$WORK/TableSeeds.dat" "$d/"
}
set_outputs () { # $1 = case dir; $2 = signed count; $3 = redshift list
    awk -v n="$2" -v list="$3" '
        /^#outputs/ { printf "#outputs = %9d\n%s\n", n, list; getline; next }
        { print }' "$WORK/$1/Init.dat" > "$WORK/$1/Init.tmp" \
        && mv "$WORK/$1/Init.tmp" "$WORK/$1/Init.dat"
}
set_setup_outputs () { # $1 = Setup.dat path; $2 = signed count; $3 = list
    awk -v n="$2" -v list="$3" '
        /Number of redshifts for analysis/ {
            printf "%5d               Number of redshifts for analysis (<0: exact)\n%s\n", n, list
            getline; next }
        { print }' "$1" > "$1.tmp" && mv "$1.tmp" "$1"
}
copy_ics () {    # $1 = case dir
    cp "$WORK/ic/Run1/PMcrd.DAT" "$WORK/ic/Run1/PMcrs0.DAT" "$WORK/ic/Run1/pt.dat" "$WORK/$1/Run1/"
}
run_init () {    # $1 = case dir
    ( cd "$WORK/$1" && "$ROOT/PMP2init.exe" > init.log 2>&1 )
}
run_main () {    # $1 = case dir;  $2 = steps (default NSTEPS)
    ( cd "$WORK/$1/Run1" && echo "${2:-$NSTEPS}" | "$ROOT/PMP2main.exe" > main.log 2>&1 )
}

# --------------------------------------------------------------- shared ICs
# Simulation cases with the default da share one IC realization.
echo "== generating shared initial conditions"
mkdir -p "$WORK/ic/Run1"
cp "$HERE/Init.base.dat" "$WORK/ic/Init.dat"
cp "$WORK/PkTable.dat" "$WORK/TableSeeds.dat" "$WORK/ic/"
( cd "$WORK/ic" && "$ROOT/PMP2init.exe" > init.log 2>&1 ) || { echo "FATAL: PMP2init failed for ICs"; exit 2; }
( cd "$WORK/ic/Run1" && echo 1 | "$ROOT/PMP2start.exe" > start.log 2>&1 ) \
    || { echo "FATAL: PMP2start failed"; tail -5 "$WORK/ic/Run1/start.log"; exit 2; }

# ------------------------------------------------------------------ sim cases
echo "== case5a_legacy: positive #outputs (legacy nearest-step schedule)"
prep_case case5a_legacy
run_init case5a_legacy || echo "  (init failed)"
copy_ics case5a_legacy && run_main case5a_legacy

echo "== case5b_trailing: Setup.dat with stale trailing lines"
prep_case case5b_trailing
{ cat "$WORK/case5a_legacy/Setup.dat"
  echo ' !------------ exact output redshifts ---------------'
  echo '    0              Number of exact output redshifts'
} > "$WORK/case5b_trailing/Setup.dat"
copy_ics case5b_trailing && run_main case5b_trailing

echo "== case1_merge: exact z on an existing output moment (#outputs = -2)"
prep_case case1_merge
set_outputs case1_merge -2 " 2.00 0.49337"
run_init case1_merge || echo "  (init failed)"
copy_ics case1_merge && run_main case1_merge

echo "== case2_insert: exact z between two scheduled steps (#outputs = -2)"
prep_case case2_insert
set_outputs case2_insert -2 " 2.00 0.512"
run_init case2_insert || echo "  (init failed)"
copy_ics case2_insert && run_main case2_insert

echo "== case3_multi: three exact z, given unsorted (#outputs = -3)"
prep_case case3_multi
set_outputs case3_multi -3 " 2.5 0.512 1.1"
run_init case3_multi || echo "  (init failed)"
copy_ics case3_multi && run_main case3_multi

echo "== case6_tailinsert: z=0 inserted after the legacy stop step (da=6.9e-4)"
# With da0=6.9e-4 the last sub-1 ladder point (a=0.99448) satisfies the legacy
# half-step exit rule, while the z=0 target must be INSERTED after it; the run
# must still execute the inserted step (regression test for the NlastX fix).
prep_case case6_tailinsert
sed -i 's/8.0000E-04/6.9000E-04/' "$WORK/case6_tailinsert/Init.dat"
set_outputs case6_tailinsert -1 " 0.00"
run_init case6_tailinsert || echo "  (init failed)"
( cd "$WORK/case6_tailinsert/Run1" && echo 1 | "$ROOT/PMP2start.exe" > start.log 2>&1 )  # own ICs: da differs
run_main case6_tailinsert 300

# ---------------------------------------------------------------- error cases
echo "== case4a_zbig / case4b_zneg: rejected by PMP2init"
prep_case case4a_zbig
set_outputs case4a_zbig -1 " 60.0"
run_init case4a_zbig
prep_case case4b_zneg
set_outputs case4b_zneg -1 " -0.5"
run_init case4b_zneg

echo "== case4c_pair / case4d_range: rejected by PMP2main (hand-edited Setup.dat)"
mkdir -p "$WORK/case4c_pair/Run1" "$WORK/case4d_range/Run1"
cp "$WORK/case5a_legacy/Setup.dat" "$WORK/case4c_pair/Setup.dat"
set_setup_outputs "$WORK/case4c_pair/Setup.dat" -2 "   0.51200   0.51000"
cp "$WORK/case5a_legacy/Setup.dat" "$WORK/case4d_range/Setup.dat"
set_setup_outputs "$WORK/case4d_range/Setup.dat" -1 "  60.00000"
copy_ics case4c_pair;  run_main case4c_pair
copy_ics case4d_range; run_main case4d_range

# ------------------------------------------------------------------ verdicts
echo
python3 "$HERE/asserts.py" "$WORK"
rc=$?
exit $rc
