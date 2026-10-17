#!/bin/bash
# ic2pm ingestion test (review R3-M1): particles are placed at array index = Gadget ID, the ID
# lattice order is detected and tagged in the PM header (IDz: 2LPTic-like int32 sequential IDs,
# IDx: FML-like int64 IDs in shuffled file order), and a non-permutation ID set is refused.
# Verified through GLAM itself: one PM step with Subsample M = 1 lists the particles in index
# order, so particle k must sit at the lattice site of ID k+1. Seconds; 16^3 particles, 32^3 mesh.
set -u
ROOT=$(cd "$(dirname "$0")/../.." && pwd); W=$ROOT/test/extras/work/ic2pm
PY="micromamba run -n cosemu python3 -I"
if [ -z "${SKIP_MODULES:-}" ]; then module purge; module load intel_comp/2024.2.0 compiler-rt tbb compiler mpi >/dev/null 2>&1; fi
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
for b in PMP2init.exe PMP2main.exe ic2pm.exe; do [ -x "$ROOT/$b" ] || { echo "FATAL: $ROOT/$b missing"; exit 2; }; done
rm -rf "$W"; mkdir -p "$W"; fail=0
for mode in z32seq:IDz x64shuf:IDx dup:FATAL; do
    m=${mode%%:*}; want=${mode##*:}; d=$W/$m; mkdir -p "$d/Run1"
    echo "== $m (expect $want)"
    $PY "$ROOT/test/extras/make_gadget_lattice.py" "$d/ic" "$m" > "$d/make.log" 2>&1 || { echo "FAIL generator"; fail=1; continue; }
    awk '/^#outputs/ { printf "#outputs = %9d\n%s\n", -1, " 48.0"; getline; next } { print }' "$ROOT/test/precise_z/Init.base.dat" \
     | sed -e 's/^Box      =  150.000/Box      =   16.000/' -e 's/^Nrow     =        36/Nrow     =        16/' \
           -e 's/^Ngrid    =        72/Ngrid    =        32/' > "$d/Init.dat"
    { echo 'Output classes            =         0  0-all class 1'; echo 'Lens slabs per axis       =         0  0-off'
      echo 'Lens plane grid           =         0  0-Ngrid'; echo 'Extra outputs             =         1  0-off'
      echo 'kc_delta                  =     1.500  h/Mpc'; echo 'kc_mom                    =     0.750  h/Mpc'
      echo 'N_dep_mom                 =         0  0-Ngrid/2'; echo 'Subsample M               =         1  every particle'
      echo 'fR at halos               =         0  0-off'; } >> "$d/Init.dat"
    cp "$ROOT/test/PkTable.dat" "$ROOT/test/TableSeeds.dat" "$d/"
    ( cd "$d" && nice -n 19 "$ROOT/PMP2init.exe" > init.log 2>&1 ) || { echo "FAIL PMP2init"; tail -3 "$d/init.log"; fail=1; continue; }
    ( cd "$d/Run1" && nice -n 19 "$ROOT/ic2pm.exe" "$d/ic." 1.0 half > ic2pm.log 2>&1 ); rc=$?
    if [ "$want" = FATAL ]; then
        if [ $rc -ne 0 ] && grep -q 'FATAL' "$d/Run1/ic2pm.log"; then echo "PASS ic2pm refused the duplicate-ID set: $(grep FATAL "$d/Run1/ic2pm.log" | head -1)"
        else echo "FAIL ic2pm accepted a non-permutation ID set (rc=$rc)"; fail=1; fi
        continue
    fi
    [ $rc -eq 0 ] || { echo "FAIL ic2pm rc=$rc"; tail -3 "$d/Run1/ic2pm.log"; fail=1; continue; }
    grep -E 'ID lattice|placed at' "$d/Run1/ic2pm.log"
    ( cd "$d/Run1" && echo 1 | nice -n 19 "$ROOT/PMP2main.exe" > main.log 2>&1 ) || { echo "FAIL PMP2main"; tail -3 "$d/Run1/main.log"; fail=1; continue; }
    $PY "$ROOT/test/extras/check_ic2pm_order.py" "$d/Run1" "$want" || fail=1
done
rm -f "$W"/*/Run1/PMcrs*.DAT "$W"/*/ic.[0-9]
if [ $fail -eq 0 ]; then echo "ALL ic2pm TESTS PASSED"; else echo "SOME ic2pm TESTS FAILED (logs in $W)"; fi
exit $fail
