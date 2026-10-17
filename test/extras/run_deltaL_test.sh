#!/bin/bash
# I7 (route 1) regression test for the FML delta_L dump.
#   V7  determinism: the 2LPT smoke-test cube is byte-identical with 1, 2 and 4 MPI ranks
#   A/B layout: check_deltaL.py on 1LPT (exact: -i k.Psi of the lattice displacements) and 2LPT ICs
#   V3  end to end: FML IC -> ic2pm (half) -> PMP2main (128^3/256^3, 40 steps, Extra outputs) ->
#       deltam/delta_L = D_code(a) at k <= 0.15, r_cross > 0.99, no coordinate offset
# Light (login-node) job: ~4 min, 4 MPI ranks / 4 OpenMP threads, < 1 GB.  All artifacts in
# test/extras/work/deltaL (gitignored, safe to delete).  Exit 0 = all checks passed.
set -u
ROOT=$(cd "$(dirname "$0")/../.." && pwd)
EX=$ROOT/2LPTIC_Gui/FML/FML/LPT/example
W=$ROOT/test/extras/work/deltaL
PY="micromamba run -n cosemu python3 -I"
if [ -z "${SKIP_MODULES:-}" ]; then module purge; module load intel_comp/2024.2.0 compiler-rt tbb compiler mpi >/dev/null 2>&1; fi
export LD_LIBRARY_PATH=/cosma/local/fftw/intel_2024.2.0_intel_mpi_2024.2.0/3.3.10-epyc/lib:/cosma/local/gsl/2.7.1/lib:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4} MPLCONFIGDIR=/tmp/matplotlib-config
fail=0
rm -rf "$W"; mkdir -p "$W"
[ -x "$EX/test" ] || ( cd "$EX" && make test > "$W/build.log" 2>&1 ) || { echo "FATAL: generator build failed"; tail -5 "$W/build.log"; exit 2; }
for b in PMP2init.exe PMP2main.exe ic2pm.exe; do [ -x "$ROOT/$b" ] || { echo "FATAL: $ROOT/$b missing (make PMP2init PMP2main ic2pm)"; exit 2; }; done

gen () {   # $1 = lpt order, $2 = MPI ranks  -> $W/fml_<lpt>lpt_np<ranks>/snap/IC*
    local d="$W/fml_${1}lpt_np$2"; mkdir -p "$d/snap"
    sed -e "s/^lpt_order = 2/lpt_order = $1/" -e "s#\"pofk_bli_z49.txt\"#\"$EX/pofk_bli_z49.txt\"#" \
        -e 's#"snap/IC_test_Np128_L256"#"snap/IC"#' "$EX/input.lua" > "$d/input.lua"
    grep -q '^write_deltaL *= *true' "$d/input.lua" || { echo 'write_deltaL = true'; echo 'deltaL_kc = 0.5'; } >> "$d/input.lua"
    ( cd "$d" && OMP_NUM_THREADS=1 nice -n 19 mpirun -np "$2" "$EX/test" input.lua > run.log 2>&1 ) \
        || { echo "FAIL generator ${1}LPT with $2 ranks"; tail -3 "$d/run.log"; fail=1; }
}

echo "== V7 determinism (2LPT, 1/2/4 ranks)"
gen 2 1; gen 2 2; gen 2 4
m1=$(md5sum "$W/fml_2lpt_np1/snap/IC_deltaL.bin" | cut -c1-32); m2=$(md5sum "$W/fml_2lpt_np2/snap/IC_deltaL.bin" | cut -c1-32); m4=$(md5sum "$W/fml_2lpt_np4/snap/IC_deltaL.bin" | cut -c1-32)
if [ "$m1" = "$m2" ] && [ "$m2" = "$m4" ]; then echo "PASS identical cubes with 1, 2 and 4 ranks ($m1)"; else echo "FAIL md5 differ: $m1 $m2 $m4"; fail=1; fi

echo "== layout A/B, 1LPT (exact)"
gen 1 4
$PY "$ROOT/test/extras/check_deltaL.py" "$W/fml_1lpt_np4/snap/IC" "$W/fml_1lpt_np4/snap/IC" > "$W/check_1lpt.log" 2>&1 || fail=1
grep -E '^(FAIL|A |B |all)' "$W/check_1lpt.log"
echo "== layout A/B, 2LPT"
$PY "$ROOT/test/extras/check_deltaL.py" "$W/fml_2lpt_np4/snap/IC" "$W/fml_2lpt_np4/snap/IC" > "$W/check_2lpt.log" 2>&1 || fail=1
grep -E '^(FAIL|A |B |all)' "$W/check_2lpt.log"

echo "== V3 end to end (GLAM 128^3 / 256^3, L = 256, z 49 -> 30, 40 steps)"
G="$W/glam"; mkdir -p "$G/Run1"
awk '/^#outputs/ { printf "#outputs = %9d\n%s\n", -2, " 48.0 30.0"; getline; next } { print }' "$ROOT/test/precise_z/Init.base.dat" \
 | sed -e 's/^Box      =  150.000/Box      =  256.000/' -e 's/^Nrow     =        36/Nrow     =       128/' \
       -e 's/^Ngrid    =        72/Ngrid    =       256/' -e 's/^step da  =  8.0000E-04/step da  =  6.0000E-04/' \
       -e 's/^DM power spectrum         =         0/DM power spectrum         =         1/' > "$G/Init.dat"
{ echo 'Output classes            =         0  0-all class 1'
  echo 'Lens slabs per axis       =         0  0-off'
  echo 'Lens plane grid           =         0  0-Ngrid'
  echo 'Extra outputs             =         1  0-off'
  echo 'kc_delta                  =     0.500  h/Mpc'
  echo 'kc_mom                    =     0.300  h/Mpc'   # below k_N(N_dep = 128)/4 = 0.39
  echo 'N_dep_mom                 =         0  0-Ngrid/2'
  echo 'Subsample M               =         0  0-off'
  echo 'fR at halos               =         0  0-off'; } >> "$G/Init.dat"
# PkTable header -> Om = 0.3, OmL = 0.7, h = sqrt((Omb+Omc)/Om) = 0.7, matching input.lua; the body is not used (ICs come from ic2pm)
{ echo '#Om_b h^2             =  0.022000'; echo '#Om_c h^2             =  0.125000'; echo '#Om_Lambda            =  0.7000'
  echo '#Om_m (1-Om_K-Om_L)   =  0.3000'; tail -n +5 "$ROOT/test/PkTable.dat"; } > "$G/PkTable.dat"
cp "$ROOT/test/TableSeeds.dat" "$G/"
( cd "$G" && nice -n 19 "$ROOT/PMP2init.exe" > init.log 2>&1 ) || { echo "FAIL PMP2init"; tail -3 "$G/init.log"; fail=1; }
( cd "$G/Run1" && nice -n 19 "$ROOT/ic2pm.exe" "$W/fml_2lpt_np4/snap/IC." 1.0 half > ic2pm.log 2>&1 ) || { echo "FAIL ic2pm"; tail -3 "$G/Run1/ic2pm.log"; fail=1; }
grep -E 'velocity epoch|rms written' "$G/Run1/ic2pm.log"
( cd "$G/Run1" && echo 40 | nice -n 19 "$ROOT/PMP2main.exe" > main.log 2>&1 ) || { echo "FAIL PMP2main"; tail -3 "$G/Run1/main.log"; fail=1; }
grep -vE '^#' "$G/Run1/OutputSchedule.dat" 2>/dev/null
$PY "$ROOT/test/extras/check_glam_deltaL.py" "$W/fml_2lpt_np4/snap/IC" "$G/Run1" > "$W/check_glam.log" 2>&1 || fail=1
grep -E '^(FAIL|==|   0\.|   wrong|   no coord|PASS   no|all)' "$W/check_glam.log" | grep -vE '^PASS   [0-9]'
rm -f "$G/Run1"/PMcrs*.DAT "$W"/fml_*/snap/IC.[0-9]*      # particle files (16 MB each), not needed after the checks
if [ $fail -eq 0 ]; then echo "ALL delta_L TESTS PASSED"; else echo "SOME delta_L TESTS FAILED (logs in $W)"; fi
exit $fail
