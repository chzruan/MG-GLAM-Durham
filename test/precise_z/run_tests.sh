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
#                      inserted, or an unpinned step within 0.5% in 1+z is
#                      moved onto it; targets may be arbitrarily close)
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
#   case4c/4d        - duplicate z / out-of-range z in a hand-edited Setup.dat
#                      rejected by PMP2main
#   case6_tailinsert - z=0 inserted after the legacy stop step (NlastX fix)
#   case7_closepairs - the four DEGRACE pairs closer than 0.5%: all scheduled
#   case8_pinned     - second member of a pair must not re-move a pinned step
#   case9_classes    - per-output class list (1 full / 2 station) reaches
#                      OutputSchedule.dat; schedule unchanged by classes
#   case4e           - class count /= #outputs rejected by PMP2init
#   case4f           - optional line whose comment contains "=" rejected by PMP2init
#   case12_frhalo    - f(R) run with BDM: fR_halo rows match the catalogue rows
#                      binary catalogue (.bin + .meta) equals the text column by column
#   case13_catbin    - BDM with CatalogueFormat = binary: .bin + .meta only, no text
#   case4g_badcatfmt - BDM.config CatalogueFormat = hdf5 is rejected at the first catalogue step
#   case11_modes     - low-k density modes from POWERfft5 (P(k) from the cube equals
#                      PowerDM.R bin by bin), momentum modes on N_dep=36, and the
#                      full quantised subsample M=1 (decoded spectra equal the dumps: V6)
#   case10_lens72/36 - lens planes (2 slabs per axis) at plane grid 72 and 36:
#                      slab means average to 1, slab average equals the
#                      ProjectDENSIT full-box projection, 36 = 2x2 block mean
#                      of 72; schedule unchanged
#   every sim case   - OutputSchedule.dat lists exactly the marked steps with
#                      their epochs, comoving distances and classes
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

echo "== case7_closepairs: four pairs closer than 0.5% in a (#outputs = -8)"
prep_case case7_closepairs
set_outputs case7_closepairs -8 " 2.0 1.98541 1.75 1.74475 1.00202 1.0 0.37562 0.375"
run_init case7_closepairs || echo "  (init failed)"
copy_ics case7_closepairs && run_main case7_closepairs 140   # reach the z=0.375 pair (steps ~108/109)

echo "== case8_pinned: merged target followed 0.09% later by an inserted one (#outputs = -2)"
prep_case case8_pinned
set_outputs case8_pinned -2 " 0.49337 0.4920"
run_init case8_pinned || echo "  (init failed)"
copy_ics case8_pinned && run_main case8_pinned

echo "== case9_classes: per-output classes + OutputSchedule.dat (#outputs = -3, classes 2 1 2)"
prep_case case9_classes
set_outputs case9_classes -3 " 2.5 0.512 1.1"
{ echo 'Output classes            =         3  0-all class 1'
  echo ' 2 1 2'
  echo 'Lens slabs per axis       =         0  0-off'
  echo 'Lens plane grid           =         0  0-Ngrid'
} >> "$WORK/case9_classes/Init.dat"
run_init case9_classes || echo "  (init failed)"
copy_ics case9_classes && run_main case9_classes

echo "== case10_lens72 / case10_lens36: lens planes at full and half plane resolution (+ ProjectDENSIT reference)"
for c in case10_lens72:72 case10_lens36:36; do
    d=${c%%:*}; g=${c##*:}
    prep_case "$d"
    set_outputs "$d" -3 " 2.5 0.512 1.1"
    sed -i 's/^Dump random fraction      =         0/Dump random fraction      =         1/' "$WORK/$d/Init.dat"
    [ "$d" = case10_lens72 ] && sed -i 's/^DM power spectrum         =         0/DM power spectrum         =         1/' "$WORK/$d/Init.dat"   # L3: P_2D vs 3-D P(k)
    { echo 'Output classes            =         0  0-all class 1'
      echo 'Lens slabs per axis       =         2  0-off'
      echo "Lens plane grid           =        $g  0-Ngrid"
    } >> "$WORK/$d/Init.dat"
    run_init "$d" || echo "  (init failed)"
    copy_ics "$d" && run_main "$d"
done

echo "== case11_modes: low-k density/momentum modes + full quantised subsample (#outputs = -3, iPower=1, kc 0.5, N_dep 36)"
prep_case case11_modes
set_outputs case11_modes -3 " 2.5 0.512 1.1"
sed -i 's/^DM power spectrum         =         0/DM power spectrum         =         1/' "$WORK/case11_modes/Init.dat"
{ echo 'Output classes            =         0  0-all class 1'
  echo 'Lens slabs per axis       =         0  0-off'
  echo 'Lens plane grid           =         0  0-Ngrid'
  echo 'Extra outputs             =         1  0-off'
  echo 'kc_delta                  =     0.500  h/Mpc'
  echo 'kc_mom                    =     0.185  h/Mpc'   # below k_N(N_dep=36)/4 = 0.1885 (no interlacing)
  echo 'N_dep_mom                 =        36  0-Ngrid/2'
  echo 'Subsample M               =         1  0-off; 1-full quantised snapshot'
  echo 'fR at halos               =         0  0-off'
} >> "$WORK/case11_modes/Init.dat"
run_init case11_modes || echo "  (init failed)"
copy_ics case11_modes && run_main case11_modes

echo "== case12_frhalo: f(R) (F5) run with BDM and f_R at halo centres (#outputs = -2)"
prep_case case12_frhalo
set_outputs case12_frhalo -2 " 1.1 0.512"
sed -i 's/^Find BDM halos            =         0/Find BDM halos            =         1/; s/^MG_flag                   =         0/MG_flag                   =         1/' "$WORK/case12_frhalo/Init.dat"
{ echo 'Output classes            =         0  0-all class 1'
  echo 'Lens slabs per axis       =         0  0-off'
  echo 'Lens plane grid           =         0  0-Ngrid'
  echo 'Extra outputs             =         1  0-off'
  echo 'kc_delta                  =     0.500  h/Mpc'
  echo 'kc_mom                    =     0.185  h/Mpc'   # below k_N(N_dep=36)/4 = 0.1885 (no interlacing)
  echo 'N_dep_mom                 =        36  0-Ngrid/2'
  echo 'Subsample M               =         0  0-off'
  echo 'fR at halos               =         1  0-off'
} >> "$WORK/case12_frhalo/Init.dat"
cp "$ROOT/test/BDM.config" "$WORK/case12_frhalo/Run1/BDM.config"
mkdir -p "$WORK/case12_frhalo/Run1/CATALOGS"     # BDM needs it
run_init case12_frhalo || echo "  (init failed)"
copy_ics case12_frhalo && run_main case12_frhalo

echo "== case13_catbin: LCDM run with BDM and CatalogueFormat = binary (no text catalogue; #outputs = -1)"
prep_case case13_catbin
set_outputs case13_catbin -1 " 0.512"
sed -i 's/^Find BDM halos            =         0/Find BDM halos            =         1/' "$WORK/case13_catbin/Init.dat"
{ cat "$ROOT/test/BDM.config"; echo; echo '          CatalogueFormat = binary     ! text, binary or both (default both)'; } > "$WORK/case13_catbin/Run1/BDM.config"
mkdir -p "$WORK/case13_catbin/Run1/CATALOGS"
run_init case13_catbin || echo "  (init failed)"
copy_ics case13_catbin && run_main case13_catbin

# ---------------------------------------------------------------- error cases
echo "== case4a_zbig / case4b_zneg: rejected by PMP2init"
prep_case case4a_zbig
set_outputs case4a_zbig -1 " 60.0"
run_init case4a_zbig
prep_case case4b_zneg
set_outputs case4b_zneg -1 " -0.5"
run_init case4b_zneg
echo "== case4f_badopt: optional line with '=' inside its comment: rejected by PMP2init (no silent default)"
prep_case case4f_badopt
set_outputs case4f_badopt -3 " 2.5 0.512 1.1"
echo 'Extra outputs             =         1  (1 = on)' >> "$WORK/case4f_badopt/Init.dat"
run_init case4f_badopt

echo "== case4g_badcatfmt: BDM.config CatalogueFormat = hdf5: rejected by PMP2main at the first catalogue step, no catalogue written"
prep_case case4g_badcatfmt
set_outputs case4g_badcatfmt -1 " 0.512"
sed -i 's/^Find BDM halos            =         0/Find BDM halos            =         1/' "$WORK/case4g_badcatfmt/Init.dat"
{ cat "$ROOT/test/BDM.config"; echo; echo '          CatalogueFormat = hdf5'; } > "$WORK/case4g_badcatfmt/Run1/BDM.config"
mkdir -p "$WORK/case4g_badcatfmt/Run1/CATALOGS"
run_init case4g_badcatfmt
copy_ics case4g_badcatfmt; run_main case4g_badcatfmt

echo "== case4e_classes: class count differs from #outputs: rejected by PMP2init"
prep_case case4e_classes
set_outputs case4e_classes -3 " 2.5 0.512 1.1"
{ echo 'Output classes            =         2  0-all class 1'
  echo ' 2 1'
} >> "$WORK/case4e_classes/Init.dat"
run_init case4e_classes

echo "== case4c_dup / case4d_range: rejected by PMP2main (hand-edited Setup.dat)"
mkdir -p "$WORK/case4c_dup/Run1" "$WORK/case4d_range/Run1"
cp "$WORK/case5a_legacy/Setup.dat" "$WORK/case4c_dup/Setup.dat"
set_setup_outputs "$WORK/case4c_dup/Setup.dat" -2 "   0.51200   0.51200"
cp "$WORK/case5a_legacy/Setup.dat" "$WORK/case4d_range/Setup.dat"
set_setup_outputs "$WORK/case4d_range/Setup.dat" -1 "  60.00000"
copy_ics case4c_dup;   run_main case4c_dup
copy_ics case4d_range; run_main case4d_range

# ------------------------------------------------------------------ verdicts
echo
python3 "$HERE/asserts.py" "$WORK"
rc=$?
exit $rc
