#!/bin/bash
# Step 1: 2LPT initial conditions at z = ZINIT with the FML generator (one Slurm job, under a minute).
# Writes $WORK/ICs: the Gadget files IC_L<box>Np<n>_seed<seed>.*, md5sum.txt, the linear modes (_deltaL), the
# input spectrum at z = ZINIT and the job log. Skipped when the ICs are complete; if only the Gadget files were
# deleted, they are regenerated and checked against the recorded md5sum.txt.
set -euo pipefail
source "$(dirname "$0")/config.sh"
trap stop_waiting INT TERM
ICS=$WORK/ICs
PREFIX=IC_L${BOX}Np${NP}_seed${SEED}
NTASKS=8                                   # FML needs Nmesh % ntasks == 0
shopt -s nullglob
n_ics() { local f=("$ICS/$PREFIX".[0-9]*); echo ${#f[@]}; }
if [ -f "$ICS/md5sum.txt" ] && [ "$(n_ics)" -eq "$NTASKS" ]; then echo "ICs already in $ICS: skipped"; exit 0; fi
[ -x "$FML_DIR/test" ] || { echo "FML generator $FML_DIR/test not found: build it first (see README)"; exit 1; }
mkdir -p "$ICS"
cd "$ICS"

# linear P(k) at z = ZINIT: the z = 0 spectrum divided by [D(0)/D(z_ini)]^2 (flat matter + Lambda growth)
$PYTHON - "$PK_Z0" "$OMEGA_M" "$ZINIT" <<'PY'
import sys
import numpy as np
from scipy.integrate import quad
pk_file, om, zini = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
E = lambda x: np.sqrt(om / x**3 + 1 - om)
D = lambda a: 2.5 * om * E(a) * quad(lambda x: 1 / (x * E(x))**3, 0, a, epsabs=0, epsrel=1e-12)[0]
ratio = D(1.0) / D(1.0 / (1.0 + zini))
k, p = np.loadtxt(pk_file, unpack=True)
np.savetxt("pofk_zini.txt", np.c_[k, p / ratio**2], fmt="%.10e")      # no header: FML rejects comment lines
print(f"pofk_zini.txt: D(0)/D({zini:g}) = {ratio:.6f}, {len(k)} rows")
PY

sed -e "s/@BOX@/$BOX/" -e "s/@NP@/$NP/g" -e "s/@SEED@/$SEED/" -e "s/@ZINIT@/$ZINIT/" -e "s/@PREFIX@/$PREFIX/" \
    -e "s/@OMEGA_M@/$OMEGA_M/" -e "s/@HUBBLE@/$HUBBLE/" "$EXAMPLE/input.lua.in" > input.lua

cat > ics_job.sh <<EOF
#!/bin/bash -l
# FML 2LPT ICs for the twin example (written by $EXAMPLE/1_ics.sh)
source "$FML_DIR/setup_env.sh"
export LD_LIBRARY_PATH=/cosma/local/fftw/intel_2024.2.0_intel_mpi_2024.2.0/3.3.10-epyc/lib:/cosma/local/gsl/2.7.1/lib:\${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=1
echo "job \$SLURM_JOB_ID on \$(hostname), \$SLURM_NTASKS tasks, generator md5 \$(md5sum "$FML_DIR/test" | cut -c1-32)"
mpirun -np \$SLURM_NTASKS "$FML_DIR/test" input.lua
EOF
submit_wait ex_ics "$ICS" -p "$PARTITION" -A "$ACCOUNT" --nodes=1 --ntasks=$NTASKS --cpus-per-task=1 --mem=16G \
            -t 00:15:00 -o "$ICS/ics_%j.out" -e "$ICS/ics_%j.err" ics_job.sh

[ "$(n_ics)" -eq "$NTASKS" ] || { echo "expected $NTASKS IC files, found $(n_ics); see $ICS/ics_*.out"; exit 1; }
if [ -f md5sum.txt ]; then
    md5sum -c --quiet md5sum.txt && echo "regenerated ICs are identical to the recorded ones (md5sum.txt)"
else
    md5sum "$PREFIX".[0-9]* > md5sum.txt
fi
# the measured P(k) of the ICs against the input spectrum (2LPT at z = 49: a ratio within ~1 % of 1)
$PYTHON - <<'PY'
import numpy as np
m = np.loadtxt("pofk_ics.txt")          # k, P_linear, P_particle_2LPT, r_cross, P_ini_over_input
w = (m[:, 0] > 0.1) & (m[:, 0] < 1.0)
r = m[w, 4]
print(f"IC power / input spectrum over 0.1 < k < 1 h/Mpc: mean {r.mean():.4f} (min {r.min():.3f}, max {r.max():.3f})")
PY
echo "ICs written to $ICS"
