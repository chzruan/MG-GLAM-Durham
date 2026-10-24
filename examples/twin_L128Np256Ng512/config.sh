# Settings of the L128Np256Ng512 twin example, sourced by every step script. To change a setting after a run,
# use a new WORK: the steps skip work that exists already.

# Required: WORK, a new directory for the ICs, the two runs and the figure on the data file system (not $HOME;
# about 0.6 GB while the ICs exist), and ACCOUNT, your Slurm project account.
: "${WORK:?set WORK first, e.g. export WORK=/cosma8/data/<project>/$USER/mgglam-example}"
: "${ACCOUNT:?set ACCOUNT to your Slurm project account, e.g. export ACCOUNT=<project>}"
WORK=$(realpath -m -- "$WORK")
EXAMPLE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO=$(cd "$EXAMPLE/../.." && pwd)
PARTITION=${PARTITION:-cosma8-serial}       # a shared COSMA8 partition

# Python with numpy, scipy, h5py and matplotlib; glamio is taken from this repository
PYTHON=${PYTHON:-python3}
export PYTHON MAMBA_USE_LOCKFILES=false
export PYTHONPATH="$REPO/python${PYTHONPATH:+:$PYTHONPATH}"

# Executables (see README, "Before you start"): GLAM from the repository root, the FML 2LPT generator
GLAM_BIN=${GLAM_BIN:-$REPO}
FML_DIR=${FML_DIR:-$REPO/2LPTIC_Gui/FML/FML/LPT/example}

# The box: the resolution of the 1 Gpc/h production runs (0.5 Mpc/h particle spacing, 0.25 Mpc/h force
# cells) in a 128 Mpc/h volume
BOX=128; NP=256; NG=512
SEED=2026                           # FML random seed
REALISATION=0                       # the run writes into Run0/ and tags its products .0000
ZINIT=49                            # 2LPT start; da0 = 6e-4 (gen_inputs default) suits 0.25 Mpc/h cells
FULL_Z=${FULL_Z:-"1.0 0.5 0.0"}     # full outputs: P(k), lens planes, BDM haloes, mode cubes, subsample
STATION_Z=${STATION_Z-"0.75 0.25"}  # stations: lens planes and P(k) only (may be empty: STATION_Z="")

# Cosmology: Planck 2015 with the CAMB z = 0 linear spectrum shipped here
OMEGA_M=0.3089; HUBBLE=0.6774; SIGMA8=0.8159; OMBH2=0.0223
PK_Z0="$EXAMPLE/PK_CAMB_z00_Planck2015.dat"

# The twins: the same ICs evolved with GR and with Hu-Sawicki f(R), n = 1, fR0 = -1e-5 (first: the LCDM reference)
TWINS="LCDM fR5"
FR0=-1e-5

export BOX NP NG ZINIT REALISATION TWINS FR0           # read by analyse.py
source "$EXAMPLE/lib.sh"
