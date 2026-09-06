#!/bin/bash
# Resources and output paths are supplied by the sealed submitter.
set -euo pipefail
export LINES=40 COLUMNS=120
module purge
module load intel_comp/2024.2.0
module load compiler-rt tbb compiler
export BDM_AUDIT_NATIVE_LIBS="$LD_LIBRARY_PATH"
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OMP_DYNAMIC=FALSE OMP_PROC_BIND=close OMP_PLACES=cores
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
ulimit -c 0
micromamba run -n cosemu python3 -B "$1" --plan-sha256 "$2" --threads "$SLURM_CPUS_PER_TASK"
