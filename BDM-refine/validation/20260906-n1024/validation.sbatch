#!/bin/bash
#SBATCH --job-name=bdm-validation
#SBATCH --partition=cosma8-serial
#SBATCH --account=dp004
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=BDM-refine/validation/20260906-n1024/work/slurm-%j.log
#SBATCH --error=BDM-refine/validation/20260906-n1024/work/slurm-%j.log
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
micromamba run -n cosemu python3 -B "$1" "${@:2}" --threads "$SLURM_CPUS_PER_TASK"
