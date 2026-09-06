#!/bin/bash
#SBATCH --job-name=bdm-repair
#SBATCH --partition=cosma8-serial
#SBATCH --account=dp004
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=768M
#SBATCH --time=00:05:00
#SBATCH --output=BDM-refine/repairs/20260906/work/native-%j.log
#SBATCH --error=BDM-refine/repairs/20260906/work/native-%j.log
set -euo pipefail
export LINES=40 COLUMNS=120
module purge
module load intel_comp/2024.2.0
module load compiler-rt tbb compiler
export BDM_AUDIT_NATIVE_LIBS="$LD_LIBRARY_PATH"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
micromamba run -n cosemu python3 -B BDM-refine/repairs/20260906/native_validation.py "$@"
