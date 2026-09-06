#!/bin/bash
#SBATCH --job-name=bdm-n1024-plots
#SBATCH --partition=cosma8-serial
#SBATCH --account=dp004
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=BDM-refine/validation/20260906-n1024/work/slurm-%j.log
#SBATCH --error=BDM-refine/validation/20260906-n1024/work/slurm-%j.log
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
ulimit -c 0
micromamba run -n cosemu python3 -B "$1"
