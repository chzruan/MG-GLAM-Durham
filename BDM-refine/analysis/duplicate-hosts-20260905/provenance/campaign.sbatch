#!/bin/bash
#SBATCH --job-name=bdm-dup-clean-z025
#SBATCH --partition=cosma8-serial
#SBATCH --account=dp004
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --array=0-7%8
#SBATCH --output=BDM-refine/analysis/duplicate-hosts-20260905/logs/campaign-%A_%a.out
#SBATCH --error=BDM-refine/analysis/duplicate-hosts-20260905/logs/campaign-%A_%a.err
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
task_root=BDM-refine/analysis/duplicate-hosts-20260905
micromamba run -n cosemu python3 "$task_root/tools/run_campaign.py" \
 "$task_root/provenance/z0p25_inventory.json" "$task_root/products/strict-v1" --workers 8
