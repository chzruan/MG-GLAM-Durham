#!/bin/bash
#SBATCH --job-name=bdm-policy-replay
#SBATCH --partition=cosma8-serial
#SBATCH --account=dp004
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:05:00
#SBATCH --output=BDM-refine/analysis/duplicate-hosts-20260905/logs/policy-replay-%j.out
#SBATCH --error=BDM-refine/analysis/duplicate-hosts-20260905/logs/policy-replay-%j.err
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
task_root=BDM-refine/analysis/duplicate-hosts-20260905
micromamba run -n cosemu python3 "$task_root/tests/policy_replay.py" \
 /cosma8/data/dp203/dc-ruan1/DESI_MGx100/data/GR/Run1/CATALOGS/CatshortV.0137.0001.DAT \
 "$task_root/products/strict-v1/DESI_MGx100/data/GR/Run1/CATALOGS/CatshortV.0137.0001.cleaning.hdf5" \
 "$task_root/provenance/policy-replay.json"
