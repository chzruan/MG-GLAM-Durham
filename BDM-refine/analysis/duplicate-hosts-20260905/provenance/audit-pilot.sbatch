#!/bin/bash
#SBATCH --job-name=bdm-dup-audit
#SBATCH --partition=cosma8-serial
#SBATCH --account=dp004
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=BDM-refine/analysis/duplicate-hosts-20260905/logs/audit-%j.out
#SBATCH --error=BDM-refine/analysis/duplicate-hosts-20260905/logs/audit-%j.err
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
micromamba run -n cosemu python3 BDM-refine/analysis/duplicate-hosts-20260905/tools/audit_catalogue.py \
 /cosma8/data/dp203/dc-ruan1/DESI_MGx100/data/GR/Run1/CATALOGS/CatshortV.0137.0001.DAT \
 BDM-refine/analysis/duplicate-hosts-20260905/provenance/box1-audit.json --velocities
