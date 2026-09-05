#!/bin/bash
#SBATCH --job-name=bdm-dup-clean-pilot
#SBATCH --partition=cosma8-serial
#SBATCH --account=dp004
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=BDM-refine/analysis/duplicate-hosts-20260905/logs/clean-pilot-%j.out
#SBATCH --error=BDM-refine/analysis/duplicate-hosts-20260905/logs/clean-pilot-%j.err
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
task_root=BDM-refine/analysis/duplicate-hosts-20260905
product_root=$task_root/products/strict-v1
micromamba run -n cosemu python3 "$task_root/tools/clean_catalogue.py" \
 /cosma8/data/dp203/dc-ruan1/DESI_MGx100/data/GR/Run1/CATALOGS/CatshortV.0137.0001.DAT \
 "$product_root/DESI_MGx100/data/GR/Run1/CATALOGS/CatshortV.0137.0001.DAT" \
 "$product_root/DESI_MGx100/data/GR/Run1/CATALOGS/CatshortV.0137.0001.cleaning.hdf5" \
 --gravity LCDM --imodel 0 --ibox 1 --redshift 0.25 --velocities
