#!/bin/bash -l
#SBATCH -J icval_pk
#SBATCH -o icval_pk_%J.out
#SBATCH -e icval_pk_%J.err
#SBATCH -p cosma8-serial
#SBATCH -A dp203
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=400G
#SBATCH -t 0:45:00

module purge
NBKT=/cosma/apps/durham/dc-ruan1/micromamba/envs/nbkt
export OMP_NUM_THREADS=1
export LD_LIBRARY_PATH=$NBKT/lib:$LD_LIBRARY_PATH

cd $SLURM_SUBMIT_DIR
echo "job $SLURM_JOB_ID, $SLURM_NTASKS tasks, host $(hostname), start $(date)"

$NBKT/bin/mpirun -n $SLURM_NTASKS $NBKT/bin/python3 measure_validate.py \
    "$PWD/snap/IC_ours_Np1024_L1024.*" \
    "/cosma7/data/dp004/bl267/Runs/DEGRACE/ICs/IC_data/L1024/Node_002/ics.*" \
    1024 "$PWD/pk_validate_lowres.npz"

echo "end $(date)"
