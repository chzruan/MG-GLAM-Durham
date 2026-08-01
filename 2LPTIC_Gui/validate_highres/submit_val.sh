#!/bin/bash -l
#SBATCH -J icval_pk2048
#SBATCH -o icval_pk_%J.out
#SBATCH -e icval_pk_%J.err
#SBATCH -p cosma8-serial
#SBATCH -A dp203
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem=900G
#SBATCH -t 1:30:00

module purge
NBKT=/cosma/apps/durham/dc-ruan1/micromamba/envs/nbkt
export OMP_NUM_THREADS=1
export LD_LIBRARY_PATH=$NBKT/lib:$LD_LIBRARY_PATH

cd $SLURM_SUBMIT_DIR
echo "job $SLURM_JOB_ID, $SLURM_NTASKS tasks, host $(hostname), start $(date)"

$NBKT/bin/mpirun -n $SLURM_NTASKS $NBKT/bin/python3 ../validate_lowres/measure_validate.py \
    "$PWD/snap/IC_ours_Np2048_L1024.*" \
    "/cosma8/data/dp203/bl267/Projects/Ongoing/HEFT/ICs/IC_highres/IC_Np1d_2048_L_1024_2LPT.*" \
    2048 "$PWD/pk_validate_highres.npz"

echo "end $(date)"
