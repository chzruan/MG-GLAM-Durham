#!/bin/bash -l
#SBATCH --ntasks 64
#SBATCH -t 72:00:00
#SBATCH -p cosma8
#SBATCH -A dp203
#SBATCH -J ICs
#SBATCH --output=run.sh.%j.out
#SBATCH --error=run.sh.%j.err
#SBATCH --mail-type=END # notifications for job done & fail
#SBATCH --mail-user=ydsp26@durham.ac.uk

# Run as:
# sbatch run.sh
code=/cosma8/data/dp203/dc-bran2/new_FML/FML/LPT/example/test
setup_env=/cosma8/data/dp203/dc-bran2/new_FML/FML/LPT/example/setup_env.sh
param_file=input.lua

source $setup_env

echo "=== Job started ==="
echo "Node list: " ${SLURM_JOB_NODELIST}
echo "Ntasks: " ${SLURM_NTASKS}

time mpirun -np $SLURM_NTASKS $code $param_file
