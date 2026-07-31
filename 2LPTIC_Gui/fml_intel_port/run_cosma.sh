#!/bin/bash -l
# 2048^3 IC generation with the FML LPT generator (Intel build).
# Memory estimate for 2048^3 particles (HEFTParticle = 88 B/particle):
#   particles ~0.76 TB x buffer_factor + FFT grids ~0.3 TB + displacement
#   buffers ~0.2 TB  =>  use >= 2 full cosma8 nodes (2 TB); 3 nodes is safe.
# Before submitting: cp input_production_L1024_Np2048.lua input.lua
#
# sbatch run_cosma.sh
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=2
#SBATCH -t 4:00:00
#SBATCH -p cosma8
#SBATCH -A dp203
#SBATCH -J ICs2048
#SBATCH --output=run_cosma.sh.%j.out
#SBATCH --error=run_cosma.sh.%j.err

code_dir=/cosma8/data/dp203/dc-ruan1/mgglam_claude/MG-GLAM/2LPTIC_Gui/FML/FML/LPT/example

cd $code_dir
source $code_dir/setup_env.sh
export LD_LIBRARY_PATH=/cosma/local/fftw/intel_2024.2.0_intel_mpi_2024.2.0/3.3.10-epyc/lib:/cosma/local/gsl/2.7.1/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

echo "=== Job started ==="
echo "Node list: " ${SLURM_JOB_NODELIST}
echo "Ntasks: " ${SLURM_NTASKS}

time mpirun -np $SLURM_NTASKS $code_dir/test input.lua
