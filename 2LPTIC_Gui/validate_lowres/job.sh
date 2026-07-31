#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=2
#SBATCH --mem=500G
#SBATCH -t 1:00:00
#SBATCH -p cosma8-serial
#SBATCH -A dp203
#SBATCH -J ICval1024
#SBATCH --output=job.%j.out
#SBATCH --error=job.%j.err

code=/cosma8/data/dp203/dc-ruan1/mgglam_claude/MG-GLAM/2LPTIC_Gui/FML/FML/LPT/example/test
setup_env=/cosma8/data/dp203/dc-ruan1/mgglam_claude/MG-GLAM/2LPTIC_Gui/FML/FML/LPT/example/setup_env.sh

cd $SLURM_SUBMIT_DIR
source $setup_env
export LD_LIBRARY_PATH=/cosma/local/fftw/intel_2024.2.0_intel_mpi_2024.2.0/3.3.10-epyc/lib:/cosma/local/gsl/2.7.1/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

echo "=== Node list: ${SLURM_JOB_NODELIST}, ntasks ${SLURM_NTASKS}"
time mpirun -np $SLURM_NTASKS $code input.lua
