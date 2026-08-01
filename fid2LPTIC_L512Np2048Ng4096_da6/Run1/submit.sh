#!/bin/bash -l
# GLAM run from the 2LPTIC (FML) fid-cosmology IC, seed 2026:
# ic2pm conversion (S_vel=1.0, standard Gadget velocities) then PMP2main.
# Config mirrors ../../fid_LCDM_L512Np2048Ng4096 except z_init=49, da=8e-4
# (same da/a=0.04 accuracy as the fid suite; 123 steps vs their 157).

#SBATCH --ntasks 1
#SBATCH --exclusive
#SBATCH -J 2lpt-L512-da6
#SBATCH -o DE2L_%J.dump
#SBATCH -e DE2L_%J.err
#SBATCH -p cosma8
#SBATCH -A dp004
#SBATCH -t 24:00:00

module unload gnu_comp
module load intel_comp/2024.2.0
export I_MPI_F90=ifx
module load compiler-rt tbb compiler mpi

echo "=== ic2pm conversion start $(date)"
../../ic2pm.exe /cosma8/data/dp203/dc-ruan1/mgglam_claude/MG-GLAM/2LPTIC_Gui/fidcosmo_L512/snap/IC_fid_Np2048_L512. 1.0
echo "=== ic2pm done, PMP2main start $(date)"

../../PMP2main.exe <<< 170

echo "=== done $(date)"
