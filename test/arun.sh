# Symmetron (MG_model=3) example run for MG-GLAM.
# Prerequisite: build PMP2start.exe and the MG binary PMP2MG.exe in the repo
# root (e.g.  make PMP2start PMP2MG ). The symmetron config lives in ../Setup.dat
# (MG flag=1, MG_model=3, a_*=0.5, xi=1e-3, beta_*=0.1). See README.md.
for box in {2..2}
do
    
    mkdir -p ./Run${box}/CATALOGS/ # the folder CATALOGS is for storing halo catalogues
    cp ./BDM.config  ./Run${box}/ # halo finder config file

    cd ./Run${box}/

cat <<EOF >submit.sh
#!/bin/bash -l

#SBATCH --ntasks 1
#SBATCH -J symA_z1
#SBATCH -o _symA_%J.dump
#SBATCH -e _symA_%J.err
#SBATCH -p cosma8
#SBATCH -A dp203
#SBATCH -t 36:00:00

module unload gnu_comp
module load intel_comp/2024.2.0
export I_MPI_F90=ifx
module load compiler-rt tbb compiler mpi

export OMP_NUM_THREADS=\${SLURM_CPUS_PER_TASK:-\${SLURM_CPUS_ON_NODE:-1}}
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export OMP_STACKSIZE=512M

../../PMP2start.exe <<< ${box}

../../PMP2MG.exe <<< 157


EOF
    sbatch submit.sh

    cd ../



done
