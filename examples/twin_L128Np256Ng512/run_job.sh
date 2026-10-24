#!/bin/bash -l
# The simulation job of one twin, submitted by 3_run.sh from the twin's run directory with IC (the Gadget file
# prefix) and REALISATION exported: PMP2init (Setup.dat), ic2pm (Gadget ICs -> GLAM particles, velocities at the
# leapfrog half step), then PMP2main or PMP2MG (bin/main_exe.txt). The particle files are deleted at the end and
# Run<N>/job_done.txt marks the finished run.
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=8G
#SBATCH -t 01:00:00
module purge
module load intel_comp/2024.2.0
export I_MPI_F90=ifx
module load compiler-rt tbb compiler mpi
set -euo pipefail
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK OMP_PROC_BIND=close OMP_PLACES=cores
EXE=$(cat bin/main_exe.txt)
echo "job $SLURM_JOB_ID on $(hostname), $OMP_NUM_THREADS threads, $EXE, $(cat bin/git_hash.txt)"
(cd bin && md5sum -c --quiet md5sum.txt) || { echo "bin/ does not match bin/md5sum.txt"; exit 1; }
[ -f "${IC}0" ] || { echo "IC file ${IC}0 missing"; exit 1; }
[ -f Setup.dat ] || ./bin/PMP2init.exe > init.log 2>&1
[ -f Setup.dat ] || { echo "PMP2init failed:"; tail -3 init.log; exit 1; }
cd "Run$REALISATION"
echo "=== ic2pm $(date)"
../bin/ic2pm.exe "$IC" 1.0 half "$REALISATION" > ic2pm.log
grep -E "velocity epoch|rms written" ic2pm.log | head -2 || true
echo "=== $EXE $(date)"
/usr/bin/time -v "../bin/$EXE" <<< 300 > main.log 2> main.time
rm -f PMcrd.DAT PMcrs*.DAT
echo "=== done $(date): $(grep -E 'Elapsed' main.time | sed 's/^\s*//'), peak $(awk '/Maximum resident/{printf "%.2f GB", $NF/1048576}' main.time)"
echo "job $SLURM_JOB_ID finished $(date -Is)" > job_done.txt
echo "outputs: $(ls PowerDM.R.*.dat | wc -l) P(k), $(ls LensPlanes.*.bin | wc -l) lens-plane sets, $(ls CATALOGS/CatshortV.*.bin | wc -l) catalogues, $(ls deltam.*.bin | wc -l) density cubes, $(ls Subsample.*.bin | wc -l) subsamples, $(ls fR_halo.*.bin 2>/dev/null | wc -l) fR_halo"
