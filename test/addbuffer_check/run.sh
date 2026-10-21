#!/bin/bash -l
# Bit-identity check of the halo finder's two particle-buffering paths (memory audit 2026-10-21, cut 4):
# RetainParticles + AddBuffer(rescaleRetained=.true.) against RescaleCoords(1) + AddBuffer, exact restoration
# and a balanced memory counter, at 1 and 4 threads. Links against the objects of a full build in the repo root
# (make PMP2main first). Light: 300,000 particles, a few seconds.
set -e
cd "$(dirname "$0")"
R=../..
module purge; module load intel_comp/2024.2.0; export I_MPI_F90=ifx; module load compiler-rt tbb compiler mpi
OBJ="PMP2mod_tools.o PMP2mod_fft5.o PMP2mod_random.o PMP2mod_density.o PMP2mod_extras.o PMP2mod_subsample.o glam_hash.o glam_kmp_alloc.o PMP2mod_power.o PMP2linker.o PMP2mod_analyze.o PMP2MG_subroutines.o PMP2mod_MGbackground.o PMP2MGsolver_fR.o PMP2extradof.o PMP2MGsolver_DGP.o PMP2MGsolver_sym.o PMP2MGsolver_kmf.o PMP2MGsolver_csf.o"
objs=""; for o in $OBJ; do objs="$objs $R/$o"; done
ifx -O3 -g -traceback -ftz -unroll -qopenmp -march=core-avx2 -mfma -fp-model fast=1 -shared-intel -mcmodel=medium -convert big_endian -I$R -w addbuffer_check.f90 $objs -o addbuffer_check.exe
for t in 1 4; do echo "== threads $t"; OMP_NUM_THREADS=$t ./addbuffer_check.exe | grep -E 'ADDBUFFER_CHECK|FAIL|originals'; done
rm -f fort.17 addbuffer_check.exe
