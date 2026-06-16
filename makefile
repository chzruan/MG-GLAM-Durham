#FC = gfortran 
#FFLAGS =   -O2  -g -fbacktrace   -fopenmp   -mcmodel=medium -fconvert=big-endian
#LDFLAGS =  -O2  -g -fbacktrace   -fopenmp   -mcmodel=medium -fconvert=big-endian
# FC = ifort 
FC = ifx
FFLAGS =   -O3 -g -traceback -ftz -unroll -qopenmp -march=core-avx2 -mfma -fp-model fast=1 -qopt-report=2 -qopt-report-phase=vec,openmp -shared-intel -mcmodel=medium -convert big_endian
LDFLAGS =  -O3 -g -traceback -ftz -unroll -qopenmp -march=core-avx2 -mfma -fp-model fast=1 -shared-intel -mcmodel=medium -convert big_endian
FMPI = mpiifort
MPIFLAGS =  -O3 -lmpi -g -traceback -ftz -unroll -qopenmp -march=core-avx2 -mfma -fp-model fast=1 -shared-intel -mcmodel=medium -convert big_endian

# Pre-Pass-4 flags: no -O3, no -march=core-avx2, no -mfma, no -fp-model fast=1, no -qopt-report*.
# Used by the PMP2main-bitmatch target to produce a binary whose FP reduction order is
# identical to the pre-OMP-pass baseline, so Pass-1/2/3 source changes can be verified
# race-free by bit-for-bit residual comparison across thread counts.
FFLAGS_BITMATCH  = -O2 -g -traceback -ftz -unroll -qopenmp -shared-intel -mcmodel=medium -convert big_endian
LDFLAGS_BITMATCH = -O2 -g -traceback -ftz -unroll -qopenmp -shared-intel -mcmodel=medium -convert big_endian

#FC     = pgf77
#FFLAGS = -O3  -mp  -byteswapio -mcmodel=medium -Mlarge_arrays -Mnoframe -Munroll -Knoieee
#LDFLAGS = -O3  -mp   -byteswapio -mcmodel=medium -Mlarge_arrays -Mnoframe -Munroll -Knoieee

OBJ = PMP2mod_tools.o  PMP2mod_fft5.o PMP2mod_random.o PMP2mod_density.o PMP2mod_power.o PMP2linker.o PMP2mod_analyze.o PMP2MG_subroutines.o PMP2mod_MGbackground.o PMP2MGsolver_fR.o  PMP2extradof.o  PMP2MGsolver_DGP.o PMP2MGsolver_sym.o PMP2MGsolver_kmf.o PMP2MGsolver_csf.o

PMP2main: $(OBJ)  PMP2main.o
	$(FC) $(LDFLAGS) -o $@.exe $^

# BDM halo finder (entry point PMP2bdm.f90). Reuse the full OBJ list so every
# module dependency (Tools, Density, Structures/LinkerList, Random/LUXURY,
# ExtradofBackgroundData) is satisfied.
PMP2BDM: $(OBJ) PMP2bdm.o
	$(FC) $(LDFLAGS) -o $@.exe $^

# Gadget-2 -> MG-GLAM PM converter (reuses Tools/WriteDataPM).
gadget2pm: PMP2mod_tools.o gadget2pm.o
	$(FC) $(LDFLAGS) -o $@.exe $^

# Build PMP2main with pre-Pass-4 FP-equivalent flags for the Task A race test.
# Stashes an existing PMP2main.exe (if present), wipes all .o/.mod so the recursive
# build uses bitmatch flags, links to PMP2main.exe, then renames. Restores the
# stashed production binary and cleans intermediate .o/.mod on exit, so this
# target leaves the tree clean regardless of invocation order.
#
# Dep-ordering quirk: the OBJ list puts PMP2mod_density.o BEFORE
# PMP2mod_MGbackground.o, but density USE's ExtradofBackgroundData which is
# defined in MGbackground. So we must pre-build two modules explicitly:
#   1. PMP2mod_tools.o     -> Tools.mod                    (needed by ~all files)
#   2. PMP2mod_MGbackground.o -> ExtradofBackgroundData.mod (needed by density,
#                              extradof, csf, kmf; USEs Tools from step 1)
# With those two .mod files present, the final $(MAKE) PMP2main processes the
# remaining OBJ files in declaration order without missing-module failures.
PMP2main-bitmatch:
	@echo "Building PMP2main-bitmatch.exe with pre-Pass-4 flags (-O2, no AVX2/FMA/fast=1)..."
	@if [ -f PMP2main.exe ]; then mv PMP2main.exe PMP2main.exe.stash; fi
	rm -f *.o *.mod
	+$(MAKE) -j1 FFLAGS="$(FFLAGS_BITMATCH)" LDFLAGS="$(LDFLAGS_BITMATCH)" PMP2mod_tools.o
	+$(MAKE) -j1 FFLAGS="$(FFLAGS_BITMATCH)" LDFLAGS="$(LDFLAGS_BITMATCH)" PMP2mod_MGbackground.o
	+$(MAKE) -j1 FFLAGS="$(FFLAGS_BITMATCH)" LDFLAGS="$(LDFLAGS_BITMATCH)" PMP2main
	mv PMP2main.exe PMP2main-bitmatch.exe
	@if [ -f PMP2main.exe.stash ]; then mv PMP2main.exe.stash PMP2main.exe; fi
	rm -f *.o *.mod

PMP2analyze: $(OBJ)  PMP2analyze.o                 
	$(FC) $(LDFLAGS) -o $@.exe $^                
PMP2spectrum: PMP2mod_tools.o  PMP2mod_fft5.o PMP2mod_random.o PMP2mod_density.o PMP2mod_power.o PMP2mod_analyze.o  PMP2spectrum.o                 
	$(FC) $(LDFLAGS) -o $@.exe $^                
PMP2init: PMP2init.o
	$(FC) $(LDFLAGS)  -o $@.exe $^
PMP2start: PMP2mod_tools.o PMP2mod_fft5.o PMP2mod_random.o PMP2MG_subroutines.o PMP2mod_MGbackground.o PMP2mod_density.o PMP2mod_power.o PMP2start.o
	$(FC) $(LDFLAGS) -o $@.exe $^
PMP2correlation: PMP2mod_tools.o  PMP2mod_random.o PMP2correlation.o
	$(FC) $(LDFLAGS) -o $@.exe  $^
PMP2density: PMP2mod_tools.o PMP2MG_subroutines.o PMP2mod_MGbackground.o PMP2density.o
	$(FC) $(LDFLAGS) -o $@.exe $^
PMP2startMPI.o: PMP2startMPI.f90
	$(FMPI) $(MPIFLAGS) -c PMP2startMPI.f90 
PMP2startMPI: $(OBJ) PMP2startMPI.o
	$(FMPI) $(MPIFLAGS)  -o $@.exe $^
PMP2mainMPI.o: PMP2mainMPI.f90
	$(FMPI) $(MPIFLAGS) -c PMP2mainMPI.f90 
PMP2mainMPI: $(OBJ) PMP2mainMPI.o
	$(FMPI) $(MPIFLAGS)  -o $@.exe $^

# Compile MG program
PMP2MG: PMP2mod_tools.o PMP2mod_fft5.o PMP2mod_random.o PMP2mod_density.o PMP2mod_power.o PMP2linker.o PMP2mod_analyze.o PMP2MG_subroutines.o PMP2mod_MGbackground.o PMP2MGsolver_fR.o  PMP2extradof.o PMP2main.o PMP2MGsolver_DGP.o PMP2MGsolver_sym.o PMP2MGsolver_kmf.o PMP2MGsolver_csf.o 
	$(FC) $(LDFLAGS) -o $@.exe $^

.f.o: 
	$(FC) -c $(FFLAGS) $<

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(FC) $(FFLAGS) -w -c $<

clean:
	rm -f *.o *.exe *.mod
