#FC = gfortran 
#FFLAGS =   -O2  -g -fbacktrace   -fopenmp   -mcmodel=medium -fconvert=big-endian
#LDFLAGS =  -O2  -g -fbacktrace   -fopenmp   -mcmodel=medium -fconvert=big-endian
# FC = ifort 
FC = ifx
FFLAGS =   -O3 -g -traceback -ftz -unroll -qopenmp -march=core-avx2 -mfma -fp-model fast=1 -shared-intel -mcmodel=medium -convert big_endian
LDFLAGS =  -O3 -g -traceback -ftz -unroll -qopenmp -march=core-avx2 -mfma -fp-model fast=1 -shared-intel -mcmodel=medium -convert big_endian
# C helper for the hashed particle subsample (glam_hash.c, called through iso_c_binding)
CC = icx
CFLAGS = -O2 -g -march=core-avx2 -mcmodel=medium
# In-line Rockstar catalogues (rockstar_inline/README.md): make ROCKSTAR_DIR=<MPI-Rockstar source> links the
# adapter (GNU GPL v3 or later, as it contains MPI-Rockstar); without it, the stub (no Rockstar code), and GLAM
# refuses Find Rockstar halos = 1 at start-up. The adapter's C++ guard needs -lstdc++.
ROCKSTAR_DIR ?=
RS_SRCS = $(wildcard rockstar_inline/*.c rockstar_inline/*.cpp rockstar_inline/*.h rockstar_inline/shadow/*.h) rockstar_inline/makefile
ifneq ($(ROCKSTAR_DIR),)
RS_LIB = rockstar_inline/build/librs_inline.a
else
RS_LIB = rockstar_inline/build/librs_stub.a
endif
RS_LIBS = -lstdc++
# The finder validates IEEE domains and uses reproducible diagnostic arithmetic.
# Append the finder safety flags after overrides, including recursive bitmatch
# builds. Intel resolves repeated FP-model options in favour of the last one.
BDM_FFLAGS = $(FFLAGS)

# Build everything in the order that satisfies the Fortran module dependencies
# from a clean tree (PMP2init first builds the shared modules; see the OBJ note).
all: PMP2init PMP2start PMP2main PMP2MG ic2pm
.PHONY: all
BDM_PRECISE_FLAGS = $(if $(filter gfortran%,$(notdir $(FC))),-fno-fast-math -ffp-contract=off,-fp-model precise)
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

OBJ = PMP2mod_tools.o  PMP2mod_fft5.o PMP2mod_random.o PMP2mod_density.o PMP2mod_extras.o PMP2mod_subsample.o glam_hash.o glam_kmp_alloc.o PMP2mod_power.o PMP2linker.o PMP2mod_rockstar.o PMP2mod_analyze.o PMP2MG_subroutines.o PMP2mod_MGbackground.o PMP2MGsolver_fR.o  PMP2extradof.o  PMP2MGsolver_DGP.o PMP2MGsolver_sym.o PMP2MGsolver_kmf.o PMP2MGsolver_csf.o

PMP2main: $(OBJ)  PMP2main.o $(RS_LIB)
	$(FC) $(LDFLAGS) -o $@.exe $^ $(RS_LIBS)

rockstar_inline/build/librs_inline.a: $(RS_SRCS)
	$(MAKE) -C rockstar_inline ROCKSTAR_DIR=$(ROCKSTAR_DIR) build/librs_inline.a
rockstar_inline/build/librs_stub.a: rockstar_inline/rs_stub.c rockstar_inline/rs_inline.h
	$(MAKE) -C rockstar_inline stub

# BDM halo finder (entry point PMP2bdm.f90). Reuse the full OBJ list so every
# module dependency (Tools, Density, Structures/LinkerList, Random/LUXURY,
# ExtradofBackgroundData) is satisfied.
PMP2BDM: $(OBJ) PMP2bdm.o $(RS_LIB)
	$(FC) $(LDFLAGS) -o $@.exe $^ $(RS_LIBS)

glam_hash.o: glam_hash.c
	$(CC) $(CFLAGS) -c $<

# Fortran ALLOCATE through glibc malloc instead of libiomp5's caching allocator (memory audit 2026-10-21, cut 1).
glam_kmp_alloc.o: glam_kmp_alloc.c
	$(CC) $(CFLAGS) -c $<

PMP2linker.o: PMP2linker.f90 PMP2mod_tools.o PMP2mod_density.o PMP2mod_extras.o
	$(FC) $(BDM_FFLAGS) $(BDM_PRECISE_FLAGS) -w -c $<

# Gadget-2 -> MG-GLAM PM converter (reuses Tools/WriteDataPM).
gadget2pm: PMP2mod_tools.o glam_kmp_alloc.o gadget2pm.o
	$(FC) $(LDFLAGS) -o $@.exe $^

# 2LPTic (Gadget format-1) IC -> MG-GLAM PM (reuses Tools::ReadSetup + WriteDataPM).
ic2pm: PMP2mod_tools.o glam_kmp_alloc.o ic2pm.o
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

PMP2analyze: $(OBJ)  PMP2analyze.o $(RS_LIB)
	$(FC) $(LDFLAGS) -o $@.exe $^ $(RS_LIBS)
PMP2spectrum: PMP2mod_tools.o  PMP2mod_fft5.o PMP2mod_random.o PMP2mod_density.o PMP2mod_extras.o PMP2mod_subsample.o glam_hash.o glam_kmp_alloc.o PMP2mod_power.o PMP2mod_rockstar.o PMP2mod_analyze.o  PMP2spectrum.o $(RS_LIB)
	$(FC) $(LDFLAGS) -o $@.exe $^ $(RS_LIBS)
PMP2init: PMP2init.o glam_kmp_alloc.o
	$(FC) $(LDFLAGS)  -o $@.exe $^
PMP2start: PMP2mod_tools.o PMP2mod_fft5.o PMP2mod_random.o PMP2MG_subroutines.o PMP2mod_MGbackground.o PMP2mod_density.o PMP2mod_extras.o PMP2mod_power.o glam_kmp_alloc.o PMP2start.o
	$(FC) $(LDFLAGS) -o $@.exe $^
PMP2correlation: PMP2mod_tools.o  PMP2mod_random.o glam_kmp_alloc.o PMP2correlation.o
	$(FC) $(LDFLAGS) -o $@.exe  $^
PMP2density: PMP2mod_tools.o PMP2MG_subroutines.o PMP2mod_MGbackground.o glam_kmp_alloc.o PMP2density.o
	$(FC) $(LDFLAGS) -o $@.exe $^
PMP2startMPI.o: PMP2startMPI.f90
	$(FMPI) $(MPIFLAGS) -c PMP2startMPI.f90 
PMP2startMPI: $(OBJ) PMP2startMPI.o $(RS_LIB)
	$(FMPI) $(MPIFLAGS)  -o $@.exe $^ $(RS_LIBS)
PMP2mainMPI.o: PMP2mainMPI.f90
	$(FMPI) $(MPIFLAGS) -c PMP2mainMPI.f90 
PMP2mainMPI: $(OBJ) PMP2mainMPI.o $(RS_LIB)
	$(FMPI) $(MPIFLAGS)  -o $@.exe $^ $(RS_LIBS)

# Compile MG program
PMP2MG: PMP2mod_tools.o PMP2mod_fft5.o PMP2mod_random.o PMP2mod_density.o PMP2mod_extras.o PMP2mod_subsample.o glam_hash.o glam_kmp_alloc.o PMP2mod_power.o PMP2linker.o PMP2mod_rockstar.o PMP2mod_analyze.o PMP2MG_subroutines.o PMP2mod_MGbackground.o PMP2MGsolver_fR.o  PMP2extradof.o PMP2main.o PMP2MGsolver_DGP.o PMP2MGsolver_sym.o PMP2MGsolver_kmf.o PMP2MGsolver_csf.o $(RS_LIB)
	$(FC) $(LDFLAGS) -o $@.exe $^ $(RS_LIBS)

.f.o: 
	$(FC) -c $(FFLAGS) $<

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(FC) $(FFLAGS) -w -c $<

# The Fortran runtime must import no allocator entry point beyond the set glam_kmp_alloc.c covers
# or knows about (memory audit 2026-10-21, cut 1): a newer libifcoremt that allocates through another
# name would split allocate/free between glibc and libiomp5. Run after a build: make check-alloc-imports
KNOWN_RT_ALLOC_IMPORTS = kmp_aligned_malloc kmp_realloc kmp_free __kmpc_aligned_alloc __kmpc_alloc __kmpc_free kmpc_sharable_malloc kmpc_sharable_free _mm_malloc _mm_free
check-alloc-imports: PMP2main.exe
	@lib=$$(ldd PMP2main.exe | awk '/libifcoremt/{print $$3}'); echo "libifcoremt: $$lib"; \
	 imports=$$(readelf -W --dyn-syms "$$lib" | awk '$$7=="UND" && $$8 ~ /kmp.*(alloc|free)|_mm_malloc|_mm_free|scalable|omp_alloc|omp_free/ {print $$8}' | sed 's/@.*//'); \
	 bad=""; for s in $$imports; do case " $(KNOWN_RT_ALLOC_IMPORTS) " in *" $$s "*) ;; *) bad="$$bad $$s";; esac; done; \
	 if [ -n "$$bad" ]; then echo "UNKNOWN allocator imports in libifcoremt:$$bad -- extend glam_kmp_alloc.c before using this runtime"; exit 1; \
	 else echo "allocator imports of libifcoremt all known: $$(echo $$imports | tr '\n' ' ')"; fi
	@nm -D PMP2main.exe | awk '$$2=="T" && ($$3=="kmp_aligned_malloc"||$$3=="kmp_free"||$$3=="kmp_realloc"||$$3=="kmp_malloc"||$$3=="kmp_calloc")' | wc -l | \
	 awk '{ if ($$1==5) print "PMP2main.exe defines the five kmp_* entry points"; else { print "PMP2main.exe lacks kmp_* definitions"; exit 1 } }'

clean:
	rm -f *.o *.exe *.mod
