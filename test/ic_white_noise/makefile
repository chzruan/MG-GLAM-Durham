# Standalone test of MG-GLAM's white-noise draw (see ranlux_lags.f90 and
# ../../docs/ic-white-noise-20261019.md).  Builds against the repository's
# own PMP2mod_random.f90 and PMP2mod_fft5.f90 with the production flags.
#   module load intel_comp/2024.2.0 compiler-rt tbb compiler
#   make            -> build/ranlux_lags.exe
#   make GLAM=<tree from before 2026-10-19> B=build_orig ORIG=1
#                   -> the same against the unmodified module, with the stand-in
#                      ranluxlevel_stub.f90 for the missing RanluxLevel()
GLAM = ../..
FC = ifx
FFLAGS = -O3 -g -traceback -ftz -unroll -qopenmp -march=core-avx2 -mfma -fp-model fast=1 -shared-intel -mcmodel=medium -convert big_endian
B = build
ifdef ORIG
STUB = $(B)/ranluxlevel_stub.o
endif

$(B)/ranlux_lags.exe: $(B)/PMP2mod_random.o $(B)/PMP2mod_fft5.o $(B)/ranlux_lags.o $(STUB)
	$(FC) $(FFLAGS) -o $@ $^

$(B)/ranluxlevel_stub.o: ranluxlevel_stub.f90 | $(B)
	$(FC) $(FFLAGS) -module $(B) -c $< -o $@

$(B)/PMP2mod_random.o: $(GLAM)/PMP2mod_random.f90 $(GLAM)/luxuryp.h | $(B)
	$(FC) $(FFLAGS) -I$(GLAM) -module $(B) -c $< -o $@

$(B)/PMP2mod_fft5.o: $(GLAM)/PMP2mod_fft5.f90 | $(B)
	$(FC) $(FFLAGS) -module $(B) -c $< -o $@

$(B)/ranlux_lags.o: ranlux_lags.f90 $(B)/PMP2mod_random.o $(B)/PMP2mod_fft5.o | $(B)
	$(FC) $(FFLAGS) -I$(GLAM) -module $(B) -c $< -o $@

$(B):
	mkdir -p $(B)

clean:
	rm -rf $(B)
.PHONY: clean
