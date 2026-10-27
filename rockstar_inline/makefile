# In-line Rockstar adapter of MG-GLAM (GNU GPL v3 or later; see rs_inline.h).
#
#   make ROCKSTAR_DIR=<MPI-Rockstar checkout>   build/librs_inline.a, the test harness and the FoF test
#   make stub                                     build/librs_stub.a (no Rockstar code)
#
# Intel toolchain as GLAM (module load intel_comp/2024.2.0; module load compiler-rt tbb compiler mpi),
# so the adapter and GLAM share one OpenMP runtime (libiomp5). Validated against MPI-Rockstar 0ae386f.
# A program linking librs_inline.a also needs -lstdc++ (the adapter's guard is C++).
CC  = icx
CXX = icpx
B   = build

# Rockstar's struct halo layout (312 bytes) and source flags (MPI-Rockstar src/Makefile)
RS_DEFS  = -DOUTPUT_RVMAX -DOUTPUT_NFW_CHI2 -DOUTPUT_INERTIA_TENSOR
RS_POSIX = -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -D_DEFAULT_SOURCE -D_POSIX_SOURCE \
           -D_POSIX_C_SOURCE=200809L -D_SVID_SOURCE
COMMON   = -O2 -g -march=core-avx2 -mcmodel=medium -fp-model=precise -fno-math-errno -fexceptions -fiopenmp $(RS_DEFS)
DEPS     = -MMD -MP
# the FoF's distance test is plain float arithmetic: no fused multiply-add (rs_fof.h)
NOFMA    = -ffp-contract=off

ifneq ($(ROCKSTAR_DIR),)
RS_SRC    = $(ROCKSTAR_DIR)/src
RS_GIT    = $(shell git -C $(ROCKSTAR_DIR) rev-parse --short=12 HEAD 2>/dev/null)
RS_DIRTY  = $(shell test -n "$$(git -C $(ROCKSTAR_DIR) status --porcelain --untracked-files=no 2>/dev/null)" && echo -dirty)
RS_COMMIT = $(if $(RS_GIT),$(RS_GIT)$(RS_DIRTY),unknown)
RS_FLAGS_STR = $(shell $(CC) --version 2>/dev/null | head -1): $(CC) $(COMMON) $(RS_POSIX); -std=c99 (Rockstar), -std=gnu11 (adapter), $(NOFMA) (FoF)
# shadow/ first: assert() that throws, an mpi.h with types only
INC       = -Ishadow -I$(RS_SRC) -I. -I$(B)
CFLAGS_RS = $(COMMON) $(DEPS) -std=c99 $(RS_POSIX) $(INC)
CFLAGS_AD = $(COMMON) $(DEPS) -std=gnu11 $(INC)
CXXFLAGS  = $(COMMON) $(DEPS) -std=c++11 $(INC)

# the per-group kernel of MPI-Rockstar and what it needs, objects rk_*.o (groupies.c #includes properties.c;
# potential.c and subhalo_metric.c are compiled through the adapter's rs_wrap_*.c files); the adapter's are rs_*.o
RS_KERNEL = groupies fof nfw jacobi hubble universe_time integrate distance config_vars check_syscalls bounds inthash
RS_IO     = read_config stringparse
RS_OBJ    = $(RS_KERNEL:%=$(B)/rk_%.o) $(RS_IO:%=$(B)/rk_io_%.o) $(B)/rk_error.o
AD_RS_OBJ = $(B)/rs_wrap_potential.o $(B)/rs_wrap_subhalo_metric.o $(B)/rs_glue.o $(B)/rs_parents.o \
            $(B)/rs_strict.o
AD_OBJ    = $(B)/rs_fof.o $(B)/rs_write.o $(B)/rs_guard.o
LIB_OBJ   = $(RS_OBJ) $(AD_RS_OBJ) $(AD_OBJ) $(B)/rs_inline.o
# the harness links the same objects, except rs_inline.c built with its failure-injection hooks
HARNESS_OBJ = $(RS_OBJ) $(AD_RS_OBJ) $(AD_OBJ) $(B)/rs_inline_hooks.o

all: $(B)/librs_inline.a harness/rs_harness harness/test_fof

$(B)/librs_inline.a: $(LIB_OBJ)
	rm -f $@
	ar rcs $@ $^

# what every object is built from (the MPI-Rockstar checkout, its commit, the compilers and all flags): rewritten
# only when that changes, and then every object is rebuilt, so the stamp below never describes stale objects
BUILD_TEXT = $(abspath $(ROCKSTAR_DIR)) $(RS_COMMIT) | $(shell $(CC) --version 2>/dev/null | head -1) | $(CC) $(CFLAGS_RS) | $(CFLAGS_AD) $(NOFMA) | $(CXX) $(CXXFLAGS)
$(B)/build.stamp: FORCE | $(B)
	@echo '$(BUILD_TEXT)' > $@.new
	@if cmp -s $@.new $@; then rm -f $@.new; else mv $@.new $@; fi

# the Rockstar commit and the compiler flags, rewritten only when they change (rs_inline.c includes it)
$(B)/rs_version.h: FORCE | $(B)
	@printf '#define RS_ROCKSTAR_COMMIT "%s"\n#define RS_BUILD_FLAGS "%s"\n' "$(RS_COMMIT)" "$(RS_FLAGS_STR)" > $@.new
	@if cmp -s $@.new $@; then rm -f $@.new; else mv $@.new $@; fi

$(B)/rk_%.o: $(RS_SRC)/%.c makefile $(B)/build.stamp | $(B)
	$(CC) $(CFLAGS_RS) -c -o $@ $<
$(B)/rk_io_%.o: $(RS_SRC)/io/%.c makefile $(B)/build.stamp | $(B)
	$(CC) $(CFLAGS_RS) -c -o $@ $<
$(B)/rk_error.o: $(RS_SRC)/error.cpp makefile $(B)/build.stamp | $(B)
	$(CXX) $(CXXFLAGS) -c -o $@ $<
# adapter files that include Rockstar sources or headers
$(B)/rs_wrap_%.o: rs_wrap_%.c makefile $(B)/build.stamp | $(B)
	$(CC) $(CFLAGS_RS) -c -o $@ $<
$(B)/rs_glue.o $(B)/rs_parents.o $(B)/rs_strict.o: $(B)/%.o: %.c makefile $(B)/build.stamp | $(B)
	$(CC) $(CFLAGS_AD) $(RS_POSIX) -c -o $@ $<
$(B)/rs_inline.o: rs_inline.c $(B)/rs_version.h makefile $(B)/build.stamp | $(B)
	$(CC) $(CFLAGS_AD) $(RS_POSIX) -c -o $@ $<
$(B)/rs_inline_hooks.o: rs_inline.c $(B)/rs_version.h makefile $(B)/build.stamp | $(B)
	$(CC) $(CFLAGS_AD) $(RS_POSIX) -DRS_INLINE_TEST_HOOKS -c -o $@ $<
$(B)/rs_fof.o: rs_fof.c makefile $(B)/build.stamp | $(B)
	$(CC) $(CFLAGS_AD) $(NOFMA) -c -o $@ $<
$(B)/rs_write.o: rs_write.c makefile $(B)/build.stamp | $(B)
	$(CC) $(CFLAGS_AD) -c -o $@ $<
$(B)/rs_guard.o: rs_guard.cpp makefile $(B)/build.stamp | $(B)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

harness/rs_harness: harness/rs_harness.c harness/pm_read.c harness/fof_ref.c harness/pm_read.h harness/fof_ref.h \
                    $(HARNESS_OBJ) makefile
	$(CC) $(COMMON) $(NOFMA) -std=gnu11 -Wall -I. -o $@ harness/rs_harness.c harness/pm_read.c harness/fof_ref.c \
	    $(HARNESS_OBJ) -lstdc++ -lm
harness/test_fof: harness/test_fof.c $(B)/rs_fof.o rs_fof.h rs_convert.h makefile
	$(CC) $(COMMON) $(NOFMA) -std=gnu11 -Wall -I. -o $@ harness/test_fof.c $(B)/rs_fof.o -lm

-include $(wildcard $(B)/*.d)
else
all:
	@echo "set ROCKSTAR_DIR=<MPI-Rockstar checkout> (validated: commit 0ae386f), or make stub"; false
endif

stub: $(B)/librs_stub.a
$(B)/librs_stub.a: rs_stub.c rs_inline.h | $(B)
	$(CC) -O2 -g -march=core-avx2 -mcmodel=medium -std=gnu11 -I. -c -o $(B)/rs_stub.o rs_stub.c
	rm -f $@
	ar rcs $@ $(B)/rs_stub.o

$(B):
	mkdir -p $(B)

clean:
	rm -rf $(B) harness/rs_harness harness/test_fof

FORCE:
.PHONY: all stub clean FORCE
