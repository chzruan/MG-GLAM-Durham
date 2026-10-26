# In-line Rockstar adapter of MG-GLAM (GNU GPL v3 or later; see rs_inline.h).
#
#   make ROCKSTAR_DIR=<MPI-Rockstar checkout>   build/librs_inline.a (+ the test harness)
#   make stub                                     build/librs_stub.a (no Rockstar code)
#
# Intel toolchain as GLAM (module load intel_comp/2024.2.0; module load compiler-rt tbb compiler mpi),
# so the adapter and GLAM share one OpenMP runtime (libiomp5). Validated against MPI-Rockstar 0ae386f.
CC  = icx
CXX = icpx
B   = build

# Rockstar's struct halo layout (312 bytes) and source flags (MPI-Rockstar src/Makefile)
RS_DEFS  = -DOUTPUT_RVMAX -DOUTPUT_NFW_CHI2 -DOUTPUT_INERTIA_TENSOR
RS_POSIX = -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -D_DEFAULT_SOURCE -D_POSIX_SOURCE \
           -D_POSIX_C_SOURCE=200809L -D_SVID_SOURCE
COMMON   = -O2 -g -march=core-avx2 -mcmodel=medium -fp-model=precise -fno-math-errno -fexceptions -fiopenmp $(RS_DEFS)

ifneq ($(ROCKSTAR_DIR),)
RS_SRC    = $(ROCKSTAR_DIR)/src
RS_COMMIT = $(shell git -C $(ROCKSTAR_DIR) rev-parse --short=12 HEAD 2>/dev/null)$(shell git -C $(ROCKSTAR_DIR) diff --quiet 2>/dev/null || echo -dirty)
RS_FLAGS_STR = $(CC) $(COMMON)
# shadow/ first: assert() that throws, an mpi.h with types only
INC       = -Ishadow -I$(RS_SRC) -I.
CFLAGS_RS = $(COMMON) -std=c99 $(RS_POSIX) $(INC)
CFLAGS_AD = $(COMMON) -std=gnu11 $(INC) -DRS_ROCKSTAR_COMMIT='"$(RS_COMMIT)"' -DRS_BUILD_FLAGS='"$(RS_FLAGS_STR)"'
CXXFLAGS  = $(COMMON) -std=c++11 $(INC)

# the per-group kernel of MPI-Rockstar and what it needs, objects rk_*.o (groupies.c #includes properties.c;
# potential.c and subhalo_metric.c are compiled through the adapter's rs_wrap_*.c files); the adapter's are rs_*.o
RS_KERNEL = groupies fof nfw jacobi hubble universe_time integrate distance config_vars check_syscalls bounds inthash
RS_IO     = read_config stringparse
RS_OBJ    = $(RS_KERNEL:%=$(B)/rk_%.o) $(RS_IO:%=$(B)/rk_io_%.o) $(B)/rk_error.o
AD_RS_OBJ = $(B)/rs_wrap_potential.o $(B)/rs_wrap_subhalo_metric.o $(B)/rs_glue.o $(B)/rs_parents.o \
            $(B)/rs_strict.o $(B)/rs_inline.o
AD_OBJ    = $(B)/rs_fof.o $(B)/rs_write.o $(B)/rs_guard.o

all: $(B)/librs_inline.a harness/rs_harness

$(B)/librs_inline.a: $(RS_OBJ) $(AD_RS_OBJ) $(AD_OBJ)
	ar rcs $@ $^

$(B)/rk_%.o: $(RS_SRC)/%.c | $(B)
	$(CC) $(CFLAGS_RS) -c -o $@ $<
$(B)/rk_io_%.o: $(RS_SRC)/io/%.c | $(B)
	$(CC) $(CFLAGS_RS) -c -o $@ $<
$(B)/rk_error.o: $(RS_SRC)/error.cpp | $(B)
	$(CXX) $(CXXFLAGS) -c -o $@ $<
# adapter files that include Rockstar sources or headers
$(B)/rs_wrap_%.o: rs_wrap_%.c | $(B)
	$(CC) $(CFLAGS_RS) -c -o $@ $<
$(B)/rs_glue.o $(B)/rs_parents.o $(B)/rs_strict.o $(B)/rs_inline.o: $(B)/%.o: %.c rs_inline_internal.h rs_record.h rs_fof.h rs_convert.h rs_inline.h | $(B)
	$(CC) $(CFLAGS_AD) $(RS_POSIX) -c -o $@ $<
$(B)/rs_fof.o: rs_fof.c rs_fof.h rs_convert.h | $(B)
	$(CC) $(CFLAGS_AD) -c -o $@ $<
$(B)/rs_write.o: rs_write.c rs_inline_internal.h rs_record.h | $(B)
	$(CC) $(CFLAGS_AD) -c -o $@ $<
$(B)/rs_guard.o: rs_guard.cpp rs_inline.h | $(B)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

harness/rs_harness: harness/rs_harness.c harness/pm_read.c harness/fof_ref.c harness/pm_read.h harness/fof_ref.h \
                    $(B)/librs_inline.a
	$(CC) $(COMMON) -std=gnu11 -Wall -I. -o $@ harness/rs_harness.c harness/pm_read.c harness/fof_ref.c \
	    $(B)/librs_inline.a -lstdc++ -lm
else
all:
	@echo "set ROCKSTAR_DIR=<MPI-Rockstar checkout> (validated: commit 0ae386f), or make stub"; false
endif

stub: $(B)/librs_stub.a
$(B)/librs_stub.a: rs_stub.c rs_inline.h | $(B)
	$(CC) -O2 -g -march=core-avx2 -mcmodel=medium -std=gnu11 -I. -c -o $(B)/rs_stub.o rs_stub.c
	ar rcs $@ $(B)/rs_stub.o

$(B):
	mkdir -p $(B)

clean:
	rm -rf $(B) harness/rs_harness

.PHONY: all stub clean
