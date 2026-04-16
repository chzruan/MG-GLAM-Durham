# OpenMP Optimization of the MG-GLAM Modified-Gravity Solvers

Applied 2026-04-15. Plan: `~/.claude/plans/rustling-wobbling-wozniak.md`.

## Summary of changes

### Pass 1 — Remove spurious `!$OMP ATOMIC`
- `PMP2MGsolver_csf.f90` (around lines 148–156)
- `PMP2MGsolver_DGP.f90` (around lines 223–230)

The red-black Gauss-Seidel sweep writes each `(M1,M2,M3)` cell at most once
per call, and different threads own disjoint `M3` slabs, so no race exists.
`fR / kmf / sym` already used the direct write — csf/DGP are now consistent.
ATOMIC was forcing a `lock cmpxchg` per inner-loop store and serializing the
cache line.

### Pass 2 + Pass 3 — Add `COLLAPSE(3) SCHEDULE(STATIC)`
Applied to all 43 active triple-nested `!$OMP PARALLEL DO` directives across:

- `PMP2MGsolver_fR.f90` (6 sites)
- `PMP2MGsolver_csf.f90` (6 sites)
- `PMP2MGsolver_DGP.f90` (6 sites)
- `PMP2MGsolver_kmf.f90` (6 sites)
- `PMP2MGsolver_sym.f90` (6 sites)
- `PMP2extradof.f90` (13 sites)

Directive change pattern:
```fortran
! before
!$OMP PARALLEL DO DEFAULT(SHARED) &
! after
!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC) DEFAULT(SHARED) &
```

`COLLAPSE(3)` flattens the `M3 / M2 / M1` iteration space so coarse multigrid
levels (where `ngrid_level` drops to 16, 8, 4 — smaller than the thread count)
keep all threads busy. `SCHEDULE(STATIC)` is correct because per-cell work is
uniform; it removes reliance on the compiler default.

All triple loops were verified perfectly nested (no statements between `DO`
headers, rectangular `1, ngrid_level` bounds). The red-black `CYCLE` skip
remains legal under COLLAPSE — collapsed iterations of the wrong color simply
become no-ops. The CYCLE was deliberately *not* rewritten as a strided start,
because that would make the `M1` bound depend on `M2/M3` and break COLLAPSE
legality.

### Pass 4 — Makefile flag upgrade

`makefile` lines 5–9:

```make
# before
FFLAGS  =   -O2  -g -traceback -ftz -unroll  -qopenmp  -shared-intel -mcmodel=medium -convert big_endian
LDFLAGS =   -O2  -g -traceback -ftz -unroll  -qopenmp  -shared-intel -mcmodel=medium -convert big_endian
MPIFLAGS =  -O2 -lmpi  -g -traceback -ftz -unroll  -qopenmp  -shared-intel -mcmodel=medium -convert big_endian

# after
FFLAGS  =   -O3 -g -traceback -ftz -unroll -qopenmp -march=core-avx2 -mfma -fp-model fast=1 -qopt-report=2 -qopt-report-phase=vec,openmp -shared-intel -mcmodel=medium -convert big_endian
LDFLAGS =   -O3 -g -traceback -ftz -unroll -qopenmp -march=core-avx2 -mfma -fp-model fast=1 -shared-intel -mcmodel=medium -convert big_endian
MPIFLAGS =  -O3 -lmpi -g -traceback -ftz -unroll -qopenmp -march=core-avx2 -mfma -fp-model fast=1 -shared-intel -mcmodel=medium -convert big_endian
```

Per-flag rationale:
- `-O3` — enables loop interchange, software pipelining, more aggressive unrolling.
- `-march=core-avx2` — targets AVX2 (4 doubles per vector). **Cosma8 is AMD
  EPYC 7542 (Zen 2 Rome), which supports AVX2 but NOT AVX512.** Using
  `-march=` instead of Intel's `-xCORE-AVX2` is critical: the `-x` forms
  embed a `GenuineIntel` vendor-ID dispatch check that refuses to run on AMD
  CPUs regardless of the ISA support, producing the runtime error
  *"Please verify that both the operating system and the processor support
  Intel(R) X87, CMOV, MMX, ..."*. The `-march=` / `-m<feature>` forms emit
  the same ISA without the vendor check.
- `-mfma` — enables 3-operand fused multiply-add (Zen 2 has native FMA).
- `-fp-model fast=1` — enables reordering of associative reductions; cheaper transcendentals (DSQRT, DEXP, DSIGN).
- `-qopt-report=2 -qopt-report-phase=vec,openmp` — writes `*.optrpt` next to each object, compile-time only.

This is allowed because the user accepted a **0.1% matter power spectrum**
reproducibility bar (not bit-for-bit), so FP reordering is in budget.

### Pass 5 — Runtime OpenMP environment

`fid_F5/arun.sh` and `fid_LCDM/arun.sh` now export, after the module loads:

```bash
export OMP_NUM_THREADS=128
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export OMP_STACKSIZE=512M
```

Cosma8 nodes are 2 sockets × 64 cores = 128 cores. `--exclusive` jobs get the
whole node. `close + cores` gives best stencil locality. `512M` stack avoids
blowing the default with the MG temporaries.

For an MPI launch (`PMP2mainMPI.exe`), prefer 2 ranks/node × 64 threads so
each rank's sub-domain sits in one socket's local memory.

## Build verification

```
module purge
module load intel_comp/2024.2.0
export I_MPI_F90=ifx
module load compiler-rt tbb compiler mpi
make PMP2mod_MGbackground.o     # build module first (existing dep ordering quirk)
make PMP2main
```

Result: `PMP2main.exe` built cleanly. `LOOP WAS VECTORIZED` counts in the
`*.optrpt` reports for the MG hot path:

| File                     | LOOP WAS VECTORIZED |
|--------------------------|---------------------|
| PMP2MGsolver_fR.optrpt   | 20                  |
| PMP2MGsolver_csf.optrpt  | 21                  |
| PMP2MGsolver_DGP.optrpt  |  3                  |
| PMP2MGsolver_kmf.optrpt  |  7                  |
| PMP2MGsolver_sym.optrpt  |  6                  |
| PMP2extradof.optrpt      |  8                  |

## Verification still TODO

- **3-step regression** for f(R), csf, DGP at `OMP_NUM_THREADS = 1, 16, 64, 128`.
  After Pass 1 + Pass 2 the residual norms should be bit-for-bit identical to
  the unoptimized baseline (no FP reordering yet at `-O2`). Any divergence at
  the bit level beyond `~1e-15` flags a lurking race.
- **fid_F5 reference run** with the new `-O3 -xCORE-AVX512 -fp-model fast=1`
  binary, comparing P(k) at z=0 to the saved baseline. Pass condition:
  `max |P_new(k)/P_ref(k) - 1| < 1e-3` over the resolved k-range.
- **Benchmark sweep**: median of 3 repeats, 5 timesteps (drop step 1), report
  `T(relax)/cycle` from the existing `CALL TimingMain(3,0/1)` brackets at
  `OMP_NUM_THREADS = 16, 32, 64, 128`. Expected: csf and DGP show 4–8× at
  128 threads (ATOMIC removal); fR shows 1.3–2× (COLLAPSE + AVX-512 only).

## Out of scope (intentionally not done)

- Refactoring the 5 solvers into a shared relaxation template.
- Loop tiling / cache blocking of the inner stencil.
- MPI changes (decomposition, halo exchange, communication overlap).
- Algorithmic multigrid changes (smoother count, V/F/W, coarsest solver).
- Single-precision conversion.
- Particle-side parallelization (CIC deposit, force interpolation).
