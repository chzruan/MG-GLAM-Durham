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

## Verification

Results live in `verify/report.md`. Scripts in `scripts/`.

- **Task A — race detection (bit-match).** `bash scripts/verify_race.sh`
  builds `PMP2main-bitmatch.exe` via `make PMP2main-bitmatch` (pre-Pass-4
  flags: `-O2`, no AVX2/FMA/`fast=1`) so FP reduction order matches the
  unoptimized baseline. Runs 3 steps at `OMP_NUM_THREADS ∈ {1,16,64,128}`
  from a Run11 IC snapshot, tees stdout to `verify/race_T{N}.log`, and
  calls `scripts/compare_residuals.py`. Pass bar: `max|Δres| ≤ 1e-14`
  (expected 0 — F20.12 print precision is the limiter at ~1e-12).
- **Task B — fid_F5 reference P(k).** Already validated in commit
  `b233ea0`: `max|P/P_ref − 1| = 7.6e-6` at z = 1.0, 0.5, 0.25, 0.0 —
  three orders of magnitude below the 0.1% bar.
- **Task C — benchmark sweep.** `bash scripts/bench_sweep.sh` uses the
  production `PMP2main.exe`, runs 5 timesteps × 3 repeats at
  `OMP_NUM_THREADS ∈ {16,32,64,128}`, preserves each run's `timing.log`
  as `verify/bench_T{N}_r{R}.log`, and calls `scripts/bench_report.py`
  to emit a speedup table. Expected: csf and DGP 4–8× at T=128 vs T=16
  (ATOMIC removal); fR 1.3–2× (COLLAPSE + AVX2 only).

### Pass 6 — Red-black stride-2 rewrite of f(R) smoother (2026-04-16)

Applied to all 6 `fr_n × ilevel` specializations of `relaxation_iterations_fR`
in `PMP2MGsolver_fR.f90` (lines 84–290). Replaces the per-cell `CYCLE`-on-parity
mask with an inner loop whose start index encodes red-black selection:

```fortran
parity_off = MERGE(0, 1, redstep)   ! hoisted, once per call

!$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(STATIC) DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3,M1start,M1u,M1l,M2u,M2l,M3u,M3l, …)
DO M3 = 1, ngrid_level
    DO M2 = 1, ngrid_level
        M1start = 1 + MOD(M2 + M3 + parity_off, 2)
    !DIR$ IVDEP
        DO M1 = M1start, ngrid_level, 2
            ! body unchanged — no CYCLE lines
        END DO
    END DO
END DO
```

Parity math: red pass processes `MOD(M1+M2+M3,2)=1`, black pass processes
`=0`. Enumerating, `M1start = 1 + MOD(M2+M3+parity_off, 2)` reproduces the
correct start for all four `{redstep, (M2+M3) mod 2}` cases, so the stride-2
inner loop hits exactly the same cells the CYCLE version touched.

Why this supersedes the Pass-2/3 note ("CYCLE deliberately not rewritten"):
after the prior Phase B hoist split `relaxation_iterations_fR` into 6 body
specializations, the branch removal became mechanical per-specialization,
and COLLAPSE(2) on the outer (M3,M2) pair still yields ≥ 262144 chunks at
`NGRID=512` (262k ÷ 128 threads = 2048 chunks/thread, comfortably above
the scheduler's static-chunk threshold). Coarsest levels (`ngrid_level=2,4`)
drop below the thread count but are already compute-negligible.

IVDEP safety: in one red (or black) pass, every M1 iteration's 6-point
stencil neighbors are opposite-parity cells (not written in this pass), so
stride-2 iterations are independent by construction. `FI2` and `FI3` are
disjoint allocatables (`PMP2mod_tools.f90:41-42`), ruling out cross-array
aliasing.

**Vectorization result** — `PMP2MGsolver_fR.optrpt` after build:

| Specialization (line) | LOOP WAS VECTORIZED | Vector length |
|-----------------------|:-------------------:|:-------------:|
| fr_n=0, ilevel=0 (98) | yes                 | 2             |
| fr_n=0, ilevel>0 (124)| yes                 | 2             |
| fr_n=1, ilevel=0 (151)| yes                 | 4             |
| fr_n=1, ilevel>0 (191)| yes                 | 4             |
| fr_n=2, ilevel=0 (232)| yes                 | 4             |
| fr_n=2, ilevel>0 (268)| yes                 | 4             |

All `#15344 vector dependence` remarks on smoother loops gone. Remaining
`#15344` in `PMP2MGsolver_fR.optrpt` (line 762) is the `restrict_residual_fR`
ilevel>1 branch — no red-black there, untouched by this pass.

`LOOP WAS VECTORIZED` count for `PMP2MGsolver_fR.optrpt`: **18** (down from
the previous 20 because the Phase-A pass dropped 3 non-smoother routines to
COLLAPSE(2), which merges some inner+remainder loops into one; net
vectorization coverage is wider — all smoother bodies now included).

`fr_n=1,2` bodies pick up length-4 SIMD via ifx's SVML-vectorized
DACOS/DCOS; `fr_n=0` hits length-2 because the simpler `DSQRT(P*P-Q4)`
body mixes the Real*4 `FI` load with Real*8 arithmetic temporaries.

## Out of scope (intentionally not done)

- Refactoring the 5 solvers into a shared relaxation template.
- Loop tiling / cache blocking of the inner stencil.
- MPI changes (decomposition, halo exchange, communication overlap).
- Algorithmic multigrid changes (smoother count, V/F/W, coarsest solver).
- Single-precision conversion.
- Particle-side parallelization (CIC deposit, force interpolation).
- Applying the same stride-2 rewrite to csf/DGP/kmf/sym smoothers (separate
  pass once fid_F5 benchmark confirms the fR wall-time reduction).
