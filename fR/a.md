# f(R) Solver Optimization Summary

Applied 2026-04-16. Plan: `optimize-the-running-speed-async-crane.md`
(also copied into this folder).

Target: reduce wall time of the f(R) modified-gravity multigrid solver
(`PMP2MGsolver_fR.f90`) without touching the multigrid algorithm, array
precision, or MPI. Baseline is the prior OpenMP pass from
`OMP_OPTIMIZATION_NOTES.md` (COLLAPSE(3), AVX2+FMA, -O3 -fp-model fast=1,
`OMP_PROC_BIND=close`, `OMP_PLACES=cores`).

Reproducibility budget: 0.1% matter power spectrum (same as prior pass),
so associative FP reordering is permitted.

## Summary of changes (`PMP2MGsolver_fR.f90`)

### Phase B — Hoist `fr_n` / `ilevel` dispatch above the triple-nested loop

For each hot subroutine, the per-cell `IF (fr_n == 0/1/2)` and
`IF (ilevel == 0)` branches were loop-invariant across a whole call but
still issued a branch in every cell, suppressing ifx's SIMD pattern
matching. Phase B moves the dispatch above the OMP parallel region, so
each specialization contains a clean straight-line cell body.

| Subroutine | Specializations produced | Notes |
|---|---|---|
| `relaxation_iterations_fR` | 6 (3 fr_n × 2 ilevel) | Red-black CYCLE preserved; no IVDEP (smoother not Phase-A-targetable) |
| `calculate_residual_fR` | 6 (3 fr_n × 2 ilevel) | See Phase E below |
| `restrict_residual_fR` | 3 for ilevel==1 branch | ilevel>1 branch has no fr_n dispatch — untouched |
| `calculate_physical_right_hand_side_fR` | 3 (fr_n only) | Single ilevel path |

### Phase A — Enable inner-M1 vectorization (IVDEP on non-smoother routines)

After Phase B, the inner M1 loops in the three non-smoother routines are
straight-line. Added `!DIR$ IVDEP` immediately above each inner M1 loop
in `calculate_residual_fR`, `restrict_residual_fR` (ilevel==1 branch),
and `calculate_physical_right_hand_side_fR`. This tells ifx that writes
to `FI3(M1,...)` do not alias reads from other offsets into `FI3`, which
is true by construction but which ifx was treating conservatively.

Implementation detail: `COLLAPSE(3)` is incompatible with `!DIR$ IVDEP`
between the `M2` and `M1` DO headers (the perfect-nesting requirement
of COLLAPSE(3) makes the inner loop disappear into the flattened
iteration space, leaving IVDEP with no target). The three non-smoother
routines therefore drop to `COLLAPSE(2)` on the outer pair — still
≥ 256 collapsed chunks for `ngrid_level ≥ 16`, well above the thread
count on typical Cosma8 runs. The smoother keeps `COLLAPSE(3)` because
the `CYCLE` guard on red-black parity blocks vectorization anyway.

Inlined `x**2` as `x*x` and `x**3` as `x*x*x` throughout the specialized
inner bodies (Intel ifx did not always recognize the integer-power
pattern for 6+ `**` per cell).

### Phase E — Fuse residual write with L² reduction (ilevel==0 only)

In the original `calculate_residual_fR`, the ilevel==0 branch wrote
`OP` to `FI3(M1,M2,M3)` in one triple loop, then a second triple loop
read `FI3(N1,N2,N3)**2` and summed it into `RES` for the L² norm. This
second pass over the full PM grid (`Real*4` FI3, ~512 MB at NGRID=512)
is pure DRAM traffic.

Phase E merges the two: the value `OP` written to `FI3` at
`(M1,M2,M3)` is exactly what the second loop later squared, so we can
do `RES = RES + OP*OP` inside the main loop with
`REDUCTION(+:RES)` on the directive. The second DO N3/N2/N1 loop and
its `IF (ilevel == 0)` guard are deleted. `RES = DSQRT(RES/ngrid³)`
and `res_PM_grid = RES` are preserved.

Correctness: the transformation is bit-identical within a call
(`OP*OP == FI3(N1,N2,N3)**2` with the same rounding). Associative
reorder of the reduction is already permitted under `-fp-model fast=1`.

### Phase C — Red-black strided M1 (APPLIED, 2026-04-16)

Rewrote all 6 specializations of `relaxation_iterations_fR` to stride-2
inner loops with `M1start = 1 + MOD(M2 + M3 + parity_off, 2)`, where
`parity_off = MERGE(0, 1, redstep)` is hoisted once per call. Dropped
`COLLAPSE(3) → COLLAPSE(2)` on the outer (M3, M2) pair and added
`!DIR$ IVDEP` on the inner M1 loop.

All 6 smoother inner M1 loops now vectorize:

| fr_n | ilevel | line | vector length |
|-----:|-------:|-----:|--------------:|
| 0    | 0      | 98   | 2             |
| 0    | >0     | 124  | 2             |
| 1    | 0      | 151  | 4             |
| 1    | >0     | 191  | 4             |
| 2    | 0      | 232  | 4             |
| 2    | >0     | 268  | 4             |

`fr_n=1,2` bodies pick up length-4 via SVML-vectorized DACOS/DCOS.
`fr_n=0` lands at length 2 because of the mixed Real*4 FI / Real*8 P
arithmetic. No `#15344 vector dependence` remarks on any smoother loop.

Coarsest-level thread underutilization (`ngrid_level ∈ {2, 4}` gives
4 or 16 collapsed chunks vs. 128 threads) is negligible because those
levels each do ≤ 4³ = 64 cells of work per call — submicrosecond.

## Build verification

```
module purge
module load intel_comp/2024.2.0
export I_MPI_F90=ifx
module load compiler-rt tbb compiler mpi
make PMP2MGsolver_fR.o   # clean build, no warnings
make PMP2main            # links cleanly
```

### Vectorization (`PMP2MGsolver_fR.optrpt`)

All 12 target inner M1 loops in the three non-smoother routines
vectorize at AVX2 `vector length 4`:

- `calculate_residual_fR` — 6 specializations, all vectorized, 3 with
  tree-reduce for `REDUCTION(+:RES)`.
- `restrict_residual_fR` — 3 specializations (ilevel==1 branch), all
  vectorized. The ilevel>1 branch is untouched (no fr_n dispatch there)
  and legitimately carries a data dependence.
- `calculate_physical_right_hand_side_fR` — 3 specializations, all
  vectorized.

`relaxation_iterations_fR`'s 6 specializations do not vectorize (the
red-black `CYCLE` fragments vector runs); this is expected and the
reason Phase A intentionally excluded the smoother.

## Benchmark status

Benchmark infrastructure was rebuilt during this pass:

- `scripts/bench_sweep.sh` — now snapshots `PMcrs*.DAT` (globbed) along
  with the small IC files, invalidates the snapshot when the live
  `PMcrs0.DAT` is newer, deletes `timing.log` before each run (so a
  stale log can't impersonate a failed run's output), and greps stdout
  for `"Error reading"` / `"Did not get all files"` to catch the
  IC-read failure case where PMP2main exits 0 but did no work.
- `scripts/bench_report.py` — now reads column index 9 (the "MG" column
  from `TimingMain(8,1)` in `PMP2mod_tools.f90:510`) instead of column 4
  (Density, which is the density-deposition bracket and ~1 ms/step on
  this config). Unit label corrected to min/step — the timing.log is
  written as `CPU/60.` so values are minutes, not seconds.

Bench results will be populated after a compute-node run of
`bash scripts/bench_sweep.sh` with the regenerated IC.

## Out of scope

- Loop tiling / cache blocking of the inner stencil.
- Converting `FI2`/`FI3` from `Real*4` to `Real*8` (would affect all
  five solvers, not just fR).
- MPI / halo-exchange changes.
- Multigrid algorithmic changes (cycle shape, smoother count, coarsest
  solver).
- Touching the other four solvers (`csf`, `DGP`, `kmf`, `sym`).
- Refactoring the five solvers into a shared relaxation template.
