# Optimize Running Speed of f(R) Model in MG-GLAM

## Context

f(R) multigrid in MG-GLAM dominates per-step wall time. A prior OpenMP pass
(see `OMP_OPTIMIZATION_NOTES.md`) established the baseline we build on:

- `COLLAPSE(3) SCHEDULE(STATIC)` on every triple-nested `!$OMP PARALLEL DO`.
- `-O3 -march=core-avx2 -mfma -fp-model fast=1` in `makefile` (AMD-safe AVX2
  for Cosma8 Zen 2; see existing rationale in `OMP_OPTIMIZATION_NOTES.md`).
- `OMP_NUM_THREADS=128 OMP_PROC_BIND=close OMP_PLACES=cores`.
- Power-spectrum reproducibility budget: 0.1 % (FP reordering allowed).

Per the prior notes, `PMP2MGsolver_fR.optrpt` already has 20 `LOOP WAS
VECTORIZED` hits, so the hot path is partially vectorized. The remaining
wall-time headroom is in the inner `M1` loops of the four hot subroutines,
where:

1. The per-iteration `IF (fr_n .EQ. 0/1/2)` and `IF (ilevel .EQ. 0)` tests
   are **loop-invariant** across a whole call but still issue branches in
   every cell, suppressing some of ifx's SIMD pattern matching.
2. `calculate_residual_fR` on the PM grid writes `FI3` then re-reads it in a
   second triple loop for the L² reduction — an avoidable DRAM pass over
   `NGRID³` entries.
3. `FI2` and `FI3` are **`Real*4`** (declared in `PMP2mod_tools.f90:41-42`),
   so AVX2 gives 8 lanes per vector — the payoff of unlocking any remaining
   unvectorized inner loops is large.

Goal: make the hot `M1` loops more amenable to SIMD by hoisting the
model/level dispatch out of the loop nest, and drop the redundant residual
pass — without touching the multigrid algorithm, array precision, or MPI.

## Structure of the change

```
current                                after Phase B/E
=======                                ===============
relaxation_iterations_fR               relaxation_iterations_fR
└─ DO M3/M2/M1 (COLLAPSE 3)              ├─ IF fr_n==0 / ilevel==0 → DO M3/M2/M1
   └─ CYCLE (red-black)                  ├─ IF fr_n==0 / ilevel>0  → DO M3/M2/M1
   └─ IF fr_n==0  ──┐                    ├─ IF fr_n==1 / ilevel==0 → DO M3/M2/M1
   └─ ELSE fr_n==1 ─┼─ IF ilevel==0/>0   ├─ IF fr_n==1 / ilevel>0  → DO M3/M2/M1
   └─ ELSE fr_n==2 ─┘                    ├─ IF fr_n==2 / ilevel==0 → DO M3/M2/M1
                                         └─ IF fr_n==2 / ilevel>0  → DO M3/M2/M1
                                            (each: !DIR$ IVDEP on M1, red-black CYCLE)

calculate_residual_fR                  calculate_residual_fR
├─ DO M3/M2/M1  writes FI3(M1,M2,M3)    ├─ 6 specialized triple loops, each with
├─ TimingMain(3,1)                      │    !DIR$ IVDEP on M1
└─ IF ilevel==0:                        │    on ilevel==0 branches:
   DO N3/N2/N1  RES += FI3(...)**2      │      REDUCTION(+:RES); RES += OP*OP  ← Phase E
                                        └─ TimingMain(3,1) kept once; reduction loop deleted

restrict_residual_fR / calculate_physical_right_hand_side_fR
└─ same hoist treatment (no red-black, no reduction merge)
```

## Files touched

- `PMP2MGsolver_fR.f90` — all four hot subroutines.
- `makefile` — not changed (flags already optimal).
- `PMP2mod_tools.f90` — not changed (shared module).
- Other solvers (`csf`, `DGP`, `kmf`, `sym`), density/particle code, FFT,
  MPI — not changed.

## Hot-path inventory (verified line numbers in current HEAD)

| Subroutine | Lines | `!$OMP PARALLEL DO` | Role |
|---|---|---|---|
| `relaxation_iterations_fR` | 20–247 | 84 | Red-black Gauss-Seidel smoother (CYCLE at 100–101) |
| `calculate_residual_fR` | 254–446 | 313 (main), 429 (reduction) | Residual + L² reduction |
| `restrict_residual_fR` | 454–686 | 493 (ilevel==1), 663 (ilevel>1) | Fine→coarse restriction |
| `calculate_physical_right_hand_side_fR` | 693–795 | 740 | Coarse-level RHS |

`ilevel>1` branch of `restrict_residual_fR` has no `fr_n` dispatch (plain
8-cell average) — skip Phase B there.

## Phases (ordered; Phase A, B, E land together, Phase C is gated)

### Phase B — Hoist `fr_n` / `ilevel` dispatch above the triple-nested loop

For each of the four subroutines, replace the in-loop
`IF (fr_n .EQ. …) … ELSE IF …` and `IF (ilevel .EQ. 0) … ELSE …` with an
outer dispatch so each specialization contains a clean triple loop with
straight-line cell bodies.

Concrete example for `calculate_residual_fR` — pattern, apply equivalently
to the other three:

```fortran
IF (fr_n == 0) THEN
    IF (ilevel == 0) THEN
        !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC) DEFAULT(SHARED) &
        !$OMP PRIVATE(M1,M2,M3,M1u,M1l,M2u,M2l,M3u,M3l,OP) &
        !$OMP REDUCTION(+:RES)                     ! Phase E (only on ilevel==0)
        DO M3 = 1, ngrid_level
            DO M2 = 1, ngrid_level
            !DIR$ IVDEP
                DO M1 = 1, ngrid_level
                    M1u = M1 + 1; IF (M1u > ngrid_level) M1u = 1
                    ! … (body from existing ilevel==0 branch, no **2 / **3, no branch)
                    FI3(M1, M2, M3) = OP
                    RES = RES + OP*OP             ! Phase E merge
                END DO
            END DO
        END DO
    ELSE
        ! body from ilevel>0 branch, no REDUCTION, no RES update
    END IF
ELSE IF (fr_n == 1) THEN
    ! same split; inline **2 as x*x
ELSE IF (fr_n == 2) THEN
    ! same split; inline **3 as x*x*x
END IF
```

Notes:

- The body is a cut-and-paste of the existing branch contents; no numerics
  change. Verify each specialization against the original by reading the
  diff side-by-side.
- `relaxation_iterations_fR` has three outer cases (`fr_n=0/1/2`) × two
  `ilevel` cases = 6 specialized loops. Keep the two `IF ((redstep) .AND. …)
  CYCLE` lines inside each specialized body for now — Phase C optionally
  rewrites them.
- `restrict_residual_fR` needs the split only inside the `ilevel==1` branch;
  the `ELSE` branch has no `fr_n` dispatch and stays as one loop.
- `calculate_physical_right_hand_side_fR` needs three specializations
  (`fr_n=0/1/2`) — there is only one `ilevel` path.

### Phase A — Enable inner-M1 vectorization where still blocked

After Phase B, inner loops are straight-line and distinct per `(fr_n,
ilevel)`. Add `!DIR$ IVDEP` immediately above each inner `M1` loop in the
**non-smoother** routines (`calculate_residual_fR`,
`restrict_residual_fR`, `calculate_physical_right_hand_side_fR`). IVDEP is
preferred over `!$OMP SIMD` because it preserves the existing `COLLAPSE(3)`
— `!$OMP SIMD` on the inner loop is silently ignored when the outer
`COLLAPSE(3)` flattens the nest.

The smoother (`relaxation_iterations_fR`) is not Phase-A-targeted: the
`CYCLE` at 100–101 fragments vector runs. Phase C (optional) addresses it.

Rationale: the module-level `Real*4` arrays `FI2` and `FI3` are declared in
`PMP2mod_tools.f90` as distinct allocatables — they cannot alias. But ifx
conservatively treats writes to one FI3 section as aliasing reads from
another FI3 section (different offsets into the same array) and falls back
to scalar. `!DIR$ IVDEP` tells ifx these reads and writes are independent
across `M1` iterations, which is true by construction (red-black smoother
excepted).

### Phase E — Merge residual write and L² reduction (ilevel==0 only)

Inside the `ilevel==0` specializations produced by Phase B for
`calculate_residual_fR`, fold the `RES = RES + FI3(N1,N2,N3)**2` loop into
the main triple loop using `RES = RES + OP*OP` with `REDUCTION(+:RES)`.
Delete the second `DO N3/N2/N1` loop and the `IF (ilevel .EQ. 0)` guard.
Keep `RES = DSQRT(RES/(DBLE(ngrid_level))**3)` and `res_PM_grid = RES`.

Correctness: `OP` is exactly the value stored to `FI3(M1,M2,M3)` in the
existing `ilevel==0` branches, so `OP*OP == FI3(N1,N2,N3)**2` bit-for-bit
within one call, and the reduction is already associative-reorderable under
`-fp-model fast=1`. Savings: one full `NGRID³` pass over `FI3` (`Real*4`,
~512 MB at `NGRID=512`) per residual call.

### Phase C — Red-black strided M1 (gated on benchmark)

Optional. Rewrite the inner of `relaxation_iterations_fR` to drop the
`CYCLE`:

```fortran
DO M3 = 1, ngrid_level
    DO M2 = 1, ngrid_level
        M1start = 1 + MOD(M2 + M3 + parity, 2)   ! parity = 0 if redstep else 1
    !DIR$ IVDEP
        DO M1 = M1start, ngrid_level, 2
            ! body, no CYCLE
        END DO
    END DO
END DO
```

Cost: `M1` bounds depend on `M2+M3`, so `COLLAPSE(3)` is no longer legal —
drop to `COLLAPSE(2)` on the outer pair. Still ≥ 256 collapsed chunks for
`ngrid_level ≥ 16`; coarsest levels (`ngrid_level = 4, 8`) give 16 or 64
chunks (below the 128-thread count) but those levels are cheap and not the
bottleneck.

Gate: merge only if the `bench_T128` run improves the relaxation bracket
`TimingMain(8, …)` by ≥ 15 %. The previous optimization pass explicitly
chose to keep `CYCLE + COLLAPSE(3)` until a benchmark justified otherwise.

## Edit order and what to run between steps

1. Phase B + A + E on `calculate_residual_fR`. Build; inspect optrpt.
2. Phase B + A on `restrict_residual_fR` (ilevel==1 branch only; leave
   `ilevel>1` branch alone). Build; inspect optrpt.
3. Phase B + A on `calculate_physical_right_hand_side_fR`. Build; inspect
   optrpt.
4. Phase B on `relaxation_iterations_fR`, `CYCLE` preserved. Build;
   inspect optrpt.
5. Benchmark. If relaxation time is still the bottleneck, do Phase C on
   `relaxation_iterations_fR` as a separate commit, benchmark again, keep
   or revert.

Each step is a separate commit so regressions are bisectable.

## Verification

**Build** (module set is non-persistent; re-source per shell):

```bash
module purge
module load intel_comp/2024.2.0
export I_MPI_F90=ifx
module load compiler-rt tbb compiler mpi
make PMP2mod_MGbackground.o   # module-ordering quirk (see OMP_OPTIMIZATION_NOTES.md)
make PMP2main
```

**Vectorization check** — after each build step:

```
grep -cE "LOOP WAS VECTORIZED" PMP2MGsolver_fR.optrpt
grep -E "vector dependence" PMP2MGsolver_fR.optrpt
```

Expected: `LOOP WAS VECTORIZED` count rises above the baseline of 20
recorded in `OMP_OPTIMIZATION_NOTES.md`. Any `vector dependence` hits on
`M1` loops in the three non-smoother routines after Phase A should be
investigated (likely a missing `PRIVATE` variable on the directive).

**Correctness** — re-run user's fid_F5 reference P(k) and confirm
`max |P_new/P_ref − 1| < 1e-3` across the resolved k-range. The 0.1 % bar
is already accepted. Optional sanity: run at `OMP_NUM_THREADS ∈ {1, 128}`
and compare the printed `res_PM_grid` across threads to the 0.1 % bar.

**Benchmark target** — ≥ 20 % wall-time reduction on a 5-timestep f(R) run
at `OMP_NUM_THREADS=128` relative to the current `HEAD`, read from the
existing `CALL TimingMain(8, …)` brackets. Phase C adds upside if gated in.

## Out of scope

- Loop tiling / cache blocking.
- Converting `FI2` / `FI3` precision (any change would affect all five
  solvers).
- MPI / halo-exchange changes.
- Multigrid algorithm changes (cycle shape, smoother count, coarsest
  solver).
- Touching the other four solvers (`csf`, `DGP`, `kmf`, `sym`).
- Refactoring the five solvers into a shared template.
