# SO normalization and catalogue version 3

This follow-up addresses P1 and clarifies P4 in the preserved
[merged review](../../../analysis/review-20260906-merged/REVIEW.md), against
`54f53aae2aac6700d1ddb8490bebb0719f71165d`. That review and its original
experiments are unchanged. Earlier version-2 validation results, including the
N1024 campaign, do not become version-3 validation by this change.

## Density and units

`SetParameters` retains its existing particle-mass constructor, including its
single-precision storage. Let `m` be that actual stored `MassOne`, `L` the
comoving `Box`, and `N = NROW^3` the number of original equal-mass particles in
the full periodic box. The comoving mean matter density is

```
rho_mean_comoving = m * (NROW / L)^3
M_enclosed       = n_enclosed * m
M_enclosed       = (4*pi/3) * Ovdens * rho_mean_comoving * Rso^3.
```

Mass is in solar masses/h, distance in comoving Mpc/h, and the density in
(solar masses/h)/(Mpc/h)^3. `GetHalo` evaluates the right-hand coefficient in
double precision with `pi = acos(-1.d0)`, promoting the stored inputs before
division and multiplication. `NROW` describes the original particle grid;
neither the force mesh `NGRID`, the buffered count `Np`, nor a local gathered
count defines the box mean. A nonpositive `NROW` is rejected on direct entry.

This is algebraic consistency between enclosed mass and background density,
not a claim of bitwise cancellation. The stored/output arrays remain real*4.
The old coefficient `1.150d12 * Om0 * Ovdens` was inconsistent with the
particle mass constructed from approximately `2.774e11 * Om0`. For the same
enclosed count it made the radius about 0.34% larger and the effective mean
overdensity about 1.03% lower. Particle membership can change near the SO
edge, so a uniform rescaling of existing catalogue radii is insufficient.

## All four overdensity modes

The existing flat matter-plus-Lambda background and empirical constants are
preserved. Define
`f = Om0 / (Om0 + (1-Om0)*a^3)` and `x = f-1`, where `a = AEXPN`.
The stored `Ovdens` is always relative to the mean matter density:

| `iVirial` | Stored mean-density threshold | Convention |
|---|---|---|
| 0 | `200/f` | 200 times critical density |
| 1 | `(178+82*x-39*x^2)/f` | Existing virial prescription |
| 2 | `200` | 200 times mean matter density |
| 3 | `(178+82*x-39*x^2)/f * (200/178)` | Existing Abacus-labelled convention |

Physical mean density is `rho_mean_comoving/a^3`; physical critical density
is `rho_mean_comoving/(a^3*f)` in this convention. Thus conversion to a
comoving SO radius introduces no additional power of `a` into the coefficient.
The normalization `178` is preserved, rather than replaced by `18*pi^2`.
No gravitational, Hubble-flow, velocity-unit, unbinding, or shape formula
changes in this follow-up.

## Published catalogue contract

The first header line now ends in `[BDM finder v3]`. There are still exactly
eight ASCII header lines and 24 data columns, with the same ordering, labels,
storage precision and output formatting. Readers that skip eight lines and
load 24 columns remain compatible. A reader that explicitly requires the v2
tag must deliberately accept the changed v3 semantics.

The outermost self-consistent discrete SO crossing within the declared search
cap defines the unextended sphere. Bound membership and `Mbound` use particles
inside that sphere after unbinding. The retained empirical `Rext` correction
defines the reported radius and the enclosed `Mtotal` aperture. Therefore
`Mbound`, `Mtotal`, and the reported radius still are not a single ordinary
SO mass/radius pair. The greater SO mass-radius bound applies to the supported
enclosed population; it does not rule out all asymmetric overlap or bridging.
The review's companion near `1.135 Rso` describes its particular equal-clump
fixture and supplies no cosmological incidence estimate. A companion can
contribute some particles even when its centre lies outside the aperture.

The effective publication floor remains **at least 20 bound particles,
inclusive at exactly 20**, together with the requested `MassMin`. For
mass-consistent rows the writer requires
`Mbound >= max(MassMin, 20*MassOne)`. The old `a8c7715` writer already imposed
this mass floor before its later `aNpart < 10` check. Consequently review P4's
description of an ordinary 10-to-20 published-sample change is incorrect;
the later explicit `<20` membership guard made an existing floor direct.
The new writer additionally checks exact membership count against stored mass.
The independent old/new writer controls below establish this clarification
without changing the historical review or its evidence. A requested mass cut
can dominate the 20-particle floor; the precise boundary uses stored mass and
the writer's floating-point comparisons.

## Targeted validation and fixture adaptation

`normalization_regression.py` compiles actual production `SetParameters`,
`SetOverdensity`, `GetHalo`, gathers, sorting, potential, eigenvalue and writer
routines under GNU and Intel checked/precise optimized modes. Its only
`GetHalo` instrumentation prints the computed double threshold before use;
it does not replace the arithmetic. It also compiles the exact v2 source at
`54f53aa` and the actual `a8c7715` writer as controls. The writer-only floor
fixture fixes concentration to a finite constant in both versions to isolate
the publication cuts.

Fixtures contain all 8^3 or 16^3 original particles, with a local cold clump
and remote matter beyond the local search cap. Their mean density, stored
particle mass, coordinates and SO reference use consistent units. An
independent sorted-interval oracle finds every discrete crossing; the expected
normalization also follows directly from the particle-number volume equation.
All four modes are checked at three matter fractions and three scale factors.
Additional controls cover box length, particle-grid density, force-mesh size,
deliberate stored-mass perturbations, multiple SO crossings, a real v2/v3
edge-membership difference, the unchanged `Rext` aperture, the `NROW` guard,
finite ordered shapes, and eight-header/24-column publication. The artificial
stored-mass perturbations are arithmetic probes after the real constructor;
their transient fixture headers are not cosmological output evidence.

The old and new writers receive 9, 10, 19, 20, 21, 25, 26 and 64-particle
mass-consistent rows. Both publish `[20,21,25,26,64]` at requested floor 0 or
20, `[25,26,64]` at 25, and `[26,64]` at 25.5 particles.

Earlier `halo_*` fixtures which independently chose `MassOne`, `Box` and
`NROW`, or built their radius from `1.150e12`, describe the old normalization.
They must not be used unchanged as v3 SO oracles. The complete-box suite checks
the normalization itself; the executable [v3 algorithm adapter](#current-source-algorithm-regressions)
below retains their algorithm tests against current source. Earlier JSON
receipts remain evidence only for their recorded source. The current peaks
test changes only its expected header suffix to v3; its earlier results are
retained.

The retained receipts record **110 successful cases for each compiler**:
[GNU](normalization-gnu.json) and [Intel](normalization-ifx.json), each split
between checked and precise optimized builds. The
[updated peaks/header controls](peaks-v3-gnu.json) pass 136 GNU cases. Both
compiler suites reproduce the same deliberate edge fixture changing from
257 bound particles under v2 to 256 under v3.

Run from the repository root, always choosing a new receipt path:

```sh
micromamba run -n cosemu python3 -B BDM-refine/repairs/20260906-review-followup/physics/normalization_regression.py --compiler gfortran --output BDM-refine/repairs/20260906-review-followup/physics/normalization-gnu.json
```

Intel requires the compiler and runtime modules loaded before the same command
with `--compiler ifx`; pass the native runtime library path through
`BDM_AUDIT_NATIVE_LIBS` when the Python environment overrides it. No simulations
or Slurm jobs are launched by this suite. Native production replay remains a
separate integration requirement.

## Current-source algorithm regressions

`halo_v3_regression.py` supplies normalization-aware fixtures and Python SO
oracles for the existing `halo_regression`, `halo_followup`, `halo_review`, and
`halo_gaps` tests. Production routines are extracted **verbatim** from the
current `PMP2linker.f90`: no production coefficient replacement or threshold
instrumentation is used. Historical scripts, Fortran fixtures, review notes,
and JSON receipts remain byte-for-byte unchanged. The wrapper records their
hashes, the adapted Python AST hashes, the generated fixture source hashes,
compiler commands, and the current production source hash.

The direct `GetHalo` fixtures provide local excerpts of a hypothetical full
periodic box. The unprovided particles are outside the tested query apertures;
these tests do not exercise snapshot loading or periodic buffer generation.
`Box=32`, `NROW=16`, `Om0=0.3` give the actual stored particle mass of about
`6.6576e11` solar masses/h through the existing real32 constructor. The 100,000
row radial-profile cases use `NROW=128`, preserving a nontrivial crossing
inside their radial range. The independent analytic profile radius and the
sorted-interval SO oracle both use that declared mean density. Shell-fixture
velocities are rescaled by `sqrt(new_mass/old_mass)` to retain meaningful
cold/hot/mixed classifications; all statistics are then checked against the
actual float32 phase data. The separate 110-case suite supplies the complete
small-box check of the global normalization.

The adapter accounts for these cases in checked and precise optimized modes:

| Preserved suite | Cases | Assertions retained or added |
|---|---:|---|
| `halo_gaps` | 166 | All 144 exact-membership cases at 1/8 threads, four energy ladders, and 18 equal-companion scans |
| `halo_followup` | 26 | Outermost crossings, sorted fallback after slow contraction, bound-only statistics, small populations, heap-index width and exact int64 sorting |
| `halo_regression` | 36 | Discrete SO/bin-spacing controls, `Rext`, central/singular/unresolved/capped cases, independent pair energy and unbinding, kinetic/spin/bulk statistics, and parallel candidates |
| Remaining `halo_review` controls | 14 | A single central survivor and all six isolated spherical-potential domain/pair controls in both modes |

The retained [GNU](halo-v3-gnu.json) and [Intel](halo-v3-ifx.json) receipts each
pass all **242 cases** against the combined source with the particle-list and
buffer follow-up integrated. The exact tested source hash and worktree commit
are recorded in each receipt. These are new v3 results, not relabelled historical
receipts.

The ladder reference still gives 91 and 191 passes for 200 and 400 initial
rows, respectively. Production diagnostics must equal those pass counts and
the independent arithmetic-series work totals, 10,010 and 40,110 active
particle rows. The surviving identities must be the final ten antipodal pairs.
The original cold cases retain 64 or 1,000 rows, the mixed/bulk shell cases
retain 900, and the deliberately hot shell loses all rows; the energy controls
cannot pass merely by deleting every halo. The central-survivor probe also
checks that a subsequent empty call clears both diagnostics.

The equal-companion scan remains a constructed fixture, not a survey of
cosmological bridging. Its rescaled absolute radii must not be interpreted as
the isolated production effect of the SO normalization repair: the old fixture
also chose a particle mass inconsistent with its `NROW`/`Box` metadata.

Run the complete current-source suite with a new receipt path:

```sh
micromamba run -n cosemu python3 -B BDM-refine/repairs/20260906-review-followup/physics/halo_v3_regression.py --compiler gfortran --output BDM-refine/repairs/20260906-review-followup/physics/halo-v3-gnu.json
```

For Intel, load the same runtime modules described above and use
`--compiler ifx`. To inspect another Git worktree without editing it, pass
`--source /absolute/worktree/PMP2linker.f90`; all compiled routine extraction,
heap-declaration checks and source hashes then use that file. The wrapper
rejects a source change during the run. The 100,000-row SO controls use sorting
and linear statistics;
they do not execute the quadratic all-pairs oracle. All compiler/run products
are temporary, and only the requested new receipt is written.
