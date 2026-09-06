# Preserve original-particle precision across periodic images

The review reproduced an additional rounding error introduced by storing a
positive periodic image in float32. For a box of 32 Mpc/h, the already-stored
float32 particle coordinate `0.1` and centre `31.9` have an exact double
minimum-image separation of `0.20000038295984268`. A radius of exactly 0.2 must
exclude that particle. The old stored image instead gave separation
`0.19999885559082031`, and the actual `AddBuffer` + `BdmHaloGather` path included
it. This is additional image arithmetic, independent of the input snapshot's
existing precision. Source hashes and the observed output are preserved in
`halo_particle_review.json`; `halo_periodic_cutoff.f90` is the small fixture.

`BdmParticleCoordinate(row, component)` now reconstructs a coordinate as the
original float32 row promoted to double plus its integer image shift times the
box size. The shift is inferred from the stored ghost; its small rounding error
cannot change the integer. `BdmParticlePosition` exposes the three coordinates.
Explicit image shifts remain distinct: minimum-imaging every ghost separately
would let multiple copies represent the same physical particle in one query.

The particle list, halo gather, unbinding/statistic offsets and centring all
use these reconstructed coordinates. List bins and query limits now evaluate
division and centre-plus/minus-radius in double precision, using the same
ceiling convention. This preserves the monotonic broad-phase relationship:
a coordinate inside the query interval cannot receive a cell outside its
bounds. GetHalo rounds its real32 broad-phase radius upward before applying
the exact real64 spherical cutoff. The existing per-iteration centring wrap
into `[0,Box)` is preserved and explicitly tested for a negative mean across a
face and a corner.

This keeps the original analysis rows in their existing float32 storage and
removes extra ghost-coordinate rounding. It does not recreate precision absent
from the stored original analysis coordinates. The integer identity map refers
to original array rows; each original row maps to itself. Test-only arbitrary
ID permutations are therefore inappropriate for that map and were replaced
with valid original/ghost mappings.

`halo_periodic_regression.py` runs checked and optimized GNU or ifx builds at
1, 2 and 4 threads. Each compiler passes 30 cases covering independent
minimum-image membership oracles, exact-radius face/corner/list boundaries,
14,336 original rows exercising the parallel list path, near-half-box searches
without double counting, canonical centring, and exact agreement of halo
properties for dyadic translations between a face and a corner. The retained
PM arrays and their restoration were not edited; all 104 particle regressions
still pass. The original halo-property regressions are rerun with the shared
helpers and the valid original-ID fixture.

```
micromamba run -n cosemu python3 -B BDM-refine/repairs/20260906/tests/halo_periodic_regression.py \
  --compiler gfortran --output BDM-refine/repairs/20260906/tests/halo_periodic_gnu_results.json
```

For ifx, initialize Intel's compiler/runtime modules and export
`BDM_AUDIT_NATIVE_LIBS="$LD_LIBRARY_PATH"` before entering `cosemu`; then use
`--compiler ifx`. The test driver keeps Python in `cosemu` and gives native
compiler/finder subprocesses that recorded Intel runtime path.
