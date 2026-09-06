# Repairs following the independent halo review

The two findings recorded in [halo_review.md](halo_review.md) are repaired in
the production routines. Original finding evidence in
`halo_review_results.json` remains unchanged.

The SO calculation gathers the **complete supported physical search cap**.
For n particles inside that cap, every valid SO root has radius no larger than
`(n*MassOne/threshold)^(1/3)`. Removing particles outside this upper bound
cannot remove a particle enclosed by any valid root. Repeating the operation
therefore converges from above to the largest self-consistent population.
If a pass would reduce the population below ten, no supported root remains.

Contraction is limited to sixteen passes; slower cases finish through a
sorted interval scan of the remaining particles. Thus an adversarial profile
that removes one particle per pass cannot make the SO stage quadratic. Most
diffuse cap-neighbourhood particles are discarded before any sorting. The
Rext aperture is still gathered and sorted separately for iterative binding.
The candidate is explicitly rejected if the search cap remains overdense.

Both radius and ID heapsorts now use **int64 local sizes and indices**. The
large-index regression compiles the actual production index declaration and
checks a valid 1,200,000,000-row heap's next child at 2,200,000,000, without
allocating that array. This child exceeds the heap and is skipped; it cannot
wrap to a negative subscript.

Run:

```sh
micromamba run -n cosemu python3 -B BDM-refine/repairs/20260906/tests/halo_followup.py
```

Twenty-six checked/optimized GNU probes pass:

- The independent two-shell interval oracle gives Rso=1.449183548 Mpc/h and
  210 enclosed/bound rows at each of Cell=0.5, 1, 2 and 4 Mpc/h.
- A deliberately slow contraction profile exercises the sorted fallback and
  returns the correct ten-row outermost root.
- Hot contaminants leave a rotating, drifting cold population. Bulk velocity,
  kinetic energy, shell pair energy, RMS radius and spin agree with separate
  compensated/pairwise calculations using only the final survivors.
- Zero/one/nine/ten initial rows return their expected populations; a later
  empty call clears previous membership and mass. Production identity maps
  use actual source-row indices.
- Both production heap index types preserve large arithmetic. The standalone
  ID sorter also sorts values exceeding 2^31 correctly.

Full source/test hashes, commands and output are in
[halo_followup_results.json](halo_followup_results.json). The broader existing halo suite completed its checked cases and printed all
optimized shell/control successes, but was stopped during the remaining Python
oracle work before an aggregate result was produced; no full-suite pass is
claimed here. The integrated full suite and native finder validation are
coordinated on the integration branch; no new performance claim is made here.
