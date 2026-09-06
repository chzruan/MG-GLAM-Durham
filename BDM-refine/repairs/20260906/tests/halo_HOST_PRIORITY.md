# Deterministic exclusion of equal-mass host candidates

The first repaired 128³ native replay exposed three close host pairs with equal
bound masses whose particle sets differed only at their edges. Exact-set
comparison correctly found them unequal, while the legacy host test required a
strictly larger mass and therefore kept both centres. The preserved native
fixture contains the original candidate indices, all properties and full bound
IDs, with the source-archive hash in `halo_equal_mass_native.json`.

| Candidate pair | Bound particles, each | Shared particles | Representative |
|---|---:|---:|---:|
| 261, 263 | 2,110 | 2,109 | 261 |
| 673, 695 | 329 | 326 | 673 |
| 1086, 1087 | 1,445 | 1,441 | 1086 |

`RemoveDuplicates` now orders candidate priority by larger stored bound mass,
then lower stable candidate index for an exact mass tie. A lower-priority
candidate is suppressed when its centre lies **strictly inside** the
higher-priority candidate's reported radius. Equality at the radius is retained.
The measurements remain immutable during the parallel pass; only afterward is
the combined mask applied. Thus a candidate may be suppressed by a host that is
itself suppressed, preserving the pre-existing host-chain convention.

This is a geometric distinct-host rule. It adds no fuzzy particle-overlap or
fixed-distance identity threshold. Disjoint low-mass hosts and hosts with only
overlapping outskirts remain separate when their centres lie outside the
higher-priority host radius. Exact bound-particle-set comparison remains the
final fallback for identical objects. The empirical Rext aperture remains the
radius used by the existing host exclusion contract.

For equal masses with different radii, priority remains the stable candidate
index. Consequently, if the lower-priority centre lies outside the
higher-priority radius, both survive this host test even when the reverse
containment holds. Swapping which radius belongs to the higher-priority index
can change that result; the tests make this convention explicit.

Host distance arithmetic and periodic shifted query bounds now use float64.
Two additional numerical regressions establish why both changes are needed:

- With a radius of `0.200000599026680`, the stored float32 positions 0.1 and
  31.9 in a box of 32 have true periodic distance `0.200000382959843`.
  Float32 subtraction can instead give `0.200000762939453`, rejecting genuine
  containment.
- A host at x=28 with radius `4.100000858306885` contains a candidate at
  x=`0.10000038892030716` across the same boundary. The float32 shifted query
  becomes `32.10000228881836`; its lower bound moves to
  `28.000001430511475` and misses the host's boundary cell for Cell=3.5.
  Direct double query bounds agree with `ListMaxima`'s double ceiling bins.

The dedicated `halo_host_ties.py` harness compiles the actual before/after
routines. Both GNU and ifx pass **160 repaired cases each**: 20 scenarios,
checked and optimized arithmetic, at 1, 2, 4 and 8 threads. Independent all-pairs
priority/containment and exact-set oracles cover the native pairs, periodic
translations, strict cutoffs, unequal-radius ties, immutable chains, disjoint
hosts, overlapping outskirts, sparse stable indices and clustered candidates.
The historical routine reproduces the equal-mass, narrow-phase precision and
shifted-query failures. Full native finder replays are performed separately
on the integration branch.

```
micromamba run -n cosemu python3 -B BDM-refine/repairs/20260906/tests/halo_host_ties.py \
  --compiler gfortran --output BDM-refine/repairs/20260906/tests/halo_host_ties_gnu_results.json
```

The ifx mode uses the recorded Intel runtime through `BDM_AUDIT_NATIVE_LIBS`;
Python remains in `cosemu`. Temporary binaries and input files are removed
on completion. Existing core test sources and historical audit artifacts are
unchanged by this repair.
