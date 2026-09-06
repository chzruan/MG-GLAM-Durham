# Focused regressions for review gaps

Added during the independent review of `a8c7715..cz`. These cover three places
where the merged suites assert a behaviour without exercising it. They add no
production change and rewrite no existing receipt; results are in
[halo_gaps_results.json](halo_gaps_results.json).

```sh
micromamba run -n cosemu python3 -B BDM-refine/repairs/20260906/tests/halo_gaps.py
```

**166 GNU experiments pass** (checked and optimized builds).

## 1. Exact-membership duplicate merging beyond a single pair

`halo_host_ties.py` exercises `MergeNumericalDuplicates` with exactly one
two-member identical set (`exact_set_fallback`). Its bottom-up merge sort, the
representative advance across larger groups, its equal-length lexicographic
comparison and its interaction with a host-removed lowest index were untested.

144 new cases at 1 and 8 threads reuse the existing immutable-priority oracle:

- three- and four-member identical sets, and two groups plus singletons;
- interleaved high candidate indices (7/42/105) where the lowest surviving
  index must win independently of linked-list traversal order;
- equal-length sets differing only in their last identity, which must both
  survive the lexicographic comparison;
- a group whose lowest index is removed by the host pass, so the representative
  must advance to the next survivor;
- a candidate-count sweep over `n` = 1–24, 31–33, 40, 63–65, crossing every
  merge width including the `width > n/2` early exit.

All pass. No defect found; the sort and the grouping are correct.

## 2. Unbinding pass count

The SO contraction caps at sixteen passes and finishes with a sorted interval
scan, precisely so an adversarial profile cannot make it quadratic. The
unbinding loop kept only the trivial `n+1` bound, and neither production nor
the suites record how many passes actually occur.

The `unbinding_ladder` fixtures place antipodal pairs on one shell, so every
survivor shares the discrete shell potential `K*(n-1)` and the survivor bulk
velocity stays zero. Speeds are set so that `0.5*(u + H*R)^2 = K*(n-1) + K`,
which leaves exactly the leading pair marginally unbound at each pass. A
Python replica of the production loop predicts both the pass count and the
final membership; the production routine must reproduce the membership exactly.

| Rows | Unbinding passes | Survivors | Production loop bound |
| ---: | ---: | ---: | ---: |
| 200 | 91 | 20 | 201 |
| 400 | 191 | 20 | 401 |

Passes grow as `n/2`, so the stage is O(n^2) on this profile. This is a
constructed ladder, not a measurement of cosmological incidence; it shows the
bound is real and unmeasured, not that production haloes reach it.

## 3. Outermost-crossing absorption of an equal companion

`two_so_crossings` proves the outermost crossing is selected but does not say
when a companion is absorbed. Two 120-row clumps of radius 0.15 Mpc/h, mass
1e12 Msun/h, Ovdens=200, are scanned in separation. `Rso` for one clump is
1.20257 Mpc/h and for both is 2^(1/3) times that, 1.51514 Mpc/h.

| Separation / single `Rso` | `Rso` (Mpc/h) | Enclosed rows | Companion absorbed |
| ---: | ---: | ---: | :--- |
| 0.6, 0.9, 1.0, 1.1 | 1.51514 | 240 | yes |
| 1.2, 1.25, 1.3, 1.5, 2.0 | 1.20257 | 120 | no |

The measured transition brackets the analytic threshold
`(2^(1/3)*Rso - r_clump)/Rso = 1.135`. The convention is therefore
self-limiting: the selected radius is always the one its own enclosed
population supports, so it can never exceed `(M_enclosed/threshold)^(1/3)`,
and a companion it absorbs necessarily has its centre inside the resulting
aperture, where the distinct-host rule already excludes it. The production
membership equals the independent interval oracle's enclosed population in
every case.
