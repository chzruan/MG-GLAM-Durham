# Independent native comparison review

No blocking comparison error was found. The retained test source is
`native_compare_review.py`; `native-compare-review.json` records its exact
command, source hash, results, and the complete earlier probe receipt. This is
a bounded checker review, not a simulation or full native replay result.

The original probe tested analyzer SHA
`e538efe8cf50722821b53fd0d41dd1c0e68f84a403a55c596af3eec1386d97d4`.
The retained script also passed against root commit
`a4cb9deb2406604659fa30c99f65413ce1926fb1`, analyzer SHA
`26c820d43e687e61fe6c42325a3ed915ce66585af9f57f81d71f2c90a7dad561`,
which includes the clarified comparison interpretation. Neither run changed
the reviewed repository, production source, or running replay code.

Four synthetic raw tapes passed the actual `check_memberships` parser. Two
comparisons, with box sizes 512 and 32 Mpc/h, exercised:

- Matching candidate IDs 5 and 9 after their published row numbers shifted.
- One exact and one changed bound set; one unmatched candidate on each side.
- A matched pair above both mass cuts and another below the percentage cut.
- Distinct sentinels for bound/total mass, aperture radius, Vmax, and Rrms.
- A periodic face crossing and bulk-velocity drifts of 50 and 20 km/s.
- Six rejected comparisons with altered raw-tape, configuration, or FI hashes.

All assertions passed. Frozen source hashes matched `build.json`, and
`FindMaxima`, `SetOverdensity`, and `IsDensityMaximum` were identical across
reference, optimized-v2, and v3. Together with the verified common snapshot,
fixed FI, background, configuration, and retained candidate ordering, these
support the candidate-ID join for this campaign. Unmatched published candidate
IDs do not establish that physical objects appeared or disappeared.

Property indices and units are consistent with the raw tape: positions and
reported radii use comoving Mpc/h, velocities use km/s, and masses use
Msun/h. The reported radius includes Rext. Percentage summaries require both
matched masses to be at least `10^12.5 Msun/h` and both quantities to be
positive. Membership overlap is `|A intersection B| / max(|A|, |B|)`; it is
not Jaccard similarity. Exact byte equality is valid for sets because the
upstream checker verifies sorted, unique int64 IDs, and the full analyzer
verifies the stage hashes linking the indexes to their raw tapes.

Two scientific wording limits remain essential:

- Quantizing the seed radius cannot guarantee global reproducibility. Adjacent
  float32 values around 1.25 round to different 0.1 bins in the retained
  counterexample. Ordinary density accumulation can also change which peaks
  pass a threshold or win a local comparison. The measured fixed-FI control
  has a narrower scope than arbitrary thread/compiler reproducibility.
- The equal-clump radius bound is fixture-specific. With SO coefficient 100,
  masses 100 at radius 0.1 and 20 at radius 0.99 give an SO radius 1.06266;
  another 200 at radius 2 cannot establish a larger crossing because even
  all 320 would give radius 1.47361. Thus a companion's nearby particles may
  contribute while its distant centre remains outside the aperture. This is
  a geometric counterexample, not a cosmological incidence measurement.

To repeat the small controls from the repository root, use a new receipt path:

```bash
micromamba run -n cosemu python3 -B BDM-refine/repairs/20260906-review-followup/integration-review/native_compare_review.py --repo "$PWD" --output /tmp/bdm-native-review-repeat.json
```

The recorded command additionally embedded the initial temporary receipt.
Temporary fixture directories are automatically removed. The initial receipt
was deleted only after its hash and complete decoded content were verified
against the retained `initial_probe` entry.
