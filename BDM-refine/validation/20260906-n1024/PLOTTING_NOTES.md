# Same-snapshot BDM property comparison

`compare_properties.py` compares the pre-audit finder with the repaired finder
on identical saved snapshots. The initial revisions were `a8c7715` and
`b6e96669`; validation followups can use later repaired revisions. Revision
labels come from metadata, not hard-coded source identities: `old_finder_commit`
and `new_finder_commit` take precedence; `source_commit` is the new-side fallback.
Its input is a consolidated NPZ containing
`old_z0`, `new_z0`, `old_z1`, `new_z1`, `old_z2`, `new_z2`: each is an `(N,24)`
native Catshort table. An empty population is `(0,24)`. A genuinely available
subset of epochs is supported; never copy an epoch's rows into another epoch.
The adjacent metadata JSON must contain `box_mpc_h` and `nrow`; main-preparation.json
is accepted directly, and its full contents/hash are preserved in the summary.
Keep actual epoch, snapshot, finder-source and binary hashes in that metadata.

Run through the required environment, using a fresh output directory:

```sh
MPLCONFIGDIR=/tmp/bdm-matplotlib-config micromamba run -n cosemu python3 -B \
  BDM-refine/validation/20260906-n1024/compare_properties.py \
  --catalogues main-comparison-catalogues.npz --metadata main-preparation.json \
  --outdir property-comparison
```

The script uses float64 diagnostic arithmetic, a one-worker periodic cKDTree,
one-core BLAS, and Matplotlib Agg with the copied publication house style.
LaTeX is the default; `--no-tex` explicitly selects mathtext if required.
Full simulation inputs should run in their allocation; only the archived N128
mechanics checks below were run on the login host.

`--skip-plots` separates matching/statistics from rendering. Replot without
repeating data construction with:

```sh
MPLCONFIGDIR=/tmp/bdm-matplotlib-config micromamba run -n cosemu python3 -B \
  BDM-refine/validation/20260906-n1024/compare_properties.py \
  --plot-ready property-comparison/comparison-plot-ready.npz \
  --summary property-comparison/comparison-summary.json --outdir replotted
```

Outputs are seven vector PDFs, compressed plot-ready arrays, a scientific summary
with input/provenance hashes, and a figure receipt with output hashes:

- `hmf.pdf`: full published populations, plus fractional abundance differences.
- `matched_primary.pdf`: catalogued mass, aperture radius, Vmax, concentration,
  spin and bulk-velocity changes.
- `matched_structure.pdf`: both axis ratios, centre offset, radial RMS, velocity
  RMS and direct Rmax if available.
- `matched_total_mass.pdf`: extended-aperture mass changes.
- `matched_energy_direction.pdf`: signed virial-ratio changes and sign-independent
  major-axis angles, with conservative direction conditioning.
- `matching.pdf`: mass-dependent match coverage, normalized separation CDF, and
  exhaustive matched/unmatched/ambiguity counts.
- `quality.pdf`: invalid-value and zero-sentinel fractions in published rows.

The figure headers state the actual particle count and box. Inputs with
`nrow < 1024` are labelled **Small integration data** by default; `sample_label`
can provide a more specific provenance label. These figures never establish
precision or resolution convergence.

Optional metadata `reference_logmass: 12.5` adds a dotted mass marker and light
grey shading below it on every mass-axis figure, including HMF and match coverage.
This is a **z=0 literature reference** supplied by the producer;
the script neither derives nor validates it. Preserve the associated citation
in metadata (for example `reference_source`) for the deck. The marker applies
no additional selection and does not establish convergence of the revised
finder. Its line at all redshifts is a common visual reference, not evidence of
z=1 or z=2 convergence. A mass range outside the measured bins can be shown solely to make the
reference visible; no data are extrapolated into that range.

Plot-ready schema 3 includes energy and direction fields with a direction cut
that is independent of transverse-axis label ordering. Schema-1 statistics lack
the energy/direction fields; schema 2 used b/a alone and can retain near-degenerate
historical objects with reversed labels. Regenerate earlier statistics from the
source catalogue NPZ into a fresh output directory before using replot mode;
historical receipts stay unchanged.

## Units and definitions

The columns follow `PMP2linker.f90::WriteFiles`, in zero-based indexing:
0–2 position [comoving Mpc/h]; 3–5 bulk velocity [km/s]; 6 catalogued/bound mass
and 7 total aperture mass [Msun/h]; 8 reported radius [comoving kpc/h]; 9 Vrms
and 10 Vmax [km/s]; 11 catalogue ID; 12 concentration; 13 particle count;
14 distinct/sub flag; 15 Xoff; 16 2K/Ep−1; 17 spin; 18 radial RMS [comoving
kpc/h]; 19–20 b/a,c/a; 21–23 major-axis direction.

**Rmax is absent from Catshort.** To compare it, supply both optional
`old_rmax_z0`/`new_rmax_z0` (and corresponding epoch) vectors, aligned exactly
with catalogue rows, and set metadata `rmax_units` to `comoving_kpc_h`.
Raw internal membership-tap Rmax is Mpc/h and must be converted by the producer.
The script requires that unit declaration. It never reconstructs Rmax from
concentration: the normal concentration branch uses a mass/radius/Vmax inversion,
so the fallback relation cannot be applied to every halo.

The repaired mass and bulk/energy/spin/shape/RMS properties use converged bound
membership. The old implementation mixed populations and normalizations. The
same numeric column names therefore do not imply identical definitions. Radius
is the extended aperture with the Rext correction, not the unextended SO crossing.
Xoff is dimensionless, normalized to the reported radius. Concentration is a
diagnostic with a fallback, not a fit to an independently measured density profile.

Native ASCII serialization limits the comparison: positions have four decimal
places, bulk velocities two, masses four significant digits and radii five.
Small differences can be hidden by that rounding. A zero catalogue residual
does not establish bitwise equality of the internal property or particle set.

## Selection and matching

HMF means dn/dlog10(Mcat) per comoving volume. It includes **every published row
with a finite positive mass**, regardless of match status or distinct/sub flag.
The summary preserves flag counts. It applies no extra mass cut, avoiding a
second threshold on rounded ASCII masses. The actual catalogue selection is
the finder-defined `max(MassMin,20*MassOne)` (MassMin=2.5e12, iVirial=1,
Rext=.15 for the intended large replay). The two finders apply selection to
their own mass definitions; new/old counts can differ because of membership,
property definitions, duplicate handling, or crossing that selection boundary.

Positions are canonicalized periodically. Matching requires finite positions,
finite positive reported radii, mutual nearest neighbours, and
`d <= 0.25 * min(Rold,Rnew)` after converting the radii from kpc/h to Mpc/h.
On **both** sides the nearest alternative must be more than twice as distant
and must differ by more than 1e-4 Mpc/h. This ambiguity veto excludes exact
duplicates and ties comparable to the ASCII coordinate precision. These choices
are configurable (`--match-fraction`, `--ambiguity-ratio`, `--tie-atol`) and
saved in the summary; comparisons across choices should use fresh output dirs.
IDs are never matched. Geometric agreement is not proof of shared particles.

Unmatched categories are exclusive in this order: invalid matching input,
no eligible target, non-mutual nearest, outside radius cutoff, ambiguous
neighbour. Thus ambiguous counts refer to otherwise acceptable mutual matches;
non-mutual candidates can also lie in a crowded region. Their sum plus matched
rows equals the original catalogue count, separately for each finder.

Matched property statistics are conditioned on both catalogues publishing an
object, this conservative geometric selection, and validity of that property.
The plot-ready NPZ retains input row indices, matching statuses, paired values,
changes, bin counts and percentiles. Unmatched objects remain in the HMF.
Trends use old catalogued mass on the horizontal axis. Percent changes exclude
zero old denominators. Vmax/concentration/Rmax zero sentinels are also excluded
on the new side. Spin, axis ratios and Xoff use **absolute differences** to
avoid unstable fractional changes near zero; bulk velocity uses the norm of
the vector difference in km/s. Component-wise bulk differences are retained.

Virial ratio is the **signed** column-16 quantity `2K/Ep-1`. Its comparison is
absolute new minus old, keeping finite negative and zero values without division
by the old value. Since energy normalizations changed, this comparison alone is
not a common-definition test of equilibrium.

Major-axis directions use columns 21–23. Each finite nonzero vector is normalized,
then the plotted angle is `acos(abs(dot(old_hat,new_hat)))` in degrees (0–90).
Taking the absolute dot product removes the arbitrary eigenvector sign; neither
non-unit input norms nor sign flips imply a directional change. Zero or nonfinite
vectors are invalid. To avoid directions that may be unstable near axis degeneracy,
the plotted trend additionally requires `max(b/a,c/a) < 0.9` in **both** catalogues
by default (`--direction-ba-max` retains its CLI name and changes this threshold).
Both reported ratios must be finite and pass the same axis-ratio domain as the
property diagnostics, [0,1] with the 5e-4 upper ASCII allowance. Taking the larger
ratio makes the cut independent of historical transverse-axis label reversals;
it never rewrites either reported ratio. This conservative proxy also excludes
oblate cases; raw eigenvalue gaps are unavailable and this is not a calibrated
orientation-reliability cut. The summary and plot-ready NPZ
retain invalid-vector counts, shape-conditioning exclusions, the condition mask,
and unconditioned angles, separately from the retained direction sample.

HMF bins with 0<N<20 have open markers; zero counts are retained in evidence but
cannot be drawn on log axes. HMF ratio bins require both counts >=20. Marginal
sqrt(N) bars provide a counting scale and are not an uncertainty estimate for
the correlated difference on the same snapshot. Matched trends require >=20
valid pairs per bin; bands show 16th–84th percentiles, not errors on the median.
Low-count choices are configurable, and all unmasked counts are retained.

Invalid domains are nonfinite/nonpositive masses or radii, negative speed/RMS,
spin/offset/concentration values, or axis ratios outside [0,1] (allowing 5e-4
ASCII rounding). A separate `axis_order` diagnostic counts c/a>b/a+5e-4. The
retained independent empirical corrections can invert the order of reported
axis ratios, even after correctly ordered eigenvalues. In-range calibrated
ratios remain in matched trends; the ordering anomaly is separately visible in
the quality PDF and JSON. Zero Vmax, concentration and auxiliary Rmax are separately
reported as unresolved sentinels; physical zeros in spin/offset/shape are not
automatically invalid. Quality denominators contain published rows only:
candidate rejection rates, SO truncation and singular-centre statuses must be
read from the finder logs. The catalogue alone cannot establish their absence.

## Bounded mechanics validation

`--self-check` checks against a tiny independent dense periodic-distance oracle,
boundary wrapping, duplicate-centre ambiguity, empty counterparts, invalid
positions, minimum-radius cutoff/unit conversion, final mass-bin inclusion,
zero-sentinel handling, empty full analysis, direct auxiliary Rmax and its
required units, and the independent repaired-axis-order assertion. Additional
controls cover signed virial-ratio subtraction through zero, normalized sign-flip
and orthogonal directions, zero/nonfinite vectors, extreme finite vector norms,
both-side direction conditioning, and reference-marker selection independence.
Further controls reject reversed near-degenerate ratios on either catalogue side
and nonfinite/negative/out-of-range transverse ratios, while checking that raw
ratios remain unchanged. These artificial values are tests, not plot inputs.

The real-data mechanics fixture uses the archived N128, L128 z=0 catalogues:
`analysis/full-audit-20260906/catalogues-n128.npz:baseline` (757 old rows) and
`repairs/20260906/native-n128-catalogues.npz:baseline_0_t1` (1215 repaired rows).
The corresponding historical receipts identify the fixed saved snapshot and
source hashes. This fixture has one real epoch, and is labelled small integration
data. The pending N1024 results must be plotted from their own consolidated NPZ.

During this plotting validation, the archived N128 sample exposed c/a>b/a in
2/757 old and 7/1215 repaired rows, beyond the 5e-4 rounding allowance. The parent
audit assigned a separate repair to order the two empirically corrected values
without changing their coefficients. These archived results remain unchanged.
For the final main replay, pass `--expect-ordered-new` to independently require
zero ordering violations in the actual published repaired catalogues. The
default keeps historical anomalies visible rather than rewriting either axis.
