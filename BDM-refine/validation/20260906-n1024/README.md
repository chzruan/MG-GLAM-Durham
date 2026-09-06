# BDM production-scale validation and property comparison

This experiment follows the completed [audit repairs](../../repairs/20260906/README.md)
merged into `cz` at `b6e96669`. It evolves a new GR realization with **1024³
particles, a 2048³ PM mesh and a 512 h⁻¹ Mpc periodic box**, and compares the
pre-audit (`a8c7715`) and refined finders on the **same saved particle snapshots**.
The experiment is developed on `validation/bdm-n1024-ng2048-20260906`.

The scientific comparison figures and [completed validation receipts](validation-summary.json)
are tracked here. The presentation PDF, LaTeX source, theme assets and generated
result macros are local, ignored files in `slides/`; they are excluded from the
current Git tree. Their original build evidence remains in
[presentation-validation.json](presentation-validation.json).

## Completed results

| Redshift | Pre-audit haloes | Refined haloes | Change | Matched pairs |
| --- | ---: | ---: | ---: | ---: |
| 2 | 49,205 | 52,370 | +6.43% | 48,998 |
| 1 | 114,454 | 119,688 | +4.57% | 113,641 |
| 0 | 141,357 | 147,042 | +4.02% | 140,598 |

The three **refined membership-checked catalogues** have zero identical bound
sets, zero repeated particle IDs within a halo, zero mass/count mismatches and
zero distinct-host exclusion violations. Their 24 published fields are finite
and their transverse axis ratios are ordered correctly. This does not turn
the pre-audit catalogues into membership-validated data.

For matched z=0 haloes with **both catalogued masses at least 10¹²·⁵ h⁻¹ M⊙**,
median changes are **+2.184% in bound mass, +2.241% in reported radius and
+0.068% in Vmax**. The mass/radius cohort contains 118,939 pairs; the Vmax
cohort contains 118,937 after excluding unresolved values. Other medians are
−2.033% in the concentration diagnostic, −0.002213 in spin, +0.0048 in b/a,
+0.0056 in c/a and −0.00532 in the dimensionless centre offset. Full
16th/50th/84th percentiles and property-specific sample sizes are retained in
`validation-summary.json`; the bands in the plots describe halo-to-halo
scatter, not the uncertainty of the median.

Seven standalone vector figures are available:

- [Mass functions and abundance ratios](slides/figs/hmf.pdf).
- [Matched bound mass, radius and Vmax](slides/figs/matched_primary.pdf).
- [Concentration, spin, shape and centre offset](slides/figs/matched_structure.pdf).
- [RMS radius/velocity, bulk velocity, energy and major-axis direction](slides/figs/matched_energy_direction.pdf).
- [Matched total aperture mass](slides/figs/matched_total_mass.pdf).
- [Matching coverage and exclusions](slides/figs/matching.pdf).
- [Invalid, unresolved and axis-order diagnostics](slides/figs/quality.pdf).

## Experiment

The main run starts at z=100 and saves z=2, 1 and 0, using the repository GLAM
initial-condition generator, realization 1, Ωm=0.3089, ΩΛ=0.6911 and σ8=0.8159.
The nominal PM cell size is 0.25 h⁻¹ Mpc and particle mass is approximately
1.07×10¹⁰ h⁻¹ M⊙. The native finder uses its existing density normalization
2.774×10¹¹ rather than the rounded 2.775×10¹¹ used in planning estimates;
membership validation always uses the actual stored particle mass.

Both finders request virial overdensity, `MassMin=2.5e12` h⁻¹ M⊙ and `Rext=0.15`.
The old first-call peak/configuration-order defect is intentionally retained in
the baseline. Each historical replay starts a fresh process, so these results
describe that defined standalone use; they are not a universal correction for
historical multi-output inline runs, whose configuration could carry between calls.
The simulator uses the refined finder inline. Replaying the identical saved
particles with the old finder isolates catalogue changes from changed dynamics.

## What this resolution tests

The large run tests duplicate membership, distinct-host exclusion, large-array
indexing and memory use, published-field validity, repeated/parallel finder
agreement and preservation of the simulation particle state. Matched-object and
population plots quantify the complete refinement's effect, including changed
property definitions. They do not assign each change to one individual repair.

This dimensional configuration was studied in [Ruan et al., MG-GLAM,
arXiv:2110.00328, Appendix A](https://arxiv.org/pdf/2110.00328), including force
and mass resolution tests. Its published halo-abundance resolution near
10¹²·⁵ h⁻¹ M⊙ at z=0 motivates a labelled reference mass in the figures; it is not
new convergence evidence for the revised finder. One GR realization cannot
establish precision mass-function convergence, inner-profile/concentration/shape
accuracy or model-specific modified-gravity binding. The catalogue cut is below
that reference mass, so lower-mass results are shown with that distinction.

The refined catalogue uses converged bound members for bound mass and internal
properties. Its reported radius and total mass retain the empirical aperture
extension, so they are not an unmodified standard SO mass/radius pair. Historical
catalogues contain different one-pass/property-population conventions. An
external analytic SO mass function is consequently not used as a correctness
oracle for this comparison.

## Resources and provenance

The authorized particle limit is strictly below 1200³. A 512³/1024³ resource
pilot in the same box precedes the main allocation. The pilot changes both
particle and force resolution and is used for sizing and diagnostics; it is
not treated as a matched-phase convergence pair.

Jobs use the shared `cosma8-serial` queue, actual core requests and explicit
memory. `resource-plan.json` records earlier pilot accounting and predictions;
`jobs.json` records submission commands, immutable submitted script hashes,
job IDs and dependencies. Main sizing follows inspection of the completed
pilot's ReqTRES, AllocTRES, billing, elapsed/CPU time and MaxRSS. No exclusive
128-core node is requested without a measured reason.

`preparation.json` freezes the source and executable hashes and reference-input
archive hashes. Only the needed archived binaries and three small IC input files
are restored into `work/`. Historical audit/repair archives and receipts remain
unchanged. `dimension-preflight.json` records the read-only large-index review;
its annotated `cell_Mpch` field refers to the BDM linked-list scale L/NROW,
not the PM mesh scale L/NGRID.

Per-process logs and atomic JSON stage receipts live in `work/`. Completed
stages can be reused only with matching binary, input selection and log hashes;
snapshot bytes are verified independently. Failed/interrupted stages are
retained for inspection. The driver refuses to silently restart an incomplete
Fortran publication or overwrite unexplained outputs. Periodic overwriting
checkpoints are disabled; the small number of IC and named snapshot files are
kept for recovery and replay. Comparison arrays and numerical summaries are
consolidated to limit inode use.

The main simulation job **11948372** completed in 3,477 s (58.0 min) on
64 shared cores, with 166.6 GiB batch MaxRSS against 192 GiB requested. The
membership/ordinary-replay job **11948467**, fixed-density control
**11948491** and final plotting job **11948494** also completed. All eight
simulation, pilot, diagnostic and plotting submissions used **110.887 billed
core-hours**, including the preserved failed pilots; the superseded plot job
was cancelled while pending and consumed none. Full ReqTRES, AllocTRES,
billing, elapsed time, TotalCPU, MaxRSS, utilization and time-limit estimates
are recorded in [accounting.json](accounting.json). Cleanup accounting is
recorded separately in `cleanup.json`.

At z=0, standard standalone replay took 45.97 s before the audit and 91.83 s
after refinement at 64 threads: **2.00 times the elapsed time**, including
snapshot reading, density calculation, finder and catalogue output. The
refined algorithm does additional SO search and iterative unbinding. In the
separate fixed-density experiment, two finder-only repeats per thread count
give median times 177.01 s at 32 threads and 108.97 s at 64 threads:
**1.62 times speedup, 81% of ideal doubling**. That instrumented experiment
holds extra density/particle arrays, so its times are not substituted for the
standard replay measurement. Future full plotting jobs should request 3 GiB
for headroom above the measured 1.99 GiB batch peak.

## Validation and plotting procedure

`run_validation.py` generates ICs and evolves the simulation under Slurm, checks
native headers, file lengths and final epochs, and saves snapshot SHA256 hashes.
Refined membership taps are read one halo at a time, keeping memory proportional
to the number of haloes plus the largest halo. Checks require strictly unique
original IDs within each halo, consistent mass/count, no identical surviving
member sets, and no lower-priority centre inside its host's radius, using the
same documented mass/index priority and periodic geometry.

Ordinary instrumented, standard and state-probe replays record byte agreement
and any row differences from the inline catalogue. Parallel float32 CIC density
accumulation is not bitwise reproducible: the N512 pilot isolated one-ulp changes
in two density-dependent centring apertures, affecting 2/71,587 published rows
while preserving every bound mass, particle count and bulk velocity. Holding
the exact density field fixed gives byte-identical finder catalogues across
16/32/16/32 threads; replacing only that field reproduces both original rows.
The original strict pilot failure remains in its receipt; the controlled
diagnosis is separate evidence, not a retroactive pass of that test.

The state probes check all six particle arrays bitwise after each finder call,
their allocation sizes and released per-call buffers. The main numerical test
therefore separates ordinary density sensitivity from the immutable-density
finder-thread comparison. Old/new standard z=0 timings use equal thread counts;
membership-dump or state-probe timings are not substituted for normal timings.

At 1024³, all six fixed-density calls and the two-call ordinary state probe
preserve the six particle arrays bitwise, their sizes and counts. Each tested
fixed field gives byte-identical finder catalogues at 32 and 64 threads.
Changing only the density field reproduces the ordinary 32/64-thread
difference: **7 of 147,042 printed rows**, including one bound-particle count
changing from 2,564 to 2,566 (about 0.080% in printed bound mass). The maximum
coordinate-component change is 0.7 h⁻¹ kpc. All 794,570 candidate indices and
positions agree between fields; density values and seed radii can differ.
These are row comparisons across different fields, not evidence for identical
particle membership in the two calls. The ordinary membership/standard
64-thread replay differs from the inline z=0 catalogue in one row while
preserving every bound mass, count and bulk velocity; z=1 and z=2 agree exactly.
See the [main control report](main-thread-control/README.md) for the isolation
experiment and the preserved [pilot diagnosis](numerical-threading/README.md).

The plots use the exact standalone refined catalogues covered by the raw
membership checks, paired with the pre-audit standalone catalogue on the same
particles. `main-verified-comparison-catalogues.npz` and `comparison-metadata.json`
bind every plotted refined table to its catalogue and raw-membership hashes.
Inline catalogues and their ordinary replay differences remain separate
evidence; matching printed rows does not transfer a raw-membership claim to
another finder call.

The plot mechanics check also exposed a retained empirical shape correction
that could invert the two transverse axis ratios. The additional
[shape repair](shape-repair/README.md) orders the corrected values without
changing their coefficients, principal direction, other properties or
memberships. It passed 574 cases under each of GNU and Intel, including full
production GetHalo fixtures. `main-preparation.json` and `native-build.json`
freeze the resulting build separately from the original b6 pilot build.

Plots use all catalogues for population changes and conservative periodic
positional matches for individual-property changes. Candidate/output IDs are
not persistent object identifiers. Match fractions, ambiguous/unmatched objects,
invalid and unresolved values, mass cuts and sparse-bin exclusions accompany
the figures. Rmax is not present in the historical 24-column catalogue and is
not reconstructed from concentration.

Direction comparisons normalize each major-axis vector and remove its arbitrary
sign. Both catalogues must have finite, valid transverse ratios with
`max(b/a,c/a)<0.9`. Checking both ratios also handles reversed historical labels;
the raw ratios are retained. This is a conservative proxy for a well-separated
major direction, not an eigenvalue-gap measurement. Plot-ready schema 3 prevents
older direction masks from acquiring the revised interpretation on replot.
The z=0 literature mass marker is displayed as a common visual reference at
all epochs; it provides no separate z=1 or z=2 convergence evidence.

The tracked vector figures are kept in `slides/figs/`. Python is always executed
with `micromamba run -n cosemu python3 -B`; native executables use the Intel
2024.2 runtime saved before entering the Python environment.

## Reproduce results and recover archived evidence

From this experiment directory, regenerate the numerical summary with:

```sh
micromamba run -n cosemu python3 -B summarize_validation.py
```

This also writes the ignored `slides/validation_results.tex` macros. Rebuilding
the optional presentation requires its local LaTeX source and theme assets;
those files are not included in a fresh checkout of the current Git tree.

`compare_properties.py --help` documents the comparison/replot commands.
The consolidated `comparison-plot-ready.npz`, comparison summaries, figure
receipts and all source catalogue arrays remain under this directory. See
[PLOTTING_NOTES.md](PLOTTING_NOTES.md) for field definitions and matching tests.
Rebuilding a deliverable changes its hashes; regenerate its presentation
verification receipt before treating it as the same checked PDF.

The final PDF was checked on all 16 pages. It has no overfull/underfull boxes,
missing assets, undefined references or visually missing glyphs. The inherited
template still reports fallback CJK families (this deck contains no CJK text),
a TU Palatino text-shape fallback, a removed translation command in an appendix
PDF bookmark and null-byte cleanup of an Inconsolata font name. These warnings
do not clip or omit slide content; the required template font declarations
are retained. Exact source, figure and PDF hashes plus build evidence are in
[presentation-validation.json](presentation-validation.json).

Completed small build files, submitted scripts, logs and intermediate catalogues
are consolidated into `work-artifacts.tar.gz`. The archive is read back and
every member checked before removal of its corresponding loose file. All raw
membership tapes, their indexes, PM headers, particle snapshots and large
density tapes remain in place. `work-archive.json` records the member hashes,
removed paths and retained data; `cleanup.json` also records removed clean
agent worktrees. Historical audit and repair archives remain unchanged.

Cleanup job **11948648** consolidated **559 files**, removed **41 empty
directories** and retained **32 data/log files**. It used one shared core for
87 s (0.0242 core-hours), with 1.32 GiB peak RSS against a 2 GiB request.
Three clean agent worktrees, four rebuildable bytecode files and the 16
temporary slide previews were also removed. The archive's SHA256 and the two
historical archive hashes were checked again after cleanup.

To recover the small files, inspect `work-archive.json` first, then restore the
archive **from this experiment directory**:

```sh
tar -xzf work-artifacts.tar.gz
```

Restore before reusing the original native launch/replay drivers. Original
numerical-control receipts retain their historical absolute worktree paths;
[relocation.json](main-thread-control/relocation.json) maps those paths to the
corresponding `numerical-threading/work` and `main-thread-control/work` paths
here. The nested `work/plotter-evidence.tar.gz` and
`plotter-evidence-archive.json` preserve the earlier plotter test files.
