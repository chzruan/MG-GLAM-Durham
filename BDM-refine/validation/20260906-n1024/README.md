# BDM production-scale validation and property comparison

This experiment follows the completed [audit repairs](../../repairs/20260906/README.md)
merged into `cz` at `b6e96669`. It evolves a new GR realization with **1024³
particles, a 2048³ PM mesh and a 512 h⁻¹ Mpc periodic box**, and compares the
pre-audit (`a8c7715`) and refined finders on the **same saved particle snapshots**.
The experiment is developed on `validation/bdm-n1024-ng2048-20260906`.

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
10¹²·⁵ h⁻¹ M⊙ motivates a labelled reference mass in the figures; it is not
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

## Validation and plotting procedure

`run_validation.py` generates ICs and evolves the simulation under Slurm, checks
native headers, file lengths and final epochs, and saves snapshot SHA256 hashes.
Refined membership taps are read one halo at a time, keeping memory proportional
to the number of haloes plus the largest halo. Checks require strictly unique
original IDs within each halo, consistent mass/count, no identical surviving
member sets, and no lower-priority centre inside its host's radius, using the
same documented mass/index priority and periodic geometry.

The instrumented, standard and two-call state-probe catalogues must equal the
inline catalogue byte for byte. The state probe checks all six particle arrays
bitwise after each finder call, their allocation sizes and released per-call
buffers. A half-thread z=0 replay checks parallel agreement. Old/new standard
z=0 timings use equal thread counts; membership-dump or state-probe timings
are not substituted for normal finder timings.

Plots use all catalogues for population changes and conservative periodic
positional matches for individual-property changes. Candidate/output IDs are
not persistent object identifiers. Match fractions, ambiguous/unmatched objects,
invalid and unresolved values, mass cuts and sparse-bin exclusions accompany
the figures. Rmax is not present in the historical 24-column catalogue and is
not reconstructed from concentration.

The editable Beamer deck and vector figures are kept in `slides/`. Python is
always executed with `micromamba run -n cosemu python3 -B`; native executables
use the Intel 2024.2 runtime saved before entering the Python environment.
