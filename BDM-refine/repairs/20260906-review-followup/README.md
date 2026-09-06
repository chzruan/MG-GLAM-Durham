# Independent-review follow-up: SO v3 and computational repairs

This follow-up addresses the [independent review](../../analysis/review-20260906-merged/REVIEW.md)
of `cz` at `0fa6a5b`. The supplied review and its 166 GNU experiments were
preserved unchanged in `54f53aa`. Repairs were developed and cross-reviewed on
separate branches, then integrated on `repair/bdm-review-followup-20260906`.

The final production source is the integrated `c68c22d` tree. Its
`PMP2linker.f90` SHA-256 is
`39f7e514d1db9fb3b9c7545fc7f3269b4e20b5ee2d1366945f8916b997e4cd5d`;
the makefile SHA-256 is
`6d4a0127d52f22a5645d08ef788440687dfe82218b852b9d81b0e24c8536ec74`.
Later commits retain validation and documentation. The
[source-scope audit](integration-source.json) lists all changed routines.

## Review dispositions

| Finding | Resolution and evidence |
|---|---|
| P1: inconsistent SO normalization | **Repaired.** The threshold is `(4*pi/3)*Ovdens*MassOne*(NROW/Box)^3`, evaluated in float64 from the actual stored particle mass and full periodic box. The catalogue is labelled **BDM finder v3**. [Physics derivation and controls](physics/README.md). |
| P2: outermost crossing and companion absorption | **Convention retained, claim narrowed.** The equal-companion bound concerns that constructed fixture; it does not exclude partial contributions, asymmetric bridging, or supply a cosmological incidence estimate. |
| P3: unmeasured unbinding work | **Measured without truncating convergence.** Per-candidate pass counts and active-row work, plus a log histogram, expose the tail. Published production haloes can exceed 32 passes. The quadratic worst case remains; a silent cap would change the scientific sample. [Diagnostics](unbinding/README.md). |
| P4: alleged 10-to-20 publication-floor change | **Review premise corrected.** The original `a8c7715` writer already enforced `MassMin >= 20*MassOne` before its redundant 10-particle check. Actual old/new writer controls publish the same mass-consistent boundary populations. |
| P5: host-check rejection not exercised at production scale | **Positive controls added.** The unchanged checker rejects one deliberately invalid host pair in each of three complete catalogue copies; all original inputs retain their hashes. [Evidence](host-control/README.md). |
| N1: precise flags lost under overrides | **Repaired.** Compiler-specific precise flags follow `FFLAGS`/`BDM_FFLAGS`, including recursive bitmatch builds. Real GNU/Intel make and runtime controls reproduce the original failure and pass the repair. [Build controls](build-config/README.md). |
| N2: density-sensitive peak/seed channel | **Residual documented.** Fixed-FI byte identity is conditional on identical density input. Ordinary density accumulation can change peak eligibility, strict ties, and the seed radius. Seed quantization adds rounding boundaries and would not guarantee peak stability; it is not introduced here. |
| N3: repeated coordinate work in `List` | **Repaired.** Cache each z cell once, using int16 when representable and int32 otherwise; preserve exact link order, guarded conversion and memory admission. Slab cache scans still cost `O(threads*Np)`. [Particle controls](particles/README.md). |
| N4: `AddBuffer` prefix scratch and serial work | **Repaired.** Byte image counts and 16,384-row blocks retain exact image/particle order while reducing prefix scratch by 7.00 GiB at N1024. |
| N5: inaccessible memory limit | **Repaired.** `BDM.config` accepts finite positive `MaxMemory` in GiB, resets its default to 500 on each read, and reports it. This is finder admission accounting, not a Slurm request or total-process RSS limit. |
| N6: unsupported legacy entry points | **Made explicit.** Dormant profile/simple-duplicate routines fail before allocating/indexing stale workspace or creating files. Unsupported rescaling flags direct restoration callers to `RemoveBuffer`. |

The change preserves the four existing overdensity conventions, particle-mass
constructor, empirical `Rext`, and eight-header/24-column catalogue layout.
`Mbound` describes the converged bound population inside the unextended SO
sphere. The reported radius and `Mtotal` still include the empirical aperture
extension; they must not be presented as an ordinary unextended SO mass/radius
pair. Correcting the approximately 1.03% normalization mismatch can change
discrete membership and selection, so rescaling old catalogue radii is not an
adequate substitute for a finder replay.

## Targeted and independent validation

All Python drivers run through `micromamba run -n cosemu python3 -B`.
Receipts identify their exact tested source; earlier isolated receipts are not
relabelled as tests of later combined source.

| Validation | Result |
|---|---|
| Complete-box SO normalization, all four modes, units, old/new floor and a real edge-membership change | 110 cases per GNU/Intel compiler; the constructed edge changes 257 to 256 bound particles. |
| Existing halo algorithms adapted to coherent v3 mass/box metadata | 242 cases per GNU/Intel compiler on the final combined finder, including all 166 supplied review experiments and 144 duplicate-merge cases. Historical fixtures/receipts are preserved. |
| Peaks/header checks | 136 GNU cases with the v3 version marker. |
| Particle list/buffer controls | 342 focused, 104 legacy and 48 benchmark cases per GNU/Intel compiler; 256 successful focused output hashes also agree across compilers. |
| Build/configuration/entry-point integration | 104 configuration/entry checks, six direct Intel precedence controls and nine actual make/runtime builds, including `FC="gfortran -m64"`. |
| Unbinding counters | 12 comparisons/24 executions per v2 and v3 receipt; all properties and exact IDs preserved, including the 91/191-pass energy ladders and empty-call resets. |
| Full linked native preflight | 21 controls covering all three replay variants, populated/empty repeated calls, exact PM restoration and workspace release, malformed density/header rejection, late-invalid publication preservation and real production/diagnostic output identity. |
| Independent integration review | Actual memory/list/init/release checks and native-build identity review; separate real-parser controls for candidate matching, membership overlap, property columns/cuts, periodic shifts and deliberate provenance mismatches. |

See [native build and protocol](native/README.md),
[memory/particle integration review](integration-review/particles.md), and
[comparison integration review](integration-review/native-compare-review.md).
Intel checked fixtures use bounds/pointer checks where appropriate; compiler
MemorySanitizer startup failed before the baseline program entered its tests,
so these receipts do not claim MemorySanitizer coverage.

## N1024 replay design

The native replay uses the completed validation simulation's **1024^3 particles,
2048^3 force mesh, 512 Mpc/h box**, and snapshots at **z=2, 1, 0**. No second
cosmological evolution is needed to isolate finder changes. The original
snapshot and saved density hashes are checked before use and fully rehashed
before group completion. The replay wrapper verifies every original particle
bit and all finder workspace release after each call, and restores the exact
same density bits before the next call.

Three executables use the same frozen common objects and precise Intel flags:
the v2 reference with output-invariant telemetry; the final computational code
with only the SO expression/version marker reverted in a diagnostic copy
(`optimized-v2`); and the final v3 production finder. This separates numerical
identity of the computational repairs from the intentional normalization change.
The real production make target also agrees byte for byte with diagnostic v3
on the native preflight. Instrumentation is limited to guarded post-publication
tapes and the replay wrapper; the selection arithmetic is unchanged.

At z=0, each variant runs twice at 32 and 64 threads on the same saved FI.
At z=2 and z=1, each v2/v3 pair shares an exact density tape generated once.
The paired comparison matches the unchanged initial density-peak candidate ID,
not the row number in the filtered published catalogue. Unmatched candidate
IDs alone do not establish the appearance or disappearance of physical objects.
The launch [plan](native/plan.json), [submissions](native/submission.json), and
[analysis submission](native/analysis-submission.json) retain hashes, resource
measurements, job IDs and dependencies.

## Native results

Both replay groups completed successfully. All ten stages report zero exact
duplicate member sets, repeated original IDs within a halo, mass/count
mismatches, and selected-host exclusion violations. Every catalogue's 23
physical columns are linked to the raw properties, with the particle mass
checked against the snapshot; the remaining column is the published row ID.
The host counts are still null results, complemented by the separate
production-size positive controls described above.

The computational changes alone give byte-identical v2 catalogues, raw member
tapes, and unbinding diagnostics. At fixed FI, all three variants' catalogue,
member and unbinding bytes also agree between 32 and 64 threads. All z=0
repeated calls produce identical catalogues, restore every original PM bit,
and release the finder workspace. The v2 z=0 catalogue reproduces the earlier
saved-density validation catalogue exactly.

| Redshift | V2 published haloes | V3 published haloes | Count change |
|---|---:|---:|---:|
| 2 | 52,370 | 51,879 | -0.938% |
| 1 | 119,688 | 118,788 | -0.752% |
| 0 | 147,042 | 146,131 | -0.620% |

The reference and v3 runs use `MassMin=2.5e12` solar masses/h, `iVirial=1` and
`Rext=0.15`. These changes isolate the SO normalization on the same snapshot
and FI; they are additional to the earlier legacy-to-v2 refinement results.
They describe this box and its three snapshots, not resolution convergence or
a general cosmological bias estimate.

### Computational cost

These z=0 medians compare reference v2 with optimized-v2, keeping the physics
and output bytes identical. Each median contains two stage samples.

| Threads | Stage | Reference seconds | Optimized seconds | Speedup |
|---|---|---:|---:|---:|
| 32 | List | 9.460 | 3.235 | 2.92x |
| 64 | List | 13.315 | 3.960 | 3.36x |
| 32 | AddBuffer | 2.825 | 1.885 | 1.50x |
| 64 | AddBuffer | 3.570 | 1.665 | 2.14x |

The new List still takes longer at 64 than at 32 threads in these measurements;
the cache reduces repeated coordinate work and memory traffic but does not
remove the `O(threads*Np)` slab scans. The two thread groups ran on different
shared nodes, so the within-group comparison is the controlled result.
`ParametersDistinct` medians are 171.460 to 157.565 s at 32 threads and 102.110
to 106.565 s at 64 threads. First-pass whole-finder times, excluding diagnostic
tape output, are 182.913 to 187.731 s and 142.860 to 133.154 s respectively.
Consequently no uniform whole-finder speedup is claimed from these few samples.
Final-pass whole-finder timings include diagnostic I/O and are unsuitable for
that claim.

At N1024 the buffer-count/prefix scratch falls from 8,589,934,592 to
1,074,266,120 bytes, saving 6.9995 GiB. The z=0 int16 z-cell cache uses
2.18610 GiB while List runs. The 12-byte-per-candidate unbinding counters and
all list-cache allocations have checked initialization, teardown and memory
accounting. Raw timing arrays and full unbinding histograms remain in
[`results-t32.json`](native/results-t32.json) and
[`results-t64.json`](native/results-t64.json).

### Halo-property response to the normalization correction

The [completed comparison](native/comparison.json) rechecks every stage's
output hashes and proves that peak selection routines and the paired FI/config
are identical. Percentage summaries use raw float32 properties promoted to
float64, require **both bound masses >= 10^12.5 solar masses/h**, and exclude
nonpositive unresolved values separately for each property.

| z | Pairs above both mass cuts | Median Mbound change | Median Mtotal change | Median reported radius change | Median Vmax change | Median Rrms change |
|---|---:|---:|---:|---:|---:|---:|
| 2 | 37,409 | -0.367% | -0.350% | -0.432% | 0.000% | -0.280% |
| 1 | 93,227 | -0.439% | -0.429% | -0.454% | 0.000% | -0.326% |
| 0 | 121,362 | -0.442% | -0.421% | -0.453% | 0.000% | -0.347% |

Vmax has two fewer valid pairs at z=1 and z=0. The table is a population
summary, not a uniform correction: at z=0 the bound-mass change among these
pairs ranges from -44.28% to +2.48%, while its 16th/84th percentiles are
-1.006%/-0.157%. Discrete aperture and subsequent unbinding changes can be much
larger than the normalization offset in individual candidates. Full counts,
means, extrema and 16th/50th/84th percentiles are retained in the JSON.

| z | Matched published candidates | V2-only candidates | V3-only candidates | Exactly unchanged bound sets | Median membership overlap |
|---|---:|---:|---:|---:|---:|
| 2 | 51,547 | 823 | 332 | 13,172 | 0.99624 |
| 1 | 117,407 | 2,281 | 1,381 | 22,054 | 0.99581 |
| 0 | 144,179 | 2,863 | 1,952 | 20,996 | 0.99569 |

Overlap here means `|A intersection B| / max(|A|,|B|)`, not Jaccard, and uses
all matched published candidates. Their peak-centre displacement is exactly
zero; median bulk-velocity changes are 0.810, 0.732 and 0.527 km/s at z=2, 1, 0.
All retained shapes are finite and ordered and their axes are normalized or
the documented zero sentinel. The compact
[`comparison-catalogues.npz`](native/comparison-catalogues.npz) contains all six
published tables and their initial candidate IDs. Raw memberships and indexes
remain available under the paths recorded in the replay receipts.

### Slurm resources

All jobs use shared `cosma8-serial` and account `dp004`; none requests an
exclusive 128-core node. The prior native pilot's 189.05 GiB Slurm MaxRSS
justified 256 GiB per replay allocation. The two independent jobs requested
their actual 32/64 thread counts and could share capacity according to the
scheduler. The serial property analysis ran separately on one core.

| Job | Work | Cores / requested memory | Elapsed | CPU utilization | Batch MaxRSS | Allocated and billed core-hours |
|---|---|---|---|---:|---|---:|
| 11949499 | 64-thread replay group | 64 / 256 GiB | 38m 30s | 52.8% | 188,094,612 KiB | 41.0667 |
| 11949500 | 32-thread replay group | 32 / 256 GiB | 25m 46s | 68.8% | 168,116,104 KiB | 13.7422 |
| 11949535 | Full paired analysis | 1 / 8 GiB | 73s | 78.1% | 8,376,344 KiB | 0.0203 |

All three jobs completed with exit `0:0`: **54.8292 allocated/billed core-hours**
against **96.25 time-limit core-hours**. The original expected replay range was
20--35 minutes per group; the 64-thread group exceeded its upper estimate but
completed within the one-hour limit. Whole-allocation utilization includes
serial hashing/validation, input and tape I/O, and the wrapper's full-bit
restoration checks. The separate host positive-control job used 0.02223
core-hours. [Accounting](native/accounting.json) preserves requested/allocated
TRES, billing, elapsed, TotalCPU, MaxRSS and available controller records.

The analysis process reports 286,620 KiB peak RSS while Slurm records nearly
the full 8 GiB request. The reason for that difference was not measured; a
repeat should request 12 GiB rather than infer headroom from process RSS alone.

## Verified cleanup and retained evidence

After replay and comparison completion, shared one-core job **11949711**
compressed and read back every archived file, checking SHA-256, mode and
symlink target before any removal. It consolidated **441 scratch files/links**
into one **234,277,483-byte** archive and removed **66 empty directories**.
The archive SHA-256 is
`37e35c6bcdd64fcddac9eb6ca0d271ef192033968b770e86d48062ba472a0247`.
The complete [manifest and cleanup receipt](work-archive.json) records each
removed and retained path. Twenty-three entries remain unpacked, including
the large density/member tapes, indexes and archive-job log. Build products,
small fixtures, frozen launch copies, catalogues and diagnostic logs can be
restored selectively from the ignored `work-artifacts.tar.gz` using those
manifest paths; the compact comparison tables remain directly accessible.

The cleanup completed in 65 s with 55.666 CPU s, 1,247,640 KiB batch MaxRSS,
and **0.01806 allocated/billed core-hours**, within its one-core/4 GiB shared
request. [Submission](cleanup-submission.json) and
[accounting](cleanup-accounting.json) preserve the measured pilot, script
hashes and resources. The initial submission was rejected because completed
parent jobs had expired from Slurm's live dependency records; its
[failure receipt](cleanup-submission-failed-dependency.json) is preserved.
The successful retry required the immutable all-completed accounting and
comparison receipts before submitting without live dependencies.

The three completed repair worktrees were removed only after their commits
were merged into this branch and each worktree was clean, including ignored
files: [two COSMA worktrees](completed-worktrees-cleanup.json) and
[the temporary physics worktree](remaining-worktree-cleanup.json). Disposable
development builds were also retired after final artifact identity checks.
Historical review files, earlier validation archives and original simulation
snapshots are preserved. This follow-up adds no presentation source or slide
deck to Git. The [final verification](final-verification.json) rechecks all
43 frozen production inputs, the comparison tables, all four current/historical
archives, retained science-file metadata and supplied review files.
An incremental, one-thread [Git repack](git-maintenance.json) also consolidated
288 loose objects (307 to 19 at that check), preserving refs and the working
tree; the repository connectivity check passed. Unreachable loose objects were
retained, and no history rewrite or immediate pruning was used.

## Remaining limits

This closes the actionable review repairs within the supported equal-mass,
periodic-box distinct-host finder. Unbinding retains a quadratic worst case;
ordinary density accumulation still permits thread-sensitive peaks/seeds;
outermost SO crossing and the empirical aperture remain modelling choices.
The List cache improves measured cost without eliminating all anti-scaling.
The N1024 snapshots test the integrated finder on production-size data but
do not establish resolution convergence or validate every possible cosmology,
input domain, shape model or compiler. Earlier v2 plots and presentations
retain their original meaning; the new receipts describe the additional v3
normalization change.
