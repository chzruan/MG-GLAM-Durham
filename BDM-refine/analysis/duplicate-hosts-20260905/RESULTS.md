# Completed BDM duplicate-host repair

Both hand-off prompts were executed on `fix/bdm-duplicate-hosts-20260905`,
without merging into `cz`. All **740 nominal z=0.25 catalogues** and their
original-row cleaning lists are complete and validated. Production catalogues,
snapshots, mirrors and downstream emulator inputs were not modified.

## Cause, repair and evidence

The original `RemoveDuplicates` in `PMP2linker.f90` removed a candidate only
when it lay inside a more massive candidate's radius and had strictly smaller
bound mass. Two candidates with the same bound particle set therefore survived
the equal-mass tie, even at zero separation, and received distinct output IDs.
The complete source chain and exact original line references are in
[DIAGNOSIS.md](DIAGNOSIS.md).

The patch makes the unequal-mass host decision read-only and periodic, then
merges connected components of the conservative numerical-duplicate graph.
The lowest original candidate index survives; its measurements are unchanged.
The postprocessor uses the same `strict-v1` graph: periodic separation
<0.2 Mpc/h, identical Mbound and Nparticles, relative speed <5 km/s, and
|dlog10 Mtot|<0.005. These are exactly Freyja's strict criteria.

The actual source routines reproduce the failure on a synthetic 128-particle
system: two candidates contain the same 64 particle IDs, and a third contains
a disjoint set of 64 IDs. The original retains three candidates; the fixed
routine retains two with unchanged survivor properties at 1, 2 and 4 OpenMP
threads. Periodic-boundary, transitive-component, merger-control and immutable
host-mask regressions pass. ASCII/HDF5 preservation, strict boundaries,
row-mask loading and overwrite protection also pass. The full finder builds
with ifx 2024.2.0 and the repository's default makefile flags.
Evidence: [regression_results.json](provenance/regression_results.json).
The compiled executable is retained as `bin/PMP2BDM.fixed.exe`.

For the reference box 1, replaying the compiled Fortran predicate and component
policy at stored catalogue precision gives exactly the postprocessor's
6,191-row mask. Every retained row is byte-identical in all 24 columns.
This is a catalogue-stage replay, not a historical particle-to-catalogue
rerun; see [policy-replay.json](provenance/policy-replay.json).

## Delivered catalogues and incidence

Products are under [products/strict-v1](products/strict-v1/), mirroring the
original fiducial `DESI_MGx100/data/GR/Run*/CATALOGS/` and training
`mg_glam/DurMun_hmfemu_*/Run*/CATALOGS/` paths. Each ASCII catalogue has one
adjacent `.cleaning.hdf5` sidecar containing masks, all group members and kept
IDs, criteria, source/output hashes, tool commit and validation receipt.
[README.md](README.md) documents the format and downstream loading interface.

| Sample | Catalogues | Original rows, all masses | Removed rows, all masses | Removed fraction at log10 Mtot >=12.4 | Per-catalogue range at >=12.4 |
|---|---:|---:|---:|---:|---:|
| Fiducial LCDM, boxes 1--100 | 100 | 333,224,063 | 599,164 | 0.309276% | 0.298619--0.317271% |
| LCDM, models 1--64, boxes 1--5 | 320 | 1,024,801,853 | 1,745,822 | 0.293975% | 0.206813--0.370065% |
| fRn1, models 1--64, boxes 1--5 | 320 | 1,144,322,737 | 2,161,149 | 0.331556% | 0.213412--0.462908% |

Across the full grid, **4,506,135 of 2,502,348,653 rows** are removed, leaving
2,497,842,518 byte-preserved rows. Cleaned ASCII files occupy 706.890 GB and
the 740 compressed sidecars occupy 216.204 MB (decimal units). The original
eight header lines are retained, with one added comment containing provenance;
original Nhalo values and survivor order are retained.

For all 100 fiducials combined, the 12.4--12.6 bin loses 67,881 of 75,508,527
rows (0.0898985%), and the 13.4--13.6 bin loses 58,973 of 6,380,821 rows
(0.924223%). Detailed deliverables:

- [catalogue_validation.csv](catalogue_validation.csv): all 740 catalogues,
  incidence, pair counts, hashes and exact product/mask paths.
- [mass_bin_validation.csv](mass_bin_validation.csv): all 8,140 catalogue/bin
  records, 0.2 dex bins from 12.4 to 14.4 plus the final 14.4--14.5 bin.
- [campaign_summary.json](campaign_summary.json): 740 expected, 740 complete,
  no missing or invalid products, and totals by gravity/model.

Every catalogue has zero remaining **strict-rule** edges. The broader
same-Mbound/Nparticles condition leaves **755,417 pairs** within 0.2 Mpc/h
across the grid, including **1,041 in box 1**. Their original row indices,
IDs and failed velocity/mass criteria are available in each sidecar's
`validation_json`. The requested zero count for this broader class is not met:
removing it would require a more aggressive scientific rule without membership
evidence, and could remove genuine neighbouring haloes or mergers.

For box 1, the mass>=12.4 pairwise radial peculiar-velocity standard deviation
over 0.5--2 Mpc/h rises from **186.950110 to 187.013296 km/s**; pair counts
fall from 382,067 to 379,805. This is dispersion around the mean. Radial RMS
instead falls from 281.956030 to 281.869118 km/s. All three individual radial
bins also have increasing standard deviation. Values and definitions are in
[box1-audit.json](provenance/box1-audit.json).

## Freyja reconciliation and remaining limitations

Selecting box 1 at log10 Mtot>=12.4 before grouping reproduces Freyja exactly:
5,967 strict removals and 484 exact-field removals from 1,882,717 rows.
The delivered mask covers the entire original source and removes 6,191 rows
(509 exact-field removals). It removes 5,969 rows above 12.4 because two groups
cross the selection boundary. Apply the full-row mask before any mass cut.
All six representative model 1/32/64 LCDM/fRn1 box 1 audits also agree exactly
with the archived downstream counts.

In Freyja's actual HOD range, with no minimum mass and log10 Mtot<14.5,
removed-row counts agree in boxes 1--5: 6,191, 5,961, 5,975, 5,972 and 6,038.
The delivered full masks additionally remove eight rows at or above 14.5 in
boxes 3--5. See [freyja-reconciliation.json](provenance/freyja-reconciliation.json)
and [representative-audit-check.json](provenance/representative-audit-check.json).

No production particle snapshots or membership files were found in the 130
searched roots; the snapshot inventory records the search scope. Consequently,
the production fraction validated as an identical particle set is unavailable.
The rule is catalogue-level and cannot settle genuine-major-merger membership.
Equal bound mass and particle count encode the same cardinality, not two
independent membership measurements.

The upstream FI peak-generation race remains outside this patch. Making the
host pass immutable and periodic can change pre-write host membership relative
to the old racy execution. Full historical finder host-set equivalence and
full-run bitwise equality across thread counts have not been established.
In-memory versus rounded catalogue thresholds can also differ near cutoffs.

If snapshots are restored, the conditional 100-box finder planning estimate is
3,200--6,400 core-hours, with 12,800 core-hours at one-hour limits on 128 CPUs;
this is not a measured standalone-finder estimate. Restoring or regenerating
snapshots adds substantial cost. The memory justification and caveats are in
[DIAGNOSIS.md](DIAGNOSIS.md). Other redshifts were inventoried but not cleaned.
Consistent downstream use still requires remeasuring HMF, halo correlations
and velocity targets, then validating/retraining the affected emulators.

## Resources and cleanup

The verified pilot preceded campaign submission. All Slurm jobs used
`cosma8-serial`, account `dp004`, one billed CPU per worker, with explicit
memory requests and eight-way capped concurrency. Campaign array **11942910**
completed in 68 min 9 s wall time, consuming **7.8864 allocated/billed
core-hours** and 4.4277 CPU-hours (56.14% CPU efficiency). Its eight two-hour
limits represented 16 core-hours. I/O stalls lowered efficiency relative to
the pilot; no worker failed. Slurm batch-step MaxRSS was 4083.95--4084.00 MiB
under the 4 GiB requests; process RSS samples and high-water marks are recorded
separately rather than equating this accounting value with application heap.

Including audit pilot **11942773**, cleaning pilot **11942890**, and compiled
reference replay **11942925**, total usage was **7.9003 core-hours**, versus
16.5 time-limit core-hours. All eleven allocations completed with exit 0:0.
Exact scripts, hashes, dependency evidence, ReqTRES/AllocTRES, MaxRSS and CPU
accounting are retained in [provenance](provenance/) and
[accounting_summary.json](provenance/accounting_summary.json).

After catalogue completion, **974 reviewed historical text logs** were
consolidated into one verified 4.75 MB archive, and their loose copies removed.
Another **139 disposable build/test/cache files or symlinks and empty error
logs**, plus **six directories**, were removed. The compiled finder, regression
evidence and all science data were retained. Recovery information and file
manifests are in [cleanup-receipt.json](provenance/cleanup-receipt.json) and
[cleanup-scratch-receipt.json](provenance/cleanup-scratch-receipt.json).
Workspace file counts excluding `.git` fell from 10,085 immediately before
cleanup to 8,977 after cleanup and final reports: a net reduction of 1,108
files. Both counts already include the required catalogue products. The
counting scope is recorded in
[cleanup-workspace-counts.json](provenance/cleanup-workspace-counts.json).
