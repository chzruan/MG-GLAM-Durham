# BDM audit repairs, September 2026

This directory contains repairs and validation following the two-stage
[physics and numerical audit](../../analysis/full-audit-20260906/AUDIT.md)
of `cz` at `a8c7715`. Development and independent reviews took place on
isolated `repair/bdm-*` branches, integrated through
`repair/bdm-audit-20260906`. The historical audit and its failing experiments
are preserved as evidence against their recorded source revisions.

The active distinct-host finder now has an explicit version-2 catalogue
contract. This changes the interpretation and values of several fields;
historical catalogues should retain their original definitions.

## Physics repairs

| Audit finding | Resulting behavior |
|---|---|
| F01, F15: configuration order, dispatch and metadata | Validate configuration before any dependent work; share the overdensity calculation between peaks and haloes; report the actual snapshot cosmology. |
| F02: translation-dependent centring | Accumulate double-precision local offsets and velocities; canonicalize the centre each iteration; define empty-neighbour behavior. |
| F03, F11, F12: mutable density and unsafe peak buffers | Read an immutable field, break equal-density ties deterministically, then use exact counts and stable offsets. Empty results publish a valid header-only catalogue and return. |
| F04: disjoint hosts merged by numerical field agreement | Exact sorted original-particle IDs establish identical membership. Distinct-host geometry is checked separately. |
| F05: shifted SO interpolation | Use the outermost self-consistent crossing of the discrete enclosed-particle profile over the supported search domain. |
| F06: one-pass unbinding | Recompute the survivor potential and bulk velocity until no particle is removed. Exclude each particle's self-potential. |
| F07, F08: inconsistent property populations and shell energy | Use the converged bound set for bound mass, bulk velocity, kinetic and potential energy, RMS radius, shape, spin and circular-velocity profile. |
| F09, F16: missing images and heuristic capacity | Count periodic images exactly, include half-open faces, and cover every supported search. Reconstruct image geometry from original coordinates in double precision. |
| F10: inline analysis changes integrator state | Retain the original six allocations and move them back after analysis; counts, array sizes and every float32 bit are preserved. |
| F13: invalid central/unresolved arithmetic | Guard domains, terminate bounded solvers, and report finite unresolved values or explicit candidate rejection. |
| F14: incorrect principal direction | Use a scaled symmetric Jacobi eigensolver, validated by independent LAPACK eigenvalues and eigenpair residuals. |

Three follow-up reviews found additional defects and produced targeted repairs:

- A sampled SO search could select different crossings when `Cell` changed.
  Starting from the complete supported cap and monotonically tightening the
  enclosed-mass bound finds the outermost crossing. After 16 contractions,
  a sorted interval scan bounds adversarial work. Both heapsorts use int64
  indices. See [the independent review](tests/halo_review.md) and
  [the repair](tests/halo_followup.md).
- Storing a periodic image as float32 could change a cutoff decision even
  when subsequent subtraction used double precision. Geometry now uses the
  original representable coordinate plus its integer image shift. See
  [the periodic precision report](tests/halo_PERIODIC_PRECISION.md).
- Writing headers and rows directly into the final catalogue could replace
  a valid result with a partial file. Same-directory staging and atomic
  replacement publish only a successfully closed, validated catalogue.
  See [publication_notes.md](publication_notes.md).

The first combined 128^3 native replay then exposed a further equal-mass
host-exclusion hole. Three surviving pairs shared 2,109/2,110, 326/329 and
1,441/1,445 particles. Their centres were separated by only 0.0035, 0.0245
and 0.0155 Mpc/h, well inside their approximately 1.216, 0.699 and
1.082 Mpc/h radii. Unequal-mass exclusion skipped them; exact membership
identity could not remove them. The full numerical evidence is in
[equal-mass-discovery.json](equal-mass-discovery.json).

Host priority is descending stored bound mass, then ascending original candidate
index for an exact mass tie. A lower-priority centre strictly inside the
higher-priority host's radius is suppressed. The mask reads immutable
measurements, including through host chains. Exact-membership removal is
an additional identity check. Distinct objects outside the relevant host
radius may share outskirts. Equal-mass objects with unequal radii follow
the same priority convention; the lower-priority radius does not reverse it.
Host displacements and shifted cell bounds use double precision, including
the periodic boundary-cell miss reproduced during this follow-up.
The [host-priority report](tests/halo_HOST_PRIORITY.md) records 160 GNU and
160 Intel cases, including the actual native pairs and both periodic failures.

## Catalogue definition and supported domain

- Input particles have equal mass and the snapshot's `NROW^3` population.
  Background formulae retain the existing flat matter-plus-Lambda convention.
- SO is defined using the code's existing density normalization and requested
  overdensity mode. The empirical `Rext` enlargement still defines the
  reported radius and `Mtotal`; bound membership lies inside the unextended
  SO sphere. Consequently these three fields are not an unmodified standard
  SO mass/radius pair. `dLogR` remains accepted for legacy profile settings
  but no longer determines the active discrete SO crossing.
- Unbinding uses an isolated, spherical Newtonian potential. It is not an
  exact aspherical pair-force calculation or model-specific modified-gravity
  binding calculation. Survivor-only energy and bulk velocity use the same
  membership. The legacy spin proxy, axis corrections and resolved Vmax
  correction remain empirical conventions.
- Searches are strictly smaller than half the periodic box. A cap that
  remains overdense, or a corrected aperture beyond the verified domain,
  rejects that candidate explicitly. Multiple particles exactly at the
  centre are rejected because their unsoftened potential is singular.
  An unresolved Vmax/Rmax has the finite zero sentinel, including zero
  concentration. Five candidate-status counts are printed in the run log.
- The output retains eight header lines and 24 columns, with a `[BDM finder v2]`
  identifier. Continuous quantities keep their storage/output dtypes.
  Invalid finite/domain checks prevent publication of corrupt catalogues.
- `HaloProfile`, `GetProfiles`, `WriteProfiles`, and `RemoveDuplicatesSimple`
  remain dormant legacy routines. They are not supported entry points or
  validated subhalo algorithms; the active BDM path does not call them.

## Validation

The test drivers extract the actual production routines. They do not replace
their arithmetic with a second implementation. Independent references include
sorted SO intervals, compensated pair-potential sums, geometric pair energies,
minimum-image memberships, bitwise state comparisons, exact integer identities,
and NumPy/LAPACK eigenpairs. Checked builds trap invalid/zero/overflow arithmetic
and check bounds; optimized builds exercise the production computation choices.

Combined-source GNU verification passed 136 peak/configuration/empty-result
experiments, 104 particle experiments, 36 halo-property experiments, 26
publication experiments, 26 SO/heap follow-ups and 30 periodic experiments.
The core suite additionally checks 129 shape matrices and six concentration
inversions per build, domain failures, 257 identity sets, writer rejection
and duplicate controls through eight threads. Peak cases extend to 16 threads.
The original focused particle and periodic branches also passed Intel tests.
Exact tested source revisions and driver hashes are retained in
[integrated-validation.json](integrated-validation.json) and the focused receipts.

After the final host and coordinate changes, ten focused suite invocations
passed again on the final production SHA, including 160 host-priority and
30 periodic cases under each compiler and five adversarial native-driver
controls. No test extraction fixes were needed. Their complete commands and
hashes are in [final-focused-validation.json](final-focused-validation.json).
The earlier full 36-case aggregate remains applicable to the unchanged
`GetHalo` implementation, with its changed coordinate helpers covered by
the final periodic and native tests.

The ten historical Python catalogue-cleaner tests also pass, including the
previous HDF5 float64 comparison and interrupted-receipt recovery repairs.
That historical numeric cleaning policy is distinct from the new finder
membership/host definition. Its archived finder-reproduction drivers describe
their original source revision; use the repair tests for the revised finder.

Native validation builds the full repository with ifx 2024.2.0, plus a checked
finder, a read-only membership tap and a state-restoration entry point.
The finder alone uses `-fp-model precise`; the simulator's existing flags are
unchanged. Frozen source and binary hashes are recorded before replays.
The tap must leave the published catalogue byte-identical. The state probe
calls BDM twice, checks every original particle bit and array size, and checks
that per-call membership buffers have been released.

The 64^3 and 128^3 z=0 GR snapshots are reused from the original audited
simulations after SHA256 verification. Replays check all 24 output fields,
unique original memberships, mass/count consistency, host exclusion, checked
versus optimized agreement, and repeat/thread agreement. These boxes share
particle resolution; they are integration checks, not a resolution-convergence
or production-scale calibration. General mass-function, concentration,
clustering and modified-gravity calibration remains outside these repairs.

Final native source SHA256:
`868bb0d6c89e5daf1ca0f631419d86577dcfe31dde7969375398cb0408140fc0`.

| Particles | Original `cz` rows | Repaired rows | Exact duplicate sets | Host-exclusion violations | Checked/optimized fields |
|---|---:|---:|---:|---:|---|
| 64^3 | 94 | 154 | 0 | 0 | Identical |
| 128^3 | 757 | 1,215 | 0 | 0 | Identical |

Both membership taps leave the catalogue byte-identical, and both two-call
state probes preserve every input particle bit and array size. The three
equal-mass secondary candidates are removed; every retained N128 membership
and all 21 raw properties are unchanged by the final host/coordinate changes.
One remaining pair shares a single outskirts particle, with both centres
outside the relevant host radius; this is permitted by the distinct-host
definition. See [native-pilot.json](native-pilot.json),
[native-n128.json](native-n128.json), and
[native-host-repair-comparison.json](native-host-repair-comparison.json).
The row-count changes combine several physics and configuration repairs and
are not a measured correction to a production halo mass function.

## Resources, efficiency and reproduction

The user authorized light login-node work when idle. A one-second resource
sample showed about 70% idle CPU across 128 logical CPUs with approximately
886 GiB available. Native builds and small replays used one core; parallel
timings use the shared `cosma8-serial` partition with explicit CPU and memory
requests. See [login-resource-check.json](login-resource-check.json).

Memory traffic is reduced by exact image allocation, one fill followed by
`move_alloc`, exact restoration of original allocations, threshold prefiltering
before neighbour comparisons and exact peak compaction. Linked-list z-slab
prefiltering avoids unnecessary x/y work; its outer scans remain O(T*Np).
The tested bucket prototype was slower at two/four threads and was discarded.
Earlier list-kernel speedups are documented in
[the particle report](tests/particles_README.md), not presented as full-finder
speedups. Complete-cap SO search and iterative unbinding perform substantially
more physical work than the old faulty finder; full-finder timings must be
reported separately.

An exact fast path avoids reconstructing an original particle from itself.
Four alternating one-core N128 pairs reduce median runtime from 5.599 to
3.970 seconds (1.410x, or 29.1% less time). All complete catalogues and the
raw membership/property streams are byte-identical. Twenty GNU/Intel periodic
controls also pass. This is a comparison against the repaired finder before
the final host-tie fix, not against the faster but incorrect historical finder.
See [coordinate_notes.md](coordinate_notes.md).

The final N128 shared-node sweep (job `11947875`) gives:

| Threads | Median wall time (s) | Speedup | Parallel efficiency |
|---|---:|---:|---:|
| 1 | 4.036 | 1.00x | 100% |
| 2 | 2.307 | 1.75x | 87% |
| 4 | 1.238 | 3.26x | 81% |
| 8 | 0.929 | 4.35x | 54% |

All twelve complete catalogues are byte-identical. These are short shared-node
measurements; they do not predict production-scale speedups. See
[scaling-summary.json](scaling-summary.json) and
[native-scaling-t8.json](native-scaling-t8.json).

The job requested and received eight CPUs and 512 MiB on `cosma8-serial`,
with a two-minute limit. It completed in 30 seconds, used 60.567 CPU seconds
and cost 0.06667 billed core-hours, against an expected 0.08889 and a limit
of 0.26667 core-hours. Allocation CPU utilization was 25.2%, reflecting the
intentional serial/thread sweep within one allocation. Finder MaxRSS peaked
at 299,540 KiB, while the whole batch reached 512,088 KiB. Future identical
sweeps now default to 768 MiB to leave headroom for the complete workflow.
The submitted script/driver versions are frozen in the local archive;
[jobs.json](jobs.json) retains hashes, dependencies, `scontrol`, and `sacct`.

Run the focused tests with the required environment, for example:

```sh
micromamba run -n cosemu python3 -B BDM-refine/repairs/20260906/tests/core_regressions.py
micromamba run -n cosemu python3 -B BDM-refine/repairs/20260906/tests/halo_followup.py
micromamba run -n cosemu python3 -B BDM-refine/repairs/20260906/tests/halo_periodic_regression.py
```

For a new native build in an empty `work/native-build` directory:

```sh
export LINES=40 COLUMNS=120
module purge
module load intel_comp/2024.2.0
module load compiler-rt tbb compiler
export BDM_AUDIT_NATIVE_LIBS="$LD_LIBRARY_PATH"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
micromamba run -n cosemu python3 -B BDM-refine/repairs/20260906/build_native.py
micromamba run -n cosemu python3 -B BDM-refine/repairs/20260906/native_validation.py --phase pilot
```

Native subprocesses use the matching Intel runtime; Python remains in cosemu.
The archived original audit snapshots are required. Parallel scaling must be
submitted with matching `--cpus-per-task` and driver `--threads`, after a measured
pilot. The driver checks allocation size and runs three interleaved repeats
through the requested maximum thread count on the same shared node.

Builds, verified snapshots, pre-fix native outputs, full regression receipts
and submission scripts are consolidated into the ignored local
`work-artifacts.tar.gz`. Each archived file is verified before its loose copy
is removed. Compact reports and final catalogues/memberships remain tracked.
Disposable agent worktrees are removed while their branches remain available.
See [cleanup.json](cleanup.json) and [validation.json](validation.json) for
archive verification, Git object packing and final source/artifact checks.
