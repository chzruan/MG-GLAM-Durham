# BDM duplicate-host diagnosis, 2026-09-05

Branch: `fix/bdm-duplicate-hosts-20260905`; base: `9294bcd`.
All new diagnosis products are confined to this directory. Production data,
snapshots and mirrors were opened read-only. The working finder and the archived
2021 finder are byte-identical before the fix, SHA256
`8e1c143fc6701ae2c9a72f3c610b01fef8db2f6e8646d6bfa3fdbb6c77651238`.
An immutable source copy is in `provenance/PMP2linker.original.f90`.

## Root cause and source chain

Line numbers below refer to that original copy of `PMP2linker.f90`.

| Stage / routine | Location | Exact behavior |
|---|---|---|
| `BDM` | 71--174 | FindMaxima, rescale/buffer particles, FindDistinctCandidates, ParametersDistinct/GetHalo, rebuild maxima list, RemoveDuplicates, WriteFiles. |
| `FindMaxima` | 1092--1272 | A cell must have `FI > Ovdens/3` and must not be strictly lower than any of its 26 periodic neighbours. Equal-density plateau cells are eligible. Accepted cells have `Dmaxim+1` added to `FI` in the same OpenMP pass that reads adjacent `FI`: a separate race. |
| `FindDistinctCandidates` | 1008--1084 | Four centroid iterations, averaging positions and velocities of all particles with `d² < Radius²`; `Radius=Cell*max(0.5,min(2,log10(Xoff+10)/2))`. `Xoff` here holds original peak density. Coordinates are wrapped into the box. `SizeList` starts `Cell=2*Box/NGRID` (0.5 Mpc/h here), enlarged by factors of 1.1 if its memory estimate requires it. |
| `GetHalo` | 722--1002 | Log shells with `dLogR=0.02`, overdensity-interpolated aperture, extended by `(Box/NGRID)*min(Rext/(aR/(Box/NGRID))**SlopeR,0.75)`. Production config: `iVirial=1`, `Rext=0.6`, `SlopeR=1.667`, `MassMin=5e11`, `NradP=100`, `dLogP=0.01`. Bound count uses `ee=-Fi(ii)+0.5*|v-vbulk+H(a)*a*dx|² <= 0` inside the **unextended** `Rvir` at that stage. This is one energy classification, not iterative removal and potential recomputation. Final `Mvir=Ncount*MassOne`, `Mtotal=Mtot` in the extended aperture, `Rvir=Radius` (extended). |
| `RemoveDuplicates` | 634--673 | Search linked cells within 3.5 Mpc/h; remove `ip` iff `d²<Rvir(jp)²` **and** `Mvir(ip)<Mvir(jp)`. Centres are directly subtracted without periodic minimum image. The pass reads and zeros shared `Mvir` concurrently. |
| `WriteFiles` | 347--389 | Require centre in output domain and `Mvir>=max(MassMin,20*MassOne)`, then assign `iHalo=iHalo+1` in candidate-array order. Write this as `Nhalo`; write `Nparticles=Mvir/MassOne` and hard-code `Distinct/Sub=0`. No particle IDs or original candidate keys are emitted. |

Two candidates with identical bound sets have equal `Mvir`. Strict `<` is false
in both directions even at zero separation; both survive with new output IDs.
Identical `Mbound` and `Nparticles` are **not independent membership evidence**:
they encode the same particle count. Distinct sets of the same cardinality also
have equal mass in this equal-particle-mass simulation.

`tests/source_reproduction.py` extracts the actual source routines and inserts
only a read-only tap at the bound-particle predicate. On a synthetic system of
128 particles, candidates 1 and 2 converge onto the same 64 particle IDs;
candidate 3 contains the other disjoint 64 IDs. All three survive the original
routine with 1, 2 and 4 OpenMP threads. See `provenance/regression_results.json`.
This proves the failure mechanism for the synthetic set, not membership of
individual production pairs.

## Box 1 reproduction

Source: `/cosma8/data/dp203/dc-ruan1/DESI_MGx100/data/GR/Run1/CATALOGS/CatshortV.0137.0001.DAT`.
SHA256: `8ba7290efaf853636825f406fd76042fe71df0362aada50a3aa88d621f123fa4`.
There are 3,385,000 data rows, of which 1,882,717 have `log10(Mtot)>=12.4`.
The nominal z=0.25 file actually records `a=0.79755`, so the physical redshift
is approximately 0.25384; downstream labels refer to the nominal snapshot.

Periodic pairs, selected mass sample, with strict upper distance boundaries:

| r [Mpc/h] | All pairs | Same Mbound and Nparticles | Fraction |
|---|---:|---:|---:|
| <0.001 | 1,135 | 1,135 | 100% |
| <0.05 | 6,400 | 6,398 | 99.9688% |
| <0.2 | 6,930 | 6,916 | 99.7980% |
| <0.5 | 7,974 | 6,949 | 87.1457% |

Exact defining-field duplicates remove 484 selected rows. Freyja's strict rule
removes 5,967 selected rows (0.31694%), matching the archived audit exactly.
Applying the rule before any mass selection removes 6,191 rows, including 5,969
above the mass threshold: two components cross the threshold and their lowest
original row lies below it. These different selection frames must not be
confused downstream. Full diagnostics and mass bins: `provenance/box1-audit.json`.

The read-only pilot was job `11942773`, cosma8-serial/dp004, 1 CPU, explicit 4 GiB,
15 minute limit. Elapsed 14 s, TotalCPU 12.768 s, MaxRSS 815.51 MiB, CPU efficiency
91.2%, allocated and billed CPUs both 1. Actual core-hours 0.00389; limit 0.25.
See saved sacct/scontrol and script hashes under `provenance/`.

## Proposed repair and limits

Use immutable candidate measurements to form a connected-component graph with
periodic separation <0.2 Mpc/h, identical bound mass/count, speed difference
<5 km/s, and |dlog10 Mtot|<0.005. Retain the lowest original candidate index;
apply the drop mask only after graph construction. Keep the existing unequal-
mass host test as a separate read-only decision pass, then apply its mask.
This removes the deletion race and avoids a blind `<` to `<=` replacement.
The numerical-duplicate graph must not alter survivor measurements.

For historical files, use exactly Freyja's strict rule and lowest original
data-row index, over the **entire** source catalogue. This maintains a common
definition instead of adding an unvalidated fractional-radius threshold.
An exact-field-only sidecar provides a conservative sensitivity bracket.
Future particle-set comparisons are preferable; no production membership
validation fraction can currently be measured. Do not turn same particle
**count** into a claim of same particle **set**.

The strict rule leaves **1,041** same-bound-mass/count pairs below 0.2 Mpc/h
in the full box 1 catalogue, but zero pairs satisfying all strict conditions.
They fail the velocity and/or total-mass conditions. A demand for zero of the
broader class conflicts with conservative merger handling; keep these pairs
explicitly listed for review instead of deleting them to meet that target.

For mass >=12.4, peculiar pairwise radial velocity dispersion at 0.5--2 Mpc/h
rises from 186.950110 to 187.013296 km/s (382,067 to 379,805 pairs).
RMS instead decreases from 281.956030 to 281.869118 km/s: state which statistic
is used. The sign at those separations is not guaranteed by removing edges
below 0.2 Mpc/h, because it changes weighting of other halo pairs.

## Snapshot availability and rerun costs

No `PMcr*`, particle snapshot or particle-membership products were found in the
entire DESI_MGx100 tree. The `degrace_pilot/snapshots` directory contains summary
correlation HDF5 files, not particle snapshots. Checked training runs also lack
the particle dumps, and their submission scripts explicitly run `rm -f PMcr*`
on completion. The completed inventory of all 130 searched roots is saved in
`provenance/snapshot_inventory.json`. Absence is limited to these searched roots, not
a claim about every archive on COSMA.

A historical-snapshot finder replay and production particle-membership test
are therefore unavailable. A fixed catalogue-stage replay can test preservation
of fields, but is not a replacement for a full particle-to-catalogue rerun.
Upstream peak-finding races also prevent promising full-run bitwise equality
across thread counts; this patch must make only its own decision deterministic.

For scale only: one 2048^3 snapshot is about 206 GB for six float32 phase-space
arrays; 100 retained snapshots would be about 20.6 TB. The density phase alone
needs about 481 GB for those arrays plus a 4096^3 float32 grid, before scratch.
Particle and linked-list allocations require a measured memory pilot before
choosing queue and thread count. Archived Run1 timing at step 137 records
11.46 minutes in the broader analysis stage on 128 threads, **not a standalone
finder measurement**. An illustrative 15--30 minutes per finder replay at
128 CPUs would be 3,200--6,400 core-hours for 100 boxes; a one-hour limit would
be 12,800 limit core-hours. These are conditional planning numbers, not verified
runtime estimates or authorization for exclusive-node jobs. Restoring or
regenerating deleted snapshots is additional work and cost. Catalogue cleaning
is feasible on a single shared CPU per file and is measured separately.
