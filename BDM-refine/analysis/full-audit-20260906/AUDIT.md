# BDM audit: changes required before a general correctness claim

Audited baseline: `a8c7715` on `cz`. Audit branch:
`audit/bdm-physics-numerics-20260905`. Completed 6 September 2026.

The two new GR pilot catalogues have no identical surviving particle sets, but
the finder does **not** pass a general physics or numerical-correctness audit.
Configuration order changes catalogue counts substantially, centring can leave
the particle cloud, and several independent boundary and numerical failures
remain. These are findings against the existing source. This branch adds audit
drivers and evidence; labelled experimental binaries are retained in a local
archive. Production finder source and historical cleaned catalogues are unchanged.

## Evidence and scope

- Read the full `PMP2linker.f90` implementation, its standalone entry point,
  the inline `Analysis` call, and the relevant density, snapshot, timing and
  allocation code in `PMP2mod_density.f90` and `PMP2mod_tools.f90`.
- Ran **46 controlled experiments** using actual extracted Fortran routines,
  including bounds/FP checks, independent particle-energy calculations,
  translation tests, sparse/dense peaks and periodic boundaries.
- Built the full finder and GR simulator with ifx 2024.2.0 and the repository's
  production flags. Ran GR simulations from z=100 to exact z=0 with **64^3 and
  128^3 particles**, meshes 128^3 and 256^3, boxes 64 and 128 Mpc/h, and fixed
  realization 1. The largest run has 2,097,152 particles, below the authorized
  exclusive upper limit of 1200^3 = 1,728,000,000.
- Recorded original snapshot-row IDs, including the identity of periodic
  copies, during the existing binding decision. This instrumentation preserves
  the two baseline ASCII catalogues **byte for byte**.
- Ran full-catalogue repeat/thread checks at 1, 2, 4 and 8 threads, plus native
  bounds/precise-build comparisons and paired performance experiments.

Machine-readable evidence: [summary.json](summary.json),
[unit-results.json](unit-results.json), [replay-results.json](replay-results.json),
[bounds-results.json](bounds-results.json), and [jobs.json](jobs.json).
Failures in the controlled experiments are retained deliberately: a successful
audit-driver exit means the experiments were recorded, not that the finder
passed them. All source-line references below refer to unchanged `a8c7715`.

## Stage 1: physics correctness

### F01 — P1: configuration is applied after peak selection

`PMP2linker.f90:119–126` calls `FindMaxima` before `ReadParameters` and
`SetParameters`. `iVirial` has no explicit initial value. On the observed first
call, the peak pass behaves as mode 0, even when `BDM.config` requests mode 1.
Later inline calls can inherit the previous configuration.

At z=0, Omega_m=0.3089, the pilot logs show peak selection using **647.4587**,
whereas the requested virial catalogue header gives **332.4786**. A trial that
only moves configuration/scales ahead of density and peaks produces:

| Particles | Existing order | Configuration first | Change |
|---|---:|---:|---:|
| 64^3 | 94 | 139 | +47.87% |
| 128^3 | 757 | 1,088 | +43.73% |

These are output-count changes in small diagnostic runs, not calibrated
corrections to the production halo mass function. Added objects still require
the other physics checks. Read and validate configuration before any dependent
work, and use one shared overdensity calculation for peaks and halo properties.

### F02 — P1: float32 centroid sums can put a centre outside its particles

`FindDistinctCandidates`, lines 1168–1174, accumulates absolute positions and
velocities in implicitly single-precision variables. A symmetric 262,144-particle
cloud has a known centre and an extent of only about 0.027 Mpc/h per axis.
Moving that same cloud changes the result:

| True centre, each coordinate | Computed centre displacement |
|---|---:|
| 5.125 Mpc/h | 0.00000302 Mpc/h |
| 1000.125 Mpc/h | 6.54893 Mpc/h |

The distant cloud is reported near (1003.906, 1003.906, 1003.906), outside its
coordinate range. Subsequent centring iterations find no nearby particles and
retain that erroneous centre. This violates translation invariance and the
basic requirement that an average lie within its input range. Accumulate local
offsets and velocities in float64, with a defined empty-neighbour outcome.
The test isolates centring; it does not assume this extreme cloud is resolved
by the PM force mesh.

### F03 — P1: peak selection reads a density field that it also modifies

At line 1297, the parallel peak pass adds `Dmaxim+1` to selected cells while
neighbouring iterations read those cells. A constant-density plateau gives
125 candidates at each tested thread count but **three different spatial peak
sets** across 1, 2, 4 and 8 threads. Sorting the coordinates does not remove the
difference. Repeats within each thread count happened to agree in this fixture;
that does not remove the read/write race.

Use an immutable density field and deterministic tie handling followed by stable
compaction. The existing atomic CIC deposit protects its updates; it does not
make this subsequent unprotected peak pass safe. See the
[OpenMP memory model](https://www.openmp.org/spec-html/5.0/openmpsu14.html) and
[atomic-access rules](https://www.openmp.org/spec-html/5.1/openmpsu105.html).

### F04 — P1 for general low-mass use: the fixed duplicate rule can merge disjoint hosts

`BdmDuplicateRules` and `MergeNumericalDuplicates`, lines 11–21 and 733–795,
apply a fixed 0.2 Mpc/h separation criterion without requiring common particles
or overlapping halo radii. The `halo_disjoint` case measures two isolated
64-particle objects with equal masses and velocities:

- Bound mass of each: 6.4e8 Msun/h; disjoint original particle IDs.
- Measured radius of each: 0.021024 Mpc/h; separation: 0.15 Mpc/h.
- Both centres lie outside the other's halo; the duplicate pass retains **one**.

This uses an explicitly lowered mass selection and concerns the general finder,
not a claim that the old production mass range has this incidence. Numerical
field agreement is insufficient evidence of object identity. For general use,
validate duplicate identity with membership or an explicitly validated physical
exclusion rule. Connected-component transitivity also means a group's endpoints
need not satisfy the original pair thresholds.

The distinct-host convention permits overlapping outskirts while requiring a
distinct centre to lie outside a larger host's virial radius; it does not imply
that all close equal-mass objects are identical.
[CosmoSim BDM definitions](https://www.cosmosim.org/cms/data/halo-finders/).

### F05 — P2: the SO interpolation associates cumulative masses with shifted radii

In `GetHalo:914–924`, cumulative `MassP(i-1)` and `MassP(i)` are paired with
`Radius*10^(i*dLogR)` and `Radius*10^((i+1)*dLogR)`, respectively. The deposited
profile's indices correspond to radii one bin smaller.

For 100,000 particles with M(<r) proportional to r and a true SO radius of
1 Mpc/h **using the code's own density normalization**, the result is:

| dLogR | Measured radius | Enclosed density / requested density |
|---:|---:|---:|
| 0.040 | 0.958014 | 1.089580 |
| 0.020 | 0.977988 | 1.045513 |
| 0.010 | 0.988719 | 1.022949 |
| 0.005 | 0.994264 | 1.011566 |

The default spacing produces a 2.20% radius bias and 4.55% excess enclosed
density in this fixture. `Rext=0` isolates interpolation from the empirical
aperture correction. Pair masses with their actual bin edges and validate
against particle-sorted SO crossings. The defining SO relation is given in
[Bryan & Norman, section 2.1](https://arxiv.org/pdf/astro-ph/9710107).

### F06 — P2: one-pass unbinding does not establish a self-bound survivor set

`GetHalo:984–1048` constructs the potential using all aperture particles, tests
binding once, and never recomputes the potential after removing particles.
In a 1,000-particle shell with a deliberately hot contaminating population,
100 particles are labelled bound. Recomputing the isolated Newtonian pair
potential from those 100 survivors, excluding self-potential and using the
same G and scale factor, finds **39 with positive energy**.

This is a failure of self-bound membership, not a small change in G. Iterate
potential, velocity and membership to a defined stopping criterion if the
published bound mass is meant to represent a self-bound system. The published
BDM procedure explicitly repeats the unbinding calculation using its survivors.
[Klypin et al., Appendix A](https://arxiv.org/pdf/1002.3660).

### F07 — P2: reported velocity dispersion mixes two particle populations

`GetHalo:1019–1042` accumulates kinetic energy and other structural statistics
from **all** particles inside the original SO radius, even after a particle
fails the binding test. `WriteFiles:389` divides that kinetic energy by the
**bound** mass. In the hot-shell fixture, the reported expression gives about
6021 km/s; the input all-particle RMS is 1904 km/s and the retained-population
RMS is 500 km/s (small Hubble terms do not explain the difference).

Define the intended population for each field and use matching numerators,
denominators and bulk velocities. Replacing this with either a bound-only or
all-particle statistic is a catalogue-definition change requiring validation.

### F08 — P2: coarse-shell potential energy overcounts shell self-interaction

The sum at `GetHalo:1062` uses the outer cumulative mass times the whole shell
mass, rather than treating the mass within a populated shell consistently.
In the same thin-shell case, `EpotM=5.41094e18`, while a direct pair sum over
the identical all-particle population gives `2.61395e18` in matching units:
a factor of **2.07**. The innermost bin also has no explicit energy term in
that loop. Validate the energy discretization separately from the potential
used for unbinding before interpreting `2K/Ep-1` as a virial diagnostic.
The thin-shell discrepancy is not an estimate of its incidence in cosmological
profiles.

### F09 — P2: periodic particle coverage is incomplete in valid boundary cases

`AddBuffer:1878–1880` uses `x>Xright` for exterior images, although the primary
domain is half-open. For 128 particles at x=0, y=z=16 in a box of size 32,
it creates **no images at x=32**. A halo near the opposite face can miss them.

Separately, the buffer is fixed at 5 Mpc/h while `GetHalo` may search to
15*Cell and apply an additional radius correction. A radius-6 search centred
at x=0.1 should find 128 particles at x=26.2 via periodic images; the buffered
representation finds **zero**. Use a consistent half-open convention and
either periodic indexing or a verified buffer that covers every actual search.

### F10 — P2: inline halo finding perturbs simulation particle arrays

`RescaleCoords` and `RemoveBuffer:1993–2007` transform the simulation's float32
positions and velocities in place and then invert the transform. A 1,000-row
test changes 19 values in each coordinate component and seven in each velocity
component. Maximum changes are 6.1035e-5 mesh units and 3.8147e-6 stored velocity
units. The transformations are mathematically inverse but not bitwise inverse.

Treat analysis as a read-only operation on the integrator state, or explicitly
document and measure this numerical intervention. No claim about its amplified
effect on a complete simulation trajectory is made here.

## Stage 2: numerical correctness and efficiency

### F11 — P1: peak-buffer allocation is unsafe; zero peaks exits successfully

At `FindMaxima:1306–1338`, `Nbuff=Nmaxima/Nthreads*5` is an average-load
heuristic, and `Mth` is allocated with extent `Nbuff` although indexed by thread.
A single valid maximum triggers bounds errors at 2, 4, 8 and 16 threads.
Uneven populations can also exceed a thread's position buffer. Replace this
with exact counts and offsets, not another oversized heuristic.

The zero-maximum case issues `STOP 'No density maxima found'` with exit status
0. In inline analysis that ends the simulation; in standalone analysis it does
not produce a valid empty catalogue. Empty results need an explicit successful
output path; fatal states need nonzero status.

### F12 — P2: float32 peak marking can lose already-counted candidates

At lines 1297 and 1333, adding a large density offset and later comparing to that
offset assumes the smaller peak survives rounding. With isolated peaks 1e10
and 100, above the chosen threshold, the code counts **two** but stores only
**one**. Arrays retain space for the missing candidate with uninitialized
contents. The plateau fix should also remove this arithmetic marking scheme.
This is an extreme dynamic-range test, not a measured production frequency.

### F13 — P1: central particles and unresolved radii cause exceptions/NaNs

`GetHalo:893,962,1015` evaluates log10(r/Radius) for r=0 before clamping the
integer bin. A particle exactly at the candidate centre raises SIGFPE in the
checked GNU driver. Clamping after log/integer conversion is too late.

At line 1110, `sqrt(1-dRvmax/R)` has no domain check and uses the loop's current
R rather than explicitly the saved Vmax radius. A compact finite cloud produces
NaN Vmax in the optimized driver. `Concentration:1644` has no iteration limit or
finite-input check; with NaN Vmax the optimized GNU experiment times out after
five seconds. Guard these domains and represent unresolved properties explicitly
rather than publishing invalid values or silently changing the physics.

### F14 — P2: the shape eigensolver can return the wrong major direction

`EigenValues:1553–1618` uses one fixed starting vector, scales by the largest
signed component, and divides by quantities that may be zero. Zero and rank-one
positive-semidefinite tensors produce NaNs. More seriously, a positive-definite
tensor with eigenvalues (3,1,0.5) returns those eigenvalues but a direction with
principal-eigenpair residual **2.017**, rather than approximately zero: the
initial vector is orthogonal to the principal eigenvector.

Use a tested symmetric 3x3 eigensolver and check residual, normalization,
ordering and degenerate cases. Sorting eigenvalues alone does not associate
the returned vector with the largest eigenvalue.

### F15 — P2: configuration dispatch, validation and metadata are inconsistent

`ReadParameters:238–263` silently ignores assignments without an inline `!`,
and accepts invalid values such as iVirial=77, dLogR=0 and MassMin=-1.
`FindMaxima:1230` dispatches Abacus mode as **23**, whereas configuration and
`SetParameters` use **3**. A mode-3 unit call leaves the previous overdensity
untouched. `SetParameters:1763` writes undefined local `Oml0`, rather than the
snapshot's `OmL`; both GR pilots advertise Omega_L=0.0000 despite using 0.6911.

Use a single validated configuration contract and actual snapshot metadata;
reject invalid input with an informative nonzero exit. These related mistakes
are currently hidden by implicit typing and the makefile's warning suppression.

### F16 — P2: buffer capacity is estimated from a uniform distribution

`AddBuffer:1822–1823` reserves `(1+4*dBuffer/Box)^3*Np` particles. A valid compact
cloud near a box corner requires eight images per original particle and exceeds
that estimate for the tested box/buffer. The routine stops with exit status 0.
Count required images or grow safely, using integer capacities and explicit
allocation-failure handling.

### Observed compiler and thread sensitivity

The production N128 catalogue's physical columns are identical after sorting
and excluding sequential row IDs in all 12 samples (three repeats at 1, 2, 4,
8 threads). This is a useful integration result, not evidence that the plateau
race is absent.

Native bounds/precise runs complete with 94 and 757 finite rows, and reciprocal
matching finds every baseline halo within 0.01 Mpc/h. The largest matched centre
change is 0.0001 Mpc/h. For N128, one bound mass changes by 0.2044%, one total
mass by 0.03198%, and three radii change slightly. The experiment changes both
optimization level and floating-point mode; it does not identify which result
is physically correct. Near-threshold membership and float32 accumulations
require further numerical control.

The first `-check all` mixed-object build enabled ifx MemorySanitizer and stopped
inside the uninstrumented Fortran runtime during output. Those warnings are
preserved in `replay-results.json` and are **not classified as finder defects**.
The follow-up uses bounds checks and precise arithmetic without MemorySanitizer.
Explicit FP exceptions in the controlled GNU executable are the direct evidence
for F13. Build commands and hashes are in `diagnostic-build.json`.

### Measured scaling and a small optimization trial

Median of three N128 full-finder samples in the same shared allocation:

| Threads | Wall time (s) | Speedup | Parallel efficiency | Reserved core-s | List build (ms) |
|---:|---:|---:|---:|---:|---:|
| 1 | 0.845 | 1.00 | 100% | 0.845 | 13.118 |
| 2 | 0.580 | 1.46 | 73% | 1.160 | 10.094 |
| 4 | 0.428 | 1.97 | 49% | 1.713 | 8.706 |
| 8 | 0.313 | 2.70 | 34% | 2.503 | 7.962 |

The list microbenchmark uses the actual `List` routine and 128^3 particles,
ten rebuilds per sample, with every particle reachable exactly once. Its outer
loop at lines 1487–1488 scans all particles once per partition/thread: O(T*Np)
work. Eight threads give only 1.65x list speedup. Consider deterministic cell
counts/offsets and particle bucketing to eliminate repeated scans.

A separate paired, alternating-order one-thread experiment adds an early
density-threshold check before examining 26 neighbours. Median full-finder
time changes from **1.134 to 0.877 s (1.29x)** and median CPU time from 0.81 to
0.63 s. The tested catalogue is identical. This is an audit-only candidate;
the immutable peak-selection fix must come first, followed by renewed tests.
Short-run timing has visible cache, I/O and shared-node noise; these numbers
are not extrapolated to production-scale catalogues.

Other reviewed costs: `AddBuffer` zeroes oversized arrays and copies all six
particle components multiple times; `RemoveBuffer` sizes its temporary arrays
to the buffered count although it only copies original particles; standalone
BDM restores particles and reallocates FI immediately before exiting. `Label`
and `Lst` contain int64 elements while `Memory` counts four-byte words, and
`Structures.TotalMemory` is never synchronized with that accounting. The fixed
500-GiB list budget is not the Slurm allocation. Use measured RSS and exact byte
counts before redesigning memory limits. The timings' shared `t0` is also reset
by individual routines; external wall/CPU measurements were used for this audit.

## Physics contracts requiring an explicit decision

These are limitations or intended legacy choices, distinct from the reproduced
implementation failures above:

- `Rext` enlarges the reported aperture, Mtot includes particles in that enlarged
  aperture, while bound membership is counted inside the original SO radius.
  The resulting Mbound, Mtot and Rvir are not interchangeable with an unmodified
  M200c/M200m or virial SO definition. Concentration and spin also use mixed legacy
  quantities and empirical resolution corrections. Version any revised schema
  and revalidate downstream users; do not silently reinterpret old products.
- Background formulae assume flat matter plus Lambda through `1-Om0`; the code
  does not use a general background expansion function. The virial fit's stated
  domain is explicit in [Bryan & Norman](https://arxiv.org/pdf/astro-ph/9710107).
  `MassOne` additionally assumes equal-mass particles with Nparticles=NROW^3.
- Unbinding and circular velocity are Newtonian even when particles came from an
  MG simulation. This can be an intentional common catalogue convention, but it
  is not a validation of dynamical binding under each screened force law. A
  concrete example of adapting halo unbinding to an extra force is
  [Hellwing et al., section 3](https://arxiv.org/pdf/1111.7257); its ReBEL potential
  is not a substitute for the f(R), DGP, symmetron or coupled-field potentials here.
- `HaloProfile`, `GetProfiles` and `WriteProfiles` are dormant in the active BDM
  path. They have uninitialized `ih` indexing and unchecked profile-array/binner
  assumptions and are not validated for activation. The unused
  `RemoveDuplicatesSimple` retains its old nonperiodic in-place race. Active
  subhalo finding is disabled; this audit certifies no subhalo catalogue.

## What passed, and what remains outside this audit

The original-row membership taps found **zero identical survivor sets**, zero
repeated periodic-image IDs inside a survivor, and zero mass/count mismatches
among 94+757 baseline haloes. Strict duplicate diagnostics also find zero rows
to remove. These checks validate the two new snapshots, not every possible
input or the historical production catalogues.

The same particle resolution and different box sizes were used for the two GR
runs. This is not a resolution-convergence experiment, independent-finder
comparison, production replay, or validation of all modified-gravity models.
Mass-function, concentration, clustering, velocity and force-resolution
calibration remains necessary after the correctness repairs. No such campaign
was launched while the hard correctness failures remained unresolved.

## Repair sequence and resource record

1. Fix configuration/metadata, immutable deterministic peaks, exact compaction,
   and float64 local centring. Require empty-output and periodic-boundary tests.
2. Define and validate SO bin edges, bound membership convergence, consistent
   property populations, potential energy, duplicate identity and MG conventions.
3. Replace unsafe scalar domains/eigensolver and add finite-output checks;
   establish precise-reference tolerances and analysis noninterference.
4. Reapply and measure the prefilter, deterministic cell bucketing and buffer
   memory reductions against the corrected implementation. Then perform the
   scientifically necessary convergence and downstream comparisons.

All five allocations used `cosma8-serial` and explicit memory/core requests:
11946517 (initial runtime-library failure), 11946524 (resumed N64), 11946543
(N128), 11946554 (membership/scaling), and 11946652 (bounds/prefilter).
The runtime issue was cosemu's older libiomp5 overriding the matching ifx
library; Python stayed in cosemu, while native subprocesses used the compiler
runtime. Completed initial conditions were reused.

Total billed usage was **0.1161 core-hours**, measured CPU usage **0.07273
core-hours**, against a sum of allocation time limits of **1.8833 core-hours**.
No exclusive COSMA8 node was used. The N128 evolution used four cores and
finished in 44.3 s; that allocation achieved about 93% CPU utilization. The
eight-core replay allocation also contained necessary serial checks, explaining
its lower utilization. `jobs.json` preserves script/driver hashes, job IDs,
dependencies, requests, allocations, MaxRSS and accounting, including the failed
pilot. Audit-owned scratch and redundant logs were archived and verified before
removing 185 loose filesystem entries. The archive and verification receipt use
two inodes, giving a net cleanup reduction of 183 inodes. See
[validation.json](validation.json) and [reproduction instructions](README.md).
