# BDM v3 convergence campaign, 8 September 2026

Branch: `validation/bdm-convergence-20260908`, based on `cz` at
`e289c3f7ea28d38390a1663e50bb2784b85abc84`. The user authorized the seven-run
campaign, simulations below 1200^3 particles, and removal/consolidation of
unnecessary working artifacts. Production finder source is frozen. This is a
convergence measurement, not an assumption that the finest run is exact.

**IC choice:** this already launched campaign uses a frozen native GLAM
first-order (Zel'dovich) generator at z_init=100. It is an existing-run exception
to the project's [default 2LPTIC workflow](../../../AGENTS.md) for new
simulations. Resume it with its recorded ICs; its convergence measurements
describe this first-order suite and do not establish convergence of a 2LPTIC
suite.

## Completed results, reviewed 9 September 2026

All seven simulations, all 21 legacy/v3 replay pairs, the full analysis and
the dependent presentation/cleanup stages completed successfully. Initial-mode
controls passed, including the full 1024^3-row comparison. Each of the 21 v3
catalogues passed the independent membership checks: 377,475 published rows
in total, zero exact duplicate member sets, repeated IDs within a halo,
mass/count mismatches or priority-ordered host-exclusion violations. The host
rule uses the higher-priority halo's extended aperture. The checker's examined
count counts neighbours already inside that aperture, so zero is expected on
these passing production catalogues and does not demonstrate test coverage.
The independent review found four reverse-orientation centre-in-aperture pairs,
all outside both unextended SO radii; they are allowed by this priority rule.
It also measured excess memberships of 0.06–0.13%: exact member sets are unique,
but different haloes can share particles. That global rate is not a bound on
individual-halo errors. Summing bound masses counts shared particles repeatedly;
the catalogue is not a partition of the particle set.

The full results contain seven comparison pairs at three redshifts, each
assessed with 100/300/1000-particle floors (63 assessments). At z=0, the
screen with a nominal 300-particle floor gives the following results for
**abundance, median bound mass, median resolved Vmax and completeness**:

| Comparison | Fixed setting | Passing log10 mass interval, Msun/h |
|---|---|---|
| C/E: 512^3 → 1024^3 particles | 2048^3 force mesh | 12.75–14.75 |
| E/F: 2048^3 → 4096^3 force mesh | 1024^3 particles | 14.00–14.75 |
| F/T: every timestep halved | 1024^3 particles, 4096^3 force mesh | 12.50–14.75 |

The screen requires abundance and median bound mass within 5%, median resolved
Vmax within 2%, at least 90% reference completeness, abundance jackknife sigma
at most 5%, and at least 30 objects per required statistic. These are chosen
relative-stability criteria, not absolute-accuracy guarantees.

The passing intervals do not certify shape. At z=0, the maximum absolute
bin-median shape shifts **inside those intervals** are:

| Pair | delta(b/a), % | delta(c/a), % |
|---|---:|---:|
| C/E | 1.31 | 1.59 |
| E/F | 4.11 | 5.86 |
| F/T | 2.92 | 4.91 |

F/T's median c/a shifts reach 4.70–5.42% across the three redshifts; the coarser
A/E particle comparison reaches the table maximum of 8.92% at z=1. All pair/redshift shape summaries
are in `CONVERGENCE.md`, and all floors in `convergence-assessment.json`.
These are the reported empirically corrected axis ratios, not raw tensor
ratios. Shape, velocity, scatter and tails need separate acceptance criteria.

The effective cut is `max(2.5e12, nominal_floor * max(particle_mass))`, followed
by the whole-bin cut. For E/F and F/T, the publication limit requires **1868
bound particles** (1867.22 particle masses), and the first whole bin requires
**2362 particles**. Nominal floors 100, 300 and 1000 therefore give identical
results for those pairs. This identifies the force sensitivity in a sample
already containing more than 2000 particles per halo.

The paired octant jackknife measures spatial resampling variation in this
realization; it can legitimately be zero for identical octant counts. It is
not a confidence interval or evidence of zero ensemble uncertainty. The report
now provides `sqrt(1/N_left + 1/N_right)` as a hypothetical independent-count
relative scale alongside it. Correlated catalogue counts require a covariance
term, so this benchmark does not replace or bound the paired-ratio uncertainty.
High-mass abundance agreement has limited statistical discrimination even when
the measured medians satisfy the screen. The acceptance decisions are unchanged.

Force resolution limits the tested lower-mass z=0 sample. In the E/F bin
12.50 <= log10(M/[Msun/h]) < 12.75, the 2048^3-mesh run has 9.2% fewer haloes;
matched median mass and Vmax are both about 7.3% lower. In contrast, C/E
maximum absolute bin-median shifts are 0.63% in mass and 0.29% in Vmax over
its eight eligible z=0 bins. The timestep conclusion is redshift dependent:
F/T passes only over log10 mass 12.50–13.00 at z=2 and 12.50–13.50 at z=1.
The force response also depends on epoch: in log10 mass 13.00–13.25,
E/F median mass shifts are +1.85%, −0.17%, −4.07% at z=2,1,0; C/D gives
+2.36%, +0.08%, −3.80%. Small z=1 medians lie near a sign transition in
the population response. These are separately matched mass-bin populations,
not tracked haloes across epochs. Both example bins satisfy the mass condition
at every epoch; abundance and Vmax cause their z=0 failures. Wider z=1 force
intervals do not establish stability across epochs. z=0 restricts lower-mass
force comparisons most strongly here; earlier epochs still restrict timesteps.

The force response in median mass agrees between C/D and E/F to 0.35 percentage
points in eight common z=0 bins (0.28 pp for Vmax), supporting approximate
force–particle separability for those statistics. Timestep interactions remain
untested. Per-halo scatter is appreciable: C/E at z=0, log10 mass 12.75–13.00
has mass-shift percentiles −8.30/−0.52/+6.74% (16th/median/84th).
Full tables, scatter and individual criteria are in [CONVERGENCE.md](CONVERGENCE.md)
and `convergence-assessment.json`. No finder-mesh, box-size or 2LPTIC
convergence claim follows from this suite.

The independent [review](../../analysis/review-20260909-convergence/REVIEW.md)
reproduced all 63 comparisons and decisions. [REVIEW-RESPONSE.md](REVIEW-RESPONSE.md)
records the repairs and clarifies the paired uncertainty and host-rule
interpretations. The subsequent [claims review](../../analysis/review-20260909-convergence-claims/CLAIMS-REVIEW.md)
is a second pass by the same reviewer, not a blind independent opinion.
[CLAIMS-RESPONSE.md](CLAIMS-RESPONSE.md) records its reproduced measurements
and necessary qualifications. The retained catalogues now support 114
abundance comparisons: smaller absolute shift in 68 bins for v3, 43 for legacy,
and 3 exact ties; median absolute shifts 2.78% and 3.21%. These correlated bins
establish neither a statistically supported advantage nor equivalence.
Equivalent membership-based legacy matching still requires member lists absent
from the current replay receipts. Candidate-count sensitivity at one finder
mesh cannot rule out attenuation or shared catalogue bias; finder-mesh tests
remain necessary. No new simulation or replay was needed for this claims response.

| Completed stage | Slurm job |
|---|---:|
| D simulation → paired replay | 11955919 → 11955978 |
| F simulation → paired replay | 11955920 → 11956205 |
| T simulation → paired replay | 11955921 → 11956206 |
| A / E paired replays | 11956146 / 11956204 |
| Full analysis | 11956207 |
| Full plots and Beamer | 11956208 |
| Final launch-file cleanup | 11956701 |

Original-campaign Slurm accounting is **2133.48 billed core-hours** and
**1891.57 consumed CPU-hours**, including pilots, replays, initial presentation
attempts and cleanup. Review IC check 11960736 adds 0.1678 billed and 0.0936
consumed CPU-hours, bringing the recorded totals to **2133.65 / 1891.66**.
The original estimate was about 2400 core-hours; original submitted wall limits
summed to 3888.75 core-hours, with another 0.5 for the review check.
D/F/T elapsed times were 7.46/7.96/13.76 hours on
64 cores each. The final analysis took 152 s on one core. The final narrative
and visual review use the completed measurements; large simulations and
membership matching were not repeated. The lightweight review render is
recorded separately in `final-review.json`.

Scientific inputs/outputs remain unpacked. Completed launch paths are retained
in the verified archives described below and can be restored for reproduction.
The branch contains code, results and receipts; the presentation and figures
remain local artifacts excluded from Git.

All runs use GR, L=256 Mpc/h, Omega_m=0.3089, Omega_Lambda=0.6911,
h=0.6774, sigma8=0.8159, z_init=100, and exact outputs z=2,1,0.
`MG_flag=0` and `MG_model=3` are retained together. The master-normalized
initial Fourier modes are shared across particle resolutions; extra resolved
short modes are allowed in finer runs. The origin and Fourier conventions
are controlled explicitly; using the same native random seed is insufficient.
Initial PM coordinates retain the native periodic upper-edge guard, which
subtracts 1e-3 mesh units when a coordinate rounds to Ngrid+1. Its physical
size changes with mesh spacing. This is distinct from unmatched Fourier
modes. The completed all-row E/F check measured a maximum periodic physical
position difference of 6.103515625e-5 Mpc/h and identical physical velocities.
Review job 11960736 reran the current source-bound validator: coordinate-component
RMS is 2.3835e-7 Mpc/h; no component exceeds the pure-roundoff bound and none
has nonlocal or unexplained excess. Maximum inferred displacement is 0.34851
Mpc/h. The receipt tests excess over the roundoff bound for a local, exact-clamp
explanation; it does not measure total clamp incidence. The conservative
source-derived E/F bound, including the native guard, is
1.8310546875e-4 Mpc/h (rounding plus the possible differential edge clamp),
and the B/D bound is 3.0517578125e-4 Mpc/h. These are numerical coordinate
effects included in the mesh comparison; do not claim bit-identical physical
positions across evolution meshes. F/T uses the same mesh and must have
identical initial position words and correctly different staggered velocities.

| ID | Particle count | Evolution mesh | Comparison |
|---|---:|---:|---|
| A | 256^3 = 16,777,216 | 2048^3 | Low particle resolution |
| B | 512^3 = 134,217,728 | 1024^3 | Coarse force mesh |
| C | 512^3 | 2048^3 | Intermediate particle resolution |
| D | 512^3 | 4096^3 | Fine force mesh |
| E | 1024^3 = 1,073,741,824 | 2048^3 | High particle resolution |
| F | 1024^3 | 4096^3 | Finest reference |
| T | 1024^3 | 4096^3 | Each F step subdivided into two |

Primary halo comparisons use a **common 2048^3 analysis mesh**, with identical
physical mass cuts and configuration, so changing the force mesh does not
simultaneously change the finder mesh. Selected additional replays can measure
finder-mesh sensitivity separately. Native Ng/Np=8 finder meshes (A and D)
can generate many discreteness seeds; their cost requires a dedicated pilot.
Legacy means `a8c7715`, not the intermediate v2 normalization.

The normal exact-redshift schedule has 158 steps; T has 316, with identical
normal endpoints and output epochs. Simply halving the adaptive input step
would instead yield 312 steps with a different schedule. T starts with its own
correctly staggered velocities. Output velocities remain at a_out-da/2 in the
native PM convention, so timestep differences include this output staggering.

The completed comparisons are A/C/E (particle mass), B/C/D and E/F (evolved
force resolution at common finder mesh), and F/T (time). The earlier
L512/N1024/Ng2048 realization is a separate volume comparison, not a matched
same-object pair. Halo abundance, matched masses/radii/Vmax/velocity/shape,
completeness and independent membership invariants determine usable
mass/redshift ranges. Agreement must be measured above particle-count floors
and with finite-volume/counting uncertainty.

## Resources and storage

The prior N1024/Ng2048 job 11948372 used 64 shared cores, 3477 s,
61.81 billed core-hours, 52.10 CPU hours and 166.62 GiB batch MaxRSS.
The initial seven-run estimate was 1600-2800 billed core-hours including
replays and pilots; the planning allowance is 3000. This is not measured
4096-grid performance. New pilot accounting determines final allocations.
Both expected and wall-limit core-hours are recorded before submissions.

Measured four-step N1024/Ng4096 pilots used 720.06 s at 64 cores and
408.81 s at 128 cores inside the simulation. Doubling the cores was 1.76
times faster but used 13.5% more simulation core-hours (22.6% more when
including the fixed job overhead). All three full 4096-mesh runs therefore
use 64 shared cores. Their 384-GiB requests cover the measured 289.7-GiB
Slurm MaxRSS with headroom; two such tasks fit concurrently on one 1-TB,
128-core node. Expected runtimes are 480, 500 and 1000 minutes for D, F and
T; wall limits are 720, 720 and 1440 minutes. These are projections from
early steps, with margin for late-time density work and snapshot I/O.

Use `cosma8-serial`, actual thread counts and explicit memory. A 4096^3
float32 field alone occupies 256 GiB; a saved-density double-mesh replay
would need over 560 GiB before finder workspace at N1024. The primary replay
adapter avoids that extra full mesh. Full-node use requires a measured reason.

Starting user quota: 241104 files / 400000 soft / 440000 hard; space
23472442148 KiB / 32212254720 KiB soft. Keep one frozen build and reference
input set, small receipts, and scientific snapshots. Verify archives before
deleting owned build/fixture scratch. Preserve unrelated working directories.
Figures, Beamer source and PDFs are local artifacts, excluded from Git.
The 21 common-analysis density tapes occupy about 672 GiB but only 21 files;
retain them for exact replay/resume. Inode cleanup targets verified compiler,
fixture and agent-worktree artifacts, not these scientific controls.

Completed cleanup consolidated 465 files/links into three verified archives:
23 root compiler outputs, 329 completed build/pilot fixtures, and 113 replay
driver fixtures. The three receipts record every original hash and each
archive hash. Executables, active data and queued-job bundles remain in place.
An incremental Git repack removed another 197 object files and 136 empty
object directories. The new pack was verified and every Git ref remained
identical; `git-cleanup.json` records the before/after inventory. No history
or unreachable objects were expired. Six temporary slide/plot inspection
images were removed after the PDFs passed visual review.
On 8 September, cleanup job 11956691 consolidated another 92 completed-job
scripts, bundles and logs into `launches-finished-launches.tar.gz`, a net
reduction of 90 files after its archive and receipt. Every archive member was
read back and hashed before any original was removed. Scientific products,
executables, configurations and active/unknown launch files remain unpacked.
The job used 6 s, 1.191 CPU s and 42,620 KiB Slurm batch MaxRSS. Final cleanup
job 11956701 completed after render job 11956208, consolidating another 27
launch files into `launches-post-render-launches.tar.gz` after checking complete
analysis and matching PDF/manifest hashes. It used 10 s and 55,552 KiB MaxRSS
on one shared core with 256 MiB requested. Across the five scratch/launch
archives, 584 files/links were consolidated, a net reduction of 574 files
after archives and receipts. Git object packing is counted separately.
After the review, `launches-review-20260909.tar.gz` consolidates six more
completed launch files (the prior final-cleanup job and the IC recheck), with
member-by-member verification before removal. The six archives now hold 590
files/links, a net reduction of 578 after archives and receipts. Restore this
additional batch with `tar -xzf launches-review-20260909.tar.gz -C .` here.
`launch-cleanup-preflight.json` records controls for active/unknown jobs,
changed originals, corrupt archives and interrupted-deletion recovery.
The workspace inventory counted about 12,564 files/directories before this
pass; the account quota also includes files outside this repository. Existing
catalogue products account for 3093 entries under the strict-v1 products tree
and are retained as scientific results.
To restore old build or fixture paths for rerunning controls, extract the
corresponding `work-artifacts.tar.gz` or `work-controls.tar.gz` in this directory.
Restore launch paths with `tar -xzf launches-finished-launches.tar.gz -C .`
here, or use `launches-post-render-launches.tar.gz` for the later jobs. Archive
receipts retain original paths and hashes even while those paths are packed.
No rebuild is needed to execute the frozen campaign.

## Analysis and continuation

`jobs.json` records every job, frozen bundle, dependencies, expected and
time-limit core-hours, plus any resource amendment. Paired B replays measured
94.5 GiB batch MaxRSS (native finder processes peaked at 35.2 GiB); pending C/D
paired replays were increased from 96 to 128 GiB through recorded `scontrol`
updates. Completed C subsequently measured 122.15 GiB, so still-pending D
was increased to 160 GiB. The E z=0 v3 pilot completed in 544 s including
validation (457.77 s inside the native finder), with 122.90 GiB batch MaxRSS.
E/F/T pairs used 192 GiB at 32 cores. F's completed batch peak reached
191.99 GiB, and full analysis reached 7.83 GiB of its 8-GiB request. For future
equivalent full reruns, allow at least 256 GiB for high-resolution paired
replays and 16 GiB for analysis, using those completed peaks plus headroom.
For an E-equivalent 1024^3-particle/2048^3-mesh evolution, request at least
128 GiB: the original 96-GiB allocation reached 95.99 GiB MaxRSS. The original
IC-check job reached 11.99 GiB; its completed review rerun reports 15.99 GiB
against a 16-GiB request. Allow 24 GiB for an equivalent future check until
process RSS and charged file cache have been measured separately.
The original submitted scripts and allocation amendments remain recorded.

`analyze.py` reads independently validated v3 memberships. It matches shared
initial-lattice IDs, records both overlap fractions, and compares bound masses
as float64 member counts times stored particle masses. Reported total aperture
mass and radius are labelled separately from the SO definition. A zero Vmax
is an unresolved sentinel, excluded from shifts and counted separately.
Entire mass bins must lie above the publication/common particle-count floor.
The analysis records 100/300/1000-particle floors, paired eight-octant abundance
uncertainty, matched scatter, completeness and source/input hashes.

`assess.py` writes `CONVERGENCE.md` and `convergence-assessment.json`. Its
adjustable working screen requires 30 objects per bin, abundance and median
bound mass within 5%, resolved median Vmax within 2%, reference completeness
of at least 90%, and abundance jackknife sigma no greater than 5%. These are
descriptive comparisons; they do not certify absolute physical accuracy.

`render.py` runs frozen plotting and local Beamer sources after analysis.
One multipage vector figure PDF contains the main and additional diagnostics;
the Beamer deck shows the design, membership checks and measured ranges.
Final visual review corrected a plotting-only defect: setting the lower y
limit before adding data had frozen the absolute velocity/centre-difference
panels at an upper limit of 1. Those panels now scale to the plotted data.
`plot-controls.json` checks coverage of every displayed 84th percentile on
all 18 affected axes; measured statistics and all 63 assessments are unchanged.
Partially completed measurements are labelled explicitly. Both plot and
presentation receipts distinguish automated checks from visual review.
Presentation source, themes, PDF, figures and launch bundles are ignored by Git.
The full render also refreshes Slurm accounting; its own few remaining seconds
are not yet included at that timestamp. `visual-review.json` binds the inspected
final PDFs by hash; a subsequent render would require a new visual review.

Run every Python command through `micromamba run -n cosemu python3 -B`.
Receipts and the final report distinguish preparation, successful process
completion, scientific validation and actual convergence; a submitted job
does not count as completed validation.
