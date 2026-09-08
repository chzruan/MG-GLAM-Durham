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

## Launch status, 8 September 2026, 01:15 UTC

All seven initial conditions are complete and their shared-mode checks passed,
including the full 1024^3-row comparison. Simulations A, B, C and E completed
all three snapshots. D, F and T are running. Independent v3 catalogue checks
have passed for every completed replay; this is not yet a completed convergence
result. The preliminary report and slides currently use eight catalogues and
five pair/redshift comparisons, frozen before later replays arrived.

| Remaining stage | Slurm job | Dependency |
|---|---:|---|
| D simulation → paired replay | 11955919 → 11955978 | Replay waits for simulation |
| F simulation → paired replay | 11955920 → 11956205 | Replay waits for simulation |
| T simulation → paired replay | 11955921 → 11956206 | Replay waits for simulation |
| A paired replay | 11956146 | Running; legacy is slower on this coarse particle load |
| E paired replay | 11956204 | Running; retains the completed z=0 v3 pilot |
| Full convergence analysis | 11956207 | All five remaining replay jobs above |
| Final plots and Beamer | 11956208 | Full analysis |
| Final launch-file cleanup | 11956701 | Successful full plots and Beamer |

Cleanup-time update (05:13 UTC): A/E paired replays have also completed at all
three redshifts. Twelve of 21 independent v3 catalogues are now published; the
saved preliminary analysis still uses its original eight inputs.

Controller dependencies, requested cores/memory, frozen script/bundle hashes,
shell syntax and bundled Python syntax have been checked. B/C replays already
completed and are verified through their receipts. Full analysis refuses
missing v3 catalogues; legacy failures or its 900-s per-stage time cap are
reported independently and do not discard validated v3 measurements.

At 01:14 UTC, accounting recorded 339 billed core-hours so far. Replacing
finished-job estimates with actual usage gives about 2400 core-hours for the
campaign. Summed submitted time limits allow 3889 core-hours; that ceiling is
not the expected charge. T remains the critical path, projected at roughly
16-17 hours of evolution plus replay and analysis; queue delays are additional.

On continuation, inspect jobs 11956207/11956208 and `accounting.json`, confirm
`convergence.json` has `completed: true` and 21 inputs, and inspect the final
plots/slides before drawing conclusions. Record a new PDF-hash-bound visual
review and commit the completed measurement receipts on this branch. The
queued jobs do not commit or merge automatically. If an upstream job fails,
inspect its frozen log and receipt before resuming; preserve successful stages.
Completed launch files may have been consolidated as described below; restore
the relevant archive when the original script, bundle or log path is needed.

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
position difference of 6.103515625e-5 Mpc/h and identical physical velocities;
it did not record edge-clamp incidence or RMS differences. The conservative
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

The planned comparisons are A/C/E (particle mass), B/C/D and E/F (evolved
force resolution at common finder mesh), and F/T (time). The earlier
L512/N1024/Ng2048 realization is a separate volume comparison, not a matched
same-object pair. Halo abundance, matched masses/radii/Vmax/velocity/shape,
completeness and independent membership invariants will determine usable
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
job 11956701 is queued after render job 11956208 on one shared core and
256 MiB (5-minute limit; 15-second expectation). It requires complete analysis
and matching PDF/manifest hashes, then archives newly finished launch files.
`launch-cleanup-preflight.json` records controls for active/unknown jobs,
changed originals, corrupt archives and interrupted-deletion recovery.
The workspace inventory counted about 12,564 files/directories before this
pass; the account quota also includes files outside this repository. Existing
catalogue products account for 3093 entries under the strict-v1 products tree
and are retained as scientific results.
To restore old build or fixture paths for rerunning controls, extract the
corresponding `work-artifacts.tar.gz` or `work-controls.tar.gz` in this directory.
Restore launch paths with `tar -xzf launches-finished-launches.tar.gz -C .`
here, or use the corresponding post-render archive once it exists. Archive
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
Remaining E/F/T pairs use 192 GiB at 32 cores, allowing additional headroom
for the finer-force density fields. The original submitted scripts remain
immutable. The final
expected campaign cost is approximately 2400 billed core-hours, with the
3000-hour planning allowance retained. Wall limits provide additional margin
and their summed maximum is larger than expected usage.

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
Partially completed measurements are labelled explicitly. Both plot and
presentation receipts distinguish automated checks from visual review.
Presentation source, themes, PDF, figures and launch bundles are ignored by Git.
The full render also refreshes Slurm accounting; its own few remaining seconds
are not yet included at that timestamp. `visual-review.json` binds the inspected
preliminary PDFs by hash and does not certify later automatic renders.

Run every Python command through `micromamba run -n cosemu python3 -B`.
Receipts and the final report distinguish preparation, successful process
completion, scientific validation and actual convergence; a submitted job
does not count as completed validation.
