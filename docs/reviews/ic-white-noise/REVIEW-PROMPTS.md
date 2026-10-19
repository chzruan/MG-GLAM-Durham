# Review prompts for branch `ic-white-noise-20261019` (fix commit e604768, `cz..ic-white-noise-20261019`)

Three reviews before the merge into `cz`, carried out in one fresh session or
subagent, in order: the Fortran fix; the mechanism and the evidence; the note,
the record and the consequences for the existing runs. The COMMON PREAMBLE is
read once; each REVIEW block is then self-contained. Findings go to the
reviewer's own folder outside the repository
(`/cosma8/data/dp203/dc-ruan1/scratch/mgglam-ic-review-20261019/findings/REVIEW-<k>-<topic>-findings.md`);
the author copies them into `docs/reviews/ic-white-noise/` with a triage
before the merge. The reviewer builds and runs only in the two detached
worktrees prepared by the author (never in the main checkout, never in the
author's worktree).

These files are instructions. No review has been carried out by writing them.

---------------------------------------------------------------------------
## COMMON PREAMBLE (read first)

You are reviewing code on the COSMA HPC system (Durham). Repository:
`/cosma8/data/dp203/dc-ruan1/mgglam_claude/MG-GLAM` (MG-GLAM, the Fortran
particle-mesh N-body code of the DEGRACE simulations, LCDM and modified
gravity, plus tools). Its mainline is `cz`, NOT `main`; `cz` is checked out
there, and that checkout holds terabytes of untracked run directories. The
branch `ic-white-noise-20261019` (the fix commit e604768 on top of `cz` =
1024c1e, plus the commit that adds this prompts file) is checked out in the
author's worktree `/cosma8/data/dp203/dc-ruan1/scratch/mgglam-ic-white-noise-20261019`:
do not touch it. You took no part in the work. Review only: find what is
wrong, unsupported or missing, and say what is sound. Fix nothing, commit
nothing, merge nothing.

### Background

On 2026-10-19 a blind review of control-variate work in the halobias
repository measured that the white noise of the DEGRACE initial conditions
(L = 1024 Mpc/h, 2048^3 particles, the native `PMP2start` generator) is not
homogeneous: the Gaussian numbers are correlated at even lags 12-28 along the
innermost index i (up to 0.015, the same for two seeds), which modulates the
variance of the linear field along x at harmonics 6-12 of the box by 0.5-2.6
per cent, the same in every box; the squared fields carry coherent power at
the axis modes (+-j, 0, 0), j = 7-12; the simulations' halo fields show it
(+20 to +67 per cent at |n|^2 = 64); the raw halo correlation function of
every box is +0.1 to +0.6 per cent high at 30-40 Mpc/h. That review did not
locate the origin in the code. Its text is read-only:
`/cosma8/data/dp203/dc-ruan1/halobias/analysis/zcv_xihh_review_20261019_v1/implementation/REPORT.md`
(finding 1, lines 11-60) and `.../REVIEW.md` (finding 1, lines 24-60); its
measurements are in `.../implementation/results/test_whitenoise_lags.json`
and `test_variance_profile.json`. The branch under review is the simulation
side's answer.

### What the commit claims (every statement is a claim to test, not evidence)

The fix commit e604768 (`git -C <wt> show --stat e604768`: 32 files, 4721
insertions, 5 deletions) says:

1. Cause. `PMP2mod_random.f90` holds RANLUX (James 1994) as a Fortran
   module. The statement `DATA ndskip / 0, 24, 73, 199, 365 /` that sets the
   table of numbers to discard per luxury level is a comment (line 167 on
   `cz`, so since the first commit 51cb6ed), and nothing else sets `ndskip`
   of the common block `/LUX2/` (`luxuryp.h`), so the table was zero and
   `rluxgo(lux, ...)` set `nskip = 0` for every `lux`: every stream ran at
   luxury level 0 (plain subtract-with-borrow, RCARRY), whatever the four
   callers asked (`lux = 2` in `SetRandomN` of `PMP2start.f90` and
   `PMP2startMPI.f90`, `getRandom` of `PMP2mod_random.f90`, `getGauss` of
   `PMP2mod_density.f90`).
2. Mechanism. With nothing discarded, the recurrence
   x(n) = x(n-10) - x(n-24) - carry (mod 1) ties Box-Muller pair p+12 to
   pairs p and p+7 (pair p+12 = pair p+7 - pair p, mod 1, component by
   component; holds for 0.99994 of the pairs at level 0, 2-3e-7 at levels
   2-3). The rejection step of `GAUSS3` (polar Box-Muller, 21.5 per cent of
   the pairs rejected) makes the acceptance of pair p+7 depend on pairs p and
   p+12 jointly, so the delivered Gaussian numbers are correlated at lags
   that depend on the acceptances in between: even lags 12-28,
   anti-correlated at 12-18, correlated at 20-28, the same for every seed,
   along i because that is the direction of the stream (one stream per plane
   k; `SPECTR` draws row j, column i, i innermost). Conditional products:
   x(p) x(p+12) given pair p+7 rejected / accepted = -0.458 / +0.125; the
   Gaussians g(p) g(p+12) = -0.548 / +0.170.
3. Fix. `rluxgo` sets `ndskip = 0, 24, 73, 199, 365` on every call; a new
   function `RanluxLevel()` returns the level for every stream (environment
   variable `GLAM_RANLUX_LUXURY`: 0-4, or 24-2000 as p itself; default 3;
   printed once); the four callers use it instead of `lux = 2`; `SetRandomN`
   (both files) and `getGauss` set `iFlag = 0` before the first `GAUSS3` (it
   was an undefined local; the note infers it behaved as 0 in the production
   builds). `GLAM_RANLUX_LUXURY=0` reproduces bit for bit the draws of every
   run made before the change.
4. Verification (`test/ic_white_noise/`: `ranlux_lags.f90` built against the
   repository's own `PMP2mod_random.f90` and `PMP2mod_fft5.f90`, `makefile`,
   `analyse_cubes.py`, `results/`). Against the unmodified module, lux 0, 2
   and 3 give the same numbers (`results/orig_*`). Against the fixed module,
   level 0 gives the same first eight values, row sum and lag table as the
   unmodified module (`fix_s1_lux0`); level 0 reproduces the halobias
   review's lag table to 3e-4; levels 1-4 give |corr| <= 4e-4 along i (noise
   1-2e-4), 1-D variance-profile harmonics <= 0.04 per cent (noise 0.02), no
   axis-mode spikes of the squared 256^3 field (ratios 0.4-2 against 40-1600
   at level 0); the profile of the two seeds correlates 0.975 at level 0 and
   -0.02 at level 3. `PMP2init`, `PMP2start`, `PMP2main` build from the
   patched sources; `PMP2startMPI` was not built (no MPI module loaded). The
   switch was exercised (unset: 3; 0: the old numbers; invalid: 3 and a
   message).
5. Cost. Level 3 generates 9.3 times more numbers than level 0; 7 s to 20 s
   for 134 million Gaussians on one core (with the statistics); about 15
   minutes of one core for a 2048^3 run, shared by the OpenMP threads of
   `SPECTR`.
6. Consequences. Every run with native initial conditions (`PMP2start`,
   `PMP2startMPI`; the DEGRACE L = 1024 and L = 1800 suites) carries the
   pattern; runs from 2LPTIC do not. The mean and the power spectrum of the
   linear field are right; the pattern sits in the fourth-order statistics.
   Bin-averaged power spectra, mass functions and the correlation function
   below about 20 Mpc/h are affected well below 1 per cent; axis modes,
   quadratic fields and the correlation function beyond 30 Mpc/h (0.1-0.6
   per cent) carry it, and it does not average out over realisations. The
   galaxy biasing (`getRandom`) and `getGauss` draws change at the new
   default. Not established: the L = 1800 runs (not measured), the evolved
   matter field, `iFlag` = 0 in the production binaries (inferred), a closed
   form of the lag correlations.
7. Texts. `docs/ic-white-noise-20261019.md` (the note: summary table,
   mechanism, fix, verification tables, consequences, what is not
   established, record); a new bullet in `AGENTS.md` ("Default initial
   conditions"); a comment in `luxuryp.h`; the commit message.

### Where things are

| item | path |
|---|---|
| your worktree of the branch (detached at its HEAD, writable, yours) | `/cosma8/data/dp203/dc-ruan1/scratch/mgglam-ic-review-20261019/wt` |
| your worktree of the unmodified code (detached at `cz` = 1024c1e, writable, yours) | `/cosma8/data/dp203/dc-ruan1/scratch/mgglam-ic-review-20261019/cz` |
| your findings and scratch | `/cosma8/data/dp203/dc-ruan1/scratch/mgglam-ic-review-20261019/findings/` and `.../work/` |
| the author's worktree (do not touch) | `/cosma8/data/dp203/dc-ruan1/scratch/mgglam-ic-white-noise-20261019` |
| the main checkout (`cz`; read with pathspecs only; change nothing) | `/cosma8/data/dp203/dc-ruan1/mgglam_claude/MG-GLAM` |
| the halobias review that measured the pattern (read-only) | `/cosma8/data/dp203/dc-ruan1/halobias/analysis/zcv_xihh_review_20261019_v1/` (`REVIEW.md`, `implementation/REPORT.md`, `science/REPORT.md`, `implementation/results/`) |
| the halobias replica of the draw (read-only) | `/cosma8/data/dp203/dc-ruan1/abcs_haloes/degrace_ic/src/lowk_ts.f90`, `glamic.py`; the built `.../degrace_ic/build/fast/lowk_ts.exe` |
| DEGRACE run folders (read-only: seed tables, launch scripts, job output) | `/cosma8/data/dp203/dc-ruan1/mg_glam/DurMun_hmfemu_LCDM_wide_sample_first_64_model1_L1024Np2048Ng4096/` (`TableSeeds.dat`, `Setup.dat`, `Run1/`); `/cosma8/data/dp203/dc-ruan1/mg_glam/F5n1_L1024Np2048Ng4096/TableSeeds.dat` (box 2 seed 2009256176); `/cosma8/data/dp203/dc-ruan1/mgglam_claude/MG-GLAM/DEGRACE_LCDM_L1800Np3200Ng6400/arun.sh` (line 26 calls `PMP2start.exe`) |
| the small example (128^3 particles, 256 mesh, 750 Mpc/h; seeds 1156811761 for box 1, 2009256176 for box 2) | `test/` of each worktree: `Setup.dat`, `TableSeeds.dat`, `PkTable.dat`, `README.md` ("Build", "Run"), `arun.sh` |

### Hard rules (shared system; the user is personally accountable for what you do)

- Change nothing in the main checkout, the author's worktree, the run
  folders, the halobias repository or `abcs_haloes`. Commit nothing
  anywhere. No git state changes (no checkout, reset, stash, branch,
  worktree add or remove, clean) except `git checkout -- <file>` in YOUR
  worktrees to restore a file you changed for a mutation check. Build and
  run only in your two worktrees and `work/`. At the end
  `git -C <your worktree> diff --stat` must be empty for both, and you say
  which built files you left. Do not start further agents, background
  processes, servers, watchers or tmux sessions. Never touch jobs, processes
  or tmux sessions you did not start (the user's own are running).
- Compute: login node, light tier only: one process at a time under
  `nice -n 19`, at most 4 threads (`OMP_NUM_THREADS=4`,
  `OMP_STACKSIZE=512M`), each step about a minute and under 2 GB. Measured
  by the author: builds 20 s serial; `ranlux_lags.exe` mode 1 at
  2048 x 256 x 256 1-20 s and 30-140 MB, mode 2 at 20 million pairs 2-5 s
  and 310 MB; `analyse_cubes.py` 53 s and 840 MB; `PMP2start.exe` on the
  128^3 example a few seconds. Check `uptime` (load well below 128) and
  `free -g` before every start. Builds with `make -j1`. No Slurm job (if you
  think one is needed, say so in the report instead). Nothing in `/tmp`:
  cubes and PM files go under `work/` and are deleted at the end (a 256^3
  cube is 67 MB; the example's PM files about 50 MB per run).
- Compiler: `module load intel_comp/2024.2.0 compiler-rt tbb compiler`
  (ifx 2024.2.0, the production compiler; in a login shell, `bash -lc`, or
  after `source /etc/profile.d/modules.sh`). For the MPI compile-only
  check: `export I_MPI_F90=ifx; module load mpi` after the above (the
  example's `arun.sh` does this); if it is not available, say so.
- Python only as
  `nice -n 19 /cosma/apps/durham/dc-ruan1/micromamba/envs/cosemu/bin/python3 -I ...`
  (numpy available). Never `micromamba run` to launch a Fortran executable
  (the environment shadows a runtime library and the executable dies at
  start-up); never bare `python`, `python3` or `pip`; install nothing.
- No recursive grep, find, du or `ls -R` over `/cosma*`, `/snap*`, the home
  directory or `/cosma8/data/dp203/dc-ruan1`. In the main checkout use
  `git grep -n <pattern> -- '*.f90' '*.h'` and `git status` or `git diff`
  ONLY with a pathspec; in your worktrees `git status` is fine. `grep` may
  be a shell function here: use `/usr/bin/grep`. Do not read or touch
  `~/.ssh`.
- Read the diff and the code around it:
  `git -C <wt> diff 1024c1e e604768 -- <paths>`, `sed -n 'a,bp' <file>`.
  Line numbers in your findings refer to the branch; say so when they refer
  to `cz`.
- Wording: plain astronomy and everyday words (the user is a cosmologist);
  code identifiers verbatim; numbers in tables; no software or
  machine-learning shorthand (not "gate", "artifact", "harness", "flag" as a
  verb, "registry", "pin", "payload", "sidecar"). State facts; do not
  describe your reasoning process.

### Decisions that belong to the user, not to you

Report their consequences and the cost of the remedies; do not count them
as defects: that the fix lives on a branch of this fork rather than
upstream; the default level 3 (the RANLUX authors' recommendation) rather
than 2 (what the callers had asked for); an environment variable as the
switch rather than a line of `Setup.dat`; whether anything is re-run;
whether and how the MG-GLAM authors are told; the merge itself.

### Deliverable, per review block

Write `findings/REVIEW-<k>-<topic>-findings.md` with these sections.
(1) Summary, 3-6 lines, ending with a merge recommendation: ready / ready
after the listed fixes / not ready. (2) Findings table with the columns ID,
Severity, File:line, Description, Evidence, Suggested fix. Severity:
Blocker = wrong numbers or a crash in ordinary use, or the old draws NOT
reproduced at level 0; Major = wrong results in some configurations or
programs, a false central claim of the note, a silent change of a product,
resource misuse on the shared node; Minor = a limited defect, a gap in the
evidence, robustness; Nit. Keep apart what the commit introduced, what
existed before, what the commit itself discloses as not established, and
the user's decisions. (3) Verified OK: the specific properties you checked
and found correct, each with how (command and result), so nobody re-reviews
them. (4) Open questions for the author. (5) The runs you made: command,
wall time, peak memory, what was deleted and what was left. Cite file:line,
quote the expression, give a failing input where you can. Assume defects
exist and look for them; do not pad with generic advice.

---------------------------------------------------------------------------
## REVIEW 1: the Fortran fix (`PMP2mod_random.f90`, its callers, `luxuryp.h`, the builds)

Findings file: `findings/REVIEW-1-fortran-fix-findings.md`

Scope: `git -C <wt> diff 1024c1e e604768 -- PMP2mod_random.f90
PMP2mod_density.f90 PMP2start.f90 PMP2startMPI.f90 luxuryp.h`, read line by
line with the unchanged code around it. On the branch: `PMP2mod_random.f90`
whole (379 lines: `GAUSS3` 59-90, `getRandom` 114-125, `RanluxLevel`
136-157, `ranlux` 188-261, `rluxgo` 265-378); `luxuryp.h` (the common
blocks: `/LUX1/` is `threadprivate`, `/LUX2/` and `/AllSeeds/` are shared);
`PMP2start.f90` `SetRandomN` 75-100, the seed from the table 140-160 and
205-240, `SPECTR` 266-366 (the parallel region at 302 restores one stream
per plane); `PMP2startMPI.f90` `SetRandomN` 76-101 and `SPECTR` 319-419;
`PMP2mod_density.f90` `getGauss` 388-403 and its caller 459-471 (inside
`!$OMP critical` at 468); `PMP2mod_analyze.f90` 232-262 (the caller of
`getRandom`); `PMP2correlation.f90` 1028-1090 (another caller of `rluxgo`
and `ranlux` with `lux = 2` at line 1032, NOT changed by the branch);
`makefile`.

Check at least:

- The cause as stated. On `cz`, is the commented `DATA ndskip` the only
  place the table could be set; is `ndskip` zero at every `rluxgo` call in
  every program (static storage of the common block; confirm by printing
  `nskip` after `rluxgo(2, ...)` against the unmodified module); did every
  caller ask for 2?
- Every caller. List every use of `rluxgo`, `ranlux`, `GAUSS3`, `RANDd`
  and `RanluxLevel` in the tree
  (`git grep -n -i 'rluxgo\|call ranlux\|GAUSS3\|RanluxLevel' -- '*.f90'`).
  Does the branch cover every program the makefile builds (`all`:
  `PMP2init PMP2start PMP2main PMP2MG ic2pm`; also `PMP2analyze`,
  `PMP2spectrum`, `PMP2BDM`, `PMP2correlation`, `PMP2density`,
  `PMP2startMPI`, `PMP2mainMPI`)? `PMP2correlation.f90` line 1032-1034 sets
  `lux = 2` and calls `rluxgo`: what does it draw, is it affected by the
  table, should it use `RanluxLevel()`, and is the omission a defect or out
  of scope (is that program built and used)?
- The table set in `rluxgo` on every call, into the SHARED common block
  `/LUX2/`. Where is `rluxgo` called from inside an OpenMP parallel region
  (the `!$OMP critical` around `getGauss`; `getRandom`'s loop in
  `PMP2mod_analyze.f90`; `SetRandomN` is serial)? Are the writes to
  `/LUX2/` (`ndskip`, `nskip`, `luxlev`, `next`, `notyet`, `twom24`,
  `twom12`, `inseed`) a new race or the same situation as before the
  branch? Is `nskip` the same for every stream of a program, so that the
  shared value is right for every plane?
- The skipping and the stream save and restore. With `nskip > 0`, `ranlux`
  discards `nskip` numbers every time `in24` reaches 24. `SetRandomN` saves
  `i24, j24, in24, kount, carry, seeds` per plane right after `rluxgo`;
  `SPECTR` restores them per plane inside the parallel loop and sets
  `mkount = 0`. Is everything the skipping needs saved and restored, so that
  a plane's stream is the same whatever thread draws it and whatever the
  order of the planes? (The end-to-end run below at 1 and 4 threads is the
  direct test.)
- `RanluxLevel()`: the `save` variable, the critical section, the first call
  from inside a parallel region, the read of `levelSaved` outside the
  critical; the parsing of `GLAM_RANLUX_LUXURY` with `Read (str, *)`: what
  happens for `3 `, `03`, `3.5`, `-1`, `5`, `23`, `2001`, an empty string, a
  non-numeric string, `0x3`; is the printed message right in each case; is
  the accepted range the same as `rluxgo`'s (`lux` 24-2000 sets p directly;
  `rluxgo` maps p = 24, 48, 97, 223, 389 back to the levels)? Run the
  built `ranlux_lags.exe` (REVIEW 2's build; mode 1 with `lux = -1` uses
  the function) with a few of these values.
- `iFlag`. On `cz`, `SetRandomN` stores a local `iFlag` it never set. Where
  is it declared (implicit typing through `luxuryp.h`; a local or a module
  variable), and is it undefined there? The note infers it behaved as 0 in
  the production builds. Test it end to end (below): the `cz` build and the
  fixed build at `GLAM_RANLUX_LUXURY=0` must give bit-identical initial
  conditions; if they do not, say whether the difference is the undefined
  `iFlag` or the fix itself.
- End-to-end reproduction of the old draws, the central claim of the
  switch. In each of your worktrees build serially
  (`module load ...; make -j1 PMP2init PMP2start`, about 20 s; the same in
  the `cz` worktree). Then, in a copy of `test/` under `work/` (one copy
  per build and setting; the executables run from `Run<box>/` and read
  `../Setup.dat`, `../TableSeeds.dat`, `../PkTable.dat`; `README.md`
  "Run"): `cd Run2 && OMP_NUM_THREADS=4 nice -n 19 <build>/PMP2start.exe <<< 2`
  (seconds; writes the PM files `PMcrd.DAT` and `PMcrs0.DAT` of 128^3
  particles, about 50 MB; check what else it writes). Compare with `cmp`:
  (a) `cz` build against the fixed build with `GLAM_RANLUX_LUXURY=0`: must
  be identical; (b) fixed build at the default against
  `GLAM_RANLUX_LUXURY=0`: must differ; (c) fixed build at the default with
  1 thread against 4 threads: must be identical (a plane's stream does not
  depend on the thread); (d) `GLAM_RANLUX_LUXURY=3` against unset:
  identical; (e) the printed level line in each case, and whether any
  output file of the run records the level or the seed. Report wall time
  and peak memory (`/usr/bin/time -v`) for levels 0 and 3, the only cost
  measurement asked for (do not run anything at 2048^3). Delete the PM
  files when done.
- The MPI generator. `PMP2startMPI.f90` has the same two-line change and was
  not built. Try a compile-only check in your branch worktree after the
  module objects exist (`export I_MPI_F90=ifx; module load mpi; mpiifort
  <MPIFLAGS of the makefile> -c PMP2startMPI.f90`); if the module is not
  available, say so. Does `PMP2startMPI` print the level once per rank, and
  does every rank see the environment variable under `mpirun`?
- Builds. Do `PMP2start`, `PMP2main`, `PMP2MG`, `ic2pm`, `PMP2analyze` and
  `PMP2spectrum` compile and link with the changed module (`make -j1
  <target>`; the makefile's `.f90.o` rule uses `-w`, so look for real
  errors only)? Compile `PMP2mod_random.f90` once with `-warn all` into
  `work/` and report any warning about the new code.
- `luxuryp.h` is unchanged except a comment: are `ndskip(0:maxlev)`,
  `maxlev = 4` and `lxdflt = 3` consistent with the function's default, the
  table and the note?
- Cost. Is the extra generation time at level 3 in `SPECTR` (parallel over
  planes) and in `getGauss` and `getRandom` (serial or inside a critical
  section in the analysis programs) what the note says (9.3 times the
  numbers generated; "about 15 minutes of one core" for 2048^3)? Scale
  from your example measurement and from REVIEW 2's test timings.

---------------------------------------------------------------------------
## REVIEW 2: the mechanism and the evidence (`test/ic_white_noise/`, the note's tables)

Findings file: `findings/REVIEW-2-mechanism-evidence-findings.md`

Scope: `test/ic_white_noise/ranlux_lags.f90` (295 lines), `makefile`,
`analyse_cubes.py`, `.gitignore`, `results/` (22 files); the sections
"Mechanism", "Reproduction and verification" and the summary table of
`docs/ic-white-noise-20261019.md`; the halobias review's measurements
(`implementation/REPORT.md` lines 11-60, `results/test_whitenoise_lags.json`,
`test_variance_profile.json`) and its replica `lowk_ts.f90`.

Check at least:

- Faithfulness of the test. Does mode 1 of `ranlux_lags.f90` reproduce
  `SetRandomN` and `SPECTR` (seeding `rluxgo(lux, Ns, 0, 0)` per plane with
  `Ns` advanced by `RANDd`; `iFlagA(k) = 0` where `SetRandomN` stores an
  undefined local; draw order j outer, i inner, one `GAUSS3` per cell; the
  restore of the stream per plane)? Does its lag estimator do what its
  output claims (the lag correlation along i over all columns; along j and
  k on the kept cube; the variance profile along x from `rfft1b` of each
  row, the basis of `VECTOR`)? Is the `.vprof` profile the variance of the
  white field along x (one row of Gaussian coefficients synthesised with
  `rfft1b` is a 1-D field along x; its square averaged over rows and planes
  is the variance profile)? Which index of GLAM's coefficient cube is x in
  real space (`PMP2start.f90` `VECTOR`, `map3`), and is "along i means
  along x" right?
- Rebuild and re-run. In your branch worktree
  (`cd test/ic_white_noise && make`, then the runs listed in the note's
  "Reproduction and verification", each seconds to 20 s under `nice -n 19`,
  outputs under `work/`): `fix_s1_lux0`, `fix_s1_lux3`, `fix_s2_lux0`, the
  mode-2 run `fix_pairs_s1_lux0`, and against the unmodified module
  (`make GLAM=/cosma8/data/dp203/dc-ruan1/scratch/mgglam-ic-review-20261019/cz B=build_cz`,
  then `orig_s1_lux0`, `orig_s1_lux2`, `orig_s1_lux3`). Compare with
  `results/` line by line (bit-identical expected for the same inputs: say
  whether they are). Compare the level-0 lag table with the halobias
  review's (`test_whitenoise_lags.json`) and state the largest difference.
- Recompute the mechanism with your own code from a mode-2 run at level 0
  (the test writes no raw pairs: add a small program of your own in `work/`
  that calls the module's `rluxgo` and `ranlux`, built against the
  unmodified module): the three-term relation between pairs p, p+7, p+12
  (positions n, n+14, n+24 of the stream; from
  x(n) = x(n-10) - x(n-24) - carry, is "pair p+12 = pair p+7 - pair p" right
  for both components, with the carry and the padding of small numbers in
  `ranlux` lines 219-224 accounted for?), the conditional products given
  the acceptance of pairs p+6, p+7, p+8, and the Gaussian lag correlations.
  Are the signs and "even lags only" explained by the relation (the same
  component of the two pairs), and is the lag range 12-28 what the
  acceptance statistics predict (12 pairs x 0.785 x 2 = 19 +- 3)? Say
  whether the mechanism is established or only consistent with the numbers.
- The field-level claims (`analyse_cubes.py`): the synthesis
  f = sum TS(i,j,k) B(x,i) B(y,j) B(z,k) with the basis the test writes; the
  variance-profile harmonics along x, y, z; the squared-field power at the
  axis modes over the neighbours. Re-run it on your cubes (`cube256_s1_lux0`,
  `cube256_s2_lux0`, `cube256_s1_lux3`; 53 s, 840 MB) and compare with
  `results/analyse_cubes.json`. Is "over the mean power of the non-axis
  modes with |n| within 0.5 of m" a fair reference (two axis modes against
  about 4 pi m^2 modes on the shell)? Is every cell of the note's three
  tables the value in `results/`?
- Noise levels. Are the quoted noise levels (lag correlation 1-2e-4; profile
  0.02 per cent; cube 0.06-0.14 per cent) right for the sample sizes (134
  million values; 256 x 256 rows of 2048; 256^3), and are the level 1-4
  numbers within them? Is 0.0149 at lag 22 with noise 1e-4 "the same for
  both seeds" to the precision claimed (3e-4)?
- Sizes. The test draws 256 planes x 256 rows x 2048 columns; the runs are
  2048^3. Does anything in the mechanism depend on the row length or the
  number of rows (the stream runs across rows within a plane), so that
  2048-column rows are the right test and a 256^3 cube a fair field-level
  stand-in? Does the lag pattern repeat identically in every row and plane
  (it should, if it is a property of the recurrence)? Test it: the lag
  table from rows 1-64 against rows 193-256.
- Seeds. 1156811761 and 2009256176 are the seeds of fiducial box 1 and box
  2 (`TableSeeds.dat` of model 1 and of F5). Does the test's seeding give
  the same per-plane seeds as a production run (how `PMP2start.f90` sets
  `Nseed` from the table row of the realisation, lines 140-160 and
  205-240)?
- Was anything the note calls verified not actually run? Every `results/`
  file against the commands in the note (file names, header lines,
  `lux_requested`, `luxlev_set`, `nskip_used`). Are the results reproducible
  from the committed sources alone (the cubes were deleted; `analyse_cubes.json`
  depends on them: can it be rebuilt, and did you rebuild it)?

---------------------------------------------------------------------------
## REVIEW 3: the note, the record and the consequences for the runs

Findings file: `findings/REVIEW-3-note-consequences-findings.md`

Scope: `docs/ic-white-noise-20261019.md` (263 lines), the `AGENTS.md`
bullet, the `luxuryp.h` comment, the commit message of e604768; the
sections "Consequences for the existing runs", "Not established" and
"Record"; the halobias review (read-only) for the measured consequences; the
launch scripts and seed tables of the run folders (read-only) for which
suites used the native generator; the halobias replica `lowk_ts.f90` and
`glamic.py` as a downstream user of the module.

Check at least:

- Every number in the note against `results/` or the halobias review (a
  number copied is not checked: give a table of note value, source value,
  match). Every statement of the summary table.
- The consequences table. Which rows are measured (by the halobias review
  on the 95 fiducial L = 1024 boxes: the quadratic spectra at the six
  shells, the halo power at |n|^2 = 64, the raw xi_hh at 29.6-40.5 Mpc/h,
  the linear field) and which are inferred ("bin-averaged power spectra,
  mass functions and the correlation function below about 20 Mpc/h affected
  well below 1 per cent"; "the evolved matter field"; "the L = 1800
  suites"; "does not average out over realisations")? Is each inferred
  statement marked as such, and is the inference sound (how many modes per
  bin at the affected shells; the pattern's amplitude; what "the same in
  every box" implies for a mean over realisations)?
- Which runs are affected. Do the launch scripts of the DEGRACE suites call
  `PMP2start.exe` (`DEGRACE_LCDM_L1800Np3200Ng6400/arun.sh` line 26; for
  the L = 1024 training and fiducial runs read one launch script or
  `Run1/*.dump` of the model-1 folder named above), and which, if any, used
  `PMP2startMPI.exe` or 2LPTIC (`AGENTS.md` makes 2LPTIC the default for new
  runs; the `fid2LPTIC_*` folders of the main checkout)? Does the note name
  the affected suites correctly? Is "runs from 2LPTIC are not affected"
  right (does `ic2pm` or the 2LPTIC path call the module's generator:
  `git grep -n -i 'ranlux\|rluxgo' -- '2LPTIC_Gui' 'ic2pm.f90'`)?
- Downstream users of the module. The halobias replica `lowk_ts.f90`
  hard-codes `lux = 2` and takes `iFlag0` as input, and is built against
  "the repository's own `PMP2mod_random.f90` and `luxuryp.h`" (which copy?
  its build folder `abcs_haloes/degrace_ic/build/fast/`; its build script
  if there is one). After the merge, a rebuild of the replica against the
  fixed module would draw at level 2 and no longer reproduce the runs
  unless it sets `GLAM_RANLUX_LUXURY=0` or `lux = 0`. Does the note or
  `AGENTS.md` warn the downstream users? Is anything else in this
  repository built against the module and used to reproduce old draws
  (`git grep -n -i 'rluxgo' -- '*.f90' 'python' 'scripts' 'test'`)?
- The switch as a record. An environment variable is not written into
  `Setup.dat` or any output of the run. Does the code print the level
  anywhere that survives (the job's standard output only)? Would a user of
  an old run's reproduction know which level was used? Is the `AGENTS.md`
  instruction ("set `GLAM_RANLUX_LUXURY=0` and record it in the campaign's
  provenance") enough? ("Provenance" is the repository's own term in
  `AGENTS.md`: say so, do not count it.)
- The "Not established" list: is it complete? Add what you found missing
  (for example the end-to-end reproduction of a production initial
  condition at level 0, the MPI build, the evolved field, the halo
  catalogues, the galaxy-biasing draws, the on-the-fly products).
- The "Record" section: are the compute statements (login node, one core,
  `nice`, times, memory, no Slurm job) consistent with the results files
  and with the timings REVIEW 2 measured?
- Wording rule: list every software-shorthand word in the note, the
  `AGENTS.md` bullet and the commit message that the user's rule forbids
  (gate, artifact, harness, flag as a verb, registry, pin, payload,
  sidecar; "provenance" only outside the repository's own term); numbers in
  sentences that belong in tables.
- Is the commit message accurate, bullet by bullet against the diff? Is
  anything in the diff not described (the `luxuryp.h` comment, the
  `.gitignore`, the `results/` files)?

---------------------------------------------------------------------------
## KICK-OFF (paste into one fresh session or subagent; the three blocks in order)

```text
Review, do not fix, do not merge. You took no part in the work under review.

Read in full, in this order:
- /cosma8/data/dp203/dc-ruan1/scratch/mgglam-ic-review-20261019/wt/docs/reviews/ic-white-noise/REVIEW-PROMPTS.md (the COMMON PREAMBLE and the three REVIEW blocks)
- /cosma8/data/dp203/dc-ruan1/scratch/mgglam-ic-review-20261019/wt/docs/ic-white-noise-20261019.md
- the diff: git -C /cosma8/data/dp203/dc-ruan1/scratch/mgglam-ic-review-20261019/wt diff 1024c1e e604768 -- . ':!test/ic_white_noise/results'

Carry out the three reviews of the branch ic-white-noise-20261019 of MG-GLAM (fix commit e604768 on cz = 1024c1e, not merged) in order: REVIEW 1 (the Fortran fix, its callers, the builds, the end-to-end reproduction of the old draws at level 0), REVIEW 2 (the mechanism and the evidence: rebuild and re-run the standalone test, recompute the note's tables with your own code), REVIEW 3 (the note, the record, the consequences for the existing runs, the downstream users of the module, the wording). Write one findings file per review under /cosma8/data/dp203/dc-ruan1/scratch/mgglam-ic-review-20261019/findings/ as the prompts say.

Follow the hard rules of the preamble: your two worktrees (wt, cz) and work/ are the only places you build, run or write; change nothing else; commit nothing; no git state changes beyond restoring a file you changed in your worktrees; login node light tier only (one process, nice -n 19, at most 4 threads, about a minute and under 2 GB per step, load checked before each start); no Slurm; nothing in /tmp; no recursive scans; Python only as nice -n 19 /cosma/apps/durham/dc-ruan1/micromamba/envs/cosemu/bin/python3 -I; never micromamba run for a Fortran executable; nothing of yours left running at the end (ps -u $(id -un)); delete the cubes and PM files you made and say which built files you left. Plain astronomy wording; numbers in tables. Do not take the decisions the preamble reserves for the user.

Finish with the three findings files and a reply of at most one screen: for each review the merge recommendation and the counts per severity, one line per Blocker or Major, the runs you made, and the full text of any findings file you could not write.
```

## RESUME (if the session was interrupted)

```text
Resume the review of the branch ic-white-noise-20261019 of MG-GLAM (fix commit e604768, not merged). Read /cosma8/data/dp203/dc-ruan1/scratch/mgglam-ic-review-20261019/wt/docs/reviews/ic-white-noise/REVIEW-PROMPTS.md, then everything under /cosma8/data/dp203/dc-ruan1/scratch/mgglam-ic-review-20261019/findings/ and work/. Work out which review blocks and items are finished, which findings are already written and checked, and whether your builds exist (git -C <worktree> diff --stat must be empty). Tell me in a few lines what you found and what you will do next, then continue under the same rules: review only, change nothing outside your worktrees and work/, commit nothing, plain astronomy wording. Do not repeat a finished run unless its output is missing or inconsistent.
```
