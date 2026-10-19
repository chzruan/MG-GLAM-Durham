# Triage of the review of `ic-white-noise-20261019` (2026-10-19)

The three reviews ran as one fresh subagent from `REVIEW-PROMPTS.md` against the
fix commit e604768 (prompts commit 54c90ff); findings in `REVIEW-<k>-*-findings.md`.
Reported totals: **Blocker 0, Major 1, Minor 9, Nit 11.** Recommendation of every
review: "ready after the listed fixes". The fixes are in the commit that adds this
file; the branch is then merged into `cz`.

| review | Blocker | Major | Minor | Nit | scope |
|---|---|---|---|---|---|
| R1 the Fortran fix, its callers, the builds | 0 | 0 | 3 | 4 | end-to-end reproduction of the old draws on the 128^3 example: bit-identical |
| R2 the mechanism and the evidence | 0 | 0 | 2 | 4 | all 22 results files and `analyse_cubes.json` reproduced bit for bit; the mechanism established for lags 12-24 |
| R3 the note, the record, the consequences | 0 | 1 | 4 | 3 | every number traced to its source; the downstream replica of the draw unwarned |

Found sound by the reviews (not re-reviewed here): the cause (the table never set,
`nskip = 0` at `lux = 2`, every caller asked for 2); the fix introduces no race
(`/LUX2/` was written on every call before as well; every call site is serial or
inside a critical section); the streams are complete under skipping (1 and 4
threads give the same draw); `RanluxLevel` is thread-safe; `iFlag` behaved as 0
in a fresh build of the unmodified code; the test is a faithful copy of
`SetRandomN` and `SPECTR`; "along i means along x"; the noise levels; the seeds
are production seeds; the consequences table is measured, not inferred; every
launch script of the run folders calls `PMP2start.exe`, none `PMP2startMPI.exe`;
the 2LPTIC runs do not touch the generator; no forbidden word in the texts.

## Major, verified and handled

| ID | finding | done |
|---|---|---|
| R3-1 | The halobias replica `abcs_haloes/degrace_ic/src/lowk_ts.f90` (line 37) hard-codes `lux = 2`; rebuilt against the merged module it would draw at level 2 and reproduce no run, silently (`zcv/icfield.py` caches by file name). Neither the note nor `AGENTS.md` warned. | Verified from the file. A section "Downstream users of the module" in the note and a sentence in the `AGENTS.md` bullet. The one-line change in `abcs_haloes` (outside this repository) is the user's decision; the executable in use (2026-10-14, recorded by checksum in halobias) was built against the old module and draws at level 0. |

## Minors and Nits: handled

| ID | finding | done |
|---|---|---|
| R1-1 | `PMP2correlation.f90:1032` kept `lux = 2` (the random particle subsample of the correlation-function program). | `lux = RanluxLevel()`; the program builds (`make PMP2correlation`); the note names it as the fifth caller. |
| R1-2 | Permissive parsing of `GLAM_RANLUX_LUXURY` (`3.5`, `3,4`, `3 4`, `1e3` accepted silently; `/` left `ival` undefined). | The value must be digits with an optional sign in front; `ival` initialised; every other string prints the message and uses the default. Checked on 23 values (table below). |
| R1-3 | The message printed p as the level for the 24-2000 form. | The message now prints the level and p (`RanluxP`); for the 24-2000 form the "level" is p itself, as `rluxgo` treats it. |
| R1-4 | `PMP2startMPI` not built. | Compiled (object only) here with Intel MPI 2021.13 as the review did; the note says how. Not linked or run: listed under "Not established". |
| R2-1 | The test cannot link against the unmodified module (`RanluxLevel` missing). | `ranluxlevel_stub.f90` and `make ... ORIG=1`; the `orig_*` results re-made with it are bit-identical. |
| R2-2 | The mechanism paragraph conflated the first-order lags 12-24 with the second-order lags 26-30. | Rewritten; the review's table of measured against predicted-from-the-pairs added. |
| R2-3 | The cube and the 2048-row profile values differ by the end-of-row factor; the note did not say which applies to the runs. | The identity amplitude(m) = 2 c(2m) (1 - 2m/N) stated; the 2048-row values named as the ones for the L = 1024 runs (factor 0.993 at m = 11 for L = 1800). |
| R2-4 | Two cells rounded up (56, 383). | 55, 382. |
| R2-5 | "to 3e-4" and "the same ranking" overstated. | "to 4 x 10^-4"; "the same six harmonics on top". |
| R2-6 | Dead code and an unbound `hdr` in `analyse_cubes.py`. | Removed; `hdr = ''`. The `i3` formats of the test's text output are left (changing them would change every results file). |
| R3-2 | No output of a run recorded the level or the seed. | `PMP2start` and `PMP2startMPI` write the level, p, the white-noise seed and the realisation into `Results.log` right after the files are opened. Writing them into the `PMcrd.DAT` header (the review's other suggestion) would change the header of every new run: left to the user. |
| R3-3 | "Fourth-order statistics" understated a position-dependent two-point defect. | Sentence replaced as the review suggested. |
| R3-4 | The affected suites were not named. | A table of the run folders from the review's search of the launch scripts. |
| R3-5 | "Not established" incomplete. | Extended with the review's items (end-to-end at 2048^3, the MPI generator at run time, `getGauss`/`getRandom`/`PMP2correlation` untested, the halo catalogues, the replica). |
| R3-6 | Four values for one modulation amplitude. | The 2048-row value (2.95 per cent at harmonic 11; 0.5-3.0 at 6-12) in the summary table and the consequences; the cube row explained. |
| R3-7 | Numbers in sentences. | The cost and the level-2/3 relation numbers moved into tables; the rest left (the mechanism paragraph reads as prose). |

## Left as they are

| ID | finding | why |
|---|---|---|
| R1-5 | `PMP2spectrum` does not build on `cz` either (the target lacks `PMP2linker.o`). | Pre-existing; out of scope. Noted in the note. |
| R1-6 | `/LUX2/` written on every `rluxgo` call. | Not a new race (the review's own conclusion); recorded. |
| R1-7 | 119 implicit-typing warnings with `-warn all`, the old style. | Nothing new; out of scope. |
| R3-8 | The commit message of e604768 does not list the `luxuryp.h` comment, the `.gitignore` and the 22 results files. | Not amended (the branch's history stays as reviewed); the merge commit lists them. |
| R3-1 (part) | The line in `abcs_haloes/degrace_ic/src/lowk_ts.f90`. | Outside this repository; the user's decision. |

## Decisions for the user

- The default level 3 (the RANLUX authors' recommendation) rather than 2 (what the
  callers had asked for); the cost is in the note's table.
- An environment variable as the switch rather than a line of `Setup.dat`; and
  whether the level and the seed should also go into the `PMcrd.DAT` header.
- The one-line change of the halobias replica (`lux = 0`) before anyone rebuilds it.
- Whether and how the MG-GLAM authors are told.

## Verification after the fixes

Login node, `nice -n 19`, one process at a time, load 22-26 of 128 before the
starts; `make -j1`; the example at `OMP_NUM_THREADS=4`.

| check | result |
|---|---|
| builds: `PMP2init`, `PMP2start`, `PMP2correlation` (ifx 2024.2.0); `PMP2startMPI.o` (Intel MPI 2021.13, compile only); the test against this branch and, with `ORIG=1`, against `cz` | all build; no error |
| the 128^3 example, box 2: the unmodified build (`cz`) against this branch at `GLAM_RANLUX_LUXURY=0` | `PMcrd.DAT`, `PMcrs0.DAT`, `PowerDM.R`, `fort.17`, `pt.dat` bit-identical; `Results.log` differs by the new line only |
| the same, this branch at the default against level 0 | the particle files and the spectra differ; `pt.dat` identical |
| `Results.log` line of the default run | `RANLUX luxury level = 3  p = 223  white-noise seed = 2009256176  realisation 2` |
| `GLAM_RANLUX_LUXURY` parsing, 23 values: unset, `0`, `3`, `4`, `24`, `389`, `2000`, `03`, `+2`, `3 `, empty | accepted (empty and unset give 3; `24` gives p = 24 = level 0) |
| `3.5`, `1e3`, `3 4`, `3,4`, `/`, `abc`, `0x3`, `-1`, `5`, `23`, `2001`, a blank | the message and the default 3, every one |
| the committed results files re-made with the fixed sources: `fix_s1_lux0`, `fix_s1_lux3`, `fix_env_default`, `fix_env_0`, `fix_env_bad`, `orig_s1_lux0`, `orig_s1_lux2`, `orig_s1_lux3`, `fix_pairs_s1_lux0`, `orig_pairs_s1_lux2` | all ten bit-identical |

The reviewer's own programs, logs and text outputs stay in
`/cosma8/data/dp203/dc-ruan1/scratch/mgglam-ic-review-20261019/work/` (6 MB) with
the author's verification outputs in `work_author/` (0.1 MB); the reviewer's
worktrees were removed after the merge.
