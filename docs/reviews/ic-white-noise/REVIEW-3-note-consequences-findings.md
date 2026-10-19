# REVIEW 3: the note, the record and the consequences for the runs

Branch `ic-white-noise-20261019`, fix commit e604768 on `cz` = 1024c1e. Scope:
`docs/ic-white-noise-20261019.md`, the `AGENTS.md` bullet, the `luxuryp.h` comment, the commit
message; the halobias review (read-only) for the measured consequences; the launch scripts
and seed tables of the run folders (read-only); the halobias replica `lowk_ts.f90` and
`glamic.py` (read-only). Reviewer: a fresh session that took no part in the work. Line
numbers refer to the branch.

## 1. Summary

Every number of the note traces to `results/` or to the halobias review (tables in section
3); the four rows of the consequences table are measured quantities of that review, correctly
quoted; the inferred statements are sound, one of them unmarked (the mass function). Every
launch script found calls `PMP2start.exe`; none calls `PMP2startMPI.exe`; the 2LPTIC runs go
through `ic2pm`, which does not use the generator. Two things are missing. The halobias
replica `lowk_ts.f90` hard-codes `lux = 2` and would silently draw at level 2 after a rebuild
against the merged module, and neither the note nor `AGENTS.md` says so. No output of a run
records the level (or the seed); the level lives in the job's standard output only. The
"fourth-order statistics" sentence understates what is affected, and the note never names
the affected suites.

Merge recommendation: ready after the listed fixes (R3-1 is one line outside this
repository plus one sentence in the note; the rest is wording and the "Not established" list).

## 2. Findings

| ID | Severity | File:line | Description | Evidence | Suggested fix |
|---|---|---|---|---|---|
| R3-1 | Major (a silent change of a downstream reproduction after the merge; the file is outside this repository, the remedy is the user's call) | `docs/ic-white-noise-20261019.md:85-89`; `AGENTS.md:27-33`; `abcs_haloes/degrace_ic/src/lowk_ts.f90:37` | The halobias replica sets `lux = 2` literally, takes `iFlag0` as input and is "built against the repository's own `PMP2mod_random.f90` and `luxuryp.h`" (no build script exists in `abcs_haloes/degrace_ic/`; only `build/fast/lowk_ts.exe` of 2026-10-14; `glamic.py:30` loads the compiler modules and runs the executable; which copy of the module it was built from is not recorded). Rebuilt against the module after the merge, it draws at level 2 (`nskip = 73`), reproduces no DEGRACE run, and `GLAM_RANLUX_LUXURY` does not reach it. The halobias package `zcv/icfield.py:56` calls `glamic.run_lowk_ts(...)` and caches the result by file name, so the error would enter the control variates silently. The note calls the replica "a correct replica" and stops there. | `lowk_ts.f90:37` `lux = 2`; `zcv/icfield.py:45-58`; `orig_s1_lux2.txt` against `fix_s1_lux2.txt` (`nskip_used 0` against `73`) | In the note and in `AGENTS.md`: name the replica as a downstream user and say that it must set `lux = 0` (or call `RanluxLevel()` with `GLAM_RANLUX_LUXURY=0`) when rebuilt; in `abcs_haloes` (the user's decision) change line 37 and record the checksum of the executable, which the halobias review's P2-1 already asks for. |
| R3-2 | Minor (introduced; the record) | `PMP2mod_random.f90:152-153`; `AGENTS.md:31-33`; `docs/ic-white-noise-20261019.md:99,198-201` | The level is printed once to standard output and nowhere else. No output of `PMP2start` records it or the seed: `PMcrd.DAT`'s header carries the `Setup.dat` line (`Random seed` there is 0; the seed is row `Nrealization` of `TableSeeds.dat`), `Results.log`, `fort.17`, `pt.dat` and the `PowerDM` files carry neither (checked on the example, REVIEW 1). A run made after the merge without the variable set draws at level 3 and leaves no durable trace unless the job's `.dump` is kept and read. The `AGENTS.md` instruction ("set `GLAM_RANLUX_LUXURY=0` and record it in the campaign's provenance"; "provenance" is the repository's own term, `AGENTS.md:23,34,48`) is enough only when followed. Whether the switch is an environment variable or a `Setup.dat` line is the user's decision; the missing record is not. | section 3 of REVIEW 1 | Write the level (and p) and the seed into `PMcrd.DAT`'s `extras` array (`extras(:)` is written in the header, `extras(100) = Box` already) or into `Results.log`. |
| R3-3 | Minor (introduced; a central sentence understated) | `docs/ic-white-noise-20261019.md:210-212` ("the mean power of every mode is right, the pattern sits in the fourth-order statistics (the squared field)") | The defect is a second-order statistic that is not translation-invariant: the coefficients of x frequencies m and m + Delta m are correlated (lag 2 Delta m along i), i.e. the covariance of the Fourier modes between different k along the x axis is not zero, and the two-point function of the field depends on the position along x (the variance profile, 1-3 per cent). The translation-averaged power spectrum and the pair-count correlation function do not see it at first order; any statistic that is not translation-averaged (a slab or sub-box variance, a local sigma_8, a light-cone cut along x, a cross-correlation of distinct modes) sees it at first order, and the squared field is the usual way it shows. The sentence lets a reader conclude that two-point statistics are safe. | section 3 of REVIEW 2 (the identity amplitude(m) = 2 c(2m) (1 - 2m/N)); the halobias review's variance profile | Replace "fourth-order" by "a position-dependent two-point statistic (the variance along x), invisible in the translation-averaged power spectrum and seen at first order in the squared field and at the axis modes". |
| R3-4 | Minor (introduced; completeness) | `docs/ic-white-noise-20261019.md:33-35,117-119,205-207,231-233` | The note says "every caller" (four listed; `PMP2correlation.f90:1032` is a fifth, REVIEW 1 R1-1) and "every run with native initial conditions", but never names the suites. Found read-only (section 3): every launch script of the run folders calls `PMP2start.exe`; none calls `PMP2startMPI.exe`; the 2LPTIC runs (`fid2LPTIC_*`) go through `ic2pm.exe`. A user of the runs has to repeat this search. | section 3 | Add a table of the affected suites (the list in section 3) and say that `PMP2startMPI` was never used in them. |
| R3-5 | Minor (introduced; the "Not established" list) | `docs/ic-white-noise-20261019.md:244-252` | Missing from the list: (a) the end-to-end reproduction of a production 2048^3 initial condition at level 0 (the runs' PM files were deleted by the submit scripts, `rm -f PMcr*`; the evidence is indirect: the halobias review's rebuilt low-k modes agree with the runs' halo fields, and this review's 128^3 identity); (b) the MPI generator at run time (compile-only here); (c) `getGauss` and `getRandom`: no test exercises them (they are off in every `Setup.dat` of the run folders: `Redshift distortions 0`, `Biasing model 0`, `rms V correction factor 0`, so no DEGRACE product used them); (d) `PMP2correlation`; (e) the downstream replica (R3-1); (f) the second-order part of the lag pattern at lags 26-30 (REVIEW 2 R2-2); (g) the halo catalogues (named only as "measurements on the evolved field"); (h) what "record it in the campaign's provenance" means in practice (R3-2). Established by this review and worth stating: 1 and 4 threads give the same draw at level 3; `iFlag` behaved as 0 in a fresh build of `cz`. | REVIEW 1 and 2 | Extend the list. |
| R3-6 | Nit (introduced; one quantity, four values) | `docs/ic-white-noise-20261019.md:17,171,208`; commit message; `AGENTS.md:29` | The headline modulation at harmonic 11 is 2.8 per cent in the summary table (the 256^3 cube), 2.95 in the rows table (the 2048-row value, the one that applies to the L = 1024 runs), "0.5-2.9" in the consequences text, "0.5-3" in the commit message, "0.5-2.6" in the halobias review (the rebuilt linear field, which the power spectrum weights). | REVIEW 2 R2-3 | State one value for the runs (2048 rows: 0.5-3.0 per cent at harmonics 6-12, 2.95 at 11) and say what the cube row is. |
| R3-7 | Nit (introduced; wording rule) | `docs/ic-white-noise-20261019.md:41-43,58-59,82-83,112-115,162-163,179-182,199-201,207-209,222-233` | No forbidden shorthand in the note, the `AGENTS.md` bullet or the commit message ("flags" at lines 197 and 259 is the noun, compiler flags; "provenance" in `AGENTS.md` is the repository's term). Numbers in sentences that belong in tables: "21.5 per cent", "about 19, with a spread of about 3", "2-3 x 10^-7", "0.0001-0.0007", "9.3 times", "7 s", "20 s", "about 15 minutes", "0.0149 (lag 22)", "0.0002-0.0004", "0.975", "-0.02", "0.5-2.6 per cent", "0.5-2.9 per cent", "well below 1 per cent", "85-170 Mpc/h", "0.1-0.6 per cent". "The same ranking (11, 8, 12, 9, 7, 10)" is the review's seed-2 order (seed 1: 11, 8, 9, 12, 7, 10). "raw halo correlation function of one box" (line 219): the review's number is the per-box effect common to every box and model. | reading | Move the numbers into the existing tables; "the same six harmonics". |
| R3-8 | Nit (introduced; commit message) | commit e604768 message | Bullet by bullet accurate against the diff. Not described: the `luxuryp.h` comment (3 lines), `test/ic_white_noise/.gitignore`, the 22 `results/` files (only "its results"). "0.5-3 per cent" against the note's 0.5-2.9. | `git show --stat e604768` | Add the two files to the last bullet. |

What the commit discloses as not established and is still open: the L = 1800 runs (not
measured; the index-space argument is sound), the evolved matter field, `iFlag` in the
production binaries, the closed form (now first-order explicit, REVIEW 2).

The user's decisions (not counted): environment variable against `Setup.dat` line; default 3
against 2; branch of the fork against upstream; whether anything is re-run; whether and how the
MG-GLAM authors are told; the change in `abcs_haloes`; the merge.

## 3. Verified OK

Every number of the note against its source ("match" = the rounded source value):

| note (line) | value | source | source value | match |
|---|---|---|---|---|
| summary (14-15) | table never set; level 0; 3 = 223 per 24 | `cz:PMP2mod_random.f90:167`; reviewer program; 24 + 199 | as stated | yes |
| summary (16) | lag 22: +0.015 every seed; 0.0000 (noise 0.0001) | `fix_s1_lux0`, `fix_s2_lux0`, `fix_s1_lux3` lag 22 | 0.01494, 0.01482, -0.00003; noise 8.6 x 10^-5 | yes |
| summary (17) | harmonic 11: 2.8 (every seed); 0.07 (noise 0.07) | `analyse_cubes.json` cube s1 lux0 / s2 lux0 / s1 lux3 along x, m = 11; noise | 2.77 / 2.75 / 0.066; 0.065 | yes (R3-6: the rows give 2.95 / 2.93) |
| summary (18) | 160-1600; 0.4-2 | axis-mode table m = 7-12 | 164-1604; 0.4-1.4 (y, z at level 0: 0.0-1.9) | yes |
| mechanism table (64-71) | 0.99994; -0.0004; -0.0011/-0.0002; -0.458/+0.125; +0.0001/-0.0006; -0.548/+0.170 | `fix_pairs_s1_lux0.txt` | 9.9994e-1; -4.28e-4; -1.107e-3/-2.42e-4; -0.45768/0.12460; 1.10e-4/-5.75e-4; -0.54784/0.16975 | yes |
| 82-83 | 2-3 x 10^-7; 0.0001-0.0007 | `fix_pairs_s1_lux2/3.txt` | 3.0e-7, 2.0e-7; -8.9e-5, -6.9e-5, -7.2e-4, 1.3e-5 | yes |
| 41-43 | 21.5 per cent rejected | `accepted_fraction 0.785279` | 21.47 | yes |
| 112-115 | 9.3 times; 7 s / 20 s per 134 M; about 15 min for 2048^3 | 223/24; reruns 6.85 / 20.28 s; (20.28 - 6.85) x 64 | 9.29; 6.9 / 20.3; 14.3 min | yes |
| lag table (151-160) | 7 rows x 11 cells | `results/fix_*.txt` | 0 mismatches (REVIEW 2) | yes |
| review's row (155) | 9 cells | `test_whitenoise_lags.json` seed 1 | -0.00252 ... +0.00253 | yes |
| 162-163 | odd lags within noise; 0.0149 at lag 22 both seeds; 0.0002-0.0004 at levels 1-4 | `results/` | odd <= 3.7e-4; 0.01494 / 0.01482; 0.00019-0.00044 | yes |
| profile table (169-177) | 4 rows x 11 cells + 3 cube cells | `analyse_cubes.json` | 1 boundary cell (REVIEW 2 R2-4) | yes |
| 179-182 | 0.975; -0.02; 0.5-2.6 per cent at 6-12; ranking | `analyse_cubes.json`; `test_variance_profile.json` | 0.9747; -0.0225; 0.0047-0.0262; (11, 8, 9, 12, 7, 10) and (11, 8, 12, 9, 7, 10) | yes (R3-7) |
| axis table (188-193) | 5 rows x 10 cells | `analyse_cubes.json` | 2 cells rounded up (REVIEW 2 R2-4) | yes |
| 195-201 | builds; MPI not built; switch unset / 0 / invalid | REVIEW 1; `fix_env_*.txt` | as stated | yes |
| consequences (217) | +37, +775, +25, +30, +94, +80; every other shell 0.0 +- 1.8 | `implementation/REPORT.md` table; `science/REPORT.md:36,53` ("all 213 other shells, n^2 <= 300: mean +0.04, rms 1.8") | +36.8, +775, +25.2, +30.0, +93.5, +80.3 | yes |
| consequences (218) | halo power at n^2 = 64: +20 to +67 per cent | `REVIEW.md:33-35`; `test_raw_xi_effect.json` `P_hh_excess_frac[1]` | 0.199-0.671 over the seven levels | yes |
| consequences (219) | raw xi_hh at 29.6-40.5: +0.1 to +0.6 per cent | `test_raw_xi_effect.json` windows; `REPORT.md:84-86` | +0.10 (lg n -2.90) to +0.55 (-5.25); level 0 (-6.03) -0.36, not quoted by the review either | yes |
| consequences (220) | linear field within 1-2 per cent | `test_anomalous_shells.json` `one_one_pct` | -1.5 to +1.8 | yes |
| record (258-263) | 1.3-20 s, 28-142 MB; 2-5 s, 312 MB; 53 s, 839 MB; no Slurm job; cubes deleted | this review's reruns | 1.5-20.3 s (0.1 s for the env runs), 25-142 MB; 2.4-5.8 s, 307-309 MB; 44.9 s, 832 MB; the cubes are indeed absent from the commit | consistent |

Other checks:

| property | how | result |
|---|---|---|
| The consequences table rows are measured, not inferred | each traced to a results file of the halobias review on the 95 fiducial L = 1024 boxes | yes, all four |
| Inferred statements | "bin-averaged power spectra ... well below 1 per cent": sound for the matter field (six axis modes among about 800 in a bin of unit width in n at n^2 = 64, excess of the linear field within 2 per cent) and for the halo field at n^2 = 64 (+20 to +67 per cent on six modes: 0.15-0.5 per cent of the bin); "mass functions": not argued in the note; a 1-3 per cent modulation of the variance in slabs along x changes the box mass function at second order only (of order the square of the slab-to-slab sigma_8 variation, 10^-4), so the statement holds but is inferred and should say so; "correlation function below about 20 Mpc/h": measured by the review at 10.4-21.3 Mpc/h (+0.04 to +0.22 per cent), so not only inferred; "does not average out over realisations": sound (the pattern is a property of the recurrence, identical for the two seeds to 3 x 10^-4 and present in the 95-box means); "L = 1800 ... the same in index space": sound (same code path, rows of 3200; the end-of-row factor is 0.993 at m = 11) | marked as inferred except the mass function |
| Which runs are affected (read-only, launch scripts) | `mg_glam/*/arun*.sh`: 30 scripts, all `PMP2start.exe` (the `DurMun_hmfemu_LCDM/fRn1_wide_sample_first_64_model*` suites, `DEGRACE2_fid_LCDM_L1024`, `DEGRACE2_fid_LCDM_L1024Np4096Ng8192`, `F5n1_L1024Np2048Ng4096`, `F6n1_L1024Np2048Ng4096`, `DurMun_hmfemu/arun_LCDM.sh` and `arun_fR.sh`; the `DurMun_hmfemu_fRn1_..._model*` folders hold no launch script of their own and were launched from the latter, which calls `PMP2start.exe`; no `arun*.sh` under `mg_glam` names `ic2pm`); main checkout `*/arun*.sh`: `DEGRACE_LCDM/F4n1/F5n1_L1800Np3200Ng6400`, `fid_F5`, `fid_F5n1_L512Np2048Ng4096`, `fid_LCDM`, `fid_LCDM_L1024Np4096Ng8192`, `fid_LCDM_L512Np2048Ng4096`, `fid_LCDM_L512Np2560Ng5120`, `LCDM_L1024Np2048Ng4096`, `test`; the supplement folders (`..._supplement_16_model68..72`): `Run1/submit.sh` calls `PMP2start.exe`; model 1 box 1: `Run1/submit.sh` `../../PMP2start.exe <<< 1` with `OMP_NUM_THREADS=128`, `.dump` of 2024-07-23 shows `SPECTR ==> SUMM= ...` | native generator everywhere; no script names `PMP2startMPI` |
| 2LPTIC runs are not affected | `git grep -n -i 'ranlux\|rluxgo' -- '2LPTIC_Gui' 'ic2pm.f90'` (no hit); `fid2LPTIC_L512Np2048Ng4096/Run1/submit.sh:22` `../../ic2pm.exe ... 1.0` then `PMP2main.exe` | right; the evolution's optional draws (`getGauss` for redshift-space spectra, `getRandom` for the biasing) are off in those `Setup.dat` too |
| Anything else in this repository reproduces old draws with the module | `git grep -n -i 'rluxgo' -- '*.f90' 'python' 'scripts' 'test'` | only the generators, `getGauss`, `getRandom`, `PMP2correlation` and the test |
| Seeds of the note's two tests are production seeds | `TableSeeds.dat` of model 1 rows 1-2; of F5 row 2 | 1156811761, 2009256176; 2009256176 |
| `AGENTS.md` bullet | against the diff and the note | accurate; placed before the recording bullet; links the note |
| `luxuryp.h` comment | `luxuryp.h:9-10` | accurate ("set in rluxgo since 2026-10-19; before that the table was never set and was zero") |
| Commit message | bullet by bullet against `git show --stat e604768` | accurate (R3-8 for omissions) |
| Forbidden words | `grep -w` over the note, the bullet, the message, the code comments of the diff | none (the two "flags" are nouns; "provenance" is the repository's term) |

## 4. Open questions for the author

1. Will `abcs_haloes/degrace_ic/src/lowk_ts.f90` be changed before anyone rebuilds it, and
   which copy of the module was `build/fast/lowk_ts.exe` built from (R3-1)?
2. Should `PMP2start` write the level and the seed into the `PMcrd.DAT` header (R3-2)?
3. Which modulation value is "the" value for the L = 1024 runs in the texts (R3-6)?
4. Is the `PMP2correlation` program in use anywhere (REVIEW 1 R1-1)?

## 5. Runs made

No compute for this review beyond REVIEW 1 and 2: read-only `grep`/`ls`/`sed` on the run
folders (`mg_glam/*/arun*.sh`, one `Run1/submit.sh` and one `.dump` of model 1, the
supplement model-68 `Run1/submit.sh`, the main checkout's `*/arun*.sh`, `fid2LPTIC_*/Run1/submit.sh`,
the `Setup.dat` of six folders), the halobias review files named above, and
`abcs_haloes/degrace_ic/{src,build/fast}`. One listing loop over `mg_glam/*/*/submit.sh` took
about three minutes on the file system (it ran in the background and finished by itself; no
process was left). Nothing written outside `findings/` and `work/`.
