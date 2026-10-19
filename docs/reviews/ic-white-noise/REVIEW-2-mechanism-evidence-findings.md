# REVIEW 2: the mechanism and the evidence (`test/ic_white_noise/`, the note's tables)

Branch `ic-white-noise-20261019`, fix commit e604768 on `cz` = 1024c1e. Reviewed in the
detached worktrees `wt` and `cz` under
`/cosma8/data/dp203/dc-ruan1/scratch/mgglam-ic-review-20261019/`; every run under `work/`.
Reviewer: a fresh session that took no part in the work. The reviewer's own code is
`work/mech/mech.f90` (Fortran, built against the unmodified and the fixed module) and
`work/py/review2.py` (its report is `work/py/out/review2_report.txt`).

## 1. Summary

The test is a faithful copy of `SetRandomN` and `SPECTR` (seeding, draw order, the restore of
the stream per plane, `iFlag = 0`), its estimators do what the output claims, and "along i
means along x" is right (`VECTOR` transforms index i first, `BLOCKS` maps `ip` to `XPAR`).
Every one of the 22 committed results files was reproduced bit for bit, and
`analyse_cubes.json` rebuilt identically from the regenerated cubes. Every cell of the note's
three tables is the value in `results/` (three cells rounded up by at most 0.6). The
mechanism is established for lags 12-24: the three-term relation holds exactly (residual 0 or
one carry) in both components for 99.93 per cent of the pairs, and the products of the pairs
(p, p+12) alone reproduce the measured Gaussian lag correlation at every even lag 12-24 to
4 x 10^-4; the correlations at lags 26-30 (up to 0.0025 at lag 28, the size of lag 12) are
not produced by that relation and the note should not list them under the same explanation.
The committed test cannot be linked against the unmodified module as the note describes
(`RanluxLevel` is missing there); the `orig_*` files need a stub.

Merge recommendation: ready after the listed fixes (R2-1 and R2-2 are wording and build
matters; nothing in the evidence is wrong).

## 2. Findings

| ID | Severity | File:line | Description | Evidence | Suggested fix |
|---|---|---|---|---|---|
| R2-1 | Minor (introduced; reproducibility of the evidence) | `test/ic_white_noise/ranlux_lags.f90:67`, `makefile:7` | `make GLAM=<cz> B=build_cz` compiles but does not link: `undefined reference to ranluxlevel_` (the test calls `RanluxLevel()` when `lux < 0`; the unmodified module has no such function; with implicit typing the compile step passes). The `orig_*` results of the note cannot be rebuilt from the committed sources against `cz`. With a two-line stub (`work/stub/ranluxlevel_stub.f90`) linked in, all four `orig_*` files reproduce bit for bit. | section 5, builds | Make the `lux < 0` branch optional (a preprocessor symbol, or a second mode), or document the stub in the makefile and the note. |
| R2-2 | Minor (introduced; the explanation, not the numbers) | `docs/ic-white-noise-20261019.md:49-83` ("Mechanism"), commit message | The relation pair p+12 = pair p+7 - pair p, with the rejection of pair p+7, explains the lags 12-24 only: the delivered lag between the first Gaussians of pairs p and p+12 is 2 x (accepted pairs among p..p+11), at most 24. Predicted from those products alone against measured (table in section 3): lags 12-24 agree to <= 4 x 10^-4; lags 26, 28, 30 (+0.0009, +0.0025, +0.0013) have no such pairs and come from chained relations (pair p+7 is shared by the triplets (p, p+7, p+12) and (p+7, p+14, p+19)). The note's "even lags 12-28" and "the mechanism is shown by the conditional products" conflate the two. The note's own "Not established" (no closed form) covers the gap only in part. | `work/mech/mech_cz_lux2.txt`; `work/py/out/review2_report.txt` section 6 | Say that the one relation accounts for lags 12-24 (the table of predicted against measured), and that lags 26-30 are second-order consequences of the same recurrence. |
| R2-3 | Nit (introduced; a reader's trap) | `docs/ic-white-noise-20261019.md:165-178` (profile table) | The 256^3 cube rows are 6-8 per cent below the "rows of 2048" rows (harmonic 11: 2.77 against 2.95) because the modulation of a row of N numbers is 2 c(2m) (1 - 2m/N): the cube loses the lag pairs that fall off the row end. The note does not say which row applies to the runs (the 2048-row values; for the L = 1800 runs with N = 3200 the loss is 0.3 per cent at m = 11). The identity amplitude(m) = 2 c(2m) (1 - 2m/N) holds for every m = 5-14 of both tables (section 3). | section 3 | One sentence under the table; state the identity, which also ties the lag table and the profile table together. |
| R2-4 | Nit (introduced) | `docs/ic-white-noise-20261019.md:188-193`, `:173` | Rounded up: 55.48 -> 56 (s1, x, m = 6), 382.47 -> 383 (s2, x, m = 7); one noise cell at the rounding boundary (rows level 3 s1, m = 9: 0.035 -> 0.04). | section 3 | Round half to even or copy the two-decimal values. |
| R2-5 | Nit (introduced) | `docs/ic-white-noise-20261019.md:155,179-182` | Against the halobias review's level-0 table the largest difference is 3.7 x 10^-4 (seed 1, lag 28) and 2.4 x 10^-4 (seed 2), i.e. within the combined noise of two different samples of the same streams (the review's cube keeps the first columns only); the note's row for the review quotes seed 1 only though both seeds are in the JSON. The note's "same ranking (11, 8, 12, 9, 7, 10)" is the review's seed-2 ranking; the review's seed-1 ranking is (11, 8, 9, 12, 7, 10) (0.0164 against 0.0161, within its noise 0.0006). | section 3 | Say "to 4 x 10^-4" and "the same six harmonics". |
| R2-6 | Nit (introduced) | `test/ic_white_noise/analyse_cubes.py:83-86`, `:29-37` | Dead code (`tags = ...; if len(tags) == 2: pass`); `read_lags` returns `hdr` unbound when the text file has no `lux_requested` line; `ranlux_lags.f90:156,186` print `***` for p >= 1000 (`i3`). Cosmetic. | reading | Remove; `hdr = ''` before the loop; `i5`. |

What the commit discloses as not established and this review confirms: no closed form of the
lag correlations (R2-2 narrows the gap: the first-order term is now explicit).

## 3. Verified OK

| property | how | result |
|---|---|---|
| Faithfulness of mode 1 to `SetRandomN` and `SPECTR` | line-by-line: `ranlux_lags.f90:69-89` against `PMP2start.f90:81-99` (seed `Ns = Nseed`, `rluxgo(lux, Ns, 0, 0)` per plane, `RANDd(Ns)` after, save of `i24, j24, in24, kount, carry, gSet, seeds`), `:98-114` against `:306-324` (restore, `mkount = 0`, j outer, i inner, one `GAUSS3` per cell, the (1,1,1) cell drawn too) | identical; `iFlagA(k) = 0` where `SetRandomN` stored the undefined local |
| The per-plane seeds equal a production run's | `Initialize` reads row `Nrealization` of `../TableSeeds.dat` into `Nseed` (`PMP2start.f90:148-153`); `TableSeeds.dat` of model 1 (L = 1024): rows 1, 2 = 1156811761, 2009256176; of F5: row 2 = 2009256176 (row 1 = 123456) | the test's two seeds are fiducial LCDM boxes 1 and 2 and F5 box 2 |
| "Along i means along x" | `SPECTR` `GRX(mapz(i), j, k) = TRX*map3(i)` (`:336`), `VECTOR` first `rfft1b` over `r(i) = GRX(i, j, k)` (`:470-479`), `BLOCKS` `XPAR(Icurrent) = Q1 - XCONS*GRX(ip, jp, kp)` (`:398-402`) | index i is the x frequency and, after the transform, the x position |
| The lag estimator along i | pooled products within rows, `sl(L)/nl(L)` minus mean squared over the variance (`:124-129, 162`); reviewer's own estimator on the raw 256^3 cube (mean removed, all three axes) | identical to 4 decimals at lags 12-28 along i; j, k within noise |
| `.vprof` is the variance of the row-synthesised white field along x | reviewer's row synthesis `B(x, i) TS(i, j, k)` squared and averaged over rows and planes | relative difference <= 1.1 x 10^-8 from `.vprof`; against the full 3-D synthesis (which weights the constant and Nyquist rows twice) the harmonics differ by <= 0.01 (x100) |
| The basis the test writes | against the analytic basis: column 1 constant, columns 2k and 2k+1 cos and sin of harmonic k, column N the Nyquist cosine | largest difference 1.7 x 10^-8; column norms 256 (constant, Nyquist) and 128 (others), off-diagonal <= 2 x 10^-6 |
| Reruns against `results/` (22 files) | `diff` of every text output after rebuilding the test (`make`; `make GLAM=<cz> B=build_cz` plus the stub) | all 22 bit-identical, including `nskip_used 0` for `orig_s1_lux2/3` and `orig_pairs_s1_lux2` |
| `analyse_cubes.json` | `analyse_cubes.py work/t2 work/t2/out` on the regenerated cubes | `cmp` identical to the committed file; 45 s, 832 MB |
| Every cell of the note's lag table (7 rows x 11 cells) | parsed from `results/fix_s*_lux*.txt`, rounded to 4 decimals | 0 mismatches; odd lags <= 3.7 x 10^-4 at every level; the largest value over lags 1-40 at level 0 is 0.01494 (s1) and 0.01482 (s2) at lag 22, at levels 1-4 0.0002-0.0004 |
| Every cell of the profile table and of the axis-mode table | parsed from `results/analyse_cubes.json` | 0 mismatches beyond R2-4; `profile_correlation` 0.9747 (level 0) and -0.0225 (level 3) as quoted; the summary table's "0.07 per cent (noise 0.07)" is harmonic 11 of the level-3 cube (0.066, noise 0.065), its "2.8 per cent" the level-0 cube (2.77; the s2 cube, not in the note, gives 2.75) |
| Level-0 lag table against the halobias review | `test_whitenoise_lags.json`, lags 1-30 | largest difference 3.7 x 10^-4 (s1) and 2.4 x 10^-4 (s2) along i, 5 x 10^-4 along j, k; seed 1 against seed 2 in the test: largest difference 3.4 x 10^-4 over 40 lags (noise of a difference 1.2 x 10^-4) |
| The mechanism's three-term relation (both components, the carry and the padding accounted for) | reviewer program, 20 million pairs, unmodified module, `lux = 2`: residual u(p+12) - u(p+7) + u(p) mod 1 classified exactly | residual 0: 49.96 per cent; residual -2^-24 (the carry): 49.97 per cent; a padded small number: 0.06 per cent; other: 0.012 per cent; the same for component 2; at level 3: 5 x 10^-8 |
| Conditional products | reviewer program | x1 and x2 components: -0.458 / +0.125 given pair p+7 rejected / accepted; given p+6 and p+8: <= 0.0011; Gaussians g_a (delivered first) and g_b (second): -0.548 / +0.170 both; all as the note |
| Lag range and sign pattern | histogram of the delivered lag between pairs p and p+12 (both accepted) | mean 19.23, spread 2.73 (binomial 2 + 2 x 11 x 0.785 = 19.28, 2.72); given p+7 rejected 17.7, accepted 19.7: the anti-correlation sits at the short lags, the correlation at the long ones; odd lags carry nothing because the relation ties the same component of the two pairs |
| Noise levels | 1/sqrt(N): 8.6 x 10^-5 (134 M), 1.7 x 10^-4 (34 M); measured rms of the level-3 lag correlations over 40 lags 7.9 x 10^-5; profile harmonic noise predicted 0.022 (rows of 2048) and 0.061 (cube), the largest of 24 about 2.2 times that; measured 0.019-0.020 and 0.065, largest 0.13-0.14 | the quoted 0.0001-0.0002, 0.02, 0.06-0.14 are right; the level 1-4 values (<= 4.4 x 10^-4 at 34 M values, 2.5 sigma for the largest of 40 lags) are noise |
| Rows 1-64 against rows 193-256 (seed 1, 256 planes, rows of 2048) | reviewer program mode 2 | largest difference over lags 1-30: 5.7 x 10^-4 at lag 24 (noise of a difference 2.4 x 10^-4); the even-lag pattern is the same in both sets (lag 22: 0.0151 / 0.0149; lag 16: -0.0116 / -0.0117) |
| Axis-mode reference | the shell "|n| within 0.5 of m" without the axis modes holds 344-2492 modes (m = 5-14) with a relative scatter of 2-5 per cent; the two-mode numerator dominates the noise | fair; at level 0 the two x-axis modes hold 20-71 per cent of the whole shell's power at m = 6-12 (0.1-0.8 per cent at level 3) |
| Sizes | the mechanism lives in the stream and is independent of the row length and the number of rows except for the end-of-row loss (R2-3); the rows test shows the same table in rows 1-64 and 193-256 | 2048-column rows are the right test; the 256^3 cube is a fair stand-in once the (1 - 2m/N) factor is allowed for |

Lag correlation of the delivered Gaussians, measured against the part predicted from the
(p, p+12) products alone (reviewer program, 20 million pairs, level 0):

| lag | 12 | 14 | 16 | 18 | 20 | 22 | 24 | 26 | 28 | 30 |
|---|---|---|---|---|---|---|---|---|---|---|
| measured | -0.0027 | -0.0071 | -0.0119 | -0.0085 | +0.0053 | +0.0152 | +0.0088 | +0.0009 | +0.0025 | +0.0013 |
| predicted from (p, p+12) | -0.0027 | -0.0071 | -0.0116 | -0.0083 | +0.0058 | +0.0157 | +0.0091 | 0 | 0 | 0 |
| pairs (p, p+12) at this lag, millions | 0.17 | 0.63 | 1.62 | 2.94 | 3.54 | 2.56 | 0.84 | 0 | 0 | 0 |

Variance-profile amplitude against the lag table, 2 c(2m) (1 - 2m/N), per cent:

| m | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 |
|---|---|---|---|---|---|---|---|---|---|
| rows of 2048, note | 0.53 | 1.43 | 2.30 | 1.65 | 1.01 | 2.95 | 1.74 | 0.22 | 0.43 |
| 2 c(2m) (1 - 2m/2048) from `fix_s1_lux0` | 0.53 | 1.43 | 2.30 | 1.65 | 1.01 | 2.96 | 1.74 | 0.22 | 0.43 |
| 256^3 cube, note | 0.53 | 1.37 | 2.13 | 1.54 | 0.90 | 2.77 | 1.54 | 0.13 | 0.43 |
| 2 c(2m) (1 - 2m/256) from the cube's own lags | 0.55 | 1.40 | 2.17 | 1.52 | 0.92 | 2.79 | 1.56 | 0.14 | 0.43 |

## 4. Open questions for the author

1. Should the note give the first-order explanation (the table above) and say that lags 26-30
   are second order (R2-2)?
2. Is the test meant to be buildable against `cz` after the merge (R2-1)? If not, the
   `orig_*` files should be marked as made with a stub.
3. Which value of the modulation does a user of the runs take: the 2048-row value (2.95 per
   cent at harmonic 11 for L = 1024) or the cube value (2.77)? (R2-3)

## 5. Runs made

Login node, `nice -n 19`, one process at a time, `OMP_NUM_THREADS=1` for the Fortran runs and
4 for Python (`MKL`/`OPENBLAS` threads 4), `uptime` checked before each batch (1-min load
8-27 of 128; one spike to 87 from other users during the reruns), available memory >= 1.4 TB.

| run | command (stdin as in the note) | wall s | peak RSS MB |
|---|---|---|---|
| test builds | `make`; `make GLAM=../../../cz B=build_cz` (link fails, R2-1); stub link `ifx <FFLAGS> -o build_cz/ranlux_lags.exe build_cz/*.o work/stub/ranluxlevel_stub.o` | 3 + 3 | - |
| `fix_s1_lux0`, `fix_s2_lux0` | mode 1, 2048 x 256 x 256, level 0 | 6.9, 6.9 | 75, 76 |
| `fix_s1_lux3`, `fix_s2_lux3` | mode 1, level 3 | 20.3, 20.3 | 74, 73 |
| `fix_s1_lux1/2/4` | mode 1, 64 rows | 2.2, 3.1, 8.0 | 25-27 |
| `orig_s1_lux0/2/3` | mode 1, 64 rows, stub build | 1.8-1.9 | 26-27 |
| `fix_pairs_s1_lux0/2/3`, `orig_pairs_s1_lux2` | mode 2, 20 million pairs | 2.4, 3.5, 5.8, 2.4 | 307-309 |
| `fix_env_default/0/bad` | mode 1, 8 x 8 x 2048, `lux = -1` | 0.1 | 9-11 |
| `cube256_s1_lux0`, `cube256_s2_lux0`, `cube256_s1_lux3` | mode 1, 256^3 with cube | 1.6, 1.5, 3.2 | 140-142 |
| `analyse_cubes.py work/t2 work/t2/out` | the branch's script | 44.9 | 832 |
| `mech_cz.exe` mode 1 (20 M pairs, lux 2), `mech_fix.exe` mode 1 (lux 3) | reviewer program | 2.6, 5.6 | 540, 541 |
| `mech_cz.exe` mode 2 (256 planes, rows of 2048, rows 1-64 and 193-256) | reviewer program | 51 | 10 |
| `review2.py` | reviewer analysis (full 3-D synthesis of two cubes, FFTs) | 20.7 | 1266 |

Deleted at the end: the three cubes and their basis files (3 x 67 MB + 3 x 0.5 MB) under
`work/t2/`. Left: the text outputs under `work/t2/` and `work/t2/out/`, `work/mech/`,
`work/py/`, `work/stub/`, the test builds `wt/test/ic_white_noise/build/` and `build_cz/`
(both ignored by git). `git diff --stat` and `git status --short` are empty in both
worktrees; no process of the reviewer is left.
