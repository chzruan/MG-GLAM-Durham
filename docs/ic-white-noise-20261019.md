# The white noise of the native initial conditions was not homogeneous along x (2026-10-19)

Branch `ic-white-noise-20261019`, from `cz`. Test and results in
[`test/ic_white_noise/`](../test/ic_white_noise/). Found by the blind review
of the halobias control-variate work on 2026-10-19
(`halobias/analysis/zcv_xihh_review_20261019_v1`, finding 1 / P1-1), which
measured the pattern on the DEGRACE runs and on the rebuilt noise but did not
locate its origin in the code.

## Summary

| | before | after (default) |
|---|---|---|
| luxury table `ndskip` of RANLUX | never set, so all zero | 0, 24, 73, 199, 365 |
| luxury level actually used, whatever the callers ask | 0 (plain subtract-with-borrow, no numbers discarded) | 3 (223 numbers generated per 24 delivered) |
| correlation of the white noise along the first index, lag 22 | +0.015 (every seed) | 0.0000 (noise 0.0001) |
| variance of the white field along x, harmonic 11 of the box (rows of 2048, the runs' row length) | 2.95 per cent (every seed) | 0.01-0.04 per cent (noise 0.02) |
| initial conditions of the 128^3 example, unmodified build against this branch at `GLAM_RANLUX_LUXURY=0` (the review) | | bit-identical |
| power of the squared field at the x-axis modes m = 7-12, over the neighbours | 160-1600 | 0.4-2 |
| switch to reproduce old runs | - | `GLAM_RANLUX_LUXURY=0` |

## What the code did

`PMP2mod_random.f90` holds the RANLUX generator of James (1994), moved into a
Fortran module. The original routine sets the table of numbers to discard per
luxury level with a `DATA ndskip / 0, 24, 73, 199, 365 /` statement. In this
repository that line is a comment (so since the first commit `51cb6ed`), and
nothing else sets the array `ndskip` of the common block `/LUX2/`
(`luxuryp.h`). The array is zero. `rluxgo(lux, ...)` then sets
`nskip = ndskip(lux) = 0` for every `lux`, and the generator delivers every
number of the subtract-with-borrow recurrence
x(n) = x(n-10) - x(n-24) - carry (mod 1): luxury level 0, the RCARRY
generator, which the header of the module itself describes as "fails many
tests". Every caller asks for level 2 (`lux = 2`: `SetRandomN` in
`PMP2start.f90` and `PMP2startMPI.f90`, `getRandom` in `PMP2mod_random.f90`,
`getGauss` in `PMP2mod_density.f90`, and `ReadDataPMcoords` in
`PMP2correlation.f90`, which draws the random particle subsample of the
correlation-function program; that program is not built by `make all` and
no launch script of the runs uses it) and got level 0.

The initial conditions draw one RANLUX stream per plane k (`SetRandomN`:
`rluxgo(2, Ns)`, Ns advanced by the congruential generator `RANDd`) and take
the Gaussian numbers of the plane in the order (j, i), i innermost, one
`GAUSS3` per cell (`SPECTR`). `GAUSS3` is the polar Box-Muller method: a pair
of uniform numbers (u1, u2) from `ranlux(RanNum, 2)` gives a point of the
square, the point is rejected if it lies outside the unit circle (21.5 per
cent of the pairs), an accepted pair gives two Gaussian numbers. The first
index i of the noise cube TS(i, j, k) is the x frequency index of GLAM's
real basis (cosine and sine pairs).

## Mechanism

With no numbers discarded, the recurrence ties every uniform number to the
ones 10 and 24 draws before it. In terms of the Box-Muller pairs (pair p is
the uniform numbers 2p-1 and 2p), pair p+12 = pair p+7 - pair p (mod 1,
component by component). This three-term relation is invisible in any
two-point statistic of the uniform numbers: given pair p, pair p+12 is
uniform over the square because pair p+7 is. The rejection step breaks this.
Whether pair p+7 is accepted depends on pairs p and p+12 jointly, and the
number of Gaussian numbers delivered between pair p and pair p+12 depends on
that acceptance. So the delivered Gaussian numbers are correlated at a lag
that depends on the acceptances in between: on average 12 pairs x 0.785
accepted x 2 numbers, i.e. about 19, with a spread of about 3. Measured on
20 million pairs of one stream (`ranlux_lags.exe`, mode 2, seed 1156811761;
x = 2u - 1 is the first component of a pair, the products are divided by the
variance of x):

| quantity, level 0 | value |
|---|---|
| fraction of pairs for which u(p+12) = u(p+7) - u(p) holds to 1.5 x 2^-24 | 0.99994 |
| mean of x(p) x(p+12), all pairs | -0.0004 |
| mean of x(p) x(p+12), pair p+6 rejected / accepted | -0.0011 / -0.0002 |
| mean of x(p) x(p+12), pair p+7 rejected / accepted | -0.458 / +0.125 |
| mean of x(p) x(p+12), pair p+8 rejected / accepted | +0.0001 / -0.0006 |
| mean of the Gaussian product g(p) g(p+12), both accepted, pair p+7 rejected / accepted | -0.548 / +0.170 |

When pair p+7 is rejected the two pairs are anti-correlated (their sum is
near 1 in each component, the most probable way to put pair p+7 in a corner
of the square), and fewer Gaussian numbers lie between them: the delivered
Gaussian numbers are anti-correlated at the shorter lags (12-18) and
correlated at the longer ones (20-24). The relation ties the same component
of the two pairs, so only even lags are affected. The delivered lag between
the first Gaussians of pairs p and p+12 is twice the number of accepted
pairs among p to p+11, at most 24, so this one relation accounts for the
lags 12-24; the smaller values at lags 26-30 come from chained relations
(pair p+7 belongs to the triplets (p, p+7, p+12) and (p+7, p+14, p+19)).
The review of 2026-10-19 measured the part of the lag correlation that the
products of the pairs (p, p+12) alone produce (20 million pairs, level 0;
`docs/reviews/ic-white-noise/REVIEW-2-mechanism-evidence-findings.md`):

| lag | 12 | 14 | 16 | 18 | 20 | 22 | 24 | 26 | 28 | 30 |
|---|---|---|---|---|---|---|---|---|---|---|
| measured | -0.0027 | -0.0071 | -0.0119 | -0.0085 | +0.0053 | +0.0152 | +0.0088 | +0.0009 | +0.0025 | +0.0013 |
| from the pairs (p, p+12) alone | -0.0027 | -0.0071 | -0.0116 | -0.0083 | +0.0058 | +0.0157 | +0.0091 | 0 | 0 | 0 |

The pattern depends on the recurrence alone, not on the seed, and it lies
along i because that is the direction of the stream. At any level above 0
the delivered numbers are separated by discarded ones and the relation no
longer connects them:

| quantity | level 0 | level 2 | level 3 |
|---|---|---|---|
| fraction of pairs for which the relation holds | 0.99994 | 3 x 10^-7 | 2 x 10^-7 |
| largest conditional product of the pairs p, p+12 | 0.458 | 0.0007 | 0.0001 |

The reproduction of the halobias review (`abcs_haloes/degrace_ic/src/lowk_ts.f90`
with `glamic.py`) copies `SetRandomN` and `SPECTR` faithfully and compiles the
repository's own `PMP2mod_random.f90` and `luxuryp.h`, so it inherits the
zero table and reproduces the level-0 draw of the runs. It is a correct
replica, not the source of the pattern.

### Downstream users of the module

That replica, `abcs_haloes/degrace_ic/src/lowk_ts.f90` (line 37), hard-codes
`lux = 2` and is built against a copy of this module. Its executable of
2026-10-14 (`abcs_haloes/degrace_ic/build/fast/lowk_ts.exe`, recorded by
checksum in halobias `analysis/zcv_xihh_20261018_v1/records/external_inputs.json`)
was built against the old module and draws at level 0, as the review of the
halobias work established from its agreement with the runs. Rebuilt against
the module of this branch it would draw at level 2 (p = 97), reproduce no
run, and `GLAM_RANLUX_LUXURY` would not reach it; the halobias package
`zcv/icfield.py` caches its output by file name, so the error would enter the
control variates silently. Before any rebuild its line 37 must ask for level
0 (`lux = 0`). No other program in this repository is built against the
module to reproduce old draws (`git grep -n -i rluxgo -- '*.f90' python
scripts test`: the two generators, `getGauss`, `getRandom`,
`PMP2correlation` and the test).

## The fix (this branch)

1. `PMP2mod_random.f90`, `rluxgo`: the table `ndskip` is set to
   0, 24, 73, 199, 365 on every call. The requested level now takes effect.
2. `PMP2mod_random.f90`: new function `RanluxLevel()`, which returns the
   level for every stream the code starts. The environment variable
   `GLAM_RANLUX_LUXURY` selects it (0-4 are Luscher's levels, 24-2000 is the
   number generated per 24 delivered); the default is 3 (223 per 24), the
   level the RANLUX authors recommend. The level is printed once at the
   first use. `GLAM_RANLUX_LUXURY=0` reproduces bit for bit the draws of
   every run made before this change.
3. The five callers use `RanluxLevel()` instead of the literal `lux = 2`
   (`SetRandomN` of `PMP2start.f90` and `PMP2startMPI.f90`, `getRandom`,
   `getGauss`, `ReadDataPMcoords` of `PMP2correlation.f90`). The variable
   must be an integer written with digits (an optional sign in front);
   `3.5`, `3,4`, `1e3`, `/` and the like print a message and give the
   default (the first version of this branch accepted some of them
   silently: review finding R1-2).
4. `SetRandomN` (both files) and `getGauss` set `iFlag = 0` before the first
   `GAUSS3`; `SetRandomN` stored an undefined local variable. In the
   production builds it behaved as 0 (a non-zero value would have made the
   first Gaussian of every plane zero and removed the modes (0, 0, n) along
   z, which the review's checks of the linear field at the axis shells rule
   out); the line only makes this explicit. In a fresh build of the
   unmodified code the variable behaved as 0: its initial conditions of the
   128^3 example are bit-identical to this branch's at level 0 (below).
5. `PMP2start` and `PMP2startMPI` write the level, p and the white-noise
   seed of the realisation into `Results.log` of the run, right after the
   files are opened. Before this branch no output of a run recorded either
   (the `PMcrd.DAT` header carries the `Setup.dat` line, whose `Random
   seed` is 0; the seed is the realisation's row of `TableSeeds.dat`); the
   level is also printed once on standard output.

Cost of the new default:

| quantity | level 0 | level 3 |
|---|---|---|
| numbers generated per 24 delivered | 24 | 223 |
| the test's draw of 134 million Gaussian numbers with its statistics, one core | 6.9 s | 20.3 s |
| `SPECTR` of the 128^3 example, 4 threads (the review) | 0.11 s | 0.16 s |
| a 2048^3 run, scaled from the test | | about 14 extra minutes of one core, shared by the OpenMP threads of `SPECTR` |

Runs from 2LPTIC (the default for new simulations, `AGENTS.md`) do not use
this generator and are not affected. Every run whose initial conditions came
from `PMP2start` or `PMP2startMPI` carries the pattern.

## Reproduction and verification

`test/ic_white_noise/ranlux_lags.f90` calls the module's own `rluxgo`,
`ranlux`, `RANDd` and `GAUSS3` with the seeding of `SetRandomN` and the draw
order of `SPECTR`. Mode 1 draws planes 1-256 of rows 1-256 (or 64) of 2048
numbers each (the row length of the DEGRACE runs) and measures the lag
correlation along i on all columns and along j and k on the first 256
columns; it also synthesises each row with `rfft1b` (the basis of `VECTOR`)
for the variance profile along x. Mode 2 is the mechanism table above.
`analyse_cubes.py` takes a 256^3 cube, builds the white field with the same
basis along the three axes, and measures the variance profile along each
axis and the power of the squared field at the axis modes. Build and run
(one core, seconds each, under `nice`):

    module load intel_comp/2024.2.0 compiler-rt tbb compiler
    cd test/ic_white_noise && make
    make GLAM=<tree from before 2026-10-19> B=build_orig ORIG=1    # against the unmodified module
    printf '1\n1156811761 2048 256 256 256 0\nnone\nout_lux0.txt\n' | build/ranlux_lags.exe
    printf '1\n1156811761 2048 256 256 256 3\nnone\nout_lux3.txt\n' | build/ranlux_lags.exe
    printf '2\n1156811761 20000000 0\npairs_lux0.txt\n' | build/ranlux_lags.exe
    printf '1\n1156811761 256 256 256 256 0\ncube.bin\ncube.txt\n' | build/ranlux_lags.exe
    micromamba run -n cosemu python3 analyse_cubes.py <dir with the cubes> results

Against the unmodified module (built with `ORIG=1`, which links a stand-in
for the missing `RanluxLevel`), `lux` = 0, 2 and 3 all give `nskip = 0` and
the same numbers (`results/orig_s1_lux*.txt`). Against the fixed module,
level 0 gives the same first eight values, the same row sum and the same
lag table as the unmodified module (`results/fix_s1_lux0.txt`).

Lag correlation of the Gaussian cube along i (rows of 2048; 134 million
values for levels 0 and 3, 34 million for 1, 2, 4; noise 0.0001-0.0002):

| level | seed | lag 12 | 14 | 16 | 18 | 20 | 22 | 24 | 26 | 28 | largest over lags 1-40 along j / k |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 (the runs) | 1156811761 | -0.0027 | -0.0072 | -0.0116 | -0.0083 | +0.0051 | +0.0149 | +0.0088 | +0.0011 | +0.0022 | 0.0008 / 0.0005 |
| 0 (the runs) | 2009256176 | -0.0026 | -0.0071 | -0.0116 | -0.0087 | +0.0052 | +0.0148 | +0.0089 | +0.0012 | +0.0023 | 0.0006 / 0.0007 |
| review's table, seed 1156811761 | | -0.0025 | -0.0072 | -0.0116 | -0.0086 | +0.0053 | +0.0150 | +0.0089 | +0.0011 | +0.0025 | |
| 1 | 1156811761 | +0.0003 | -0.0002 | -0.0001 | -0.0002 | -0.0002 | +0.0001 | -0.0002 | +0.0002 | -0.0001 | 0.0016 / 0.0011 |
| 2 | 1156811761 | +0.0000 | -0.0004 | -0.0001 | +0.0001 | +0.0000 | +0.0000 | -0.0000 | -0.0001 | -0.0002 | 0.0021 / 0.0010 |
| 3 (new default) | 1156811761 | -0.0002 | +0.0000 | -0.0000 | +0.0001 | -0.0000 | -0.0000 | -0.0000 | +0.0000 | -0.0000 | 0.0005 / 0.0006 |
| 3 (new default) | 2009256176 | -0.0000 | +0.0001 | -0.0001 | -0.0000 | -0.0000 | -0.0002 | -0.0001 | +0.0001 | -0.0001 | 0.0006 / 0.0005 |
| 4 | 1156811761 | -0.0003 | -0.0002 | -0.0000 | +0.0001 | -0.0003 | +0.0001 | +0.0002 | +0.0001 | +0.0002 | 0.0016 / 0.0008 |

Odd lags are within noise at every level. At level 0 the largest value over
lags 1-40 is 0.0149 (lag 22) for both seeds; at levels 1-4 it is 0.0002-0.0004.
The review's table (both seeds in its JSON) agrees with the level-0 rows to
4 x 10^-4, the noise of two samples of the same streams.

Variance of the white field along x, fractional amplitude of the harmonics
of the box (per cent). Rows: the 1-D synthesis of the 2048-long rows; the
256^3 cube synthesised along the three axes (`results/analyse_cubes.json`):

| harmonic m | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | noise (m = 40-120) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| level 0, seed 1156811761, rows of 2048 | 0.14 | 0.53 | 1.43 | 2.30 | 1.65 | 1.01 | 2.95 | 1.74 | 0.22 | 0.43 | 0.02 |
| level 0, seed 2009256176, rows of 2048 | 0.14 | 0.52 | 1.41 | 2.31 | 1.71 | 1.03 | 2.93 | 1.75 | 0.23 | 0.46 | 0.02 |
| level 3, seed 1156811761, rows of 2048 | 0.01 | 0.04 | 0.01 | 0.02 | 0.04 | 0.01 | 0.01 | 0.02 | 0.02 | 0.01 | 0.02 |
| level 3, seed 2009256176, rows of 2048 | 0.02 | 0.01 | 0.03 | 0.02 | 0.02 | 0.03 | 0.04 | 0.03 | 0.02 | 0.04 | 0.02 |
| level 0, seed 1156811761, 256^3 cube, along x | 0.06 | 0.53 | 1.37 | 2.13 | 1.54 | 0.90 | 2.77 | 1.54 | 0.13 | 0.43 | 0.06 |
| level 0, seed 1156811761, 256^3 cube, along y / z (largest of m = 1-24) | | | | | | | | | | | 0.14 / 0.14 |
| level 3, seed 1156811761, 256^3 cube, along x (largest of m = 1-24) | | | | | | | | | | | 0.13 |

The profile along x is the same for the two seeds at level 0 (correlation
0.975 between the 2048-point profiles) and unrelated at level 3 (-0.02).
The review measured 0.5-2.6 per cent at harmonics 6-12 on the rebuilt
linear field of boxes 1 and 2 with the same six harmonics on top. The
amplitude of harmonic m is 2 c(2m) (1 - 2m/N), with c the lag correlation
of the table above and N the row length (the review checked m = 5-14 of
both rows): the 256^3 cube rows are 6-8 per cent below the 2048 rows
because the lag pairs that fall off the end of a shorter row are lost. The
2048-row values are the ones for the L = 1024 runs; for the L = 1800 runs
(rows of 3200) the factor is 0.993 at m = 11.

Power of the squared white field f^2 at the axis modes (+-m, 0, 0), over the
mean power of the non-axis modes with |n| within 0.5 of m, 256^3 cube
(homogeneous noise: about 1, with the scatter of a two-mode mean):

| m | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 |
|---|---|---|---|---|---|---|---|---|---|---|
| level 0, seed 1156811761, x axis | 0.6 | 55 | 403 | 942 | 534 | 164 | 1604 | 521 | 3.3 | 39 |
| level 0, seed 2009256176, x axis | 3.5 | 41 | 382 | 1061 | 512 | 188 | 1506 | 536 | 5.7 | 40 |
| level 0, seed 1156811761, y axis / z axis | 3.4 / 3.5 | 1.2 / 0.1 | 0.9 / 0.0 | 1.0 / 1.0 | 1.9 / 0.0 | 0.1 / 0.0 | 0.6 / 0.3 | 0.9 / 1.2 | 0.4 / 3.0 | 1.4 / 0.2 |
| level 3, seed 1156811761, x axis | 0.2 | 1.9 | 0.7 | 1.0 | 1.4 | 0.8 | 0.9 | 0.4 | 0.0 | 0.1 |

The production executables `PMP2init`, `PMP2start`, `PMP2main`, `PMP2MG`,
`ic2pm` and `PMP2analyze` build from the patched sources (serial `make`,
ifx 2024.2.0, the makefile's flags; `PMP2spectrum` does not build on `cz`
either, a target that lacks `PMP2linker.o`). `PMP2startMPI` compiles with
Intel MPI 2021.13 (`source /cosma/local/intel/oneAPI_2024.2.0/mpi/latest/env/vars.sh;
I_MPI_F90=ifx mpiifort <the makefile's MPIFLAGS> -c PMP2startMPI.f90`; the
`mpi` module of `intel_comp/2024.2.0` conflicts with the login shell's
OpenMPI) and was not linked or run. The switch was exercised: with the
variable unset the level is 3, with `GLAM_RANLUX_LUXURY=0` it is 0 and the
numbers are the old ones, with an invalid value the default is used and a
message printed.

End to end, on the initial conditions of the 128^3 example in `test/`
(box 2, seed 2009256176; the review of 2026-10-19, `cmp` of every output
file of the run):

| comparison | result |
|---|---|
| unmodified build against this branch at `GLAM_RANLUX_LUXURY=0` | every file bit-identical (`PMcrd.DAT`, `PMcrs0.DAT`, the power spectra, `Results.log`, `fort.17`, `pt.dat`) |
| this branch, default against `GLAM_RANLUX_LUXURY=0` | the particle files and the spectra differ |
| this branch, default, 4 threads against 1 thread | bit-identical |
| this branch, `GLAM_RANLUX_LUXURY=3` against unset | bit-identical |

## Consequences for the existing runs

Nothing can be re-run now. Every run with native initial conditions
(`PMP2start`, `PMP2startMPI`), whatever the seed, box or resolution, has the
same pattern: the variance of the linear field is modulated along x at
harmonics 6-12 of the box, by 0.5-3.0 per cent of the variance (the
2048-row values of the table above, 2.95 at harmonic 11; the same in every
box), and not along y or z. The mean of the field and its shell-averaged
power spectrum are not affected: the mean power of every mode is right.
The pattern is a two-point statistic that depends on position: the
coefficients of the x frequencies m and m + Delta m are correlated, so the
variance of the field changes along x. It is invisible in anything averaged
over translations and appears at first order in the squared field, at the
axis modes and in any statistic that is not so averaged (the variance of a
slab or a sub-box, a light-cone cut along x, the covariance of distinct
modes). The suites, from the launch scripts of the run folders (the
review): every one calls `PMP2start.exe`; none calls `PMP2startMPI.exe`;
the `fid2LPTIC_*` runs go through `ic2pm.exe` and are not affected.

| run folders (native generator, affected) | box, particles |
|---|---|
| `mg_glam/DurMun_hmfemu_{LCDM,fRn1}_wide_sample_first_64_model1-64` and the supplement models 65-72 of this checkout | 1024, 2048^3 |
| `mg_glam/DEGRACE2_fid_LCDM_L1024` (the fiducial boxes), `DEGRACE2_fid_LCDM_L1024Np4096Ng8192` | 1024, 2048^3 and 4096^3 |
| `mg_glam/F5n1_L1024Np2048Ng4096`, `F6n1_L1024Np2048Ng4096` | 1024, 2048^3 |
| `DEGRACE_{LCDM,F4n1,F5n1}_L1800Np3200Ng6400` | 1800, 3200^3 |
| `fid_F5`, `fid_F5n1_L512Np2048Ng4096`, `fid_LCDM`, `fid_LCDM_L1024Np4096Ng8192`, `fid_LCDM_L512Np2048Ng4096`, `fid_LCDM_L512Np2560Ng5120`, `LCDM_L1024Np2048Ng4096`, `test` | 512-1024 |
| `fid2LPTIC_L512Np2048Ng4096` and its `_da4`, `_da6` (2LPTIC through `ic2pm`): not affected | 512, 2048^3 |

From the review of 2026-10-19 on the 95 fiducial boxes (L = 1024 Mpc/h,
z = 0.25):

| product | effect |
|---|---|
| spectra of delta^2 and s^2 of the initial conditions at the lattice shells n^2 = 49, 64, 81, 100, 121, 144 (the axis modes (+-j, 0, 0), j = 7-12) | +37, +775, +25, +30, +94, +80 per cent over the homogeneous expectation; every other shell 0.0 +- 1.8 |
| halo power at n^2 = 64 (six axis modes), every density level | +20 to +67 per cent over the neighbouring shells |
| raw halo correlation function of one box at 29.6-40.5 Mpc/h | +0.1 to +0.6 per cent from these modes |
| linear field (1-1, 1-delta, delta-delta, L-L) at these shells | within 1-2 per cent |

A user of the runs should know: statistics that average over all modes of
a shell or over the box are affected well below 1 per cent (power spectra
in bins: six axis modes among about 800 of a unit shell at |n|^2 = 64,
measured within 2 per cent on the linear field and +20 to +67 per cent on
the halo field at those six modes, 0.15-0.5 per cent of the bin; the
correlation function at 10.4-21.3 Mpc/h, measured +0.04 to +0.22 per cent;
mass functions, inferred, not measured: a 1-3 per cent modulation of the
variance in slabs along x changes the box's mass function at second order,
of order 10^-4);
statistics that are sensitive to the variance of the field at wavelengths
L/12 to L/6 along x (85-170 Mpc/h for L = 1024), to the axis modes, or to
the quadratic fields (control variates, bispectra at the axis modes, the
correlation function beyond 30 Mpc/h at the 0.1-0.6 per cent level) carry
it; the pattern is the same in every box, so it does not average out over
realisations. The L = 1800 DEGRACE suites were not measured here; the
generator and the draw order are the same, so the pattern is the same in
index space (harmonics 6-12 of the box).

Other products of the code: the on-the-fly power spectra and the halo
catalogues are measurements on the evolved field, so they carry what the
field carries (above); the measurement routines themselves do not use the
generator. The galaxy biasing of `PMP2mod_analyze.f90` (`getRandom`,
uniform numbers per row) and `getGauss` of `PMP2mod_density.f90` used level
0 as well; uniform numbers at level 0 have no two-point defect, Gaussian
numbers from `getGauss` have the one above. With this branch both use level
3 unless `GLAM_RANLUX_LUXURY=0` is set, so their draws change.

## Not established

- The lag correlations were measured, not derived in closed form; the
  first-order part (lags 12-24) is explained by the one relation, the
  second-order part (lags 26-30) only named.
- The size of the effect on the evolved matter field of the runs (the
  review measured the initial quadratic fields and the halo field) and on
  the halo catalogues beyond the halo power at |n|^2 = 64.
- The L = 1800 runs were not measured directly.
- That `iFlag` was 0 in the production binaries is inferred from the
  review's checks and from a fresh build of the unmodified code (bit-identical
  to level 0 on the 128^3 example), not read from the binaries.
- The end-to-end reproduction of a production 2048^3 initial condition at
  level 0: the runs' particle files were deleted by their submit scripts;
  the evidence is the 128^3 identity above and the agreement of the
  halobias review's rebuilt low-k modes with the runs' halo fields.
- The MPI generator at run time (compiled only, not linked or run).
- `getGauss`, `getRandom` and `PMP2correlation` are exercised by no test;
  the first two are off in every `Setup.dat` of the run folders (no
  redshift-space spectra, no biasing model), so no DEGRACE product used
  them.
- The downstream replica (above) has not been rebuilt; its executable of
  2026-10-14 is the one in use.

## Record

| item | value |
|---|---|
| worktree | `/cosma8/data/dp203/dc-ruan1/scratch/mgglam-ic-white-noise-20261019` |
| compiler | ifx 2024.2.0 (`intel_comp/2024.2.0`, `compiler-rt tbb compiler`), the makefile's flags |
| test runs | login node, one core, `nice -n 19`; 1.3-20 s and 28-142 MB each (mode 1), 2-5 s and 312 MB (mode 2), 53 s and 839 MB (`analyse_cubes.py`) |
| no Slurm job | |
| review | blind review of 2026-10-19 (`docs/reviews/ic-white-noise/`: three findings files, `TRIAGE.md`); its fixes are in the second commit of the branch |
| files kept | `test/ic_white_noise/results/`: the text outputs of every run above and `analyse_cubes.json`; the cubes (67 MB each) were deleted |
| file names in `results/` | `orig_*`: against the unmodified module; `fix_*`, `cube256_*`: against the fixed module; `s1` seed 1156811761, `s2` seed 2009256176, `luxN` the requested level, `pairs` mode 2, `env_*` the switch (unset, 0, invalid); the `lag` lines hold lag, correlation along i, j, k |
