# REVIEW 1: the Fortran fix (`PMP2mod_random.f90`, its callers, `luxuryp.h`, the builds)

Branch `ic-white-noise-20261019`, fix commit e604768 on `cz` = 1024c1e. Reviewed in the
detached worktrees `wt` (branch, HEAD 54c90ff) and `cz` (1024c1e) under
`/cosma8/data/dp203/dc-ruan1/scratch/mgglam-ic-review-20261019/`; builds and runs only there
and under `work/`. Line numbers refer to the branch unless marked `cz`. Reviewer: a fresh
session that took no part in the work.

## 1. Summary

The cause is as stated: on `cz` the table `ndskip` is set nowhere (the `DATA` statement at
`cz:PMP2mod_random.f90:167` is a comment), the unmodified module reports `nskip = 0` after
`rluxgo(2, ...)`, and all four callers asked for level 2. The fix is small, correct and
complete for the programs the makefile's `all` builds; it introduces no new race. The central
claim of the switch holds end to end: the `cz` build and the fixed build with
`GLAM_RANLUX_LUXURY=0` write bit-identical initial conditions and bit-identical secondary
outputs for the 128^3 example, the default differs, 1 and 4 threads agree at level 3, and
`=3` equals unset. One caller of `rluxgo` is left on the literal `lux = 2`
(`PMP2correlation.f90`), the environment-variable parsing is permissive, and `PMP2startMPI`
was compile-checked here but not linked or run.

Merge recommendation: ready after the listed fixes (none is a Blocker; R1-1 and R1-2 are
the ones worth doing before the merge).

## 2. Findings

| ID | Severity | File:line | Description | Evidence | Suggested fix |
|---|---|---|---|---|---|
| R1-1 | Minor (pre-existing caller, behaviour changed by the commit) | `PMP2correlation.f90:1031-1034`, `:1087-1089` | The only other caller of `rluxgo` keeps the literal `lux = 2` and draws its random particle subsample with `ranlux(Gg,Nrecord)`, `Gg(ip) < fraction`. With the table now set it draws at level 2; before the commit it drew at level 0. `GLAM_RANLUX_LUXURY` does not reach it, so a correlation measured with this program on an old run selects a different random subsample than the same program did before the branch, and the old subsample cannot be recovered. The program is not in `all`; no launch script of the run folders names it (`grep` of `mg_glam/*/*.sh` and the main checkout's `*/arun*.sh`). | `git grep -n -i 'rluxgo' -- '*.f90'` lists five call sites; four changed, this one not. | `lux = RanluxLevel()` at line 1032, or a sentence in the note that this program is out of scope. |
| R1-2 | Minor (introduced) | `PMP2mod_random.f90:138,145-146` | `ival` is not initialised before `Read (str, *, iostat=istat) ival`; a list-directed read is permissive. Measured with `ranlux_lags.exe` mode 1, `lux = -1` (table in section 3): `3.5`, `3,4`, `3 4` give 3 and `1e3` gives 1000 silently; `/` terminates the read with `istat = 0` and leaves `ival` undefined (here the message appeared; elsewhere it may not). The note says an invalid value prints a message; for these inputs it does not. | section 3, table "parsing" | `ival = -1` before the read; accept only strings of digits (optional sign), or document the permissive read. |
| R1-3 | Nit (introduced) | `PMP2mod_random.f90:152-153`; `test/ic_white_noise/ranlux_lags.f90:156,186` | The message prints p, not the level, for the 24-2000 form: `GLAM_RANLUX_LUXURY=24` prints `RANLUX luxury level = 24` although it is level 0 and reproduces the old draws (`luxlev_set 0 nskip_used 0`). The test's `i3` formats print `***` for p >= 1000. | section 3 | Print both the level and p (`luxlev`, `nskip + 24`) after the first `rluxgo`, or say "p =" in the message. |
| R1-4 | Minor (gap in the evidence, disclosed by the note; partly closed here) | `PMP2startMPI.f90:83,85` | Not built by the author. Compile-only check here passes (object written, no diagnostics) with Intel MPI 2021.13: `source /cosma/local/intel/oneAPI_2024.2.0/mpi/latest/env/vars.sh; I_MPI_F90=ifx mpiifort <MPIFLAGS> -c PMP2startMPI.f90`. The `mpi` module of `intel_comp/2024.2.0` cannot be loaded in a login shell here (conflict with the profile's `openmpi`), which is why the author did not build it. Not linked, not run. Each rank prints the level once (`levelSaved` is per process); the variable reaches every rank only if the launcher exports the environment (Hydra and Slurm do by default; not tested). | `work/warn/PMP2startMPI.o` | Link and run the MPI generator once on a compute node (2 ranks, the example) before relying on it; say in the note how the MPI module is obtained. |
| R1-5 | Nit (pre-existing, not the branch's) | `makefile:87` | `PMP2spectrum` does not build on `cz` (`PMP2mod_analyze.f90(22): error #7002 ... [LINKERLIST]`, the target lacks `PMP2linker.o`) nor on the branch (after `PMP2main` has made the modules: undefined `linkerlist_mp_bdm_`, `extradofbackgrounddata_mp_backgroundevolution_`). | `work/warn/cz_spectrum_build.log` | Out of scope; add the missing objects to the target when the program is next needed. |
| R1-6 | Nit (introduced) | `PMP2mod_random.f90:282-286`, `luxuryp.h:18` | `ndskip` is written into the shared common block `/LUX2/` on every `rluxgo` call. Not a new race: every call site is serial (`SetRandomN`, `BiasDensity`'s loop) or inside `!$OMP critical` (`getGauss` in `DENSITrsd`), every call writes the same constants, and `/LUX2/` was already written on every call before the branch (`nskip`, `next`, `twom24`, `twom12`, `inseed`, `luxlev`, `notyet`). Recorded so nobody re-reviews it. | code reading; `PMP2mod_density.f90:468-470`; `PMP2mod_analyze.f90:243-262` has no OpenMP directive | none |
| R1-7 | Nit (introduced) | `PMP2mod_random.f90` (whole), `work/warn/warn_all.log` | `-warn all` gives 119 warnings, all #6717 (implicit typing) of the old style; the new lines inherit it (`lux` at 122, `ndskip` at 282, `luxlev` at 288). No warning of another kind. | `grep -o 'warning #[0-9]*' work/warn/warn_all.log | sort | uniq -c` -> `119 warning #6717` | none |

What the commit discloses as not established and this review confirms as still open: `iFlag` = 0
in the production binaries (section 3 shows it was 0 in this build of `cz`); the MPI generator
at run time.

The user's decisions (not counted): environment variable rather than a `Setup.dat` line;
default 3 rather than 2; branch of the fork rather than upstream.

## 3. Verified OK

| property | how | result |
|---|---|---|
| The commented `DATA ndskip` is the only place the table could be set on `cz` | `git -C cz grep -n -i ndskip -- '*.f90' '*.h'` | three hits: the comment at `PMP2mod_random.f90:167`, the read at 254 and 258, the declaration in `luxuryp.h:17`; no assignment |
| `ndskip` and `nskip` are zero after `rluxgo(2, ...)` against the unmodified module | reviewer program `work/mech/mech_cz.exe` (built against `cz`), prints the common block right after `rluxgo(2, 1156811761, 0, 0)` | `luxlev = 2  nskip = 0  ndskip = 0 0 0 0 0`; the committed test gives `orig_s1_lux2.txt: luxlev_set 2 nskip_used 0`, `orig_s1_lux3.txt: nskip_used 0` |
| Every caller asked for 2 on `cz` | `cz:PMP2start.f90:83`, `cz:PMP2startMPI.f90:84`, `cz:PMP2mod_random.f90:122`, `cz:PMP2mod_density.f90:397`, `PMP2correlation.f90:1032` | all `lux = 2` (the last unchanged, R1-1) |
| Every use of `rluxgo`, `ranlux`, `GAUSS3`, `RanluxLevel` in the tree | `git grep -n -i 'rluxgo\|call ranlux\|GAUSS3\|RanluxLevel' -- '*.f90'` | `rluxgo`: `PMP2start.f90:86`, `PMP2startMPI.f90:87`, `PMP2mod_random.f90:123`, `PMP2mod_density.f90:398`, `PMP2correlation.f90:1034` (a commented one at 770); `ranlux`: `GAUSS3`/`GAUSS3fake`, `getRandom`, `PMP2correlation.f90:1087`; `GAUSS3`: `SPECTR` of both generators, `getGauss`, the test; `RanluxLevel`: the four callers and the test |
| Programs the makefile builds with the changed module | `make -j1` in `wt` | `PMP2init`, `PMP2start`, `PMP2main`, `PMP2MG`, `ic2pm`, `PMP2analyze` build and link (9 s for the first two, 48 s for the rest); `PMP2spectrum` fails on both trees (R1-5); `PMP2startMPI` compiles (R1-4); `PMP2BDM`, `PMP2correlation`, `PMP2density`, `PMP2mainMPI` not tried |
| Streams are complete under skipping (what `ranlux` needs with `nskip > 0`: `in24`, `kount`, `carry`, `i24`, `j24`, `seeds` are saved per plane in `/LUX1/` (`threadprivate`); `nskip`, `next`, `twom24`, `twom12` are shared constants set identically by every `rluxgo`) | code reading of `SetRandomN` 85-99, `SPECTR` 306-317, `luxuryp.h:15-23`; direct test below (1 against 4 threads) | bit-identical |
| End-to-end reproduction of the old draws (128^3 example, box 2, seed 2009256176, `OMP_NUM_THREADS=4` unless stated, `OMP_STACKSIZE=512M`, `nice -n 19`) | `cmp` of every output file of `Run2/` | see table below |
| `GLAM_RANLUX_LUXURY` parsing (`ranlux_lags.exe` mode 1, `lux = -1`, 2 planes x 2 rows of 64) | one run per value | see table below |
| `luxuryp.h` constants | `maxlev = 4`, `lxdflt = 3`, `ndskip(0:maxlev)`, table 0, 24, 73, 199, 365 | consistent with `RanluxLevel`'s default 3, with `rluxgo`'s mapping of p = 24, 48, 97, 223, 389 back to levels 0-4 (measured: 24 -> level 0, 389 -> level 4), and with the note |
| `RanluxLevel` thread safety | named critical `ranlux_level` around the one-time read; the only call from inside a parallel region is `getGauss`, itself inside an unnamed critical (nesting of differently named criticals is allowed); the read of `levelSaved` after `END CRITICAL` follows the implied flush | sound |
| `iFlag` on `cz` is an undefined implicit local of `SetRandomN` (not in any common block of `luxuryp.h`) | `cz:PMP2start.f90:75-100`, `luxuryp.h` | confirmed; in this build it behaved as 0 (the `cz` run equals the fixed run at level 0, which sets it) |
| Cost of level 3 | 223/24 = 9.3 numbers generated per number delivered; the example's `Spectr` wall time; the test's 134 million Gaussians (REVIEW 2) | table below |

End-to-end comparisons (`work/ex_*/Run2`):

| comparison | `PMcrd.DAT` | `PMcrs0.DAT` (48 MB) | `PowerDM.R`, `PowerDM.log`, `Results.log`, `fort.17`, `pt.dat` |
|---|---|---|---|
| (a) `cz` build vs fixed build, `GLAM_RANLUX_LUXURY=0` | identical | identical | all identical |
| (b) fixed build default vs `=0` | differ (byte 79) | differ (byte 2) | differ except `pt.dat` |
| (c) fixed build default, 4 threads vs 1 thread | identical | identical | - |
| (d) fixed build `=3` vs unset | identical | identical | - |
| (e) level line printed | `RANLUX luxury level = 0/3 (GLAM_RANLUX_LUXURY; ...)` once, on standard output only | no output file records the level or the seed (`PMcrd.DAT` header carries the `Setup.dat` line; `Setup.dat` has `Random seed = 0`; the seed comes from the realisation's row of `TableSeeds.dat`) | |

Wall time and peak memory of the example (`/usr/bin/time -v`):

| run | threads | level | wall s | `Spectr` s | peak RSS MB |
|---|---|---|---|---|---|
| `cz` build | 4 | 0 | 1.26 | 0.16 | 147 |
| fixed, `=0` | 4 | 0 | 1.18 | 0.11 | 147 |
| fixed, unset | 4 | 3 | 1.32 | 0.16 | 147 |
| fixed, unset | 1 | 3 | 3.68 | 0.61 | 146 |
| fixed, `=3` | 4 | 3 | 1.39 | 0.16 | 147 |

The standalone test (one core) draws 134 million Gaussians in 6.85 s at level 0 and 20.3 s at
level 3 including its statistics; scaled by 64 to 2048^3 that is 14 extra minutes of one core
at level 3, shared by the threads of `SPECTR`, as the note says.

`GLAM_RANLUX_LUXURY` parsing (value -> level used; "msg" = the "not 0-4 or 24-2000" message):

| value | result | value | result |
|---|---|---|---|
| unset | 3 | `23` | msg, 3 |
| `3 ` | 3 | `24` | p = 24 = level 0, `nskip 0`; prints "level = 24" |
| `03` | 3 | `389` | p = 389 = level 4 |
| `3.5` | 3, no msg | `2000` | p = 2000, `nskip 1976` |
| `-1` | msg, 3 | `2001` | msg, 3 |
| `5` | msg, 3 | empty | 3, no msg |
| ` ` (blank) | msg, 3 | `abc`, `0x3`, `T`, `4abc` | msg, 3 |
| `3,4`, `3 4` | 3, no msg | `1e3` | p = 1000, no msg |
| `+2` | 2 | `/` | msg, 3 here (`ival` undefined) |

## 4. Open questions for the author

1. Is `PMP2correlation` in use anywhere? If yes, should it take `RanluxLevel()` (R1-1)?
2. Should the level (and p) be written into a run output that survives (`PMcrd.DAT` header
   extras, `Results.log`), given that nothing but standard output records it? (REVIEW 3 takes
   this up.)
3. Was `iFlag` ever observed non-zero? A non-zero value would zero the first Gaussian of every
   plane, i.e. the (0, 0, kz) modes; the halobias review's linear-field checks argue against it
   for the L = 1024 runs, but nothing has been checked for the L = 1800 runs.
4. The MPI generator: is it still meant to be used, and can it be linked and run once on a
   compute node?

## 5. Runs made

All on the login node, `nice -n 19`, one process at a time, `uptime` 1-min load 8-27 on
128 cores before each start (one spike to 87 from other users during the REVIEW 2 reruns),
`free -g` available >= 1.4 TB.

| step | command | wall | peak RSS |
|---|---|---|---|
| builds in `wt` | `make -j1 PMP2init PMP2start`; `make -j1 PMP2main PMP2MG ic2pm PMP2analyze PMP2spectrum` | 9 s; 48 s (last target fails) | - |
| builds in `cz` | `make -j1 PMP2init PMP2start`; `make -j1 PMP2spectrum` (fails) | 9 s | - |
| warnings | `ifx <FFLAGS> -warn all -module work/warn -c PMP2mod_random.f90 -o work/warn/PMP2mod_random.o` | 2 s | - |
| MPI compile-only | `mpiifort <MPIFLAGS> -c PMP2startMPI.f90 -o work/warn/PMP2startMPI.o` (Intel MPI 2021.13 env sourced) | 3 s | - |
| five example runs | `OMP_NUM_THREADS=N nice -n 19 /usr/bin/time -v <build>/PMP2start.exe <<< 2` in `work/ex_*/Run2` | 1.2-3.7 s each | 146-147 MB |
| parsing tests | 22 runs of `ranlux_lags.exe` mode 1 at 64 x 2 x 2 | 0.1 s each | 10 MB |

Deleted at the end: the PM files `PMcrs0.DAT` (48 MB each, five copies) and `PMcrd.DAT`
of `work/ex_*/Run2`. Left in place: the executables and `.o`/`.mod` files of both worktrees
(`PMP2init.exe`, `PMP2start.exe` in `cz`; those plus `PMP2main.exe`, `PMP2MG.exe`,
`ic2pm.exe`, `PMP2analyze.exe` in `wt`), `wt/test/ic_white_noise/build/` and `build_cz/`
(both ignored by git: `git status --short` is empty), `work/` (2.5 MB: logs, text outputs,
the reviewer's programs). `git diff --stat` is empty in both worktrees. No process of the
reviewer is left (`ps -u`); a VS Code `agent host` process of the user's own session is
running (not touched).
