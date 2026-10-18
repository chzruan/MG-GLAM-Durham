# Triage of the five reviews of `glamio-readers` (2026-10-18)

Reviews run as independent subagents from `REVIEW-PROMPTS.md` against commit
0a6ecff; findings in `REVIEW-<k>-*-findings.md`. Reported totals: **Blocker 0,
Major 9 (seven distinct issues), Minor 45, Nit 24.** Recommendations: R1
"ready", R2, R3, R4, R5 "ready after the listed fixes". Every Major and most
Minors were fixed in the follow-up commit on the branch (this file is part of
it); the rest are listed below as deferred.

| review | Blocker | Major | Minor | Nit | scope |
|---|---|---|---|---|---|
| R1 snapshot, subsample, schedule, PowerDM vs the Fortran | 0 | 0 | 9 | 3 | bit-exact re-quantisation of a snapshot into the subsample |
| R2 lens planes, mode cubes, fR_halo vs the Fortran | 0 | 3 | 4 | 3 | independent reconstruction from the snapshots |
| R3 the package as software | 0 | 3 | 13 | 6 | API, robustness, login-node resources |
| R4 tests, executed in a worktree | 0 | 2 | 8 | 3 | 27-mutant battery, suite 221/221 reproduced |
| R5 docs, conventions, downstream | 0 | 1 | 11 | 9 | README, plan, deck, halocat needs |

What the reviews confirmed (no re-review needed): the 529-byte header
layout and byte order; the 24 MiB page records, their padding and the file
split; `a_v = AEXPN - ASTEP/2` at every GLAM snapshot and `ic2pm half` IC
(R1 re-quantised all 6 x 46656 values of a subsample from the snapshot with
zero mismatches); the lens-plane orientation on all three axes (R2 rebuilt
the planes from the particles); the mode-cube layout and the `exp(-i k.x)`
sign convention (conjugate and transpose alternatives are off by 8-13 rms);
`power()` = `PowerDM.R` structurally (one packed entry per full-sphere mode,
`iCIC = 0`, complete shells for `b <= n`); the fR_halo row order and the
catalogue join; splitmix64 against a pure-Python implementation and the
SplitMix64 reference vectors; the subsample decode to the last int16.

## Majors, verified and fixed

| ID | What | Fix |
|---|---|---|
| R2-M1 = R3-M1 (+ R2-m4) | `Modes.power()` held 9x the cube in temporaries (4.7 GB and ~12 s on the production n = 244 cube) | slab-by-slab `np.bincount` over the `ix` index; a few MB of temporaries, memory-mapped cubes read one slab at a time; `index_grids(sparse=True)`; results equal the old ones to 1e-13 |
| R2-M2 = R3-m7 | `power(nbins > n)` returned truncated shells 5-19 % low without warning | `ValueError` for `nbins > n`; docstring states that only shells `b <= n` are complete |
| R2-M3 = R4-1 = R3-m1 = R1-5 | `load_array(count=0, mmap=True)` crashed on a legitimate empty product (`rows = 0` fR_halo, `n_sel = 0` subsample) | returns an empty array for `count = 0`; unit tests for both empty products with and without `mmap` |
| R3-M2 = R5-m8 | `Snapshot.read()` held every page plus the concatenation and had no cap (412 GB for 2048^3 in physical units) | preallocated arrays filled page by page; default `max_particles = 2**27` (3.2 / 6.4 GB) with a message pointing to `iter_pages()`; README states 24 / 48 B per particle |
| R3-M3 = R1-6 (+ R3-m9) | `Subsample.ids()` allocated `arange(N_p)` (69 GB for 2048^3 at M = 1) before any check; `M = 0` kept every ID | count check first, `max_ids = 2**28` default cap, `M >= 1` validated in `keep_mask`, `select_ids` and `ids()`; cost documented |
| R4-3 (+ R1-3, R4-12) | the PMcrd unit-test fixture used the reader's own dtype, so a wrong field order was invisible; the suite pinned five header words | fixture packed with `struct`; every packed value is read back in the test; the suite pins Nrow, realization, ASTEP (= 2(a - a_v)), Om, OmL, h, PARTW and `AU0 = 0` against the subsample sidecar; `AEXP0`/`Nseed` documented as meaningless |
| R5-M1 = R2-m1 = R3-m2 | docstring, deck and plan claimed every lens plane averages to 1 | corrected to "per-axis mean 1; a plane averages to 1 + <delta>_slab" in the docstring, README, plan (I2 and L1) and deck (two lines); `plane_means()` added |

## Minors and Nits: handled

R1-1 `sync`-mode IC files (`ic2pm` header without the `v at a-da/2` tag,
`ISTEP = 0`) get `a_v = a` (`PMHeader.ic_velocity_epoch`); R1-2 the
records-per-file split is detected from the size of `PMcrs0` (512 for GLAM,
256 for PMP2start ICs); R1-3/R4-12 `AEXP0` and `Nseed` documented; R1-4/R1
Q5/R3-n4 the window helper renamed `cic_window_jing_sq` (Jing's alias-summed
1-D form, not a shell average), refused on `Z` files, `read_powerdm(box=)`
gives an exact `kNyq`; R1-7/R3-m8a the caption is searched after the epoch
line; R1-8/R3-n6/R5-m5 `lattice_tag` matches only the trailing `, IDx`
tag, `None` documented, `lattice` kept as an alias; R1-9/R3-m13/R5-m7
`particle_mass()` uses the finder's 2.774e11 (`RHO_CRIT_GLAM`); R1-11
`resolve_pair` on a directory raises `IsADirectoryError`; R1-12 stale
trailing records documented; R2-m2/R3-m6 `power()` without `Box` raises
`ValueError`; R2-m3/R3-m5 `catalogue(prefer=)` validated; R2-n1/R3-n2
boolean axes rejected; R2-n2 lens-plane size note; R2-n3 `fR_mean` fallback
wording; R3-m3 `read_modes` accepts only `MODE_PRODUCTS` by default; R3-m4
`product` returns the first token; R3-m8b schedule values validated
(`exact` 0/1, `class >= 1`, `a > 0`); R3-m10 `PowerDM.n_grid` alias;
R3-m12/R5-m11 helpers exported, README says to re-run the editable install
after a version bump (done, 0.2.1); R3-n1 unused import; R3-n3 empty
`read()` dtype; R4-2 README states the test invocations (script mode is the
contract; pytest without `-I`); R4-4 the schedule check covers case9
(classes 2 2 1) and case5a (exact 0 0); R4-5 `P.a`/`P.z` pinned; R4-6 a
non-default `SALT` fixture; R4-7 the wrap fixture sits at `NGRID + 1`; R4-8
exact-multiple and single-page fixtures, a 3-records-per-file split; R4-9/
R5-n7 `GLAMIO_TEST_TMPDIR` created on demand, `python/tmp/` ignored;
R4-10 the failing frame is appended to a `glamio` verdict; R4-11 asymmetric
Hermitian probes and the direct branch; R4-13/R5-n6 the extras README count
and the row placement; R5-m1 README lead paragraph qualified; R5-m2 the
example reads lens planes and P(k) at every output; R5-m3 deck and plan V2
wording (raw CIC); R5-m4 plan section 3 keys and Q9; R5-m6 `product_path()`
/ `products_at()` and `Schedule.row_at_z()`; R5-m9 the "still unreadable"
list; R5-m10 check_run.py stays an independent decoder (plan says so);
R5-n2/n3/n4/n8 README byte order, deck header-line count and sidecar epochs,
`ProjectionXY/XZ/YZ`.

## Minors and Nits: left as they are

R1-10 asterisk overflow of `i4` fields (unreachable: NstepM = 2000, Ng <=
6400); R3-m10 the mixed `__getitem__` conventions (`Schedule` by name or
index, `PowerDM` by name; documented); R3-m11 the catalogue module's own
size/layout code (reviewed in the previous round, left verbatim); R3-n5 no
`.bin`/`.meta` mtime check; R3-n6 `MetaProduct._int` rounding; R5-n1 fixed
in the README only; R5-n5 the deck's `h^{-1} Mpc` vs the package's `Mpc/h`
(same units); R5-n9 plural/singular column names (documented). Open: the
halocat ingest (a halocat session), an `M > 1` subsample on real suite
output (`test/extras/run_subsample_test.sh` covers the Fortran side and R4
confirmed glamio reproduces its M = 100/7 files), `*_deltaL` and `PowerDM.Z`
on real files, a `plane_means()`-style QA in `check_run.py`.

## Verification after the fixes

- `python/tests`: 16 tests (9 readers + 7 catalogue), script mode and pytest.
- `test/precise_z/run_tests.sh`: 221/221 with the strengthened `glamio` checks.
- `Modes.power()` on a synthetic n = 100 cube: tracemalloc peak a few MB
  beyond the data (was 9x); results equal the previous implementation.
- Deck rebuilt; the changed frames inspected.
