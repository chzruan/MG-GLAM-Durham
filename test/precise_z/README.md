# Tests: exact-redshift snapshot outputs (`#outputs < 0`)

Automated test suite for the precise-redshift snapshot feature. A tiny LCDM
box (36³ particles, 72³ mesh, 150 Mpc/h, z_init = 49) is evolved for 100
steps under several output configurations, and every produced
`PMcrd.NNNN.DAT` header is checked against an independent Real*4 replication
of the schedule arithmetic (`asserts.py`, pure-stdlib Python).

## Config convention under test

The existing `#outputs` line + redshift list in `Init.dat`/`Setup.dat` —
no new config lines:

```
#outputs =         3        (<0: hit these redshifts exactly)
 1.0 0.5 0.0
```

- `#outputs = 3` — legacy: analyze at the scheduled step *closest* to each z.
- `#outputs = -3` — exact: the same three redshifts are hit *exactly*; each
  becomes a new timestep, or a scheduled step within 0.5% (in 1+z) that is
  not already pinned to another target is moved onto it. Targets may be
  arbitrarily close to each other (only bitwise duplicates are rejected);
  the later member of a close pair becomes a short inserted step.
  `PMP2init` then writes the list to `Setup.dat` in full Real*4 precision
  (`es16.8`) instead of the legacy `f8.3`.

## Run

```bash
cd test/precise_z
bash run_tests.sh            # loads the Intel modules itself; ~1-2 min
# or, with the environment already loaded:
SKIP_MODULES=1 bash run_tests.sh
```

Exit code 0 = all checks passed. All artifacts land in `work/` (safe to
delete). The executables are built automatically if missing.

## Cases

| case | configuration | expectation |
|------|---------------|-------------|
| `case1_merge` | `-2` / `2.00 0.49337` (0.49337 is within 0.5% of a scheduled step) | merged: the step is *moved* onto the target — one snapshot, no duplicate step; header bit-exact at 1/(1+z) |
| `case2_insert` | `-2` / `2.00 0.512` (both between scheduled steps) | two new timesteps inserted; snapshots land bit-exactly on 1/(1+z) |
| `case3_multi` | `-3` / `2.5 0.512 1.1` (deliberately unsorted) | all three produced, snapshot numbers ascending as z decreases; z=2.5 exercises the merge path, the others insert |
| `case4a/4b` | `-1` / z = 60 (≥ z_init) resp. z = −0.5 | `PMP2init` stops with a clear message; no `Setup.dat` written |
| `case4c/4d` | duplicate z resp. z = 60, hand-edited into `Setup.dat` | `PMP2main` stops with a clear message; no snapshots written |
| `case5a_legacy` | `+2` / `2.00 0.50` | legacy nearest-step schedule, output moments match the ladder oracle |
| `case5b_trailing` | `case5a` `Setup.dat` with stale trailing lines appended | accepted and ignored; expansion-factor sequence identical to `case5a` step by step |
| `case7_closepairs` | `-8` / the four DEGRACE pairs closer than 0.5% in a (`2.0 1.98541 1.75 1.74475 1.00202 1.0 0.37562 0.375`) | all 8 pinned on distinct steps, headers bit-exact, the run executes every short step in the oracle's order |
| `case8_pinned` | `-2` / `0.49337 0.4920` (first merges onto a ladder step, second is 0.09% later) | the second target is inserted; the step pinned to the first is not moved again |
| `case9_classes` | `-3` / `2.5 0.512 1.1` with `Output classes = 3` / `2 1 2` appended to `Init.dat` | schedule identical to `case3_multi`; `OutputSchedule.dat` lists the three pinned steps with classes 2, 1, 2 attached to the right redshifts after sorting |
| `case4e_classes` | class count 2 for 3 outputs | `PMP2init` stops ("one output class per listed redshift"); no `Setup.dat` |
| `case4f_badopt` | `Extra outputs = 1  (1 = on)`: an optional line whose comment contains `=` | `PMP2init` stops ("malformed optional line") instead of silently taking the default; no `Setup.dat` |
| every simulation case | `Run1/OutputSchedule.dat` | exactly the marked steps; epochs equal the headers; exact flags; χ(a) matches flat matter+Λ to 1e-5; classes 1 unless requested |
| `case10_lens72` / `case10_lens36` | `-3` / `2.5 0.512 1.1`, `Lens slabs per axis = 2`, plane grid 72 resp. 36, `Dump random fraction = 1` for the `ProjectDENSIT` reference, `DM power spectrum = 1` at grid 72 | planes and `.meta` for every output; per axis the slab means average to 1; at grid 72 the slab average equals `ProjectionYZ/XZ/XY` to 1e-5; the 36 planes are the 2×2 block means of the 72 planes; L3: the plane power P_2D(k_⊥) equals the exact periodic-box projection of `PowerDM.R` (Σ_kz P_3D \|W_slab\|²/L) within 3σ per ring and 10 % on average for 0.08 < k_⊥ < 0.5 |
| `case11_modes` | `-3` / `2.5 0.512 1.1`, `DM power spectrum = 1`, `Extra outputs = 1`, `kc_delta = 0.5`, `kc_mom = 0.185` (below k_N(N_dep)/4), `N_dep_mom = 36`, `Subsample M = 1` | `deltam` (n = 11), `mom_x/y/z` (n = 4) cubes and `.meta` for every output; Hermitian iz = 0 plane; k = 0 zero; P(k) from the `deltam` cube equals `PowerDM.R` bin by bin to 1e-4 (V2); `Subsample` files in full mode (n_sel = N), whose decoded particles reproduce the `deltam` (5 bins) and `mom_x` (3 bins) spectra to 1e-3 (V6) and equal both dumps mode by mode to 1e-2 of the rms mode, while transposed or conjugated cubes fail (layout, axis order and sign of Im of the production dumps). Unit tests: `test/extras/run_modes_test.sh` (V1 layout, f_R interpolation) and `test/extras/run_subsample_test.sh` (V5: hash, order, quantisation, clipping) |
| `case12_frhalo` | `-2` / `1.1 0.512`, `MG_flag = 1` (F5), `Find BDM halos = 1`, `Extra outputs = 1`, `fR at halos = 1`, `test/BDM.config` | f(R) run with the halo finder: `fR_halo.<step>.0001.{bin,meta}` with one (fR_centre, fR_mean) pair per `CatshortV` row, finite, negative like f_R0. The interpolation itself is unit-tested in `test/extras/run_modes_test.sh` (`check_frhalo.py`) |
| `case12_frhalo` (binary catalogue) | same run, `BDM.config` with `CatalogueFormat = both` appended | at every BDM output `CatshortV.<step>.0001.{bin,meta}` sit beside the text catalogue: meta rows = .bin size/96 = text rows, meta `a` = header A and `a_v = a − da/2`, all 24 binary columns equal the text at its printed precision, `row` = 1..N, `distinct` = 0, no staged `.Catshort*` temporaries; read through `python/glamio` (validates its reader against the writer) |
| `case13_catbin` | `-1` / `0.512`, `Find BDM halos = 1`, `BDM.config` with `CatalogueFormat = binary` | `.bin` + `.meta` only (no `.DAT`), readable from the stem, `text_catalogue = none`, finite positive masses, `row` = 1..N |
| `case15_cattext` | the `case13` run with `CatalogueFormat = text` | `.DAT` only, no `.bin`/`.meta`; "BDM.config validated at start-up" in the log; the text rows equal `case13`'s binary rows column by column (bitwise-reproducible runs) |
| `case14_catdefault` | `-2` / `1.1 0.512`, `Find BDM halos = 1`, **no** `BDM.config` | the finder writes its default file at the first catalogue step with the full key (`CatalogueFormat = binary`, the compiled default since 2026-10-18) and accepts it at the second; `.bin` + `.meta` at both outputs, readable with glamio, `text_catalogue = none`, `row` = 1..N, no `.DAT`, no staged temporaries; "BDM.config absent" in the log, no configuration error (regression test for the `T20` tab stop that truncated the key to `Catalogue`) |
| `case4g_badcatfmt` | `BDM.config` with `CatalogueFormat = hdf5` | `PMP2main` stops at start-up (`CheckBdmConfig`, "CatalogueFormat must be text, binary or both") before any time step; no `CATALOGS/Catshort*`, no hidden `.Catshort*`, no snapshot written |
| `case6_tailinsert` | `-1` / `0.00` with `da = 6.9e-4`: the z=0 step must be **inserted after** the step where the legacy half-step exit rule (`a >= 1 - da/2`) already fires | the run keeps stepping to the inserted moment (`NlastX` stop index) and writes the snapshot at AEXPN exactly 1.0 |
| `glamio` cross-checks | the outputs of `case10_lens72`, `case11_modes`, `case12_frhalo` | every `python/glamio` reader agrees with this file's independent decoders: schedule rows and Om/OmL; lens planes equal, `plane('x', 1)` orientation, axis means 1; the 12 mode cubes equal, Hermitian accessor, `power()` = `binned_power`; PowerDM tables and epoch lines; subsample decode and ids, PMcrd header fields, every PMcrs particle equal to the subsample within the quanta; fR_halo pairs and `catalogue()` rows |

Checks common to every simulation case: the *set* of numbered snapshots
equals the oracle's prediction exactly (marked moments + the end-of-run
`iSave` dump), requested epochs are **bit-exact** `float32(1/(1+z))`, and all
other epochs match the ladder to ≤ 2 ulp.

## Notes

- The oracle in `asserts.py` mirrors `Initialize`/`SetExactSteps` in
  `PMP2main.f90` including the float32 rounding of every operation; if the
  scheduling code changes, the oracle must change with it (a mismatch shows
  up as a hard FAIL, not a silent pass). The bitwise assertions have real
  discriminating power: they caught the legacy `f8.3` truncation of the
  `Setup.dat` redshift list, which silently degraded exact targets.
- `ZINIT`, `DA0`, `NSTEPS`, and the per-case redshift lists are mirrored
  between `run_tests.sh` and `asserts.py` — keep them in sync.
- Halo catalogs share the same trigger path as particle snapshots (the
  `Nlist` block in the main loop); BDM is disabled here (`Find BDM halos = 0`)
  because the box is too small to host halos at the test redshifts. The halo
  path was validated separately on a larger box (see feature commit message).
- P(k) and RNG-seed tables are borrowed from `../PkTable.dat` /
  `../TableSeeds.dat` (fallback: `fid2LPTIC_L512Np2048Ng4096/`, or a
  synthetic seed table); the test only exercises scheduling, so the P(k)
  normalization is irrelevant.
