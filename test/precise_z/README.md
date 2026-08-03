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
  becomes a new timestep, or a scheduled step within 0.5% (in 1+z) is moved
  onto it. `PMP2init` then writes the list to `Setup.dat` in full Real*4
  precision (`es16.8`) instead of the legacy `f8.3`.

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
| `case4c/4d` | pair closer than 0.5% resp. z = 60, hand-edited into `Setup.dat` | `PMP2main` stops with a clear message; no snapshots written |
| `case5a_legacy` | `+2` / `2.00 0.50` | legacy nearest-step schedule, output moments match the ladder oracle |
| `case5b_trailing` | `case5a` `Setup.dat` with stale trailing lines appended | accepted and ignored; expansion-factor sequence identical to `case5a` step by step |
| `case6_tailinsert` | `-1` / `0.00` with `da = 6.9e-4`: the z=0 step must be **inserted after** the step where the legacy half-step exit rule (`a >= 1 - da/2`) already fires | the run keeps stepping to the inserted moment (`NlastX` stop index) and writes the snapshot at AEXPN exactly 1.0 |

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
