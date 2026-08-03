# Tests: exact-redshift snapshot outputs (`Nexact`/`zexact`)

Automated test suite for the precise-redshift snapshot feature. A tiny LCDM
box (36³ particles, 72³ mesh, 150 Mpc/h, z_init = 49) is evolved for 100
steps under several output configurations, and every produced
`PMcrd.NNNN.DAT` header is checked against an independent Real*4 replication
of the schedule arithmetic (`asserts.py`, pure-stdlib Python).

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
| `case1_merge` | `zexact` equal to an existing output moment (z = 0.49337, the legacy `zout=0.5` nearest-step moment) | merged: **one** snapshot, no duplicate, no index shift; header bit-exact at 1/(1+z) |
| `case2_insert` | `zexact = 0.512`, between two scheduled steps | a **new** timestep is inserted; snapshot lands bit-exactly on 1/(1+z); legacy moments preserved |
| `case3_multi` | `zexact = 2.5 0.512 1.1` (deliberately unsorted) | all three produced, snapshot numbers ascending as z decreases; z = 2.5 exercises the merge path, the others insert |
| `case4a/4b` | z = 60 (≥ z_init) / z = −0.5 in `Init.dat` | `PMP2init` stops with a clear message; no `Setup.dat` written |
| `case4c/4d` | pair closer than 0.5% / z = 60, hand-edited into `Setup.dat` | `PMP2main` stops with a clear message; no snapshots written |
| `case5a_legacy` | no `Nexact` anywhere | legacy nearest-step schedule, output moments match the ladder oracle |
| `case5b_stripped` | `Setup.dat` with the trailing block **removed** (pre-feature format) | accepted; expansion-factor sequence identical to `case5a` step by step |
| `case6_tailinsert` | `zexact = 0`, `da = 6.9e-4`: the z=0 step must be **inserted after** the step where the legacy half-step exit rule (`a >= 1 - da/2`) already fires | the run keeps stepping to the inserted moment (`NlastX` stop index) and writes the snapshot at AEXPN exactly 1.0, ending there without a duplicate final dump |

Checks common to every simulation case: the *set* of numbered snapshots
equals the oracle's prediction exactly (marked moments + the end-of-run
`iSave` dump), requested epochs are **bit-exact** `float32(1/(1+z))`, and all
other epochs match the ladder to ≤ 2 ulp.

## Notes

- The oracle in `asserts.py` mirrors `Initialize`/`SetExactSteps` in
  `PMP2main.f90` including the float32 rounding of every operation; if the
  scheduling code changes, the oracle must change with it (a mismatch shows
  up as a hard FAIL, not a silent pass).
- `ZOUT_LEGACY`, `NSTEPS`, `z_init`, `da` are mirrored between
  `Init.base.dat`, `run_tests.sh`, and `asserts.py` — keep them in sync.
- Halo catalogs share the same trigger path as particle snapshots (the
  `Nlist` block in the main loop); BDM is disabled here (`Find BDM halos = 0`)
  because the box is too small to host halos at the test redshifts. The halo
  path was validated separately on a larger box (see feature commit message).
- P(k) and RNG-seed tables are borrowed from `../PkTable.dat` /
  `../TableSeeds.dat` (fallback: `fid2LPTIC_L512Np2048Ng4096/`, or a
  synthetic seed table); the test only exercises scheduling, so the P(k)
  normalization is irrelevant.
