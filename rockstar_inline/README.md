# In-line Rockstar halo finding for MG-GLAM

This directory holds an adapter. At an output step it finds haloes with the Rockstar
halo finder (Behroozi et al. 2013), working directly on MG-GLAM's particle arrays in
memory, and writes a catalogue next to BDM's. The adapter is licensed under the GNU
GPL, version 3 or later, because it links GPLv3 code from
[MPI-Rockstar](https://github.com/Tomoaki-Ishiyama/mpi-rockstar) (Tokuue et al.
2024). A GLAM build without the adapter links `rs_stub.c` instead, which contains no
Rockstar code.

## What it does

Each call (`rs_inline_run`, declared in `rs_inline.h`) runs these steps:

1. **FoF.** An exact, periodic friends-of-friends with linking length b = 0.28 over
   all particles (`rs_fof.c`): parallel, deterministic, int64 throughout. Particles are
   friends when d² < r² in float, as in Rockstar's own test. Groups smaller than
   `MIN_HALO_PARTICLES` are dropped.
2. **Gather.** Each group's members are copied into Rockstar's particle array, sorted
   by particle id and unwrapped around the first member. Positions are converted to
   comoving Mpc/h and velocities to physical peculiar km/s (`rs_convert.h`). The
   velocities use the epoch a_v = AEXPN − ASTEP/2, because GLAM's momenta lag its
   positions by half a leapfrog step.
3. **Kernel.** MPI-Rockstar's per-group kernel (`find_subs`) runs on every group in an
   OpenMP loop, largest group first, as MPI-Rockstar's `do_halo_finding` does. The
   haloes are merged in group order, so the catalogue does not depend on the thread
   count.
4. **Output.**
   - Rockstar's output filter is applied (`MIN_HALO_OUTPUT_SIZE`, `UNBOUND_THRESHOLD`).
   - Strict spherical-overdensity masses are computed from every particle
     (`rs_strict.c`, as Rockstar's `STRICT_SO_MASSES`).
   - Parents are assigned as Rockstar's `find_parents` does (`rs_parents.c`).
   - The catalogue is published atomically (`rs_write.c`).

The particle arrays are only read. Nothing in the adapter ends the process:
- Rockstar's `exit()` and failed assertions become C++ exceptions (`shadow/assert.h`)
  and are caught for each group (`rs_guard.cpp`).
- Allocation failures and failed memory pre-checks become a status.
- In every case a catalogue is still published, and its `.meta` says what happened.

## Build

Use the Intel toolchain, as for GLAM, so that GLAM and the adapter share one OpenMP
runtime:

```
module load intel_comp/2024.2.0; module load compiler-rt tbb compiler mpi
make -C rockstar_inline ROCKSTAR_DIR=<MPI-Rockstar checkout>   # build/librs_inline.a + harness/rs_harness
make -C rockstar_inline stub                                    # build/librs_stub.a
```

The adapter has been validated against MPI-Rockstar commit 0ae386f. The commit and the
compiler flags are recorded in every catalogue's `.meta`.

## Configuration

The configuration is read from `Rockstar.config` in the run directory, with Rockstar's
own parser. `FILE_FORMAT` must be present, although it is not used.

**Values the adapter sets from the run** (they are not read from the file):
- `BOX_SIZE`, `TOTAL_PARTICLES`, `Om`, `Ol`, `h0`, `SCALE_NOW`;
- `PARTICLE_MASS` = ρ_c Ω_m L³/N, and the mean particle spacing.

**Values the adapter forces:**
- `TEMPORAL_HALO_FINDING = 0`, `PERIODIC = 1`, `LIGHTCONE = 0`, `OUTPUT_LEVELS = 0`.

**`FORCE_RES`.** It is capped at `FORCE_RES_PHYS_MAX`/a, as Rockstar does. Its value
matters at GLAM's resolution. At FORCE_RES = 0.25 Mpc/h (one PM cell at 0.5 Mpc/h mean
particle spacing), Rockstar removes every host below ~2×10¹² Msun/h. Host counts are
converged for FORCE_RES ≤ 0.1 Mpc/h.

**Keys for GLAM:**
- `GLAM_MEMORY_CEILING_GIB` (default 900): the output is skipped when the memory
  pre-checks would exceed it.
- `GLAM_STRICT_SO` (default 1).

## Catalogue

`CATALOGS/Rockstar.<step>.<box>.bin` holds raw little-endian records of 376 bytes.
Each record is:
- Rockstar's `struct halo`, 312 bytes, with `OUTPUT_RVMAX`, `OUTPUT_NFW_CHI2` and
  `OUTPUT_INERTIA_TENSOR`;
- then `pid`, `upid`, the phase-space parent, the FoF group, and the strict SO masses
  with flags saying which thresholds were reached.

The sidecar `.meta` (key = value) gives:
- every column with its type, byte offset and unit;
- the effective configuration and the velocity epoch;
- FoF statistics, timings and memory;
- the status: `ok`, `skipped` or `failed`.

Catalogues that are skipped or failed have zero rows. `harness/rsbin.py` reads the
catalogue into a numpy structured array.

## Tests (`harness/`)

`harness/rs_harness` reads a GLAM PM snapshot (from a run with `Save snapshots = 1`).
Its commands:

| command | what it does |
|---|---|
| `inline DIR STEP CONFIG OUTDIR [--shift F] [--fail kernel\|assert\|alloc] [--repeat N]` | runs the in-line finder exactly as GLAM calls it |
| `fofcmp DIR STEP` | test T2: the parallel FoF against the serial reference FoF (`fof_ref.c`), checking that the two partitions are exactly equal |
| `internal DIR STEP OUTDIR NFILES [...]` | writes Rockstar INTERNAL particle files plus a config for standalone MPI-Rockstar (the reference) |

The Python scripts compare the in-line catalogue with standalone MPI-Rockstar run on
the same particles with the same configuration (`t1_compare.py`, `t1_props.py`). They
also check units against BDM (`t4_units.py`).

**Validation** (L128 and L256 boxes at 0.5 Mpc/h particle spacing, ΛCDM and f(R),
z = 0 and 1):
- **FoF:** identical partitions to the serial reference.
- **Determinism:** catalogues are bitwise identical between runs with 3 and 4 threads,
  and on reruns.
- **Error paths:** each injected failure is reported, and a clean call afterwards
  gives the same catalogue.
- **Mass functions and properties against standalone MPI-Rockstar:**
  - host and subhalo mass functions agree within Poisson errors (0.1-1 %);
  - median Vmax, Rs, spin, b/a and T/|U| agree within 0.3 %;
  - subhalo fractions agree as well.

  Standalone runs that differ only in rank count or input order differ from each
  other by the same amounts.
- **Halo-by-halo identity.** 92-93 % of hosts with ≥ 100 particles match a standalone
  host within 0.1 R200m, and 89 % match within |Δln M| < 0.02. Between standalone runs
  the figures are 94-95 % and 92 %.
  - The difference comes from the FoF. Rockstar's own FoF merges groups within 2r of
    particles with more than 16 neighbours, and does so in an order that depends on
    the run.
  - A deterministic emulation of that merging lowered the agreement, so the exact FoF
    is used.
- **Velocities against BDM:** bulk velocities agree with BDM's to the expected epoch
  factor AEXPN/a_v.
