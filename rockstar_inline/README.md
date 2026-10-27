# In-line Rockstar halo finding for MG-GLAM

This directory holds an adapter. At an output step it finds haloes with the Rockstar
halo finder (Behroozi et al. 2013), working directly on MG-GLAM's particle arrays in
memory, and writes a catalogue of distinct haloes next to BDM's. Subhaloes are not
written: GLAM's fixed mesh resolution does not resolve substructure. The adapter is
licensed under the GNU
GPL, version 3 or later, because it links GPLv3 code from
[MPI-Rockstar](https://github.com/Tomoaki-Ishiyama/mpi-rockstar) (Tokuue et al.
2024). A GLAM build without the adapter links `rs_stub.c` instead, which contains no
Rockstar code.

## What it does

Each call (`rs_inline_run`, declared in `rs_inline.h`) runs these steps:

1. **FoF.** An exact, periodic friends-of-friends with linking length b = 0.28 over
   all particles (`rs_fof.c`): parallel and deterministic.
   - Particles are friends when d² < r² in float, as in Rockstar's own test, with
     d² = (dx² + dy²) + dz² and no fused multiply-add.
   - Indices and group sizes are int64; group labels are int32 (fewer than 2³¹ groups,
     checked).
   - Groups smaller than `MIN_HALO_PARTICLES` are dropped.
   - Particle coordinates that are not finite, or more than one box outside GLAM's
     range, fail the output instead of being moved.
2. **Gather.** Each group's members are copied into Rockstar's particle array, sorted
   by particle id and unwrapped around the first member. Positions are converted to
   comoving Mpc/h and velocities to physical peculiar km/s (`rs_convert.h`). The
   velocities use the epoch a_v = AEXPN − ASTEP/2, because GLAM's momenta lag its
   positions by half a leapfrog step. A group spanning half the box or more fails the
   output (it could not be unwrapped).
3. **Kernel.** MPI-Rockstar's per-group kernel (`find_subs`) runs on every group in an
   OpenMP loop, largest group first, as MPI-Rockstar's `do_halo_finding` does. The
   haloes are merged in group order, so the catalogue does not depend on the thread
   count. One failing group abandons the output: after an error inside `find_subs`,
   Rockstar's per-thread buffers cannot be reused until they are freed.
4. **Output.**
   - Halo centres are wrapped into [0, Box), and Rockstar's output filter is applied
     (`MIN_HALO_OUTPUT_SIZE`, `UNBOUND_THRESHOLD`).
   - Only distinct haloes are kept: those without a parent in Rockstar's
     `find_parents` (`rs_parents.c`), with the kernel radius `r`. Halos of equal vmax
     are taken in row order, so the selection does not depend on the C library's sort.
     The strict radius is not used here: a subhalo's strict sphere also counts its
     host's particles, which at L256 turned about 3 % of the distinct haloes into
     subhaloes.
   - Strict spherical-overdensity masses of the distinct haloes are computed from every
     particle (`rs_strict.c`, as Rockstar's `STRICT_SO_MASSES`).
   - The catalogue is published atomically (`rs_write.c`).

   The kernel itself is unchanged, and it still finds substructure: Rockstar separates
   neighbouring haloes in one FoF group through the same phase-space hierarchy, and it
   computes a halo's masses and properties from its own particles and those of its
   substructure (`num_child_particles`).

The particle arrays are only read. Failures inside the adapter do not end the process:
- Rockstar's `exit()` and failed assertions become C++ exceptions (`shadow/assert.h`)
  and are caught (`rs_guard.cpp`): per group in the kernel, and around the
  configuration parser and the search for distinct haloes.
- Allocation failures and failed memory pre-checks become a status.
- A catalogue is then still published, and its `.meta` says what happened. If a
  complete catalogue of the same step already exists (a restart that repeats the
  output), it is kept. If the file system refuses the write, only the run's standard
  output records the reason.

Signals (a segmentation fault) and the kernel's or Slurm's out-of-memory killer still
end the process. Under Linux's memory overcommit an allocation rarely fails, so the
memory pre-checks are the real protection against running out of memory: each output
is skipped when VmRSS plus the estimated need would exceed the ceiling, the smaller of
`GLAM_MEMORY_CEILING_GIB` and 97 % of the process's memory limit (its cgroups, such as
a Slurm job's, or MemTotal).

## In GLAM

- **Switch:** the optional `Init.dat` line `Find Rockstar halos = 1`, after
  `fR at halos` (`gen_inputs.py --rockstar` writes it); `PMP2init` copies it into
  the extra-outputs block of `Setup.dat`. Older `Setup.dat` files without the line
  read as 0.
- **Configuration:** `Run<box>/Rockstar.config`, next to `BDM.config`
  (`gen_inputs.py --rockstar` writes the validated settings).
- **Build:** GLAM's makefile links `build/librs_inline.a` with
  `make ROCKSTAR_DIR=<MPI-Rockstar source>`, and the stub otherwise.
- **Start-up check:** `PMP2main`/`PMP2MG` stop before any compute when the switch is
  on but the build has only the stub, `Rockstar.config` is missing, or it is not
  usable (`CheckRockstarConfig` in `PMP2mod_rockstar.f90`).
- **When it runs:** at class-1 outputs, after BDM (`Analysis` in
  `PMP2mod_analyze.f90`). The PM density array `FI` is freed for the call and
  allocated again afterwards: after BDM it is not resident, and `DENSIT` rebuilds it
  at the next step. The scalar field of f(R) runs (`FI2`) is not touched.
- **Particle arrays:** only read. With `GLAM_ROCKSTAR_CHECKSUM=1` in the environment,
  GLAM compares checksums of the six particle arrays before and after each call and
  stops if they differ (a test, not for production).
- **Not supported:** the MPI build switches the catalogues off (each process holds
  only part of the particles).
- **Post-processing:** `check_run.py` and glamio read and pack the catalogue (see
  `python/README.md`).

## Build

Use the Intel toolchain, as for GLAM, so that GLAM and the adapter share one OpenMP
runtime:

```
module load intel_comp/2024.2.0; module load compiler-rt tbb compiler mpi
make -C rockstar_inline ROCKSTAR_DIR=<MPI-Rockstar checkout>   # build/librs_inline.a, harness/rs_harness, harness/test_fof
make -C rockstar_inline stub                                    # build/librs_stub.a
```

A program that links `librs_inline.a` also needs `-lstdc++`. The adapter has been
validated against MPI-Rockstar commit 0ae386f. The commit (with `-dirty` for local
changes), the compiler version and the flags are recorded in every catalogue's
`.meta`. The failure-injection hooks of the tests are compiled into the harness only.

## Configuration

The configuration is read from `Rockstar.config` in the run directory, with Rockstar's
own parser. `FILE_FORMAT` must be present, although it is not used. The file is read
once per process (at the start-up check or the first output), so every output of a run
uses the same settings; a file that cannot be read or holds no entries is an error. The
`.meta` records the file's checksum (`config_fnv1a64`) and every effective Rockstar
setting (`config.<KEY>`).

**Values the adapter sets from the run** (they are not read from the file):
- `BOX_SIZE`, `TOTAL_PARTICLES`, `Om`, `Ol`, `h0`, `SCALE_NOW`;
- `PARTICLE_MASS` = ρ_c Ω_m L³/N, and the mean particle spacing.

**Values the adapter forces:**
- `TEMPORAL_HALO_FINDING = 0`, `PERIODIC = 1`, `LIGHTCONE = 0`, `OUTPUT_LEVELS = 0`,
  `PARALLEL_IO = 1`, `NON_COSMOLOGICAL = 0`, `W0 = -1`, `WA = 0`.

**Keys that have no effect in-line:** `STRICT_SO_MASSES` and `UNFILTERED_HALO_OUTPUT`
(a warning is printed when they are set; the strict pass is `GLAM_STRICT_SO`, and the
output filter always applies), and `OUTPUT_FORMAT`.

**Checked:** FORCE_RES > 0; 0 < FOF_LINKING_LENGTH ≤ 0.5; 0 < FOF_FRACTION < 1;
MIN_HALO_PARTICLES ≥ 2; MIN_HALO_OUTPUT_SIZE ≥ 1; 1 ≤ MIN_PART_PER_BIN ≤
MIN_SCALE_PART (otherwise Rockstar's NFW fit divides by zero); 1 ≤ MAX_SCALE_BINS ≤
1000.

**`FORCE_RES`.** It is capped at `FORCE_RES_PHYS_MAX`/a, as Rockstar does. Its value
matters at GLAM's resolution. At FORCE_RES = 0.25 Mpc/h (one PM cell at 0.5 Mpc/h mean
particle spacing), Rockstar removes every host below ~2×10¹² Msun/h. Host counts are
converged for FORCE_RES ≤ 0.1 Mpc/h.

**Keys for GLAM:**
- `GLAM_MEMORY_CEILING_GIB` (default 900): the output is skipped when the memory
  pre-checks would exceed it (or 97 % of the process's memory limit, if smaller).
- `GLAM_STRICT_SO` (default 1).

## Catalogue

`CATALOGS/Rockstar.<step>.<box>.bin` holds raw little-endian records of 352 bytes, one
per distinct halo. Each record is:
- Rockstar's `struct halo`, 312 bytes, with `OUTPUT_RVMAX`, `OUTPUT_NFW_CHI2` and
  `OUTPUT_INERTIA_TENSOR`;
- then the FoF group, and the strict SO masses with flags saying which thresholds were
  reached.

The sidecar `.meta` (key = value) follows the other GLAM sidecars:
- `product = Rockstar`, `status` (`ok`, `skipped` or `failed`) and `reason`;
- the column lists `columns`, `dtypes`, `offsets` and `units`, and one
  `column_NN = name: description` line per column;
- the epochs (`a`, `a_x`, `z`, `a_v`, `da`, `velocity_scale`), `step`,
  `realization`, the run (`Box`, `N_grid`, `N_row`, `N_p`, `particle_mass`), the
  cosmology (`Omega_m`, `Omega_L`, `hubble`), the gravity model (`MG_flag`,
  `MG_model`, `fR0`) and the run's `header`;
- the effective Rockstar configuration (`config.<KEY>`) and its file's checksum;
- FoF statistics, timings and memory (VmRSS sampled at each phase's peak, VmHWM).

Column notes:
- With `BOUND_PROPS = 1` (the default here), Rockstar computes every property except
  `m` from the bound particles: `r` is the radius of the bound mass `mgrav`, and the
  `alt_m_*` are bound masses.
- A mass of 0 in `alt_m_*` or `alt_m_strict_*` means that the threshold was not
  reached at this `FORCE_RES` (`strict_crossed` says which were); the `*2` shape and
  inertia columns (`MASS_DEFINITION4`, 500c by default) are then 0 too.

Catalogues that are skipped or failed have zero rows. An interrupted write leaves a
`.bin` without `.meta` and possibly the hidden temporaries
`.Rockstar.<step>.<box>.{bin,meta}.tmp`; a rerun of the step replaces them.
`glamio.read_rockstar` reads the catalogue, native or packed, into a numpy structured
array; it refuses a skipped or failed one unless asked (`python/README.md`, "In-line
Rockstar catalogues"). `harness/rsbin.py` is the test harness's minimal reader.

## Tests (`harness/`)

`harness/rs_harness` reads a GLAM PM snapshot (from a run with `Save snapshots = 1`).
Its commands:

| command | what it does |
|---|---|
| `inline DIR STEP CONFIG OUTDIR [--shift F] [--fail kernel\|assert\|alloc] [--repeat N]` | runs the in-line finder exactly as GLAM calls it |
| `fofcmp DIR STEP` | test T2: the parallel FoF against the serial reference FoF (`fof_ref.c`), checking that the two partitions are exactly equal |
| `internal DIR STEP OUTDIR NFILES [...]` | writes Rockstar INTERNAL particle files plus a config for standalone MPI-Rockstar (the reference) |

`harness/test_fof` needs no snapshot: it checks the FoF against an O(n²) brute force
with its own conversion and distance test, on adversarial particle sets at the
production geometry (pairs at r(1 ± 2×10⁻⁴) across the periodic boundary and across
cell boundaries, GLAM's coordinate edges, a corner cloud, a blob, a chain). It also
checks the canonical order, bitwise identical output with 1 to 4 threads, and the
error codes.

The Python scripts compare the in-line catalogue with standalone MPI-Rockstar run on
the same particles with the same configuration:
- `t1_compare.py` compares the mass functions of distinct haloes ("hosts") and
  matches them halo by halo.
  - The mass-function ratios get their errors from a paired bootstrap over
    subvolumes, because both catalogues hold mostly the same haloes.
  - Halo by halo, the matched fraction and the fraction within |Δln M| < 0.02 must
    each be within 3 points of the same figures between two standalone runs of the
    same output (`--baseline`; `--max-gap`).
  - The reference's hosts are its rows without a parent in a run with
    `STRICT_SO_MASSES = 0`, whose parents use the kernel radius as the in-line
    selection does.
  - The strict masses are compared with a `STRICT_SO_MASSES = 1` reference. Its hosts
    are taken from its non-strict twin, which holds the same haloes.
- `t1_props.py` lists property medians by mass (descriptive).
- `t4_units.py` compares each host's bulk velocity with the mean velocity of the
  snapshot particles inside its radius, converted independently by glamio; it also
  reports the ratio to BDM's velocities.

**Validation** (ΛCDM and f(R) boxes of 128 and 256 Mpc/h at 0.5 Mpc/h particle spacing;
L128 at z = 0 and 1, L256 at z = 0):
- **FoF:** `test_fof` passes. The partitions are identical to the serial reference
  (T2) at all six outputs, up to 534 098 groups.
- **Determinism:** the `.bin` is bitwise identical with 3, 4 and 16 threads, and on
  reruns.
- **Error paths:** an injected Rockstar `exit()`, failed assertion or allocation
  failure, an unreadable configuration and an invalid one each give status `failed`
  with the reason. A clean call afterwards in the same process gives the identical
  catalogue, and an existing complete catalogue survives a failing rerun.
- **Against standalone MPI-Rockstar** (the same particles and configuration; its
  distinct haloes, "hosts", are those without a parent in a `STRICT_SO_MASSES = 0` run;
  the strict masses are compared with the `STRICT_SO_MASSES = 1` run of the same haloes):
  - The mass functions of distinct haloes (bound mass, `m` and the strict mass) agree
    within the paired-bootstrap errors in 125 of 126 bins. The exception is L128 f(R)
    z = 0, 10¹³-3×10¹³: 542 against 551, 3.06σ.
  - Median Vmax, Rs, spin, b/a and T/|U| agree within the reference's statistical
    error in every mass bin: within 1.5 % at L128, 0.6 % at L256.
  - Halo by halo, the reference's distinct haloes with ≥ 100 particles matched within
    0.1 R200m, and of those the fraction within |Δln M| < 0.02:

    | output | in-line vs standalone | standalone vs standalone |
    |---|---|---|
    | L128 ΛCDM z = 0 | 92.6 %, 89.2 % | 94.5 %, 91.7 % |
    | L128 ΛCDM z = 1 | 90.2 %, 90.0 % | 92.8 %, 92.5 % |
    | L128 f(R) z = 0 | 93.2 %, 89.2 % | 94.9 %, 91.8 % |
    | L128 f(R) z = 1 | 91.3 %, 90.6 % | 93.4 %, 93.1 % |
    | L256 ΛCDM z = 0 | 92.4 %, 89.3 % | 94.3 %, 91.6 % |
    | L256 f(R) z = 0 | 92.9 %, 89.2 % | 94.7 %, 91.7 % |
    | L256 ΛCDM z = 0, positions shifted by L/2 | 92.0 %, 88.3 % | 93.8 %, 90.5 % (shifted standalone) |
    | strict masses, L128 / L256 ΛCDM z = 0 | 92.7 %, 89.1 % / 92.6 %, 89.2 % | 94.7 %, 91.6 % / 94.4 %, 91.5 % |

    The in-line catalogue is 1.7-2.7 points below the standalone baseline of the
    same output in every row. The criterion of `t1_compare.py` allows 3 points; a
    fixed threshold would not fit, because the baseline itself is lower at z = 1.
  - The difference comes from the FoF. Rockstar's own FoF merges groups within 2r of
    particles with more than 16 neighbours, in an order that depends on the run. A
    deterministic emulation of that merging lowered the agreement, so the exact FoF is
    used.
- **Units (T4):**
  - The haloes' bulk velocity divided by the mean velocity of the snapshot particles
    inside r is 1.0000-1.0001. Those particles were converted independently by glamio
    at a_v.
  - A conversion at AEXPN is excluded at more than 40σ (z = 0 and 1).
  - Against BDM, which converts with AEXPN, the ratio is 1.0014-1.0064, against the
    expected 1.0027-1.0078. BDM's halo velocity is a different estimator.
- **In GLAM** (test runs from 2LPT initial conditions: L128 ΛCDM and f(R) with outputs
  at z = 1, 0.5 and 0; L256 ΛCDM and f(R) at z = 1 and 0):
  - catalogues at every full output, `check_run.py` without failures, and packing
    and unpacking byte for byte;
  - the particle checksums are identical before and after every call;
  - with Rockstar on and off, P(k) differs by at most 7×10⁻⁶ (two runs with it off:
    6×10⁻⁶, GLAM's round-off), and the BDM catalogues hold the same number of haloes
    above 10¹², 10¹³ and 10¹⁴;
  - a run stopped after an output and restarted from the checkpoint before it
    replaces that output's catalogue, with nothing duplicated or left over;
  - at L256 on 32 threads a call takes 12-19 s, and its memory stays below the peak
    BDM has already set. VmRSS before the call grows by 0.4 GiB between the two
    outputs; `MALLOC_ARENA_MAX=1` removes that but makes the kernel 2.6-5 times
    slower, so it is not set.
- **Cost at L256** (134 M particles, 16 threads): 27-31 s per output; the peak is
  4.3-4.4 GiB above the entry VmRSS (34-35 B per particle), during the gather. The
  catalogue holds 223 287 (ΛCDM) and 234 230 (f(R)) distinct haloes, 79-82 MB.
