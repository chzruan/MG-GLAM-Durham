# glamio

Numpy-only Python readers for GLAM / MG-GLAM simulation outputs: the BDM halo
catalogues (legacy text and binary), the in-line Rockstar halo catalogues, the
output schedule, the lens planes, the
low-k mode cubes, the quantised particle subsample, the f(R) halo values, the
power spectra and the PM particle snapshots. The same readers read packed runs
(one HDF5 file per output step, see [Packed runs](#packed-runs-hdf5)), which
need h5py in addition. Every reader validates the file
against its sidecar (product name, required keys, byte count) and raises
`ValueError` instead of misreading. The precise_z suite cross-checks each reader
against independent decoders on real GLAM output
(`test/precise_z/asserts.py`, the `glamio` checks).

## Install

From the repository root, into the `cosemu` environment:

```
micromamba run -n cosemu python3 -m pip install -e python/
```

Only `numpy>=1.24` and Python >= 3.10 are required; packed runs also need
`h5py>=3` (the `hdf5` extra; cosemu has h5py 3.15). h5py is imported only
when a packed file is opened. Re-run the line after a version bump: the
editable install records the version at install time (`pip show glamio` lags
`glamio.__version__` otherwise).

Tests: `python3 tests/test_readers.py`, `python3 tests/test_catalogue.py`,
`python3 tests/test_pack.py` and `python3 tests/test_rockstar.py` (script mode;
`test_pack.py` and the packing tests of `test_rockstar.py` need h5py and skip
without it) or `python3 -m pytest -q tests` when pytest is available (it is not
part of cosemu; `python3 -I` hides a user-site pytest).

## Halo catalogues

GLAM's BDM finder writes, per output step and box, `CATALOGS/Catshort?.<step>.<box>.*`
(step and box 4-digit zero-padded; the letter after `Catshort` is `V`, `W`, `M`
or `A` for the overdensity definition, `V` = virial being the usual one):

* `.DAT`: legacy text, 8 header lines then one row per halo with 24 columns,
  Fortran format `(3f11.4,3x,3f10.2,1p,2g12.4,g12.5,26g12.4)`.
* `.bin` + `.meta`: headerless little-endian binary, one 96-byte record per
  halo (24 four-byte fields, same rows and order as the text file), plus a
  `key = value` text sidecar (`columns`, `dtypes`, `rows`, `record_bytes`, `a`,
  `z`, `a_v`, `velocity_scale`, `step`, `Box`, `N_grid`, `N_row`, `Omega_m`,
  `Omega_L`, `hubble`, `particle_mass`, `overdensity`, ...).

Since 2026-10-18 the finder writes the binary form by default
(`CatalogueFormat = binary` in `BDM.config`); `text` and `both` are opt-in.

Record layout (`glamio.CATALOGUE_DTYPE`, itemsize 96):

| # | name | dtype | unit / meaning |
|---|---|---|---|
| 1-3 | x, y, z | `<f4` | Mpc/h, comoving, in [0, Box) |
| 4-6 | vx, vy, vz | `<f4` | km/s, peculiar velocity as written (converted with 1/a) |
| 7 | Mbound | `<f4` | Msun/h, bound virial mass |
| 8 | Mtot | `<f4` | Msun/h, total mass inside Rvir |
| 9 | Rvir | `<f4` | kpc/h |
| 10 | Vrms | `<f4` | km/s |
| 11 | Vcirc | `<f4` | km/s, maximum circular velocity |
| 12 | row | `<i4` | 1-based row index in this file |
| 13 | Cvir | `<f4` | concentration |
| 14 | Npart | `<f4` | number of particles (a float in the code) |
| 15 | distinct | `<i4` | 0 = distinct halo |
| 16 | Xoff | `<f4` | centre-of-mass offset of the bound particles, in units of Rvir |
| 17 | virial_ratio | `<f4` | 2K/\|Ep\| - 1 |
| 18 | lambda | `<f4` | spin parameter |
| 19 | RadRMS | `<f4` | kpc/h |
| 20 | b_a | `<f4` | axis ratio b/a |
| 21 | c_a | `<f4` | axis ratio c/a |
| 22-24 | ax_x, ax_y, ax_z | `<f4` | major-axis direction |

### Usage

```python
import glamio

# .bin, .meta, .DAT or the bare stem all work; .DAT selects the text reader
cat = glamio.read_catalogue("Run1/CATALOGS/CatshortV.0181.0001")
len(cat), cat.columns              # number of haloes, column names
cat["Mbound"]                      # a column view (structured array field)
cat.a, cat.z, cat.a_v, cat.box, cat.step   # parsed from the meta (None if absent)
pos = cat.positions()              # (N, 3) float64, Mpc/h
vel = cat.velocities_at_av()       # (N, 3) float64, km/s at the velocity epoch
meta = glamio.read_meta("Run1/CATALOGS/CatshortV.0181.0001.meta")
```

The binary reader requires both `.bin` and `.meta`, checks the file size
against the 96-byte record (and `rows` when present) and refuses a sidecar whose
`columns`/`dtypes` differ from `CATALOGUE_DTYPE`, so a format change is noticed
rather than misread.

### Velocity epoch

The `vx, vy, vz` columns are the velocities as the finder wrote them: GLAM's
leapfrog velocities live at the half step `a_v = a - da/2`, but they were
converted to km/s with `1/a` of the output epoch `a`. `Catalogue.velocities_at_av()`
returns `v_written * a / a_v`, i.e. the velocities converted with `1/a_v`, their
own epoch. For text catalogues `a_v` is derived from the header (`A - Step/2`).

## In-line Rockstar catalogues

With `Find Rockstar halos = 1` in `Init.dat` (and `Run<box>/Rockstar.config`),
GLAM runs the in-line Rockstar adapter (`rockstar_inline/`, see its README) at
every class-1 output, after BDM, and writes the distinct haloes (no subhaloes)
to `CATALOGS/Rockstar.<step>.<box>.{bin,meta}`:

* `.bin`: headerless little-endian records of `record_bytes` (352) bytes,
  `rows` of them: Rockstar's `struct halo` (312 bytes; MPI-Rockstar built with
  `OUTPUT_RVMAX`, `OUTPUT_NFW_CHI2`, `OUTPUT_INERTIA_TENSOR`) followed by the
  adapter's `fof_group`, strict SO masses (`m_strict`, `r_strict`,
  `alt_m_strict_1..4`), `strict_crossed` and `pad`.
* `.meta`: a `key = value` sidecar that describes the layout itself:
  `columns`, `dtypes` (`i8`, `f4`, `i4`, little-endian), `offsets` (bytes into
  the record) and `units` are space-separated lists with one entry per column,
  and `column_NN = name: description` lines explain each column. It also holds
  `status`, `reason`, `rows`, `record_bytes`, the epochs (`a`, `a_x`, `z`,
  `a_v`, `da`), the run (`step`, `realization`, `Box`, `N_grid`, `N_row`,
  `N_p`, `Omega_m`, `Omega_L`, `hubble`, `MG_flag`, `MG_model`, and `fR0` in
  f(R) runs), `particle_mass`, `header`, `finder` (MPI-Rockstar commit,
  compiler and flags), the configuration file and its checksum, every effective
  Rockstar setting (`config.<KEY>`), FoF statistics, timings and memory.

The reader builds the record type from the sidecar, so columns are found by
name rather than through a layout fixed in glamio. It checks that
the lists have one entry per column, that every type is a little-endian
integer or float and that the fields lie inside the record without
overlapping.

Units (the `.meta` column descriptions are authoritative): positions comoving
Mpc/h in [0, Box); velocities (`vx, vy, vz`, `corevel_*`, `bulkvel_*`, `vmax`,
`vrms`) physical peculiar km/s at the velocity epoch `a_v = a - da/2` (the
adapter converts GLAM's momenta at `a_v`, so no `a/a_v` correction is needed,
unlike the BDM catalogue); radii comoving kpc/h; masses Msun/h; `id` is the
0-based row number. With `BOUND_PROPS = 1` (the default) every property except
`m` comes from the bound particles; a mass of 0 in `alt_m_*` or
`alt_m_strict_*` means the threshold was not reached (`strict_crossed` holds
the bits). `rockstar_inline/README.md` explains the columns and the validation.

Status: every class-1 output publishes a pair, also when the finder does not
succeed (unless the file system refuses the write, or a restart repeats an
output whose complete catalogue exists, which is kept). `status = ok` is the
catalogue. `status = skipped` (a memory pre-check would exceed the ceiling) or
`status = failed` (an error inside the adapter, e.g. an invalid
`Rockstar.config`) comes with a 0-byte `.bin`, `rows = 0`, the same column
lists and a `reason`. These are legitimate products: they pack, verify, unpack
and prune like any other. They are not empty volumes, so the reader refuses
them unless asked. An interrupted write leaves a `.bin` without `.meta` and
possibly the hidden temporaries `CATALOGS/.Rockstar.<step>.<box>.{bin,meta}.tmp`;
`glamio.pack` and `DEGRACE-companion/check_run.py` refuse both.

```python
import glamio

stem = glamio.product_path("Run0", "rockstar", 165, box=0)   # Run0/CATALOGS/Rockstar.0165.0000 (or .../Output.0165.0000.h5/Rockstar)
cat = glamio.read_rockstar(stem)            # .bin, .meta, stem or packed path; mmap=True maps the records
len(cat), cat.columns                       # number of distinct haloes, column names
cat["m"], cat["mgrav"], cat["vmax"]         # columns of the structured array cat.data
cat.units["r"], cat.descriptions["r"]       # 'kpc/h' and the column description of the .meta
pos, vel = cat.positions(), cat.velocities()   # (N, 3) float64: Mpc/h, and km/s at a_v
cat.a, cat.a_v, cat.box, cat.particle_mass, cat.meta["config.FORCE_RES"]

c = glamio.read_rockstar(glamio.product_path("Run0", "rockstar", 114, box=0), allow_status=True)
if not c.ok:
    print(c.status, c.reason)               # e.g. 'skipped', 'memory pre-check: ...'; len(c) == 0
```

Without `allow_status`, `read_rockstar` raises `ValueError` naming the status
and the reason of a skipped or failed catalogue. It always raises when the
`.bin` does not hold `rows x record_bytes` bytes, when the `id` column is not
`0..N-1` (byte-swapped, misaligned or corrupted file), for a status other than
`ok`, `skipped` or `failed`, and for a skipped or failed sidecar with rows.
`products_at` lists the catalogue as `'rockstar'`, and
`run_file(rundir, "Rockstar.config")` finds the configuration natively or in
`RunInfo`. `glamio.rockstar_dtype(meta)` returns the record type of a parsed
sidecar.

## Other products

The sidecar products (`LensPlanes`, `deltam`, `mom_x/y/z`, `Subsample`,
`fR_halo`, and the catalogues, these under `CATALOGS/`) share the file scheme
`<product>.<step>.<box>.{bin,meta}` in the run directory (`Run1/`): raw
little-endian binary without record markers plus a `key = value` sidecar.
Their readers accept the `.bin`, the `.meta` or the bare stem, take
`mmap=True` to memory-map the binary, and return objects that expose the
sidecar as `.meta` and the common scalars as properties (`a`, `z`, `a_v`,
`box`, `step`, `realization`, `n_grid`, `n_row`, `n_p`, `omega_m`,
`omega_l`, `hubble`, `header`; `None` when the key is absent). The three
other products differ: `OutputSchedule.dat` and `PowerDM.*.dat` are text
(no sidecar, no `mmap`), and the PM snapshot is a pair of big-endian Fortran
files whose header is returned as `Snapshot.header` (a `PMHeader`).
`glamio.product_path(rundir, product, step)` builds the file names and
`glamio.products_at(rundir, step)` lists what a step has.

| reader | file | returns |
|---|---|---|
| `read_schedule(path)` | `OutputSchedule.dat` | `Schedule`: structured `rows` (`step a z chi exact class`), `omega_m/omega_l` from the comment, `row(step)`, `primary()`, `stations()` |
| `read_lensplanes(stem)` | `LensPlanes.*` | `LensPlanes`: `planes[plane, j, i]` float32 (`3 N_slab` planes X1..XN, Y1..YN, Z1..ZN), `plane(axis, slab)` as `[i, j]` (`slab` 1-based), `slabs(axis)`, `axis_means()` (1 by construction), `plane_means()` (`1 + <delta>_slab`), `pixel_edges()`, `slab_edges()` |
| `read_modes(stem)` | `deltam.*`, `mom_x/y/z.*`, `*_deltaL.*` | `Modes`: `data[ix+n, iy+n, iz]` complex64, `mode(ix, iy, iz)` with Hermitian completion for `iz < 0`, `wavenumbers()`, `k_magnitude()`, `power()` (shell-averaged `Box^3 <|delta_k|^2>` for the complete shells `1..n`, slab by slab, a few MB of temporaries) |
| `read_subsample(stem)` | `Subsample.*` | `Subsample`: `q[6, n_sel]` int16, `positions()` Mpc/h, `velocities()` km/s at `a_v`, `ids()` (0-based IDs from the splitmix64 rule; 8 B per particle, capped at 2 GB by default), `lattice_tag` |
| `read_frhalo(stem)` | `fR_halo.*` | `FRHalo`: `data` (`fR_centre`, `fR_mean`) float32 per catalogue row, `catalogue()` reads the catalogue of the same step |
| `read_powerdm(path, box=None)` | `PowerDM.R/log.*.dat` (and `Z`, only written by the redshift-space branch) | `PowerDM`: `k`, `nmodes`, `pk` (raw CIC), `a`, `z`, `step`, `ngrid`, `k_nyq` (exact when `box` is given), `binning`; `pk_deconvolved()` with Jing's alias-summed 1-D CIC window |
| `read_pm_header(path)` | `PMcrd.<step>.DAT` | `PMHeader`: the 529-byte header record (byte order detected; `aexpn`, `astep`, `istep`, `nrow`, `ngrid`, `nparticles`, `extras`, `box`, `a_v`, `ic_velocity_epoch`, `velocity_scale()`, `particle_mass()`) |
| `read_snapshot(path)` | `PMcrd` + `PMcrs<i>` | `Snapshot`: `iter_pages()` streams 1024^2-particle pages, `read()` loads all into preallocated arrays (24 B per particle in `'code'` units, 48 B in `'physical'`; refuses more than 2^27 particles unless `max_particles` is raised), `units='code'` (grid units, code momenta) or `'physical'` (Mpc/h in [0, Box), km/s at `a_v`); the records-per-file split (512 for GLAM, 256 for PMP2start ICs) is detected |

```python
import glamio

sched = glamio.read_schedule("Run1/OutputSchedule.dat")
for step in sched.steps:                              # lens planes and P(k) exist at every output (classes 1 and 2)
    planes = glamio.read_lensplanes(glamio.product_path("Run1", "LensPlanes", step), mmap=True)
    sigma_x1 = planes.plane("x", 1)                   # [iy, iz], Sigma/Sigma_bar of the first X slab
    pk_file = glamio.read_powerdm(glamio.product_path("Run1", "PowerDM.R", step), box=planes.box)
for step in sched.primary()["step"]:                 # the extras and the catalogue only at class-1 outputs
    modes = glamio.read_modes(glamio.product_path("Run1", "deltam", step))
    k, pk, nmodes = modes.power()                     # equals PowerDM.R bin by bin (raw CIC)
    sub = glamio.read_subsample(glamio.product_path("Run1", "Subsample", step))
    pos, vel, ids = sub.positions(), sub.velocities(), sub.ids()
    print(glamio.products_at("Run1", step))           # {'LensPlanes': ..., 'deltam': ..., 'catalogue': ..., ...}

snap = glamio.read_snapshot("Run1/PMcrd.0181.DAT")    # header only; particles stream on demand
for x, y, z, vx, vy, vz in snap.iter_pages(units="physical"):
    ...                                               # one 1024^2-particle page at a time
```

Conventions the readers encode (from the Fortran writers):

* Lens planes: `(N_lens, N_lens)` Fortran column-major per plane, X planes
  indexed `(y, z)`, Y planes `(x, z)`, Z planes `(x, y)`; with 1-based
  indices as in the sidecar, pixel `i` covers `[(i-1), i) Box/N_lens` and
  slab `s` covers `[(s-1), s) Box/N_slab` (0-based numpy index `i` covers
  `[i, i+1)`). `planes[p]` is the raw `[j, i]` layout; `plane(axis, s)`
  transposes to `[i, j]`. Only the per-axis mean is 1 by construction; a
  single plane averages to `1 + <delta>_slab` (`plane_means()`).
* Mode cubes: `delta_k = N^-3 sum_x f(x) exp(-i k.x)`, `k = (2 pi/Box)(ix, iy, iz)`,
  cube `|ix|, |iy| <= n`, `0 <= iz <= n`, `k = 0` zeroed; `P(k) = Box^3 <|delta_k|^2>`
  over the full sphere (`iz > 0` modes count twice). Momentum cubes hold
  `(1 + delta) v_i` with `v` in km/s at `a_v`.
* Subsample: ID = array index - 1; kept when `splitmix64(ID xor SALT) mod M == 0`
  (`SALT = 0x5EED00DE6C0A2026`); `x = (q + 32768 + 0.5) Box/65536`,
  `v = s dv_quant` (`0.25` km/s), velocities at `a_v = a - da/2`.
* Snapshots: GLAM is built with `-convert big_endian` (the reader detects the
  byte order from the record marker); coordinates in grid units `[1, N_grid+1)`,
  `x = (X - 1) Box/N_grid`; `v = V 100 Box/N_grid / a_v` with `a_v = AEXPN - ASTEP/2`
  (the subsample's convention; `velocity_epoch="a"` gives the legacy catalogue
  conversion; an IC file from `ic2pm` in its legacy `sync` mode holds velocities
  at `AEXPN` and `PMHeader.a_v` says so). The last page of a `PMcrs` file is
  padded to 1024^2 particles; the reader drops the padding. `AEXP0` and `Nseed`
  in the header carry no information (an uninitialised word and the running
  RNG state). `particle_mass()` uses the finder's constant `2.774e11`.
* `PowerDM`: raw CIC density on the `N_grid` mesh, window not deconvolved, shot
  noise not subtracted; `pk_deconvolved()` divides by Jing's alias-summed 1-D
  CIC window `1 - (2/3) sin^2(pi k/2 k_Nyq)` evaluated at the bin `|k|` (not a
  shell average of the 3-D window; the two differ by up to 37 % near `k_Nyq`).

Helpers that are not readers (`keep_mask`, `select_ids`, `splitmix64`,
`product_path`, `products_at`, `run_file`, `rockstar_dtype`, `load_meta`,
`check_size`, `load_array`, `resolve_pair`, `meta_int`) are exported from the
package too. Still without a
reader: `ProjectionXY/XZ/YZ`, `timing.log`, `main.log`, the BDM `outputB`
logs.

## Packed runs (HDF5)

A finished run directory holds about 400-430 files (19 per full output in an
f(R) run, 17 in ΛCDM, 4 per station; about 40 % of them 1-2.5 kB `.meta`
sidecars). `glamio.pack` packs them, after the run, into one HDF5 file per
output step:

```
Run0/
├── Output.<step>.<box>.h5     one per output step (42 in a DEGRACE run)
├── RunInfo.<box>.h5           schedule, BDM.config, Rockstar.config, logs
└── packed.md5                 md5 of the .h5 files, to check a copy
```

Inside `Output.<step>.<box>.h5` every product is a group named after it
(`CatshortV`, `Rockstar`, `fR_halo`, `Subsample`, `deltam`, `mom_x/y/z`,
`LensPlanes`, `PowerDM/R`, `PowerDM/log`):

* `<group>/data` holds the `.bin` bytes unchanged: same little-endian type,
  one contiguous block (no chunking or compression), aligned to 4 KiB in files
  written since the 2026-10-24 review, with the shape the reader gives the file
  (`CatshortV/data` a 24-column compound table, `Rockstar/data` a `(rows,)`
  compound table with the sidecar's record type, offsets and record size kept
  (`(0,)` for a skipped or failed output), `Subsample/data`
  `(6, n_sel)` int16, `deltam/data` `(2n+1, 2n+1, n+1)` complex64 (an HDF5
  compound `{r, i}`), `LensPlanes/data` `(3 N_slab, N_lens, N_lens)` float32, ...);
* the `.meta` is kept verbatim in the group attribute `meta_text` (the
  authoritative copy); the other group attributes repeat each `key = value`
  line as `glamio.read_meta` splits it, comments included (a line with two
  keys, such as `MG_flag = 1   MG_model = 1`, stays one attribute);
* text files are verbatim `uint8` datasets (`PowerDM/R/text`,
  `CatshortV/outputB`, `logs/<name>`, `BDM.config` and `Rockstar.config` in
  RunInfo);
  `PowerDM/<kind>/table` and `schedule/table` repeat the parsed columns;
* root attributes: `format = glamio-pack`, `format_version`, `layout` (a
  pointer to this section), the schedule row (`step a z chi exact class`),
  `products`, the run-wide sidecar values (`Box`, `N_grid`, `N_row`, `N_p`,
  `Omega_m`, `Omega_L`, `hubble`, `MG_flag`, `MG_model`, `fR0`, which is the
  string `none` in ΛCDM runs, and `header`), `realization`, `run`, `box_dir`,
  `source_dir`, `glamio_version`, `glamio_code_sha256`, `packed_utc` and, when
  `<run dir>/bin/git_hash.txt` exists, `glam_build`; RunInfo has `steps` (the
  output steps with products, one Output file each);
* `manifest` lists every original file: `path` (relative to `Run<box>/`),
  `bytes`, `md5` and `object`, the dataset holding the bytes or, when no
  dataset has that name, `<group>@<attribute>` (the sidecars). `unpack`
  rebuilds the files byte for byte from it. File times and permissions are
  not kept.

Files are written with HDF5 objects readable by HDF5 1.10 and later (the
2026-10-24 files also read with 1.8.20). `format_version` changes only when an
older reader would misread a new file; glamio refuses newer versions.

```
micromamba run -n cosemu python3 -m glamio.pack pack   Run0 [--out DIR]   # write the packed files (originals untouched)
micromamba run -n cosemu python3 -m glamio.pack verify Run0 [--against Run0]   # --against only while the originals exist
micromamba run -n cosemu python3 -m glamio.pack unpack Run0 --out DIR     # rebuild the original files (never overwrites)
micromamba run -n cosemu python3 -m glamio.pack prune  Run0 [--yes]       # delete the originals after re-checking them
micromamba run -n cosemu python3 -m glamio.pack status Run0               # native, incomplete, packed, partially-pruned, pruned, ...
```

`pack` refuses a run directory with anything it cannot place: another
subdirectory, a symbolic link, PM particle or checkpoint files, an unknown
step-tagged file, a text catalogue (`CATALOGS/*.DAT`), a `.bin` without its
`.meta`, the temporaries of an interrupted Rockstar write
(`CATALOGS/.Rockstar.*.tmp`), a wrong box tag, a step missing from
`OutputSchedule.dat`, or a `packed.md5` it did not write. Each file is written as `<name>.partial`,
verified (each dataset's bytes in the file against the md5 of its original,
plus the parsed tables and attributes), flushed to disk and renamed; a failure
removes the `.partial`. A rerun keeps an existing packed file only if its
manifest lists exactly the current originals of its step with the recorded
md5s, and refuses otherwise. `prune` refuses while a `.partial` exists or
while any original is in no packed file, re-verifies every packed file and
every original still present, deletes nothing if anything differs, and
without `--yes` only reports. One core and at most about 1 GB of memory; an
L1024 Np2048 box (75 GB) takes about 20 min to pack and 30 min to prune, so
run it under Slurm: `DEGRACE-companion/pack_run.sh` does pack, a full byte
check and then prune (unless `PRUNE=no`), and resumes from `status`.

Reading: the sidecar-product, binary-catalogue (BDM and Rockstar), P(k) and
schedule readers take `<file>.h5/<group>` where they take a native path (`mmap=True` maps the
dataset at its offset in the HDF5 file; the text catalogue and snapshot
readers are native only). `product_path` / `products_at` return such paths
when the native product is absent or incomplete; the schedule is found with
any `box` when the directory holds a single RunInfo file. `run_file(rundir,
name, box)` does the same for the run-level files (`timing.log`, `main.log`,
`ic2pm.log`, `BDM.config`, `Rockstar.config`, `OutputSchedule.dat`, ...), which
`read_text` reads in either form. So code written
against the native layout runs unchanged, and native-only runs never import
h5py:

```python
cat = glamio.read_catalogue(glamio.product_path("Run0", "catalogue", 184, box=0))   # Run0/Output.0184.0000.h5/CatshortV
rs = glamio.read_rockstar(glamio.product_path("Run0", "rockstar", 184, box=0))      # Run0/Output.0184.0000.h5/Rockstar
planes = glamio.read_lensplanes("Run0/Output.0184.0000.h5/LensPlanes", mmap=True)
sched = glamio.read_schedule(glamio.product_path("Run0", "schedule", 0))            # Run0/RunInfo.0000.h5/schedule
```

Without glamio:

```python
import h5py, numpy as np
with h5py.File("Run0/Output.0184.0000.h5") as f:
    mbound = f["CatshortV/data"]["Mbound"]                     # one column of the compound table
    mvir = f["Rockstar/data"]["mgrav"]                         # the Rockstar catalogue (status: f["Rockstar"].attrs["status"])
    pk = f["PowerDM/R/table"][:]
    z = f.attrs["z"]
    meta = f["deltam"].attrs["meta_text"].decode()             # the sidecar text
    log = f["CatshortV/outputB"][()].tobytes().decode()        # a text file (uint8 dataset)
    d = f["deltam/data"]
    cube = np.memmap("Run0/Output.0184.0000.h5", dtype=d.dtype, mode="r",
                     offset=d.id.get_offset(), shape=d.shape)  # memory map, no copy
```

`h5dump -d /deltam/data -b LE -o out.bin` writes nothing for compound
datasets (the catalogue, the complex cubes). To extract a dataset's bytes
with shell tools, take `OFFSET` and `SIZE` from `h5dump -H -p -d /deltam/data
Output.0184.0000.h5` and run `tail -c +$((OFFSET + 1)) Output.0184.0000.h5 |
head -c SIZE`; the result has the manifest's md5.

## Definitions and caveats (from the finder's code, 2026-10-18 review)

* `Rvir` is the spherical-overdensity radius `R_SO` (enclosed density =
  `overdensity` x mean matter density) plus the finder's resolution correction
  `min(Rext/(R_SO/cell)^SlopeR, 0.75)` cells; `Mtot` counts every particle
  inside `Rvir`, `Mbound` the particles inside `R_SO` that survive unbinding.
  The binary sidecar states this in its `definitions` line and records `Rext`,
  `SlopeR`, `cell`, `MG_flag`, `MG_model` and `a_x` as well.
* `Vcirc` and `Cvir` are 0 when the circular-velocity maximum is unresolved
  (0.11 % of the pilot's rows).
* Velocities (`vx, vy, vz`, and the derived `Vrms`, `virial_ratio`, `lambda`)
  are written with the legacy `1/a` conversion: they are LOW by `a_v/a`
  (0.8 % at z = 2, 0.3 % at z = 0 in the pilots). `velocities_at_av()` multiplies
  `vx, vy, vz` by `a/a_v`.
* The text reader returns the same float32 layout as the binary reader; the
  cast adds up to one float32 ulp (6e-5 Mpc/h for x > 256) to positions the
  text printed with 4 decimals. Use `np.loadtxt` directly if you need the
  printed decimals exactly. Text and binary `meta` share the key names `a`,
  `da`, `a_v`, `z`, `step`, `realization`, `rows`, `product`; the binary sidecar
  has many more (see above) and its `format` value is a description, while the
  text reader sets `format = text`.
* `read_catalogue(stem, mmap=True)` memory-maps the `.bin` instead of loading
  it (about 1 GB per output for 3200^3 runs). The binary reader refuses a file
  whose `row` column is not `1..N` (byte-swapped, misaligned or corrupted) and a
  sidecar without `columns`/`dtypes` lines.

## Examples

`examples/` holds two scripts that read every product of the node-2 runs and
make the figures of the status deck (`measure_products.py`, `plot_products.py`;
see `examples/README.md`). They are the reference for the reading patterns and
the physics checks (P(k) vs linear theory, cube vs table, momentum continuity,
plane growth, velocity rms, mass functions, chameleon screening).
