# glamio

Numpy-only Python readers for GLAM / MG-GLAM simulation outputs: the BDM halo
catalogues (legacy text and binary), the output schedule, the lens planes, the
low-k mode cubes, the quantised particle subsample, the f(R) halo values, the
power spectra and the PM particle snapshots. Every reader validates the file
against its sidecar (product name, required keys, byte count) and raises
`ValueError` instead of misreading. The precise_z suite cross-checks each reader
against independent decoders on real GLAM output
(`test/precise_z/asserts.py`, the `glamio` checks).

## Install

From the repository root, into the `cosemu` environment:

```
micromamba run -n cosemu python3 -m pip install -e python/
```

Only `numpy>=1.24` and Python >= 3.10 are required.

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

## Other products

All per-output products share the file scheme `<product>.<step>.<box>.{bin,meta}`
in the run directory (`Run1/`), raw little-endian binary without record
markers plus a `key = value` sidecar; each reader accepts the `.bin`, the
`.meta` or the bare stem, and `mmap=True` memory-maps the binary. The returned
objects expose the sidecar as `.meta` and the common scalars as properties
(`a`, `z`, `a_v`, `box`, `step`, `realization`, `n_grid`, `n_row`, `n_p`,
`omega_m`, `omega_l`, `hubble`, `header`; `None` when the key is absent).

| reader | file | returns |
|---|---|---|
| `read_schedule(path)` | `OutputSchedule.dat` | `Schedule`: structured `rows` (`step a z chi exact class`), `omega_m/omega_l` from the comment, `row(step)`, `primary()`, `stations()` |
| `read_lensplanes(stem)` | `LensPlanes.*` | `LensPlanes`: `planes[plane, j, i]` float32 (`3 N_slab` planes X1..XN, Y1..YN, Z1..ZN), `plane(axis, slab)` as `[i, j]`, `slabs(axis)`, `axis_means()`, `pixel_edges()`, `slab_edges()` |
| `read_modes(stem)` | `deltam.*`, `mom_x/y/z.*`, `*_deltaL.*` | `Modes`: `data[ix+n, iy+n, iz]` complex64, `mode(ix, iy, iz)` with Hermitian completion for `iz < 0`, `wavenumbers()`, `k_magnitude()`, `power()` (shell-averaged `Box^3 <|delta_k|^2>`) |
| `read_subsample(stem)` | `Subsample.*` | `Subsample`: `q[6, n_sel]` int16, `positions()` Mpc/h, `velocities()` km/s at `a_v`, `ids()` (0-based IDs from the splitmix64 rule), `lattice` tag |
| `read_frhalo(stem)` | `fR_halo.*` | `FRHalo`: `data` (`fR_centre`, `fR_mean`) float32 per catalogue row, `catalogue()` reads the catalogue of the same step |
| `read_powerdm(path)` | `PowerDM.R/log/Z.*.dat` | `PowerDM`: `k`, `nmodes`, `pk` (raw CIC), `a`, `z`, `step`, `ngrid`, `k_nyq`, `binning`; `pk_deconvolved()` with the approximate shell-averaged CIC window |
| `read_pm_header(path)` | `PMcrd.<step>.DAT` | `PMHeader`: the 529-byte big-endian header (`aexpn`, `astep`, `istep`, `nrow`, `ngrid`, `nparticles`, `extras`, `box`, `a_v`, `velocity_scale()`) |
| `read_snapshot(path)` | `PMcrd` + `PMcrs<i>` | `Snapshot`: `iter_pages()` streams 1024^2-particle pages, `read()` loads all (24 B per particle), `units='code'` (grid units, code momenta) or `'physical'` (Mpc/h in [0, Box), km/s at `a_v`) |

```python
import glamio

sched = glamio.read_schedule("Run1/OutputSchedule.dat")
for step in sched.primary()["step"]:
    planes = glamio.read_lensplanes(f"Run1/LensPlanes.{step:04d}.0001", mmap=True)
    sigma_x1 = planes.plane("x", 1)                   # [iy, iz], Sigma/Sigma_bar of the first X slab
    modes = glamio.read_modes(f"Run1/deltam.{step:04d}.0001")
    k, pk, nmodes = modes.power()                     # equals PowerDM.R bin by bin (raw CIC)
    sub = glamio.read_subsample(f"Run1/Subsample.{step:04d}.0001")
    pos, vel, ids = sub.positions(), sub.velocities(), sub.ids()

snap = glamio.read_snapshot("Run1/PMcrd.0181.DAT")    # header only; particles stream on demand
for x, y, z, vx, vy, vz in snap.iter_pages(units="physical"):
    ...                                               # one 1024^2-particle page at a time
```

Conventions the readers encode (from the Fortran writers):

* Lens planes: `(N_lens, N_lens)` Fortran column-major per plane, X planes
  indexed `(y, z)`, Y planes `(x, z)`, Z planes `(x, y)`; pixel `i` covers
  `[i, i+1) Box/N_lens`, slab `s` (1-based) covers `[(s-1), s) Box/N_slab`.
  `planes[p]` is the raw `[j, i]` layout; `plane(axis, s)` transposes to `[i, j]`.
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
  conversion). The last page of a `PMcrs` file is padded to 1024^2 particles;
  the reader drops the padding.
* `PowerDM`: raw CIC density on the `N_grid` mesh, window not deconvolved, shot
  noise not subtracted; `pk_deconvolved()` divides by
  `1 - (2/3) sin^2(pi k/2 k_Nyq)`, an approximation of the shell-averaged window.

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
