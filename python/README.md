# glamio

Numpy-only Python readers for GLAM / MG-GLAM simulation outputs. The first
module reads the BDM halo catalogues written to `CATALOGS/`, in both the legacy
text form and the new binary form.

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

## Planned modules (not implemented yet)

* `schedule`: `OutputSchedule.dat`
* `lensplanes`: lens-plane outputs
* `modes`: `deltam` / `mom_*` / `deltaL` cubes
* `subsample`: particle subsamples
* `frhalo`: f(R) halo scalar-field outputs
* `powerdm`: matter power spectra
* `pmcrd`: big-endian PM snapshot headers

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
