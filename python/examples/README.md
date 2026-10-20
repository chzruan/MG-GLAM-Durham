# glamio examples: the on-the-fly products as data

Two scripts that read every product of a pair of GLAM runs through `glamio`
and turn them into six summary figures (the "In the data" frames of the
DEGRACE companion status deck). They double as worked examples of the reader
API and of the physics checks each product allows.

| script | does | runtime (login node, 4 threads) |
|---|---|---|
| `measure_products.py {fast,cat,planes} [--lcdm RUN1] [--fr RUN1] [--out DIR]` | reads the products and caches plot-ready arrays in `products_<stage>.npz` (JSON provenance under the key `provenance`) | `fast` 60 s, `cat` 11 s, `planes` 15 s; < 4 GB |
| `plot_products.py [--data DIR] [--out DIR] [fig ...]` | draws `products_{pk,modes,fields,planes,velocities,haloes}.pdf/.png` in the house style (`house_style.py`, `chz-paper.mplstyle` beside it) | 20 s |
| `halo_environment.py RUN1 [--z Z | --step N] [--R R ...] [--filter gaussian|tophat] [--ngrid N] [--interp cubic|linear] [--out FILE.npz]` | smoothed overdensity δ_R, tidal eigenvalues λ₁ ≥ λ₂ ≥ λ₃, anisotropy α_R and T-web class at every halo of one output, from the `deltam` cube and the catalogue | 40 s and 5 GB per R at 512³; 100 s and 15 GB at 768³ |

The defaults point at the node-2 runs under `DEGRACE-companion/runs/` (the ΛCDM
confirmation run, job 12185817, and the f(R) pilot, job 12180110); any two
`Run1` directories with the same output schedule work (the f(R) one needs
`fR_halo` files for the `cat` stage). Run from the repository root with the
cosemu environment:

```sh
micromamba run -n cosemu python3 python/examples/measure_products.py fast   --out /tmp/figs
micromamba run -n cosemu python3 python/examples/measure_products.py cat    --out /tmp/figs
micromamba run -n cosemu python3 python/examples/measure_products.py planes --out /tmp/figs
MPLCONFIGDIR=/tmp/matplotlib-config micromamba run -n cosemu python3 python/examples/plot_products.py --data /tmp/figs
```

What each figure checks (node 2, L = 1024 Mpc/h, 1024^3 particles, N_g = 2048):

* `pk`: `PowerDM.R` at the 17 primary outputs against the linear input spectrum
  (`PkTable.dat` times D^2(a); ±2 % for 0.02 < k < 0.1 h/Mpc) and the f(R)/ΛCDM
  ratio of the two runs (+23 % at k = 5 at z = 0, +4 % at z = 2).
* `modes`: the `deltam` cube's own `power()` equals the table shell by shell
  (4e-6); θ_k = i k·p_k from the three `mom_*` cubes gives P_θθ = (aHf)^2 P_δδ
  to 2 % for k < 0.06 with no free parameter, and the density-weighted term of
  p = (1+δ)v dominates beyond k ≈ 0.1.
* `fields`: a lens plane, the inverse transform of the `deltam` cube over the
  same 256 Mpc/h slab (Parseval-checked) and the particle subsample coloured by
  v_z, at z = 0.
* `planes`: rms of Σ/Σ̄ on 8 Mpc/h pixels for all 12 planes at all 42 outputs
  against D(a) (follows it to z ≈ 1, +7 % at z = 0; the native 0.5 Mpc/h pixels
  carry 0.125 of particle discreteness in quadrature) and the pixel PDF.
* `velocities`: the subsample's 1-D velocity rms at the 17 primaries against
  a H f D [∫P_lin dk / 6π^2]^(1/2) (+4/+8/+15 % at z = 2/1/0) and the velocity
  PDF.
* `haloes`: cumulative mass functions of both runs with their ratio, and
  |f_R(centre)| / |f_R,bg(a)| against mass from `fR_halo` (screening
  transition near 3e13 M_sun/h at z = 0, everything above 2e12 screened at z = 2).

The shell-averaging in `measure_products.shell_power` is the same integer-shell
binning as GLAM's `PowerDM` (modes with `round(|k|/k_f) = s`), with Hermitian
pairs weighted twice, so the cube spectra can be compared with the table bin by
bin. The momentum cubes sit on the N_dep = 1024 mesh; the plot removes each
field's own CIC window with Jing's shell-averaged form before taking ratios.

## Halo environment for assembly-bias work: `halo_environment.py`

```sh
nice -n 19 micromamba run -n cosemu python3 python/examples/halo_environment.py \
    DEGRACE-companion/runs/Node_002_LCDM_L1024Np1024Ng2048_confirm20261018/Run1 \
    --z 0 --R 2 5 --out /tmp/env_z0.npz
```

Recipe: the `deltam` cube (all modes to k_c = 1.5 h/Mpc) is placed on the
rfft layout of an N³ grid, divided by the CIC window, multiplied by the
smoothing kernel W(kR) (Gaussian by default, `--filter tophat` optional) and
inverse-transformed; T_ij = (k_i k_j / k²) δ W gives six more inverse FFTs.
Each field is read at the halo positions by periodic cubic-spline
interpolation (the B-spline prefilter is applied in Fourier space, where it
is an exact division for periodic data; scipy's recursive prefilter costs
60 s per field). Per halo the output holds `delta_R<R>`, `lambda_R<R>`
(n_halo × 3, descending, summing to δ_R), `alpha_R<R>` = √q²/(1+δ_R) with
q² = ½ Σ(λ_i − λ_j)² (Paranjape, Hahn & Sheth 2018), `nweb_R<R>` (number of
eigenvalues above `--lambda-th`, default 0: void/sheet/filament/knot), the
catalogue columns `row x y z Mtot Mbound Rvir Vcirc Vrms Cvir Npart Xoff
virial_ratio lambda RadRMS b_a c_a`, and JSON `provenance` with the grid-level
checks.

Validation (node-2 ΛCDM confirmation run, z = 0, 777,663 haloes):

* a direct Fourier sum over all 58 million modes at six halo positions
  (including the most massive one, δ_2 = 222) agrees with the 512³ cubic
  values to 2 % at R = 2 and 0.02 % at R = 5; a 768³ grid brings R = 2 to
  0.25 %. Over all haloes the 512³ vs 768³ difference at R = 2 is 1 % in the
  median and 2.8 % at the 99th percentile (δ_R, λ₁ and α_R alike, T-web class
  identical for 99.2 %); at R = 5 it is 0.02 %. The script warns when the
  cells are coarser than 0.67 R. Linear interpolation (`--interp linear`) is
  10–40 % low at R = 2 and 1–2 % low at R = 5: haloes sit at peaks.
* the sum of the eigenvalues equals δ_R to 2e-5 (float32), the grid variance
  equals the Parseval sum of the filtered cube, σ(δ_2) = 1.56 and σ(δ_5) = 0.67.
* mass trend at z = 0: ⟨δ_5⟩ = 0.29, 0.57, 1.0, 1.9, 3.7, 7.7 and ⟨α_5⟩ =
  0.55, 0.50, 0.42, 0.32, 0.23, 0.15 for half-dex bins from 10^12.2 to
  10^15.2 Msun/h: massive haloes sit in denser, more isotropic (knot)
  surroundings, as in the literature.
* the f(R) pilot (text catalogue, 856,309 haloes) runs through the same
  path: σ(δ_2) = 1.62, σ(δ_5) = 0.69.

Limits: the cube cannot resolve halo-scale tides (Gaussian R below about
1.5 Mpc/h loses more than 10 % of the kernel to k > k_c); the density-built
T_ij is the geometric tidal field, identical in definition for both gravity
models (the f(R) fifth force is not in it).

