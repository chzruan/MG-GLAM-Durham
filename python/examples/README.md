# glamio examples: the on-the-fly products as data

Two scripts that read every product of a pair of GLAM runs through `glamio`
and turn them into six summary figures (the "In the data" frames of the
DEGRACE companion status deck). They double as worked examples of the reader
API and of the physics checks each product allows.

| script | does | runtime (login node, 4 threads) |
|---|---|---|
| `measure_products.py {fast,cat,planes} [--lcdm RUN1] [--fr RUN1] [--out DIR]` | reads the products and caches plot-ready arrays in `products_<stage>.npz` (JSON provenance under the key `provenance`) | `fast` 60 s, `cat` 11 s, `planes` 15 s; < 4 GB |
| `plot_products.py [--data DIR] [--out DIR] [fig ...]` | draws `products_{pk,modes,fields,planes,velocities,haloes}.pdf/.png` in the house style (`house_style.py`, `chz-paper.mplstyle` beside it) | 20 s |

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
