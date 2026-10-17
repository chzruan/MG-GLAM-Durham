# DEGRACE ECOSMOG suite: pinned metadata for MG-GLAM companion runs

Baojiu Li's DEGRACE simulations (ECOSMOG = RAMSES + f(R) solver) and their
ΛCDM twins. Surveyed 2026-10-16 from the files on disk; every number below
was read from the generator scripts, the parameter files or the run
directories, and cross-checked by `extract_metadata.py` (all 96 nodes pass
every check). The goal is to run MG-GLAM companions with identical initial
conditions, cosmologies and output epochs.

Files here:

| file | content |
|---|---|
| `nodes.csv` | one row per node: design parameters, derived cosmology, the values actually found in `2LPTic.param`, the CAMB inputs actually run, ECOSMOG f(R) parameters, run status, output-list checks |
| `aout.json` | the 42 output expansion factors of every node (from the ECOSMOG namelists) |
| `extract_metadata.py` | regenerates both files and prints the consistency report (reads bl267's directories on cosma8 and cosma7) |

## 1. Where everything is

| item | path | notes |
|---|---|---|
| design table | `/cosma7/data/dp004/bl267/Runs/DEGRACE/ICs/CAMB/Parameter_files/Python_codes/DEGRACE_L512.csv` | 512 Sobol rows × 12 columns; rows 1–96 are the nodes |
| generators | same directory: `Run_CAMB.py`, `Run_2LPTic.py`, `Gen_ECOSMOG_Params.py`, `Gen_POWMES_Params.py`, `Gen_ECOSMOG_Params_single.py` (Node_900) | the pipeline design → CAMB → 2LPTic → ECOSMOG |
| CAMB inputs | same directory: `params.ini` (template), `params_n001..096.ini` | the `.ini` actually run per node |
| CAMB linear P(k) | `ICs/IC_data/L1024/Node_XXX/PK_CAMB_z00_XXX.dat` (cosma8 copy; original `test_PK_CAMB_z00_XXX.dat` on cosma7, byte-identical) | `k/h`, `P(k) [(Mpc/h)³]`, 617 rows, z = 0, CAMB `do_nonlinear = 0` |
| initial conditions | `/cosma8/data/dp203/bl267/Projects/Ongoing/DEGRACE/ICs/IC_data/L1024/Node_XXX/ics.0 … ics.27` | copy of cosma7 `ICs/IC_data/L1024` (generated as `L1024_new`); 96 nodes |
| ΛCDM runs | `/cosma8/data/dp203/bl267/Projects/Ongoing/DEGRACE/Runs/LCDM/Node_XXX/` | `namelist/cosmo.nml`, `data/output_000NN/info_000NN.txt`, `gadget_data/` |
| f(R) runs | `/cosma7/data/dp004/bl267/Runs/DEGRACE/Runs/L1024/FR/Node_XXX/` | same layout; `nrestart` set (checkpoint restarts) |
| f(R) power spectra | `/cosma7/data/dp004/bl267/Runs/DEGRACE/Runs/L1024/FR/power_spectra/Node_XXX/PK_FR_NXXX_SNN.dat` | POWMES, 43 snapshots per node |
| IC power spectra | `ICs/IC_data/L1024/Node_XXX/PK_LCDM_NXXX_S01.dat`, `inputspec_ics.txt`, `log_2LPTic.txt`, `log_POWMES.txt` | POWMES on the ICs (z = 49); 2LPTic echo of the input spectrum |
| fiducial node | `/cosma7/data/dp004/bl267/Runs/DEGRACE/Runs/L1024/FR/Node_900/` | Ω_m = 0.3089, h = 0.6774, Ω_b = 0.0486: only `run_cosma7.bash` (8 Oct 2026); no ICs, not run |
| tests | `Runs/L512/Node_002` (ΛCDM `data/` and f(R) `data_fR/`, 43 outputs each), `Runs/L256/LCDM/Node_001`, `ICs/IC_data/L512`, `L256` | small boxes, node 2 cosmology |
| not run | `/cosma7/.../Runs/L1024/GR/Np2048/Node_002` | 2048³ ΛCDM namelist pointing at non-existent `IC_data/L1024h` |
| obsolete | `ICs/IC_data/L1024_old/` (101 nodes) | an earlier design with different cosmologies per node (e.g. Node_001: Ω_m = 0.3328, h = 0.7306). Do not use |

Nodes 097–100 exist only as `Runs/LCDM/Node_097..100` with a namelist and no
ICs or outputs (the generators loop over 1–96). The suite is 96 nodes.

## 2. The design

`DEGRACE_L512.csv` is a Sobol sequence: columns 7–12 are the unit-cube
coordinates, columns 1–6 the mapped parameters. Row 1 is the lower corner of
every range, row 2 the centre.

| column | parameter | range (rows 1–96) |
|---|---|---|
| 1 | Ω_m | 0.25 – 0.3492 |
| 2 | h | 0.65 – 0.7492 |
| 3 | S_8 | 0.65 – 0.9453 |
| 4 | n_s | 0.90 – 1.0488 |
| 5 | log10 \|f_R0\| | −6.5 – −4.02 |
| 6 | unused by the scripts | −0.45 – 0.98 (probably a second model's parameter) |

Derived per node (all verified against the files on disk):

- Ω_b = 0.02237 / h², Ω_Λ = 1 − Ω_m (flat), Ω_c = Ω_m − Ω_b in the scripts' intent
- σ_8 = S_8 √(0.3 / Ω_m); the ICs carry this to within 1e-4 (the CAMB A_s tuning tolerance), column `sigma8_ic`
- f_R0 = −10^(col 5), Hu–Sawicki n = 1; ECOSMOG `fr_c1oc22 = |f_R0| (12/Ω_m − 9)²` (= `|f_R0| [3(1 + 4(1−Ω_m)/Ω_m)]²`)
- ECOSMOG `dens_V_on` = 80 for |f_R0| < 1e-6, 120 for < 1e-5, 160 otherwise (19 / 39 / 38 nodes)

## 3. Linear power spectra (CAMB) and a pitfall

Template values (`params.ini`, Planck 2018 base): ω_b = 0.02237, ω_ν = 0.0006442
(one massive neutrino, 0.06 eV, N_eff 3.046, `share_delta_neff`), Ω_k = 0,
w = −1, T_CMB 2.7255, Y_He 0.24, τ = 0.0544, pivot 0.05 Mpc⁻¹, no running,
`transfer_kmax = 20`, `transfer_k_per_logint = 5`, high-precision transfer, one
output redshift z = 0. Per node the script sets `hubble`, `omch2`,
`scalar_spectral_index(1)` and iterates `scalar_amp(1)` until CAMB's σ_8(z=0)
equals the design σ_8 within 1e-4 (column `A_s`).

**Pitfall.** `Run_CAMB.py` computes `Oc = Om − Ob` but then writes
`omch2 = Om * h0**2` (the total matter density) into the CDM slot. The CAMB
cosmology therefore has Ω_m^CAMB = Ω_m + Ω_b + Ω_ν, 0.041–0.055 above the Ω_m
used by 2LPTic and ECOSMOG (column `Omega_m_camb`). The *shape* of every
linear spectrum corresponds to this higher ω_m, while the N-body background,
the growth used by 2LPTic to scale to z = 49, and the σ_8 normalisation
target use the design Ω_m. For companion runs this does not matter as long as
MG-GLAM ingests the same ICs, but it does matter for anything that labels the
simulations with a linear theory (emulator training, P(k) ratios to "linear",
halo-model inputs). Raise it with Baojiu before building anything on the
`PK_CAMB` files; `nodes.csv` records both Ω_m values.

## 4. Initial conditions

2LPTic (N-GenIC lineage, MPI, 28 ranks), identical settings for every node:

| setting | value |
|---|---|
| Box | 1024 Mpc/h |
| particles | 1024³ (`Nsample = Nmesh = 1024`, `GlassTileFac = 1024` of a one-particle glass, i.e. a regular grid) |
| start | z = 49 (a = 0.02), 2LPT |
| seed | **2026 for all nodes** (fixed phases across the design) |
| spectrum | `WhichSpectrum = 2`, file `PK_2LPTic_z00_XXX.dat` = (log10 k, log10 4π k³ P) made from the CAMB file; `PrimordialIndex = 1` (no extra tilt); renormalised to `Sigma8` (log line `Normalization adjusted to Sigma8=… (Normfac=…)`, column `normfac_2lptic`) |
| units | Mpc/h, 1e10 M_sun/h, km/s |
| output | Gadget-1, little-endian, 28 files, 1.09 GB each (30 GB per node); header Box = 1024, Ω_m/Ω_Λ/h of the node, mass 6.94 (node 1) |
| growth | `inputspec_ics.txt` first line = (49, D_+ ratio), column `Dplus_z49` |

The same files feed the ΛCDM and f(R) runs of a node, so the pair is a
matched-phase twin.

## 5. ECOSMOG configuration

| namelist item | value |
|---|---|
| `levelmin / levelmax` | 10 / 18 (1024³ root grid, 8 refinement levels, finest cell 3.9 kpc/h) |
| `m_refine` | 4,4,4,8,… ; `nsubcycle = 1,2,2,…`; `nremap = 10`; `ngridtot = 3e9`, `nparttot = 1.5e9` |
| ΛCDM twin | identical namelist with `extradof = .false.` (f(R) parameters still present) |
| job | 560 MPI ranks, `cosma8` partition (ΛCDM, account dp203) / cosma7 (f(R), dp004), 72 h, restarts |
| cosmology | read from the Gadget IC header (RAMSES `info` files agree with `2LPTic.param` for all nodes) |

Status: 96/96 ΛCDM and 96/96 f(R) runs have all 43 outputs and end at a = 1.

## 6. Output epochs

Each run writes 43 outputs: `output_00001` at a = 0.02 (the ICs) plus the 42
`aout` values. The generator builds them exactly as the Slack list for the
fiducial cosmology: 16 fixed redshifts 0.125 … 2.0, a = 1, and 25 stations
where the comoving distance χ(a) = 256 k Mpc/h, k = 1 … 25, computed with
`c/H0 = 2997.92458 Mpc/h` and E²(a) = Ω_m a⁻³ + 1 − Ω_m (matter + Λ only, the
node's own Ω_m). Recomputing the stations reproduces every list to 2.3e-7 in
a. The lists are in `aout.json`; they differ per node, and every node has
between 1 and 4 neighbouring targets closer than 0.5 % in a (handled by the
exact-z scheduler since `cz` commit `ae4fa9d`). ECOSMOG dumps at the first
step past each target, so the actual snapshot epochs carry small offsets (the
Node_082 headers in the earlier memory note show this).

## 7. Snapshot products

- The raw RAMSES particle files were deleted after conversion; each
  `data/output_000NN/` keeps only `info_000NN.txt`.
- Converted snapshots: `gadget_data/output_000NN/gadget.0 … gadget.559` (one
  file per rank, r2g conversion; `r2g.config.snap_NN` and `log_r2g_snap_NN.txt`
  beside them). Little-endian Gadget-1, positions in **kpc/h** (header
  BoxSize = 1023999.98), masses 6.94 × 1e10 M_sun/h, IDs 1 … 1024³, 26–29 GB
  per snapshot, about 1.2 TB per node and gravity model. The velocity
  convention of r2g has not been verified here.
- Nodes 001–080 also carry `gadget_data_old` (conversion of 15 Mar 2026,
  byte-different from the 27 Mar `gadget_data`); nodes 081–096 were converted
  once, later. Use `gadget_data` and ask Baojiu what changed in r2g.
- f(R): POWMES P(k) for all 43 snapshots per node (`PK_FR_NXXX_SNN.dat`,
  1024³ grid, folding; `plot.py` in the IC directory shows the column use:
  k = col1 × 2π/L, P = col4 × L³). No ΛCDM snapshot P(k) was found.

## 8. What an MG-GLAM companion run needs per node

| MG-GLAM input | source |
|---|---|
| ICs | `ICs/IC_data/L1024/Node_XXX/ics.*` through `ic2pm` (default `half` velocity epoch, see `AGENTS.md`); Box 1024, 1024³ particles, z_init 49 |
| cosmology (`Init.dat`, `PkTable.dat` header) | `nodes.csv`: Ω_m, h, Ω_b = 0.02237/h², σ_8 = `sigma8_ic`, n_s |
| linear P(k) table | `PK_CAMB_z00_XXX.dat` normalised to `sigma8_ic`, with the Ω_m^CAMB caveat of section 3 |
| f(R) | `MG_model` Hu–Sawicki n = 1, `fR0` from `nodes.csv` (negative); ΛCDM twin with `MG_flag = 0` |
| outputs | `#outputs = -42` with z = 1/a − 1 from `aout.json` (exact-z mode) |
| resolution | ECOSMOG resolves to 3.9 kpc/h; GLAM at Ng = 2048 (0.5 Mpc/h cells) or 4096 (0.25 Mpc/h) with the same 1024³ particles. A 2048³ GLAM run would need new ICs; 2LPTic with the same seed keeps the large-scale modes only if its seed-table scheme is resolution-independent, which has not been checked |
| comparison data | f(R) P(k) at all 43 epochs; Gadget snapshots of both gravities for halo catalogues and lens planes |

## 9. Open questions for Baojiu

1. Was `omch2 = Ω_m h²` in the CAMB inputs intended (section 3)?
2. What changed between the 15 Mar and 27 Mar r2g conversions, and is `gadget_data_old` disposable?
3. r2g velocity convention (Gadget √a factor or peculiar km/s)?
4. Is Node_900 (fiducial) going to be run as f(R) or ΛCDM, and with which f_R0? Its output list is the Slack list.
5. Are the Node_097–100 directories and `L1024_old` disposable?
6. Will ΛCDM snapshot P(k) be measured, or should the GLAM side measure both?
