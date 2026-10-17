# MG-GLAM implementation plan for the DEGRACE companion runs

One plan for everything MG-GLAM needs before the DEGRACE companion campaign:
the output schedule and lens planes agreed with Baojiu Li on Slack, and the
extra on-the-fly outputs for single-index halo statistics (initial-field
modes, late-time low-k modes, particle subsample, scalar field at haloes).

Merged 2026-10-17 from `docs/degrace_lensplanes_implementation_plan.md`
(Slack-derived plan, 2026-10-16, now this file) and
`/cosma8/data/dp203/dc-ruan1/glam-on-the-fly-output/opus_glam_extra_outputs_plan.md`
(extra-outputs brief). This file supersedes both. Read it whole before
touching the code; answer the Phase 0 questions of section 13 from the code
before any edit and record the answers in section 16.

Status (2026-10-17): the cross-check of the 42 output epochs is done
(section 2); **I1 is on `cz`** (`ae4fa9d`); **Phase 0 is answered**
(section 16); **I2, I3, I4, I5, I8, I9, I10 are implemented on branch
`degrace-track-a`** (175/175 checks plus two unit-test scripts); the
node-2 ΛCDM pilot is running (Slurm job 12180082, section 16). I6 is
encoded in the input generator and launcher; I7 and I11 are open.

Related references: `docs/degrace_ecosmog/README.md` (Baojiu's ECOSMOG
suite: design, ICs, output lists, snapshots), `2LPTIC_Gui/VALIDATION.md`
(2LPTIC and time-stepping validation), review decks under `docs/slides/`.

## 0. Ground rules

1. **Initial conditions are always 2LPTIC** (the FML generator in
   `2LPTIC_Gui/`, or Baojiu's original 2LPTic `ics.*` for matched-phase
   companions), ingested with `ic2pm` (default `half` velocity epoch).
   GLAM's own `PMP2start` Zel'dovich generator is used only to reproduce
   legacy runs. Start at z_init = 49.
2. **Time stepping**: da0 = 6e-4 for force resolution Δx ≤ 0.25 Mpc/h
   (163 steps, the legacy late-time ladder), 8e-4 acceptable at Δx = 0.5,
   4e-4 for a converged reference (see `2LPTIC_Gui/VALIDATION.md` and the
   memory note on stepping). Pipe a step count ≥ the ladder length.
3. **Do not change physics or numerics.** PM solver, multigrid, stepping,
   BDM and the existing output formats stay as they are.
4. **Every addition is behind a runtime flag that defaults to off.** With
   all flags off the code must be bit-identical to the current `cz`.
5. Branch from `cz` (the project mainline, not `main`). New code goes in
   new source files where possible, wired in through minimal call sites.
6. **No new N_g³ allocations** unless the memory budget (section 13, Q8)
   allows it. f(R) at N_g = 4096 is memory-tight (FI, FI2, FI3 plus the
   multigrid scratch). Reuse `FI` at the output step (section 4) or use
   coarser deposit grids.
7. Compute: test runs through Slurm on COSMA8 compute nodes
   (`cosma8-serial`/`cosma8-shm`), only light jobs on login nodes. Ask
   before any job above ~50 node-hours and before any production-size run.
8. Test configurations: the 36³/72³ box in `test/precise_z/` for schedule
   logic; a 256³ particle, L = 256 Mpc/h, N_g = 512 box with a few outputs
   for the output products, one ΛCDM and one f(R) (F5) with the same seed.
9. After each phase (section 15) append what was done, file names and
   test numbers to section 16.

## 1. Agreed output design (Slack, Baojiu, Oct 2026)

Scope: the **L1024 Np2048 Ng4096** production (100 ΛCDM boxes in the fiducial
Planck 2015 cosmology, later the design cosmologies) and the companions of
Baojiu's 96-node ECOSMOG suite (same box, 1024³ particles, his ICs).

| Item | Decision |
|---|---|
| Particle snapshots | **None** (`iSave=0`). 200 GB per snapshot at 2048³. Checkpoints are still written for restart and deleted at the end. |
| Velocities | Kept in the BDM halo catalogues, and in the 1% quantised subsample (I9). |
| Lens planes | At every output: project the density along each of the 3 axes in **4 slabs of 256 Mpc/h** → **12 planes of 4096² float32** = 0.8 GB per output. |
| Outputs per box | **42** = 17 *primary* redshifts common to all cosmologies (z = 0, 0.125, …, 2.0) + 25 *secondary* at comoving distance χ = 256 k Mpc/h, k = 1…25. FORGE / DEGRACE-ECOSMOG convention. |
| Secondary set | Cosmology-dependent: recompute a(χ = 256 k) per cosmology (section 5, I3). |
| Plane resolution | 4096² (= N_g); keep a block-average option. |

## 2. Verified numbers

**Cross-check (2026-10-16).** Baojiu's 42 scale factors reproduce flat ΛCDM
with Ω_m = 0.3089, Ω_Λ = 0.6911, **no radiation**, c/H₀ = 2997.92458 Mpc/h:

| comparison | max \|Δa/a\| over the 25 stations |
|---|---|
| Ω_m = 0.3089, matter + Λ only | 7.7e-8 |
| same + radiation (h = 0.6774, N_eff = 3.046) | 1.7e-3 (2.8 Mpc/h at χ = 6400) |
| Ω_m = 0.3156 | 2.9e-2 |

The same convention generates every ECOSMOG node's list (verified to 2.3e-7
for all 96 nodes, `docs/degrace_ecosmog/aout.json`); the lightcone code must
use the same matter + Λ χ(a).

**The 42 fiducial targets, ascending in a** (`P` primary, `S k` station k):

```
0.09885636873988636 S25   0.11436339923360468 S24   0.13101192462219735 S23
0.14880957763634076 S22   0.16776699787816512 S21   0.18789864408246054 S20
0.2092238939643731  S19   0.2317683231209623  S18   0.25556525683790143 S17
0.2806576241072646  S16   0.30710027529154543 S15   0.3333333333333333  P z=2.000
0.33496251467163335 S14   0.34782608695652173 P z=1.875   0.36363636363636365 P z=1.750
0.3643315849360659  S13   0.38095238095238093 P z=1.625   0.3953165138959508  S12
0.4                 P z=1.500   0.42105263157894735 P z=1.375   0.42805303395318656 S11
0.4444444444444444  P z=1.250   0.4627096181121927  S10   0.47058823529411764 P z=1.125
0.49949502290382997 S9    0.5                 P z=1.000   0.5333333333333333  P z=0.875
0.5386678234812446  S8    0.5714285714285714  P z=0.750   0.5805485755215236  S7
0.6153846153846154  P z=0.625   0.6255358719332695  S6    0.6666666666666666  P z=0.500
0.6741267910946328  S5    0.7269451849577723  S4    0.7272727272727273  P z=0.375
0.7847794661808203  S3    0.8                 P z=0.250   0.8486352893278191  S2
0.8888888888888888  P z=0.125   0.9198100312830662  S1    1.0                 P z=0 (= S0)
```

Four pairs are closer than the 0.5 % merge tolerance the exact-z code had
before `ae4fa9d` (z = 2.0/1.98541, 1.75/1.74475, 1.00202/1.0, 0.37562/0.375;
the closest 0.045 %, Δa = 3.3e-4). Every ECOSMOG node has 1–4 such pairs.
Since I1 they are accepted: the later member becomes a short inserted step,
which the variable-step leapfrog integrates without loss of accuracy
(velocity error ~2e-6 per short step, below float32 noise).

Ladder for z_init = 49, da0 = 4e-4: 244 steps, late step 6.8e-3 (0.68 % of a
at a = 1); da0 = 6e-4: 163 steps, late step 1.025e-2 (identical to the legacy
z = 100 ladder at z < 2).

Reproduce the station cross-check (cosemu env):

```python
from scipy.integrate import quad; from scipy.optimize import brentq
cH0, Om = 2997.92458, 0.3089
chi = lambda a: cH0*quad(lambda x: 1/(x*x*((Om/x**3 + 1-Om)**0.5)), a, 1, epsrel=1e-12)[0]
a_station = [brentq(lambda x: chi(x) - 256.0*k, 0.01, 0.9999999) for k in range(1, 26)]
```

## 3. Shared conventions for every new product

**Parameters.** New flags are appended to the **end** of `Init.dat` and
mirrored in `Setup.dat` by `PMP2init` (writer ≈ `PMP2init.f90:330`, parser ≈
line 444, template ≈ line 540) and read by `ReadSetup` in
`PMP2mod_tools.f90` with `iostat`, defaulting to off, so every existing run
directory stays valid. Section 7 lists the flags.

**Fourier convention** (modes products I7, I8).
δ̃(k) = (1/N³) Σ_x δ(x) e^{−i k·x}, k = (2π/L)(ix, iy, iz), P(k) = L³⟨|δ̃|²⟩.
Convert from FFTPACK5's native layout and normalisation in one copy routine
used by every dump, because the dumps are cross-correlated mode by mode.

**Mode selection.** Store the cube |ix|, |iy| ≤ n, 0 ≤ iz ≤ n with
n = floor(k_c L / 2π); keep the Hermitian-redundant iz = 0 modes; set k = 0
to zero. (2n+1)²(n+1) complex(4) values.

**Array layout.** Fortran `C(0:n, -n:n, -n:n)` (first index iz, then iy, ix);
numpy `np.fromfile(f, dtype='<c8').reshape(2n+1, 2n+1, n+1)` → `[ix+n, iy+n, iz]`.
Index ix must be the physical x axis of the halo catalogue; record `axis_map`.

**Mass assignment.** Never deconvolve inside GLAM; record `assignment` (CIC),
`N_dep` and `interlaced`; offline code deconvolves.

**File format, all products (lens planes included).** Data: raw
little-endian binary, Fortran `access='stream', form='unformatted'`, explicit
`convert='little_endian'` on the OPEN (the build is globally big-endian),
no record markers. Metadata: a plain-text sidecar `.meta` with `key = value`
lines. Required keys: `git_hash`, `model`, `fR0` (or `none`), `L`, `N_p`,
`N_g`, `N_dep`, `a`, `z`, `a_x`, `a_v`, `step`, `box`, `units`, `axis_map`,
`D_code(a)` (the code's own ΛCDM linear growth at this a, normalised to
D(1) = 1), plus the item-specific keys below. Names:
`extra/<run>_a<a:.4f>_<product>.{bin,meta}`, `extra/<run>_deltaL.{bin,meta}`,
`extra/<run>_a<a:.4f>_lensplanes.{bin,meta}` (adapt to the existing
`<product>.<step>.<box>.dat` scheme if that is preferred; decide once).

**Velocity epoch.** GLAM's leapfrog leaves velocities half a step behind the
positions at every output (v at a − da/2 with the step just taken; in
exact-z mode that step can be as short as 3.3e-4). Every product that
carries velocities (I8 momentum, I9 subsample, the BDM catalogue) uses the
same (x, v) pair and records `a_x` and `a_v`.

## 4. Sequence at an output step

At a marked step the main loop (`PMP2main.f90` ≈ line 78) calls
`WriteDataPM` (if `iSave`) and `Analysis(Path, mDENSIT)` in
`PMP2mod_analyze.f90`. Today `Analysis` does: PDF (`iDensityDistr`) → power
spectrum (`iPower`: `DENSIT` if needed, `GetPower` FFTs `FI` in place and sets
`mDENSIT = 1`) → dump-random-fraction (`ProjectDENSIT`) → bias models →
`BDM(mDENSIT)` (re-deposits if needed). Proposed order with the new products:

1. `DENSIT` → `FI` holds the CIC overdensity δ on N_g³ at the output positions.
2. **Lens planes (I2)** from `FI` (sum over slabs; no FFT needed).
3. `GetPower` forward FFT of `FI` → **I8 density modes**: copy |k_i| ≤ k_c,δ
   right after the FFT (N_dep = N_g, no extra memory), then the usual P(k).
4. **I8 momentum**: reuse `FI` three times (deposit p_x, FFT, copy; p_y; p_z).
5. **I9 subsample** from the particle arrays (no grid).
6. `BDM` (re-deposits density into `FI` as it does now).
7. **I10 f_R at halo centres** while the converged scalar field of this step
   is still in memory, after BDM wrote the catalogue.

Everything between the P(k) FFT and `BDM` may overwrite `FI`; nothing else
reads it until the next step's `DENSIT` (verify in Phase 0, Q5). Lens planes
and the I8 density share the same deposit, so they are consistent by
construction.

## 5. Work items

### I1. Exact-z schedule: close targets — **done** (`cz` commit `ae4fa9d`)

`SetExactSteps` (`PMP2main.f90`) now rejects only bitwise-duplicate targets
(relative separation < 1e-6 in a); the 0.5 % merge still moves a scheduled
step onto a target but never a step already pinned (`MarkX = 1`); a
post-check asserts a strictly increasing schedule with consistent `dAlist`
and exactly `Nexact` pinned steps. Tests: `test/precise_z` (65 checks,
`case7_closepairs` runs the four DEGRACE pairs for 140 steps, `case8_pinned`,
`case4c_dup`). Accuracy argument: the kick spans (da_old + da_new)/2 with the
background prefactor at the current a; a 3.3e-4 step after a 6.8e-3 one puts
the kick midpoint 0.2 % of a off, 0.05 % on that one kick, ~2e-6 in velocity;
the legacy 1.5× step growth does the same seven times per run. Cost: one PM
step per inserted target.

### I2. Lens planes

Where: new routine next to `ProjectDENSIT` (`PMP2mod_density.f90` ≈ 308),
called from `Analysis` after `DENSIT` and before `GetPower`.

- Parameters: `Lens slabs per axis` (4 → depth = Box/N_slab = 256 Mpc/h;
  refuse if N_g mod N_slab ≠ 0) and `Lens plane grid` (4096, must divide N_g).
- Algorithm: per axis, slab s = ⌊(i−1)/(N_g/N_slab)⌋ + 1; accumulate `FI + 1`
  over the slab's cells of each column and divide by their number → mean
  (1+δ) of the slab column (Σ/Σ̄); optional block average to the plane grid;
  OpenMP over the outer loop as in `ProjectDENSIT`. Extra memory
  12 × 4096² × 4 B = 0.8 GB.
- Checks logged per output: each plane averages to 1 within float error; the
  4 slabs sum to the full-box projection; timing.
- File: one `.bin` per output holding the 12 planes in fixed order X₁…₄,
  Y₁…₄, Z₁…₄, each N_lens × N_lens float32 column-major, plus `.meta` with
  `N_slab`, `depth`, `N_lens`, `normalisation = mean(1+delta)`.
- Keep `ProjectDENSIT` untouched.

### I3. Output schedule input and generator

- `Init.dat`: `#outputs = -42` followed by the 42 redshifts (exact mode).
  **Never list z_init** (exact mode enforces 0 ≤ z < z_init; the IC file is
  the z = 49 snapshot). `zout(1000)`, `NstepM = 2000` are ample; Real*4
  `zout` (≈ 6e-8) suffices; `PMP2init` writes `es16.8` to `Setup.dat`.
- Generator script (next to the supplement generators in `DEGRACE-pilot/`
  or in `docs/degrace_ecosmog/`): for an ECOSMOG companion take the node's
  list from `docs/degrace_ecosmog/aout.json` and its cosmology from
  `nodes.csv`; for a new cosmology compute the 25 stations from χ(a)
  (matter + Λ, the run's Ω_m), merge the 17 primary redshifts, assert no
  duplicates, and write the `#outputs` block, `Init.dat` cosmology lines and
  the `PkTable.dat` header. Verify it reproduces Baojiu's fiducial list to
  < 1e-7 and every ECOSMOG list to < 1e-6.

### I4. Step → epoch index file

Snapshot files are named by `ISTEP`, irregular once targets are inserted. At
run start, after the `SetExactSteps` schedule print, write
`OutputSchedule.dat`: one line per marked step with `ISTEP a z chi exact
class` (χ matter + Λ). halocat loaders map snapshot numbers from it.

### I5. Per-output content control (one mechanism for everything)

`iBDM`, `iPower`, the lens-plane flag and `EXTRA_OUTPUT` must not all fire at
all 42 outputs by default (budget, section 9). Add one per-output class list
(a second list after the redshift list in `Init.dat`, or classes derived from
the primary/secondary distinction) that the output block consults: e.g.
class 1 = primary (BDM + P(k) + lens planes + extras), class 2 = station
(lens planes + P(k) only). `EXTRA_OUTPUT_LIST` of the brief becomes this
list. Decide the default with Baojiu (section 11).

### I6. Run configuration and provenance

- `Init.dat`: `iSave=0`, `iPower=1`, `iBDM=1`, lens planes on, `Pk tune =
  1.0` (inert for ic2pm runs), z_init 49, da0 by rule 2 of section 0,
  `Ncheckpoint` small enough that losing a 72 h job does not restart from
  z_init (~250 PM steps per box).
- Record per box: `Init.dat`, `OutputSchedule.dat`, the ic2pm log line
  `velocity epoch =`, the IC provenance (generator, seed, z_init) as
  `AGENTS.md` requires, and the `.meta` sidecars.

### I7. Item A: initial Gaussian field dump δ_L(k)

Purpose: a cosmic-variance-free halo propagator Γ_h(k|M) = Σ Re[δ_h δ_L*] /
Σ|δ_L|² and an exact phase reference for every other dump.

Spec (from the brief): dump δ_L(k) for |k_i| ≤ k_c,L (default 1.5 h/Mpc),
normalised to linear theory at a = 1 with the code's own ΛCDM growth; record
`D_code(a_ini)` and the input P_lin table name; f(R) and ΛCDM runs with the
same seed must give identical files (ICs are GR; f(R) linear growth is applied
offline); a regenerate-only mode (`IC_DUMP_ONLY`) that writes the dump and
exits before allocating particles.

**Conflict with rule 1 and how to resolve it.** The brief places the dump in
`PMP2start`, but new runs never use `PMP2start`. Three routes, to decide:

1. FML-generated runs (new boxes): add the δ_L dump to our FML generator in
   `2LPTIC_Gui/` (we control it; same file format, same copy convention).
2. Baojiu's `ics.*` (ECOSMOG companions): either (a) derive δ_L from the IC
   displacement field at ingest in `ic2pm` (first-order inversion
   δ_L = −∇·Ψ on the lattice, with the 2LPT second-order contamination,
   a few % at k ~ 1 and negligible at k < 0.1 at z = 49, quantified by V3),
   or (b) patch his 2LPTic to dump the Gaussian field in regenerate-only
   mode with his `2LPTic.param` and seed 2026 (the code and parameters are
   on cosma7).
3. Legacy pilots (PMP2start ZA at z = 100): the brief's `IC_DUMP_ONLY` mode
   in `PMP2start`, which needs the Phase 0 answers on its RNG (Q1).

### I8. Item B: low-k Fourier dumps of the late-time density and momentum

Purpose: large-scale P_hm(k|M) and halo × momentum cross-spectra for any
marks with no subsample noise; the matter propagator G_m = P_mL / P_L; a D²
calibration against I7; offline environment and tidal fields at R ≥ 3/k_c,δ.

Spec, at each output of the extras class, at the (x, v) pair of the BDM
catalogue:

1. Density: CIC on N_dep³, forward FFT, copy |k_i| ≤ k_c,δ (default 1.5
   h/Mpc → R ≥ 2 Mpc/h; 3.0 if tidal/environment definitions need R ≈ 1,
   e.g. adaptive R = 4 R_200 for low-mass haloes).
2. Momentum: for i = x, y, z deposit p_i(x) = (1/n̄_cell) Σ_p v_{i,p}
   W_CIC(x − x_p) with v the peculiar velocity in km/s and
   n̄_cell = N_p/N_dep³, so p_i = (1+δ) v_i; FFT; copy |k_i| ≤ k_c,p
   (default 1.0 h/Mpc).
3. Buffer: preferred N_dep = N_g reusing `FI` sequentially as in section 4
   (no extra memory); otherwise a separate N_g/2 buffer (34 GB at 2048³) only
   if Q8 shows headroom. If N_dep < N_g and k_c > k_N,dep/4, interlace (two
   deposits shifted by half a cell, combined with the phase factor) and set
   `interlaced = true`.
4. Files: `deltam`, `mom_x`, `mom_y`, `mom_z`, each with `.meta`.

### I9. Item C: ID-hashed, quantised particle subsample

Purpose: small-scale (r ≲ 10 Mpc/h) halo–matter and halo–momentum pair
statistics, validation of I8; Lagrangian positions enable hybrid-EFT fits
with the true displacements.

- Particle ID: ID = i + N(j + N k), 0-based, lattice indices in the
  generator's order. **Known (2026-10-17):** Baojiu's `ics.*` carry IDs
  1…N³ sequential across files 0…27 (checked at every file boundary),
  i.e. file order = lattice order with k fastest; `ic2pm` skips the ID block
  and fills the arrays in file order (verify it preserves that order, and
  that nothing in GLAM permutes the arrays: Q2). The FML generator's order
  must be checked the same way.
- Selection: keep if `splitmix64(ID xor SALT) mod M == 0`, M = 100,
  SALT = `0x5EED00DE6C0A2026` (never change). splitmix64 in a small C file
  via `iso_c_binding`:

```c
#include <stdint.h>
uint64_t splitmix64(uint64_t x) {
    uint64_t z = x + 0x9E3779B97F4A7C15ULL;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}
int glam_select(uint64_t id, uint64_t salt, uint64_t m) {
    return (splitmix64(id ^ salt) % m) == 0;
}
```

  The Python reference (numpy uint64, `np.errstate(over='ignore')`) must
  reproduce it exactly.
- Ordering: write selected particles in increasing ID order with a
  deterministic OpenMP two-pass compaction (count per chunk, prefix sum,
  fill). IDs are implicit: the k-th stored particle is the k-th ID passing
  the rule.
- Quantisation: q = floor(x/L × 65536), x wrapped into [0, L), stored as
  integer(2) holding q − 32768, decode x = (q + 0.5) L/65536 (15.6 kpc/h at
  L = 1024); s = nint(v/Δv), Δv = 0.25 km/s, clipped to ±32767, clipped
  count in `.meta`.
- Layout: one `.bin` with qx, qy, qz, vx, vy, vz (each n_sel × int16,
  contiguous, in that order; 12 B per particle). Meta keys: `N_side`, `M`,
  `SALT` (hex), `n_sel`, `dx_quant`, `dv_quant`, `n_clip_v`, `a_x`, `a_v`,
  `v_def = peculiar km/s`.
- Full quantised snapshot option: `FULL_QUANT_OUTPUTS` lists output indices
  written with M = 1 (~103 GB at 2048³), for one validation node at the main
  emulation redshift only.

### I10. Item D: scalar field at halo centres (f(R) only)

The one quantity not reconstructible offline. After BDM, while the converged
scalar field of this step is in memory, write per halo in catalogue order
(with the halo ID if one exists): `fR_centre` (CIC-interpolated f_R at the
BDM centre) and `fR_mean` (mean over cells within the BDM radius), float32;
record the code-variable → f_R conversion and f_R0 in `.meta`. Off for ΛCDM.

### I11. Optional fallback E: on-the-fly environment and tidal ladder

Only if the I8 density dump at the required k_c,δ is judged too large:
truncate δ(k) to an N_sub³ sub-grid (1024³ = 4.3 GB, or 2048³ = 34 GB for
R = 1), for R ∈ {1, 2, 4, 8} Mpc/h (Gaussian window) and the 6 components
T_ij(k) = (k_i k_j/k²) δ_R(k): inverse FFT, CIC-interpolate to halo centres,
store float32, discard; one sub-grid resident at a time; no
eigen-decomposition in Fortran.

## 6. What NOT to implement in GLAM

Marked correlation functions or spectra, mass binning, rank or Legendre
weighting; eigen-decomposition or any environment definition; Python
readers (they belong in halocat; list reader work as follow-ups in
section 16, do not edit halocat from a GLAM session).

## 7. Runtime flags (all default off; appended to `Init.dat`/`Setup.dat`)

| Flag | Default | Meaning |
|---|---|---|
| `LENS_SLABS` | 0 (off) | slabs per axis for lens planes (4 for 256 Mpc/h at L = 1024) |
| `LENS_GRID` | N_g | plane grid (must divide N_g) |
| `OUTPUT_CLASS` | all outputs class 1 | per-output class list (I5) |
| `EXTRA_OUTPUT` | off | master switch for I8–I10 |
| `KC_LIN` | 1.5 | k_c for I7 (h/Mpc) |
| `KC_DELTA` | 1.5 | k_c for `deltam` (h/Mpc) |
| `KC_MOM` | 1.0 | k_c for the momentum dumps (h/Mpc) |
| `N_DEP` | N_g | deposit grid for I8 (N_g/2 fallback) |
| `SUBSAMPLE_M` | 100 | keep 1 in M particles |
| `SUBSAMPLE_SALT` | 0x5EED00DE6C0A2026 | hash salt, never change between runs |
| `FULL_QUANT_OUTPUTS` | none | output indices for full quantised snapshots |
| `FR_AT_HALOS` | on if f(R), else off | I10 |
| `IC_DUMP` | on when `EXTRA_OUTPUT` is on | I7 during normal runs (where the generator allows) |
| `IC_DUMP_ONLY` | off | regenerate-only mode for the legacy pilots |

Names may be adapted to the `Init.dat` style (`key = value` with comments).

## 8. Validation tests

Schedule and lens planes:

| ID | Test | Pass criterion |
|---|---|---|
| S1 | `test/precise_z` suite | 65/65 (already passing) |
| S2 | generator reproduces Baojiu's fiducial list and all 96 ECOSMOG lists | < 1e-7 / < 1e-6 in a |
| L1 | plane means and slab sums on the test box | mean = 1 to float error; 4 slabs = full projection |
| L2 | planes vs. offline projection of a full snapshot (test box, `iSave=1`) | identical to float rounding |
| L3 | plane power spectrum vs. the 3-D P(k) at k ≪ k_N (Limber-free check on the test box) | within shot noise |

Extra outputs (from the brief):

| ID | Test | Pass criterion |
|---|---|---|
| V1 | layout unit test: cos(2π m·x/L), sin(…) for a few integer m through the copy routine | amplitude 1/2 at ±m only, correct signs; zero elsewhere to round-off |
| V2 | P(k) from `deltam` (CIC deconvolved offline) vs the code's own P_mm(k), same binning, k < k_c | < 1e-3 |
| V3 | `deltam` at the first step vs `deltaL` | Re(ratio) = D_code(a_ini) to < 0.5 % at k < 0.1; \|Im/Re\| < 1e-3; r > 0.999 |
| V4 | halo axes/offset: BDM haloes (M ≳ 10¹³) deposited offline, crossed with `deltam` and `deltaL` | r ≳ 0.9 in the lowest bins; no linear-in-k phase |
| V5 | subsample on the test box with a full-precision snapshot | n_sel within 3σ of N/M; reconstructed IDs match to dx_quant/2; n_clip_v = 0 or logged; C and Python hashes agree on 10⁷ IDs |
| V6 | momentum spectrum from the subsample vs `mom_*`, k < 0.3 | within subsample shot noise |
| V7 | determinism of the regenerate-only mode, twice and with another thread count | identical md5, or thread dependence documented and enforced |
| V8 | f(R) vs ΛCDM, same seed: `deltaL` | identical md5 |
| V9 | overhead at production size (one short job, approval first) | extra wall time < 5 %; peak RSS within node memory with margin |

## 9. Storage budget (L = 1024 Mpc/h)

Mode count (2n+1)²(n+1), n = floor(k_c L/2π), 8 B per mode.

| product | per output, 2048³ | per output, 1024³ |
|---|---|---|
| lens planes (12 × 4096²) | 0.8 GB | 0.8 GB (or 0.2 GB at 2048²) |
| `CatshortV` catalogue | 0.3–0.4 GB | ~0.1 GB |
| `deltam`, k_c = 1.5 | 0.47 GB (3.7 GB at k_c = 3.0) | same |
| `mom_x,y,z`, k_c = 1.0 | 0.41 GB | same |
| subsample, M = 100 | 1.03 GB | 0.13 GB |
| f_R at haloes | ~8 B per halo | same |
| `deltaL` (once per run) | 0.47 GB | same |
| reference: float32 snapshot | 206 GB | 26 GB |

Per box of 42 outputs at 2048³: lens planes 34 GB, catalogues ~15 GB,
extras at all 42 outputs ~80 GB (≈130 GB per box, 13 TB for 100 boxes);
extras at the 17 primary outputs only ~32 GB (≈80 GB per box, 8 TB for
100 boxes). Baojiu's Slack budget covered lens planes only (3.4 TB).
Decision needed (section 11). Optional offline `zstd` gains are modest.

## 10. Phasing

Two tracks; order to be confirmed by CZ. Each phase ends with code on the
branch, a short test report with numbers, and an entry in section 16.

Track A, campaign readiness (cheap, needed for any run):
1. Phase 0: reconnaissance (section 13), no code edits; report.
2. I3 generator + S2; I4 index file; I5 class list.
3. I2 lens planes + L1–L3.
4. Pilot node 2 (ECOSMOG design centre): ingest its ICs, run ΛCDM and f(R)
   with the 42-output schedule, compare P(k) with the ECOSMOG twins at the
   common epochs, measure runtime, `sacct` CPU efficiency and file sizes.

Track B, extra outputs (the brief's phases):
5. I7 (route 1 or 2 decided) + V1, V7, V8; then the first science product:
   δ_L for one existing pilot via the regenerate-only route and V4 against
   its halo catalogues (unlocks Γ_h for the pilots).
6. I8 + V2, V3, V4 on the test box.
7. I9 + V5, V6.
8. I10 on the f(R) test box.
9. V9 production-size overhead run, after approval.

## 11. Open questions

For Baojiu:
1. BDM catalogues at all 42 outputs or only the 17 primary redshifts?
2. Lens planes at the primary redshifts too, or only at the 26 stations
   (z = 0 and k = 1…25)? His 34 GB budget assumes all 42.
3. Plane normalisation: mean (1+δ) of the slab column (proposed) or surface
   density; must match his ECOSMOG planes.
4. Keep both members of each close pair (now supported) or drop one? They
   differ by ≤ 1.5 Mpc/h in χ and 0.03 % in growth: a science choice.
5. The questions in `docs/degrace_ecosmog/README.md` section 9 (CAMB ω_c,
   r2g conversions, Node_900, ΛCDM P(k)).

For CZ:
6. Extras cadence: all 42 outputs or primary only (section 9).
7. k_c,δ = 1.5 or 3.0 (do tidal/environment definitions need R ≈ 1 Mpc/h)?
8. I7 route for the FML generator and for Baojiu's ICs (section 5, I7).
9. File naming: the brief's `extra/<run>_a<a>_<product>` or GLAM's
   `<product>.<step>.<box>.dat`.
10. Companion resolution: N_g = 2048 or 4096 on Baojiu's 1024³ particles.

## 12. Downstream use (orientation only; not implemented in GLAM)

- Halo propagator Γ_h(k|M) = Σ_shell Re[δ_h δ_L*] / Σ_shell |δ_L|², with δ_L
  scaled to the output redshift (f(R): scale-dependent D(k, z)); divide by
  the realised Σ|δ_L|², never by a theory P_L; fit b1 + a k² at
  k ≲ 0.05–0.08 h/Mpc.
- Growth calibration D²_eff = Σ|δ_m|² / Σ|δ_L|² in the lowest shells.
- Marked statistics: halo weights (2a+1) P̃_a(u), u = rank within the mass
  bin, a ≤ 1–2; large scales from I8 (k-space), small scales from I9 pair
  counts (Corrfunc), spliced at ~10–30 Mpc/h.
- Environment and tidal fields offline from `deltam`, R ≥ 3/k_c,δ,
  T_ij = (k_i k_j/k²) δ_R, interpolated to halo centres.
- k-space factorisation P_hh(k|M1, M2) = P_hL(M1) P_hL(M2)/P_L + P_⊥, the
  residual emulated only at r ≲ 10 Mpc/h.

## 13. Phase 0 reconnaissance questions (answer from the code, record in section 16)

1. **IC generator.** Only for the legacy-pilot regenerate route (I7 route 3):
   where `PMP2start` generates the Gaussian field, in which array, δ(k) or
   potential, normalisation and growth factor, RNG seeding, thread
   dependence (record the pilots' thread count). For the FML generator: the
   same questions in `2LPTIC_Gui/`.
2. **Particle ordering.** Does anything in GLAM sort or permute the particle
   arrays after `ic2pm` fills them (BDM, linker lists, `ParticlesRandom`)?
   Which loop order maps (i, j, k) to the array index and to x, y, z in
   `ic2pm` and in the FML generator?
3. **Units.** Positions in grid units [1, N_g+1) (`ic2pm`: X = pos·N_g/Box
   + 1); the momentum code units and the conversion to peculiar km/s used by
   the snapshot or halo writer (reuse it).
4. **Time synchronisation.** Confirm v at a − da/2 at outputs and that BDM
   uses that pair; record `ASTEP` per output (`PMcrd` header).
5. **Density buffer.** Confirm the `Analysis` order of section 4, that
   `GetPower` FFTs `FI` in place, FFTPACK5's half-complex layout and
   normalisation, and that `FI` is free between the P(k) FFT and `BDM`.
6. **Scalar field.** Which array holds the converged scalar field at level 0
   after the multigrid solve (`FI2`/`FI3`), the conversion to f_R, and
   whether it is still valid when BDM runs (note `MG_scratch_release`
   around outputs).
7. **BDM.** Call site (end of `Analysis`), catalogue order, halo ID column,
   which centre is written (needed for I10).
8. **Memory.** Peak RSS of a production f(R) run at N_p = 2048³, N_g = 4096
   (`sacct` MaxRSS of a past job). Slack: "a 4096 modified-gravity run takes
   950+ GB"; ΛCDM pilot ~206 GB; L1800/N_g 6400 1099 GB. Headroom on a
   1 TB node decides N_dep.
9. **I/O and parameters.** `Init.dat` sequential parser + `Setup.dat`
   (section 3); output locations (pilot data live under the model
   directories in `MG-GLAM/`, e.g.
   `DurMun_hmfemu_*_model65_L1024Np2048Ng4096/Run{box}/`, not under
   `/cosma8/data/dp203/dc-ruan1/degrace_pilot/`; verify).

## 14. Checklist

- [x] Cross-check of the 42 scale factors (2026-10-16)
- [x] I1 close-pair fix + tests (`cz` `ae4fa9d`, 2026-10-16)
- [x] ECOSMOG suite metadata pinned (`docs/degrace_ecosmog/`, 2026-10-16)
- [ ] Phase 0 answers in section 16
- [x] I3 generator + S2; I4 `OutputSchedule.dat`; I5 class list (`degrace-track-a` `cb09e3c`)
- [x] I2 lens planes + L1, L2 (`degrace-track-a` `b774f8d`); L3 (plane P(k) vs 3-D P(k)) pending
- [ ] Pilot node 2 vs the ECOSMOG twins (ΛCDM job 12180066 submitted 2026-10-17; f(R) next)
- [ ] I7 δ_L dump (route decided) + V7, V8; δ_L for one legacy pilot + V4 (V1 layout test exists)
- [x] I8 + V1, V2 (`degrace-track-a` `a41ae37`); V3 needs I7, V4 needs the pilot catalogues
- [x] I9 + V5, V6 (`degrace-track-a`, commit after `09c4aa4`)
- [x] I10 (`degrace-track-a`, commit after `a41ae37`); unit test + f(R) smoke case
- [ ] V9 overhead run
- [ ] Reply to Baojiu with the cross-check result and questions 1–5

## 15. Related memory

`project_degrace_lensplanes_plan`, `project_degrace_ecosmog_suite`,
`project_exact_z_snapshots`, `project_stepping_for_2lptic`,
`feedback_ic_generator`, `project_2lptic_ic_integration`.

## 16. Notes (append findings here, newest last)

- 2026-10-17: plans merged. Facts already established in the
  2026-10-16/17 sessions, to be re-verified at implementation time:
  (Q2) Baojiu's `ics.*` are in lattice order with sequential IDs 1…N³
  (checked at all 28 file boundaries); `ic2pm` reads POS/VEL only.
  (Q4) outputs are not synchronous: x at a, v at a − ASTEP/2
  (`project_exact_z_snapshots`); exact-z mode makes ASTEP vary per output.
  (Q5) `DENSIT` writes δ into `FI` (initialised to −1, CIC adds);
  `GetPower` FFTs `FI` in place and sets `mDENSIT = 1`; `Analysis` runs
  PDF → P(k) → dump/project → bias → `BDM(mDENSIT)`; `ProjectDENSIT`
  (full-box projections) exists under `iDumpPart`. (Q6) `FI2`, `FI3` are
  allocated for MG; `FI2` enters `MOVE` as the fifth-force factor; the
  f_R mapping is unverified. (Q8) Slack: f(R) N_g = 4096 needs 950+ GB.
  (Q9) `Init.dat` is parsed sequentially (`iParseLine`); `Setup.dat` is
  written by `PMP2init` and read by `ReadSetup`; stale trailing lines in
  `Setup.dat` are ignored. I7 conflicts with the 2LPTIC rule; see I7.
- 2026-10-17, Phase 0 answers (four read-only reconnaissance passes over the
  root sources; line numbers refer to `cz`):
  - **Q1 (legacy ICs)** `PMP2start` builds Zel'dovich displacements in
    `SPECTR` (`GRX/GRY/GRZ = TS k_i/k²` in a real cos/sin basis, not δ or the
    potential), RANLUX lux 2 seeded from `TableSeeds.dat` row `Nrealization`,
    one saved state per k-plane restored inside the OpenMP loop: the random
    numbers do not depend on the thread count (only the real*8 reduction
    of the normalisation sum, at ~1e-15). Normalisation uses the realised
    total power (`ALPHA = AMPLT/sqrt(SUMM)`), not the table.
  - **Q2 (ordering)** Nothing permutes the particle arrays: `MOVE`, density
    deposits, BDM (`move_alloc` copies, ghost rows appended after `Np`, index
    sorts only, originals restored), `ParticlesRandom`, `WriteDataPM`/
    `ReadDataPM` (pages in index order). `ic2pm` fills the arrays in Gadget
    file order and never reads the ID block; Baojiu's `ics.*` carry
    sequential IDs 1…N³ across files 0…27 (lattice order, k fastest), so
    file order = lattice order for the ECOSMOG companions. `PMP2start`
    lattice: `Icurrent = ip + (jp-1)·NROW + (kp-1)·NROW²`, ip → x.
    Known bug: `ParticlesRandom` sets `Nparticles = NROW³` in int4 (overflows
    for NROW ≥ 1291; only with `iDumpPart = 1`).
  - **Q3 (units)** Positions in grid units [1, NGRID+1), node i at
    (i−1)·Box/NGRID; `VX` is the comoving momentum with x in grid units; every
    existing writer (BDM, `ParticlesRandom`, `DumpParticles`) uses
    v[km/s] = VX·(100·Box/NGRID)/AEXPN although VX sits at a − ASTEP/2, so
    their km/s are low by a_v/a; the new products must divide by a_v and
    record `a_v`. `ic2pm`: V = v_pec·a·NGRID/(100·Box) via
    `vfac = S_vel·a^1.5·NGRID/(100·Box)`.
  - **Q4 (epochs)** Outputs are not synchronous: x at a, v at a − ASTEP/2
    (ASTEP = the step just taken, in the `PMcrd` header and the BDM catalogue
    header `Step=`). BDM uses the arrays as they stand.
  - **Q5 (density buffer)** `DENSIT`: CIC into `FI` initialised to −1 with
    weight NGRID³/Np, so FI = δ with zero mean; `FI(ix,iy,iz)`, ix fastest.
    `GetPower` → `POWERfft5`: three real 1-D `rfft1f` passes (x, y, transpose,
    z), √2 scaling of indices 1 and N per pass, then **squares the
    coefficients in place** (`FI = FI²/8`) and zeroes k = 0, so after
    `GetPower` FI holds neither δ nor the coefficients. Low-k mode dumps (I8)
    must be copied inside `POWERfft5` between the FFT and the squaring, or
    use their own transform. Packing per dimension (even N, phase index
    j−1): p=1 → (1/N)Σx; p=2m → (2/N)Σx cos; p=2m+1 → +(2/N)Σx sin; p=N →
    Nyquist; after the transpose the coefficient sits at `FI(pz,py,px)`.
    `Analysis` order: PDF → P(k) (`mDENSIT = 1` after) → dump/project →
    bias → `BDM(mDENSIT)` (re-deposits if needed). Lens planes now run
    first, from the fresh deposit, and leave `mDENSIT = 0` for the power
    spectrum. `MG_scratch_release` frees only `FI3`; `FI2` is kept.
  - **Q6 (scalar field)** `FI2` holds u = |f_R|^{1/(n+1)} normalised by
    u_bg after `scalar_field`; f_R(x) = fR_bg(a_old)·FI2^{n+1} with
    fR_bg = fR0·(R0_bg/R_bg)^{n+1}, R_bg = 3(Ω_m/a³ + 4Ω_Λ), evaluated at the
    pre-step a_old (FI2 was solved before MOVE/ADDTIME). FI2 is intact at
    output time (only the MG solve writes it). The `fric_fac·FI2` term in
    MOVE is symmetron-only; the f(R) force enters through `POTENT_FOR_MG`.
  - **Q7 (BDM)** Called at the end of `Analysis`; `CatshortV.<step>.<box>.DAT`
    has an 8-line header (A, Step=ASTEP, cosmology, mass, Ovdens) and 24
    columns: x,y,z [Mpc/h], Vx,Vy,Vz [km/s], Mvir, Mtotal, Rvir [kpc/h],
    Vrms, Vmax, iHalo (row index, per file only), Cvir, Npart, 0, Xoff,
    2K/Ep−1, λ, RadRMS, b/a, c/a, axis direction. Rows in `FindMaxima` peak
    order (z, y, x), not by mass. Centre = density peak refined by four
    centre-of-mass iterations; radius = Bryan–Norman virial aperture.
    Per-halo arrays (`xMaxx…`, `Rvir`, `Mvir`) live in module `Structures`
    and are freed by `ReleaseMaxima` right after `WriteFiles`, so I10 must
    hook in between (inside `BDM`, `PMP2linker.f90` ≈ 208–217).
  - **Q8 (memory)** `sacct` MaxRSS: f(R) L1024 Np2048 Ng4096 (`fR_m65_b1`,
    job 11448712) **985 GB**, 16 h 46 m on 128 cores; ΛCDM twin 489 GB,
    6 h 20 m; the L512 2LPTIC reference 655 GB. On a 1005 GB node f(R) at
    N_g = 4096 has ~20 GB headroom: no new N_g³ buffer, and not even the
    N_g/2 fallback; FI reuse is mandatory there. During `Analysis` the MG
    run has `FI3` freed (one N_g³ block) but `FI2` must stay.
  - **Q9 (parameters)** `Init.dat` is parsed sequentially by
    `ParseLine`/`iParseLine` (last '=' on the line, fixed 12/10-char field,
    no `iostat`); `Setup.dat` is written by `PMP2init` and read by
    `ReadSetup` in the same order, which stops after `csf_beta`. The new
    optional block ("extra outputs") sits after `csf_beta` in both files
    and is read only when its header line is present. Outputs land in the
    run directory (`Run1/`); executables are addressed by absolute path;
    `PMP2main` takes only the step count on stdin. The makefile has no
    dependency rules (full rebuild after a module change) and no C rule
    (needed for I9's splitmix64).
- 2026-10-17, Track A implemented on `degrace-track-a`:
  - `cb09e3c` I3 `DEGRACE-companion/gen_inputs.py` (S2: fiducial list to
    7.7e-8, all 96 ECOSMOG lists to 1.3e-7); I4 `OutputSchedule.dat`
    (module `ExtraOutputs`, `PMP2mod_extras.f90`); I5 per-output classes
    (optional trailing `Init.dat`/`Setup.dat` block, class carried through
    the sorted exact-z targets, BDM only for class 1).
  - `b774f8d` I2 lens planes (`LensPlanes`, same module): real*8
    accumulation, threads own plane rows, one deposit shared with P(k);
    `.bin` + `.meta` per output. Tests `test/precise_z`: 146/146 (new
    cases 9, 10_lens72, 10_lens36, 4e; `OutputSchedule.dat` checked for
    every case). L3 (plane P(k) vs 3-D P(k)) not yet written.
  - Pilot: `DEGRACE-companion/runs/Node_002_LCDM_L1024Np1024Ng2048/`
    (node 2 = design centre, Baojiu's ICs via `ic2pm half`, Ng 2048, da0
    6e-4, classes primary = 1 / station = 2, lens planes 4 × 2048²,
    frozen binaries from `b774f8d`), Slurm job **12180066** on
    `cosma8-serial`, 64 cores, 200 GB, 8 h, account dp203; pending at
    submission (nodes reserved). Expected ~3 h once running. f(R) twin to
    follow once it completes.
  - Cosmetic bug fixed on the branch: `PMP2init` printed the exact-z list
    with a recycled `(a,10f10.5)` format, which emitted raw float bytes
    into the log for more than 10 redshifts.
- 2026-10-17, Track B on `degrace-track-a`:
  - `a41ae37` I8: `FFT3Dforward` (copy of `POWERfft5`'s passes), `DumpModes`
    (packed → complex cube, V1 layout test passes: cos → +½ real at ±m,
    sin → −¼i at +m), `MomentumDeposit`, `ExtraModes`; density modes
    hooked into `POWERfft5` before the squaring (no extra FFT), momentum
    on a transient N_dep³ grid (default N_g/2). V2: the `deltam` cube
    reproduces `PowerDM.R` bin by bin to 1e-4 on the test box. New flags
    `Extra outputs`, `kc_delta`, `kc_mom`, `N_dep_mom`, `Subsample M`,
    `fR at halos` (optional `Init.dat` lines, `Setup.dat` block extended;
    older blocks keep the defaults). `PMP2start`/`PMP2spectrum` now link
    the extras module because `Power` uses it.
  - I10: `fR_halo.<step>.<box>.bin` written from inside BDM's row loop
    (`PMP2linker.f90`, 3 lines) with `FI2` and the pre-step epoch; unit
    test on a synthetic ramp; f(R) smoke case on the 36³ box (BDM found
    1 and 5 haloes at the two outputs, rows match).
  - Pilot job 12180066 failed at once (`$0` resolves to Slurm's spool
    copy; now `SLURM_SUBMIT_DIR`); resubmitted as **12180082** on m8504,
    running at ~0.64 min per step on 64 cores (163 steps ≈ 1.8 h plus
    outputs). `Run1/CATALOGS` had to be created by hand while it ran.
  - Known: `PowerDM`/extras velocities use a_v = a − da/2 (recorded in
    every `.meta`); the legacy BDM catalogue divides by a instead.
  - I9 integrated: `PMP2mod_subsample.f90` + `glam_hash.c` (written by a
    subagent to the brief's spec, new files only), wired into `Analysis`
    for class-1 outputs; V5 unit test (hash agreement on 1e7 IDs, order,
    quantisation, clipping) and V6 in the suite (full-mode subsample
    reproduces the `deltam`/`mom_x` spectra to 1e-3). Lesson recorded
    for the parameter files: a comment containing `=` after a value
    breaks `ParseLine` (it keys on the last `=`); the generator's
    comments avoid it.
