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
(section 16); **I2–I5 and I7–I10 are implemented and merged into `cz`**
(2026-10-17; 180/180 checks plus four unit-test scripts); both node-2
pilots are complete and verified (ΛCDM job 12180082, f(R) job 12180110,
section 16). I6 is encoded in the input generator and launcher; I7 is
done (route 1); I11 is open.

Related references: `docs/degrace_ecosmog/README.md` (Baojiu's ECOSMOG
suite: design, ICs, output lists, snapshots), `2LPTIC_Gui/VALIDATION.md`
(2LPTIC and time-stepping validation), review decks kept locally under `docs/slides/` (not tracked).

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
lines. Keys as implemented (the brief's list, amended 2026-10-18 to what
the writers emit and `glamio` reads): `product`, `format`, `step`,
`realization`, `a`, `z`, `a_x`, `a_v` (products that carry velocities;
`fR_halo` records `a_old` instead), `Box`, `units`, `axis_map`, `N_grid`,
`N_row`, `N_p`, `N_dep` (cubes), `D_code(a)` (cubes; the code's own ΛCDM
linear growth at this a, normalised to D(1) = 1), `Omega_m`, `Omega_L`,
`hubble`, `MG_flag`, `MG_model`, `fR0` (or `none`), `header` (the 45-char
run header), plus the item-specific keys below. No `git_hash`: the
campaign's provenance file records the build. Names: GLAM's
`<product>.<step:04d>.<box:04d>.{bin,meta}` in the run directory (the
brief's `extra/<run>_a<a>_<product>` was not adopted; Q9 decided).

**Velocity epoch.** GLAM's leapfrog leaves velocities half a step behind the
positions at every output (v at a − da/2 with the step just taken; in
exact-z mode that step can be as short as 3.3e-4). The new products that
carry velocities (I8 momentum, I9 subsample) use the same (x, v) pair and
record `a_x` and `a_v`. The legacy BDM catalogue is unchanged on this
branch and converts its velocities with 1/a instead of 1/a_v
(`PMP2linker.f90`, `Vscale = 100 Xscale/AEXPN`), so its velocities are low
by a_v/a: 0.8 % at z = 2 (the first catalogue output, Δa = 0.0051), 0.3 % at
z = 0 (Δa = 0.0054); V_rms, the virial ratio and λ carry the same offset
(corrected 2026-10-18 after review: earlier text had the direction inverted).
Correcting it would change every GLAM catalogue, not only the companions;
the choice is listed in section 11 and the catalogue readers must multiply
by a/a_v meanwhile (`glamio` `velocities_at_av()` does).

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
- Checks logged per output: the mean over the slabs of each axis is 1 within
  float error (a single plane averages to 1 + ⟨δ⟩_slab, its slab's mean
  overdensity: up to 1e-3 on the pilots); the
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

**Decision (CZ, 2026-10-17): route 1.** The dump is implemented in our FML
generator (`2LPTIC_Gui/`); Baojiu's ECOSMOG ICs and the legacy pilots are
out of scope for I7. The other routes stay documented for reference:

1. FML-generated runs (new boxes): add the δ_L dump to our FML generator in
   `2LPTIC_Gui/` (we control it; same file format, same copy convention).
2. Baojiu's `ics.*` (ECOSMOG companions): either (a) derive δ_L from the
   ICs at ingest in `ic2pm`. The naive version takes δ_1 = −∇·Ψ with
   Ψ = x − q on the lattice and carries the 2LPT term Ψ_2 as a
   contamination (negligible at k < 0.1, a few per cent near k ~ 1 at
   z = 49). The exact version uses positions and velocities together:
   2LPTic writes x − q = Ψ_1 + Ψ_2 and v/(aH) = f_1 Ψ_1 + f_2 Ψ_2 with its
   own f_1 = Ω_m(a)^{5/9}, f_2 = 2 Ω_m(a)^{6/11}, so
   Ψ_1 = (f_2 (x − q) − v/(aH)) / (f_2 − f_1) per particle, then
   δ_L(a = 1) = −∇·Ψ_1 / D(a_ini) by one FFT of the lattice field.
   Or (b) patch his 2LPTic (source on cosma7) to dump its Gaussian field
   in a regenerate-only mode with his `2LPTic.param` and seed 2026.
3. Legacy pilots (PMP2start ZA at z = 100): the brief's `IC_DUMP_ONLY` mode
   in `PMP2start`, which needs the Phase 0 answers on its RNG (Q1).

**Implemented (2026-10-17, `degrace-track-a`).** `DeltaLDump.h` + three
hunks in `Main.cpp` of the FML LPT example (tracked copies in
`2LPTIC_Gui/fml_intel_port/`, see its `APPLY.md`):

- Lua: `write_deltaL` (default false; on in the smoke-test and production
  inputs), `deltaL_kc` (default 1.5 h/Mpc), `deltaL_output_prefix`
  (default `output_prefix`). Output `<prefix>_deltaL.bin` + `.meta` next
  to the Gadget files, written right after the Gaussian field is drawn,
  before any LPT displacement.
- Content: δ_k = N⁻³ Σ_x δ(x) e^{−ik·x} of the field at z_ini divided by
  D(a_ini) (flat matter+Λ Simpson integral, D(1) = 1, the formula of
  GLAM's `GrowthFactor`), i.e. linear theory at a = 1; cube |ix|,|iy| ≤ n,
  0 ≤ iz ≤ n, n = ⌊k_c L/2π⌋ (n = 244 at k_c = 1.5, L = 1024: 0.47 GB),
  complex64 `C(0:n,−n:n,−n:n)` (iz fastest), k = 0 zeroed: the `deltam`
  layout exactly. FML's FFTWGrid already uses the N⁻³ normalisation and
  node j at j L/N, so the copy is index-exact. The `.meta` records n,
  k_c, N_mesh, N_part_1D, Box, z_ini, D_gen(a_ini), Ω_m, Ω_Λ, h, seed,
  lpt_order, fix_amplitude, reverse_phases, the P_lin file, the lattice
  convention (particles at i L/Np, ID = 1 + ix + Np(iy + Np iz)), the
  axis map and the cube rms.
- MPI: each task fills its slab's modes into a full float cube already
  divided by D (D is identical on all tasks), then a chunked `MPI_Reduce`
  (exact: one non-zero contribution per element) gathers on task 0; 8 B
  per mode per task (0.47 GB per rank at production size, 30 GB per node
  at 64 ranks).
- Tests (`test/extras/run_deltaL_test.sh`, 75 s on a login node):
  V7 identical md5 with 1, 2 and 4 ranks; layout check A on Zel'dovich
  ICs (δ_L = −ik·Ψ(k)/D from the lattice displacements, exact: residual
  rms 8e-6 of the cube rms, every mode within 1e-3, x↔y / conjugate /
  half-cell hypotheses off by 10⁴×); A on 2LPT ICs (residual 8e-3, the
  second-order displacement); B (CIC density of the particles, half-cell-
  shifted deposit, deconvolved: mean ratio 1.0002); V3 end to end (FML IC
  → `ic2pm half` → PMP2main 128³/256³ L256 with extras, outputs at
  z = 48 and 30): ⟨deltam/W_CIC/δ_L⟩/D_code(a) = 0.9997–1.0012 in the
  bins k ≤ 0.15, |Im| ≤ 0.2 %, r_cross ≥ 0.9995 (0.997 at 0.2–0.3), and
  the per-axis coordinate offset inferred from the mean phase is
  ≤ 0.03 Mpc/h (half a cell = 0.5), so `ic2pm`'s `x = pos N/L + 1` and
  GLAM's node i at (i−1) L/N introduce no shift. V8 holds by
  construction: the generator has no gravity-model input, and f(R)/ΛCDM
  twins share one IC file.
- Caveats found while testing: the CIC density of particles that sit
  exactly on deposit-grid nodes is not linear in the displacement (the
  weight 1 − |Ψ|/Δ gives one-sided differences), which leaves per-mode
  phase scatter ~ kΔ/2 (10 % at k = 0.5 on the 256 grid) but no bias in
  the mean to O((kΔ)²); every GLAM output from a lattice IC carries this at
  early times. Offline cross-spectra should therefore use binned means or
  a TSC/half-cell-shifted deposit, never single-mode ratios. The
  regenerate-only mode of the brief is not needed: the dump costs ~1 s and
  runs inside the normal IC job; a dump-only run is the generator with
  `write_gadget = false`.

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
   `interlaced = true`. **Status (2026-10-17): interlacing is not
   implemented.** `PMP2init` refuses k_c,mom > k_N,dep/4 and the default is
   k_c,mom = 0.75 h/Mpc (k_N,dep/4 = 0.785 at N_dep = 1024, L = 1024); every
   `.meta` writes `interlaced = false` truthfully. The node-2 f(R) pilot was
   launched with k_c,mom = 1.0 before this rule existed, so its `mom_*` modes
   above k = 0.78 h/Mpc carry uncorrected aliasing. The transient N_dep³
   buffer (not FI reuse) is what the code does; at N_g = 4096 f(R) it would
   not fit (Q8) and must be revisited before such a run.
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
weighting; eigen-decomposition or any environment definition. Python
readers of the raw formats live in `python/glamio` of this repository
(decision 2026-10-17, amending the earlier "readers belong in halocat"):
they track the Fortran writers and the test suite validates them against
real output. halocat keeps its catalogue layer (hdf5 cache, snapshot maps)
and should import glamio; that integration is a halocat-session task, do
not edit halocat from a GLAM session.

## 7. Runtime flags (as implemented: optional trailing `Init.dat` lines, copied by `PMP2init` into a `Setup.dat` block)

All default off; a missing line means the default. The value is the first
token after the last `=` of the line, so comments must not contain `=`.
A present but unreadable line stops `PMP2init` (since the review fixes of
2026-10-17; before, it silently took the default).

| `Init.dat` key | Default | Meaning |
|---|---|---|
| `Output classes` + list line | 0 (all outputs class 1) | number of classes = `#outputs`, then one class per redshift (1 full, 2 station), I5 |
| `Lens slabs per axis` | 0 (off) | slabs per axis for the lens planes (4 for 256 Mpc/h at L = 1024), I2 |
| `Lens plane grid` | 0 (= N_g) | plane grid, must divide N_g |
| `Extra outputs` | 0 | master switch for I8–I10 |
| `kc_delta` | 1.5 | k_c of `deltam` (h/Mpc), deposit on N_g |
| `kc_mom` | 0.75 | k_c of `mom_*` (h/Mpc); must be below k_N,dep/4 (no interlacing) |
| `N_dep_mom` | 0 (= N_g/2) | momentum deposit grid, 8 … N_g |
| `Subsample M` | 0 (off) | keep 1 in M particles (1 = every particle, quantised) |
| `fR at halos` | 0 | f_R at BDM halo centres; forced to 0 unless `MG_flag = 1` and `MG_model = 1` |

One flag lives in `BDM.config` instead (the halo finder's own file, parsed
by `ReadParameters` in `PMP2linker.f90`; branch `bdm-binary-catalogue`,
2026-10-17):

| `BDM.config` key | Default | Meaning |
|---|---|---|
| `CatalogueFormat` | `binary` (compiled default since 2026-10-18, user decision; `gen_inputs.py --catalogue-format both` or `text` while a consumer still needs the `.DAT`) | `text`: legacy `CatshortV.<step>.<box>.DAT`; `binary`: `CatshortV.<step>.<box>.bin` (96-byte little-endian records, 24 four-byte fields = the 24 text columns, 22 float32 + 2 int32, same rows and order) plus a `key = value` sidecar `.meta` (columns, dtypes, units, rows, a, z, a_v, da, cosmology, particle mass, MassMin, overdensity, iVirial, finder, header); `both`: all three files. Any other value stops the run at the first catalogue step. All files are staged as hidden temporaries and renamed into place text → bin → meta, so a published `.meta` marks a complete `.bin`. |

The companion input generator writes `Run1/BDM.config`; the key is written
only when `--catalogue-format` is given (a build older than 2026-10-18
rejects the key, so the generated inputs stay usable with the frozen pilot
binaries, which write the text). Without the key a current build writes
the binary only; pass `both` or `text` while a consumer (halocat, until it
imports glamio) still needs the `.DAT`. Since 2026-10-18 `PMP2main` validates
`BDM.config` at start-up when `Find BDM halos = 1` (`CheckBdmConfig`), so a
typo or an unknown key stops the run before any compute instead of at the
first catalogue step. A finder killed between staging and publication
leaves hidden `CATALOGS/.Catshort*.tmp.*` files (up to three per output);
delete them when no finder is running (`check_run.py` reports them).
Readers: `python/glamio` (`read_catalogue(stem)`, text or binary,
`velocities_at_av()` for the a/a_v correction below).

Not implemented as flags: the hash salt (constant `0x5EED00DE6C0A2026` in
`PMP2mod_subsample.f90`), per-output full quantised snapshots
(`FULL_QUANT_OUTPUTS` of the brief: only the global `Subsample M`), the
regenerate-only IC dump (`IC_DUMP_ONLY`: not needed, the FML dump costs ~1 s
inside the IC job). The δ_L dump is a Lua switch of the FML generator
(`write_deltaL`, default false; `deltaL_kc`, default 1.5), independent of
the GLAM flags.

## 8. Validation tests

Schedule and lens planes:

| ID | Test | Pass criterion |
|---|---|---|
| S1 | `test/precise_z` suite | 65/65 (already passing) |
| S2 | generator reproduces Baojiu's fiducial list and all 96 ECOSMOG lists | < 1e-7 / < 1e-6 in a |
| L1 | axis means and slab sums on the test box | per-axis mean = 1 to float error (a plane's own mean is 1 + ⟨δ⟩_slab); 4 slabs = full projection |
| L2 | planes vs. offline projection of a full snapshot (test box, `iSave=1`) | identical to float rounding |
| L3 | plane power spectrum vs. the 3-D P(k) at k ≪ k_N (Limber-free check on the test box) | within shot noise. **Passed 2026-10-17** (`case10_lens72`): P_2D(k_⊥) of the six planes vs the exact periodic-box projection of `PowerDM.R`, (1/L) Σ_kz P_3D \|W_slab(kz)\|², agrees within 3σ in every ring 0.08 < k_⊥ < 0.5 and on average to 7/4/3 % at the three outputs (sample variance of a few modes per ring) |

Extra outputs (from the brief):

| ID | Test | Pass criterion |
|---|---|---|
| V1 | layout unit test: cos(2π m·x/L), sin(…) for a few integer m through the copy routine | amplitude 1/2 at ±m only, correct signs; zero elsewhere to round-off |
| V2 | P(k) from `deltam` (raw CIC, no deconvolution: the code's P_mm(k) is raw too), same binning, k < k_c | < 1e-3 (as run: 4e-6, `Modes.power()`) |
| V3 | `deltam` at the first step vs `deltaL` | Re(ratio) = D_code(a_ini) to < 0.5 % at k < 0.1; \|Im/Re\| < 1e-3; r > 0.999. **Passed 2026-10-17** on the 128³/256³ test box at z = 48 and 30 (ratio/D = 0.9997–1.0012 for k ≤ 0.15, ⟨Im⟩ ≤ 2e-3, r = 0.99998 at k ≤ 0.05, 0.9995 at 0.1–0.15; see I7) |
| V4 | halo axes/offset: BDM haloes (M ≳ 10¹³) deposited offline, crossed with `deltam` and `deltaL` | r ≳ 0.9 in the lowest bins; no linear-in-k phase |
| V5 | subsample on the test box with a full-precision snapshot | n_sel within 3σ of N/M; reconstructed IDs match to dx_quant/2; n_clip_v = 0 or logged; C and Python hashes agree on 10⁷ IDs |
| V6 | momentum spectrum from the subsample vs `mom_*`, k < 0.3 | within subsample shot noise |
| V7 | determinism of the regenerate-only mode, twice and with another thread count | identical md5, or thread dependence documented and enforced. **Passed 2026-10-17**: `deltaL` md5 identical with 1, 2 and 4 MPI ranks (`run_deltaL_test.sh`) |
| V8 | f(R) vs ΛCDM, same seed: `deltaL` | identical md5. Holds by construction (2026-10-17): the FML generator has no gravity-model input; twins share one IC file |
| V9 | overhead at production size (one short job, approval first) | extra wall time < 5 %; peak RSS within node memory with margin. **Measured on the node-2 pilots (1024³, N_g = 2048; section 16, 2026-10-17)**: the extras cost at most 1.4 min per primary output (≤ 6 % of the f(R) run, ≤ 15 % of a ΛCDM run; an upper bound until `ExtraOutputs` has per-product timers); MaxRSS 118 GB (ΛCDM) and 177 GB (f(R)) on a 1 TB node. Production size (N_g = 4096) not yet timed |

## 9. Storage budget (L = 1024 Mpc/h)

Mode count (2n+1)²(n+1), n = floor(k_c L/2π), 8 B per mode.

| product | per output, 2048³ | per output, 1024³ |
|---|---|---|
| lens planes (12 × 4096²) | 0.8 GB | 0.8 GB (or 0.2 GB at 2048²) |
| `CatshortV` catalogue, text | 0.3–0.4 GB apparent, ≈ 0.5× on disk (this Lustre compresses the text 1.9×) | 0.22 GB apparent, 0.11 GB on disk (z = 0: 777,651 haloes × 283 B) |
| `CatshortV` catalogue, binary (`.bin`, 96 B per halo) | ≈ 0.12 GB apparent, ≈ 0.11 GB on disk (float32 compresses 1.09×) | 0.075 GB apparent, ≈ 0.07 GB on disk |
| `deltam`, k_c = 1.5 | 0.47 GB (3.7 GB at k_c = 3.0) | same |
| `mom_x,y,z`, k_c = 1.0 | 0.41 GB | same |
| subsample, M = 100 | 1.03 GB | 0.13 GB |
| f_R at haloes | ~8 B per halo | same |
| `deltaL` (once per run) | 0.47 GB | same |
| reference: float32 snapshot | 206 GB | 26 GB |

Per box of 42 outputs at 2048³: lens planes 34 GB, catalogues ~15 GB
apparent as text (≈ 7.5 GB on disk; ≈ 5 GB as binary; `both` adds ≈ 60 %
to the on-disk catalogue space),
extras at all 42 outputs ~80 GB (≈130 GB per box, 13 TB for 100 boxes).
**Decision (CZ, 2026-10-17): the extras are written at all 42 outputs**
(every output is class 1; BDM at all 42 as well unless Baojiu asks for
primary-only catalogues, which the class list can still express). For
the 96 ECOSMOG companions at 1024³ the extras cost ~1.2 GB per output
(~50 GB per run). Baojiu's Slack budget covered lens planes only
(3.4 TB). Optional offline `zstd` gains are modest.

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
6. ~~Extras cadence~~ decided 2026-10-17: all 42 outputs (section 9).
7. k_c,δ = 1.5 or 3.0 (do tidal/environment definitions need R ≈ 1 Mpc/h)?
8. I7 route for the FML generator and for Baojiu's ICs (section 5, I7).
9. File naming: decided for GLAM's `<product>.<step>.<box>.{bin,meta}` in
   the run directory (as implemented; `glamio.product_path` builds the names).
10. Companion resolution: N_g = 2048 or 4096 on Baojiu's 1024³ particles.
11. BDM catalogue velocities: keep the legacy 1/a conversion (text and
    binary agree; the binary `.meta` records `a`, `a_v`, `da` and
    `velocity_scale`, glamio's `velocities_at_av()` applies a/a_v) or
    switch the writer to 1/a_v. A switch changes every GLAM catalogue, so
    it is a separate commit after this decision.
12. Production catalogue files: `binary` (the compiled default since
    2026-10-18, user decision); `both` only for a run whose catalogues
    halocat must ingest before it reads the binary. On disk the binary saves ≈ 40 %
    of the catalogue space (≈ 1.3 GB per 2048³ box with 17 catalogues,
    ≈ 0.25 TB over the suite), not the 3× of the apparent sizes, because
    this Lustre compresses the text 1.9× and float32 only 1.09×.

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
- [x] I2 lens planes + L1, L2 (`degrace-track-a` `b774f8d`); L3 (plane P(k) vs 3-D P(k)) in the suite since 2026-10-17
- [x] Pilot node 2 vs the ECOSMOG twins: both pilots complete and verified (ΛCDM job 12180082, f(R) job 12180110, 2026-10-17, section 16); the ≈ +0.5 % linear-scale offset is the same in both gravity models (code/ingestion; origin open, da0 test)
- [x] I7 δ_L dump, route 1 in the FML generator + V3, V7, V8 (`degrace-track-a`, 2026-10-17; `test/extras/run_deltaL_test.sh`); legacy pilots out of scope (decision)
- [x] I8 + V1, V2 (`degrace-track-a` `a41ae37`); V3 passed with I7; V4 needs the pilot catalogues
- [x] I9 + V5, V6 (`degrace-track-a`, commit after `09c4aa4`); ic2pm places particles at index = Gadget ID and tags the lattice order since the 2026-10-17 review fixes
- [x] Review of the branch (six subagents, `docs/reviews/degrace-track-a/`, 040f011) and fixes of the eight Majors (commit after 040f011)
- [x] I10 (`degrace-track-a`, commit after `a41ae37`); unit test + f(R) smoke case
- [ ] V9 overhead run: ΛCDM baseline measured (planes < 1 %, P(k) at 42 outputs 11 %, BDM at 17 primaries 23 % of the wall time); f(R) pilot: the extras (modes, subsample, f_R at haloes) ≤ 1.4 min per primary = ≤ 6 % of the run (upper bound), MaxRSS 177 GB; still open: per-product timers in `ExtraOutputs`, and a production-size (N_g = 4096) timing
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
- 2026-10-17, I7 route 1 on `degrace-track-a`: δ_L dump in the FML
  generator (`DeltaLDump.h`, Lua `write_deltaL`/`deltaL_kc`/
  `deltaL_output_prefix`; tracked in `2LPTIC_Gui/fml_intel_port/`, the
  build directory `2LPTIC_Gui/FML/` is gitignored). The port copy of
  `Main.cpp` was 7 weeks behind the build directory (an ID index rounding
  fix of 2026-08-27, `llround` instead of `floor`, with a
  `test-id-round-v1` Makefile target); both are now synced and the
  provenance is in `APPLY.md`. Tests: `test/extras/run_deltaL_test.sh`
  (V7, layout A/B on 1LPT and 2LPT ICs, V3 end to end through `ic2pm` and
  PMP2main) all pass; numbers in I7. Lesson: a test that divides mode by
  mode is ill-conditioned where δ_L is small and the CIC of particles on
  grid nodes is non-analytic in Ψ, so the assertions use residual rms,
  binned means, cross-correlation coefficients and the mean phase per axis
  (a coordinate-offset estimator) instead of per-mode ratios.
- 2026-10-17, L3 and run verification: `case10_lens72` now also writes
  `PowerDM.R`; `asserts.py` compares the plane power with the exact
  projection formula (177/177 checks). `DEGRACE-companion/check_run.py`
  verifies a run directory (schedule vs `Init.dat` and the pinned ECOSMOG
  list, 42 PowerDM + LensPlanes with sizes/epochs/means, catalogues at the
  17 primaries, ic2pm epoch, `timing.log` output overhead, sacct). Dry run
  on the running ΛCDM pilot after 140 steps: non-output step 0.61 min;
  station outputs 0.99 min (Analysis 0.43 = 2048³ P(k) + planes 1.8 s,
  IO 0.21 for the 201 MB plane file); primary outputs 2.53 min (BDM, logged
  in the Analysis column); output steps so far 25 % of the wall time, of
  which the lens planes themselves are < 1 %.
- 2026-10-17, ΛCDM pilot (node 2, job 12180082) **completed and verified**
  (`check_run.py`: 0 FAIL): 42 exact outputs on the pinned ECOSMOG epochs
  (max rel diff 1e-7), 42 `PowerDM.R`, 42 `LensPlanes` (201 MB each, axis
  means 1 to 1e-7, Σ/Σ̄ from 0.33 to 1.8 at z = 8.7 and from 0.11 to 50 at
  z = 0), 17 `CatshortV` at the primaries, 9.1 GB per box in all.
  Resources: 2 h 32 min wall on 64 cores of one `cosma8-serial` node,
  CPU efficiency 83 %, MaxRSS 118 GB (sacct; 91 GB for `PMP2main` alone),
  181 steps (163 ladder + 18 inserted), median step 0.60 min; output steps
  cost 43 min = 29 % of the wall time (primaries 2.6 min each, almost all
  BDM; stations 0.96 min: 2048³ P(k) 25 s, planes 1.8 s, 201 MB write
  12 s). The lens planes are < 1 % of the run. Production estimate per
  ΛCDM box at this size: ~165 core-hours (plus ~10 % for the extras, to
  be measured).
  Twin comparison (`compare_ecosmog_pk.py`; ECOSMOG has only f(R)
  spectra, so this is ΛCDM vs the f(R) twin from the same IC files):
  P_GLAM/P_ECOSMOG at k = 0.01–0.03 is +0.47 % at z = 2 rising to +0.69 %
  at z = 0 and is already +0.70 % at z = 8.7 where f(R) is negligible at
  all k; at 0.03–0.1 the ratio goes 1.006 → 0.995, at 0.1–0.3 1.006 →
  0.953 and at 0.3–1 1.05 → 0.897 (z = 2 → 0): the f(R) enhancement for
  |f_R0| = 5.6e-6 with the expected scale and time dependence. The
  constant ~+0.5 % at linear scales (mode-exact, same IC) is unexplained:
  candidates are the half-step velocity convention of `ic2pm` vs
  ECOSMOG's start, POWMES vs raw-CIC normalisation, or a 0.25 % growth
  difference between the codes. The f(R) pilot against the f(R) twin
  decides whether it is a code/ingestion offset (ratio stays 1.005) or
  not; a da0 = 3e-4 restart of the first 20 steps would isolate the
  half-step term (expected to move linear P by ≈ 0.3 % per the
  `2LPTIC_Gui/README.md` numbers).
  At k = 0.3–1 and z ≳ 3 our raw spectra exceed POWMES by ~10 % where
  P ~ 30 (Mpc/h)³: discreteness of the perturbed lattice (white shot
  noise alone is 1 (Mpc/h)³), not a pipeline error; compare below
  0.3 h/Mpc at high z.
- 2026-10-17, review fixes (eight Majors of `docs/reviews/degrace-track-a/TRIAGE.md`):
  - R3-M1: `ic2pm` now reads the Gadget ID block (4- or 8-byte, stream
    access), places every particle at array index = ID, verifies that the
    IDs are a permutation of 1..N (sentinel fill), detects the lattice
    order from the first file (Baojiu's 2LPTic files: IDs sequential in
    file order, z fastest, verified on `Node_002/ics.0`; FML: x fastest)
    and appends `IDx`/`IDz`/`ID?` to the PM header, which every product's
    `.meta` carries; the subsample `.meta` explains the tags. For 2LPTic
    files the particle order is unchanged (sequential IDs); multi-rank
    FML files are re-ordered. `PMcrd.DAT` differs from older conversions
    by the header tag.
  - R2-M1: the optional `Init.dat` parsers read the first token after the
    last `=` (any comment after a space is fine, `kc_delta = 1` is 1.0)
    and stop with the offending line when the token is not a number;
    EOF or a line without `=` still means "absent". `ReadSetup` stops on
    a malformed extras line and defaults only at end of file.
  - R2-M2: `fR at halos` needs `MG_flag = 1` and `MG_model = 1` (PMP2init
    demotes with a message; the BDM hook checks both); `MG_model` is
    written into the LensPlanes, modes and Subsample sidecars.
  - R5-M1: `compare_ecosmog_pk.py` takes the epoch from
    `OutputSchedule.dat` by step (header fallback with the f8.5 rounding
    tolerance).
  - R5-M2: case 11 compares the decoded full subsample's CIC cubes with
    the `deltam` and `mom_x` dumps mode by mode (rms residual < 1e-2 of
    the rms mode) and requires the transposed and conjugated cubes to
    fail; `kc_mom` of the test cases is 0.185 (below k_N(36)/4 = 0.1885).
  - R6-M1/M2/M3: section 7 rewritten with the implemented keys and
    defaults; section 3 states the BDM catalogue epoch offset; I8 item 3
    records the interlacing status and the new `kc_mom` rule/default
    (0.75; `gen_inputs.py` default likewise).
  - Tests after the fixes: suite 180/180 (new case4f: a malformed
    optional line stops PMP2init; mode-by-mode dump test), modes and
    subsample unit tests, `run_deltaL_test.sh` (its `kc_mom` lowered to
    0.3 for the new rule) and the new `test/extras/run_ic2pm_test.sh`
    (synthetic two-file Gadget-1 ICs: int32 sequential z-fastest IDs →
    `IDz`, int64 shuffled x-fastest → `IDx`, duplicate ID → refused;
    verified through a one-step GLAM run with `Subsample M = 1`). The
    `compare_ecosmog_pk.py` rerun matches all 42 ΛCDM-pilot epochs.
- 2026-10-17, f(R) pilot (node 2, job 12180110, c353887 binaries, extras
  on: `kc_delta = 1.5`, `kc_mom = 1.0`, `N_dep_mom = 0` → 1024,
  `Subsample M = 100`, `fR at halos = 1`) **completed and verified**
  (`check_run.py`: 0 FAIL, 0 WARN; the outputs are saved beside the run
  as `check_run_12180110.txt` and `compare_ecosmog_pk_12180110.txt`):
  42 exact outputs on the pinned epochs, 42 `PowerDM.R`, 42 `LensPlanes`,
  17 `CatshortV`, and at the 17 primaries `deltam` (n = 244, 469 MB),
  `mom_x/y/z` (n = 162, 138 MB each), `Subsample` (10.74 M particles =
  N/100 to 0.03 %, n_clip_v = 0, 129 MB) and `fR_halo` (one pair per
  catalogue row, 1.12 M rows at z = 1, 9 MB); 25 GB per f(R) box
  (ΛCDM 9.1 GB: the extras add 1.0 GB per primary).
  Resources: 6 h 08 min wall on 64 cores (`cosma8-serial`, 400 GB
  requested), CPU efficiency 83 %, MaxRSS 177 GB (ΛCDM 118 GB), 181
  steps, median non-output step 1.72 min of which the multigrid solver
  1.1 min; output steps 58.7 min = 16 % of the wall time. Station outputs
  2.09 min (+0.37 min: P(k) and planes, as in ΛCDM); primary outputs
  4.62 min (+2.9 min: BDM, modes, subsample and f_R at haloes, all logged
  in the Analysis column). Against the ΛCDM primaries (+2.0 min, BDM
  only) the extras cost at most 1.4 min per primary, ≤ 23 min = 6 % of
  this run and ≤ 15 % of a ΛCDM run (V9; an upper bound, because the BDM
  cost itself differs between the twins; per-product timers in
  `ExtraOutputs` would settle the split). ≈ 390 core-hours per f(R) box
  at this size.
  Twin comparison, f(R) vs the f(R) POWMES spectra (the real twin test):
  P_GLAM/P_ECOSMOG at k = 0.01–0.03 is +0.70 % at z = 8.7, +0.48 % at
  z = 2 and +0.82 % at z = 0; at 0.03–0.1 +0.96 % → +0.65 % → +0.75 %;
  at 0.1–0.3 +1.35 % → +0.97 % → +0.39 %; at 0.3–1 +10 % (discreteness)
  → +7 % → −0.6 %. The ≈ +0.5 % at linear scales is therefore the same
  in both gravity models (ΛCDM-vs-f(R) gave +0.47 % → +0.69 % at the
  same bins): a code or ingestion offset common to both runs, not a
  gravity-model difference; the differences between the two pilots'
  low-k ratios (≤ 0.15 %) are at the level of the f(R) enhancement and
  the 3-modes-per-bin scatter. Within ±1 % up to k = 0.3 at every epoch;
  at z = 0 our f(R) spectrum is 0.6 % below ECOSMOG's at k = 0.3–1 where
  the ΛCDM one sits 10 % below, i.e. the f(R) enhancement is reproduced
  to ≈ 0.6 % of P there (PM force resolution at dx = 0.5 Mpc/h; the
  deficit grows to ≈ 9 % at k = 1.3). The da0 = 3e-4 restart test of 20
  steps remains the way to isolate the half-step term of the offset; a
  POWMES normalisation check against a raw spectrum of one ECOSMOG
  snapshot is the other candidate test.
- 2026-10-17, repository and data: `degrace-track-a` merged into `cz`
  (fast-forward, 17 commits), status deck in
  `docs/slides/degrace_companion_status/`, finished runs and notes
  consolidated into `archive/` (README with restore commands). At the
  user's request the 2048³ validation IC (`2LPTIC_Gui/validate_highres/
  snap`, 275 GB, regenerable in 8 min) and the step-100 checkpoint of the
  finished L1800 ΛCDM Run1 (786 GB) were deleted; `archive/README.md`
  records both.
- 2026-10-17, binary halo catalogue (branch `bdm-binary-catalogue` from
  `cz`): `BDM.config` gained `CatalogueFormat = text | binary | both`
  (default `both`). The finder writes, beside or instead of the text
  `CatshortV.<step>.<box>.DAT`, a headerless little-endian `.bin` of
  96-byte records (the 24 text columns as 22 float32 + 2 int32, same rows
  and order) and a `key = value` `.meta` sidecar (columns, dtypes, units,
  rows, a, z, a_v, da, velocity_scale, step, Box, grids, cosmology,
  particle mass, MassMin, overdensity, iVirial, buffer width, finder,
  header). Text and binary rows come from the same variables, so they
  agree to the text's printed precision; the binary keeps the full real*4
  values. (Corrected 2026-10-18: the compiled default is `binary`, see the
  entries below; the
  size figures below are apparent bytes; on disk the saving is ≈ 40 %, see
  section 9.) The staged-file publication was generalised to the three files
  (renamed text → bin → meta). HDF5 was rejected for the Fortran side: the
  cosma modules are gnu-built, so ifx would need its own HDF5 build, and
  the sidecar already covers self-description; HDF5 stays on the Python
  side (halocat). Measured on the ΛCDM pilot at z = 0: 220 MB text → 75
  MB binary per catalogue (1.9 → 0.65 GB per box at 1024³; ≈ 6 → 2 GB per
  box at 2048³; ≈ 0.75 TB over the suite); the write takes 7 of 88 s of
  the finder at z = 0. Readers: new package `python/glamio` (numpy only,
  editable install into cosemu; `read_catalogue` for text and binary,
  `velocities_at_av`), its own unit tests, and the precise_z suite now
  checks binary = text column by column on the halo-finder case, a
  binary-only case and a bad-value rejection.
- 2026-10-18, review of `bdm-binary-catalogue` before the merge: five
  subagent reviews (`docs/reviews/bdm-binary-catalogue/`, prompts, five
  findings files, `TRIAGE.md`). One Blocker (the default `BDM.config` the
  finder writes when none exists truncated the new key to `Catalogue` via a
  `T20` tab stop, so every run without a `BDM.config` and with ≥ 2 catalogue
  outputs died at the second one; found by two reviewers, one by running
  it) and six Majors, all fixed the same day: format fixed and a
  default-config test added (case14); `CheckBdmConfig` validates
  `BDM.config` at start-up; the compiled default was set to `text` in this
  round (`gen_inputs.py` writes the key only when asked); the sidecar's
  `definitions` line now states Rvir = R_SO + resolution correction and
  Mbound = bound mass inside R_SO; `check_run.py` verifies `.bin`/`.meta`
  and reports stale staged temporaries; the storage figures distinguish
  apparent from on-disk bytes; the velocity-epoch prose had the direction
  inverted (catalogue velocities are low by a_v/a; `glamio` multiplies by
  a/a_v, which was already right). glamio gained the row self-check,
  required layout lines, `mmap`, unified meta keys and two tests; the suite
  gained case15 (text mode equals case13's binary rows) and catches
  glamio's own exceptions as failures; 215/215 checks pass.
- 2026-10-18, compiled default `CatalogueFormat = binary` (user decision
  before the merge into `cz`): a run without `BDM.config`, or with one that
  lacks the key, writes `.bin` + `.meta` only. The review round had set the
  default to `text` to spare the text consumers; the user chose the
  binary-only default instead (≈ 40 % less catalogue space on disk, no
  duplicate files). Consequences: halocat and any script that opens
  `CatshortV.*.DAT` from a run made with a build from this commit on needs
  the glamio reader or `CatalogueFormat = text|both` in that run's
  `BDM.config`; the frozen pilot binaries and every existing catalogue are
  unaffected. The fR_halo sidecar now says for each `catalogue` line when
  the file exists. case14 checks the binary-only default (readable, rows
  > 0, `text_catalogue = none`); the suite passes 215/215 (36 s, 4 threads).
- 2026-10-18, glamio readers for the other products (branch `glamio-readers`
  from `cz`, python/glamio 0.2.0): `read_schedule` (OutputSchedule.dat),
  `read_lensplanes`, `read_modes` (deltam, mom_x/y/z and the companion
  `*_deltaL` cubes; `mode()` with Hermitian completion, `power()` = the
  shell average that equals PowerDM.R bin by bin), `read_subsample`
  (decode, `ids()` from the splitmix64 rule), `read_frhalo` (with
  `catalogue()`), `read_powerdm` (R/log/Z tables, epoch line) and the PM
  snapshot (`read_pm_header`: the 529-byte big-endian record, byte order
  detected from the marker; `Snapshot.iter_pages()` streams the 24 MiB
  `PMcrs` pages, `read()` loads all; `units='physical'` gives Mpc/h and
  km/s at a_v = AEXPN - ASTEP/2). A shared `MetaProduct` base now holds the
  sidecar scalars (`Catalogue` inherits it). Validation: the precise_z
  suite gained six `glamio` checks that compare every reader with the
  suite's own decoders on real output (case10/11/12), including
  snapshot particles = subsample particles within the quanta (which also
  ties the header's a_v to the subsample's); 221/221. Unit tests with
  synthetic files (`python/tests/test_readers.py`, 8 tests incl. the
  SplitMix64 reference vector and a little-endian snapshot): 15/15 with
  the catalogue tests. The pilot's 2048^2 lens planes read memory-mapped
  with axis means 1. Not done here: the halocat ingest (halocat session),
  `ProjectionXY/XZ/YZ` (legacy product, not in the DEGRACE output list),
  `timing.log`, `main.log` and the BDM `outputB` logs (text, parsed ad hoc by
  `check_run.py`), `PowerDM.Z` (never written by the DEGRACE runs).
  `DEGRACE-companion/check_run.py` keeps its own decoders on purpose, as an
  independent verification of the same files.
- 2026-10-18, review of `glamio-readers` before the merge: five subagent
  reviews (`docs/reviews/glamio-readers/`: prompts, five findings files,
  `TRIAGE.md`). No Blocker; seven distinct Majors, all fixed the same day:
  `Modes.power()` now works slab by slab (was 9x the cube in temporaries,
  4.7 GB on the production cube) and refuses `nbins > n` (truncated shells
  were returned 5-19 % low); `load_array` returns an empty array for a
  0-row product under `mmap=True` (was a numpy crash on a legitimate
  `rows = 0` fR_halo); `Snapshot.read()` fills preallocated arrays and
  refuses more than 2^27 particles by default (was two copies and no cap:
  412 GB for 2048^3 in physical units); `Subsample.ids()` checks the count
  first and caps the result (was 69 GB of `arange` at M = 1 for 2048^3);
  the PMcrd test fixture is packed with `struct` instead of the reader's
  own dtype and the suite pins Nrow, realization, ASTEP, Om, OmL, h and
  PARTW against the subsample sidecar (field-order errors were invisible);
  the per-plane mean claim was corrected here, in the deck and in the
  `lensplanes` docstring (only the per-axis mean is 1). Minors fixed:
  `sync`-mode IC files get `a_v = a`, the `PMcrs` split (512 or 256
  records) is detected from the file size, `particle_mass()` uses the
  finder's 2.774e11, the caption of `PowerDM` is searched after the epoch
  line, `read_powerdm(box=)` gives an exact `kNyq`, the CIC window helper
  is named for what it is (Jing's alias-summed form), `read_modes`
  accepts only the mode products by default, `M >= 1` and
  `prefer in (binary, text)` are validated, the schedule rejects bad
  `exact`/`class` values and gained `row_at_z()`, `product_path()` /
  `products_at()` build and list the files of a step, the README example
  reads lens planes at every output. glamio 0.2.1; suite 221/221 with
  the strengthened checks; unit tests 16/16.
- 2026-10-18, confirmation pilot of the current `cz` build submitted:
  `DEGRACE-companion/runs/Node_002_LCDM_L1024Np1024Ng2048_confirm20261018`,
  Slurm job 12185817 (cosma8-serial, dp203, 64 threads, 220 GB, 8 h limit).
  A first submission (12185806) stopped in PMP2init after 2 s: it had
  copied the f(R) pilot's `kc_mom = 1.0`, which the rule added in the
  review fixes rejects at N_dep = 1024 (k_N/4 = 0.785 h/Mpc; the f(R)
  pilot's build predates the rule); the job script also missed the
  failure because PMP2init's `stop` returns exit code 0 and ran ic2pm
  into the missing `Setup.dat`. Fixed: `kc_mom = 0.75` (the production
  value), the script now requires `Setup.dat`, and PMP2init was run on
  the login node before resubmitting. Follow-up for the Fortran:
  configuration refusals in PMP2init should use `error stop` so that a
  job script sees a non-zero exit code (done later the same day, see the
  entry below).
  Executables built from `cz` 88fa70a (Fortran sources unchanged since
  69bd63b; md5s in `bin/md5sum.txt`). Same ICs, schedule, classes
  (primary), planes (4 x 2048^2) and BDM as the 17 Oct pilot (job
  12180082, build b774f8d), plus `CatalogueFormat = both` and the extras
  at the 17 primaries (deltam k_c = 1.5, momentum k_c = 0.75 on N_dep =
  1024, subsample M = 100). Purpose: the production recipe with the
  merged build. Expected: P(k), lens planes and the text catalogue rows
  identical to the 17 Oct pilot (the dynamics are unchanged), the binary
  catalogue equal to the text column by column, every product readable
  with glamio, `check_run.py` 0 FAIL, ~3 h wall.
  **Result (2026-10-18): job 12185817 COMPLETED, exit 0, 2 h 39 m on
  m8407 (PMP2main 2 h 37 m), CPU efficiency 81 %.** `check_run.py` 0 FAIL,
  0 WARN; 42 `PowerDM.R`, 42 lens-plane pairs, 17 catalogues in both
  formats, 17 x 4 mode cubes, 17 subsamples, no staged temporaries.
  glamio 0.2.1 read every product: the 17 binary catalogues equal the
  text column by column at the text precision (13,656,386 rows; caveat:
  glamio's own text reader stores float32 columns, which hold a position
  in a 1024 box only to one ulp = 6.1e-5, so compare the text through a
  float64 parse or use the binary); 68 cubes finite with a consistent
  Hermitian accessor (`deltam` n = 244, k_max = 1.497 h/Mpc; `mom_*`
  n = 122, k_max = 0.749 on N_dep = 1024); the z = 0 `deltam` cube's
  `power()` equals the `PowerDM.R` table to four digits at k = 0.02-1.0;
  the z = 0 subsample keeps 10,741,150 of 1024^3 particles (1 in 99.97
  at M = 100; lattice tag IDz; unique IDs < 1024^3; positions in [0, L);
  1-D velocity rms 321 km/s); planes with unit axis means and
  `MG_model = 1` in the sidecar.
  **Against the 17 Oct pilot the dynamics agree to float32 round-off,
  not bitwise.** Planes identical through step 98; one pixel differs by
  one ulp at step 102; 10 % of the pixels differ at step 105 (max 6e-6
  relative) and 98 % at z = 0 (rms 4e-6, max 1e-4 in the densest
  pixels; slab means equal to 1e-11). 58 of ~43,000 P(k) bins differ,
  all in the last printed digit (max |P_new/P_old - 1| = 9e-6). The
  catalogues have +1, +1 and +12 rows at z = 2, 0.625 and 0 (of
  507,388 / 958,104 / 777,651) with the counts above 1e13 and 1e14
  M_sun/h equal within 1-3; 48 % of the z = 2 rows and 1.3 % of the
  z = 0 rows print identically. The cause is not the code: the ICs
  store the particles in Gadget-ID order, so the new `ic2pm` (array
  index = ID) builds the same particle array as the old file-order
  placement; `OutputSchedule.dat` is identical; and every new line in
  the time-step path is guarded by `OutClass == 1`, first true at step
  105, after the first visible difference. GLAM's CIC deposition
  (`DENSIT`, `PMP2mod_density.f90` 235-280) adds each particle's eight
  weights into the shared float32 mesh with `!$OMP ATOMIC`, so the order
  of the adds into a cell follows thread timing: at 1024^3 with 64
  threads the code is not bitwise reproducible run to run, and the
  non-linear evolution amplifies the ulp differences until they show
  (the suite's 36^3 case with 4 threads does reproduce bitwise: two
  fresh runs and the suite's own output are identical). Runs of the same
  configuration therefore compare at round-off tolerance (P(k) 1e-5,
  planes 1e-4 point-wise, halo counts 1e-5), never byte by byte.
  **Cost of the new outputs at 1024^3:** `timing.log` 157.0 vs 149.9 min
  (Analysis 57.7 vs 45.8 min = +0.70 min per primary output for the
  cubes, subsample and binary catalogue; station outputs unchanged at
  0.95 min; non-output steps 0.56 vs 0.60 min, node scatter). Process
  peak RSS (`/usr/bin/time`) 86.65 vs 86.60 GiB: +54 MB, because the
  extras' buffers (`G(N_dep^3)` = 4 GiB, the `C` cubes, the subsample
  chunks) are allocated and freed inside Analysis and never coexist
  with the finder's mesh, so the peak is unchanged. Slurm MaxRSS 121.9
  vs 112.7 GiB: the difference is the page cache of the written
  products, which the cgroup counts. Scaling the process footprint by 8
  gives the full-size need: LCDM with BDM ~0.69 TB, f(R) (123.0 GiB in
  the f(R) pilot) ~0.98 TB, equal to the 985 GB of job 11448712 in Q8;
  the earlier 1.4 TB figure scaled sacct and overstated it. Products
  22 GiB per box at this resolution: planes 7.3, `deltam` 6.9, `mom_*`
  2.6, subsamples 1.9, catalogues 1.2 binary + 2.3 text, P(k) 2 MB.
  Run directory kept (checkpoints deleted by the script); outputs of
  the checks in the session scratchpad only.
- 2026-10-18, PMP2init refusals exit with status 1. Every `stop` in
  `PMP2init.f90` (22: the Init.dat validation rules, the P(k) table
  checks, the missing-file cases and the "template written" branch) is
  now `call Refuse(message)`, a new subroutine that prints
  `PMP2init refused: <message>` on stdout and stderr and executes
  `error stop 1`. Measured with ifx 2024.2: a plain `stop 'msg'` returns
  status 0 with the message on stderr; `error stop 'msg'` returns 128;
  `error stop 1` returns 1, which is what job scripts and `set -e`
  expect. Checked with the confirmation run's own inputs: `kc_mom = 1.0`
  gives status 1 and the interlacing message, `kc_mom = 0.75` gives
  status 0 and a `Setup.dat` identical to the run's. The suite's
  `run_init` records the status in `init.rc` and a new check requires 1
  for the four refused cases (4a, 4b, 4e, 4f) and 0 for the sixteen
  accepted ones: 222/222. `submit.sh` keeps its `Setup.dat` guard. Not
  changed (same class, PMP2main's start-up): the `Setup.dat` parse
  refusals in `PMP2mod_tools.f90` (five plain `stop`s) and
  `CheckBdmConfig`, whose job-script symptom is a missing first
  catalogue rather than a silent continuation.
- 2026-10-18, campaign clean-up (file count). Removed, all regenerable:
  the Python caches (`DEGRACE-companion/__pycache__`, `python/.pytest_cache`,
  `python/glamio/__pycache__`, `python/tests/__pycache__`,
  `test/precise_z/__pycache__`, the status deck's `figs/__pycache__`),
  `python/glamio.egg-info` (a by-product of `pip install -e`; the
  editable install lives in cosemu's `site-packages` as a PEP 660 finder
  and still imports), the suite's `test/precise_z/work` and the session
  scratch. Saved beside the confirmation run, following the f(R) pilot's
  convention: `check_run_12185817.txt` and `verify_vs_pilot_12185817.txt`
  (the plane / P(k) / catalogue comparison with job 12180082 and the
  glamio reads). Nothing scientific was removed. Candidates left for a
  decision: the 17 Oct ΛCDM pilot run (226 files, 9.1 GB; superseded at
  round-off level by the confirmation run, which holds every product it
  has and more; archive as one tarball in `archive/` or delete), the
  untracked `DEGRACE-pilot/` (July supplement launch scripts, 17 files)
  and `verify/` (April OpenMP race logs, 6 files), the three BDM review
  decks of August/September under `docs/slides/`, and the 38 `.o`/`.mod`
  build intermediates in the repository root.
  Removed on the user's instruction later the same day: `DEGRACE-pilot/`
  (17 untracked files, 764 KB, July supplement launch scripts and a
  background test; not archived, gone), `verify/` (6 untracked files, 172 KB,
  April OpenMP race logs; gone), the three BDM decks
  (`bdm_duplicate_survival`, `bdm_session3_findings`,
  `bdm_session4_validated_algorithm`, 24 files; the three decks (24 files) are still in the pushed history before the
  untrack commit: `git archive 9db3aac docs/slides/bdm_session3_findings | tar -x`
  restores one), and the 22 `.o` + 16 `.mod` build
  intermediates in the repository root (`make` rebuilds them; the five
  `.exe` files stay). Left: the ΛCDM pilot run (decision pending),
  `DEGRACE-highres/` and `resume.md`.
- 2026-10-18, first look at the products as data (local review deck,
  frames "In the data (1)–(6)"; scripts `measure_products.py` and
  `plot_products.py` beside the deck, glamio reads of both node-2 runs,
  about one minute each on a login node). Findings, all from the
  on-the-fly outputs alone: `PowerDM.R` follows linear theory to ±2 % for
  0.02 < k < 0.1 at every primary; the f(R)/ΛCDM ratio reaches +23 % at
  k = 5 at z = 0 and +4 % at z = 2; the `deltam` cube's own P(k) equals
  the table to 4e-6 shell by shell; θ = i k·p from the `mom_*` cubes
  gives P_θθ = (aHf)² P_δδ to 2 % for k < 0.06 at z = 0, 1, 2 with no
  free parameter, and the density-weighted term dominates beyond k ≈ 0.1
  (15× at k = 0.7, z = 0); the lens-plane rms on 8 Mpc/h pixels follows
  D(a) from z = 8.7 to z ≈ 1 and exceeds it by 7 % at z = 0 (on the native
  0.5 Mpc/h pixels particle discreteness adds 0.125 in quadrature); the
  subsample's 1-D velocity rms is 4/8/15 % above the linear prediction at
  z = 2/1/0 with exponential tails; the f(R) mass function has 25 % more
  haloes above 3e12 at z = 1 and above 1e14 at z = 0; `fR_halo` shows the
  chameleon transition near 3e13 at z = 0 and full screening above 2e12
  at z = 2. The inverse transform of the cube, the lens plane and the
  subsample at z = 0 show the same web. The scripts live in
  `python/examples/` (with `README.md`; `--lcdm/--fr/--out/--data`
  options, defaults = the node-2 runs); the deck's `figs/` holds only
  their outputs.
- 2026-10-20, halo environment and tides from the products
  (`python/examples/halo_environment.py`, glamio reads of the `deltam`
  cube and the catalogue of one output): δ_R, the tidal eigenvalues,
  the anisotropy α_R and the T-web class per halo, cubic-spline read-off
  with the B-spline prefilter applied in Fourier space (exact for periodic
  grids; scipy's recursive prefilter took 60 s per 512³ field). Feasibility
  recorded in the deck ("Use case (1)/(2)" frames): the cube supports
  Gaussian R ≥ 2 Mpc/h (1 % kernel leakage past k_c = 1.5), the 1-in-100
  subsample is too sparse (shot noise σ = 0.44 at R = 5), halo-scale
  tides (R_G < 1) and formation times are out of reach, the catalogue has
  C_vir, λ, b/a, c/a, X_off. Validation at z = 0 of the node-2 ΛCDM run:
  direct Fourier sums at six haloes agree with the 512³ values to 2 %
  (R = 2) and 0.02 % (R = 5), 768³ brings R = 2 to 0.25 %; eigenvalue sum
  = δ_R to 2e-5; ⟨δ_5⟩ rises from 0.29 to 7.7 and ⟨α_5⟩ falls from 0.55 to
  0.15 between 10^12.2 and 10^15.2 Msun/h. 40 s and 5 GB per R at 512³
  (100 s, 15 GB at 768³) on a login core; the f(R) pilot runs through the
  same path (856,309 haloes). README section in `python/examples/`.

