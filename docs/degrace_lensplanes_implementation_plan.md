# DEGRACE L1024 production: snapshot schedule and lens-plane outputs

Implementation plan for MG-GLAM, derived from the Slack thread with Baojiu Li in
`slack-chat-degrace-prep.txt` (analysed 2026-10-16). Written for agent sessions
and for the user. Status: **nothing implemented yet**; the numerical cross-check
Baojiu asked for is done (section 2).

Review deck: `docs/slides/degrace_lensplanes_plan/degrace_lensplanes_plan.pdf`
(source `.tex` beside it; the schedule figure is made by
`docs/slides/degrace_lensplanes_plan/figs/make_output_schedule_fig.py`).

Scope: the **L1024 Np2048 Ng4096** production (100 LCDM boxes in the fiducial
Planck 2015 cosmology, later the design cosmologies). It is *not* the L1800
high-res suite: 256 Mpc/h does not divide 1800.

## 1. Agreed output design (from the chat)

| Item | Decision |
|---|---|
| Particle snapshots | **None** (`iSave=0`). 200 GB per snapshot, 9 TB per box is not feasible for 100 boxes. Checkpoints (`PMcrd/PMcrs`) are still written for restart and deleted at the end. |
| Velocities | Kept only in the BDM halo catalogues (`CatshortV`). |
| Lens planes | At every output: project the density along each of the 3 axes in **4 slabs of 256 Mpc/h** depth → **12 planes of 4096² float32** = 0.8 GB per snapshot, 34 GB per box, 3.4 TB for 100 boxes. |
| Outputs per box | **42** = 17 *primary* redshifts common to all cosmologies (z = 0, 0.125, 0.25, …, 2.0) + 25 *secondary* at comoving distance χ = 256 k Mpc/h, k = 1…25 (to 6400 Mpc/h). FORGE/DEGRACE-ECOSMOG convention. |
| Secondary set | Cosmology-dependent: recompute a(χ = 256 k) per cosmology for design runs. |
| Plane resolution | 4096² (= Ngrid). Baojiu first floated 2048²; keep a downsample option. |

## 2. Verified numbers

**Cross-check (done 2026-10-16).** Baojiu's 42 scale factors reproduce flat ΛCDM
with Ω_m = 0.3089, Ω_Λ = 0.6911, **no radiation**, c/H₀ = 2997.92458 Mpc/h:

| comparison | max \|Δa/a\| over the 25 stations |
|---|---|
| Ω_m = 0.3089, matter + Λ only | 7.7e-8 (his integrator's noise) |
| same + radiation (h = 0.6774, N_eff = 3.046) | 1.7e-3 (2.8 Mpc/h at χ = 6400) |
| Ω_m = 0.3156 | 2.9e-2 (so the list *is* cosmology-specific) |

Conclusion: the list is correct; the distance convention is matter + Λ only and
the lightcone/lensing code must use the same χ(a).

**The 42 targets, ascending in a** (17 primary marked `P`, 25 stations `S k`):

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

**Four pairs are closer than the 0.5 % merge tolerance of the exact-z code:**

| pair (z) | a | separation in a | Δa |
|---|---|---|---|
| 2.00000 / 1.98541 | 0.333333 / 0.334963 | 0.486 % | 1.6e-3 |
| 1.75000 / 1.74475 | 0.363636 / 0.364332 | 0.191 % | 7.0e-4 |
| 1.00202 / 1.00000 | 0.499495 / 0.500000 | 0.101 % | 5.1e-4 |
| 0.37562 / 0.37500 | 0.726945 / 0.727273 | 0.045 % | 3.3e-4 |

All other neighbours are ≥ 0.99 % apart. For z_init = 49 and `step da = 4e-4`
the plain ladder has 244 steps and a late-time step of 6.8e-3 (0.68 % of a at
a = 1, 2.1 % at a = 0.33), so the closest pair is ~1/20 of a step.

The older 43-entry list (memory `project_degrace_highres_output_redshifts`, from
bl267's Node_082 ECOSMOG run) is for the **Node_082 cosmology** and must not be
used for the fiducial suite. `ic2pm_val_L1024/Init.dat` used that list in legacy
nearest-step mode (`#outputs = 43`, including z_init = 49).

Reproduce the cross-check (cosemu env):

```python
from scipy.integrate import quad; from scipy.optimize import brentq
cH0, Om = 2997.92458, 0.3089
chi = lambda a: cH0*quad(lambda x: 1/(x*x*((Om/x**3 + 1-Om)**0.5)), a, 1, epsrel=1e-12)[0]
a_station = [brentq(lambda x: chi(x) - 256.0*k, 0.01, 0.9999999) for k in range(1, 26)]
```

## 3. Implementation items in GLAM (branch off `cz`)

### I1. Exact-z schedule: allow targets closer than 0.5 % — **blocker**

Where: `SetExactSteps`, `PMP2main.f90` (≈ line 232); tests in `test/precise_z/`.

Today: `ZTOL = 0.005`; a target within ZTOL of a scheduled step *moves* that step
onto the target, and any two targets closer than ZTOL are rejected with
"closer than the 0.5% merge tolerance". The rejection exists because the second
target of a close pair would otherwise re-move the step already pinned to the
first target.

Change:

1. Never move a step that is already pinned (`MarkX == 1`); fall through to the
   insertion branch instead. Keep ZTOL = 0.5 % for moving *unpinned* steps, so
   most targets still merge (35 of 42 for the ladder above) instead of adding
   PM steps.
2. Relax the pair rejection to "distinct in Real*4" (e.g. relative separation
   > 1e-5); keep the clear error message for true duplicates.
3. Assert after insertion: all `dAlist > 0`, `Alist` strictly increasing,
   exactly `Nexact` steps have `MarkX = 1`, each pinned `Alist` equals the
   Real*4 of 1/(1+z).
4. Restart check unchanged (checkpoint `AEXPN` must match `Alist(ISTEP)`).

Tests to add in `test/precise_z/` (small box, `Nrow=128`): the four Baojiu pairs
→ 42 marked steps, headers bit-exact at 1/(1+z), no duplicate step, one restart
across a pinned pair; existing cases unchanged. Note the generalised leapfrog
`FactV = 2 da_new/(da_old+da_new)` already handles the short inserted steps
(FactV ≈ 0.09 then ≈ 1.9 around the 3.3e-4 step); the known limitation that the
kick factor is evaluated at `AEXPN` rather than the kick mid-point remains.

### I2. Lens-plane output

Where: new routine next to `ProjectDENSIT` in `PMP2mod_density.f90` (≈ line
308), called from `Analysis` in `PMP2mod_analyze.f90`; parameters through
`PMP2init.f90` (parse ≈ line 444, template ≈ line 540, Setup.dat writer ≈ line
330) and `ReadSetup` in `PMP2mod_tools.f90`.

`ProjectDENSIT` already projects the full box onto three NGRID² planes
(sum of `FI+1` along the line of sight, divided by NGRID) but only runs under
`iDumpPart = 1` after `ParticlesRandom`. The lens-plane routine is its slab-wise
generalisation:

- Parameters (append to the **end** of `Init.dat` and `Setup.dat`, read with
  `iostat`, default 0 = off, so every existing run directory stays valid):
  `Lens slabs per axis = 4` and `Lens plane grid = 4096` (must divide NGRID).
  Slab depth = Box / Nslab; refuse if NGRID mod Nslab ≠ 0.
- Algorithm: after `DENSIT` (recompute if `mDENSIT = 1`) and **before**
  `GetPower`, which overwrites `FI` with its FFT. For each axis, slab
  s = (index−1)/(NGRID/Nslab) + 1, accumulate `FI+1` over the NGRID/Nslab cells
  of the slab per column, divide by NGRID/Nslab → mean (1+δ) of the slab column
  (= Σ/Σ̄). Downsample to the plane grid by block averaging. OpenMP over the
  outer loop as in `ProjectDENSIT`. Extra memory: 12 × 4096² × 4 B = 0.8 GB.
- Checks logged per output: each plane averages to 1 within float error; the
  sum over the 4 slabs equals the full-box projection; timing.
- File: one per output, `LensPlanes.<ISTEP>.<realization>.dat`,
  `access='stream'` with an explicit `convert='little_endian'` on the OPEN
  (independent of the build's global endian flags), header
  `(int32 version, real*4 AEXPN, int32 ISTEP, NGRID, Nlens, Nslab, real*4 Box,
  depth, Om, OmL, hubble)`, then the 12 planes in the fixed order axis X slabs
  1…4, Y 1…4, Z 1…4, each Nlens×Nlens real*4 in Fortran column-major order.
  Provide a numpy reader in the halocat tooling.
- Keep `ProjectDENSIT` untouched (still used by the dump-random-fraction path).

### I3. Output schedule for the production

- `Init.dat`: `#outputs = -42` followed by the 42 redshifts (exact mode).
  **Do not list z_init = 49**: exact mode enforces 0 ≤ z < z_init; the IC file is
  the z = 49 snapshot. `zout(1000)` and `NstepM = 2000` are ample.
- Precision: `zout` is Real*4 (≈ 6e-8 relative), sufficient; `PMP2init` writes
  the exact list to `Setup.dat` as `es16.8`.
- Add a pipeline script (next to the supplement generators in `DEGRACE-pilot/`)
  that computes the 25 stations from χ(a) of the run's cosmology (Ω_m from the
  model's `Init.dat`/`PkTable.dat` header, matter + Λ only), merges the 17
  primary redshifts, asserts no exact duplicates, and writes the `#outputs`
  block. Verify it reproduces Baojiu's list for Ω_m = 0.3089 to < 1e-7.

### I4. Step → redshift index file

Snapshot files are named by `ISTEP`, which becomes irregular once targets are
inserted. At run start (after the `SetExactSteps` schedule print) write
`OutputSchedule.dat` with one line per marked step: `ISTEP  a  z  chi[Mpc/h]
exact(0/1)`, χ from matter + Λ. The halocat loaders map snapshot numbers from
this file instead of hard-coded indices.

### I5. Per-output content control (only if Baojiu wants BDM at the primary redshifts only)

`iBDM`, `iPower` and the lens-plane flag are global and apply at every marked
step. If the catalogues are wanted only at the 17 primary redshifts, add a
per-redshift class (e.g. a second list of flags after the redshift list, or a
sign/offset convention) that the output block in `PMP2main.f90` consults.
Defer until section 4 is answered.

### I6. Run configuration and provenance

- `Init.dat`: `iSave=0`, `iPower=1`, `iBDM=1`, lens planes on, `Ncheckpoint`
  small enough that a 72 h job loss does not restart from z_init (the 244-step
  ladder plus ~7 inserted steps).
- Record per box: Init.dat, the `OutputSchedule.dat`, the 2LPTIC/ic2pm
  provenance (velocity epoch line) as required by `AGENTS.md`.

## 4. Open questions for Baojiu / the user

1. BDM halo catalogues at all 42 outputs, or only at the 17 primary redshifts?
   Pilot catalogues are 0.3–0.4 GB per snapshot (MassMin 5e11), so all 42 adds
   ~15 GB per box, ~1.5 TB for 100 boxes.
2. Lens planes at the primary redshifts too, or only at the 26 stations
   (z = 0 plus k = 1…25)? Baojiu's 34 GB budget assumes all 42.
3. Plane normalisation: mean (1+δ) over the slab column (proposed) versus
   surface density; must match his ECOSMOG planes.
4. Do the four close pairs stay (I1 fix) or does he drop one of each pair?
   Differences are ≤ 1.5 Mpc/h in χ.

## 5. Disk budget per box (42 outputs)

| product | per snapshot | per box | 100 boxes |
|---|---|---|---|
| lens planes (12 × 4096² float32) | 0.8 GB | 34 GB | 3.4 TB |
| `CatshortV` halo catalogues | 0.3–0.4 GB | ~15 GB | ~1.5 TB |
| `PowerDM` | KB | negligible | negligible |
| particle snapshots (`iSave=1`) | 200 GB | 8.4 TB | not an option |
| checkpoints `PMcrs` (transient) | 200 GB | 200 GB live | deleted at end |

## 6. Checklist

- [x] Cross-check of the 42 scale factors (2026-10-16)
- [ ] I1 ZTOL / pinned-step fix + tests
- [ ] I2 lens-plane routine, parameters, file format, numpy reader
- [ ] I3 schedule generator script; production `Init.dat`
- [ ] I4 `OutputSchedule.dat`
- [ ] I5 per-output content control (pending question 1)
- [ ] Pilot: one L1024 box with the 42-output schedule, lens planes and BDM;
      check runtime, CPU efficiency (`sacct`), file sizes, plane means
- [ ] Reply to Baojiu with the cross-check result and questions 1–4

Related memory files: `project_degrace_lensplanes_plan`,
`project_exact_z_snapshots`, `project_degrace_highres_output_redshifts`.
