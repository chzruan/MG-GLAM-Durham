# REVIEW 1: Fortran extras module (lens planes, Fourier dumps, momentum, f_R at haloes)

Branch `degrace-track-a` vs `cz` (merge base 69c9261). Scope: `PMP2mod_extras.f90`
(module `ExtraOutputs`, 668 lines, new) and its hooks in `PMP2mod_power.f90:341-343`,
`PMP2mod_analyze.f90:46-54,126-132,202`, `PMP2linker.f90:119,209-211,624`,
`PMP2main.f90:82,122` (`WriteOutputSchedule` call). Read-only review, 2026-10-17;
nothing was built or run except a few single-threaded Python arithmetic checks.

## 1. Summary

The module is correct in every physics-bearing path I could check: the packed
half-complex layout and the sign of the imaginary part are derived independently
from FFTPACK5 (`rfftf1`'s `tsnm = -tsn` makes r(2m+1) = +(2/N) Σ x sin, so
`(r_cos − i r_sin)/2` is the e^{−ik·x} coefficient), `FFT3Dforward` is a
line-for-line copy of `POWERfft5`'s passes, the lens-plane loops match the
documented orientation and conserve mass, the momentum deposit reproduces
DENSIT's CIC and wrap on the coarser grid with the leapfrog-consistent velocity
epoch, the f_R conversion `fR_bg(a_old)·FI2^(n+1)` agrees with `POTENT_FOR_MG`
and `scalar_field`, FI2 is intact at BDM time, and no unit number or file name
collides with live GLAM I/O. No Blocker or Major. The Minor items are gating
holes (non-f(R) MG models), a hard-coded catalogue label in the fR_halo meta, a
missing fR_halo file for empty catalogues, a documented-but-unresolved velocity
epoch mismatch with the BDM catalogue, a spec deviation on momentum interlacing,
missing provenance keys, and an untimed atomic deposit at 2048³.
**Recommendation: ready after the listed fixes** (R1-04, R1-05, R1-07 are
one-line changes; R1-02, R1-03, R1-09 need a documented decision).

## 2. Findings

| ID | Severity | File:line | Description | Evidence | Suggested fix |
|---|---|---|---|---|---|
| R1-01 | Minor | `PMP2mod_extras.f90:465-514` | `MomentumDeposit` does 8 `!$OMP ATOMIC` read-modify-writes per particle into the shared `G(Nd,Nd,Nd)` with random addresses; correct, but 3 × 8 × N_p = 2.1e11 atomics per output at 2048³. Estimated 1-3 min per output (DRAM-latency bound, 128 threads), i.e. 1-2 h over 42 outputs, 15-30 % of the 6 h ΛCDM-twin wall time quoted in plan Q8; V9 (< 5 %) is still open. At 1024³ about 10-20 s per output. | Lines 498-513 (eight `!$OMP ATOMIC`); `DENSIT` (`PMP2mod_density.f90:17-301`) uses a split decomposition (`iSplit = 8`) without atomics. | Reuse DENSIT's decomposition or give each thread ownership of a range of K-slabs (two passes or a private boundary layer); measure in V9 before the 2048³ campaign. |
| R1-02 | Minor | `PMP2mod_extras.f90:542-548`, `PMP2init.f90:276-279` | Default N_dep = N_g/2 with k_c,mom = 1.0 h/Mpc at L = 1024 (k_N,dep = 3.14 h/Mpc) gives k_c/k_N = 0.32 > 1/4, exactly the regime where plan I8 item 3 says "interlace ... and set interlaced = true"; the code never interlaces and writes `interlaced = false`. PMP2init checks only k_c < k_N. Alias estimate at k → k_c, axis-aligned worst case: first image at 2k_N − k = 1.68 k_N weighted by W_CIC = sinc²(0.84π) ≈ 0.033 times the amplitude ratio of the field, i.e. ≲ 1 % of the mode amplitude with random phase (noise, not bias, in cross spectra). | `write (30,'(a)') 'interlaced = false'` (line 396); plan `docs/degrace_companion_implementation_plan.md` I8 item 3. | Either document the deviation in the plan (the k_N/4 rule is conservative for CIC modes) or lower the default `kc_mom` to ≤ 0.78 or default N_dep = N_g with FI reuse. |
| R1-03 | Minor | `PMP2mod_extras.f90:461,531`; `PMP2linker.f90:537` | Momentum products convert with a_v = a − da/2 (correct: `MOVE` drifts with `XCONST` at `Ahalf`, `PMP2main.f90:524-560`, so VX is the half-step momentum), but the BDM catalogue converts with `Vscale = 100 Xscale/AEXPN` (post-ADDTIME a), as do `ParticlesRandom` (`PMP2mod_analyze.f90:312`) and `DENSITrsd` (`PMP2mod_density.f90:604`). Plan section 3 requires I8, I9 and the BDM catalogue to use the same (x, v) pair. Mismatch factor a/(a − da/2) ≈ 1 + da/(2a) ≈ 0.3-0.6 % on the late ladder. Legacy behaviour, not introduced here, but the new products' docs claim consistency. Recoverable from the catalogue header (A, Step=ASTEP). | `vfac = 100.*Box/NGRID/a_v*...` vs `Vscale=100.d0*Xscale/dble(AEXPN)`. | Record `catalogue_velocity_epoch = a` (and the factor) in the mom_*/fR_halo meta and in `DEGRACE-companion/README.md`; changing BDM's Vscale would alter legacy catalogues, so leave the code. |
| R1-04 | Minor | `PMP2linker.f90:209`; `PMP2init.f90:282` | fR_halo gate is `iExtra==1 .and. FRhalo==1 .and. MG_flag==1 .and. OutClass==1`; PMP2init only zeroes `FRhalo` when `MG_flag /= 1`. For DGP/symmetron/k-mouflage/CSF (`MG_model` 2-5) with `fR at halos = 1` the code silently writes an `fR_halo` file from the wrong field (FI2 is that model's scalar) scaled by `fR0*(R0/R)**(fr_n+1)` with the f(R) parameters of a non-f(R) run. | `fRHaloOpen` lines 564-577 use `fR0`, `fr_n` unconditionally. | Add `.and. MG_model == 1` to both gates (and to the ReadSetup/PMP2init validation). |
| R1-05 | Minor | `PMP2mod_extras.f90:644-645` | fR_halo `.meta` hard-codes `catalogue = CATALOGS/CatshortV.<step>.<box>.DAT`, but the label is `W.`/`V.`/`M.`/`A.` by `iVirial` (`PMP2linker.f90:355-363`). Wrong reference for any non-virial catalogue. | `write (Name,'(2(a,i4.4),a)') 'CATALOGS/CatshortV.', fRmoment, ...` | Pass `CatalogueFinalPath` (module LinkerList) to `fRHaloClose`, or record `iVirial`. |
| R1-06 | Minor | `PMP2mod_extras.f90:185-217,386-418,642-664` | Plan section 3 lists `git_hash` and `model` as required meta keys for every product; none of the three sidecars writes them (`header`, `MG_flag`, `fR0` only), and the LensPlanes meta also omits `D_code(a)`. Provenance of which binary produced a file is therefore only in the Slurm log. | Compare the `write (30,...)` key lists with plan lines 155-160. | Inject a build-time hash (`-DGIT_HASH=...` in the makefile, or a generated `version.inc`) and write `git_hash`, `MG_model`, `D_code(a)` in all sidecars; or amend the plan's required-key list. |
| R1-07 | Minor | `PMP2linker.f90:152-156` vs `209-211` | BDM's empty-catalogue path (`Nmaxima == 0`) calls `WriteFiles` and returns before `fRHaloOpen` is reached, so a header-only `CatshortV` is published with no `fR_halo.<step>.<box>.{bin,meta}`; `fRHaloClose` is never called (harmless, flag false). Downstream code that expects one fR_halo per catalogue will fail on such a step (plausible at the earliest outputs of small test boxes, cf. "BDM found 1 and 5 haloes" in plan section 16). | Early return precedes the open at line 209. | Call `fRHaloOpen` before the `Nmaxima == 0` block (and `fRHaloClose(0)` inside it) so an empty .bin and a `rows = 0` .meta are written. |
| R1-08 | Minor | `PMP2mod_extras.f90:573,629,637` | `fR_halo.bin` is written directly under its final name while the catalogue uses staged-then-rename (`PublishCatalogue`, `PMP2linker.f90:416-441`); any `error stop` inside `WriteFiles` leaves a truncated fR_halo.bin with no .meta. Also unit 31 stays open across the whole catalogue loop (fine, no other unit 31 in GLAM). | Compare with `BeginCataloguePublication`. | Write to `<name>.part` and rename in `fRHaloClose`, or state in the README that the .meta is the completion marker. |
| R1-09 | Minor | `PMP2mod_extras.f90:542-544` | `ExtraModes` allocates a separate transient `G(Nd,Nd,Nd)` (4.3 GB at N_dep = 1024; 34 GB at the N_g = 4096 default N_g/2) instead of reusing `FI` as plan section 4 / I8 item 3 prefer; plan ground rule 6 and Q8 say f(R) at N_g = 4096 has ~20 GB headroom and "not even the N_g/2 fallback" fits. No guard: that configuration would OOM at the first class-1 output after hours. Companion configs (N_g = 2048) are fine. The allocation is also not reported through `Memory()`. | `allocate (G(Nd, Nd, Nd))`; plan lines 43-45 and Q8. | Refuse at init (`iExtra==1 .and. MG_flag==1 .and. NdepMom**3*4 > headroom`) or document the limit; call `Memory()` around the allocation. |
| R1-10 | Nit | `PMP2mod_extras.f90:14,267-268` | `If (N > Nlenmax) Stop` does not cover N = Nlenmax = 16384, where `lensav = N + int(log2 N) + 4 = 16402 > size(wsave)`; `rfft1i` would then write past `wsave`. Same latent limit in `POWERfft5` (`Nlensav = 16384`). Irrelevant for N ≤ 8192. | Checked: lensav(16384) = 16402. | Test `lensav > Nlenmax` instead. |
| R1-11 | Nit | `PMP2mod_extras.f90:188` | Meta says pixel (i,j) covers `[(i-1), i) * Box/N_lens`. With block factor f the pixel is the mean over nodes at ((i−1)f + m) Box/N_g, m = 0..f−1, whose centre is [(i−1) + (f−1)/(2f)] Box/N_lens, not the pixel centre; with f = 1 the sample point is the pixel's left edge (node). Offset Box/(2N_g) = 0.25 Mpc/h at N_g = 2048, a linear-in-k phase in plane × halo cross-correlations if the user assumes pixel centres. The node statement (line 84 comment, `fourier_norm` in the modes meta) is present but implicit here. | Loop `ip = (iy-1)/f + 1` (line 125). | Add `pixel_centre = ((i-1) + (f-1)/(2f)) * Box/N_lens` (or `sample = node`) to the meta. |
| R1-12 | Nit | `PMP2mod_extras.f90:64,115,184,381,385,573,641` | No `iostat` on any OPEN/WRITE; a full disk or a missing run directory aborts the simulation with a runtime error (legacy GLAM behaves the same). The `Stop` at line 108 fires only for hand-edited `Setup.dat` since PMP2init validates divisibility (`PMP2init.f90:264-268`). | — | Optional: `iostat` + a clear message. |
| R1-13 | Nit | `PMP2mod_extras.f90:498-513` | Float atomics make `mom_*` summation order thread-dependent, so mom_* are not bit-reproducible across runs or thread counts (restart replays overwrite with slightly different bits). Fine for science. | — | Exclude mom_* from any md5-based determinism test; note in README. |
| R1-14 | Nit | `PMP2mod_extras.f90:211-216,412-417` | Sidecars write `MG_flag` and `fR0` but not `MG_model`; for MG_model 2-5 the `fR0 = <value>` line is misleading. | — | Write `MG_model` and print `fR0` only when `MG_model == 1`. |

## 3. Verified OK

- **FFTPACK5 packing and sign (independent derivation).** `rfftf1`
  (`PMP2mod_fft5.f90:1937-1965`) scales entry 1 and N by 1/N, even entries by
  2/N and odd entries by `tsnm = -2/N`. The radix-4 kernel `r1f4kf`
  (`PMP2mod_fft5.f90:1382-1476`, ido = 1 block) gives raw r(3) = x4 − x2 =
  −Σ x sin(2π(j−1)/4) for N = 4, so after scaling r(2m+1) = +(2/N) Σ x sin,
  r(2m) = (2/N) Σ x cos, r(1) = (1/N) Σ x, r(N) = (1/N) Σ (−1)^(j−1) x. Hence
  X_m = (1/N) Σ x e^{−2πi m(j−1)/N} = r(2m)/2 − i r(2m+1)/2 exactly as
  `PackedBasis` encodes (w = 0.5, −0.5 i sign(m)); the √2 factors on p = 1 and
  p = N are undone by w = 1/√2, giving X_0 = (1/N) Σ x and the Nyquist term.
  X_{−m} = conj(X_m) follows from w(2) flipping sign with real r. The V1
  expectation in `test/extras/check_modes.py` (sin field → −0.25 i at +m) is
  the correct value for this convention, so that test does pin the sign.
- **3-D separability with packed real intermediates.** The y-pass transforms
  each real sequence `G(px,:,k)`; by linearity X_{mx,my}(k) = Σ_i Σ_j
  w_x(i) w_y(j) G(px(i),py(j),k), and likewise after the z-pass. Transposition
  swaps indices 1 and 3 (`G(i,j,k) ↔ G(k,j,i)`), the z-pass runs along the new
  first index, so the coefficient sits at `G(pz,py,px)`, matching
  `val = val + wx(i)*wy(j)*wz(k)*G(pz(k), py(j), px(i))`. Overall
  normalisation N^−3; `fourier_norm` meta line is right.
- **FFT3Dforward vs POWERfft5** (`PMP2mod_extras.f90:261-337` vs
  `PMP2mod_power.f90:264-339`): identical loop nests, √2 scalings,
  transposition and `lensav` formula; only `NGRID → N`. Works for N ≠ NGRID:
  `rfft1i` re-initialises the module `wsave` per call, `COPYIN(wsave,work)`
  on each parallel region, `r(Nlenmax)` private, lenwrk = N ≤ Nlenmax.
- **DumpModes bounds and k = 0.** n = floor(k_c L/2π) (244 for k_c = 1.5,
  162 for 1.0 at L = 1024), clamped to N/2 − 1 with a message; px ≤ 2n+1 ≤
  N−1 in bounds; `C(0,0,0) = 0` after the loop; iz = 0 modes stored for both
  ±(mx,my) and are exact conjugates. Fortran order `C(0:n,-n:n,-n:n)` →
  numpy `reshape(2n+1,2n+1,n+1)[ix+n, iy+n, iz]` as the meta says. Cube size
  (2n+1)²(n+1)·8 B = 0.469 GB (k_c 1.5), 0.138 GB (k_c 1.0), matches plan §9.
  Hook in `POWERfft5` sits after the third √2 scaling and before `FI**2/8`,
  gated identically to the `iPowerDone` argument in `Analysis:129`
  (`iPower==1 .and. iPowerRSD==0`); `GetPowerGal` has no callers, so no
  galaxy field can be dumped as `deltam`.
- **ExtraModes.** Runs only for `iExtra==1 .and. OutClass==1`
  (`PMP2mod_analyze.f90:128`); re-deposits only when `mDENSIT==1`, sets
  `mDENSIT=1` after destroying FI; `a_x = AEXPN`, `a_v = AEXPN − ASTEP/2`
  written for every product; `D_code(a)` present in the modes meta.
- **LensPlanes** (`PMP2mod_extras.f90:94-219`). FI holds δ at the call
  (`DENSIT` sets FI = −1 then adds CIC with `Wpar = NGRID³/N_p`,
  `PMP2mod_density.f90:36,25`); `Analysis:46-54` deposits if `mDENSIT==1`
  before anything else touches FI. Orientation: case 1 → acc(ip←iy, jp←iz,
  s←ix) = X-planes (y,z); case 2 → (ip←ix, jp←iz, s←iy) = Y-planes (x,z);
  case 3 → (ip←ix, jp←iy, s←iz) = Z-planes (x,y); the `write(29) plane`
  column-major order makes the first index fastest, as the meta states.
  `wnorm = 1/(nc f²)` is the cell count per pixel; `1 + acc*wnorm` is
  mean(1+δ); Σ over all pixels of (plane−1) = wnorm Σ FI = 0 up to rounding,
  so mass is conserved and the printed axis means test it. Threads own
  disjoint `jp` rows in every case (no race); accumulation real*8. Cost at
  N_g = 2048: three streaming passes over FI (~100 GB read) and a 268 MB
  `acc` (Ns = 4); negligible. File size 12 × 2048² × 4 B = 201 MB.
  Divisibility is validated at init (`PMP2init.f90:264-268`); `LensGrid = 0`
  → N_g handled in both places.
- **MomentumDeposit** (`454-515`). `X = (XPAR−1)·sc + 1` maps [1, NGRID+1)
  onto [1, Nd+1); the `XN = Nd+1−1e-7` wrap and `I1 > Nd → 1` are DENSIT's
  (`PMP2mod_density.f90:23,88-93`); CIC weights identical to DENSIT.
  `vfac = (100 Box/NGRID)/a_v · Nd³/N_p` gives (1+δ) v in km/s: GLAM's drift
  `XCONST = ASTEP/(Ahalf·sqrt(Ahalf(Om+OmL Ahalf³)))` (`MOVE`) shows VX is
  the momentum at a_half, so dividing by a − da/2 is the consistent epoch
  (the legacy `1/AEXPN` is the da/(2a) bias of the ic2pm erratum class).
  `Integer*8` loop variable; atomics make the reduction correct. Transient
  memory 4.3 GB (G) + 0.14 GB (C) at N_dep = 1024.
- **f_R at haloes.** `scalar_field` (`PMP2extradof.f90`) multiplies FI2 by
  |fR_bg|^{−1/(n+1)} after convergence (and back by |fR_bg|^{1/(n+1)} as the
  next step's guess), so FI2 = u/u_bg; `POTENT_FOR_MG` (`PMP2main.f90:1105-1155`)
  adds −½ c̃² fR_bg FI2^(n+1) = −(c̃²/2) f_R, confirming f_R = fR_bg·FI2^(n+1)
  with fR_bg = fR0 (R0/R)^(n+1), R = 3(Om/a³ + 4 OmL) (same formulas in the
  solver, `PMP2MGsolver_fR.f90:59-61`). Sign consistent (FI2 ≥ 0, fR0 < 0).
  Epoch: DENSIT and `scalar_field` run at the top of the step
  (`PMP2main.f90:44-47`) before `MOVE`/`ADDTIME` (69-71), and `ASTEP` is
  updated before `MOVE`, so a_old = AEXPN − ASTEP (exact-z: Alist(ISTEP) −
  dAlist(ISTEP) = Alist(ISTEP−1) to real*4 rounding). FI2 survives to BDM:
  `MG_scratch_release` frees only FI3 (`PMP2mod_tools.f90:605-616`); no
  FI2 reference in `PMP2linker.f90`, `PMP2mod_density.f90`, `WriteDataPM`.
  Trilinear interpolation: X1 = modulo(x N/Box, N) + 1 matches the catalogue
  convention x_cat = (XPAR−1) Box/NGRID (`RescaleCoords`,
  `PMP2linker.f90:537-539`); I clamp at N+1 → weight on the periodic image;
  wraps on I1/J1/K1 and in the sphere loop; cube ±nr covers the sphere;
  ncell = 0 falls back to the centre value. Catalogue rows pass x,y,z in
  Mpc/h in [0, Box) (Xleft = 0, Xright = Box, `PMP2linker.f90:1583`) and
  `Rvir(ip)` in Mpc/h (real*4; the catalogue column is `1e3*Rvir`). Cost
  negligible (≤ 15³ cells per cluster-size halo). Row count: `fRHaloRow` is
  called right after every `write(12,...)` in the single row loop
  (`PMP2linker.f90:620-624`), `WriteFiles` runs once per step, `Nhalo = iHalo`
  compared in `fRHaloClose`. `fRHaloActive` is set only in `fRHaloOpen` and
  cleared in `fRHaloClose`, which `BDM` calls unconditionally and which
  returns early when inactive. `ReadParameters(ISTEP)` names the catalogue by
  `ISTEP`, matching `fRHaloOpen(ISTEP)`.
- **GrowthFactor / ComovingDistance.** Substitution x = ln a' turns
  D ∝ E(a) ∫ da'/(a'E)³ into E(a) ∫ e^{2.5x}/(Om + OmL e^{3x})^{1.5} dx as coded;
  χ = (c/H0) ∫ dx/(a E) from ln a to 0 with h = −x0/4000 > 0; Simpson weights
  and even interval count correct. Numerical check (Python, Om = 0.3): D(0.02),
  D(0.5) agree with a 2e5-interval reference from a = 1e-12 to 2e-11 relative;
  χ(0.02, 0.5, 0.9) to 1e-14. Lower limit 1e-8 contributes O(1e-20). All
  arithmetic in real*8 with `dble()` on the real*4 inputs; c/H0 = 2997.92458.
  Radiation and neutrinos omitted, as stated in the header.
- **WriteOutputSchedule.** Legacy mode: `Nclass(jm)` = class of the zout
  moment, fuller class wins on shared steps. Exact mode: `iord` is moved
  with `atgt` in the insertion sort, inserted steps shift `Nclass` with the
  other arrays, and `Nclass(j/jnear) = zclass(iord(i))`; pinned steps are
  never re-moved (`MarkX` guard), so the class lands on the right epoch.
  `max(1, Nclass)` covers Nclasses = 0. `status='unknown'` with no `append`
  rewrites the file on restart. Format `(i6,2es16.8,f13.4,2i6)` fits
  i ≤ 2000, χ < 1e5.
- **Units and files.** Unit 29: only here (plus commented-out code in
  `PMP2main.f90:1094-1100`). Unit 30: `Timing`/`TimingMain` open and close
  it only at `isign == 0` (`PMP2mod_tools.f90:562-587`), never while an
  extras sidecar is open; `PMP2start*` use 30+ifile in a separate program.
  Unit 31: only here. Units 12/13/17/18 open during the hooks are untouched.
  File names `LensPlanes.`, `deltam.`, `mom_[xyz].`, `fR_halo.`,
  `OutputSchedule.dat` collide with nothing in GLAM; `<step i4.4>.<box i4.4>`
  matches the legacy `PowerDM.R.<step>.<box>.dat` widths. `status='replace'`
  everywhere, so a restart replay overwrites rather than appends.
  `convert='little_endian'` on every stream OPEN (build is big-endian).
- **Hygiene.** `implicit none` at module level; every dummy has `intent`;
  `acc`, `plane`, `C`, `G` deallocated on every path; `PackedBasis` is pure;
  all privates listed in the OpenMP clauses; `Complex*8` = complex64 as the
  meta states; `Integer*8` for `Nparticles`, `system_clock` and the particle
  loop; no `Nclass`/`OutClass` default changes legacy behaviour (`OutClass`
  = 1 at start, set only on marked steps; the BDM guard `OutClass==1` is
  always true when no class list is given).

## 4. Open questions for the author

1. Will any companion configuration run f(R) at N_g = 4096 (R1-09)? If yes,
   the separate N_dep³ buffer must become FI reuse as the plan requires; if
   no, record the limit in the plan and README.
2. R1-03: do you want to state explicitly that the BDM catalogue velocities
   are at a (legacy `1/AEXPN`) while mom_* and the subsample are at a − da/2,
   and give the correction factor, or change the catalogue for `iExtra = 1`
   runs only (which would make catalogue content flag-dependent)?
3. R1-02: is k_c,mom = 1.0 with N_dep = N_g/2 and no interlacing the intended
   default, given plan I8 item 3? If so, amend the plan with the alias
   estimate above.
4. R1-01: has `MomentumDeposit` been timed at 2048³ (V9)? The atomic version
   is probably the dominant extras cost at that size.
5. R1-06: how will `git_hash` get into the sidecars (makefile `-D`, generated
   include, or dropping the requirement from plan section 3)?
6. R1-07: does `halocat` expect an `fR_halo` file for every published
   catalogue, including header-only ones?
