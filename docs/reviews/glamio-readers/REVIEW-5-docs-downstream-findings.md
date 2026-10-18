# REVIEW 5: documentation, conventions and downstream use — findings

Branch `glamio-readers` (commit 0a6ecff vs `cz`), reviewed 2026-10-18 on COSMA
(login node, read-only; Python checks with the cosemu interpreter against
`test/precise_z/work/case10_lens72`, `case11_modes`, `case12_frhalo` and the
LCDM pilot's `LensPlanes.0057.0001`). Scope: python/README.md, the plan's
section 16 entry and sections 3, 4, 6, 9, 11, test/precise_z/README.md, the
docstrings of the eight new modules, the status deck's definition frames, and
the consumers (DEGRACE-companion/check_run.py, test/extras/check_*.py, the
halocat ingest as the plan describes it). Line numbers refer to the branch.

## 1. Summary

The README, the module docstrings and the plan entry describe the readers
accurately in almost every statement I could test: the usage example runs on
real output, `Modes.power()` reproduces `PowerDM.R` to 4e-6 in all 11 bins of
case11, the conventions bullets match the Fortran sidecars and the deck, the
section 16 counts (6 suite checks, 8 + 7 unit tests, 529-byte record, pilot
planes memory-mapped with axis means 1) are verifiable from the diff and the
files. One factual error is repeated in the new `lensplanes.py` docstring, the
deck and the plan: "the mean of every plane is 1 up to float32 rounding". Only
the per-axis mean is 1; individual planes deviate by 1-6 % on the suite's
150 Mpc/h box and by up to 1e-3 on the pilot at z = 8.7 (the Fortran itself
logs per-plane means without asserting them and asserts only the axis means).
The remaining findings are documentation precision (the README's generic
paragraph does not hold for `Schedule`, `PowerDM`, `Snapshot`; the example
skips the lens planes of the 25 station outputs; stale plan/deck wording on
the V2 deconvolution and the sidecar key names) and API gaps the halocat ingest
will hit (no step/redshift-to-path helper, no ID-to-lattice mapping).
Recommendation: **ready after the listed fixes** (M-1 and the Minor doc items;
the API additions can follow in a later commit).

## 2. Findings

| ID | Severity | File:line | Description | Evidence | Suggested fix |
|---|---|---|---|---|---|
| M-1 | Major | python/glamio/lensplanes.py:5-7; docs/slides/.../degrace_companion_status.tex:284, 303; plan :228, :523 | The module docstring states "the mean of every plane is 1 up to float32 rounding"; the deck says "Each plane averages to 1 by construction" and tabulates "pilots: plane means 1 ± 1e-7"; plan I2 says "each plane averages to 1 within float error". Only the mean over the N_slab planes of one axis is 1 by construction (mass conservation); a single slab's plane has mean 1 + ⟨δ⟩_slab. A user who normalises each plane to its own mean, or asserts per-plane mean = 1 as QA, gets a wrong Σ/Σ̄ or rejects every valid file. | Measured with the reader: case10 per-plane means − 1 = [−0.01, 0.01, 0.03, −0.03, 0.02, −0.02] at z = 2.5 and [−0.02, 0.02, 0.06, −0.06, 0.05, −0.05] at z = 0.512, while axis means − 1 ≈ 5e-10; pilot (2048², z = 8.73) per-plane means − 1 up to 9.5e-4, axis means − 1 ≈ 1e-11. PMP2mod_extras.f90:171-174 logs `mean(1+delta)` per plane without a check and :180 asserts only `axis means (must be 1)`. `axis_means()` (lensplanes.py:101) and check_run.py:98 (per-axis mean over the slabs) are correct; the deck's table row labels the per-axis number "plane means". | Docstring: "the mean over the N_slab planes of an axis is 1 up to rounding (`axis_means()`); a single plane has mean 1 + ⟨δ⟩_slab". Fix the deck (:284, :303 → "axis means") and plan I2/L1 wording. |
| m-1 | Minor | python/README.md:92-98 | The "Other products" lead paragraph says *all* per-output products share the `<product>.<step>.<box>.{bin,meta}` scheme, accept `.bin`/`.meta`/stem, take `mmap=True`, and expose `.meta` plus the common scalar properties. None of this holds for `read_schedule` (text, no sidecar), `read_powerdm` (text, no `mmap`, no `.meta`, no `box`/`realization`/`n_grid`), `read_pm_header`/`read_snapshot` (big-endian Fortran records with markers, `Snapshot.header` instead of `.meta`). | `issubclass(PowerDM|Schedule|Snapshot, MetaProduct)` is False; `hasattr(PowerDM, "box")` is False; `read_powerdm(path)` has no `mmap` argument (powerdm.py:77). | State the scheme for the four sidecar products (`LensPlanes`, `Modes`, `Subsample`, `FRHalo`) and describe the three text/legacy readers separately. |
| m-2 | Minor | python/README.md:114-121 | The usage example reads `LensPlanes` inside `for step in sched.primary()["step"]`. Lens planes (and `PowerDM`) are written at every output, including the 25 class-2 stations; the extras (`deltam`, `Subsample`, catalogues) only at class 1. A reader following the example silently skips 25 of the 42 production planes; on case11 the loop raises `FileNotFoundError` at the first iteration (no planes in that case). | Pilot `OutputSchedule.dat`: steps 57-86 are class 2 and have `LensPlanes.*`; schedule.py:112-118 docstrings say class 2 = "lens planes and P(k) only". Running the example on case11: `FileNotFoundError: sidecar not found: .../LensPlanes.0070.0001.meta`. | Loop over `sched.steps` for planes/PowerDM and over `sched.primary()["step"]` for the extras; say so in a comment. |
| m-3 | Minor | deck :347 (and :323); plan :532 (V2 row) | The deck's definitions frame states "Offline P_mm(k) = L³⟨\|δ̃_m/W_CIC\|²⟩ equals PowerDM.R bin by bin to 1e-4 (V2)"; the plan's V2 row says "P(k) from deltam (CIC deconvolved offline) vs the code's own P_mm(k)". Both contradict the readers' documentation (powerdm.py:10-11 and README:147-149: PowerDM is raw CIC, not deconvolved) and the as-run checks. | asserts.py:470-485 (V2) and :316-325 (`binned_power`) compare the raw cube power with `PowerDM.R`, no window; `Modes.power()` (no deconvolution) agrees with `PowerDM.R` to 4.1e-6 in all 11 bins of case11 (`nmodes` equal, `k` to 5e-6). Dividing by W_CIC would break the agreement at the 10-30 % level near k_c. | Drop `/W_CIC` from the deck sentence and "CIC deconvolved offline" from the plan's V2 row (both sides raw). |
| m-4 | Minor | plan :152-157 (section 3), :607 (Q9) | Section 3 lists the "required" sidecar keys as `git_hash, model, L, N_g, box, ...` and file names `extra/<run>_a<a>_<product>`; the writers (and therefore `MetaProduct`, README:96-98) use `Box`, `N_grid`, `N_row`, `N_p`, `MG_flag`/`MG_model`, no `git_hash`, and the `<product>.<step>.<box>` scheme. Q9 (file naming) is still listed as open although decided by the implementation. | `product.py:76-103` reads `Box`, `N_grid`, `N_row`, `N_p`; the real sidecars (case10/11/12, pilot) have exactly those keys; the README documents `<product>.<step>.<box>.{bin,meta}`. | Add an "as implemented" note to section 3 (key names, no `git_hash`, naming scheme) or point to the README; close Q9. |
| m-5 | Minor | python/glamio/subsample.py:119-125; deck :336 | `Subsample.lattice` returns `'IDx'`, `'IDz'`, `'ID?'` or `None` and the docstring does not say what `None` means. Both pilots and every native-IC suite run return `None` (their headers carry no tag), and the sidecar's own `lattice` line distinguishes "no tag = pre-fix ic2pm, IC file order" from "native PMP2start: ID = 1 + i + Np (j + Np k)". There is also no function mapping IDs to lattice indices / Lagrangian positions although the deck advertises "the subsample index gives the Lagrangian position" and the plan (I9) names it a purpose of the product. | Pilot header `ic2pm: 2LPTic ingest, v at a-da/2` (no tag; the tag was added at ic2pm.f90:340-346 later); case11 `S.lattice` is `None`. The formulas exist only in test/extras/check_ic2pm_order.py:24-29. | Document `None` (pre-tag ic2pm or native IC) and add `Subsample.lattice_indices(ids)` / `lagrangian_positions()` using `N_side`, `Box` and the tag (raise when the order is undetermined). |
| m-6 | Minor | API (python/glamio/__init__.py, schedule.py) vs plan :255, :460-464 | What the halocat ingest needs and the API lacks: (a) a run-directory helper that builds `<product>.<step:04d>.<box:04d>` paths and lists the products present at a step (the README example hand-builds the names); (b) a redshift-to-step lookup (halocat's loaders are keyed by redshift; `Schedule` offers only `row(step)`/`by_step()`); (c) joins by step between catalogue, `fR_halo` and `Subsample` (only `FRHalo.catalogue()` exists). | README:116-120 repeats `f"Run1/<product>.{step:04d}.0001"` four times; schedule.py:101-118. | Add `Schedule.step_at(z, atol=...)` and a small `glamio.run` module (`product_path(rundir, product, step, box=1)`, `products(rundir, step)`); mention it in the plan's halocat note. |
| m-7 | Minor | python/glamio/snapshot.py:42-43, :111-113 | `RHO_CRIT = 2.7754e11` differs from the finder's `2.774e11` (PMP2linker.f90, `Dscale`), so `PMHeader.particle_mass()` disagrees with the catalogue sidecar's `particle_mass` for the same run by 5.0e-4 and `Npart * particle_mass()` does not reproduce `Mbound`. | case12: header `particle_mass()` = 6.2017e12, catalogue sidecar `particle_mass = 6.19855688E+12` (ratio 1.000505). | Use the code's constant (so both numbers agree) and name the constant `RHO_CRIT_GLAM` with a comment, or state the 0.05 % difference in the docstring and README. |
| m-8 | Minor | python/README.md:109; python/glamio/snapshot.py:26-27 | "`read()` loads all (24 B per particle)" / "a full read costs 24 bytes per particle (float32)" holds for `units='code'` only; `units='physical'` returns six float64 arrays (48 B per particle, plus the float32 page temporaries). A user sizing a login-node read from the README underestimates by 2x (412 GB vs 206 GB for 2048³). | Measured: `snap.read(units="physical")` → six `float64` arrays, 48.0 B per particle; `units="code"` → `float32`, 24.0 B. | Say "24 B (`code`, float32) or 48 B (`physical`, float64) per particle". |
| m-9 | Minor | plan :1055-1056 ("Not done here"); deck :227-253 | The entry lists only the halocat ingest and `ProjectionXY/XZ` as not covered. Still unreadable through glamio and worth listing: `timing.log` (parsed by check_run.py:143-162 and listed by the deck as a per-box file), `Run.log`/`main.log`, the BDM `outputB.<step>.dat` logs, `ProjectionYZ` (the run writes XY, XZ and YZ), `PowerDM.Z` (never produced in DEGRACE runs: PMP2mod_power.f90:77 writes it only in the redshift-space branch with the name `PowerDM.Z.<moment>.<real>.V<sigv>.d<thr>.dat`, so the README's `PowerDM.R/log/Z.*.dat` row and the `Z` column mapping are untested on a real file). | `ls test/precise_z/work/case10_lens72/Run1`: ProjectionXY/XZ/YZ, timing.log, Run.log, main.log, pt.dat; case12 `CATALOGS/outputB.0084.dat`; no `PowerDM.Z.*` anywhere. | List them in the plan entry with a one-word verdict each (timing.log: next reader, used by check_run.py; logs: not needed; PowerDM.Z: untested). |
| m-10 | Minor | DEGRACE-companion/check_run.py:19-24, 55-58, 86-100, 111-116 | check_run.py re-implements `read_meta`, the schedule parser (no column-count or ordering check), the lens-plane decode and per-axis means, the PowerDM epoch regex and the catalogue sidecar handling. The plan/README do not say whether it is meant to stay an independent decoder (like the suite's `asserts.py`, which must) or should now call glamio and gain the byte-count/product validation. | The suite and test/extras oracles must stay independent (preamble); check_run.py is a production verification tool, not a test oracle, and is not mentioned in section 6 or 16. | Decide and record it in the plan: either `try: import glamio` with the current code as fallback, or keep it independent on purpose ("second decoder on every production box") and say so in its docstring. |
| m-11 | Minor | python/README.md:17; installed metadata | The install line is still valid, but the editable install in cosemu reports the pre-bump metadata (`pip show glamio` → 0.1.0) while `glamio.__version__` is 0.2.0, because `dynamic = ["version"]` is evaluated at install time only. Anyone checking the version with pip will be misled until the author re-runs the install line. | `nice -n 19 .../cosemu/bin/python3 -I -m pip show glamio` → `Version: 0.1.0`; `glamio.__version__` → `0.2.0`. | Re-run `pip install -e python/` after the merge and add "re-run after a version bump" to the install section. |
| n-1 | Nit | python/README.md:130-133 | The lens-plane bullet mixes a 0-based pixel (`pixel i covers [i, i+1) Box/N_lens`) with a 1-based slab (`slab s covers [(s-1), s)`) in one sentence; the module docstring (:13-14) and the sidecar use 1-based `[(i-1), i)` for pixels. | lensplanes.py:13-14 vs README:131-132. | Write "0-based numpy index i" or use the 1-based form in both places. |
| n-2 | Nit | python/README.md:108 | The `read_pm_header` row says "the 529-byte big-endian header" although the reader also accepts little-endian files (snapshot.py:123-128) and the conventions bullet (:141-142) says so. | snapshot.py:116-133. | "529-byte header (byte order detected)". |
| n-3 | Nit | deck :457 vs python/glamio/powerdm.py:3-6 | The deck says `PowerDM.R` has "2 header lines (a, step, z, N_g, k_Ny)"; the file has three header lines (45-char run header, epoch line, caption) as the docstring and the real files show. | `head -3 PowerDM.R.0070.0001.dat`. | Deck: "3 header lines (run header, epoch line, caption)". |
| n-4 | Nit | deck :257 vs case12 `fR_halo.*.meta` | "Every sidecar records step, a, z, a_x, a_v" — the fR_halo sidecar carries neither `a_x` nor `a_v` (it has `a`, `z`, `step`, `a_old`), so `FRHalo.a_v`/`a_x` are `None`. The reader is right; the deck sentence (or the writer) is not. | `F.a_v is None`, `F.a_x is None` on case12; PMP2mod_extras.f90:645-658 writes no `a_v`. | Deck: "every sidecar of a product that carries velocities records a_v"; or add `a_x`/`a_v` to `fRHaloClose`. |
| n-5 | Nit | deck (13 × `h^{-1}Mpc`, 0 × `Mpc/h`) vs README/docstrings/sidecars (`Mpc/h`, `kpc/h`, `Msun/h`) | Notation differs between the deck and the package; same units. | `grep -c` on the deck. | Pick one notation per document and note the equivalence once in the README. |
| n-6 | Nit | test/precise_z/README.md:64; test/extras/README.md:15-16 | The new `glamio` row sits in the middle of the case table (between case14 and case4g) and is not a case; the extras README still says the precise_z suite has "180 checks" while the commit message and plan say 221. | Diff of test/precise_z/README.md; extras README unchanged by the commit. | Move the row to the end (or a "cross-checks" paragraph) and update the extras README count. |
| n-7 | Nit | .gitignore:92-98; python/tests/test_readers.py:35 | `python/tmp` (the `GLAMIO_TEST_TMPDIR` used in the review recipe) is not ignored; fixtures are removed on success, but a failing test leaves `glamio-readers-*` directories as untracked clutter. `__pycache__/`, `.pytest_cache/`, `python/glamio.egg-info/` and the work dirs are ignored correctly. | `git check-ignore python/tmp` prints nothing; the others resolve. | Add `python/tmp/` to .gitignore. |
| n-8 | Nit | plan :1056 | "`ProjectionXY/XZ` (legacy product...)" — the run writes three projections (XY, XZ, YZ). | case10 `Run1/ProjectionYZ.*.dat`. | "`ProjectionXY/XZ/YZ`". |
| n-9 | Nit | python/glamio/schedule.py:50-72 | Column properties are `steps`, `classes` (plural) but `a`, `z`, `chi`, `exact` (singular); `__getitem__` accepts both names and row indices while `PowerDM.__getitem__` accepts names only (powerdm.py:45). | Code. | Either all plural or document that `steps`/`classes` are plural to avoid the `step` property name and the `class` keyword. |

Counts: Blocker 0, Major 1, Minor 11, Nit 9.

## 3. Verified OK

- README install section: `pip install -e python/` from the repository root
  matches pyproject (`setuptools>=64`, `dynamic version = glamio.__version__`,
  `packages.find include = ["glamio*"]`); "Only numpy>=1.24 and Python >= 3.10"
  matches `dependencies`/`requires-python`; `glamio.__version__` is `0.2.0`
  and `test_catalogue.py` now accepts any `x.y.z`.
- Readers table: every function name, return type and listed attribute/method
  exists with that name (`read_schedule` → `Schedule.rows/omega_m/omega_l/row/primary/stations`;
  `read_lensplanes` → `planes[plane, j, i]`, `plane`, `slabs`, `axis_means`,
  `pixel_edges`, `slab_edges`; `read_modes` → `data`, `mode`, `wavenumbers`,
  `k_magnitude`, `power`; `read_subsample` → `q[6, n_sel]`, `positions`,
  `velocities`, `ids`, `lattice`; `read_frhalo` → `data` with `fR_centre`,
  `fR_mean`, `catalogue()`; `read_powerdm` → `k`, `nmodes`, `pk`, `a`, `z`,
  `step`, `ngrid`, `k_nyq`, `binning`, `pk_deconvolved`; `read_pm_header` →
  the listed fields, `box`, `a_v`, `velocity_scale`; `read_snapshot` →
  `iter_pages`, `read`, `units` values).
- Usage example, adapted paths: runs on case10 (`plane("x", 1)` returns a
  `(72, 72)` memmap, `axis_means()` = [1, 1, 1] to 5e-10) and on case11 for
  schedule, modes, subsample and snapshot; `iter_pages(units="physical")`
  yields one page of 46656 particles with x in [0.0013, 149.9993) and
  |v| ≤ 578 km/s; the checkpoint `PMcrd.DAT` reads too (step 100, a_v
  0.67646).
- "`power()` equals PowerDM.R bin by bin (raw CIC)": on case11 step 70
  (n = 11, 72³ mesh) max |pk/PowerDM − 1| = 4.1e-6 over bins 1..11 (the
  printed precision of `es13.5`), `nmodes` identical, `k` equal to 5e-6; bins
  1..n are complete because any mode with a component |i| = n+1 has
  nint|m| ≥ n+1. `nbins` defaults to `n`.
- Conventions bullets vs code, sidecars and deck: mode normalisation and cube
  (`fourier_norm`, `modes`, `layout` lines of the real sidecars; deck :261-266);
  `n = floor(k_c Box/2π)` gives 11 (k_c = 0.5) and 4 (k_c = 0.185) on case11;
  "iz > 0 modes count twice" matches `power()` and the suite's
  `binned_power`; subsample decode (`x_decode`, `v_decode`, `id_rule`, `SALT`
  of the real sidecar; deck :374-380); snapshot units `[1, N_grid+1)`,
  `(X−1) Box/N_grid`, `a_v = AEXPN − ASTEP/2` (PMP2mod_subsample.f90:106,
  :165; the reader's `h.a_v` equals the subsample's `a_v` to 1e-7 on case11
  and the suite ties them to the velocity quanta, 0.125 km/s against a 1.4 %
  a vs a_v difference); PowerDM window formula identical to the deck's :459.
- Module docstrings vs writers: `schedule.py` (comment lines, `exact`, `class`
  semantics, written at start-up and restart: PMP2mod_extras.f90:53-81);
  `frhalo.py` (pair order, `fR_mean` definition, `a_old = a − da`,
  `fR_bg(a)`, `catalogue`/`catalogue_bin` lines: PMP2mod_extras.f90:559-658
  and the case12 sidecar; `catalogue()` resolves `CATALOGS/...` against the
  run directory and returns the 5-row catalogue); `powerdm.py` (3 header
  lines, `N_grid/2` rows, `kNyq = π N_grid/Box` = 1.508 → "1.51", `log`
  skips empty bins: real files); `modes.py` momentum definition and
  `N_dep = N_g/2` default vs plan §7 and deck :344-346; `subsample.py`
  `dx_quant` 15.6 kpc/h at Box = 1024, clipping range.
- Property names map to the right keys: `n_fr` → `n`, `fR_bg_old` →
  `fR_bg(a_old)`, `d_code` → `D_code(a)`, `k_fund` = 2π/Box, `classes` →
  `class`, `a_x` → `a_x`; `plane(axis, slab)` is documented 1-based in the
  README, the docstring and the sidecar, raising `IndexError` for 0.
- Section 16 entry: six `glamio` checks (asserts.py labels), 8 tests in
  `test_readers.py` + 7 in `test_catalogue.py` = 15, SplitMix64 reference
  vectors (lines 184-185) and a little-endian snapshot fixture (line 331),
  529-byte record (files are 537 bytes), "pilot's 2048² planes read
  memory-mapped with axis means 1" (|mean − 1| ≤ 1.1e-11, 201 MB memmap),
  "snapshot particles = subsample particles within the quanta" (case11 has
  M = 1, so the subsample is the full snapshot), `Catalogue` behaviour
  unchanged by the `MetaProduct` refactor (same `a`, `z`, `a_v`, `box`,
  `step` semantics; it only gains properties). "221/221" is left to Review 4.
- Plan §9 sizes vs the reader constants and real files: 12 × 4096² × 4 B =
  805 MB (README/docstring), 12 × 2048² × 4 B = 201 326 592 B (pilot file),
  deltam k_c = 1.5 at L = 1024 → n = 244, 0.47 GB; subsample M = 100 at 2048³
  → 1.03 GB; snapshot 24 B × 2048³ = 206 GB.
- Plan §3/§4/§6 consistency: velocity-epoch text (legacy catalogue 1/a, new
  products a_v, `velocities_at_av()` multiplies by a/a_v) matches
  catalogue.py and the README's "Velocity epoch" section; §6 "readers live in
  python/glamio, halocat should import glamio" matches the entry.
- Deck definition frames vs docstrings (where not listed above): lens-plane
  layout and pixel convention (:285-287), subsample ID/hash/quantisation
  (:371-380), fR_halo definitions (:383-392), catalogue columns and units
  (:440-446) agree with lensplanes.py, subsample.py, frhalo.py, catalogue.py.
- Hygiene: `git ls-files python/` holds only sources, README, pyproject and
  tests; no untracked files under python/, docs/reviews or test/precise_z in
  the main checkout; `__pycache__/`, `.pytest_cache/`, `python/glamio.egg-info/`,
  `test/precise_z/work/`, `test/extras/work/` are ignored; the 16 files of the
  commit are all legitimate.
- test/extras/check_*.py: independent decoders of `deltaL`/`deltam` cubes and
  the subsample (same reshape by necessity); they are Fortran-unit-test
  oracles and should keep their own code.

## 4. Open questions for the author

1. M-1: do you want `axis_means()` to stay the only built-in check, or add
   `plane_means()` (per plane, informational) so users do not re-derive the
   wrong expectation? Also confirm the Fortran's logged per-plane means were
   never meant as an assertion.
2. m-3: is the deck's `δ̃_m/W_CIC` a leftover of the original V2 wording, or
   was a deconvolved comparison run at some point? The as-run V2 and
   `power()` are both raw.
3. m-5/m-6: which of the halocat needs (redshift → step, product paths per
   step, ID → lattice / Lagrangian q) should land in glamio before the
   halocat session, and which can halocat implement on top?
4. m-7: GLAM's `2.774e11` or the physical `2.7754e11` for
   `PMHeader.particle_mass()`? If the latter stays, should the catalogue
   sidecar's `particle_mass` be the documented source of truth?
5. m-9/m-10: should check_run.py become a glamio consumer (and be listed in
   the plan), and should a `read_timing` reader be planned since check_run.py
   already parses `timing.log`?
6. `read_modes` on a real `*_deltaL` file and `read_powerdm` on a real
   `PowerDM.Z` file are untested here (no such files among the allowed paths);
   Review 4 may cover the former from `test/extras/run_deltaL_test.sh`.
7. m-11: please re-run the editable install after the merge so `pip show`
   reports 0.2.0.
