# Independent review of the 20260908 BDM v3 convergence campaign

Reviewer session: 9 September 2026. Read-only review of the frozen campaign.
No production source, analysis code, receipt, catalogue, snapshot, density tape,
plot or slide was modified. Nothing was merged, pushed, committed or deleted.

## 0. Pins, actual state and divergence

| Item | Pinned in the prompt | Actually observed |
|---|---|---|
| Branch | `validation/bdm-convergence-20260908` | same |
| Campaign/result tip | `cb87021d140fd2ea72841b6cea32abf090c255dd` | present; **HEAD is `d485bc7cc9c43842c6bc7f685b5be0c0213d6ad5`**, one commit ahead |
| Base (`cz`) | `e289c3f7ea28d38390a1663e50bb2784b85abc84` | `cz` resolves to exactly this; it is the merge-base with HEAD |
| Working tree | — | **no modified tracked files**; 24 untracked paths, all unrelated to this campaign |

The single commit beyond the pin, `d485bc7 "Document independent BDM convergence
review prompts"`, adds only `BDM-refine/validation/20260908-convergence/REVIEW-PROMPTS.md`.
That is the "later addition of review instructions" the prompt excludes, so I
reviewed `e289c3f..cb87021` (18 commits, 129 files, 118 126 insertions, **zero
deletions or modifications outside `2LPTIC_Gui/README.md`**). I did not reset or
switch the shared checkout.

**Production finder is unchanged, verified two ways.** `git diff --name-status
e289c3f cb87021` touches no `*.f90` outside the new campaign directory, and
`sha256(PMP2linker.f90)` in the working tree equals `sha256` of the same file at
`e289c3f` and equals `common.py:15 FINDER_SHA` and `executables.json
production_finder_sha256`: `39f7e514…e4cd5d`. The six new Fortran files are all
campaign-only adapters under `.../20260908-convergence/{ic,replay}/`.

---

## 1. Findings, ranked

No P0 finding. I found **no defect that invalidates a stated numerical result**:
every measurement I recomputed from the retained raw data reproduced exactly.
The findings below are, in order, one interpretation gap that materially affects
how the headline should be used, two methodological/provenance defects, and a set
of smaller quantified corrections and latent hazards.

### P1-1 — The "passing mass intervals" are not property-general: median axis ratios differ by up to 5.9 % inside them, and this magnitude is never reported

*Status: confirmed, quantified. Disclosed qualitatively; magnitude never stated.*

* **Where.** `assess.py:9-11` (`CRITERIA`) and `assess.py:33-45` (`evaluate`).
  The working screen tests abundance, abundance jackknife sigma, median bound
  mass, median resolved `Vmax` and reference completeness. It does **not** test
  `axis_ba`, `axis_ca`, `aperture_total_mass`, `aperture_radius`,
  `bulk_velocity_km_s` or `centre_distance_mpc_h`, all of which
  `analyze.py:153-160` measures and stores.
* **Trigger.** Read `README.md` ("At z=0 and a 300-particle floor, passing
  log10(Mbound) intervals are C/E 12.75–14.75, E/F 14.00–14.75, F/T 12.50–14.75")
  or `CONVERGENCE.md`'s table and treat the interval as "the catalogues agree
  there".
* **Expected vs observed.** A reader would expect the unscreened properties to
  agree at a comparable level. Restricting to the bins that actually pass:

  | Pair, z | Passing interval | max abs median ΔM | max abs median ΔVmax | **max abs median Δ(b/a)** | **max abs median Δ(c/a)** |
  |---|---|---:|---:|---:|---:|
  | C/E z=0 | 12.75–14.75 | 0.63 % | 0.29 % | 1.31 % | 1.59 % |
  | C/E z=1 | 12.75–14.25 | 0.80 % | 0.18 % | 1.67 % | 1.19 % |
  | C/E z=2 | 12.75–13.75 | 0.57 % | 0.21 % | 1.25 % | 1.94 % |
  | **E/F z=0** | 14.00–14.75 | 0.54 % | 1.30 % | **4.11 %** | **5.86 %** |
  | **E/F z=1** | 12.75–14.25 | 2.08 % | 1.82 % | 2.11 % | **5.45 %** |
  | E/F z=2 | 12.50–13.25 | 1.85 % | 0.91 % | 1.30 % | 1.67 % |
  | **F/T z=0** | 12.50–14.75 | 1.09 % | 1.53 % | 2.92 % | **4.91 %** |
  | **F/T z=1** | 12.50–13.50 | 2.55 % | 1.98 % | 1.94 % | **4.70 %** |
  | **F/T z=2** | 12.50–13.00 | 2.49 % | 1.99 % | 2.53 % | **5.42 %** |

  Applying the report's own 5 % mass tolerance to the reported median `c/a`
  would remove the E/F z=0 pass at 14.00–14.25 and 14.25–14.50, and the F/T z=2
  pass at 12.50–12.75. The particle comparison C/E is the only one whose
  shapes agree at the ≤2 % level. Note also that `Axba`/`Axca` are not raw
  inertia-tensor ratios: `PMP2linker.f90:1105-1117` applies an empirical
  exponent correction whose exponent depends on `RadRms/aperture`, so a
  resolution change enters the reported shape nonlinearly. `convergence.json`'s
  `property_definitions` documents `bound_mass`, `aperture_total_mass`,
  `aperture_radius` and `vmax`, but not the axis ratios.
* **Evidence.** `review_results.json` →
  `unscreened_properties_inside_passing_intervals`, recomputed by me from
  `convergence.json` and independently reproduced from the retained membership
  index and match arrays (`review_checks.py recompute`, 0 disagreements). The
  effect is visible in figure page 15 (z=0 force, Δ(c/a) panel), where E/F sits
  outside the grey ±5 % band across the whole displayed range — so this is a
  narrative omission, not a plotting error.
* **Affected products.** `README.md` "Completed results" section,
  `CONVERGENCE.md` table and "Limits" list, slides 10–11.
* **Already disclosed?** Partly. `CONVERGENCE.md` says "Shape and velocity
  differences are reported separately; the displayed mass/Vmax screen does not
  certify every halo property", and slide 5 says "shapes can require stricter
  cuts". Neither gives the size, and neither warns that the shape disagreement
  exceeds the mass tolerance inside intervals declared to pass.
* **Minimal repair.** Add one sentence and one column to `CONVERGENCE.md` /
  `README.md` stating the maximum absolute median `b/a` and `c/a` shift inside
  each passing interval (numbers already in `convergence.json`; no rerun), and
  add `axis_ba`/`axis_ca` to `analyze.py`'s `property_definitions` noting the
  empirical correction. Optionally add a shape criterion to `CRITERIA` and
  report the resulting (shorter) intervals alongside the current ones.

### P2-1 — The paired eight-octant jackknife is not a usable uncertainty, and the criterion built on it can neither fail on perfect agreement nor discriminate at the high-mass end

*Status: confirmed, quantified. Correlation caveat disclosed; the failure modes are not.*

* **Where.** `analyze.py:144-146` (delete-one-octant paired ratio, correct
  `sqrt(7/8·Σ)` normalisation) and `assess.py:40`
  (`abundance_uncertainty = 100·sigma <= 5.0`).
* **Trigger 1 — vacuous pass.** When the two catalogues have identical octant
  counts in a bin, every leave-one-out ratio is exactly 1, the sample variance
  is 0 and `sigma = 0`. This happens in **3 eligible bins** (C/E z=1,
  14.00–14.25, at all three particle floors; octant counts
  `[10,9,7,5,2,6,4,4]` in both C and E). A 47-vs-47 count comparison is then
  reported with zero uncertainty, where a counting error alone is ≈15 %.
* **Trigger 2 — no power where it matters.** At floor 300, **21 bins pass the
  5 % abundance criterion with a jackknife sigma above 2.5 %**, i.e. above half
  the tolerance. Examples (shift ± sigma): E/F z=0 14.50–14.75 = 0.00 ± 4.16 %;
  A/E z=0 14.50–14.75 = 0.00 ± 7.15 %; B/C z=2 13.50–13.75 = +4.35 ± 6.57 %.
  These bins are consistent with a true shift anywhere up to the tolerance.
  This is precisely the high-mass region where E/F and C/D "pass".
* **Expected vs observed.** The screen reads as "abundance agrees to 5 % and the
  measurement is good to 5 %". Observed: the estimator measures how *consistently*
  the ratio differs between octants of one realisation, not the uncertainty of
  the ratio; it is exactly zero on perfect agreement and largest where counts are
  smallest.
* **Evidence.** `review_results.json` →
  `degenerate_jackknife_in_eligible_bins`, `abundance_passes_with_sigma_over_2p5pct`.
  I reproduced `jk` for all 63 comparisons from the stored octant histograms to
  1e-10 relative.
* **Mitigating fact.** The *median bound mass* and *median Vmax* criteria are
  paired statistics over 54–265 matched objects in the same bins and are far
  tighter (e.g. E/F z=0 14.00–14.25: ΔM = −0.24 %, ΔVmax = −1.30 %). The
  headline conclusion "force resolution limits lower-mass z=0 results" is driven
  by those, and survives. It is the *abundance* leg of the screen at high mass
  that carries little information.
* **Affected products.** `CONVERGENCE.md` screen description,
  `convergence-assessment.json:criteria`, slide 10.
* **Minimal repair.** State in `CONVERGENCE.md` that `sigma` is a paired
  octant-consistency estimate that can be exactly zero and is not a confidence
  interval, and report alongside it a simple counting uncertainty
  `sqrt(1/N_left + 1/N_right)` per bin so the reader can see where the abundance
  test has no power. No rerun needed; counts are in `convergence.json`.

### P2-2 — `ic-production-validation.json` was not produced by the committed `ic_validate.py`; three of its validity conditions have never been evaluated on the production ICs

*Status: confirmed defect (provenance/reproducibility). Partly disclosed.*

* **Where.** `ic_validate.py:66-81` versus the retained
  `ic-production-validation.json`.
* **Trigger.** Compare the keys the committed code writes with the keys on disk.
  The receipt contains `E_F_float32_position_bound_mpc_h`, which the committed
  code **never writes**, and is **missing nine keys the committed code always
  writes**: `E_F_pure_float32_position_bound_mpc_h`,
  `E_F_native_periodic_mapper_bound_mpc_h`, `E_F_position_component_rms_mpc_h`,
  `E_F_components_exceeding_pure_roundoff_bound`, `E_F_nonedge_excess_components`,
  `E_F_excess_without_exact_clamp_sentinel`, `max_inferred_displacement_mpc_h`,
  `coordinate_caveat`, `boundary_source`.
  Timestamps agree: the receipt is `2026-09-08T00:20:01Z`; `ic_validate.py` was
  last modified an hour later and committed once, in `c3d7dd7` at 01:20 local.
* **Expected vs observed.** A frozen campaign receipt should be reproducible by
  the committed producer. Observed: the retained receipt is from a superseded
  version, and the receipt records **no source hash**, so nothing in the frozen
  tree ties it to any code. Consequently the committed script's conditions
  `nonlocal_excess == 0`, `unexplained_excess == 0` and `max_displacement < 256`
  (`ic_validate.py:66`) have never run on the production ICs, and neither has the
  RMS position difference.
* **Impact on conclusions.** Low. The retained receipt reports
  `E_F_max_position_difference_mpc_h = 6.103515625e-05`, exactly equal to the
  pure float32 bound `2·spacing(float32(4097))·256/4096`, which means **no row
  exceeded pure round-off** and the differential periodic edge clamp had no
  measurable effect. The conservative source-derived bound of
  `1.8310546875e-4 Mpc/h` quoted in `README.md` was therefore never approached.
  `README.md` already states that the retained check "did not record edge-clamp
  incidence or RMS differences" — the undisclosed part is that the *committed
  code would now record them and adds three assertions that were never exercised*.
* **Affected products.** `ic-production-validation.json`, the `README.md`
  paragraph on initial coordinates, slide 4.
* **Minimal repair.** Re-run `micromamba run -n cosemu python3 -B ic_validate.py`
  and commit the refreshed receipt, adding a `source_sha256` field. Cost: reads
  the three 1024³ IC particle sets (~77 GB), one core, a few GiB; submit to
  `cosma8-serial` with `--cpus-per-task=1 --mem=8G`, ~15 min wall, ≈0.25
  core-hours. This resolves a material provenance gap for the single most
  important assumption of the E/F comparison.

### P2-3 — The frozen membership checker's host-exclusion test is one-sided and, as instrumented, cannot exercise its own violation branch; four real centre-in-aperture pairs exist that it cannot reach

*Status: confirmed blind spot, quantified. No scientific impact. "Never exercised" is disclosed; the one-sidedness and the four instances are new.*

* **Where.** `BDM-refine/validation/20260906-n1024/run_validation.py:498-511`.
  Candidates come from `query_ball_point(centre[i], props[i,8])`, so every `j`
  returned already satisfies `dist <= R_i`; the loop then skips unless `i` is the
  **more massive** halo, increments `examined`, and flags a violation when
  `dist² < R_i²`. Every pair that reaches `examined += 1` is therefore already a
  violation, so `examined == 0` *is* the pass condition rather than evidence of
  coverage, and the discriminating branch is never taken.
* **Trigger.** Test the symmetric condition `dist < max(R_i, R_j)` over all pairs
  instead of the priority-ordered one.
* **Expected vs observed.** `CONVERGENCE.md` reports "Host violations 0" for all
  21 catalogues, which reads as "no published halo centre lies inside another
  halo". My all-pairs test over all 21 catalogues finds **4 pairs** where it does:
  A/z1 (2), F/z1 (1), C/z0 (1). In every one, the **more massive** halo's centre
  lies inside the **less massive** halo's reported aperture — the exact
  orientation the checker's priority filter excludes. There are zero cases of the
  orientation the checker does test.
* **Quantitative evidence** (`review_results.json` → `host_exclusion_disputed_pairs`).
  Inverting the aperture formula `PMP2linker.f90:1003`
  (`aperture = rso + Cell·min(Rext/(rso/Cell)^SlopeR, 0.75)`, `Cell = 0.125`,
  `Rext = 0.15`, `SlopeR = 0.20`):

  | Catalogue | M_more / M_less | separation | / larger aperture | / larger unextended SO |
  |---|---:|---:|---:|---:|
  | A/z1 | 1.000 | 0.36835 | 0.986 | **1.028** |
  | A/z1 | 1.115 | 0.45126 | 0.981 | **1.013** |
  | F/z1 | 1.114 | 0.69849 | 0.990 | **1.009** |
  | C/z0 | 1.029 | 0.37486 | 0.998 | **1.040** |

  All four separations **exceed** the unextended SO radius of both members. These
  are `Rext`-aperture artefacts, not SO host-exclusion failures. Four pairs out
  of 377 475 published rows; no measurable effect on any comparison.
* **Affected products.** `CONVERGENCE.md` "Membership checks" table caption,
  `final-review.json:counts.host_exclusion_violations`, slide 6.
* **Minimal repair.** In the checker, query with `max(R_i, R_j)` (or query at the
  catalogue's maximum radius and test both orientations) and report the count of
  centre-in-aperture and centre-in-SO pairs separately; state in the table
  caption that exclusion is enforced on the SO radius, not the reported aperture.
  Note the checker is *frozen* (`common.py:16 CHECKER_SHA`) — changing it
  invalidates every replay receipt, so this belongs in a follow-up campaign, not
  a patch to this one.

### P3-1 — Published memberships are not disjoint; "zero exact duplicate member sets" does not mean exclusive membership

*Status: confirmed, new measurement. Not a defect; a missing quantification the report invites.*

Every one of the 21 catalogues contains particles belonging to more than one
published halo. Excess memberships (occurrences beyond the first), from a
complete pass over all 1 129 million memberships:

| Catalogue | memberships | shared IDs | rate |
|---|---:|---:|---:|
| A/z0 | 4 019 463 | 2 565 | 0.064 % |
| C/z0 | 32 586 813 | 27 457 | 0.084 % |
| E/z0 | 261 452 588 | 327 722 | **0.125 %** |
| F/z0 | 273 883 096 | 238 127 | **0.087 %** |
| T/z0 | 274 844 546 | 200 152 | 0.073 % |
| … (all 21 in `review_results.json`) | | | 0.06 – 0.13 % |

This is intrinsic to a spherical-overdensity finder with **centre** exclusion:
apertures may overlap in volume (the closest E/z0 pair is separated by
0.519 × the sum of radii) and a particle in the lens can be bound to both. The
rate differs between the two runs being compared for force resolution
(E 0.125 % vs F 0.087 %, a 44 % relative difference), but the absolute effect on
bound mass is ≤0.04 %, two orders below the 5 % tolerance, so **no conclusion
changes**. Repair: one sentence in `CONVERGENCE.md` stating that membership is
centre-exclusive but not particle-exclusive, with the measured rate.

### P3-2 — The "300-particle floor" does not bind for E/F or F/T; the effective floor there is ≈1867 particles

`analyze.py:161`: `mass_floor = max(2.5e12, floor · max(m_p))`. For the 1024³
runs `m_p = 1.3389e9`, so `2.5e12` corresponds to **1867 particles** and the
nominal floor never binds: the E/F and F/T rows of `convergence-assessment.json`
are *identical* at floors 100, 300 and 1000. The lowest eligible E/F bin
(12.50–12.75) contains only haloes with ≥2362 particles. The floor binds only
for pairs containing A (29 particles at `MassMin`) or B/C/D (233). Describing
the E/F and F/T results as "at a 300-particle floor" understates the effective
requirement ~6×. This *strengthens* the scientific claim — E/F's 9.2 % abundance
deficit at 12.50–12.75 is measured on >2000-particle haloes, so it is force
resolution, not particle discreteness — but the README wording obscures it.
Repair: quote `mass_floor_msun_h` and the implied particle count per pair.

### P3-3 — `evolve-E` peaked at 95.99 GiB of a 96 GiB request (99.99 %); no simulation-rerun memory guidance is given

`accounting.json`, job `evolve-E-t64`: `ReqTRES mem=100663296K` (96 GiB),
`MaxRSS = 95.99 GiB`. The README gives rerun headroom for replays (≥256 GiB) and
analysis (≥16 GiB) but nothing for the simulations. A rerun of E at the recorded
allocation is at serious OOM risk. Repair: record ≥128 GiB for an E-equivalent
1024³/2048³ run (D/F/T already had 384 GiB against a 304.7 GiB peak).

### P3-4 — A clean checkout cannot reproduce the launch provenance or the presentation

`.gitignore` in the campaign directory excludes `work/`, `launches-*.tar.gz`,
`work-artifacts.tar.gz`, `work-controls.tar.gz`, `root-compiler-artifacts.tar.gz`,
`figures/` and `slides/`. The five archives and the LaTeX sources therefore exist
only in this working directory. I verified all five archives independently
(recomputed archive sha256, extracted and re-hashed **every** manifest entry:
92 + 27 + 329 + 113 + 23 = 584 entries, 0 missing, 0 hash mismatches, 0 originals
left behind; the 47 entries that are symlinks are present as symlink members).
Provenance is therefore *sound where the archives exist*, and `jobs.json` records
`script_sha256`, `bundle_sha256`, `driver_sha256` and `submit_command` for all 41
jobs — but a fresh clone has the hashes and not the content. This is disclosed
("presentation and figures remain local artifacts excluded from Git"). Repair, if
the branch is to be self-contained: commit the two small `launches-*.tar.gz`
(3.9 MB total) or record the exact `submit.py` argument vector per job so an
equivalent script can be regenerated.

### P3-5 — Smaller provenance and code observations

* `convergence-assessment.json` records `input_sha256` of `convergence.json`
  (verified: `31d96e95…`) but **no `assess.py` source hash**, unlike
  `analyze.py`'s `analysis_source_sha256` (verified: `c1fdbb62…` equals the
  on-disk file). `final-review.json:artifact_sha256.assess.py` supplies the link
  indirectly. Repair: add `assess_source_sha256` in `assess.py:75`.
* `E-ic.json` and `F-ic.json` lack `sampled_modes_sha256`, present in A, B, C, D
  and T — they predate `campaign.py:161-163`, which added both the field and the
  `assert sha(matched_modes.bin) == sha(master/...)` guard. Nothing consumes the
  field, and `ic-production-validation.json` records all seven hashes as
  identical (`40f4c0e4…`), so this is receipt completeness only.
* `analyze.py:91` and `:94` increment `spatial_candidates` and `examined` in the
  same branch, so the two receipt fields are always equal (e.g. A/C z=2: both
  7046). The receipt suggests two independent diagnostics but reports one.
* `common.py:50` stages every JSON write at a fixed `<name>.tmp`. Two concurrent
  writers of the same path would corrupt each other. Not triggered here (paths
  are per-tag and jobs are dependency-gated), but latent.
* `campaign.py:45` silently leaves a line unchanged when a key in `changes` is
  absent from the reference `Init.dat`. I verified all 12 keys are present in the
  generated files, so nothing was dropped; a future reference-input change could
  silently disable e.g. `MG_flag` with no error. Repair: assert every key in
  `changes` was consumed.
* Absolute mass scale: the finder's `MassOne` uses ρ_crit = 2.774e11 (verified
  against the tape header to 7.8e-9 relative), while GLAM's generated `Setup.dat`
  prints "Particle Mass" from 2.77467e11 — a 2.4e-4 relative offset. Inherited,
  common to all seven runs, cancels in every ratio, shifts absolute log10 M by
  1e-4 dex. No action.

### Unresolved hypotheses (stated, not established)

* Whether a **fixed 2048³ finder mesh suppresses sensitivity** in all seven
  comparisons is not settled by anything in this campaign. It is plausible: the
  finder's peak detection, cell size (0.125 Mpc/h) and search scales are held
  fixed while the evolution mesh varies, so any halo whose identification is
  mesh-limited would be limited identically in every run. Two facts bound it
  loosely — the measured matched centre distances are sub-cell (E/F 0.016–0.046
  Mpc/h; F/T 0.003–0.029) and completeness is ≥0.97 everywhere — but neither
  tests the mesh. The report says so explicitly. Resolving it needs the
  finder-mesh replay described in §5.
* Whether the **empirical `Axba`/`Axca` correction** amplifies or damps the true
  resolution dependence of halo shape. The correction exponent depends on
  `RadRms/aperture`, which itself changes with resolution. Untested.

---

## 2. What I independently tested, what I relied on, what I did not test

### Independently tested (my own code, not the campaign's)

| Control | Result |
|---|---|
| Re-derived all **63 assessments** from `convergence.json` with independently written criteria, without calling `assess.evaluate` | **0 mismatches** in eligibility, all five per-bin checks, verdicts and merged intervals |
| Rebuilt **all 63 comparisons** from the retained membership indices and match arrays only — histograms, octants, abundance ratios, jackknife, mass floors, valid bins, completeness, left matched fraction, and the 8×3 quantile statistics | **0 disagreements** (rtol 1e-9) |
| Re-hashed **all 21 membership index `.npz`** | 21/21 match receipts |
| Re-hashed **all 21 raw membership tapes** (`repair-members.bin`, 4 MB – 2.09 GB, ~30 GB total) | 21/21 match receipts |
| Re-hashed **all 144 science files < 1 GiB** across the 21 v3 receipts | 144/144 match; 24 files ≥1 GiB (density tapes, PMcrs) size-verified only |
| Parsed the E/z0 raw tape from scratch: header, 403 probed rows against the index, sorted/in-range IDs, `props[6] == count·m_p` | all consistent; `m_p` matches `2.774e11·Ω_m·(L/N)³` to 7.8e-9 |
| **Complete** membership disjointness over all 21 catalogues (1.129 × 10⁹ memberships) | 0.06–0.13 % shared IDs → **P3-1** |
| **Symmetric all-pairs host exclusion** over all 21 catalogues | 4 centre-in-aperture pairs → **P2-3** |
| `project_ids` against explicit lattice construction for ratios 2, 4, 2 and 4 (incl. non-power-of-two 12/3) | bijective onto the coarse lattice; fine index `r·i` ↔ coarse `i` at identical Lagrangian position |
| `analyze.match` against **exhaustive** all-pairs matching on 4 fixtures: clean nested, periodic wrap + split object + disjoint pair, two coarse haloes competing for one fine halo, empty reference catalogue | identical matches in all four |
| **Wider-radius rematch on real data**: E/F z=0, 64³ Mpc/h subvolume, candidate radius ×6 → 9525 pairs examined vs 420 within the production filter | **identical 411 matches**, none changed, none lost; max matched centre distance 0.66 Mpc/h vs a 1.58 Mpc/h filter radius |
| Slurm accounting recomputed from a fresh `sacct` query | billed 2133.4772 ch and consumed 1891.5669 CPU-h, both matching `accounting.json`/README; no parent/step double counting; all 41 jobs on `cosma8-serial` |
| All 18 `plot-controls.json` axis-coverage values recomputed from `convergence.json` | exact match; every displayed upper limit covers the highest displayed 84th percentile |
| All five cleanup archives: archive sha256 + every manifest member extracted and re-hashed | 584/584 entries accounted for, 0 mismatches, 0 originals surviving |
| `final-review.json` `artifact_sha256` (13 entries) and `simulation_and_pair_receipt_sha256` (28 entries) against current files | 41/41 match |
| Visual inspection of figure pages 15–18 and slides 1–2, plus full text of all 29 slides | see §4 |
| Cross-check of the seven `Init.dat` and `Setup.dat` files actually used | differ **only** in `Nrow`, `Ngrid` and (T) `step da`; everything else byte-identical |

### Relied on prior evidence (not re-derived)

* The 21 density tapes (32 GiB each, 672 GiB) — existence and size only. I relied
  on `replay-*.json` `density.sha256` and on `shared_bit_identical_density`, which
  I confirmed is `true` for **all 21** pairs, for the claim that legacy and v3
  consumed bit-identical density fields.
* The PM snapshot files (`PMcrs*.DAT`, ~26 GB per 1024³ epoch) — I relied on the
  `*-simulation.json` output manifests and on `replays.py:65-68`, which verifies
  them before each replay.
* `ic-production-validation.json`'s all-row E/F and F/T comparison — I read and
  checked the code and reproduced its arithmetic (the `F_T_staggered_velocity_ratio`
  of 1.015502 matches `((a_T)/(a_F))^1.5` = 1.015506 from the recorded steps), but
  I did not re-read the 77 GB of IC particle files. See **P2-2**.
* The campaign's own fixture controls (`analysis-controls.json`,
  `driver-controls.json`, `replay-resume-controls.json`,
  `launch-cleanup-preflight.json`, `replay/preflight*.json`, `ic/checks-*.json`,
  `schedule-preflight.json`). I read their scope and confirmed they cover the
  resume/failure-injection, corrupt-archive, active-job-guard and interrupted-
  deletion cases in Stage-2 items 10–11; I did not re-run them. In particular
  `schedule-preflight.json` records a real 16³/32³ paired run in which the
  table-driven schedule reproduced the native schedule **bit-identically** at all
  three output epochs — a strong control I did not attempt to reproduce.
* Compiler flags and build linkage (`build.json`, `ic-build.json`,
  `replay-build.json`, `executables.json`): I verified the recorded
  `production_finder_sha256` and the binary inventory, and that `replays.py:38`
  binds each stage receipt to `sha(PMP2replay.{variant}.exe)`; I did not rebuild
  anything or re-verify object/module hashes.

### Not tested

* No simulation, replay, finder run or IC generation was executed. The
  legacy-vs-v3 finder difference itself was not analysed (see §5 Q4).
* Absolute physical accuracy against any external reference (theory, another
  code, another box). Out of scope by construction.
* Finder-mesh sensitivity, box size, cosmology, 2LPTIC ICs — none is measured by
  this suite, and the report says so.
* The 29-slide deck was inspected in full as text and on 2 rendered pages; I
  rendered and examined 4 of the 18 figure pages. The campaign's own
  `visual-review.json` inspected 6 slide pages and 4 figure pages of the final
  PDFs and states honestly that its inspection is "representative, not a claim to
  have inspected every page".

---

## 3. Claim / evidence / verdict

| Claim | Evidence I checked | Verdict |
|---|---|---|
| Production `PMP2linker.f90` unchanged from the pinned `cz` base | `git diff` name-status over the full base→target range; sha256 of the working tree file == sha at `e289c3f` == `FINDER_SHA` == `executables.json` | **Confirmed** |
| Seven simulations at the stated (N_p, N_g), GR, L=256, z_init=100, exact outputs z=2,1,0, common 2048³ analysis mesh, `iVirial=1`, `Rext=0.15`, `MassMin=2.5e12` | Seven generated `Init.dat`/`Setup.dat` differ only in `Nrow`/`Ngrid`/(T)`step da`; snapshot headers give a = 1/3, 1/2, 1 exactly and the correct N/G/particle counts; `BDM.config` frozen by `common.py:20-23` and re-checked by `replays.py:34` | **Confirmed** |
| All 7 simulations, 21 paired replays, 21 validated v3 catalogues, 63 assessments | 21 `replay-*.json` all `pair_completed: true`; 21 `v3-validation-*.json` all `completed`; 63 comparisons in `convergence.json`; 377 475 v3 rows | **Confirmed** |
| z=0, 300-particle floor: C/E 12.75–14.75, E/F 14.00–14.75, F/T 12.50–14.75 | Re-derived independently from `convergence.json` and from raw arrays | **Confirmed** (see P1-1, P3-2 for scope) |
| E/F lowest eligible z=0 bin: −9.2 % haloes, −7.3 % median mass and Vmax | −9.186 %, −7.256 %, −7.288 % in bin 12.50–12.75 | **Confirmed** |
| C/E max absolute z=0 bin-median shifts 0.63 % (mass), 0.29 % (Vmax) | 0.6259 %, 0.2862 % over its 8 eligible bins | **Confirmed** |
| "Force resolution limits lower-mass z=0 results; higher-redshift timestep sensitivity remains" | E/F z=0 fails below 14.00 on a smooth, monotone, high-significance trend (−9.19 %, 12.0–22.2σ in the lowest bins); F/T passes only 12.50–13.00 at z=2 and 12.50–13.50 at z=1; C/D (512³ particles) fails below 14.00 at z=0 too, independently corroborating a particle-count-independent force-mesh limit | **Confirmed** |
| Matched IC modes shared across resolutions; not merely the same seed | `ic/generate.py:45-51` fixes a 1024-draw master stride and per-plane luxury seeds; `:56-68` reuses a single master ALPHA; `:70` fixes the physical origin at `Box/4096`; `ic-production-validation.json` records identical 33³ mode samples for **all seven** and one shared ALPHA; `Setup.dat` `DRho/rho` rises with N_p (0.0387 → 0.0488 → 0.0579) as expected for added short modes at fixed normalisation | **Confirmed** |
| Nested initial-lattice IDs give an exact Lagrangian correspondence | `PMP2start.f90:397,399`: `Q1 = (ip−1)·N_g/N_row + 1`, `Icurrent = ip + (jp−1)N + (kp−1)N²`, so physical Lagrangian position is `(ip−1)·Box/N_row` (corner lattice, not cell-centred) and `analyze.py:56-63` maps fine `r·i ↔ coarse i` exactly; verified numerically for 4 ratios | **Confirmed** |
| 158/316 steps, normal endpoints and output epochs preserved; T starts with correctly staggered velocities | `timetable.py:110-111` asserts endpoint and flag identity; `timetable.json` records outputs 92/109/158 → 184/218/316; `T-simulation.json` snapshots at 184/218/316 with a = 1/3, 1/2, 1; `ReadCampaignSchedule` (`timetable.py:75`) refuses an IC whose `ASTEP` disagrees with `dAlist(1)`, and `T`'s `Setup.dat` carries the exact half first step | **Confirmed** |
| F/T initial positions identical; E/F physical velocities identical; E/F max position difference 6.1035e−5 Mpc/h | Recorded over all 1024³ rows; equals exactly the pure float32 bound, so the edge clamp contributed nothing | **Confirmed**, but see **P2-2** on the producing code |
| Matching is mutual-best shared-lattice overlap with ≥50 % in each object; not nearest-centre | `analyze.py:96-107`; verified against exhaustive matching on 4 fixtures and against a 6× wider real-data search | **Confirmed** |
| Bound mass is float64 `count × m_p`, distinct from the `Rext` aperture | `analyze.py:71-72,153-154`; `replay/membership.f90:43` writes `Mvir` and `Mtotal` as separate fields; checker asserts `Mvir ≤ Mtotal` | **Confirmed** |
| Legacy and v3 consumed bit-identical density fields | `shared_bit_identical_density: true` in **all 21** pair receipts; enforced by `replays.py:157-159` | **Confirmed** (relied on recorded hashes) |
| 377 475 rows, 0 exact duplicate member sets / repeated IDs / mass-count mismatches / host violations | Structurally true; but the four "0" fields are **hard-coded constants** reached only if the preceding asserts pass (`run_validation.py:521-522`), so they are pass/fail flags, not counts. My independent tests confirm zero exact duplicates and zero SO-level exclusion failures, but find non-disjoint memberships (**P3-1**) and 4 aperture-level centre-in-radius pairs (**P2-3**) | **Confirmed as stated; incomplete as a characterisation** |
| 2133.48 billed core-hours, 1891.57 consumed CPU-hours, estimate ≈2400, ceiling 3888.75 | Independently recomputed from `sacct`: 2133.4772 and 1891.5669. 38 COMPLETED, 1 CANCELLED, 2 FAILED — all three non-completions are superseded 1–5 s render pilots totalling 0.001 core-hours, each with a recorded reason. 95.9 % of billed core-hours ran at ≥75 % CPU efficiency; only 0.9 % below 25 % | **Confirmed** |
| Render receipt "predates" the visual review — stale or contradictory? | `render-validation.json` completed 00:12:36Z saying "newly completed results still need visual scientific review"; `visual-review.json` at 00:19:58Z binds the **same** PDF hashes (`58ebad41…`, `5c56ab1e…`) that the current files still have. Sequence is correct; the earlier note is accurate at its own timestamp, not stale | **Not a contradiction** |
| Repaired absolute velocity/centre-distance axes | `plots/plot_convergence.py:93-94` autoscales **then** sets only `bottom=0`; the 18 recorded coverage values reproduce exactly; visually, figure page 16 shows \|Δv\| on 0–30 km/s and \|Δx\| on 0–0.33 Mpc/h, not clipped at 1 | **Confirmed** |
| Archives are verified and restorable | All five re-verified member by member by me | **Confirmed** (but gitignored — **P3-4**) |
| "These are conditional catalogue-stability measurements at a fixed finder mesh; they do not demonstrate finder-mesh convergence, convergence of every halo property, or absolute accuracy" | Correct, and the strongest thing the suite supports | **Confirmed** |

---

## 4. Presentation and figures

`bdm_convergence_n300.pdf` (18 pages, `5c56ab1e…`) and `bdm_convergence.pdf`
(29 pages, `58ebad41…`) match the hashes bound by `visual-review.json`,
`plot-manifest-n300.json` and `presentation-validation.json`, and those three
receipts' own hashes match `final-review.json`. Pages 15–18 render correctly:
axes labelled with units, grey reference bands explicitly captioned "reference
scales, not accuracy guarantees", 16–84 percentile shading visible, jackknife
error bars on Δn, omitted bins genuinely absent rather than plotted at zero,
completeness and left-matched-fraction panels carry the 0.9 threshold line, and
the repaired \|Δv\| / \|Δx\| axes scale to the data.

Scientific wording is careful and I found nothing overclaimed. Slide 2 separates
correctness from convergence; slide 4 states the first-order IC exception; slide
10's caption says the criteria "do not certify absolute accuracy or every
reported property"; slide 12 quantifies the z=2 timestep failure; slide 13 lists
one-realisation, finder-mesh and 2LPTIC limits explicitly. The one gap is the
one in **P1-1**: the deck shows the Δ(c/a) panels but never states that inside
the intervals it declares passing, the median `c/a` disagreement reaches 5.9 %.

Reproducibility: a clean checkout **plus the retained `work/` tree** reproduces
the measurements (`analyze.py`, `assess.py`, `common.py`, the 21 v3 receipts and
`convergence.json` are all committed, and `analyze.py`'s cached-match receipts
bind left/right membership hashes, row counts, box, analysis mesh and the
analysis source hash, so changed inputs are rejected). A clean checkout **alone**
cannot reproduce either the measurements (the `work/` data is gitignored, which
is correct at 672 GiB + 26 GB/epoch) or the presentation (`slides/` is
gitignored) — see **P3-4**.

---

## 5. Answers

### Q1. Is this branch suitable to merge into `cz`, and what must be fixed first?

**Yes, with two documentation fixes and one cheap re-run first.** The branch adds
only a self-contained campaign directory plus `AGENTS.md`; it modifies no
production source, and every measurement I could recompute reproduced exactly.
The code quality is high: frozen source/binary/config identity is enforced at
every stage, resume is guarded against stale reuse, and the failure-injection
controls are real.

Before merging:

1. **P1-1** — add the maximum median `b/a` and `c/a` shift inside each passing
   interval to `README.md` and `CONVERGENCE.md`, and define the axis-ratio
   columns in `analyze.py`'s `property_definitions`. Zero compute cost; the
   numbers are already in `convergence.json`.
2. **P2-2** — re-run `ic_validate.py` and commit the refreshed
   `ic-production-validation.json` with a source hash (≈0.25 core-hours on
   `cosma8-serial`, 1 core, 8 GiB). Without this the branch carries a receipt no
   committed code produces, for the campaign's most load-bearing assumption.
3. **P2-1** and **P3-2** — one paragraph each in `CONVERGENCE.md` on what the
   octant jackknife is and is not, and on the effective particle floor per pair.

Not merge blockers, but worth doing in the same pass: **P2-3** and **P3-1**
wording in the membership-checks section, **P3-3** memory guidance, and the
`assess.py` source hash (**P3-5**).

### Q2. Precisely which properties, mass/redshift ranges and settings support saying the refined BDM catalogues passed a convergence test?

Only this, and only with all qualifiers attached:

> For a single 256 Mpc/h GR realisation with **first-order Zel'dovich ICs at
> z_init = 100**, with the **BDM v3 finder run at a fixed 2048³ analysis mesh**,
> `iVirial = 1`, `Rext = 0.15`, `MassMin = 2.5 × 10¹² M⊙/h`, the **abundance,
> median bound mass (within 5 %), median resolved Vmax (within 2 %) and
> reference completeness (≥90 %)** of matched haloes are stable to those
> tolerances over the following ranges in log10(M_bound/[M⊙/h]), each with the
> other two numerical parameters held fixed:
>
> | Refinement | Held fixed | z=2 | z=1 | z=0 |
> |---|---|---|---|---|
> | Particles 512³→1024³ (C/E) | 2048³ force mesh | 12.75–13.75 | 12.75–14.25 | 12.75–14.75 |
> | Force mesh 2048³→4096³ (E/F) | 1024³ particles | 12.50–13.25 | 12.75–14.25 | 14.00–14.75 |
> | Timestep halved (F/T) | 1024³ particles, 4096³ mesh | 12.50–13.00 | 12.50–13.50 | 12.50–14.75 |
>
> Every bin in these ranges contains ≥30 objects per required statistic and lies
> entirely above `max(2.5e12, floor·m_p)`; for E/F and F/T that is ≥1867
> particles, not 300.

Explicitly **not** covered by that statement: halo shape (`b/a`, `c/a` differ by
up to 5.9 % inside the E/F and F/T ranges — **P1-1**), bulk velocity (which also
carries the native half-step output staggering), the tails and halo-to-halo
scatter (the 16–84 spread is many times the median shift for every property),
absolute physical accuracy, and any combination of settings not in the table.

### Q3. Which stronger claims remain unsupported, and what is the smallest additional measurement for each?

| Unsupported claim | Smallest measurement |
|---|---|
| "The catalogues are converged for log10 M ≥ 14.0 at z=0" (joint over all three refinements) | The intersection is untested at any single joint point. The cheapest genuine test is **replay-only**: none exists — it needs one new simulation at 1024³/8192³ (or 512³ particles with halved steps at 4096³) to check that the refinements are additive. C/D vs E/F already shows the force-mesh limit is roughly independent of particle count at z=0 (both fail below 14.00), which is supporting but not decisive. **Cost: one D-class run ≈500–900 core-hours.** |
| "The halo finder is converged" | Finder-mesh sensitivity is completely untested — **the single largest gap**. It is **replay-only**: re-run the existing `replay.exe` on the retained F z=0 and E z=0 snapshots at analysis meshes 1024³ and 4096³ and compare to the 2048³ catalogues. No new simulation. The adapter already supports any power-of-two target. **Cost: ≈6 replays × ~20 core-hours ≈ 120 core-hours**; the 4096³ analysis mesh needs ≥400 GiB (a 4096³ float32 field alone is 256 GiB). |
| "Halo shapes/velocities are converged" | No new run needed. Re-assess the existing `convergence.json` with a shape criterion and report the resulting intervals (**P1-1**). To go further, a 1000-particle-floor shape analysis on the existing F/T and E/F matches — also free. |
| "The v3 catalogues are converged with the production 2LPTIC ICs" | A 2LPTIC control pair at reduced cost: C-class (512³/2048³) and E-class (1024³/2048³) with `lpt_order = 2`, same box/seed, comparing the C/E particle-refinement interval to the first-order result. **Cost: ≈100 core-hours** (C ≈ 40, E ≈ 64). |
| "The measured intervals apply to other volumes or cosmologies" | Untestable without a second realisation; at minimum a second seed at C-class to separate cosmic variance from numerical shift. **Cost: ≈50 core-hours.** |
| "Abundance agrees to 5 % with 5 % uncertainty at the high-mass end" | Free: report counting uncertainties alongside the octant jackknife (**P2-1**). |
| "Published memberships are exclusive" | Already measured by this review: they are not (0.06–0.13 %). |

### Q4. Does this campaign establish any convergence advantage of v3 over legacy BDM?

**No.** This is unambiguous from the code, not an inference:

* `analyze.py:47` asserts `'v3' in report['variants']` and reads
  `report['variants']['v3']['membership']` only. The string `legacy` does not
  appear in `analyze.py`, `assess.py`, `convergence.json` or
  `convergence-assessment.json`.
* The 21 legacy catalogues **were** produced and are retained (385 631 published
  rows against v3's 377 475, i.e. legacy publishes 2.2 % more haloes), on
  bit-identical density fields, but they were never matched, binned or assessed.

The campaign therefore establishes **conditional stability of the v3 catalogues**
under particle, force and timestep refinement — nothing about whether v3
converges *better* than legacy. That question is answerable **replay-free**:
every input already exists on disk, and running the existing `analyze.py`
pipeline against `variants['legacy']` would produce a directly comparable set of
63 assessments for a few core-hours of analysis. The earlier correctness repairs
(unique hosts, consistent membership, well-defined discrete SO) are a separate
and independently established result; nothing here adds to or subtracts from
them. The presentation is correct to keep the two questions apart (slide 2).

---

## 6. Files created, and cleanup

Created, all under `BDM-refine/analysis/review-20260909-convergence/` (new
directory; `review-20260906-merged`, `full-audit-20260906` and
`duplicate-hosts-20260905` already existed and were not touched):

| File | Size | Purpose |
|---|---:|---|
| `REVIEW.md` | this file | the review |
| `review_checks.py` | 28 KB | six independent controls: `assess`, `recompute`, `tape`, `disjoint`, `hashes`, `matchfixture`, `matchwiden` |
| `review_results.json` | 32 KB | machine-readable results of all of them |

Reproduce any control with
`micromamba run -n cosemu python3 -B review_checks.py <stage>` from that
directory. `review_checks.py` writes only into its own directory; it opens the
campaign's data read-only, and `python3 -B` prevents `__pycache__` creation when
it imports `analyze.py`.

Cleanup performed: the temporary `pdftoppm` previews (6 PNGs) and two
intermediate JSON files were folded into `review_results.json` and removed; all
scratch lived in the session scratchpad, not in the repository. No campaign file
was created, modified or deleted — confirmed by `git status --porcelain` on the
campaign path (empty) and by a modification-time scan of the campaign directory.
The archives, figures, slides, LaTeX sources and all untracked working
directories are untouched and remain uncommitted. No live or queued Slurm jobs
exist for this account, so no running work was disturbed. This review applied no
repairs, performed no merge and changed nothing in the frozen campaign.
