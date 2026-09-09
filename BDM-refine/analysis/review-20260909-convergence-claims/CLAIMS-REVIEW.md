# What the 20260908 BDM campaign does and does not establish about convergence

Focused scientific review (Prompt 2 of `REVIEW-PROMPTS.md`), 9 September 2026,
against the overridden target below. Read-only: no production source, analysis
code, receipt, catalogue, snapshot, density tape, figure or slide was modified;
nothing was merged, committed, pushed or launched.

## 0. Pins, and a disclosure about independence

| Item | Requested | Observed |
|---|---|---|
| Branch | `fix/bdm-convergence-review-20260909` | same, and is HEAD |
| Target | `ceb7e8c` | `ceb7e8c182f9f0395d1876ebab15f1f2be6724e6`, is HEAD |
| Base | `e289c3f…` | `cz` resolves to exactly this; it is the merge-base |
| Working tree | — | no modified tracked files |

`e289c3f..ceb7e8c` is 22 commits. The first 19 are the original campaign
(`..d485bc7`); the last three are the repair response. **Prompt 2 asks the
reviewer not to read another new review's verdict first. I cannot satisfy that:
I wrote the Prompt-1 review, and `9625a51` committed it into this branch.** This
report is therefore a second, differently-aimed pass by the same reviewer, not a
blind second opinion. I have tried to compensate by (a) attacking the conclusion
with measurements neither the campaign nor my earlier review made, and (b)
auditing my own earlier report as adversarially as the campaign's, which turned
up five errors of my own (§1.2).

**Preservation verified.** Across `cb87021..ceb7e8c`, `convergence.json`,
`analyze.py`, `common.py`, `campaign.py`, `replays.py`, `timetable.py`,
`build.py`, `submit.py`, `cleanup.py`, `archive_scratch.py`, `preflight.py`,
`driver_checks.py`, `plots/plot_convergence.py` and all 56
`*-ic.json` / `*-simulation.json` / `replay-*.json` / `v3-validation-*.json`
receipts are **byte-identical**. `sha256(PMP2linker.f90)` still equals its value
at `e289c3f`. I re-derived all 63 acceptance decisions from `convergence.json`
with my own criteria implementation against the *new* assessment: **0
mismatches**, `criteria` unchanged. The repair changed documentation,
diagnostics and provenance only.

---

## 1. Confirmed errors

### 1.1 In the campaign as it now stands

**C-1 (interpretation, material). The broad z=1 passing intervals for the two
force-mesh comparisons are partly a zero crossing, not convergence.**

The matched median bound-mass shift for force refinement *changes sign* between
z=2 and z=0 in every bin where all three epochs are measurable:

| Pair | log10 M | ΔM at z=2 | ΔM at z=1 | ΔM at z=0 | z=1 suppression |
|---|---:|---:|---:|---:|---:|
| E/F | 12.75 | +0.69 % | −1.50 % | −5.53 % | 0.5× |
| E/F | 13.00 | +1.85 % | **−0.17 %** | −4.07 % | **11×** |
| E/F | 13.25 | +2.25 % | +0.90 % | −2.83 % | 2.5× |
| E/F | 13.50 | +2.70 % | +1.32 % | −1.88 % | 1.4× |
| C/D | 12.75 | +1.42 % | −0.89 % | −5.18 % | 1.6× |
| C/D | 13.00 | +2.36 % | **+0.08 %** | −3.81 % | **31×** |
| C/D | 13.25 | +2.64 % | +1.00 % | −2.55 % | 2.6× |
| C/D | 13.50 | +2.99 % | +1.62 % | −1.72 % | 1.1× |

All eight bins **pass** the screen at z=1. Two of them (E/F and C/D at
13.00–13.25) agree to better than 0.2 % at z=1 while the same haloes disagree by
1.9–2.4 % at z=2 and 3.8–4.1 % at z=0 — suppression by factors of 11 and 31
relative to their own neighbours in redshift. E/F's z=1 passing interval
(12.75–14.25) is consequently much wider than its z=0 interval (14.00–14.75) not
because the simulation is better converged at z=1, but because the force-resolution
error passes through zero near that epoch.

Nothing in `CONVERGENCE.md`, `README.md`, `REVIEW-RESPONSE.md` or the deck says
this. `CONVERGENCE.md` says only "The z=0 timestep result does not extend to all
higher-redshift masses; inspect the separate z=1 and z=2 ranges", which points at
F/T (where there is *no* sign change) and not at the force comparisons (where
there is). Reporting the widest of three intervals, or interpolating between
epochs, would be wrong here.

*Evidence:* `claims_results.json` → `z1_sign_change_cancellation`,
`redshift_behaviour`. *Repair:* one row in the "Limits" list stating that the
force-mesh error changes sign between z=2 and z=0, that the z=1 intervals are
therefore not a conservative summary, and that z=0 is the binding epoch.

**C-2 (presentation, minor). `CONVERGENCE.md`'s shape paragraph quotes a
non-maximal case.** Its summary sentence cites "the coarser A/C particle
comparison reaches 6.24 % at z=0", but the table immediately above it contains
A/E z=1 at **8.92 %** and A/C z=1 at **8.11 %**. `REVIEW-RESPONSE.md` does quote
the 8.92 % figure, so this is a mismatch between the response and the committed
document, not a suppressed number. *Repair:* quote the table maximum.

### 1.2 In my own Prompt-1 review (`BDM-refine/analysis/review-20260909-convergence/REVIEW.md`)

The response challenges four of my points. Three challenges are correct, one is
correct in substance but understates what remains. I also found a fifth error
myself.

**E-1 (my error, confirmed). "21 bins pass the 5 % abundance criterion with a
jackknife sigma above 2.5 %" conflated one condition with the whole screen.**
21 bins meet the abundance-*difference* condition with σ > 2.5 %; only **15**
pass the full screen. My worked example, B/C z=2 at 13.50–13.75 (+4.35 ± 6.57 %),
*fails* the σ ≤ 5 % condition and is not a passing bin at all. The response's
correction is right. **The substance survives with the corrected count**: 15
full-screen passes have σ > 2.5 %, including the three bins that carry E/F's and
F/T's z=0 high-mass intervals — E/F 14.25–14.50 (−1.59 ± 2.97 %), E/F
14.50–14.75 (0.00 ± 4.16 %), F/T 14.50–14.75 (+1.89 ± 4.73 %). The abundance leg
of the screen still has little discriminating power exactly where those
intervals end.

**E-2 (my error, confirmed). My host-exclusion framing was wrong.**
`PMP2linker.f90:694-760` (`RemoveDuplicates`) implements

```fortran
if (Mvir(ip) < Mvir(jp) .or. (Mvir(ip) == Mvir(jp) .and. jp < ip)) then
   ...
   if (dd < dble(Rvir(jp))**2) removeHost(ip) = .true.
```

with `Rvir(ip)` assigned the **Rext-extended** aperture at line 1014. That is a
priority-ordered, one-directional test on the extended aperture, and
`run_validation.py:498-511` implements the same orientation, the same radius and
the same tie-break. My claim that the checker has a "blind spot" and that a
symmetric `max(R_i, R_j)` test is the repair was wrong: a symmetric test would
check a *different* rule. My four reverse-orientation pairs are **allowed** by
the implemented rule, which is what the response says. They remain a useful
characterisation — and my measurement that all four separations exceed both
unextended SO radii still shows the catalogues are cleaner than the rule
requires — but they are not a defect. Withdraw the "minimal repair".

**E-3 (my error, confirmed). "Every pair that reaches `examined += 1` is
therefore already a violation" is false at exact equality.** The test is
`dd < R**2`, so a pair at exactly `dist == R_i` increments `examined` without a
violation. The response's `exact_aperture_boundary` control demonstrates it
(`examined: 1, violations: []`). My conclusion that `examined == 0` was the
production pass condition is still right for the 21 catalogues, and the campaign
had already disclosed the coverage limitation.

**E-4 (my error, confirmed). My Q4 cost estimate for a legacy comparison was
wrong.** I wrote that "running the existing `analyze.py` pipeline against
`variants['legacy']` would produce a directly comparable set of 63 assessments
for a few core-hours". I checked all 21 paired receipts: the `legacy` variant
carries only `catalogue`, `elapsed_seconds`, `maxrss_kib`, `stage_receipt` —
**no membership tape or index anywhere**. The legacy replay binary is the exact
`a8c7715` finder with no membership hook, so matched-property, completeness and
shared-tracer statistics are impossible without new legacy replays. The response
is right. What *is* free is catalogue-level statistics, and I have now done that
(§2, N-4).

**E-5 (my error, found here). My suggested joint-refinement run is
infeasible.** I proposed "one new simulation at 1024³/8192³ … ≈500–900
core-hours". An 8192³ float32 field is 2 TiB, which does not fit a 1 TB COSMA8
node at all, let alone alongside 1024³ particles. The response flags this. §5
replaces it with a feasible design.

---

## 2. What the experiment establishes — claims, measurements, limits

Six measurements below (N-1 … N-6) are new to this review; they use only
retained data and no new simulation, replay or matching.

| # | Claim | Supporting measurement | Precise limit |
|---|---|---|---|
| 1 | Particle refinement 512³→1024³ leaves abundance, median bound mass, median resolved Vmax and completeness stable to the screen over 12.75–14.75 at z=0 (12.75–14.25 at z=1, 12.75–13.75 at z=2) | 63 assessments, reproduced twice independently from raw arrays; C/E max median shifts 0.63 % (mass), 0.29 % (Vmax) over 8 eligible z=0 bins | At a fixed 2048³ finder mesh, first-order ICs, one realisation, whole bins above 3.213e12 M⊙/h. Says nothing about shape, tails, scatter or absolute accuracy |
| 2 | Force refinement 2048³→4096³ is **not** converged below ~10^14 M⊙/h at z=0 | E/F: −9.19 % abundance, −7.26 % median mass, −7.29 % median Vmax in 12.50–12.75, at 12–22σ of the paired octant scatter; corroborated by C/D | This is the campaign's most robust negative result. The lowest E/F bin contains only haloes with ≥2362 bound particles, so it is force resolution, not discreteness |
| 3 | Timestep halving is converged at z=0 over 12.50–14.75 but not at z=2 | F/T z=0 max median shifts 1.09 % (mass), 1.53 % (Vmax); z=2 passes only 12.50–13.00 | The F/T error decays **monotonically and without sign change** toward z=0 (N-1), so the z=0 result is genuine convergence, unlike the force case |
| **N-1** | **The force-mesh error changes sign with redshift; z=1 intervals are inflated by cancellation** | §1.1 table: E/F and C/D shifts go +2.7 %→+1.3 %→−1.9 % from z=2 to z=0; 8 passing z=1 bins straddle the zero crossing, two suppressed 11× and 31× | Demonstrated for matched median bound mass at 12.75–13.50. z=0 is the binding epoch; no interpolation in z is safe |
| **N-2** | **The refinements are approximately separable, so the intersection reading is better supported than either the campaign or my first review claimed** | (a) Force refinement measured at two particle loads agrees closely: C/D (512³) vs E/F (1024³) matched median bound mass differ by **≤0.35 pp** in every z=0 bin (≤0.72 pp at z=2), Vmax by ≤0.27 pp. (b) The particle chain composes: measured Δ(A/E) equals Δ(A/C)⊕Δ(C/E) to a median \|residual\| of **0.16 pp**, max 0.94 pp, over 18 cases | Tests force×particle and particle×particle only. Force×time and particle×time are untested. Both tests are at the same fixed analysis mesh and ICs |
| **N-3** | **The fixed 2048³ finder mesh does not suppress sensitivity** — the "shared finder limitation" hypothesis is largely refuted | At identical analysis mesh and z=0, candidate counts are 7,969,396 (A, 256³) vs 958,777 (C, 512³) vs 736,241 (E, 1024³) — a factor of ~10. At fixed particle load, refining the evolution mesh *raises* candidates by 21–30 % (C→D +21 %, E→F +23 %, B→C +30 %). F vs T differ by 0.005 % | The mesh does set a **common absolute scale**: the median published halo aperture is only 2.97–3.04 analysis cells at z=0 in *all seven* runs. Any absolute error at that scale is shared and cancels in every ratio, and remains unmeasured |
| **N-4** | **No convergence advantage of v3 over legacy exists in abundance** (stronger than "is not established") | Legacy abundance convergence computed from the retained `catalogues-pair.npz` for all 21 catalogues, zero new compute. 114 comparable bins: v3 has the smaller absolute shift in 68, legacy in 46; median \|shift\| 2.78 % (v3) vs 3.21 % (legacy); mean 4.12 % vs 4.39 %. For the well-sampled pairs (C/E, E/F, F/T, C/D, B/C) legacy−v3 is ≤0.8 pp in almost every bin | Abundance only. The 68/46 split is ≈2σ under an independence assumption that the correlated bins do not satisfy, so it is not significant. Method validated: v3 counts from the ASCII catalogue reproduce the membership-derived counts with **zero** bin disagreement across all 21 |
| **N-5** | **Median agreement coexists with large halo-to-halo scatter; the screen constrains neither tails nor per-object accuracy** | Inside the passing intervals, the 16–84 half-width of the per-halo bound-mass difference is **9–25× the median shift** for the particle comparisons (C/E z=0: ±7.5 % scatter around a 0.63 % median) and ≈1× for F/T | A downstream use needing per-halo masses, mass-function width, or halo–galaxy assignment is not covered at all by the published intervals |
| **N-6** | **Match-selection bias is real but bounded** | Inside the bins the screen uses at z=0, unmatched eligible reference haloes are 1.47 % (C/E), 0.77 % (E/F), 0.14 % (F/T), 1.57 % (A/E). The lost objects are systematically lower-Vmax: E/F 257.7 vs 302.6 km/s (−15 %), C/E 311.7 vs 347.6 (−10 %) | Matching does preferentially keep the better-resolved objects, but at ≤1.6 % of the sample it cannot move a bin median by more than a few tenths of a percent. Independently, a 6× wider candidate search on E/F z=0 returned identical matches, so the loss comes from the 50 % mutual-overlap rule, not the spatial filter |
| **N-7** | The E/F difference is a force-solver effect, not chaotic amplification of its IC seed | E/F starts from an initial-position component RMS of 2.38e-7 Mpc/h (9.5e-7 of the 0.25 Mpc/h interparticle spacing) with **identical** velocities, and ends at 7.26 % median mass difference. F/T starts from **bit-identical positions** and a 1.55 % coherent velocity difference and ends at ≤1.09 % | No pair in the suite differs by numerical noise alone, so the irreducible noise floor of the whole pipeline is unmeasured |

### Answers to the specific traps Prompt 2 lists

* **Shared finder limitation suppressing sensitivity** — largely refuted (N-3);
  but the common ~3-cell aperture scale is a shared absolute limitation that no
  ratio in this suite can see.
* **Intersection of separate passing ranges** — better supported than claimed
  (N-2), and still not a proof: force×time and particle×time are untested.
* **Selection bias from projection / spatial candidates / mutual-best** —
  quantified and bounded at ≤1.6 % (N-6); the spatial filter loses nothing.
* **Unresolved Vmax, publication limits, bin migration, both-halo cuts** —
  unresolved Vmax is negligible (≤3 objects in any bin, over all 63
  comparisons); the whole-bin rule removes straddling bins; the binding cut for
  E/F and F/T is the 2.5e12 publication mass, i.e. ≥1868 particles, not the
  nominal floor.
* **Medians vs tails, shape and velocity** — quantified (N-5); the repair now
  publishes the shape maxima, which reach 8.92 % in median c/a (A/E z=1).
* **Eight correlated octants; sparse bins and NaNs** — NaNs are propagated
  honestly and make a bin ineligible; the paired jackknife is exactly zero in 3
  eligible bins and exceeds half the tolerance in 15 passing bins.
* **Redshift dependence and cancellation** — a demonstrated cancellation (N-1).
  No generalisation across epochs is justified.
* **Common first-order ICs and native staggering** — every run shares the same
  master-normalised modes by construction, so a common IC bias would cancel in
  all seven comparisons and is invisible here. This is the one listed trap that
  the suite genuinely cannot address; it needs the 2LPTIC control.
* **Duplicate sets vs physical uniqueness** — different questions, as suspected:
  member sets are unique, but 0.06–0.13 % of memberships are shared between
  haloes.
* **v3 vs legacy** — now measured for abundance (N-4): indistinguishable.

---

## 3. Assessment of `REVIEW-RESPONSE.md`

**Overall: justified.** Every factual and operational claim I could check is
true, the numerical results are untouched, and the three challenges it raises
against my Prompt-1 review are correct (§1.2). It does not overstate the
repairs, and it declines to change the acceptance decisions — correctly, since
none of my findings implied a decision was wrong.

| Response item | Verdict | Evidence |
|---|---|---|
| "All 63 original acceptance decisions and criteria are unchanged" | **Verified** | Re-derived independently against the new assessment: 0 mismatches; `criteria` identical; `input_sha256` still equals `sha256(convergence.json)` |
| P1-1: shape maxima added for all 21 pair/redshift cases | **Verified** | I recomputed all 21 rows of the new table from `convergence.json`: **0 disagreements**. Table maximum is A/E z=1 at 8.92 %, as the response states |
| P2-1: jackknife scope clarified; independent-count scale added as a labelled benchmark | **Justified, and its statistical criticism of me is correct** | `sqrt(1/N_L+1/N_R)` does drop the covariance of two catalogues sharing a realisation; 20.63 % for 47+47 is not this ratio's uncertainty. Retaining it as a labelled scale rather than an estimator is the right call |
| P2-1: "the review's 21 bins … some fail other conditions" | **Correct; my error** | 21 meet the abundance condition, 15 pass the full screen (§1.2 E-1) |
| P2-2: IC validator rerun, source-bound | **Verified** | Job 11960736 COMPLETED, 604 s, 0.16778 billed core-hours, MaxRSS 15.99 GiB, `cosma8-serial`. New receipt carries `source_sha256`, `support_source_sha256` and all seven `input_receipt_sha256`. All previously missing fields present; `E_F_components_exceeding_pure_roundoff_bound = 0` confirms the edge clamp contributed nothing |
| P2-3: host rule is priority-ordered extended-aperture; symmetric test would test another rule | **Correct; my error** | `PMP2linker.f90:745-750` + line 1014; checker matches orientation, radius and tie-break (§1.2 E-2) |
| P2-3: violation branch *can* be exercised; six controls run the frozen block | **Verified, and genuinely rigorous** | `review_followup_checks.py:57-87` extracts the checker's host block by AST from `run_validation.py` and `exec`s it verbatim. Recorded outcomes are exactly as claimed, including `exact_aperture_boundary` with `examined: 1, violations: []` |
| P3-1/P3-2/P3-3/P3-4/P3-5 dispositions | **Verified** | Effective-cut table reproduces (1868 particles at the publication cut, 2362 at the first whole bin); `assess_source_sha256` present and correct; memory guidance added |
| "Cannot obtain legacy matched properties by swapping `variants['v3']`" | **Correct; my error** | No legacy receipt contains a `membership` key (§1.2 E-4) |
| "An 8192³ float32 field alone is 2 TiB" | **Correct; my error** | 8192³ × 4 B = 2.2e12 B |
| Cost/cleanup claims (2133.645 total, 31 slides, archive verified) | **Verified** | Independent `sacct`: 2133.64500 billed / 1891.66047 CPU-h, all 41 jobs `cosma8-serial`. Slides 31 pages, no warnings, PDF hashes bound by `visual-review.json` and matching on disk. New 6-member archive re-verified member-by-member, 0 originals surviving; all 13 `final-review.json` artifact hashes match |

**Where the response is incomplete, not wrong:**

1. It does not address the z=1 cancellation (N-1) — neither did I; it is new here
   and is the most consequential remaining interpretation gap.
2. Its own committed `CONVERGENCE.md` quotes a non-maximal shape figure while
   the response text quotes the maximum (§1.1 C-2).
3. "No superiority in numerical convergence over legacy is established by the
   current v3-only analysis" is true but weaker than the data allow: for
   abundance, none *exists* (N-4), and that measurement was free.
4. The IC recheck again ran at 99.9 % of its memory request (15.99 of 16 GiB).
   The response notices and recommends 24 GiB, which is the right disposition,
   but this is the second job in the campaign to land within 0.01 GiB of its
   ceiling (after `evolve-E` at 95.99 of 96 GiB).

**One disagreement.** The response says the reverse-orientation pairs "do not
… require disjoint memberships". Agreed — but the two are separate points, and
the 0.06–0.13 % shared-membership rate does matter for one specific downstream
use: summing bound masses over a catalogue double-counts by that amount. The
new `README.md` wording ("not a bound on individual-halo errors") is correct but
could add that catalogue-summed mass is not a partition of the particle set.

---

## 4. A defensible paragraph

> Using seven matched GR simulations of a 256 h⁻¹ Mpc box sharing
> master-normalised first-order initial conditions at z_init = 100, and running
> the refined (v3) BDM finder on a common 2048³ analysis mesh for every run, we
> measure how the published halo catalogues respond to particle number
> (256³–1024³), force mesh (1024³–4096³) and timestep (158 vs 316 steps).
> Against a descriptive working screen — abundance and median bound mass within
> 5 %, median resolved V_max within 2 %, reference completeness ≥ 90 %, ≥ 30
> objects per statistic — the catalogues are stable at z = 0 over
> log10(M_bound/[M⊙ h⁻¹]) = 12.75–14.75 for particle refinement, 12.50–14.75 for
> timestep refinement, and only 14.00–14.75 for force refinement; the 4096³
> force mesh is the limiting tested setting, with 9.2 % fewer haloes and 7.3 %
> lower matched median mass and V_max at 12.50–12.75 on the 2048³ mesh even
> though those haloes contain more than 2000 bound particles. The force-mesh
> error changes sign between z = 2 and z = 0, so the wider z = 1 intervals
> reflect a zero crossing rather than convergence and z = 0 is the binding
> epoch; the timestep error, by contrast, decays monotonically to below 1.1 % by
> z = 0. Refinements are approximately separable: the force-mesh response
> measured at 512³ and at 1024³ particles agrees to better than 0.35 percentage
> points at z = 0, and the particle-refinement chain composes to within 0.2
> percentage points. These are conditional statements about catalogue stability
> at a fixed finder mesh: they do not certify absolute accuracy, and they do not
> extend to halo shape (median c/a differs by up to 5.9 % inside the passing
> force and timestep intervals, and by 8.9 % for the coarsest particle pair), to
> halo-to-halo scatter (the 16–84 spread of the per-halo mass difference is 9–25
> times the median shift for particle refinement), to a different finder mesh,
> to other volumes or cosmologies, or to a 2LPTIC suite. Comparing the retained
> legacy catalogues over the same 114 mass bins, the refined finder shows no
> measurable abundance-convergence advantage; the v3 repairs improved
> correctness, not resolution response.

---

## 5. Prioritised minimal follow-up

Each entry names the uncertainty it resolves. Costs assume `cosma8-serial`,
explicit `--cpus-per-task` and `--mem`, never a whole-node allocation.
Nothing here is a recommendation to run without the stated question.

**Tier 1 — replay-only, no new simulation.**

1. **Finder-mesh sensitivity, minimal form: 2 replays.** Re-run the existing
   `replay.exe` on the retained E and F z=0 snapshots at analysis mesh **1024³**
   only, and compare the E/F force-resolution signal at 12.50–12.75 with the
   2048³ result. *Resolves:* whether the shared analysis mesh is setting the
   force-resolution signal. N-3 already shows the mesh does not suppress
   sensitivity; this tests the complementary question of whether it biases the
   measured amplitude. Coarsening is cheap — a 1024³ float32 field is 4 GiB —
   and the retained 2048³ results supply the control. *Estimate:* ~2 × 15
   core-hours at 32 cores; request ≥128 GiB from the measured 2048³ peak of
   191.99 GiB scaled down, and confirm with a pilot.
2. **Finder-mesh sensitivity, full form: +2 replays at 4096³.** Only if step 1
   shows a mesh dependence. *Caution:* a 4096³ float32 field alone is 256 GiB
   before particles and finder workspace, so this needs its own pilot and
   probably ≥512 GiB; it is not a routine job.
3. **Legacy catalogue-level convergence — already done here, free.** N-4
   answers the abundance question from `catalogues-pair.npz` with no compute.
   Do **not** commission 21 legacy membership replays unless a specific question
   needs matched legacy properties: the free measurement already shows no
   abundance-convergence difference, so the expensive version has low prior
   value.

**Tier 2 — small new simulations.**

4. **Complete the separability matrix: one 512³/2048³ half-timestep run.**
   N-2 establishes force×particle separability; the untested products are
   force×time and particle×time. A "C with 316 steps" run, compared with F/T,
   tests whether the timestep response depends on particle load. *Resolves:*
   whether the intersection reading extends to the time axis. *Estimate:* C took
   40.3 billed core-hours at 64 cores; doubling the steps gives ≈80–90
   core-hours, plus 3 paired replays (~20 core-hours each at 32 cores, 128 GiB).
   This is the cheapest genuine joint-refinement test and replaces the
   infeasible 8192³ proposal in my Prompt-1 review.
5. **2LPTIC control: one C-class and one E-class pair.** Same box, same seed,
   `lpt_order = 2` via the `2LPTIC_Gui` workflow, compared against the
   first-order C/E particle-refinement interval. *Resolves:* whether the
   measured intervals transfer to the project's default IC workflow — the one
   listed trap this suite structurally cannot address, since all seven runs
   share the same modes. *Estimate:* ≈40 + 64 ≈ 105 core-hours plus 6 replays.
6. **Noise floor: repeat C with a different thread count.** Identical physics,
   different OpenMP force-accumulation order. *Resolves:* the irreducible
   pipeline noise, which no pair in the present suite measures (N-7). *Priority:
   low* — N-7 already shows E/F is not seeded-noise-dominated, so this would
   calibrate rather than change any conclusion. *Estimate:* ≈40 core-hours plus
   3 replays.
7. **Second seed at C-class.** *Resolves:* separation of numerical shift from
   realisation variance, which one box cannot do. Only worth it if a published
   claim needs an error bar rather than a stability range. *Estimate:* ≈40
   core-hours plus 3 replays.

**Not recommended.** A joint 1024³/8192³ run (infeasible on a 1 TB node); a
volume study (needs a campaign, not a control); 21 legacy membership replays
(see 3).

---

## 6. Files and cleanup

Created, all in `BDM-refine/analysis/review-20260909-convergence-claims/` (new
directory; no existing review directory was reused or touched):

| File | Purpose |
|---|---|
| `CLAIMS-REVIEW.md` | this report |
| `claims_checks.py` | six controls: `verify`, `additivity`, `finder`, `selection`, `tails`, `legacy`, `redshift` |
| `claims_results.json` | machine-readable results of all of them |

Reproduce with `micromamba run -n cosemu python3 -B claims_checks.py <stage>`
from that directory. The script opens campaign data read-only and writes only
into its own directory; `python3 -B` prevents bytecode files. All work ran on an
idle login node (load 2.3/128, ~970 GB free) with BLAS/OpenMP threads pinned to
1; total runtime under one minute, peak memory a few GB. No Slurm job was
submitted and no live job existed for this account.

No campaign file was created, modified or deleted; `git status` on the campaign
path is empty, and the latest modification time anywhere under it is 03:55:43,
from the repair session that produced `ceb7e8c` — nothing at or after 04:00,
when this review's controls ran. The earlier review directory
`review-20260909-convergence/` is likewise untouched (latest mtime 02:33). The
archives, figures, slides, LaTeX sources, retained density tapes, snapshots and
membership tapes are untouched. No repair was applied, no merge performed, and
no simulation launched.

---

## 7. Erratum, added 9 September 2026 after `CLAIMS-RESPONSE.md`

Appended after the campaign response at `78bb887`. **Nothing above this line has
been altered**; the report stands as it was reviewed, and the finding
identifiers (C-1, N-1 … N-7) are kept so the response's cross-references still
resolve. Four errors in §1–§5 are confirmed. All four are mine; none of them
changes a campaign measurement or any of the 63 acceptance decisions. Each was
re-verified against `convergence.json` and `convergence-assessment.json` before
being recorded here.

### E-1. §2 N-4 — legacy win counting, and an over-claimed conclusion

*As published:* "114 comparable bins: v3 has the smaller absolute shift in 68,
legacy in 46", under the heading "**No convergence advantage of v3 over legacy
exists in abundance** (stronger than 'is not established')".

*Correct:* **68 v3 / 43 legacy / 3 exact ties.** `claims_checks.py`
`stage_legacy` computes `legacy_smaller_absolute_shift = len(allbins) - wins`
from a strict `abs(v3) < abs(legacy)` test, so the three ties were silently
credited to legacy. The ties are degenerate integer-ratio bins: (13.75,
0.0/0.0), (14.50, 0.0/0.0), (14.50, 1.887/1.887). Median absolute shifts are
unaffected (2.7820 % v3, 3.2130 % legacy).

*The heading is withdrawn.* Asserting that no advantage **exists** is an
assertion of the null. Correlated bins from a single realization, with no power
analysis, support neither superiority nor equivalence. The defensible statement
is the one already in the "precise limit" column and in `CLAIMS-RESPONSE.md`:
**no advantage is established, in either direction**, on abundance only. The
corrected split is slightly more favourable to v3 than the published one, which
does not change that.

### E-2. §1.1 C-1 and §2 N-1 — mass-bin populations, not tracked haloes; and the mass shift is not the binding condition

*As published:* "Two of them (E/F and C/D at 13.00–13.25) agree to better than
0.2 % at z=1 while **the same haloes** disagree by 1.9–2.4 % at z=2 and
3.8–4.1 % at z=0."

*Correct, first point:* these are separately matched populations occupying the
same mass bin at each epoch, **not the same objects followed through time**. No
merger-tree or progenitor measurement was made, and haloes change mass bin
between z=2 and z=0. The phrase "the same haloes" is wrong and should be read as
"the catalogue statistic in the same mass bin".

*Correct, second point:* in all four example bins the `bound_mass` condition is
**True at every epoch**. The z=0 failures are `abundance` and `resolved_vmax`.
The median-mass sign transition therefore cannot by itself explain the width of
the passing intervals, as `CLAIMS-RESPONSE.md` states.

*What survives, on the correct statistic.* Re-running the test on the conditions
that actually bind, the abundance ratio and the median resolved Vmax **also**
change sign between z=2 and z=0 in **7 of the 8** bins — and those are precisely
the conditions that fail at z=2 (abundance) and at z=0 (Vmax):

| Pair | log10 M | Δn % (z2/z1/z0) | median ΔVmax % (z2/z1/z0) |
|---|---:|---|---|
| E/F | 13.00 | +0.90 / −2.81 / −6.43 | +0.57 / −0.96 / −5.19 |
| E/F | 13.25 | +7.41 / +0.74 / −4.48 | +0.93 / −0.23 / −4.15 |
| E/F | 13.50 | +14.63 / +1.05 / −2.84 | +1.47 / +0.43 / −3.13 |
| C/D | 12.75 | +0.68 / −3.00 / −6.89 | +0.13 / −1.54 / −6.14 |
| C/D | 13.00 | +4.67 / −2.95 / −5.82 | +0.72 / −0.79 / −4.92 |
| C/D | 13.25 | +6.88 / −0.21 / −3.97 | +1.19 / −0.07 / −3.99 |
| C/D | 13.50 | +9.52 / +2.34 / −2.18 | +1.55 / +0.46 / −3.07 |

The single exception is E/F 12.75, where abundance and Vmax worsen monotonically
and z=1 has simply not yet crossed tolerance. So the substantive caution — that
the wider z=1 force-mesh intervals reflect a sign transition in the population
response rather than stability across epochs, and that no interpolation in
redshift is safe — holds, but it must be argued from abundance and Vmax rather
than from median mass. Three outputs do not locate a continuous zero crossing or
identify its cause.

*§4 inherits this.* The sentence "the wider z = 1 intervals reflect a zero
crossing rather than convergence and z = 0 is the binding epoch" should read:
the wider z = 1 intervals reflect a sign transition in the binned abundance and
Vmax response, and z = 0 is the **more restrictive** epoch for the lower-mass
force comparison — not the binding epoch for every property, and not for the
timestep comparison, whose error decays monotonically.

### E-3. §2 N-5 — the scatter/median ratio mixes bins

*As published:* the 16–84 half-width is "**9–25× the median shift** for the
particle comparisons".

*Correct:* `claims_checks.py` `stage_tails` reports
`max(half_width) / max(|median|)` over the passing bins, so numerator and
denominator may come from **different bins**. The per-bin ratios for C/E z=0 are
2.7, 9.0, 10.3, 13.7, 14.4, 18.7, 28.0 — and one degenerate bin at 14.00–14.25
whose median is −0.001 % gives a meaningless 2236. The published range is an
artefact of the aggregate. The qualitative claim survives on per-bin values:
halo-to-halo scatter exceeds the median shift by roughly one order of magnitude
for particle refinement (C/E z=0, 12.75–13.00: percentiles −8.303 / −0.524 /
+6.737 %, half-width 7.520 %) and by ≈1× for F/T. Per-bin values, not the
aggregate ratio, should be quoted.

### E-4. §2 N-6 — pooled unmatched fractions understate the worst bin, and the median inference does not follow

*As published:* "unmatched eligible reference haloes are 1.47 % (C/E), 0.77 %
(E/F), 0.14 % (F/T), 1.57 % (A/E) … at ≤1.6 % of the sample it cannot move a bin
median by more than a few tenths of a percent."

*Correct, first point:* those are **pooled** fractions over all used bins. Under
the frozen completeness definition — which applies the particle floor to **both**
matched haloes — the largest single-bin z=0 loss among those four pairs is A/E
13.50–13.75 at ≈2.73 % (`CLAIMS-RESPONSE.md` quotes 2.72 %), and across all seven
pairs it reaches **2.90 %** at B/C 12.75–13.00 (n = 5352). Per-bin losses, not
pooled ones, are the relevant bound.

*Correct, second point:* the inference "≤1.6 % … cannot move a bin median by more
than a few tenths of a percent" does not follow. A small lost fraction bounds
**ranks**, not a median expressed in percentage units; with a gap in the
distribution, removing one object from 101 can move the median arbitrarily far.
The response's counterexample is valid. What N-6 does establish is that the loss
is small and rank-bounded, that the lost objects have lower Vmax **and** lower
mass (so −15 % is not an isolated Vmax bias at fixed mass), and that a 6× wider
candidate search reproduced the matches exactly in one E/F z=0 subvolume.

### Separately withdrawn: §2 N-7

N-7 argued that E/F cannot be chaotic amplification of its IC seed because F/T
begins from a "1.55 % coherent velocity difference" and ends smaller. That
1.55 % is a difference in **stored** velocities arising from different leapfrog
half-step output epochs, not a perturbation of the same physical state at the
same time, so its amplitude cannot be set against E/F's position difference. The
argument does not stand and is withdrawn. The accompanying observation is
unaffected and remains the useful part: **no pair in this suite differs by
numerical noise alone, so the irreducible noise floor of the pipeline is
unmeasured**, and §5 item 6 remains the way to measure it.

### Not affected

The 63 acceptance decisions and their criteria (re-derived twice, 0 mismatches);
the preservation checks; N-2 separability (independently reconfirmed: C/D vs E/F
median differences of 0.3483 pp in mass and 0.2726 pp in Vmax over eight common
z=0 bins, 0.7240 and 0.3334 pp across all epochs); N-3 finder-mesh candidate
counts; the §1.2 audit of the Prompt-1 review; and the §5 follow-up design as
amended by `CLAIMS-RESPONSE.md`'s resource corrections.
