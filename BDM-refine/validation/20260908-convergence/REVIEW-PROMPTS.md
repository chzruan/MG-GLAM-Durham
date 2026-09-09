# Independent review prompts for the completed BDM convergence campaign

Use Prompt 1 for the full review. Prompt 2 is an optional second, independent
scientific review focused on the convergence claim. Each can be run in a fresh
session by asking the agent to read this file and execute that prompt only.
These prompts specify review work; they do not authorize repairs or a merge.

## Prompt 1 — Full independent review before merging

```text
Independently review the completed BDM convergence branch. Read the source and
test the claims; do not treat existing reports, receipts, successful jobs or
earlier agent reviews as proof that the analysis is correct.

Workspace:
/cosma8/data/dp203/dc-ruan1/mgglam_claude/MG-GLAM

Pinned review target:
- Branch: validation/bdm-convergence-20260908
- Campaign/result tip: cb87021d140fd2ea72841b6cea32abf090c255dd
- Base, true main cz: e289c3f7ea28d38390a1663e50bb2784b85abc84
- Review the full base-to-target difference, not just the final two commits.
- Later additions of review instructions/reports are outside that frozen diff.
- Record the actual HEAD, cz, working-tree status and any divergence from these
  pins. Do not reset or switch the shared checkout to make it match.

C below means BDM-refine/validation/20260908-convergence.
Production PMP2linker.f90 is reported unchanged from the pinned cz base. Verify
this. The new branch chiefly adds matched ICs, a controlled evolution schedule,
common-mesh replay adapters, analysis, reporting, submission/resume and cleanup.
Consult the existing finder implementation and earlier repair tests wherever
their semantics affect this campaign. Distinguish inherited behavior from a
defect introduced by the reviewed changes.

Execution and preservation:
- Read AGENTS.md. Every Python execution must use
  micromamba run -n cosemu python3 -B
  and package installation, if necessary, must use
  micromamba run -n cosemu pip
- Review first. Do not modify production source, existing analysis code,
  receipts, catalogues, snapshots, density tapes, plots or slides. Do not merge,
  push, delete branches or commit review outputs as part of this prompt.
- Write only new review reports and small independent controls. Check whether
  the proposed review destination already exists before choosing it.
- Light, bounded controls may run on an idle login node after inspecting load
  and memory. Keep BLAS/OpenMP thread counts explicit. No new full simulations
  are needed for the initial review. Identify any additional simulation needed
  to resolve a material uncertainty, with its purpose and resource estimate.
- A justified substantial reanalysis/replay may use Slurm: default to shared
  cosma8-serial with actual CPUs and explicit memory. Inspect existing pilot
  accounting first. A previous full analysis peaked at 7.83 GiB: allow at least
  16 GiB for an equivalent rerun. High-resolution paired replays reached
  191.99 GiB: allow at least 256 GiB for an equivalent rerun. Do not put such
  work on a login node. cosma8 is an exclusive 128-core node allocation;
  one-core submission there can bill the whole node. Record job IDs, scripts,
  dependencies, ReqTRES, AllocTRES, TotalCPU, elapsed time, MaxRSS and both
  expected and wall-limit core-hours for any new substantial job.
- Before executing a campaign driver, inspect its side effects. Several write
  to fixed C paths even when imported or given an alternate input. Use an
  isolated review harness, explicit output paths or verified redirection.
  BDM_CONVERGENCE_ROOT also affects REPO derivation; do not assume changing it
  alone safely isolates every path. Never rerun analyze.py, assess.py or
  render.py against the original campaign merely to compare their output.
- Avoid unnecessary large-file hashing. State exactly which raw files you
  rechecked, which hashes you only verified through receipts, and which earlier
  checks you relied on. Receipt agreement is not independent data validation.
- File quota matters. Use one small control script and a compact results file
  where practical; do not create a file per test. Read needed archive members
  directly or restore selected members into owned scratch. Remove your own
  temporary builds/previews after retaining evidence. Preserve unrelated
  untracked files and all existing scientific products. Slides/LaTeX/PDFs must
  remain uncommitted.

Start with:
- C/README.md, C/CONVERGENCE.md and C/common.py.
- C/ic/README.md, C/replay/README.md and their implementation files.
- C/analyze.py, C/assess.py and C/plots/plot_convergence.py.
- C/convergence.json, C/convergence-assessment.json and C/final-review.json.
- C/jobs.json, C/accounting.json, C/build.json, C/ic-build.json,
  C/replay-build.json and C/executables.json.
- C/*-ic.json, C/*-simulation.json, C/replay-*-ng2048.json and
  C/v3-validation-*-ng2048.json, read programmatically rather than dumping all
  large JSON documents into context.
- The frozen independent membership checker is
  BDM-refine/validation/20260906-n1024/run_validation.py; examine its assumptions
  as well as those of the producer.

Stage 1: physics and statistical validity

1. Reconstruct the actual experiment from code, inputs and receipts:
   A: 256^3 particles / 2048^3 evolution mesh
   B: 512^3 / 1024^3
   C: 512^3 / 2048^3
   D: 512^3 / 4096^3
   E: 1024^3 / 2048^3
   F: 1024^3 / 4096^3
   T: 1024^3 / 4096^3, every normal timestep subdivided into two.
   Claimed common settings: GR, L=256 Mpc/h, z_init=100, outputs z=2,1,0,
   common finder mesh 2048^3, iVirial=1, Rext=0.15, MassMin=2.5e12 Msun/h.
   Check actual settings, output epochs, particle masses and configuration
   parsing, including legacy first-call behavior.

2. Review matched IC construction: shared Fourier amplitudes and phases,
   master normalization, mode truncation/Nyquist treatment, physical origin,
   nested particle IDs, FFT conventions and velocity normalization. Shared
   seed values alone are insufficient. Separate complete small controls from
   sampled production modes and all-row production coordinate comparisons.
   This existing campaign used native first-order Zel'dovich ICs, although new
   simulations now default to corrected FML 2LPTIC. Verify that the exception
   is truthful and that the conclusions do not silently extend to 2LPTIC.

3. Review the 158/316-step construction and leapfrog velocity epochs. Check
   whether normal endpoints are preserved, how T's initial velocities change,
   and which measured differences include output-time staggering. Examine
   periodic coordinate conversion and its rounding/edge guard. E/F physical
   initial positions are reported close, not bit-identical.

4. Review common-mesh adapters and publication membership extraction. Check
   physical units, one-based coordinates, periodic images, particle ordering,
   density-tape metadata, identical legacy/v3 density inputs, compiler flags,
   and whether instrumentation changes finder state. Legacy means a8c7715;
   the v3 executable must match the frozen audited source.

5. Independently challenge halo matching. Test ID projection between lattices,
   both overlap denominators, mutual-best decisions/ties and the spatial
   candidate filter. Compare against exhaustive matching in small fixtures;
   a nearest-centre check alone is not sufficient. Consider empty catalogues,
   zero shared tracers, split/merged objects and periodic boundary cases.
   Assess selection bias from losing difficult matches, including its effect
   on apparently stable median properties and reported completeness.

6. Check bound mass from float64 member counts times stored particle mass;
   distinguish it from the Rext-expanded aperture mass/radius. Check bin edges,
   whole-bin mass cuts, publication completeness, both members' particle floors,
   positive resolved Vmax filtering and unresolved incidence. Verify reference
   mass binning, units, ratio direction, match completeness and left matched
   fraction. Examine shape and velocity results without assuming the mass/Vmax
   acceptance screen certifies those properties.

7. Independently recompute the 63 assessments from measured numbers, without
   calling assess.evaluate as the oracle. Validate minimum counts, finite-value
   handling, paired eight-octant jackknife normalization/denominators and merged
   passing intervals. Distinguish quantile scatter, median uncertainty, spatial
   jackknife uncertainty and absolute accuracy. Ask whether the chosen cuts or
   post-hoc mass-range selection make the conclusion misleading.

8. Challenge membership validation with positive and negative controls. The
   reported 377,475 rows have zero exact duplicate member sets, repeated IDs or
   mass/count inconsistencies. All 21 production host checks also examined zero
   eligible higher-priority neighbours. Do not count this null host check as a
   demonstrated production positive control. Distinguish exact member-set
   duplication from broader physical duplication and overlapping memberships.

Stage 2: numerical behavior, reproducibility and resource use

9. Examine dtype promotion, integer ranges, byte order, periodic distances,
   finite/empty inputs, thread-dependent arithmetic and threshold boundaries.
   Inspect meaningful existing controls, then add small independent tests for
   uncovered failure modes rather than duplicating implementation logic.

10. Review common.py, campaign.py, timetable.py, build.py, replays.py,
    submit.py, accounting.py, cleanup.py and archive_scratch.py. Verify frozen
    source/configuration/binary linkage, cached-match identity and rejection of
    changed inputs. In isolated fixtures, inject failures during output/receipt
    publication and resume. Look for false completion, stale reuse, missing
    receipts, concurrent writers, partial files and destructive retry behavior.

11. Review archive verification and cleanup recovery: active/unknown job guards,
    changed originals, corrupt archives, interrupted deletion, selected-member
    restoration and preservation of scientific inputs/outputs. Some original
    build and launch paths now exist only in documented archives; a missing
    unpacked path alone is not proof that provenance was lost.

12. Reconcile achieved Slurm cost and memory against allocation records. The
    report claims 2133.48 billed core-hours and 1891.57 consumed CPU-hours,
    including superseded render attempts. Check for parent/step double counting
    and confusion between requested, billed and consumed CPU time. Distinguish
    expected cost from the summed wall-limit ceiling. Identify material scaling
    or memory hazards with measured evidence.

13. Inspect reporting and local artifacts if available:
    C/figures/bdm_convergence_n300.pdf (18 pages),
    C/slides/bdm_convergence.pdf (29 pages), its .tex and make_slides.py,
    C/plot-controls.json, C/plot-manifest-n300.json,
    C/presentation-validation.json, C/render-validation.json and
    C/visual-review.json. Slides and figure artifacts are intentionally ignored
    by Git. Match hashes and visually inspect representative scientific panels.
    Specifically check the repaired absolute velocity/centre-distance axes,
    scatter visibility, labels, cuts, omitted bins and scientific wording.
    An automated render receipt predates its separate visual review; compare
    the actual hashes and sequence before calling it stale or contradictory.
    Assess whether a clean checkout plus retained provenance can reproduce the
    measurements, and separately whether it can reproduce the presentation.

Claims to verify, not expected answers:
- All seven simulations and 21 paired replays completed, with 21 validated v3
  catalogues and 63 pair/redshift/particle-floor assessments.
- At z=0 and a 300-particle floor, passing log10(Mbound/[Msun/h]) intervals are
  C/E: 12.75-14.75, E/F: 14.00-14.75, F/T: 12.50-14.75.
- E/F's lowest eligible z=0 bin has about 9.2% fewer haloes and 7.3% lower
  matched median mass and Vmax on the 2048^3 evolution mesh.
- Force resolution limits lower-mass z=0 results under the adopted screen;
  higher-redshift timestep sensitivity remains.
- These are conditional catalogue-stability measurements at a fixed finder
  mesh. They do not by themselves demonstrate finder-mesh convergence,
  convergence of every halo property, or absolute physical accuracy.

Deliverables:
Write BDM-refine/analysis/review-20260909-convergence/REVIEW.md, choosing a new
suffix if that directory already contains another review. Retain at most a
small review_checks.py and review_results.json alongside it when needed.

Lead with actionable findings ranked P0-P3. For each finding give the pinned
file/line, concrete trigger, expected versus observed behavior, reproduction
or quantitative evidence, affected conclusions/products, and a minimal repair
direction. Clearly mark confirmed defects, unresolved hypotheses and already
disclosed limitations. Do not manufacture findings to fill a quota.

Then state what was independently tested, what relied on prior evidence, and
what was not tested. Include a claim/evidence/verdict table and answer:
1. Is this branch suitable to merge into cz, and what must be fixed first?
2. Precisely which properties, mass/redshift ranges and numerical settings
   support saying the refined BDM catalogues passed a convergence test?
3. Which stronger claims remain unsupported, and what smallest additional
   measurement would address each material gap?
4. Does this campaign establish any convergence advantage over legacy BDM,
   beyond the separate earlier correctness repairs? Do not assume it does.

If no actionable defect is found, say so explicitly and describe the remaining
limits. Report your new files and cleanup. Finish the review without applying
repairs, merging or changing the frozen campaign.
```

## Prompt 2 — Independent challenge to the convergence conclusion

```text
Read BDM-refine/validation/20260908-convergence/REVIEW-PROMPTS.md. Use the pinned
scope, execution rules and preservation requirements in Prompt 1, but execute
only this focused scientific review. Do not read another new review's verdict
until you have recorded your own conclusions.

The question is: what exactly does this completed experiment establish about
convergence of the refined BDM halo finder and its catalogues?

Independently reconstruct the seven simulation settings and the 21 comparisons
at each of the three particle floors from retained data and code. Derive the
acceptance decisions yourself; do not use assess.py's decisions as the oracle.
Use small controls and selected retained data to challenge the interpretation,
preserving every original product. Do not launch a new simulation campaign.

Concentrate on these possible ways the apparent agreement could mislead:
- The finder analysis mesh stays fixed while evolution resolution changes.
  Can a shared finder limitation suppress sensitivity in every comparison?
- Each refinement is measured at particular settings of the other parameters.
  Does taking the intersection of separate passing ranges justify any joint
  convergence statement at an untested particle/force/timestep combination?
- Shared-lattice projection, restricted spatial candidates and mutual-best
  selection may preferentially keep stable haloes. Can the completeness
  definition and minimum counts reveal the resulting selection bias?
- Unresolved Vmax filtering, publication limits, bin migration and selection
  on both haloes' particle counts may affect apparent convergence near cuts.
- Median agreement does not constrain tails, halo-to-halo scatter, shape or
  velocity accuracy. Does the scientific wording preserve that distinction?
- Eight correlated spatial octants in one realization do not provide all
  relevant uncertainty estimates. Are sparse bins or NaNs treated honestly?
- Force and timestep sensitivity may vary with redshift; resolution errors
  may be nonmonotonic or cancel between settings. Is any generalization beyond
  the measured comparisons justified?
- All runs share first-order ICs, native velocity staggering and related
  approximations. Could common biases survive every passing comparison?
- Zero exact duplicate sets and null host-exclusion checks answer different
  questions from physical uniqueness and numerical convergence.
- Legacy/v3 paired outputs exist. Has anyone actually demonstrated superior
  convergence of v3, or only established conditional stability of v3 outputs?

Do not label every missing experiment a code defect. Quantify the impact of
each confirmed issue where possible, and distinguish a failed tolerance from
an untested configuration or a limitation already disclosed in the report.

Write a separate new report under
BDM-refine/analysis/review-20260909-convergence-claims/ (use a fresh suffix if
occupied). Keep independent control code/results compact. Include:
1. Any confirmed calculation or interpretation errors, with file/line evidence.
2. A table of supported claims, supporting measurements and precise limits.
3. One defensible paragraph for a paper or presentation answering whether the
   refined BDM passed convergence testing.
4. A prioritized minimal follow-up design, separating finder-mesh sensitivity,
   additional force/time refinement and a 2LPTIC control. Reuse existing
   snapshots where scientifically appropriate; distinguish replay-only work
   from genuinely required new simulations. Recommend no run without stating
   which uncertainty it resolves.

Do not alter slides, measurements or production source. Report what you
actually tested and clean up only scratch created by this review.
```
