# Response to the convergence claims review, 9 September 2026

Branch: `fix/bdm-convergence-claims-20260909`, from `ceb7e8c`.
The [Prompt-2 report](../../analysis/review-20260909-convergence-claims/CLAIMS-REVIEW.md)
and its two supporting files are preserved as received. Its disclosure matters:
this was a second pass by the Prompt-1 reviewer, not a blind second opinion.

All seven executable review controls reproduce exactly. All 63 assessments and
their criteria are unchanged. The 56 production IC/simulation/replay/checker
receipts, `convergence.json`, `analyze.py`, the production finder and frozen
checker match `ceb7e8c`. No simulation, replay, matching or large-tape scan was
needed. `claims_response.py` reproduces these checks and records the input/source
hashes and additional measurements in `claims-response.json`.

## Measurements adopted, with their actual scope

**Redshift dependence is now explicit in the report and slides.** Eight E/F and
C/D fixed mass bins change sign between z=2 and z=0. At log10 mass 13.00–13.25,
E/F median mass shifts are +1.851%, −0.166%, −4.072%; C/D gives +2.363%,
+0.077%, −3.805%. Small z=1 medians are consistent with a sign transition in
the population response. They are not evidence of stability across epochs.

Two qualifications to C-1/N-1 are necessary. These are separately matched
populations in the same mass bin at each epoch, **not the same haloes followed
through time**; no merger-tree measurement was made. Three outputs do not
locate a continuous zero crossing or identify its cause. Also, both example
bins pass the 5% mass condition at all three epochs: their z=0 failures are
**abundance and resolved Vmax**. The mass sign transition alone cannot explain
the full-screen interval widths. The recorded z=1 passes remain valid under
the stated screen. z=0 is more restrictive for the lower-mass force comparison;
it is not the binding epoch for every property or for timesteps.

**The shape summary now quotes the table maximum:** A/E z=1, median c/a 8.92%,
alongside E/F z=0 at 5.86%. The abundance uncertainty paragraph now states that
15 of the 21 bins satisfying the difference-only cut also pass the full screen.
Summed bound mass explicitly counts shared particles repeatedly; a catalogue
is not a partition of the particle set.

**There is evidence for approximate force–particle separability of medians.**
C/D and E/F differ by at most 0.3483 pp in median mass and 0.2726 pp in median
Vmax in their eight common usable z=0 bins. Across all three epochs the maxima
are 0.7240 and 0.3334 pp. This does not certify abundance, shapes or interactions
with time resolution. The particle-chain residual is a median 0.162 pp and
maximum 0.942 pp over 18 **bin/property cases combining mass and Vmax**; 0.2 pp
is not an upper bound, and multiplying population medians is not a test of
individual matched triplets.

**Legacy abundance has now been measured, but absence of an advantage is not
proved.** Exact integer-ratio comparisons give 68 v3 wins, **43 legacy wins
and 3 ties**, over 114 usable bins. The review's 46 legacy wins included its
three ties. Median absolute shifts reproduce as 2.7817% versus 3.2130%.
Correlated bins from one realization do not establish statistical superiority,
equivalence, or the categorical claim that no advantage exists. Legacy matched
properties remain unavailable without legacy member lists. Correctness evidence
for the repairs and numerical convergence are distinct claims.

**Scatter is reported beside its own bin median.** For C/E z=0,
12.75–13.00, the mass-shift percentiles are −8.303/−0.524/+6.737%; the
16–84 half-width is 7.520%. The review's aggregate scatter/median ratio divides
the largest width by the largest absolute median, potentially in different
bins. The stored per-bin values avoid that ambiguity.

## Stronger interpretations the measurements do not support

- **Finder mesh:** varying candidate counts show sensitivity at one analysis
  mesh; they do not rule out attenuation of the response of published haloes.
  Median apertures of 2.97–3.04 cells are observations in samples sharing a
  publication cut. They do not establish a mesh-induced radius floor, an
  absolute error, or exact cancellation of that error in ratios.
- **Selection:** 0.14–1.57% are pooled unmatched fractions above whole-bin lower
  cuts, not per-bin bounds or a calculation of the complete property selection.
  The frozen completeness definition also applies particle floors to both
  matched haloes. The largest loss in the four pairs' eligible z=0 bins is
  2.72%. E/F's pooled unmatched sample has lower Vmax, but it also has lower
  mass, so −15% is not an isolated Vmax bias at fixed mass. A small lost fraction
  alone bounds ranks, not a median shift in percentage units: the reproducible
  101-object counterexample loses one object yet shifts its median by 25%.
  The widened spatial search established equality in one E/F z=0 subvolume;
  it is not a global proof for every pair/epoch.
- **Noise and ICs:** F/T's 1.55% initial stored velocity difference represents
  different leapfrog half-step epochs. It is not a physical velocity
  perturbation at the same time, and its amplitude cannot be compared with an
  E/F position perturbation to exclude chaotic amplification. No same-physics,
  noise-only evolution pair was run. The source-bound IC check rules out
  excess above its roundoff bound; it does not prove the periodic clamp never
  acted. A common first-order IC prescription can interact with resolution;
  its bias need not cancel exactly. The 2LPTIC control remains useful.

## Follow-up design and resource corrections

No additional jobs were launched. Retain shared `cosma8-serial`, measured
thread counts and explicit memory with headroom. The successful IC recheck's
15.99/16-GiB peak remains a reason to request at least 24 GiB next time;
rerunning it only to obtain a lower percentage of requested memory adds no
validation evidence. E-equivalent evolution needs at least 128 GiB.

| Question | Smallest informative addition | Limitation / cost basis |
|---|---|---|
| Does finder coarsening alter the E/F force signal? | E and F z=0 at a 1024^3 finder mesh: two replays | Useful cheap diagnostic; a null result does not establish stability under refinement to 4096^3. Pilot runtime and memory; do not scale total memory by mesh size alone. |
| Does finder refinement alter that signal? | E and F z=0 at a 4096^3 finder mesh, with retained 2048^3 controls | The field alone is 256 GiB; particle and finder workspace need a measured allocation. A coarsening pass does not eliminate this question. |
| Does the timestep response change at a cheaper combined setting? | C with 316 steps | About 80–90 evolution core-hours from C's measured 40.28; it changes both particle load and force mesh relative to F/T. It tests a combined response, not the complete interaction matrix. |
| Isolate force–time interaction | E with 316 steps, compared with E and F/T | Particle load stays 1024^3; rough evolution estimate 110–130 core-hours from E's 63.57 and measured F/T scaling. Replays and a resource pilot are additional. |
| Isolate particle–time interaction at the fine force mesh | D with 316 steps, compared with D and F/T | Force mesh stays 4096^3; D already cost 477.40 core-hours. This is substantially more expensive than C with half steps. |
| Does the C/E particle result transfer to default ICs? | Matched 2LPTIC C and E runs, with verified common modes | About 104 evolution core-hours from existing C+E; IC generation/conversion, validation and six epoch replays are additional, not included in that total. |
| How variable is a refinement shift across realizations? | A second-seed matched refinement pair, e.g. C/E | A single unpaired second-seed C cannot measure variation of the paired refinement shift. |

An 8192^3 float32 field alone needs 2 TiB and remains infeasible on a 1-TB node.
The current pass intervals certify only the specified conditional screen, with
separate mass and redshift limits; none of these interpretations changes its
63 decisions.

The updated Beamer PDF and LaTeX remain local and Git-ignored as requested.
All 35 pages pass the layout checks. The template retains a harmless hyperref
warning when it removes the translation command from the appendix bookmark;
this affects bookmark metadata, not the rendered scientific text or layout.
The final verification records the render, retained hashes and owned temporary
file cleanup. Earlier `review-followup.json` is historical evidence for the
Prompt-1 repair at `ceb7e8c`; its recorded artifact hashes refer to that version.
