# Response to the independent convergence review, 9 September 2026

Repair branch: `fix/bdm-convergence-review-20260909`, from `d485bc7`.
The [independent review](../../analysis/review-20260909-convergence/REVIEW.md)
covered the original `e289c3f..cb87021` campaign. Its original report and
controls are retained unchanged. Production finder source, the frozen
membership checker, raw scientific products, `analyze.py` and
`convergence.json` remain unchanged.

The review's requested documentation changes and production IC validation
rerun are addressed. All 63 original acceptance decisions and criteria are
unchanged. The revised statement is conditional stability of the measured
abundance, median bound mass, median resolved Vmax and completeness, at the
tested settings. It does not certify shape, tails, velocities, absolute
accuracy, a different finder mesh or a 2LPTIC suite.

## Findings and disposition

| Finding | Action and evidence |
|---|---|
| P1-1: unscreened shapes | `assess.py` now derives shape maxima inside the passing intervals, with counts, for every comparison and particle floor. `CONVERGENCE.md` reports all 21 pair/redshift cases at nominal floor 300; the README and slides state the principal results. E/F z=0 reaches 5.86% in c/a; F/T reaches 4.70–5.42%; the coarser A/E comparison reaches 8.92% at z=1. Axis definitions explicitly include empirical correction and transverse-axis reordering. |
| P2-1: abundance uncertainty | Clarified the scope of the paired octant jackknife. Added both relative and ratio-unit hypothetical independent-count scales per bin. They are labelled as benchmarks, with covariance omitted, and are not substituted for the paired uncertainty or used to change decisions. See the statistical qualification below. |
| P2-2: superseded IC validation | Reran the current validator as Slurm job 11960736. All 1024^3 rows passed. The new receipt binds the executed source, support modules and all seven input receipts, including the frozen zipapp case. The old receipt is retained inside `review-followup.json`. |
| P2-3: host check | Corrected the documentation to the actual priority-ordered **extended-aperture** rule. Executed the frozen checker's host block on six controls, including three cases that exercise its rejection branch. The four reverse-orientation pairs found by the review are allowed by this rule. No checker or physics change is warranted by those pairs. |
| P3-1: shared memberships | State the measured 0.06–0.13% rate of excess memberships. This is a global rate of repeated occurrences across haloes, not a bound on any individual halo's mass error. |
| P3-2: nominal floor | Report the common mass cut and implied integer particle thresholds per pair. E/F and F/T require at least 1868 particles at the publication cut and 2362 at the first whole bin. Nominal floors 100/300/1000 do not supply independent tests for these pairs. |
| P3-3: memory | Added at least 128 GiB for E-equivalent evolution. The IC recheck used one shared core and 16 GiB, based on the previous 11.99-GiB peak; its new accounting peak is 15.99 GiB. Recommend 24 GiB for an equivalent future check pending separate process-RSS/cache measurements. |
| P3-4: ignored artifacts | Preserve the documented local-data dependency. Scientific data and verified launch archives remain available in this workspace with hashes/restoration instructions. Slides and LaTeX remain uncommitted as requested. A clean Git checkout alone is not claimed to reproduce the data or presentation. |
| P3-5: assessment source | Added `assess_source_sha256`, verified for ordinary execution and import from a frozen zipapp. The new property definitions live in the assessment so that the frozen measurement producer and cached-match identity remain valid. |

The pre-existing E/F IC receipts remain historical records; their omitted
sampled-mode fields are supplied by the new source-bound all-seven validation.
No old receipt has been made to claim an assertion it did not originally run.

The redundant matching counters, concurrent same-path JSON writer hazard and
missing-key handling for a future reference `Init.dat` remain disclosed latent
items. The reviewed campaign used dependency-gated writers and contained all
requested configuration keys. No failing production case was found. This
repair does not change those frozen execution paths or present the branch as a
general-purpose concurrent campaign service. The inherited small absolute mass
normalization offset is unchanged.

## Statistical qualification of P2-1

Identical paired octant counts can yield zero jackknife scatter without an
implementation error. For the reviewed C/E bin there are 47 haloes in each
catalogue and the octant count vectors are identical. Omitting any common
octant leaves the ratio at one. This estimates zero resolved spatial variation
for that realization; it does not establish zero ensemble uncertainty.

For a ratio R = N_L/N_R, a first-order variance calculation contains

`Var(log R) ~= Var(N_L)/N_L^2 + Var(N_R)/N_R^2 - 2 Cov(N_L,N_R)/(N_L N_R)`.

The suggested `sqrt(1/N_L + 1/N_R)` drops the covariance and assumes Poisson
variances. It is 20.63% for **two independent** counts of 47, whereas these
catalogues share their initial realization. It is neither the measured paired
ratio uncertainty nor a general bound on it. We retain it as an explicitly
labelled comparison scale and preserve the existing descriptive screen.

The review's 21 floor-300 bins pass the abundance-difference condition alone;
some fail the uncertainty or other conditions of the full screen. We state
that distinction. Tight paired mass/Vmax medians do not establish abundance
precision. A confidence statement about the underlying abundance difference
would require an appropriate covariance/error analysis and sufficient data.

## Host-rule qualification of P2-3

`PMP2linker.f90:696–750` specifies and implements the following condition: a
lower-priority centre is removed if it lies inside a higher-priority halo's
**extended** aperture. `Rvir` is assigned that aperture at line 1014. Priority
is bound mass, then stable candidate index. The frozen checker implements the
same orientation and radius. Replacing it with symmetric exclusion on
unextended SO radii, as suggested in the review, would test another rule.

The checker's examined count counts pairs already inside the query aperture;
zero therefore provides no broad coverage evidence. It can exercise its
violation branch: the follow-up controls execute the actual checker block and
reject an ordinary violation, an equal-mass priority violation and a periodic
violation. They accept a separated pair, the allowed reverse orientation and
an exact-radius boundary pair. The last has an examined count of one and no
violation, consistent with the strict inequality.

The four reverse-orientation pairs and the overlap scan are useful additional
measurements. They do not establish a defect in the implemented host rule or
require disjoint memberships.

## Production IC rerun

Job 11960736 completed in 604 seconds on one `cosma8-serial` core with a 16-GiB
request. Expected cost was 0.25 core-hours; the wall-limit ceiling was 0.5.
Actual cost was 0.16778 billed core-hours and 0.09357 CPU-hours, with Slurm
MaxRSS 15.99 GiB. The actual request/allocation and immutable launch hashes are
recorded in `jobs.json` and `accounting.json`.

All 1,073,741,824 fine particle rows were examined. E/F physical velocities
are identical; F/T position words are identical and velocity staggering passes.
E/F maximum coordinate difference remains 6.103515625e-5 Mpc/h; component RMS
is 2.3834971611433514e-7 Mpc/h. There are zero components exceeding the pure
roundoff bound, zero nonlocal excesses and zero excesses lacking the exact
clamp sentinel. Maximum inferred displacement is 0.34850568142395566 Mpc/h.
These excess diagnostics do not count all uses of the native edge clamp.

No simulation, finder replay or membership matching was rerun for this repair.
Original-campaign cost remains 2133.47722 billed core-hours; with the IC recheck
the total is 2133.64500. Local presentation work is recorded separately.

## Follow-up experiments remain separate

- Finder-mesh sensitivity can use existing snapshots. For E and F at z=0,
  analysis meshes 1024^3 and 4096^3 require **four new replays** if the retained
  2048^3 results are reused. Pilot the 4096^3 case before assigning resources:
  its float32 field alone is 256 GiB, before particles and finder workspace.
  The review's cost estimate is provisional. These four replays would address
  this limited panel, not all redshifts or configurations.
- Existing legacy catalogues allow abundance and distribution comparisons
  without new simulations. Equivalent shared-member matched-property tests
  cannot be obtained simply by changing `variants['v3']` to `variants['legacy']`:
  none of the 21 legacy receipts contains a membership tape/index. They require
  validated legacy membership extraction and finder replays, or a clearly
  different matching method with its own controls. No superiority in numerical
  convergence over legacy is established by the current v3-only analysis.
- A joint-refinement claim, a 2LPTIC claim or a volume/cosmology claim needs
  controls designed for that question. A second seed should repeat a matched
  refinement pair if the aim is to separate numerical shifts from realization
  variation. An 8192^3 float32 field alone is 2 TiB, so that proposal is not a
  D-class single-node rerun on a 1-TB COSMA8 node.

`review-followup-controls.json` records the independent reassessment, source
identity, host-rule, shape-summary and effective-cut controls.
`review-followup.json` records the original review/receipt identities and the
completed repair evidence. No merge into `cz` is performed by this response.

The revised deck has 31 slides and no build warnings. Representative scientific
and layout inspection is bound to the current PDFs in `visual-review.json`.
The six finished launch files for the prior final-cleanup job and the IC
recheck were consolidated into `launches-review-20260909.tar.gz`; every member
was verified before removing the originals. Owned temporary slide previews
were removed after review. Scientific data and unrelated working files remain.
