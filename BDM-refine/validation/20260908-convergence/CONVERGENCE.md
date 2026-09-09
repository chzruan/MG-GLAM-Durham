# BDM convergence measurements

**Complete seven-run measurement.**

All primary catalogues use the same 2048^3 analysis mesh. The simulation box is 256 Mpc/h; outputs are z=2,1,0. Particle resolution, evolved force resolution and timestep size are compared separately.

**IC scope:** all seven runs use native GLAM first-order (Zel'dovich) initial conditions at z_init=100. This completed suite is an existing-run exception to the default 2LPTIC workflow for new simulations; it does not establish convergence with 2LPTIC initial conditions.

Bound mass is the original-member count times the stored particle mass, evaluated in float64. Reported aperture mass and radius include the empirical Rext expansion. Positive Vmax values are compared only when both haloes resolve them; unresolved values are counted separately. Halo matches require mutual best shared-lattice membership overlap, with at least 50% in each object.

The table applies a descriptive working screen **only to abundance, median bound mass, median resolved Vmax and reference completeness**: at least 30 objects in each required bin; abundance within 5% with paired octant sigma at most 5%; median bound mass within 5%; median resolved Vmax within 2%; reference completeness at least 90%. These choices are adjustable. They do not certify those absolute accuracies. Scatter, individual criteria and 100/300/1000-particle floors are retained in the JSON results.

| Comparison | z | Nominal particle floor | Eligible bins | log10 mass intervals meeting the stated screen |
|---|---:|---:|---:|---|
| A/C (particle) | 2 | 300 | 1 | None |
| C/E (particle) | 2 | 300 | 4 | 12.75–13.75 |
| A/E (particle) | 2 | 300 | 1 | None |
| B/C (force) | 2 | 300 | 4 | 13.25–13.50 |
| C/D (force) | 2 | 300 | 4 | 12.75–13.25 |
| E/F (force) | 2 | 300 | 5 | 12.50–13.25 |
| F/T (time) | 2 | 300 | 5 | 12.50–13.00 |
| A/C (particle) | 1 | 300 | 3 | 13.50–14.25 |
| C/E (particle) | 1 | 300 | 6 | 12.75–14.25 |
| A/E (particle) | 1 | 300 | 3 | 13.50–14.25 |
| B/C (force) | 1 | 300 | 6 | 13.75–14.00 |
| C/D (force) | 1 | 300 | 6 | 12.75–14.00 |
| E/F (force) | 1 | 300 | 7 | 12.75–14.25 |
| F/T (time) | 1 | 300 | 7 | 12.50–13.50 |
| A/C (particle) | 0 | 300 | 5 | 13.50–14.00; 14.25–14.50 |
| C/E (particle) | 0 | 300 | 8 | 12.75–14.75 |
| A/E (particle) | 0 | 300 | 5 | 13.50–14.00 |
| B/C (force) | 0 | 300 | 8 | None |
| C/D (force) | 0 | 300 | 8 | 14.00–14.75 |
| E/F (force) | 0 | 300 | 9 | 14.00–14.75 |
| F/T (time) | 0 | 300 | 9 | 12.50–14.75 |

## Shape shifts inside the passing intervals

The screen does not test shape. The following are maximum absolute bin-median percentage shifts **inside its passing intervals**, with at least 30 objects for the shape statistic. They are not halo-to-halo scatter or shape acceptance limits. Both reported axis ratios include the inherited concentration-dependent empirical correction and transverse-axis reordering (PMP2linker.f90:1105–1117); they are not the uncorrected tensor ratios.

| Pair | z | max abs median delta(b/a), % | max abs median delta(c/a), % |
|---|---:|---:|---:|
| A/C | 2 | Unmeasured | Unmeasured |
| C/E | 2 | 1.25 | 1.94 |
| A/E | 2 | Unmeasured | Unmeasured |
| B/C | 2 | 1.10 | 0.15 |
| C/D | 2 | 1.62 | 3.49 |
| E/F | 2 | 1.30 | 1.67 |
| F/T | 2 | 2.53 | 5.42 |
| A/C | 1 | 2.56 | 8.11 |
| C/E | 1 | 1.67 | 1.19 |
| A/E | 1 | 2.52 | 8.92 |
| B/C | 1 | 3.12 | 5.08 |
| C/D | 1 | 1.89 | 4.87 |
| E/F | 1 | 2.11 | 5.45 |
| F/T | 1 | 1.94 | 4.70 |
| A/C | 0 | 2.09 | 6.24 |
| C/E | 0 | 1.31 | 1.59 |
| A/E | 0 | 2.71 | 6.12 |
| B/C | 0 | Unmeasured | Unmeasured |
| C/D | 0 | 3.91 | 5.97 |
| E/F | 0 | 4.11 | 5.86 |
| F/T | 0 | 2.92 | 4.91 |

In the completed suite, E/F at z=0 reaches 5.86% in median c/a inside its mass/Vmax passing interval. F/T reaches 4.70–5.42% across the three outputs; the table maximum is 8.92% for A/E at z=1. Shape, velocity, tails and scatter require their own criteria.

## Effective particle and publication cuts

The common mass cut is max(2.5e12 Msun/h, nominal_floor × max(m_particle)). Only whole bins above that cut are used. The table below uses the nominal 300-particle floor; the values do not depend on redshift in this campaign.

| Pair | Common mass cut, Msun/h | Minimum bound particles, left / reference | First whole-bin lower log10 mass |
|---|---:|---:|---:|
| A/C | 2.570666e+13 | 300 / 2400 | 13.50 |
| C/E | 3.213332e+12 | 300 / 2400 | 12.75 |
| A/E | 2.570666e+13 | 300 / 19200 | 13.50 |
| B/C | 3.213332e+12 | 300 / 300 | 12.75 |
| C/D | 3.213332e+12 | 300 / 300 | 12.75 |
| E/F | 2.5e+12 | 1868 / 1868 | 12.50 |
| F/T | 2.5e+12 | 1868 / 1868 | 12.50 |

For E/F and F/T, the publication mass cut corresponds to 1867.22 particle masses, hence at least **1868 bound particles**. Their first whole bin begins at log10 mass 12.50, requiring at least **2362 particles**. Nominal floors of 100, 300 and 1000 therefore give the same results for these pairs; they do not constitute three independent resolution tests.

## What the abundance diagnostic measures

The paired delete-one-octant jackknife measures variation of the abundance ratio across spatial omissions in this realization. It can legitimately be zero when both catalogues have identical octant counts. This does not establish zero ensemble uncertainty. The screen uses measured shifts and this diagnostic; it is not a confidence statement that the true abundance shift is below 5%.

For comparison, supplemental_diagnostics in convergence-assessment.json reports sqrt(1/Nleft + 1/Nright) as a **hypothetical independent-count relative scale**, and also that scale multiplied by the absolute ratio for comparison with delta n. The actual catalogues are paired and correlated, so their covariance must be included for a sampling-error interpretation. The independent-count scale is not used to replace the jackknife, bound the actual uncertainty, or change the acceptance decisions.

C/E at z=1, log10 mass 14.00–14.25 has identical octant counts [10,9,7,5,2,6,4,4]: 47 objects in each catalogue, paired sigma 0, and independent-count relative scale 20.63%. This is one unique bin repeated at three nominal floors. At floor 300, 21 bins meet the abundance-difference cut while paired sigma exceeds 2.5%; 15 of those bins pass the full screen. High-mass abundance agreement has limited statistical discrimination.


## Redshift dependence and additional claims checks

Matched median bound-mass shifts change sign in eight E/F and C/D mass bins between z=2 and z=0. The following example uses the convention 100 × (coarse/reference − 1). These are separately matched populations at each epoch, not the same haloes tracked through time.

| Pair | log10 mass interval | z=2 median shift, % | z=1, % | z=0, % |
|---|---|---:|---:|---:|
| E/F | 13.00–13.25 | +1.85 | -0.17 | -4.07 |
| C/D | 13.00–13.25 | +2.36 | +0.08 | -3.80 |

The small z=1 medians lie near a sign transition in the population response. They do not establish stability across epochs or locate a continuous zero crossing. Both example bins pass the 5% mass condition at every epoch; their z=0 full-screen failures are abundance and resolved Vmax. The wider z=1 force intervals must therefore be read with all criteria, not attributed solely to mass cancellation. z=0 is the more restrictive lower-mass force comparison here; higher redshifts remain restrictive for timesteps.

The C/D and E/F force responses agree to 0.35 percentage points in median mass and 0.28 in median Vmax in the eight common usable z=0 bins. This supports approximate force–particle separability for those medians. The particle chain has median absolute residual 0.162 pp, maximum 0.942 pp over 18 bin/property cases (mass and Vmax combined); this is not a 0.2-pp bound. Interactions with timesteps and all-property separability remain untested.

The retained legacy/v3 catalogues allow 114 comparable abundance bins: v3 has the smaller absolute shift in 68, legacy in 43, with 3 exact ties. Median absolute shifts are 2.78% and 3.21%, respectively. These correlated measurements establish neither a statistically supported convergence advantage nor equivalence. Matched legacy properties still require missing membership information.

In C/E at z=0, log10 mass 12.75–13.00, the per-halo mass-shift 16th/median/84th percentiles are −8.30/−0.52/+6.74%. Median agreement does not bound individual errors. The review's 0.14–1.57% unmatched fractions pool broad mass ranges; per-bin losses and the both-particle-floor selection are recorded separately. Small selection fractions alone do not bound a median shift in percentage units.

Candidate counts respond strongly to particle load and force refinement at the fixed finder mesh. This demonstrates sensitivity but cannot exclude attenuation or a shared bias in catalogue statistics. The z=0 median apertures span 2.97–3.04 finder cells; similar radii are also affected by the common publication cut. Neither an absolute error nor its exact cancellation is measured.

See [CLAIMS-RESPONSE.md](CLAIMS-RESPONSE.md) and the source-bound `claims-response.json` for the reproduced controls, population definitions and follow-up design.

## Membership checks

The frozen checker applies the production **priority-ordered extended-aperture rule**: a lower-priority centre must lie outside the higher-priority halo's reported aperture. Priority is bound mass, then stable candidate index. This is the rule in PMP2linker.f90:696–750; it is not symmetric centre exclusion or a test using unextended SO radii. The zero fields are assertions that passed; the examined count counts neighbours already inside the queried aperture and is not a broad coverage measure.

| Catalogue | Published haloes | Exact duplicate sets | Repeated IDs within a halo | Mass/count mismatches | Priority-aperture violations | In-aperture priority pairs examined |
|---|---:|---:|---:|---:|---:|---:|
| A/z2 | 8198 | 0 | 0 | 0 | 0 | 0 |
| B/z2 | 7420 | 0 | 0 | 0 | 0 | 0 |
| C/z2 | 8395 | 0 | 0 | 0 | 0 | 0 |
| D/z2 | 8551 | 0 | 0 | 0 | 0 | 0 |
| E/z2 | 8435 | 0 | 0 | 0 | 0 | 0 |
| F/z2 | 8754 | 0 | 0 | 0 | 0 | 0 |
| T/z2 | 9019 | 0 | 0 | 0 | 0 | 0 |
| A/z1 | 19085 | 0 | 0 | 0 | 0 | 0 |
| B/z1 | 16626 | 0 | 0 | 0 | 0 | 0 |
| C/z1 | 19480 | 0 | 0 | 0 | 0 | 0 |
| D/z1 | 20394 | 0 | 0 | 0 | 0 | 0 |
| E/z1 | 19635 | 0 | 0 | 0 | 0 | 0 |
| F/z1 | 20777 | 0 | 0 | 0 | 0 | 0 |
| T/z1 | 20986 | 0 | 0 | 0 | 0 | 0 |
| A/z0 | 24898 | 0 | 0 | 0 | 0 | 0 |
| B/z0 | 20793 | 0 | 0 | 0 | 0 | 0 |
| C/z0 | 25811 | 0 | 0 | 0 | 0 | 0 |
| D/z0 | 27855 | 0 | 0 | 0 | 0 | 0 |
| E/z0 | 25913 | 0 | 0 | 0 | 0 | 0 |
| F/z0 | 28226 | 0 | 0 | 0 | 0 | 0 |
| T/z0 | 28224 | 0 | 0 | 0 | 0 | 0 |

The independent review found four reverse-orientation centre-in-aperture pairs: the higher-priority centre lies inside the lower-priority aperture. This orientation is allowed by the production rule. All four separations exceed both unextended SO radii. They are not host-rule failures. The production data provide no positive control of the violation branch; fixtures remain necessary.

Exact member-set uniqueness does not require disjoint memberships. The review measured excess memberships (occurrences after the first appearance of a particle ID) at 0.06–0.13% of total memberships. Summed bound mass counts shared particles repeatedly; it is not a partition of the particle set. This global rate does not bound the fractional mass error of an individual halo. See [the independent review](../../analysis/review-20260909-convergence/REVIEW.md) for the full scan.

## Limits of the measurement

- The finest simulation is a comparison reference, not an independent physical truth.
- One matched realization isolates numerical changes but does not measure box-size or cosmology dependence.
- At z=0, particle refinement C/E and timestep refinement F/T meet the working screen over broader mass ranges than force refinement E/F. The force mesh is the limiting tested setting for lower-mass z=0 haloes in this suite.
- The z=0 timestep result does not extend to all higher-redshift masses; inspect the separate z=1 and z=2 ranges.
- The force response changes sign with redshift in several fixed mass bins; the wider z=1 intervals are epoch-specific. Small mass medians near this transition do not establish convergence across epochs.
- The initial E/F positions differ by at most 6.103515625e-5 Mpc/h; physical velocities match exactly. The native periodic edge guard is retained. F/T initial positions match exactly.
- Native output velocities are staggered by half a timestep. The F/T velocity difference includes this output-time effect.
- The normal schedule has 158 steps and T has 316. All normal endpoints and output epochs are retained.
- Shape and velocity differences are reported separately; the displayed mass/Vmax screen does not certify every halo property.

Source and input receipt, membership, index and density hashes are recorded in convergence.json. Completed job allocations and billing are recorded in accounting.json.
