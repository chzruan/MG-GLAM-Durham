# BDM convergence measurements

**Partial campaign: conclusions await the missing runs.**

All primary catalogues use the same 2048^3 analysis mesh. The simulation box is 256 Mpc/h; outputs are z=2,1,0. Particle resolution, evolved force resolution and timestep size are compared separately.

Bound mass is the original-member count times the stored particle mass, evaluated in float64. Reported aperture mass and radius include the empirical Rext expansion. Positive Vmax values are compared only when both haloes resolve them; unresolved values are counted separately. Halo matches require mutual best shared-lattice membership overlap, with at least 50% in each object.

The table applies a descriptive working screen: at least 30 objects in each required bin; abundance within 5% with jackknife sigma at most 5%; median bound mass within 5%; median resolved Vmax within 2%; reference completeness at least 90%. These choices are adjustable. They do not certify those absolute accuracies. Scatter, individual criteria and 100/300/1000-particle floors are retained in the JSON results.

| Comparison | z | Particle floor | Eligible bins | log10 mass intervals meeting working criteria |
|---|---:|---:|---:|---|
| A/C (particle) | 2 | 300 | 1 | None |
| B/C (force) | 2 | 300 | 4 | 13.25–13.50 |
| B/C (force) | 1 | 300 | 6 | 13.75–14.00 |
| A/C (particle) | 0 | 300 | 5 | 13.50–14.00; 14.25–14.50 |
| B/C (force) | 0 | 300 | 8 | None |

## Membership checks

| Catalogue | Published haloes | Exact duplicate sets | Repeated particle IDs | Mass/count mismatches | Host violations | Higher-priority neighbours examined |
|---|---:|---:|---:|---:|---:|---:|
| A/z2 | 8198 | 0 | 0 | 0 | 0 | 0 |
| B/z2 | 7420 | 0 | 0 | 0 | 0 | 0 |
| C/z2 | 8395 | 0 | 0 | 0 | 0 | 0 |
| B/z1 | 16626 | 0 | 0 | 0 | 0 | 0 |
| C/z1 | 19480 | 0 | 0 | 0 | 0 | 0 |
| A/z0 | 24898 | 0 | 0 | 0 | 0 | 0 |
| B/z0 | 20793 | 0 | 0 | 0 | 0 | 0 |
| C/z0 | 25811 | 0 | 0 | 0 | 0 | 0 |

A zero host-violation count is not a production positive control when no eligible higher-priority neighbours were examined. The independent small fixtures exercise the failure branch.

## Limits of the measurement

- The finest simulation is a comparison reference, not an independent physical truth.
- One matched realization isolates numerical changes but does not measure box-size or cosmology dependence.
- The initial E/F positions differ by at most 6.103515625e-5 Mpc/h; physical velocities match exactly. The native periodic edge guard is retained. F/T initial positions match exactly.
- Native output velocities are staggered by half a timestep. The F/T velocity difference includes this output-time effect.
- The normal schedule has 158 steps and T has 316. All normal endpoints and output epochs are retained.
- Shape and velocity differences are reported separately; the displayed mass/Vmax screen does not certify every halo property.

Source and input receipt, membership, index and density hashes are recorded in convergence.json. Completed job allocations and billing are recorded in accounting.json.
