# BDM convergence measurements

**Complete seven-run measurement.**

All primary catalogues use the same 2048^3 analysis mesh. The simulation box is 256 Mpc/h; outputs are z=2,1,0. Particle resolution, evolved force resolution and timestep size are compared separately.

**IC scope:** all seven runs use native GLAM first-order (Zel'dovich) initial conditions at z_init=100. This completed suite is an existing-run exception to the default 2LPTIC workflow for new simulations; it does not establish convergence with 2LPTIC initial conditions.

Bound mass is the original-member count times the stored particle mass, evaluated in float64. Reported aperture mass and radius include the empirical Rext expansion. Positive Vmax values are compared only when both haloes resolve them; unresolved values are counted separately. Halo matches require mutual best shared-lattice membership overlap, with at least 50% in each object.

The table applies a descriptive working screen: at least 30 objects in each required bin; abundance within 5% with jackknife sigma at most 5%; median bound mass within 5%; median resolved Vmax within 2%; reference completeness at least 90%. These choices are adjustable. They do not certify those absolute accuracies. Scatter, individual criteria and 100/300/1000-particle floors are retained in the JSON results.

| Comparison | z | Particle floor | Eligible bins | log10 mass intervals meeting working criteria |
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

## Membership checks

| Catalogue | Published haloes | Exact duplicate sets | Repeated particle IDs | Mass/count mismatches | Host violations | Higher-priority neighbours examined |
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

A zero host-violation count is not a production positive control when no eligible higher-priority neighbours were examined. The independent small fixtures exercise the failure branch.

## Limits of the measurement

- The finest simulation is a comparison reference, not an independent physical truth.
- One matched realization isolates numerical changes but does not measure box-size or cosmology dependence.
- At z=0, particle refinement C/E and timestep refinement F/T meet the working screen over broader mass ranges than force refinement E/F. The force mesh is the limiting tested setting for lower-mass z=0 haloes in this suite.
- The z=0 timestep result does not extend to all higher-redshift masses; inspect the separate z=1 and z=2 ranges.
- The initial E/F positions differ by at most 6.103515625e-5 Mpc/h; physical velocities match exactly. The native periodic edge guard is retained. F/T initial positions match exactly.
- Native output velocities are staggered by half a timestep. The F/T velocity difference includes this output-time effect.
- The normal schedule has 158 steps and T has 316. All normal endpoints and output epochs are retained.
- Shape and velocity differences are reported separately; the displayed mass/Vmax screen does not certify every halo property.

Source and input receipt, membership, index and density hashes are recorded in convergence.json. Completed job allocations and billing are recorded in accounting.json.
