# Halo property repairs: F05–F08 and GetHalo portions of F13

The halo-property path now determines spherical-overdensity (SO) radii from the
sorted particle radii. It chooses the outermost constant-enclosed-mass interval
containing the exact SO root, with at least ten enclosed particles. This removes
the shifted logarithmic-bin interpolation and makes the answer independent of
`dLogR`. The normalization remains `1.150e12 * Om0 * Ovdens`, now evaluated in
double precision.

Unbinding starts with particles inside the unextended SO radius. Each pass
recomputes the peculiar bulk velocity and the isolated spherical Newtonian
potential from its current survivors, excludes each particle's self term, and
removes particles with positive energy. Membership only decreases; each
non-final pass removes at least one particle, so the process terminates within
the initial population plus one passes. A successful final pass checks every
survivor against the same final potential and velocity. The peculiar velocity
and Hubble term retain the existing units/convention:
`v - <v>bound + 100*a*sqrt(Om0/a^3 + 1 - Om0)*(x - centre)`.

The potential uses the discrete spherical pair kernel
`G*m/max(ri,rj)/a`, excluding self. Its positive potential-energy magnitude is
`G*m*m/a * sum((i-1)/ri)` in ascending radius order. This includes the inner
shell and counts every distinct pair once. It is the exact discrete **spherical
approximation**, not an exact potential for a general asymmetric particle set.
The direct geometric pair oracle in the tests measures that approximation
separately from the unbinding and energy-discretization fixes.

The following catalogue definitions change deliberately:

- `Mvir` is final bound mass inside the unextended SO sphere. The enlarged Rext
  aperture no longer supplies potential from particles excluded from this
  population.
- Bulk velocity, kinetic/potential energy, reduced inertia tensor, RMS radius,
  centre offset, angular momentum/spin proxy, and circular-speed profile all
  use those same bound particles. The legacy spin and empirical shape proxy
  forms remain, with the bound mass in the spin denominator.
- The Rext correction itself is preserved. Reported `Rvir` remains the enlarged
  aperture and `Mtotal` counts every aperture particle; these fields should not
  be relabelled as uncorrected SO radius/mass.
- Vmax is the maximum of the final bound cumulative particle profile. The
  existing resolution correction now uses that maximum's saved radius and the
  same `G=4.333e-9` as the binding calculation. If the maximum lies at or below
  `0.1*Box/NGRID`, both Vmax and Rmax are explicitly zero (unresolved).

`HaloStatus` records insufficient SO population (1), an SO/aperture search that
exceeds the supported domain (2), unresolved Vmax (4), an empty final bound set
(8), or multiple particles exactly at the centre (16). The last case has an
undefined unsoftened spherical potential, so the candidate is rejected without
inventing a softening length. Rejected candidates have zero bound mass; all
catalogue-facing arrays receive finite initial values on every call.

`BdmHaloMembership` and `BoundParticleIds(:)` retain the ascending int64 original
particle IDs of final survivors. Candidate-local arrays contain only gathered
neighbours. No full-Np workspace is allocated per halo. The outer membership
array is initialized before the parallel `ParametersDistinct` loop; each
candidate owns its component. The integrated duplicate pass must free this
workspace after exact membership comparisons. The particle repair supplies
`OriginalParticleId`, `HaloSearchRadius` and `ParticleSearchRadius`; local tests
inject only these declarations if that branch has not yet been integrated.

Run the regressions with:

```
micromamba run -n cosemu python3 -B BDM-refine/repairs/20260906/tests/halo_regression.py \
  --output BDM-refine/repairs/20260906/tests/halo_results.json
```

The tests extract the actual production routines and run checked/FP-trapping
and optimized GNU builds. Independent oracles cover SO spacing and Rext,
iterative spherical membership, direct geometric pair energies for cold/mixed
controls, bound-only kinetic/bulk/RMS/spin/offset statistics, centre particles,
unresolved Vmax, singular and clipped-domain rejection, reversed original-ID
mapping and 1/2/4-thread candidate equality. Temporary build/run artifacts are
removed automatically. `halo_results.json` records source/test hashes and all
checks. Full native finder replays are coordinated on the integration branch.
