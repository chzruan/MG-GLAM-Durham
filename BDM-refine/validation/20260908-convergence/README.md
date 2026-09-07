# BDM v3 convergence campaign, 8 September 2026

Branch: `validation/bdm-convergence-20260908`, based on `cz` at
`e289c3f7ea28d38390a1663e50bb2784b85abc84`. The user authorized the seven-run
campaign, simulations below 1200^3 particles, and removal/consolidation of
unnecessary working artifacts. Production finder source is frozen. This is a
convergence measurement, not an assumption that the finest run is exact.

All runs use GR, L=256 Mpc/h, Omega_m=0.3089, Omega_Lambda=0.6911,
h=0.6774, sigma8=0.8159, z_init=100, and exact outputs z=2,1,0.
`MG_flag=0` and `MG_model=3` are retained together. The master-normalized
initial Fourier modes are shared across particle resolutions; extra resolved
short modes are allowed in finer runs. The origin and Fourier conventions
are controlled explicitly; using the same native random seed is insufficient.

| ID | Particle count | Evolution mesh | Comparison |
|---|---:|---:|---|
| A | 256^3 = 16,777,216 | 2048^3 | Low particle resolution |
| B | 512^3 = 134,217,728 | 1024^3 | Coarse force mesh |
| C | 512^3 | 2048^3 | Intermediate particle resolution |
| D | 512^3 | 4096^3 | Fine force mesh |
| E | 1024^3 = 1,073,741,824 | 2048^3 | High particle resolution |
| F | 1024^3 | 4096^3 | Finest reference |
| T | 1024^3 | 4096^3 | Each F step subdivided into two |

Primary halo comparisons use a **common 2048^3 analysis mesh**, with identical
physical mass cuts and configuration, so changing the force mesh does not
simultaneously change the finder mesh. Selected additional replays can measure
finder-mesh sensitivity separately. Native Ng/Np=8 finder meshes (A and D)
can generate many discreteness seeds; their cost requires a dedicated pilot.
Legacy means `a8c7715`, not the intermediate v2 normalization.

The normal exact-redshift schedule has 158 steps; T has 316, with identical
normal endpoints and output epochs. Simply halving the adaptive input step
would instead yield 312 steps with a different schedule. T starts with its own
correctly staggered velocities. Output velocities remain at a_out-da/2 in the
native PM convention, so timestep differences include this output staggering.

The planned comparisons are A/C/E (particle mass), B/C/D and E/F (evolved
force resolution at common finder mesh), and F/T (time). The earlier
L512/N1024/Ng2048 realization is a separate volume comparison, not a matched
same-object pair. Halo abundance, matched masses/radii/Vmax/velocity/shape,
completeness and independent membership invariants will determine usable
mass/redshift ranges. Agreement must be measured above particle-count floors
and with finite-volume/counting uncertainty.

## Resources and storage

The prior N1024/Ng2048 job 11948372 used 64 shared cores, 3477 s,
61.81 billed core-hours, 52.10 CPU hours and 166.62 GiB batch MaxRSS.
The initial seven-run estimate was 1600-2800 billed core-hours including
replays and pilots; the planning allowance is 3000. This is not measured
4096-grid performance. New pilot accounting determines final allocations.
Both expected and wall-limit core-hours are recorded before submissions.

Use `cosma8-serial`, actual thread counts and explicit memory. A 4096^3
float32 field alone occupies 256 GiB; a saved-density double-mesh replay
would need over 560 GiB before finder workspace at N1024. The primary replay
adapter avoids that extra full mesh. Full-node use requires a measured reason.

Starting user quota: 241104 files / 400000 soft / 440000 hard; space
23472442148 KiB / 32212254720 KiB soft. Keep one frozen build and reference
input set, small receipts, and scientific snapshots. Verify archives before
deleting owned build/fixture scratch. Preserve unrelated working directories.
Figures, Beamer source and PDFs are local artifacts, excluded from Git.

Run every Python command through `micromamba run -n cosemu python3 -B`.
Receipts and the final report distinguish preparation, successful process
completion, scientific validation and actual convergence; a submitted job
does not count as completed validation.
