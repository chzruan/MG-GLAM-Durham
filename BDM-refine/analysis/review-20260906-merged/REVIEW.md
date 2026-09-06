# Independent review of the merged BDM changes, `a8c7715..cz`

Reviewed 2026-09-06 against `cz` at `0fa6a5b`. Scope: `PMP2linker.f90` (+1456/-1234),
`makefile`, and the audit, repair and 1024^3 validation evidence merged with them.
Root checkout was read-only except for the new tests recorded in
[tests/halo_gaps.md](../../repairs/20260906/tests/halo_gaps.md).

The repairs are sound. The catalogue-v2 contract, the SO/unbinding/population
rework, exact periodic image reconstruction, exact PM state restoration and
atomic publication all hold up under reading, and the accompanying evidence is
unusually complete and honestly caveated. Everything below is either a residual
defect, an unmeasured cost, or a claim whose strength differs from how it reads.

## Findings

### P1 — SO threshold normalisation disagrees with the particle-mass normalisation (medium)

`GetHalo` uses `threshold = 1.150d12*Om0*Ovdens` (`PMP2linker.f90:1126`), which
implies `rho_crit,0 = 1.150e12/(4*pi/3) = 2.7454e11`. `SetParameters` builds the
particle mass from a different constant, `Dscale = 2.774e+11*(Box/NROW)**3`
(`:2020`). The simulation's own mean matter density is `NROW^3*MassOne/Box^3 =
Om0*2.774e11`, so the finder applies

    Delta_effective = 0.98970 * Ovdens

relative to the box it is measuring, about **1.04% below the overdensity printed
in the catalogue header**. Propagated through the crossing condition this is a
few per mille in `Rso` and in bound mass, depending on the local profile slope.

Both constants are inherited, so this does not contaminate the old/new
comparison. It matters now because the crossing itself has been made exact: this
is the largest remaining systematic in the SO definition, and it is the only one
that is free to remove. The threshold can be derived rather than hard-coded:

```fortran
threshold = 4.1887902d0*dble(Ovdens)*dble(MassOne)*dble(NROW)**3/dble(Box)**3
```

which is exact by construction. This changes catalogue values, so it belongs in a
deliberate, versioned change with its own validation, not in a cleanup.

### P2 — the outermost-crossing convention is self-limiting (informational)

`halo_review.md` records that the incidence of multi-crossing profiles in
cosmological catalogues was not estimated. The new pair scan closes that with a
bound rather than a survey. For two equal clumps, the outer root is
`2^(1/3)*Rso`, and the companion is absorbed only while its centre lies within
`1.135*Rso`; measured transition between 1.1 and 1.2 (see the test note).

Generally, a second crossing at `R2` requires `N(<R2)*m = threshold*R2^3`, so the
selected radius is always the one its own enclosed population supports and can
never exceed `(M_enclosed/threshold)^(1/3)`. Any companion it absorbs therefore
has its centre inside the resulting aperture, where the distinct-host rule
already excludes it. The convention cannot inflate a radius arbitrarily. This
belongs in the catalogue contract next to the existing SO wording.

### P3 — unbinding has no pass cap and no measurement (medium)

The SO contraction caps at sixteen passes and finishes with a sorted interval
scan (`:1151-1185`), explicitly so that an adversarial profile cannot make it
quadratic. The unbinding loop (`:1209`) kept only the trivial `n+1` bound and
nothing records how many passes actually occur.

A constructed energy ladder needs **91 passes for 200 rows and 191 for 400**
(one antipodal pair removed per pass), i.e. passes ~ n/2 and O(n^2) work. At
1024^3 the largest bound sets reach ~10^5 rows, so one pathological candidate can
dominate `ParametersDistinct` with no diagnostic. Recommend a pass cap (32 is
generous next to the SO stage's 16), a sixth `HaloStatus` bit, and printing the
maximum pass count alongside the five existing status counts. The fixture is not
a claim about cosmological incidence — it shows the bound is real and unmeasured.

### P4 — the published-sample floor moved from 10 to 20 bound particles (low, documentation)

`WriteFiles` now cuts at `aNpart < 20.` (`:590`); the pre-audit writer cut at 10.
Immaterial for the validated configuration, where `MassMin = 2.5e12` corresponds
to 233 particles, but material for any run with a small `MassMin`, and it is not
in the repairs table or the catalogue contract. `publication_regression.py` pins
the new boundary (`Mvir = 20*MassOne` must publish) without stating the change.

### P5 — the production-scale host check is a null assertion (informational)

`validation-summary.json` reports `host_exclusion_violations: 0` at all three
redshifts, but also `higher_priority_neighbours_examined: 0`. The KD-tree found
no lower-priority centre inside any published `Rvir`, so the violation branch in
`run_validation.py:509` never executed. That is the expected outcome of a correct
implementation and the check is definitionally aligned with production, but the
1024^3 run demonstrates no ability to fail. The positive controls are at unit
scale, where `halo_host_ties.py` reproduces the three native equal-mass failures.
The headline should be read as "no published pair violates the rule", not as an
independent production-scale test with a demonstrated failure mode. Injecting a
known violating pair into a copy of the membership tape would cost minutes.

### N1 — `BDM_FFLAGS` silently loses `-fp-model precise` under an FFLAGS override (medium)

```makefile
BDM_FFLAGS = $(subst -fp-model fast=1,-fp-model precise,$(FFLAGS))
```

is a recursive variable, so `$(FFLAGS)` re-expands at use. The repository's own
`PMP2main-bitmatch` target overrides `FFLAGS` with `FFLAGS_BITMATCH`, which
contains no `-fp-model` at all, so the substitution matches nothing. Verified:

```
$ make -n FFLAGS="-O2 -g -qopenmp" PMP2linker.o
ifx -O2 -g -qopenmp -c PMP2linker.f90
```

The finder is then built with ifx's default fast FP model, which is exactly what
the rule exists to prevent. Appending instead of substituting,
`BDM_FFLAGS = $(FFLAGS) -fp-model precise`, is override-proof provided ifx
resolves a repeated `-fp-model` in favour of the last one — documented
behaviour, but not verified here, since ifx was not loadable in the review
shell. If that is not wanted as an assumption, list the finder's flags
explicitly instead of deriving them from `FFLAGS`. Low impact today (bitmatch
targets the PM solver), but the guarantee is silent when it lapses. The explicit
rule also drops the `-w` used by the `.f90.o` suffix rule; harmless, but
gratuitous.

### N2 — the seven density-sensitive rows are fully explained, with one stated residual

The explanation can be closed structurally, not just empirically. `FI` is read
only inside `FindMaxima` (`:1626-1706`) and is deallocated immediately after
(`:154`), so exactly two channels carry the density field into the catalogue:
the peak set, and `Xoff(slot)=FI(...)` at `:1666`, whose only consumer is the
centring seed radius at `:1548`. With the control's 794,570 peak indices and
positions identical between the two fields, the seed radius is the *only*
remaining path, and the seven changed rows are what that path produces. That is
a complete causal account, and stronger than the log's row-matching evidence.

The residual is that peak-set identity is itself density-dependent, through two
exact float comparisons: `FI > Ovdens/3` (`:1631`, `:1663`) and the strict
`>`/index tie rule (`:1705-1706`). The observed `max|dFI| = 0.0664` is far larger
than a ULP near the cut (`Ovdens/3 ~ 113` for this cosmology), so peak-set
stability was *observed here*, not guaranteed in general. That should be said
explicitly wherever the seven rows are described as fully understood.

The cheapest way to remove the last channel is to quantise the seed radius — it
only interpolates between `0.5*Cell` and `2*Cell`, so rounding
`log10(max(Xoff,0)+10)/2` to a coarse grid makes ULP-level `FI` noise unable to
move any centre. Given a reproducible peak set the finder would then be
reproducible, and these seven rows would not exist. It changes catalogue values,
so it is a versioned change.

### N3 — `List` anti-scales; measured, and ~10% of the finder at 64 threads (medium)

`particles_README.md` discloses that the z-slab prefilter keeps O(T*Np) outer
scans and makes no linear-work claim, but the production cost was never measured.
The fixed-density control ran six calls on the same immutable field, three at
each thread count, which isolates it exactly:

| Stage | 32 threads (s) | 64 threads (s) | Ratio |
| --- | --- | --- | ---: |
| `ParametersDistinct` | 142.04, 144.39, 144.51 | 82.38, 80.11, 78.36 | 1.80x faster |
| `List` | 8.66, 8.80, 9.23 | 10.94, 11.11, 10.05 | **0.80x, i.e. 24% slower** |

Doubling the threads makes `List` slower, which is the O(T*Np) signature: every
thread reads every row's z coordinate through `BdmParticleCoordinate`, including
an `OriginalParticleId` indirection, and discards ~63/64 of that work. It is
10.9 s of a ~109 s finder now, and on a full 128-core node it would keep growing
while everything else halves — roughly a quarter of the finder.

An output-identical fix: one parallel pass writing each row's z-cell index into
an `integer*2` (or `integer*4`) scratch array, then the existing slab loop reads
only that. Traffic per row per thread drops from an indirection plus two float32
reads to 2-4 contiguous bytes; the ascending-`jp` link order is untouched. Given
the COSMA CPU-efficiency rules this is worth doing before the next production
campaign.

### N4 — `AddBuffer` peak memory and serial prefix (low)

`imageEnd` is `integer*8` of length `originalCount` — 8.6 GB at 1024^3 — and its
prefix sum is serial over 1.07e9 entries. Per-row image counts never exceed 8, so
a 1-byte count array plus a blocked parallel prefix removes ~8 GB of the peak and
part of the ~4 s. Correctness is not affected; the exact-count design is right.

### N5 — `MaxMemory` is unreachable from configuration (low)

Refusing to coarsen `Cell` in `SizeList` is correct now that `Cell` sets the
physical search radii — silently changing the physics to fit an allocation
estimate was the worse behaviour. But `MaxMemory = 500` is a source parameter and
`ReadParameters` rejects unknown keys with `error stop`, so an operator who hits
`'BDM linked-list allocation exceeds configured memory limit'` has no
configuration path and must edit and rebuild. Add `maxmemory` to the accepted
keys and validate it like the others.

### N6 — dead and dormant code (low)

`RescaleCoords(iFlag /= 1)` is unreachable: the only call site passes 1 and `BDM`
calls `RemoveBuffer` directly. `GetProfiles`, `HaloProfile`, `WriteProfiles` and
`RemoveDuplicatesSimple` are documented as dormant, but would fault if revived —
`GetProfiles` indexes with `ih`, which is never set, and the `MassH1`-family
arrays are never allocated. An `error stop 'unsupported entry point'` at the top
of each would make the documented status enforced rather than advisory.

## Assumptions behind the 1024^3 evidence, tested

Deliberately trying to break the validation's own claims, three hold and two
should be stated more narrowly.

- **Bitwise finder determinism at fixed density**: holds. Six calls, two fields,
  strict byte equality within each group, plus bitwise particle-state
  restoration. The `List` link order is thread-count independent by construction
  (each cell has one writer; ascending `jp` traversal), which is consistent.
- **Membership integrity**: holds, and is the strongest part of the evidence —
  147,042 published rows tied to raw sorted original IDs, zero repeats, zero
  identical sets, mass/count consistency, 23 of 24 columns joined to raw values.
- **Host exclusion at production scale**: weaker than it reads (P5).
- **"Seven rows fully explained"**: true, and now provable from the code, but
  conditional on peak-set identity, which was observed rather than guaranteed (N2).
- **The 2.00x replay factor**: correctly scoped in the summary. The per-call
  breakdown attributes it: at 64 threads `ParametersDistinct` 78-82 s,
  `List` 10-11 s, `AddBuffer` 2.8-4.6 s, `FindMaxima` ~2 s,
  `RemoveDuplicates` 1.4-2 s, `WriteFiles` ~1.8 s. The cost is genuine new
  physics work (full-cap SO gather, sort, iterative unbinding); about 10% is
  recoverable in `List` (N3) and a few percent in `AddBuffer` (N4).

One cost that deserves a line in the operating notes: the periodic buffer is now
sized as `ParticleSearchRadius = 30.75*Box/NGRID` rather than the legacy 5 Mpc/h,
so the extra-row fraction is roughly `184.5/NGRID` — about 9% at NGRID=2048, 19%
at 1024, 40% at 512. Coarse-mesh configurations should check that before scaling.

## Tests added

[tests/halo_gaps.md](../../repairs/20260906/tests/halo_gaps.md) — 166 GNU
experiments, checked and optimized, no production change:

1. 144 duplicate-merge cases. `MergeNumericalDuplicates` previously had exactly
   one two-member identical set; now groups of three and four, two groups plus
   singletons, high/interleaved candidate indices, equal-length sets differing
   only in the last identity, a group whose lowest index is host-removed, and a
   candidate-count sweep crossing every bottom-up merge width including the
   `width > n/2` exit. All pass against the existing oracle. No defect.
2. Unbinding pass counts (P3).
3. Companion-absorption scan (P2).

Still uncovered, and cheap: a positive control for the production-scale host
check, and any production instrumentation of unbinding pass counts.
