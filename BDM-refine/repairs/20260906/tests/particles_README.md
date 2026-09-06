# Particle preparation and centring repairs

This repair covers audit F02, F09, F10 and F16 in `PMP2linker.f90`. It also
corrects the linked-list z partition and removes repeated x/y index work.
The historical audit fixtures and results remain unchanged.

- Centring sums double-precision **local displacements** and velocities. The
  262,144-row symmetric cloud stays exactly centred at both 5.125 and 1000.125
  Mpc/h; the old far-cloud displacement was 6.55 Mpc/h. A one-row aperture has
  its measured centre/velocity; an empty aperture retains its last position
  and clears its mass, radius and velocity estimate.
- Periodic storage counts actual images before allocating, including x=0
  images at x=Box. Original rows remain first and every image carries the
  same `integer*8 OriginalParticleId` as its source. The output arrays are
  filled once and moved into place; oversized zeroing and a second final
  particle-array copy are eliminated.
- `RescaleCoords(1)` moves the six original PM allocations into retained
  storage, then creates the analysis workspace. `RemoveBuffer` frees that
  workspace and moves the originals back. Restoration is bitwise exact,
  including signed zero and subnormal input values, over repeated calls.
- Particle `Cell` remains `2*Box/NGRID`; it is never increased after preparing
  periodic coverage to satisfy a memory estimate. List storage estimates use
  exact int64 byte counts. The configured `MaxMemory` remains a coarse overall
  memory limit, not a measurement of Slurm's allowed RSS.
- `List` partitions **z** bounds with local thread counts and stable descending
  row-ID links. It computes x/y cell indices only for rows in the current slab.
  The z scans still perform O(T*Np) work; no linear-work scaling claim is made.

## Integration contract

`PrepareParticleSearch()` must run after `ReadParameters` and before
`SetParameters` writes catalogue headers; `AddBuffer` and `SizeList` also call
it defensively. It sets:

```
Cell                 = 2 * Box / NGRID
halfBox              = nearest(0.5 * Box, -1.0)
HaloSearchRadius     = min(15 * Cell, halfBox)
ParticleSearchRadius = min(HaloSearchRadius + 0.75 * Box/NGRID, halfBox)
dBuffer              = min(max(dBuffer, ParticleSearchRadius), halfBox)
```

`GetHalo` must use these shared limits and explicitly reject unresolved SO
crossings or expanded apertures beyond them. Every accepted spherical search
has radius strictly below half the periodic box, so two images of the same
original particle cannot belong to it. This is a supported-domain restriction,
not a claim to measure haloes larger than half the box. `OriginalParticleId`
exists from `AddBuffer` until `RemoveBuffer`; per-halo membership storage must
retain the **IDs**, independently of that workspace's lifetime.

The BDM allocation/free hooks for int64 `Lst` and `Label` must pass **twice**
the element count to `Memory`, whose argument counts four-byte words. The
particle preparation routines account for their own allocations with int64
word counts.

## Verification

Run from the repository root:

```sh
micromamba run -n cosemu python3 -B BDM-refine/repairs/20260906/tests/particles_tests.py
```

Native ifx verification uses the matching runtime outside Python's libraries:

```sh
export LINES=40 COLUMNS=120
module purge
module load intel_comp/2024.2.0
module load compiler-rt tbb compiler
export BDM_AUDIT_NATIVE_LIBS="$LD_LIBRARY_PATH"
micromamba run -n cosemu python3 -B BDM-refine/repairs/20260906/tests/particles_tests.py \
  --compiler ifx --output BDM-refine/repairs/20260906/tests/particles_results_ifx.json
```

Both GNU and ifx pass **104 experiments each**: 13 cases, checked and optimized
builds, and 1/2/4/8 threads. The ifx optimized build uses the repository's
AVX2/FMA, fast=1, FTZ and unroll choices; this is an extracted-routine test,
not a complete native finder replay. GNU checks include full bounds checks and
invalid/zero/overflow FP traps. ifx uses bounds checks and precise arithmetic
in its checked variant.

The cases include independent double-precision centring and velocity oracles,
exact face/corner image counts, minimum-image membership oracles at all eight
box corners, a radius-six search with a five-Mpc default buffer, a box smaller
than that default, exact full-array restoration through two cycles, and a serial
row-order oracle for every link and head in an asymmetric list domain.
See [particles_results.json](particles_results.json) and
[particles_results_ifx.json](particles_results_ifx.json) for compiler commands,
source/test hashes and complete process output.

The paired 128^3-particle list benchmark performs ten rebuilds per sample and
three alternating baseline/repaired samples per thread count. Every sample
also checks every link against the independent oracle. The latest median
baseline/repaired ratios are approximately 1.08, 1.10, 1.27 and 1.85 at 1, 2,
4 and 8 threads; short shared-host timings vary, and these are kernel results,
not production-finder speedups. An earlier tested linear-work bucket prototype
regressed 2/4-thread timings and was removed. Reproduce the retained simple
prefilter comparison with:

```sh
micromamba run -n cosemu python3 -B BDM-refine/repairs/20260906/tests/particles_benchmark.py
```

Full samples and commands: [particles_benchmark_results.json](particles_benchmark_results.json).
The scripts remove their temporary builds and process directories automatically.
