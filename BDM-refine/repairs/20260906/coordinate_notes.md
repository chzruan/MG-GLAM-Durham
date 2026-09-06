# Original-row coordinate shortcut

`BdmParticleCoordinate` now returns its stored coordinate when the original
particle ID equals the row ID. `BdmParticlePosition` loads the three original
coordinates directly in that case. Original rows occupy the leading block and
map to themselves; periodic images retain the existing double-precision integer
shift reconstruction. The change does not alter search radii, neighbour tests,
halo definitions or membership algorithms.

Four alternating baseline/optimized and optimized/baseline pairs were run after
one warmup of each executable, using the verified 128^3 GR snapshot and a single
pinned CPU. The native baseline source is
`e89b30162715478a8172bd23875cd13dfaa588a16b0c2eb2ea0bb309cddc10c6`.
Only the two coordinate helpers were replaced in the paired optimized build;
all other source, objects, snapshot input, configuration and runtime libraries
were shared.

| Metric | Baseline median | Shortcut median |
|---|---:|---:|
| Elapsed time | 5.599 s | 3.970 s |
| CPU time | 5.525 s | 3.920 s |

The measured median speedup is **1.410x**, or **29.1% less elapsed time**. Every
pair favoured the shortcut; paired speedups were 1.73x, 1.42x, 1.41x and 1.41x.
The first pair was slower overall, so the medians and three later pairs are the
useful estimate. This is a one-core result for this snapshot, not a general
scaling claim. It does not replace the required full-cap SO search.

All warmup, timed and instrumented runs wrote the identical 1,218-row catalogue
byte for byte. Separate baseline/optimized membership taps also wrote identical
raw original-ID sets and all 21 recorded raw properties. Twenty one-core GNU
and ifx checked/optimized periodic controls passed: original and ghost positions,
face/corner/list-boundary membership oracles, half-box coverage, centring,
Hubble kinetic energy and RMS radius. The existing periodic oracle driver was
used without changing its assertions.

Evidence: `coordinate-results.json`, `coordinate-periodic-results.json` and
`coordinate-replays.log.gz`. `coordinate_benchmark.py` records compiler commands,
binary/source/input hashes, per-run wall/CPU/MaxRSS and paired results. The
benchmark requires a copied native-build workspace and `copied-manifest.json`;
initialize Intel modules and `BDM_AUDIT_NATIVE_LIBS`, then invoke it with
`micromamba run -n cosemu python3 -B` and `--workspace <copied-workspace>`.
It checks copied hashes and pins all timed processes to one CPU. Failed timed
process groups are terminated rather than leaving the time wrapper's child.
The copied build and snapshot scratch were removed after evidence was retained.

To reconstruct the copied workspace from the verified repair archive, map
`work/equal-mass-discovery/native-build/` to `build/` and
`work/snapshots/n128/` to `snapshot/`. Write the `copied_manifest` object in
`coordinate-results.json` to `copied-manifest.json` in that workspace. The
driver verifies these original input hashes before compiling either variant.
