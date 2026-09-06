# Configuration and deterministic peak repairs

Addresses audit findings F01, F03, F11, F12 and F15. The isolated branch changes
configuration, peak selection, maxima cleanup and catalogue metadata; it does
not repair halo membership or particle buffering.

`BDM` now reads and validates configuration and establishes scales before any
density/peak work. `SetOverdensity` supplies the same density threshold to
header creation and peak selection for modes 0, 1, 2 and 3. It retains the
existing flat-Lambda background convention and historical normalization 178.
The ASCII header now records the actual `OmL` and appends `[BDM finder v2]` to
its first line. There are still eight header lines and 24 data columns.

Configuration accepts compact assignments, mixed-case names, tabs, optional
leading/trailing comments and the existing parameter aliases. Invalid modes,
nonfinite values, nonpositive bin widths/profile counts, negative mass/radius
correction values, unknown names, malformed values and truncated long lines
fail with nonzero status before catalogue outputs are opened. The obsolete
`Nne` option remains accepted with an explicit notice that it has no effect.

Peak selection reads an immutable density field. A cell must exceed one third
of the selected overdensity and be no lower than all 26 periodic neighbours.
For equal densities, a neighbour with a lower global `(z,y,x)` index rejects
the cell. This local deterministic rule removes adjacent equal-density seeds;
it is not a watershed or a claim that every arbitrarily shaped connected
plateau has only one seed. A rectangular uniform plateau yields its lowest
index. Stable per-plane counts and prefix offsets give exact storage and
thread-independent catalogue candidate order with O(NGRID) scratch memory.
There are no arithmetic density markers or average thread-capacity estimates.
The threshold check precedes neighbour comparisons in both scans.

Zero peaks produce the normal header-only catalogue and return from BDM.
`ReleaseMaxima` cleans up the common zero/nonzero allocations. The empty path
leaves the density array allocated and does not rescale/buffer particles.

Validation:

```sh
micromamba run -n cosemu python3 -B BDM-refine/repairs/20260906/tests/peaks_regression.py
```

136 asserted experiments passed in GNU bounds/FP-check and O3 builds. Fixtures
cover sparse and empty fields, large dynamic range, plateaus, periodic ties,
load skew, a uniform mesh and random quantized density. Coordinates and original
peak density match an independent NumPy 26-neighbour oracle exactly and are
byte-identical at 1, 2, 4, 8 and 16 threads. All four density modes agree between
headers and peak selection; invalid configuration fails before output creation.
Repeated empty BDM calls preserve particle arrays exactly and retain FI.
`tests/peaks_results.json` records source hashes, commands and experiment output.
Temporary build/run directories are removed automatically.

Integration: the focused driver extracts `BDM`, `ReleaseMaxima`,
`ReadParameters`, `ConfigurationError`, `ValidateParameters`, `SetOverdensity`,
`SetParameters`, `FindMaxima`, `IsDensityMaximum` and `WriteFiles`. Empty-path
particle/halo routines are fail-on-use stubs, as documented in the test. Adapt
those stubs for any new BDM entry/exit helpers from the independent particle
repair. Add membership cleanup to `ReleaseMaxima` when integrating that work.
Native ifx builds, combined repaired-finder tests and small simulation replay
remain integration checks; no Slurm job was submitted by this branch.
