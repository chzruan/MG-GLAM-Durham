# Common analysis mesh for the convergence campaign

These campaign-only adapters run the exact legacy `a8c7715` finder or the
audited v3 finder on particles from any power-of-two evolution mesh, using a
common analysis mesh. The intended production target is 2048. Simulation and
production finder files are not changed. The legacy source is compiled with
its historical `-fp-model fast=1`; v3 and the adapter use `-fp-model precise`.
Both link the caller's same frozen runtime and density objects.

## Build integration

After building the 14 common objects listed in `build_adapter.py`, with their
Fortran module files alongside them, run:

```bash
module load intel_comp/2024.2.0 compiler-rt tbb compiler
export BDM_AUDIT_NATIVE_LIBS="$LD_LIBRARY_PATH"
micromamba run -n cosemu python3 -B replay/build_adapter.py \
  --common-dir /absolute/frozen/production \
  --v3-source /absolute/frozen/production/PMP2linker.f90 \
  --legacy-source /absolute/exact-a8c7715-linker.f90 \
  --output-dir /absolute/new-adapter-build
```

Alternatively import `build_adapter.build(common_dir, v3_source, legacy_source,
output_dir, compiler='ifx', env=None)`. The output directory must be new.
`build.json` records source, object, module, adapter and binary hashes and the
exact compiler commands. Source hashes are pinned to the audited v3 and exact
legacy files. Object/module hashes are checked again after the build. Module
files are copied into separate variant directories, so neither variant can
overwrite the caller's production modules.

The resulting executables are `legacy/replay.exe` and `v3/replay.exe`.
Execute in a fresh directory containing `BDM.config`, `CATALOGS/` and links to
the native `PMcrd.STEP.DAT` and `PMcrs*.STEP.DAT` files:

```text
replay.exe STEP TARGET_NGRID THREADS [write|read DENSITY_PATH]
```

With three arguments the entry calculates density normally. For a controlled
pair, run the first variant with `write /absolute/shared-fi.bin`, then the
second with `read /absolute/shared-fi.bin`. The second reads exactly the first
density field, avoiding OpenMP atomic accumulation differences. There is only
one full FI allocation per process, and no converted particle files or second
in-memory density copy. At analysis Ng2048 the temporary tape is 32 GiB plus
28 bytes.

The density tape is big endian: `>qqiff` stores analysis mesh, particle count,
snapshot step, expansion factor and box length, followed by float32 Fortran
planes. Reads validate this metadata and exact byte length. The calling
campaign must hash the original snapshot/configuration, density tape and
outputs; a completed writer exit and successful checks are required before
reuse. A killed writer may leave an incomplete, unreceipted tape; the adapter
never overwrites it. After both validated outputs and their receipts are
durable, the campaign can unlink the tape. This entry does not implement a
campaign resume or cleanup policy.

## Physical meaning and rounding

PM positions use a one-based mesh origin. The unit change is

```text
x_analysis = 1 + (x_evolution - 1) * Ng_analysis / Ng_evolution
v_analysis = v_evolution * Ng_analysis / Ng_evolution
```

Thus `(x-1)*Box/Ng` and `v*100*Box/(Ng*a)` describe the same physical particle.
Coordinates are canonicalized periodically, including an exact upper-image
coordinate `Ng+1` and a converted coordinate rounded onto that image. The
original files, particle order, NROW, box, cosmology and epoch stay unchanged.

Power-of-two velocity scaling is exact for ordinary normal float32 values.
The one-based position addition can introduce float32 rounding, including when
coarsening the mesh. Each call measures the maximum periodic physical position
error and velocity error and prints them. The position error must stay below
one float32 ULP at the analysis grid's upper boundary, converted to physical
units. This conservative bound is 3.052e-5 Mpc/h for L256/Ng2048. The adapter
does not claim bitwise invariant physical coordinates. It rejects invalid
coordinates, nonfinite values and overflow; underflow follows the documented
native compiler flags, including `-ftz`.

This common mesh holds finder peak detection, centring and search scales fixed
while the evolution mesh changes. It does not remove finite-particle effects.
Snapshots contain leapfrog velocities staggered by half their last timestep;
this adapter preserves that convention and does not synchronize velocities.

## Legacy behavior and diagnostics

No `ReadParameters` or `SetParameters` call is inserted before legacy BDM.
Its first-call configuration-order behavior is therefore retained. Supply a
legacy-compatible configuration with an initial comment line, spaces around
assignments and inline comments. Configuration differences must not be silently
removed when interpreting the legacy-to-v3 comparison.

V3 has only two added calls, immediately after successful catalogue writing
(ordinary and empty paths). The added routine writes `repair-members.bin` for
published haloes only, in the existing repair format:

- Header: int64 selected count, int64 candidate count, float32 particle mass.
- Each row: int64 candidate ID and member count, the existing 21 float32
  properties, then the sorted int64 original particle IDs.

The writer verifies its selected count equals `Nhalo`, stages to `.part`, and
publishes by rename. It creates no all-candidate raw unbinding tape. Existing
v3 stdout retains aggregate status, unbinding-pass and active-row diagnostics.

## Small controls

```bash
micromamba run -n cosemu python3 -B replay/test_adapter.py \
  --build-dir /absolute/new-adapter-build \
  --repo /absolute/MG-GLAM \
  --receipt /absolute/preflight.json
```

Controls cover power-of-two grids from 8 through 4096, exponent transitions,
both periodic images, invalid inputs, and physical velocity preservation. A
64-particle halo stored on evolution grids 8, 16 and 32 is analyzed on grid 16
using one shared density tape. Catalogues and v3 membership must agree byte for
byte within each variant; original snapshot hashes must stay unchanged. Both
variant outputs on the identity-grid fixture must also match their original,
unmodified standalone entries linked with the same runtime objects. All
temporary fixture files are removed after a successful receipt is written.
