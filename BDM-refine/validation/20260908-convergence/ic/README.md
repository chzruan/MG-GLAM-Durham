# Matched native initial conditions for the BDM convergence campaign

This is a **campaign-only, first-order (Zel'dovich) IC variant** generated from
`PMP2start.f90`. It leaves every production source and default untouched. Both
legacy and audited finder replays must consume the resulting identical PM
snapshots; this IC change is not attributed to BDM refinement.

## Why the same native seed is insufficient

Native `SPECTR` gives each packed Fourier k plane its own luxury RNG state, but
consumes NROW random draws per j row. Changing NROW therefore changes shared
modes. It also computes `ALPHA = AMPLT/sqrt(SUMM)*sqrt(8)` independently for each
realization and particle resolution: even matching Gaussian draws would not
match Fourier amplitudes. Finally, native `BLOCKS` shifts particles by half a PM
cell, which changes the physical origin when NGRID changes.

The campaign variant makes these controls explicit:

1. Every j row consumes **1024** Gaussian draws (the configured master NROW),
   retaining the native luxury seed stream and Box-Muller implementation. The
   first Box-Muller flag is explicitly zero. Smaller runs keep their own native
   spherical cutoff `|k| < NROW/2`; newly available shorter modes are added at
   higher particle resolution. Discarded coefficients still consume draws.
2. One NROW=1024 master determines ALPHA from its own native Setup amplitude and
   spectrum. Every consumer reuses that exact float32 value. Per-plane sums are
   accumulated in a fixed serial order, so the normalizer is thread independent.
   This retains native master realization normalization; it does not impose an
   independent realization normalization on coarser runs. The header AMPLT still
   records that run's Setup value, while the receipt records the actual ALPHA.
3. All terminal packed Fourier planes are zero. FFT5's final packed element is
   a Nyquist cosine, but the original map assigns it wave number zero. Allowing
   those planes creates spurious modes and inconsistent resolution cutoffs.
   Uniform modes are also zero. Nonterminal common modes preserve the master
   stream exactly; no attempt is made to reproduce the native Nyquist artifacts.
4. The fixed physical shift is `Box/(2*2048) = 0.0625 Mpc/h` in the L256 box.
   `BLOCKS` uses `xShift = 0.5*NGRID/2048`, so physical particle coordinates
   `(PMx-1)*Box/NGRID` equal lattice position minus displacement **plus** that
   common shift. Positions remain stored in native float32 PM coordinates.

The native FFT5 inverse is a cosine/sine series with **no NROW normalization
factor**. The independent analytic series and NumPy Fourier reconstruction in
`checks.py` verify this, including zero-coordinate planes. Native sqrt(2)
weights for zero-frequency coordinates remain in place.

The variant uses the native first-order velocity formula evaluated at
`a_velocity = AEXPN0 - ASTEP0/2`. For timestep repeat T, supply its actual stored
first step in Setup.dat; positions and displacement modes remain the same,
while initial velocity momenta use the correct different staggered time. This
retains the native high-redshift growth approximation, rather than changing
IC order or growth prescriptions during the convergence experiment.

## Build and run API

All commands below assume the repository root as working directory. Use an
isolated build directory outside tracked source. On COSMA, load Intel compiler
modules first, for example `intel_comp/2025.3.0 compiler-rt tbb umf compiler`.

```
micromamba run -n cosemu python3 -B BDM-refine/validation/20260908-convergence/ic/build.py \
  --repo . --work /path/to/campaign/build/ic --compiler ifx
```

The result is `PMP2start.matched.exe`, generated source, build.log and build.json
with exact compiler flags, input/generated-source hashes and executable hash.
The generated source sets `nbyteword=1`; Intel **must** use `-assume byterecl`
(the supplied builder does). GNU uses its default byte record lengths. Both
write the native big-endian PM format. Do not relink with mismatched RECL units.

Create native Setup.dat, PkTable.dat and TableSeeds.dat in the run's parent.
For the master (NROW=1024, normally NGRID=2048):

```
micromamba run -n cosemu python3 -B BDM-refine/validation/20260908-convergence/ic/configure.py \
  --run-dir /path/to/master/Run1 --master-nrow 1024 --origin-ngrid 2048 --realization 1
```

Run the executable **inside Run1**, passing realization `1` on stdin. Set
`OMP_NUM_THREADS` to the requested cores and `OMP_STACKSIZE=128M` or larger.
Missing explicit seed tables stop immediately; the slow native seed-table
creation fallback is disabled. The required generated namelist is:

```
&matched_ic
 ic_master_nrow=1024, ic_origin_ngrid=2048,
 ic_alpha=-1, ic_normalize_only=.false.
/
```

`ic_alpha=-1` is accepted only at the master NROW. After the master finishes its
spectrum, it writes `matched_ic_receipt.txt`; this is a normalization receipt,
**not a simulation/IC completion marker**. Confirm executable success and valid
PM outputs separately. For a consumer, including another NGRID or timestep:

```
micromamba run -n cosemu python3 -B BDM-refine/validation/20260908-convergence/ic/configure.py \
  --run-dir /path/to/consumer/Run1 --master-run-dir /path/to/master/Run1 \
  --master-nrow 1024 --origin-ngrid 2048 --realization 1
```

This verifies the master input/control hashes, realization/seed, epoch, box,
cosmology and byte-identical Pk/seed tables before copying ALPHA into the new
namelist. NROW, NGRID, native Setup amplitude and timestep may differ as intended.
A new realization requires its own master normalization. `--normalize-only`
creates a master scalar-normalization pass without particles or FFT arrays;
its receipt is otherwise compatible. It still evaluates the entire master
spectrum, so schedule that work appropriately.

Additional artifacts are deliberately small:

- `matched_ic_inputs.json`: input hashes and master receipt/control linkage.
- `matched_ic_receipt.txt`: exact round-trip ALPHA, spectrum sum, requested
  amplitude, seed, particle/PM sizes, origin, position/velocity epochs and factors.
- `matched_modes.bin`: one big-endian int32 extent, followed by three float32
  packed displacement-component cubes of that extent in Fortran order, before
  FFT. Extent is min(33,NROW-1). At production NROW>=256 every sampled mode lies
  inside every spherical cutoff, so these files should be byte identical across
  the suite. A small sample supplements the complete-mode small controls;
  it is not a full-volume production-mode survey.

Configuration publication never overwrites different existing controls or input
provenance. Identical files can be reused; missing files from an interrupted
configuration are published atomically without overwriting concurrent writers.

There is no `GetPower` call or allocation of the NGRID^3 density grid during IC
creation. Live main arrays are nine float32 NROW^3 arrays, **36 GiB** at 1024^3,
plus about 24 MiB of native output buffers per thread and small FFT/RNG scratch.
A 32-core, 64-GiB master IC allocation has reasonable headroom; measure its actual
MaxRSS and CPU use. Once evidence is retained, archive/remove isolated build
objects, compiler module files and temporary controls. No job is launched by
these utilities.

Two incidental native IC bookkeeping problems are made explicit in the
campaign-generated module: AEXP0/AU0 are initialized shared values for the
output header (the original scopes leave local header values undefined), and
the displacement log divides by Nparticles rather than a private loop index.
Neither changes particle physics. The FFT argument length is NROW, independent
of NGRID. No production source edit is required for any of these controls.

## Validation and scope

```
micromamba run -n cosemu python3 -B BDM-refine/validation/20260908-convergence/ic/checks.py \
  --repo . --compiler gfortran --output /path/to/checks-gnu.json
micromamba run -n cosemu python3 -B BDM-refine/validation/20260908-convergence/ic/checks.py \
  --repo . --compiler ifx --output /path/to/checks-ifx.json
```

Each compiler runs 13 small cases plus the actual standalone-entry PM I/O smoke
and input-provenance rejection controls. Tests check all packed modes and all
real-space particles at 8^3, 16^3, 32^3, three PM grids, 1/2/4 threads, a second seed,
1024-draw stride, no-FFT normalization, and half-step initial velocities. Real
fields agree with independent analytic trigonometric synthesis; shared complex
Fourier coefficients agree within float32 transform roundoff. Matching packed
coefficients and same-resolution threaded fields/particles are byte identical.
The PM I/O smoke verifies the native 529-byte header, initial epoch, big-endian
record layout and particle values after reading actual Setup/Pk inputs.

The test driver removes its complete temporary directory after success and
stores a compact JSON result. It does not assert cross-compiler bit identity,
resolution convergence, Gaussian-distribution quality beyond the inherited
RNG, or a completed full 1024^3 physical-mode check. Production master and
consumer receipts should be checked before evolving expensive runs.
