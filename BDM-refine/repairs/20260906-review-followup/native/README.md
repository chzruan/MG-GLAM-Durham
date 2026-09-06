# Controlled native follow-up replays

`build_native.py` freezes all 43 top-level Fortran/header/makefile inputs from
`--repo`, builds the actual production `make -j1 PMP2BDM` target, and builds
three diagnostic executables from the same 14 common objects. It checks those
common sources against commit `800eaac76ef7fe78e6a392c85c3a88822b1c9135`.
Source, object, module, executable and diagnostic source hashes, every compiler
command, and small native preflight results are recorded in `native/build.json`
under `--root`. Existing receipts and build directories are rejected, so a
failed build needs a new output root. Build artifacts remain in
`work/native-build` for provenance and later verified cleanup.

The variants are:

* `reference`: commit `800eaac` (v2 physics with unbinding work counters).
* `optimized-v2`: final source with exactly the SO threshold expression and
  catalogue version marker changed back to v2 in a diagnostic copy.
* `v3`: final source with its actual-box mean-density SO normalization.

The generated finder copies add a diagnostic flag, two guarded calls immediately
after successful `WriteFiles`, and the routines in `diagnostics.f90`. The flag
is true only on the final replay pass. Neither the finder selection routines nor
their arithmetic are instrumented. All diagnostic finders and wrappers use the
same explicit ifx precise-O3 options. The unmodified production executable is
also run on the small PM fixture and compared byte for byte with diagnostic v3.

Example build, after loading Intel 2024.2 `compiler-rt tbb compiler` and saving
the module runtime `LD_LIBRARY_PATH` as `BDM_AUDIT_NATIVE_LIBS`:

```bash
micromamba run -n cosemu python3 -B native/build_native.py --repo /absolute/MG-GLAM --root /absolute/follow-up
```

The native replay command is:

```text
replay.exe STEP THREADS PASSES DENSITY_FILE read|write
```

It reads `PMcrd.STEP.DAT` and the matching PM particle pages in its working
directory. The PM header is checked before particle allocation: its actual
step must match, `0 < NROW < 1200`, and `Nparticles == NROW**3`. This wrapper
supports the complete equal-mass PM snapshot used by this audit. It never
evolves the simulation. Density paths have a 4096-character argument buffer;
oversized arguments fail explicitly.

In `write` mode, the wrapper computes DENSIT once at the requested thread count
and stages/publishes that exact density tape without replacing an existing
file. In `read` mode it checks tape dimensions and byte length and reads it
without mutation. The caller verifies snapshot and density SHA256 identities.
Both modes keep one immutable FI copy and six arrays of original PM float bits.
Before every `BDM(0)` call they restore FI and verify every density bit; after
every call they verify all six original particle arrays, counts, FI allocation,
publication completion, and all finder workspace arrays, including both new
telemetry arrays. Diagnostic copies need 24 additional bytes per particle and
4 bytes per mesh cell beyond the normal finder. These extra control copies
are outside production `Memory` admission accounting; the Slurm allocation
must include them. Initial PM FI accounting matches the production entry.

Successful catalogues are preserved as `p1.DAT`, `p2.DAT`, etc. Existing pass
outputs are rejected. Each pass emits:

```text
REPLAY FIXED FI pass=1
REPLAY FINDER pass=1 seconds=      123.456789
REPLAY RESTORED pass=1
```

Only the final pass writes both diagnostic tapes and emits
`REPLAY DIAGNOSTICS candidates=N selected=M`. Full success ends with
`REPLAY COMPLETE`; exit status alone is insufficient. Final-pass finder timing
includes tape output, so use the first pass for whole-finder timing and stage
logs for List/AddBuffer comparisons. The parent driver records individual
measurements and does not infer a general scaling law from this replay.

All tapes are unformatted streams with big-endian numeric fields and no record
markers or alignment padding:

* Density: int64 mesh dimension, int64 original particle count, int32 density
  thread count, followed by float32 FI in Fortran `(x,y,z)` order. Size is
  `20 + 4*NGRID**3` bytes and matches the earlier fixed-field control tapes.
* `repair-members.bin`: the established membership format, int64 selected
  count, int64 candidate count, float32 MassOne; each selected row has int64
  candidate ID, int64 member count, 21 float32 raw properties in the existing
  validation order, and that many int64 original particle IDs. The selected
  count must equal the production writer's `Nhalo`.
* `unbinding.bin`: int64 candidate count, then exactly 28 bytes per implicit
  candidate `1..Nmaxima`: int32 passes, int32 status, int64 sum of active rows,
  int64 retained bound-ID count, float32 post-selection Mvir. A duplicate/host
  rejection can leave retained IDs while setting Mvir to zero. These counters
  cover all attempted candidates, including unpublished ones. A singular
  potential attempt counts as one pass. No pass cap is introduced.

The small one-core preflight suite exercises every variant with two calls on
both populated and empty meshes, verifies exact 64-particle memberships,
checks the final-pass dump count and tape layout, and tests publication/cleanup
against the actual linked modules. A late-invalid second candidate must leave
a prior valid catalogue unchanged. Positive failure controls inject retained
pass/work arrays, malformed or truncated density tapes, a header/step mismatch,
and a 1200-cubed header (rejected before large allocation). Inputs are rehashed
after all runs. All synthetic data and binaries stay under the ignored work
directory; no historical audit or repair receipts are changed.
