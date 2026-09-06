# Build/configuration review follow-up: N1, N5, N6

This repair starts from `54f53aa` and addresses the corresponding findings in
`analysis/review-20260906-merged/REVIEW.md`. It changes the finder build flags,
memory-limit configuration and unsupported entry-point handling. It does not
change the SO/unbinding algorithms, particle list/buffer algorithms, peak
selection, catalogue writer or the eight-line catalogue header.

## N1: append precise flags after overrides

The finder rule appends `-fp-model precise` after its other Intel compiler flags,
including an overridden `FFLAGS` or `BDM_FFLAGS`. The recursive
`PMP2main-bitmatch` target therefore receives the same guarantee. Other object
and link commands keep their existing flags. GNU `FC=gfortran` overrides remain
supported, with trailing `-fno-fast-math -ffp-contract=off`. The explicit finder
rule again includes the `-w` present in the suffix rule. `BDM_PRECISE_FLAGS` is
the explicit compiler-specific safety setting; ordinary flag overrides do not
silently replace it.

Actual **ifx 2024.2.0** compilation and execution establishes the ordering. For
runtime input `a=1+2^-23`, `b=1-2^-23`, `c=-1`, the real32 expression `a*b+c`
produces zero with precise rounding and `-2^-46` (bits `A8800000`) with fast FMA
contraction on this AVX2 machine. `fast=1` followed by `precise` gives zero;
reversing the flags gives `-2^-46`. The compiler's warning 10121 explicitly
identifies which earlier option was overridden. `-Ofast` and `-ffast-math`
followed by precise also produce zero.

Eight **actual make/compile/link/runtime** fixtures exercise the old rule's
failure, the repaired default, no-model and conflicting-model `FFLAGS`, an
overridden `BDM_FFLAGS`, the canonical recursive bitmatch build, an AVX2/FMA
bitmatch override, and GNU flags. The fixtures use the actual makefile with tiny
stand-in module sources so that the compiler and recursive-target behavior can
be tested without a simulation. Every repaired fixture produces zero; the old
`FFLAGS` override reproducer produces `-2^-46`. This is compiler/runtime evidence,
not solely a printed make command.

## N5: configurable MaxMemory

`BDM.config` accepts `MaxMemory`, case-insensitively, as a finite positive real
number in **GiB**. The per-call default remains **500 GiB**. Generated default
configuration and stdout report the value and units; omitted keys revert to
the default on each `ReadParameters` call. Zero, negative, nonfinite,
unrepresentable, malformed and multi-value entries fail validation before
opening the analysis log or staged catalogue. Existing output files are
preserved on invalid configuration.

For example, `MaxMemory = 512 ! GiB` changes the existing linked-list admission
gate without rebuilding. This remains an estimate of tracked finder storage
plus the prospective int64 linked lists, **not a process RSS limit or a Slurm
memory request**. The unchanged gate rejects an insufficient limit; it does not
coarsen `Cell` or change physical search radii. The admission tests use a large
integer particle count but allocate no particle/list/mesh arrays.

## N6: enforce unsupported entry points

The production Fortran call graph before this repair contains only
`BDM -> RescaleCoords(1)` and the internal dormant
`GetProfiles -> HaloProfile` edge. It has no callers of `GetProfiles`,
`WriteProfiles` or `RemoveDuplicatesSimple`; the old simple-duplicate call in
`BDM` is commented out. The recorded call graph checks all top-level production
Fortran sources, not the test or archived copies.

The four dormant routines are explicit error stubs with their existing names
and argument types. This removes `HaloProfile`'s automatic profile array, whose
extent could otherwise be evaluated before an executable guard. They fail with
named unsupported-entry errors before indexing unavailable workspace or
creating files. `RescaleCoords` rejects every flag other than 1 and directs
restoration callers to `RemoveBuffer`. The active forward path is unchanged;
an exact particle-bit roundtrip through the real forward/restore routines is
included in the tests.

## Validation and reproduction

`results.json` records **104 configuration/entrypoint experiments**: 26 each
under GNU 11.1 checked and optimized builds and ifx 2024.2.0 checked and optimized
builds. Cases include default/annotated/scientific configuration, reset between
calls, invalid input with prior-output preservation, memory admission above and
below the same fixed geometry, all unsupported calls, and active rescaling with
exact restoration. Six direct compiler precedence controls and eight real make
builds also pass. Source, generated-fixture, driver and binary hashes and exact
commands/stdout/stderr are retained.

GNU checked mode uses `-fcheck=all`; ifx checked mode uses bounds and pointer
checks. An initial ifx `-check all` attempt terminated in MemorySanitizer in the
Intel `fast_mem_ops.c` startup constructor before test entry. No memory-sanitizer
coverage is claimed. No Slurm allocation or simulation was required for these
small one-core checks; the parent audit runs the integrated saved N1024 replay.

Load `intel_comp/2024.2.0`, `compiler-rt`, `tbb`, `compiler`, capture
`BDM_AUDIT_NATIVE_LIBS="$LD_LIBRARY_PATH"`, then run:

```bash
micromamba run -n cosemu python3 -B \
  BDM-refine/repairs/20260906-review-followup/build-config/regression.py \
  --output /path/to/new-build-config-results.json
```

Use a new output path to preserve the committed receipt. Scratch modules,
objects, executables and run files are removed automatically. Existing audit,
review and regression receipts are unchanged.

## Integration follow-up

Cross-review found that inspecting only the last word of `FC` misclassified
`FC="gfortran -m64"`. The family check now recognizes a `gfortran` command
anywhere in the command words, including paths and ordinary command wrappers.
Opaque wrapper names require an explicit `BDM_PRECISE_FLAGS` appropriate to
the underlying compiler. `integration-results.json` adds a real GNU build
with an option in `FC`: all 104 configuration/entry checks, six Intel
precedence controls, and nine actual make/runtime builds passed. The original
`results.json` remains evidence for the earlier recorded source.
