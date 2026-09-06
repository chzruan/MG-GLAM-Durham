# N3/N4 particle-workspace optimizations

Base: `54f53aae2aac6700d1ddb8490bebb0719f71165d`. This change addresses N3 and N4
in the immutable [merged review](../../../analysis/review-20260906-merged/REVIEW.md).
Production edits are confined to `List` and `AddBuffer`. SO normalization,
unbinding, halo populations, periodic coordinate reconstruction, the PM restore
routines, configuration, output serialization and build flags are unchanged.

## Linked-list cache

The parallel list path computes each exact periodic z-cell index once. It stores
a signed int16 index when both declared z bounds fit [-32768,32767], otherwise
int32. Every slab still scans rows in ascending order; each cell has one writer,
so its linked traversal remains the same descending row-ID sequence as the base.
The scan remains O(T*Np), but its repeated reads are contiguous two- or four-byte
cell indices. Coordinate reconstruction and ceiling arithmetic are O(Np).

Plane counts, slab boundaries, empty-slab endpoints and label-zeroing loop
indices use int64. This avoids an overflowing int32 bound difference or a final
loop increment at an int32 extremum. Cache conversion is performed only after
checking the signed storage range. Scratch is charged and released through the
existing four-byte-word `Memory` API, rounding the int16 allocation up by at most
two bytes, and checked against `MaxMemory` before allocation.

For a validated `AddBuffer` workspace, original coordinates are in [0,Box) and
image shifts are -1,0,1. If twice Box is representable in float32 and
`2*Box/Cell < huge(int32)`, the original fast ceiling/subtraction is provably
representable. Retaining it avoids the measured ifx cost of unnecessary guarded
coordinate conversions. Unbuffered callers or other scale ranges use an explicit
clamp before conversion: an interior value satisfies
`lower+1 < coordinate/Cell <= upper`, which makes the int32 ceiling and subtraction
safe. This gives the same cell as `clamp(ceiling(coordinate/Cell)-1)` at exact cell
faces and signed bound extremes. Invalid cell sizes, reversed bounds, negative
particle counts and nonfinite coordinates on the guarded path fail explicitly.

## Periodic image counts and prefix

`PrepareParticleSearch` retains its strict half-box limit, so each coordinate
has at most two allowed shifts. Each original therefore has 0..7 ghost rows.
`AddBuffer` stores that count in one byte and records an int64 total per 16,384
originals. Counting and filling run over independent blocks; only the block
prefix remains serial. Within a block the original rows and each row's k/j/i
image order are unchanged. Originals remain first and every replica preserves
its original int64 ID and all six float32 fields.

Each row still checks its actual fill against its counted number; each block
also checks its endpoint. Coordinate-domain, active-buffer and original-count
checks are retained. Additional checks prevent invalid byte counts, negative or
unrepresentable original counts, overflowing prefixes and overflowing allocation
accounting. Empty input and partial final blocks are supported. Scratch is freed
before the old analysis arrays are freed and the populated arrays moved into
place; no PM originals are altered.

For N original particles the old scratch was `8*N` bytes. The new scratch is
`N + 8*(ceil(N/16384)+1)` bytes, up to three extra accounted bytes from the word
interface. At 1024^3 this is 8,589,934,592 versus 1,074,266,120 bytes, a reduction
of 7,515,668,472 bytes (6.9995 GiB). The serial prefix has 65,536 entries instead
of 1,073,741,824. The List cache is a separate, later transient allocation;
these scratch figures are not an estimate of the complete finder's peak RSS.

## Verification and reproduction

All new evidence is in this directory. The historical review and old receipts
were not edited. `run_checks.py` extracts actual production routines, compiles
the base and updated variants with matching flags, and checks their complete
particle/identity/list-array byte streams. Its independent image oracle enumerates
all -1,0,1 shifts and filters by the physical buffer boundaries. Its list oracle
uses the wide-integer ceiling definition and builds each row's reference link.

The focused matrix has 342 experiments per compiler: checked and optimized
builds, 1/2/4/8 threads, counts 0/1/16383/16384/16385/32769, exact faces and their
adjacent float32 values, fractional box size, the half-box cap, periodic replicas,
full int16 range and both fallback boundaries, more threads than z slabs,
int32 extrema, and positive failure controls. Int32-extreme cases use only the
updated routine and the independent oracle: the base's default-integer loop or
conversion was undefined there. Failure checks reject a missing expected error
rather than accepting the test driver's own failure as a production rejection.

The unchanged historical particle harness additionally runs its 104 experiments
per compiler, including exact PM restoration over two cycles with signed zeros
and subnormals, centring, minimum-image membership and asymmetric list bounds.

```sh
micromamba run -n cosemu python3 -B \
  BDM-refine/repairs/20260906-review-followup/particles/run_checks.py \
  --compiler gfortran --output /tmp/particles-focused-gnu.json
micromamba run -n cosemu python3 -B \
  BDM-refine/repairs/20260906/tests/particles_tests.py \
  --output /tmp/particles-legacy-gnu.json
```

For ifx use the native runtime before entering `cosemu`:

```sh
export LINES=40 COLUMNS=120
module purge
module load intel_comp/2024.2.0
module load compiler-rt tbb compiler
export BDM_AUDIT_NATIVE_LIBS="$LD_LIBRARY_PATH"
micromamba run -n cosemu python3 -B \
  BDM-refine/repairs/20260906-review-followup/particles/run_checks.py \
  --compiler ifx --output /tmp/particles-focused-ifx.json
```

GNU checked builds include bounds, integer-overflow and FP traps; ifx uses
`-check bounds -fpe0`, with precise arithmetic in both focused modes. An initial
`ifx -check all` harness attempt failed in MemorySanitizer's
`_GLOBAL__sub_I_fast_mem_ops.c` before even the baseline entered the program.
The established bounds-check configuration avoids that compiler-runtime startup
failure. The historical ifx optimized harness separately retains its previous
fast=1/AVX2/FTZ settings and reports them in its new receipt.

## Kernel measurements and integration

The bounded pilot uses 262,144 original rows and 850,823 buffered rows, with
Box=32 and NGRID=128. Each sample repeats a kernel six times; three alternating
base/updated samples are taken at each of 1/2/4/8 threads. Every sample also
passes the independent oracles and complete-array byte comparison. Timed regions
exclude oracle construction, file hashing and the untimed setup. GNU-time
resource receipts cover the complete script including compilation and checking.

These are small synthetic kernels on a lightly loaded login host, not measured
N1024 finder speedups. The coarse mesh gives a larger replica fraction than the
production mesh. Shared-host timing, allocation and cache noise remain; buffer
runtime gains in particular are small and variable. The deterministic memory
reduction is the main established benefit of N4. Final timings and resource
measurements are summarized in `summary.json`, with full samples in the two
`benchmark-*.json` receipts.

To isolate these output-identical optimizations from the separately versioned
SO normalization and unbinding changes, transplant only `List` and `AddBuffer`
from this commit onto the base, retaining the same compiler flags and immutable
particle/density inputs. The receipts preserve both routine hashes. Require
catalogue/membership/PM equality for this variant first; assess changed halo
values from the normalization variant separately. Parent integration owns the
saved N1024 replays and full finder stage timings.

Final median base/updated ratios (larger means faster):

| Kernel/compiler | 1 thread | 2 threads | 4 threads | 8 threads |
| --- | ---: | ---: | ---: | ---: |
| list/gnu | 0.941 | 1.002 | 1.433 | 1.685 |
| list/ifx | 0.993 | 0.986 | 1.193 | 1.664 |
| buffer/gnu | 1.114 | 1.082 | 1.081 | 1.088 |
| buffer/ifx | 0.969 | 1.069 | 1.050 | 1.217 |

All 256 successful focused-case byte streams also agree across GNU and ifx.
