# N512 fixed-density thread control

The two catalogue differences from pilot job 11948326 are caused by differences
in the upstream `DENSIT` field. Holding that field fixed makes the unchanged BDM
finder's output byte-identical across 16 and 32 threads. Swapping only the field
reproduces the two original catalogues byte for byte, including exactly the two
changed rows. This identifies the cause of this pilot discrepancy; it does not
claim universal reproducibility for every input or thread count.

The original pilot validation remains failed under its byte-equality contract.
Its receipts and thresholds have not been altered. This separate diagnostic adds
controlled evidence explaining that failure. No production source was edited.

## Controlled experiment and results

Job **11948371**, on `cosma8-serial`, read the SHA-verified L512 N512 Ng1024
snapshot at step 157, z=0. The binary links a new diagnostic module/entry to the
unchanged final-b6 native objects restored from the repair archive. The finder
source SHA is `868bb0d6c89e5daf1ca0f631419d86577dcfe31dde7969375398cb0408140fc0`.
The archived build records its pre-merge build commit as `32087c8`; those source
bytes are the production code merged into `cz` at `b6e9666`.

1. Compute `DENSIT` once at 32 threads. Save the exact float32 field to a stream
   tape and an independent immutable array.
2. Restore that same saved array before each `BDM(0)` call, at 16, 32, 16 and 32
   threads. Verify every density bit before each call and all six original
   particle arrays after each call. Preserve the complete published catalogue.
3. Compute a second `DENSIT` field once at 16 threads, compare it with the first,
   then restore this same second field before calls at 16 and 32 threads.

| Input field | Finder threads | Full catalogue SHA-256 |
| --- | --- | --- |
| One saved 32-thread field | 16, 32, 16, 32 | `4f042274a341651771545f51bd5eb2b2c895194e6fa7f7ddf180b6f07c154313` |
| One saved 16-thread field | 16, 32 | `bd6c779cb9e45e043934a4daa73dd8c935522578aa7855a5831422ffef2bed07` |

Every catalogue contains **71,587 rows**. The first hash exactly matches the
original inline-32/member-32 output. The second exactly matches the original
normal baseline-16 output. All six density-restoration and particle-restoration
checks pass. The experiment therefore separates density accumulation order from
finder execution order without using independently recomputed fields as a
purported fixed-input control.

The two independently accumulated fields differ in **4,394,702 / 1,073,741,824
cells**, with maximum absolute difference **0.01611328125**. All **121,835 peak
positions and candidate indices agree**. Peak field values differ for 23,446
candidates; the corresponding initial centering-radius formula differs for
7,375 candidates.

| Published row | Candidate | Peak FI at 32 threads | Peak FI at 16 threads | Initial radius at 32 threads (Mpc/h) | Initial radius at 16 threads (Mpc/h) |
| --- | --- | --- | --- | --- | --- |
| 42511 | 72204 | 567.1731567382812 | 567.1732788085938 | 1.3806530237197876 | 1.3806531429290771 |
| 55520 | 94460 | 3097.221923828125 | 3097.220458984375 | 1.7461861371994019 | 1.7461860179901123 |

Each radius changes by one float32 step, about 0.119 pc/h. The centres then
change by up to 1.9 kpc/h, and a small set of aperture/property values changes.
Every published bound mass, particle count and bulk velocity agrees across the
two catalogues. Row-by-row values and all changed columns are in
[results.json](results.json). This uninstrumented probe did not dump bound member
IDs, so agreement of membership sets across the two different fields is not
claimed from equality of their counts.

The causal path is visible in the unchanged production code:

- `PMP2mod_density.f90:237–281` adds CIC contributions through OpenMP atomic
  updates into float32 `FI`. Atomicity prevents lost updates; it does not impose
  one floating-point addition order across thread counts.
- `PMP2linker.f90:1660` temporarily stores the selected peak's field value in
  `Xoff`. The peak locations and ordering are identical in this experiment.
- `PMP2linker.f90:1542` uses that value to set the centering aperture, and
  `:1554` applies a strict particle-distance cut during recentering. These are
  the density-dependent inputs that remain after the mesh is deallocated.

The exact catalogue reproduction plus fixed-field controls establish the source
of the observed difference. Individual centering iterations and particles on
the changing aperture boundary were not separately traced.

## Resource accounting

The successful job allocated **32 CPUs and 32 GiB for 136 seconds**, with
`TotalCPU=33:47.275`, about 46.6% allocation CPU utilisation. Actual billing was
1.209 core-hours; the five-minute limit was 2.667 core-hours. The native process
took 115.554 seconds and reported MaxRSS **20,359,012 KiB (19.416 GiB)**. Slurm
batch MaxRSS was **31,325,680 KiB (29.874 GiB)**. Raw `sacct` and `scontrol`,
request/allocation TRES, estimates and calculations are in
[accounting.json](accounting.json) and [submission.json](submission.json).

The probe's additional resident arrays are one immutable FI copy,
`4 * Ngrid**3` bytes, and the saved particle bits, `24 * Nrow**3` bytes. These add
4+3=7 GiB at N512/Ng1024, and 32+24=56 GiB at N1024/Ng2048. The normal active FI
is separate; BDM deallocates it during halo processing. The two persistent
density tapes are sequential disk outputs, not two simultaneous immutable
arrays. Each contains a 20-byte header followed by `4 * Ngrid**3` bytes.

Resident-array counts alone underestimate the observed batch accounting.
Scaling measured batch MaxRSS by eight gives about 239.0 GiB for the larger
case; adding 20% gives 286.8 GiB. Root selected a separate **288 GiB** shared
allocation for the later N1024 control. This is a sizing estimate, not a measured
N1024 result; its file-cache behaviour and finder allocations must be checked.

The first attempt, **11948363**, stopped after its first successful finder call
because the diagnostic tried to move `CatalogueFinalPath`, which publication
deliberately clears. The corrected probe uses the retained `outputName` and
passes a tiny allocation-free publication preflight. The failed job's source,
build, submission, output and accounting remain in
[attempt-11948363.json](attempt-11948363.json). It supplies no cross-thread result.

## Reuse and retained evidence

`thread_probe.f90` contains `RunBdmThreadControls(high_threads, low_threads)`;
`thread_probe_entry.f90` reads `snapshot_step high_threads low_threads` from
standard input, loads the native snapshot and allocates FI. Use distinct thread
counts, with the higher first. This diagnostic is configured for `iVirial=1`.
It can be linked against another frozen production finder without changing that
finder's source or objects. `run_probe.py` deliberately describes this one N512
experiment; use a separate copied/parameterized runner for N1024 and preserve
these build/results receipts.

In a fresh copied diagnostic directory, load the same Intel environment and
link against the final native build and its matching receipt:

```bash
export LINES=40 COLUMNS=120
module purge
module load intel_comp/2024.2.0
module load compiler-rt tbb compiler
export BDM_AUDIT_NATIVE_LIBS="$LD_LIBRARY_PATH"
micromamba run -n cosemu python3 -B build_probe.py \
  --native-dir /absolute/path/to/final-native-build \
  --build-receipt /absolute/path/to/native-build.json
```

The builder verifies the frozen source receipt, records every linked object and
module hash, compiles only the diagnostic code, and runs the tiny publication
preflight on one core. Large native runs require Slurm. The present Slurm
wrapper and exact submission command are retained, with their hashes.

The complete native log is also in [successful-probe.log.gz](successful-probe.log.gz).
`work/run/` retains both exact successful density tapes, both peak tapes and
logs, all six catalogues, timings and snapshot symlinks. `work/run-11948363/`
retains the failed attempt's artifacts. Both corresponding native build
directories are retained. Root must verify their archive or relocation before
removing this worktree. No heavy archive or cleanup has been performed here.
