# N1024 fixed-density thread control

This independent experiment isolates finder threading from the upstream
floating-point accumulation order in `DENSIT`. It uses the final shape-repaired
production objects recorded in `../native-build.json`, with the diagnostic
Fortran from `../numerical-threading/` copied byte-for-byte. Production code is
not recompiled or edited. The earlier N512 build, failures, results and density
tapes remain untouched.

The job depends on successful main simulation **11948372**. It reads that job's
completed `main-simulation.json`, verifies the final z=0 snapshot, and derives
the actual snapshot step from the receipt and native big-endian PM header.
The N1024 particle count is 1,073,741,824, below the authorized 1200³ limit.
All large allocations and full snapshot/tape hashing run inside Slurm.

The unchanged native diagnostic runs this sequence:

1. Compute `FI` once with 64 density threads and save the exact field.
2. Restore that immutable field for finder calls at 32, 64, 32 and 64 threads.
3. Compute `FI` once with 32 density threads and measure exact cell differences.
4. Restore the second immutable field for finder calls at 32 and 64 threads.

Every call verifies exact restored density bits and original particle position
and velocity bits. All six calls must finish and preserve their complete
published catalogues. Catalogue bytes must agree within each fixed-field group.
Comparisons between the two fields and against normal simulation/replay
catalogues are reported without a tolerance or a requirement of equality.
Such differences are density sensitivity, and do not by themselves establish a
finder race. Published row matching across distinct fields is descriptive;
candidate identity is not assumed. This diagnostic does not save raw bound
particle membership sets; those are checked in the separate normal validation.

`results.json` retains the pilot-compatible `completed`,
`fixed_density_controls_passed`, `fixed_density_comparisons`,
`between_density_comparison`, `density_difference`, `native_finder_timings`,
`peaks`, native CPU/time/RSS and `outputs_sha256` fields. It additionally records
`fixed_density_groups`, `density_field_difference`, `restoration`,
`normal_comparisons` and `density_tapes`. A failed byte comparison leaves a
failed receipt and preserved output, not a successful validation. Existing
results/run directories cause a refusal to restart; inspection and a separate
attempt directory are required. Normal validation job 11948467 is optional at
collection time, not a dependency. Completed native replay receipts are marked
as process-level evidence if their membership check is still in progress.

The request is **cosma8-serial, dp004, 64 CPUs, 288 GiB, 45 minutes**. The measured
N512 control used 31,325,680 KiB batch MaxRSS (29.8745 GiB). Eightfold volume and
20% headroom gives 286.80 GiB, rounded to 288 GiB. At N1024, diagnostic storage
adds one immutable 32 GiB field and 24 GiB saved particle bits. Active `FI` is
separate and BDM releases/reallocates it. The two retained 32 GiB density tapes
are disk outputs, not two simultaneous immutable in-memory fields. Expected
runtime is 10–20 minutes (10.7–21.3 allocated core-hours); the limit is 48
core-hours. Shared scheduling leaves the remaining 64 cores available for the
independent normal validation.

Reproduction uses `micromamba run -n cosemu python3 -B` throughout. Load
`intel_comp/2024.2.0`, `compiler-rt`, `tbb`, `compiler`, and capture
`BDM_AUDIT_NATIVE_LIBS="$LD_LIBRARY_PATH"` before entering the Python environment.
`prepare_control.py --original-repo /path/to/MG-GLAM` copies the proven builder,
Fortran, modules, objects and production sources into new ignored `work/`
scratch and runs the tiny one-core publication preflight. `test_control.py`
checks synthetic headers, source/snapshot mismatches and false completion.
`submit_control.py` verifies the sealed build, sources, runner and test receipt,
then records the exact Slurm command, hashes, dependency, resources and job ID
in `submission.json`. The plan hash is passed explicitly to the frozen runner.

Retain `work/run` logs, catalogues, peak tapes, both exact density tapes and
all receipts until archive/relocation has been verified. This directory's
evidence is independent of the earlier pilot evidence.

## Completed result: job 11948491

The job completed on 2026-09-06 with both fixed-field groups passing strict
catalogue byte equality. Every catalogue contains **147,042 rows**. All six
restored fields matched exactly, all six calls restored the original particle
arrays bit-for-bit, and the snapshot/config hashes still matched after the
experiment. Actual final snapshot step was **158**, selected from the completed
main simulation receipt rather than assumed from the pilot.

The two density fields differ in **46,743,773 cells**, with maximum absolute
`FI` difference **0.06640625**. Their **794,570 peak indices and positions are
identical**; 204,646 peak density values and 67,074 empirical seed radii differ.
Changing only the saved field changes seven published rows: 23225, 26108,
62218, 65341, 66851, 94799 and 99155. The largest coordinate-component change is
0.0007 Mpc/h. Six rows keep the same bound mass, count and bulk velocity. Row
99155 changes bound count from 2564 to 2566 and printed bound mass from
2.7463e13 to 2.7485e13 Msun/h (about 0.080%). These differences remain recorded
as numerical density sensitivity; no acceptance thresholds were relaxed.

The controlled `d64` catalogue exactly matches the main inline 64-thread
catalogue. The controlled `d32` catalogue exactly matches the normal 32-thread
baseline. Separate normal 64-thread baseline and membership runs differ from
`d64`/inline in one additional row, 7383, in total mass, offset and virial ratio
columns. This illustrates variation between separate normal density
realizations even at the same configured thread count. The state replay was
still pending when optional comparisons were collected, as explicitly recorded
in the immutable result. Changing finder threads with a frozen field changes
no catalogue bytes in this experiment; this is a bounded control rather than a
claim that every possible finder input is free of races.

The six finder times were 174.35, 110.76, 179.66, 107.18, 178.20 and 104.17 s.
Native elapsed time was 1007.35 s and process MaxRSS was 159,517,752 KiB
(152.13 GiB). The complete allocation lasted **1151 s**, billed **20.4622
core-hours**, and consumed **29,429 CPU seconds** (39.95% CPU utilization).
Measured batch MaxRSS was **198,231,384 KiB (189.05 GiB)**, leaving 34.36% of the
288 GiB memory request unused. Requested and allocated TRES both record
64 CPUs, 288 GiB, one shared node and billing=64. The run stayed inside its
10–20 minute estimate and 48 core-hour time limit.

`results.json` contains full row differences and the hashes of both exact
32 GiB density tapes. `accounting.json` preserves final `sacct` output and cost
calculation. `successful-probe.log.gz` is the complete native log compressed
without changing its contents. All raw output, inputs and frozen build files
are preserved through the relocation map below.

## Consolidation

`relocation.json` maps the original worktree paths for both
`numerical-threading/work` and `main-thread-control/work` to the matching
directories under the main MG-GLAM repository. Each directory was moved by
same-filesystem `os.rename`, with every entry's device, inode, size and mtime
verified unchanged and every input symlink still resolving. Existing validated
SHA manifests are linked in the relocation receipt. Historical execution
receipts retain their original paths; use this map to locate their artifacts.
