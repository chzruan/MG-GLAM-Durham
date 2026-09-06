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
