# BDM finder audit: physics, then numerics and efficiency

Branch: `audit/bdm-physics-numerics-20260905`.
Source baseline: `a8c7715`, the BDM repair merged into `cz`.

This is an audit of the actual source in that revision, including the standalone
entry point and the simulation's analysis call. Production catalogues and finder
source are not modified by the audit. Confirmed defects are reproduced in small
Fortran cases extracted verbatim from the source, with independent physical
oracles where appropriate. Read-only instrumentation is labelled explicitly.

Stage 1 covers configuration and background conventions; density peaks and
centring; periodic neighbours; spherical-overdensity radii and mass semantics;
unbinding, velocities, spin and shape; host exclusion and duplicate survival;
output validity; and effects on the simulation particle arrays.

Stage 2 covers bounds, undefined values and floating-point exceptions; compiler
and thread sensitivity; memory requirements; linked-list and finder scaling;
and reproducible resource accounting. Failed correctness checks are not treated
as passing simply because an executable exits successfully.

Slurm runs use the shared `cosma8-serial` partition, explicit cores and memory,
and a pilot before further runs. The particle limit is strictly below `1200^3`.
Logs and results are consolidated to limit inode use. A small GR simulation is
an integration diagnostic, not a mass-function convergence study.

The audit is complete with **changes required**: see [AUDIT.md](AUDIT.md) for
16 prioritized finding groups, evidence, scope limits and the recommended
repair order. The two new baseline catalogues contain no identical surviving
particle sets; controlled cases expose other physics and numerical failures.
Experimental configuration and performance changes remain separate from the
production source.

## Results and provenance

| Artifact | Contents |
|---|---|
| [summary.json](summary.json) | Finding index, membership checks, catalogue comparisons, timings and total resource use |
| [unit-results.json](unit-results.json) | All 46 controlled experiments, including intentional failures, source hashes and compiler commands |
| `simulation-n64.json`, `simulation-n128.json` | Simulation inputs, snapshot/executable hashes, process measurements and completion state |
| `replay-results.json`, `bounds-results.json` | Instrumented, configuration-order, thread, bounds and prefilter experiments |
| `diagnostic-build.json`, `build.log.gz` | Native build commands, diagnostics and variant hashes |
| `*.npz`, `*.log.gz` | Consolidated numerical catalogues and full process logs |
| [jobs.json](jobs.json) | Five allocations, submission provenance, dependencies and final Slurm accounting |
| [validation.json](validation.json) | Final integrity checks, artifact hashes and verified cleanup inventory |

The local, Git-ignored `work-artifacts.tar.gz` preserves scratch builds,
executables, simulation inputs, initial conditions, final snapshots and batch
logs. Every archived regular file was checked by SHA-256 and every symlink by
target before the loose scratch tree was removed. Original external IC input
files are also included under `reference-inputs/`. Build symlinks retain their
original absolute targets; rebuild them if moving to a different checkout.
The archive is local evidence and is not included in a clone of this branch.

## Reproduction

Use a separate checkout to retain the original result files: the drivers write
new results at the same paths. Commands below start at the repository root.
Python always runs in `cosemu`; the controlled cases use GNU Fortran 14.1.0,
and the native runs use ifx 2024.2.0.

```bash
audit_dir=BDM-refine/analysis/full-audit-20260906
micromamba run -n cosemu python3 -B "$audit_dir/run_audit.py"
bash "$audit_dir/build_full.sh"
export LINES=40 COLUMNS=120
module purge
module load intel_comp/2024.2.0
module load compiler-rt tbb compiler
export BDM_AUDIT_NATIVE_LIBS="$LD_LIBRARY_PATH"
micromamba run -n cosemu python3 -B "$audit_dir/build_diagnostics.py"
micromamba run -n cosemu python3 -B "$audit_dir/build_diagnostics.py" --bounds-only
```

For simulation reproduction, supply the recorded `test/Init.dat`,
`test/PkTable.dat` and `test/TableSeeds.dat` inputs; copies are in the local
archive. Submit the small pilot first, inspect its completion and accounting,
then submit the N128 case:

```bash
sbatch "$audit_dir/pilot.sbatch"
# After the pilot completes and its output and resource use are checked:
sbatch --cpus-per-task=4 --time=00:02:00 "$audit_dir/pilot.sbatch" --nrow 128 --threads 4
# After both final snapshots are verified:
sbatch "$audit_dir/replays.sbatch"
sbatch --cpus-per-task=1 --time=00:02:00 "$audit_dir/replays.sbatch" --bounds-only
```

Record new job IDs and script hashes separately from the historical `jobs.json`.
The `--resume` simulation option applies only to an incomplete run with verified
initial conditions. The normal driver refuses to replace an existing run.
Allocation sizes above reproduce these small experiments; they are not sizing
advice for a larger campaign.

For replay using the retained snapshots, first restore `work/` from the archive:

```bash
tar -xzf "$audit_dir/work-artifacts.tar.gz" -C "$audit_dir" work
```

To regenerate the summary from recorded evidence without new simulations or
Slurm queries:

```bash
micromamba run -n cosemu python3 -B "$audit_dir/summarize_audit.py"
```

Add `--refresh-accounting` to query `sacct` again for the recorded job IDs.
The final validation manifest describes the committed evidence; intentional
reruns produce different timings, logs and hashes.
