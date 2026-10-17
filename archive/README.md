# archive/ — consolidated finished runs and notes (tarballs are gitignored)

Created 2026-10-17 to reduce the file count of the repository root
(AGENTS.md: consolidate completed test scratch and finished-job files).
Every tarball was verified against the originals (`tar -d`) before they were
removed. Restore with `tar -xzf archive/<name>.tar.gz` from the repository
root (the notes tarball extracts into the current directory).

| tarball | contents | why archived |
|---|---|---|
| `conv_da4.tar.gz`, `conv_da8.tar.gz`, `conv_da16.tar.gz` | July 2026 time-step convergence runs (L512, 36 steps per `da`; Init/Setup/PkTable, `Run1/` P(k) outputs, logs; Slurm jobs 11555300-11555302) | pre-fix `sync` velocity runs; their "+1-5 % convergence" was the half-step offset (see `2LPTIC_Gui/README.md`, VALIDATION.md 2026-09-17 erratum) |
| `ic2pm_val_L1024.tar.gz` | July 2026 `ic2pm` validation run inputs/logs (8 files) | superseded by the half-step fix and `test/extras/run_ic2pm_test.sh` |
| `halfstep_ab.tar.gz` | September 2026 half-step A/B test: `da4/da8/da16 × half/sync` runs, serial diagnostics, Slurm logs, frozen `bin/` (602 files, 22 MB) | finished; conclusions in `halfstep_ab/README.md` (kept in place with `compare_halfstep.py`) and VALIDATION.md |
| `review-and-session-notes-2026-04-to-09.tar.gz` | ic2pm half-step review/handoff prompts and findings (`2LPTIC_Gui/REVIEW-*`, `FIX-PROMPTS`, `HANDOFF`, `VERIFY-6b629f1`), BDM duplicated-haloes debug prompts, `op.md` (OMP summary, superseded by `OMP_OPTIMIZATION_NOTES.md`), `thread_num.md`, `conv_jobids.txt` | finished tasks; results live in the tracked docs |

Not archived (still in use or awaiting a decision): `DEGRACE-companion/runs/`
(pilots; the f(R) pilot was running), `DEGRACE-pilot/`, `DEGRACE-highres/`,
`verify/`, the production run directories (`DurMun_*`, `DEGRACE_*`, `fid_*`),
`fid2LPTIC_*` (superseded pre-fix fiducial runs, 15 GB of catalogues),
`2LPTIC_Gui/validate_highres/snap` (213 GB regenerable IC),
`BDM-refine/validation/*/work` (459 GB of PM files).
