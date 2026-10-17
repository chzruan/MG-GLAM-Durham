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
(pilots), `DEGRACE-pilot/`, `DEGRACE-highres/`, `verify/`, the production run
directories (`DurMun_*`, `DEGRACE_*`, `fid_*`), `fid2LPTIC_*` (superseded
pre-fix fiducial runs, 15 GB of catalogues), `BDM-refine/validation/*/work`
(459 GB of PM files).

## Removed without an archive (user decision 2026-10-17)

| removed | what it was | how to recover |
|---|---|---|
| `2LPTIC_Gui/validate_highres/snap/IC_ours_Np2048_L1024.{0..63}` (64 Gadget-1 files, 275 GB, written 2026-07-28 by job 11658977) | the regenerated 2048^3 L1024 seed-2026 validation IC (HEFT normalisation). Its validation record is `2LPTIC_Gui/VALIDATION.md` ("High-res validation") with `pk_validate_highres.npz`; `snap/pofks_ours_2048.txt` (the generator's own P(k)) is kept | deterministic: `cd 2LPTIC_Gui/validate_highres && sbatch job.sh` (64 ranks, cosma8-shm, 3000G, 8 min; `pofk_bli_z49.txt` is beside `input.lua`) |
| `DEGRACE_LCDM_L1800Np3200Ng6400/Run1/PMcrd.DAT` + `PMcrs{0..61}.DAT` (62 particle files, 786 GB = 732 GiB, written 2026-02-12) | the mid-run restart checkpoint of the finished L1800 LCDM Run1: header `N=3200x6400L=1800zi100da=4.00E-04f0.040`, ISTEP = 100, AEXPN = 0.41824 (z = 1.39), Nseed = 13979821. The run went on to step 157 (z = 0); all 11 `PowerDM.*` spectra and `CATALOGS/CatshortV.{0092..0157}` (19 GB) are kept. Run2..Run12 never had PMcrs files | not recoverable as such (no full particle snapshot of Run1 remains); a rerun of box 1 (`PMP2start.exe <<< 1; PMP2main.exe <<< 157`, 72 h cosma8-shm slots with checkpoint restarts) reproduces it |
