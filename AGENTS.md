# Repository instructions for agents

## Default initial conditions

- **Use 2LPTIC by default for new simulations**, including production,
  validation and convergence runs, unless the task explicitly requires another
  IC method.
- The project workflow uses the corrected FML LPT generator in
  [`2LPTIC_Gui/`](2LPTIC_Gui/README.md), configured with `lpt_order = 2`.
  Read its [validation record](2LPTIC_Gui/VALIDATION.md) when preparing a run.
- Generate the Gadget-format ICs with 2LPTIC, convert them to GLAM PM files
  with [`ic2pm`](ic2pm.f90), then evolve with GLAM. Match the IC epoch,
  cosmology, spectrum normalization and velocity conventions to the run.
- Use `ic2pm` built from commit 6b629f1 or later (on `cz` since
  2026-09-17). Its default velocity epoch is `half` (velocities
  shifted to a_init − ASTEP/2, GLAM's leapfrog convention); never shift
  velocities by editing the Gadget IC files. Use `sync` only to reproduce
  runs converted before ee44100, and pass it explicitly
  (`ic2pm.exe <IC> <S_vel> sync`): the legacy submit scripts omit the
  argument and would now get `half`. Pre-fix runs have late-time P(k) high
  by ≈0.6 da/a_init at linear k (more at 0.3<k<1); see the 2026-09-17
  erratum in `2LPTIC_Gui/VALIDATION.md`. Record the epoch (the
  `velocity epoch =` line of the ic2pm log) in each campaign's provenance.
- GLAM's native `PMP2start` generators use first-order Zel'dovich ICs. Use
  them for explicitly requested legacy reproduction or IC-method controls;
  do not silently substitute them for the default 2LPTIC workflow.
- The native generators' white noise was drawn at RANLUX luxury level 0
  until 2026-10-19 (the luxury table was never set), which modulates the
  variance of the field along x at harmonics 6-12 of the box; see
  [`docs/ic-white-noise-20261019.md`](docs/ic-white-noise-20261019.md). The
  default is now level 3. To reproduce the initial conditions of a run made
  before that date, set `GLAM_RANLUX_LUXURY=0` and record it in the
  campaign's provenance; `Results.log` of the run records the level and
  the seed. Anything built against `PMP2mod_random.f90` to reproduce old
  draws (the halobias replica `abcs_haloes/degrace_ic/src/lowk_ts.f90`,
  which hard-codes `lux = 2`) must ask for level 0 explicitly when it is
  rebuilt; see the note's "Downstream users of the module".
- Record the actual generator/version, LPT order, starting redshift, seed
  and input configuration in each campaign's provenance. Resolution
  comparisons also require verified shared modes and consistent normalization.
- Preserve the recorded IC method of existing runs when resuming them.
  The already launched
  [20260908 BDM convergence campaign](BDM-refine/validation/20260908-convergence/README.md)
  uses a frozen native first-order variant and is an exception, not a template
  for the default IC choice in new campaigns.

## After a simulation finishes: verify, pack, delete the originals

- A finished box directory `Run<box>/` holds about 400-430 files (17-19 per
  full output, 4 per lens-plane-only output, about 40 % of them small `.meta`
  sidecars). Post-process every finished box that is kept (production, pilot
  or test): pack it into one HDF5 file per output
  step with `glamio.pack` (format in [`python/README.md`](python/README.md),
  "Packed runs"): 42 outputs give `Output.<step>.<box>.h5` x42,
  `RunInfo.<box>.h5` and `packed.md5`, 44 files. The bytes are stored
  unchanged, and glamio's readers, `product_path` and `products_at` read
  packed runs transparently. Then delete the native originals, to cut the
  file count, but only after the byte check below has passed.
- Order: wait for the simulation job to finish, run
  `DEGRACE-companion/check_run.py` (it reads native and packed runs alike, but
  run it before packing so a failed check is caught while the originals
  exist) and any script that still opens native files by name, then pack.
  Never pack a directory a job is still writing, and never run two pack/prune
  jobs on the same box.
- Use [`DEGRACE-companion/pack_run.sh`](DEGRACE-companion/pack_run.sh)
  (cosma8-serial, one core). Submit it from the run directory as
  `sbatch -J pack_<name> <repo>/DEGRACE-companion/pack_run.sh "$PWD/Run<box>"`.
  It packs in place and proves byte identity: it rebuilds every original from
  the packed files in scratch, requires an identical file list and `cmp`
  equality of every file, then removes the rebuilt copy (also on failure).
  It then deletes the originals in the same job, via `glamio.pack prune`,
  which re-verifies every md5 first, refuses while any original is in no
  packed file, and deletes nothing if anything differs. This deletion is the
  default end of the post-processing and needs no extra confirmation; submit
  with `--export=ALL,PRUNE=no` only when the originals must be kept (for
  example while a native-only script still has to run). Never delete
  originals without the byte check, and never by hand.
- Reruns resume from `python -m glamio.pack status Run<box>`. A finished pack
  is re-checked against the originals; an interrupted prune
  (`partially-pruned`) is finished by a resubmission (not with `PRUNE=no`),
  which first byte-checks the originals still present;
  `busy` (a `.partial` left by a running or crashed pack) and `inconsistent`
  stop the job without changing anything. Investigate those with the user
  instead of deleting files to get past them.
- `pack` refuses another subdirectory, symbolic links, PM particle or
  checkpoint files, step-tagged files of an unknown product, text catalogues
  (`CATALOGS/*.DAT`, written when `CatalogueFormat` is `text` or `both`), a
  `.bin` without its `.meta`, a wrong box tag, a step missing from
  `OutputSchedule.dat` and a `packed.md5` it did not write. Resolve those with
  the user first, never by deleting them to get past it.
- Restore native files with `python -m glamio.pack unpack Run<box> --out DIR`
  (byte-identical; `micromamba run -n cosemu python3 -I`, with
  `MAMBA_USE_LOCKFILES=false` in jobs). For an L1024 box this writes 75 GB in
  about 15 min, so run it under Slurm, not on a login node. To copy a packed
  box elsewhere (e.g. to a collaborator's account), copy the whole run directory:
  `Init.dat`, `Setup.dat`, `PkTable.dat`, `TableSeeds.dat`,
  `inputs_provenance.json` and `bin/` carry the provenance and are not packed.
  Then run `md5sum -c packed.md5` in each `Run<box>/` at the destination.
- Cost for an L1024 Np2048 box (75 GB): pack 20 min, check 15-30 min, prune
  30 min, peak RSS 0.8 GB (sacct MaxRSS reports the 8 GB limit because it
  counts page cache). The check needs free scratch equal to the box. An L256
  pilot takes 30-100 s (Lustre load varies) and may run on the login node (light tier).
- Record the pack/prune job IDs, file counts and bytes in the campaign's
  provenance (e.g. the `DEGRACE-MGGLAM/README.md` row of the box).

## Cleanup and file quota

- At the end of work, reduce unnecessary file and directory counts. Consolidate
  completed build/test scratch and finished-job launch files when they are no
  longer needed unpacked; verify archived contents before removing originals.
- Preserve scientific outputs, input configurations, provenance and everything
  required by running or queued jobs. Untracked files are not automatically
  disposable. Record what was removed and how retained archives can be restored.
