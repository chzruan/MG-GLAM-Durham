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
- Use `ic2pm` built from commit ee44100 or later (branch `ic2pm-halfstep`
  until merged into `cz`). Its default velocity epoch is `half` (velocities
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
- Record the actual generator/version, LPT order, starting redshift, seed
  and input configuration in each campaign's provenance. Resolution
  comparisons also require verified shared modes and consistent normalization.
- Preserve the recorded IC method of existing runs when resuming them.
  The already launched
  [20260908 BDM convergence campaign](BDM-refine/validation/20260908-convergence/README.md)
  uses a frozen native first-order variant and is an exception, not a template
  for the default IC choice in new campaigns.

## Cleanup and file quota

- At the end of work, reduce unnecessary file and directory counts. Consolidate
  completed build/test scratch and finished-job launch files when they are no
  longer needed unpacked; verify archived contents before removing originals.
- Preserve scientific outputs, input configurations, provenance and everything
  required by running or queued jobs. Untracked files are not automatically
  disposable. Record what was removed and how retained archives can be restored.
