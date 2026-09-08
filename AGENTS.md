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
