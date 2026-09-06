# Atomic catalogue publication

Independent review of core repair `cd728d6c4da3fce6794d04e0973572a075cd13e6`
and peak integration `0720cd42b2a97f5f6a8353cf07d798a150f8cf27` found a concrete
publication hazard. A valid candidate followed by a nonfinite kinetic energy
made `WriteFiles` fail, but the final catalogue already contained eight header
lines plus a valid 24-column row prefix. An existing catalogue had been replaced.
Configuration/header creation also exposed an incomplete final pathname before
any finder work. No further actionable defect was identified in the exact-ID
merge sorting, Jacobi eigenpair calculation, concentration inversion or default
finder-specific precise compilation flags during this focused review.

`BeginCataloguePublication(final_path)` now exclusively creates a hidden staged
file in the final catalogue's own directory. Process ID, clock and a collision
counter distinguish concurrent writers; Fortran `STATUS='NEW'` prevents an
existing stage from being truncated and respects the process umask. Headers
and rows use this stage on unit 12.

After the complete writer succeeds, `PublishCatalogue` verifies the open stage,
checks the close result, then uses C `rename` through `ISO_C_BINDING` to replace
the final path atomically on the supported POSIX filesystems. The prior catalogue
remains visible until replacement. A failed or interrupted run can leave a
hidden `.Catshort*.DAT.tmp.*` stage; it cannot expose that stage under the final
catalogue name. A leftover stage can be removed after its writer has exited.
Successful publication removes the staged pathname. No receipt is required.

The final filenames, eight header lines and 24 columns are unchanged. A valid
empty run publishes the eight-line header. Direct `WriteFiles` callers must
start with `BeginCataloguePublication`; header-only test callers explicitly
finish with `PublishCatalogue`. The focused peaks/core drivers have minimal
adapters for this contract and extract both helpers.

Validation in `publication-results.json` includes:

- 26 GNU checked/O3 publication experiments: valid and empty replacement; late
  invalid candidate, post-header interruption, prematurely closed stage and
  rename failure; preservation of previously generated valid catalogue bytes;
  concurrent writers with unique stages; normal umask permissions.
- All 136 existing peak/configuration/empty-BDM assertions and the core identity,
  eigenpair, concentration and writer suite, after the caller adapters.
- Native ifx O3 precise compilation and valid, late-invalid and empty controls.

Run the publication regressions with:

```sh
micromamba run -n cosemu python3 -B BDM-refine/repairs/20260906/tests/publication_regression.py
```

Tests build and run in temporary directories and remove those directories on
exit. No root checkout or historical audit result was changed during review.
