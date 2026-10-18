# Review prompts for branch `bdm-binary-catalogue` (1 commit, `cz..bdm-binary-catalogue`)

Five independent reviews before the merge into `cz`. Each prompt is
self-contained: paste the COMMON PREAMBLE followed by one REVIEW block into a
fresh session or subagent. Findings go to
`docs/reviews/bdm-binary-catalogue/REVIEW-<k>-<topic>-findings.md`. Reviews 1,
2, 3 and 5 are read-only; review 4 builds and runs the tests in a detached
worktree prepared by the author (never in the main checkout).

---------------------------------------------------------------------------
## COMMON PREAMBLE (prepend to every review)

You are reviewing code on the COSMA HPC system (Durham). Repository:
/cosma8/data/dp203/dc-ruan1/mgglam_claude/MG-GLAM (MG-GLAM, a Fortran
particle-mesh N-body code for LCDM and modified gravity, plus tools). Branch
`bdm-binary-catalogue` is checked out; the mainline is `cz` (NOT `main`).
Review range: `git diff cz..bdm-binary-catalogue` (one commit, 167fb1c,
17 files; `git show --stat 167fb1c`). The commit adds a binary form of the
BDM halo-finder catalogue and a Python reader package:

- `BDM.config` (the halo finder's own key = value file, parsed by
  `ReadParameters` in PMP2linker.f90) gains `CatalogueFormat = text | binary
  | both`, default `both`. Beside or instead of the legacy text
  `CATALOGS/CatshortV.<step>.<box>.DAT` (8 header lines, 24 columns, Fortran
  format `(3f11.4,3x,3f10.2,1p,2g12.4,g12.5,26g12.4)`), `WriteFiles` writes
  `CatshortV.<step>.<box>.bin`: headerless little-endian stream, one 96-byte
  record per halo = the same 24 values, 22 as real*4 and 2 as integer*4 (the
  row index and the distinct flag), same rows and order as the text; and a
  `key = value` sidecar `CatshortV.<step>.<box>.meta` written by
  `WriteCatalogueMeta`. The staged-file publication (hidden temporary created
  with STATUS=NEW, C rename into place) was generalised from one file to up
  to three (`BeginCataloguePublication`, `StageCatalogueFile`,
  `PublishCatalogue`), renamed in the order text, bin, meta.
- `python/glamio` (numpy only; `pip install -e python/` into the cosemu env):
  `read_catalogue(path)` for the text or the binary catalogue, a `Catalogue`
  class with column access and the epochs from the sidecar,
  `velocities_at_av()`; unit tests in `python/tests/`.
- `test/precise_z` (36^3 particles, 72^3 mesh, 150 Mpc/h, 100 steps):
  case12_frhalo now also checks the binary catalogue against the text one
  through glamio; new case13_catbin (binary only) and case4g_badcatfmt
  (bad value rejected). The author reports 194/194 checks passing.
- Docs: docs/degrace_companion_implementation_plan.md (sections 6, 7, 9, 11,
  16), DEGRACE-companion/README.md and gen_inputs.py (`--catalogue-format`,
  writes Run1/BDM.config), test/precise_z/README.md, python/README.md, two
  frames of docs/slides/degrace_companion_status/degrace_companion_status.tex.

Treat the commit message, the plan's section 16 entry and the READMEs as
claims to check, not as evidence.

Hard rules (shared system; the user is personally accountable for what you do):
- Do not edit tracked files in the main checkout, do not commit, do not
  rebuild binaries there, do not run simulations there, do not submit Slurm
  jobs, do not start background processes, do not use `micromamba run` to
  launch any Fortran executable (the env shadows a runtime library and the
  executable dies at start-up). Write exactly one file: your findings (path
  given below). Small Python checks are fine, run as
  `nice -n 19 /cosma/apps/durham/dc-ruan1/micromamba/envs/cosemu/bin/python3 -I -c '...'`
  (seconds, small). The `glamio` package is importable from that interpreter
  (editable install pointing at python/ of the main checkout).
- Never run recursive grep/find/du/ls -R over /cosma*, /snap*, the home
  directory or the data root /cosma8/data/dp203/dc-ruan1. Search only inside
  the repository paths named below with explicit pathspecs
  (`git grep -n <pattern> -- '*.f90'`, `grep -n` on a named file). Do not run
  `git status` or `git grep` without a pathspec: the checkout holds terabytes
  of untracked run directories. Do not read or touch ~/.ssh.
- Use `git diff cz..bdm-binary-catalogue -- <paths>` and `git show 167fb1c --
  <path>` to see exactly what changed, and read the surrounding unchanged
  code when correctness depends on it (`sed -n 'a,bp' <file>`).
- Real data you may read (read-only, do not copy): the LCDM pilot text
  catalogue
  DEGRACE-companion/runs/Node_002_LCDM_L1024Np1024Ng2048/Run1/CATALOGS/CatshortV.0181.0001.DAT
  (220 MB, 777,651 haloes, z = 0) and its directory's `outputB.0181.dat`,
  `../timing.log`, `../BDM.config`. There is no binary catalogue of a real
  run yet; the only binary catalogues are the ones the test suite produces.

Deliverable: write the findings file with these sections. (1) Summary, 3-6
lines, ending with a merge recommendation: ready / ready after the listed
fixes / not ready. (2) Findings table with columns ID, Severity, File:line,
Description, Evidence, Suggested fix. Severity: Blocker = wrong results, data
loss or a crash in production use; Major = wrong results in some
configurations, silent misuse, resource blow-up; Minor = robustness or
maintainability; Nit. (3) Verified OK: a bullet list of the specific
properties you checked and found correct, so nobody re-reviews them.
(4) Open questions for the author. Be concrete: cite file:line, quote the
offending expression, give a failing input where you can. Assume bugs exist
and look for them; do not pad with generic advice. In your final message
return only the findings file path, the counts per severity and one line per
Blocker or Major.

---------------------------------------------------------------------------
## REVIEW 1: Fortran writer, config parser and staged publication

Findings file: docs/reviews/bdm-binary-catalogue/REVIEW-1-fortran-findings.md

Scope: the diff of PMP2linker.f90 (module Structures declarations,
ReadParameters, BeginCataloguePublication, StageCatalogueFile,
PublishCatalogue, SetParameters header loop, WriteFiles, WriteCatalogueMeta)
and PMP2mod_extras.f90 (fRHaloClose sidecar line). Read the unchanged
callers: subroutine BDM (PMP2linker.f90, the sequence ReadParameters ->
SetParameters -> ... -> WriteFiles -> fRHaloClose -> close(12), including
the early `Nmaxima == 0` return), RescaleCoords, the Tools module
(PMP2mod_tools.f90: AEXPN, ASTEP, ISTEP, Nrealization, Box, NGRID, NROW,
Om, OmL, hubble, HEADER) and the BDM call site in PMP2mod_analyze.f90.

Check at least:
- Record layout: the binary `write(CatUnitBin)` item list is exactly the 24
  values of the text write in the same order; each item is 4 bytes
  (`int(iHalo,4)`, `0_4`, every real*4 expression such as `1.e3*RadRms(ip)`
  and the real*4 locals x, y, z, rr, Vrms, Cvir, aNpart, VirRat); no default
  real*8 promotion, no integer*8, no derived-type padding; stream unformatted
  writes carry no record markers. Confirm `convert='little_endian'` is
  honoured for stream access by ifx and is consistent with the fR_halo writer.
- Parser: the new `case('catalogueformat', ...)` branch (io reset, lowercase
  conversion, length guard, accepted values), interaction with the generic
  `if (io /= 0)` check after the select, the "one scalar value" guard, the
  default `both` reset at every ReadParameters call, the default-file writer
  (format 30) and the parameter echo. Would an old binary reading a new
  BDM.config (with the key) fail? Would a new binary reading an old file
  change behaviour? Is `ValidateParameters` still complete?
- Staged publication: CatPubN/CatPubUnit/CatPubFinal/CatPubStaged bookkeeping
  across repeated BDM calls in one run (every class-1 output), the
  `CataloguePublicationPending` state machine, the `close(unit)` before
  STATUS=NEW, unit numbers 12/32/33 free of conflicts elsewhere
  (`git grep -n -E 'open *\( *(12|32|33)' -- '*.f90'`), the inquire/close/
  rename loops in PublishCatalogue, behaviour when the text is off (unit 12
  never opened: the header loop guard in SetParameters, `close(12)` in BDM,
  anything else writing to unit 12 during BDM), path length limits
  (character(256), suffix 96), the `+5` in the length check, and what is
  left on disk if the run dies between staging and publication.
- WriteCatalogueMeta: every value comes from the right variable at the right
  time (MassMin after the `max(MassMin, 20*MassOne)` adjustment; Ovdens
  relative to the mean density per SetOverdensity; a_v = AEXPN - ASTEP/2 at
  the moment BDM runs, i.e. after Move/AddTime of this step; `da = ASTEP`),
  format widths cannot overflow (asterisks) for production values (ISTEP up
  to 9999, Box 1024-1800, NGRID 8192, MassMin 1e11, Ovdens 100-400), the
  `columns`/`dtypes`/`units` lines agree with the record layout, and the
  `format` and `centre` free-text lines are correct statements about the
  code (Mbound = Mvir(ip) the bound mass; Mtot = Mtotal(ip); Xoff in units
  of the aperture; Rvir = 1e3*Rvir(ip) in kpc/h; velocities from
  RescaleCoords with Vscale = 100 Xscale/AEXPN).
- Error paths: every new `error stop` / ConfigurationError is reachable and
  the message matches what test case4g asserts ("CatalogueFormat must be
  text, binary or both"); no path can publish a `.meta` without its `.bin`.
- Regression safety: with `CatalogueFormat = text` the text catalogue is
  byte-identical to the pre-change output (same unit, same format, same
  header) and the timing/outputB logs are unchanged.
- Threading: WriteFiles runs serially after the OpenMP regions; confirm no
  new shared state is touched inside a parallel region.

---------------------------------------------------------------------------
## REVIEW 2: format specification and sidecar semantics

Findings file: docs/reviews/bdm-binary-catalogue/REVIEW-2-format-findings.md

Scope: the meaning and consistency of every documented quantity, across the
Fortran writer (WriteCatalogueMeta and the text header writer in
PMP2linker.f90), the reader (python/glamio/catalogue.py: CATALOGUE_DTYPE,
COLUMN_UNITS, read_catalogue_text's header parsing, velocities_at_av), and
the prose: python/README.md, DEGRACE-companion/README.md, test/precise_z/
README.md, the plan (docs/degrace_companion_implementation_plan.md sections
3, 6, 7, 9, 11 and the 2026-10-17 entries of section 16) and the deck source
docs/slides/degrace_companion_status/degrace_companion_status.tex (the
"What a box writes" table and the "Definitions: BDM halo catalogue" frame).

Check at least, from the code and not from the prose:
- Units and definitions of all 24 columns as computed in WriteFiles / GetHalo
  / ParametersDistinct (PMP2linker.f90): positions in [0, Box) Mpc/h after
  RescaleCoords; velocities in km/s with the 1/AEXPN scale; Mbound vs Mtot;
  Rvir in kpc/h; Vrms = sqrt(2 Ekin/Mvir); Vcirc = VmaxM; Cvir from
  Concentration(); Npart = Mvir/MassOne or the membership count; distinct
  always 0; Xoff = |offset|/aperture; VirRat = 2K/Ep - 1; LambdaM; RadRms in
  kpc/h; axis ratios and axis direction. Are the `units` line, COLUMN_UNITS
  and the README table right? Is "virial_ratio = 2K/|Ep| - 1" the actual
  expression (sign of EpotM)?
- The velocity epoch claim. In the PM main loop (PMP2main.f90, the Nlist
  block near line 78 and the Move/AddTime order), at the moment BDM is called
  are positions at AEXPN and velocities at AEXPN - ASTEP/2 with ASTEP the
  step just taken? Compare with how PMP2mod_extras.f90 records a_v for the
  lens planes and momentum dumps (`a_v = AEXPN - ASTEP/2.`). State precisely
  what `velocities_at_av()` returns and whether the factor a/a_v is the
  correct first-order correction (what does the 1/a conversion actually
  assume about the velocity variable?).
- Overdensity: SetOverdensity's `density` for iVirial = 0..3; is the sidecar
  text "(relative to the mean matter density at this epoch)" true for every
  iVirial, and is the deck's "338 at z = 0, 203 at z = 1" consistent with
  the formula for Omega_m = 0.3089?
- MassMin in the sidecar vs the selection actually applied (`Mvir < MassMin`
  and `aNpart < 20` cuts, and the domain cut) and the text "rows also need
  Npart >= 20".
- Text-vs-binary precision statements: derive the printed precision of each
  text column from the Fortran format (F vs E form of G editing with the 1p
  scale factor, significant digits) and check the tolerances used in
  test/precise_z/asserts.py (TEXT_TOL) are tight enough to catch a wrong
  column yet loose enough never to fail on a legitimate value, including
  values that switch between F and E form (|x| < 0.1 or >= 1e4 for g12.4)
  and negative values.
- Sizes quoted in the plan (section 9 rows for the catalogue), the deck
  table and the commit message against the record size and the pilot's row
  count (777,651 rows: 283 B/row text, 96 B/row binary).
- Naming: the sidecar `product = CatshortV` vs the file label for iVirial
  0/2/3 (W/M/A); does glamio detect the label? Does anything assume `V`?
- Every key the README lists as written by the Fortran is actually written,
  and vice versa; the `fR_halo` sidecar's new `catalogue_bin` line.

---------------------------------------------------------------------------
## REVIEW 3: the glamio Python package

Findings file: docs/reviews/bdm-binary-catalogue/REVIEW-3-glamio-findings.md

Scope: python/pyproject.toml, python/glamio/__init__.py, meta.py,
catalogue.py, python/tests/test_catalogue.py, python/README.md, the
.gitignore additions. Cross-check against the Fortran writer (PMP2linker.f90
WriteFiles and WriteCatalogueMeta) and the text format.

Check at least:
- CATALOGUE_DTYPE field order, names and kinds against the Fortran write
  statement; itemsize 96; endianness explicit (`<`); the layout check
  `_check_layout` (what mismatches does it catch, which does it miss, e.g. a
  reordered `columns` line with the same names, an extra column, `>f4`).
- Binary reader: size checks (not a multiple of 96, rows mismatch,
  record_bytes), zero-row files, the stem resolution rules (bin/meta
  preferred, .DAT/.dat fallback, missing files), `np.fromfile` vs memory
  mapping for 300 MB files, no copies hidden in `__getitem__`.
- Text reader: the 8-line header parse against the real pilot file and the
  test fixture (`A =`, `Step =`, `I =`, `Nrow`, `Ngrid`, `Omega_0`,
  `Omega_L`, `hubble`, `buffer width`, mass of the smallest particle,
  overdensity limit); a_v = a - Step/2 formatting (6 decimals: enough?);
  np.loadtxt with `ndmin=2` on a header-only file; rows with != 24 columns;
  Fortran `g12.4` outputs such as `1.2649E+13`, `-0.3133`, `4.0263E-02`;
  does a `Catshort?` label other than V parse; is the float32 cast of
  positions lossless relative to the 4 printed decimals?
- `velocities_at_av()` and the properties `a`, `z`, `a_v`, `step`, `box`
  (`meta_float` semantics: first token, D exponents, non-numeric -> error).
- Packaging: pyproject validity for setuptools >= 64 editable installs,
  `requires-python`, package discovery, whether tests are excluded from the
  wheel, version string in one place.
- Tests: do they exercise the real failure modes (wrong endianness, swapped
  fields, truncated file, text/binary mismatch) or only the happy path? Are
  fixtures written to a sane temporary location (GLAMIO_TEST_TMPDIR /
  TMPDIR, never /tmp on a login node with data)? Can they run under pytest
  and as a script? Run them once: `cd python && nice -n 19
  /cosma/apps/durham/dc-ruan1/micromamba/envs/cosemu/bin/python3 -m pytest -q
  -p no:cacheprovider tests` and delete any __pycache__ you create under
  python/.
- Code quality: typing, docstrings that state formulas, error messages that
  name the file, no prints in library code, no dead code, numpy 2.x
  compatibility (the env has numpy 2.3 / Python 3.14).

---------------------------------------------------------------------------
## REVIEW 4: test suite changes, executed in a worktree

Findings file: docs/reviews/bdm-binary-catalogue/REVIEW-4-tests-findings.md

This review may build and run. Use ONLY the detached worktree the author
prepared at
/cosma8/data/dp203/dc-ruan1/scratch/review-bdm-binary-catalogue-20261018/wt
(a checkout of commit 167fb1c; tracked files only, nothing else). Do not
build or run in the main checkout. Environment for building and running:

    module purge; module load intel_comp/2024.2.0; export I_MPI_F90=ifx
    module load compiler-rt tbb compiler mpi
    export PATH=/cosma/apps/durham/dc-ruan1/micromamba/envs/cosemu/bin:$PATH
    export OMP_NUM_THREADS=4 SKIP_MODULES=1
    cd <worktree> && nice -n 19 make PMP2init PMP2start PMP2main   # ~10 s, serial (make -j breaks the module order)
    cd <worktree>/test/precise_z && nice -n 19 bash run_tests.sh   # ~45 s, 4 threads, < 100 MB

Everything you run must stay within this light tier (one process, <= 4
threads, minutes, well under 16 GB) and must be finished, with no background
processes, before you report. `run_tests.sh` deletes and recreates its
`work/` directory; `asserts.py work` can be rerun standalone on an existing
`work/` directory. The worktree's `python/glamio` is what `asserts.py`
imports (it inserts `../../python` at the front of sys.path).

Scope: the diff of test/precise_z/asserts.py, run_tests.sh and README.md,
plus the harness pieces they depend on (prep_case, set_outputs, run_init,
run_main, `oracle_full`, `Report`, `ERR_CASES`).

Check at least:
- Run the suite once and record the verdict lines for case12_frhalo,
  case13_catbin and case4g_badcatfmt and the totals (the author reports
  194/194).
- Discriminating power, by mutation on the produced files (no rebuild
  needed): after a run, corrupt `work/case12_frhalo/Run1/CATALOGS/
  CatshortV.<step>.0001.bin` in several ways (flip one byte in a float
  column, swap two records, truncate by one record, append 96 zero bytes,
  byte-swap one field), edit the `.meta` (rows off by one, `columns` with
  two names exchanged, `a_v` wrong), leave a fake staged file
  `.CatshortV.0098.0001.bin.tmp.1.2.3` in CATALOGS/, and rerun
  `nice -n 19 python3 asserts.py work` each time. Report which mutations
  are caught and by which check, and which are not. Restore the files (or
  rerun the suite) between mutations.
- TEXT_TOL: derive the printed precision of each column from the Fortran
  format and judge whether the tolerances could mask a real column error or
  could fail spuriously (values crossing between F and E form of g12.4,
  negatives, zeros, |x| < 1e-12).
- case13: does it really prove the text-off path (no `.DAT`, header loop
  skipped, publication of two files)? Does `text_catalogue = none` get
  checked strictly? Would a run that silently wrote the `.DAT` anyway be
  caught?
- case4g: is the run really stopped by the config error at the first
  catalogue step (check main.log) and could a false pass arise if the run
  died earlier for another reason? Is the forbidden glob
  `Run1/CATALOGS/Catshort*` the right one?
- Harness hygiene: the `echo` inserted before appending to the copied
  test/BDM.config (which has no trailing newline): is that fix robust, and
  should test/BDM.config itself get a final newline? Are the Init.dat `sed`
  edits anchored exactly (they match fixed-width lines)? Any case that could
  pass vacuously (zero haloes)? The totals in README.md.
- Coverage gaps worth a new check (e.g. a run that sets `CatalogueFormat =
  text` and asserts byte-identical output to the pre-change catalogue; a
  `both` run compared against an independent reader without glamio).
- At the end delete `<worktree>/test/precise_z/work`, any `__pycache__`
  under the worktree, and the built objects are fine to leave; do NOT
  remove the worktree itself (the author will).

---------------------------------------------------------------------------
## REVIEW 5: downstream impact, operations and repository hygiene

Findings file: docs/reviews/bdm-binary-catalogue/REVIEW-5-impact-findings.md

Scope: everything that consumes or configures the catalogue, and the
documentation of the change. Read-only.

Check at least:
- Consumers of the text catalogue inside the repository and whether the
  default `both` or a later switch to `binary` breaks them: the Fortran
  tools (`git grep -n Catshort -- '*.f90'`: PMP2catalReform.f90,
  PMP2correlation.f90, Massfunction.f90, PMP2catalogCorrect.f90,
  PMP2haloStat.f90), the Python scripts (`git grep -ln CatshortV --
  '*.py'`), the companion checker DEGRACE-companion/check_run.py (does it
  count catalogue files by glob and would the new `.bin`/`.meta` confuse
  its counts or sizes?), and the external halocat package at
  /cosma8/data/dp203/dc-ruan1/halocat (read only halocat/halocat/io.py,
  pipeline.py, config.py: how it finds and parses `CatshortV.*.DAT`, what
  an added `.bin`/`.meta` next to it would do, what it would need to read
  the binary).
- The default: `both` for every run whose BDM.config lacks the key, i.e.
  every existing campaign directory and every run that lets the code write
  its default BDM.config. Quantify the extra disk for the running or
  planned campaigns you can identify from the repository (DEGRACE-companion
  pilots, the L1800 DEGRACE suites, fid_*), and say whether `both` or
  `text` is the safer default for the merge.
- gen_inputs.py: the new `bdm_config()` content versus the pilots' actual
  Run1/BDM.config (DEGRACE-companion/runs/Node_002_*/Run1/BDM.config) and
  test/BDM.config; provenance (is `catalogue_format` recorded in
  inputs_provenance.json?); the `--catalogue-format` help text; does the
  existing submit.sh or README still describe how BDM.config reaches a run?
- Restart and crash behaviour: a run killed after staging and before
  publication leaves hidden `.Catshort*.tmp.*` files; on restart from a
  checkpoint the output step is redone and new temporaries are created
  (different pid/stamp). Is there any cleanup, and does the plan or README
  tell the operator to remove stale temporaries? Does `check_run.py` notice
  them?
- Documentation accuracy: plan sections 6, 7 (the new BDM.config table),
  9 (sizes), 11 (items 11 and 12), 16 (the entry); DEGRACE-companion/
  README.md; test/precise_z/README.md; python/README.md; the deck's two
  frames. Flag every statement that the code does not support.
- Repository rules (AGENTS.md, the user's CLAUDE.md): nothing committed
  that should not be (egg-info, caches, PDFs are intentionally tracked for
  the deck), .gitignore additions, file naming, no stray scratch left in
  the tree, commit message accuracy (e.g. "194/194 checks pass").
- Anything the commit changes for runs that do not use BDM at all
  (`Find BDM halos = 0`): confirm nothing in the diff executes then.
