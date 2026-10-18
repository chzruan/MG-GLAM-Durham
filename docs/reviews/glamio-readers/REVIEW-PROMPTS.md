# Review prompts for branch `glamio-readers` (1 commit, `cz..glamio-readers`)

Five independent reviews before the merge into `cz`. Each prompt is
self-contained: paste the COMMON PREAMBLE followed by one REVIEW block into a
fresh session or subagent. Findings go to
`docs/reviews/glamio-readers/REVIEW-<k>-<topic>-findings.md`. Reviews 1, 2,
3 and 5 are read-only (small Python checks allowed); review 4 builds and runs
the tests in a detached worktree prepared by the author (never in the main
checkout).

---------------------------------------------------------------------------
## COMMON PREAMBLE (prepend to every review)

You are reviewing code on the COSMA HPC system (Durham). Repository:
/cosma8/data/dp203/dc-ruan1/mgglam_claude/MG-GLAM (MG-GLAM, a Fortran
particle-mesh N-body code for LCDM and modified gravity, plus tools). Branch
`glamio-readers` is checked out; the mainline is `cz` (NOT `main`). Review
range: `git diff cz..glamio-readers` (one commit, 0a6ecff, 16 files;
`git show --stat 0a6ecff`). The commit extends the numpy-only Python package
`python/glamio` (version 0.1.0 -> 0.2.0; it previously read only the BDM halo
catalogues) with one reader per remaining GLAM output product. Every product
except the schedule and the power spectra is a headerless little-endian
`.bin` file beside a `key = value` text sidecar `.meta` written by the
Fortran code; the readers must reproduce the Fortran layouts, units and
epochs exactly:

- `glamio/product.py`: `MetaProduct` base (sidecar scalars `a`, `z`, `a_v`,
  `a_x`, `box`, `step`, `realization`, `n_grid`, `n_row`, `n_p`, `omega_m`,
  `omega_l`, `hubble`, `header`) and helpers `resolve_pair`, `load_meta`
  (product line + required keys), `check_size` (exact byte count),
  `load_array` (np.fromfile or read-only np.memmap). `Catalogue` now inherits
  `MetaProduct` (glamio/catalogue.py; behaviour claimed unchanged).
- `glamio/schedule.py`: `read_schedule` for `OutputSchedule.dat` (writer:
  `WriteOutputSchedule`, PMP2mod_extras.f90).
- `glamio/lensplanes.py`: `read_lensplanes` for `LensPlanes.<step>.<box>`
  (writer: `LensPlanes`, PMP2mod_extras.f90): `planes[plane, j, i]`,
  `plane(axis, slab)` transposed to `[i, j]`, `slabs`, `axis_means`.
- `glamio/modes.py`: `read_modes` for `deltam`, `mom_x/y/z` and the
  companion `*_deltaL` cubes (writer: `DumpModes`, PMP2mod_extras.f90):
  `data[ix+n, iy+n, iz]` complex64, `mode(ix, iy, iz)` with Hermitian
  completion for iz < 0, `wavenumbers`, `k_magnitude`, `power()` (the
  shell average claimed to equal `PowerDM.R` bin by bin).
- `glamio/subsample.py`: `read_subsample` for `Subsample.<step>.<box>`
  (writer: `WriteSubsample`, PMP2mod_subsample.f90, with glam_hash.c):
  `positions()`, `velocities()`, `ids()` from the splitmix64 rule,
  `splitmix64`, `select_ids`, `lattice` tag.
- `glamio/frhalo.py`: `read_frhalo` for `fR_halo.<step>.<box>` (writer:
  `fRHaloOpen/Row/Close`, PMP2mod_extras.f90) and `FRHalo.catalogue()`.
- `glamio/powerdm.py`: `read_powerdm` for `PowerDM.R/log/Z.*.dat` (writer:
  `GetPower`, PMP2mod_power.f90): table, epoch line, `pk_deconvolved()`.
- `glamio/snapshot.py`: `read_pm_header` (the `PMcrd` record) and
  `Snapshot` (`iter_pages`, `read`, units `code`/`physical`) for the PM
  particle files `PMcrs<i>` (writer/reader: `WriteDataPM`/`ReadDataPM`,
  PMP2mod_tools.f90; the build uses `-convert big_endian`, makefile line 6).
- Tests: `python/tests/test_readers.py` (8 tests with synthetic files;
  `python/tests/test_catalogue.py` only had its version assertion relaxed)
  and six new `glamio` checks in `test/precise_z/asserts.py` that compare
  every reader with the suite's own independent decoders (`read_schedule`,
  `read_planes`, `read_modes`, `read_subsample`, `read_powerdm`,
  `read_aexpn`, `binned_power` in the same file) on the real outputs of
  case10_lens72, case11_modes and case12_frhalo (36^3 particles, 72^3 mesh,
  150 Mpc/h). The author reports 221/221 suite checks and 15/15 unit tests.
- Docs: python/README.md ("Other products" section), test/precise_z/README.md
  (one row), docs/degrace_companion_implementation_plan.md (section 16
  entry dated 2026-10-18).

Treat the commit message, the docstrings, the README and the plan entry as
claims to check, not as evidence.

Hard rules (shared system; the user is personally accountable for what you do):
- Do not edit tracked files in the main checkout, do not commit, do not
  rebuild binaries there, do not run simulations there, do not submit Slurm
  jobs, do not start background processes, do not use `micromamba run` to
  launch any Fortran executable (the env shadows a runtime library and the
  executable dies at start-up). Write exactly one file: your findings (path
  given below). Small Python checks are fine, run as
  `nice -n 19 /cosma/apps/durham/dc-ruan1/micromamba/envs/cosemu/bin/python3 -I -c '...'`
  (seconds, well under 1 GB). The `glamio` package of the main checkout is
  importable from that interpreter (editable install pointing at python/).
- Never run recursive grep/find/du/ls -R over /cosma*, /snap*, the home
  directory or the data root /cosma8/data/dp203/dc-ruan1. Search only inside
  the repository paths named below with explicit pathspecs
  (`git grep -n <pattern> -- '*.f90'`, `grep -n` on a named file). Do not run
  `git status` or `git grep` without a pathspec: the checkout holds terabytes
  of untracked run directories. Do not read or touch ~/.ssh.
- Use `git diff cz..glamio-readers -- <paths>` and `git show 0a6ecff --
  <path>` to see exactly what changed, and read the surrounding unchanged
  code when correctness depends on it (`sed -n 'a,bp' <file>`).
- Real files you may read (read-only, do not copy, do not delete): the
  suite's outputs of the author's last run in the main checkout,
  test/precise_z/work/case10_lens72/Run1 (LensPlanes.*, PMcrd.*/PMcrs0.*,
  PowerDM.*, OutputSchedule.dat), test/precise_z/work/case11_modes/Run1
  (deltam.*, mom_?.*, Subsample.*, PMcrd.*/PMcrs0.*, PowerDM.*,
  OutputSchedule.dat) and test/precise_z/work/case12_frhalo/Run1 (fR_halo.*,
  CATALOGS/CatshortV.*); and the LCDM pilot's full-size lens planes
  DEGRACE-companion/runs/Node_002_LCDM_L1024Np1024Ng2048/Run1/LensPlanes.0057.0001.{bin,meta}
  (201 MB, 12 planes of 2048^2) and its OutputSchedule.dat and PowerDM.*
  files. The pilot has no snapshots, subsamples or mode cubes.

Deliverable: write the findings file with these sections. (1) Summary, 3-6
lines, ending with a merge recommendation: ready / ready after the listed
fixes / not ready. (2) Findings table with columns ID, Severity, File:line,
Description, Evidence, Suggested fix. Severity: Blocker = wrong results or a
crash in ordinary use; Major = wrong results in some configurations, silent
misuse, resource blow-up on a shared login node; Minor = robustness or
maintainability; Nit. (3) Verified OK: a bullet list of the specific
properties you checked and found correct, so nobody re-reviews them.
(4) Open questions for the author. Be concrete: cite file:line, quote the
offending expression, give a failing input where you can. Assume bugs exist
and look for them; do not pad with generic advice. In your final message
return only the findings file path, the counts per severity and one line per
Blocker or Major.

---------------------------------------------------------------------------
## REVIEW 1: snapshot, subsample, schedule and PowerDM readers against the Fortran writers

Findings file: docs/reviews/glamio-readers/REVIEW-1-snapshot-subsample-findings.md

Scope: glamio/snapshot.py, glamio/subsample.py, glamio/schedule.py,
glamio/powerdm.py and the parts of glamio/product.py they use, checked
line by line against the Fortran that writes the files:
PMP2mod_tools.f90 (the Tools module declarations, lines 1-70: `HEADER`
character*45, `extras(100)`, the integer kinds of `Nparticles`, `NROW`,
`NGRID`, `ISTEP`; `ReadDataPM` lines 283-427; `WriteDataPM` lines 433-530),
PMP2main.f90 (where `WriteDataPM` is called relative to the `ASTEP` update at
line 405 and the `AEXPN = AEXPN + ASTEP` at line 482; the checkpoint call
with iFlag = 0 at line 75 and the output calls at 83 and 101),
PMP2mod_subsample.f90 (whole module, 274 lines; plus glam_hash.c),
PMP2mod_extras.f90 lines 53-81 (`WriteOutputSchedule`) and
PMP2mod_power.f90 lines 1-110 (`GetPower`), ic2pm.f90 lines 190-205 and
335-350 (the header tags `v at a-da/2`, `IDx`/`IDz`), makefile line 6
(`-convert big_endian`, no `-assume byterecl`).

Check at least:
- The `PMcrd` record: field order, kinds and byte count of the Fortran
  `write (4) HEADER, AEXPN, ..., Nparticles, extras` against
  `snapshot._FIELDS` (45 + 4x4 + 4 + 7x4 + 4x4 + 3x4 + 8 + 400 = 529?); is
  `Nparticles` integer*8 in the Tools module? Is the record really a single
  sequential record with 4-byte markers in this build (ifx, `-convert
  big_endian`, default record markers)? Is the byte-order detection from the
  marker sound, and the `extras(100) = Box` claim (ic2pm.f90:197,
  PMP2start.f90:156)?
- `PMcrs` pages: `RECL=NACCESS` with `Naccess = Nrecord*6` and no
  `-assume byterecl` means RECL in 4-byte units, i.e. 24 MiB records of
  `Xb, Yb, Zb, VXb, VYb, VZb` each `Nrecord = 1024**2` real*4: confirm the
  reader's `record_bytes = 6*page*4`, the file split at `Nrecpage = 512`
  records (`ifile = (i-1)/Nrecpage`), the file names `PMcrs<i>.<step>.DAT`
  with `i0` (no zero padding) and `PMcrs<i>.DAT` for checkpoints, and that
  only the first `NinPage` entries of the last page are valid (the writer
  copies only `jfirst..jlast` into the buffers; the rest holds the previous
  page's values). Does the Python reader drop exactly that padding, and does
  it handle `Nparticles` being an exact multiple of the page?
- Coordinates: the writer stores grid units in `[1, NGRID+1)`; `ReadDataPM`
  clamps `INT(Xb)==Ngrid+1` to `Xb - 1e-3`. The Python `units='physical'`
  wraps `(X-1) Box/NGRID` modulo Box instead; is that a correct treatment of
  a particle at exactly `NGRID+1`, and is `units='code'` documented as raw?
- Velocities: the Python claims `v [km/s] = V_code (100 Box/NGRID)/a_v` with
  `a_v = AEXPN - ASTEP/2` at the time of writing. Verify against
  `WriteSubsample` (`av = AEXPN - ASTEP/2`, `vscl = 100 Box/NGRID/av`) and
  against the order of operations in PMP2main.f90: when `WriteDataPM(1,...)`
  runs, does `ASTEP` hold the step just taken (line 405) so that the header's
  `ASTEP` gives the same `a_v`? Also for the checkpoint (`WriteDataPM(0,...)`
  at line 75) and for the initial `PMcrd.DAT` from ic2pm (`half` mode: the
  header says `v at a-da/2`; what are `ASTEP` and `AEXPN` there, and is the
  reader's `a_v` right for an IC file?). The legacy catalogue divides by
  `AEXPN`; the reader offers `velocity_epoch='a'` for that.
- `PMHeader.particle_mass()` uses `RHO_CRIT = 2.7754e11`; the finder uses
  `Dscale = 2.774e+11*(Box/NROW)**3` (PMP2linker.f90:2066) and the catalogue
  sidecar records the finder's value. Is a 0.05 % discrepancy acceptable in a
  reader, or should the reader use the code's constant / drop the method?
- Subsample: decode formulas (`x = (q + 32768 + 0.5) Box/65536`,
  `v = s dv_quant`), the int16 ranges (`q - 32768` so that q = 65535 maps to
  32767; velocities clipped to +-32767), the file order (six contiguous
  arrays), `n_sel` from the sidecar vs the byte count, the `SALT` parse from
  `0x5EED00DE6C0A2026` (`int(tok, 0)`), `splitmix64` against glam_hash.c
  (constants, shift amounts, uint64 wrap-around in numpy), `select_ids`
  blockwise equality, the claim "ID = index - 1 in lattice order" and the
  `lattice` tag parse (`IDx`/`IDz`/`ID?` substrings of the 45-character
  header; could `ID?` match something else?). Memory of `ids()` for a
  2048^3 run at M = 1 (`np.arange(8.6e9, dtype=uint64)` = 69 GB) and at
  M = 100 (`select_ids` over 512 blocks of 16.7M).
- Schedule: the Fortran format `(i6,2es16.8,f13.4,2i6)` vs the whitespace
  split (any field that can run into its neighbour: step >= 100000? chi
  >= 1e8?), the `Om= ... OmL=` regex against the `f8.5` writes, `class`
  `max(1, Nclass)`, strict step ordering.
- PowerDM: the row format `(i6,f9.5,es11.4,7es13.5)` and the epoch line
  `(a,f8.5,a,i4,a,f8.3,a,i4,a,f8.2)` vs the regex `_EPOCH_RE` (what happens
  for `Ngrid >= 10000` or `Step >= 10000`, printed as asterisks; negative
  redshift at a > 1; `kNyq >= 100`); the `log` file skips empty bins; the
  `Z` file's caption and 6 columns; whether the first line (the run header,
  which can contain any characters) can be mistaken for the caption; the
  `cic_window_sq` formula vs what the plan/deck call the shell-averaged CIC
  window (docs/degrace_companion_implementation_plan.md, search `W_CIC` or
  `sin^2`); `binning` inferred from the file name.
- `product.py`: `resolve_pair` on a path ending in `.bin.meta` or on a
  directory, `load_meta` with a `product` value that carries a comment,
  `check_size` on a 0-row file, `load_array(mmap=True)` with count 0.

---------------------------------------------------------------------------
## REVIEW 2: lens planes, mode cubes and fR_halo readers against the Fortran writers

Findings file: docs/reviews/glamio-readers/REVIEW-2-planes-modes-findings.md

Scope: glamio/lensplanes.py, glamio/modes.py, glamio/frhalo.py and the parts
of glamio/product.py they use, checked line by line against
PMP2mod_extras.f90: `LensPlanes` (lines 85-220: the three accumulation
loops, `write (29) plane` and the sidecar lines), `FFT3Dforward` (lines
256-340: the packed layout and the sqrt(2) factors), `DumpModes` and
`PackedBasis` (lines 343-450), `MomentumDeposit` and `ExtraModes` (lines
453-545), `fRHaloOpen/Row/Close` (lines 553-665); PMP2mod_power.f90 /
POWERfft5 (PMP2mod_fft5.f90 or wherever `POWERfft5` bins `|delta_k|^2`:
`git grep -n "POWERfft5" -- '*.f90'`) for the binning that `Modes.power()`
claims to reproduce; and the `deltaL` cube writer of the companion tools
(`git grep -n -i "deltaL" -- DEGRACE-companion test/extras '*.py'`).

Check at least:
- Lens planes orientation: from the loops, X-planes accumulate
  `acc(ip, jp, s)` with `ip` from `iy` and `jp` from `iz`, so the Fortran
  `plane(i, j)` is `(y, z)`; `write (29) plane` streams column-major (i
  fastest). The reader reshapes the file to `[plane, j, i]` and `plane(axis,
  s)` returns `.T` = `[i, j]`. Confirm for all three axes (Y: `(x, z)`, Z:
  `(x, y)`), the plane order X1..XN, Y1..YN, Z1..ZN, the 1-based `slab`
  argument vs the sidecar's `slab s covers [(s-1), s) depth`, and the pixel
  convention `pixel (i,j) covers [(i-1), i) Box/N_lens` (grid node i sits at
  `(i-1) Box/NGRID`, block-averaged by `f = NGRID/N_lens`). Is `pixel_edges`
  / `slab_edges` consistent with that?
- Lens plane values: `plane = 1 + acc*wnorm` with `wnorm = 1/(nc f f)`; is
  `axis_means()` expected to be exactly 1 up to float32 rounding (mean of
  delta over the mesh is 0 by construction of DENSIT)? What does the
  reader do with `N_lens` absent / `N_lens = NGRID` (Fortran `Nl = NGRID`
  when `LensGrid <= 0`)?
- Mode cubes: `C(0:nm, -nm:nm, -nm:nm)` with the first index iz, `write (29)
  C` streams iz fastest, then iy, then ix; the reader's
  `reshape(2n+1, 2n+1, n+1)` -> `[ix+n, iy+n, iz]`. Check the sign/conjugate
  convention of `delta_k = N^-3 sum f(x) exp(-i k.x)` through `PackedBasis`
  (`w(2) = (0, -0.5 sign(m))`, "X_{-m} = conj") so that `Modes.mode(ix, iy,
  iz < 0) = conj(data[-ix, -iy, -iz])` is the right completion for the
  iz < 0 half (the field is real for deltam and for each momentum
  component); the `nm` clamp to `N/2 - 1`; `k = 0` zeroed; the Nyquist
  branch of `PackedBasis` (`abs(m) == N/2`) cannot be reached in the dump
  because of the clamp, right?
- `Modes.power()`: bins `nint(|m|)`, `iz > 0` weight 2 and `iz = 0` weight
  1, `Box^3 |delta_k|^2`, `k_mean` weighted. Does POWERfft5 bin
  `PowerDM.R` the same way (bin index = nint(k/k_fund)? which modes,
  which weights, does it include the k = 0 mode or the Nyquist plane,
  is `kharm` the mean k of the bin)? The suite's V2 check and the author's
  validation found `power()` = `PowerDM.R` to 4e-6 for n = 11 on a 72^3
  mesh: judge whether that agreement is structural or coincidental for
  small n (e.g. the Nyquist-plane treatment only matters at the mesh
  edge). Memory of `power()` for the production cube (k_c = 1.5, Box =
  1024: n = 244, 58.6M modes): count the temporaries (`index_grids` int64
  x3, `mm`, `w`, `b`, `p2` complex128 copy...) and estimate the peak;
  a login-node user calling it on one cube should not need > 1-2 GB.
- Momentum cubes: units `(1 + delta) v_i` with `v` in km/s at `a_v`
  (`vfac = 100 Box/NGRID/a_v * Nd^3/N_p`); `N_dep` vs `N_grid`; nothing in
  the reader applies a factor, confirm the docstrings say so.
- `deltaL` cubes: which keys does their sidecar carry (`n`, `Box`, `product
  = deltaL`?), does `read_modes` work with only `n` and `Box`, and does the
  layout match (same `[ix+n, iy+n, iz]`, same normalisation, k = 0 zeroed?).
- fR_halo: pairs written per `fRHaloRow` call in catalogue row order
  (`WriteFiles` in PMP2linker.f90 calls it once per written row, after the
  mass and particle-count cuts: `git grep -n fRHaloRow -- PMP2linker.f90`);
  `rows` in the sidecar = `nfRrows`; the `catalogue`/`catalogue_bin` lines
  are relative to the run directory (`CATALOGS/...`) and the reader resolves
  them against `dirname(fR_halo file)`: correct when the file is read via a
  path from elsewhere? `FRHalo.catalogue(prefer=...)` with an invalid
  `prefer`; `rows = 0` (a box without haloes).
- Byte order: every product here is `convert='little_endian'` in its
  `open`; confirm each (the build default is big-endian) and that the
  reader's dtypes are explicitly `<`.

---------------------------------------------------------------------------
## REVIEW 3: the Python package as software (API, robustness, resources)

Findings file: docs/reviews/glamio-readers/REVIEW-3-package-findings.md

Scope: every file under python/ (glamio/*.py, tests/*.py, pyproject.toml,
README.md) as a library its users will call from analysis scripts and from
halocat; not the Fortran. Run small checks with the cosemu interpreter
against the real files listed in the preamble (read-only).

Check at least:
- API consistency across readers: argument forms (`.bin`, `.meta`, stem),
  `mmap` semantics (np.memmap through `reshape` and `.T` stays a memmap;
  which methods silently load everything, e.g. `axis_means`, `power`,
  `positions`), 0-based vs 1-based indices (`LensPlanes.plane(axis, slab)`
  is 1-based while numpy users expect 0-based: is it documented loudly
  enough, and is `IndexError` for 0 right?), property names (`n_fr`,
  `fR_bg_old`, `classes`, `d_code`, `k_fund`), `__len__`/`__getitem__`
  conventions (`Schedule.__getitem__` accepts both a column name and a row
  index; `PowerDM.__getitem__` only names), `repr`s.
- Robustness: every `ValueError`/`FileNotFoundError` path; what happens with
  an empty product (`rows = 0`, `n_sel = 0`, `n = 0`), a sidecar value with
  a trailing comment (e.g. `a_v = 0.28   (velocity epoch)`), `D` exponents,
  missing optional keys (`Box` absent -> `k_fund` raises, `positions`
  raises: are the messages clear?), a `.bin` newer than its `.meta`, a path
  that is a directory, non-ASCII bytes in the 45-char header, the
  `_EPOCH_RE` on asterisks.
- Resource use on a shared login node (the user's rules: <= 16 GB, light
  jobs): `Snapshot.read()` has no default `max_particles` (a 2048^3 run is
  206 GB); `Subsample.ids()` at M = 1 allocates `N_p` uint64; `Modes.power()`
  temporaries for the production cube (n = 244); `LensPlanes.axis_means()`
  on an 805 MB memmap; `iter_pages` reads with `np.fromfile(offset=...)`
  (one 24 MiB read per page; fine?) Suggest guards or chunking where a
  careless call would exceed a few GB.
- Correctness of pure-Python logic you can test without Fortran: the
  Hermitian completion, `power()` against a brute-force double loop with
  random cubes of several n (do it), `splitmix64` against the two reference
  vectors in the tests and against a pure-Python integer implementation for
  random 64-bit inputs (do it), `select_ids` block independence, the
  `_header_dtype` itemsize and field offsets (print them), the PowerDM
  caption-to-name mapping, `resolve_pair` edge cases, `Schedule.row` on an
  absent step, `FRHalo.catalogue_paths` when the sidecar value has no
  comment.
- Package metadata: `__init__` exports vs `__all__`, the version bump and
  `tests/test_catalogue.py`'s relaxed version assertion, pyproject
  (`packages.find include = ["glamio*"]`), Python >= 3.10 features used
  (`str | None`, `match`?), numpy >= 1.24 features used, any accidental
  dependency; `from .catalogue import ...` in frhalo.py creating an import
  cycle risk.
- Code quality: dead code, duplicated helpers between product.py and
  catalogue.py (the catalogue's own `_check_layout`, `_read_catalogue_binary`
  could use `check_size`/`load_array`), inconsistent float32/float64
  returns, docstrings that disagree with the code (quote them).

---------------------------------------------------------------------------
## REVIEW 4: tests, executed in a worktree, with a mutation battery

Findings file: docs/reviews/glamio-readers/REVIEW-4-tests-findings.md

This review may build and run. Use ONLY the detached worktree the author
prepared at
/cosma8/data/dp203/dc-ruan1/scratch/review-glamio-readers-20261018/wt
(a checkout of commit 0a6ecff; tracked files only, nothing else). Do not
build or run in the main checkout. Environment for building and running:

    module purge; module load intel_comp/2024.2.0; export I_MPI_F90=ifx
    module load compiler-rt tbb compiler mpi
    export PATH=/cosma/apps/durham/dc-ruan1/micromamba/envs/cosemu/bin:$PATH
    export OMP_NUM_THREADS=4 SKIP_MODULES=1
    cd <worktree> && nice -n 19 make PMP2init PMP2start PMP2main   # ~20 s, serial (make -j breaks the module order)
    cd <worktree>/test/precise_z && nice -n 19 bash run_tests.sh   # ~40 s, 4 threads, < 100 MB
    cd <worktree>/python && GLAMIO_TEST_TMPDIR=<worktree>/python/tmp nice -n 19 python3 -I -m pytest -q tests
    cd <worktree>/python && nice -n 19 python3 -I tests/test_readers.py     # script mode

Everything you run must stay within this light tier (one process, <= 4
threads, minutes, well under 16 GB) and must be finished, with no background
processes, before you report. `run_tests.sh` deletes and recreates its
`work/` directory; `asserts.py work` can be rerun standalone on an existing
`work/` directory (seconds). The worktree's `python/glamio` is what the
worktree's `asserts.py` and `tests/` import (both insert the worktree's
`python/` at the front of sys.path), so a mutation of the worktree's
`python/glamio/*.py` is exercised by both without reinstalling anything;
the main checkout's editable install is not involved as long as you run
from the worktree with `-I`.

Scope: the diff of test/precise_z/asserts.py (the `glamio` block and its
helper `glamio_check`), test/precise_z/README.md, python/tests/test_readers.py
and the one-line change in python/tests/test_catalogue.py.

Check at least:
- Run the suite once and record the six `glamio:` verdict lines and the
  totals (the author reports 221/221, ~40 s), and run the unit tests both
  ways (15 passed / 8 passed + 7 passed).
- Discriminating power by mutation of the worktree's reader code (restore
  the file after each; `git -C <worktree> checkout -- python/glamio` restores
  everything): at least these, and report for each whether the suite
  cross-check (`asserts.py work`), the unit tests, both or neither catch it:
  lensplanes: return `planes[...]` instead of `.T` in `plane()`; drop the
  `3 * ns` size check. modes: drop the conj in `mode()`; use weight 1 for
  iz > 0 in `power()`; reshape as `(n+1, 2n+1, 2n+1)`. subsample: decode
  with `+ 32768` only (no `+ 0.5`); velocities without `dv`; `ids()` with a
  different salt. snapshot: force `bo = "<"`; `extras[98]` for the box;
  swap `aexp0` and `aexpn` in `_FIELDS`; keep the padding of the last
  page; `a_v = aexpn + 0.5*astep`; `velocity_scale` with `/ aexpn`.
  powerdm: map `Preal` to `nmodes`; use `float(epoch.group(2))` for `a`.
  schedule: read `class` from column 4. product: `check_size` with `<=`.
- Independence of the oracles: do the suite's decoders share any code or
  assumption with glamio (same reshape, same formula) such that a shared
  mistake would pass? Where they are the same formula by necessity, say so.
- Coverage gaps: multi-file `PMcrs` with the production constants (page
  1024^2, 512 records) is only exercised with `page=4, records_per_file=2`
  in the unit test; `rows = 0`/`n_sel = 0` products; `mmap=True` paths in
  the suite; `*_deltaL` cubes on real files; a `PowerDM.Z` file; the
  checkpoint `PMcrd.DAT` (the suite's work dirs have one); a big-endian vs
  little-endian snapshot pair; `Subsample` with `M > 1` on real output
  (test/extras/run_subsample_test.sh produces such files from a Fortran
  unit test: does it still run in the worktree, and could its outputs feed
  a glamio check?).
- Test hygiene: fixtures under `GLAMIO_TEST_TMPDIR`, removed afterwards;
  random fixtures seeded; pytest/script-mode parity; the `_raises` helper
  (does it mask a wrong exception type?); the assertion
  `abs(S.dx_quant - 150.0 / NQ) < 1e-10`; the SplitMix64 reference vectors
  (verify the two constants independently); the asserts block's
  `glamio_check` catching every exception as a failure (good) but also
  hiding tracebacks (is the message enough to debug?).
- At the end delete `<worktree>/test/precise_z/work`, `<worktree>/python/tmp`,
  any `__pycache__` and `.pytest_cache` under the worktree; built objects
  are fine to leave; restore any mutated file; do NOT remove the worktree
  itself (the author will).

---------------------------------------------------------------------------
## REVIEW 5: documentation, conventions and downstream use

Findings file: docs/reviews/glamio-readers/REVIEW-5-docs-downstream-findings.md

Scope: python/README.md (whole file, especially "Other products" and the
conventions list), docs/degrace_companion_implementation_plan.md (the new
section 16 entry and, for consistency, sections 3, 4, 6, 9 and the item
list in 11), test/precise_z/README.md, the docstrings of every new module,
and the consumers that will use the readers: DEGRACE-companion/check_run.py
(its own decoders of the same files), test/extras/check_*.py, and the
halocat ingest described in the plan (halocat is a separate repository,
not available here; judge the API from the plan's description of what
halocat needs: catalogue + fR_halo joins, subsample IDs, streamed
snapshots). Read-only; small Python checks allowed to verify README
examples against the real files listed in the preamble.

Check at least:
- Every README statement against the code: the readers table, the usage
  example (run it against test/precise_z/work/case11_modes/Run1 and
  case10_lens72/Run1 with the paths adapted; the `for step in
  sched.primary()["step"]` loop will not find lens planes in case11: does
  the example mislead?), the conventions bullets (lens plane orientation,
  mode normalisation and the "`iz > 0` modes count twice" statement,
  subsample decode, snapshot units and `a_v`, PowerDM window) and the
  "equals PowerDM.R bin by bin (raw CIC)" claim.
- Consistency with the plan's conventions (section 3: epochs `a`, `a_v`,
  velocity conversions; section 4: what is written when; section 9: file
  sizes) and with the status deck's definitions frame
  (docs/slides/degrace_companion_status/degrace_companion_status.tex, the
  "Definitions" frames: the deck's statement of each product's layout and
  units must match the readers' docstrings).
- The section 16 entry: is every claim in it verifiable from the diff and
  the suite (counts, "snapshot particles = subsample particles within the
  quanta", "ties the header's a_v to the subsample's")?
- Downstream: what halocat needs that the API lacks (e.g. reading a
  catalogue + fR_halo + subsample of the same step by step number from a
  run directory; iterating the schedule's primary outputs; mapping subsample
  IDs to lattice positions; a function that lists the products present for
  a step); duplicated decoders in check_run.py and test/extras that could
  now call glamio (should they, given the suite must stay an independent
  oracle?); what is still unreadable (ProjectionXY/XZ, timing.log,
  PowerDM.Z on a production run, the BDM `outputB` logs) and whether the
  plan should list it.
- Naming and units in docstrings: `n_fr`, `fR_bg_old`, `classes`, the
  1-based `slab`, `k_fund`, `d_code`; `Mpc/h` vs `h^-1 Mpc` consistency
  with the deck; the `RHO_CRIT` constant vs the finder's 2.774e11.
- Repository hygiene of the commit: files that should not be tracked, the
  .gitignore entries for `python/tmp`, `.pytest_cache`, `__pycache__`;
  the README's install line still valid after the version bump; the
  precise_z README row.
