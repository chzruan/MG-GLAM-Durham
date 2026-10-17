# REVIEW 4: test suite changes, executed in a worktree

Branch `bdm-binary-catalogue` (commit 167fb1c, `cz..bdm-binary-catalogue`).
Scope: `test/precise_z/{asserts.py,run_tests.sh,README.md}` and the harness
pieces they use. Built and run only in the detached worktree
`/cosma8/data/dp203/dc-ruan1/scratch/review-bdm-binary-catalogue-20261018/wt`
(ifx 2024.2.0, serial `make`, `OMP_NUM_THREADS=4`, `nice -n 19`; build 21 s /
110 MB, suite 41 s / 67 MB). Reviewer: Claude (Fable 5.1), 2026-10-18.

## 1. Summary

The suite reproduces the author's figure: **194/194 checks pass** in a clean
build of 167fb1c, and the three new verdicts (case12 binary catalogue, case13,
case4g) are genuine: 28 of 33 deliberate corruptions of the produced files are
detected (see Appendix A), the two undetected data mutations are below the
text's printed precision by construction, and the TEXT_TOL tolerances are
derived correctly from the Fortran format with ~20 % headroom. However, a
probe of a configuration the suite does not cover found a production defect:
with no `BDM.config` present, the finder writes its default file with the key
truncated to `Catalogue = both` (format 30's `T20` overwrites the 15-character
key) and the run dies at the **second** catalogue output with `unrecognized
parameter: catalogue` (F1, Blocker; Fortran, i.e. Review 1 territory, but
confirmed here by execution). The suite also has no `text`-mode case, no guard
against a vacuous pass with zero haloes, and aborts with a traceback (losing
~19 later verdicts) whenever glamio rejects a file.
**Recommendation: not ready.** Fix F1 and add the regression test F2; F3-F5 are
cheap and recommended before merge.

## 2. Findings

| ID | Severity | File:line | Description | Evidence | Suggested fix |
|----|----------|-----------|-------------|----------|---------------|
| F1 | **Blocker** | `PMP2linker.f90:372` (`30 format(10x,a,T20,' = ',a,T40,a)`), used at `:357` (default `BDM.config` writer) and `:369` (parameter echo) | `10x` puts `CatalogueFormat` in columns 11-25; `T20` then writes `' = '` over columns 20-22 and the value from column 23, so the default file contains `Catalogue = both`. On the next `ReadParameters` (every later catalogue step) `select case(trim(key))` hits `case default` (`:338`) and the run stops. Every run without a `BDM.config` and with >= 2 catalogue outputs dies at the second one. New in this commit: on `cz` the longest key written through these formats is `MaxMemory` (9 chars), so `T20` was safe. | Probe `x_default` (LCDM, `Find BDM halos = 1`, `#outputs = -2` at z = 1.1, 0.512, no `BDM.config`): step 84 writes `Run1/BDM.config` whose line 10 reads `          Catalogue = both             ! text (CatshortV.*.DAT), ...`; `main.log:963`: ` BDM configuration error, line 10: unrecognized parameter: catalogue`, exit 1, no step-98 catalogue. Also visible in every suite log: `main.log:927  Catalogue = binary` (case13), `:1151  Catalogue = both` (case12). | Drop the tab or widen it, e.g. `30 format(10x,a,' = ',a,T44,a)` (or `T27`/`T44` for all three formats so the columns stay aligned). Add F2 as the regression test. `gen_inputs.py:219` writes the key itself and is unaffected; runs that let the code write the default are the ones hit. |
| F2 | Major | `test/precise_z/run_tests.sh` (no case), `asserts.py` (no check) | The path "no `BDM.config` -> code writes its default -> re-reads it at the next catalogue step" is not exercised, although it is the default for every run directory without the file. A green 194/194 therefore says nothing about the default configuration, which is broken (F1). | `x_default` probe above; the suite only ever supplies a `BDM.config` (`test/BDM.config` without the key for case12, with the key for case13/case4g). | Add `case14_catdefault`: LCDM, `Find BDM halos = 1`, `#outputs = -2 " 1.1 0.512"`, `rm -f Run1/BDM.config`; assert the run completes (`Step = 100` or the z = 0.512 snapshot), the written `BDM.config` re-parses (`grep -c '^ *CatalogueFormat *= *both' == 1`), and `CatshortV.0084/0098.0001.{DAT,bin,meta}` all exist. Runs in ~1.4 s. |
| F3 | Minor | `asserts.py:596` (`cat = read_catalogue(stem)`), `:627`, and the direct `cat.meta["rows"]`/`["text_catalogue"]`/`["da"]` accesses at `:599-603`, `:629` | glamio raises `ValueError`/`KeyError` on a truncated or extended `.bin`, a wrong `rows`, a reordered `columns` line or a `>f4` dtype; asserts.py does not catch it, so the script dies with a traceback: no `N/M checks passed` line, and every later check (case13, lens36 block means, L3, cases 9/1/2/3/6/7/8/5b and all 14 error-case checks) is never evaluated. Exit code is 1, so CI would still fail, but the report is lost. | Mutations m07, m08, m13, m15, m16, m19 in Appendix A: `rc=1`, `(no totals line)`, `traceback=True`. | Wrap the reads: `try: cat = read_catalogue(stem) except Exception as exc: rep.check(False, case, f"glamio refused {stem}: {exc}"); continue`; use `cat.meta.get(...)` with a FAIL on `None`. |
| F4 | Minor | `asserts.py:599-611` (case12), `:629-631` (case13) | The sidecar's layout descriptors are optional for both glamio (`_check_layout` only tests keys that are present) and the test. A `.meta` without `columns`, or without `dtypes` and `record_bytes`, still passes 194/194, so the lines the README advertises as the reader's self-description are unverified by the suite. | m20 (`columns` removed) and m21 (`dtypes` + `record_bytes` removed): `194/194 checks passed`. | In case12 and case13 assert `cat.meta["columns"].split() == list(CATALOGUE_DTYPE.names)`, `cat.meta["dtypes"].split() == [...]` and `int(cat.meta["record_bytes"]) == 96` explicitly. |
| F5 | Minor | `asserts.py:598` (`nb_total`), `:628` (`n13`), `:615`, `:633` | The row totals are printed but never asserted positive. All column checks are `np.all` over the rows, so a catalogue with zero haloes passes vacuously. The margin is thin: the f(R) run has 1 halo at z = 1.1 and 5 at z = 0.512 (6 in total), the LCDM run 5 at z = 0.512 and **0** at z = 1.1 (`x_default` step 84: header-only `.DAT`, 0-byte `.bin`, `rows = 0`). A different seed table, a higher `MassMin` or dropping the z = 0.512 output turns the suite green with nothing compared. | `PASS case12_frhalo: ... (rows in total: 6)`, `PASS case13_catbin: ... (rows: 5)`; m33 shows a consistent empty catalogue is only caught because the text still has 5 rows. | `rep.check(nb_total >= 1, "case12_frhalo", ...)` and `rep.check(n13 >= 1, "case13_catbin", ...)`. Optionally add an explicit zero-halo check on an LCDM z = 1.1 output (0-byte `.bin`, `rows = 0`, 8-line `.DAT`): glamio reads it (verified), and `text_rows` handles it (with a numpy UserWarning, F11). |
| F6 | Minor | `run_tests.sh:247-254` (only `binary` and the implicit `both` are run) | No case runs `CatalogueFormat = text`, so the suite never shows that the text catalogue is unaffected by the binary machinery, nor that `.bin` is identical between `binary` and `both`. Both properties hold and are cheap to assert: the runs are bitwise reproducible at 4 threads (1.2-1.4 s each). The pre-change byte identity of the `.DAT` (Review 1's regression question) cannot be shown inside the suite without a `cz` build; the probe at least shows text-mode output == both-mode output. | Probes: `x_text` (`text`) writes `CatshortV.0098.0001.DAT` only, no `.bin`/`.meta`, no `fort.*`; `cmp x_text/.DAT x_both/.DAT` -> identical; `cmp x_both/.bin case13/.bin` -> identical; `.meta` differ only in `text_catalogue`. Mixed-case `Both` accepted. | Add an LCDM `both` twin of case13 plus a `text` twin; assert `.DAT(text) == .DAT(both)` and `.bin(both) == .bin(binary)` byte for byte (`filecmp.cmp(..., shallow=False)`), and that the `text` run leaves no `.bin`/`.meta`. Use `Both`/`TEXT` in one of them to cover the lowercase conversion. |
| F7 | Minor | `asserts.py:72` (`"Run1/CATALOGS/Catshort*"`), `run_tests.sh:270-277`, README `:62` | (a) The forbidden glob does not match hidden staged temporaries `.Catshort*` nor `outputB.*`; today nothing is left because `ReadParameters` fails before `BeginCataloguePublication` and before `outputB` is opened, but the check would not notice a regression that stages first. (b) README says the run "stops at the first catalogue step"; the test only requires the message somewhere in `main.log`, not that step 98 was reached. No false pass is possible today (the message is unique to `ConfigurationError`), but an earlier validation (e.g. at start-up) would pass silently. | case4g `main.log`: `:956 Step = 98`, `:963 BDM configuration error, line 10: CatalogueFormat must be text, binary or both`, then `1`; `CATALOGS/` empty; `PMcrd.0098.DAT`/`PMcrs0.0098.DAT` written before BDM. | Assert `os.listdir(Run1/CATALOGS) == []` (covers hidden files and `outputB`), and `re.search(r"Step =\s+98 ", text)` if "at the first catalogue step" is the intended contract (see Q1). |
| F8 | Nit | `PMP2linker.f90:369` | Parameter echo in every BDM log reads `Catalogue = both/binary/text`. Same root cause as F1; cosmetic on its own. | `work/*/Run1/main.log:927`. | Fixed by F1. |
| F9 | Nit | `test/BDM.config` (last line `dLogP = ...` has no newline), `run_tests.sh:251,274` | The `{ cat ...; echo; echo '...'; }` construction is robust: if the file later gains a final newline the extra blank line is skipped by the parser (`PMP2linker.f90:291 if (len_trim(line) == 0) cycle`), and the `! ...` comment on the appended line is stripped at `:289-290`. case12 copies the unterminated file verbatim and ifx reads the last record fine. Still, the file should end with a newline (POSIX text file; `tail`, `wc -l`, diff tools). | `cat -A test/BDM.config` ends `...for profiles` without `$`. | Add the trailing newline to `test/BDM.config`; keep the `echo`. |
| F10 | Nit | `asserts.py:630`, README `:61` ("finite positive masses") | case13 checks only `Mbound`; `Mtot` and the other 21 columns are unchecked there (m32: a flipped exponent byte in `Rvir` passes). Acceptable since the `.bin` is byte-identical to the `both` run's (F6) and case12 checks all 24 columns, but the README overstates. | m32 -> 194/194. | Say "finite positive `Mbound`" in the README, or check `Mtot` too. |
| F11 | Nit | `asserts.py:581` (`np.loadtxt(cat, skiprows=8, ndmin=2)`) | On a header-only catalogue numpy emits `UserWarning: loadtxt: input contained no data` (glamio suppresses the same warning; asserts.py does not). Harmless noise; relevant once F5's zero-halo check exists. | Verified on `x_default/.../CatshortV.0084.0001.DAT`: shape `(0, 1)`, one warning; `text_rows` correctly returns `zeros((0, 24))`. | Wrap in `warnings.catch_warnings()` as glamio does. |
| F12 | Nit | `asserts.py:543-566` (fR_halo block) | The `fR_halo` sidecar's new `catalogue_bin` line (`PMP2mod_extras.f90:648-649`) is not checked. | Diff of `PMP2mod_extras.f90`; no `catalogue_bin` in asserts.py. | One line: `os.path.exists(os.path.join(rundir, meta["catalogue_bin"].split()[0]))` and rows equal. |

## 3. Verified OK

Execution (worktree `wt`, commit 167fb1c, clean build):
- `nice -n 19 make PMP2init PMP2start PMP2main` serial: OK, 21.3 s, 110 MB.
- `nice -n 19 bash run_tests.sh`: `194/194 checks passed`, exit 0, 41.1 s wall,
  67 MB max RSS, no `FAIL` line. Verdicts observed:
  - `PASS case12_frhalo: binary catalogue .bin + .meta beside every text catalogue; no staged temporaries`
  - `PASS case12_frhalo: meta rows == .bin size/96 == text rows; meta a == header A, a_v = a - da/2; text_catalogue named`
  - `PASS case12_frhalo: all 24 binary columns equal the text at its printed precision; row = 1..N, distinct = 0 (rows in total: 6)`
  - `PASS case13_catbin: CatalogueFormat = binary: .bin + .meta only, readable, text_catalogue = none (rows: 5)`
  - `PASS case4g_badcatfmt: error message "CatalogueFormat must be text, binary or both" reported`
  - `PASS case4g_badcatfmt: no Run1/CATALOGS/Catshort* produced`
  - plus the unchanged schedule/OutputSchedule checks for case12 (snapshots 84, 98, 100) and case13 (98, 100).
- `asserts.py work` standalone reruns in 0.6 s; after all mutations were
  restored the final run was again 194/194.
- asserts.py imports the worktree's glamio (`sys.path.insert` verified via
  `glamio.__file__`), not the editable install of the main checkout;
  numpy 2.3.5, Python 3.14.2.

Discriminating power (Appendix A):
- Record layout and order: swapped records (m06), swapped `Mbound`/`Mtot`
  fields (m10), `distinct = 1` (m11), `row = 0` (m12) are all caught by the
  column check.
- Endianness: a byte-swapped field (m09) is caught; combined with the
  byte-level equality through `<f4` this proves the file is little-endian
  although the whole code is compiled with `-convert big_endian`
  (`convert='little_endian'` on the stream `open` is honoured by ifx).
- Size/rows consistency: `rows` off by one (m13), truncation (m07),
  extension (m08) are refused by glamio; a *consistent* lie (`rows = 4` with
  a 4-record file, m14) and a consistent empty pair (m33) are caught through
  the text row count. Missing `.bin` or `.meta` (m25, m26) caught.
- Epochs: `a_v = a` (m17), `a_v` off by 1e-6 (m18), `da` doubled (m23)
  caught; `a` agrees with the text header `A` to 1e-5 (f8.5 rounding is
  5e-6).
- `text_catalogue`: wrong in `both` (m22) and in `binary` (m28) caught.
- Staged temporaries `.Catshort*` detected in both cases (m24, m29); a stray
  `.DAT` in binary mode detected (m27); negative `Mbound` (m30) and reordered
  rows (m31) detected in case13.
- Not detectable by design: 1 ulp in a position (m02) and up to 255 ulp in a
  g12.4 column (m04) lie below the text's printed precision; the text is the
  only reference, so this is inherent, not a weakness.

TEXT_TOL (`asserts.py:576-578`), derived from
`(3f11.4,3x,3f10.2,1p,2g12.4,g12.5,26g12.4)`:
- `f11.4` positions: rounding <= 5e-5 absolute; tol 6e-5 atol. Correct for
  any Box (the text is a rounding of the same real*4 value).
- `f10.2` velocities: <= 5e-3; tol 6e-3 atol. Overflow (asterisks) needs
  |v| >= 1e7 km/s, impossible.
- `1p,g12.4` masses: always E form (>= 1e4), scale factor 1 gives 5
  significant digits -> relative error <= 5e-5; tol rtol 6e-5. (A mass < 1e4
  would fall into F form with 4 digits and fail spuriously; impossible since
  rows need >= 20 particles.)
- `g12.5` Rvir: F form for 0.1 <= R < 1e5 with 5 significant digits (the
  scale factor is ignored in F form, as the output `2802.3` for 2802.3098
  confirms) -> <= 5e-5 relative; tol 6e-5.
- `g12.4` for the 15 remaining columns: F form 4 significant digits ->
  <= 5e-4 relative (worst case mantissa 1.000, e.g. Vcirc `1119.` for
  1118.87), E form (|x| < 0.1 or >= 1e4) 5 digits -> <= 5e-5; tol rtol 6e-4
  plus atol 1e-12 for exact zeros (`0.000`). Negative values (`-0.5741`) and
  F/E switches (`6.8017E-02`) behave identically. The two integer columns
  (`row`, `distinct`) are printed with I editing and compared exactly within
  the same tolerance. No legitimate value can fail; a wrong column, unit
  (1e3 on RadRMS) or epoch factor (1/a) differs by O(1) and is caught.

case13 (`text` off):
- `.DAT` absence is checked strictly (`not os.path.exists`), `text_catalogue`
  must equal `none` exactly; a run that silently wrote the `.DAT` fails (m27).
- No `fort.12`, `fort.32` or `fort.33` appears in any `Run1/` (the header
  loop guard `If(kf==1.and..not.CatText) cycle` at `PMP2linker.f90:2056` is
  effective; `close(12)` on a never-opened unit is harmless). The only
  `fort.*` litter is `fort.17` in `work/ic/Run1` and `work/case6_tailinsert/
  Run1`, written by `PMP2start.exe` and pre-existing (not this commit).
- Publication of two files: both exist with the right sizes (480 B = 5 x 96,
  1660 B meta); the rename order cannot be tested post hoc.

case4g:
- The run really executes 98 steps, writes the step-98 snapshot, then dies in
  `ReadParameters` (line 10 of the config) before `outputB` or any catalogue
  file is created: `CATALOGS/` is empty, including hidden files. The message
  reaches `main.log` via `2>&1` (`error_unit`). A run that died earlier for
  another reason cannot false-pass: the message text is unique to this check.

Harness hygiene:
- The Init.dat `sed` anchors (`^Find BDM halos            =         0`,
  `^MG_flag                   =         0`) each match exactly one line of
  `Init.base.dat` (grep -c = 1); a failed match would turn BDM off and the
  catalogue checks would FAIL, not pass.
- `set_outputs` awk edit and `copy_ics` unchanged; case13/case4g follow the
  existing case pattern (`run_init ... || echo`, `copy_ics && run_main`).
- The `echo` fix for the unterminated `test/BDM.config` is robust (F9).
- `test/BDM.config` has no `CatalogueFormat` key, so case12 covers the
  "old config file read by the new binary" path -> `both`.
- README rows 60-62 describe what asserts.py checks (except the "masses"
  wording, F10); the README quotes no check total; the commit message's
  194/194 is reproduced.
- Suite runtime and memory fit the login-node light tier.

Probes outside the suite (worktree only, each ~1.3 s):
- `text`: `.DAT` + `outputB` only; byte-identical `.DAT` to the `both` run.
- `Both` (mixed case): accepted; `.bin` byte-identical to case13's `binary`
  `.bin`; `.meta` differs only in `text_catalogue`.
- Runs are bitwise reproducible with 4 OpenMP threads, so byte-identity
  regression checks are feasible.
- Zero-halo output (`x_default` step 84): 0-byte `.bin`, `rows = 0`, 8-line
  `.DAT`; `glamio.read_catalogue` returns 0 rows with `a`/`a_v` intact.

## 4. Open questions for the author

1. case4g spends 98 steps before the configuration error (the snapshot is
   written first). Should `BDM.config` be validated at `PMP2main` start-up
   when `Find BDM halos = 1`, and should the test then assert the early stop
   instead of "at the first catalogue step"?
2. After fixing F1, is the default-file writer worth keeping at all, or
   should a missing `BDM.config` be an error when BDM is on? Either way F2
   should pin the chosen behaviour.
3. Is a 0-byte `.bin` with `rows = 0` the intended on-disk form of an empty
   catalogue (what the `Nmaxima == 0` path produces)? If so, python/README
   should say so and F5's zero-halo check should assert it.
4. The column-by-column check rests on 6 haloes from an f(R) run. Would you
   accept a dedicated LCDM `both` twin (F6) as the primary layout test, with
   case12 keeping the fR_halo cross-check?

## Appendix A: mutation battery (asserts.py rerun after each, files restored)

Baseline: `rc=0 194/194`. "Exception" = glamio raised, asserts.py aborted
with a traceback, no totals line, later checks skipped (F3).

| Mutation | Result | Detected by |
|---|---|---|
| m01 `.bin`: flip exponent/sign byte of `x`, record 0 | FAIL (193/194) | case12 "all 24 binary columns equal the text" |
| m02 `.bin`: flip bit 0 of `x`, record 0 (1 ulp = 7.6e-6) | **pass** 194/194 | below printed precision (5e-5) - inherent |
| m03 `.bin`: flip bit 7 of low byte of `x` (128 ulp = 9.8e-4) | FAIL | columns |
| m04 `.bin`: flip whole low byte of `Cvir` (<= 255 ulp ~ 3e-5 rel) | **pass** | below g12.4 precision (5e-4) - inherent |
| m05 `.bin`: flip bit 7 of byte 1 of `Cvir` (2^15 ulp ~ 3.9e-3) | FAIL | columns |
| m06 `.bin`: swap records 0 and 1 | FAIL | columns (+ row index) |
| m07 `.bin`: truncate by one record | Exception | glamio: size != rows x 96 |
| m08 `.bin`: append 96 zero bytes | Exception | glamio: size != rows x 96 |
| m09 `.bin`: byte-swap field `x` of record 2 | FAIL | columns |
| m10 `.bin`: swap `Mbound` and `Mtot`, record 1 | FAIL | columns |
| m11 `.bin`: `distinct = 1`, record 3 | FAIL | columns (`distinct == 0`) |
| m12 `.bin`: `row = 0`, record 0 | FAIL | columns (`row == 1..N`) |
| m13 `.meta`: `rows = 6` (bin untouched) | Exception | glamio: size != rows x 96 |
| m14 `.meta` `rows = 4` + `.bin` truncated to 4 records | FAIL x2 | meta (`!= text rows`) and columns |
| m15 `.meta`: `columns` with `x`/`y` exchanged | Exception | glamio `_check_layout` |
| m16 `.meta`: `columns` with `Mbound`/`Mtot` exchanged | Exception | glamio `_check_layout` |
| m17 `.meta`: `a_v = a` | FAIL | meta (`a_v = a - da/2`) |
| m18 `.meta`: `a_v` off by 1e-6 | FAIL | meta |
| m19 `.meta`: `dtypes` with `>f4` for `x` | Exception | glamio `_check_layout` |
| m20 `.meta`: `columns` line removed | **pass** 194/194 | nothing (F4) |
| m21 `.meta`: `dtypes` and `record_bytes` removed | **pass** 194/194 | nothing (F4) |
| m22 `.meta`: `text_catalogue = none` while `.DAT` exists | FAIL | meta (`endswith`) |
| m23 `.meta`: `da` doubled | FAIL | meta |
| m24 dir: fake `.CatshortV.0098.0001.bin.tmp.1.2.3` | FAIL | files ("no staged temporaries") |
| m25 dir: delete step-84 `.meta` | FAIL | files |
| m26 dir: delete step-84 `.bin` | FAIL | files |
| m27 case13: fake empty `CatshortV.0098.0001.DAT` | FAIL | case13 (`.DAT` must not exist) |
| m28 case13: `text_catalogue = CATALOGS/...DAT` | FAIL | case13 (`== none`) |
| m29 case13: fake `.CatshortV.0098.0001.meta.tmp.9.9.9` | FAIL | case13 glob |
| m30 case13: `Mbound = -1`, record 0 | FAIL | case13 (`Mbound > 0`) |
| m31 case13: swap records 0 and 4 | FAIL | case13 (`row == 1..N`) |
| m32 case13: flip exponent byte of `Rvir`, record 1 | **pass** 194/194 | not checked in case13 (F10) |
| m33 `.bin` emptied + `.meta` `rows = 0` (consistent) | FAIL x2 | meta and columns, via the text row count |

## Appendix B: cleanup

Deleted after the review: `wt/test/precise_z/work` (incl. the probe
directories `x_text`, `x_both`, `x_default`), `wt/python/glamio/__pycache__`,
and the reviewer's `wt/build.log`, `wt/suite_run1.log` (copies kept in the
session scratchpad). Built objects and executables left in `wt`; the worktree
itself left for the author. No tracked file in `wt` or the main checkout was
modified; the only file written outside `wt` is this one.
