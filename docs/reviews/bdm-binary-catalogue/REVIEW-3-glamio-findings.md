# REVIEW 3: the glamio Python package — findings

Branch `bdm-binary-catalogue`, commit 167fb1c, range `cz..bdm-binary-catalogue`.
Scope: `python/pyproject.toml`, `python/glamio/{__init__,meta,catalogue}.py`,
`python/tests/test_catalogue.py`, `python/README.md`, the `.gitignore` additions,
cross-checked against `PMP2linker.f90` (`WriteFiles`, `WriteCatalogueMeta`, the text
header loop in `SetParameters`) and the LCDM pilot text catalogue
`DEGRACE-companion/runs/Node_002_LCDM_L1024Np1024Ng2048/Run1/CATALOGS/CatshortV.0181.0001.DAT`.
Read-only review; the only executions were the unit tests (pytest and script mode),
one read of the pilot text catalogue through glamio (1.8 s, 410 MB RSS) and small
numpy probes on fixtures in the session scratchpad (removed afterwards, `__pycache__`
deleted). Reviewed 2026-10-18.

## 1. Summary

The package is small, numpy-only, typed and documented, and it reads what the Fortran
writes: the 24-field `CATALOGUE_DTYPE` matches the binary `write(CatUnitBin)` item list
field by field (all items are 4-byte real*4 / integer*4, itemsize 96, explicit `<`), the
sidecar `columns`/`dtypes` lines written by `WriteCatalogueMeta` equal the dtype names and
kinds, and the text reader parsed the real 777,651-row pilot catalogue with every header
key and every column agreeing with an independent parse. The six unit tests pass under
pytest and as a script. No Blocker or Major was found. The Minor items concern
robustness and contract consistency: no self-consistency check on decoded records (a
byte-swapped `.bin` is returned as valid data), inconsistent meta key names between the
text and binary forms (`step_da` vs `da`, `format`), a float32 downcast in the text
reader that can exceed the 6e-5 position tolerance used by `test/precise_z/asserts.py`
on boxes larger than ~256 Mpc/h (safe in the 150 Mpc/h test box), a text fixture that
does not reproduce Fortran G editing, and the version string kept in three places.
Merge recommendation: **ready**. None of the findings produces wrong results from a
correctly written catalogue; R3-01, R3-02 and R3-03 should follow in a small commit before
the binary catalogue is consumed in production or the case12 check is reused on a
production box.

## 2. Findings

Counts: Blocker 0, Major 0, Minor 8, Nit 10.

| ID | Severity | File:line | Description | Evidence | Suggested fix |
|---|---|---|---|---|---|
| R3-01 | Minor | `python/glamio/catalogue.py:234-255` | The binary reader validates the sidecar and the file size but never the decoded values. A byte-swapped, field-shifted or misaligned `.bin` whose sidecar says `f4` (the Fortran line carries no byte order; `little-endian` is only prose in `format`) is returned as a valid `Catalogue`. The `row` column is 1..N by construction (`PMP2linker.f90:695,704`) and `distinct` is always `0_4`, so a self-check is free. | Probe: `arr.byteswap().tofile(stem+".bin")` with the Fortran-style `.meta` → `read_catalogue(stem)` succeeds, `vx = [7.18e-41 2.61e-41 5.39e-41]`, no error. | After `np.fromfile`: `if len(data) and not (data["row"][0] == 1 and data["row"][-1] == len(data)): raise ValueError(f"{bin_path}: row column is not 1..N (byte order or record misalignment)")`; optionally a `check: bool = True` kwarg that also asserts `(data["distinct"] == 0).all()`. Add a byte-swapped-file test. |
| R3-02 | Minor | `python/glamio/catalogue.py:277,286-294` vs `PMP2linker.f90` `WriteCatalogueMeta` | The two readers return different meta keys for the same quantities, so code written against one form breaks on the other. Text: `step_da`, `format = text`, `header` = HEADER + ` [BDM finder v3]`, no `rows`/`realization`/`Box`. Binary (Fortran): `da`, `format = raw little-endian records, 96 bytes per halo: ...` (prose; the unit-test fixture writes `format = binary`, which the Fortran never does), `header` = raw HEADER plus a separate `finder` key, `rows`, `realization`, `Box`. | `test/precise_z/asserts.py` (diff line 61) reads `cat.meta["da"]`; on a text `Catalogue` that is a `KeyError`. `test_catalogue.py:97` writes `"format": "binary"`; `WriteCatalogueMeta` writes the prose sentence; `meta["format"] == "binary"` is therefore never true for a real file. | In `_parse_text_header` also set `da` (keep `step_da` as alias), `rows = <n>` (after loading), and `step`/`realization` from the file name (`Catshort?.<step>.<box>`); strip the ` [BDM finder v3]` suffix into `finder`. For the Fortran (Review 2): write `format = binary` and move the prose to `layout = ...`; make the fixture copy the Fortran text verbatim. |
| R3-03 | Minor | `python/glamio/catalogue.py:323-328` (cast to `<f4`) together with `test/precise_z/asserts.py` `TEXT_TOL` positions atol 6e-5 | The text reader downcasts the printed decimals to float32. The writer's float32 `b` is printed with 4 decimals (`|t-b| <= 5e-5`), then cast back (`|c-t| <= ulp/2`). For 256 <= x < 1024 Mpc/h the result lands 1-2 float32 ulps = 6.1035e-5 from the writer's value, which is *above* the 6e-5 atol that case12 uses to compare text (via glamio) with binary. It passes today only because the test box is 150 Mpc/h (max error 5.34e-5). | Measured over 60k random float32 positions per band: max\|cast-writer\| = 5.34e-5 for x<150; 4.58e-5 for 150-256; **6.104e-5 for 256-512 (8.5 % of values > 6e-5)**; **6.104e-5 for 512-1024 (39 % > 6e-5)**; 0 for 1024-1800 (there the 4 decimals pin the float32 exactly). | Either compare the text as float64 (`np.loadtxt` output, error <= 5e-5) against the binary in asserts.py, or widen the position atol to >= 1.3e-4 (half printed unit + one ulp at 2048), or let `read_catalogue_text(..., float64=True)` return the raw decimals. Cross-reference Review 4 (TEXT_TOL). |
| R3-04 | Minor | `python/tests/test_catalogue.py:114-119,122-129` | `_row_text` is documented as the "Python equivalent of the Fortran format" but prints every `g12.4`/`g12.5` column in E form (`7.7100E+02`), whereas Fortran G editing prints F form for 0.1 <= \|x\| < 1e4 (`771.00      `, `191.1       `, `-0.3133     `, `88.00`), which is what every Rvir/Vrms/Vcirc/Cvir/Npart/b_a/c_a value and most Xoff/lambda values look like. The text-reader unit test therefore never parses a real Fortran row; `TEXT_ROWS[0:2]` are the pilot's first two rows, re-formatted. | Pilot row 1 (verbatim): `   394.8253     5.1722     0.0920      -361.42   -355.91      8.70  7.3234E+12  1.2649E+13  771.00       191.1       303.5               1   22.96       88.00               0  4.0263E-02 -0.3133      2.6562E-02   275.4      0.8410      0.8140     -9.9848E-02  0.9850      0.1407    `; fixture output for the same row: `7.3234E+12  1.2649E+13  7.7100E+02  1.9110E+02 ...`. (The real F form does parse — verified on the pilot — so this is test fidelity, not a reader bug.) | Embed the pilot's first two rows (and one row with an E-form `Xoff` and a negative `virial_ratio`, as above) as literal strings in the fixture; keep `_row_text` only if it implements the G editing F/E switch, otherwise rename its docstring. |
| R3-05 | Minor | `python/glamio/catalogue.py:208-231`; `python/README.md:69-72` | `_check_layout` is a no-op when the sidecar has none of `columns`/`dtypes`/`record_bytes`, so a writer that drops the keys is not noticed, contrary to the README claim that the reader "refuses a sidecar whose columns/dtypes differ ... so a format change is noticed rather than misread". `record_bytes = 96.0` and `float32` spellings pass (harmless). | Probe: `_check_layout({"rows": "3"}, "m")` returns without error. | Require `columns` and `dtypes` when `product` starts with `Catshort` (or always for files produced by GLAM); raise `ValueError(f"{meta_path}: no 'columns'/'dtypes' line")`, or at least `warnings.warn`. Adjust the README sentence accordingly. |
| R3-06 | Minor | `python/tests/test_catalogue.py` | Failure modes not exercised: byte order (a `>f4` dtypes line — the check does catch it, verified — or a byte-swapped file, not caught, see R3-01); `columns` with the same names in a different order (caught, verified, untested); a sidecar with the Fortran's actual content (prose `format`, `velocity_scale` containing a second `=`, `a_v` with trailing comment, `rows =` with i10 padding, `es16.8` numbers); a `CatshortW/M/A` file name (the `_PRODUCT_RE` label); `pathlib.Path` input; a `.bin` with extra whole records when `rows` is absent (passes by design — document). | Only `Mbound`→`Mvir` rename, `i4`→`f4` swap, `record_bytes = 100`, two truncations and a 23-column row are tested. | Add one test per item; a Fortran-exact `.meta` string (transcribed from `WriteCatalogueMeta`) is the most valuable, because it is the real contract. |
| R3-07 | Minor | `python/pyproject.toml:7`, `python/glamio/__init__.py:6`, `python/tests/test_catalogue.py:348` | The version string lives in three places (`version = "0.1.0"`, `__version__ = "0.1.0"`, `assert glamio.__version__ == "0.1.0"`), so every bump touches three files and an editable install can report a stale dist version. | `importlib.metadata.version("glamio")` = 0.1.0 today; it would diverge after a bump of `__init__.py` alone. | `dynamic = ["version"]` in `[project]`, `[tool.setuptools.dynamic] version = {attr = "glamio.__version__"}`; in the test compare `glamio.__version__` with `importlib.metadata.version("glamio")` instead of a literal. |
| R3-08 | Minor | `python/glamio/catalogue.py:254` | `np.fromfile` loads the whole `.bin`. Fine for the pilot (75 MB) but the L1800 Np3200 catalogues will be ~1 GB per output and users often want one or two columns. | 96 B/row; pilot 777,651 rows = 75 MB; an L1800 run with ~10 M haloes = 0.96 GB resident per catalogue. | Add `mmap: bool = False`; when true use `np.memmap(bin_path, dtype=CATALOGUE_DTYPE, mode="r")` (structured field access stays lazy). Document that `__getitem__` then returns a memmap view. |
| R3-09 | Nit | `python/glamio/catalogue.py:65-90` | `COLUMN_UNITS["Xoff"] = "centre offset"` omits the unit; the code computes `Xoff(ip)=real(sqrt(sum(centre**2))/aperture)` with `Rvir(ip)=real(aperture)` (`PMP2linker.f90:1133,1214`), i.e. in units of Rvir, as the Fortran `units` line says. `COLUMN_UNITS` is also unused, not in `__all__`, and never compared with the sidecar `units` line. | `git grep COLUMN_UNITS -- '*.py'` → only its definition. | Set `"Xoff": "centre-of-mass offset / Rvir"`; export it and (optionally) check `meta["units"].split()` has 24 tokens in `_check_layout`, or drop the dict. |
| R3-10 | Nit | `python/glamio/catalogue.py:119-120` | `__repr__` evaluates `self.a`, which raises `ValueError` for a non-numeric `a` value; a `repr` should never raise. | `repr(Catalogue(arr, {"a": "unknown"}, "p"))` → `ValueError` (probe). | Use `self.meta.get("a")` in the repr. |
| R3-11 | Nit | `python/glamio/catalogue.py:319-321` | With `ndmin=2` the branch `raw.ndim != 2` is unreachable, and `ncol = raw.ndim` would report a dimension count as a column count. | `np.loadtxt(..., ndmin=2)` always returns a 2-D array (header-only → shape `(0, 1)`). | Reduce to `if raw.shape[1] != _NCOL: raise ValueError(...)`. |
| R3-12 | Nit | `python/glamio/catalogue.py:62` | Module-level `assert CATALOGUE_DTYPE.itemsize == RECORD_BYTES` is stripped under `python -O`. | — | `if CATALOGUE_DTYPE.itemsize != RECORD_BYTES: raise RuntimeError(...)`. |
| R3-13 | Nit | `python/glamio/catalogue.py:193` | Extension match is case-sensitive to exactly `.DAT`/`.dat`; `foo.Dat` is treated as a stem and fails with `FileNotFoundError: ... foo.Dat.bin`. | `"a.Dat".endswith((".DAT", ".dat"))` → `False`. | `p.lower().endswith(".dat")`. |
| R3-14 | Nit | `python/glamio/catalogue.py:200-205` | A stem with `.meta` and `.DAT` present but `.bin` absent (e.g. `.bin` files deleted to save space, `.meta` kept) raises `FileNotFoundError` instead of falling back to the text. Documented, but surprising. | Probe: stem with `.meta` + `.DAT` only → `FileNotFoundError ... .bin`. | Prefer binary only when `<stem>.bin` exists; otherwise fall back to `.DAT` and mention the orphan `.meta` in a warning. |
| R3-15 | Nit | `python/glamio/meta.py:49-53`; `catalogue.py:169-177` | `meta_float` accepts `nan`/`inf` (`float("nan")`), so `a = nan` yields NaN velocities from `velocities_at_av()` without an error; only `a_v == 0` is rejected. | `meta_float({"k": "nan"}, "k")` → `nan` (probe). | Reject non-finite `a`/`a_v` in `velocities_at_av()` (`math.isfinite`). |
| R3-16 | Nit | `python/glamio/catalogue.py:217-228` | The `dtypes` comparison is byte-order-naive: `np.dtype("f4")` is *native*, so the Fortran's `f4` line equals `<f4` only on little-endian hosts; on a big-endian host the check would fail spuriously (and `>f4` would pass). Not a practical problem on COSMA. | `np.dtype(">f4") != np.dtype("<f4")` → True; `np.dtype("f4") == np.dtype("<f4")` → True on x86. | Compare `(d.kind, d.itemsize)` pairs, or `np.dtype(c).newbyteorder("<")`. |
| R3-17 | Nit | `python/README.md:63,80`; `python/tests/test_catalogue.py:47-48` | README implies `cat.box` works for any catalogue; it is `None` for every text catalogue (Box is not in the text header), and `cat.a` has 5 decimals from the text (`f8.5`) but 9 significant digits from the sidecar (`es16.8`), so `a` differs between the two forms of one output in the 6th decimal. Also: when `GLAMIO_TEST_TMPDIR` is unset the fixtures go to `tempfile.gettempdir()` (= `$TMPDIR` or `/tmp`); they are < 10 kB and removed in `finally`, which satisfies the /tmp rule, but the README does not mention the variable. | Pilot header: `A    =    1.00000`; sidecar: `a =   1.00000000E+00`. | Say in the README that `box`, `rows`, `realization` come only from the sidecar and that text `a`/`da` carry 5 decimals; document `GLAMIO_TEST_TMPDIR` in the README's test instructions. |
| R3-18 | Nit | `python/pyproject.toml`; `python/glamio/__init__.py:8` | No `license`, `authors`, `[project.optional-dependencies] test = ["pytest"]`; `__all__` omits `RECORD_BYTES` and `COLUMN_UNITS` although `asserts.py`/tests import from the submodule. | — | Add the metadata and export the constants (or state that `glamio.catalogue` is the stable import path). |

## 3. Verified OK

- **Record layout.** `CATALOGUE_DTYPE` has 24 fields in the order of the Fortran binary write
  (`PMP2linker.f90:703-706`): `x y z | VxMaxx VyMaxx VzMaxx | Mvir Mtotal | rr | Vrms | VmaxM |
  int(iHalo,4) | Cvir | aNpart | 0_4 | Xoff | VirRat | LambdaM | 1.e3*RadRms | Axba Axca | Xax Yax Zax`
  → `x y z vx vy vz Mbound Mtot Rvir Vrms Vcirc row Cvir Npart distinct Xoff virial_ratio lambda
  RadRMS b_a c_a ax_x ax_y ax_z`. Every item is 4 bytes: the locals are `real*4` (`:634`), the
  module arrays `Real*4, ALLOCATABLE` (`:97-102`), `1.e3*RadRms(ip)` is default real,
  `int(iHalo,4)` and `0_4` are integer*4; `Nhalo` is `Integer*4` (`:111`). `itemsize == 96`, byte
  order explicit `<`. No real*8 promotion flag is involved (no Makefile `-r8`; the test suite's
  size/96 check also pins this).
- **Sidecar layout lines.** The literal `columns` and `dtypes` strings in `WriteCatalogueMeta`
  equal `list(CATALOGUE_DTYPE.names)` and the per-field kinds (programmatic comparison); the
  `units` line has exactly 24 tokens.
- **`_check_layout` catches** (probed): two names exchanged (`Mbound`/`Mtot`), an extra column, a
  missing column, `>f4`, `f8`, a non-numeric `record_bytes`, and the tests' `Mvir` rename,
  `i4→f4` swap and `record_bytes = 100`. Unparsable dtype codes raise `TypeError` in numpy 2.4 and
  are converted to `ValueError`.
- **Binary reader.** Size not a multiple of 96 → `ValueError` naming the file; `rows` mismatch
  (fewer or more whole records) → `ValueError` naming both files; missing `.bin` or `.meta` →
  `FileNotFoundError` naming the missing file; zero-byte file with `rows = 0` → empty `Catalogue`
  with the right dtype and `velocities_at_av().shape == (0, 3)`; `__getitem__` returns a field
  view of the structured array (no copy); `positions()`/`velocities_at_av()` are the documented
  float64 copies; `Catalogue.path` is the `.bin` (binary) or the `.DAT` (text).
- **Fortran-exact sidecar parses.** A `.meta` transcribed from the `WriteCatalogueMeta` formats
  (`i10` rows, `i6` step, `f12.4` Box, `es16.8` a/z/a_v/da with the trailing `a_v` comment, the
  `velocity_scale` value containing a second `=`, the `centre` line with four `=`) gives
  `a = 1.0`, `z = 0.0`, `a_v = 0.9973`, `step = 181`, `box = 1024.0`, intact `header` and
  `text_catalogue`; `velocities_at_av()` = `v × 1.0027073`.
- **Text reader on the pilot** (777,651 rows, 220 MB): 1.8 s, 410 MB peak RSS; rows = `wc -l` − 8;
  all 15 header keys match the file (`a=1.00000`, `step_da=0.00540`, `step=181`, `N_row=1024`,
  `N_grid=2048`, `Omega_m=0.3000`, `Omega_L=0.7000`, `hubble=0.7000`, `buffer_width=15.3750`,
  `particle_mass=8.3220E+10`, `overdensity=338.3000`, `product=CatshortV`, `format=text`,
  `z=0.000000`, `a_v=0.997300`); first and last rows equal an independent `float()` parse in all
  24 columns; `row` = 1..N and `distinct` = 0 throughout; `x` in [0.0023, 1023.9985] (consistent
  with `[0, Box)`); F-form (`771.00`, `-0.3133`, `0.1190`) and E-form (`4.0263E-02`, `-9.9848E-02`,
  594,206 small `Xoff` values, 693,537 negative virial ratios) all parsed.
- **Header regexes vs the Fortran formats** (`PMP2linker.f90:2084-2104`): `2(a,f8.5)` for A/Step,
  `3(a,i5)` for I/Nrow/Ngrid, `6(a,f8.4)` for the cosmology line, `(a,i4)`, `(a,1p,g12.4)` for the
  particle mass, `(a,f8.4)` for the overdensity; the padded `character` labels are absorbed by
  `\s*=\s*`. `a_v = A − Step/2` formatted with 6 decimals is exact because A and Step carry 5
  decimals (f8.5). A real*4 value never needs a 3-digit exponent, so Fortran's `1.0000+100` form
  (unparsable by `float`) cannot occur in a catalogue; `np.loadtxt` handles every form the
  `(3f11.4,3x,3f10.2,1p,2g12.4,g12.5,26g12.4)` format can emit for finite real*4 values.
- **`np.loadtxt(..., ndmin=2)`** on a header-only file returns shape `(0, 1)`, handled by the
  `raw.size == 0` branch (test passes, warning suppressed); a 23-column row raises `ValueError`
  (from numpy when mixed with 24-column rows, from glamio when uniform).
- **Properties.** `a`, `z` (explicit key wins over `1/a − 1`), `a_v`, `box`, `step` (`int(round())`)
  behave as documented; `meta_float` converts the first token only, accepts `D` exponents,
  returns the default for absent/empty values and raises `ValueError` naming the key for a
  non-numeric token; `read_meta` splits on the first `=`, strips, skips blank/`=`-less lines and
  keeps the last duplicate.
- **`velocities_at_av()`** returns `v_written × a/a_v` as a float64 `(N, 3)` copy and raises when
  `a` or `a_v` is missing or `a_v == 0`. This is the algebraic consequence of the sidecar's own
  statement (`vx = code velocity × 100 Box/N_grid / a`): if the leapfrog momentum at `a_v` was
  divided by `a`, dividing by `a_v` instead multiplies by `a/a_v` (exact under that convention,
  not a first-order approximation). Whether the convention statement itself is right is
  Review 2's question.
- **Packaging.** `pyproject.toml` is valid for setuptools ≥ 64 editable installs (the cosemu
  install resolves `glamio.__file__` to `python/glamio/__init__.py`, dist version 0.1.0);
  `[tool.setuptools.packages.find] include = ["glamio*"]` excludes `tests/` and `glamio.egg-info`
  from the wheel; `requires-python >= 3.10` is consistent with the syntax used (`X | None` only in
  annotations under `from __future__ import annotations`, `dict[str, str]`); runs on the env's
  Python 3.14.2 / numpy 2.4.3; `readme = "README.md"` exists; no `print(` in library code; no dead
  functions other than R3-11.
- **Repository hygiene.** `git ls-files -- python` shows exactly the six intended files;
  `python/glamio.egg-info/` (present, untracked) and `.pytest_cache/` are ignored by the new
  `.gitignore` lines, `__pycache__/` by the existing one (`git check-ignore -v`).
- **Tests.** `python -m pytest -q -p no:cacheprovider tests` → 6 passed in 0.32 s; script mode
  (`python -I tests/test_catalogue.py`) → 6 passed; fixtures honour `GLAMIO_TEST_TMPDIR` and are
  removed (0 leftovers). Both `__pycache__` directories created by the runs were deleted.

## 4. Open questions for the author

1. Should the `row == 1..N` / `distinct == 0` self-check (R3-01) be mandatory in
   `_read_catalogue_binary`? It costs one pass over 4 bytes per row and would have caught every
   mutation Review 4 is asked to try except a flipped byte in a float column.
2. Will the text and binary meta keys be unified in glamio (`da`, `rows`, `realization`,
   `format = text|binary`) so that downstream code such as `asserts.py` (`cat.meta["da"]`) works
   on either form? And should the Fortran write `format = binary` with the prose under another
   key (Review 2)?
3. Is the float32 downcast of the text columns intended as part of the API? It makes the text
   `Catalogue` layout-identical to the binary one, but it adds up to one float32 ulp to values
   that the text had pinned to 1e-4 (R3-03); a `float64` option or returning the raw `loadtxt`
   array alongside would let callers choose.
4. For the L1800 catalogues (~1 GB per `.bin`), is a `mmap` mode wanted (R3-08)?
5. Should `Catalogue.box`/`step`/`realization` for text catalogues be filled from the file name
   or from a sibling `.meta` when present, so that the text and binary objects expose the same
   properties?
6. Note for Review 4 (TEXT_TOL): the position atol of 6e-5 in `asserts.py` holds only for
   |x| < 256 Mpc/h once the text is read through glamio's float32 cast (measured worst case
   6.104e-5 for 256 <= x < 1024, affecting 8-39 % of values); the velocity atol 6e-3 (f10.2,
   half unit 5e-3, ulp ≤ 1.2e-4 for |v| < 2048 km/s) is safe.
