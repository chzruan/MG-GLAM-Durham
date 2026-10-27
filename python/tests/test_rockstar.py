"""Tests for the in-line Rockstar catalogue: reader, run paths and packing, with synthetic files.

The synthetic catalogues follow ``rockstar_inline/rs_write.c`` (raw
little-endian records, the layout in the sidecar's ``columns`` / ``dtypes`` /
``offsets`` / ``record_bytes``, ``status`` ok / skipped / failed) with a
smaller 64-byte record that leaves a 4-byte gap between two columns, so that
packing is shown to keep bytes no column covers.  The pack tests need h5py
(skipped without it).  Runs under pytest or as a plain script; fixtures go to
a fresh ``tempfile.mkdtemp`` directory under ``$GLAMIO_TEST_TMPDIR`` (or the
default temp dir) and are removed afterwards.
"""

from __future__ import annotations

import filecmp
import os
import shutil
import subprocess
import sys
import tempfile
import textwrap

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import glamio  # noqa: E402

try:
    import h5py
    from glamio import pack as gpack
except ImportError:  # no h5py: the pack tests are skipped
    h5py = None

try:
    import pytest
except ImportError:  # plain-script mode
    pytest = None

BOX = 150.0
HEADER = "ic2pm: 2LPTic ingest, v at a-da/2, IDz"
SCHEDULE = ("# MG-GLAM output schedule: one line per analysis output (snapshot number = step)\n"
            "# chi: comoving distance [Mpc/h], flat matter+Lambda, Om= 0.30890 OmL= 0.69110, c/H0 = 2997.92458 Mpc/h\n"
            "#  step               a               z          chi exact class\n"
            "    70  2.85714298E-01  2.49999976E+00    4045.9986     1     1\n"
            "    84  4.76190507E-01  1.09999990E+00    2464.3611     1     2\n"
            "    98  6.61375642E-01  5.12000084E-01    1345.2087     1     1\n")
A = {70: 0.285714298, 84: 0.476190507, 98: 0.661375642, 112: 0.8}       # 112: not in the schedule
#: (name, dtype, offset, unit) of the synthetic record; bytes 56..59 belong to no column.
COLUMNS = [("id", "i8", 0, "-"), ("x", "f4", 8, "Mpc/h"), ("y", "f4", 12, "Mpc/h"), ("z", "f4", 16, "Mpc/h"),
           ("vx", "f4", 20, "km/s"), ("vy", "f4", 24, "km/s"), ("vz", "f4", 28, "km/s"), ("m", "f4", 32, "Msun/h"),
           ("r", "f4", 36, "kpc/h"), ("num_p", "i8", 40, "-"), ("fof_group", "i8", 48, "-"),
           ("strict_crossed", "i4", 60, "-")]
RECORD_BYTES = 64


def _tmpdir() -> str:
    base = os.environ.get("GLAMIO_TEST_TMPDIR") or None
    if base:
        os.makedirs(base, exist_ok=True)
    return tempfile.mkdtemp(prefix="glamio-rockstar-", dir=base)


def _write(path: str, text: str) -> str:
    with open(path, "w") as fh:
        fh.write(text)
    return path


def _raises(exc, fn, *args, **kwargs):
    try:
        fn(*args, **kwargs)
    except exc as e:
        return str(e)
    raise AssertionError(f"{fn.__name__} did not raise {exc.__name__}")


def _dtype() -> np.dtype:
    return np.dtype({"names": [c[0] for c in COLUMNS], "formats": ["<" + c[1] for c in COLUMNS],
                     "offsets": [c[2] for c in COLUMNS], "itemsize": RECORD_BYTES})


def _meta(step: int, box: int, rows: int, status: str, reason: str, columns=COLUMNS,
          record_bytes: int = RECORD_BYTES, product: str = "Rockstar") -> str:
    a, da = A[step], 0.008
    text = (f"product = {product}\nformat = raw little-endian records, {record_bytes} bytes per distinct halo; no header\n"
            f"record_bytes = {record_bytes}\nstatus = {status}\nreason = {reason}\n"
            f"a = {a:.8f}\na_x = {a:.8f}   (epoch of the positions; same as a)\nz = {1 / a - 1:.8f}\n"
            f"a_v = {a - da / 2:.8f}   (epoch of the velocities, a - da/2)\nda = {da:.8f}\n"
            f"velocity_scale = 1/a_v\nstep = {step}\nrealization = {box}\nBox = {BOX:.4f}\nN_grid = 72\nN_row = 36\n"
            "N_p = 46656\nOmega_m = 0.308900\nOmega_L = 0.691100\nhubble = 0.677400\n"
            "particle_mass = 6.03780000e+09   (Msun/h, Rockstar's: CRITICAL_DENSITY Omega_m Box^3 / N_p)\n"
            "MG_flag = 0\nMG_model = 1\nfinder = in-line Rockstar of MG-GLAM (rockstar_inline/), test\n"
            f"config_file = Rockstar.config\nconfig_fnv1a64 = 0123456789abcdef\nheader = {HEADER}\nrows = {rows}\n")
    if status == "ok":
        text += "fof_groups = 7\nconfig.FORCE_RES = 0.1\nconfig.PARALLEL_IO_SERVER_INTERFACE = \n"   # an empty value
    text += ("columns = " + " ".join(c[0] for c in columns) + "\ndtypes = " + " ".join(c[1] for c in columns)
             + "\noffsets = " + " ".join(str(c[2]) for c in columns) + "\nunits = " + " ".join(c[3] for c in columns) + "\n")
    text += "".join(f"column_{i + 1:02d} = {c[0]}: {c[3]} (column {c[0]})\n" for i, c in enumerate(columns))
    return text


def _records(rows: int, seed: int = 2026) -> bytes:
    """``rows`` synthetic records; the bytes outside the columns are random (non-zero), as padding may be."""
    rng = np.random.default_rng(seed)
    raw = rng.integers(1, 255, size=rows * RECORD_BYTES, dtype=np.uint8)
    rec = raw.view(_dtype())
    rec["id"] = np.arange(rows)
    for c in ("x", "y", "z"):
        rec[c] = rng.uniform(0, BOX, rows)
    for c in ("vx", "vy", "vz"):
        rec[c] = rng.normal(0, 300, rows)
    rec["m"] = np.geomspace(1e11, 1e14, rows) if rows else []
    rec["r"] = rng.uniform(50, 900, rows)
    rec["num_p"] = rng.integers(20, 10000, rows)
    rec["fof_group"] = np.arange(rows) // 2
    rec["strict_crossed"] = 31
    return raw.tobytes()


def _rockstar(catdir: str, step: int, box: int = 1, rows: int = 5, status: str = "ok", reason: str = "-",
              **meta_kw) -> str:
    """Write ``Rockstar.<step>.<box>.{bin,meta}`` into ``catdir``; return the stem."""
    os.makedirs(catdir, exist_ok=True)
    stem = os.path.join(catdir, f"Rockstar.{step:04d}.{box:04d}")
    with open(stem + ".bin", "wb") as fh:
        fh.write(_records(rows, seed=step) if status == "ok" else b"")
    _write(stem + ".meta", _meta(step, box, rows if status == "ok" else 0, status, reason, **meta_kw))
    return stem


def _lens(rundir: str, step: int, box: int = 1, ns: int = 2, nl: int = 4) -> None:
    stem = os.path.join(rundir, f"LensPlanes.{step:04d}.{box:04d}")
    a = A[step]
    (1.0 + 0.01 * np.random.default_rng(step).standard_normal((3 * ns, nl, nl))).astype("<f4").tofile(stem + ".bin")
    _write(stem + ".meta", "product = LensPlanes\nformat = raw little-endian float32, no record markers\n"
                           f"step = {step}\nrealization = {box}\na = {a:.8E}\nz = {1 / a - 1:.8E}\na_v = {a - 0.004:.8E}\n"
                           f"Box = {BOX:.4f}\nN_grid = 72\nN_row = 36\nN_p = 46656\nOmega_m = 0.308900\n"
                           f"Omega_L = 0.691100\nhubble = 0.677400\nMG_flag = 0\nMG_model = 1\nfR0 = none\n"
                           f"header = {HEADER}\nN_slab = {ns}\ndepth = {BOX / ns:.4f}\nN_lens = {nl}\nN_planes = {3 * ns}\n")


def _run_fixture(d: str, box: int = 1) -> str:
    """``Run<box>/``: lens planes at the three outputs, Rockstar ok (5 rows) at step 70 and skipped (0 rows) at 98."""
    rundir = os.path.join(d, f"Run{box}")
    os.makedirs(rundir)
    _write(os.path.join(rundir, "OutputSchedule.dat"), SCHEDULE)
    _write(os.path.join(rundir, "BDM.config"), "MaxMemory = 900\n")
    _write(os.path.join(rundir, "Rockstar.config"), 'FILE_FORMAT = "GLAM_INLINE"\nFORCE_RES = 0.1\n')
    _write(os.path.join(rundir, "main.log"), "step 70\nstep 84\nstep 98\n")
    for step in (70, 84, 98):
        _lens(rundir, step, box)
    cat = os.path.join(rundir, "CATALOGS")
    _rockstar(cat, 70, box, rows=5)
    _rockstar(cat, 98, box, status="skipped", reason="memory pre-check: VmRSS 0.4 GiB + FoF 0.6 GiB > ceiling 0.5 GiB")
    return rundir


def _tree(root: str) -> list[str]:
    out = []
    for dirpath, _, names in os.walk(root):
        out += [os.path.relpath(os.path.join(dirpath, n), root) for n in names]
    return sorted(out)


def _same_tree(a: str, b: str) -> bool:
    ta, tb = _tree(a), _tree(b)
    return ta == tb and all(filecmp.cmp(os.path.join(a, f), os.path.join(b, f), shallow=False) for f in ta)


def _quiet(*_args, **_kwargs):
    pass


def _skip_without_h5py() -> bool:
    if h5py is not None:
        return False
    if pytest is not None:
        pytest.skip("h5py not installed")
    return True


# ----------------------------------------------------------------- reader
def test_read_ok():
    d = _tmpdir()
    try:
        stem = _rockstar(d, 70, rows=5)
        want = np.frombuffer(_records(5, seed=70), dtype=_dtype())
        for path in (stem, stem + ".bin", stem + ".meta"):
            for mmap in (False, True):
                cat = glamio.read_rockstar(path, mmap=mmap)
                assert cat.data.dtype == _dtype() and cat.data.tobytes() == want.tobytes(), (path, mmap)
                assert isinstance(cat.data, np.memmap) == mmap
        assert len(cat) == 5 and cat.columns == [c[0] for c in COLUMNS] and cat.path == stem + ".bin"
        assert cat.status == "ok" and cat.ok and cat.reason is None and cat.record_bytes == RECORD_BYTES
        assert cat.product == "Rockstar" and cat.step == 70 and cat.realization == 1 and cat.box == BOX
        assert abs(cat.a - A[70]) < 1e-8 and abs(cat.a_v - (A[70] - 0.004)) < 1e-8 and cat.n_p == 46656
        assert cat.particle_mass == 6.0378e9 and cat.header == HEADER
        assert cat.units["m"] == "Msun/h" and cat.units["id"] == "-" and len(cat.units) == len(COLUMNS)
        assert cat.descriptions["vx"] == "km/s (column vx)" and len(cat.descriptions) == len(COLUMNS)
        assert np.array_equal(cat["id"], np.arange(5)) and np.array_equal(cat["m"], want["m"])
        pos, vel = cat.positions(), cat.velocities()
        assert pos.shape == (5, 3) and pos.dtype == np.float64 and np.array_equal(pos[:, 1], want["y"].astype(np.float64))
        assert vel.shape == (5, 3) and np.array_equal(vel[:, 2], want["vz"].astype(np.float64))
        assert cat.meta["config.PARALLEL_IO_SERVER_INTERFACE"] == "" and "RockstarCatalogue(5 haloes" in repr(cat)
        dt = glamio.rockstar_dtype(glamio.read_meta(stem + ".meta"))
        assert dt == _dtype() and dt.itemsize == RECORD_BYTES
        empty = _rockstar(os.path.join(d, "empty"), 70, rows=0)                # an ok catalogue may hold no halo
        cat = glamio.read_rockstar(empty, mmap=True)
        assert len(cat) == 0 and cat.ok and cat.data.dtype == _dtype()
    finally:
        shutil.rmtree(d, ignore_errors=True)


def test_status_skipped_and_failed():
    d = _tmpdir()
    try:
        for status, reason in (("skipped", "memory pre-check: VmRSS 0.4 GiB + FoF 0.6 GiB > ceiling 0.5 GiB"),
                               ("failed", "Rockstar.config: FORCE_RES > 0 violated")):
            stem = _rockstar(os.path.join(d, status), 98, status=status, reason=reason)
            assert os.path.getsize(stem + ".bin") == 0
            msg = _raises(ValueError, glamio.read_rockstar, stem)
            assert f"status = {status}" in msg and reason in msg and "allow_status" in msg
            for mmap in (False, True):
                cat = glamio.read_rockstar(stem, allow_status=True, mmap=mmap)
                assert len(cat) == 0 and cat.status == status and not cat.ok and cat.reason == reason
                assert cat.data.dtype == _dtype() and cat.columns == [c[0] for c in COLUMNS] and abs(cat.a - A[98]) < 1e-8
        stem = os.path.join(d, "skipped", "Rockstar.0098.0001")
        meta = open(stem + ".meta").read()
        _write(stem + ".meta", meta.replace("rows = 0", "rows = 2"))       # skipped but rows: inconsistent
        with open(stem + ".bin", "wb") as fh:
            fh.write(_records(2))
        assert "no rows" in _raises(ValueError, glamio.read_rockstar, stem, allow_status=True)
        _write(stem + ".meta", meta.replace("status = skipped", "status = partial"))
        open(stem + ".bin", "wb").close()
        assert "expected one of" in _raises(ValueError, glamio.read_rockstar, stem, allow_status=True)
        _write(stem + ".meta", meta.replace("status = skipped\n", ""))
        assert "status" in _raises(ValueError, glamio.read_rockstar, stem, allow_status=True)
    finally:
        shutil.rmtree(d, ignore_errors=True)


def test_layout_errors():
    d = _tmpdir()
    try:
        stem = _rockstar(d, 70, rows=4)
        good = open(stem + ".meta").read()
        raw = open(stem + ".bin", "rb").read()
        with open(stem + ".bin", "wb") as fh:                                 # truncated
            fh.write(raw[:-1])
        assert "rows = 4" in _raises(ValueError, glamio.read_rockstar, stem)
        with open(stem + ".bin", "wb") as fh:
            fh.write(raw)
        bad_metas = {
            "product": good.replace("product = Rockstar", "product = CatshortV"),
            "big-endian": good.replace("dtypes = i8 f4", "dtypes = i8 >f4"),
            "complex": good.replace("dtypes = i8 f4", "dtypes = i8 c8"),
            "count": good.replace("offsets = 0 8 ", "offsets = 8 "),
            "overlap": good.replace("offsets = 0 8 12", "offsets = 0 6 12"),
            "outside": good.replace(" 60\nunits", " 62\nunits"),
            "record": good.replace("record_bytes = 64", "record_bytes = 0"),
            "units": good.replace("units = - ", "units = "),
            "names": good.replace("columns = id x y", "columns = id x x"),
            "rows": good.replace("rows = 4", "rows = -4"),
        }
        for why, text in bad_metas.items():
            assert text != good, why
            _write(stem + ".meta", text)
            try:
                glamio.read_rockstar(stem)
            except ValueError:
                continue
            raise AssertionError(f"layout variant {why!r} was not refused")
        _write(stem + ".meta", good)
        assert len(glamio.read_rockstar(stem)) == 4
        swapped = np.frombuffer(raw, dtype=_dtype()).copy()                   # a byte-swapped id column
        swapped["id"] = swapped["id"].byteswap()
        swapped.tofile(stem + ".bin")
        assert "0..3" in _raises(ValueError, glamio.read_rockstar, stem)
        os.remove(stem + ".meta")
        assert "missing" in _raises(FileNotFoundError, glamio.read_rockstar, stem)
    finally:
        shutil.rmtree(d, ignore_errors=True)


def test_run_paths_native():
    d = _tmpdir()
    try:
        rundir = _run_fixture(d)
        stem = os.path.join(rundir, "CATALOGS", "Rockstar.0070.0001")
        assert glamio.product_path(rundir, "rockstar", 70) == stem
        assert glamio.products_at(rundir, 70)["rockstar"] == stem
        assert glamio.products_at(rundir, 98)["rockstar"].endswith("Rockstar.0098.0001")
        assert "rockstar" not in glamio.products_at(rundir, 84)
        assert glamio.read_rockstar(glamio.product_path(rundir, "rockstar", 98), allow_status=True).status == "skipped"
        assert glamio.run_file(rundir, "Rockstar.config") == os.path.join(rundir, "Rockstar.config")
        os.remove(stem + ".meta")                                            # an incomplete pair is not a product
        assert "rockstar" not in glamio.products_at(rundir, 70) and glamio.product_path(rundir, "rockstar", 70) == stem
        _raises(FileNotFoundError, glamio.read_rockstar, glamio.product_path(rundir, "rockstar", 70))
        prelude = f"import sys\nsys.path.insert(0, {os.path.dirname(os.path.dirname(os.path.abspath(__file__)))!r})\n"
        reads = textwrap.dedent(f"""
            sys.modules['h5py'] = None
            import glamio
            r = {rundir!r}
            p = glamio.products_at(r, 98)["rockstar"]
            assert glamio.read_rockstar(p, allow_status=True).status == "skipped"
            assert ".h5" not in glamio.product_path(r, "rockstar", 70)
            """)
        subprocess.run([sys.executable, "-I", "-c", prelude + reads], check=True)   # native reads need no h5py
    finally:
        shutil.rmtree(d, ignore_errors=True)


# ----------------------------------------------------------------- packing
def test_pack_roundtrip_with_rockstar():
    if _skip_without_h5py():
        return
    d = _tmpdir()
    old = gpack.BLOCK_BYTES
    try:
        rundir = _run_fixture(d)
        out = os.path.join(d, "packed", "Run1")
        gpack.pack(rundir, out, log=_quiet)
        for path in gpack._packed_files([out]):
            problems, stats = gpack.verify_file(path, against=rundir)
            assert problems == [] and stats["originals_checked"] == stats["entries"], (path, problems)
        with h5py.File(os.path.join(out, "Output.0070.0001.h5"), "r") as h:
            g = h["Rockstar"]
            assert h.attrs["products"].split() == ["Rockstar", "LensPlanes"] and h.attrs["Box"] == BOX
            assert g["data"].dtype == _dtype() and g["data"].shape == (5,) and g["data"].chunks is None
            assert g["data"].id.get_offset() % 4096 == 0
            assert g.attrs["status"] == "ok" and g.attrs["rows"] == "5" and g.attrs["config.PARALLEL_IO_SERVER_INTERFACE"] == ""
            assert g.attrs["meta_text"].decode() == open(os.path.join(rundir, "CATALOGS", "Rockstar.0070.0001.meta")).read()
            paths = sorted(e["path"].decode() for e in h["manifest"][()])
            assert "CATALOGS/Rockstar.0070.0001.bin" in paths and "CATALOGS/Rockstar.0070.0001.meta" in paths
        with h5py.File(os.path.join(out, "Output.0098.0001.h5"), "r") as h:
            assert h["Rockstar/data"].shape == (0,) and h["Rockstar"].attrs["status"] == "skipped"
        with h5py.File(os.path.join(out, "RunInfo.0001.h5"), "r") as h:
            assert h["Rockstar.config"][()].tobytes().startswith(b"FILE_FORMAT") and "logs/Rockstar.config" not in h

        for step in (70, 84, 98):                                         # packed reads equal the native ones
            found = glamio.products_at(out, step)
            assert ("rockstar" in found) == (step != 84)
            if step == 84:
                continue
            pn, pp = glamio.product_path(rundir, "rockstar", step), glamio.product_path(out, "rockstar", step)
            assert pp == found["rockstar"] and pp.endswith(f"Output.{step:04d}.0001.h5/Rockstar")
            for mmap in (False, True):
                a, b = glamio.read_rockstar(pn, allow_status=True, mmap=mmap), glamio.read_rockstar(pp, allow_status=True, mmap=mmap)
                assert a.data.dtype == b.data.dtype and a.data.tobytes() == b.data.tobytes() and a.meta == b.meta
                assert (a.status, a.reason, a.a, a.box) == (b.status, b.reason, b.a, b.box)
                if mmap and len(b):
                    assert isinstance(b.data, np.memmap)
        assert "skipped" in _raises(ValueError, glamio.read_rockstar, glamio.product_path(out, "rockstar", 98))
        rc = glamio.run_file(out, "Rockstar.config")
        assert rc.endswith("RunInfo.0001.h5/Rockstar.config") and glamio.read_text(rc) == glamio.read_text(glamio.run_file(rundir, "Rockstar.config"))

        rebuilt = os.path.join(d, "rebuilt")                              # byte-identical, the 0-row .bin and the gaps included
        assert gpack.unpack([out], rebuilt, log=_quiet) == len(_tree(rundir))
        assert _same_tree(rundir, rebuilt)

        gpack.BLOCK_BYTES = 40                                            # one record per block
        small = os.path.join(d, "small")
        gpack.pack(rundir, small, box=1, log=_quiet)
        gpack.unpack([small], os.path.join(d, "rebuilt_small"), log=_quiet)
        assert _same_tree(rundir, os.path.join(d, "rebuilt_small"))
        gpack.BLOCK_BYTES = old

        copy = os.path.join(d, "copy", "Run1")                            # pack in place and prune, then read from the packed files
        shutil.copytree(rundir, copy)
        gpack.pack(copy, log=_quiet)
        assert gpack.status(copy)[0] == "packed"
        assert gpack.prune(copy, yes=True, log=_quiet) == len(_tree(rundir))
        assert gpack.status(copy)[0] == "pruned" and not os.path.exists(os.path.join(copy, "CATALOGS"))
        assert len(glamio.read_rockstar(glamio.product_path(copy, "rockstar", 70))) == 5
        assert glamio.read_rockstar(glamio.product_path(copy, "rockstar", 98), allow_status=True).status == "skipped"
        assert gpack.main(["verify", copy]) == 0
    finally:
        gpack.BLOCK_BYTES = old
        shutil.rmtree(d, ignore_errors=True)


def test_pack_refuses_rockstar_problems():
    if _skip_without_h5py():
        return
    d = _tmpdir()
    try:
        rundir = _run_fixture(d)
        cat = os.path.join(rundir, "CATALOGS")
        for bad, why in ((os.path.join(cat, ".Rockstar.0070.0001.bin.tmp"), "interrupted Rockstar catalogue write"),
                         (os.path.join(cat, ".Rockstar.0070.0001.meta.tmp"), "interrupted Rockstar catalogue write"),
                         (os.path.join(cat, "Rockstar.0084.0001.bin"), "needs both"),
                         (os.path.join(cat, "Rockstar.0070.0002.meta"), "box tag"),
                         (os.path.join(cat, "Rockstar.0070.0001.txt"), "Rockstar.<step>.<box>.bin/.meta")):
            _write(bad, "x")
            try:
                gpack.pack(rundir, os.path.join(d, "out"), log=_quiet)
                raise AssertionError(f"{bad} was not refused")
            except gpack.PackError as exc:
                assert why in str(exc), (why, str(exc))
            os.remove(bad)
        _rockstar(cat, 112)                                               # a step outside the schedule
        try:
            gpack.pack(rundir, os.path.join(d, "out"), log=_quiet)
            raise AssertionError("a Rockstar catalogue outside the schedule was not refused")
        except gpack.PackError as exc:
            assert "step 112" in str(exc) and "not in OutputSchedule.dat" in str(exc)
        os.remove(os.path.join(cat, "Rockstar.0112.0001.bin"))
        os.remove(os.path.join(cat, "Rockstar.0112.0001.meta"))
        with open(os.path.join(cat, "Rockstar.0070.0001.bin"), "ab") as fh:  # a size that does not match rows
            fh.write(b"\0" * RECORD_BYTES)
        _raises(ValueError, gpack.pack, rundir, os.path.join(d, "out"), log=_quiet)
        assert not [n for n in os.listdir(os.path.join(d, "out")) if n.startswith("Output.0070")]   # nothing left
    finally:
        shutil.rmtree(d, ignore_errors=True)


def test_rockstar_config_packed_as_log():
    """A RunInfo that holds Rockstar.config under logs/ (glamio 0.3) is still found by run_file."""
    if _skip_without_h5py():
        return
    d = _tmpdir()
    try:
        rundir = _run_fixture(d)
        out = os.path.join(d, "packed", "Run1")
        obj = glamio.run.RUN_OBJECTS.pop("Rockstar.config")       # the placement before glamio 0.4
        try:
            gpack.pack(rundir, out, log=_quiet)
        finally:
            glamio.run.RUN_OBJECTS["Rockstar.config"] = obj
        with h5py.File(os.path.join(out, "RunInfo.0001.h5"), "r") as h:
            assert "logs/Rockstar.config" in h and "Rockstar.config" not in h
        p = glamio.run_file(out, "Rockstar.config")
        assert p.endswith("RunInfo.0001.h5/logs/Rockstar.config") and glamio.read_text(p).startswith("FILE_FORMAT")
    finally:
        shutil.rmtree(d, ignore_errors=True)


TESTS = [test_read_ok, test_status_skipped_and_failed, test_layout_errors, test_run_paths_native,
         test_pack_roundtrip_with_rockstar, test_pack_refuses_rockstar_problems, test_rockstar_config_packed_as_log]

if __name__ == "__main__":
    failed = 0
    for t in TESTS:
        try:
            t()
            print(f"PASS {t.__name__}")
        except Exception as exc:                         # noqa: BLE001
            failed += 1
            print(f"FAIL {t.__name__}: {type(exc).__name__}: {exc}")
    print(f"{len(TESTS) - failed}/{len(TESTS)} passed")
    sys.exit(1 if failed else 0)
