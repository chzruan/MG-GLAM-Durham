"""Tests for glamio.catalogue: binary roundtrip, layout checks, text reader.

Runs under pytest (``python -m pytest -q tests``) or as a plain script
(``python tests/test_catalogue.py``).  Fixtures are written to a fresh
``tempfile.mkdtemp`` directory under ``$GLAMIO_TEST_TMPDIR`` (or the default
temp dir when unset) and removed afterwards.
"""

from __future__ import annotations

import os
import shutil
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import glamio  # noqa: E402
from glamio import CATALOGUE_DTYPE, Catalogue, read_catalogue, read_meta  # noqa: E402
from glamio.catalogue import RECORD_BYTES, read_catalogue_text  # noqa: E402
from glamio.meta import meta_float  # noqa: E402

try:
    import pytest
except ImportError:  # plain-script mode
    pytest = None

NAMES = list(CATALOGUE_DTYPE.names)
DTYPE_CODES = " ".join("i4" if CATALOGUE_DTYPE[n].kind == "i" else "f4" for n in NAMES)

HEADER_8 = (
    "ic2pm: 2LPTic ingest, v at a-da/2 [BDM finder v3]\n"
    " A    =    1.00000 Step =    0.00540\n"
    " I    =     181 Nrow =    1024 Ngrid=    2048\n"
    " Omega_0=   0.3000 Omega_L=   0.7000 hubble =   0.7000 buffer width (Mpch) =                   15.3750\n"
    " Number of radial bins                =  100\n"
    " Mass of smallest particle (Msunh)    =   8.3220E+10\n"
    " Overdensity limit                    = 338.3000\n"
    "     XYZ(Mpch)                          Vxyz(km/s)                  Mbound     Mtot/Msunh    Rvir(kpch)"
    " Vrms(km/s) Vcirc(km/s)       Nhalo  Cvir    Nparticles  Distinct/Sub   Xoff      2K/Ep-1   Lambda"
    "   RadRMS/k  b/a  c/a MajorAxis:  x      y      z  \n"
)


def _mkdir() -> str:
    return tempfile.mkdtemp(prefix="glamio-test-", dir=os.environ.get("GLAMIO_TEST_TMPDIR"))


if pytest is not None:

    @pytest.fixture
    def workdir():
        d = _mkdir()
        try:
            yield d
        finally:
            shutil.rmtree(d, ignore_errors=True)


def synthetic_catalogue(n: int = 50, seed: int = 2026) -> np.ndarray:
    """Structured array with random but realistic halo values."""
    rng = np.random.default_rng(seed)
    arr = np.empty(n, dtype=CATALOGUE_DTYPE)
    for c in ("x", "y", "z"):
        arr[c] = rng.uniform(0.0, 1024.0, n)
    for c in ("vx", "vy", "vz"):
        arr[c] = rng.normal(0.0, 300.0, n)
    arr["Mbound"] = 10.0 ** rng.uniform(12.0, 15.0, n)
    arr["Mtot"] = arr["Mbound"] * rng.uniform(1.0, 2.0, n)
    arr["Rvir"] = rng.uniform(200.0, 2000.0, n)
    arr["Vrms"] = rng.uniform(100.0, 1500.0, n)
    arr["Vcirc"] = rng.uniform(100.0, 1800.0, n)
    arr["row"] = np.arange(1, n + 1)
    arr["Cvir"] = rng.uniform(3.0, 25.0, n)
    arr["Npart"] = np.floor(arr["Mbound"] / 8.322e10)
    arr["distinct"] = 0
    arr["Xoff"] = rng.uniform(0.0, 0.2, n)
    arr["virial_ratio"] = rng.normal(-0.3, 0.1, n)
    arr["lambda"] = rng.uniform(0.0, 0.1, n)
    arr["RadRMS"] = arr["Rvir"] * rng.uniform(0.2, 0.5, n)
    arr["b_a"] = rng.uniform(0.5, 1.0, n)
    arr["c_a"] = rng.uniform(0.3, 0.9, n)
    ax = rng.normal(size=(n, 3))
    ax /= np.linalg.norm(ax, axis=1)[:, None]
    arr["ax_x"], arr["ax_y"], arr["ax_z"] = ax[:, 0], ax[:, 1], ax[:, 2]
    return arr


def write_binary(stem: str, arr: np.ndarray, extra: dict[str, str] | None = None,
                 columns: str | None = None) -> None:
    """Write ``<stem>.bin`` with ``.tofile`` and a hand-written ``<stem>.meta``."""
    arr.tofile(stem + ".bin")
    meta = {
        "product": "CatshortV",
        "format": "binary",
        "columns": " ".join(NAMES) if columns is None else columns,
        "dtypes": DTYPE_CODES,
        "rows": str(len(arr)),
        "record_bytes": str(RECORD_BYTES),
        "a": "0.5",
        "a_v": "0.4973",
        "Box": "1024.0",
        "step": "100",
    }
    if extra:
        meta.update(extra)
    with open(stem + ".meta", "w") as fh:
        for k, v in meta.items():
            fh.write(f"{k} = {v}\n")


def _row_text(values: list[float]) -> str:
    """Python equivalent of the Fortran format (3f11.4,3x,3f10.2,1p,2g12.4,g12.5,26g12.4)."""
    s = "".join(f"{v:11.4f}" for v in values[0:3]) + "   "
    s += "".join(f"{v:10.2f}" for v in values[3:6])
    s += "".join(f"{v:12.4E}" for v in values[6:])
    return s + "\n"


TEXT_ROWS = [
    [394.8253, 5.1722, 0.0920, -361.42, -355.91, 8.70, 7.3234e12, 1.2649e13, 771.00, 191.1, 303.5,
     1, 22.96, 88.0, 0, 4.0263e-2, -0.3133, 2.6562e-2, 275.4, 0.8410, 0.8140, -9.9848e-2, 0.9850, 0.1407],
    [163.9297, 5.8444, 1023.9327, -215.26, 80.88, 19.00, 1.0902e13, 1.7976e13, 806.92, 220.9, 345.7,
     2, 19.56, 131.0, 0, 7.0409e-2, -0.2776, 4.2018e-2, 323.3, 0.8255, 0.7728, 0.5978, 0.3271, 0.7319],
    [512.0001, 1000.25, 17.5, 1234.56, -0.01, 999.99, 3.5e14, 4.1e14, 1950.3, 1100.2, 1402.7,
     3, 4.21, 4205.0, 0, 0.1234, 0.0512, 0.0101, 801.5, 0.6102, 0.4411, 0.0, -1.0, 0.0],
]


def write_text(path: str, rows: list[list[float]] = TEXT_ROWS) -> list[str]:
    """Write a .DAT fixture; return the formatted row strings."""
    lines = [_row_text(r) for r in rows]
    with open(path, "w") as fh:
        fh.write(HEADER_8)
        fh.writelines(lines)
    return lines


# ---------------------------------------------------------------------------


def test_binary_roundtrip(workdir: str) -> None:
    arr = synthetic_catalogue()
    stem = os.path.join(workdir, "CatshortV.0100.0001")
    write_binary(stem, arr)

    cats = [read_catalogue(stem), read_catalogue(stem + ".bin"), read_catalogue(stem + ".meta")]
    for cat in cats:
        assert isinstance(cat, Catalogue)
        assert len(cat) == 50
        assert cat.columns == NAMES
        assert cat.data.dtype == CATALOGUE_DTYPE
        np.testing.assert_array_equal(cat.data, arr)
        np.testing.assert_array_equal(cat["Mbound"], arr["Mbound"])
        assert cat.a == 0.5
        assert cat.z == 1.0  # fallback 1/a - 1 (no 'z' key written)
        assert cat.a_v == 0.4973
        assert cat.box == 1024.0
        assert cat.step == 100
        assert cat.path == stem + ".bin"

    cat = cats[0]
    pos = cat.positions()
    assert pos.shape == (50, 3) and pos.dtype == np.float64
    np.testing.assert_array_equal(pos[:, 1], arr["y"].astype(np.float64))
    vel = cat.velocities_at_av()
    assert vel.shape == (50, 3) and vel.dtype == np.float64
    expected = np.stack([arr[c].astype(np.float64) for c in ("vx", "vy", "vz")], axis=1) * (0.5 / 0.4973)
    np.testing.assert_allclose(vel, expected, rtol=1e-12)

    # missing a / a_v -> ValueError
    bare = Catalogue(arr, {"a": "0.5"}, stem)
    try:
        bare.velocities_at_av()
    except ValueError as exc:
        assert "a_v" in str(exc)
    else:
        raise AssertionError("velocities_at_av without a_v should raise")

    # explicit z wins over the fallback
    assert Catalogue(arr, {"a": "0.5", "z": "0.9  # comment"}, stem).z == 0.9
    assert Catalogue(arr, {}, stem).a is None and Catalogue(arr, {}, stem).z is None

    # missing sidecar -> FileNotFoundError naming it
    os.remove(stem + ".meta")
    try:
        read_catalogue(stem)
    except FileNotFoundError as exc:
        assert stem + ".meta" in str(exc)
    else:
        raise AssertionError("missing .meta should raise FileNotFoundError")


def test_mismatched_columns_raises(workdir: str) -> None:
    arr = synthetic_catalogue(10)
    stem = os.path.join(workdir, "CatshortV.0100.0002")
    bad = " ".join(NAMES[:7] + ["Mvir"] + NAMES[8:])
    write_binary(stem, arr, columns=bad)
    try:
        read_catalogue(stem)
    except ValueError as exc:
        assert "columns" in str(exc)
    else:
        raise AssertionError("mismatched columns should raise ValueError")

    write_binary(stem, arr, extra={"dtypes": DTYPE_CODES.replace("i4", "f4", 1)})
    try:
        read_catalogue(stem)
    except ValueError as exc:
        assert "dtypes" in str(exc)
    else:
        raise AssertionError("mismatched dtypes should raise ValueError")

    write_binary(stem, arr, extra={"record_bytes": "100"})
    try:
        read_catalogue(stem)
    except ValueError as exc:
        assert "record_bytes" in str(exc)
    else:
        raise AssertionError("wrong record_bytes should raise ValueError")


def test_truncated_binary_raises(workdir: str) -> None:
    arr = synthetic_catalogue(10)
    stem = os.path.join(workdir, "CatshortV.0100.0003")
    write_binary(stem, arr)
    with open(stem + ".bin", "r+b") as fh:
        fh.truncate(10 * RECORD_BYTES - 7)
    try:
        read_catalogue(stem)
    except ValueError as exc:
        assert "multiple" in str(exc)
    else:
        raise AssertionError("truncated .bin should raise ValueError")

    # whole records missing but rows says 10 -> ValueError too
    with open(stem + ".bin", "r+b") as fh:
        fh.truncate(8 * RECORD_BYTES)
    try:
        read_catalogue(stem)
    except ValueError as exc:
        assert "rows" in str(exc)
    else:
        raise AssertionError("rows mismatch should raise ValueError")


def test_text_reader(workdir: str) -> None:
    path = os.path.join(workdir, "CatshortV.0181.0001.DAT")
    lines = write_text(path)
    cat = read_catalogue(path)
    assert isinstance(cat, Catalogue)
    assert cat.meta["format"] == "text"
    assert cat.meta["product"] == "CatshortV"
    assert cat.meta["header"].startswith("ic2pm: 2LPTic ingest")
    assert cat.a == 1.0
    assert cat.z == 0.0
    assert abs(cat.a_v - (1.0 - 0.0054 / 2)) < 1e-9
    assert meta_float(cat.meta, "da") == 0.0054
    assert cat.meta["rows"] == str(len(cat)) and cat.meta["realization"] == "1"
    assert cat.step == 181
    assert meta_float(cat.meta, "N_row") == 1024
    assert meta_float(cat.meta, "N_grid") == 2048
    assert meta_float(cat.meta, "Omega_m") == 0.3
    assert meta_float(cat.meta, "Omega_L") == 0.7
    assert meta_float(cat.meta, "hubble") == 0.7
    assert meta_float(cat.meta, "buffer_width") == 15.375
    assert meta_float(cat.meta, "particle_mass") == 8.322e10
    assert meta_float(cat.meta, "overdensity") == 338.3
    assert cat.box is None

    assert len(cat) == 3 and cat.data.dtype == CATALOGUE_DTYPE
    for i, line in enumerate(lines):
        tokens = [float(t) for t in line.split()]
        assert len(tokens) == 24
        for name, tok in zip(NAMES, tokens):
            if name in ("row", "distinct"):
                assert cat[name][i] == int(round(tok))
            else:
                assert cat[name][i] == np.float32(tok)
    assert cat["x"][0] == np.float32(394.8253)
    assert list(cat["row"]) == [1, 2, 3]

    # stem resolution falls back to the .DAT when no .bin/.meta exist
    assert len(read_catalogue(path[: -len(".DAT")])) == 3

    # header-only file -> empty catalogue with parsed meta
    empty = os.path.join(workdir, "CatshortV.0181.0002.DAT")
    with open(empty, "w") as fh:
        fh.write(HEADER_8)
    ecat = read_catalogue_text(empty)
    assert len(ecat) == 0 and ecat.data.dtype == CATALOGUE_DTYPE and ecat.a == 1.0

    # a short row -> ValueError
    short = os.path.join(workdir, "CatshortV.0181.0003.DAT")
    with open(short, "w") as fh:
        fh.write(HEADER_8)
        fh.write(lines[0])
        fh.write(" ".join(lines[1].split()[:23]) + "\n")
    try:
        read_catalogue(short)
    except ValueError:
        pass
    else:
        raise AssertionError("row with 23 columns should raise ValueError")

    # truncated header -> ValueError
    stub = os.path.join(workdir, "CatshortV.0181.0004.DAT")
    with open(stub, "w") as fh:
        fh.write("".join(HEADER_8.splitlines(True)[:5]))
    try:
        read_catalogue(stub)
    except ValueError:
        pass
    else:
        raise AssertionError("short header should raise ValueError")


def test_text_binary_equivalence(workdir: str) -> None:
    path = os.path.join(workdir, "CatshortV.0181.0001.DAT")
    write_text(path)
    tcat = read_catalogue_text(path)
    stem = os.path.join(workdir, "CatshortV.0181.0001")
    write_binary(stem, tcat.data, extra={"a": tcat.meta["a"], "a_v": tcat.meta["a_v"]})
    bcat = read_catalogue(stem)  # stem now prefers the binary
    assert bcat.path == stem + ".bin"
    np.testing.assert_array_equal(bcat.data, tcat.data)
    assert bcat.a == tcat.a and bcat.a_v == tcat.a_v
    np.testing.assert_array_equal(bcat.velocities_at_av(), tcat.velocities_at_av())


def test_read_meta(workdir: str) -> None:
    p = os.path.join(workdir, "x.meta")
    with open(p, "w") as fh:
        fh.write("a = 0.5\n\nno equals here\n  key with spaces  =  v = w \nempty =\nd = 1.0D-5 comment\n")
    m = read_meta(p)
    assert m == {"a": "0.5", "key with spaces": "v = w", "empty": "", "d": "1.0D-5 comment"}
    assert meta_float(m, "a") == 0.5
    assert meta_float(m, "d") == 1e-5
    assert meta_float(m, "empty") is None
    assert meta_float(m, "missing", 7.0) == 7.0
    try:
        meta_float(m, "key with spaces")
    except ValueError:
        pass
    else:
        raise AssertionError("non-numeric value should raise ValueError")
    assert glamio.__version__ == "0.1.0"


TESTS = [
    test_binary_roundtrip,
    test_mismatched_columns_raises,
    test_truncated_binary_raises,
    test_text_reader,
    test_text_binary_equivalence,
    test_read_meta,
]

def test_byte_swapped_and_missing_layout(workdir: str) -> None:
    """A byte-swapped .bin fails the row self-check; a sidecar without 'columns' is refused."""
    arr = synthetic_catalogue(12)
    stem = os.path.join(workdir, "CatshortV.0100.0001")
    write_binary(stem, arr)
    cat = read_catalogue(stem, mmap=True)                 # memory-mapped read is equivalent
    assert np.array_equal(cat.data, arr) and len(cat) == 12
    swapped = arr.byteswap()                               # same bytes reinterpreted big-endian
    swapped.view(CATALOGUE_DTYPE).tofile(stem + ".bin")
    try:
        read_catalogue(stem)
    except ValueError as exc:
        assert "'row' column is not 1..12" in str(exc)
    else:
        raise AssertionError("byte-swapped .bin was accepted")
    arr.tofile(stem + ".bin")
    lines = [l for l in open(stem + ".meta") if not l.startswith("columns")]
    open(stem + ".meta", "w").writelines(lines)
    try:
        read_catalogue(stem)
    except ValueError as exc:
        assert "no 'columns' line" in str(exc)
    else:
        raise AssertionError("sidecar without 'columns' was accepted")


TESTS = TESTS + [test_byte_swapped_and_missing_layout] if "TESTS" in globals() else None  # noqa: F821


if __name__ == "__main__":
    failures = 0
    for fn in TESTS:
        d = _mkdir()
        try:
            fn(d)
            print(f"OK    {fn.__name__}")
        except Exception as exc:  # noqa: BLE001
            failures += 1
            print(f"FAIL  {fn.__name__}: {exc!r}")
        finally:
            shutil.rmtree(d, ignore_errors=True)
    print(f"{len(TESTS) - failures} passed, {failures} failed")
    sys.exit(1 if failures else 0)
