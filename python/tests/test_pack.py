"""Tests for glamio.pack (pack, verify, unpack, prune) and for reading packed runs.

A synthetic run directory ``Run1/`` with two full outputs (one with an empty
catalogue) and one station is packed; every product read back through
:func:`glamio.product_path` from the packed files must equal the native read,
unpacking must rebuild every file byte for byte, and verify / prune must
refuse damaged files.  Needs h5py (skipped without it).  Runs under pytest or
as a plain script; fixtures go to a fresh ``tempfile.mkdtemp`` directory under
``$GLAMIO_TEST_TMPDIR`` (or the default temp dir) and are removed afterwards.
"""

from __future__ import annotations

import filecmp
import hashlib
import os
import shutil
import subprocess
import sys
import tempfile
import textwrap

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import glamio  # noqa: E402
from glamio import CATALOGUE_DTYPE, FRHALO_DTYPE  # noqa: E402

try:
    import h5py
    from glamio import pack as gpack
except ImportError:  # no h5py: nothing to test
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


def _tmpdir() -> str:
    base = os.environ.get("GLAMIO_TEST_TMPDIR") or None
    if base:
        os.makedirs(base, exist_ok=True)
    return tempfile.mkdtemp(prefix="glamio-pack-", dir=base)


def _write(path: str, text: str) -> str:
    with open(path, "w") as fh:
        fh.write(text)
    return path


def _raises(exc, fn, *args, **kwargs):
    try:
        fn(*args, **kwargs)
    except exc:
        return True
    raise AssertionError(f"{fn.__name__} did not raise {exc.__name__}")


def _common(step: int) -> str:
    a = A[step]
    return (f"step =   {step}\nrealization =      1\na =   {a:.8E}\nz =   {1 / a - 1:.8E}\n"
            f"a_v =   {a - 0.004:.8E}   (velocity epoch of this output)\nBox =    {BOX:.4f}\nunits = Mpc/h\n"
            "N_grid =     72\nN_row =     36\nN_p =        46656\nOmega_m =   0.308900\nOmega_L =   0.691100\n"
            f"hubble =   0.677400\nMG_flag =   0\nMG_model =   1\nfR0 = none\nheader = {HEADER}\n")


def _tag(rundir: str) -> str:
    """The 4-digit box tag of a ``Run<box>`` fixture directory."""
    return f"{int(os.path.basename(rundir)[3:]):04d}"


def _lens(rundir: str, step: int, ns: int = 2, nl: int = 4) -> None:
    stem = os.path.join(rundir, f"LensPlanes.{step:04d}.{_tag(rundir)}")
    (1.0 + 0.01 * np.random.default_rng(step).standard_normal((3 * ns, nl, nl))).astype("<f4").tofile(stem + ".bin")
    _write(stem + ".meta", "product = LensPlanes\nformat = raw little-endian float32, no record markers\n"
                           + _common(step) + f"N_slab =      {ns}\ndepth = {BOX / ns:.4f}\nN_lens = {nl}\nN_planes = {3 * ns}\n")


def _modes(rundir: str, step: int, product: str, n: int) -> None:
    rng = np.random.default_rng(sum(map(ord, product)) + step)
    shape = (2 * n + 1, 2 * n + 1, n + 1)
    stem = os.path.join(rundir, f"{product}.{step:04d}.{_tag(rundir)}")
    (rng.standard_normal(shape) + 1j * rng.standard_normal(shape)).astype("<c8").tofile(stem + ".bin")
    _write(stem + ".meta", f"product = {product}\nformat = raw little-endian complex64\nn =    {n}\nk_c = 0.5\nN_dep = 72\n"
                           + _common(step) + "D_code(a) = 0.55\n")


def _subsample(rundir: str, step: int, n_sel: int = 11) -> None:
    q = np.random.default_rng(step).integers(-32768, 32767, size=(6, n_sel)).astype("<i2")
    stem = os.path.join(rundir, f"Subsample.{step:04d}.{_tag(rundir)}")
    q.tofile(stem + ".bin")
    _write(stem + ".meta", "product = Subsample\nformat = raw little-endian int16\n" + _common(step)
           + f"M =        100\nSALT = 0x5EED00DE6C0A2026\nn_sel = {n_sel}\ndv_quant =   0.2500   (km/s)\n")


def _catalogue(rundir: str, step: int, rows: int) -> None:
    os.makedirs(os.path.join(rundir, "CATALOGS"), exist_ok=True)
    stem = os.path.join(rundir, "CATALOGS", f"CatshortV.{step:04d}.{_tag(rundir)}")
    data = np.zeros(rows, dtype=CATALOGUE_DTYPE)
    data["row"] = np.arange(1, rows + 1)
    data["Mbound"] = np.linspace(1e12, 3e12, rows)
    data["x"] = np.linspace(1.0, 140.0, rows)
    data.tofile(stem + ".bin")
    names = list(CATALOGUE_DTYPE.names)
    codes = " ".join("i4" if CATALOGUE_DTYPE[c].kind == "i" else "f4" for c in names)
    _write(stem + ".meta", f"product = CatshortV\nformat = raw little-endian records\ncolumns = {' '.join(names)}\n"
                           f"dtypes = {codes}\nrows =    {rows}\nrecord_bytes = 96\ntext_catalogue = none\n" + _common(step))
    _write(os.path.join(rundir, "CATALOGS", f"outputB.{step:04d}.dat"), f"     Boundaries(Mpch) =  0.0  {BOX}\n time = {step}\n")
    stem = os.path.join(rundir, f"fR_halo.{step:04d}.{_tag(rundir)}")
    np.array([(-1e-6 * (i + 1), -2e-6) for i in range(rows)], dtype=FRHALO_DTYPE).tofile(stem + ".bin")
    _write(stem + ".meta", "product = fR_halo\n"
                           f"catalogue = CATALOGS/CatshortV.{step:04d}.{_tag(rundir)}.DAT   (text catalogue, when CatalogueFormat is text or both)\n"
                           f"catalogue_bin = CATALOGS/CatshortV.{step:04d}.{_tag(rundir)}.bin   (binary catalogue)\nrows = {rows}\n"
                           "fR_bg(a_old) =  -2.7E-06\n" + _common(step))


def _power(rundir: str, step: int) -> None:
    a = A[step]
    for kind, nrow in (("R", 4), ("log", 3)):
        rows = "".join(f"  {i + 1:4d}  {0.05 * (i + 1):.5f} {18.0 * (i + 1):.4E}  {1000.0 / (i + 1):.5E}\n" for i in range(nrow))
        _write(os.path.join(rundir, f"PowerDM.{kind}.{step:04d}.{_tag(rundir)}.dat"),
               f"{HEADER:45s}\n Aexpn = {a:.5f} Step= {step:4d} Redshift= {1 / a - 1:8.3f} Ngrid=   72 kNyq=     1.51\n"
               f"   bin  k/Mpch    N          Preal\n{rows}")


def _run_fixture(d: str, box: int = 1) -> str:
    rundir = os.path.join(d, f"Run{box}")
    os.makedirs(rundir)
    _write(os.path.join(rundir, "OutputSchedule.dat"), SCHEDULE)
    _write(os.path.join(rundir, "BDM.config"), "MaxMemory = 900\nCatalogueFormat = binary\n")
    _write(os.path.join(rundir, "main.log"), "step 70\r\nstep 84\nstep 98\n")      # a CR LF line on purpose
    _write(os.path.join(rundir, "fort.17"), "")                                   # an empty file
    for step in (70, 98):
        _catalogue(rundir, step, rows=3 if step == 70 else 0)
    for step in (70,):
        _subsample(rundir, step)
        _modes(rundir, step, "deltam", 3)
        for c in "xyz":
            _modes(rundir, step, f"mom_{c}", 2)
    for step in (70, 84, 98):
        _lens(rundir, step)
        _power(rundir, step)
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


def _skip() -> bool:
    if h5py is not None:
        return False
    if pytest is not None:
        pytest.skip("h5py not installed")
    return True


def _compare_reads(native: str, packed: str, box: int = 1) -> None:
    """Every product read through product_path from ``packed`` equals the native read."""
    readers = {"LensPlanes": (glamio.read_lensplanes, "planes"), "deltam": (glamio.read_modes, "data"),
               "mom_x": (glamio.read_modes, "data"), "mom_y": (glamio.read_modes, "data"), "mom_z": (glamio.read_modes, "data"),
               "Subsample": (glamio.read_subsample, "q"), "fR_halo": (glamio.read_frhalo, "data"),
               "catalogue": (glamio.read_catalogue, "data")}
    for step in (70, 84, 98):
        found_n, found_p = glamio.products_at(native, step, box), glamio.products_at(packed, step, box)
        assert set(found_n) == set(found_p), (step, found_n, found_p)
        for name in found_n:
            pn, pp = glamio.product_path(native, name, step, box), glamio.product_path(packed, name, step, box)
            assert pp == found_p[name] and ".h5/" in pp and ".h5/" not in pn
            if name.startswith("PowerDM"):
                a, b = glamio.read_powerdm(pn, box=BOX), glamio.read_powerdm(pp, box=BOX)
                assert np.array_equal(a.data, b.data) and a.columns == b.columns and a.binning == b.binning
                assert (a.a, a.z, a.step, a.ngrid, a.k_nyq, a.header) == (b.a, b.z, b.step, b.ngrid, b.k_nyq, b.header)
                continue
            reader, attr = readers[name]
            for mmap in (False, True):
                a, b = reader(pn, mmap=mmap), reader(pp, mmap=mmap)
                xa, xb = getattr(a, attr), getattr(b, attr)
                assert xa.dtype == xb.dtype and xa.shape == xb.shape and np.array_equal(xa, xb), (step, name, mmap)
                assert a.meta == b.meta and a.z == b.z and a.box == b.box
                if mmap and xb.size:
                    assert isinstance(xb.base if not isinstance(xb, np.memmap) else xb, np.memmap), (name, type(xb))
            if name == "fR_halo":
                ca, cb = glamio.read_frhalo(pn).catalogue(), glamio.read_frhalo(pp).catalogue()
                assert np.array_equal(ca.data, cb.data) and ".h5/CatshortV" in cb.path
    sn, sp = glamio.read_schedule(glamio.product_path(native, "schedule", 0)), glamio.read_schedule(glamio.product_path(packed, "schedule", 0))
    assert np.array_equal(sn.rows, sp.rows) and sn.comments == sp.comments and sp.path.endswith(f"RunInfo.{box:04d}.h5/schedule")
    assert glamio.read_meta(glamio.product_path(packed, "deltam", 70, box)) == glamio.read_meta(glamio.product_path(native, "deltam", 70, box) + ".meta")


# ----------------------------------------------------------------- tests
def test_pack_roundtrip():
    if _skip():
        return
    d = _tmpdir()
    try:
        rundir = _run_fixture(d)
        out = os.path.join(d, "packed", "Run1")
        written = gpack.pack(rundir, out, log=_quiet)
        assert sorted(os.listdir(out)) == ["Output.0070.0001.h5", "Output.0084.0001.h5", "Output.0098.0001.h5",
                                           "RunInfo.0001.h5", "packed.md5"]
        assert len(written) == 4
        _compare_reads(rundir, out)

        with h5py.File(os.path.join(out, "Output.0070.0001.h5"), "r") as h:          # plain h5py view of the layout
            assert h.attrs["format"] == "glamio-pack" and h.attrs["step"] == 70 and h.attrs["class"] == 1
            assert abs(h.attrs["a"] - A[70]) < 1e-12 and h.attrs["N_p"] == 46656 and h.attrs["fR0"] == "none"
            assert h.attrs["header"] == HEADER and h.attrs["run"] == os.path.basename(d) and h.attrs["box_dir"] == "Run1"
            assert np.array_equal(h["CatshortV/data"]["Mbound"], np.linspace(1e12, 3e12, 3).astype("<f4"))
            assert h["CatshortV"].attrs["rows"] == "3" and h["deltam"].attrs["n"] == "3"
            assert h["deltam/data"].shape == (7, 7, 4) and h["deltam/data"].chunks is None
            assert h["LensPlanes/data"].shape == (6, 4, 4) and h["Subsample/data"].shape == (6, 11)
            assert list(h["PowerDM/R/table"].dtype.names) == ["bin", "k", "nmodes", "pk"] and h["PowerDM/R/table"].shape == (4,)
            assert h["CatshortV/outputB"][()].tobytes().startswith(b"     Boundaries")
            assert len(h["manifest"]) == 2 * 8 + 1 + 2         # 8 sidecar pairs, outputB, two PowerDM tables
        with h5py.File(os.path.join(out, "RunInfo.0001.h5"), "r") as h:
            assert list(h["schedule/table"]["step"]) == [70, 84, 98] and list(h.attrs["steps"]) == [70, 84, 98]
            assert h["logs/main.log"][()].tobytes() == b"step 70\r\nstep 84\nstep 98\n" and h["logs/fort.17"].shape == (0,)
            assert h["BDM.config"][()].tobytes().startswith(b"MaxMemory")
        with open(os.path.join(out, "packed.md5")) as fh:
            for line in fh:
                md5, name = line.split()
                assert hashlib.md5(open(os.path.join(out, name), "rb").read()).hexdigest() == md5

        for path in written:
            problems, stats = gpack.verify_file(path, against=rundir)
            assert problems == [] and stats["originals_checked"] == stats["entries"]

        rebuilt = os.path.join(d, "rebuilt")
        assert gpack.unpack([out], rebuilt, log=_quiet) == len(_tree(rundir))
        assert _same_tree(rundir, rebuilt)
        assert _raises(gpack.PackError, gpack.unpack, [out], rebuilt, log=_quiet)      # never overwrites

        assert gpack.pack(rundir, out, log=_quiet) == written                          # existing files verified and kept
        assert gpack.main(["verify", out, "--against", rundir]) == 0
    finally:
        shutil.rmtree(d, ignore_errors=True)


def test_pack_refuses_unknown_files():
    if _skip():
        return
    d = _tmpdir()
    try:
        rundir = _run_fixture(d)
        for bad, why in ((os.path.join(rundir, "Projection.0070.0001.dat"), "step-tagged"),
                         (os.path.join(rundir, "PMcrd.0070.DAT"), "PM particle"),
                         (os.path.join(rundir, "CATALOGS", "CatshortV.0070.0001.DAT"), "unknown file in CATALOGS"),
                         (os.path.join(rundir, "LensPlanes.0112.0001.bin"), "needs both"),
                         (os.path.join(rundir, "Subsample.0070.0002.bin"), "box tag")):
            _write(bad, "x")
            try:
                gpack.pack(rundir, os.path.join(d, "out"), log=_quiet)
                raise AssertionError(f"{bad} was not refused")
            except gpack.PackError as exc:
                assert why in str(exc), (why, str(exc))
            os.remove(bad)
        _lens(rundir, 112)                                                             # a step outside the schedule
        assert _raises(gpack.PackError, gpack.pack, rundir, os.path.join(d, "out"), log=_quiet)
        os.remove(os.path.join(rundir, "LensPlanes.0112.0001.bin"))
        os.remove(os.path.join(rundir, "LensPlanes.0112.0001.meta"))
        os.makedirs(os.path.join(rundir, "extra"))
        assert _raises(gpack.PackError, gpack.pack, rundir, os.path.join(d, "out"), log=_quiet)
        assert not os.path.exists(os.path.join(d, "out"))                              # refused before writing anything
    finally:
        shutil.rmtree(d, ignore_errors=True)


def test_verify_detects_damage():
    if _skip():
        return
    d = _tmpdir()
    try:
        rundir = _run_fixture(d)
        out = os.path.join(d, "packed")
        gpack.pack(rundir, out, box=1, log=_quiet)
        f = os.path.join(out, "Output.0070.0001.h5")
        with h5py.File(f, "r+") as h:
            x = h["Subsample/data"]
            x[0, 0] = x[0, 0] + 1
        problems, _ = gpack.verify_file(f)
        assert len(problems) == 1 and "Subsample/data" in problems[0]
        assert gpack.main(["verify", f]) == 1
        assert _raises(gpack.PackError, gpack.pack, rundir, out, box=1, log=_quiet)    # a damaged existing file is not kept
        g = os.path.join(out, "Output.0084.0001.h5")
        with h5py.File(g, "r+") as h:
            h["LensPlanes"].attrs["meta_text"] = np.bytes_(b"product = LensPlanes\n")
        problems, _ = gpack.verify_file(g)
        assert len(problems) == 1 and "LensPlanes@meta_text" in problems[0]
        with open(os.path.join(rundir, "PowerDM.R.0098.0001.dat"), "a") as fh:     # an original changed after packing
            fh.write("\n")
        problems, _ = gpack.verify_file(os.path.join(out, "Output.0098.0001.h5"), against=rundir)
        assert len(problems) == 1 and "differs" in problems[0]
    finally:
        shutil.rmtree(d, ignore_errors=True)


def test_prune_in_place():
    if _skip():
        return
    d = _tmpdir()
    try:
        rundir = _run_fixture(d)
        copy = os.path.join(d, "copy")
        shutil.copytree(rundir, copy)
        assert _raises(gpack.PackError, gpack.prune, rundir, yes=True, log=_quiet)      # nothing packed yet
        gpack.pack(rundir, log=_quiet)                                                 # in place, beside the originals
        assert ".h5/" not in glamio.product_path(rundir, "deltam", 70)                # natives preferred while present
        assert gpack.prune(rundir, log=_quiet) == 0 and len(_tree(rundir)) == len(_tree(copy)) + 5
        with open(os.path.join(rundir, "main.log"), "a") as fh:                        # an original that differs blocks it all
            fh.write("late line\n")
        assert _raises(gpack.PackError, gpack.prune, rundir, yes=True, log=_quiet)
        assert len(_tree(rundir)) == len(_tree(copy)) + 5
        shutil.copy(os.path.join(copy, "main.log"), os.path.join(rundir, "main.log"))
        assert gpack.prune(rundir, yes=True, log=_quiet) == len(_tree(copy))
        assert sorted(os.listdir(rundir)) == ["Output.0070.0001.h5", "Output.0084.0001.h5", "Output.0098.0001.h5",
                                              "RunInfo.0001.h5", "packed.md5"]
        _compare_reads(copy, rundir)                                                   # now from the packed files
        gpack.unpack([rundir], os.path.join(d, "rebuilt"), log=_quiet)
        assert _same_tree(copy, os.path.join(d, "rebuilt"))
    finally:
        shutil.rmtree(d, ignore_errors=True)


def test_box0_and_path_rules():
    if _skip():
        return
    d = _tmpdir()
    try:
        rundir = _run_fixture(d, box=0)
        out = os.path.join(d, "packed", "Run0")
        gpack.pack(rundir, out, log=_quiet)
        _compare_reads(rundir, out, box=0)
        assert glamio.product_path(out, "schedule", 0).endswith("RunInfo.0000.h5/schedule")    # default box=1 still finds it
        for name, obj in (("main.log", "logs/main.log"), ("BDM.config", "BDM.config"), ("OutputSchedule.dat", "schedule")):
            pn, pp = glamio.run_file(rundir, name, 0), glamio.run_file(out, name)
            assert pn == os.path.join(rundir, name) and pp.endswith("RunInfo.0000.h5/" + obj)
            assert glamio.read_text(pn) == glamio.read_text(pp)
        assert glamio.run_file(out, "nothing.log") == os.path.join(out, "nothing.log")
        assert len(glamio.read_schedule(glamio.product_path(out, "schedule", 0))) == 3
        mixed = os.path.join(d, "mixed", "Run0")                       # an incomplete native pair falls back to the packed copy
        shutil.copytree(rundir, mixed)
        shutil.copy(os.path.join(out, "Output.0070.0000.h5"), mixed)
        os.remove(os.path.join(mixed, "fR_halo.0070.0000.bin"))
        p = glamio.product_path(mixed, "fR_halo", 70, 0)
        assert p.endswith("Output.0070.0000.h5/fR_halo") and glamio.products_at(mixed, 70, 0)["fR_halo"] == p
        assert len(glamio.read_frhalo(p)) == 3
        try:
            glamio.read_modes(os.path.join(out, "Output.0999.0000.h5", "deltam"))
            raise AssertionError("a missing packed file was not reported")
        except FileNotFoundError as exc:
            assert "packed file not found" in str(exc) and "Output.0999.0000.h5" in str(exc)
        pp = glamio.product_path(out, "deltam", 70, 0)
        m, a = glamio.read_modes(pp, mmap=True).data, glamio.read_modes(pp).data
        f = os.path.abspath(os.path.join(out, "Output.0070.0000.h5"))
        with h5py.File(f, "r") as h:
            off = h["deltam/data"].id.get_offset()
        assert isinstance(m, np.memmap) and m.filename == f and m.offset == off and m.offset % 4096 == 0
        assert not isinstance(a, np.memmap) and not isinstance(a.base, np.memmap) and np.array_equal(a, m)
        os.makedirs(os.path.join(rundir, "CATALOGS", "CatshortV.0070.0000"))     # a directory named like the stem
        assert len(glamio.read_catalogue(os.path.join(rundir, "CATALOGS", "CatshortV.0070.0000"))) == 3
    finally:
        shutil.rmtree(d, ignore_errors=True)


def test_native_reads_need_no_h5py():
    if _skip():
        return
    d = _tmpdir()
    try:
        rundir = _run_fixture(d)
        prelude = f"import sys\nsys.path.insert(0, {os.path.dirname(os.path.dirname(os.path.abspath(__file__)))!r})\n"
        reads = textwrap.dedent(f"""
            import glamio
            r = {rundir!r}
            R = {{"LensPlanes": glamio.read_lensplanes, "deltam": glamio.read_modes, "mom_x": glamio.read_modes,
                  "mom_y": glamio.read_modes, "mom_z": glamio.read_modes, "Subsample": glamio.read_subsample,
                  "fR_halo": glamio.read_frhalo, "catalogue": glamio.read_catalogue,
                  "PowerDM.R": glamio.read_powerdm, "PowerDM.log": glamio.read_powerdm}}
            for step in (70, 84, 98):
                found = glamio.products_at(r, step)
                assert found and all(".h5" not in p for p in found.values()), found
                for name, path in found.items():
                    R[name](path)
                assert ".h5" not in glamio.product_path(r, "deltam", 84)      # absent product: native path
            glamio.read_frhalo(glamio.product_path(r, "fR_halo", 70)).catalogue()
            glamio.read_schedule(glamio.product_path(r, "schedule", 0))
            """)
        native_only = prelude + reads + "assert 'h5py' not in sys.modules, 'h5py was imported'\n"
        subprocess.run([sys.executable, "-I", "-c", native_only], check=True)
        gpack.pack(rundir, log=_quiet)                                 # packed files beside the natives
        blocked = prelude + "sys.modules['h5py'] = None\n" + reads
        subprocess.run([sys.executable, "-I", "-c", blocked], check=True)
    finally:
        shutil.rmtree(d, ignore_errors=True)


def test_resume_status_and_prune_rules():
    if _skip():
        return
    d = _tmpdir()
    try:
        rundir = _run_fixture(d)
        assert gpack.status(rundir)[0] == "native"
        gpack.pack(rundir, log=_quiet)
        assert gpack.status(rundir)[0] == "packed"
        path = os.path.join(rundir, "deltam.0070.0001.bin")            # an original changed after packing
        raw = open(path, "rb").read()
        with open(path, "r+b") as fh:
            fh.write(bytes([raw[0] ^ 1]))
        try:
            gpack.pack(rundir, log=_quiet)
            raise AssertionError("a stale packed file was kept")
        except gpack.PackError as exc:
            assert "does not match the originals" in str(exc)
        with open(path, "wb") as fh:
            fh.write(raw)
        late = _write(os.path.join(rundir, "check_run.log"), "0 FAIL\n")   # a file added after packing
        assert gpack.status(rundir)[0] == "incomplete"
        assert _raises(gpack.PackError, gpack.pack, rundir, log=_quiet)
        try:
            gpack.prune(rundir, yes=True, log=_quiet)
            raise AssertionError("prune deleted with an unpacked file present")
        except gpack.PackError as exc:
            assert "in no packed file" in str(exc)
        os.remove(late)
        part = _run_fixture(os.path.join(d, "part"))                   # packed one step only: prune refuses
        gpack.pack(part, steps=[70], log=_quiet)
        n_before = len(_tree(part))
        assert _raises(gpack.PackError, gpack.prune, part, yes=True, log=_quiet) and len(_tree(part)) == n_before
        assert gpack.status(part)[0] == "incomplete"
        for f in ("LensPlanes.0084.0001.meta", "deltam.0070.0001.meta"):   # an interrupted prune
            os.remove(os.path.join(rundir, f))
        assert gpack.status(rundir)[0] == "partially-pruned"
        assert _raises(gpack.PackError, gpack.pack, rundir, log=_quiet)
        _write(os.path.join(rundir, "Output.0084.0001.h5.partial"), "")
        assert gpack.status(rundir)[0] == "busy" and _raises(gpack.PackError, gpack.prune, rundir, yes=True, log=_quiet)
        os.remove(os.path.join(rundir, "Output.0084.0001.h5.partial"))
        assert gpack.prune(rundir, yes=True, log=_quiet) > 0 and gpack.status(rundir)[0] == "pruned"
        assert gpack.main(["verify", rundir]) == 0
    finally:
        shutil.rmtree(d, ignore_errors=True)


def test_unpack_ownership_and_verify_gaps():
    if _skip():
        return
    d = _tmpdir()
    try:
        rundir = _run_fixture(d)
        _write(os.path.join(rundir, "slurm@m8380.out"), "a log whose name holds an @\n")
        out = os.path.join(d, "packed")
        gpack.pack(rundir, out, box=1, log=_quiet)
        u = os.path.join(d, "u")
        os.makedirs(u)
        precious = _write(os.path.join(u, "main.log.partial"), "precious")
        gpack.unpack([os.path.join(out, "RunInfo.0001.h5")], u, log=_quiet)
        assert open(precious).read() == "precious" and not [n for n in os.listdir(u) if n.startswith(".")]
        umask = os.umask(0)
        os.umask(umask)
        assert os.stat(os.path.join(u, "main.log")).st_mode & 0o777 == 0o666 & ~umask
        assert open(os.path.join(u, "slurm@m8380.out")).read().startswith("a log")
        _write(os.path.join(rundir, "packed.md5"), "my notes\n")         # a user's file with the reserved name
        try:
            gpack.pack(rundir, os.path.join(d, "o2"), log=_quiet)
            raise AssertionError("a foreign packed.md5 was accepted")
        except gpack.PackError as exc:
            assert "not a checksum list" in str(exc)
        os.remove(os.path.join(rundir, "packed.md5"))
        empty = os.path.join(d, "empty")
        os.makedirs(empty)
        assert gpack.main(["verify", empty]) == 1                       # nothing to verify is a failure
        moved = os.path.join(d, "moved.h5")
        os.rename(os.path.join(out, "Output.0084.0001.h5"), moved)
        assert gpack.main(["verify", out]) == 1                         # a missing output file is noticed
        os.rename(moved, os.path.join(out, "Output.0084.0001.h5"))
        assert gpack.main(["verify", out]) == 0
    finally:
        shutil.rmtree(d, ignore_errors=True)


def test_format_details_and_small_blocks():
    if _skip():
        return
    d = _tmpdir()
    old = gpack.BLOCK_BYTES
    try:
        rundir = _run_fixture(d)
        gpack.BLOCK_BYTES = 40                                          # rows split across blocks, multi-block datasets
        out = os.path.join(d, "packed")
        gpack.pack(rundir, out, box=1, log=_quiet)
        gpack.unpack([out], os.path.join(d, "rebuilt"), log=_quiet)
        assert _same_tree(rundir, os.path.join(d, "rebuilt"))
        gpack.BLOCK_BYTES = old
        f = os.path.join(out, "Output.0070.0001.h5")
        with h5py.File(f, "r") as h:
            offsets = [h[g]["data"].id.get_offset() for g in h if isinstance(h[g], h5py.Group) and "data" in h[g] and h[g]["data"].size]
            assert len(offsets) == 8 and all(o % 4096 == 0 for o in offsets)
            assert "Packed runs" in h.attrs["layout"] and len(h.attrs["glamio_code_sha256"]) == 64
        newer = os.path.join(d, "newer", "Output.0070.0001.h5")
        os.makedirs(os.path.dirname(newer))
        shutil.copy(f, newer)
        with h5py.File(newer, "r+") as h:
            h.attrs["format_version"] = np.int32(99)
        assert _raises(ValueError, glamio.read_lensplanes, newer + "/LensPlanes")
        assert "format_version 99" in gpack.verify_file(newer)[0][0]
    finally:
        gpack.BLOCK_BYTES = old
        shutil.rmtree(d, ignore_errors=True)


def test_round2_guards():
    if _skip():
        return
    d = _tmpdir()
    try:
        rundir = _run_fixture(d, box=0)
        other = _run_fixture(os.path.join(d, "other"), box=1)
        gpack.pack(other, os.path.join(d, "other_packed"), box=1, log=_quiet)
        gpack.pack(rundir, log=_quiet)
        stray = os.path.join(rundir, "Output.0001.0001.h5")             # a packed file of another box sorts first
        shutil.copy(os.path.join(d, "other_packed", "Output.0070.0001.h5"), stray)
        assert gpack.status(rundir)[0] == "inconsistent"
        n_before = len(_tree(rundir))
        try:
            gpack.prune(rundir, yes=True, log=_quiet)
            raise AssertionError("prune ran with a packed file of another box present")
        except gpack.PackError as exc:
            assert "boxes [0, 1]" in str(exc)
        assert len(_tree(rundir)) == n_before
        os.remove(stray)
        extra = os.path.join(rundir, "Output.0099.0000.h5")             # an Output file RunInfo does not list
        shutil.copy(os.path.join(rundir, "Output.0084.0000.h5"), extra)
        assert gpack.status(rundir)[0] == "inconsistent"
        assert _raises(gpack.PackError, gpack.prune, rundir, yes=True, log=_quiet) and len(_tree(rundir)) == n_before   # stray out, extra in
        os.remove(extra)
        assert gpack.status(rundir)[0] == "packed"
        os.remove(os.path.join(rundir, "Output.0084.0000.h5"))          # a stale .partial from another writer is kept
        busy = _write(os.path.join(rundir, "Output.0084.0000.h5.partial"), "another job's file")
        try:
            gpack.pack(rundir, log=_quiet)
            raise AssertionError("pack wrote over another writer's .partial")
        except gpack.PackError as exc:
            assert "exists" in str(exc)
        assert open(busy).read() == "another job's file"
        os.remove(busy)
        gpack.pack(rundir, log=_quiet)                                   # rewrites the missing output, keeps the rest
        assert gpack.status(rundir)[0] == "packed"
        out = os.path.join(d, "out")                                      # a user's packed.md5 in --out is not overwritten
        os.makedirs(out)
        _write(os.path.join(out, "packed.md5"), "my precious notes\n")
        assert _raises(gpack.PackError, gpack.pack, rundir, out, log=_quiet)
        assert open(os.path.join(out, "packed.md5")).read() == "my precious notes\n"
        with open(os.path.join(rundir, "RunInfo.0000.h5"), "r+b") as fh:   # an unreadable RunInfo: reported, no traceback
            fh.truncate(100)
        assert gpack.main(["verify", rundir]) == 1
    finally:
        shutil.rmtree(d, ignore_errors=True)


TESTS = [test_pack_roundtrip, test_pack_refuses_unknown_files, test_verify_detects_damage, test_prune_in_place,
         test_box0_and_path_rules, test_native_reads_need_no_h5py, test_resume_status_and_prune_rules,
         test_unpack_ownership_and_verify_gaps, test_format_details_and_small_blocks, test_round2_guards]

if __name__ == "__main__":
    if h5py is None:
        print("SKIP: h5py not installed")
        sys.exit(0)
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
