"""Tests for the non-catalogue glamio readers with synthetic files.

schedule, lens planes, mode cubes, subsample (and the splitmix64 selection),
fR_halo, PowerDM and the PM snapshot.  Runs under pytest or as a plain script
(``python tests/test_readers.py``); fixtures go to a fresh ``tempfile.mkdtemp``
directory under ``$GLAMIO_TEST_TMPDIR`` (or the default temp dir) and are
removed afterwards.  The precise_z suite cross-checks the same readers against
GLAM's real output (``test/precise_z/asserts.py``).
"""

from __future__ import annotations

import os
import shutil
import struct
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import glamio  # noqa: E402
from glamio import CATALOGUE_DTYPE  # noqa: E402
from glamio.snapshot import HEADER_BYTES  # noqa: E402
from glamio.subsample import NQ, SALT_DEFAULT, keep_mask, select_ids, splitmix64  # noqa: E402

try:
    import pytest
except ImportError:  # plain-script mode
    pytest = None


def _tmpdir() -> str:
    base = os.environ.get("GLAMIO_TEST_TMPDIR") or None
    if base:
        os.makedirs(base, exist_ok=True)
    return tempfile.mkdtemp(prefix="glamio-readers-", dir=base)


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


# ----------------------------------------------------------------- schedule
def test_schedule():
    d = _tmpdir()
    try:
        p = _write(os.path.join(d, "OutputSchedule.dat"),
                   "# MG-GLAM output schedule\n"
                   "# chi: comoving distance [Mpc/h], flat matter+Lambda, Om= 0.30890 OmL= 0.69110, c/H0 = 2997.92458 Mpc/h\n"
                   "#  step  a  z  chi exact class\n"
                   "    70  2.85714298E-01  2.49999976E+00    4045.9986     1     1\n"
                   "    84  4.76190507E-01  1.09999990E+00    2464.3611     0     2\n"
                   "\n"
                   "    98  6.61375642E-01  5.12000084E-01    1345.2087     1     1\n")
        s = glamio.read_schedule(p)
        assert len(s) == 3 and list(s.steps) == [70, 84, 98]
        assert s.omega_m == 0.3089 and s.omega_l == 0.6911 and s.c_over_h0 == 2997.92458
        assert abs(s.row(84)["chi"] - 2464.3611) < 1e-9 and s.row(84)["class"] == 2 and s.row(84)["exact"] == 0
        assert list(s.primary()["step"]) == [70, 98] and list(s.stations()["step"]) == [84]
        assert set(s.by_step()) == {70, 84, 98} and s["a"][0] == float("2.85714298E-01")
        assert _raises(KeyError, s.row, 85)
        assert s.row_at_z(1.1)["step"] == 84 and _raises(KeyError, s.row_at_z, 1.2)
        _write(os.path.join(d, "bad.dat"), "70 0.5 1.0 100.0 1\n")
        assert _raises(ValueError, glamio.read_schedule, os.path.join(d, "bad.dat"))
        _write(os.path.join(d, "order.dat"), "84 0.5 1.0 100.0 1 1\n70 0.3 2.3 200.0 1 1\n")
        assert _raises(ValueError, glamio.read_schedule, os.path.join(d, "order.dat"))
        _write(os.path.join(d, "class0.dat"), "70 0.5 1.0 100.0 1 0\n")
        assert _raises(ValueError, glamio.read_schedule, os.path.join(d, "class0.dat"))
    finally:
        shutil.rmtree(d, ignore_errors=True)


# ----------------------------------------------------------------- lens planes
def _lens_fixture(d, ns=2, nl=3, box=150.0):
    planes = np.zeros((3 * ns, nl, nl), dtype="<f4")      # [plane, j, i] = plane*100 + j*10 + i
    for p in range(3 * ns):
        for j in range(nl):
            for i in range(nl):
                planes[p, j, i] = p * 100 + j * 10 + i
    stem = os.path.join(d, "LensPlanes.0070.0001")
    planes.tofile(stem + ".bin")
    _write(stem + ".meta",
           "product = LensPlanes\nformat = raw little-endian float32\nstep =     70\nrealization = 1\n"
           f"a =   2.85714298E-01\nz = 2.5\na_v = 0.2817\nBox =     {box:.4f}\nN_grid = 72\nN_slab = {ns}\n"
           f"depth = {box / ns:.4f}\nN_lens = {nl}\nN_planes = {3 * ns}\nheader = N=0036x0072\n")
    return stem, planes


def test_lensplanes():
    d = _tmpdir()
    try:
        stem, planes = _lens_fixture(d)
        for arg, mmap in ((stem, False), (stem + ".bin", False), (stem + ".meta", True)):
            L = glamio.read_lensplanes(arg, mmap=mmap)
            assert np.array_equal(L.planes, planes) and L.n_slab == 2 and L.n_lens == 3 and L.n_planes == 6
        L = glamio.read_lensplanes(stem)
        assert L.a == float("2.85714298E-01") and L.box == 150.0 and L.depth == 75.0 and L.pixel_size == 50.0 and L.step == 70
        # plane('x', 1) is [i, j]: value i + 10 j of the first plane
        px = L.plane("x", 1)
        assert px[2, 1] == 12.0 and px[1, 2] == 21.0 and np.array_equal(px, planes[0].T)
        assert np.array_equal(L.plane("z", 2), planes[5].T) and L.slabs("y").shape == (2, 3, 3) and L.slabs(1)[0, 0, 0] == 200.0
        assert np.allclose(L.axis_means(), [np.mean(planes[0:2]), np.mean(planes[2:4]), np.mean(planes[4:6])])
        pm = L.plane_means()
        assert pm.shape == (3, 2) and np.allclose(pm, [[np.mean(planes[p]) for p in (0, 1)], [np.mean(planes[p]) for p in (2, 3)], [np.mean(planes[p]) for p in (4, 5)]])
        assert _raises(ValueError, L.plane, True, 1) and _raises(ValueError, L.slabs, 1.0)
        assert np.array_equal(L.pixel_edges(), [0, 50, 100, 150]) and np.array_equal(L.slab_edges(), [0, 75, 150])
        assert _raises(IndexError, L.plane, "x", 0) and _raises(ValueError, L.plane, "w", 1)
        # wrong size, wrong product, missing key
        with open(stem + ".bin", "ab") as fh:
            fh.write(b"\0\0\0\0")
        assert _raises(ValueError, glamio.read_lensplanes, stem)
        planes.tofile(stem + ".bin")
        meta = open(stem + ".meta").read()
        _write(stem + ".meta", meta.replace("product = LensPlanes", "product = deltam"))
        assert _raises(ValueError, glamio.read_lensplanes, stem)
        _write(stem + ".meta", meta.replace("N_lens = 3\n", ""))
        assert _raises(ValueError, glamio.read_lensplanes, stem)
        os.remove(stem + ".meta")
        assert _raises(FileNotFoundError, glamio.read_lensplanes, stem)
    finally:
        shutil.rmtree(d, ignore_errors=True)


# ----------------------------------------------------------------- mode cubes
def _brute_power(M, nbins):
    """Reference: mean of Box^3 |delta_k|^2 over the full cube -n..n^3 in bins nint(|m|), through mode()."""
    n = M.n
    sums = np.zeros(nbins)
    counts = np.zeros(nbins)
    for ix in range(-n, n + 1):
        for iy in range(-n, n + 1):
            for iz in range(-n, n + 1):
                b = int(np.rint(np.sqrt(ix * ix + iy * iy + iz * iz)))
                if 1 <= b <= nbins:
                    sums[b - 1] += abs(M.mode(ix, iy, iz)) ** 2 * M.box**3
                    counts[b - 1] += 1
    with np.errstate(invalid="ignore"):
        return sums / counts, counts


def test_modes():
    d = _tmpdir()
    try:
        n, box = 2, 100.0
        rng = np.random.default_rng(3)
        data = (rng.standard_normal((2 * n + 1, 2 * n + 1, n + 1)) + 1j * rng.standard_normal((2 * n + 1, 2 * n + 1, n + 1))).astype("<c8")
        data[n, n, 0] = 0
        stem = os.path.join(d, "deltam.0084.0001")
        data.tofile(stem + ".bin")
        _write(stem + ".meta", f"product = deltam\nn = {n}\nk_c = 0.5\nN_dep = 72\nN_grid = 72\nstep = 84\n"
                               f"a = 0.476\na_v = 0.470\nD_code(a) = 0.55\nBox = {box:.4f}\n")
        M = glamio.read_modes(stem)
        assert np.array_equal(M.data, data) and M.n == n and M.k_c == 0.5 and M.n_dep == 72 and M.d_code == 0.55
        assert abs(M.k_fund - 2 * np.pi / box) < 1e-15 and M.product == "deltam"
        assert M.mode(1, -2, 1) == complex(data[1 + n, -2 + n, 1])
        assert M.mode(1, -2, -1) == complex(np.conj(data[-1 + n, 2 + n, 1]))
        assert M.mode(0, 0, 0) == 0
        assert _raises(IndexError, M.mode, 3, 0, 0)
        kx, ky, kz = M.wavenumbers()
        assert kx.shape == data.shape and kx[0, 0, 0] == -n * M.k_fund and kz[0, 0, n] == n * M.k_fund
        assert abs(M.k_magnitude()[n + 1, n, 0] - M.k_fund) < 1e-12
        km, pk, w = M.power()
        ref, cnt = _brute_power(M, n)
        assert np.allclose(pk, ref, rtol=1e-12) and np.array_equal(w, cnt)
        assert _raises(ValueError, M.power, n + 1)
        sx, sy, sz = M.index_grids(sparse=True)
        assert sx.shape == (2 * n + 1, 1, 1) and sz.shape == (1, 1, n + 1) and np.array_equal(sx + 0 * sy + 0 * sz, M.index_grids()[0])
        assert abs(km[0] - M.k_fund * np.sum([1] * 6 + [np.sqrt(2)] * 12) / 18) < 1e-12   # bin 1: 6 axis + 12 face-diagonal modes
        M2 = glamio.read_modes(stem, mmap=True, expect_product=("deltam", "mom_x"))
        assert np.array_equal(M2.data, data)
        assert _raises(ValueError, glamio.read_modes, stem, expect_product="mom_x")
        data[:, :, :-1].tofile(stem + ".bin")            # wrong size
        assert _raises(ValueError, glamio.read_modes, stem)
        # a deltaL-style sidecar with only n and Box
        data.tofile(stem + ".bin")
        _write(stem + ".meta", f"product = deltaL\nn = {n}\nBox = {box}\nD_gen(a_ini) = 0.02\n")
        assert glamio.read_modes(stem).k_c is None and glamio.read_modes(stem).d_gen == 0.02
        _write(stem + ".meta", f"product = fR_halo\nn = {n}\nBox = {box}\n")
        assert _raises(ValueError, glamio.read_modes, stem)                       # not a mode product
        assert glamio.read_modes(stem, expect_product=None).product == "fR_halo"
        _write(stem + ".meta", f"product = deltam   (cube)\nn = {n}\n")
        M3 = glamio.read_modes(stem)
        assert M3.product == "deltam" and _raises(ValueError, M3.power) and _raises(ValueError, M3.k_magnitude)
    finally:
        shutil.rmtree(d, ignore_errors=True)


# ----------------------------------------------------------------- subsample
def test_splitmix64_and_selection():
    assert int(splitmix64(np.uint64(0))) == 0xE220A8397B1DCDAF           # reference vector of SplitMix64 for seed 0
    assert int(splitmix64(np.uint64(1))) == 0x910A2DEC89025CC1
    ids = select_ids(1000, SALT_DEFAULT, 7)
    assert ids.dtype == np.uint64 and np.all(np.diff(ids.astype(np.int64)) > 0)
    assert np.all(splitmix64(ids ^ np.uint64(SALT_DEFAULT)) % np.uint64(7) == 0)
    assert np.array_equal(select_ids(1000, SALT_DEFAULT, 7, block=64), ids)
    assert np.array_equal(select_ids(5, SALT_DEFAULT, 1), np.arange(5, dtype=np.uint64))
    assert 1000 / 7 * 0.5 < ids.size < 1000 / 7 * 1.5


def _subsample_fixture(d, m, n_p, box=150.0, dv=0.25, step=70, salt=SALT_DEFAULT, header="ic2pm: 2LPTic ingest, v at a-da/2, IDz"):
    ids = select_ids(n_p, salt, m)
    n_sel = ids.size
    q = np.zeros((6, n_sel), dtype="<i2")
    q[0] = (np.arange(n_sel) * 7919 % NQ) - 32768
    q[1] = -32768
    q[2] = 32767
    q[3] = np.arange(n_sel) - 3
    q[4] = -32767
    q[5] = 32767
    stem = os.path.join(d, f"Subsample.{step:04d}.0001")
    q.tofile(stem + ".bin")
    _write(stem + ".meta",
           f"product = Subsample\nstep = {step}\nrealization = 1\na = 0.2857\nz = 2.5\na_v = 0.2817   (velocity epoch)\n"
           f"Box = {box:.4f}\nN_grid = 72\nN_row = 36\nN_p = {n_p}\nM = {m}\nSALT = 0x{salt:016X}\n"
           f"n_sel = {n_sel}\ndx_quant = {box / NQ:.8E}   (Box/65536 Mpc/h)\ndv_quant = {dv:.4f}   (km/s)\n"
           f"n_clip_v = 0\nv_scale = 739.4   (km/s per code unit)\nheader = {header}\n")
    return stem, q, ids


def test_subsample():
    d = _tmpdir()
    try:
        stem, q, ids = _subsample_fixture(d, m=5, n_p=2000)
        S = glamio.read_subsample(stem)
        assert len(S) == S.n_sel == ids.size and S.m == 5 and S.n_p == 2000 and S.salt == SALT_DEFAULT and S.lattice_tag == "IDz" == S.lattice
        assert np.array_equal(S.q, q) and abs(S.dx_quant - 150.0 / NQ) < 1e-10 and S.dv_quant == 0.25 and S.a_v == 0.2817
        pos = S.positions()
        assert pos.shape == (ids.size, 3) and np.all(pos >= 0) and np.all(pos < 150.0)
        assert pos[0, 1] == 0.5 * 150.0 / NQ and abs(pos[0, 2] - (65535.5 * 150.0 / NQ)) < 1e-12
        vel = S.velocities()
        assert vel[0, 0] == -0.75 and vel[0, 1] == -32767 * 0.25 and vel[0, 2] == 32767 * 0.25
        assert np.array_equal(S.ids(), ids)
        assert np.array_equal(glamio.read_subsample(stem + ".meta", mmap=True).q, q)
        meta = open(stem + ".meta").read()
        _write(stem + ".meta", meta.replace("N_p = 2000", "N_p = 4000"))          # the rule then keeps about twice as many IDs
        assert _raises(ValueError, glamio.read_subsample(stem).ids)
        _write(stem + ".meta", meta.replace(f"n_sel = {ids.size}", f"n_sel = {ids.size + 1}"))
        assert _raises(ValueError, glamio.read_subsample, stem)
        stem1, q1, ids1 = _subsample_fixture(d, m=1, n_p=12, step=84, header="N=0036x0072L=0150zi049da=8.00E-04f0.040")
        S1 = glamio.read_subsample(stem1)
        assert np.array_equal(S1.ids(), np.arange(12, dtype=np.uint64)) and len(S1) == 12 and S1.lattice_tag is None
        assert _raises(ValueError, S1.ids, max_ids=11) and S1.ids(max_ids=12).size == 12
        # a non-default salt from the sidecar selects a different set; M = 0 is rejected
        stem2, q2, ids2 = _subsample_fixture(d, m=3, n_p=3000, step=98, salt=0x0123456789ABCDEF)
        S2 = glamio.read_subsample(stem2)
        assert S2.salt == 0x0123456789ABCDEF and np.array_equal(S2.ids(), ids2) and not np.array_equal(ids2, select_ids(3000, SALT_DEFAULT, 3))
        assert _raises(ValueError, select_ids, 10, SALT_DEFAULT, 0) and _raises(ValueError, keep_mask, [1, 2], SALT_DEFAULT, 0)
        # an empty subsample (n_sel = 0) reads with and without mmap
        stem3 = os.path.join(d, "Subsample.0099.0001")
        open(stem3 + ".bin", "wb").close()
        _write(stem3 + ".meta", "product = Subsample\nn_sel = 0\nBox = 150.0\ndv_quant = 0.25\nN_p = 5\nM = 7\n")
        for mm in (False, True):
            S3 = glamio.read_subsample(stem3, mmap=mm)
            assert len(S3) == 0 and S3.positions().shape == (0, 3) and S3.velocities().shape == (0, 3)
    finally:
        shutil.rmtree(d, ignore_errors=True)


# ----------------------------------------------------------------- fR_halo
def _catalogue_fixture(rundir, step, rows):
    os.makedirs(os.path.join(rundir, "CATALOGS"), exist_ok=True)
    stem = os.path.join(rundir, "CATALOGS", f"CatshortV.{step:04d}.0001")
    data = np.zeros(rows, dtype=CATALOGUE_DTYPE)
    data["row"] = np.arange(1, rows + 1)
    data["Mbound"] = 1e12
    data.tofile(stem + ".bin")
    names = list(CATALOGUE_DTYPE.names)
    codes = " ".join("i4" if CATALOGUE_DTYPE[n].kind == "i" else "f4" for n in names)
    _write(stem + ".meta", "product = CatshortV\nformat = raw little-endian\n"
                           f"columns = {' '.join(names)}\ndtypes = {codes}\nrows = {rows}\nrecord_bytes = 96\n"
                           f"a = 0.4762\na_v = 0.4704\nz = 1.1\nstep = {step}\nBox = 150.0\ntext_catalogue = none\n")
    return stem


def test_frhalo():
    d = _tmpdir()
    try:
        rundir = os.path.join(d, "Run1")
        os.makedirs(rundir)
        pairs = np.array([(-1e-6, -1.1e-6), (-2e-6, -2e-6), (-3e-7, -4e-7)], dtype="<f4")
        stem = os.path.join(rundir, "fR_halo.0084.0001")
        pairs.tofile(stem + ".bin")
        _write(stem + ".meta",
               "product = fR_halo\nformat = raw little-endian float32 pairs\n"
               "catalogue = CATALOGS/CatshortV.0084.0001.DAT   (text catalogue, when CatalogueFormat is text or both)\n"
               "catalogue_bin = CATALOGS/CatshortV.0084.0001.bin   (binary catalogue)\nrows = 3\n"
               "fR0 = -1.00000000E-05\nn = 1\na_old = 0.4646\nfR_bg(a_old) = -2.76490598E-06\na = 0.4762\nz = 1.1\n"
               "step = 84\nrealization = 1\nBox = 150.0\n")
        F = glamio.read_frhalo(stem)
        assert len(F) == F.rows == 3 and np.array_equal(F["fR_centre"], pairs[:, 0]) and np.array_equal(F["fR_mean"], pairs[:, 1])
        assert F.fR0 == -1e-5 and F.n_fr == 1 and F.a_old == 0.4646 and F.fR_bg_old == float("-2.76490598E-06") and F.step == 84
        paths = F.catalogue_paths()
        assert paths["binary"] == os.path.join(rundir, "CATALOGS", "CatshortV.0084.0001.bin") and paths["text"].endswith(".DAT")
        assert _raises(FileNotFoundError, F.catalogue) and _raises(ValueError, F.catalogue, "bogus")
        _catalogue_fixture(rundir, 84, 3)
        cat = F.catalogue()
        assert len(cat) == 3 and cat.path.endswith(".bin")
        _catalogue_fixture(rundir, 84, 4)
        assert _raises(ValueError, F.catalogue)
        pairs[:2].tofile(stem + ".bin")
        assert _raises(ValueError, glamio.read_frhalo, stem)
        # a box without haloes: empty .bin and rows = 0, readable with and without mmap
        stem0 = os.path.join(rundir, "fR_halo.0070.0001")
        open(stem0 + ".bin", "wb").close()
        _write(stem0 + ".meta", "product = fR_halo\nrows = 0\nfR0 = -1e-5\nstep = 70\n")
        for mm in (False, True):
            F0 = glamio.read_frhalo(stem0, mmap=mm)
            assert len(F0) == 0 and F0.data.dtype == glamio.FRHALO_DTYPE
    finally:
        shutil.rmtree(d, ignore_errors=True)


# ----------------------------------------------------------------- PowerDM
def test_powerdm():
    d = _tmpdir()
    try:
        p = _write(os.path.join(d, "PowerDM.R.0070.0001.dat"),
                   "N=0036x0072L=0150zi049da=8.00E-04f0.040      \n"
                   " Aexpn = 0.28571 Step=  70 Redshift=    2.500 Ngrid=   72 kNyq=     1.51\n"
                   "   bin  k/Mpch    N          Preal\n"
                   "     1  0.05345 1.8000E+01  1.37139E+03\n"
                   "     2  0.09344 6.2000E+01  7.31680E+02\n")
        P = glamio.read_powerdm(p)
        assert len(P) == 2 and P.columns == ["bin", "k", "nmodes", "pk"] and P.binning == "linear"
        assert P.a == 0.28571 and P.z == 2.5 and P.step == 70 and P.ngrid == 72 and P.k_nyq == 1.51
        assert P.header.startswith("N=0036x0072") and P.k[1] == 0.09344 and P.nmodes[0] == 18 and P.pk[0] == 1371.39
        w2 = P.cic_window_jing_sq()
        assert np.all((0 < w2) & (w2 <= 1)) and np.allclose(P.pk_deconvolved() * w2, P.pk)
        assert _raises(KeyError, P.__getitem__, "p0") and P.n_grid == 72 and P.box is None
        Pb = glamio.read_powerdm(p, box=150.0)
        assert abs(Pb.k_nyq - np.pi * 72 / 150.0) < 1e-12 and Pb.box == 150.0
        hb = _write(os.path.join(d, "PowerDM.R.0071.0001.dat"), open(p).read().replace("N=0036x0072L=0150zi049da=8.00E-04f0.040      ", "LCDM bin 3 test header"))
        Ph = glamio.read_powerdm(hb)
        assert Ph.header == "LCDM bin 3 test header" and len(Ph) == 2 and Ph.columns == P.columns
        z = _write(os.path.join(d, "PowerDM.Z.0250.0001.V300.d2.dat"),
                   "hdr\n Aexpn = 0.2857 Step=  70 Redshift=    2.500 Ngrid=   36 kNyq=     0.75\n"
                   "   bin  k/Mpch    N          Preal        P0_redshift  P2_redshift\n"
                   "     1  0.05345 1.8000E+01  1.37139E+03  1.5E+03  2.0E+02\n")
        Z = glamio.read_powerdm(z)
        assert Z.columns[-2:] == ["p0", "p2"] and Z["p2"][0] == 200.0 and Z.binning == "redshift"
        assert _raises(ValueError, Z.cic_window_jing_sq) and glamio.read_powerdm(z, box=150.0).k_nyq == 0.75
        _write(os.path.join(d, "bad.dat"), "hdr\n   bin  k/Mpch    N          Preal\n     1  0.05 1.8E+01  1.3E+03\n")
        assert _raises(ValueError, glamio.read_powerdm, os.path.join(d, "bad.dat"))
        _write(os.path.join(d, "bad2.dat"), "hdr\n Aexpn = 0.2857 Step=  70 Redshift=    2.500 Ngrid=   36 kNyq=     0.75\n"
                                            "   bin  k/Mpch    N          Preal\n     1  0.05 1.8E+01\n")
        assert _raises(ValueError, glamio.read_powerdm, os.path.join(d, "bad2.dat"))
    finally:
        shutil.rmtree(d, ignore_errors=True)


# ----------------------------------------------------------------- PM snapshot
HEADER_VALUES = dict(aexpn=0.5, aexp0=0.02, amplt=0.8, astep=0.02, istep=70, partw=56.888889, tintg=0.0, ekin=1.5,
                     ekin1=2.5, ekin2=3.5, au0=0.0, aeu0=4.5, nrow=2, ngrid=8, nrealization=3, nseed=2026, om=0.3, oml=0.7,
                     hubble=0.7)


def _pack_header(byteorder, header, n, box, **overrides):
    """The 529-byte PMcrd record packed with struct (independent of the reader's dtype)."""
    v = dict(HEADER_VALUES, **overrides)
    extras = [0.0] * 100
    extras[99] = box
    return struct.pack(byteorder + "45s4fi7f4i3fq100f", header.encode().ljust(45),
                       v["aexpn"], v["aexp0"], v["amplt"], v["astep"], v["istep"],
                       v["partw"], v["tintg"], v["ekin"], v["ekin1"], v["ekin2"], v["au0"], v["aeu0"],
                       v["nrow"], v["ngrid"], v["nrealization"], v["nseed"], v["om"], v["oml"], v["hubble"], n, *extras)


def _snapshot_fixture(d, byteorder=">", n=9, page=4, rpf=2, ngrid=8, box=100.0, step=70, a=0.5, da=0.02,
                      header="synthetic PM snapshot, v at a-da/2, IDx", istep=None):
    rec = _pack_header(byteorder, header, n, box, aexpn=a, astep=da, istep=step if istep is None else istep, ngrid=ngrid)
    assert len(rec) == HEADER_BYTES
    end = "big" if byteorder == ">" else "little"
    marker = HEADER_BYTES.to_bytes(4, end)
    tag = f"{step:04d}"
    crd = os.path.join(d, f"PMcrd.{tag}.DAT")
    with open(crd, "wb") as fh:
        fh.write(marker + rec + marker)
    # particles: x = 1 + i*0.5 (grid units), y = 1 + i, z = NGRID + 1 (the periodic image of 0), v = i - 4
    i = np.arange(n, dtype=np.float32)
    comps = [1 + 0.5 * i, 1 + (i % ngrid), np.full(n, ngrid + 1, np.float32), i - 4, 2 * i, -i]   # z = NGRID + 1 wraps to 0
    n_pages = (n + page - 1) // page
    files = {}
    for ip in range(n_pages):
        fi = ip // rpf
        block = np.zeros((6, page), dtype=byteorder + "f4")
        lo, hi = ip * page, min(n, (ip + 1) * page)
        for c in range(6):
            block[c, :hi - lo] = comps[c][lo:hi]
        block[:, hi - lo:] = -99.0                              # padding the reader must ignore
        files.setdefault(fi, []).append(block.tobytes())
    for fi, blocks in files.items():
        with open(os.path.join(d, f"PMcrs{fi}.{tag}.DAT"), "wb") as fh:
            fh.write(b"".join(blocks))
    return crd, comps


def test_snapshot():
    d = _tmpdir()
    try:
        for bo in (">", "<"):
            sub = os.path.join(d, "be" if bo == ">" else "le")
            os.makedirs(sub)
            crd, comps = _snapshot_fixture(sub, byteorder=bo)
            h = glamio.read_pm_header(crd)
            assert h.byteorder == bo and h.header.startswith("synthetic") and h.nparticles == 9 and h.ngrid == 8 and h.box == 100.0
            for name, value in HEADER_VALUES.items():          # every packed field comes back (order and kinds)
                if name in ("aexpn", "astep", "istep", "ngrid"):
                    continue
                got = getattr(h, name)
                assert abs(got - value) < 1e-6 * max(1.0, abs(value)), (name, got, value)
            assert abs(h.aexpn - 0.5) < 1e-7 and h.istep == 70 and abs(h.a_v - 0.49) < 1e-7 and h.cell == 12.5 and h.ic_velocity_epoch is None
            assert abs(h.velocity_scale() - 100 * 100.0 / 8 / h.a_v) < 1e-9 and abs(h.velocity_scale("a") - 100 * 100.0 / 8 / h.aexpn) < 1e-9
            assert abs(h.particle_mass() - glamio.RHO_CRIT_GLAM * h.om * 50.0**3) < 1.0   # h.om is the float32 of 0.3
            snap = glamio.read_snapshot(crd, page=4, records_per_file=2)
            assert snap.n_pages == 3 and len(snap.files) == 2 and snap.tag == "0070" and snap.n_particles == 9
            auto = glamio.read_snapshot(crd, page=4)                 # split detected from the size of PMcrs0
            assert auto.records_per_file == 2 and len(auto.files) == 2
            x6 = snap.read(max_particles=None)
            assert all(c.dtype == np.float32 for c in x6)
            for c in range(6):
                assert np.array_equal(x6[c], np.asarray(comps[c], dtype=np.float32)), c
            pages = list(snap.iter_pages(velocities=False))
            assert [len(p[0]) for p in pages] == [4, 4, 1] and len(pages[0]) == 3
            x, y, z, vx, vy, vz = snap.read(units="physical")
            assert x.dtype == np.float64 and np.allclose(x, (comps[0] - 1) * 12.5) and np.allclose(y, (comps[1] - 1) * 12.5)
            assert np.all(z == 0.0)                                   # NGRID + 1 is the periodic image of 0
            assert np.allclose(vx, comps[3] * h.velocity_scale()) and np.allclose(vz, comps[5] * h.velocity_scale())
            _, _, _, vxa, _, _ = snap.read(units="physical", velocity_epoch="a")
            assert np.allclose(vxa, comps[3] * h.velocity_scale("a"))
            assert _raises(ValueError, snap.read, max_particles=8)
            assert _raises(ValueError, snap.iter_pages(units="grid").__next__)
        # exact multiple of the page (no padding, one file) and a single page
        for n_part, exp_pages, exp_files in ((8, 2, 1), (4, 1, 1), (12, 3, 2)):
            sub = os.path.join(d, f"n{n_part}")
            os.makedirs(sub)
            crd, comps = _snapshot_fixture(sub, n=n_part)
            snap = glamio.read_snapshot(crd, page=4)
            assert snap.n_pages == exp_pages and len(snap.files) == exp_files and snap.records_per_file == (2 if exp_files == 2 else exp_pages)
            x6 = snap.read()
            assert len(x6[0]) == n_part and np.array_equal(x6[0], np.asarray(comps[0], dtype=np.float32))
        # a 256-records-per-file split (native PMP2start ICs) is detected without options
        sub = os.path.join(d, "split3")
        os.makedirs(sub)
        crd, comps = _snapshot_fixture(sub, n=11, rpf=3)
        snap = glamio.read_snapshot(crd, page=4)
        assert snap.records_per_file == 3 and len(snap.files) == 1 and len(snap.read()[0]) == 11
        # IC files converted by ic2pm: half-step tag vs legacy sync mode
        sub = os.path.join(d, "ic")
        os.makedirs(sub)
        crd, _ = _snapshot_fixture(sub, header="ic2pm: 2LPTic ingest, v at a-da/2, IDz", istep=0, step=0)
        h = glamio.read_pm_header(crd)
        assert h.ic_velocity_epoch == "half" and abs(h.a_v - 0.49) < 1e-7
        crd, _ = _snapshot_fixture(sub, header="ic2pm: 2LPTic ingest", istep=0, step=1)
        h = glamio.read_pm_header(crd)
        assert h.ic_velocity_epoch == "sync" and abs(h.a_v - 0.5) < 1e-7
        crd, _ = _snapshot_fixture(sub, header="ic2pm: 2LPTic ingest", istep=70, step=2)
        assert glamio.read_pm_header(crd).ic_velocity_epoch is None     # a run started from a sync IC: later snapshots are half-step
        # damaged marker, missing particle file, short particle file, wrong file name
        crd, _ = _snapshot_fixture(os.path.join(d, "be"))
        with open(crd, "r+b") as fh:
            fh.write(b"\0\0\0\7")
        assert _raises(ValueError, glamio.read_pm_header, crd)
        crd, _ = _snapshot_fixture(os.path.join(d, "be"))
        os.remove(os.path.join(d, "be", "PMcrs1.0070.DAT"))
        assert _raises(FileNotFoundError, glamio.read_snapshot, crd, page=4, records_per_file=2)
        crd, _ = _snapshot_fixture(os.path.join(d, "be"))
        with open(os.path.join(d, "be", "PMcrs1.0070.DAT"), "r+b") as fh:
            fh.truncate(10)
        assert _raises(ValueError, glamio.read_snapshot, crd, page=4, records_per_file=2)
        assert _raises(ValueError, glamio.read_snapshot, os.path.join(d, "be", "PMcrs0.0070.DAT"))
    finally:
        shutil.rmtree(d, ignore_errors=True)


# ----------------------------------------------------------------- run paths
def test_run_paths():
    d = _tmpdir()
    try:
        rundir = os.path.join(d, "Run1")
        os.makedirs(os.path.join(rundir, "CATALOGS"))
        assert glamio.product_path(rundir, "deltam", 84) == os.path.join(rundir, "deltam.0084.0001")
        assert glamio.product_path(rundir, "catalogue", 84, 2) == os.path.join(rundir, "CATALOGS", "CatshortV.0084.0002")
        assert glamio.product_path(rundir, "PowerDM.R", 84).endswith("PowerDM.R.0084.0001.dat")
        assert glamio.product_path(rundir, "snapshot", 84).endswith("PMcrd.0084.DAT") and glamio.product_path(rundir, "schedule", 0).endswith("OutputSchedule.dat")
        assert _raises(ValueError, glamio.product_path, rundir, "bogus", 84)
        assert glamio.products_at(rundir, 84) == {}
        _lens_fixture(rundir)                                      # LensPlanes.0070.0001
        _write(os.path.join(rundir, "PowerDM.R.0070.0001.dat"), "x\n")
        _catalogue_fixture(rundir, 70, 2)
        found = glamio.products_at(rundir, 70)
        assert set(found) == {"LensPlanes", "PowerDM.R", "catalogue"} and found["LensPlanes"].endswith("LensPlanes.0070.0001")
    finally:
        shutil.rmtree(d, ignore_errors=True)


TESTS = [test_schedule, test_lensplanes, test_modes, test_splitmix64_and_selection, test_subsample,
         test_frhalo, test_powerdm, test_snapshot, test_run_paths]

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
