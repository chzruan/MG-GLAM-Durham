"""Readers for GLAM BDM halo catalogues (``CATALOGS/Catshort?.<step>.<box>.*``).

Two on-disk forms carry the same 24 columns in the same row order:

* Legacy text ``Catshort?.<step>.<box>.DAT``: an 8-line header (run header
  string; ``A`` and ``Step``; ``I``, ``Nrow``, ``Ngrid``; cosmology and buffer
  width; number of radial bins; smallest particle mass; overdensity limit; the
  column-caption line) followed by one row per halo of 24 whitespace-separated
  numbers, written with the Fortran format
  ``(3f11.4,3x,3f10.2,1p,2g12.4,g12.5,26g12.4)``.
* Binary ``Catshort?.<step>.<box>.bin`` + sidecar ``.meta``: headerless,
  little-endian, row-major, one 96-byte record per halo made of 24 four-byte
  fields (:data:`CATALOGUE_DTYPE`).  The sidecar is a ``key = value`` text file
  (see :mod:`glamio.meta`) describing the layout (``columns``, ``dtypes``,
  ``rows``, ``record_bytes``) and the run (``a``, ``a_v``, ``Box``, ...).

Both readers return a :class:`Catalogue` whose ``data`` is a structured numpy
array with :data:`CATALOGUE_DTYPE`.  The letter after ``Catshort`` (``V``,
``W``, ``M`` or ``A``) names the overdensity definition; ``V`` (virial) is the
usual one.

Velocity convention
-------------------
The ``vx, vy, vz`` columns are the velocities *as written* by the finder: the
leapfrog velocities, which live at the half step ``a_v = a - da/2``, converted
to km/s with the factor ``1/a`` of the output epoch ``a``.  The conversion
proper to their own epoch uses ``1/a_v`` instead, so

    v(a_v) = v_written * a / a_v,

which is what :meth:`Catalogue.velocities_at_av` returns (``a`` and ``a_v`` are
taken from the meta).  Positions are comoving Mpc/h in ``[0, Box)``; masses are
Msun/h; ``Rvir`` and ``RadRMS`` are kpc/h.
"""

from __future__ import annotations

import os
import re
import warnings

import numpy as np

from .meta import meta_float, read_meta

#: Record layout of the binary catalogue (and of the arrays returned by both readers).
CATALOGUE_DTYPE = np.dtype(
    [
        ("x", "<f4"), ("y", "<f4"), ("z", "<f4"),
        ("vx", "<f4"), ("vy", "<f4"), ("vz", "<f4"),
        ("Mbound", "<f4"), ("Mtot", "<f4"), ("Rvir", "<f4"),
        ("Vrms", "<f4"), ("Vcirc", "<f4"), ("row", "<i4"),
        ("Cvir", "<f4"), ("Npart", "<f4"), ("distinct", "<i4"),
        ("Xoff", "<f4"), ("virial_ratio", "<f4"), ("lambda", "<f4"),
        ("RadRMS", "<f4"), ("b_a", "<f4"), ("c_a", "<f4"),
        ("ax_x", "<f4"), ("ax_y", "<f4"), ("ax_z", "<f4"),
    ]
)

#: Bytes per halo record in the binary catalogue.
RECORD_BYTES = 96
assert CATALOGUE_DTYPE.itemsize == RECORD_BYTES

#: Physical unit (or meaning) of every column.
COLUMN_UNITS: dict[str, str] = {
    "x": "Mpc/h (comoving, in [0, Box))",
    "y": "Mpc/h (comoving, in [0, Box))",
    "z": "Mpc/h (comoving, in [0, Box))",
    "vx": "km/s (peculiar, as written: converted with 1/a)",
    "vy": "km/s (peculiar, as written: converted with 1/a)",
    "vz": "km/s (peculiar, as written: converted with 1/a)",
    "Mbound": "Msun/h (bound virial mass)",
    "Mtot": "Msun/h (total mass inside Rvir)",
    "Rvir": "kpc/h",
    "Vrms": "km/s",
    "Vcirc": "km/s (maximum circular velocity)",
    "row": "1-based row index in this file",
    "Cvir": "concentration",
    "Npart": "number of particles",
    "distinct": "0 = distinct halo",
    "Xoff": "centre offset",
    "virial_ratio": "2K/|Ep| - 1",
    "lambda": "spin parameter",
    "RadRMS": "kpc/h",
    "b_a": "axis ratio b/a",
    "c_a": "axis ratio c/a",
    "ax_x": "major-axis direction",
    "ax_y": "major-axis direction",
    "ax_z": "major-axis direction",
}

#: Number of header lines in the legacy text catalogue.
TEXT_HEADER_LINES = 8

_INT_COLUMNS = ("row", "distinct")
_NCOL = len(CATALOGUE_DTYPE.names)
_PRODUCT_RE = re.compile(r"^(Catshort[A-Za-z])\.")


class Catalogue:
    """A halo catalogue: structured array ``data``, sidecar/header ``meta``, source ``path``.

    ``cat['Mbound']`` returns the column view, ``len(cat)`` the number of
    haloes.  The scalar properties are parsed from ``meta`` on access and are
    ``None`` when the corresponding key is missing.
    """

    def __init__(self, data: np.ndarray, meta: dict[str, str], path: str):
        self.data = data
        self.meta = meta
        self.path = path

    def __len__(self) -> int:
        return len(self.data)

    def __getitem__(self, name: str) -> np.ndarray:
        return self.data[name]

    def __repr__(self) -> str:
        return f"Catalogue({len(self)} haloes, a={self.a}, path={self.path!r})"

    @property
    def columns(self) -> list[str]:
        """Column names in file order."""
        return list(self.data.dtype.names)

    @property
    def a(self) -> float | None:
        """Expansion factor of the output (meta key ``a``)."""
        return meta_float(self.meta, "a")

    @property
    def z(self) -> float | None:
        """Redshift: meta key ``z``, else ``1/a - 1``, else ``None``."""
        z = meta_float(self.meta, "z")
        if z is None:
            a = self.a
            if a is not None and a != 0.0:
                z = 1.0 / a - 1.0
        return z

    @property
    def a_v(self) -> float | None:
        """Velocity epoch ``a_v = a - da/2`` (meta key ``a_v``)."""
        return meta_float(self.meta, "a_v")

    @property
    def box(self) -> float | None:
        """Box side in Mpc/h (meta key ``Box``)."""
        return meta_float(self.meta, "Box")

    @property
    def step(self) -> int | None:
        """Time-step number of the output (meta key ``step``)."""
        value = meta_float(self.meta, "step")
        return None if value is None else int(round(value))

    def positions(self) -> np.ndarray:
        """Return ``(N, 3)`` float64 comoving positions ``[x, y, z]`` in Mpc/h."""
        return np.stack([self.data[c].astype(np.float64) for c in ("x", "y", "z")], axis=1)

    def velocities_at_av(self) -> np.ndarray:
        """Return ``(N, 3)`` float64 velocities at their own epoch ``a_v``, in km/s.

        The written velocities were converted with ``1/a``; the leapfrog
        velocities live at ``a_v = a - da/2``, so ``v(a_v) = v_written * a / a_v``.
        Raises ``ValueError`` when ``a`` or ``a_v`` is not in the meta.
        """
        a, a_v = self.a, self.a_v
        if a is None or a_v is None:
            raise ValueError(
                "velocities_at_av() needs both 'a' and 'a_v' in the catalogue meta; "
                f"got a={a!r}, a_v={a_v!r} for {self.path}"
            )
        if a_v == 0.0:
            raise ValueError(f"velocities_at_av(): a_v is zero in the meta of {self.path}")
        factor = a / a_v
        vel = np.stack([self.data[c].astype(np.float64) for c in ("vx", "vy", "vz")], axis=1)
        vel *= factor
        return vel


def read_catalogue(path: str | os.PathLike) -> Catalogue:
    """Read a BDM catalogue from its ``.bin``/``.meta`` pair or its legacy ``.DAT``.

    ``path`` may be the ``.bin``, the ``.meta``, the ``.DAT`` or the stem without
    extension.  A ``.DAT`` (or ``.dat``) path goes to :func:`read_catalogue_text`.
    A bare stem is read as binary when ``<stem>.bin`` or ``<stem>.meta`` exists,
    as text when only ``<stem>.DAT`` exists, and as binary otherwise (which
    raises ``FileNotFoundError`` naming the missing file).
    """
    p = os.fspath(path)
    if p.endswith((".DAT", ".dat")):
        return read_catalogue_text(p)
    if p.endswith(".bin"):
        stem = p[: -len(".bin")]
    elif p.endswith(".meta"):
        stem = p[: -len(".meta")]
    else:
        stem = p
        if not (os.path.isfile(stem + ".bin") or os.path.isfile(stem + ".meta")):
            for ext in (".DAT", ".dat"):
                if os.path.isfile(stem + ext):
                    return read_catalogue_text(stem + ext)
    return _read_catalogue_binary(stem)


def _check_layout(meta: dict[str, str], meta_path: str) -> None:
    """Raise ``ValueError`` if the sidecar's ``columns``/``dtypes``/``record_bytes`` differ from CATALOGUE_DTYPE."""
    if "columns" in meta:
        names = meta["columns"].split()
        if names != list(CATALOGUE_DTYPE.names):
            raise ValueError(
                f"{meta_path}: 'columns' does not match glamio.CATALOGUE_DTYPE "
                f"(file: {names}; expected: {list(CATALOGUE_DTYPE.names)})"
            )
    if "dtypes" in meta:
        codes = meta["dtypes"].split()
        try:
            got = [np.dtype(c) for c in codes]
        except TypeError as exc:
            raise ValueError(f"{meta_path}: unparsable 'dtypes' entry: {exc}") from exc
        expected = [CATALOGUE_DTYPE[i] for i in range(_NCOL)]
        if len(got) != _NCOL or any(g != e for g, e in zip(got, expected)):
            raise ValueError(
                f"{meta_path}: 'dtypes' does not match glamio.CATALOGUE_DTYPE "
                f"(file: {codes}; expected: {[e.str for e in expected]})"
            )
    rb = meta_float(meta, "record_bytes")
    if rb is not None and int(rb) != RECORD_BYTES:
        raise ValueError(f"{meta_path}: record_bytes = {int(rb)}, expected {RECORD_BYTES}")


def _read_catalogue_binary(stem: str) -> Catalogue:
    """Read ``<stem>.bin`` with the layout checked against ``<stem>.meta``."""
    bin_path, meta_path = stem + ".bin", stem + ".meta"
    for p in (bin_path, meta_path):
        if not os.path.isfile(p):
            raise FileNotFoundError(f"binary catalogue needs both .bin and .meta; missing {p}")
    meta = read_meta(meta_path)
    _check_layout(meta, meta_path)

    size = os.path.getsize(bin_path)
    if size % RECORD_BYTES:
        raise ValueError(
            f"{bin_path}: size {size} bytes is not a multiple of the {RECORD_BYTES}-byte record"
        )
    rows = meta_float(meta, "rows")
    if rows is not None and size != int(rows) * RECORD_BYTES:
        raise ValueError(
            f"{bin_path}: size {size} bytes, but {meta_path} says rows = {int(rows)} "
            f"({int(rows) * RECORD_BYTES} bytes)"
        )
    data = np.fromfile(bin_path, dtype=CATALOGUE_DTYPE)
    return Catalogue(data, meta, bin_path)


def _grab(pattern: str, line: str, meta: dict[str, str], *keys: str) -> None:
    """Store the regex groups of ``pattern`` found in ``line`` under ``keys`` (no-op if no match)."""
    m = re.search(pattern, line)
    if m is None:
        return
    for key, value in zip(keys, m.groups()):
        meta[key] = value


def _parse_text_header(lines: list[str], path: str) -> dict[str, str]:
    """Turn the 8 header lines of a text catalogue into a meta dict (strings).

    Keys: ``header``, ``a``, ``step_da``, ``a_v = a - step_da/2``, ``z = 1/a - 1``,
    ``step`` (the ``I`` value), ``N_row``, ``N_grid``, ``Omega_m``, ``Omega_L``,
    ``hubble``, ``buffer_width``, ``particle_mass``, ``overdensity``, ``product``
    (the ``Catshort?`` label from the file name) and ``format = text``.
    """
    num = r"([-+]?[0-9.]+(?:[eEdD][-+]?[0-9]+)?)"
    meta: dict[str, str] = {"header": lines[0].strip()}
    _grab(rf"\bA\s*=\s*{num}\s+Step\s*=\s*{num}", lines[1], meta, "a", "step_da")
    _grab(rf"\bI\s*=\s*{num}\s+Nrow\s*=\s*{num}\s+Ngrid\s*=\s*{num}", lines[2], meta, "step", "N_row", "N_grid")
    _grab(rf"Omega_0\s*=\s*{num}", lines[3], meta, "Omega_m")
    _grab(rf"Omega_L\s*=\s*{num}", lines[3], meta, "Omega_L")
    _grab(rf"hubble\s*=\s*{num}", lines[3], meta, "hubble")
    _grab(rf"buffer width\s*\(Mpch\)\s*=\s*{num}", lines[3], meta, "buffer_width")
    _grab(rf"Mass of smallest particle\s*\(Msunh\)\s*=\s*{num}", lines[5], meta, "particle_mass")
    _grab(rf"Overdensity limit\s*=\s*{num}", lines[6], meta, "overdensity")

    a, da = meta_float(meta, "a"), meta_float(meta, "step_da")
    if a is not None:
        meta["z"] = f"{1.0 / a - 1.0:.6f}"
        if da is not None:
            meta["a_v"] = f"{a - 0.5 * da:.6f}"
    m = _PRODUCT_RE.match(os.path.basename(path))
    if m is not None:
        meta["product"] = m.group(1)
    meta["format"] = "text"
    return meta


def read_catalogue_text(path: str | os.PathLike) -> Catalogue:
    """Read a legacy text catalogue ``Catshort?.<step>.<box>.DAT``.

    The 8-line header becomes ``meta`` (see :func:`_parse_text_header`); the rows
    are parsed with ``np.loadtxt`` and cast to :data:`CATALOGUE_DTYPE` (float
    columns to ``<f4``, ``row`` and ``distinct`` to ``<i4``).  Raises ``ValueError``
    when the header is shorter than 8 lines or a row does not have exactly 24
    columns.  A header-only file gives a zero-length array.
    """
    p = os.fspath(path)
    with open(p, encoding="utf-8", errors="replace") as fh:
        header = [fh.readline() for _ in range(TEXT_HEADER_LINES)]
    if any(line == "" for line in header):
        raise ValueError(f"{p}: text catalogue header shorter than {TEXT_HEADER_LINES} lines")
    meta = _parse_text_header(header, p)

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message=".*input contained no data.*", category=UserWarning)
        raw = np.loadtxt(p, skiprows=TEXT_HEADER_LINES, ndmin=2)
    if raw.size == 0:
        return Catalogue(np.empty(0, dtype=CATALOGUE_DTYPE), meta, p)
    if raw.ndim != 2 or raw.shape[1] != _NCOL:
        ncol = raw.shape[1] if raw.ndim == 2 else raw.ndim
        raise ValueError(f"{p}: expected {_NCOL} columns per row, found {ncol}")

    data = np.empty(raw.shape[0], dtype=CATALOGUE_DTYPE)
    for i, name in enumerate(CATALOGUE_DTYPE.names):
        if name in _INT_COLUMNS:
            data[name] = np.rint(raw[:, i]).astype("<i4")
        else:
            data[name] = raw[:, i]
    return Catalogue(data, meta, p)
