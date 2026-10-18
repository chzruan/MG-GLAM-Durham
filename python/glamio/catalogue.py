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

Column definitions follow the finder: ``Rvir`` is the spherical-overdensity
radius ``R_SO`` plus a resolution correction of up to 0.75 mesh cells,
``Mbound`` the mass of the particles inside ``R_SO`` that survive unbinding,
``Mtot`` all mass inside ``Rvir`` (the sidecar's ``definitions`` line states
this).  Velocities are written with the legacy ``1/a`` conversion, i.e. LOW by
``a_v/a`` relative to the peculiar velocity at the leapfrog epoch ``a_v``;
:meth:`Catalogue.velocities_at_av` multiplies by ``a/a_v``.  The text reader
casts the printed decimals to float32 so that both forms share one layout; the
cast adds up to one float32 ulp (6e-5 Mpc/h at x > 256) to positions the text
pinned to 1e-4.
"""

from __future__ import annotations

import os
import re
import warnings

import numpy as np

import math

from .meta import meta_float, read_meta
from .product import MetaProduct

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
if CATALOGUE_DTYPE.itemsize != RECORD_BYTES:  # not an assert: must survive python -O
    raise RuntimeError(f"glamio.CATALOGUE_DTYPE itemsize {CATALOGUE_DTYPE.itemsize} != {RECORD_BYTES}")

#: Physical unit (or meaning) of every column.
COLUMN_UNITS: dict[str, str] = {
    "x": "Mpc/h (comoving, in [0, Box))",
    "y": "Mpc/h (comoving, in [0, Box))",
    "z": "Mpc/h (comoving, in [0, Box))",
    "vx": "km/s (peculiar, as written: converted with 1/a)",
    "vy": "km/s (peculiar, as written: converted with 1/a)",
    "vz": "km/s (peculiar, as written: converted with 1/a)",
    "Mbound": "Msun/h (mass of the particles inside the spherical-overdensity radius R_SO that remain bound after unbinding)",
    "Mtot": "Msun/h (all particles inside Rvir)",
    "Rvir": "kpc/h (R_SO plus the finder's resolution correction of up to 0.75 mesh cells; see the sidecar 'definitions')",
    "Vrms": "km/s (3-D rms velocity of the bound particles; same 1/a scale as vx, vy, vz)",
    "Vcirc": "km/s (maximum circular velocity of the bound profile; 0 when unresolved)",
    "row": "1-based row index in this file",
    "Cvir": "concentration (0 when Vcirc is unresolved)",
    "Npart": "number of particles",
    "distinct": "0 = distinct halo",
    "Xoff": "centre-of-mass offset of the bound particles, in units of Rvir",
    "virial_ratio": "2K/|Ep| - 1 (K from the 1/a-scaled velocities)",
    "lambda": "spin parameter (built from the 1/a-scaled velocities)",
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
_PRODUCT_RE = re.compile(r"^(Catshort[A-Za-z])\.(\d+)\.(\d+)\.")


class Catalogue(MetaProduct):
    """A halo catalogue: structured array ``data``, sidecar/header ``meta``, source ``path``.

    ``cat['Mbound']`` returns the column view, ``len(cat)`` the number of
    haloes.  The scalar properties (``a``, ``z``, ``a_v``, ``box``, ``step``,
    ...) come from :class:`glamio.product.MetaProduct`: parsed from ``meta``
    on access, ``None`` when the corresponding key is missing.
    """

    def __init__(self, data: np.ndarray, meta: dict[str, str], path: str):
        super().__init__(meta, path)
        self.data = data

    def __len__(self) -> int:
        return len(self.data)

    def __getitem__(self, name: str) -> np.ndarray:
        return self.data[name]

    def __repr__(self) -> str:
        return f"Catalogue({len(self)} haloes, a={self.meta.get('a', '?').split()[0] if self.meta.get('a') else None}, path={self.path!r})"

    @property
    def columns(self) -> list[str]:
        """Column names in file order."""
        return list(self.data.dtype.names)

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
        if not (math.isfinite(a) and math.isfinite(a_v)) or a_v == 0.0:
            raise ValueError(f"velocities_at_av(): a = {a}, a_v = {a_v} in the meta of {self.path} are not usable")
        factor = a / a_v
        vel = np.stack([self.data[c].astype(np.float64) for c in ("vx", "vy", "vz")], axis=1)
        vel *= factor
        return vel


def read_catalogue(path: str | os.PathLike, mmap: bool = False) -> Catalogue:
    """Read a BDM catalogue from its ``.bin``/``.meta`` pair or its legacy ``.DAT``.

    ``path`` may be the ``.bin``, the ``.meta``, the ``.DAT`` or the stem without
    extension.  A ``.DAT`` path (any letter case) goes to :func:`read_catalogue_text`.
    A bare stem is read as binary when ``<stem>.bin`` exists, as text when only a
    ``.DAT`` exists, and as binary otherwise (which raises ``FileNotFoundError``
    naming the missing file).  ``mmap=True`` maps the ``.bin`` read-only instead of
    loading it (useful for the ~1 GB catalogues of 3200^3 runs); the text reader
    ignores it.
    """
    p = os.fspath(path)
    if p.lower().endswith(".dat"):
        return read_catalogue_text(p)
    if p.endswith(".bin"):
        stem = p[: -len(".bin")]
    elif p.endswith(".meta"):
        stem = p[: -len(".meta")]
    else:
        stem = p
        if not os.path.isfile(stem + ".bin"):
            for ext in (".DAT", ".dat"):
                if os.path.isfile(stem + ext):
                    return read_catalogue_text(stem + ext)
    return _read_catalogue_binary(stem, mmap=mmap)


def _check_layout(meta: dict[str, str], meta_path: str) -> None:
    """Raise ``ValueError`` if the sidecar's ``columns``/``dtypes``/``record_bytes`` differ from CATALOGUE_DTYPE.

    ``columns`` and ``dtypes`` are required: a sidecar without them is not a
    GLAM binary-catalogue sidecar (the Fortran always writes both).
    """
    for key in ("columns", "dtypes"):
        if key not in meta:
            raise ValueError(f"{meta_path}: no '{key}' line; not a binary-catalogue sidecar")
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


def _read_catalogue_binary(stem: str, mmap: bool = False) -> Catalogue:
    """Read ``<stem>.bin`` with the layout checked against ``<stem>.meta``.

    After reading, the ``row`` column must be ``1..N``: a byte-swapped,
    misaligned or truncated-in-the-middle file fails this (a flipped byte in a
    float column cannot be detected here).
    """
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
    if mmap:
        data = np.memmap(bin_path, dtype=CATALOGUE_DTYPE, mode="r")
    else:
        data = np.fromfile(bin_path, dtype=CATALOGUE_DTYPE)
    n = len(data)
    if n and not (data["row"][0] == 1 and data["row"][-1] == n and bool(np.all(np.diff(data["row"]) == 1))):
        raise ValueError(
            f"{bin_path}: the 'row' column is not 1..{n}; wrong byte order, misaligned or corrupted file"
        )
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

    Keys (the same names the binary sidecar uses where both exist): ``header``,
    ``a``, ``da`` (the ``Step`` value), ``a_v = a - da/2``, ``z = 1/a - 1``,
    ``step`` (the ``I`` value), ``realization`` (from the file name), ``N_row``,
    ``N_grid``, ``Omega_m``, ``Omega_L``, ``hubble``, ``buffer_width``,
    ``particle_mass``, ``overdensity``, ``product`` (the ``Catshort?`` label from
    the file name) and ``format = text``; ``rows`` is added by the reader.
    """
    num = r"([-+]?[0-9.]+(?:[eEdD][-+]?[0-9]+)?)"
    meta: dict[str, str] = {"header": lines[0].strip()}
    _grab(rf"\bA\s*=\s*{num}\s+Step\s*=\s*{num}", lines[1], meta, "a", "da")
    _grab(rf"\bI\s*=\s*{num}\s+Nrow\s*=\s*{num}\s+Ngrid\s*=\s*{num}", lines[2], meta, "step", "N_row", "N_grid")
    _grab(rf"Omega_0\s*=\s*{num}", lines[3], meta, "Omega_m")
    _grab(rf"Omega_L\s*=\s*{num}", lines[3], meta, "Omega_L")
    _grab(rf"hubble\s*=\s*{num}", lines[3], meta, "hubble")
    _grab(rf"buffer width\s*\(Mpch\)\s*=\s*{num}", lines[3], meta, "buffer_width")
    _grab(rf"Mass of smallest particle\s*\(Msunh\)\s*=\s*{num}", lines[5], meta, "particle_mass")
    _grab(rf"Overdensity limit\s*=\s*{num}", lines[6], meta, "overdensity")

    a, da = meta_float(meta, "a"), meta_float(meta, "da")
    if a is not None:
        meta["z"] = f"{1.0 / a - 1.0:.6f}"
        if da is not None:
            meta["a_v"] = f"{a - 0.5 * da:.6f}"
    m = _PRODUCT_RE.match(os.path.basename(path))
    if m is not None:
        meta["product"] = m.group(1)
        meta["realization"] = m.group(3).lstrip("0") or "0"
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
        meta["rows"] = "0"
        return Catalogue(np.empty(0, dtype=CATALOGUE_DTYPE), meta, p)
    if raw.shape[1] != _NCOL:          # ndmin=2 guarantees a 2-D array
        raise ValueError(f"{p}: expected {_NCOL} columns per row, found {raw.shape[1]}")
    meta["rows"] = str(raw.shape[0])

    data = np.empty(raw.shape[0], dtype=CATALOGUE_DTYPE)
    for i, name in enumerate(CATALOGUE_DTYPE.names):
        if name in _INT_COLUMNS:
            data[name] = np.rint(raw[:, i]).astype("<i4")
        else:
            data[name] = raw[:, i]
    return Catalogue(data, meta, p)
