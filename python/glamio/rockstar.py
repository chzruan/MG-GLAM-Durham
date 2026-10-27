"""Reader for the in-line Rockstar halo catalogue (``CATALOGS/Rockstar.<step>.<box>.*``).

With ``Find Rockstar halos = 1`` in ``Init.dat``, GLAM runs the in-line
Rockstar adapter (``rockstar_inline/``) at every class-1 (full) output and
publishes ``CATALOGS/Rockstar.<step>.<box>.bin`` + ``.meta``:

* ``.bin``: headerless, little-endian, one ``record_bytes``-byte record per
  distinct halo (352 bytes today: Rockstar's ``struct halo`` followed by the
  adapter's FoF-group and strict-SO columns), ``rows`` records;
* ``.meta``: a ``key = value`` sidecar (:mod:`glamio.meta`) that describes the
  layout itself: ``columns``, ``dtypes`` (``f4``, ``i8``, ``i4``, ...,
  little-endian), ``offsets`` (bytes into the record) and ``units`` are
  space-separated lists with one entry per column, and ``column_NN = name:
  description`` lines explain each column.  It also carries ``status``
  (``ok``, ``skipped`` or ``failed``), ``reason``, ``rows``, the epochs
  (``a``, ``a_x``, ``z``, ``a_v``, ``da``), the run (``Box``, ``N_grid``,
  ``N_row``, ``N_p``, cosmology, ``MG_flag``, ``MG_model``, ``fR0`` for f(R)),
  ``particle_mass``, ``finder``, the effective configuration (``config.*``),
  FoF statistics, timings and memory.

A skipped output (a memory pre-check) or a failed one (an error inside the
adapter) still publishes a pair: an empty ``.bin`` and a ``.meta`` with
``status = skipped|failed``, the ``reason`` and ``rows = 0``.  These are
legitimate products; :func:`read_rockstar` refuses them unless
``allow_status=True``, so that an analysis does not take an empty catalogue
for an empty volume.  An interrupted write can leave the hidden temporaries
``CATALOGS/.Rockstar.<step>.<box>.{bin,meta}.tmp``; they are not products.

Units and epochs (from the column descriptions): positions comoving Mpc/h in
``[0, Box)``; velocities physical peculiar km/s at the velocity epoch
``a_v = a - da/2`` (no ``1/a`` correction is needed, unlike the BDM
catalogue); radii comoving kpc/h; masses Msun/h.  ``id`` is the 0-based row
number.
"""

from __future__ import annotations

import os
import re

import numpy as np

from .meta import meta_float
from .product import MetaProduct, load_array, load_meta, nbytes, product_exists, resolve_pair

#: Value of the sidecar key ``product``.
PRODUCT = "Rockstar"
#: The values of ``status`` the writer uses; ``ok`` is the only one with rows.
STATUSES = ("ok", "skipped", "failed")
#: Sidecar keys the reader needs.
_REQUIRED = ("columns", "dtypes", "offsets", "record_bytes", "rows", "status")
_TYPE_RE = re.compile(r"^[<|]?([iuf])([1248])$")
_COLUMN_LINE_RE = re.compile(r"^column_\d+$")


class RockstarCatalogue(MetaProduct):
    """An in-line Rockstar catalogue: structured array ``data``, sidecar ``meta``, source ``path``.

    ``cat['m']`` returns a column, ``len(cat)`` the number of distinct haloes
    (0 for a skipped or failed output).  ``status`` and ``reason`` say why a
    catalogue is empty; the scalar properties (``a``, ``z``, ``a_v``, ``box``,
    ``step``, ...) come from :class:`glamio.product.MetaProduct`.
    """

    def __init__(self, data: np.ndarray, meta: dict[str, str], path: str):
        super().__init__(meta, path)
        self.data = data

    def __len__(self) -> int:
        return len(self.data)

    def __getitem__(self, name: str) -> np.ndarray:
        return self.data[name]

    def __repr__(self) -> str:
        return f"RockstarCatalogue({len(self)} haloes, status={self.status}, a={self.a}, path={self.path!r})"

    @property
    def columns(self) -> list[str]:
        """Column names in file order."""
        return list(self.data.dtype.names)

    @property
    def status(self) -> str:
        """``ok``, ``skipped`` or ``failed`` (first token of the meta key ``status``)."""
        return self.meta["status"].split()[0]

    @property
    def ok(self) -> bool:
        return self.status == "ok"

    @property
    def reason(self) -> str | None:
        """Why the output was skipped or failed (meta key ``reason``); ``None`` for the writer's ``-``."""
        value = self.meta.get("reason")
        return None if value is None or value.strip() in ("", "-") else value

    @property
    def record_bytes(self) -> int:
        return self.data.dtype.itemsize

    @property
    def particle_mass(self) -> float | None:
        """Particle mass in Msun/h (meta key ``particle_mass``, Rockstar's value)."""
        return self._float("particle_mass")

    @property
    def units(self) -> dict[str, str]:
        """``{column: unit token}`` from the meta key ``units`` (``-`` = dimensionless or internal); empty when absent."""
        tokens = self.meta.get("units", "").split()
        return dict(zip(self.columns, tokens)) if len(tokens) == len(self.columns) else {}

    @property
    def descriptions(self) -> dict[str, str]:
        """``{column: description}`` from the ``column_NN = name: description`` lines of the meta."""
        out: dict[str, str] = {}
        for key, value in self.meta.items():
            if _COLUMN_LINE_RE.match(key) and ":" in value:
                name, desc = value.split(":", 1)
                out[name.strip()] = desc.strip()
        return out

    def positions(self) -> np.ndarray:
        """Return ``(N, 3)`` float64 comoving positions ``[x, y, z]`` in Mpc/h."""
        return np.stack([self.data[c].astype(np.float64) for c in ("x", "y", "z")], axis=1)

    def velocities(self) -> np.ndarray:
        """Return ``(N, 3)`` float64 halo velocities ``[vx, vy, vz]``: physical peculiar km/s at ``a_v``, as written."""
        return np.stack([self.data[c].astype(np.float64) for c in ("vx", "vy", "vz")], axis=1)


def rockstar_dtype(meta: dict[str, str], meta_path: str = "<meta>") -> np.dtype:
    """The record type described by a Rockstar sidecar's ``columns``, ``dtypes``, ``offsets`` and ``record_bytes``.

    Every type must be a little-endian integer or float (``i4``, ``i8``,
    ``f4``, ...; a leading ``<`` is accepted), the lists must have one entry
    per column (``units`` too, when present), the names must be distinct and
    the fields must lie inside the record without overlapping.  Raises
    ``ValueError`` otherwise.
    """
    names = meta.get("columns", "").split()
    codes = meta.get("dtypes", "").split()
    offs = meta.get("offsets", "").split()
    if not names:
        raise ValueError(f"{meta_path}: empty 'columns'")
    if len(codes) != len(names) or len(offs) != len(names):
        raise ValueError(f"{meta_path}: {len(names)} columns but {len(codes)} dtypes and {len(offs)} offsets")
    units = meta.get("units")
    if units is not None and len(units.split()) != len(names):
        raise ValueError(f"{meta_path}: {len(names)} columns but {len(units.split())} units")
    if len(set(names)) != len(names):
        raise ValueError(f"{meta_path}: repeated column names in 'columns'")
    formats = []
    for name, code in zip(names, codes):
        m = _TYPE_RE.match(code)
        if m is None:
            raise ValueError(f"{meta_path}: column {name!r} has dtype {code!r}; expected a little-endian "
                             "integer or float such as f4, i4, i8")
        formats.append(np.dtype("<" + m.group(1) + m.group(2)))
    try:
        offsets = [int(o) for o in offs]
    except ValueError as exc:
        raise ValueError(f"{meta_path}: unparsable 'offsets': {exc}") from exc
    rb = meta_float(meta, "record_bytes")
    if rb is None or rb != int(rb) or rb <= 0:
        raise ValueError(f"{meta_path}: record_bytes = {meta.get('record_bytes')!r}, expected a positive integer")
    record_bytes = int(rb)
    end = 0
    for off, fmt, name in sorted(zip(offsets, formats, names), key=lambda t: t[0]):
        if off < 0 or off + fmt.itemsize > record_bytes:
            raise ValueError(f"{meta_path}: column {name!r} at offset {off} ({fmt.itemsize} bytes) lies outside "
                             f"the {record_bytes}-byte record")
        if off < end:
            raise ValueError(f"{meta_path}: column {name!r} at offset {off} overlaps the previous column (ends at {end})")
        end = off + fmt.itemsize
    return np.dtype({"names": names, "formats": formats, "offsets": offsets, "itemsize": record_bytes})


def read_rockstar(path: str | os.PathLike, allow_status: bool = False, mmap: bool = False) -> RockstarCatalogue:
    """Read an in-line Rockstar catalogue from its ``.bin``/``.meta`` pair or a packed file.

    ``path`` may be the ``.bin``, the ``.meta``, the stem
    ``CATALOGS/Rockstar.<step>.<box>`` or the packed product
    ``Output.<step>.<box>.h5/Rockstar`` (see :mod:`glamio.packed`).  The
    record type comes from the sidecar (:func:`rockstar_dtype`); the ``.bin``
    must hold exactly ``rows * record_bytes`` bytes and, when an integer
    ``id`` column exists, it must be the row number ``0..N-1`` (a byte-swapped,
    misaligned or corrupted file fails this).  ``mmap=True`` maps the ``.bin``
    read-only instead of loading it.

    A catalogue whose ``status`` is not ``ok`` (``skipped`` or ``failed``,
    always with zero rows) raises ``ValueError`` naming the status and the
    reason, unless ``allow_status=True``: then it is returned with zero rows
    and ``.status`` / ``.reason`` say why.  An unknown status value, or a
    skipped or failed catalogue with rows, always raises.
    """
    bin_path, meta_path, _ = resolve_pair(path)
    for p in (bin_path, meta_path):
        if not product_exists(p):
            raise FileNotFoundError(f"Rockstar catalogue needs both .bin and .meta; missing {p}")
    meta = load_meta(meta_path, expect_product=PRODUCT, required=_REQUIRED)
    tokens = meta["status"].split()
    status = tokens[0] if tokens else ""
    if status not in STATUSES:
        raise ValueError(f"{meta_path}: status = {meta['status']!r}, expected one of {STATUSES}")
    reason = meta.get("reason", "-")
    dt = rockstar_dtype(meta, meta_path)
    rows = meta_float(meta, "rows")
    if rows is None or rows != int(rows) or rows < 0:
        raise ValueError(f"{meta_path}: rows = {meta['rows']!r}, expected a non-negative integer")
    rows = int(rows)
    if status != "ok":
        if rows != 0:
            raise ValueError(f"{meta_path}: status = {status} with rows = {rows}; a {status} catalogue has no rows")
        if not allow_status:
            raise ValueError(f"{meta_path}: Rockstar catalogue status = {status} ({reason}); it has no haloes. "
                             "Pass allow_status=True to read it anyway")
    size = nbytes(bin_path)
    if size != rows * dt.itemsize:
        raise ValueError(f"{bin_path}: size {size} bytes, but {meta_path} says rows = {rows} "
                         f"x record_bytes = {dt.itemsize} ({rows * dt.itemsize} bytes)")
    data = load_array(bin_path, dt, rows, mmap=mmap)
    if rows and "id" in dt.names and dt["id"].kind in "iu":
        ids = data["id"]
        if not (ids[0] == 0 and ids[-1] == rows - 1 and bool(np.all(np.diff(ids) == 1))):
            raise ValueError(f"{bin_path}: the 'id' column is not the row number 0..{rows - 1}; "
                             "wrong byte order, misaligned or corrupted file")
    return RockstarCatalogue(data, meta, bin_path)
