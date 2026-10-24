"""Reading packed runs: one HDF5 file per output step instead of ~19 native files.

:mod:`glamio.pack` turns a finished run directory (``Run<box>/``) into
``Output.<step>.<box>.h5`` (one per output step), ``RunInfo.<box>.h5`` (the
schedule, ``BDM.config`` and the logs) and ``packed.md5``.  Layout
(``format = glamio-pack``, ``format_version = 1``):

* every ``.bin`` becomes the dataset ``<group>/data``: the file's bytes
  unchanged (same little-endian type, one contiguous block, no chunking,
  compression or filters), shaped as the native reader shapes the file, so
  it can be memory-mapped at its offset in the HDF5 file;
* every ``.meta`` becomes attributes of ``<group>``: the verbatim sidecar in
  ``meta_text`` (the authoritative copy) plus, for convenience, one string per
  ``key = value`` line split as :func:`glamio.read_meta` splits them (a line
  holding two keys stays one attribute, as in ``read_meta``);
* text files are kept verbatim as ``uint8`` datasets (``PowerDM/R/text``,
  ``CatshortV/outputB``, ``logs/main.log``, ...); ``PowerDM/<kind>/table`` and
  ``schedule/table`` repeat the parsed columns for plain-h5py users;
* ``manifest`` lists every original file: ``path`` (relative to the run
  directory), ``bytes``, ``md5`` and ``object``, the HDF5 dataset holding the
  bytes or, when no dataset of that name exists, ``<group>@<attribute>``.
  :mod:`glamio.pack` rebuilds the files byte for byte from it.

Group names are the native product names (``LensPlanes``, ``deltam``,
``mom_x``, ``Subsample``, ``fR_halo``, ``CatshortV``) and ``PowerDM/R``,
``PowerDM/log``.  A product inside a packed file is named by the path
``<file>.h5/<group>``, e.g. ``Run0/Output.0184.0000.h5/deltam``; the
sidecar-product, binary-catalogue, P(k) and schedule readers accept such a
path where they accept a native one, and :func:`glamio.product_path` returns
one when the native files are absent.  h5py is imported on first use only,
so reading native files still needs numpy alone.

Versioning: ``format_version`` is raised when a reader of the previous
version would misread a new file (an object renamed, retyped or removed, or
its meaning changed); added attributes or objects and the physical placement
(alignment, HDF5 object versions) do not raise it.  Readers refuse files of a
newer version than :data:`FORMAT_VERSION` and files without ``format =
glamio-pack``.
"""

from __future__ import annotations

import functools
import io
import os
import warnings

import numpy as np

#: Suffix of packed files.
PACKED_SUFFIX = ".h5"
#: Value of the root attribute ``format`` of every packed file.
FORMAT = "glamio-pack"
#: Layout version written by :mod:`glamio.pack` (and the newest one the readers accept).
FORMAT_VERSION = 1
#: Group of the halo catalogue inside an output file.
CATALOGUE_GROUP = "CatshortV"


def output_file_name(step: int, box: int) -> str:
    """``Output.<step:04d>.<box:04d>.h5``: the packed file of one output step."""
    return f"Output.{step:04d}.{box:04d}{PACKED_SUFFIX}"


def runinfo_file_name(box: int) -> str:
    """``RunInfo.<box:04d>.h5``: schedule, ``BDM.config`` and logs of one box."""
    return f"RunInfo.{box:04d}{PACKED_SUFFIX}"


def h5py_module():
    """Import h5py, with an installation hint when it is missing."""
    try:
        import h5py
    except ImportError as exc:                                  # pragma: no cover - depends on the environment
        raise ImportError("reading a packed (HDF5) run needs h5py: pip install h5py, "
                          "or install glamio with the 'hdf5' extra") from exc
    return h5py


def split_packed(path: str | os.PathLike) -> tuple[str, str] | None:
    """``(file, object)`` when ``path`` names an object inside an existing ``.h5`` file, else ``None``.

    ``path`` is ``<file>.h5/<object>``; the file part must be a regular file
    (a directory whose name ends in ``.h5`` does not count) and ``object``
    comes back normalised (no empty components, no leading or trailing ``/``).
    """
    p = os.fspath(path)
    marker = PACKED_SUFFIX + "/"
    i = p.find(marker)
    while i >= 0:
        f = p[: i + len(PACKED_SUFFIX)]
        if os.path.isfile(f):
            obj = "/".join(c for c in p[i + len(marker):].split("/") if c)
            return (f, obj) if obj else None
        i = p.find(marker, i + 1)
    return None


def missing_packed_file(path: str | os.PathLike) -> str | None:
    """The ``<file>.h5`` part of ``path`` when ``path`` points inside a packed file that does not exist, else ``None``."""
    p = os.fspath(path)
    marker = PACKED_SUFFIX + "/"
    i = p.find(marker)
    while i >= 0:
        f = p[: i + len(PACKED_SUFFIX)]
        if os.path.isfile(f):
            return None
        if not os.path.exists(f):
            return f
        i = p.find(marker, i + 1)
    return None


def _require(path: str | os.PathLike) -> tuple[str, str]:
    packed = split_packed(path)
    if packed is None:
        missing = missing_packed_file(path)
        raise FileNotFoundError(f"packed file not found: {missing}" if missing else
                                f"{os.fspath(path)} does not name an object inside an existing {PACKED_SUFFIX} file")
    return packed


def check_format(h, name: str) -> None:
    """Raise ``ValueError`` unless the open file ``h`` is a glamio-pack file of a version this glamio reads."""
    fmt = h.attrs.get("format")
    if isinstance(fmt, bytes):
        fmt = fmt.decode()
    if fmt != FORMAT:
        raise ValueError(f"{name}: not a glamio-pack file (format = {fmt!r})")
    version = int(h.attrs.get("format_version", 0))
    if not 1 <= version <= FORMAT_VERSION:
        raise ValueError(f"{name}: format_version {version}; this glamio reads 1..{FORMAT_VERSION} (upgrade glamio)")


def _open(file: str):
    h = h5py_module().File(file, "r")
    try:
        check_format(h, file)
    except Exception:
        h.close()
        raise
    return h


@functools.lru_cache(maxsize=256)
def _index_cached(file: str, mtime_ns: int, size: int, ino: int) -> tuple[frozenset, frozenset]:
    h5py = h5py_module()
    names: set[str] = set()
    products: set[str] = set()
    with _open(file) as h:
        def visit(name, obj):
            names.add(name)
            if isinstance(obj, h5py.Dataset) and "/" in name and name.rsplit("/", 1)[-1] in ("data", "text"):
                products.add(name.rsplit("/", 1)[0])
        h.visititems(visit)
    return frozenset(names), frozenset(products)


def _index(file: str) -> tuple[frozenset, frozenset]:
    """``(all object names, groups holding a product)`` of a packed file, cached while the file is unchanged."""
    st = os.stat(file)
    return _index_cached(os.path.abspath(file), st.st_mtime_ns, st.st_size, st.st_ino)


def as_bytes(value) -> bytes:
    """An HDF5 string attribute (fixed-length bytes or variable-length str) as bytes."""
    if isinstance(value, bytes):
        return bytes(value)
    if isinstance(value, str):
        return value.encode("utf-8")
    if isinstance(value, np.ndarray) and value.dtype == np.uint8:
        return value.tobytes()
    raise TypeError(f"cannot read {type(value).__name__} as text")


def _text(raw: bytes) -> str:
    """Decode as ``open(..., encoding='utf-8', errors='replace')`` would, newlines included."""
    return io.StringIO(raw.decode("utf-8", errors="replace"), newline=None).read()


def exists(path: str | os.PathLike) -> bool:
    """True when ``path`` names an object inside an existing packed file."""
    packed = split_packed(path)
    if packed is None:
        return False
    f, obj = packed
    return obj in _index(f)[0]


def read_meta_text(path: str | os.PathLike) -> str:
    """The verbatim sidecar of the packed product ``<file>.h5/<group>`` (attribute ``meta_text``)."""
    f, obj = _require(path)
    if obj.endswith("/data") or obj == "data":
        obj = obj[: -len("data")].rstrip("/")
    with _open(f) as h:
        g = h.get(obj)
        if g is None:
            raise FileNotFoundError(f"{os.fspath(path)}: no object {obj!r} in {f}")
        if "meta_text" not in g.attrs:
            raise FileNotFoundError(f"{os.fspath(path)}: {obj!r} carries no sidecar (attribute meta_text)")
        return _text(as_bytes(g.attrs["meta_text"]))


def read_text(path: str | os.PathLike) -> str:
    """A verbatim text file stored in a packed file: the ``uint8`` dataset ``path``, or ``path/text`` for a group."""
    f, obj = _require(path)
    h5py = h5py_module()
    with _open(f) as h:
        o = h.get(obj)
        if o is None:
            raise FileNotFoundError(f"{os.fspath(path)}: no object {obj!r} in {f}")
        if isinstance(o, h5py.Group):
            if "text" not in o:
                raise FileNotFoundError(f"{os.fspath(path)}: group {obj!r} has no 'text' dataset")
            o = o["text"]
        if o.dtype != np.uint8 or o.ndim != 1:
            raise ValueError(f"{os.fspath(path)}: {o.name} is {o.dtype} {o.shape}, not a verbatim text (uint8) dataset")
        return _text(o[()].tobytes())


def dataset_nbytes(path: str | os.PathLike) -> int:
    """Byte count of the dataset ``path`` (``<file>.h5/<group>/data``)."""
    f, obj = _require(path)
    h5py = h5py_module()
    with _open(f) as h:
        d = h.get(obj)
        if not isinstance(d, h5py.Dataset):
            raise FileNotFoundError(f"{os.fspath(path)}: no dataset {obj!r} in {f}")
        return int(d.size) * d.dtype.itemsize


def load_dataset(path: str | os.PathLike, dtype: np.dtype | str, count: int, mmap: bool = False) -> np.ndarray:
    """``count`` items of ``dtype`` from the dataset ``path``, flattened as the native readers' ``np.fromfile``.

    The dataset's type must equal ``dtype`` (byte order included) and its
    size ``count``.  ``mmap=True`` maps the dataset read-only at its offset
    in the HDF5 file, which needs the contiguous layout :mod:`glamio.pack`
    writes; any other layout is read into memory instead, with a warning.
    """
    f, obj = _require(path)
    h5py = h5py_module()
    dt = np.dtype(dtype)
    offset = None
    with _open(f) as h:
        d = h.get(obj)
        if not isinstance(d, h5py.Dataset):
            raise FileNotFoundError(f"{os.fspath(path)}: no dataset {obj!r} in {f}")
        if d.dtype != dt:
            raise ValueError(f"{os.fspath(path)}: dataset type {d.dtype}, expected {dt}")
        if int(d.size) != count:
            raise ValueError(f"{os.fspath(path)}: {d.size} items, expected {count}")
        if count == 0:
            return np.zeros(0, dtype=dt)
        if mmap and d.chunks is None:
            offset = d.id.get_offset()                          # None when the storage is not one block
        if offset is None:
            if mmap:
                warnings.warn(f"{os.fspath(path)}: dataset is not one contiguous block; read into memory instead of mapped",
                              stacklevel=2)
            return np.asarray(d[()]).reshape(-1)
    return np.memmap(f, dtype=dt, mode="r", offset=offset, shape=(count,))


def objects_with_data(file: str) -> set[str]:
    """Groups of a packed output file that hold a product (a ``data`` or ``text`` dataset), e.g. ``{'deltam', 'PowerDM/R'}``."""
    return set(_index(file)[1])
