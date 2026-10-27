"""Pack a finished run directory into one HDF5 file per output step; verify, unpack, prune.

usage: python -m glamio.pack pack   RUNDIR [--out DIR] [--box N] [--steps S [S ...]]
       python -m glamio.pack verify FILE_OR_DIR [...] [--against RUNDIR]
       python -m glamio.pack unpack FILE_OR_DIR [...] --out DIR
       python -m glamio.pack prune  RUNDIR [--yes]
       python -m glamio.pack status RUNDIR

``pack`` reads ``RUNDIR`` (the ``Run<box>/`` directory of a finished GLAM
run), checks every binary product with its glamio reader and writes, into
``--out`` (default ``RUNDIR``), ``Output.<step>.<box>.h5`` for each output
step, ``RunInfo.<box>.h5`` (schedule, ``BDM.config``, ``Rockstar.config``, logs) and
``packed.md5`` (md5 of those files, to check a later copy).  The layout is
described in :mod:`glamio.packed`.  Each file is written as
``<name>.partial``, verified (every dataset's bytes in the file against the
md5 of its original, plus the derived tables and attributes), flushed to
disk and only then renamed; a failure removes the ``.partial``.  On a rerun
an existing packed file is kept only if it still matches the run directory:
its manifest must list exactly the current original files of its step, and
every one of them must have the recorded size and md5.  The original files
are never modified.

``verify`` re-checks packed files against their manifests and, for a
directory, that no output file is missing; with ``--against RUNDIR`` it also
re-reads the original files (use it only while they exist).  ``unpack``
rebuilds the original files byte for byte under ``--out`` and never
overwrites.  ``prune`` deletes the original files of a packed ``RUNDIR``:
it refuses while a ``.partial`` exists or while any original is in no
packed file, verifies every packed file and compares every original still
present with its manifest entry (size and md5), deletes nothing if anything
differs, and without ``--yes`` only reports what it would delete.
``status`` prints the state of a run directory (``native``, ``incomplete``,
``packed``, ``partially-pruned``, ``pruned``, ``busy`` or ``inconsistent``),
which ``DEGRACE-companion/pack_run.sh`` uses to resume.

``pack`` refuses a run directory holding anything it cannot place (another
subdirectory, symbolic links, PM particle or checkpoint files, a step-tagged
file of an unknown product, a text catalogue, a ``.bin`` without its
``.meta``, the hidden temporaries ``CATALOGS/.Rockstar.*.tmp`` of an
interrupted Rockstar catalogue write, a wrong box tag, a step missing from
``OutputSchedule.dat``, a ``packed.md5`` it did not write) rather than leave
it behind silently.  The in-line Rockstar catalogues
(``CATALOGS/Rockstar.<step>.<box>.{bin,meta}``) are packed as the group
``Rockstar``, skipped and failed ones (an empty ``.bin``) included.  One
process, one core; memory stays below about 1 GB (blocks of at most 256 MiB).
File modification times and permissions are not recorded.
"""

from __future__ import annotations

import argparse
import datetime
import hashlib
import os
import re
import signal
import sys
import tempfile
import time
from dataclasses import dataclass, field

import numpy as np

from . import __version__
from .catalogue import read_catalogue
from .frhalo import read_frhalo
from .lensplanes import read_lensplanes
from .meta import parse_meta
from .modes import read_modes
from .packed import (CATALOGUE_GROUP, FORMAT, FORMAT_VERSION, ROCKSTAR_GROUP, _text, as_bytes, check_format,
                     h5py_module, output_file_name, runinfo_file_name)
from .powerdm import _parse_powerdm, read_powerdm
from .rockstar import read_rockstar
from .run import RUN_OBJECTS, SIDECAR_PRODUCTS
from .schedule import _parse_schedule, read_schedule
from .subsample import read_subsample

#: Bytes per block when copying or hashing.
BLOCK_BYTES = 1 << 28
#: Name of the checksum list written beside the packed files.
MD5_FILE = "packed.md5"
#: Root attribute ``layout`` of every packed file: where the format is described.
LAYOUT_DOC = ("glamio-pack format 1: one group per product; <group>/data holds the .bin bytes unchanged (contiguous, "
              "memory-mappable at its file offset); <group> attribute meta_text holds the .meta verbatim (the other "
              "attributes are a convenience); text files are uint8 datasets; manifest lists path, bytes, md5 and "
              "object of every original file. See MG-GLAM python/README.md, 'Packed runs (HDF5)'.")
#: h5py.File arguments for writing: objects readable by HDF5 >= 1.10 (and 1.8 for what is used today); every
#: allocation on a 4 KiB boundary, so the memory-mapped arrays are aligned.
_WRITE_KW = {"libver": ("earliest", "v110"), "alignment_threshold": 1, "alignment_interval": 4096}
#: Largest sidecar stored as one string attribute (HDF5's compact attribute limit is 64 KiB).
_MAX_META_BYTES = 60000

#: Reader and array attribute that validate and shape each sidecar product.
_READERS = {"LensPlanes": (read_lensplanes, "planes"), "deltam": (read_modes, "data"), "mom_x": (read_modes, "data"),
            "mom_y": (read_modes, "data"), "mom_z": (read_modes, "data"), "Subsample": (read_subsample, "q"),
            "fR_halo": (read_frhalo, "data")}
#: Order of the groups in an output file.
_ORDER = ("Catshort", "outputB", "Rockstar", "fR_halo", "Subsample", "deltam", "mom_x", "mom_y", "mom_z", "LensPlanes",
          "PowerDM")
#: Run-wide sidecar keys copied, typed, to the root attributes of each output file.
_COMMON_KEYS = (("Box", float), ("N_grid", int), ("N_row", int), ("N_p", int), ("Omega_m", float), ("Omega_L", float),
                ("hubble", float), ("MG_flag", int), ("MG_model", int), ("fR0", float), ("header", str))

_PACKED_RE = re.compile(r"^(Output\.\d{4}\.\d{4}|RunInfo\.\d{4})\.h5$")
_OWN_RE = re.compile(r"^((Output\.\d{4}\.\d{4}|RunInfo\.\d{4})\.h5|packed\.md5)(\.partial)?$")
_MD5_LINE_RE = re.compile(r"^[0-9a-f]{32}  (Output\.\d{4}\.\d{4}|RunInfo\.\d{4})\.h5$")
_SIDECAR_RE = re.compile(r"^(%s)\.(\d{4})\.(\d{4})\.(bin|meta)$" % "|".join(SIDECAR_PRODUCTS))
_POWER_RE = re.compile(r"^PowerDM\.([A-Za-z]+)\.(\d{4})\.(\d{4})\.dat$")
_CAT_RE = re.compile(r"^(Catshort[A-Za-z])\.(\d{4})\.(\d{4})\.(bin|meta)$")
_ROCKSTAR_RE = re.compile(r"^(%s)\.(\d{4})\.(\d{4})\.(bin|meta)$" % ROCKSTAR_GROUP)
_ROCKSTAR_TMP_RE = re.compile(r"^\.%s\." % ROCKSTAR_GROUP)
_OUTPUTB_RE = re.compile(r"^outputB\.(\d{4})\.dat$")
_TAGGED_RE = re.compile(r"\.\d{4}\.")
_PM_RE = re.compile(r"^PMcr[ds]")
_RUNDIR_RE = re.compile(r"^Run(\d+)$")


class PackError(RuntimeError):
    """A run directory or a packed file that cannot be packed, verified or pruned safely."""


@dataclass
class Item:
    """One product of one output step: its HDF5 group and its original files (relative to the run directory)."""
    kind: str                                   # 'sidecar', 'catalogue', 'power' or 'outputB'
    group: str
    files: dict[str, str] = field(default_factory=dict)     # role ('bin', 'meta', 'text') -> relative path

    def rank(self) -> int:
        key = "outputB" if self.kind == "outputB" else self.group
        return next(i for i, prefix in enumerate(_ORDER) if key.startswith(prefix))


@dataclass
class Inventory:
    box: int
    outputs: dict[int, list[Item]]
    runfiles: list[str]
    problems: list[str]

    def step_files(self, step: int) -> list[str]:
        """Original files (relative paths) of one output step, sorted."""
        return sorted(f for item in self.outputs.get(step, []) for f in item.files.values())

    def all_files(self) -> list[str]:
        return sorted([f for step in self.outputs for f in self.step_files(step)] + self.runfiles)


# ----------------------------------------------------------------- inventory
def _box_of(rundir: str, box: int | None) -> int:
    if box is not None:
        return box
    m = _RUNDIR_RE.match(os.path.basename(os.path.abspath(rundir)))
    if m is None:
        raise PackError(f"{rundir}: cannot infer the box from the directory name (expected Run<box>); pass --box")
    return int(m.group(1))


def _is_own_md5(path: str) -> bool:
    """True when ``path`` reads as a checksum list written by :func:`pack` (or is empty)."""
    try:
        with open(path, encoding="ascii") as fh:
            lines = fh.read().splitlines()
    except (OSError, UnicodeDecodeError):
        return False
    return all(_MD5_LINE_RE.match(line) for line in lines if line)


def inventory(rundir: str, box: int | None = None) -> Inventory:
    """Classify every entry of ``rundir`` (and of ``CATALOGS/``) into output products, run files and problems."""
    box = _box_of(rundir, box)
    items: dict[tuple[int, str], Item] = {}
    runfiles: list[str] = []
    problems: list[str] = []

    def add(step: int, kind: str, group: str, role: str, rel: str, tag_box: int | None):
        if tag_box is not None and tag_box != box:
            problems.append(f"{rel}: box tag {tag_box:04d}, expected {box:04d}")
            return
        items.setdefault((step, group), Item(kind, group)).files[role] = rel

    for entry in sorted(os.scandir(rundir), key=lambda e: e.name):
        name = entry.name
        if entry.is_symlink():
            problems.append(f"{name}: symbolic link")
        elif entry.is_dir(follow_symlinks=False):
            if name != "CATALOGS":
                problems.append(f"{name}/: unexpected subdirectory")
                continue
            for sub in sorted(os.scandir(entry.path), key=lambda e: e.name):
                rel = os.path.join("CATALOGS", sub.name)
                m, b = _CAT_RE.match(sub.name) or _ROCKSTAR_RE.match(sub.name), _OUTPUTB_RE.match(sub.name)
                if sub.is_symlink() or not sub.is_file(follow_symlinks=False):
                    problems.append(f"{rel}: not a regular file")
                elif m:
                    add(int(m.group(2)), "catalogue", m.group(1), m.group(4), rel, int(m.group(3)))
                elif b:
                    add(int(b.group(1)), "outputB", f"{CATALOGUE_GROUP}/outputB", "text", rel, None)
                elif _ROCKSTAR_TMP_RE.match(sub.name):
                    problems.append(f"{rel}: hidden temporary of an interrupted Rockstar catalogue write (GLAM "
                                    "stopped while publishing, or an unpack into this directory was killed; the step's "
                                    "Rockstar catalogue may be missing or older): check the run with check_run.py and "
                                    "delete the temporary once no job writes to this run")
                else:
                    problems.append(f"{rel}: unknown file in CATALOGS/ (only Catshort?.<step>.<box>.bin/.meta, "
                                    "Rockstar.<step>.<box>.bin/.meta and outputB.<step>.dat are packed)")
        elif not entry.is_file(follow_symlinks=False):
            problems.append(f"{name}: not a regular file")
        elif name == MD5_FILE and not _is_own_md5(entry.path):
            problems.append(f"{name}: not a checksum list written by glamio.pack (a user file?); move it away first")
        elif _OWN_RE.match(name):
            continue
        elif (m := _SIDECAR_RE.match(name)):
            add(int(m.group(2)), "sidecar", m.group(1), m.group(4), name, int(m.group(3)))
        elif (m := _POWER_RE.match(name)):
            add(int(m.group(2)), "power", f"PowerDM/{m.group(1)}", "text", name, int(m.group(3)))
        elif _PM_RE.match(name):
            problems.append(f"{name}: PM particle or checkpoint file (delete or move the particle files before packing)")
        elif _TAGGED_RE.search(name):
            problems.append(f"{name}: step-tagged file of an unknown product")
        else:
            runfiles.append(name)

    for (step, group), item in items.items():
        if item.kind in ("sidecar", "catalogue") and set(item.files) != {"bin", "meta"}:
            problems.append(f"step {step} {group}: needs both .bin and .meta, found {sorted(item.files.values())}")
    outputs: dict[int, list[Item]] = {}
    for (step, _), item in sorted(items.items()):
        outputs.setdefault(step, []).append(item)
    for step in outputs:
        outputs[step].sort(key=Item.rank)
    if "OutputSchedule.dat" not in runfiles:
        problems.append("OutputSchedule.dat missing")
    else:
        scheduled = set(int(s) for s in read_schedule(os.path.join(rundir, "OutputSchedule.dat")).steps)
        for step in sorted(set(outputs) - scheduled):
            problems.append(f"step {step}: products present but the step is not in OutputSchedule.dat")
    return Inventory(box, outputs, runfiles, problems)


# ----------------------------------------------------------------- byte helpers
def _md5_bytes(raw: bytes) -> str:
    return hashlib.md5(raw).hexdigest()


def _md5_file(path: str) -> tuple[int, str]:
    h, n = hashlib.md5(), 0
    with open(path, "rb") as fh:
        while chunk := fh.read(1 << 24):
            h.update(chunk)
            n += len(chunk)
    return n, h.hexdigest()


def _fsync(path: str) -> None:
    """Flush a file (or a directory entry) to stable storage."""
    fd = os.open(path, os.O_RDONLY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def _blocks(shape: tuple[int, ...], itemsize: int):
    """``(index, block_shape)`` covering an array of ``shape`` in C order, each block at most ``BLOCK_BYTES``.

    Rows along axis 0 are grouped while they fit; a single row larger than
    ``BLOCK_BYTES`` is split along the next axis (recursively), so the blocks
    stay bounded whatever the shape.  Concatenating the blocks in order gives
    the array's bytes in C order.
    """
    def rec(prefix, dims):
        sub = itemsize * int(np.prod(dims[1:], dtype=np.int64))
        if len(dims) == 1 or sub <= BLOCK_BYTES:
            step = max(1, BLOCK_BYTES // max(sub, 1))
            for i0 in range(0, dims[0], step):
                i1 = min(dims[0], i0 + step)
                yield prefix + (slice(i0, i1),), (i1 - i0,) + tuple(dims[1:])
        else:
            for i in range(dims[0]):
                yield from rec(prefix + (i,), dims[1:])

    if len(shape) == 0 or 0 in shape:
        return
    yield from rec((), tuple(int(n) for n in shape))


def _copy_into(dset, path: str) -> tuple[int, str]:
    """Fill ``dset`` with the bytes of ``path`` block by block; return (bytes, md5 of the file)."""
    shape, dt = dset.shape, dset.dtype
    total = dt.itemsize * int(np.prod(shape, dtype=np.int64))
    size = os.path.getsize(path)
    if size != total:
        raise PackError(f"{path}: {size} bytes, the dataset needs {total}")
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for idx, bshape in _blocks(shape, dt.itemsize):
            n = dt.itemsize * int(np.prod(bshape, dtype=np.int64))
            buf = fh.read(n)
            if len(buf) != n:
                raise PackError(f"{path}: short read at block {idx}")
            h.update(buf)
            dset[idx] = np.frombuffer(buf, dtype=dt).reshape(bshape)
        if fh.read(1):
            raise PackError(f"{path}: grew while it was being packed")
    return total, h.hexdigest()


def _iter_dataset_bytes(dset):
    """The bytes of ``dset`` in C order, block by block, as uint8 arrays read through h5py."""
    for idx, _ in _blocks(dset.shape, dset.dtype.itemsize):
        yield np.ascontiguousarray(dset[idx]).reshape(-1).view(np.uint8)


def _iter_raw_bytes(path: str, offset: int, nbytes: int):
    """The bytes ``[offset, offset + nbytes)`` of ``path``, block by block (what a memory map of the dataset sees)."""
    with open(path, "rb") as fh:
        fh.seek(offset)
        left = nbytes
        while left > 0:
            buf = fh.read(min(left, BLOCK_BYTES))
            if not buf:
                raise PackError(f"{path}: ends inside a dataset at byte {offset + nbytes - left}")
            left -= len(buf)
            yield buf


def _base_alignment(dt: np.dtype) -> int:
    """Alignment the native memory reads want: the largest field alignment of a compound type."""
    if dt.fields:
        return max(_base_alignment(f[0]) for f in dt.fields.values())
    return dt.alignment


def _entry_blocks(h, obj: str):
    """Byte blocks of a manifest object: the dataset ``obj`` or, when no dataset has that name, ``<group>@<attr>``."""
    h5py = h5py_module()
    if obj in h and isinstance(h[obj], h5py.Dataset):
        return _iter_dataset_bytes(h[obj])
    if "@" in obj:
        grp, attr = obj.rsplit("@", 1)
        if grp in h and attr in h[grp].attrs:
            return iter([np.frombuffer(as_bytes(h[grp].attrs[attr]), dtype=np.uint8)])
    raise PackError(f"{h.filename}: manifest object {obj} missing")


def _put_file(h, obj: str, path: str) -> tuple[int, str]:
    """Store the file ``path`` verbatim as the uint8 dataset ``obj``."""
    d = h.create_dataset(obj, shape=(os.path.getsize(path),), dtype=np.uint8)
    return _copy_into(d, path)


def _typed(value: str, kind):
    if kind is str:
        return value
    token = value.split()[0] if value.split() else ""
    try:
        x = float(token.replace("D", "E").replace("d", "e"))
        return np.int64(round(x)) if kind is int else np.float64(x)
    except ValueError:
        return token                                # e.g. fR0 = none in a LCDM run


def _manifest_array(entries: list[tuple[str, int, str, str]]) -> np.ndarray:
    lp = max([len(e[0]) for e in entries] + [1])
    lo = max([len(e[3]) for e in entries] + [1])
    dt = np.dtype([("path", f"S{lp}"), ("bytes", "<i8"), ("md5", "S32"), ("object", f"S{lo}")])
    return np.array([(p.encode(), n, m.encode(), o.encode()) for p, n, m, o in entries], dtype=dt)


def _manifest_paths(path: str) -> list[str]:
    with h5py_module().File(path, "r") as h:
        return sorted(e["path"].decode() for e in h["manifest"][()])


def _code_sha256() -> str:
    """sha256 of the glamio sources in use (provenance of the packing code)."""
    pkg = os.path.dirname(os.path.abspath(__file__))
    h = hashlib.sha256()
    for name in sorted(n for n in os.listdir(pkg) if n.endswith(".py")):
        h.update(name.encode() + b"\0")
        with open(os.path.join(pkg, name), "rb") as fh:
            h.update(fh.read())
    return h.hexdigest()


def _run_attrs(rundir: str, box: int) -> dict:
    """Root attributes shared by the output files and RunInfo of one box."""
    absdir = os.path.abspath(rundir)
    attrs = {"format": FORMAT, "format_version": np.int32(FORMAT_VERSION), "layout": LAYOUT_DOC,
             "glamio_version": __version__, "glamio_code_sha256": _code_sha256(),
             "realization": np.int32(box), "run": os.path.basename(os.path.dirname(absdir)),
             "box_dir": os.path.basename(absdir), "source_dir": absdir,
             "packed_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")}
    build = os.path.join(os.path.dirname(absdir), "bin", "git_hash.txt")
    if os.path.isfile(build):
        attrs["glam_build"] = open(build, encoding="utf-8", errors="replace").read().strip()
    return attrs


# ----------------------------------------------------------------- writing
def _write_sidecar(h, rundir: str, item: Item, manifest: list) -> dict[str, str]:
    bin_path = os.path.join(rundir, item.files["bin"])
    meta_path = os.path.join(rundir, item.files["meta"])
    if item.kind == "catalogue" and item.group == ROCKSTAR_GROUP:
        obj = read_rockstar(bin_path, allow_status=True, mmap=True)     # skipped and failed (empty) catalogues pack too
        arr = obj.data
    elif item.kind == "catalogue":
        obj = read_catalogue(bin_path, mmap=True)
        arr = obj.data
    else:
        reader, attr = _READERS[item.group]
        obj = reader(bin_path, mmap=True)
        arr = getattr(obj, attr)
    shape, dtype, meta = arr.shape, arr.dtype, obj.meta
    del obj, arr                                    # release the memory map
    with open(meta_path, "rb") as fh:
        raw = fh.read()
    if not raw or raw.endswith(b"\0") or len(raw) > _MAX_META_BYTES:
        raise PackError(f"{meta_path}: a sidecar must be non-empty, below {_MAX_META_BYTES} bytes and not end in NUL "
                        f"to be stored as a string attribute ({len(raw)} bytes)")
    g = h.create_group(item.group)
    g.attrs["meta_text"] = np.bytes_(raw)
    for key, value in meta.items():
        if key not in ("meta_text",):
            try:
                g.attrs[key] = value
            except (KeyError, TypeError, ValueError):   # a name HDF5 refuses: meta_text still holds the line
                pass
    manifest.append((item.files["meta"], len(raw), _md5_bytes(raw), f"{item.group}@meta_text"))
    d = g.create_dataset("data", shape=shape, dtype=dtype)     # contiguous: no chunks, no filters
    n, md5 = _copy_into(d, bin_path)
    manifest.append((item.files["bin"], n, md5, f"{item.group}/data"))
    return meta


def _power_table(p) -> np.ndarray:
    table = np.zeros(len(p), dtype=[(c, "<i4" if c == "bin" else "<f8") for c in p.columns])
    for i, c in enumerate(p.columns):
        table[c] = np.rint(p.data[:, i]) if c == "bin" else p.data[:, i]
    return table


def _write_power(h, rundir: str, item: Item, manifest: list) -> None:
    path = os.path.join(rundir, item.files["text"])
    p = read_powerdm(path)                          # refuses a malformed table
    g = h.create_group(item.group)
    n, md5 = _put_file(g, "text", path)
    manifest.append((item.files["text"], n, md5, f"{item.group}/text"))
    g.create_dataset("table", data=_power_table(p))
    g.attrs.update({"header": p.header, "a": p.a, "z": p.z, "step": np.int32(p.step), "ngrid": np.int32(p.ngrid),
                    "kNyq_printed": p.k_nyq, "binning": p.binning or "", "columns": " ".join(p.columns)})


def _write_output(path: str, rundir: str, inv: Inventory, step: int, row, run_attrs: dict) -> None:
    h5py = h5py_module()
    manifest: list[tuple[str, int, str, str]] = []
    metas: dict[str, dict[str, str]] = {}
    with h5py.File(path, "w", **_WRITE_KW) as h:
        for item in inv.outputs[step]:
            if item.kind in ("sidecar", "catalogue"):
                metas[item.group] = _write_sidecar(h, rundir, item, manifest)
            elif item.kind == "power":
                _write_power(h, rundir, item, manifest)
            else:                                   # outputB: the finder's log of this step
                g = h.require_group(CATALOGUE_GROUP)
                n, md5 = _put_file(g, "outputB", os.path.join(rundir, item.files["text"]))
                manifest.append((item.files["text"], n, md5, item.group))
        h.attrs.update(run_attrs)
        h.attrs.update({"step": np.int32(step), "a": np.float64(row["a"]), "z": np.float64(row["z"]),
                        "chi": np.float64(row["chi"]), "exact": np.int32(row["exact"]), "class": np.int32(row["class"]),
                        "products": " ".join(item.group for item in inv.outputs[step] if item.kind != "outputB")})
        ref = metas.get("LensPlanes") or next(iter(metas.values()), {})
        for key, kind in _COMMON_KEYS:
            if key in ref:
                h.attrs[key] = _typed(ref[key], kind)
        h.create_dataset("manifest", data=_manifest_array(manifest))


def _write_runinfo(path: str, rundir: str, inv: Inventory, run_attrs: dict) -> None:
    h5py = h5py_module()
    manifest: list[tuple[str, int, str, str]] = []
    with h5py.File(path, "w", **_WRITE_KW) as h:
        for name in inv.runfiles:
            src = os.path.join(rundir, name)
            if name == "OutputSchedule.dat":
                obj = "schedule/text"
                h.create_dataset("schedule/table", data=read_schedule(src).rows)
            else:
                obj = RUN_OBJECTS.get(name, f"logs/{name}")     # BDM.config, Rockstar.config at the root
            n, md5 = _put_file(h, obj, src)
            manifest.append((name, n, md5, obj))
        h.attrs.update(run_attrs)
        h.attrs["steps"] = np.array(sorted(inv.outputs), dtype=np.int32)   # steps with products: one Output file each
        h.create_dataset("manifest", data=_manifest_array(manifest))


# ----------------------------------------------------------------- verification
def _check_derived(h, path: str) -> list[str]:
    """The convenience copies (per-key attributes, parsed tables, typed root attributes) against the verbatim texts."""
    h5py = h5py_module()
    name = os.path.basename(path)
    problems: list[str] = []
    metas: dict[str, dict[str, str]] = {}
    for grp, g in h.items():
        if isinstance(g, h5py.Group) and "meta_text" in g.attrs:
            meta = parse_meta(as_bytes(g.attrs["meta_text"]).decode("utf-8", errors="replace"))
            metas[grp] = meta
            for key, value in meta.items():
                if key in g.attrs and g.attrs[key] != value:
                    problems.append(f"{name}: {grp} attribute {key!r} = {g.attrs[key]!r}, meta_text says {value!r}")
    if "PowerDM" in h:
        for kind in h["PowerDM"]:
            g = h["PowerDM"][kind]
            if "table" in g and "text" in g:
                p = _parse_powerdm(_text(g["text"][()].tobytes()), f"{name}/PowerDM/{kind}", f"PowerDM.{kind}.", None)
                t = g["table"][()]
                if list(t.dtype.names) != p.columns or len(t) != len(p) or not all(
                        np.array_equal(t[c].astype(np.float64), p.data[:, i]) for i, c in enumerate(p.columns)):
                    problems.append(f"{name}: PowerDM/{kind}/table differs from its text")
    if "schedule" in h and "table" in h["schedule"] and "text" in h["schedule"]:
        sched = _parse_schedule(_text(h["schedule/text"][()].tobytes()), f"{name}/schedule")
        if not np.array_equal(h["schedule/table"][()], sched.rows):
            problems.append(f"{name}: schedule/table differs from its text")
    ref = metas.get("LensPlanes") or next(iter(metas.values()), {})
    for key, kind in _COMMON_KEYS:
        if key in h.attrs and key in ref:
            want, got = _typed(ref[key], kind), h.attrs[key]
            if (got != want) if isinstance(want, str) else (isinstance(got, str) or float(got) != float(want)):
                problems.append(f"{name}: root attribute {key} = {got!r}, sidecar says {want!r}")
    return problems


def verify_file(path: str, against: str | None = None, missing_ok: bool = False) -> tuple[list[str], dict]:
    """Check a packed file: every manifest entry has the recorded size and md5, the convenience copies agree.

    Each dataset is hashed over its bytes in the file (the region a memory
    map sees); its first and last blocks are also read through h5py and must
    match.  Binary datasets must be single contiguous blocks.  With
    ``against`` (the run directory) each original file still present is
    re-read and compared too; a missing original is a problem unless
    ``missing_ok``.  Returns ``(problems, stats)``; ``stats['unaligned']``
    lists data datasets whose offset does not suit their type (still
    readable, slower to map).
    """
    h5py = h5py_module()
    problems: list[str] = []
    stats = {"entries": 0, "bytes": 0, "originals_checked": 0, "originals": [], "unaligned": []}
    name = os.path.basename(path)
    try:
        h = h5py.File(path, "r")
    except OSError as exc:
        return [f"{name}: cannot open ({exc})"], stats
    with h:
        try:
            check_format(h, name)
        except ValueError as exc:
            return [str(exc)], stats
        if "manifest" not in h:
            return [f"{name}: no manifest"], stats
        for e in h["manifest"][()]:
            rel, size, md5, obj = e["path"].decode(), int(e["bytes"]), e["md5"].decode(), e["object"].decode()
            hh, got_size = hashlib.md5(), 0
            if obj in h and isinstance(h[obj], h5py.Dataset):
                d = h[obj]
                nbytes = int(d.size) * d.dtype.itemsize
                offset = d.id.get_offset() if d.size and d.chunks is None else None
                if d.size and offset is None:
                    problems.append(f"{name}: {obj} is not one contiguous block")
                    blocks = _iter_dataset_bytes(d)
                else:
                    blocks = _iter_raw_bytes(path, offset, nbytes) if nbytes else iter(())
                    if nbytes:
                        if offset % _base_alignment(d.dtype):
                            stats["unaligned"].append(obj)
                        spots = list(_blocks(d.shape, d.dtype.itemsize))
                        pos, edges = 0, []
                        for idx, bshape in spots:
                            n = d.dtype.itemsize * int(np.prod(bshape, dtype=np.int64))
                            edges.append((idx, pos, n))
                            pos += n
                        with open(path, "rb") as fh:
                            for idx, start, n in {edges[0][1]: edges[0], edges[-1][1]: edges[-1]}.values():
                                fh.seek(offset + start)
                                if fh.read(n) != np.ascontiguousarray(d[idx]).tobytes():
                                    problems.append(f"{name}: {obj} reads differently through h5py and from the file")
            elif "@" in obj:
                grp, attr = obj.rsplit("@", 1)
                if grp not in h or attr not in h[grp].attrs:
                    problems.append(f"{name}: {obj} missing")
                    continue
                blocks = iter([as_bytes(h[grp].attrs[attr])])
            else:
                problems.append(f"{name}: dataset {obj} missing")
                continue
            for block in blocks:
                hh.update(block)
                got_size += len(block)
            got_md5 = hh.hexdigest()
            stats["entries"] += 1
            stats["bytes"] += got_size
            if (got_size, got_md5) != (size, md5):
                problems.append(f"{name}: {obj} ({rel}) has {got_size} bytes md5 {got_md5}, manifest {size} {md5}")
            if against is not None:
                native = os.path.join(against, rel)
                if os.path.isfile(native):
                    if _md5_file(native) != (size, md5):
                        problems.append(f"{name}: original {native} differs from the manifest")
                    stats["originals_checked"] += 1
                    stats["originals"].append(native)
                elif not missing_ok:
                    problems.append(f"{name}: original {native} not found")
        problems += _check_derived(h, path)
    return problems, stats


def _completeness(directory: str) -> list[str]:
    """Problems with the set of packed files in ``directory``: missing or extra outputs, a stale packed.md5, .partial files."""
    names = sorted(os.listdir(directory))
    problems = [f"{n}: left by a pack that is running or crashed" for n in names if n.endswith(".partial")]
    runinfos = [n for n in names if re.match(r"^RunInfo\.\d{4}\.h5$", n)]
    outputs = {n for n in names if re.match(r"^Output\.\d{4}\.\d{4}\.h5$", n)}
    if len(runinfos) != 1:
        problems.append(f"{directory}: {len(runinfos)} RunInfo files, expected 1")
    else:
        try:
            with h5py_module().File(os.path.join(directory, runinfos[0]), "r") as h:
                check_format(h, runinfos[0])
                box = int(h.attrs["realization"])
                steps = [int(s) for s in h.attrs["steps"]]
        except (OSError, KeyError, ValueError) as exc:
            return problems + [f"{runinfos[0]}: unreadable ({exc})"]
        want = {output_file_name(s, box) for s in steps}
        problems += [f"{n}: missing (RunInfo lists step {n[7:11]})" for n in sorted(want - outputs)]
        problems += [f"{n}: not listed in {runinfos[0]}" for n in sorted(outputs - want)]
    if MD5_FILE in names:
        with open(os.path.join(directory, MD5_FILE)) as fh:
            listed = {line.split()[-1] for line in fh if line.strip()}
        present = outputs | set(runinfos)
        if listed != present:
            problems.append(f"{MD5_FILE}: lists {len(listed)} files, the directory holds {len(present)} packed files")
    return problems


def _verified_or_raise(path: str) -> dict:
    problems, stats = verify_file(path)
    if problems:
        raise PackError(f"{path} fails verification:\n  " + "\n  ".join(problems))
    return stats


def _check_existing(final: str, rundir: str, expected: list[str], log) -> None:
    """Keep an existing packed file only if it holds exactly the current originals, byte for byte."""
    problems, _ = verify_file(final, against=rundir)
    try:
        listed = _manifest_paths(final)
    except (OSError, KeyError) as exc:
        listed = []
        problems.append(f"no readable manifest ({exc})")
    if listed != expected:
        problems.append(f"its manifest lists {len(listed)} files, the run directory has {len(expected)} for it: "
                        f"only in the file {sorted(set(listed) - set(expected))[:5]}, "
                        f"only in the directory {sorted(set(expected) - set(listed))[:5]}")
    if problems:
        raise PackError(f"{final} exists but does not match the originals in {rundir}:\n  " + "\n  ".join(problems)
                        + "\n  delete it (or the whole packed set) to repack")
    with h5py_module().File(final, "r") as h:
        src = h.attrs.get("source_dir")
    if src != os.path.abspath(rundir):
        log(f"note: {os.path.basename(final)} was packed from {src}; its content matches {rundir}")


def _packed_files(paths: list[str]) -> list[str]:
    out: list[str] = []
    for p in paths:
        if os.path.isdir(p):
            out += [os.path.join(p, n) for n in sorted(os.listdir(p)) if _PACKED_RE.match(n)]
        else:
            out.append(p)
    return out


def _write_verified(final: str, writer, log) -> None:
    """Write ``final`` through ``writer(tmp)`` as ``<final>.partial``, verify, flush, rename; remove the partial on failure."""
    tmp = final + ".partial"
    try:
        os.close(os.open(tmp, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644))
    except FileExistsError:
        raise PackError(f"{tmp} exists: another pack is running, or one crashed (check, then delete it)") from None
    try:
        writer(tmp)
        _verified_or_raise(tmp)
        _fsync(tmp)
        os.replace(tmp, final)
        _fsync(os.path.dirname(os.path.abspath(final)))
    except BaseException:
        if os.path.exists(tmp):
            os.remove(tmp)
        raise


# ----------------------------------------------------------------- commands
def pack(rundir: str, out: str | None = None, box: int | None = None, steps: list[int] | None = None,
         log=print) -> list[str]:
    """Pack ``rundir`` into ``out`` (default ``rundir``); return the packed files, raise :class:`PackError` on any problem."""
    out = out or rundir
    inv = inventory(rundir, box)
    if inv.problems:
        raise PackError(f"{rundir}: cannot pack:\n  " + "\n  ".join(inv.problems))
    os.makedirs(out, exist_ok=True)
    sched = read_schedule(os.path.join(rundir, "OutputSchedule.dat"))
    todo = sorted(inv.outputs) if steps is None else sorted(set(steps))
    missing = [s for s in todo if s not in inv.outputs]
    if missing:
        raise PackError(f"steps {missing} have no products in {rundir}")
    run_attrs = _run_attrs(rundir, inv.box)
    written: list[str] = []
    t_start = time.time()
    new_in = new_files = kept = 0
    for step in todo:
        final = os.path.join(out, output_file_name(step, inv.box))
        files = inv.step_files(step)
        n_bytes = sum(os.path.getsize(os.path.join(rundir, f)) for f in files)
        if os.path.exists(final):
            _check_existing(final, rundir, files, log)
            log(f"step {step:4d}: {os.path.basename(final)} exists and matches the originals, kept")
            written.append(final)
            kept += 1
            continue
        t0 = time.time()
        _write_verified(final, lambda tmp: _write_output(tmp, rundir, inv, step, sched.row(step), run_attrs), log)
        log(f"step {step:4d} class {int(sched.row(step)['class'])}: {len(files):2d} files {n_bytes / 1e9:7.3f} GB -> "
            f"{os.path.basename(final)} {os.path.getsize(final) / 1e9:7.3f} GB  {time.time() - t0:6.1f} s (write + verify)")
        written.append(final)
        new_in += n_bytes
        new_files += len(files)
    final = os.path.join(out, runinfo_file_name(inv.box))
    if os.path.exists(final):
        _check_existing(final, rundir, sorted(inv.runfiles), log)
        log(f"{os.path.basename(final)} exists and matches the originals, kept")
        kept += 1
    else:
        _write_verified(final, lambda tmp: _write_runinfo(tmp, rundir, inv, run_attrs), log)
        log(f"{os.path.basename(final)}: {len(inv.runfiles)} run files")
        new_files += len(inv.runfiles)
        new_in += sum(os.path.getsize(os.path.join(rundir, f)) for f in inv.runfiles)
    written.append(final)
    t0 = time.time()
    md5_path = os.path.join(out, MD5_FILE)
    if os.path.exists(md5_path) and not _is_own_md5(md5_path):
        raise PackError(f"{md5_path} is not a checksum list written by glamio.pack (a user file?); move it away first")
    present = _packed_files([out])
    try:
        fh = open(md5_path + ".partial", "x")
    except FileExistsError:
        raise PackError(f"{md5_path}.partial exists: another pack is running, or one crashed (check, then delete it)") from None
    try:
        with fh:
            for p in present:
                fh.write(f"{_md5_file(p)[1]}  {os.path.basename(p)}\n")
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(md5_path + ".partial", md5_path)
    except BaseException:
        if os.path.exists(md5_path + ".partial"):
            os.remove(md5_path + ".partial")
        raise
    _fsync(os.path.abspath(out))
    log(f"{MD5_FILE}: {len(present)} files, {time.time() - t0:.1f} s")
    log(f"packed {new_files} files ({new_in / 1e9:.3f} GB) of {rundir} into {len(written) - kept} new files in {out}"
        + (f"; {kept} existing files matched the originals and were kept" if kept else "")
        + f"; {time.time() - t_start:.0f} s")
    return written


def unpack(paths: list[str], out: str, log=print) -> int:
    """Rebuild the original files of the packed files ``paths`` under ``out``; return the number of files written.

    Each file is written to a fresh temporary name, checked against the
    manifest md5 and then linked to its name, which fails if the name exists:
    nothing is ever overwritten.
    """
    h5py = h5py_module()
    count = 0
    root = os.path.abspath(out)
    umask = os.umask(0)
    os.umask(umask)
    files = _packed_files(paths)
    if not files:
        raise PackError(f"no packed files in {paths}")
    for path in files:
        with h5py.File(path, "r") as h:
            check_format(h, path)
            for e in h["manifest"][()]:
                rel, size, md5, obj = e["path"].decode(), int(e["bytes"]), e["md5"].decode(), e["object"].decode()
                target = os.path.abspath(os.path.join(root, rel))
                if os.path.commonpath([root, target]) != root:
                    raise PackError(f"{path}: manifest path {rel!r} leaves the output directory")
                if os.path.lexists(target):
                    raise PackError(f"{target} exists; unpack does not overwrite")
                os.makedirs(os.path.dirname(target), exist_ok=True)
                fd, tmp = tempfile.mkstemp(dir=os.path.dirname(target), prefix=f".{os.path.basename(target)}.", suffix=".unpack")
                try:
                    hh, n = hashlib.md5(), 0
                    with os.fdopen(fd, "wb") as fh:
                        for block in _entry_blocks(h, obj):
                            fh.write(block)
                            hh.update(block)
                            n += len(block)
                    if (n, hh.hexdigest()) != (size, md5):
                        raise PackError(f"{path}: {obj} rebuilds to {n} bytes md5 {hh.hexdigest()}, manifest {size} {md5}")
                    os.chmod(tmp, 0o666 & ~umask)
                    try:
                        os.link(tmp, target)
                    except FileExistsError:
                        raise PackError(f"{target} appeared during unpack; not overwritten") from None
                finally:
                    if os.path.lexists(tmp):
                        os.remove(tmp)
                count += 1
        log(f"{os.path.basename(path)}: unpacked")
    return count


def status(rundir: str) -> tuple[str, list[str]]:
    """State of a run directory and a few example files.

    ``native`` (no packed files), ``incomplete`` (some originals are in no
    packed file and none is missing: pack resumes), ``packed`` (every
    original is in a packed file and present), ``partially-pruned`` (every
    remaining original is in a packed file, some are gone: prune resumes),
    ``pruned`` (no originals left), ``busy`` (a ``.partial`` exists) or
    ``inconsistent`` (originals missing and others unpacked, unreadable
    packed files, ...).  Only names and manifests are read, no data.
    """
    partial = [n for n in sorted(os.listdir(rundir)) if n.endswith(".partial")]
    if partial:
        return "busy", partial
    files = _packed_files([rundir])
    natives = _native_names(rundir)
    if not files:
        return "native", []
    listed: set[str] = set()
    try:
        for f in files:
            listed |= set(_manifest_paths(f))
    except (OSError, KeyError) as exc:
        return "inconsistent", [f"unreadable packed file: {exc}"]
    uncovered, gone = sorted(natives - listed), sorted(listed - natives)
    if uncovered:
        return ("incomplete" if not gone else "inconsistent"), uncovered[:5] + (["... and originals missing: "] + gone[:5] if gone else [])
    problems = _completeness(rundir)                # every original is packed: the packed set must be whole
    if problems:
        return "inconsistent", problems[:5]
    if not natives:
        return "pruned", []
    return ("packed", []) if not gone else ("partially-pruned", gone[:5])


def _native_names(rundir: str) -> set[str]:
    """Every entry of ``rundir`` (and of ``CATALOGS/``) that is not one of pack's own files, relative to ``rundir``."""
    natives: set[str] = set()
    for n in os.listdir(rundir):
        p = os.path.join(rundir, n)
        if n == "CATALOGS" and os.path.isdir(p) and not os.path.islink(p):
            natives |= {os.path.join("CATALOGS", m) for m in os.listdir(p)}
        elif not _OWN_RE.match(n):
            natives.add(n)
    return natives


def prune(rundir: str, yes: bool = False, log=print) -> int:
    """Delete the original files of a packed ``rundir`` after full verification; return the number deleted (0 without ``yes``)."""
    files = _packed_files([rundir])
    if not files:
        raise PackError(f"{rundir}: no packed files")
    partial = [n for n in os.listdir(rundir) if n.endswith(".partial")]
    if partial:
        raise PackError(f"{rundir}: nothing deleted: {partial[:3]} exist (a pack is running or crashed)")
    boxes = set()
    for path in files:
        try:
            with h5py_module().File(path, "r") as h:
                check_format(h, path)
                boxes.add(int(h.attrs["realization"]))
        except (OSError, KeyError, ValueError) as exc:
            raise PackError(f"{rundir}: nothing deleted: {os.path.basename(path)} unreadable ({exc})") from None
    m = _RUNDIR_RE.match(os.path.basename(os.path.abspath(rundir)))
    if len(boxes) != 1 or (m is not None and int(m.group(1)) not in boxes):
        raise PackError(f"{rundir}: nothing deleted: packed files of boxes {sorted(boxes)} in a directory named "
                        f"{os.path.basename(os.path.abspath(rundir))}")
    box = boxes.pop()
    incomplete = _completeness(rundir)
    if incomplete:
        raise PackError(f"{rundir}: nothing deleted: the packed set is not whole:\n  " + "\n  ".join(incomplete))
    inv = inventory(rundir, box)
    originals: list[str] = []
    problems: list[str] = []
    for path in files:
        p, stats = verify_file(path, against=rundir, missing_ok=True)
        problems += p
        originals += stats["originals"]
    if problems:
        raise PackError(f"{rundir}: nothing deleted:\n  " + "\n  ".join(problems))
    covered = {os.path.abspath(p) for p in originals}
    left = [f for f in inv.all_files() if os.path.abspath(os.path.join(rundir, f)) not in covered]
    if left:
        raise PackError(f"{rundir}: nothing deleted: {len(left)} original files are in no packed file "
                        f"(pack them first): {left[:5]}")
    originals = list(dict.fromkeys(originals))      # manifest order: the run files (RunInfo) go last
    n_bytes = sum(os.path.getsize(p) for p in originals)
    log(f"{len(originals)} original files ({n_bytes / 1e9:.3f} GB) match the manifests of {len(files)} packed files")
    stay = sorted(_native_names(rundir) - {os.path.relpath(p, rundir) for p in originals} - {"CATALOGS"})
    if stay:
        log(f"note: {len(stay)} entries are in no packed file and stay untouched: {stay[:5]}")
    if not yes:
        log("dry run: pass --yes to delete them")
        return 0
    for path in files:
        _fsync(path)
    _fsync(os.path.abspath(rundir))
    for p in originals:
        os.remove(p)
    cat = os.path.join(rundir, "CATALOGS")
    if os.path.isdir(cat) and not os.listdir(cat):
        os.rmdir(cat)
    _fsync(os.path.abspath(rundir))
    log(f"deleted {len(originals)} files")
    return len(originals)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="python -m glamio.pack", description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("pack", help="pack a run directory")
    p.add_argument("rundir")
    p.add_argument("--out", help="directory for the packed files (default: RUNDIR)")
    p.add_argument("--box", type=int, help="box number (default: from the Run<box> directory name)")
    p.add_argument("--steps", type=int, nargs="+", help="pack only these output steps")
    v = sub.add_parser("verify", help="verify packed files")
    v.add_argument("paths", nargs="+", help="packed files or directories holding them")
    v.add_argument("--against", help="run directory with the original files to compare as well (while they exist)")
    u = sub.add_parser("unpack", help="rebuild the original files")
    u.add_argument("paths", nargs="+")
    u.add_argument("--out", required=True)
    r = sub.add_parser("prune", help="delete the original files of a packed run directory")
    r.add_argument("rundir")
    r.add_argument("--yes", action="store_true", help="really delete (default: report only)")
    s = sub.add_parser("status", help="print the state of a run directory")
    s.add_argument("rundir")
    a = ap.parse_args(argv)
    signal.signal(signal.SIGTERM, lambda signum, frame: sys.exit(128 + signum))   # time limit: clean up like an error
    try:
        if a.cmd == "pack":
            pack(a.rundir, a.out, a.box, a.steps)
        elif a.cmd == "verify":
            files = _packed_files(a.paths)
            bad = 0 if files else 1
            if not files:
                print(f"FAIL no packed files in {a.paths}")
            for path in files:
                problems, stats = verify_file(path, against=a.against)
                bad += bool(problems)
                print(("FAIL " if problems else "OK   ") + f"{os.path.basename(path)}: {stats['entries']} entries, "
                      f"{stats['bytes'] / 1e9:.3f} GB" + (f", {stats['originals_checked']} originals compared" if a.against else "")
                      + (f", {len(stats['unaligned'])} unaligned datasets" if stats["unaligned"] else ""))
                for q in problems:
                    print("     " + q)
            for d in a.paths:
                if os.path.isdir(d):
                    problems = _completeness(d)
                    bad += bool(problems)
                    print(("FAIL " if problems else "OK   ") + f"{d}: completeness")
                    for q in problems:
                        print("     " + q)
            return 1 if bad else 0
        elif a.cmd == "unpack":
            print(f"{unpack(a.paths, a.out)} files written")
        elif a.cmd == "prune":
            prune(a.rundir, a.yes)
        else:
            state, examples = status(a.rundir)
            print(f"status: {state}")
            for e in examples:
                print("     " + e)
            return 1 if state in ("inconsistent", "busy") else 0
    except PackError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
