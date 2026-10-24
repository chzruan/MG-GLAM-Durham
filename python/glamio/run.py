"""Paths of the products of one output step in a run directory.

GLAM writes every product of step ``s`` and box ``b`` as
``<product>.<s:04d>.<b:04d>.*`` in the run directory (``Run1/``), the halo
catalogue under ``CATALOGS/`` and the snapshot as ``PMcrd.<s:04d>.DAT``.
These helpers build those names so that callers (halocat, analysis scripts)
do not hand-assemble them, and list which products a step actually has.
A packed run (:mod:`glamio.packed`) keeps the products of step ``s`` in
``Output.<s:04d>.<b:04d>.h5`` and the schedule in ``RunInfo.<b:04d>.h5``;
the helpers return paths inside those files when the native files are absent.
"""

from __future__ import annotations

import os
import re

from . import packed

_RUNINFO_RE = re.compile(r"^RunInfo\.\d{4}\.h5$")

#: Products with a ``.bin`` + ``.meta`` pair in the run directory.
SIDECAR_PRODUCTS = ("LensPlanes", "deltam", "mom_x", "mom_y", "mom_z", "Subsample", "fR_halo")


def _native_path(r: str, product: str, step: int, box: int) -> str:
    tag = f"{step:04d}.{box:04d}"
    if product in SIDECAR_PRODUCTS:
        return os.path.join(r, f"{product}.{tag}")
    if product == "catalogue":
        return os.path.join(r, "CATALOGS", f"CatshortV.{tag}")
    if product in ("PowerDM.R", "PowerDM.log"):
        return os.path.join(r, f"{product}.{tag}.dat")
    if product == "schedule":
        return os.path.join(r, "OutputSchedule.dat")
    if product == "snapshot":
        return os.path.join(r, f"PMcrd.{step:04d}.DAT")
    raise ValueError(f"unknown product {product!r}")


def _native_complete(product: str, path: str) -> bool:
    """True when the native form of ``product`` is complete (both files of a pair; a catalogue in either form)."""
    if product in SIDECAR_PRODUCTS:
        return os.path.isfile(path + ".bin") and os.path.isfile(path + ".meta")
    if product == "catalogue":
        return ((os.path.isfile(path + ".bin") and os.path.isfile(path + ".meta"))
                or os.path.isfile(path + ".DAT") or os.path.isfile(path + ".dat"))
    return os.path.isfile(path)


def _runinfo(r: str, box: int) -> str:
    """``RunInfo.<box>.h5`` in ``r``, else the only ``RunInfo.*.h5`` there (the native schedule does not depend on the box)."""
    f = os.path.join(r, packed.runinfo_file_name(box))
    if os.path.isfile(f) or not os.path.isdir(r):
        return f
    found = sorted(n for n in os.listdir(r) if _RUNINFO_RE.match(n))
    if len(found) > 1:
        raise ValueError(f"{r} holds several RunInfo files {found}; pass the box")
    return os.path.join(r, found[0]) if found else f


def _packed_path(r: str, product: str, step: int, box: int) -> str | None:
    """Where ``product`` lives in a packed run (``None`` for the snapshot, which is never packed)."""
    if product == "schedule":
        return os.path.join(_runinfo(r, box), "schedule")
    if product == "snapshot":
        return None
    group = {"catalogue": packed.CATALOGUE_GROUP, "PowerDM.R": "PowerDM/R", "PowerDM.log": "PowerDM/log"}.get(product, product)
    return os.path.join(r, packed.output_file_name(step, box), group)


def product_path(rundir: str | os.PathLike, product: str, step: int, box: int = 1) -> str:
    """The path (stem for the sidecar products) of ``product`` at ``step`` in ``rundir``.

    ``product`` is one of :data:`SIDECAR_PRODUCTS` (stem of the ``.bin`` +
    ``.meta`` pair), ``'catalogue'`` (stem of ``CATALOGS/CatshortV.<step>.<box>``,
    for :func:`glamio.read_catalogue`), ``'PowerDM.R'`` / ``'PowerDM.log'``
    (the ``.dat`` file), ``'schedule'`` (``OutputSchedule.dat``) or
    ``'snapshot'`` (``PMcrd.<step>.DAT``).  The native path is returned when
    the native product is complete; otherwise, when the run is packed and holds
    the product, the path inside the packed file (``Output.<step>.<box>.h5/<group>``,
    ``RunInfo.<box>.h5/schedule``, where a run with a single ``RunInfo`` file
    is found whatever ``box`` says, as the native schedule is); otherwise the
    native path, so that a reader names the missing file.  Without h5py the
    packed files are not looked into.
    """
    r = os.fspath(rundir)
    native = _native_path(r, product, step, box)
    if _native_complete(product, native):
        return native
    inside = _packed_path(r, product, step, box)
    try:
        if inside is not None and packed.exists(inside):
            return inside
    except ImportError:
        pass
    return native


def run_file(rundir: str | os.PathLike, name: str, box: int = 1) -> str:
    """Path of a run-level file of ``rundir`` (``OutputSchedule.dat``, ``BDM.config``, ``timing.log``, ``ic2pm.log``, ...).

    The native file when it exists; otherwise, in a packed run, its copy in
    ``RunInfo.<box>.h5`` (``schedule`` for ``OutputSchedule.dat``,
    ``BDM.config``, ``logs/<name>`` for the others; a run with a single
    RunInfo file is found whatever ``box`` says); otherwise the native path.
    Read it with :func:`glamio.read_text`, which takes either form.
    """
    r = os.fspath(rundir)
    native = os.path.join(r, name)
    if os.path.isfile(native):
        return native
    obj = {"OutputSchedule.dat": "schedule", "BDM.config": "BDM.config"}.get(name, f"logs/{name}")
    try:
        inside = os.path.join(_runinfo(r, box), obj)
        if packed.exists(inside):
            return inside
    except ImportError:
        pass
    return native


def products_at(rundir: str | os.PathLike, step: int, box: int = 1) -> dict[str, str]:
    """``{product: path}`` of the products present for ``step`` in ``rundir`` (sidecar products need both files; the catalogue either form).

    Products absent natively but present in ``Output.<step>.<box>.h5`` are
    listed with their path inside that file (not looked for without h5py).
    """
    r = os.fspath(rundir)
    out: dict[str, str] = {}
    for name in SIDECAR_PRODUCTS:
        stem = _native_path(r, name, step, box)
        if os.path.isfile(stem + ".bin") and os.path.isfile(stem + ".meta"):
            out[name] = stem
    stem = _native_path(r, "catalogue", step, box)
    if (os.path.isfile(stem + ".bin") and os.path.isfile(stem + ".meta")) or os.path.isfile(stem + ".DAT"):
        out["catalogue"] = stem
    for name in ("PowerDM.R", "PowerDM.log", "snapshot"):
        p = _native_path(r, name, step, box)
        if os.path.isfile(p):
            out[name] = p
    h5 = os.path.join(r, packed.output_file_name(step, box))
    absent = [name for name in SIDECAR_PRODUCTS + ("catalogue", "PowerDM.R", "PowerDM.log") if name not in out]
    if absent and os.path.isfile(h5):
        try:
            groups = packed.objects_with_data(h5)
        except ImportError:
            return out
        for name in absent:
            inside = _packed_path(r, name, step, box)
            if inside[len(h5) + 1:] in groups:
                out[name] = inside
    return out
