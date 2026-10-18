"""Paths of the products of one output step in a run directory.

GLAM writes every product of step ``s`` and box ``b`` as
``<product>.<s:04d>.<b:04d>.*`` in the run directory (``Run1/``), the halo
catalogue under ``CATALOGS/`` and the snapshot as ``PMcrd.<s:04d>.DAT``.
These helpers build those names so that callers (halocat, analysis scripts)
do not hand-assemble them, and list which products a step actually has.
"""

from __future__ import annotations

import os

#: Products with a ``.bin`` + ``.meta`` pair in the run directory.
SIDECAR_PRODUCTS = ("LensPlanes", "deltam", "mom_x", "mom_y", "mom_z", "Subsample", "fR_halo")


def product_path(rundir: str | os.PathLike, product: str, step: int, box: int = 1) -> str:
    """The path (stem for the sidecar products) of ``product`` at ``step`` in ``rundir``.

    ``product`` is one of :data:`SIDECAR_PRODUCTS` (stem of the ``.bin`` +
    ``.meta`` pair), ``'catalogue'`` (stem of ``CATALOGS/CatshortV.<step>.<box>``,
    for :func:`glamio.read_catalogue`), ``'PowerDM.R'`` / ``'PowerDM.log'``
    (the ``.dat`` file), ``'schedule'`` (``OutputSchedule.dat``) or
    ``'snapshot'`` (``PMcrd.<step>.DAT``).
    """
    r = os.fspath(rundir)
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


def products_at(rundir: str | os.PathLike, step: int, box: int = 1) -> dict[str, str]:
    """``{product: path}`` of the products present for ``step`` in ``rundir`` (sidecar products need both files; the catalogue either form)."""
    out: dict[str, str] = {}
    for name in SIDECAR_PRODUCTS:
        stem = product_path(rundir, name, step, box)
        if os.path.isfile(stem + ".bin") and os.path.isfile(stem + ".meta"):
            out[name] = stem
    stem = product_path(rundir, "catalogue", step, box)
    if (os.path.isfile(stem + ".bin") and os.path.isfile(stem + ".meta")) or os.path.isfile(stem + ".DAT"):
        out["catalogue"] = stem
    for name in ("PowerDM.R", "PowerDM.log", "snapshot"):
        p = product_path(rundir, name, step, box)
        if os.path.isfile(p):
            out[name] = p
    return out
