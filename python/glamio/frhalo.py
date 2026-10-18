"""Reader for ``fR_halo.<step>.<box>.bin`` + ``.meta``: f_R at the BDM halo centres (f(R) runs).

One float32 pair per catalogue row, in the row order of the ``CatshortV``
catalogue of the same step: ``fR_centre`` (f_R interpolated with CIC at the
halo centre) and ``fR_mean`` (the mean of f_R over the mesh cells whose
centres lie within the catalogue radius ``Rvir`` of the centre; equal to
``fR_centre`` when the radius is below half a cell).  The field is the
converged scalar field of the pre-step epoch ``a_old = a - da``:
``f_R(x) = fR_bg(a_old) u(x)^(n+1)`` with ``fR_bg(a) = fR0 (R0/R(a))^(n+1)``,
``R(a) = 3 (Omega_m/a^3 + 4 Omega_L)``.  The sidecar names the text and the
binary catalogue (``catalogue``, ``catalogue_bin``); which of them exists
depends on ``CatalogueFormat`` in ``BDM.config``.
"""

from __future__ import annotations

import os

import numpy as np

from .catalogue import Catalogue, read_catalogue
from .product import MetaProduct, check_size, load_array, load_meta, meta_int, resolve_pair

#: Record layout of the fR_halo file.
FRHALO_DTYPE = np.dtype([("fR_centre", "<f4"), ("fR_mean", "<f4")])


class FRHalo(MetaProduct):
    """f_R at halo centres: structured ``data`` (:data:`FRHALO_DTYPE`) plus the sidecar ``meta``."""

    def __init__(self, data: np.ndarray, meta: dict[str, str], path: str):
        super().__init__(meta, path)
        self.data = data

    def __len__(self) -> int:
        return len(self.data)

    def __getitem__(self, name: str) -> np.ndarray:
        return self.data[name]

    def __repr__(self) -> str:
        return f"FRHalo({len(self)} rows, fR0={self.fR0}, a={self.a}, path={self.path!r})"

    @property
    def rows(self) -> int:
        return len(self)

    @property
    def fR0(self) -> float | None:
        return self._float("fR0")

    @property
    def n_fr(self) -> int | None:
        """The f(R) index ``n`` of the Hu-Sawicki model (meta key ``n``)."""
        return self._int("n")

    @property
    def a_old(self) -> float | None:
        """Pre-step epoch at which the scalar field was solved (meta key ``a_old``)."""
        return self._float("a_old")

    @property
    def fR_bg_old(self) -> float | None:
        """Background ``f_R`` at ``a_old`` (meta key ``fR_bg(a_old)``)."""
        return self._float("fR_bg(a_old)")

    def catalogue_paths(self) -> dict[str, str]:
        """``{'text': ..., 'binary': ...}`` paths of the catalogue files named by the sidecar, relative to the run directory of the file."""
        rundir = os.path.dirname(os.path.abspath(self.path))
        out = {}
        for key, name in (("catalogue", "text"), ("catalogue_bin", "binary")):
            value = self.meta.get(key, "")
            if value.split():
                out[name] = os.path.join(rundir, value.split()[0])
        return out

    def catalogue(self, prefer: str = "binary") -> Catalogue:
        """Read the catalogue of the same step (``prefer`` ``'binary'`` or ``'text'``; falls back to the other form).

        Raises ``FileNotFoundError`` when neither exists and ``ValueError``
        when its row count differs from this file's.
        """
        paths = self.catalogue_paths()
        order = ("binary", "text") if prefer == "binary" else ("text", "binary")
        for form in order:
            p = paths.get(form)
            if p and os.path.isfile(p):
                cat = read_catalogue(p)
                if len(cat) != len(self):
                    raise ValueError(f"{self.path}: {len(self)} rows but the catalogue {p} has {len(cat)}")
                return cat
        raise FileNotFoundError(f"no catalogue found for {self.path}: {paths}")


def read_frhalo(path: str | os.PathLike, mmap: bool = False) -> FRHalo:
    """Read ``fR_halo.<step>.<box>.bin`` + ``.meta`` from either file or the bare stem.

    The sidecar must carry ``product = fR_halo`` and ``rows``; the ``.bin``
    must hold exactly ``rows`` float32 pairs.
    """
    bin_path, meta_path, _ = resolve_pair(path)
    meta = load_meta(meta_path, expect_product="fR_halo", required=("rows",))
    rows = meta_int(meta, "rows")
    if rows is None or rows < 0:
        raise ValueError(f"{meta_path}: rows = {meta.get('rows')!r} is not usable")
    check_size(bin_path, 8 * rows, f"rows={rows} float32 pairs")
    arr = load_array(bin_path, FRHALO_DTYPE, rows, mmap=mmap)
    return FRHalo(arr, meta, bin_path)
