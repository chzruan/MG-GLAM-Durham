"""Reader for the lens-plane product ``LensPlanes.<step>.<box>.bin`` + ``.meta``.

At every output step GLAM projects the CIC overdensity on the ``N_grid^3``
mesh along each axis in ``N_slab`` slabs of depth ``Box/N_slab``.  Every
plane is the mean of ``1 + delta`` over the slab column, block-averaged to
``N_lens x N_lens`` pixels (``Sigma/Sigma_bar``).  Only the mean over the
``N_slab`` planes of one axis is 1 by construction (mass conservation);
a single plane averages to ``1 + <delta>_slab``, the mean overdensity of its
slab, which is physical signal (a few percent on a 150 Mpc/h box, 1e-3 on
the 1024 Mpc/h pilots).  The CIC window is not deconvolved.

File layout: raw little-endian float32, no record markers, ``3 N_slab``
planes in the order X1..XN, Y1..YN, Z1..ZN (projection axis, slab).  Each
plane is ``(N_lens, N_lens)`` in Fortran column-major order, indexed
``(y, z)`` for the X planes, ``(x, z)`` for the Y planes and ``(x, y)`` for
the Z planes; pixel ``(i, j)`` (1-based) covers ``[(i-1), i) Box/N_lens``
along the first plane axis and slab ``s`` covers ``[(s-1), s) depth`` along
the projection axis.  The axes follow the particle ``x, y, z`` of the halo
catalogue.

:func:`read_lensplanes` returns the planes as a numpy array of shape
``(3 N_slab, N_lens, N_lens)`` in the natural reshape of the file, i.e.
indexed ``[plane, j, i]`` with ``i`` the FIRST (fastest) plane axis.
:meth:`LensPlanes.plane` returns one plane transposed to ``[i, j]`` so that
``planes.plane('x', 1)[iy, iz]`` reads as the sidecar describes.
"""

from __future__ import annotations

import os

import numpy as np

from .product import MetaProduct, check_size, load_array, load_meta, meta_int, resolve_pair

AXES = ("x", "y", "z")


def _axis_index(axis) -> int:
    if isinstance(axis, (bool, np.bool_)):
        raise ValueError(f"axis must be one of {AXES} or 0..2, got {axis!r}")
    if isinstance(axis, str):
        if axis.lower() not in AXES:
            raise ValueError(f"axis must be one of {AXES} or 0..2, got {axis!r}")
        return AXES.index(axis.lower())
    if not isinstance(axis, (int, np.integer)) or axis not in (0, 1, 2):
        raise ValueError(f"axis must be one of {AXES} or 0..2, got {axis!r}")
    return int(axis)


class LensPlanes(MetaProduct):
    """Lens planes of one output: ``planes[plane, j, i]`` (float32) plus the sidecar ``meta``."""

    def __init__(self, planes: np.ndarray, meta: dict[str, str], path: str):
        super().__init__(meta, path)
        self.planes = planes

    def __repr__(self) -> str:
        return f"LensPlanes({self.n_planes} planes of {self.n_lens}x{self.n_lens}, {self.n_slab} slabs per axis, a={self.a}, path={self.path!r})"

    @property
    def n_planes(self) -> int:
        return int(self.planes.shape[0])

    @property
    def n_slab(self) -> int:
        """Slabs per projection axis."""
        return self.n_planes // 3

    @property
    def n_lens(self) -> int:
        """Pixels per plane side."""
        return int(self.planes.shape[1])

    @property
    def depth(self) -> float | None:
        """Slab depth in Mpc/h (meta key ``depth``, else ``Box/N_slab``)."""
        d = self._float("depth")
        if d is None and self.box is not None:
            d = self.box / self.n_slab
        return d

    @property
    def pixel_size(self) -> float | None:
        """Pixel side in Mpc/h (``Box/N_lens``)."""
        return None if self.box is None else self.box / self.n_lens

    def slabs(self, axis) -> np.ndarray:
        """All planes projected along ``axis`` (``'x'``/``'y'``/``'z'`` or 0..2): ``(N_slab, N_lens, N_lens)``, ``[slab, j, i]``."""
        iax = _axis_index(axis)
        return self.planes[iax * self.n_slab:(iax + 1) * self.n_slab]

    def plane(self, axis, slab: int) -> np.ndarray:
        """One plane as ``[i, j]`` (first plane axis first), ``slab`` 1-based as in the sidecar (X1..XN).

        For ``axis='x'`` the result is indexed ``[iy, iz]``, for ``'y'``
        ``[ix, iz]``, for ``'z'`` ``[ix, iy]``; index ``i`` covers
        ``[i, i+1) Box/N_lens`` (0-based here).
        """
        iax = _axis_index(axis)
        if not 1 <= slab <= self.n_slab:
            raise IndexError(f"slab must be 1..{self.n_slab}, got {slab}")
        return self.planes[iax * self.n_slab + slab - 1].T

    def axis_means(self) -> np.ndarray:
        """Mean of ``Sigma/Sigma_bar`` over the planes of each axis, ``(3,)``; each must be 1 (loads every plane once)."""
        return np.array([float(np.mean(self.slabs(iax), dtype=np.float64)) for iax in range(3)])

    def plane_means(self) -> np.ndarray:
        """Mean of each plane, ``(3, N_slab)`` indexed ``[axis, slab-1]``: ``1 + <delta>_slab``, the numbers the Fortran logs per plane."""
        ns = self.n_slab
        return np.array([[float(np.mean(self.planes[iax * ns + s], dtype=np.float64)) for s in range(ns)] for iax in range(3)])

    def pixel_edges(self) -> np.ndarray:
        """Pixel edges along a plane axis, ``(N_lens + 1,)`` in Mpc/h."""
        if self.box is None:
            raise ValueError(f"Box is not in the sidecar of {self.path}")
        return np.linspace(0.0, self.box, self.n_lens + 1)

    def slab_edges(self) -> np.ndarray:
        """Slab edges along the projection axis, ``(N_slab + 1,)`` in Mpc/h."""
        if self.box is None:
            raise ValueError(f"Box is not in the sidecar of {self.path}")
        return np.linspace(0.0, self.box, self.n_slab + 1)


def read_lensplanes(path: str | os.PathLike, mmap: bool = False) -> LensPlanes:
    """Read ``LensPlanes.<step>.<box>.bin`` + ``.meta`` from either file or the bare stem.

    The sidecar must carry ``product = LensPlanes``, ``N_slab``, ``N_lens`` and
    ``Box``; the ``.bin`` must hold exactly ``3 N_slab N_lens^2`` float32
    values.  ``mmap=True`` memory-maps the file (``3 N_slab N_lens^2`` float32:
    805 MB for ``N_slab = 4``, ``N_lens = 4096``) instead of loading it.
    """
    bin_path, meta_path, _ = resolve_pair(path)
    meta = load_meta(meta_path, expect_product="LensPlanes", required=("N_slab", "N_lens", "Box"))
    ns, nl = meta_int(meta, "N_slab"), meta_int(meta, "N_lens")
    if ns is None or nl is None or ns <= 0 or nl <= 0:
        raise ValueError(f"{meta_path}: N_slab = {meta.get('N_slab')!r}, N_lens = {meta.get('N_lens')!r} are not usable")
    n_planes = meta_int(meta, "N_planes")
    if n_planes is not None and n_planes != 3 * ns:
        raise ValueError(f"{meta_path}: N_planes = {n_planes} but 3 N_slab = {3 * ns}")
    count = 3 * ns * nl * nl
    check_size(bin_path, 4 * count, f"3 x N_slab={ns} x N_lens^2={nl * nl} float32")
    arr = load_array(bin_path, "<f4", count, mmap=mmap).reshape(3 * ns, nl, nl)
    return LensPlanes(arr, meta, bin_path)
