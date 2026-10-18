"""Reader for the low-k Fourier mode cubes: ``deltam``, ``mom_x/y/z`` and ``*_deltaL``.

GLAM dumps, for a field ``f(x)`` on an ``N_dep^3`` mesh (CIC, not
interlaced), the modes ``delta_k = N^-3 sum_x f(x) exp(-i k.x)`` with
``k = (2 pi/Box)(ix, iy, iz)`` for the cube ``|ix|, |iy| <= n``,
``0 <= iz <= n`` (``n = floor(k_c Box/2 pi)``, the Hermitian-redundant
``iz = 0`` modes are kept, ``k = 0`` is set to zero).  ``P(k) = Box^3
<|delta_k|^2>``; grid node ``i`` sits at ``x = (i-1) Box/N``.  Products:

* ``deltam``: the matter overdensity on the ``N_grid`` mesh (``k_c`` from
  ``kc_delta``);
* ``mom_x``, ``mom_y``, ``mom_z``: ``(1 + delta) v_i`` with ``v`` the
  peculiar velocity in km/s at the velocity epoch ``a_v``, deposited on the
  ``N_dep`` momentum mesh (``k_c`` from ``kc_mom``);
* ``<prefix>_deltaL``: the linear density modes of the ICs produced by the
  companion tools, same layout (only ``n`` and ``Box`` are guaranteed in
  its sidecar).

File layout: raw little-endian complex64 (re, im float32), no record
markers, Fortran ``C(0:n, -n:n, -n:n)`` with the first index ``iz``, so
``np.fromfile(dtype=complex64).reshape(2n+1, 2n+1, n+1)`` is indexed
``[ix+n, iy+n, iz]``.  The axes follow the particle ``x, y, z`` of the halo
catalogue.
"""

from __future__ import annotations

import math
import os

import numpy as np

from .product import MetaProduct, check_size, load_array, load_meta, meta_int, resolve_pair


class Modes(MetaProduct):
    """A mode cube: ``data[ix+n, iy+n, iz]`` complex64 plus the sidecar ``meta``."""

    def __init__(self, data: np.ndarray, meta: dict[str, str], path: str):
        super().__init__(meta, path)
        self.data = data

    def __repr__(self) -> str:
        return f"Modes({self.product}, n={self.n}, k_c={self.k_c}, N_dep={self.n_dep}, a={self.a}, path={self.path!r})"

    @property
    def n(self) -> int:
        """Half-width of the cube in mode units (``|ix|, |iy|, iz <= n``)."""
        return int(self.data.shape[2]) - 1

    @property
    def k_c(self) -> float | None:
        """Requested cut-off wavenumber in h/Mpc (meta key ``k_c``)."""
        return self._float("k_c")

    @property
    def n_dep(self) -> int | None:
        """Mesh the field was deposited on (meta key ``N_dep``)."""
        return self._int("N_dep")

    @property
    def d_code(self) -> float | None:
        """Linear growth factor ``D(a)/D(1)`` of the code's background (meta key ``D_code(a)``)."""
        return self._float("D_code(a)")

    @property
    def k_fund(self) -> float:
        """Fundamental wavenumber ``2 pi/Box`` in h/Mpc."""
        if self.box is None:
            raise ValueError(f"Box is not in the sidecar of {self.path}")
        return 2.0 * math.pi / self.box

    def mode(self, ix: int, iy: int, iz: int) -> complex:
        """``delta_k`` at integer mode ``(ix, iy, iz)``; ``iz < 0`` uses the Hermitian symmetry ``C(-k) = conj(C(k))``."""
        n = self.n
        if max(abs(ix), abs(iy), abs(iz)) > n:
            raise IndexError(f"mode ({ix}, {iy}, {iz}) outside the cube |i| <= {n}")
        if iz >= 0:
            return complex(self.data[ix + n, iy + n, iz])
        return complex(np.conj(self.data[-ix + n, -iy + n, -iz]))

    def index_grids(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Integer mode indices ``(ix, iy, iz)`` broadcast to the cube shape (``indexing='ij'``)."""
        n = self.n
        return np.meshgrid(np.arange(-n, n + 1), np.arange(-n, n + 1), np.arange(0, n + 1), indexing="ij")

    def wavenumbers(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """``(kx, ky, kz)`` in h/Mpc, cube shape."""
        kf = self.k_fund
        return tuple(kf * g for g in self.index_grids())

    def k_magnitude(self) -> np.ndarray:
        """``|k|`` in h/Mpc, cube shape."""
        ix, iy, iz = self.index_grids()
        return self.k_fund * np.sqrt(ix**2 + iy**2 + iz**2)

    def power(self, nbins: int | None = None) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Shell-averaged ``P(k) = Box^3 <|delta_k|^2>`` in bins ``b = nint(|m|)``, ``b = 1..nbins``.

        The average runs over the full sphere: every ``iz > 0`` mode counts
        twice (for its Hermitian partner), ``iz = 0`` modes once.  Returns
        ``(k_mean, pk, weight)`` where ``k_mean`` is the weighted mean ``|k|``
        of the bin and ``weight`` the number of modes counted; empty bins give
        ``nan``.  ``nbins`` defaults to ``n``.  With the ``deltam`` cube this
        matches ``PowerDM.R`` bin by bin (CIC window not deconvolved).
        """
        if nbins is None:
            nbins = self.n
        ix, iy, iz = self.index_grids()
        mm = np.sqrt(ix**2 + iy**2 + iz**2)
        w = np.where(iz > 0, 2.0, 1.0)
        b = np.rint(mm).astype(int)
        p2 = np.abs(self.data.astype(np.complex128)) ** 2 * self.box**3
        kmean = np.full(nbins, np.nan)
        pk = np.full(nbins, np.nan)
        weight = np.zeros(nbins)
        for k in range(1, nbins + 1):
            sel = b == k
            ws = np.sum(w[sel])
            weight[k - 1] = ws
            if ws > 0:
                pk[k - 1] = np.sum(w[sel] * p2[sel]) / ws
                kmean[k - 1] = self.k_fund * np.sum(w[sel] * mm[sel]) / ws
        return kmean, pk, weight


def read_modes(path: str | os.PathLike, mmap: bool = False, expect_product: str | tuple[str, ...] | None = None) -> Modes:
    """Read a mode cube from its ``.bin``, its ``.meta`` or the bare stem.

    The sidecar must carry ``n`` (and ``Box`` for the k-space helpers); the
    ``.bin`` must hold exactly ``(2n+1)^2 (n+1)`` complex64 values.
    ``expect_product`` restricts the accepted ``product`` line.
    """
    bin_path, meta_path, _ = resolve_pair(path)
    meta = load_meta(meta_path, expect_product=expect_product, required=("n",))
    n = meta_int(meta, "n")
    if n is None or n < 0:
        raise ValueError(f"{meta_path}: n = {meta.get('n')!r} is not usable")
    count = (2 * n + 1) ** 2 * (n + 1)
    check_size(bin_path, 8 * count, f"(2n+1)^2 (n+1) complex64 with n={n}")
    arr = load_array(bin_path, "<c8", count, mmap=mmap).reshape(2 * n + 1, 2 * n + 1, n + 1)
    return Modes(arr, meta, bin_path)
