"""Reader for the low-k Fourier mode cubes: ``deltam``, ``mom_x/y/z`` and ``*_deltaL``.

GLAM dumps, for a field ``f(x)`` on an ``N_dep^3`` mesh (CIC, not
interlaced), the modes ``delta_k = N^-3 sum_x f(x) exp(-i k.x)`` with
``k = (2 pi/Box)(ix, iy, iz)`` for the cube ``|ix|, |iy| <= n``,
``0 <= iz <= n`` (``n = floor(k_c Box/2 pi)``, clamped to ``N/2 - 1``; the
Hermitian-redundant ``iz = 0`` modes are kept, ``k = 0`` is set to zero).
``P(k) = Box^3 <|delta_k|^2>``; grid node ``i`` sits at ``x = (i-1) Box/N``.
Products:

* ``deltam``: the matter overdensity on the ``N_grid`` mesh (``k_c`` from
  ``kc_delta``);
* ``mom_x``, ``mom_y``, ``mom_z``: ``(1 + delta) v_i`` with ``v`` the
  peculiar velocity in km/s at the velocity epoch ``a_v``, deposited on the
  ``N_dep`` momentum mesh (``k_c`` from ``kc_mom``); the reader applies no
  factor;
* ``<prefix>_deltaL``: the linear density modes of the ICs produced by the
  companion tools (``product = deltaL``), same layout and normalisation,
  divided by the generator's growth factor (its sidecar carries ``n``,
  ``k_c``, ``N_mesh``, ``N_dep``, ``Box``, ``D_gen(a_ini)``; no ``D_code(a)``).

File layout: raw little-endian complex64 (re, im float32), no record
markers, Fortran ``C(0:n, -n:n, -n:n)`` with the first index ``iz``, so
``np.fromfile(dtype=complex64).reshape(2n+1, 2n+1, n+1)`` is indexed
``[ix+n, iy+n, iz]``.  The axes follow the particle ``x, y, z`` of the halo
catalogue.

Memory: the production ``deltam`` cube (``k_c = 1.5``, ``Box = 1024``,
``n = 244``) holds 58.6 M modes (0.47 GB).  :meth:`Modes.power` works slab
by slab and needs only a few MB beyond the data; :meth:`Modes.k_magnitude`
and the dense form of :meth:`Modes.wavenumbers` allocate one or three
float64 cubes (0.47 GB each) and are meant for small cubes or sparse use.
"""

from __future__ import annotations

import math
import os

import numpy as np

from .product import MetaProduct, check_size, load_array, load_meta, meta_int, resolve_pair

#: Products written in this layout by GLAM and the companion tools.
MODE_PRODUCTS = ("deltam", "mom_x", "mom_y", "mom_z", "deltaL")


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
        """Linear growth factor ``D(a)/D(1)`` of the code's background (meta key ``D_code(a)``; ``None`` for ``deltaL`` cubes, which carry ``D_gen(a_ini)``)."""
        return self._float("D_code(a)")

    @property
    def d_gen(self) -> float | None:
        """Growth factor the IC generator divided a ``deltaL`` cube by (meta key ``D_gen(a_ini)``)."""
        return self._float("D_gen(a_ini)")

    @property
    def k_fund(self) -> float:
        """Fundamental wavenumber ``2 pi/Box`` in h/Mpc (``ValueError`` when ``Box`` is not in the sidecar)."""
        if self.box is None:
            raise ValueError(f"Box is not in the sidecar of {self.path}")
        return 2.0 * math.pi / self.box

    def mode(self, ix: int, iy: int, iz: int) -> complex:
        """``delta_k`` at integer mode ``(ix, iy, iz)``; ``iz < 0`` uses the Hermitian symmetry ``C(-k) = conj(C(k))`` of a real field."""
        n = self.n
        if max(abs(ix), abs(iy), abs(iz)) > n:
            raise IndexError(f"mode ({ix}, {iy}, {iz}) outside the cube |i| <= {n}")
        if iz >= 0:
            return complex(self.data[ix + n, iy + n, iz])
        return complex(np.conj(self.data[-ix + n, -iy + n, -iz]))

    def index_grids(self, sparse: bool = False) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Integer mode indices ``(ix, iy, iz)`` broadcastable to the cube shape (``indexing='ij'``).

        ``sparse=True`` returns the three 1-D vectors shaped ``(2n+1, 1, 1)``,
        ``(1, 2n+1, 1)`` and ``(1, 1, n+1)`` (no memory); ``sparse=False`` the
        dense int64 cubes (24 bytes per mode).
        """
        n = self.n
        return np.meshgrid(np.arange(-n, n + 1), np.arange(-n, n + 1), np.arange(0, n + 1), indexing="ij", sparse=sparse)

    def wavenumbers(self, sparse: bool = False) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """``(kx, ky, kz)`` in h/Mpc, broadcastable to the cube shape (dense float64 cubes unless ``sparse=True``)."""
        kf = self.k_fund
        return tuple(kf * g for g in self.index_grids(sparse=sparse))

    def k_magnitude(self) -> np.ndarray:
        """``|k|`` in h/Mpc as a dense float64 cube (8 bytes per mode)."""
        kf = self.k_fund
        ix, iy, iz = self.index_grids(sparse=True)
        return kf * np.sqrt(ix**2 + iy**2 + iz**2)

    def power(self, nbins: int | None = None) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Shell-averaged ``P(k) = Box^3 <|delta_k|^2>`` in bins ``b = nint(|m|)``, ``b = 1..nbins``.

        The average runs over the full sphere: every ``iz > 0`` mode counts
        twice (for its Hermitian partner), ``iz = 0`` modes once.  Returns
        ``(k_mean, pk, weight)`` where ``k_mean`` is the weighted mean ``|k|``
        of the bin and ``weight`` the number of modes counted; empty bins give
        ``nan``.  Only the shells ``b <= n`` lie entirely inside the cube, so
        ``nbins`` defaults to ``n`` and larger values raise ``ValueError``
        (the truncated shells would be biased low).  With the ``deltam`` cube
        this matches ``PowerDM.R`` bin by bin (CIC window not deconvolved,
        the same ``nint`` binning as POWERfft5).  The cube is processed one
        ``ix`` slab at a time, so the temporaries stay at a few MB whatever
        the cube size (memory-mapped data is read slab by slab).
        """
        kf = self.k_fund
        n = self.n
        if nbins is None:
            nbins = n
        if nbins > n:
            raise ValueError(f"nbins = {nbins} exceeds n = {n}: shells above n are cut by the cube and would be biased")
        if nbins < 1:
            return np.zeros(0), np.zeros(0), np.zeros(0)
        box3 = float(self.box) ** 3
        iy = np.arange(-n, n + 1, dtype=np.float64)[:, None]
        iz = np.arange(0, n + 1, dtype=np.float64)[None, :]
        yz2 = iy**2 + iz**2                                   # (2n+1, n+1)
        w = np.broadcast_to(np.where(iz > 0, 2.0, 1.0), yz2.shape)
        sum_w = np.zeros(nbins + 1)
        sum_p = np.zeros(nbins + 1)
        sum_k = np.zeros(nbins + 1)
        for ix in range(-n, n + 1):
            mm = np.sqrt(ix * ix + yz2)
            b = np.rint(mm).astype(np.intp)
            sel = (b >= 1) & (b <= nbins)
            if not np.any(sel):
                continue
            slab = np.asarray(self.data[ix + n])
            p2 = slab.real.astype(np.float64) ** 2 + slab.imag.astype(np.float64) ** 2
            bb = b[sel]
            ws = w[sel]
            sum_w += np.bincount(bb, weights=ws, minlength=nbins + 1)
            sum_p += np.bincount(bb, weights=ws * p2[sel], minlength=nbins + 1)
            sum_k += np.bincount(bb, weights=ws * mm[sel], minlength=nbins + 1)
        weight = sum_w[1:]
        with np.errstate(invalid="ignore", divide="ignore"):
            pk = np.where(weight > 0, box3 * sum_p[1:] / weight, np.nan)
            kmean = np.where(weight > 0, kf * sum_k[1:] / weight, np.nan)
        return kmean, pk, weight


def read_modes(path: str | os.PathLike, mmap: bool = False,
               expect_product: str | tuple[str, ...] | None = MODE_PRODUCTS) -> Modes:
    """Read a mode cube from its ``.bin``, its ``.meta`` or the bare stem.

    The sidecar must carry ``n`` (and ``Box`` for the k-space helpers) and a
    ``product`` line among ``expect_product`` (default
    :data:`MODE_PRODUCTS`; pass ``None`` to accept any product); the
    ``.bin`` must hold exactly ``(2n+1)^2 (n+1)`` complex64 values.
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
