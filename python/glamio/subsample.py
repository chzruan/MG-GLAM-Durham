"""Reader for the quantised particle subsample ``Subsample.<step>.<box>.bin`` + ``.meta``.

Selection: particle ID = array index - 1 (0-based, lattice order; see the
sidecar's ``lattice`` line and :attr:`Subsample.lattice`).  ID is kept when
``splitmix64(ID xor SALT) mod M == 0`` with ``SALT = 0x5EED00DE6C0A2026``;
``M = 1`` keeps every particle.  IDs are not stored: the k-th stored particle
is the k-th ID that passes the rule, so :meth:`Subsample.ids` recovers them
by re-evaluating the rule over ``0 .. N_p - 1``.

Quantisation: positions ``q = floor(x/Box * 65536)`` with ``x`` wrapped into
``[0, Box)``, stored as int16 ``q - 32768``; decoded as ``x = (q + 0.5)
Box/65536`` (``dx_quant``, 15.6 kpc/h at ``Box = 1024``).  Velocities:
peculiar ``v [km/s] = V_code (100 Box/N_grid)/a_v`` at the velocity epoch
``a_v = a - da/2``, stored as int16 ``s = nint(v/dv_quant)`` clipped to
``[-32767, 32767]`` (``dv_quant = 0.25`` km/s; the clipped count is in the
sidecar).

File layout: raw little-endian int16, no record markers, the six arrays
``qx, qy, qz, vx, vy, vz`` of ``n_sel`` values each, contiguous, in that
order (12 bytes per particle).  The axes follow the particle ``x, y, z`` of
the halo catalogue.
"""

from __future__ import annotations

import os

import numpy as np

from .product import MetaProduct, check_size, load_array, load_meta, meta_int, resolve_pair

#: Position quanta per box side.
NQ = 65536
#: The selection salt (never changes).
SALT_DEFAULT = 0x5EED00DE6C0A2026


def splitmix64(x) -> np.ndarray:
    """The splitmix64 finaliser on uint64 with wrap-around arithmetic (as in C)."""
    x = np.asarray(x, dtype=np.uint64)
    with np.errstate(over="ignore"):
        z = x + np.uint64(0x9E3779B97F4A7C15)
        z = (z ^ (z >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
        z = (z ^ (z >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
        return z ^ (z >> np.uint64(31))


def keep_mask(ids, salt: int = SALT_DEFAULT, m: int = 1) -> np.ndarray:
    """Boolean mask of the kept IDs: ``splitmix64(ID xor SALT) mod M == 0``."""
    ids = np.asarray(ids, dtype=np.uint64)
    return (splitmix64(ids ^ np.uint64(salt)) % np.uint64(m)) == 0


def select_ids(n_p: int, salt: int = SALT_DEFAULT, m: int = 1, block: int = 1 << 24) -> np.ndarray:
    """IDs (0-based, increasing, uint64) of the kept particles among ``0 .. n_p - 1``."""
    if m == 1:
        return np.arange(n_p, dtype=np.uint64)
    out = []
    for i0 in range(0, n_p, block):
        ids = np.arange(i0, min(n_p, i0 + block), dtype=np.uint64)
        out.append(ids[keep_mask(ids, salt, m)])
    return np.concatenate(out) if out else np.zeros(0, np.uint64)


class Subsample(MetaProduct):
    """A quantised subsample: ``q[6, n_sel]`` int16 (``qx, qy, qz, vx, vy, vz``) plus the sidecar ``meta``."""

    def __init__(self, q: np.ndarray, meta: dict[str, str], path: str):
        super().__init__(meta, path)
        self.q = q

    def __len__(self) -> int:
        return int(self.q.shape[1])

    def __repr__(self) -> str:
        return f"Subsample({len(self)} particles, M={self.m}, N_p={self.n_p}, a={self.a}, path={self.path!r})"

    @property
    def n_sel(self) -> int:
        return len(self)

    @property
    def m(self) -> int | None:
        """Subsampling factor ``M`` (meta key ``M``)."""
        return self._int("M")

    @property
    def salt(self) -> int:
        """The selection salt (meta key ``SALT``, hexadecimal), :data:`SALT_DEFAULT` when absent."""
        value = self.meta.get("SALT")
        if not value:
            return SALT_DEFAULT
        return int(value.split()[0], 0)

    @property
    def dx_quant(self) -> float | None:
        """Position quantum ``Box/65536`` in Mpc/h."""
        d = self._float("dx_quant")
        if d is None and self.box is not None:
            d = self.box / NQ
        return d

    @property
    def dv_quant(self) -> float | None:
        """Velocity quantum in km/s (meta key ``dv_quant``)."""
        return self._float("dv_quant")

    @property
    def v_scale(self) -> float | None:
        """km/s per code momentum unit, ``100 Box/N_grid/a_v`` (meta key ``v_scale``)."""
        return self._float("v_scale")

    @property
    def n_clip_v(self) -> int | None:
        """Velocity components clipped to +-32767, summed over x, y, z (meta key ``n_clip_v``)."""
        return self._int("n_clip_v")

    @property
    def lattice(self) -> str | None:
        """Lattice-order tag of the particle IDs from the run header: ``'IDx'``, ``'IDz'``, ``'ID?'`` or ``None``."""
        h = self.header or ""
        for tag in ("IDx", "IDz", "ID?"):
            if tag in h:
                return tag
        return None

    def positions(self) -> np.ndarray:
        """``(n_sel, 3)`` float64 comoving positions in Mpc/h, in ``[0, Box)`` (decoded cell centres)."""
        if self.box is None:
            raise ValueError(f"Box is not in the sidecar of {self.path}")
        return (self.q[:3].T.astype(np.float64) + (32768.0 + 0.5)) * (self.box / NQ)

    def velocities(self) -> np.ndarray:
        """``(n_sel, 3)`` float64 peculiar velocities in km/s at ``a_v``."""
        dv = self.dv_quant
        if dv is None:
            raise ValueError(f"dv_quant is not in the sidecar of {self.path}")
        return self.q[3:].T.astype(np.float64) * dv

    def ids(self, block: int = 1 << 24) -> np.ndarray:
        """The 0-based particle IDs (uint64, increasing) of the stored particles, recomputed from the selection rule."""
        n_p, m = self.n_p, self.m
        if n_p is None or m is None:
            raise ValueError(f"N_p and M must be in the sidecar of {self.path} to recover the IDs")
        ids = select_ids(n_p, self.salt, m, block)
        if ids.size != len(self):
            raise ValueError(f"{self.path}: the selection rule keeps {ids.size} of {n_p} IDs but the file holds {len(self)}")
        return ids


def read_subsample(path: str | os.PathLike, mmap: bool = False) -> Subsample:
    """Read ``Subsample.<step>.<box>.bin`` + ``.meta`` from either file or the bare stem.

    The sidecar must carry ``product = Subsample``, ``n_sel``, ``Box`` and
    ``dv_quant``; the ``.bin`` must hold exactly ``6 n_sel`` int16 values.
    ``mmap=True`` memory-maps the file (``12 N_p/M`` bytes).
    """
    bin_path, meta_path, _ = resolve_pair(path)
    meta = load_meta(meta_path, expect_product="Subsample", required=("n_sel", "Box", "dv_quant"))
    n_sel = meta_int(meta, "n_sel")
    if n_sel is None or n_sel < 0:
        raise ValueError(f"{meta_path}: n_sel = {meta.get('n_sel')!r} is not usable")
    check_size(bin_path, 12 * n_sel, f"6 x n_sel={n_sel} int16")
    arr = load_array(bin_path, "<i2", 6 * n_sel, mmap=mmap).reshape(6, n_sel)
    return Subsample(arr, meta, bin_path)
