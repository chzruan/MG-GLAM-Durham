"""Shared base class and helpers for the sidecar-described products.

Every binary product GLAM writes (``LensPlanes``, ``deltam``/``mom_*``
cubes, ``Subsample``, ``fR_halo``, the binary halo catalogue) comes as a
headerless ``.bin`` file beside a ``key = value`` text sidecar ``.meta``
(:mod:`glamio.meta`).  :class:`MetaProduct` gives the readers one place for
the scalar properties every sidecar carries (``a``, ``z``, ``a_v``, ``Box``,
``step``, ...); :func:`resolve_pair`, :func:`load_meta`, :func:`check_size`
and :func:`load_array` implement the common "find the pair, validate the
sidecar, check the byte count, read or memory-map" sequence.
"""

from __future__ import annotations

import os

import numpy as np

from .meta import meta_float, read_meta


class MetaProduct:
    """A product with a ``meta`` dict (sidecar or parsed header) and a source ``path``.

    The scalar properties are parsed from ``meta`` on access and are ``None``
    when the corresponding key is missing; the readers validate the keys they
    need, so a ``None`` only appears for optional keys.
    """

    def __init__(self, meta: dict[str, str], path: str):
        self.meta = meta
        self.path = path

    def _float(self, key: str) -> float | None:
        return meta_float(self.meta, key)

    def _int(self, key: str) -> int | None:
        value = meta_float(self.meta, key)
        return None if value is None else int(round(value))

    def _str(self, key: str) -> str | None:
        value = self.meta.get(key)
        return None if value is None else value

    @property
    def product(self) -> str | None:
        """Product name: the first token of the meta key ``product`` (a trailing comment is dropped)."""
        value = self.meta.get("product")
        if value is None:
            return None
        tokens = value.split()
        return tokens[0] if tokens else ""

    @property
    def a(self) -> float | None:
        """Expansion factor of the output (meta key ``a``)."""
        return self._float("a")

    @property
    def z(self) -> float | None:
        """Redshift: meta key ``z``, else ``1/a - 1``, else ``None``."""
        z = self._float("z")
        if z is None:
            a = self.a
            if a is not None and a != 0.0:
                z = 1.0 / a - 1.0
        return z

    @property
    def a_v(self) -> float | None:
        """Velocity epoch ``a_v = a - da/2`` of the output (meta key ``a_v``)."""
        return self._float("a_v")

    @property
    def a_x(self) -> float | None:
        """Position epoch (meta key ``a_x``, equal to ``a``)."""
        return self._float("a_x")

    @property
    def box(self) -> float | None:
        """Box side in Mpc/h (meta key ``Box``)."""
        return self._float("Box")

    @property
    def step(self) -> int | None:
        """Time-step number of the output (meta key ``step``)."""
        return self._int("step")

    @property
    def realization(self) -> int | None:
        """Box / realization number (meta key ``realization``)."""
        return self._int("realization")

    @property
    def n_grid(self) -> int | None:
        """PM mesh size per side (meta key ``N_grid``)."""
        return self._int("N_grid")

    @property
    def n_row(self) -> int | None:
        """Particles per side (meta key ``N_row``)."""
        return self._int("N_row")

    @property
    def n_p(self) -> int | None:
        """Total number of particles (meta key ``N_p``)."""
        return self._int("N_p")

    @property
    def omega_m(self) -> float | None:
        return self._float("Omega_m")

    @property
    def omega_l(self) -> float | None:
        return self._float("Omega_L")

    @property
    def hubble(self) -> float | None:
        return self._float("hubble")

    @property
    def header(self) -> str | None:
        """The 45-character GLAM run header string (meta key ``header``)."""
        return self._str("header")


def resolve_pair(path: str | os.PathLike, bin_suffix: str = ".bin", meta_suffix: str = ".meta") -> tuple[str, str, str]:
    """Return ``(bin_path, meta_path, stem)`` for a ``.bin``, a ``.meta`` or a bare stem."""
    p = os.fspath(path)
    if os.path.isdir(p):
        raise IsADirectoryError(f"{p} is a directory; pass the product's .bin, .meta or stem")
    for suffix in (bin_suffix, meta_suffix):
        if p.endswith(suffix):
            stem = p[: -len(suffix)]
            break
    else:
        stem = p
    return stem + bin_suffix, stem + meta_suffix, stem


def load_meta(meta_path: str, expect_product: str | tuple[str, ...] | None = None,
              required: tuple[str, ...] = ()) -> dict[str, str]:
    """Read a sidecar and check its ``product`` line and the keys a reader needs.

    Raises ``FileNotFoundError`` when the sidecar is missing and ``ValueError``
    when ``product`` is not one of ``expect_product`` (a name or a tuple of
    names; ``None`` skips the check) or a ``required`` key is absent.
    """
    if not os.path.isfile(meta_path):
        raise FileNotFoundError(f"sidecar not found: {meta_path}")
    meta = read_meta(meta_path)
    if expect_product is not None:
        expected = (expect_product,) if isinstance(expect_product, str) else tuple(expect_product)
        value = meta.get("product", "")
        name = value.split()[0] if value.split() else ""
        if name not in expected:
            raise ValueError(f"{meta_path}: product = {value!r}, expected one of {expected}")
    missing = [key for key in required if key not in meta]
    if missing:
        raise ValueError(f"{meta_path}: missing sidecar key(s) {missing}")
    return meta


def meta_int(meta: dict[str, str], key: str, default: int | None = None) -> int | None:
    """Integer form of :func:`glamio.meta.meta_float` (first token, rounded)."""
    value = meta_float(meta, key)
    return default if value is None else int(round(value))


def check_size(path: str, expected_bytes: int, what: str) -> None:
    """Raise ``ValueError`` unless ``path`` holds exactly ``expected_bytes``."""
    if not os.path.isfile(path):
        raise FileNotFoundError(f"binary file not found: {path}")
    size = os.path.getsize(path)
    if size != expected_bytes:
        raise ValueError(f"{path}: {size} bytes, expected {expected_bytes} ({what})")


def load_array(path: str, dtype: np.dtype | str, count: int, mmap: bool = False) -> np.ndarray:
    """Read (or memory-map, read-only) ``count`` items of ``dtype`` from ``path``; ``count = 0`` gives an empty array either way."""
    dt = np.dtype(dtype)
    if count == 0:
        return np.zeros(0, dtype=dt)              # np.memmap cannot map an empty file; an empty product is legitimate
    if mmap:
        arr = np.memmap(path, dtype=dt, mode="r", shape=(count,))
    else:
        arr = np.fromfile(path, dtype=dt, count=count)
    if arr.size != count:
        raise ValueError(f"{path}: read {arr.size} items of {dt}, expected {count}")
    return arr
