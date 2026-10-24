"""Reader for the power-spectrum text files ``PowerDM.R/log/Z.<step>.<box>.dat``.

Layout: line 1 the 45-character run header, line 2 the epoch line
``Aexpn = ... Step= ... Redshift= ... Ngrid= ... kNyq= ...`` (``kNyq = pi
N_grid/Box`` in h/Mpc), line 3 the caption ``bin k/Mpch N Preal`` (plus
``P0_redshift P2_redshift`` for the ``Z`` files), then one row per bin:
the bin index, the mean ``k`` of the bin in h/Mpc, the number of modes and
``P`` in ``(Mpc/h)^3``.  ``R`` files use linear bins of width ``2 pi/Box``
(``N_grid/2`` rows), ``log`` files logarithmic bins (empty bins skipped).
The estimate is the raw CIC density on the ``N_grid`` mesh: the CIC window
is not deconvolved and the shot noise not subtracted.  ``Z`` files are only
written by the redshift-space branch (never in the DEGRACE runs); their
``Ngrid`` is the halved mesh while ``kNyq`` is the full-grid value.
"""

from __future__ import annotations

import math
import os
import re

import numpy as np

from .packed import split_packed
from .product import read_text

_EPOCH_RE = re.compile(r"Aexpn\s*=\s*([0-9.]+)\s+Step=\s*(\d+)\s+Redshift=\s*([-0-9.]+)\s+Ngrid=\s*(\d+)\s+kNyq=\s*([0-9.]+)")
_COLUMN_NAMES = {"bin": "bin", "k/Mpch": "k", "N": "nmodes", "Preal": "pk",
                 "P0_redshift": "p0", "P2_redshift": "p2"}


class PowerDM:
    """A power-spectrum table: ``data[row, column]`` float64 with named ``columns``."""

    def __init__(self, data: np.ndarray, columns: list[str], header: str, a: float, z: float,
                 step: int, ngrid: int, k_nyq: float, binning: str | None, path: str, box: float | None = None):
        self.data = data
        self.columns = columns
        self.header = header
        self.a = a
        self.z = z
        self.step = step
        self.ngrid = ngrid
        self.k_nyq = k_nyq
        self.binning = binning
        self.path = path
        self.box = box
        if box is not None and binning != "redshift":
            self.k_nyq = math.pi * ngrid / box              # exact, instead of the two-decimal print

    def __len__(self) -> int:
        return int(self.data.shape[0])

    def __getitem__(self, name: str) -> np.ndarray:
        if name not in self.columns:
            raise KeyError(f"{name!r} is not a column of {self.path}: {self.columns}")
        return self.data[:, self.columns.index(name)]

    def __repr__(self) -> str:
        return f"PowerDM({len(self)} bins, {self.binning}, a={self.a}, step={self.step}, columns={self.columns}, path={self.path!r})"

    @property
    def n_grid(self) -> int:
        """Alias of ``ngrid`` (the other readers' name)."""
        return self.ngrid

    @property
    def k(self) -> np.ndarray:
        """Mean wavenumber of each bin, h/Mpc."""
        return self["k"]

    @property
    def nmodes(self) -> np.ndarray:
        """Number of modes in each bin."""
        return self["nmodes"]

    @property
    def pk(self) -> np.ndarray:
        """Real-space ``P(k)`` in ``(Mpc/h)^3`` (raw CIC)."""
        return self["pk"]

    def cic_window_jing_sq(self) -> np.ndarray:
        """Jing (2005) alias-summed 1-D CIC window ``1 - (2/3) sin^2(pi k/2 k_Nyq)`` evaluated at the bin ``|k|``.

        This is the usual quick correction of a CIC power spectrum, not a
        shell average of the 3-D window ``prod_i sinc^4`` (the two differ by
        up to 37 % near ``k_Nyq``).  ``k_nyq`` is the printed two-decimal
        value unless ``read_powerdm(..., box=...)`` was given.  Not defined
        for ``Z`` files (their ``ngrid`` and ``kNyq`` refer to different meshes).
        """
        if self.binning == "redshift":
            raise ValueError(f"{self.path}: the CIC window of a Z file is ambiguous (ngrid is the halved mesh, kNyq the full one)")
        return 1.0 - (2.0 / 3.0) * np.sin(0.5 * np.pi * self.k / self.k_nyq) ** 2

    def pk_deconvolved(self) -> np.ndarray:
        """``P(k)`` divided by :meth:`cic_window_jing_sq`."""
        return self.pk / self.cic_window_jing_sq()


def read_powerdm(path: str | os.PathLike, box: float | None = None) -> PowerDM:
    """Parse a ``PowerDM.*.dat`` file.

    ``box`` (Mpc/h) replaces the printed two-decimal ``kNyq`` by
    ``pi Ngrid/box``.  Raises ``ValueError`` when the epoch line or the
    caption is missing or a data row does not have the caption's number of
    columns.  ``path`` may be a packed table (``Output.<step>.<box>.h5/PowerDM/R``,
    see :mod:`glamio.packed`), whose verbatim text is parsed the same way.
    """
    p = os.fspath(path)
    inside = split_packed(p)
    name = os.path.basename(p) if inside is None else inside[1].replace("/", ".") + "."   # PowerDM/R -> PowerDM.R.
    return _parse_powerdm(read_text(p), p, name, box)


def _parse_powerdm(text: str, p: str, name: str, box: float | None) -> PowerDM:
    """Parse the text of a PowerDM table; ``name`` (a file name) selects the binning, ``p`` labels errors."""
    lines = text.splitlines()
    if not lines:
        raise ValueError(f"{p}: empty file")
    header = lines[0].strip()
    epoch = None
    epoch_index = None
    for i, line in enumerate(lines[1:6], 1):                  # line 0 is the free-text run header
        m = _EPOCH_RE.search(line)
        if m:
            epoch, epoch_index = m, i
            break
    if epoch is None:
        raise ValueError(f"{p}: epoch line 'Aexpn = ... kNyq= ...' not found")
    caption_index = None
    for i in range(epoch_index + 1, min(epoch_index + 4, len(lines))):
        if "bin" in lines[i].split():
            caption_index = i
            break
    if caption_index is None:
        raise ValueError(f"{p}: caption line 'bin k/Mpch N ...' not found after the epoch line")
    raw_columns = lines[caption_index].split()
    columns = [_COLUMN_NAMES.get(c, c) for c in raw_columns]
    ncol = len(columns)
    rows = []
    for lineno, line in enumerate(lines[caption_index + 1:], caption_index + 2):
        tokens = line.split()
        if not tokens:
            continue
        if len(tokens) != ncol:
            raise ValueError(f"{p}:{lineno}: {len(tokens)} columns, caption has {ncol}")
        try:
            rows.append([float(t) for t in tokens])
        except ValueError as exc:
            raise ValueError(f"{p}:{lineno}: cannot parse {line!r}") from exc
    data = np.array(rows, dtype=np.float64).reshape(-1, ncol)
    binning = "linear" if ".R." in name else "log" if ".log." in name else "redshift" if ".Z." in name else None
    return PowerDM(data, columns, header, float(epoch.group(1)), float(epoch.group(3)), int(epoch.group(2)),
                   int(epoch.group(4)), float(epoch.group(5)), binning, p, box)
