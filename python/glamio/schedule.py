"""Reader for ``OutputSchedule.dat``: one line per analysis output.

``PMP2main`` writes the file at start-up (and identically at a restart) in
the run directory.  Comment lines start with ``#``; one of them records the
cosmology used for the comoving distance (``Om=``, ``OmL=``, ``c/H0 =
2997.92458 Mpc/h``).  Each data line holds ``step a z chi exact class``:

* ``step``: the time step that names the snapshot and product files;
* ``a``, ``z``: the expansion factor and redshift of the output;
* ``chi``: comoving distance in Mpc/h, flat matter + Lambda;
* ``exact``: 1 when the step is pinned on a requested redshift (``#outputs < 0``
  in ``Init.dat``), 0 for the nearest scheduled step;
* ``class``: 1 = full output (halo finder, extras), 2 = station (lens planes
  and P(k) only).
"""

from __future__ import annotations

import os
import re

import numpy as np

from .product import read_text

#: Row layout returned by :func:`read_schedule`.
SCHEDULE_DTYPE = np.dtype([("step", "<i4"), ("a", "<f8"), ("z", "<f8"), ("chi", "<f8"),
                           ("exact", "<i4"), ("class", "<i4")])

_COSMO_RE = re.compile(r"Om=\s*([0-9.eEdD+-]+)\s+OmL=\s*([0-9.eEdD+-]+)")
_C_H0_RE = re.compile(r"c/H0\s*=\s*([0-9.eEdD+-]+)")


class Schedule:
    """The output schedule: structured ``rows`` (:data:`SCHEDULE_DTYPE`) plus the header comments."""

    def __init__(self, rows: np.ndarray, comments: list[str], path: str):
        self.rows = rows
        self.comments = comments
        self.path = path

    def __len__(self) -> int:
        return len(self.rows)

    def __getitem__(self, key):
        """``sched['a']`` gives a column, ``sched[i]`` the i-th row."""
        return self.rows[key]

    def __repr__(self) -> str:
        return f"Schedule({len(self)} outputs, steps {self.steps.min() if len(self) else None}..{self.steps.max() if len(self) else None}, path={self.path!r})"

    @property
    def steps(self) -> np.ndarray:
        return self.rows["step"]

    @property
    def a(self) -> np.ndarray:
        return self.rows["a"]

    @property
    def z(self) -> np.ndarray:
        return self.rows["z"]

    @property
    def chi(self) -> np.ndarray:
        return self.rows["chi"]

    @property
    def exact(self) -> np.ndarray:
        return self.rows["exact"]

    @property
    def classes(self) -> np.ndarray:
        return self.rows["class"]

    def _cosmo(self) -> tuple[float, float] | None:
        for line in self.comments:
            m = _COSMO_RE.search(line)
            if m:
                return float(m.group(1).replace("D", "E")), float(m.group(2).replace("D", "E"))
        return None

    @property
    def omega_m(self) -> float | None:
        """Omega_m of the comoving-distance cosmology (from the header comment)."""
        c = self._cosmo()
        return None if c is None else c[0]

    @property
    def omega_l(self) -> float | None:
        c = self._cosmo()
        return None if c is None else c[1]

    @property
    def c_over_h0(self) -> float | None:
        """The ``c/H0`` the writer used (Mpc/h), from the header comment."""
        for line in self.comments:
            m = _C_H0_RE.search(line)
            if m:
                return float(m.group(1).replace("D", "E"))
        return None

    def by_step(self) -> dict[int, np.void]:
        """``{step: row}`` for every output."""
        return {int(r["step"]): r for r in self.rows}

    def row(self, step: int) -> np.void:
        """The row of output ``step``; ``KeyError`` when the step is not an output."""
        hit = np.nonzero(self.rows["step"] == step)[0]
        if hit.size == 0:
            raise KeyError(f"step {step} is not in {self.path}")
        return self.rows[hit[0]]

    def row_at_z(self, z: float, atol_a: float = 1e-4) -> np.void:
        """The row whose ``a`` is nearest to ``1/(1+z)``; ``KeyError`` when the nearest is farther than ``atol_a``."""
        a = 1.0 / (1.0 + z)
        i = int(np.argmin(np.abs(self.rows["a"] - a)))
        if len(self.rows) == 0 or abs(self.rows["a"][i] - a) > atol_a:
            raise KeyError(f"no output within {atol_a} in a of z = {z} in {self.path}")
        return self.rows[i]

    def primary(self) -> np.ndarray:
        """The class-1 rows (full outputs: halo finder and extras); classes other than 1 and 2 are not written today."""
        return self.rows[self.rows["class"] == 1]

    def stations(self) -> np.ndarray:
        """The class-2 rows (lens planes and P(k) only)."""
        return self.rows[self.rows["class"] == 2]


def read_schedule(path: str | os.PathLike) -> Schedule:
    """Parse ``OutputSchedule.dat``.

    Raises ``ValueError`` for a data line without exactly six tokens or when
    the steps are not strictly increasing (the writer emits them in step
    order).  ``path`` may be the packed copy (``RunInfo.<box>.h5/schedule``,
    see :mod:`glamio.packed`).
    """
    p = os.fspath(path)
    return _parse_schedule(read_text(p), p)


def _parse_schedule(text: str, p: str) -> Schedule:
    """Parse the text of ``OutputSchedule.dat``; ``p`` labels errors and becomes ``Schedule.path``."""
    comments: list[str] = []
    rows: list[tuple] = []
    for lineno, line in enumerate(text.split("\n"), 1):
        s = line.strip()
        if not s:
            continue
        if s.startswith("#"):
            comments.append(s)
            continue
        tokens = s.split()
        if len(tokens) != 6:
            raise ValueError(f"{p}:{lineno}: expected 6 columns (step a z chi exact class), got {len(tokens)}")
        try:
            rows.append((int(tokens[0]), float(tokens[1]), float(tokens[2]), float(tokens[3]),
                         int(tokens[4]), int(tokens[5])))
        except ValueError as exc:
            raise ValueError(f"{p}:{lineno}: cannot parse {s!r}") from exc
    arr = np.array(rows, dtype=SCHEDULE_DTYPE)
    if len(arr) > 1 and not np.all(np.diff(arr["step"]) > 0):
        raise ValueError(f"{p}: steps are not strictly increasing")
    if len(arr) and (not np.all(np.isin(arr["exact"], (0, 1))) or np.any(arr["class"] < 1) or np.any(arr["a"] <= 0)):
        raise ValueError(f"{p}: exact must be 0/1, class >= 1 and a > 0 on every row")
    return Schedule(arr, comments, p)
