"""Readers for the PM particle snapshots ``PMcrd.<step>.DAT`` + ``PMcrs<i>.<step>.DAT``.

GLAM is built with ``-convert big_endian``, so both files are big-endian
Fortran unformatted files (a little-endian build is detected from the
record marker and handled too).

``PMcrd`` holds one sequential record of 529 bytes (4-byte length markers
around it): the 45-character run ``HEADER``, then ``AEXPN, AEXP0, AMPLT,
ASTEP`` (float32), ``ISTEP`` (int32), ``PARTW, TINTG, EKIN, EKIN1, EKIN2,
AU0, AEU0`` (float32), ``NROW, NGRID, Nrealization, Nseed`` (int32), ``Om,
OmL, hubble`` (float32), ``Nparticles`` (int64) and ``extras(100)`` (float32;
``extras(100)`` is the box side in Mpc/h).  The checkpoint pair
``PMcrd.DAT`` / ``PMcrs<i>.DAT`` (no step in the name) has the same layout.

``PMcrs<i>`` are direct-access files of fixed records ("pages") of
``NPAGE = 1024^2`` particles: each record is ``X(NPAGE), Y(NPAGE), Z(NPAGE),
VX(NPAGE), VY(NPAGE), VZ(NPAGE)`` float32 (24 MiB), no markers; only the
first ``NinPage`` entries of the last page are meaningful (the rest is
whatever the writer's buffer held).  512 records go into one file; the next
file is ``PMcrs<i+1>``.  Coordinates are in grid units ``[1, NGRID+1)``
(``x_Mpc/h = (X - 1) Box/NGRID``); velocities are code momenta,
``v [km/s] = V (100 Box/NGRID)/a_v`` with the leapfrog epoch ``a_v = AEXPN -
ASTEP/2`` (the convention of the quantised subsample; the legacy catalogue
divides by ``AEXPN`` instead).

Memory: a full read costs 24 bytes per particle (``float32``); use
:meth:`Snapshot.iter_pages` to stream large snapshots page by page.
"""

from __future__ import annotations

import math
import os
import re

import numpy as np

#: Particles per direct-access record ("page") and records per ``PMcrs`` file.
PAGE = 1024 ** 2
RECORDS_PER_FILE = 512
HEADER_BYTES = 529
#: Critical density in Msun/h per (Mpc/h)^3 used for :meth:`PMHeader.particle_mass`.
RHO_CRIT = 2.7754e11

_FIELDS = (("header", "S45"), ("aexpn", "f4"), ("aexp0", "f4"), ("amplt", "f4"), ("astep", "f4"),
           ("istep", "i4"), ("partw", "f4"), ("tintg", "f4"), ("ekin", "f4"), ("ekin1", "f4"),
           ("ekin2", "f4"), ("au0", "f4"), ("aeu0", "f4"), ("nrow", "i4"), ("ngrid", "i4"),
           ("nrealization", "i4"), ("nseed", "i4"), ("om", "f4"), ("oml", "f4"), ("hubble", "f4"),
           ("nparticles", "i8"), ("extras", "f4", (100,)))
_CRD_RE = re.compile(r"^PMcrd(?:\.(\d+))?\.DAT$")


def _header_dtype(byteorder: str) -> np.dtype:
    fields = []
    for name, code, *shape in _FIELDS:
        dt = code if code.startswith("S") else byteorder + code
        fields.append((name, dt, tuple(shape[0])) if shape else (name, dt))
    dt = np.dtype(fields)
    if dt.itemsize != HEADER_BYTES:
        raise RuntimeError(f"PMcrd header dtype is {dt.itemsize} bytes, expected {HEADER_BYTES}")
    return dt


class PMHeader:
    """The ``PMcrd`` header record as attributes (``header``, ``aexpn``, ``astep``, ``istep``, ``nrow``, ``ngrid``, ``nparticles``, ``extras``, ...)."""

    def __init__(self, record: np.void, byteorder: str, path: str):
        self.byteorder = byteorder
        self.path = path
        for name, *_ in _FIELDS:
            value = record[name]
            if name == "header":
                value = value.decode("ascii", errors="replace").rstrip()
            elif name == "extras":
                value = np.array(value, dtype=np.float32)
            else:
                value = value.item()
            setattr(self, name, value)

    def __repr__(self) -> str:
        return (f"PMHeader(a={self.aexpn:.6g}, step={self.istep}, Nrow={self.nrow}, Ngrid={self.ngrid}, "
                f"Np={self.nparticles}, Box={self.box:g}, header={self.header!r})")

    @property
    def box(self) -> float:
        """Box side in Mpc/h (``extras(100)``)."""
        return float(self.extras[99])

    @property
    def z(self) -> float:
        return 1.0 / self.aexpn - 1.0

    @property
    def a_v(self) -> float:
        """Leapfrog velocity epoch ``AEXPN - ASTEP/2`` (``ASTEP`` = the step just taken; the ICs from ``ic2pm half`` follow the same rule)."""
        return self.aexpn - 0.5 * self.astep

    @property
    def cell(self) -> float:
        """Mesh cell in Mpc/h."""
        return self.box / self.ngrid

    def velocity_scale(self, epoch: str = "a_v") -> float:
        """km/s per code momentum unit: ``100 Box/NGRID`` divided by ``a_v`` (default) or by ``a`` (legacy catalogue convention)."""
        if epoch == "a_v":
            return 100.0 * self.box / self.ngrid / self.a_v
        if epoch == "a":
            return 100.0 * self.box / self.ngrid / self.aexpn
        raise ValueError("epoch must be 'a_v' or 'a'")

    def particle_mass(self) -> float:
        """``rho_crit Omega_m (Box/Nrow)^3`` in Msun/h with :data:`RHO_CRIT`; the finder's own value is in the catalogue sidecar."""
        return RHO_CRIT * self.om * (self.box / self.nrow) ** 3


def read_pm_header(path: str | os.PathLike) -> PMHeader:
    """Read a ``PMcrd`` file; the byte order is taken from the record marker (529)."""
    p = os.fspath(path)
    with open(p, "rb") as fh:
        raw = fh.read(4 + HEADER_BYTES + 4)
    if len(raw) < 4 + HEADER_BYTES + 4:
        raise ValueError(f"{p}: {len(raw)} bytes, a PMcrd header record needs {4 + HEADER_BYTES + 4}")
    if int.from_bytes(raw[:4], "big") == HEADER_BYTES:
        bo = ">"
    elif int.from_bytes(raw[:4], "little") == HEADER_BYTES:
        bo = "<"
    else:
        raise ValueError(f"{p}: first record marker {raw[:4].hex()} is not {HEADER_BYTES} in either byte order")
    end = "big" if bo == ">" else "little"
    if int.from_bytes(raw[4 + HEADER_BYTES:], end) != HEADER_BYTES:
        raise ValueError(f"{p}: trailing record marker does not match the header length")
    record = np.frombuffer(raw, dtype=_header_dtype(bo), count=1, offset=4)[0]
    return PMHeader(record, bo, p)


class Snapshot:
    """A PM snapshot: the header plus streaming access to the ``PMcrs`` particle pages.

    ``path`` is the ``PMcrd.<step>.DAT`` (or ``PMcrd.DAT``) file; the particle
    files are located beside it.  ``page`` and ``records_per_file`` are the
    writer's constants and only need changing for files from another GLAM
    build.
    """

    def __init__(self, path: str | os.PathLike, page: int = PAGE, records_per_file: int = RECORDS_PER_FILE):
        self.path = os.fspath(path)
        self.header = read_pm_header(self.path)
        m = _CRD_RE.match(os.path.basename(self.path))
        if not m:
            raise ValueError(f"{self.path}: not a PMcrd.<step>.DAT or PMcrd.DAT file name")
        self.tag = m.group(1)
        self.page = int(page)
        self.records_per_file = int(records_per_file)
        self.record_bytes = 6 * self.page * 4
        n = int(self.header.nparticles)
        self.n_pages = (n + self.page - 1) // self.page
        n_files = (self.n_pages + self.records_per_file - 1) // self.records_per_file
        rundir = os.path.dirname(os.path.abspath(self.path))
        self.files = []
        for i in range(n_files):
            name = f"PMcrs{i}.{self.tag}.DAT" if self.tag is not None else f"PMcrs{i}.DAT"
            f = os.path.join(rundir, name)
            if not os.path.isfile(f):
                raise FileNotFoundError(f"particle file {f} of {self.path} is missing")
            records = min(self.records_per_file, self.n_pages - i * self.records_per_file)
            size = os.path.getsize(f)
            if size < records * self.record_bytes or size % self.record_bytes:
                raise ValueError(f"{f}: {size} bytes, expected a multiple of {self.record_bytes} holding at least {records} records")
            self.files.append(f)

    @property
    def n_particles(self) -> int:
        return int(self.header.nparticles)

    def __repr__(self) -> str:
        return f"Snapshot({self.header!r}, {self.n_pages} pages in {len(self.files)} file(s))"

    def iter_pages(self, velocities: bool = True, units: str = "code", velocity_epoch: str = "a_v"):
        """Yield ``(x, y, z[, vx, vy, vz])`` arrays page by page (up to ``page`` particles each).

        ``units='code'``: float32 grid coordinates in ``[1, NGRID+1)`` and code
        momenta, as stored.  ``units='physical'``: float64 Mpc/h in ``[0, Box)``
        (``(X-1) Box/NGRID``, wrapped) and km/s with
        :meth:`PMHeader.velocity_scale` (``velocity_epoch`` ``'a_v'`` or ``'a'``).
        """
        if units not in ("code", "physical"):
            raise ValueError("units must be 'code' or 'physical'")
        h = self.header
        dt = np.dtype(h.byteorder + "f4")
        ncomp = 6 if velocities else 3
        remaining = self.n_particles
        for ip in range(self.n_pages):
            f = self.files[ip // self.records_per_file]
            offset = (ip % self.records_per_file) * self.record_bytes
            n_in = min(self.page, remaining)
            remaining -= n_in
            block = np.fromfile(f, dtype=dt, count=ncomp * self.page, offset=offset)
            if block.size != ncomp * self.page:
                raise ValueError(f"{f}: short read at page {ip}")
            block = block.reshape(ncomp, self.page)[:, :n_in]
            if units == "code":
                yield tuple(np.ascontiguousarray(block[i]) for i in range(ncomp))
            else:
                scale = h.box / h.ngrid
                pos = [np.mod((block[i].astype(np.float64) - 1.0) * scale, h.box) for i in range(3)]
                if velocities:
                    vs = h.velocity_scale(velocity_epoch)
                    pos += [block[i].astype(np.float64) * vs for i in range(3, 6)]
                yield tuple(pos)

    def read(self, velocities: bool = True, units: str = "code", velocity_epoch: str = "a_v",
             max_particles: int | None = None) -> tuple[np.ndarray, ...]:
        """Read every particle: ``(x, y, z[, vx, vy, vz])`` as in :meth:`iter_pages`, concatenated.

        ``max_particles`` guards against loading a snapshot larger than
        intended (``ValueError`` when ``nparticles`` exceeds it).
        """
        if max_particles is not None and self.n_particles > max_particles:
            raise ValueError(f"{self.path}: {self.n_particles} particles exceed max_particles={max_particles}")
        parts = list(self.iter_pages(velocities, units, velocity_epoch))
        ncomp = 6 if velocities else 3
        return tuple(np.concatenate([p[i] for p in parts]) if parts else np.zeros(0) for i in range(ncomp))


def read_snapshot(path: str | os.PathLike, **kwargs) -> Snapshot:
    """Open a PM snapshot from its ``PMcrd`` file (see :class:`Snapshot`)."""
    return Snapshot(path, **kwargs)
