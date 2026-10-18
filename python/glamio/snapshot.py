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
Two words carry no information in GLAM's own snapshots: ``AEXP0`` is an
uninitialised local of the writer (only the native IC generator sets it)
and ``Nseed`` is the running random-number state, not the IC seed; ``AU0``
is always 0.

``PMcrs<i>`` are direct-access files of fixed records ("pages") of
``NPAGE = 1024^2`` particles: each record is ``X(NPAGE), Y(NPAGE), Z(NPAGE),
VX(NPAGE), VY(NPAGE), VZ(NPAGE)`` float32 (24 MiB), no markers; only the
first ``NinPage`` entries of the last page are meaningful (the rest is
whatever the writer's buffer held).  GLAM's writer puts 512 records in one
file before opening ``PMcrs<i+1>``; the native IC generator ``PMP2start``
uses 256.  The reader takes the split from the size of ``PMcrs0`` (a
writer fills every file but the last), so both are read without options.
A direct-access file is never truncated by Fortran, so a rerun with fewer
particles in the same directory leaves stale trailing records; the reader
ignores anything beyond the ``n_pages`` it needs.

Coordinates are in grid units ``[1, NGRID+1)`` (``x_Mpc/h = (X - 1)
Box/NGRID``; a particle at exactly ``NGRID + 1`` is the periodic image of
0); velocities are code momenta, ``v [km/s] = V (100 Box/NGRID)/a_v`` with
the leapfrog epoch ``a_v = AEXPN - ASTEP/2`` (``ASTEP`` = the step just
taken; this is the convention of the quantised subsample, and the legacy
catalogue divides by ``AEXPN`` instead).  Exception: an IC file converted
by ``ic2pm`` in its legacy ``sync`` mode holds velocities at ``AEXPN``
itself; such a file is recognised by its header (``ic2pm: ...`` without the
``v at a-da/2`` tag, and ``ISTEP = 0``) and :attr:`PMHeader.a_v` then
returns ``AEXPN``.

Memory: :meth:`Snapshot.read` costs 24 bytes per particle in code units
(float32) and 48 in physical units (float64), with no second copy; it
refuses more than ``max_particles`` (default ``2**27``, 3.2 / 6.4 GB) so
that a production snapshot (1024^3 = 26 / 52 GB, 2048^3 = 206 / 412 GB) is
not loaded by accident on a login node.  :meth:`Snapshot.iter_pages`
streams one 24 MiB page at a time.
"""

from __future__ import annotations

import os
import re

import numpy as np

#: Particles per direct-access record ("page") and GLAM's records per ``PMcrs`` file.
PAGE = 1024 ** 2
RECORDS_PER_FILE = 512
HEADER_BYTES = 529
#: The finder's critical density, Msun/h per (Mpc/h)^3 (``Dscale`` in PMP2linker.f90), for :meth:`PMHeader.particle_mass`.
RHO_CRIT_GLAM = 2.774e11
#: Default cap of :meth:`Snapshot.read` in particles (2**27 = 3.2 GB in code units).
MAX_PARTICLES_DEFAULT = 2 ** 27

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
    """The ``PMcrd`` header record as attributes (``header``, ``aexpn``, ``astep``, ``istep``, ``nrow``, ``ngrid``, ``nparticles``, ``extras``, ...).

    ``aexp0`` and ``nseed`` are exposed because they are in the record, but
    carry no information in GLAM's snapshots (see the module docstring).
    """

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
    def ic_velocity_epoch(self) -> str | None:
        """For an ``ic2pm`` IC file (``ISTEP = 0``): ``'half'`` (velocities at ``a - da/2``) or ``'sync'`` (at ``a``); ``None`` otherwise."""
        if self.istep == 0 and self.header.startswith("ic2pm"):
            return "half" if "v at a-da/2" in self.header else "sync"
        return None

    @property
    def a_v(self) -> float:
        """Velocity epoch: ``AEXPN - ASTEP/2`` (``ASTEP`` = the step just taken; also the ``ic2pm half`` and native ICs), or ``AEXPN`` for an ``ic2pm sync`` IC."""
        if self.ic_velocity_epoch == "sync":
            return self.aexpn
        return self.aexpn - 0.5 * self.astep

    @property
    def cell(self) -> float:
        """Mesh cell in Mpc/h."""
        return self.box / self.ngrid

    def velocity_scale(self, epoch: str = "a_v") -> float:
        """km/s per code momentum unit: ``100 Box/NGRID`` divided by :attr:`a_v` (default) or by ``AEXPN`` (``epoch='a'``, the legacy catalogue convention)."""
        if epoch == "a_v":
            return 100.0 * self.box / self.ngrid / self.a_v
        if epoch == "a":
            return 100.0 * self.box / self.ngrid / self.aexpn
        raise ValueError("epoch must be 'a_v' or 'a'")

    def particle_mass(self) -> float:
        """``RHO_CRIT_GLAM Omega_m (Box/Nrow)^3`` in Msun/h, the finder's own value (the catalogue sidecar's ``particle_mass``)."""
        return RHO_CRIT_GLAM * self.om * (self.box / self.nrow) ** 3


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
    files are located beside it.  ``page`` is the writer's record length in
    particles (``1024**2`` for every GLAM build); ``records_per_file`` is
    taken from the size of ``PMcrs0`` unless given.
    """

    def __init__(self, path: str | os.PathLike, page: int = PAGE, records_per_file: int | None = None):
        self.path = os.fspath(path)
        self.header = read_pm_header(self.path)
        m = _CRD_RE.match(os.path.basename(self.path))
        if not m:
            raise ValueError(f"{self.path}: not a PMcrd.<step>.DAT or PMcrd.DAT file name")
        self.tag = m.group(1)
        self.page = int(page)
        self.record_bytes = 6 * self.page * 4
        n = int(self.header.nparticles)
        self.n_pages = (n + self.page - 1) // self.page
        rundir = os.path.dirname(os.path.abspath(self.path))

        def name(i: int) -> str:
            return os.path.join(rundir, f"PMcrs{i}.{self.tag}.DAT" if self.tag is not None else f"PMcrs{i}.DAT")

        first = name(0)
        if not os.path.isfile(first):
            raise FileNotFoundError(f"particle file {first} of {self.path} is missing")
        size0 = os.path.getsize(first)
        if size0 % self.record_bytes:
            raise ValueError(f"{first}: {size0} bytes is not a multiple of the {self.record_bytes}-byte page record")
        if records_per_file is None:
            records0 = size0 // self.record_bytes
            records_per_file = self.n_pages if records0 >= self.n_pages else records0
            if records_per_file < 1:
                raise ValueError(f"{first}: empty particle file for {n} particles")
        self.records_per_file = int(records_per_file)
        n_files = max(1, (self.n_pages + self.records_per_file - 1) // self.records_per_file)
        self.files = []
        for i in range(n_files):
            f = name(i)
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
             max_particles: int | None = MAX_PARTICLES_DEFAULT) -> tuple[np.ndarray, ...]:
        """Read every particle: ``(x, y, z[, vx, vy, vz])`` as in :meth:`iter_pages`, each a full-length array.

        The arrays are allocated once (24 bytes per particle in code units,
        48 in physical units) and filled page by page.  ``max_particles``
        (default :data:`MAX_PARTICLES_DEFAULT`) refuses larger snapshots with
        ``ValueError``; pass ``None`` to lift the cap, or stream with
        :meth:`iter_pages`.
        """
        n = self.n_particles
        if max_particles is not None and n > max_particles:
            raise ValueError(f"{self.path}: {n} particles exceed max_particles={max_particles} "
                             f"({24 if units == 'code' else 48} bytes each); stream with iter_pages() or pass max_particles=None")
        ncomp = 6 if velocities else 3
        dtype = np.float32 if units == "code" else np.float64
        out = tuple(np.empty(n, dtype=dtype) for _ in range(ncomp))
        pos = 0
        for parts in self.iter_pages(velocities, units, velocity_epoch):
            m = len(parts[0])
            for c in range(ncomp):
                out[c][pos:pos + m] = parts[c]
            pos += m
        if pos != n:
            raise ValueError(f"{self.path}: read {pos} particles, header says {n}")
        return out


def read_snapshot(path: str | os.PathLike, **kwargs) -> Snapshot:
    """Open a PM snapshot from its ``PMcrd`` file (see :class:`Snapshot`)."""
    return Snapshot(path, **kwargs)
