#!/usr/bin/env python3
"""Reference decoder and checker for the GLAM particle subsample
(item I9, module Subsample in PMP2mod_subsample.f90). Pure numpy.

Usage (driven by run_subsample_test.sh):
    python3 check_subsample.py <workdir> [hash_ref.bin]

The work dir holds the files written by test_subsample.exe:
    subsample_test_particles.bin   raw XPAR, YPAR, ZPAR, VX, VY, VZ (float32 LE)
    Subsample.<step>.<box>.{bin,meta}
and optionally hash_ref.bin from hash_ref.c.

Checks (PASS/FAIL per line, exit 1 on any FAIL):
  a  numpy splitmix64 / selection == C routines on the hash_ref IDs
  b  every Subsample file: byte size, n_sel == number of IDs passing the
     rule, decoded positions within dx_quant/2 (periodic) and decoded
     velocities within dv_quant/2 of the raw arrays, clipped components
     stored as +-32767 and counted in n_clip_v(x,y,z)
  c  n_sel within 3 sigma of N_p/M (M > 1)
  d  the M = 1 file holds all N_p particles in order

The module-level functions (splitmix64, select_ids, read_meta,
read_subsample, decode) are the reference implementation for readers
of production files.
"""
import glob
import os
import sys

import numpy as np

SALT_DEFAULT = 0x5EED00DE6C0A2026
NQ = 65536


# ----------------------------------------------------------------- selection
def splitmix64(x):
    """splitmix64 finaliser on uint64 (wrap-around arithmetic, as in C)."""
    x = np.asarray(x, dtype=np.uint64)
    with np.errstate(over='ignore'):
        z = x + np.uint64(0x9E3779B97F4A7C15)
        z = (z ^ (z >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
        z = (z ^ (z >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
        return z ^ (z >> np.uint64(31))


def keep(ids, salt, m):
    """Boolean mask: splitmix64(ID xor SALT) mod M == 0."""
    ids = np.asarray(ids, dtype=np.uint64)
    return (splitmix64(ids ^ np.uint64(salt)) % np.uint64(m)) == 0


def select_ids(n_p, salt, m, block=1 << 24):
    """IDs (0-based, increasing) of the kept particles among 0 .. n_p-1."""
    out = []
    for i0 in range(0, n_p, block):
        ids = np.arange(i0, min(n_p, i0 + block), dtype=np.uint64)
        out.append(ids[keep(ids, salt, m)])
    return np.concatenate(out) if out else np.zeros(0, np.uint64)


# ----------------------------------------------------------------- file I/O
def read_meta(path):
    """key -> value string (first token of the value, comments in
    parentheses stripped); 'raw' keeps the full lines."""
    meta = {}
    with open(path) as f:
        for line in f:
            if ' = ' not in line:
                continue
            key, val = line.split(' = ', 1)
            meta[key.strip()] = val.strip()
    return meta


def mnum(meta, key, typ=float):
    tok = meta[key].split()[0]
    if typ is int:
        return int(tok, 0)
    return typ(tok)


def read_subsample(binpath, metapath=None):
    """Return (meta, qx, qy, qz, sx, sy, sz) int16 arrays of n_sel."""
    if metapath is None:
        metapath = binpath[:-4] + '.meta'
    meta = read_meta(metapath)
    n_sel = mnum(meta, 'n_sel', int)
    data = np.fromfile(binpath, dtype='<i2')
    if data.size != 6 * n_sel:
        raise ValueError('%s: %d int16 values, expected 6 x n_sel = %d' % (binpath, data.size, 6 * n_sel))
    return (meta,) + tuple(data[k * n_sel:(k + 1) * n_sel] for k in range(6))


def decode(meta, qx, qy, qz, sx, sy, sz):
    """Positions [Mpc/h] in [0, Box) and peculiar velocities [km/s]."""
    box = mnum(meta, 'Box')
    dv = mnum(meta, 'dv_quant')
    pos = [(q.astype(np.float64) + 32768.0 + 0.5) * box / NQ for q in (qx, qy, qz)]
    vel = [s.astype(np.float64) * dv for s in (sx, sy, sz)]
    return pos, vel


# ----------------------------------------------------------------- checks
class Report:
    def __init__(self):
        self.nfail = 0

    def check(self, ok, label, detail=''):
        tag = 'PASS' if ok else 'FAIL'
        if not ok:
            self.nfail += 1
        print('%s  %s%s' % (tag, label, ('  [' + detail + ']') if detail else ''))
        return ok


def check_hash_ref(rep, path):
    with open(path, 'rb') as f:
        n = int(np.fromfile(f, dtype='<i8', count=1)[0])
        ids = np.fromfile(f, dtype='<u8', count=n)
        h = np.fromfile(f, dtype='<u8', count=n)
        s100 = np.fromfile(f, dtype='u1', count=n)
        s7 = np.fromfile(f, dtype='u1', count=n)
    salt = np.uint64(SALT_DEFAULT)
    hp = splitmix64(ids ^ salt)
    rep.check(bool(np.array_equal(hp, h)), 'a  splitmix64(ID xor SALT): numpy == C on %d IDs' % n,
              'mismatches %d' % int(np.count_nonzero(hp != h)))
    k100 = keep(ids, SALT_DEFAULT, 100).astype(np.uint8)
    k7 = keep(ids, SALT_DEFAULT, 7).astype(np.uint8)
    rep.check(bool(np.array_equal(k100, s100)) and bool(np.array_equal(k7, s7)),
              'a  selection mod 100 and mod 7: numpy == C',
              'kept M=100: %d, M=7: %d' % (int(k100.sum()), int(k7.sum())))


def check_file(rep, binpath, raw, n_p_raw):
    name = os.path.basename(binpath)
    try:
        meta, qx, qy, qz, sx, sy, sz = read_subsample(binpath)
    except ValueError as e:
        rep.check(False, 'b  %s: file size' % name, str(e))
        return
    n_p = mnum(meta, 'N_p', int)
    m = mnum(meta, 'M', int)
    salt = mnum(meta, 'SALT', int)
    n_sel = mnum(meta, 'n_sel', int)
    box = mnum(meta, 'Box')
    ngrid = mnum(meta, 'N_grid', int)
    a_v = mnum(meta, 'a_v')
    dx = mnum(meta, 'dx_quant')
    dv = mnum(meta, 'dv_quant')
    rep.check(n_p == n_p_raw and salt == SALT_DEFAULT, 'b  %s: N_p and SALT as expected' % name,
              'N_p %d, SALT 0x%016X, M %d' % (n_p, salt, m))
    rep.check(abs(dx - box / NQ) <= 1e-12 * box and abs(dv - 0.25) < 1e-12,
              'b  %s: dx_quant = Box/65536, dv_quant = 0.25' % name)

    ids = select_ids(n_p, salt, m)
    rep.check(ids.size == n_sel, 'b  %s: n_sel == number of IDs passing the rule' % name,
              'n_sel %d, rule %d' % (n_sel, ids.size))
    if ids.size != n_sel:
        return
    idx = ids.astype(np.int64)

    # positions: the Fortran path in float64 (grid units from 0, wrapped, floor)
    pos, vel = decode(meta, qx, qy, qz, sx, sy, sz)
    worst = 0.0
    qexact = True
    for comp, q in zip(range(3), (qx, qy, qz)):
        s = raw[comp][idx].astype(np.float64) - 1.0
        s = np.mod(s, float(ngrid))
        x_true = s * box / ngrid                                  # [0, Box)
        d = pos[comp] - x_true
        d = (d + box / 2) % box - box / 2                         # periodic difference
        worst = max(worst, float(np.max(np.abs(d))) if d.size else 0.0)
        qref = np.minimum(np.floor(s * (NQ / ngrid)).astype(np.int64), NQ - 1) - 32768
        qexact = qexact and bool(np.array_equal(qref, q.astype(np.int64)))
    rep.check(worst <= 0.5 * dx * (1 + 1e-9), 'b  %s: decoded positions within dx_quant/2' % name,
              'max |dx| = %.3e Mpc/h, dx_quant/2 = %.3e, exact q match: %s' % (worst, 0.5 * dx, qexact))

    # velocities
    vfac = 100.0 * box / ngrid / a_v
    nclip_ref = []
    worst = 0.0
    sign_ok = True
    for comp, s in zip(range(3), (sx, sy, sz)):
        v_true = raw[3 + comp][idx].astype(np.float64) * vfac
        r = np.rint(v_true / dv)
        clipped = np.abs(r) > 32767
        nclip_ref.append(int(np.count_nonzero(clipped)))
        if np.any(~clipped):
            worst = max(worst, float(np.max(np.abs(vel[comp][~clipped] - v_true[~clipped]))))
        if np.any(clipped):
            sign_ok = sign_ok and bool(np.all(s[clipped].astype(np.int64) == np.sign(r[clipped]) * 32767))
    rep.check(worst <= 0.5 * dv * (1 + 1e-6) + 1e-6, 'b  %s: decoded velocities within dv_quant/2' % name,
              'max |dv| = %.4e km/s, dv_quant/2 = %.4e' % (worst, 0.5 * dv))
    nclip_meta = [mnum(meta, 'n_clip_v' + c, int) for c in 'xyz']
    rep.check(nclip_meta == nclip_ref and mnum(meta, 'n_clip_v', int) == sum(nclip_ref) and sign_ok,
              'b  %s: clipped components stored as +-32767 and counted' % name,
              'meta %s, reference %s' % (nclip_meta, nclip_ref))

    # c: count within 3 sigma of the binomial expectation
    if m > 1:
        mean = n_p / m
        sig = np.sqrt(n_p / m * (1 - 1.0 / m))
        rep.check(abs(n_sel - mean) <= 3 * sig, 'c  %s: n_sel within 3 sigma of N_p/M' % name,
                  'n_sel %d, N_p/M %.1f, sigma %.1f, pull %+.2f' % (n_sel, mean, sig, (n_sel - mean) / sig))
    # d: full mode
    if m == 1:
        rep.check(n_sel == n_p and bool(np.array_equal(idx, np.arange(n_p))),
                  'd  %s: M = 1 holds all N_p particles in ID order' % name, 'n_sel %d' % n_sel)


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(2)
    work = sys.argv[1]
    rep = Report()
    if len(sys.argv) > 2 and os.path.exists(sys.argv[2]):
        check_hash_ref(rep, sys.argv[2])
    else:
        rep.check(False, 'a  hash reference file missing')

    rawpath = os.path.join(work, 'subsample_test_particles.bin')
    rawall = np.fromfile(rawpath, dtype='<f4')
    n_p_raw = rawall.size // 6
    raw = [rawall[k * n_p_raw:(k + 1) * n_p_raw] for k in range(6)]
    files = sorted(glob.glob(os.path.join(work, 'Subsample.*.bin')))
    rep.check(len(files) > 0, 'b  Subsample files present', '%d files' % len(files))
    rep.check(not os.path.exists(os.path.join(work, 'Subsample.0010.0001.bin')) and
              not os.path.exists(os.path.join(work, 'Subsample.0010.0001.meta')),
              'b  M = 0 (step 10) wrote nothing')
    ms = []
    for f in files:
        check_file(rep, f, raw, n_p_raw)
        ms.append(mnum(read_meta(f[:-4] + '.meta'), 'M', int))
    rep.check(1 in ms, 'd  an M = 1 file was tested', 'M values %s' % ms)
    print('%s: %d check(s) failed' % ('FAIL' if rep.nfail else 'PASS', rep.nfail))
    sys.exit(1 if rep.nfail else 0)


if __name__ == '__main__':
    main()
