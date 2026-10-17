#!/usr/bin/env python3
"""After ic2pm + one PMP2main step with Subsample M = 1: particle k of the subsample (index order)
must sit at the lattice site of Gadget ID k+1 under the lattice formula of the header tag.
usage: check_ic2pm_order.py <Run1 dir> <expected tag IDx|IDz> [Np=16] [Box=16]"""
import glob, os, sys
import numpy as np
run, tag = sys.argv[1], sys.argv[2]
Np = int(sys.argv[3]) if len(sys.argv) > 3 else 16; Box = float(sys.argv[4]) if len(sys.argv) > 4 else 16.0
fb = sorted(glob.glob(os.path.join(run, 'Subsample.*.0001.bin')))[0]; fm = fb[:-4] + '.meta'
meta = {}
for l in open(fm):
    if '=' in l: k, v = l.split('=', 1); meta[k.strip()] = v.strip()
nsel = int(meta['n_sel'].split()[0]) if 'n_sel' in meta else None
raw = np.fromfile(fb, dtype='<i2')
n = raw.size // 6; assert raw.size == 6 * n, raw.size
q = raw.reshape(6, n)[:3].astype(np.int64)
x = (q + 32768 + 0.5) * Box / 65536.0                 # x_decode of the meta
fails = []
def check(ok, msg):
    print(('PASS ' if ok else 'FAIL ') + msg)
    if not ok: fails.append(msg)
check(n == Np ** 3, 'subsample holds all %d particles (n = %d)' % (Np ** 3, n))
check(tag in meta.get('header', ''), 'header carries the expected tag %s: "%s"' % (tag, meta.get('header')))
ID = np.arange(1, n + 1) - 1
if tag == 'IDx':
    ix, iy, iz = ID % Np, (ID // Np) % Np, ID // (Np * Np)
else:
    iz, iy, ix = ID % Np, (ID // Np) % Np, ID // (Np * Np)
lat = np.stack([ix, iy, iz], 0) * (Box / Np)
d = x - lat; d -= Box * np.round(d / Box)
dmax = np.abs(d).max()
check(dmax < 0.5 * Box / Np, 'every particle k sits at the lattice site of ID k+1 (%s): max |x - q| = %.3f, half cell %.3f'
      % (tag, dmax, 0.5 * Box / Np))
# a wrong formula must fail: swap x and z
lat2 = np.stack([iz, iy, ix], 0) * (Box / Np); d2 = x - lat2; d2 -= Box * np.round(d2 / Box)
check(np.abs(d2).max() > 0.5 * Box / Np, 'the other lattice order is rejected (max |x - q| = %.3f)' % np.abs(d2).max())
sys.exit(1 if fails else 0)
