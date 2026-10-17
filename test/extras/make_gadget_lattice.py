#!/usr/bin/env python3
"""Write a tiny synthetic Gadget-1 IC (two files) for the ic2pm ingestion test.

usage: make_gadget_lattice.py <out_prefix> <mode> [Np=16] [Box=16]
modes: z32seq  int32 IDs, sequential in file order, z-fastest lattice (2LPTic-like)
       x64shuf int64 IDs, x-fastest lattice (FML-like), particles in random file order
       dup     as z32seq but one ID duplicated (ic2pm must refuse the file set)
Cosmology matches test/PkTable.dat (Om 0.3089, OmL 0.6911, h = sqrt((0.022299+0.119445)/0.3089)).
"""
import struct, sys
import numpy as np

out, mode = sys.argv[1], sys.argv[2]
Np = int(sys.argv[3]) if len(sys.argv) > 3 else 16
Box = float(sys.argv[4]) if len(sys.argv) > 4 else 16.0
N = Np ** 3; dx = Box / Np
rng = np.random.default_rng(12345)
ix, iy, iz = np.meshgrid(np.arange(Np), np.arange(Np), np.arange(Np), indexing='ij')
ix, iy, iz = ix.ravel(), iy.ravel(), iz.ravel()           # iz fastest in this flattening
if mode in ('z32seq', 'dup'):
    ids = 1 + iz + Np * (iy + Np * ix); order = np.arange(N); idt = '<i4'
else:
    ids = 1 + ix + Np * (iy + Np * iz); order = rng.permutation(N); idt = '<i8'
assert np.array_equal(np.sort(ids), np.arange(1, N + 1))
if mode == 'dup':
    ids = ids.copy(); ids[6] = ids[5]
q = np.stack([ix, iy, iz], 1) * dx
pos = (q + 0.1 * dx * rng.uniform(-1, 1, q.shape)) % Box
vel = rng.normal(0, 50.0, q.shape) / np.sqrt(0.02)      # Gadget u = v_pec / sqrt(a)
pos, vel, ids = pos[order].astype('<f4'), vel[order].astype('<f4'), ids[order].astype(idt)
h = np.sqrt((0.022299 + 0.119445) / 0.3089)
nf = 2; cuts = [0, N // 2, N]
for f in range(nf):
    sl = slice(cuts[f], cuts[f + 1]); n = cuts[f + 1] - cuts[f]
    hdr = struct.pack('<6I', 0, n, 0, 0, 0, 0) + struct.pack('<6d', 0, 1.0, 0, 0, 0, 0) \
        + struct.pack('<2d', 0.02, 49.0) + struct.pack('<2i', 0, 0) + struct.pack('<6I', 0, N, 0, 0, 0, 0) \
        + struct.pack('<2i', 0, nf) + struct.pack('<4d', Box, 0.3089, 0.6911, h) + struct.pack('<2i', 0, 0) \
        + struct.pack('<6I', 0, 0, 0, 0, 0, 0)
    hdr += b'\0' * (256 - len(hdr))
    def rec(b): return struct.pack('<i', len(b)) + b + struct.pack('<i', len(b))
    with open('%s.%d' % (out, f), 'wb') as fh:
        fh.write(rec(hdr) + rec(pos[sl].tobytes()) + rec(vel[sl].tobytes()) + rec(ids[sl].tobytes()))
print('wrote %s.0/.1: %d particles, mode %s, %s IDs' % (out, N, mode, idt))
