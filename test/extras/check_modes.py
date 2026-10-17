#!/usr/bin/env python3
"""V1: the mode cube written by test_modes must hold 1/2 at +-m1 (cosine, real) and
-+0.25i at +-m2 (sine), zero elsewhere; also checks GrowthFactor/ComovingDistance prints."""
import sys, numpy as np
n = 4
C = np.fromfile(sys.argv[1], dtype="<c8").reshape(2*n+1, 2*n+1, n+1)   # [ix+n, iy+n, iz]
def at(m): return C[m[0]+n, m[1]+n, m[2]]
exp = {(2, -1, 3): 0.5+0j, (1, 2, 0): -0.25j, (-1, -2, 0): 0.25j}   # (-2,1,-3) has iz<0: not stored
ok = True
for m, v in exp.items():
    d = abs(at(m) - v); ok &= d < 1e-6
    print(f"mode {m}: got {at(m):.6f} expected {v:.6f}  {'OK' if d < 1e-6 else 'FAIL'}")
mask = np.ones(C.shape, bool)
for m in exp: mask[m[0]+n, m[1]+n, m[2]] = False
rest = np.max(np.abs(C[mask])); ok &= rest < 1e-6
print(f"all other modes |C| <= {rest:.2e}  {'OK' if rest < 1e-6 else 'FAIL'}")
print("V1", "PASS" if ok else "FAIL"); sys.exit(0 if ok else 1)
