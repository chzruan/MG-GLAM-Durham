#!/usr/bin/env python3
"""fR_halo unit check: reproduce the three rows written by test_modes from the synthetic ramp."""
import sys, numpy as np
N, L = 16, 100.0
Om, OmL, a, da, fR0, n = 0.3, 0.7, 0.5, 0.01, -1e-5, 1
u = 0.5 + 0.02 * np.arange(N)                    # u(ix), ix = 0..N-1 (Fortran index ix+1), independent of y,z
U = np.broadcast_to(u[:, None, None], (N, N, N)).astype(np.float32)
a_old = a - da
R = lambda aa: 3 * (Om / aa**3 + 4 * OmL)
fRbg = fR0 * (R(1.0) / R(a_old)) ** (n + 1)
def cic(x, y, z):
    X, Y, Z = [(v * N / L) % N + 1 for v in (x, y, z)]
    I, J, K = int(X), int(Y), int(Z)
    D = [X - I, Y - J, Z - K]
    val = 0.0
    for di, wi in ((0, 1 - D[0]), (1, D[0])):
        for dj, wj in ((0, 1 - D[1]), (1, D[1])):
            for dk, wk in ((0, 1 - D[2]), (1, D[2])):
                ii, jj, kk = [(q - 1 + d) % N for q, d in ((I, di), (J, dj), (K, dk))]
                val += wi * wj * wk * float(U[ii, jj, kk]) ** (n + 1)
    return val
def mean_in_sphere(x, y, z, rvir):
    X, Y, Z = [(v * N / L) % N + 1 for v in (x, y, z)]
    r = rvir * N / L
    if r < 0.5: return cic(x, y, z)
    I, J, K = int(X), int(Y), int(Z); nr = int(r) + 1
    tot, cnt = 0.0, 0
    for kk in range(K - nr, K + nr + 1):
        for jj in range(J - nr, J + nr + 1):
            for ii in range(I - nr, I + nr + 1):
                if (ii - X)**2 + (jj - Y)**2 + (kk - Z)**2 <= r * r:
                    tot += float(U[(ii - 1) % N, (jj - 1) % N, (kk - 1) % N]) ** (n + 1); cnt += 1
    return tot / cnt if cnt else cic(x, y, z)
rows = [(10.0, 0.0, 0.0, 0.0), (10.0, 0.0, 0.0, 20.0), (99.5, 50.0, 25.0, 0.0)]
got = np.fromfile(sys.argv[1], dtype="<f4").reshape(-1, 2)
ok = got.shape[0] == 3
for (x, y, z, rv), (gc, gm) in zip(rows, got):
    ec, em = fRbg * cic(x, y, z), fRbg * mean_in_sphere(x, y, z, rv)
    dc, dm = abs(gc / ec - 1), abs(gm / em - 1)
    ok &= dc < 1e-5 and dm < 1e-5
    print(f"row ({x},{y},{z},r={rv}): centre {gc:.6e} vs {ec:.6e} ({dc:.1e}); mean {gm:.6e} vs {em:.6e} ({dm:.1e})")
meta = {k.strip(): v.strip() for k, v in (l.split("=", 1) for l in open(sys.argv[2]) if "=" in l)}
ok &= abs(float(meta["fR_bg(a_old)"]) / fRbg - 1) < 1e-6 and int(meta["rows"]) == 3
print("fR_halo", "PASS" if ok else "FAIL"); sys.exit(0 if ok else 1)
