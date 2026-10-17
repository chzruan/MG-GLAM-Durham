#!/usr/bin/env python3
"""End-to-end check (V3) of the delta_L dump against MG-GLAM's deltam cubes.

Pipeline under test: FML generator (delta_L cube + Gadget IC) -> ic2pm (half)
-> PMP2main with `Extra outputs = 1` -> deltam.<step>.<box>.{bin,meta}.

For every deltam cube in the run directory the complex ratio
    r(k) = deltam(k) / W_CIC(k; N_dep) / delta_L(k)
is formed mode by mode.  In the linear regime r(k) = D_code(a) (real,
k-independent), where D_code(a) is the growth factor MG-GLAM wrote into the
.meta file (same Simpson integral as the generator uses for D(a_ini)).  The
particle lattice sits on the deposit-grid nodes, where the CIC response is
non-analytic in the displacement (one-sided differences), which leaves
per-mode phase scatter ~ k dx/2 but no bias in the mean to O((k dx)^2).

Asserts, per deltam cube:
  - bin means of Re r / D_code(a) within 1 per cent of 1 for k <= 0.15 h/Mpc,
  - mean Im r / D_code(a) below 1 per cent in those bins,
  - cross-correlation coefficient |<deltam delta_L*>| / sqrt(<|deltam|^2><|delta_L|^2>) > 0.99 for k <= 0.3,
  - wrong hypotheses (x<->y, conjugate, half-cell shift) are far off.
usage: check_glam_deltaL.py <deltaL_prefix> <run_dir> [box_index=1]
"""
import glob, math, os, sys
import numpy as np

def read_meta(f):
    m = {}
    for line in open(f):
        if '=' in line:
            k, v = line.split('=', 1); m[k.strip()] = v.strip()
    return m

def read_cube(fb, n):
    arr = np.fromfile(fb, dtype='<c8')
    assert arr.size == (2 * n + 1) ** 2 * (n + 1), (fb, arr.size, n)
    return arr.reshape(2 * n + 1, 2 * n + 1, n + 1)

def main():
    dlp, run = sys.argv[1], sys.argv[2]
    ibox = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    ml = read_meta(dlp + '_deltaL.meta')
    nL = int(ml['n']); box = float(ml['Box'].split()[0])
    L = read_cube(dlp + '_deltaL.bin', nL)
    fails = []
    def check(ok, msg):
        print(('PASS ' if ok else 'FAIL ') + msg)
        if not ok: fails.append(msg)
    metas = sorted(glob.glob(os.path.join(run, 'deltam.*.%04d.meta' % ibox)))
    check(len(metas) > 0, '%d deltam cubes found in %s' % (len(metas), run))
    for fm in metas:
        m = read_meta(fm)
        n = int(m['n']); Nd = int(m['N_dep']); a = float(m['a']); D = float(m['D_code(a)'])
        boxg = float(m['Box'])
        check(abs(boxg - box) < 1e-3 * box, '%s: Box %g matches the generator' % (os.path.basename(fm), boxg))
        nn = min(n, nL)
        G = read_cube(fm[:-5] + '.bin', n)
        # common cube |ix|,|iy| <= nn, 0 <= iz <= nn
        Gc = G[n - nn:n + nn + 1, n - nn:n + nn + 1, :nn + 1]
        Lc = L[nL - nn:nL + nn + 1, nL - nn:nL + nn + 1, :nn + 1]
        kk = np.arange(-nn, nn + 1)
        KX, KY, KZ = np.meshgrid(kk, kk, np.arange(0, nn + 1), indexing='ij')
        kf = 2 * math.pi / box
        kmag = kf * np.sqrt(KX ** 2 + KY ** 2 + KZ ** 2)
        def W1(ki):
            s = np.pi * ki / Nd
            return np.where(ki == 0, 1.0, (np.sin(s) / np.where(s == 0, 1, s)) ** 2)
        Gd = Gc / (W1(KX) * W1(KY) * W1(KZ))
        sel = kmag > 0
        print('\n== %s: a = %.6g (z = %.3f), D_code(a) = %.6g, n = %d, N_dep = %d, %d modes compared'
              % (os.path.basename(fm), a, 1 / a - 1, D, n, Nd, sel.sum()))
        edges = [0, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.55]
        print('   k-bin [h/Mpc]      modes   <Re r>/D   <Im r>/D   rms(Re r/D - 1)   r_cross   P_m/P_L/D^2')
        for lo, hi in zip(edges[:-1], edges[1:]):
            s = sel & (kmag > lo) & (kmag <= hi)
            if not s.any(): continue
            r = Gd[s] / Lc[s] / D
            rc = abs(np.sum(Gd[s] * np.conj(Lc[s]))) / math.sqrt(np.sum(abs(Gd[s]) ** 2) * np.sum(abs(Lc[s]) ** 2))
            pr = np.sum(abs(Gd[s]) ** 2) / np.sum(abs(Lc[s]) ** 2) / D ** 2
            print('   %.2f < k <= %.2f  %6d   %8.5f   %8.5f   %8.4f          %7.5f   %8.5f'
                  % (lo, hi, s.sum(), r.real.mean(), r.imag.mean(), math.sqrt(((r.real - 1) ** 2).mean()), rc, pr))
            if hi <= 0.15:
                check(abs(r.real.mean() - 1) < 1e-2, '  %.2f<k<=%.2f: <Re r>/D = %.5f within 1 per cent of 1' % (lo, hi, r.real.mean()))
                check(abs(r.imag.mean()) < 1e-2, '  %.2f<k<=%.2f: |<Im r>|/D = %.5f below 1 per cent' % (lo, hi, abs(r.imag.mean())))
            if hi <= 0.30:
                check(rc > 0.99, '  %.2f<k<=%.2f: cross-correlation coefficient %.5f > 0.99' % (lo, hi, rc))
        # layout discriminators: transposed / conjugated cubes decorrelate completely
        s = sel & (kmag <= 0.3)
        def rcross(A, B): return abs(np.sum(A * np.conj(B))) / math.sqrt(np.sum(abs(A) ** 2) * np.sum(abs(B) ** 2))
        r_true = rcross(Gd[s], Lc[s])
        for k, v in {'x<->y transposed': np.transpose(Gd, (1, 0, 2)), 'conjugate (exp sign)': np.conj(Gd)}.items():
            rv = rcross(v[s], Lc[s])
            check(rv < 0.5 and r_true > 0.99, '  wrong hypothesis "%s" rejected (r_cross %.3f vs %.5f for the true layout, k <= 0.3)' % (k, rv, r_true))
        # coordinate offset per axis: a shift s_i multiplies r by exp(i k_i s_i); over the modes with
        # k_i > 0 the mean phase of r gives s_i = arg<r> / <k_i>.  (The lattice-on-nodes scatter has
        # zero mean phase.)  ic2pm's x = pos N/Box + 1 and GLAM's node i at (i-1) Box/N must give s = 0.
        for ax, K in (('x', KX), ('y', KY), ('z', KZ)):
            s = sel & (kmag > 0.05) & (kmag <= 0.15) & (K > 0)
            r = Gd[s] / Lc[s] / D
            shift = math.atan2(r.imag.mean(), r.real.mean()) / (kf * K[s].mean())
            check(abs(shift) < 0.05 * box / Nd, '  no coordinate offset along %s: inferred shift %+.4f Mpc/h (half a deposit cell = %.3f, %d modes)'
                  % (ax, shift, 0.5 * box / Nd, s.sum()))
    if fails:
        print('\n%d check(s) FAILED' % len(fails)); sys.exit(1)
    print('\nall GLAM deltam/delta_L checks passed')

if __name__ == '__main__':
    main()
