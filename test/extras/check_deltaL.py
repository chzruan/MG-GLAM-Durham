#!/usr/bin/env python3
"""Layout / convention check of the FML delta_L dump (I7, route 1).

Two independent reconstructions of the Gaussian field from the Gadget-1 IC
particles written by the FML generator are compared mode by mode with the
delta_L cube (both in MG-GLAM's convention
delta_k = N^-3 sum_x delta(x) exp(-i k.x), k = (2 pi/Box)(ix,iy,iz), grid
value j at x = j Box/N):

 A. Displacement field.  The particles sit on the lattice q = i Box/Np (no
    half-cell offset); Psi(q) = pos - q is sampled exactly on that lattice,
    so its FFT is exact and linear theory gives
        delta(k) = -i k . Psi(k),   delta_L = delta / D(a_ini).
    For 1LPT (Zel'dovich) ICs this is exact to float32 rounding; for 2LPT the
    second-order displacement adds O(delta(z_ini)) ~ 1e-3 per mode.
 B. CIC density of the particles, deconvolved, divided by D(a_ini).  The
    particles sit on the nodes of the natural grid, where the CIC weight
    1 - |Psi|/dx is not linear in Psi, so the deposit grid is shifted by half
    a cell (weights 1/2 -+ Psi/dx, linear for |Psi| < dx/2) and the phase is
    corrected.  Residual second-order terms (Psi_x Psi_y cross terms, the
    real second-order density) leave per-mode scatter of a few per cent at
    k ~ 0.5 h/Mpc, so only the mean ratio is asserted.

Wrong hypotheses (x<->y transposed, opposite exponent sign, half-cell shift)
are evaluated too, to show the test discriminates.

usage: check_deltaL.py <IC_prefix> <deltaL_prefix>
"""
import glob, struct, sys, math
import numpy as np

def read_gadget1(prefix):
    files = sorted(glob.glob(prefix + '.[0-9]*'), key=lambda f: int(f.rsplit('.', 1)[1]))
    assert files, 'no IC files ' + prefix + '.*'
    pos, ids, hdr = [], [], None
    for f in files:
        with open(f, 'rb') as fh:
            def block():
                n, = struct.unpack('<i', fh.read(4)); d = fh.read(n)
                m, = struct.unpack('<i', fh.read(4)); assert n == m, (f, n, m); return d
            h = block(); assert len(h) == 256
            npart = struct.unpack('<6I', h[:24])
            time, z = struct.unpack('<2d', h[72:88]); nall = struct.unpack('<6I', h[96:120])
            nfiles, box = struct.unpack('<id', h[124:136]); om, ol, hub = struct.unpack('<3d', h[136:160])
            n = sum(npart)
            p = np.frombuffer(block(), dtype='<f4').reshape(n, 3); block()  # vel
            idb = block(); ids.append(np.frombuffer(idb, dtype='<i8' if len(idb) == 8 * n else '<u4'))
            pos.append(p); hdr = dict(a=time, z=z, box=box, Om=om, OL=ol, h=hub, ntot=sum(nall), nfiles=nfiles)
    return np.concatenate(pos).astype(np.float64), np.concatenate(ids), hdr

def cic_halfshift(pos, box, N):
    x = pos / box * N - 0.5                      # deposit-grid nodes at (j + 1/2) dx
    i0 = np.floor(x).astype(np.int64); d = x - i0
    rho = np.zeros((N, N, N))
    for dx in (0, 1):
        wx = d[:, 0] if dx else 1 - d[:, 0]
        for dy in (0, 1):
            wy = d[:, 1] if dy else 1 - d[:, 1]
            for dz in (0, 1):
                wz = d[:, 2] if dz else 1 - d[:, 2]
                np.add.at(rho, ((i0[:, 0] + dx) % N, (i0[:, 1] + dy) % N, (i0[:, 2] + dz) % N), wx * wy * wz)
    return rho / rho.mean() - 1

def growth_simpson(a, Om, OL, nint=4000):
    x0 = math.log(1e-8); h = (math.log(a) - x0) / nint
    xs = x0 + h * np.arange(nint + 1); f = np.exp(2.5 * xs) / (Om + OL * np.exp(3 * xs)) ** 1.5
    w = np.full(nint + 1, 2.0); w[1::2] = 4.0; w[0] = w[-1] = 1.0
    return math.sqrt(Om / a ** 3 + OL) * (f * w).sum() * h / 3
def D_simpson(a, Om, OL): return growth_simpson(a, Om, OL) / growth_simpson(1.0, Om, OL)
def D_quad(a, Om, OL):
    from scipy.integrate import quad
    E = lambda x: math.sqrt(Om / x ** 3 + OL)
    g = lambda x: E(x) * quad(lambda y: 1.0 / (y * E(y)) ** 3, 0.0, x, epsabs=0, epsrel=1e-10, limit=200)[0]
    return g(a) / g(1.0)

def read_meta(f):
    m = {}
    for line in open(f):
        if '=' in line:
            k, v = line.split('=', 1); m[k.strip()] = v.strip()
    return m

def main():
    icp, dlp = sys.argv[1], sys.argv[2]
    meta = read_meta(dlp + '_deltaL.meta')
    n = int(meta['n']); N = int(meta['N_mesh']); box = float(meta['Box'].split()[0])
    lpt = int(meta.get('lpt_order', 2))
    cube = np.fromfile(dlp + '_deltaL.bin', dtype='<c8')
    assert cube.size == (2 * n + 1) ** 2 * (n + 1), (cube.size, n)
    cube = cube.reshape(2 * n + 1, 2 * n + 1, n + 1)
    pos, ids, hdr = read_gadget1(icp)
    fails = []
    def check(ok, msg):
        print(('PASS ' if ok else 'FAIL ') + msg)
        if not ok: fails.append(msg)
    print('IC %s: %d particles, a = %.6g, box = %g; cube n = %d (k_c = %s), lpt_order = %d'
          % (icp, pos.shape[0], hdr['a'], hdr['box'], n, meta['k_c'], lpt))
    check(hdr['ntot'] == pos.shape[0] == int(meta['N_part_1D']) ** 3, 'particle count %d' % pos.shape[0])
    check(abs(hdr['box'] - box) < 1e-6 * box and abs(1 / hdr['a'] - 1 - float(meta['z_ini'])) < 1e-6,
          'IC header box/z match meta (box %g, z %g)' % (hdr['box'], hdr['z']))
    Dm = float(meta['D_gen(a_ini)']); Ds = D_simpson(hdr['a'], hdr['Om'], hdr['OL'])
    check(abs(Dm / Ds - 1) < 1e-6, 'D(a_ini): meta %.8g, Simpson %.8g' % (Dm, Ds))
    try:
        Dq = D_quad(hdr['a'], hdr['Om'], hdr['OL'])
        check(abs(Dm / Dq - 1) < 1e-5, 'D(a_ini) vs adaptive quadrature of the Heath integral %.8g' % Dq)
    except ImportError:
        print('skip scipy cross-check of D')
    check(cube[n, n, 0] == 0, 'k = 0 mode zero')
    herm = np.abs(cube[::-1, ::-1, 0] - np.conj(cube[:, :, 0])).max() / np.abs(cube[:, :, 0]).max()
    check(herm < 1e-6, 'Hermitian iz = 0 plane, max rel diff %.2e' % herm)

    # ---- lattice from IDs: ID = 1 + ix + Np (iy + Np iz), q = i Box/Np
    Np = int(meta['N_part_1D']); i = (ids - 1); ixl, iyl, izl = i % Np, (i // Np) % Np, i // (Np * Np)
    check(len(np.unique(ids)) == ids.size and ids.min() == 1 and ids.max() == Np ** 3, 'IDs unique 1..Np^3')
    q = np.stack([ixl, iyl, izl], 1) * (box / Np)
    psi = pos - q; psi -= box * np.round(psi / box)
    check(np.abs(psi).max() < 0.5 * box / Np, 'ID lattice: max |pos - q| = %.3f Mpc/h (< half a cell %.3f)'
          % (np.abs(psi).max(), 0.5 * box / Np))
    kf = 2 * math.pi / box
    kk = np.arange(-n, n + 1)
    KX, KY, KZ = np.meshgrid(kk, kk, np.arange(0, n + 1), indexing='ij')
    kmag = kf * np.sqrt(KX ** 2 + KY ** 2 + KZ ** 2)
    sel = kmag > 0
    def report(name, Fc, tol_mean, tol_resid, tol_mode=None, alt_factor=10.0):
        ratio = Fc[sel] / cube[sel]
        rr, ri = ratio.real, ratio.imag
        resid = np.sqrt((np.abs(Fc[sel] - cube[sel]) ** 2).mean() / (np.abs(cube[sel]) ** 2).mean())
        print('%s: %d modes, ratio to delta_L: <Re> = %.6f, rms(Re-1) = %.2e, rms Im = %.2e, max |ratio-1| = %.2e, rms residual / rms delta_L = %.2e'
              % (name, sel.sum(), rr.mean(), np.sqrt(((rr - 1) ** 2).mean()), np.sqrt((ri ** 2).mean()), np.abs(ratio - 1).max(), resid))
        kb = np.linspace(0, kmag[sel].max() * 1.0001, 6)
        for lo, hi in zip(kb[:-1], kb[1:]):
            s2 = sel & (kmag > lo) & (kmag <= hi)
            r2 = Fc[s2] / cube[s2]
            print('    %.3f < k <= %.3f : %5d modes  <Re> = %.6f  rms(Re-1) = %.2e  rms Im = %.2e'
                  % (lo, hi, s2.sum(), r2.real.mean(), np.sqrt(((r2.real - 1) ** 2).mean()), np.sqrt((r2.imag ** 2).mean())))
        check(abs(rr.mean() - 1) < tol_mean, '%s: mean ratio within %.0e of 1' % (name, tol_mean))
        check(resid < tol_resid, '%s: rms residual below %.0e of rms delta_L' % (name, tol_resid))
        if tol_mode is not None:
            check(np.abs(ratio - 1).max() < tol_mode, '%s: every mode within %.0e of 1' % (name, tol_mode))
        alts = {'x<->y transposed': np.transpose(Fc, (1, 0, 2)),
                'conjugate (exp sign)': np.conj(Fc),
                'half-cell shift': Fc * np.exp(1j * np.pi * (KX + KY + KZ) / N)}
        base = np.abs(ratio - 1).mean()
        for k, v in alts.items():
            dev = np.abs(v[sel] / cube[sel] - 1).mean()
            check(dev > alt_factor * base, '%s: wrong hypothesis "%s" rejected (mean |ratio-1| %.2e vs %.2e, factor %.0f required)' % (name, k, dev, base, alt_factor))

    # ---- A: exact displacement field on the lattice (Np = Nmesh assumed for the index mapping)
    if Np == N:
        PSI = np.zeros((N, N, N, 3))
        PSI[ixl, iyl, izl, :] = psi
        Pk = np.stack([np.fft.fftn(PSI[..., c]) / N ** 3 for c in range(3)], -1)
        Pk = Pk[KX % N, KY % N, KZ % N, :]
        delta_lin = -1j * kf * (KX * Pk[..., 0] + KY * Pk[..., 1] + KZ * Pk[..., 2])
        if lpt == 1:   # exact up to float32 position rounding
            report('A displacement (-i k.Psi, 1LPT)', delta_lin / Ds, 1e-5, 1e-4, tol_mode=2e-3)
        else:          # second-order displacement ~ delta(z_ini) ~ 1e-3 .. 1e-2 of the modes
            report('A displacement (-i k.Psi, %dLPT)' % lpt, delta_lin / Ds, 2e-3, 2e-2)
    else:
        print('skip A: Np != Nmesh')
    # ---- B: CIC density, half-cell-shifted deposit, deconvolved
    F = np.fft.fftn(cic_halfshift(pos, box, N)) / N ** 3
    f = np.fft.fftfreq(N) * N
    F = F * np.exp(-1j * np.pi * (f[:, None, None] + f[None, :, None] + f[None, None, :]) / N)
    def W1(ki):
        s = np.pi * ki / N
        return np.where(ki == 0, 1.0, (np.sin(s) / np.where(s == 0, 1, s)) ** 2)
    Fc = F[KX % N, KY % N, KZ % N] / (W1(KX) * W1(KY) * W1(KZ)) / Ds
    report('B CIC density', Fc, 2e-3, 1e-1, alt_factor=3.0)   # second-order CIC terms leave per-mode scatter
    print('rms |delta_L| cube = %.6g (meta %s)' % (np.sqrt((np.abs(cube) ** 2).mean()), meta.get('rms_deltaL_cube')))
    if fails:
        print('\n%d check(s) FAILED' % len(fails)); sys.exit(1)
    print('\nall delta_L layout checks passed')

if __name__ == '__main__':
    main()
