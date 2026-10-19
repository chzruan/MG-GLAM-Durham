"""Field-level checks of the white noise drawn by ranlux_lags.exe (mode 1, with a cube written).

For each cube (NROW^3, GLAM's draw order i, j, k) the white field f(x, y, z) = sum_ijk TS(i, j, k) B(x, i) B(y, j) B(z, k)
is synthesised with the basis of VECTOR (B(:, i) = rfft1b of the i-th unit vector, written by the test as cube.basis).
Reported, per cube:
  - the variance profile along each axis, v(x) = <f^2>_(y, z): the fractional amplitude 2 |V_m| / v0 of its harmonics m = 1..24,
    and the noise level (median over m = 40..120);
  - the power of delta^2 = f^2 - <f^2> at the axis modes (+-m, 0, 0), (0, +-m, 0), (0, 0, +-m) over the mean power of the
    non-axis modes with |n| within 0.5 of m, for m = 1..20 (homogeneous noise: about 1, with the scatter of a two-mode mean);
  - the lag correlation of TS along i, j, k from the test's own text output.
Also the 1-D variance profile of the NROW = 2048 runs (outfile.vprof): harmonics m = 1..30.
Usage: micromamba run -n cosemu python3 analyse_cubes.py <workdir> <outdir>
Files read: <workdir>/cube256_*.bin(.basis, .txt, .txt.vprof) and <workdir>/fix_s*_lux*.txt(.vprof); results to <outdir>.
"""
import glob, json, os, sys
import numpy as np

work, outdir = sys.argv[1], sys.argv[2]
os.makedirs(outdir, exist_ok=True)
out = {}


def harmonics(v, mmax):
    V = np.fft.rfft(v - v.mean()) / v.size
    amp = 2 * np.abs(V) / v.mean()
    return amp[1:mmax + 1], float(np.median(amp[40:121])) if amp.size > 121 else float(np.median(amp[mmax + 1:]))


def read_lags(txt):
    lags, hdr = {}, ''
    for ln in open(txt):
        if ln.startswith('lag '):
            p = ln.split()
            lags[int(p[1])] = [float(p[2]), float(p[3]), float(p[4])]
        if ln.startswith('lux_requested'):
            hdr = ln.strip()
    return hdr, lags


# ---- the 256^3 cubes
for cube in sorted(glob.glob(os.path.join(work, 'cube256_*.bin'))):
    tag = os.path.basename(cube)[:-4]
    hdr, lags = read_lags(cube[:-4] + '.txt')
    n = int(round((os.path.getsize(cube) / 4) ** (1 / 3)))
    ts = np.fromfile(cube, '>f4').reshape(n, n, n, order='F').astype(np.float64)       # TS(i, j, k)
    B = np.fromfile(cube + '.basis', '>f8').reshape(n, n, order='F')                     # B(x, i)
    f = np.einsum('xi,ijk->xjk', B, ts)
    f = np.einsum('yj,xjk->xyk', B, f)
    f = np.einsum('zk,xyk->xyz', B, f)
    del ts
    res = {'header': hdr, 'n': n, 'field_variance': float(np.mean(f**2)), 'axes': {}}
    for ax, name in enumerate('xyz'):
        other = tuple(a for a in range(3) if a != ax)
        v = np.mean(f**2, axis=other)
        amp, noise = harmonics(v, 24)
        res['axes'][name] = {'rel_modulation_amp_m1_24': amp.round(5).tolist(), 'noise_m40_120': round(noise, 5),
                             'top_harmonics': [int(i) for i in np.argsort(-amp)[:6] + 1]}
    d2 = f**2 - np.mean(f**2)
    del f
    P = np.abs(np.fft.fftn(d2))**2
    del d2
    nn = np.fft.fftfreq(n, 1 / n)
    n2 = nn[:, None, None]**2 + nn[None, :, None]**2 + nn[None, None, :]**2
    axis = ((nn[:, None, None] != 0).astype(int) + (nn[None, :, None] != 0) + (nn[None, None, :] != 0)) == 1
    ratios = {name: [] for name in 'xyz'}
    for m in range(1, 21):
        ref = P[(~axis) & (np.abs(np.sqrt(n2) - m) <= 0.5)].mean()
        ratios['x'].append(float((P[m, 0, 0] + P[-m, 0, 0]) / 2 / ref))
        ratios['y'].append(float((P[0, m, 0] + P[0, -m, 0]) / 2 / ref))
        ratios['z'].append(float((P[0, 0, m] + P[0, 0, -m]) / 2 / ref))
    res['d2_axis_mode_power_over_neighbours_m1_20'] = {k: np.round(v, 2).tolist() for k, v in ratios.items()}
    res['lag_corr_ijk_1_40'] = lags
    out[tag] = res
    print(f'{tag}: {hdr}')
    for name in 'xyz':
        a = res['axes'][name]
        print(f'  variance profile along {name}: harmonics 1-24 (x100): ' + ' '.join(f'{100 * x:.2f}' for x in a['rel_modulation_amp_m1_24'])
              + f' | noise {100 * a["noise_m40_120"]:.2f} | top {a["top_harmonics"]}')
        print(f'  delta^2 power at (+-m,0,0)-type modes along {name} over neighbours, m = 1..20: '
              + ' '.join(f'{x:.1f}' for x in ratios[name]))
    del P

# ---- the 1-D profiles of the NROW = 2048 runs
prof = {}
for txt in sorted(glob.glob(os.path.join(work, 'fix_s*_lux*.txt'))):
    tag = os.path.basename(txt)[:-4]
    hdr, lags = read_lags(txt)
    v = np.loadtxt(txt + '.vprof')[:, 1]
    amp, noise = harmonics(v, 30)
    prof[tag] = v / v.mean() - 1
    out[tag] = {'header': hdr, 'rows': None, 'vprof_rel_modulation_amp_m1_30': amp.round(5).tolist(), 'vprof_noise_m40_120': round(noise, 5),
                'top_harmonics': [int(i) for i in np.argsort(-amp)[:6] + 1], 'lag_corr_ijk_1_40': lags,
                'max_abs_lag_corr_i': max(abs(c[0]) for c in lags.values()), 'max_abs_lag_corr_j': max(abs(c[1]) for c in lags.values()),
                'max_abs_lag_corr_k': max(abs(c[2]) for c in lags.values())}
    print(f'{tag}: {hdr}: 1-D variance profile harmonics 1-30 (x100): ' + ' '.join(f'{100 * x:.2f}' for x in amp)
          + f' | noise {100 * noise:.2f} | top {out[tag]["top_harmonics"]}')
for a, b in (('fix_s1_lux0', 'fix_s2_lux0'), ('fix_s1_lux3', 'fix_s2_lux3'), ('fix_s1_lux0', 'fix_s1_lux3')):
    if a in prof and b in prof:
        c = float(np.corrcoef(prof[a], prof[b])[0, 1])
        out[f'profile_correlation_{a}_{b}'] = c
        print(f'correlation of the 1-D variance profiles {a} vs {b}: {c:+.3f}')
json.dump(out, open(os.path.join(outdir, 'analyse_cubes.json'), 'w'), indent=1)
