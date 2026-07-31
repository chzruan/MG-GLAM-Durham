#!/usr/bin/env python3
"""Compare the z~0 P(k) of the 2LPTIC-IC GLAM run (Run1 here) with the five
fiducial GLAM boxes (fid_LCDM_L512Np2048Ng4096/Run1-5, GLAM ZA ICs @ z=100).

Same box/Np/Ngrid/cosmology -> identical k bins. Differences expected:
realization scatter (seed 2026 vs seeds 1-5) at low k, and any IC-method
systematic (2LPT@z49 vs tuned-ZA@z100) at all k. Small a-offsets of the
final outputs are corrected with linear growth.
"""
import numpy as np, glob, re, sys
from scipy.integrate import quad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE = '/cosma8/data/dp203/dc-ruan1/mgglam_claude/MG-GLAM'
FID = f'{BASE}/fid_LCDM_L512Np2048Ng4096'
OURS = f'{BASE}/fid2LPTIC_L512Np2048Ng4096'
Om, OL = 0.3089, 0.6911

BLUE, ORANGE = '#2a78d6', '#eb6834'
INK, INK2, MUTED, GRID = '#0b0b0b', '#52514e', '#898781', '#e1e0d9'
SURF = '#fcfcfb'

def growth(a):
    E = lambda x: np.sqrt(Om/x**3 + OL)
    I, _ = quad(lambda x: 1.0/(x*E(x))**3, 1e-8, a, limit=200)
    return 2.5*Om*E(a)*I

def read_powerdm(fn):
    with open(fn) as f:
        f.readline()
        hdr = f.readline()
    a = float(re.search(r'Aexpn\s*=\s*([\d.]+)', hdr).group(1))
    d = np.loadtxt(fn, skiprows=3)
    return a, d[:, 1], d[:, 2], d[:, 3]     # a, k, Nmodes, P

def last_output(rundir, target_a=1.0):
    files = sorted(glob.glob(f'{rundir}/PowerDM.log.*.dat'))
    best, besta = None, None
    for fn in files:
        a, k, n, P = read_powerdm(fn)
        if besta is None or abs(a - target_a) < abs(besta - target_a):
            best, besta = (a, k, n, P), a
    return best

# ---- load ----
a_o, k_o, n_o, P_o = last_output(f'{OURS}/Run1')
a_c, k_c, n_c, P_c = last_output(f'{OURS}_da4/Run1')   # 244-step converged run
a_d, k_d, n_d, P_d = last_output(f'{OURS}_da6/Run1')   # 163-step default schedule
fid_data = [last_output(f'{FID}/Run{i}') for i in range(1, 6)]
a_ref = 1.0
D_t = growth(a_ref)

# growth-correct every P to a=1
Pc_o = P_o * (D_t / growth(a_o))**2
Pc_c = P_c * (D_t / growth(a_c))**2
Pc_d = P_d * (D_t / growth(a_d))**2
print(f"ours 123-step: a={a_o:.5f}; 163-step: a={a_d:.5f}; 244-step: a={a_c:.5f}")
Pfid = []
for i, (a, k, n, P) in enumerate(fid_data):
    assert np.allclose(k, k_o), f"k-bin mismatch Run{i+1}"
    c = (D_t / growth(a))**2
    print(f"fid Run{i+1}: a={a:.5f} (growth corr {c:.4f})")
    Pfid.append(P * c)
Pfid = np.array(Pfid)
Pmean, Pstd = Pfid.mean(axis=0), Pfid.std(axis=0, ddof=1)

ratio = Pc_o / Pmean
ratio_c = Pc_c / Pmean
ratio_d = Pc_d / Pmean
band = Pstd / Pmean / np.sqrt(5)          # error on the 5-box mean
scat = Pstd / Pmean                        # single-box scatter

# ---- plot ----
plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 9.5,
    'text.color': INK, 'axes.labelcolor': INK2,
    'axes.edgecolor': MUTED, 'xtick.color': MUTED, 'ytick.color': MUTED,
    'xtick.labelcolor': INK2, 'ytick.labelcolor': INK2,
    'axes.facecolor': SURF, 'figure.facecolor': SURF,
    'grid.color': GRID, 'grid.linewidth': 0.6,
    'axes.grid': True, 'axes.axisbelow': True, 'legend.frameon': False,
})
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.0, 7.2), sharex=True,
                               gridspec_kw={'height_ratios': [1.2, 1.0], 'hspace': 0.08})

for i in range(5):
    ax1.plot(k_o, Pfid[i], '-', color=MUTED, lw=0.7, alpha=0.6,
             label='fid GLAM boxes (5 seeds, ZA IC z=100)' if i == 0 else None)
ax1.plot(k_o, Pc_o, '-', color=BLUE, lw=1.6, label='2LPTIC IC run (seed 2026, 2LPT z=49)')
ax1.set_xscale('log'); ax1.set_yscale('log')
ax1.set_ylabel(r'$P(k,z{=}0)\ \ [(\mathrm{Mpc}/h)^3]$')
ax1.set_title('GLAM $z=0$ power spectrum: 2LPTIC IC vs the five fiducial boxes\n'
              r'$L=512\,\mathrm{Mpc}/h$, $2048^3$, Ng=4096, fid cosmology',
              fontsize=10.5, color=INK, pad=10)
ax1.legend(loc='lower left', fontsize=8.8)

ax2.axhline(1.0, color=MUTED, lw=0.8)
# measured fid step-0 amplitude: P = 0.989 x intended table (incl. Pk tune x1.010).
# If fid boxes only had that IC offset and grew perfectly, ours/fid = 1/0.989.
ax2.axhline(1/0.989, color=INK2, lw=1.0, ls='--',
            label='fid IC amplitude offset alone (step-0 measured, incl. Pk tune)')
ax2.fill_between(k_o, 1-scat, 1+scat, color=MUTED, alpha=0.25,
                 label='single-box scatter (5 fid boxes)')
ax2.fill_between(k_o, 1-band, 1+band, color=MUTED, alpha=0.45,
                 label='error on 5-box mean')
ax2.plot(k_o, ratio, '-', color=BLUE, lw=0.9, alpha=0.4,
         label='2LPTIC run, 123 steps / fid mean')
ax2.plot(k_d, ratio_d, '-', color=ORANGE, lw=1.3,
         label='2LPTIC run, 163 steps (default schedule) / fid mean')
ax2.plot(k_c, ratio_c, '-', color=BLUE, lw=1.7,
         label='2LPTIC run, 244 steps (converged) / fid mean')
ax2.set_xlabel(r'$k\ \ [h/\mathrm{Mpc}]$')
ax2.set_ylabel(r'$P_{\rm 2LPTIC}\,/\,\langle P_{\rm fid}\rangle$')
ax2.set_ylim(0.9, 1.3)
ax2.legend(loc='upper left', fontsize=8.8)
ax2.text(0.985, 0.04,
         'ratio above the dashed line = pure ZA@z=100 transient under-growth\n'
         'in the fid boxes (2LPT@z=49 avoids it); dashed line = their measured\n'
         'step-0 amplitude deficit (0.989 in P, ALREADY including Pk tune = 1.005,\n'
         'i.e. x1.010 in P — without the tune the gap would be ~1% larger)',
         transform=ax2.transAxes, ha='right', va='bottom', fontsize=7.6, color=MUTED)

fig.savefig(f'{OURS}/pk_fid_compare.png', dpi=200, bbox_inches='tight')
fig.savefig(f'{OURS}/pk_fid_compare.pdf', bbox_inches='tight')

sel = (k_o > 0.1) & (k_o < 1.0)
print(f"\n123-step ratio (0.1<k<1): median={np.median(ratio[sel]):.4f}")
print(f"244-step ratio (0.1<k<1): median={np.median(ratio_c[sel]):.4f}")
sel2 = (k_o > 1.0) & (k_o < 5.0)
print(f"244-step ratio (1<k<5):   median={np.median(ratio_c[sel2]):.4f}")
print("wrote pk_fid_compare.png/.pdf")
