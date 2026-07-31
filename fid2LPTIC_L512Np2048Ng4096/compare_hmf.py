#!/usr/bin/env python3
"""Halo mass function comparison at z=0: the three 2LPTIC-IC GLAM runs
(123/163/244 steps, same IC) vs the five fiducial ZA-IC boxes.

Masses: BDM Mtot [Msun/h] from CatshortV (overdensity ~333.5, mp=1.339e9).
Each run's dn/dlog10M is interpolated linearly in (ln n, ln a) to a=1
using its last two outputs (final a's differ by up to 0.7% between
schedules, which matters at the exponential tail).
"""
import numpy as np, glob, re
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE = '/cosma8/data/dp203/dc-ruan1/mgglam_claude/MG-GLAM'
BOX = 512.0
VOL = BOX**3
MP = 1.3389e9

BLUE, ORANGE, INK, INK2, MUTED, GRID = '#2a78d6', '#eb6834', '#0b0b0b', '#52514e', '#898781', '#e1e0d9'
SURF = '#fcfcfb'

bins = np.arange(11.3, 15.61, 0.1)
mid = 0.5*(bins[1:]+bins[:-1])

def read_cat_hmf(fn):
    a = float(re.search(r'A\s*=\s*([\d.]+)', open(fn).readlines()[1]).group(1))
    m = pd.read_csv(fn, sep=r'\s+', skiprows=8, usecols=[7], header=None,
                    dtype=np.float64).values[:, 0]
    cnt, _ = np.histogram(np.log10(m), bins=bins)
    return a, cnt

def hmf_at_a1(rundir):
    files = sorted(glob.glob(f'{rundir}/CATALOGS/CatshortV.*.DAT'))[-2:]
    (a1, c1), (a2, c2) = read_cat_hmf(files[0]), read_cat_hmf(files[1])
    n1, n2 = c1/VOL/0.1, c2/VOL/0.1
    ok = (c1 > 0) & (c2 > 0)
    lnn = np.full(len(mid), -np.inf)
    # interpolate ln n in ln a to a=1; where either count is 0, keep final
    lnn[ok] = np.log(n2[ok]) + (np.log(n1[ok]) - np.log(n2[ok])) * \
              (np.log(1.0) - np.log(a2)) / (np.log(a1) - np.log(a2))
    lnn[~ok & (c2 > 0)] = np.log(n2[~ok & (c2 > 0)])
    print(f"{rundir.split('/')[-2]:42s} a=({a1:.5f},{a2:.5f}) Nhalo(final)={c2.sum():,}")
    return np.exp(lnn), c2

runs = {
    '123': f'{BASE}/fid2LPTIC_L512Np2048Ng4096/Run1',
    '163': f'{BASE}/fid2LPTIC_L512Np2048Ng4096_da6/Run1',
    '244': f'{BASE}/fid2LPTIC_L512Np2048Ng4096_da4/Run1',
}
ours, cnts = {}, {}
for k_, d in runs.items():
    ours[k_], cnts[k_] = hmf_at_a1(d)

fid_n, fid_c = [], []
for r in range(1, 6):
    n, c = hmf_at_a1(f'{BASE}/fid_LCDM_L512Np2048Ng4096/Run{r}')
    fid_n.append(n); fid_c.append(c)
fid_n = np.array(fid_n)
nmean = fid_n.mean(axis=0)
nscat = fid_n.std(axis=0, ddof=1)
fid_ctot = np.sum(fid_c, axis=0)

plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 9.5,
    'text.color': INK, 'axes.labelcolor': INK2,
    'axes.edgecolor': MUTED, 'xtick.color': MUTED, 'ytick.color': MUTED,
    'xtick.labelcolor': INK2, 'ytick.labelcolor': INK2,
    'axes.facecolor': SURF, 'figure.facecolor': SURF,
    'grid.color': GRID, 'grid.linewidth': 0.6,
    'axes.grid': True, 'axes.axisbelow': True, 'legend.frameon': False,
})
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.0, 7.4), sharex=True,
                               gridspec_kw={'height_ratios': [1.15, 1.0], 'hspace': 0.08})

for i in range(5):
    ax1.plot(mid, fid_n[i], '-', color=MUTED, lw=0.7, alpha=0.6,
             label='fid GLAM boxes (5 seeds, ZA IC z=100)' if i == 0 else None)
ax1.plot(mid, ours['123'], '-', color=BLUE, lw=0.9, alpha=0.4, label='2LPTIC, 123 steps')
ax1.plot(mid, ours['163'], '-', color=ORANGE, lw=1.3, label='2LPTIC, 163 steps (default)')
ax1.plot(mid, ours['244'], '-', color=BLUE, lw=1.7, label='2LPTIC, 244 steps (converged)')
ax1.set_yscale('log')
ax1.set_ylabel(r'$dn/d\log_{10}M\ \ [(h/\mathrm{Mpc})^3]$')
ax1.set_title('GLAM $z=0$ halo mass function: 2LPTIC IC vs the five fiducial boxes\n'
              r'BDM $M_{\rm tot}$, $L=512\,\mathrm{Mpc}/h$, $2048^3$, Ng=4096; all interpolated to $a=1$',
              fontsize=10.5, color=INK, pad=10)
ax1.legend(loc='lower left', fontsize=8.6)

sel = fid_ctot >= 100
ax2.axhline(1.0, color=MUTED, lw=0.8)
ax2.fill_between(mid[sel], (1-nscat/nmean)[sel], (1+nscat/nmean)[sel],
                 color=MUTED, alpha=0.25, label='single-box scatter (5 fid boxes)')
for key, col, lw, al in (('123', BLUE, 0.9, 0.4), ('163', ORANGE, 1.3, 1.0), ('244', BLUE, 1.7, 1.0)):
    r = ours[key]/nmean
    err = r/np.sqrt(np.maximum(cnts[key], 1))
    ax2.plot(mid[sel], r[sel], '-', color=col, lw=lw, alpha=al,
             label=f'2LPTIC {key} steps / fid mean')
    if key == '244':
        ax2.fill_between(mid[sel], (r-err)[sel], (r+err)[sel], color=col, alpha=0.15)
ax2.set_xlabel(r'$\log_{10}\,M_{\rm tot}\ \ [M_\odot/h]$')
ax2.set_ylabel(r'$n_{\rm 2LPTIC}\,/\,\langle n_{\rm fid}\rangle$')
ax2.set_ylim(0.85, 1.42)
ax2.legend(loc='upper left', fontsize=8.6)
ax2.text(0.985, 0.97,
         'low-mass excess is stepping-independent (all 3 curves overlap):\n'
         'ZA@z=100 starts under-produce marginally-resolved haloes\n'
         r'($\lesssim$300 particles); well-resolved $10^{13}\,M_\odot/h$ agree to ~1%',
         transform=ax2.transAxes, ha='right', va='top', fontsize=7.6, color=MUTED)

fig.savefig(f'{BASE}/fid2LPTIC_L512Np2048Ng4096/hmf_fid_compare.png', dpi=200, bbox_inches='tight')
fig.savefig(f'{BASE}/fid2LPTIC_L512Np2048Ng4096/hmf_fid_compare.pdf', bbox_inches='tight')

print('\nratio n_2LPTIC / <n_fid> at z=0 (a=1):')
print(f"{'log10M':>8} {'123':>8} {'163':>8} {'244':>8} {'fid scat':>9}")
for lm in (11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 14.9):
    i = np.argmin(np.abs(mid-lm))
    print(f"{mid[i]:8.2f} {ours['123'][i]/nmean[i]:8.3f} {ours['163'][i]/nmean[i]:8.3f} "
          f"{ours['244'][i]/nmean[i]:8.3f} {nscat[i]/nmean[i]:9.3f}")
print('\nwrote hmf_fid_compare.png/.pdf')
