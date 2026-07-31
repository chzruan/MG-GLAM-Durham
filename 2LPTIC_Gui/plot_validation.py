#!/usr/bin/env python3
"""P(k) validation figure: our Intel-built FML LPT ICs vs the existing
DEGRACE low-res (original 2LPTic) and HEFT high-res (Gui's FML) ICs,
plus the input linear spectra."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE = '/cosma8/data/dp203/dc-ruan1/mgglam_claude/MG-GLAM/2LPTIC_Gui'

# palette (dataviz reference, light mode)
BLUE, ORANGE = '#2a78d6', '#eb6834'
INK, INK2, MUTED, GRID = '#0b0b0b', '#52514e', '#898781', '#e1e0d9'
SURF = '#fcfcfb'

def load_npz(path):
    d = np.load(path)
    k, Po, Pr, Px, modes = d['k'], d['Pours'], d['Pref'], d['Pcross'], d['modes']
    sel = (k > 0) & (modes > 0) & np.isfinite(Px) & np.isfinite(Po) & np.isfinite(Pr)
    return k[sel], Po[sel], Pr[sel], Px[sel]

k_lo, Po_lo, Pr_lo, Px_lo = load_npz(f'{BASE}/validate_lowres/pk_validate_lowres.npz')
k_hi, Po_hi, Pr_hi, Px_hi = load_npz(f'{BASE}/validate_highres/pk_validate_highres.npz')

lin_lo = np.loadtxt(f'{BASE}/validate_lowres/pofk_lowres_exact_z49.txt')   # low-res norm
lin_hi = np.loadtxt(f'{BASE}/validate_highres/pofk_bli_z49.txt')           # HEFT norm (x1.0181)

def lin_interp(lin, k):
    return np.exp(np.interp(np.log(k), np.log(lin[:, 0]), np.log(lin[:, 1])))

Plin_lo = lin_interp(lin_lo, k_lo)
Plin_hi = lin_interp(lin_hi, k_hi)

plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 9.5,
    'text.color': INK, 'axes.labelcolor': INK2,
    'axes.edgecolor': MUTED, 'xtick.color': MUTED, 'ytick.color': MUTED,
    'xtick.labelcolor': INK2, 'ytick.labelcolor': INK2,
    'axes.facecolor': SURF, 'figure.facecolor': SURF,
    'grid.color': GRID, 'grid.linewidth': 0.6,
    'axes.grid': True, 'axes.axisbelow': True,
    'legend.frameon': False,
})

fig, (ax1, ax2, ax3) = plt.subplots(
    3, 1, figsize=(7.0, 9.2), sharex=True,
    gridspec_kw={'height_ratios': [1.35, 1.0, 1.0], 'hspace': 0.08})

# ---------------- (a) absolute P(k) ----------------
kk = lin_lo[(lin_lo[:, 0] > 5e-3) & (lin_lo[:, 0] < 8), 0]
ax1.plot(kk, lin_interp(lin_lo, kk), '-', color=INK2, lw=1.2, zorder=3,
         label='input linear $P(k,z=49)$')
s = slice(3, None, 12)   # subsample reference markers
ax1.plot(k_lo[s], Pr_lo[s], 'o', ms=4.5, mfc='none', mec=BLUE, mew=1.2, zorder=4,
         label=r'2LPTic $1024^3$ (DEGRACE, ics.*)')
ax1.plot(k_lo, Po_lo, '-', color=BLUE, lw=1.6, zorder=5,
         label=r'ours $1024^3$ (FML, Intel)')
ax1.plot(k_hi[s], Pr_hi[s], 's', ms=4.0, mfc='none', mec=ORANGE, mew=1.2, zorder=4,
         label=r'FML $2048^3$ (Gui, HEFT)')
ax1.plot(k_hi, Po_hi, '-', color=ORANGE, lw=1.6, zorder=5,
         label=r'ours $2048^3$ (FML, Intel)')
ax1.set_xscale('log'); ax1.set_yscale('log')
ax1.set_ylabel(r'$P(k)\ \ [(\mathrm{Mpc}/h)^3]$')
ax1.set_title('IC power spectra: our Intel-built 2LPTIC (FML) vs existing DEGRACE / HEFT ICs\n'
              r'$L=1024\,\mathrm{Mpc}/h$, $z_{\rm ini}=49$, seed 2026',
              fontsize=10.5, color=INK, pad=10)
ax1.legend(loc='lower left', fontsize=8.5, handlelength=2.2)
ax1.set_ylim(4e-3, 40)

# ---------------- (b) measured / input linear ----------------
ax2.axhline(1.0, color=MUTED, lw=0.8)
ax2.plot(k_lo, Pr_lo/Plin_lo, 'o', ms=3.6, mfc='none', mec=BLUE, mew=1.0,
         markevery=4, label=r'2LPTic $1024^3$ / input')
ax2.plot(k_lo, Po_lo/Plin_lo, '-', color=BLUE, lw=1.5, label=r'ours $1024^3$ / input')
ax2.plot(k_hi, Pr_hi/Plin_hi, 's', ms=3.2, mfc='none', mec=ORANGE, mew=1.0,
         markevery=4, label=r'FML $2048^3$ / input')
ax2.plot(k_hi, Po_hi/Plin_hi, '-', color=ORANGE, lw=1.5, label=r'ours $2048^3$ / input')
ax2.set_ylabel(r'$P_{\rm measured}\,/\,P_{\rm lin,\,input}$')
ax2.set_ylim(0.88, 1.14)
leg2 = ax2.legend(loc='upper left', fontsize=8.5, ncol=2, handlelength=2.0,
                  frameon=True, facecolor=SURF, edgecolor=GRID, framealpha=0.95)
leg2.set_zorder(10)
ax2.text(0.985, 0.05,
         'each set divided by its own input;\n'
         r'HEFT norm = 1.0181$\times$ low-res ($\sigma_8$: 0.805 vs 0.8001)',
         transform=ax2.transAxes, ha='right', va='bottom', fontsize=7.8, color=MUTED)

# ---------------- (c) ours vs reference, machine level ----------------
eps = 1e-12
ax3.plot(k_lo, np.abs(Po_lo/Pr_lo - 1) + eps, '-', color=BLUE, lw=1.5,
         label=r'$|P_{\rm ours}/P_{\rm ref}-1|$, $1024^3$')
ax3.plot(k_hi, np.abs(Po_hi/Pr_hi - 1) + eps, '-', color=ORANGE, lw=1.5,
         label=r'$|P_{\rm ours}/P_{\rm ref}-1|$, $2048^3$')
r_lo = Px_lo/np.sqrt(Po_lo*Pr_lo)
r_hi = Px_hi/np.sqrt(Po_hi*Pr_hi)
ax3.plot(k_lo, 1 - r_lo + eps, ':', color=BLUE, lw=1.5,
         label=r'$1-r(k)$, $1024^3$')
ax3.plot(k_hi, 1 - r_hi + eps, ':', color=ORANGE, lw=1.5,
         label=r'$1-r(k)$, $2048^3$')
ax3.set_yscale('log')
ax3.set_ylim(1e-9, 1e-2)
ax3.set_ylabel('deviation from reference IC')
ax3.set_xlabel(r'$k\ \ [h/\mathrm{Mpc}]$')
leg3 = ax3.legend(loc='upper left', fontsize=8.5, ncol=1, handlelength=2.0,
                  frameon=True, facecolor=SURF, edgecolor=GRID, framealpha=0.95)
leg3.set_zorder(10)
ax3.text(0.985, 0.97,
         'phases and amplitudes match the existing ICs\n'
         r'to $\lesssim 3\times10^{-4}$ ($1024^3$) and $\lesssim 10^{-6}$ ($2048^3$)',
         transform=ax3.transAxes, ha='right', va='top', fontsize=7.8, color=MUTED)

ax3.set_xlim(7e-3, 8)
for ax in (ax1, ax2, ax3):
    ax.tick_params(which='both', direction='out', length=3)

fig.savefig(f'{BASE}/pk_validation.png', dpi=200, bbox_inches='tight')
fig.savefig(f'{BASE}/pk_validation.pdf', bbox_inches='tight')
print('wrote pk_validation.png / .pdf')
