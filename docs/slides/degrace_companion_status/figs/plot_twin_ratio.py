#!/usr/bin/env python3
"""P_GLAM / P_ECOSMOG(f(R), POWMES) against k for the two node-2 pilots: LCDM (left) and f(R) (right)."""
import json, math, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(__file__))
from house_style import use_house_style

ROOT = "/cosma8/data/dp203/dc-ruan1/mgglam_claude/MG-GLAM"
RUNS = [(ROOT + "/DEGRACE-companion/runs/Node_002_LCDM_L1024Np1024Ng2048/Run1", r"$\Lambda$CDM pilot (job 12180082) vs the f(R) twin"),
        (ROOT + "/DEGRACE-companion/runs/Node_002_FR_L1024Np1024Ng2048/Run1", r"f(R) pilot (job 12180110) vs the f(R) twin")]
PK = "/cosma7/data/dp004/bl267/Runs/DEGRACE/Runs/L1024/FR/power_spectra/Node_002"
BOX, NGRID = 1024.0, 2048
aout = sorted(json.load(open(ROOT + "/docs/degrace_ecosmog/aout.json"))["Node_002"])
steps = [57, 90, 105, 127, 146, 161, 181]
kf, kny = 2 * math.pi / BOX, math.pi * NGRID / BOX


def ratio(run, st):
    sched = {int(l.split()[0]): float(l.split()[1]) for l in open(run + "/OutputSchedule.dat") if not l.startswith("#")}
    rows = [l.split() for l in open(run + "/PowerDM.R.%04d.0001.dat" % st) if len(l.split()) == 4 and l.split()[0].isdigit()]
    kg = np.array([float(r[1]) for r in rows]); pg = np.array([float(r[3]) for r in rows])
    pg = pg / (1 - 2.0 / 3.0 * np.sin(math.pi * kg / (2 * kny)) ** 2)   # shell-averaged CIC deconvolution
    a = sched[st]; i = int(np.argmin(np.abs(np.array(aout) - a)))
    E = np.loadtxt(PK + "/PK_FR_N002_S%02d.dat" % (i + 2)); ke = E[:, 0] * kf; pe = E[:, 3] * BOX ** 3
    s = (ke > 0) & (pe > 0)
    pe_at = np.exp(np.interp(np.log(kg), np.log(ke[s]), np.log(pe[s])))
    m = (kg >= 0.008) & (kg <= 1.5)
    return a, kg[m], pg[m] / pe_at[m]


use_house_style()
fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.4), sharey=True)
cmap = plt.get_cmap("viridis")
for ax, (run, title) in zip(axes, RUNS):
    for st in steps:
        a, k, r = ratio(run, st)
        ax.plot(k, r, color=cmap((a - 0.1) / 0.9), lw=1.4, label=r"$z=%.2f$" % (1 / a - 1))
    ax.axhspan(0.99, 1.01, color="0.85", zorder=0)
    ax.axhline(1.0, color="k", lw=0.8)
    ax.set_xscale("log"); ax.set_xlim(0.008, 1.5); ax.set_ylim(0.86, 1.14)
    ax.set_xlabel(r"$k\ [h\,\mathrm{Mpc}^{-1}]$")
    ax.text(0.0095, 1.12, title, fontsize=8.5)
axes[0].set_ylabel(r"$P_{\rm GLAM}(k)\,/\,P_{\rm ECOSMOG}(k)$")
axes[0].text(0.0095, 0.875, r"node 2, same ICs, $|f_{R0}| = 5.6\times10^{-6}$; grey band $\pm1\%$", fontsize=8)
axes[1].legend(fontsize=7.5, ncol=2, loc="lower left", frameon=False)
fig.tight_layout(w_pad=0.6)
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "twin_ratio")
fig.savefig(out + ".pdf"); fig.savefig(out + ".png", dpi=200)
print("wrote", out + ".pdf")
