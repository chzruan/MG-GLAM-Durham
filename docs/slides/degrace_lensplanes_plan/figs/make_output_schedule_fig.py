#!/usr/bin/env python3
"""DEGRACE L1024 output schedule: the 42 targets on chi(a) and their spacing.

Top panel: comoving distance chi(a) for the fiducial Planck 2015 background
(flat LCDM, Omega_m = 0.3089, matter + Lambda only), with the 17 primary
redshifts and the 25 secondary stations chi = 256 k Mpc/h (k = 1..25).
Bottom panel: separation of neighbouring targets in a (per cent) against the
0.5 per cent merge tolerance of MG-GLAM's SetExactSteps, and the local PM
step size of the ladder (z_init = 49, da0 = 4e-4) for scale.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
import numpy as np  # noqa: E402
from scipy.integrate import quad  # noqa: E402

from house_style import PALETTE, save_figure, use_house_style  # noqa: E402

BLUE, TEAL, GOLD, ROSE = PALETTE
C_H0 = 2997.92458  # Mpc/h
OM = 0.3089
ZTOL = 0.5  # per cent, SetExactSteps merge tolerance

# Baojiu's 42 scale factors (Slack, slack-chat-degrace-prep.txt)
A_TARGETS = np.sort(np.array([
    0.09885636873988636, 0.11436339923360468, 0.13101192462219735,
    0.14880957763634076, 0.16776699787816512, 0.18789864408246054,
    0.2092238939643731, 0.2317683231209623, 0.25556525683790143,
    0.2806576241072646, 0.30710027529154543, 0.3333333333333333,
    0.33496251467163335, 0.34782608695652173, 0.36363636363636365,
    0.3643315849360659, 0.38095238095238093, 0.3953165138959508, 0.4,
    0.42105263157894735, 0.42805303395318656, 0.4444444444444444,
    0.4627096181121927, 0.47058823529411764, 0.49949502290382997, 0.5,
    0.5333333333333333, 0.5386678234812446, 0.5714285714285714,
    0.5805485755215236, 0.6153846153846154, 0.6255358719332695,
    0.6666666666666666, 0.6741267910946328, 0.7269451849577723,
    0.7272727272727273, 0.7847794661808203, 0.8, 0.8486352893278191,
    0.8888888888888888, 0.9198100312830662, 1.0,
]))
Z_PRIMARY = np.arange(0.0, 2.0 + 1e-9, 0.125)
A_PRIMARY = 1.0 / (1.0 + Z_PRIMARY)


def chi(a: float) -> float:
    f = lambda x: 1.0 / (x * x * np.sqrt(OM / x**3 + 1.0 - OM))  # noqa: E731
    return C_H0 * quad(f, a, 1.0, epsabs=0, epsrel=1e-11)[0]


def pm_ladder(zinit: float = 49.0, da0: float = 4e-4):
    a0 = 1.0 / (1.0 + zinit)
    sf = da0 / a0
    da, a, alist, dlist = da0, a0, [], []
    while True:
        if da < sf / 1.25 * a and a < 0.3:
            da *= 1.5
        a += da
        alist.append(a)
        dlist.append(da)
        if a >= 1.0:
            break
    return np.array(alist), np.array(dlist)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("output_schedule.pdf"))
    parser.add_argument("--no-tex", action="store_true")
    args = parser.parse_args()
    use_house_style(usetex=not args.no_tex)

    is_primary = np.array([np.any(np.isclose(a, A_PRIMARY, rtol=0, atol=1e-9)) for a in A_TARGETS])
    a_sec = A_TARGETS[~is_primary]
    a_pri = A_TARGETS[is_primary]
    chi_t = np.array([chi(a) for a in A_TARGETS])

    a_grid = np.linspace(0.085, 1.0, 400)
    chi_grid = np.array([chi(a) for a in a_grid])

    sep = 100.0 * (A_TARGETS[1:] - A_TARGETS[:-1]) / A_TARGETS[1:]
    a_mid = A_TARGETS[1:]
    bad = sep < ZTOL

    al, dl = pm_ladder()

    fig = plt.figure(figsize=(6.6, 4.9))
    gs = fig.add_gridspec(2, 1, height_ratios=(2.6, 1.6), hspace=0.0)
    ax = fig.add_subplot(gs[0])
    axb = fig.add_subplot(gs[1], sharex=ax)
    ax.tick_params(labelbottom=False)

    # --- top: chi(a) with the 42 targets
    for k in range(1, 26):
        ax.axhline(256.0 * k, color="0.85", lw=0.5, zorder=1)
    ax.plot(a_grid, chi_grid, color=BLUE, lw=1.6, zorder=3)
    ax.plot(a_sec, chi_t[~is_primary], "o", ms=4.2, mfc=GOLD, mec="k", mew=0.35, lw=0, zorder=5)
    ax.plot(a_pri, chi_t[is_primary], "s", ms=4.2, mfc=TEAL, mec="k", mew=0.35, lw=0, zorder=6)
    ax.set_ylabel(r"$\chi(a)\ [h^{-1}\,\mathrm{Mpc}]$")
    ax.set_ylim(-150, 6800)
    ax.set_xlim(0.07, 1.02)
    ax.text(0.98, 0.95, r"flat $\Lambda$CDM, $\Omega_m=0.3089$, matter $+\Lambda$ only",
            transform=ax.transAxes, ha="right", va="top", fontsize=9, color="0.3")
    ax.text(0.98, 0.86, r"grey lines: $\chi = 256\,k\ h^{-1}\mathrm{Mpc}$",
            transform=ax.transAxes, ha="right", va="top", fontsize=9, color="0.3")
    handles = [
        Line2D([], [], color=BLUE, lw=1.6, label=r"$\chi(a)$"),
        Line2D([], [], marker="s", lw=0, ms=4.2, mfc=TEAL, mec="k", mew=0.35,
               label=r"17 primary: $z = 0, 0.125, \ldots, 2$"),
        Line2D([], [], marker="o", lw=0, ms=4.2, mfc=GOLD, mec="k", mew=0.35,
               label=r"25 secondary: $\chi = 256\,k$, $k=1\ldots25$"),
    ]
    ax.legend(handles=handles, loc="lower left", fontsize=9)

    # secondary x axis in redshift (top)
    zt = np.array([9, 5, 3, 2, 1.5, 1, 0.5, 0.25, 0])
    axz = ax.secondary_xaxis("top", functions=(lambda a: 1.0 / a - 1.0, lambda z: 1.0 / (1.0 + z)))
    axz.set_xticks(zt)
    axz.set_xticklabels([f"{z:g}" for z in zt])
    axz.set_xlabel(r"$z$", fontsize=12)
    axz.tick_params(direction="in")

    # --- bottom: neighbour separation vs ZTOL
    axb.plot(al, 100.0 * dl / al, color="0.55", lw=0.9, ls="--", zorder=2)
    axb.axhspan(0.01, ZTOL, color="0.88", alpha=0.7, zorder=1)
    axb.axhline(ZTOL, color="k", lw=0.8, zorder=2)
    axb.plot(a_mid[~bad], sep[~bad], "o", ms=3.4, mfc="0.35", mec="k", mew=0.3, lw=0, zorder=4)
    axb.plot(a_mid[bad], sep[bad], "o", ms=5.0, mfc=ROSE, mec="k", mew=0.35, lw=0, zorder=5)
    for i, (a_b, s_b) in enumerate(zip(a_mid[bad], sep[bad])):
        left = i % 2 == 0          # alternate: left-below, right-above
        axb.annotate(f"{s_b:.3g}\\%" if not args.no_tex else f"{s_b:.3g}%",
                     (a_b, s_b), xytext=(-7, -7) if left else (7, 3),
                     textcoords="offset points", ha="right" if left else "left",
                     va="top" if left else "bottom", fontsize=8, color=ROSE)
    axb.set_yscale("log")
    axb.set_ylim(0.012, 600)
    axb.set_yticks([0.1, 1, 10, 100])
    axb.set_yticklabels(["0.1", "1", "10", "100"])
    axb.set_xlabel(r"scale factor $a$")
    axb.set_ylabel(r"$\Delta a/a$ to next target $[\%]$", fontsize=11)
    hb = [
        Line2D([], [], color="k", lw=0.8,
               label=(r"\texttt{ZTOL} $= 0.5\%$: pairs below are rejected" if not args.no_tex
                      else "ZTOL = 0.5%: pairs below are rejected")),
        Line2D([], [], color="0.55", lw=0.9, ls="--",
               label=r"PM step $\mathrm{d}a/a$ ($z_{\rm init}=49$, $\mathrm{d}a_0=4\times10^{-4}$)"),
        Line2D([], [], marker="o", lw=0, ms=3.4, mfc="0.35", mec="k", mew=0.3,
               label="neighbour separation"),
        Line2D([], [], marker="o", lw=0, ms=5.0, mfc=ROSE, mec="k", mew=0.35,
               label="the four close pairs"),
    ]
    axb.legend(handles=hb, loc="upper right", fontsize=7.5, ncol=2, handlelength=1.8,
               columnspacing=1.0, borderaxespad=0.3)
    axb.grid(True, ls=":", alpha=0.35)

    save_figure(fig, args.out)
    print(f"[plot] wrote {args.out}; {bad.sum()} pairs below ZTOL")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
