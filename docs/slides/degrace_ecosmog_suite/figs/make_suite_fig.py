#!/usr/bin/env python3
"""DEGRACE ECOSMOG suite: design coverage of the 96 Sobol nodes and the
per-node output stations chi = 256 k Mpc/h (data: docs/degrace_ecosmog/)."""
from __future__ import annotations
import argparse, csv, json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.cm import ScalarMappable  # noqa: E402
from matplotlib.colors import Normalize  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
import numpy as np  # noqa: E402
from house_style import PALETTE, save_figure, use_house_style  # noqa: E402

BLUE, TEAL, GOLD, ROSE = PALETTE
DOC = Path(__file__).resolve().parents[3] / "degrace_ecosmog"
PLANCK = np.sort(np.array([0.09885636873988636,0.11436339923360468,0.13101192462219735,0.14880957763634076,0.16776699787816512,0.18789864408246054,0.2092238939643731,0.2317683231209623,0.25556525683790143,0.2806576241072646,0.30710027529154543,0.33496251467163335,0.3643315849360659,0.3953165138959508,0.42805303395318656,0.4627096181121927,0.49949502290382997,0.5386678234812446,0.5805485755215236,0.6255358719332695,0.6741267910946328,0.7269451849577723,0.7847794661808203,0.8486352893278191,0.9198100312830662]))[::-1]
PRIM = 1.0 / (1.0 + np.arange(0.0, 2.0001, 0.125))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=Path("suite_overview.pdf"), help="directory anchor; two files are written beside it")
    ap.add_argument("--no-tex", action="store_true")
    a = ap.parse_args()
    use_house_style(usetex=not a.no_tex)
    rows = list(csv.DictReader(open(DOC / "nodes.csv")))
    aout = json.load(open(DOC / "aout.json"))
    om = np.array([float(r["Omega_m"]) for r in rows]); h = np.array([float(r["h"]) for r in rows])
    s8 = np.array([float(r["S_8"]) for r in rows]); ns = np.array([float(r["n_s"]) for r in rows])
    lf = np.array([float(r["log10_fR0"]) for r in rows])
    norm = Normalize(-6.5, -4.0); cmap = plt.get_cmap("viridis")

    # --- design coverage: two panels + a colour bar in its own axis
    fig = plt.figure(figsize=(6.6, 3.0))
    axes = [fig.add_axes([0.08, 0.17, 0.36, 0.78]), fig.add_axes([0.55, 0.17, 0.36, 0.78])]
    cax = fig.add_axes([0.925, 0.17, 0.018, 0.78])
    for ax, (x, y, xl, yl) in zip(axes, [(om, h, r"$\Omega_m$", r"$h$"), (s8, ns, r"$S_8$", r"$n_s$")]):
        ax.scatter(x, y, c=lf, cmap=cmap, norm=norm, s=22, edgecolors="k", linewidths=0.35, zorder=4)
        ax.scatter(x[:2], y[:2], s=90, facecolors="none", edgecolors=ROSE, linewidths=1.2, zorder=5)
        ax.annotate("1", (x[0], y[0]), xytext=(5, -9), textcoords="offset points", fontsize=8, color=ROSE)
        ax.annotate("2", (x[1], y[1]), xytext=(5, 4), textcoords="offset points", fontsize=8, color=ROSE)
        ax.set_xlabel(xl); ax.set_ylabel(yl); ax.grid(True, ls=":", alpha=0.35)
    axes[0].set_xlim(0.24, 0.36); axes[0].set_ylim(0.64, 0.76)
    axes[1].set_xlim(0.63, 0.97); axes[1].set_ylim(0.88, 1.07)
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax)
    cb.set_label(r"$\log_{10}|f_{R0}|$", fontsize=11)
    axes[0].text(0.03, 0.95, "96 Sobol nodes", transform=axes[0].transAxes, va="top", fontsize=9, color="0.3")
    save_figure(fig, a.out.with_name("design_coverage.pdf"))

    # --- output stations per node
    fig2, ax = plt.subplots(figsize=(3.9, 3.0))
    k = np.arange(1, 26)
    for node, al in aout.items():
        al = np.array(al)
        isprim = np.array([np.any(np.isclose(v, PRIM, rtol=0, atol=1e-9)) for v in al])
        st = np.sort(al[~isprim])[::-1]
        ax.plot(k, st, color="0.6", lw=0.6, alpha=0.6, zorder=2)
    ax.plot(k, PLANCK, color=BLUE, lw=1.8, zorder=4)
    ax.plot(k, PLANCK, "o", ms=3.4, mfc=GOLD, mec="k", mew=0.35, lw=0, zorder=5)
    ax.set_xlabel(r"station $k$ ($\chi = 256\,k\ h^{-1}$Mpc)"); ax.set_ylabel(r"$a$ at the station")
    ax.set_xlim(0.5, 25.5); ax.set_ylim(0.08, 0.95); ax.grid(True, ls=":", alpha=0.35)
    ax.legend(handles=[Line2D([], [], color="0.6", lw=0.8, label="96 nodes (own $\\Omega_m$)"),
                       Line2D([], [], color=BLUE, marker="o", mfc=GOLD, mec="k", mew=0.35, ms=3.4, lw=1.8,
                              label="Planck fiducial (Slack list)")], loc="upper right", fontsize=8.5)
    fig2.subplots_adjust(left=0.17, right=0.97, bottom=0.17, top=0.97)
    save_figure(fig2, a.out.with_name("station_epochs.pdf"))
    print("[plot] wrote design_coverage.pdf and station_epochs.pdf")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
