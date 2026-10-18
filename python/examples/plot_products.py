#!/usr/bin/env python3
"""Figures of GLAM's on-the-fly products from the npz files written by measure_products.py.

usage: plot_products.py [--data DIR] [--out DIR] [fig ...]
with fig in pk, modes, fields, planes, velocities, haloes (default: every figure whose
data file exists).  --data (default: the current directory) holds products_fast.npz,
products_cat.npz and products_planes.npz; the figures go to --out (default: --data) as
products_<fig>.pdf and .png.  Uses the house style of house_style.py beside this script.
"""
import argparse
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from house_style import PALETTE, save_figure, style_residual_axis, use_house_style  # noqa: E402

ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
ap.add_argument("figs", nargs="*")
ap.add_argument("--data", default=".", help="directory with products_*.npz")
ap.add_argument("--out", default=None, help="directory for the figures (default: --data)")
ARGS = ap.parse_args()
DATA = os.path.abspath(ARGS.data)
HERE = os.path.abspath(ARGS.out or DATA)
F = np.load(DATA + "/products_fast.npz")
C = np.load(DATA + "/products_cat.npz") if os.path.exists(DATA + "/products_cat.npz") else None
OM, OL = float(F["omega_m"]), float(F["omega_l"])
BOX = float(F["box"]) if "box" in F else 1024.0
NGRID, NDEP = 2048, 1024
cmap = plt.get_cmap("viridis")
trapz = getattr(np, "trapezoid", None) or np.trapz


def acol(a):
    """Colour of an epoch: viridis over a = 0.30 (z = 2.3) .. 1.02."""
    return cmap(float(np.clip((a - 0.30) / 0.72, 0.0, 1.0)))


def growth(a):
    """Linear growth factor D(a)/D(1) and rate f = dlnD/dlna of the flat LCDM background."""
    a = np.atleast_1d(np.asarray(a, float))
    ai = np.linspace(1e-3, 1.0, 40000)
    E = np.sqrt(OM / ai ** 3 + OL)
    g = 1.0 / (ai * E) ** 3
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (g[1:] + g[:-1]) * np.diff(ai))])
    cum += ai[0] ** 2.5 / 2.5 / OM ** 1.5 * 0  # the a < 1e-3 part is negligible (matter era: D ~ a)
    D = E * cum
    D /= D[-1]
    f = np.gradient(np.log(D[1:]), np.log(ai[1:]))
    return np.interp(a, ai, D), np.interp(a, ai[1:], f)


def H_kms(a):
    """H(a) in km/s per Mpc/h."""
    return 100.0 * np.sqrt(OM / a ** 3 + OL)


def cic_shell(k, kny):
    """Shell-averaged CIC window |W|^2 (Jing's alias-summed form for p = 2)."""
    return 1.0 - 2.0 / 3.0 * np.sin(np.pi * k / (2.0 * kny)) ** 2


def index_at_z(z, tag="lcdm"):
    zs = F[f"sched_{tag}_z"]
    i = int(np.argmin(np.abs(zs - z)))
    assert abs(zs[i] - z) < 2e-3, (z, zs[i])
    return i


def epoch_colorbar(fig, ax, zticks=(2.0, 1.0, 0.5, 0.0)):
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0.30, 1.02))
    cb = fig.colorbar(sm, ax=ax, pad=0.02, aspect=28)
    cb.set_ticks([1.0 / (1.0 + z) for z in zticks])
    cb.set_ticklabels([rf"$z={z:g}$" for z in zticks])
    cb.ax.tick_params(labelsize=8.5)
    return cb


# ----------------------------------------------------------------------------- figures
def fig_pk():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.8, 3.7))
    k, kny, P, A, cls = F["pk_lcdm_k"], float(F["pk_lcdm_knyq"]), F["pk_lcdm_p"], F["pk_lcdm_a"], F["sched_lcdm_class"]
    w = cic_shell(k, kny)
    m = (k > 0) & (k < 0.8 * kny)
    for i in np.where(cls == 1)[0]:
        ax1.plot(k[m], (P[i] / w)[m], color=acol(A[i]), lw=1.0)
    kl, pl = F["pk_lin_k"], F["pk_lin_p"]
    for z, ls in ((0.0, "--"), (2.0, ":")):
        D, _ = growth(1.0 / (1.0 + z))
        ax1.plot(kl, pl * D[0] ** 2, color="k", ls=ls, lw=1.0, label=rf"linear theory, $z = {z:g}$")
    ax1.set_xscale("log"); ax1.set_yscale("log")
    ax1.set_xlim(0.006, 6.0); ax1.set_ylim(2.0, 1.6e5)
    ax1.set_xlabel(r"$k\ [h\,\mathrm{Mpc}^{-1}]$"); ax1.set_ylabel(r"$P(k)\ [(h^{-1}\mathrm{Mpc})^3]$")
    ax1.legend(loc="lower left", fontsize=9)
    ax1.text(0.0075, 8.5e4, r"\texttt{PowerDM.R}, $\Lambda$CDM, the 17 primary outputs", fontsize=8.5)
    ax1.text(0.0075, 4.6e4, r"(CIC window removed; colour = epoch)", fontsize=8.5)
    epoch_colorbar(fig, ax1)
    for z in (2.0, 1.0, 0.5, 0.0):
        i, j = index_at_z(z, "lcdm"), index_at_z(z, "fr")
        assert F["sched_lcdm_step"][i] == F["sched_fr_step"][j]
        kf_, pf = F["pk_fr_k"], F["pk_fr_p"][j]
        assert np.allclose(kf_, k)
        ax2.plot(k[m], (pf / P[i])[m], color=acol(A[i]), lw=1.5, label=rf"$z = {z:g}$")
    ax2.axhline(1.0, color="k", lw=0.8)
    ax2.set_xscale("log"); ax2.set_xlim(0.006, 6.0); ax2.set_ylim(0.97, 1.33)
    ax2.set_xlabel(r"$k\ [h\,\mathrm{Mpc}^{-1}]$"); ax2.set_ylabel(r"$P_{f(R)}(k)\,/\,P_{\Lambda\mathrm{CDM}}(k)$")
    ax2.legend(loc="upper left", fontsize=9)
    ax2.text(0.0075, 0.985, r"f(R) pilot / $\Lambda$CDM run, same ICs; $|f_{R0}| = 5.6\times10^{-6}$, $n = 1$", fontsize=8.5)
    fig.tight_layout(w_pad=1.0)
    save_figure(fig, [HERE + "/products_pk.pdf", HERE + "/products_pk.png"])


def fig_modes():
    fig = plt.figure(figsize=(9.8, 4.0))
    gs = fig.add_gridspec(2, 2, height_ratios=(3.0, 1.3), hspace=0.0, wspace=0.30)
    ax = fig.add_subplot(gs[0, 0]); axr = fig.add_subplot(gs[1, 0], sharex=ax); ax2 = fig.add_subplot(gs[:, 1])
    ax.tick_params(labelbottom=False)
    kt, pt, zs = F["pk_lcdm_k"], F["pk_lcdm_p"], F["sched_lcdm_z"]
    for z in (2.0, 0.0):
        kc, pc = F[f"cube_dd_k_z{z:g}"], F[f"cube_dd_p_z{z:g}"]
        i = index_at_z(z)
        n = len(kc)
        assert np.allclose(kc, kt[:n], rtol=5e-3)   # the table prints k with 5 decimals
        col = acol(1.0 / (1.0 + z))
        ax.plot(kc, pc, color=col, lw=1.8, label=rf"\texttt{{deltam}} cube, $z = {z:g}$")
        ax.plot(kt[:n], pt[i][:n], color="k", ls="--", lw=0.9, label=r"\texttt{PowerDM.R} table" if z == 0 else None)
        axr.plot(kc, 1e4 * (pc / pt[i][:n] - 1.0), color=col, lw=1.2)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(0.006, 1.6); ax.set_ylim(20, 1.6e5)
    ax.set_ylabel(r"$P(k)\ [(h^{-1}\mathrm{Mpc})^3]$")
    ax.legend(loc="lower left", fontsize=9)
    ax.text(0.0075, 8.5e4, r"$n = 244$ shells up to $k_c = 1.5\,h\,\mathrm{Mpc}^{-1}$", fontsize=8.5)
    ax.text(0.0075, 4.6e4, r"raw CIC on $2048^3$; the table's own shells", fontsize=8.5)
    style_residual_axis(axr, band=(-1.0, 1.0), ylim=(-3.0, 3.0))
    axr.set_ylabel(r"$10^4\,(P_{\rm cube}/P_{\rm table} - 1)$", fontsize=9)
    axr.set_xlabel(r"$k\ [h\,\mathrm{Mpc}^{-1}]$")
    rows = F["cube_rows"]
    for z in (2.0, 1.0, 0.0):
        k, ptt, pdd = F[f"cube_tt_k_z{z:g}"], F[f"cube_tt_p_z{z:g}"], F[f"cube_tt_pdd_z{z:g}"]
        av = float(rows[np.argmin(np.abs(rows[:, 0] - z)), 3])
        D, f = growth(av)
        fac = (av * H_kms(av) * f[0]) ** 2
        ratio = (ptt / cic_shell(k, np.pi * NDEP / BOX)) / (fac * pdd / cic_shell(k, np.pi * NGRID / BOX))
        ax2.plot(k, ratio, color=acol(av), lw=1.7, label=rf"$z = {z:g}$ ($f = {f[0]:.2f}$)")
    ax2.axhline(1.0, color="k", lw=0.8)
    ax2.set_xscale("log"); ax2.set_yscale("log"); ax2.set_xlim(0.006, 0.8); ax2.set_ylim(0.5, 40.0)
    ax2.set_xlabel(r"$k\ [h\,\mathrm{Mpc}^{-1}]$")
    ax2.set_ylabel(r"$P_{\theta\theta}(k)\,/\,[a H f]^2 P_{\delta\delta}(k)$")
    ax2.legend(loc="lower right", fontsize=9)
    ax2.text(0.03, 0.95, r"$\theta_{\mathbf k} = i\,\mathbf{k}\cdot\mathbf{p}_{\mathbf k}$ from the \texttt{mom\_x,y,z} cubes ($N_{\rm dep} = 1024$, $k_c = 0.75$)", fontsize=8.5, transform=ax2.transAxes, va="top")
    ax2.text(0.03, 0.88, r"linear theory $= 1$; the rise beyond $k \approx 0.1$ is the density-weighted", fontsize=8.5, transform=ax2.transAxes, va="top")
    ax2.text(0.03, 0.81, r"term $\nabla\!\cdot(\delta\mathbf{v})$ of the momentum $\mathbf{p} = (1+\delta)\mathbf{v}$", fontsize=8.5, transform=ax2.transAxes, va="top")
    fig.tight_layout(w_pad=1.0)
    save_figure(fig, [HERE + "/products_modes.pdf", HERE + "/products_modes.png"])


def fig_fields():
    fig, axes = plt.subplots(1, 3, figsize=(12.4, 4.5), sharey=True, constrained_layout=True)
    ext = [0.0, BOX, 0.0, BOX]
    norm = LogNorm(0.25, 25.0)
    plane = F["plane_z0"]
    axes[0].imshow(plane.T, origin="lower", extent=ext, cmap="viridis", norm=norm, interpolation="nearest")
    proj = 1.0 + F["field_z0_proj"]
    im1 = axes[1].imshow(np.clip(proj, 0.25, None).T, origin="lower", extent=ext, cmap="viridis", norm=norm, interpolation="nearest")
    xy, vz = F["part_z0_xy"], F["part_z0_vz"]
    sc = axes[2].scatter(xy[:, 0], xy[:, 1], c=vz, s=0.12, cmap="RdBu_r", vmin=-600, vmax=600, lw=0, rasterized=True)
    axes[2].set_xlim(0, BOX); axes[2].set_ylim(0, BOX); axes[2].set_aspect("equal")
    depth = float(F["plane_depth"])
    titles = [r"\texttt{LensPlanes}: $\Sigma/\bar\Sigma$" + "\n" + rf"$z$-axis slab 1, {depth:.0f} $h^{{-1}}$Mpc deep ($2048^2$ pixels, shown at $512^2$)",
              r"\texttt{deltam} cube: $1 + \langle\delta\rangle$ over the same slab" + "\n" + r"band-limited to $k \le 1.5\,h\,\mathrm{Mpc}^{-1}$ ($489^3$ grid)",
              rf"\texttt{{Subsample}}: 1 in {int(F['sub_m'])} particles" + "\n" + rf"slab $z < {float(F['part_slab']):.0f}\,h^{{-1}}$Mpc, colour $= v_z$"]
    for ax, t in zip(axes, titles):
        ax.set_title(t, fontsize=9.5)
        ax.set_xlabel(r"$x\ [h^{-1}\mathrm{Mpc}]$")
        ax.set_xticks([0, 256, 512, 768, 1024]); ax.set_yticks([0, 256, 512, 768, 1024])
    axes[0].set_ylabel(r"$y\ [h^{-1}\mathrm{Mpc}]$")
    cb = fig.colorbar(im1, ax=axes[:2], pad=0.01, aspect=30, shrink=0.85)
    cb.set_label(r"$\Sigma/\bar\Sigma$ \quad $1 + \langle\delta\rangle$")
    cb2 = fig.colorbar(sc, ax=axes[2], pad=0.01, aspect=30, shrink=0.85)
    cb2.set_label(r"$v_z\ [\mathrm{km\,s^{-1}}]$")
    save_figure(fig, [HERE + "/products_fields.pdf", HERE + "/products_fields.png"])


def fig_planes():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.8, 3.7))
    z, cls, a = F["sched_lcdm_z"], F["sched_lcdm_class"], F["sched_lcdm_a"]
    P8 = np.load(DATA + "/products_planes.npz") if os.path.exists(DATA + "/products_planes.npz") else None
    sd, pix = (P8["plane_std8"], float(P8["pixel8"])) if P8 is not None else (F["plane_std"], 0.5)
    x = 1.0 + z
    for ax_i, (col, lab) in enumerate(zip(PALETTE[:3], ("x", "y", "z"))):
        for s_ in range(4):
            ax1.plot(x, sd[:, 4 * ax_i + s_], color=col, lw=0.8, alpha=0.7, label=rf"{lab}-axis planes (4 slabs)" if s_ == 0 else None)
    mean = sd.mean(axis=1)
    ax1.plot(x[cls == 1], mean[cls == 1], "o", ms=4.2, color="k", mfc="k", lw=0, label="primary outputs (17)")
    ax1.plot(x[cls == 2], mean[cls == 2], "o", ms=4.2, color="k", mfc="w", lw=0, label="stations (25)")
    D, _ = growth(a)
    ref = index_at_z(2.0)
    ax1.plot(x, mean[ref] * D / D[ref], "k--", lw=0.9, label=r"linear growth $\propto D(a)$, matched at $z = 2$")
    ax1.set_xscale("log"); ax1.set_yscale("log"); ax1.set_xlim(10.5, 0.97)
    zt = [0, 0.5, 1, 2, 4, 8]
    ax1.set_xticks([1 + t for t in zt]); ax1.set_xticklabels([f"{t:g}" for t in zt]); ax1.minorticks_off()
    ax1.set_xlabel(r"$z$"); ax1.set_ylabel(rf"$\sigma(\Sigma/\bar\Sigma)$ per {pix:g} $h^{{-1}}$Mpc pixel")
    ax1.legend(loc="lower right", fontsize=8.5)
    ax1.text(0.03, 0.95, r"12 planes per output: 4 slabs of 256 $h^{-1}$Mpc per axis", fontsize=8.5, transform=ax1.transAxes, va="top")
    e = F["plane_pdf_edges"]; xc = 0.5 * (e[1:] + e[:-1])
    for zz in (2.0, 1.0, 0.0):
        pdf = F[f"plane_pdf_z{zz:g}"]
        ax2.plot(xc[pdf > 0], pdf[pdf > 0], color=acol(1.0 / (1.0 + zz)), lw=1.7, label=rf"$z = {zz:g}$")
    ax2.set_yscale("log"); ax2.set_xlim(-1.2, 1.8); ax2.set_ylim(1e-5, 10)
    ax2.set_xlabel(r"$\log_{10}\,\Sigma/\bar\Sigma$"); ax2.set_ylabel(r"pixel PDF [per dex]")
    ax2.legend(loc="upper right", fontsize=9)
    ax2.text(0.03, 0.95, r"$z$-axis planes, $4\times2048^2$ pixels of 0.5 $h^{-1}$Mpc", fontsize=8.5, transform=ax2.transAxes, va="top")
    fig.tight_layout(w_pad=1.0)
    save_figure(fig, [HERE + "/products_planes.pdf", HERE + "/products_planes.png"])


def fig_velocities():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.8, 3.7))
    av, vr = F["sub_av"], F["sub_vrms"]
    zv = 1.0 / av - 1.0
    for j, lab in enumerate("xyz"):
        ax1.plot(zv, vr[:, j], "o", ms=3.8, mec="k", mew=0.35, color=PALETTE[j], lw=0, label=rf"$\sigma_{{v_{lab}}}$ (1 in 100 particles)")
    kl, pl = F["pk_lin_k"], F["pk_lin_p"]
    I = trapz(pl, kl)
    zg = np.linspace(0.0, 2.2, 120); ag = 1.0 / (1.0 + zg)
    D, f = growth(ag)
    sig = ag * H_kms(ag) * f * D * np.sqrt(I / (6.0 * np.pi ** 2))
    ax1.plot(zg, sig, "k--", lw=1.0, label=r"linear theory: $a H f D\,[\int P_{\rm lin}\,{\rm d}k/6\pi^2]^{1/2}$")
    ax1.set_xlim(2.1, -0.05)
    ax1.set_xlabel(r"$z$ (velocity epoch $a_v$)"); ax1.set_ylabel(r"1-D peculiar velocity rms [km s$^{-1}$]")
    ax1.legend(loc="upper left", fontsize=8.5)
    e = F["vel_pdf_edges"]; xc = 0.5 * (e[1:] + e[:-1])
    for zz in (2.0, 1.0, 0.0):
        pdf = F[f"vel_pdf_z{zz:g}"]
        ax2.plot(xc[pdf > 0], pdf[pdf > 0], color=acol(1.0 / (1.0 + zz)), lw=1.7, label=rf"$z = {zz:g}$")
    s0 = float(vr[np.argmin(np.abs(zv)), 0])
    ax2.plot(xc, np.exp(-0.5 * (xc / s0) ** 2) / (s0 * np.sqrt(2 * np.pi)), "k--", lw=0.9, label=rf"Gaussian, $\sigma = {s0:.0f}$ km/s")
    ax2.set_yscale("log"); ax2.set_xlim(-2200, 2200); ax2.set_ylim(1e-8, 1e-2)
    ax2.set_xlabel(r"$v_x\ [\mathrm{km\,s^{-1}}]$"); ax2.set_ylabel(r"PDF [s km$^{-1}$]")
    ax2.legend(loc="upper right", fontsize=8.5)
    ax2.text(0.03, 0.04, r"\texttt{Subsample}: $10.7$ M particles, int16 at 0.25 km/s", fontsize=8.5, transform=ax2.transAxes)
    fig.tight_layout(w_pad=1.0)
    save_figure(fig, [HERE + "/products_velocities.pdf", HERE + "/products_velocities.png"])


def fig_haloes():
    fig = plt.figure(figsize=(9.8, 4.2))
    gs = fig.add_gridspec(2, 2, height_ratios=(3.0, 1.3), hspace=0.0, wspace=0.30)
    ax = fig.add_subplot(gs[0, 0]); axr = fig.add_subplot(gs[1, 0], sharex=ax); ax2 = fig.add_subplot(gs[:, 1])
    ax.tick_params(labelbottom=False)
    M = C["hmf_edges"]
    for z in (2.0, 1.0, 0.5, 0.0):
        nl, nf = C[f"hmf_lcdm_z{z:g}"], C[f"hmf_fr_z{z:g}"]
        col = acol(1.0 / (1.0 + z))
        ok = nl * BOX ** 3 >= 20
        ax.plot(M[ok], nl[ok], color=col, lw=1.6, label=rf"$z = {z:g}$")
        ax.plot(M[ok], nf[ok], color=col, lw=1.6, ls="--")
        good = nl * BOX ** 3 >= 200
        axr.plot(M[good], 100.0 * (nf[good] / nl[good] - 1.0), color=col, lw=1.3)
    ax.plot([], [], color="0.3", lw=1.6, label=r"$\Lambda$CDM (solid)")
    ax.plot([], [], color="0.3", lw=1.6, ls="--", label=r"f(R) (dashed)")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(1e12, 3e15); ax.set_ylim(1e-8, 1e-2)
    ax.set_ylabel(r"$n(>M)\ [(h^{-1}\mathrm{Mpc})^{-3}]$")
    ax.legend(loc="upper right", fontsize=8.5, ncol=2, title=r"BDM \texttt{CatshortV}, $M_{\rm tot}$ within $R_{\rm vir}$", title_fontsize=8.5)
    style_residual_axis(axr, band=(-2.0, 2.0), ylim=(-12, 45))
    axr.set_ylabel(r"$n_{f(R)}/n_{\Lambda} - 1$ [\%]", fontsize=9)
    axr.set_xlabel(r"$M_{\rm vir}\ [h^{-1}M_\odot]$")
    e = C["fr_mass_edges"]; mc = 10 ** (0.5 * (e[1:] + e[:-1]))
    for z in (2.0, 1.0, 0.5, 0.0):
        R = C[f"fr_fR_centre_z{z:g}"]
        ok = R[:, 0] >= 20
        col = acol(1.0 / (1.0 + z))
        ax2.plot(mc[ok], 10 ** R[ok, 2], color=col, lw=1.8, label=rf"$z = {z:g}$")
        ax2.fill_between(mc[ok], 10 ** R[ok, 1], 10 ** R[ok, 3], color=col, alpha=0.18, lw=0)
    ax2.axhline(1.0, color="k", lw=0.8)
    ax2.set_xscale("log"); ax2.set_yscale("log"); ax2.set_xlim(1e12, 3e15); ax2.set_ylim(2e-4, 2.0)
    ax2.set_xlabel(r"$M_{\rm vir}\ [h^{-1}M_\odot]$")
    ax2.set_ylabel(r"$|f_R(\mathrm{centre})|\,/\,|\bar f_R(a)|$")
    ax2.legend(loc="upper right", fontsize=9, title=rf"\texttt{{fR\_halo}}, $f_{{R0}} = {float(C['fr_fR0']) * 1e6:.1f}\times10^{{-6}}$, $n = {int(C['fr_n'])}$" + "\n" + r"median and 16--84\,\% band per 0.2 dex", title_fontsize=8.5)
    fig.tight_layout(w_pad=1.0)
    save_figure(fig, [HERE + "/products_haloes.pdf", HERE + "/products_haloes.png"])


FIGS = {"pk": fig_pk, "modes": fig_modes, "fields": fig_fields, "planes": fig_planes,
        "velocities": fig_velocities, "haloes": fig_haloes}

if __name__ == "__main__":
    use_house_style()
    os.makedirs(HERE, exist_ok=True)
    names = ARGS.figs or [n for n in FIGS if n != "haloes" or C is not None]
    for name in names:
        if name == "haloes" and C is None:
            print("haloes: products_cat.npz missing, skipped"); continue
        FIGS[name]()
        plt.close("all")
        print("wrote", f"products_{name}.pdf/.png")
