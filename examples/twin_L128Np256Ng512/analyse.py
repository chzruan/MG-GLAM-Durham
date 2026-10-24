"""Read both twins with glamio and make the example's figure.

Works on native and packed runs alike: every file is found with glamio.product_path, so after step 5 the data
come straight out of Output.<step>.<tag>.h5. Panels: (a) P(k) of both twins at the full outputs, (b) the f(R)
enhancement P_fR/P_LCDM - 1 at every output, (c) the halo mass function ratio, (d) a lens plane at the lowest
station redshift, (e) f_R at the halo centres against halo mass (chameleon screening). Same ICs, so the ratios
carry almost no cosmic variance. Writes $WORK/figures/showcase.png and .pdf and prints a short summary.
Run: source config.sh && $PYTHON analyse.py  (or through run_all.sh); config.sh exports the settings read here.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

import glamio  # noqa: E402

WORK = os.environ["WORK"]
BOX = float(os.environ.get("BOX", 128))
TAG = int(os.environ.get("REALISATION", 0))
REF, MG = os.environ.get("TWINS", "LCDM fR5").split()[:2]   # the LCDM reference and the f(R) twin
RUN = {twin: os.path.join(WORK, twin, f"Run{TAG}") for twin in (REF, MG)}

try:                                                       # the repository's plotting style, when available
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "python", "examples"))
    from house_style import PALETTE, use_house_style
    use_house_style(usetex=False)
except Exception:                                          # noqa: BLE001
    PALETTE = ("#1d4e89", "#20b2aa", "#edae49", "#d1495b")

sched = glamio.read_schedule(glamio.product_path(RUN[REF], "schedule", 0, TAG))
full = [int(r["step"]) for r in sched.primary()]
stations = [int(r["step"]) for r in sched.stations()]
zof = {int(r["step"]): float(r["z"]) for r in sched.rows}
order = sorted(zof, key=zof.get)                           # one colour per redshift, low z yellow to high z purple
colour = {s: plt.cm.viridis(0.9 - 0.85 * i / max(len(order) - 1, 1)) for i, s in enumerate(order)}


def pk(twin, step):
    p = glamio.read_powerdm(glamio.product_path(RUN[twin], "PowerDM.R", step, TAG), box=BOX)
    ok = (p.nmodes > 0) & (p.k < 0.5 * p.k_nyq)            # above, the 1-D CIC window correction grows unreliable
    return p.k[ok], p.pk_deconvolved()[ok]


fig, ax = plt.subplots(2, 3, figsize=(13.5, 8.2))
ax = ax.ravel()
print(f"outputs: full z = {', '.join(f'{zof[s]:g}' for s in full)}; stations z = {', '.join(f'{zof[s]:g}' for s in stations)}")

# (a) P(k) at the full outputs, (b) the enhancement at every output
for s in full:
    k, pl = pk(REF, s)
    ax[0].loglog(k, pl, color=colour[s], label=f"z = {zof[s]:g}")
    ax[0].loglog(*pk(MG, s), color=colour[s], ls="--")
for s in sorted(zof, key=zof.get):
    k, pl = pk(REF, s)
    kf, pf = pk(MG, s)
    boost = pf / pl - 1.0
    ax[1].semilogx(k, 100 * boost, color=colour[s], ls="-" if s in full else ":", label=f"z = {zof[s]:g}")
    sel = [np.argmin(np.abs(k - kk)) for kk in (0.5, 1.0, 3.0)]
    print(f"z = {zof[s]:4.2f}: P_fR/P_LCDM - 1 = " + ", ".join(f"{100 * boost[i]:+.1f} % at k = {k[i]:.2f}" for i in sel))
ax[0].set(xlabel=r"$k\ [h\,\mathrm{Mpc}^{-1}]$", ylabel=r"$P(k)\ [(\mathrm{Mpc}/h)^3]$", title="(a) P(k): LCDM solid, f(R) dashed")
ax[1].axhline(0, color="k", lw=0.8)
ax[1].set(xlabel=r"$k\ [h\,\mathrm{Mpc}^{-1}]$", ylabel=r"$P_{f(R)}/P_{\Lambda\mathrm{CDM}} - 1\ [\%]$",
          title="(b) f(R) enhancement (dotted: stations)")
ax[0].legend(frameon=False); ax[1].legend(frameon=False)

# (c) halo mass function ratio at the full outputs
bins = np.arange(11.5, 15.01, 0.3)
mid = 0.5 * (bins[1:] + bins[:-1])
for s in full:
    lm = {t: np.log10(np.asarray(glamio.read_catalogue(glamio.product_path(RUN[t], "catalogue", s, TAG), mmap=True)["Mbound"],
                                 dtype=np.float64)) for t in RUN}
    n = {t: np.histogram(lm[t], bins=bins)[0] for t in RUN}
    ok = n[REF] >= 30                                    # Poisson noise below 20 %
    ax[2].plot(mid[ok], n[MG][ok] / n[REF][ok], "o-", color=colour[s], ms=4, label=f"z = {zof[s]:g}")
    print(f"z = {zof[s]:4.2f}: haloes above 10^{bins[0]} Msun/h: LCDM {np.sum(lm[REF] >= bins[0])}, "
          f"f(R) {np.sum(lm[MG] >= bins[0])}")
ax[2].axhline(1, color="k", lw=0.8)
ax[2].set(xlabel=r"$\log_{10} M_\mathrm{bound}\ [h^{-1}M_\odot]$", ylabel=r"$n_{f(R)}/n_{\Lambda\mathrm{CDM}}$",
          title="(c) halo mass function ratio")
ax[2].legend(frameon=False)

# (d) a lens plane at the lowest-z station (Sigma/Sigma_bar of the Z slab)
s = min(stations or full, key=zof.get)
lp = glamio.read_lensplanes(glamio.product_path(RUN[REF], "LensPlanes", s, TAG), mmap=True)
img = np.log10(np.asarray(lp.plane("z", 1), dtype=np.float64))
im = ax[3].imshow(img.T, origin="lower", extent=[0, lp.box, 0, lp.box], cmap="magma", vmin=-0.6, vmax=1.2)
fig.colorbar(im, ax=ax[3], label=r"$\log_{10}\Sigma/\bar\Sigma$", shrink=0.9)
ax[3].set(xlabel=r"$x\ [h^{-1}\mathrm{Mpc}]$", ylabel=r"$y\ [h^{-1}\mathrm{Mpc}]$",
          title=f"(d) lens plane, {REF}, z = {zof[s]:g}")
depth = lp.depth

# (e) chameleon screening: f_R at halo centres over the background value, at z = 0
s = min(full, key=zof.get)
fr = glamio.read_frhalo(glamio.product_path(RUN[MG], "fR_halo", s, TAG))
cat = fr.catalogue()                                       # the catalogue of the same step, same rows
fbg = fr.fR_bg_old                                         # background f_R at the epoch the field was solved
lm, u = np.log10(np.asarray(cat["Mbound"], dtype=np.float64)), np.asarray(fr["fR_centre"], dtype=np.float64) / fbg
ax[4].plot(lm, u, ",", color=PALETTE[0], alpha=0.3)
med = [np.median(u[(lm >= lo) & (lm < hi)]) if np.sum((lm >= lo) & (lm < hi)) >= 20 else np.nan for lo, hi in zip(bins[:-1], bins[1:])]
ax[4].plot(mid, med, "o-", color=PALETTE[3], ms=4, label="median")
ax[4].set(xlabel=r"$\log_{10} M_\mathrm{bound}\ [h^{-1}M_\odot]$", ylabel=r"$f_R(\mathrm{centre})/\bar f_R$",
          ylim=(0, 1.2), title=f"(e) screening in f(R), z = {zof[s]:g}")
ax[4].legend(frameon=False)

# (f) what was read, and from where
ax[5].axis("off")
src = glamio.product_path(RUN[MG], "catalogue", s, TAG)
ax[5].text(0, 1, "\n".join([
    f"L = {BOX:g} Mpc/h, {os.environ.get('NP', '256')}^3 particles, {os.environ.get('NG', '512')}^3 grid",
    f"LCDM and f(R) (n = 1, fR0 = {os.environ.get('FR0', '-1e-5')})",
    f"from the same 2LPT ICs at z = {os.environ.get('ZINIT', '49')}",
    f"full outputs: z = {', '.join(f'{zof[x]:g}' for x in full)}",
    f"stations: z = {', '.join(f'{zof[x]:g}' for x in stations) or 'none'}",
    f"lens planes {depth:g} Mpc/h deep",
    "",
    "read with glamio.product_path from",
    f"{'packed HDF5 files' if '.h5/' in src else 'native files'}, e.g.",
    f"{os.path.relpath(glamio.split_packed(src)[0] if '.h5/' in src else src, WORK)}",
    f"  /{glamio.split_packed(src)[1] if '.h5/' in src else ''}",
]), va="top", family="monospace", fontsize=8.5, transform=ax[5].transAxes)

fig.tight_layout()
out = os.path.join(WORK, "figures")
os.makedirs(out, exist_ok=True)
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(out, f"showcase.{ext}"), dpi=150)
print(f"figure: {out}/showcase.png")
