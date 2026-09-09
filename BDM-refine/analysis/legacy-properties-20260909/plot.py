"""Render measured halo statistics; no catalogue/density measurement here."""
from __future__ import annotations

import os
from pathlib import Path
import sys
import tempfile

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
STYLE = REPO/'BDM-refine/validation/20260908-convergence/plots'
sys.path.insert(0, str(STYLE))


def main():
    # TeX and Matplotlib caches are task-owned and removed after rendering.
    with tempfile.TemporaryDirectory(prefix='bdm-properties-mpl-') as cache:
        os.environ['MPLCONFIGDIR'] = cache
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
        from matplotlib.lines import Line2D
        from house_style import use_house_style, PALETTE
        from measure import sha, write_json
        use_house_style()
        plt.rcParams.update({'font.size': 10, 'axes.labelsize': 11,
                             'xtick.labelsize': 9, 'ytick.labelsize': 9,
                             'legend.fontsize': 9})
        figures = HERE/'figures'
        figures.mkdir(exist_ok=True)
        catalogue = np.load(HERE/'catalogue-statistics.npz')
        epochs = [np.load(HERE/f'F-z{z}-statistics.npz') for z in (0, 1, 2)]
        selections = ('mass12p5', 'mass13')
        colors = PALETTE[:2]
        outputs = []
        page_names = []
        limits_checks = []
        method_handles = [Line2D([], [], color='.25', lw=1.1, marker='o', ms=3,
                                 mec='k', mew=.35, label='Refined (v3)'),
                          Line2D([], [], color='.25', lw=1.25, ls='--', label='Legacy')]
        sample_handles = [Line2D([], [], color=c, lw=1.5, label=label) for c, label in zip(colors,
            [r'$M_{\rm bound}\geq 10^{12.5}\,h^{-1}M_\odot$',
             r'$M_{\rm bound}\geq 10^{13}\,h^{-1}M_\odot$'])]

        def line(ax, x, y, color, variant, mask=None):
            y = np.asarray(y, dtype=float).copy()
            if mask is not None:
                y[~mask] = np.nan
            if variant == 'v3':
                ax.plot(x, y, color=color, lw=.95, marker='o', ms=3,
                        mec='k', mew=.35, zorder=3)
            else:
                ax.plot(x, y, color=color, lw=1.35, ls='--', zorder=2)

        def residual(ax, x, y, color, mask=None, percent=True):
            values = np.asarray(y, dtype=float).copy()
            if mask is not None:
                values[~mask] = np.nan
            ax.plot(x, values, color=color, lw=.9, marker='o', ms=2.6, mec='k', mew=.3)
            ax.axhline(0, color='k', lw=.65, zorder=0)
            if percent:
                ax.axhspan(-1, 1, color='.9', zorder=-1)
                ax._bdm_reference_band = (-1., 1.)
            ax.grid(True, ls=':', alpha=.3)

        def ratio(new, old):
            return np.divide(new, old, out=np.full_like(np.asarray(new, dtype=float), np.nan),
                             where=np.asarray(old) != 0)*100-100

        def axes(title, pairs=1, samples=True):
            fig, axs = plt.subplots(2*pairs, 3, figsize=(11.4, 4.8 if pairs == 1 else 8.2),
                                    gridspec_kw={'height_ratios': [2.3, 1.15]*pairs}, squeeze=False)
            fig.subplots_adjust(left=.085, right=.99, bottom=.16 if pairs == 1 else .11,
                                top=.78 if pairs == 1 else .86, hspace=.11, wspace=.30)
            fig.suptitle(title, y=.99, fontsize=13)
            fig.legend(handles=(sample_handles if samples else [])+method_handles,
                       loc='upper center', bbox_to_anchor=(.5, .948), ncol=4, frameon=False)
            for z in (0, 1, 2):
                text = rf'$z={z}$'
                if samples:
                    counts = [len(epochs[z][f'v3_{s}_rows']) for s in selections]
                    text += '\n'+r'$N_{\rm v3}=$'+f'{counts[0]:,} / {counts[1]:,}'
                axs[0, z].set_title(text, fontsize=10, pad=8)
                for j in range(pairs):
                    axs[2*j, z].tick_params(labelbottom=False)
            return fig, axs

        def finish(fig, name, note, book):
            # Explicitly cover every plotted datum, including legacy-only
            # sparse tails. Setting a scale late can leave stale limits in
            # Matplotlib; the figure must not crop those extrema.
            for index, ax in enumerate(fig.axes):
                values = []
                for curve in ax.lines:
                    xx = np.asarray(curve.get_xdata(), dtype=float)
                    yy = np.asarray(curve.get_ydata(), dtype=float)
                    lo, hi = sorted(ax.get_xlim())
                    visible = np.isfinite(yy)
                    if len(xx) > 2:
                        visible &= (xx >= lo) & (xx <= hi)
                    values.extend(yy[visible].tolist())
                values.extend(getattr(ax, '_bdm_reference_band', ()))
                values = np.asarray(values)
                if ax.get_yscale() == 'log':
                    values = values[values > 0]
                if not len(values):
                    continue
                before = ax.get_ylim()
                transform = ax.yaxis.get_transform()
                yy = transform.transform(values)
                span = float(np.max(yy)-np.min(yy))
                pad = .065*span if span else max(abs(float(yy[0]))*.05, .1)
                lower, upper = transform.inverted().transform([np.min(yy)-pad, np.max(yy)+pad])
                if ax.get_yscale() == 'symlog':
                    lower = min(-1.5, values.min()-.1)
                ax.set_ylim(lower, upper)
                limits_checks.append(dict(page=name, axis=index, plotted_min=float(values.min()),
                    plotted_max=float(values.max()), ylim=[float(lower), float(upper)],
                    out_of_bounds_before=int(np.count_nonzero((values < before[0]) | (values > before[1]))),
                    out_of_bounds_after=int(np.count_nonzero((values < lower) | (values > upper)))))
                assert limits_checks[-1]['out_of_bounds_after'] == 0
            fig.text(.5, .018, note, ha='center', va='bottom', fontsize=8)
            book.savefig(fig, bbox_inches='tight', pad_inches=.055)
            page_names.append(name)
            plt.close(fig)

        common = ('Run F: 1024$^3$ particles, 4096$^3$ evolution mesh, '
                  '256 $h^{-1}$ Mpc; same 2048$^3$ finder field in each legacy/v3 pair.')
        radii = np.sqrt(epochs[0]['r_edges'][:-1]*epochs[0]['r_edges'][1:])
        masses = (catalogue['mass_edges'][:-1]+catalogue['mass_edges'][1:])/2
        bookpath = figures/'bdm_legacy_comparison.pdf'
        with PdfPages(bookpath) as book:
            fig, axs = axes('Halo mass function', samples=False)
            for z in (0, 1, 2):
                old = catalogue[f'F_z{z}_legacy_counts']
                new = catalogue[f'F_z{z}_v3_counts']
                for v in ('legacy', 'v3'):
                    line(axs[0, z], masses, catalogue[f'F_z{z}_{v}_hmf'], PALETTE[z], v,
                         catalogue[f'F_z{z}_{v}_counts'] > 0)
                residual(axs[1, z], masses, ratio(new, old), PALETTE[z], (old >= 30) & (new >= 30))
                axs[0, z].set_yscale('log')
                axs[0, z].set_xlim(12.45, 15.25)
                axs[1, z].set_xlim(12.45, 15.25)
                axs[1, z].set_xlabel(r'$\log_{10}(M_{\rm bound}/[h^{-1}M_\odot])$')
            axs[0, 0].set_ylabel(r'$dn/d\log_{10}M\ [(h^{-1}{\rm Mpc})^{-3}]$')
            axs[1, 0].set_ylabel('New / legacy\n$-1$ [per cent]')
            finish(fig, 'hmf', common+'\nResiduals require 30 haloes/bin in both catalogues; grey band: $\\pm1$ per cent reference, not an acceptance test.', book)

            fig, axs = axes(r'Halo--matter cross bias: $b_{hm}(k)=P_{hm}(k)/P_{mm}(k)$')
            for z, data in enumerate(epochs):
                for s, color in zip(selections, colors):
                    old = data[f'legacy_{s}_bias_k']
                    new = data[f'v3_{s}_bias_k']
                    for v in ('legacy', 'v3'):
                        vals = data[f'{v}_{s}_bias_k']
                        line(axs[0, z], vals[:, 0], vals[:, 4], color, v)
                    residual(axs[1, z], new[:, 0], ratio(new[:, 4], old[:, 4]), color)
                for ax in axs[:, z]:
                    ax.axvspan(.05, .15, color='lightsteelblue', alpha=.18, zorder=-2)
                    ax.set_xlim(.022, .201)
                axs[1, z].set_xlabel(r'$k\ [h\,{\rm Mpc}^{-1}]$')
            axs[0, 0].set_ylabel(r'$b_{hm}(k)$')
            axs[1, 0].set_ylabel('New / legacy\n$-1$ [per cent]')
            finish(fig, 'bias', common+'\nShaded $k$ interval: reported effective bias band (0.05--0.15 $h$ Mpc$^{-1}$). No independent-volume error estimate.', book)

            fig, axs = axes(r'Real-space halo autocorrelation $\xi_{hh}$')
            for z, data in enumerate(epochs):
                for s, color in zip(selections, colors):
                    old = data[f'legacy_{s}_xi']
                    new = data[f'v3_{s}_xi']
                    ok = (data[f'legacy_{s}_count'] >= 50) & (data[f'v3_{s}_count'] >= 50)
                    for v in ('legacy', 'v3'):
                        line(axs[0, z], radii, data[f'{v}_{s}_xi'], color, v)
                    residual(axs[1, z], radii, ratio(1+new, 1+old), color, ok)
                for ax in axs[:, z]:
                    ax.set_xscale('log')
                    ax.set_xlim(.01, 50)
                axs[0, z].set_yscale('symlog', linthresh=1.)
                axs[0, z].axhline(-1, color='.7', lw=.6, zorder=0)
                axs[1, z].set_xlabel(r'$r\ [h^{-1}{\rm Mpc}]$')
            axs[0, 0].set_ylabel(r'$\xi_{hh}(r)$')
            axs[1, 0].set_ylabel(r'$100\,\Delta\xi/(1+\xi_{\rm legacy})$')
            finish(fig, 'xi_hh', common+'\nPeriodic analytic RR; self-pairs excluded. $\\xi=-1$ denotes zero pairs. Residuals require 50 unordered pairs/bin in both catalogues.', book)

            fig, axs = axes(r'Mean radial pairwise peculiar velocity $v_{12}$')
            for z, data in enumerate(epochs):
                for s, color in zip(selections, colors):
                    for v in ('legacy', 'v3'):
                        count = data[f'{v}_{s}_count']
                        line(axs[0, z], radii, data[f'{v}_{s}_moments'][:, 0], color, v, count >= 50)
                    ok = (data[f'legacy_{s}_count'] >= 50) & (data[f'v3_{s}_count'] >= 50)
                    residual(axs[1, z], radii, data[f'v3_{s}_moments'][:, 0]-data[f'legacy_{s}_moments'][:, 0],
                             color, ok, percent=False)
                for ax in axs[:, z]:
                    ax.set_xscale('log')
                    ax.set_xlim(.15, 50)
                axs[1, z].set_xlabel(r'$r\ [h^{-1}{\rm Mpc}]$')
            axs[0, 0].set_ylabel(r'$v_{12}\ [{\rm km\,s}^{-1}]$')
            axs[1, 0].set_ylabel(r'$\Delta v_{12}\ [{\rm km\,s}^{-1}]$')
            finish(fig, 'velocity_mean', common+'\n$v_r=(\\mathbf v_j-\\mathbf v_i)\\cdot\\hat{\\mathbf r}_{ij}$: negative means infall; no Hubble velocity added. At least 50 pairs/bin.', book)

            for name, title, inds, labels, units in (
                ('velocity_dispersion', 'Pairwise velocity dispersions', [1, 2],
                 [r'$\sigma_r$', r'$\sigma_{t,\,1d}$'], True),
                ('velocity_higher_moments', 'Radial pairwise velocity: third and fourth moments', [3, 4],
                 [r'$\gamma_1=\mu_3/\sigma_r^3$', r'$\gamma_2=\mu_4/\sigma_r^4-3$'], False)):
                fig, axs = axes(title, pairs=2)
                for j, (index, label) in enumerate(zip(inds, labels)):
                    for z, data in enumerate(epochs):
                        for s, color in zip(selections, colors):
                            old = data[f'legacy_{s}_moments'][:, index]
                            new = data[f'v3_{s}_moments'][:, index]
                            ok = (data[f'legacy_{s}_count'] >= 50) & (data[f'v3_{s}_count'] >= 50)
                            for v in ('legacy', 'v3'):
                                line(axs[2*j, z], radii, data[f'{v}_{s}_moments'][:, index], color, v,
                                     data[f'{v}_{s}_count'] >= 50)
                            residual(axs[2*j+1, z], radii, ratio(new, old) if units else new-old,
                                     color, ok, percent=units)
                        for ax in axs[2*j:2*j+2, z]:
                            ax.set_xscale('log')
                            ax.set_xlim(.15, 50)
                        if j == 1:
                            axs[2*j+1, z].set_xlabel(r'$r\ [h^{-1}{\rm Mpc}]$')
                        else:
                            axs[2*j+1, z].tick_params(labelbottom=False)
                    axs[2*j, 0].set_ylabel(label+(r'$\ [{\rm km\,s}^{-1}]$' if units else ''))
                    axs[2*j+1, 0].set_ylabel('New / legacy\n$-1$ [per cent]' if units else r'New $-$ legacy')
                definition = (r'$\sigma_r^2=\langle(v_r-v_{12})^2\rangle$; '
                              r'$\sigma_{t,1d}^2=\langle|\Delta\mathbf v|^2-v_r^2\rangle/2$.') if units else (
                              'Central moments; excess kurtosis subtracts 3. High moments are sensitive to rare pairs.')
                finish(fig, name, common+'\n'+definition+' At least 50 pairs/bin; no independent-volume error estimate.', book)

            fig, axs = axes('Internal halo velocities at fixed bound mass', pairs=2, samples=False)
            for j, (ind, label) in enumerate([(1, r'Median $V_{\rm rms}$'), (2, r'Median $V_{\max}$')]):
                for z in (0, 1, 2):
                    old = catalogue[f'F_z{z}_legacy_profiles'][:, ind]
                    new = catalogue[f'F_z{z}_v3_profiles'][:, ind]
                    ok = (catalogue[f'F_z{z}_legacy_counts'] >= 30) & (catalogue[f'F_z{z}_v3_counts'] >= 30)
                    for v in ('legacy', 'v3'):
                        line(axs[2*j, z], masses, catalogue[f'F_z{z}_{v}_profiles'][:, ind], PALETTE[z], v,
                             catalogue[f'F_z{z}_{v}_counts'] >= 30)
                    residual(axs[2*j+1, z], masses, ratio(new, old), PALETTE[z], ok)
                    for ax in axs[2*j:2*j+2, z]:
                        ax.set_xlim(12.45, 14.75)
                    if j:
                        axs[2*j+1, z].set_xlabel(r'$\log_{10}(M_{\rm bound}/[h^{-1}M_\odot])$')
                    else:
                        axs[2*j+1, z].tick_params(labelbottom=False)
                axs[2*j, 0].set_ylabel(label+r'$\ [{\rm km\,s}^{-1}]$')
                axs[2*j+1, 0].set_ylabel('New / legacy\n$-1$ [per cent]')
            finish(fig, 'internal_velocities', common+'\nPublished internal $V_{\\rm rms}$ and $V_{\\max}$; fixed mass bins, not matched individual haloes. At least 30 haloes/bin.', book)

            fig, axs = axes('Selection check: fixed mass versus equal number density')
            # Method key differs on this page, so replace the default legend.
            for legend in fig.legends:
                legend.remove()
            fig.legend(handles=sample_handles+[
                Line2D([], [], color='.25', marker='o', ms=3, label='Fixed mass cut'),
                Line2D([], [], color='.25', ls='--', label='Equal number density')],
                loc='upper center', bbox_to_anchor=(.5, .948), ncol=4, frameon=False)
            for z, data in enumerate(epochs):
                for si, color in enumerate(colors):
                    for stem, ls, marker in [('mass', '-', 'o'), ('rank', '--', None)]:
                        s = stem+('12p5' if si == 0 else '13')
                        old = data[f'legacy_{s}_bias_k']
                        new = data[f'v3_{s}_bias_k']
                        axs[0, z].plot(new[:, 0], ratio(new[:, 4], old[:, 4]), color=color,
                                       ls=ls, marker=marker, ms=3, mec='k', mew=.35, lw=1)
                        ok = (data[f'legacy_{s}_count'] >= 50) & (data[f'v3_{s}_count'] >= 50)
                        change = ratio(1+data[f'v3_{s}_xi'], 1+data[f'legacy_{s}_xi'])
                        change[~ok] = np.nan
                        axs[1, z].plot(radii, change, color=color, ls=ls, marker=marker, ms=3, lw=1)
                    for ax in axs[:, z]:
                        ax.axhline(0, color='k', lw=.6)
                axs[0, z].tick_params(labelbottom=True)
                axs[0, z].set_xlabel(r'$k\ [h\,{\rm Mpc}^{-1}]$')
                axs[1, z].set_xscale('log')
                axs[1, z].set_xlabel(r'$r\ [h^{-1}{\rm Mpc}]$')
            fig.subplots_adjust(hspace=.75)
            axs[0, 0].set_ylabel(r'$100\,(b_{\rm new}/b_{\rm old}-1)$')
            axs[1, 0].set_ylabel(r'$100\,\Delta\xi/(1+\xi_{\rm legacy})$')
            finish(fig, 'selection_check', common+'\nEqual density: retain the top $N$ masses in each catalogue, with $N=\\min(N_{\\rm new},N_{\\rm legacy})$ above each cut; not object matching.', book)

            fig, axs = plt.subplots(1, 3, figsize=(11.4, 4.0))
            fig.subplots_adjust(left=.085, right=.99, top=.80, bottom=.25, wspace=.28)
            fig.suptitle('Supporting HMF comparison across all seven simulations', fontsize=13)
            allcolors = [PALETTE[0]] + [PALETTE[1]]*3 + [PALETTE[2]]*3
            allstyles = ['-', ':', '-', '--', '-', '--', '-.']
            for z, ax in enumerate(axs):
                for case, color, ls in zip('ABCDEFT', allcolors, allstyles):
                    old = catalogue[f'{case}_z{z}_legacy_counts']
                    new = catalogue[f'{case}_z{z}_v3_counts']
                    ok = (old >= 30) & (new >= 30)
                    change = ratio(new, old)
                    change[~ok] = np.nan
                    ax.plot(masses, change, color=color, ls=ls, marker='o', ms=3, lw=1, label=case)
                ax.axhline(0, color='k', lw=.6)
                ax.set_title(rf'$z={z}$')
                ax.set_xlabel(r'$\log_{10}(M_{\rm bound}/[h^{-1}M_\odot])$')
                ax.set_xlim(12.45, 14.75)
            axs[0].set_ylabel('HMF: new / legacy $-1$ [per cent]')
            fig.legend(*axs[0].get_legend_handles_labels(), loc='upper center', ncol=7,
                       bbox_to_anchor=(.5, .93))
            finish(fig, 'hmf_all_cases', 'Common mass cut corresponds to about 37 / 296 / 2362 particles in A / B,C,D / E,F,T. At least 30 haloes/bin.\nA: 256$^3$ particles; B/C/D: 512$^3$; E/F/T: 1024$^3$. Evolution meshes: A/C/E 2048$^3$, B 1024$^3$, D/F/T 4096$^3$; T halves the timestep.', book)

        outputs.append({'path': str(bookpath.relative_to(REPO)), 'sha256': sha(bookpath)})
        write_json(HERE/'visual-validation.json', dict(passed=True, pdf_sha256=sha(bookpath),
            axes=limits_checks, out_of_bounds_after=sum(x['out_of_bounds_after'] for x in limits_checks)))
        write_json(HERE/'plots.json', dict(completed=True, pages=len(page_names), page_names=page_names, outputs=outputs,
            script_sha256=sha(HERE/'plot.py'),
            styles={str(p.relative_to(REPO)): sha(p) for p in [STYLE/'house_style.py', STYLE/'chz-paper.mplstyle']},
            measurements={p.name: sha(p) for p in [HERE/'catalogue-statistics.npz']+
                          [HERE/f'F-z{z}-statistics.npz' for z in (0, 1, 2)]},
            cache='Temporary Matplotlib/TeX directory removed after rendering',
            figures_committed=False))
        for a in [catalogue]+epochs:
            a.close()
        print(bookpath)


if __name__ == '__main__':
    main()
