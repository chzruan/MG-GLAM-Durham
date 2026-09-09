"""Re-derive summary tables and verify retained measurement products."""
from __future__ import annotations

import json
import numpy as np

import measure as m


def rel(new, old):
    return 100*(np.asarray(new)/np.asarray(old)-1)


def main():
    h = np.load(m.HERE/'catalogue-statistics.npz')
    totals, hmf, bias, clustering, velocities, internal, numerical = [], [], [], [], [], [], []
    support, close_pairs = [], []
    verification = []
    with (m.HERE/'catalogue-statistics.json').open() as f:
        catalogue_receipt = json.load(f)
    assert catalogue_receipt['output_sha256'] == m.sha(m.HERE/'catalogue-statistics.npz')
    for input_record in catalogue_receipt['inputs']:
        assert m.sha(m.REPO/input_record['receipt']) == input_record['receipt_sha256']
        assert m.sha(m.REPO/input_record['catalogue_arrays']) == input_record['catalogue_arrays_sha256']
    for path, expected in catalogue_receipt['source_hashes'].items():
        assert m.sha(m.REPO/path) == expected
    for case in 'ABCDEFT':
        np_side = 256 if case == 'A' else 512 if case in 'BCD' else 1024
        minimum_particles = int(np.ceil(10.**12.5 / (1338888448*(1024/np_side)**3)))
        row = dict(case=case, nparticles_side=np_side, approximate_particles_at_mass_cut=minimum_particles,
                   lowest_bin_log_mass=[12.5, 12.75], epochs=[])
        for z in (0, 1, 2):
            old, new = (int(h[f'{case}_z{z}_{v}_counts'][0]) for v in m.VARIANTS)
            row['epochs'].append(dict(z=z, legacy=old, refined=new, change_percent=float(rel(new, old))))
        support.append(row)
    for z in (0, 1, 2):
        path = m.HERE/f'F-z{z}-statistics.npz'
        receipt = json.loads(path.with_suffix('.json').read_text())
        assert receipt['completed'] and m.sha(path) == receipt['output_sha256']
        for source, expected in receipt['source_hashes'].items():
            assert m.sha(m.REPO/source) == expected
        replay = json.loads((m.REPO/receipt['input']['receipt']).read_text())
        assert receipt['density']['density_sha256'] == replay['density']['sha256']
        assert receipt['density']['nparticles'] == 1024**3
        cats, _, _ = m.load_pair('F', z)
        for s, cut in [('mass12p5', 12.5), ('mass13', 13.)]:
            for v in m.VARIANTS:
                pos = cats[v][cats[v][:, 6] >= 10.**cut, :3]
                pairs = m.cKDTree(pos, boxsize=m.BOX).query_pairs(.1, output_type='ndarray')
                dr = pos[pairs[:, 1]]-pos[pairs[:, 0]]
                dr -= m.BOX*np.rint(dr/m.BOX)
                separation = np.linalg.norm(dr, axis=1)
                pairs = pairs[separation < .1]
                close_pairs.append(dict(z=z, selection=s, variant=v,
                    r_max_exclusive_mpc_h=.1, unordered_pairs=len(pairs),
                    separations_mpc_h=separation[separation < .1].tolist()))
        with np.load(path) as a:
            k = np.linalg.norm(a['wavevectors'], axis=1)*2*np.pi/m.BOX
            mm = abs(a['matter_modes'])**2
            radius = np.sqrt(a['r_edges'][1:]*a['r_edges'][:-1])
            rows_old, rows_new = (int(h[f'F_z{z}_{v}_total_rows']) for v in m.VARIANTS)
            totals.append(dict(z=z, legacy=rows_old, refined=rows_new,
                               change=rows_new-rows_old, change_percent=float(rel(rows_new, rows_old))))
            for b, lo in enumerate(h['mass_edges'][:-1]):
                old, new = (int(h[f'F_z{z}_{v}_counts'][b]) for v in m.VARIANTS)
                if min(old, new) >= 30:
                    hmf.append(dict(z=z, log_mass_lo=float(lo), log_mass_hi=float(h['mass_edges'][b+1]),
                                    legacy=old, refined=new, change_percent=float(rel(new, old))))
            for s in a['selections']:
                small_bias_changes, all_k_changes = [], []
                for v in m.VARIANTS:
                    key = f'{v}_{s}'
                    n = len(a[key+'_rows'])
                    expected_xi = (2*a[key+'_count']*m.BOX**3 /
                                   (n*(n-1)*4*np.pi/3*np.diff(a['r_edges']**3))) - 1
                    np.testing.assert_allclose(a[key+'_xi'], expected_xi, rtol=1e-10, atol=1e-10)
                    hm = np.real(a[key+'_halo_modes']*a['matter_modes'].conj())
                    use = (k >= .05) & (k < .15)
                    np.testing.assert_allclose(a[key+'_bias'][0], hm[use].sum()/mm[use].sum(), atol=1e-13)
                    small_bias_changes.append(float(rel(a[key+'_bias_ng256'][0], a[key+'_bias'][0])))
                    all_k_changes.extend(rel(a[key+'_bias_k_ng256'][:, -1], a[key+'_bias_k'][:, -1]))
                oldb, newb = (float(a[f'{v}_{s}_bias'][0]) for v in m.VARIANTS)
                delta = float(rel(newb, oldb))
                delta_low = float(rel(a[f'v3_{s}_bias_ng256'][0], a[f'legacy_{s}_bias_ng256'][0]))
                bias.append(dict(z=z, selection=str(s),
                    n_legacy=len(a[f'legacy_{s}_rows']), n_refined=len(a[f'v3_{s}_rows']),
                    minimum_mass_legacy=float(a[f'legacy_{s}_minimum_mass']),
                    minimum_mass_refined=float(a[f'v3_{s}_minimum_mass']),
                    bias_legacy=oldb, bias_refined=newb, change_percent=delta,
                    change_percent_ng256=delta_low,
                    alternative_band_changes_percent=rel(a[f'v3_{s}_bias'][1:], a[f'legacy_{s}_bias'][1:]).tolist()))
                numerical.append(dict(z=z, selection=str(s),
                    max_absolute_effective_bias_grid_change_percent=max(map(abs, small_bias_changes)),
                    paired_change_grid_difference_pp=abs(delta-delta_low),
                    max_absolute_kbin_bias_grid_change_percent=float(np.max(np.abs(all_k_changes)))))
                for rlo, rhi in ((.5, 1.), (1., 5.), (5., 30.)):
                    use = ((radius >= rlo) & (radius < rhi) &
                           (a[f'legacy_{s}_count'] >= 50) & (a[f'v3_{s}_count'] >= 50))
                    oldxi, newxi = a[f'legacy_{s}_xi'][use], a[f'v3_{s}_xi'][use]
                    if not use.any():
                        continue
                    dx = rel(1+newxi, 1+oldxi)
                    meaningful = np.abs(oldxi) > .1
                    clustering.append(dict(z=z, selection=str(s), rlo=rlo, rhi=rhi, bins=int(use.sum()),
                        max_absolute_pair_probability_change_percent=float(np.max(np.abs(dx))),
                        median_pair_probability_change_percent=float(np.median(dx)),
                        max_absolute_xi_change_percent_where_abs_legacy_xi_gt_point1=(
                            float(np.max(np.abs(rel(newxi[meaningful], oldxi[meaningful])))) if meaningful.any() else None),
                        min_unordered_pairs=int(min(a[f'legacy_{s}_count'][use].min(), a[f'v3_{s}_count'][use].min()))))
                    oldmom, newmom = a[f'legacy_{s}_moments'][use], a[f'v3_{s}_moments'][use]
                    velocities.append(dict(z=z, selection=str(s), rlo=rlo, rhi=rhi, bins=int(use.sum()),
                        max_absolute_mean_change_km_s=float(np.max(abs(newmom[:, 0]-oldmom[:, 0]))),
                        max_absolute_radial_sigma_change_percent=float(np.max(abs(rel(newmom[:, 1], oldmom[:, 1])))),
                        max_absolute_tangential_sigma_change_percent=float(np.max(abs(rel(newmom[:, 2], oldmom[:, 2])))),
                        max_absolute_skewness_change=float(np.max(abs(newmom[:, 3]-oldmom[:, 3]))),
                        max_absolute_excess_kurtosis_change=float(np.max(abs(newmom[:, 4]-oldmom[:, 4])))))
            use = (h[f'F_z{z}_legacy_counts'] >= 30) & (h[f'F_z{z}_v3_counts'] >= 30)
            for index, name in [(1, 'internal_vrms'), (2, 'vmax')]:
                delta = rel(h[f'F_z{z}_v3_profiles'][use, index], h[f'F_z{z}_legacy_profiles'][use, index])
                internal.append(dict(z=z, property=name, eligible_mass_bins=int(use.sum()),
                                     max_absolute_median_change_percent=float(np.max(abs(delta)))))
            verification.append(dict(z=z, output_sha256=receipt['output_sha256'],
                                     normalisation_and_bias_rederived=True,
                                     input_receipt_and_density_hash_binding=True))
    summary = dict(completed=True, selection_scope='F, bound-mass cuts and rank controls; z=0,1,2',
        total_rows=totals, hmf=hmf, bias=bias, clustering=clustering, velocity_moments=velocities,
        internal_profiles=internal, numerical_grid_checks=numerical, verification=verification,
        supporting_lowest_mass_bin=support, close_pairs=close_pairs,
        report_script_sha256=m.sha(m.HERE/'report.py'),
        uncertainties='Paired single volume; no independent-volume uncertainties or significance test')
    m.write_json(m.HERE/'results.json', summary)
    lines = ['# Measured legacy versus refined BDM differences', '',
        'Main sample: F (1024³ particles, 4096³ evolution mesh, L=256 Mpc/h), '
        'z=0, 1, 2. Both finders act on the same 2048³ density field. '
        'The [nine-page PDF](figures/bdm_legacy_comparison.pdf) contains all four requested '
        'statistics and supporting selection/internal-velocity checks.', '',
        '**These plots quantify the effect of the repairs; they do not demonstrate a '
        'new convergence advantage or supply independent-volume uncertainties.** '
        'Full definitions and display cuts are in [README.md](README.md).', '',
        '## Catalogue abundance', '',
        '| z | Legacy rows | Refined rows | Difference | Change (%) |',
        '|---|---:|---:|---:|---:|']
    for r in totals:
        lines.append(f"| {r['z']} | {r['legacy']:,} | {r['refined']:,} | {r['change']:+,} | {r['change_percent']:+.3f} |")
    lines += ['', 'Total rows use the native publication limit (2.5e12 Msun/h); the plots and '
              'clustering samples start at 10^12.5 Msun/h. Changes in the HMF are mass-dependent:', '',
              '| z | Mass bin, log10(Mbound/[Msun/h]) | Legacy | Refined | HMF change (%) |',
              '|---|---|---:|---:|---:|']
    for r in hmf:
        lines.append(f"| {r['z']} | {r['log_mass_lo']:.2f}–{r['log_mass_hi']:.2f} | {r['legacy']} | {r['refined']} | {r['change_percent']:+.3f} |")
    lines += ['', 'Only bins with at least 30 haloes in both catalogues appear above. '
              'Sparse-bin fractional changes should be read alongside their counts.', '',
              'The common mass cut does not represent the same particle count in the supporting simulations. '
              'At 10^12.5 Msun/h, A has about 37 particles per halo, B/C/D about 296 and E/F/T about 2362. '
              'The lowest-bin (12.50–12.75) HMF changes are:', '',
              '| Case | Particle floor at common mass cut | z=0 change (%) | z=1 change (%) | z=2 change (%) |',
              '|---|---:|---:|---:|---:|']
    for r in support:
        changes = ' | '.join(f"{x['change_percent']:+.3f}" for x in r['epochs'])
        lines.append(f"| {r['case']} | {r['approximate_particles_at_mass_cut']} | {changes} |")
    lines += ['', 'The much larger A response is in a poorly resolved mass range. '
              'The small F differences should not be generalised to this sample.', '',
              '## Matter-referenced effective bias', '',
              'The band estimator is sum P_hm / sum P_mm over **0.05 ≤ k < 0.15 h/Mpc**. '
              'This is a finite-band cross bias, not a fitted k→0 limit.', '',
              '| z | Selection | N legacy / refined | b legacy | b refined | Change (%) |',
              '|---|---|---:|---:|---:|---:|']
    for r in bias:
        label = {'mass12p5':'M ≥ 10^12.5','mass13':'M ≥ 10^13',
                 'rank12p5':'equal n, 12.5 cut','rank13':'equal n, 13 cut'}[r['selection']]
        lines.append(f"| {r['z']} | {label} | {r['n_legacy']:,} / {r['n_refined']:,} | {r['bias_legacy']:.5f} | {r['bias_refined']:.5f} | {r['change_percent']:+.4f} |")
    lines += ['', f"Changing the matter measurement grid from 512³ to 256³ changes the effective bias by at most "
              f"{max(r['max_absolute_effective_bias_grid_change_percent'] for r in numerical):.5f}% and the paired bias "
              f"change by at most {max(r['paired_change_grid_difference_pp'] for r in numerical):.5f} percentage points. "
              'Direct halo Fourier sums use no halo assignment grid. This validates measurement coarsening, not finder-mesh independence.', '',
              '## Clustering and pairwise velocities', '',
              'The following are **maximum absolute per-bin changes at 5–30 Mpc/h**, '
              'using geometric bin centres and at least 50 unordered pairs/bin in both catalogues. '
              'They are descriptive extrema, not error bars. The correlation column is the '
              'change in **1+xi**, not the percentage change in xi.', '',
              '| z | log10 mass cut | Δ(1+xi)/(1+xi) (%) | Δv12 (km/s) | Δsigma_r (%) | Δsigma_t (%) | Δskewness | Δexcess kurtosis |',
              '|---|---|---:|---:|---:|---:|---:|---:|']
    for r in velocities:
        if r['rlo'] != 5 or not r['selection'].startswith('mass'):
            continue
        c = next(x for x in clustering if x['z']==r['z'] and x['selection']==r['selection'] and x['rlo']==5)
        cut = '12.5' if r['selection']=='mass12p5' else '13.0'
        lines.append(f"| {r['z']} | {cut} | {c['max_absolute_pair_probability_change_percent']:.3f} | "
            f"{r['max_absolute_mean_change_km_s']:.3f} | {r['max_absolute_radial_sigma_change_percent']:.3f} | "
            f"{r['max_absolute_tangential_sigma_change_percent']:.3f} | {r['max_absolute_skewness_change']:.4f} | "
            f"{r['max_absolute_excess_kurtosis_change']:.4f} |")
    lines += ['', 'The plots retain the radial dependence, smaller-scale changes and gaps from the pair-count floor. '
              '`results.json` also reports 0.5–1 and 1–5 Mpc/h extrema, both equal-density controls, '
              'and ordinary Δxi/xi where |xi_legacy| > 0.1. All unmasked bins, raw pair counts, '
              'moments and low-k Fourier coefficients are in the three `F-z*-statistics.npz` files.', '',
              '**Below 1 Mpc/h the velocity changes are substantially larger.** '
              'For the 10^12.5 Msun/h cut, the 0.5–1 Mpc/h bins meeting the same 50-pair floor give:', '',
              '| z | Bins | Minimum pairs/catalogue/bin | Max abs Δv12 (km/s) | Max abs Δsigma_r (%) | Max abs Δsigma_t (%) |',
              '|---|---:|---:|---:|---:|---:|']
    for r in velocities:
        if r['rlo'] != .5 or r['selection'] != 'mass12p5':
            continue
        c = next(x for x in clustering if x['z']==r['z'] and x['selection']==r['selection'] and x['rlo']==.5)
        lines.append(f"| {r['z']} | {r['bins']} | {c['min_unordered_pairs']} | {r['max_absolute_mean_change_km_s']:.3f} | "
                     f"{r['max_absolute_radial_sigma_change_percent']:.3f} | {r['max_absolute_tangential_sigma_change_percent']:.3f} |")
    lines += ['',
              '### Very close halo pairs', '',
              'The spikes at the smallest separations in the legacy xi curve are produced by just a few pairs. '
              'The table uses exact periodic centre distances below 0.1 Mpc/h. '
              'A plotted bin centred below 0.1 can include pairs above 0.1; these counts are independent of the plotting bins:', '',
              '| z | Mass cut | Legacy pairs with r < 0.1 Mpc/h | Refined pairs |',
              '|---|---|---:|---:|']
    for z in (0, 1, 2):
        for s, label in [('mass12p5', '10^12.5'), ('mass13', '10^13')]:
            old, new = [next(x['unordered_pairs'] for x in close_pairs if x['z']==z and
                            x['selection']==s and x['variant']==v) for v in m.VARIANTS]
            lines.append(f'| {z} | {label} | {old} | {new} |')
    lines += ['', 'The refinement removes these very close pairs and produces a clear exclusion region. '
              'Centre distances alone do not establish that a pair had exactly identical particle membership; '
              'the legacy memberships were not retained. The spikes are not a well-sampled clustering signal.', '',
              '## Internal halo velocities', '',
              '| z | Profile | Maximum absolute change of bin median (%) |',
              '|---|---|---:|']
    for r in internal:
        lines.append(f"| {r['z']} | {r['property']} | {r['max_absolute_median_change_percent']:.3f} |")
    lines += ['', 'These use fixed mass bins with 30 haloes per catalogue; they do not match individual haloes.', '',
              '## Validation and provenance', '',
              '- All 21 catalogue-pair files were checked against the frozen replay hashes.',
              '- All three 32-GiB matter tapes were hashed during streaming and matched their replay receipts.',
              '- Seven independent small controls passed; all 24 measured selections have exact agreement '
              'between Corrfunc ordered and KDTree unordered pair counts.',
              '- The summary independently reconstructs the periodic RR normalisation and band bias from '
              'retained counts/Fourier modes; input, output and source hashes are checked.',
              '- The single shared Slurm pilot was cancelled while pending, with zero allocated/billed time. '
              'Sequential low-priority login-node execution used the user’s idle-node authorization; '
              'actual time/memory are in `execution.json`.',
              '- No production source, simulation input or frozen convergence result was changed.', '']
    (m.HERE/'RESULTS.md').write_text('\n'.join(lines))
    h.close()


if __name__ == '__main__':
    main()
