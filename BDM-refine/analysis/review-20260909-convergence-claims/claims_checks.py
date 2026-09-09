"""Independent challenge to the BDM convergence claims, and audit of REVIEW-RESPONSE.md.

micromamba run -n cosemu python3 -B claims_checks.py <stage>
Stages: verify | additivity | finder | selection | tails | legacy | response
Read-only with respect to the campaign; results go to claims_results.json.
"""
import argparse, hashlib, json, math, os, struct
from pathlib import Path

os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('MKL_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np

REPO = Path('/cosma8/data/dp203/dc-ruan1/mgglam_claude/MG-GLAM')
C = REPO / 'BDM-refine/validation/20260908-convergence'
HERE = Path(__file__).resolve().parent
RESULTS = HERE / 'claims_results.json'
EDGES = np.arange(12.25, 15.76, .25)
PUB_MIN = 2.5e12
PAIRS = [('A','C','particle'),('C','E','particle'),('A','E','particle'),
         ('B','C','force'),('C','D','force'),('E','F','force'),('F','T','time')]


def sha(p, chunk=8*1024**2):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        while b := f.read(chunk): h.update(b)
    return h.hexdigest()


def store(key, value):
    d = json.loads(RESULTS.read_text()) if RESULTS.exists() else {}
    d[key] = value
    RESULTS.write_text(json.dumps(d, indent=1, default=str) + '\n')
    print(f'[{key}] stored')


def load():
    conv = json.loads((C/'convergence.json').read_text())
    asmt = json.loads((C/'convergence-assessment.json').read_text())
    return conv, asmt, {(r['coarse'],r['reference'],r['redshift'],r['particle_floor']): r
                        for r in conv['comparisons']}


def index_of(name, z):
    rep = json.loads((C/f'v3-validation-{name}-z{z}-ng2048.json').read_text())
    ev = rep['variants']['v3']['membership']
    raw = REPO / ev['retained_raw']
    with np.load(raw.with_suffix('.index.npz')) as d:
        a = {k: d[k] for k in d.files}
    return dict(raw=raw, mass_one=ev['mass_one'], box=rep['spec']['box_mpc_h'],
                nrow=rep['spec']['nrow'], **a)


# --------------------------------------------------------------------------- #
def stage_verify():
    """Re-derive the 63 decisions from the NEW assessment; audit the response's numbers."""
    conv, asmt, rows = load()
    out = {}
    ok = lambda v: v is not None and np.isfinite(v)

    # (a) decisions unchanged, re-derived with independent code
    mism = []
    for row in conv['comparisons']:
        k = (row['coarse'],row['reference'],row['redshift'],row['particle_floor'])
        p = [x for x in asmt['comparisons']
             if (x['coarse'],x['reference'],x['redshift'],x['particle_floor']) == k][0]
        for i, valid in enumerate(row['valid_mass_bins']):
            m = row['matched_statistics']['bound_mass'][i]; v = row['matched_statistics']['vmax'][i]
            n = [row['left_counts'][i], row['right_counts'][i],
                 row['reference_eligible_counts'][i], m['count'], v['count']]
            r_, s_ = row['abundance_ratio'][i], row['abundance_ratio_jackknife8_sigma'][i]
            mm, mv, cp = m['q16_median_q84'][1], v['q16_median_q84'][1], row['reference_completeness'][i]
            elig = bool(valid) and min(n) >= 30 and all(ok(x) for x in (r_,s_,mm,mv,cp))
            ch = dict(abundance=abs(100*(r_-1)) <= 5., abundance_uncertainty=100*s_ <= 5.,
                      bound_mass=abs(mm) <= 5., resolved_vmax=abs(mv) <= 2.,
                      completeness=cp >= .9) if elig else {}
            b = p['bins'][i]
            if elig != b['eligible'] or ch != b['checks'] or bool(elig and all(ch.values())) != b['meets_working_criteria']:
                mism.append([str(k), i])
    out['decisions_reproduced'] = dict(n=len(conv['comparisons']), mismatches=mism,
                                       criteria=asmt['criteria'])
    out['assess_source_hash_now_recorded'] = (asmt.get('assess_source_sha256') == sha(C/'assess.py'))
    out['input_sha_still_bound'] = (asmt['input_sha256'] == sha(C/'convergence.json'))

    # (b) shape table: recompute maxima inside passing intervals for all 21 floor-300 rows
    shape = {}
    for l, r, _ in PAIRS:
        for z in (2,1,0):
            row = rows[(l,r,z,300)]
            p = [x for x in asmt['comparisons']
                 if (x['coarse'],x['reference'],x['redshift'],x['particle_floor']) == (l,r,z,300)][0]
            idx = [i for i,b in enumerate(p['bins']) if b['meets_working_criteria']]
            def mx(name):
                vals = [abs(row['matched_statistics'][name][i]['q16_median_q84'][1]) for i in idx
                        if row['matched_statistics'][name][i]['count'] >= 30
                        and ok(row['matched_statistics'][name][i]['q16_median_q84'][1])]
                return max(vals) if vals else None
            shape[f'{l}/{r} z{z}'] = dict(ba=mx('axis_ba'), ca=mx('axis_ca'))
    published = {}
    for d in asmt['supplemental_diagnostics']:
        if d['nominal_particle_floor'] != 300: continue
        published[f"{d['coarse']}/{d['reference']} z{d['redshift']}"] = dict(
            ba=d['unscreened_shapes_inside_passing_intervals']['axis_ba']['maximum_absolute_bin_median_shift_percent'],
            ca=d['unscreened_shapes_inside_passing_intervals']['axis_ca']['maximum_absolute_bin_median_shift_percent'])
    diff = {k: [shape[k], published[k]] for k in shape
            if not all((a is None and b is None) or (a is not None and b is not None and abs(a-b) < 1e-9)
                       for a, b in zip(shape[k].values(), published[k].values()))}
    out['shape_table'] = dict(mine=shape, disagreements=diff,
                              max_ca_overall=max((v['ca'] for v in shape.values() if v['ca']), default=None))

    # (c) my Prompt-1 "21 bins" claim: how many actually pass the FULL screen?
    full, abundance_only = [], []
    for l, r, _ in PAIRS:
        for z in (2,1,0):
            row = rows[(l,r,z,300)]
            p = [x for x in asmt['comparisons']
                 if (x['coarse'],x['reference'],x['redshift'],x['particle_floor']) == (l,r,z,300)][0]
            for i, b in enumerate(p['bins']):
                if not b['eligible'] or not b['checks'].get('abundance'): continue
                s_ = row['abundance_ratio_jackknife8_sigma'][i]
                if 100*s_ <= 2.5: continue
                rec = [f'{l}/{r}', z, round(float(EDGES[i]),2), round(100*(row['abundance_ratio'][i]-1),2),
                       round(100*s_,2), b['meets_working_criteria']]
                abundance_only.append(rec)
                if b['meets_working_criteria']: full.append(rec)
    out['sigma_over_2p5_percent'] = dict(
        meet_abundance_difference_only=len(abundance_only),
        also_pass_full_screen=len(full),
        passing_bins=full,
        note='Prompt-1 REVIEW.md cited 21 bins as "pass the 5% abundance criterion"; only the '
             'subset listed here passes the whole screen. Its B/C z=2 example fails the sigma condition.')

    # (d) effective particle cuts table
    cuts = {}
    for l, r, _ in PAIRS:
        row = rows[(l,r,0,300)]
        ml = conv['abundances'][f'{l}/z0']['mass_one']; mr = conv['abundances'][f'{r}/z0']['mass_one']
        first = next((EDGES[i] for i,v in enumerate(row['valid_mass_bins']) if v), None)
        cuts[f'{l}/{r}'] = dict(mass_floor=row['mass_floor_msun_h'],
            min_particles_left=math.ceil(row['mass_floor_msun_h']/ml),
            min_particles_reference=math.ceil(row['mass_floor_msun_h']/mr),
            first_whole_bin=float(first) if first is not None else None,
            min_particles_first_bin_left=math.ceil(10**first/ml) if first is not None else None,
            min_particles_first_bin_reference=math.ceil(10**first/mr) if first is not None else None)
    out['effective_cuts'] = cuts

    # (e) legacy receipts: is a membership tape present anywhere?
    legacy = []
    for z in (2,1,0):
        for n in 'ABCDEFT':
            rec = json.loads((C/f'replay-{n}-z{z}-ng2048.json').read_text())
            legacy.append([f'{n}/z{z}', sorted(rec['variants']['legacy'].keys()),
                           rec['variants']['legacy']['catalogue']['rows'],
                           rec['variants']['v3']['catalogue']['rows']])
    out['legacy_variant_contents'] = dict(
        any_with_membership=[r[0] for r in legacy if 'membership' in r[1]],
        keys_example=legacy[0][1],
        total_legacy_rows=sum(r[2] for r in legacy), total_v3_rows=sum(r[3] for r in legacy),
        per_catalogue=legacy)
    # do the retained pair arrays hold legacy catalogue columns?
    arrays = C/'work/replays/E-z0-ng2048/catalogues-pair.npz'
    if arrays.is_file():
        with np.load(arrays) as d:
            out['legacy_variant_contents']['pair_npz_keys'] = list(d.files)
            out['legacy_variant_contents']['pair_npz_shapes'] = {k: list(d[k].shape) for k in d.files}
    store('verification_and_response_audit', out)


# --------------------------------------------------------------------------- #
def _hist(mass, sel=None):
    return np.histogram(np.log10(mass if sel is None else mass[sel]), EDGES)[0]


def _octant(pos, box):
    o = np.floor((pos % box) / (box/2)).astype(int)
    return o[:,0] + 2*o[:,1] + 4*o[:,2]


def _jk(cl, ol, cr, orr, valid):
    with np.errstate(invalid='ignore'):
        rr = np.divide(cl-ol, cr-orr, out=np.full((8,14), np.nan), where=(cr-orr) > 0)
    s = np.sqrt(7/8*np.sum((rr-np.mean(rr, axis=0))**2, axis=0))
    s[~valid] = np.nan
    return s


def stage_additivity():
    """Are the three refinements separable? Test with the chains already measured."""
    conv, asmt, rows = load()
    out = {}
    ab = conv['abundances']

    # (1) abundance chain closure is exact by construction (ratios of the same counts)
    #     -> the informative test is on MATCHED median shifts, which use different samples.
    chain = []
    for z in (2,1,0):
        for i in range(14):
            ac = rows[('A','C',z,300)]; ce = rows[('C','E',z,300)]; ae = rows[('A','E',z,300)]
            if not (ac['valid_mass_bins'][i] and ce['valid_mass_bins'][i] and ae['valid_mass_bins'][i]):
                continue
            g = lambda r, k: r['matched_statistics'][k][i]['q16_median_q84'][1]
            n = lambda r, k: r['matched_statistics'][k][i]['count']
            if min(n(ac,'bound_mass'), n(ce,'bound_mass'), n(ae,'bound_mass')) < 30: continue
            for key in ('bound_mass','vmax'):
                a, b, d = g(ac,key), g(ce,key), g(ae,key)
                if None in (a,b,d): continue
                composed = 100*((1+a/100)*(1+b/100) - 1)
                chain.append(dict(z=z, bin=round(float(EDGES[i]),2), prop=key,
                                  d_AC=round(a,3), d_CE=round(b,3), composed_AE=round(composed,3),
                                  measured_AE=round(d,3), residual=round(d-composed,3)))
    res = [c['residual'] for c in chain]
    out['particle_chain_A_C_E'] = dict(
        n=len(chain), max_abs_residual_percent=round(max(map(abs,res)), 3) if res else None,
        median_abs_residual_percent=round(float(np.median(np.abs(res))), 3) if res else None,
        rows=chain,
        meaning='If particle refinement composes, measured A/E should equal A/C composed with C/E.')

    # (2) is the force-mesh effect independent of particle load?  C/D (512^3) vs E/F (1024^3)
    force = []
    for z in (2,1,0):
        cd = rows[('C','D',z,300)]; ef = rows[('E','F',z,300)]
        for i in range(14):
            if not (cd['valid_mass_bins'][i] and ef['valid_mass_bins'][i]): continue
            g = lambda r,k: r['matched_statistics'][k][i]['q16_median_q84'][1]
            n = lambda r,k: r['matched_statistics'][k][i]['count']
            row = dict(z=z, bin=round(float(EDGES[i]),2))
            keep = False
            for key in ('bound_mass','vmax','axis_ca'):
                a, b = g(cd,key), g(ef,key)
                if a is None or b is None or min(n(cd,key), n(ef,key)) < 30: continue
                row[f'CD_{key}'] = round(a,3); row[f'EF_{key}'] = round(b,3)
                row[f'diff_{key}'] = round(b-a,3); keep = True
            ra, rb = cd['abundance_ratio'][i], ef['abundance_ratio'][i]
            if ra is not None and rb is not None:
                row['CD_abundance'] = round(100*(ra-1),3); row['EF_abundance'] = round(100*(rb-1),3)
                row['diff_abundance'] = round(100*(rb-ra),3); keep = True
            if keep: force.append(row)
    out['force_effect_at_two_particle_loads'] = dict(
        n=len(force), rows=force,
        meaning='2048->4096 force refinement measured at 512^3 (C/D) and at 1024^3 (E/F) particles. '
                'Agreement supports treating the force term as separable from particle count.')
    store('separability', out)


def stage_finder():
    """Does the fixed 2048^3 finder mesh limit what is detected in every run?"""
    conv, _, _ = load()
    rowsout = []
    for z in (2,1,0):
        for n in 'ABCDEFT':
            rep = json.loads((C/f'v3-validation-{n}-z{z}-ng2048.json').read_text())
            ev = rep['variants']['v3']['membership']
            idx = index_of(n, z)
            mass = idx['counts'].astype(np.float64)*idx['mass_one']
            rowsout.append(dict(case=n, z=z, evolution_ngrid=rep['spec']['evolution_ngrid'],
                nrow=rep['spec']['nrow'], candidates=ev['candidates'], published=ev['selected'],
                published_over_candidates=round(ev['selected']/ev['candidates'], 5),
                particles=rep['spec']['nrow']**3,
                candidates_per_particle=ev['candidates']/rep['spec']['nrow']**3,
                min_published_mass=float(mass.min()), max_published_mass=float(mass.max()),
                median_aperture_mpc_h=float(np.median(idx['properties'][:,8])),
                min_aperture_mpc_h=float(idx['properties'][:,8].min()),
                aperture_in_analysis_cells=float(np.median(idx['properties'][:,8])/(256./2048))))
    store('finder_mesh_indicators', dict(
        analysis_cell_mpc_h=256./2048, rows=rowsout,
        meaning='Candidate count is set by peak detection on the shared 2048^3 analysis mesh. '
                'If it were mesh-limited it would be insensitive to the particle load.'))


def stage_selection():
    """Characterise the reference haloes that matching loses, at z=0."""
    conv, asmt, rows = load()
    out = {}
    for l, r in [('C','E'), ('E','F'), ('F','T'), ('A','E')]:
        L, R = index_of(l, 0), index_of(r, 0)
        with np.load(C/f'work/analysis/matches-{l}-{r}-z0.npz') as d:
            P = d['pairs']
        ii, jj = P[:,0].astype(int), P[:,1].astype(int)
        mL = L['counts'].astype(float)*L['mass_one']; mR = R['counts'].astype(float)*R['mass_one']
        row = rows[(l,r,0,300)]
        mf = row['mass_floor_msun_h']
        elig = mR >= mf
        matched = np.zeros(len(mR), bool); matched[jj] = True
        # restrict to bins the screen actually uses
        lo = next((EDGES[i] for i,v in enumerate(row['valid_mass_bins']) if v), None)
        use = elig & (np.log10(mR) >= lo)
        miss = use & ~matched
        got = use & matched
        # what are the lost ones like?
        out[f'{l}/{r}'] = dict(
            reference_eligible=int(use.sum()), unmatched=int(miss.sum()),
            unmatched_fraction=round(float(miss.sum()/max(use.sum(),1)), 5),
            median_log10M_matched=round(float(np.median(np.log10(mR[got]))), 4),
            median_log10M_unmatched=round(float(np.median(np.log10(mR[miss]))), 4) if miss.sum() else None,
            median_aperture_matched=round(float(np.median(R['properties'][got,8])), 4),
            median_aperture_unmatched=round(float(np.median(R['properties'][miss,8])), 4) if miss.sum() else None,
            median_vmax_matched=round(float(np.median(R['properties'][got,11])), 3),
            median_vmax_unmatched=round(float(np.median(R['properties'][miss,11])), 3) if miss.sum() else None,
            unmatched_in_lowest_used_bin=int((miss & (np.log10(mR) < lo+0.25)).sum()),
            eligible_in_lowest_used_bin=int((use & (np.log10(mR) < lo+0.25)).sum()),
            # what would the abundance ratio be if the unmatched reference haloes were dropped?
            hypothetical_bias_if_unmatched_excluded_percent=round(
                100*(1/(1-float(miss.sum()/max(use.sum(),1))) - 1), 4))
    store('matching_selection', dict(rows=out,
        meaning='Matched medians are computed on matched objects only. This quantifies how many '
                'eligible reference haloes are lost and whether they differ systematically.'))


def stage_tails():
    """Median shift versus halo-to-halo scatter inside the passing intervals."""
    conv, asmt, rows = load()
    out = []
    for l, r, _ in PAIRS:
        for z in (2,1,0):
            p = [x for x in asmt['comparisons']
                 if (x['coarse'],x['reference'],x['redshift'],x['particle_floor']) == (l,r,z,300)][0]
            row = rows[(l,r,z,300)]
            idx = [i for i,b in enumerate(p['bins']) if b['meets_working_criteria']]
            if not idx: continue
            for key in ('bound_mass','vmax'):
                w, m = [], []
                for i in idx:
                    q = row['matched_statistics'][key][i]['q16_median_q84']
                    if None in q: continue
                    w.append((q[2]-q[0])/2); m.append(abs(q[1]))
                if not w: continue
                out.append(dict(pair=f'{l}/{r}', z=z, prop=key,
                    max_abs_median_percent=round(max(m),3),
                    max_half_1684_width_percent=round(max(w),3),
                    min_half_1684_width_percent=round(min(w),3),
                    scatter_over_median_at_widest=round(max(w)/max(max(m),1e-9),1)))
    store('tails_vs_medians', dict(rows=out,
        meaning='The screen constrains bin medians. The 16-84 half-width is the halo-to-halo scatter '
                'that the screen never constrains.'))


# --------------------------------------------------------------------------- #
def stage_legacy():
    """Abundance convergence of the LEGACY catalogues, from retained pair arrays only.

    No new simulation, replay or matching. Legacy receipts carry no membership tape,
    so only catalogue-level (abundance / distribution) statistics are available.
    """
    conv, asmt, rows = load()
    cat = {}
    for z in (2,1,0):
        for n in 'ABCDEFT':
            rec = json.loads((C/f'replay-{n}-z{z}-ng2048.json').read_text())
            arr = Path(rec['catalogue_arrays'])
            assert sha(arr) == rec['catalogue_arrays_sha256'], arr
            with np.load(arr) as d:
                cat[n,z] = {k: d[k] for k in ('v3','legacy')}
    box = 256.
    out = {'analysis': 'Bound mass is catalogue column 6 (float32, 4 significant figures in ASCII); '
                       'bins are 0.25 dex so edge rounding is immaterial. Octant jackknife as in analyze.py.'}
    # sanity: v3 abundance from the catalogue must reproduce the membership-derived counts
    check = []
    for z in (2,1,0):
        for n in 'ABCDEFT':
            a = cat[n,z]['v3']
            c1 = _hist(a[:,6])
            c2 = np.array(conv['abundances'][f'{n}/z{z}']['counts'])
            check.append([f'{n}/z{z}', int(np.abs(c1-c2).sum()), int(c2.sum()), len(a)])
    out['v3_catalogue_vs_membership_counts'] = dict(
        rows=check, total_bin_disagreement=int(sum(r[1] for r in check)))
    table = []
    for l, r, kind in PAIRS:
        for z in (2,1,0):
            row = rows[(l,r,z,300)]
            valid = np.array(row['valid_mass_bins'])
            entry = dict(pair=f'{l}/{r}', kind=kind, z=z,
                         mass_floor=row['mass_floor_msun_h'], bins=[])
            for variant in ('v3','legacy'):
                A, B = cat[l,z][variant], cat[r,z][variant]
                cl, cr = _hist(A[:,6]), _hist(B[:,6])
                ol = np.array([_hist(A[:,6], _octant(A[:,:3], box) == q) for q in range(8)])
                orr = np.array([_hist(B[:,6], _octant(B[:,:3], box) == q) for q in range(8)])
                ratio = np.divide(cl, cr, out=np.full(14, np.nan), where=cr > 0)
                ratio[~valid] = np.nan
                entry[variant] = dict(rows=int(len(A)), reference_rows=int(len(B)),
                    shift_percent=[None if not np.isfinite(v) else round(100*(v-1),3) for v in ratio],
                    sigma_percent=[None if not np.isfinite(v) else round(100*v,3)
                                   for v in _jk(cl, ol, cr, orr, valid)],
                    left_counts=cl.tolist(), right_counts=cr.tolist())
            # per-bin comparison where both are measurable and counts suffice
            for i in range(14):
                v3s, lgs = entry['v3']['shift_percent'][i], entry['legacy']['shift_percent'][i]
                if v3s is None or lgs is None: continue
                if min(entry['v3']['left_counts'][i], entry['v3']['right_counts'][i],
                       entry['legacy']['left_counts'][i], entry['legacy']['right_counts'][i]) < 30: continue
                entry['bins'].append(dict(log10_mass=round(float(EDGES[i]),2),
                    v3_shift=v3s, legacy_shift=lgs, legacy_minus_v3=round(lgs-v3s,3),
                    v3_abs_smaller=bool(abs(v3s) < abs(lgs)),
                    v3_sigma=entry['v3']['sigma_percent'][i], legacy_sigma=entry['legacy']['sigma_percent'][i]))
            table.append(entry)
    allbins = [b for e in table for b in e['bins']]
    wins = sum(b['v3_abs_smaller'] for b in allbins)
    out['comparison'] = table
    out['summary'] = dict(
        comparable_bins=len(allbins), v3_smaller_absolute_shift=wins,
        legacy_smaller_absolute_shift=len(allbins)-wins,
        median_abs_v3=round(float(np.median([abs(b['v3_shift']) for b in allbins])),3),
        median_abs_legacy=round(float(np.median([abs(b['legacy_shift']) for b in allbins])),3),
        mean_abs_v3=round(float(np.mean([abs(b['v3_shift']) for b in allbins])),3),
        mean_abs_legacy=round(float(np.mean([abs(b['legacy_shift']) for b in allbins])),3),
        max_abs_legacy_minus_v3=round(max(abs(b['legacy_minus_v3']) for b in allbins),3),
        caveat='Abundance only. Legacy has no membership tape, so matched-property, completeness and '
               'shared-tracer statistics cannot be produced without new legacy replays.')
    store('legacy_vs_v3_abundance_convergence', out)


def stage_redshift():
    """Is the force / timestep sensitivity monotonic in mass and redshift?"""
    conv, asmt, rows = load()
    out = []
    for l, r, kind in [('E','F','force'), ('F','T','time'), ('C','D','force'), ('C','E','particle')]:
        for i in range(14):
            rec = dict(pair=f'{l}/{r}', kind=kind, log10_mass=round(float(EDGES[i]),2))
            keep = True
            for z in (2,1,0):
                row = rows[(l,r,z,300)]
                if not row['valid_mass_bins'][i]: keep = False; break
                s = row['matched_statistics']['bound_mass'][i]
                v = row['matched_statistics']['vmax'][i]
                if s['count'] < 30 or s['q16_median_q84'][1] is None: keep = False; break
                rec[f'dM_z{z}'] = round(s['q16_median_q84'][1], 3)
                rec[f'dV_z{z}'] = round(v['q16_median_q84'][1], 3) if v['q16_median_q84'][1] is not None else None
            if not keep: continue
            m = [rec['dM_z2'], rec['dM_z1'], rec['dM_z0']]
            rec['monotonic_in_z'] = bool(all(a <= b for a,b in zip(m, m[1:])) or all(a >= b for a,b in zip(m, m[1:])))
            rec['sign_changes_with_z'] = bool(min(m) < 0 < max(m))
            out.append(rec)
    store('redshift_behaviour', dict(rows=out,
        meaning='Matched median bound-mass and Vmax shifts at the same mass bin across the three outputs. '
                'Non-monotonic or sign-changing behaviour forbids extrapolating a z=0 range to other epochs.'))


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('stage')
    a = p.parse_args()
    globals()['stage_' + a.stage]()
