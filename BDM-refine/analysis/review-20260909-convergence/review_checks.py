"""Independent controls for the 20260908 BDM convergence campaign review.

Run every stage with: micromamba run -n cosemu python3 -B review_checks.py <stage>
Stages: assess | recompute | tape | matchfixture | matchwiden | report
Nothing here writes to the campaign directory; results go to review_results.json.
"""
import argparse, hashlib, json, os, struct, sys, time
from pathlib import Path

os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('MKL_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
from scipy.spatial import cKDTree

REPO = Path('/cosma8/data/dp203/dc-ruan1/mgglam_claude/MG-GLAM')
C = REPO / 'BDM-refine/validation/20260908-convergence'
HERE = Path(__file__).resolve().parent
RESULTS = HERE / 'review_results.json'
EDGES = np.arange(12.25, 15.76, .25)
PUB_MIN = 2.5e12


def sha(path, chunk=8 * 1024**2):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        while b := f.read(chunk):
            h.update(b)
    return h.hexdigest()


def store(key, value):
    data = json.loads(RESULTS.read_text()) if RESULTS.exists() else {}
    data[key] = value
    RESULTS.write_text(json.dumps(data, indent=1, default=str) + '\n')
    print(f'[{key}] stored')


def receipt(name, z):
    return json.loads((C / f'v3-validation-{name}-z{z}-ng2048.json').read_text())


def load_index(name, z, verify_index=True):
    """Load the retained membership index; optionally re-hash the small .npz itself."""
    rep = receipt(name, z)
    ev = rep['variants']['v3']['membership']
    raw = REPO / ev['retained_raw']
    idx = raw.with_suffix('.index.npz')
    got = sha(idx) if verify_index else None
    with np.load(idx) as d:
        arrays = {k: d[k] for k in d.files}
    return dict(name=name, z=z, rep=rep, raw=raw, index=idx, index_sha=got,
                index_sha_expected=ev['index_sha256'], mass_one=ev['mass_one'],
                box=rep['spec']['box_mpc_h'], nrow=rep['spec']['nrow'], **arrays)


# --------------------------------------------------------------------------- #
def stage_assess():
    """Re-derive all 63 assessments from convergence.json with independent code."""
    conv = json.loads((C / 'convergence.json').read_text())
    asmt = json.loads((C / 'convergence-assessment.json').read_text())
    ok = lambda v: v is not None and np.isfinite(v)
    mism, table = [], []
    prod = {(r['coarse'], r['reference'], r['redshift'], r['particle_floor']): r for r in asmt['comparisons']}
    for row in conv['comparisons']:
        k = (row['coarse'], row['reference'], row['redshift'], row['particle_floor'])
        passes = []
        for i, valid in enumerate(row['valid_mass_bins']):
            m = row['matched_statistics']['bound_mass'][i]
            v = row['matched_statistics']['vmax'][i]
            n = [row['left_counts'][i], row['right_counts'][i],
                 row['reference_eligible_counts'][i], m['count'], v['count']]
            r_, s_ = row['abundance_ratio'][i], row['abundance_ratio_jackknife8_sigma'][i]
            mm, mv, cp = m['q16_median_q84'][1], v['q16_median_q84'][1], row['reference_completeness'][i]
            elig = bool(valid) and min(n) >= 30 and all(ok(x) for x in (r_, s_, mm, mv, cp))
            ch = dict(abundance=abs(100 * (r_ - 1)) <= 5., abundance_uncertainty=100 * s_ <= 5.,
                      bound_mass=abs(mm) <= 5., resolved_vmax=abs(mv) <= 2.,
                      completeness=cp >= .9) if elig else {}
            p = prod[k]['bins'][i]
            if elig != p['eligible'] or ch != p['checks'] or bool(elig and all(ch.values())) != p['meets_working_criteria']:
                mism.append([str(k), i])
            passes.append(bool(elig and all(ch.values())))
        iv, start = [], None
        for i, s in enumerate([*passes, False]):
            if s and start is None: start = i
            if not s and start is not None:
                iv.append([float(EDGES[start]), float(EDGES[i])]); start = None
        if iv != prod[k]['meets_working_criteria_log10_mass_intervals']:
            mism.append([str(k), 'intervals'])
        table.append(dict(pair=f'{k[0]}/{k[1]}', z=k[2], floor=k[3], intervals=iv,
                          eligible=sum(1 for i, b in enumerate(prod[k]['bins']) if b['eligible'])))
    store('independent_assessment', dict(reproduced=len(conv['comparisons']), mismatches=mism,
                                         assess_py_source_hash_recorded_in_assessment=False,
                                         intervals_by_floor=table))


# --------------------------------------------------------------------------- #
def stage_recompute():
    """Rebuild compare() outputs for every pair from retained index+match arrays only."""
    conv = json.loads((C / 'convergence.json').read_text())
    got = {(r['coarse'], r['reference'], r['redshift'], r['particle_floor']): r for r in conv['comparisons']}
    pairs = [('A','C'),('C','E'),('A','E'),('B','C'),('C','D'),('E','F'),('F','T')]
    diffs, index_hash = [], []
    cache = {}
    for z in (2, 1, 0):
        for l, r in pairs:
            for nm in (l, r):
                if (nm, z) not in cache:
                    s = load_index(nm, z)
                    index_hash.append([f'{nm}/z{z}', s['index_sha'] == s['index_sha_expected']])
                    cache[nm, z] = s
            L, R = cache[l, z], cache[r, z]
            mp = C / f'work/analysis/matches-{l}-{r}-z{z}.npz'
            rec = json.loads(mp.with_suffix('.json').read_text())
            assert sha(mp) == rec['sha256'], mp
            with np.load(mp) as d:
                P = d['pairs']
            ii, jj = P[:, 0].astype(int), P[:, 1].astype(int)
            mL = L['counts'].astype(np.float64) * L['mass_one']
            mR = R['counts'].astype(np.float64) * R['mass_one']
            box = L['box']
            def hist(m, sel=None):
                x = np.log10(m if sel is None else m[sel])
                return np.histogram(x, EDGES)[0]
            cL, cR = hist(mL), hist(mR)
            oc = np.floor((L['properties'][:, :3] % box) / (box / 2)).astype(int)
            octL = oc[:, 0] + 2 * oc[:, 1] + 4 * oc[:, 2]
            oc = np.floor((R['properties'][:, :3] % box) / (box / 2)).astype(int)
            octR = oc[:, 0] + 2 * oc[:, 1] + 4 * oc[:, 2]
            OL = np.array([hist(mL, octL == q) for q in range(8)])
            OR = np.array([hist(mR, octR == q) for q in range(8)])
            for floor in (100, 300, 1000):
                key = (l, r, z, floor)
                g = got[key]
                mf = max(PUB_MIN, floor * max(L['mass_one'], R['mass_one']))
                vb = 10**EDGES[:-1] >= mf
                ratio = np.divide(cL, cR, out=np.full(14, np.nan), where=cR > 0); ratio[~vb] = np.nan
                with np.errstate(invalid='ignore'):
                    rr = np.divide(cL - OL, cR - OR, out=np.full((8, 14), np.nan), where=(cR - OR) > 0)
                jk = np.sqrt(7 / 8 * np.sum((rr - np.mean(rr, axis=0))**2, axis=0)); jk[~vb] = np.nan
                res = (L['counts'][ii] >= floor) & (R['counts'][jj] >= floor)
                x = np.log10(mR[jj])
                a, b = L['properties'][ii], R['properties'][jj]
                out = {}
                for lab, col in [('bound_mass', 6), ('aperture_total_mass', 7), ('aperture_radius', 8),
                                 ('vmax', 11), ('axis_ba', 16), ('axis_ca', 17)]:
                    av, bv = (mL[ii], mR[jj]) if lab == 'bound_mass' else (a[:, col], b[:, col])
                    val = 100 * np.divide(av - bv, bv, out=np.full(len(a), np.nan), where=bv > 0)
                    use = res.copy()
                    if lab == 'vmax': use &= (a[:, col] > 0) & (b[:, col] > 0)
                    out[lab] = [(int(np.sum(sel := use & (x >= lo) & (x < hi) & np.isfinite(val))),
                                 np.quantile(val[sel], [.16, .5, .84]).tolist() if sel.any() else [None] * 3)
                                for lo, hi in zip(EDGES[:-1], EDGES[1:])]
                elig = mR >= mf
                refc = hist(mR, elig)
                selc = np.histogram(x[res & elig[jj]], EDGES)[0]
                comp = np.divide(selc, refc, out=np.full(14, np.nan), where=refc > 0)
                lec = hist(mL, L['counts'] >= floor)
                mlf = np.divide(np.histogram(np.log10(mL[ii[res]]), EDGES)[0], lec,
                                out=np.full(14, np.nan), where=lec > 0)
                def cmp(name, mine, theirs):
                    t = np.array([np.nan if v is None else v for v in theirs], dtype=float)
                    if not np.allclose(np.asarray(mine, dtype=float), t, equal_nan=True, rtol=1e-9, atol=1e-12):
                        diffs.append([str(key), name])
                cmp('left_counts', cL, g['left_counts']); cmp('right_counts', cR, g['right_counts'])
                cmp('abundance_ratio', ratio, g['abundance_ratio'])
                cmp('jackknife', jk, g['abundance_ratio_jackknife8_sigma'])
                cmp('completeness', comp, g['reference_completeness'])
                cmp('left_matched_fraction', mlf, g['left_matched_fraction'])
                cmp('reference_eligible', refc, g['reference_eligible_counts'])
                if abs(mf - g['mass_floor_msun_h']) > 1e-3: diffs.append([str(key), 'mass_floor'])
                for lab in out:
                    cmp(lab + '.count', [c for c, _ in out[lab]],
                        [d['count'] for d in g['matched_statistics'][lab]])
                    for q in range(3):
                        cmp(f'{lab}.q{q}', [v[q] for _, v in out[lab]],
                            [d['q16_median_q84'][q] for d in g['matched_statistics'][lab]])
            print(f'  recomputed {l}/{r} z={z}', flush=True)
    store('independent_recompute', dict(pairs=len(pairs) * 3 * 3, disagreements=diffs,
                                        index_npz_rehashed=index_hash))


# --------------------------------------------------------------------------- #
def stage_tape():
    """Independently parse the raw membership tape for one catalogue; no receipt trust."""
    s = load_index('E', 0)
    rep = s['rep']['variants']['v3']['membership']
    t0 = time.time()
    raw_sha = sha(s['raw'])
    out = dict(catalogue='E/z0', raw=str(s['raw']), bytes=s['raw'].stat().st_size,
               raw_sha256_recomputed=raw_sha, raw_sha256_matches_receipt=raw_sha == rep['raw_sha256'],
               index_sha256_matches_receipt=s['index_sha'] == s['index_sha_expected'],
               hash_seconds=round(time.time() - t0, 1))
    with s['raw'].open('rb') as f:
        sel, cand, m1 = struct.unpack('>qqf', f.read(20))
        out.update(header_selected=sel, header_candidates=cand, header_mass_one=m1,
                   mass_one_matches_receipt=m1 == rep['mass_one'],
                   selected_matches_index=sel == len(s['counts']))
        expected = 2.774e11 * 0.3089 * (256. / 1024)**3
        out['mass_one_vs_2p774e11_formula'] = dict(expected=expected, relative=abs(m1 / expected - 1))
        rng = np.random.default_rng(20260909)
        probe = rng.choice(sel, 400, replace=False)
        probe = np.append(probe, [0, sel - 1, int(np.argmax(s['counts']))])
        bad, sorted_ok, range_ok, propok = [], True, True, True
        allids = []
        for k in probe:
            f.seek(int(s['offsets'][k]) - 100)
            cid, cnt = struct.unpack('>qq', f.read(16))
            props = np.frombuffer(f.read(84), dtype='>f4').astype(np.float64)
            ids = np.frombuffer(f.read(int(s['counts'][k]) * 8), dtype='>i8')
            if cid != s['candidates'][k] or cnt != s['counts'][k]: bad.append(int(k))
            if not np.allclose(props, s['properties'][k], rtol=0, atol=0): propok = False
            if not np.all(ids[1:] > ids[:-1]): sorted_ok = False
            if ids.min() < 1 or ids.max() > 1024**3: range_ok = False
            if abs(props[6] - cnt * m1) > 2 * abs(float(np.spacing(np.float32(props[6])))): bad.append(f'mass{k}')
            allids.append(ids)
        out.update(rows_probed=len(probe), row_header_mismatches=bad,
                   properties_match_index=propok, ids_strictly_sorted=sorted_ok, ids_in_range=range_ok)
        # cross-halo membership disjointness on the probe (physical uniqueness, not exact-set duplication)
        cat = np.concatenate(allids)
        uniq = len(np.unique(cat))
        out['probe_membership_overlap'] = dict(total_ids=int(len(cat)), unique_ids=int(uniq),
                                               shared_ids=int(len(cat) - uniq))
        # host exclusion, independently, over the full catalogue
        box = s['box']; p = s['properties']
        tree = cKDTree(p[:, :3] % box, boxsize=box)
        pairs = tree.query_pairs(float(np.max(p[:, 8])), output_type='ndarray')
        d = p[pairs[:, 0], :3] - p[pairs[:, 1], :3]
        d -= box * np.rint(d / box)
        dist = np.linalg.norm(d, axis=1)
        r0, r1 = p[pairs[:, 0], 8], p[pairs[:, 1], 8]
        inside = (dist < np.maximum(r0, r1))
        out['host_exclusion'] = dict(
            pairs_within_max_radius=int(np.sum(inside)),
            centre_inside_larger_radius=int(np.sum(dist < np.maximum(r0, r1))),
            centre_inside_smaller_radius=int(np.sum(dist < np.minimum(r0, r1))),
            closest_pair_over_sum_radii=float(np.min(dist / (r0 + r1))) if len(dist) else None,
            note='checker asserts no centre lies within a MORE massive halo radius; both counts here are over all pairs')
    store('raw_tape_independent_parse', out)


# --------------------------------------------------------------------------- #
def _load_analyze():
    sys.path.insert(0, str(C))
    import importlib
    return importlib.import_module('analyze')


def stage_matchfixture():
    """Exercise analyze.match / project_ids against exhaustive brute force on fixtures."""
    an = _load_analyze()
    rng = np.random.default_rng(7)
    results = []

    # -- project_ids: exact nested-lattice identity against an explicit construction
    for nf, nc in [(8, 4), (16, 4), (16, 8), (12, 3)]:
        if nf % nc: continue
        r = nf // nc
        ids = np.arange(1, nf**3 + 1)
        proj = an.project_ids(ids, nf, nc)
        xs, ys, zs = np.meshgrid(np.arange(nc), np.arange(nc), np.arange(nc), indexing='ij')
        want = np.sort(1 + xs.ravel() + nc * ys.ravel() + nc**2 * zs.ravel())
        results.append(dict(check='project_ids_full_lattice', nfine=nf, ncoarse=nc,
                            bijective=bool(np.array_equal(np.sort(proj), want)),
                            kept=len(proj), expected=nc**3))
        # Lagrangian equivalence: fine index r*i maps to coarse index i at identical position
        i0 = rng.integers(0, nc, 50); j0 = rng.integers(0, nc, 50); k0 = rng.integers(0, nc, 50)
        fine_id = 1 + (r * i0) + nf * (r * j0) + nf**2 * (r * k0)
        coarse_id = 1 + i0 + nc * j0 + nc**2 * k0
        same_pos = np.allclose((r * i0) / nf, i0 / nc)
        results.append(dict(check='project_ids_lagrangian', nfine=nf, ncoarse=nc,
                            ids_match=bool(np.array_equal(an.project_ids(fine_id, nf, nc), coarse_id)),
                            positions_identical=bool(same_pos)))

    # -- synthetic membership tapes -> match() vs exhaustive matching
    def build(tmp, tag, nrow, centres, radii, memberlists, box=256.):
        raw = tmp / f'members-{tag}.bin'
        offs, cnts, props = [], [], []
        with raw.open('wb') as f:
            f.write(struct.pack('>qqf', len(centres), len(centres), 1.e10))
            for n, (c, rad, ids) in enumerate(zip(centres, radii, memberlists)):
                f.write(struct.pack('>qq', n + 1, len(ids)))
                p = np.zeros(21, dtype='>f4')
                p[0:3] = c; p[6] = len(ids) * 1.e10; p[7] = p[6] * 1.1; p[8] = rad
                p[11] = 100. + n; p[16] = .8; p[17] = .6
                f.write(p.tobytes())
                offs.append(f.tell())
                f.write(np.sort(np.asarray(ids, dtype='>i8')).tobytes())
                cnts.append(len(ids)); props.append(p.astype(np.float64))
        # mirror the checker: properties are always (-1, 21)
        return dict(name=f'N{nrow}', z=0, raw=raw,
                    report=dict(spec=dict(nrow=nrow, box_mpc_h=box),
                                variants=dict(v3=dict(membership=dict(mass_one=1.e10)))),
                    counts=np.array(cnts), offsets=np.array(offs),
                    properties=np.asarray(props, dtype=np.float64).reshape(-1, 21),
                    candidates=np.arange(1, len(centres) + 1))

    def exhaustive(L, R, box=256.):
        nl, nr = L['report']['spec']['nrow'], R['report']['spec']['nrow']
        setsL, setsR = [], []
        with L['raw'].open('rb') as f:
            for i in range(len(L['counts'])):
                f.seek(int(L['offsets'][i])); setsL.append(np.frombuffer(f.read(int(L['counts'][i]) * 8), dtype='>i8'))
        with R['raw'].open('rb') as f:
            for j in range(len(R['counts'])):
                f.seek(int(R['offsets'][j])); setsR.append(an.project_ids(
                    np.frombuffer(f.read(int(R['counts'][j]) * 8), dtype='>i8').astype(np.int64), nr, nl))
        bl, br = {}, {}
        for i in range(len(setsL)):
            for j in range(len(setsR)):
                if not len(setsR[j]): continue
                ov = len(np.intersect1d(setsL[i], setsR[j], assume_unique=True))
                fa, fb = ov / len(setsL[i]), ov / len(setsR[j])
                if min(fa, fb) < .5: continue
                d = L['properties'][i, :3] - R['properties'][j, :3]
                d -= box * np.rint(d / box)
                dist = float(np.linalg.norm(d))
                sc = fa * fb
                if i not in bl or (sc, -dist, -j) > bl[i][:3]: bl[i] = (sc, -dist, -j, j)
                if j not in br or (sc, -dist, -i) > br[j][:3]: br[j] = (sc, -dist, -i, i)
        return sorted((i, v[3]) for i, v in bl.items() if br[v[3]][3] == i)

    import tempfile
    cases = []
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        # case 1: clean one-to-one, coarse 8^3 / fine 16^3
        nl, nf = 8, 16
        cen = np.array([[10., 10, 10], [50, 50, 50], [200, 30, 90], [250, 250, 250]])
        cenR = cen + np.array([[.01, 0, 0], [0, .02, 0], [-.01, 0, .01], [0, 0, -.02]])
        memL = [np.arange(1, 121) + 200 * n for n in range(4)]
        memR = [an_ids for an_ids in [
            np.array([1 + (2 * ((v - 1) % nl)) + nf * (2 * (((v - 1) // nl) % nl)) + nf**2 * (2 * ((v - 1) // nl**2))
                      for v in m]) for m in memL]]
        L = build(tmp, 'c1L', nl, cen, [2., 2, 2, 2], memL); R = build(tmp, 'c1R', nf, cenR, [2., 2, 2, 2], memR)
        cases.append(('clean_nested', L, R))
        # case 2: periodic wrap + split object + empty overlap
        cen2 = np.array([[.05, .05, .05], [128., 128, 128], [60., 60, 60]])
        cenR2 = np.array([[255.98, 255.98, 255.98], [128.4, 128, 128], [180., 60, 60]])
        memL2 = [np.arange(1, 101), np.arange(200, 300), np.arange(400, 500)]
        conv = lambda m: np.array([1 + (2 * ((v - 1) % nl)) + nf * (2 * (((v - 1) // nl) % nl)) + nf**2 * (2 * ((v - 1) // nl**2)) for v in m])
        memR2 = [conv(np.arange(1, 101)), conv(np.arange(200, 260)), conv(np.arange(900, 1000))]
        L2 = build(tmp, 'c2L', nl, cen2, [1., 3., 1.], memL2); R2 = build(tmp, 'c2R', nf, cenR2, [1., 3., 1.], memR2)
        cases.append(('periodic_split_disjoint', L2, R2))
        # case 3: two coarse haloes competing for one fine halo (tie/mutual-best)
        cen3 = np.array([[100., 100, 100], [100.3, 100, 100]])
        cenR3 = np.array([[100.15, 100, 100]])
        memL3 = [np.arange(1, 101), np.arange(1, 101)]
        memR3 = [conv(np.arange(1, 101))]
        L3 = build(tmp, 'c3L', nl, cen3, [1., 1.], memL3); R3 = build(tmp, 'c3R', nf, cenR3, [1.], memR3)
        cases.append(('competing_mutual_best', L3, R3))
        # case 4: empty right catalogue
        L4 = build(tmp, 'c4L', nl, cen[:2], [2., 2.], memL[:2])
        R4 = build(tmp, 'c4R', nf, np.zeros((0, 3)), [], [])
        cases.append(('empty_reference', L4, R4))
        for label, A, B in cases:
            try:
                pairs, desc = an.match(A, B)
                mine = sorted((int(p[0]), int(p[1])) for p in pairs)
            except Exception as e:
                results.append(dict(check='match_vs_exhaustive', case=label, error=repr(e))); continue
            ref = exhaustive(A, B)
            results.append(dict(check='match_vs_exhaustive', case=label, production=mine,
                                exhaustive=ref, agree=mine == ref,
                                spatial_candidates=desc['spatial_candidates'], matches=desc['matches']))
    store('match_fixture_controls', results)


# --------------------------------------------------------------------------- #
def stage_matchwiden(pair='E,F', z=0, half=32.0, factor=6.0):
    """On real data, redo matching in a subvolume with a much wider candidate radius."""
    an = _load_analyze()
    l, r = pair.split(',')
    L, R = load_index(l, z, verify_index=False), load_index(r, z, verify_index=False)
    box = L['box']; nl, nr = L['nrow'], R['nrow']
    a, b = L['properties'], R['properties']
    centre = np.array([box / 2] * 3)
    d = a[:, :3] - centre; d -= box * np.rint(d / box)
    sub = np.where(np.all(np.abs(d) <= half, axis=1))[0]
    tree = cKDTree(b[:, :3] % box, boxsize=box)
    wide = factor * (np.maximum(a[sub, 8], np.max(b[:, 8])) + box / 2048)
    cand = tree.query_ball_point(a[sub, :3] % box, wide, workers=8)
    cacheL, cacheR = {}, {}
    def idsL(i):
        if i not in cacheL:
            with L['raw'].open('rb') as f:
                f.seek(int(L['offsets'][i])); cacheL[i] = np.frombuffer(f.read(int(L['counts'][i]) * 8), dtype='>i8').astype(np.int64)
        return cacheL[i]
    def idsR(j):
        if j not in cacheR:
            with R['raw'].open('rb') as f:
                f.seek(int(R['offsets'][j])); v = np.frombuffer(f.read(int(R['counts'][j]) * 8), dtype='>i8').astype(np.int64)
            cacheR[j] = an.project_ids(v, nr, nl)
        return cacheR[j]
    best_l, best_r, examined, narrow_examined = {}, {}, 0, 0
    for n, i in enumerate(sub):
        for j in cand[n]:
            dd = a[i, :3] - b[j, :3]; dd -= box * np.rint(dd / box)
            dist = float(np.linalg.norm(dd))
            narrow = dist <= max(a[i, 8], b[j, 8]) + box / 2048
            IL, IR = idsL(i), idsR(j)
            if not len(IR): continue
            examined += 1; narrow_examined += narrow
            ov = len(np.intersect1d(IL, IR, assume_unique=True))
            fa, fb = ov / len(IL), ov / len(IR)
            if min(fa, fb) < .5: continue
            sc = fa * fb
            if i not in best_l or (sc, -dist, -j) > best_l[i][:3]: best_l[i] = (sc, -dist, -j, j, dist, narrow)
            if j not in best_r or (sc, -dist, -i) > best_r[j][:3]: best_r[j] = (sc, -dist, -i, i)
    wide_pairs = {i: v for i, v in best_l.items() if best_r[v[3]][3] == i}
    with np.load(C / f'work/analysis/matches-{l}-{r}-z{z}.npz') as dz:
        P = dz['pairs']
    prod = {int(p[0]): int(p[1]) for p in P}
    subset = set(sub.tolist())
    prod_sub = {i: j for i, j in prod.items() if i in subset}
    only_wide = {int(i): int(v[3]) for i, v in wide_pairs.items() if prod_sub.get(i) != v[3]}
    only_prod = {int(i): int(j) for i, j in prod_sub.items() if i not in wide_pairs}
    mL = L['counts'].astype(float) * L['mass_one']
    store('wider_radius_rematch', dict(
        pair=f'{l}/{r}', z=z, subvolume_half_mpc_h=half, radius_factor=factor,
        left_haloes_in_subvolume=len(sub), production_matches_in_subvolume=len(prod_sub),
        wide_matches=len(wide_pairs), pairs_examined_wide=examined,
        pairs_examined_within_production_filter=int(narrow_examined),
        matches_changed_by_widening=only_wide,
        production_matches_lost_when_widened=only_prod,
        max_matched_centre_distance_mpc_h=float(max((v[4] for v in wide_pairs.values()), default=0.)),
        production_filter_radius_example_mpc_h=float(max(a[sub, 8]) + box / 2048),
        unmatched_left_above_pub_cut=int(np.sum((mL[sub] >= PUB_MIN) & ~np.isin(sub, list(wide_pairs))))))


# --------------------------------------------------------------------------- #
def stage_disjoint():
    """Complete physical-uniqueness control: no particle in two published haloes.

    This is stronger than the campaign's exact-duplicate-member-set check and than
    its null host-exclusion count, and it is run over every membership of all 21
    published v3 catalogues.
    """
    out = []
    for z in (2, 1, 0):
        for name in 'ABCDEFT':
            s = load_index(name, z, verify_index=False)
            buf = s['raw'].read_bytes()
            n = int(s['counts'].sum())
            allids = np.empty(n, dtype=np.int64)
            at = 0
            for off, cnt in zip(s['offsets'], s['counts']):
                c = int(cnt)
                allids[at:at + c] = np.frombuffer(buf, dtype='>i8', count=c, offset=int(off))
                at += c
            del buf
            assert at == n
            allids.sort()
            dup = int(np.count_nonzero(allids[1:] == allids[:-1]))
            lo, hi = int(allids[0]), int(allids[-1])
            # host exclusion over every pair, independent of the checker's priority rule
            p = s['properties']; box = s['box']
            tree = cKDTree(p[:, :3] % box, boxsize=box)
            pr = tree.query_pairs(float(np.max(p[:, 8])), output_type='ndarray')
            if len(pr):
                d = p[pr[:, 0], :3] - p[pr[:, 1], :3]; d -= box * np.rint(d / box)
                dist = np.linalg.norm(d, axis=1)
                r0, r1 = p[pr[:, 0], 8], p[pr[:, 1], 8]
                inside = int(np.sum(dist < np.maximum(r0, r1)))
                closest = float(np.min(dist / (r0 + r1)))
            else:
                inside, closest = 0, None
            out.append(dict(catalogue=f'{name}/z{z}', haloes=len(s['counts']), memberships=n,
                            shared_particle_ids=dup, id_min=lo, id_max=hi,
                            id_range_valid=bool(lo >= 1 and hi <= s['nrow']**3),
                            centre_inside_any_other_radius=inside,
                            min_separation_over_sum_radii=closest))
            print(f"  {name}/z{z}: {len(s['counts'])} haloes, {n} memberships, shared={dup}, "
                  f"centre-in-radius={inside}", flush=True)
            del allids
    store('full_membership_disjointness', out)



# --------------------------------------------------------------------------- #
def stage_hashes():
    """Re-hash every raw membership tape and every science file below 1 GiB."""
    tapes, files, big = [], [], []
    for z in (2, 1, 0):
        for name in 'ABCDEFT':
            rep = receipt(name, z)
            ev = rep['variants']['v3']['membership']
            raw = REPO / ev['retained_raw']
            got = sha(raw)
            tapes.append([f'{name}/z{z}', raw.stat().st_size, got == ev['raw_sha256']])
            for path, meta in rep['science_files'].items():
                p = Path(path)
                if meta['bytes'] >= 1024**3:
                    big.append([f'{name}/z{z}', p.name, meta['bytes'],
                                p.is_file() and p.stat().st_size == meta['bytes']])
                    continue
                files.append([f'{name}/z{z}', p.name, p.is_file() and sha(p) == meta['sha256']])
            d = rep['density']
            big.append([f'{name}/z{z}', 'density.bin', d['bytes'],
                        Path(d['path']).is_file() and Path(d['path']).stat().st_size == d['bytes']])
            print(f"  {name}/z{z}: tape re-hashed={tapes[-1][2]}", flush=True)
    store('independent_hash_verification', dict(
        raw_membership_tapes=len(tapes), raw_tapes_matching=sum(t[2] for t in tapes), tapes=tapes,
        small_science_files=len(files), small_files_matching=sum(f[2] for f in files),
        small_file_failures=[f for f in files if not f[2]],
        large_files_size_only=len(big), large_files_size_ok=sum(b[3] for b in big),
        large_file_note='Density tapes (32 GiB each) and the 1024^3 raw tapes: size and existence only where >=1 GiB'))

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('stage', choices=['assess', 'recompute', 'tape', 'matchfixture', 'matchwiden', 'disjoint', 'hashes'])
    p.add_argument('--pair', default='E,F'); p.add_argument('--z', type=int, default=0)
    a = p.parse_args()
    if a.stage == 'assess': stage_assess()
    elif a.stage == 'recompute': stage_recompute()
    elif a.stage == 'tape': stage_tape()
    elif a.stage == 'matchfixture': stage_matchfixture()
    elif a.stage == 'disjoint': stage_disjoint()
    elif a.stage == 'hashes': stage_hashes()
    else: stage_matchwiden(a.pair, a.z)
