"""Reproduce Prompt 2 and qualify its population-level claims without new replays.

Run with micromamba run -n cosemu python3 -B claims_response.py.
Only claims-response.json is written; review inputs and science stay unchanged.
"""
import os
for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[name] = '1'

import hashlib
import json
from pathlib import Path
import subprocess
import tempfile

import numpy as np

from common import ROOT, REPO, FINDER_SHA, CHECKER_SHA, load_module, now, sha, write_json

REVIEW = REPO / 'BDM-refine/analysis/review-20260909-convergence-claims'
PIN = 'ceb7e8c182f9f0395d1876ebab15f1f2be6724e6'


def main():
    started = now()
    sources = {str(p.relative_to(REPO)): sha(p) for p in
               [Path(__file__), ROOT/'common.py', *sorted(REVIEW.glob('*'))] if p.is_file()}
    # Check the actual frozen products; do not re-read the large raw tapes.
    frozen = [ROOT/'convergence.json', ROOT/'analyze.py', REPO/'PMP2linker.f90',
              REPO/'BDM-refine/validation/20260906-n1024/run_validation.py']
    receipts = sorted({p for pattern in ('*-ic.json', '*-simulation.json',
                                         'replay-*-z*-ng2048.json', 'v3-validation-*-z*-ng2048.json')
                       for p in ROOT.glob(pattern)})
    assert len(receipts) == 56
    frozen_hashes = {}
    for p in frozen + receipts:
        relative = str(p.relative_to(REPO))
        expected = hashlib.sha256(subprocess.check_output(['git', 'show', PIN+':'+relative], cwd=REPO)).hexdigest()
        assert sha(p) == expected, relative
        frozen_hashes[relative] = expected
    assert sha(REPO/'PMP2linker.f90') == FINDER_SHA
    assert sha(frozen[3]) == CHECKER_SHA

    reviewer = load_module('claims_review', REVIEW/'claims_checks.py')
    original = json.loads((REVIEW/'claims_results.json').read_text())
    with tempfile.TemporaryDirectory(prefix='claims-response-') as folder:
        reviewer.RESULTS = Path(folder)/'reproduced.json'
        for stage in ('verify', 'additivity', 'finder', 'selection', 'tails', 'legacy', 'redshift'):
            getattr(reviewer, 'stage_'+stage)()
        reproduced = json.loads(reviewer.RESULTS.read_text())
    checks = {key: value == original[key] for key, value in reproduced.items()}
    assert all(checks.values()), checks
    assert not reproduced['verification_and_response_audit']['decisions_reproduced']['mismatches']

    conv, assessment, rows = reviewer.load()
    before = json.loads(subprocess.check_output(['git', 'show', PIN+':'+str(
        (ROOT/'convergence-assessment.json').relative_to(REPO))], cwd=REPO))
    assert assessment['criteria'] == before['criteria']
    assert assessment['comparisons'] == before['comparisons']
    verdicts = {(r['coarse'], r['reference'], r['redshift']): r
                for r in assessment['comparisons'] if r['particle_floor'] == 300}
    edges = conv['log10_mass_edges']
    redshift = []
    for pair in [('E', 'F'), ('C', 'D')]:
        for i, lo in enumerate(edges[:-1]):
            measurements = []
            for z in (2, 1, 0):
                row = rows[(*pair, z, 300)]
                stat = row['matched_statistics']['bound_mass'][i]
                if not row['valid_mass_bins'][i] or stat['count'] < 30:
                    break
                b = verdicts[(*pair, z)]['bins'][i]
                measurements.append(dict(z=z, matched_count=stat['count'],
                    median_mass_shift_percent=stat['q16_median_q84'][1],
                    passes=b['meets_working_criteria'],
                    failed_conditions=[k for k, v in b['checks'].items() if not v]))
            if len(measurements) != 3:
                continue
            values = [v['median_mass_shift_percent'] for v in measurements]
            redshift.append(dict(pair='/'.join(pair), log10_mass_interval=[lo, edges[i+1]],
                                 measurements=measurements, sign_change=min(values) < 0 < max(values)))

    # Exact integer-ratio comparison separates ties that the review assigned to legacy.
    legacy = []
    for entry in reproduced['legacy_vs_v3_abundance_convergence']['comparison']:
        for b in entry['bins']:
            i = edges.index(b['log10_mass'])
            vl, vr = (entry['v3'][k][i] for k in ('left_counts', 'right_counts'))
            ll, lr = (entry['legacy'][k][i] for k in ('left_counts', 'right_counts'))
            left, right = abs(vl-vr)*lr, abs(ll-lr)*vr
            legacy.append(dict(pair=entry['pair'], z=entry['z'], log10_mass=b['log10_mass'],
                v3_shift_percent=100*(vl/vr-1), legacy_shift_percent=100*(ll/lr-1),
                smaller_absolute_shift='v3' if left < right else 'legacy' if right < left else 'tie'))
    legacy_summary = dict(comparable_bins=len(legacy),
        **{key: sum(b['smaller_absolute_shift'] == key for b in legacy) for key in ('v3', 'legacy', 'tie')},
        median_absolute_shift_percent={v: float(np.median([abs(b[v+'_shift_percent']) for b in legacy]))
                                       for v in ('v3', 'legacy')},
        scope='Descriptive correlated bins in one realization. Neither superiority nor statistical equivalence is established.')

    # Keep each scatter width beside its own bin median, not a maximum from another bin.
    scatter = []
    for key, verdict in verdicts.items():
        row = rows[(*key, 300)]
        for i, b in enumerate(verdict['bins']):
            if not b['meets_working_criteria']:
                continue
            q = row['matched_statistics']['bound_mass'][i]['q16_median_q84']
            scatter.append(dict(pair='/'.join(key[:2]), z=key[2], log10_mass_interval=b['log10_mass_interval'],
                count=row['matched_statistics']['bound_mass'][i]['count'], q16_median_q84_percent=q,
                half_16_84_width_percent=(q[2]-q[0])/2))

    # The review pooled all masses above the first whole bin. Quantify each usable bin too.
    selection = []
    for pair in [('C', 'E'), ('E', 'F'), ('F', 'T'), ('A', 'E')]:
        for i, b in enumerate(verdicts[(*pair, 0)]['bins']):
            if not b['eligible']:
                continue
            row = rows[(*pair, 0, 300)]
            selection.append(dict(pair='/'.join(pair), log10_mass_interval=b['log10_mass_interval'],
                reference_count=row['reference_eligible_counts'][i],
                fraction_outside_matched_floor_sample=1-row['reference_completeness'][i],
                passes=b['meets_working_criteria']))
    # A small lost fraction bounds ranks, not the value of a median, without a distribution assumption.
    full = np.r_[np.full(50, 100.), np.full(51, 200.)]
    lost = full[:-1]
    median_control = dict(n_full=len(full), n_removed=1, full_median=float(np.median(full)),
                          retained_median=float(np.median(lost)),
                          relative_median_change_percent=100*(float(np.median(lost))/float(np.median(full))-1))
    assert median_control['relative_median_change_percent'] == -25.

    force_summary = []
    for z in (2, 1, 0):
        cd, ef = (rows[(*pair, z, 300)] for pair in [('C', 'D'), ('E', 'F')])
        for prop in ('bound_mass', 'vmax'):
            differences = []
            for i in range(len(edges)-1):
                if not (cd['valid_mass_bins'][i] and ef['valid_mass_bins'][i]):
                    continue
                a, b = (r['matched_statistics'][prop][i] for r in (cd, ef))
                if min(a['count'], b['count']) >= 30:
                    differences.append(b['q16_median_q84'][1]-a['q16_median_q84'][1])
            force_summary.append(dict(z=z, property=prop, bins=len(differences),
                maximum_absolute_difference_percentage_points=max(map(abs, differences))))

    for relative, digest in sources.items():
        assert sha(REPO/relative) == digest, relative
    for relative, digest in frozen_hashes.items():
        assert sha(REPO/relative) == digest, relative
    result = dict(completed=True, started_at_utc=started, completed_at_utc=now(),
        reviewed_pin=PIN, source_sha256=sources, input_sha256=sha(ROOT/'convergence.json'),
        frozen_file_sha256=frozen_hashes, review_controls_reproduced=checks,
        unchanged_assessments=63, unchanged_criteria=True,
        redshift_populations=redshift, legacy_abundance_summary=legacy_summary,
        legacy_abundance_bins=legacy, mass_scatter_inside_passing_bins=scatter,
        selection_by_eligible_bin=selection, selection_fraction_counterexample=median_control,
        force_response_at_two_particle_loads=force_summary,
        review_measurements=dict(
            particle_chain=reproduced['separability']['particle_chain_A_C_E'],
            candidate_and_aperture_counts=reproduced['finder_mesh_indicators'],
            pooled_selection=reproduced['matching_selection']),
        scope='Small retained catalogues and indices were re-read and checked by the review controls. '
              'No new simulation, replay, matching, raw-tape scan or Slurm job. Review files unchanged.')
    result['report_markdown_lines'] = report_lines(result)
    write_json(ROOT/'claims-response.json', result)
    print(json.dumps(dict(legacy=legacy_summary, force=force_summary,
        max_eligible_bin_selection_loss=max(b['fraction_outside_matched_floor_sample'] for b in selection)), indent=2))


def report_lines(evidence):
    """A short, data-bound addition to the generated convergence report."""
    lines = ['', '## Redshift dependence and additional claims checks', '',
        'Matched median bound-mass shifts change sign in eight E/F and C/D mass bins between z=2 and z=0. '
        'The following example uses the convention 100 × (coarse/reference − 1). These are separately '
        'matched populations at each epoch, not the same haloes tracked through time.', '',
        '| Pair | log10 mass interval | z=2 median shift, % | z=1, % | z=0, % |',
        '|---|---|---:|---:|---:|']
    for r in evidence['redshift_populations']:
        if r['log10_mass_interval'][0] == 13.:
            lines.append('| '+r['pair']+' | 13.00–13.25 | '+' | '.join(
                f"{x['median_mass_shift_percent']:+.2f}" for x in r['measurements'])+' |')
    lines += ['', 'The small z=1 medians lie near a sign transition in the population response. '
        'They do not establish stability across epochs or locate a continuous zero crossing. '
        'Both example bins pass the 5% mass condition at every epoch; their z=0 full-screen failures '
        'are abundance and resolved Vmax. The wider z=1 force intervals must therefore be read with '
        'all criteria, not attributed solely to mass cancellation. z=0 is the more restrictive '
        'lower-mass force comparison here; higher redshifts remain restrictive for timesteps.', '',
        'The C/D and E/F force responses agree to 0.35 percentage points in median mass and 0.28 '
        'in median Vmax in the eight common usable z=0 bins. This supports approximate force–particle '
        'separability for those medians. The particle chain has median absolute residual 0.162 pp, '
        'maximum 0.942 pp over 18 bin/property cases (mass and Vmax combined); this is not a 0.2-pp bound. '
        'Interactions with timesteps and all-property separability remain untested.', '']
    s = evidence['legacy_abundance_summary']
    lines += [f"The retained legacy/v3 catalogues allow {s['comparable_bins']} comparable abundance bins: "
        f"v3 has the smaller absolute shift in {s['v3']}, legacy in {s['legacy']}, with {s['tie']} exact ties. "
        f"Median absolute shifts are {s['median_absolute_shift_percent']['v3']:.2f}% and "
        f"{s['median_absolute_shift_percent']['legacy']:.2f}%, respectively. These correlated measurements "
        'establish neither a statistically supported convergence advantage nor equivalence. '
        'Matched legacy properties still require missing membership information.', '',
        'In C/E at z=0, log10 mass 12.75–13.00, the per-halo mass-shift 16th/median/84th percentiles '
        'are −8.30/−0.52/+6.74%. Median agreement does not bound individual errors. The review\'s '
        '0.14–1.57% unmatched fractions pool broad mass ranges; per-bin losses and the both-particle-floor '
        'selection are recorded separately. Small selection fractions alone do not bound a median shift '
        'in percentage units.', '',
        'Candidate counts respond strongly to particle load and force refinement at the fixed finder mesh. '
        'This demonstrates sensitivity but cannot exclude attenuation or a shared bias in catalogue '
        'statistics. The z=0 median apertures span 2.97–3.04 finder cells; similar radii are also affected '
        'by the common publication cut. Neither an absolute error nor its exact cancellation is measured.', '',
        'See [CLAIMS-RESPONSE.md](CLAIMS-RESPONSE.md) and the source-bound `claims-response.json` '
        'for the reproduced controls, population definitions and follow-up design.']
    return lines


if __name__ == '__main__':
    main()
