"""Small independent controls for native/analyze_replays.py; no simulations.

Execute with micromamba run -n cosemu python3 -B. The repository being reviewed
is read only; temporary raw tapes and indexes are removed on exit. A receipt
must have a new output path, so earlier evidence is never overwritten.
"""
import argparse
import copy
import hashlib
import importlib
import json
from pathlib import Path
import shlex
import struct
import subprocess
import sys
import tempfile
import time

import numpy as np


def properties(candidate, count, box, mass, changed=False):
    p = np.zeros(21, dtype=np.float32)
    p[:3] = {
        1: [.1 * box] * 3,
        5: [box - .1, .2 * box, .2 * box],
        9: [.4 * box] * 3,
        11: [.6 * box] * 3,
    }[candidate]
    p[3:6] = [100, 200, 300] if candidate == 5 else [50, 50, 50]
    p[6] = count * mass
    p[7] = 1.e13
    p[8] = .4
    p[9:11] = [5.e16, 1.e17]
    p[11:18] = [300, .07, .1, .035, .15, .8, .6]
    p[18:21] = [1, 0, 0]
    if changed and candidate == 5:
        p[0] = .1
        p[3:6] = [130, 240, 300]
        p[7:9] = [1.1e13, .44]
        p[11:16] = [330, .08, .2, .04, .12]
    if changed and candidate == 9:
        p[4] += 20
    return p


def catalogue(props, counts):
    """Construct the published format from the documented raw-tape fields."""
    p = np.asarray(props, dtype=np.float64)
    out = np.zeros((len(p), 24), dtype=np.float64)
    out[:, :8] = p[:, :8]
    out[:, 8] = np.float32(1000 * p[:, 8])
    out[:, 9] = np.float32(np.sqrt(2 * p[:, 9] / p[:, 6]))
    out[:, 10] = p[:, 11]
    out[:, 11] = np.arange(1, len(p) + 1)
    out[:, 12] = 5
    out[:, 13] = counts
    out[:, 15] = p[:, 13]
    out[:, 16] = np.float32(2 * p[:, 9] / p[:, 10] - 1)
    out[:, 17] = p[:, 14]
    out[:, 18] = np.float32(1000 * p[:, 15])
    out[:, 19:] = p[:, 16:]
    return out


def build_pair(tmp, box, checker):
    spec = dict(nrow=int(2 * box), box_mpc_h=box,
                cosmology={'Omega_m': .3089})
    mass = np.float32(2.774e11 * .3089 * (box / spec['nrow'])**3)
    lists = [
        [(1, 2001, 350), (5, 1, 400), (9, 1001, 280)],
        [(5, 1, 450), (9, 1001, 280), (11, 3001, 360)],
    ]
    pair = []
    for side, rows in enumerate(lists):
        path = tmp / f'box{int(box)}-side{side}.bin'
        props = []
        with path.open('wb') as stream:
            stream.write(struct.pack('>qqf', 3, 12, mass))
            for candidate, start, count in rows:
                p = properties(candidate, count, box, mass, bool(side))
                props.append(p)
                stream.write(struct.pack('>qq', candidate, count))
                stream.write(p.astype('>f4').tobytes())
                stream.write(np.arange(start, start + count, dtype='>i8').tobytes())
        m = checker.check_memberships(
            path, spec, catalogue(props, [row[2] for row in rows]))
        m['retained_raw'] = str(path)
        pair.append(dict(membership=m, density={'sha256': 'd' * 64},
                         config_sha256='c' * 64,
                         catalogue={'fixture': path.stem}))
    return pair


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--initial-receipt', type=Path,
                        help='Optional earlier probe to retain inside this receipt')
    args = parser.parse_args()
    assert not args.output.exists(), 'Preserve the earlier receipt'
    repo = args.repo.resolve()
    native = repo / 'BDM-refine/repairs/20260906-review-followup/native'
    sys.path.insert(0, str(native))
    a = importlib.import_module('analyze_replays')
    r = a.r
    checker = r.checker()
    watched = [native / 'analyze_replays.py', native / 'run_replays.py',
               r.VALIDATION / 'run_validation.py']
    source_hashes = {str(p.relative_to(repo)): r.sha(p) for p in watched}
    started = time.monotonic()
    results, negative = [], []
    with tempfile.TemporaryDirectory(prefix='bdm-native-compare-review-') as temporary:
        tmp = Path(temporary)
        checker.ROOT = tmp
        for box in [512., 32.]:
            first, second = build_pair(tmp, box, checker)
            result = {}
            x, y = a.compare(first, second, result, box)
            assert [result[k] for k in ['reference_rows', 'v3_rows',
                    'matched_candidates', 'reference_only', 'v3_only']] == [3, 3, 2, 1, 1]
            assert result['identical_bound_sets_for_matched_candidates'] == 1
            assert result['changed_bound_sets_for_matched_candidates'] == 1
            assert result['matched_both_above_log10_mass_12p5'] == 1
            assert result['finite_ordered_shapes_and_normalized_axes']
            xp, yp = x['properties'][1], y['properties'][0]
            expected_fields = {'Mbound': 6, 'Mtotal': 7,
                               'reported_aperture_radius': 8,
                               'Vmax': 11, 'Rrms': 15}
            field_changes = {}
            for field, column in expected_fields.items():
                stats = result['percentage_changes_for_both_above_mass_cut'][field]
                expected = float(100 * (yp[column] / xp[column] - 1))
                assert stats['count'] == 1
                assert np.isclose(stats['mean'], expected, rtol=0, atol=1.e-12)
                field_changes[field] = stats['mean']
            overlap = result['all_matched_membership_overlap_fraction']
            assert np.isclose(overlap['mean'], (400 / 450 + 1) / 2,
                              rtol=0, atol=1.e-15)
            assert np.isclose(overlap['percentiles']['0'], 400 / 450,
                              rtol=0, atol=1.e-15)
            drift = result['all_matched_bulk_velocity_shift_km_s']
            assert drift['count'] == 2 and drift['mean'] == 35
            assert drift['percentiles']['100'] == 50
            expected_position = float(np.float32(.1)) - float(np.float32(box - .1)) + box
            position = result['all_matched_position_shift_mpc_h']
            assert np.isclose(position['percentiles']['100'], expected_position,
                              rtol=0, atol=1.e-12)
            assert np.isclose(position['mean'], expected_position / 2,
                              rtol=0, atol=1.e-12)
            results.append(dict(
                box_mpc_h=box, independently_parsed_catalogues=2,
                matching='Candidate IDs 5 and 9 despite shifted published row indices',
                matched=2, only_each_side=1, identical_sets=1, changed_sets=1,
                mass_cut_cohort=1, percentage_changes=field_changes,
                overlap_mean=overlap['mean'],
                maximum_periodic_shift_mpc_h=position['percentiles']['100'],
                bulk_drift_mean_km_s=drift['mean']))
            for key in ['config', 'density', 'raw']:
                broken = copy.deepcopy(second)
                if key == 'config':
                    broken['config_sha256'] = 'e' * 64
                elif key == 'density':
                    broken['density']['sha256'] = 'e' * 64
                else:
                    broken['membership']['raw_sha256'] = 'e' * 64
                try:
                    a.compare(first, broken, {}, box)
                except AssertionError:
                    negative.append(dict(box_mpc_h=box,
                                         rejected=key + ' hash mismatch'))
                else:
                    raise AssertionError('Mismatch was accepted: ' + key)

        build = json.loads((native / 'build.json').read_text())
        frozen = {}
        for variant, entry in build['variants'].items():
            path = Path(entry['source_path'])
            actual = r.sha(path)
            assert actual == entry['source_sha256']
            frozen[variant] = dict(source_sha256=actual, routines={})
            for name, kind in [('FindMaxima', 'subroutine'),
                               ('SetOverdensity', 'subroutine'),
                               ('IsDensityMaximum', 'function')]:
                definition = a.routine(path.read_text(), name, kind)
                frozen[variant]['routines'][name] = hashlib.sha256(
                    definition.encode()).hexdigest()
        for name in ['FindMaxima', 'SetOverdensity', 'IsDensityMaximum']:
            assert len({v['routines'][name] for v in frozen.values()}) == 1
        boundary = np.float32(1.25)
        seeds = np.array([np.nextafter(boundary, np.float32(-np.inf)),
                          np.nextafter(boundary, np.float32(np.inf))], dtype=np.float32)
        quantized = np.rint(seeds.astype(np.float64) / .1) * .1
        assert quantized[0] != quantized[1]
        partial_radius = (120 / 100)**(1 / 3)
        all_radius = (320 / 100)**(1 / 3)
        assert .99 < partial_radius < 2 and all_radius < 2
        assert all(r.sha(p) == source_hashes[str(p.relative_to(repo))] for p in watched)

    receipt = dict(
        completed=True,
        scope='Bounded independent read-only review; actual membership checker and '
              'compare(...,box), synthetic tapes only; no simulations or full native analysis',
        root_head=subprocess.check_output(
            ['git', '-C', str(repo), 'rev-parse', 'HEAD'], text=True).strip(),
        command=shlex.join(['micromamba', 'run', '-n', 'cosemu', 'python3', '-B', *sys.argv]),
        test_source_sha256=r.sha(__file__), source_sha256=source_hashes,
        elapsed_seconds=time.monotonic() - started,
        positive_controls=results, negative_controls=negative,
        frozen_source_proof=frozen,
        interpretation=dict(
            membership_overlap_formula='Intersection cardinality divided by max(left '
                                       'cardinality,right cardinality); not Jaccard',
            fractional_property_cohort='Matched initial candidate IDs with both Mbound '
                                       '>= 10^12.5 Msun/h and both reported quantities positive',
            candidate_identity_scope='Same fixed FI and identical peak routines, background, '
                'snapshot and paired config. Unmatched published candidate IDs are not '
                'necessarily new or lost physical objects.',
            quantization_counterexample=dict(
                input_float32_seed_factors=[float(x) for x in seeds],
                round_to_point_one=[float(x) for x in quantized],
                conclusion='Finite quantization has boundaries and cannot guarantee global '
                           'reproducibility; ordinary density changes can also change the candidate set.'),
            partial_companion_counterexample=dict(
                normalized_mass_over_so_coefficient=dict(primary=1.0, companion_arm=.2,
                                                          companion_core=2.0),
                distances_from_primary=dict(primary_rows=.1, arm_rows=.99, core_rows=2.0),
                selected_so_radius=partial_radius,
                hypothetical_all_mass_radius=all_radius,
                conclusion='Arm particles can contribute with a distant companion centre '
                           'outside the selected primary radius; the equal-clump fixture '
                           'is no general bridging/incidence theorem.')))
    if args.initial_receipt:
        receipt['initial_probe'] = dict(
            sha256=r.sha(args.initial_receipt),
            result=json.loads(args.initial_receipt.read_text()))
    with args.output.open('x') as stream:
        json.dump(receipt, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps(dict(completed=True, receipt=str(args.output),
                          positive_box_fixtures=len(results),
                          valid_tapes_parsed=2 * len(results),
                          negative_controls=len(negative),
                          elapsed_seconds=receipt['elapsed_seconds'],
                          source_sha256=source_hashes), indent=2))


if __name__ == '__main__':
    main()
