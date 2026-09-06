"""Focused regressions for review gaps in the merged BDM repairs.

Three areas that the existing suites leave uncovered:

1. `MergeNumericalDuplicates` is exercised by exactly one two-member identical
   set (`exact_set_fallback` in halo_host_ties.py). Its merge sort, its
   representative advance across groups of three or more, its equal-length
   lexicographic comparison and its interaction with a host-removed lowest
   index are untested. A candidate-count sweep also stresses the `width>n/2`
   early exit of the bottom-up merge.
2. The unbinding loop in `GetHalo` is bounded only by its initial population
   plus one. Unlike the SO contraction, which caps at sixteen passes and then
   finishes with a sorted interval scan, nothing caps or measures the number
   of unbinding passes. These fixtures build an energy ladder that removes one
   antipodal pair per pass and report the resulting pass count.
3. `two_so_crossings` proves that the outermost crossing is selected but does
   not quantify when a companion is absorbed. The pair scan measures the
   separation at which an equal companion joins the primary SO sphere, and
   checks the self-limiting bound `Rso <= (M_enclosed/threshold)^(1/3)`.

Run with micromamba run -n cosemu python3 -B halo_gaps.py.
Only --output (default halo_gaps_results.json) is written; every build and
run directory is temporary. No production source or existing receipt changes.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import tempfile

os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('MKL_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np

import halo_host_ties as ties
import halo_review as review
import halo_followup as followup

HERE = Path(__file__).resolve().parent
GRAVITY = 4.333e-9
OM0 = float(np.float32(.3))
AEXPN = float(np.float32(.8))
OVDENS = 200.
THRESHOLD = 1.150e12 * OM0 * OVDENS


# --------------------------------------------------------------------------
# 1. Exact-membership duplicate merging beyond a single pair
# --------------------------------------------------------------------------

def lattice(index, spacing=4.):
    """Well-separated site so no default 0.1 Mpc/h radius can contain another."""
    a, b, c = index % 8, (index // 8) % 8, (index // 64) % 8
    return [1. + spacing * a, 1. + spacing * b, 1. + spacing * c]


def duplicate_cases():
    group = list(range(5000, 5020))
    other = list(range(7000, 7020))
    cases = []

    def add(name, records):
        cases.append(dict(name=name, box=32., mass_one=1., records=records))

    add('triplet_identical_sets',
        [ties.record(i, lattice(i), ids=group) for i in (1, 2, 3)])
    add('two_groups_plus_singletons',
        [ties.record(i, lattice(i), ids=group if i in (1, 4, 6) else
                     other if i in (2, 5) else range(9000 + 100 * i, 9020 + 100 * i))
         for i in range(1, 8)])
    # The lowest surviving index must win regardless of linked-list order.
    add('interleaved_groups_high_indices',
        [ties.record(i, lattice(i), ids=group if i in (7, 42, 105) else other)
         for i in (105, 42, 12, 7, 3)])
    # Equal-length sets differing only in their last identity must both survive.
    add('equal_length_last_id_differs',
        [ties.record(1, lattice(1), ids=group[:19] + [5019]),
         ties.record(2, lattice(2), ids=group[:19] + [5020])])
    # A host-removed lowest index must not become the duplicate representative.
    add('duplicate_representative_after_host_removal',
        [ties.record(1, [5., 8., 8.], radius=2., count=40),
         ties.record(2, [6., 8., 8.], ids=group),
         ties.record(3, [20., 8., 8.], ids=group),
         ties.record(4, [26., 8., 8.], ids=group)])
    # Candidate-count sweep across the bottom-up merge widths.
    for total in list(range(1, 25)) + [31, 32, 33, 40, 63, 64, 65]:
        add(f'count_sweep_{total:03d}',
            [ties.record(i, lattice(i),
                         ids=range(6000 + 100 * (i % 7), 6020 + 100 * (i % 7)) if i > 3
                         else range(1000 * i, 1000 * i + 20))
             for i in range(1, total + 1)])
    return cases


def run_duplicates(work, compiler, records):
    passed = 0
    for case in duplicate_cases():
        expected, edges = ties.oracle(case)
        for mode in ('checked', 'optimized'):
            for threads in (1, 8):
                result = ties.run(work, case, 'repaired', mode, threads, compiler)
                assert result['keep'] == expected, (result, expected, edges)
                result.update(expected_keep=expected, host_edges=edges)
                records.append(result)
                passed += 1
    return passed


# --------------------------------------------------------------------------
# 2. Unbinding pass count
# --------------------------------------------------------------------------

def antipodal_directions(pairs):
    """Antipodal unit vectors: the survivor bulk velocity stays exactly zero."""
    half = review.sphere(2 * pairs, 1.)[:pairs]
    half = half / np.linalg.norm(half, axis=1)[:, None]
    return np.r_[half, -half]


def energy_ladder(pairs, removed_pairs, radius=1., mass=1.e12):
    """One antipodal pair unbinds per pass, with a one-K margin either side.

    All rows share the shell radius, so the discrete shell potential is
    K*(n-1) for every survivor.  Setting 0.5*(u+H*R)^2 = K*(n_pass-1)+K makes
    exactly the leading pair marginally unbound at each pass.
    """
    total = 2 * pairs
    directions = antipodal_directions(pairs)
    position = radius * directions + 5.
    hubble_a = 100. * math.sqrt(OM0 / AEXPN ** 3 + 1. - OM0) * AEXPN
    k = (GRAVITY / AEXPN) * float(np.float32(mass)) / radius
    speed = np.zeros(pairs)
    for p in range(1, removed_pairs + 1):
        speed[p - 1] = math.sqrt(2. * k * (total - 2 * p + 2)) - hubble_a * radius
    assert np.all(speed >= 0.)
    velocity = np.r_[speed, speed][:, None] * directions
    return np.c_[position, velocity]


def reference_unbinding(phase, mass):
    """Replica of the production loop; returns its pass count and survivors."""
    data = np.asarray(phase, dtype=np.float32).astype(float)
    offset = data[:, :3] - 5.
    radius = np.linalg.norm(offset, axis=1)
    order = np.argsort(radius, kind='stable')
    offset, radius, raw = offset[order], radius[order], data[order, 3:]
    mass = float(np.float32(mass))
    hubble_a = 100. * math.sqrt(OM0 / AEXPN ** 3 + 1. - OM0) * AEXPN
    keep = np.arange(len(radius))
    passes = 0
    while len(keep):
        passes += 1
        r = radius[keep]
        bulk = raw[keep].mean(axis=0)
        inverse = np.where(r > 0., 1. / np.where(r > 0., r, 1.), 0.)
        suffix = np.r_[np.cumsum(inverse[::-1])[::-1][1:], 0.]
        interior = np.arange(len(keep)) * inverse
        potential = (GRAVITY / AEXPN) * mass * (interior + suffix)
        velocity = raw[keep] - bulk + hubble_a * offset[keep]
        survivors = keep[.5 * np.sum(velocity * velocity, axis=1) - potential <= 0.]
        if len(survivors) == len(keep):
            break
        keep = survivors
    return passes, sorted(order[keep] + 1)


def run_unbinding(work, records):
    measured = []
    for pairs, removed in ((100, 90), (200, 190)):
        phase = energy_ladder(pairs, removed)
        passes, survivors = reference_unbinding(phase, 1.e12)
        assert passes == removed + 1, (pairs, removed, passes)
        for mode in ('checked', 'optimized'):
            result = review.run_halo(work, f'unbinding_ladder_{2 * pairs}', phase,
                                     mode, cell=1., mass=1.e12, identity=True)
            assert result['ids'] == survivors, (mode, len(result['ids']), len(survivors))
            result.update(reference_passes=passes, reference_survivors=len(survivors),
                          population=2 * pairs)
            records.append(result)
        measured.append(dict(population=2 * pairs, unbinding_passes=passes,
                             survivors=len(survivors),
                             production_loop_bound=2 * pairs + 1))
    return measured


# --------------------------------------------------------------------------
# 3. Outermost-crossing absorption of an equal companion
# --------------------------------------------------------------------------

def pair_fixture(count, separation, clump=.15):
    primary = review.sphere(count, clump) + 5.
    companion = review.sphere(count, clump) + [5. + separation, 5., 5.]
    return np.c_[np.r_[primary, companion], np.zeros((2 * count, 3))]


def run_pairs(work, records):
    count, mass = 120, 1.e12
    single = (count * float(np.float32(mass)) / THRESHOLD) ** (1. / 3.)
    merged = (2 * count * float(np.float32(mass)) / THRESHOLD) ** (1. / 3.)
    scan = []
    for ratio in (.6, .9, 1.0, 1.1, 1.2, 1.25, 1.3, 1.5, 2.0):
        separation = ratio * single
        phase = pair_fixture(count, separation)
        expected, enclosed = followup.exact_outer_root(phase, mass, 15.)
        for mode in ('checked', 'optimized'):
            result = review.run_halo(work, f'pair_{ratio:.2f}', phase, mode,
                                     cell=1., mass=mass, identity=True)
            assert math.isclose(result['values'][3], expected, rel_tol=6e-8), \
                (ratio, result['values'][3], expected)
            # Every particle inside the aperture is bound here, so the reported
            # membership must reproduce the oracle's enclosed population. The
            # outermost crossing is therefore self-limiting: its radius is the
            # one the enclosed particles themselves support, never larger.
            assert len(result['ids']) == enclosed, (ratio, len(result['ids']), enclosed)
            result.update(separation=separation, separation_over_single=ratio,
                          expected_so=expected, enclosed=enclosed)
            records.append(result)
        scan.append(dict(separation_over_single_rso=ratio, separation=separation,
                         so_radius=expected, enclosed_particles=enclosed,
                         absorbed_companion=enclosed > count))
    absorbed = [s['separation_over_single_rso'] for s in scan if s['absorbed_companion']]
    return dict(single_clump_rso=single, merged_pair_rso=merged,
                largest_absorbing_separation_over_single_rso=max(absorbed, default=None),
                smallest_isolated_separation_over_single_rso=
                min((s['separation_over_single_rso'] for s in scan
                     if not s['absorbed_companion']), default=None),
                scan=scan)


# --------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--compiler', choices=['gfortran', 'ifx'], default='gfortran')
    parser.add_argument('--output', type=Path, default=HERE / 'halo_gaps_results.json')
    args = parser.parse_args()

    records, summary = [], {}
    with tempfile.TemporaryDirectory(prefix='bdm-gaps-ties-') as directory:
        work = Path(directory)
        ties_build, source_sha, before_sha = ties.build(work, args.compiler)
        summary['duplicate_experiments'] = run_duplicates(work, args.compiler, records)
        print(f"{summary['duplicate_experiments']} duplicate-merge experiments passed", flush=True)

    with tempfile.TemporaryDirectory(prefix='bdm-gaps-halo-') as directory:
        work = Path(directory)
        halo_build = review.build(work)
        summary['unbinding'] = run_unbinding(work, records)
        print('unbinding pass counts:', summary['unbinding'], flush=True)
        summary['companion_absorption'] = run_pairs(work, records)
        print('companion absorption:',
              json.dumps({k: v for k, v in summary['companion_absorption'].items()
                          if k != 'scan'}), flush=True)

    report = dict(checked_at_utc=datetime.now(timezone.utc).isoformat(),
                  compiler=args.compiler, source_sha256=source_sha,
                  host_ties_baseline_sha256=before_sha,
                  compilation=dict(host_ties=ties_build, halo=halo_build),
                  test_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                               for p in [Path(__file__), HERE / 'halo_host_ties.py',
                                         HERE / 'halo_host_ties_cases.f90',
                                         HERE / 'halo_review.py',
                                         HERE / 'halo_review_cases.f90',
                                         HERE / 'halo_followup.py']},
                  summary=summary, passed=len(records), results=records)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    print(f'{len(records)} focused gap experiments passed; wrote {args.output.name}')


if __name__ == '__main__':
    main()
