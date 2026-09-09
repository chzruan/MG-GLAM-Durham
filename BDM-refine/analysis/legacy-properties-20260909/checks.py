"""Independent small controls for periodic pairs and low-mode bias."""
import itertools
import json
import time

import numpy as np

import measure as m


def brute(pos, vel, edges, box):
    radii, radial, tangential = [], [], []
    for i, j in itertools.combinations(range(len(pos)), 2):
        dr = pos[j] - pos[i]
        # Independent minimum-image expression.
        dr = (dr + box/2) % box - box/2
        r = np.linalg.norm(dr)
        dv = vel[j] - vel[i]
        u = np.dot(dv, dr)/r if r else np.nan
        radii.append(r)
        radial.append(u)
        tangential.append((np.dot(dv, dv)-u*u)/2)
    radii, radial, tangential = map(np.asarray, (radii, radial, tangential))
    count = np.histogram(radii, edges)[0]
    moments = np.full((len(count), 5), np.nan)
    for b in range(len(count)):
        use = (radii >= edges[b]) & (radii < edges[b+1])
        u = radial[use]
        if len(u):
            mean = u.mean()
            mu2 = np.mean((u-mean)**2)
            moments[b, :3] = mean, np.sqrt(mu2), np.sqrt(np.mean(tangential[use]))
            if mu2 > 0:
                moments[b, 3:] = np.mean((u-mean)**3)/mu2**1.5, np.mean((u-mean)**4)/mu2**2 - 3
    return count, moments


def check_pairs():
    rng = np.random.default_rng(32467)
    pos = rng.uniform(0, 16, (94, 3))
    pos = np.vstack([pos, [0.01, 1, 1], [15.99, 1, 1], pos[3]])
    vel = rng.normal(0, 270, pos.shape)
    edges = np.geomspace(.01, 6, 12)
    got = m.pair_statistics(pos, vel, edges=edges, box=16.)
    count, moments = brute(pos, vel, edges, 16.)
    np.testing.assert_array_equal(got['count'], count)
    np.testing.assert_allclose(got['moments'], moments, rtol=1e-12, atol=1e-11)
    assert int(got['zero']) == 1
    return {'haloes': len(pos), 'unordered_pairs': int(count.sum()), 'zero_pairs': 1}


def check_invariance():
    rng = np.random.default_rng(74)
    pos = rng.uniform(0, 20, (85, 3))
    vel = rng.normal(0, 300, pos.shape)
    edges = np.geomspace(.07, 8, 10)
    first = m.pair_statistics(pos, vel, edges=edges, box=20.)
    changed = m.pair_statistics((pos+[3.11, 6.25, 9.73]) % 20,
                               vel+[901., -1342., 14.], edges=edges, box=20.)
    for key in ('count', 'xi', 'moments'):
        np.testing.assert_allclose(first[key], changed[key], rtol=1e-11, atol=1e-10)
    return {'periodic_translation_and_Galilean_boost': True}


def check_known_velocities():
    pos, vel = [], []
    speeds = np.array([-300., -100., 200., 400.])
    for i, u in enumerate(speeds):
        pos.extend([[1+5*i, 1, 1], [1.5+5*i, 1, 1]])
        vel.extend([[0., 0., 0.], [u, 60., 80.]])
    got = m.pair_statistics(np.array(pos), np.array(vel), edges=np.array([.1, 1.]), box=32.)
    d = speeds - speeds.mean()
    target = [speeds.mean(), np.sqrt(np.mean(d*d)), np.sqrt(5000.),
              np.mean(d**3)/np.mean(d*d)**1.5, np.mean(d**4)/np.mean(d*d)**2-3]
    np.testing.assert_allclose(got['moments'][0], target, atol=1e-12)
    assert got['count'][0] == 4
    return dict(zip(m.MOMENT_NAMES, map(float, target)))


def check_density_phase_window():
    ng = 64
    vectors = np.array([[1, -2, 1], [0, 1, -1], [2, 0, 0]], dtype=int)
    amplitudes = np.array([.3, .11, .2])
    phases = np.array([.4, -1.3, .9])
    zz, yy, xx = np.ogrid[:ng, :ng, :ng]
    fine = np.zeros((ng, ng, ng), dtype=np.float64)
    for vec, amp, phase in zip(vectors, amplitudes, phases):
        win = np.prod(np.sinc(vec/ng)**2)
        fine += amp*win*np.cos(2*np.pi/ng*(vec[0]*xx+vec[1]*yy+vec[2]*zz)+phase)
    target = .5*amplitudes*np.exp(1j*phases)
    errors = []
    for factor in (2, 4):
        actual = m.matter_modes(m.block_mean(fine, factor), vectors, ng)
        np.testing.assert_allclose(actual, target, rtol=1e-12, atol=1e-13)
        errors.append(float(np.max(np.abs(actual-target))))
    return {'max_coefficient_error': max(errors), 'axis_order_phase_CIC_and_average_windows': True}


def check_halo_fourier():
    rng = np.random.default_rng(9483)
    pos = rng.uniform(0, m.BOX, (45, 3))
    vec = m.wavevectors()
    got = m.halo_modes(pos, vec)
    target = np.exp(-2j*np.pi*(vec @ pos.T)/m.BOX).mean(axis=1)
    np.testing.assert_allclose(got, target, rtol=1e-12, atol=1e-13)
    shift = np.array([17., 23., 47.])
    changed = m.halo_modes((pos+shift) % m.BOX, vec)
    np.testing.assert_allclose(changed, got*np.exp(-2j*np.pi*(vec @ shift)/m.BOX), atol=2e-14)
    # Every member of a +/- pair appears once, no DC or Nyquist modes.
    tuples = set(map(tuple, vec))
    assert (0, 0, 0) not in tuples
    assert all(tuple(-x) not in tuples for x in vec)
    return {'independent_modes': len(vec), 'direct_DFT_max_error': float(np.max(np.abs(got-target)))}


def check_bias():
    rng = np.random.default_rng(387)
    vec = m.wavevectors()
    matter = rng.normal(size=len(vec)) + 1j*rng.normal(size=len(vec))
    halo = (1.7 + .31j)*matter
    table, bands = m.bias_statistics(halo, matter, vec)
    np.testing.assert_allclose(table[:, -1], 1.7, atol=1e-14)
    np.testing.assert_allclose(bands, 1.7, atol=1e-14)
    return {'injected_bias': 1.7, 'recovered_band_bias': bands.tolist()}


def check_selections():
    records = []
    for z in (0, 2):
        cats, _, _ = m.load_pair('F', z)
        selections = m.selections(cats)
        for j, threshold in enumerate(m.THRESHOLDS):
            n = min(len(selections[v][j]) for v in m.VARIANTS)
            for v in m.VARIANTS:
                assert len(selections[v][j+2]) == n
                ids = selections[v][j]
                np.testing.assert_array_equal(ids, np.flatnonzero(cats[v][:, 6] >= 10.**threshold))
                assert len(set(selections[v][j+2])) == n
        records.append({'z': z, 'counts': {v: list(map(len, a)) for v, a in selections.items()}})
    return {'real_selection_controls': records}


def main():
    started = time.monotonic()
    m.set_num_threads(1)
    results = []
    for func in (check_pairs, check_invariance, check_known_velocities,
                 check_density_phase_window, check_halo_fourier, check_bias, check_selections):
        t = time.monotonic()
        detail = func()
        results.append(dict(name=func.__name__, passed=True,
                            elapsed_seconds=time.monotonic()-t, detail=detail))
        print(func.__name__, 'PASS', flush=True)
    m.write_json(m.HERE/'checks.json', dict(passed=True, controls=results,
                                          elapsed_seconds=time.monotonic()-started,
                                          source_hashes=m.source_hashes()))


if __name__ == '__main__':
    main()
