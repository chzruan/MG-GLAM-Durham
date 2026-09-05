"""Periodic catalogue diagnostics; all indices refer to original data rows."""
from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

COLUMNS = ['x', 'y', 'z', 'vx', 'vy', 'vz', 'Mbound', 'Mtot', 'Rvir',
           'Vrms', 'Vcirc', 'Nhalo', 'Cvir', 'Nparticles', 'DistinctSub',
           'Xoff', '2KEp_m1', 'Lambda', 'RadRMS_k', 'b_a', 'c_a',
           'Major_x', 'Major_y', 'Major_z']
FIELDS = ['x', 'y', 'z', 'vx', 'vy', 'vz', 'Mbound', 'Mtot', 'Rvir',
          'Nhalo', 'Nparticles']
RULE = dict(name='strict-v1', separation_mpc_h=0.2, relative_speed_km_s=5.0,
            abs_dlog10_mtot=0.005, same_mbound=True, same_nparticles=True,
            geometry='periodic minimum image', grouping='connected components',
            survivor='lowest original zero-based data-row index',
            membership_validation='unavailable; catalogue-level criterion')


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(8 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def read_ascii(path):
    with open(path, 'rb') as f:
        header = [f.readline() for _ in range(8)]
    if b'Nhalo' not in header[7]:
        raise ValueError(f'Not an eight-line Catshort header: {path}')
    arr = np.loadtxt(path, skiprows=8, usecols=[COLUMNS.index(k) for k in FIELDS],
                     ndmin=2)
    data = {k: arr[:, i] for i, k in enumerate(FIELDS)}
    if not all(np.isfinite(data[k]).all() for k in FIELDS):
        raise ValueError(f'Non-finite defining fields: {path}')
    if np.any(data['Mtot'] <= 0) or np.any(data['Mbound'] <= 0):
        raise ValueError(f'Nonpositive mass: {path}')
    ids = data['Nhalo']
    if np.any(ids != np.rint(ids)) or np.any(np.diff(ids) <= 0):
        raise ValueError('Nhalo must be strictly increasing integer-valued IDs')
    return data, header


def pairs_within(data, radius, box=1024.0, rows=None):
    # Promote before wrapping/subtraction; HDF5 storage may be float32 or integer.
    xyz = np.column_stack([np.asarray(data[k], dtype=np.float64) for k in ['x', 'y', 'z']])
    if rows is not None:
        xyz = xyz[rows]
    xyz %= box
    pairs = cKDTree(xyz, boxsize=box).query_pairs(radius, output_type='ndarray')
    dr = xyz[pairs[:, 1]] - xyz[pairs[:, 0]]
    dr -= box * np.rint(dr / box)
    distance = np.linalg.norm(dr, axis=1)
    use = distance < radius  # query_pairs includes the boundary; our rule does not.
    if rows is not None:
        pairs = rows[pairs]
    return pairs[use], distance[use], dr[use]


def edge_flags(data, pairs, distance):
    i, j = pairs.T
    same = ((data['Mbound'][i] == data['Mbound'][j]) &
            (data['Nparticles'][i] == data['Nparticles'][j]))
    # Equality uses the original values; only diagnostic arithmetic is promoted.
    speed = np.sqrt(sum((np.asarray(data[k][i], dtype=np.float64) -
                         np.asarray(data[k][j], dtype=np.float64))**2
                        for k in ['vx', 'vy', 'vz']))
    dm = np.abs(np.log10(np.asarray(data['Mtot'][i], dtype=np.float64)) -
                np.log10(np.asarray(data['Mtot'][j], dtype=np.float64)))
    strict = same & (distance < 0.2) & (speed < 5) & (dm < 0.005)
    exact = same.copy()
    for k in ['x', 'y', 'z', 'vx', 'vy', 'vz', 'Mtot', 'Rvir']:
        exact &= data[k][i] == data[k][j]
    return dict(same=same, strict=strict, exact=exact, speed=speed, dlogmass=dm)


def components(nrows, edges):
    """Union by smallest source index, independent of pair enumeration order."""
    parent = np.arange(nrows, dtype=np.int64)

    def root(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for a, b in edges:
        a, b = root(a), root(b)
        if a != b:
            parent[max(a, b)] = min(a, b)
    involved = np.unique(edges)
    for i in involved:
        parent[i] = root(i)
    drop = parent != np.arange(nrows)
    return parent, drop


def mass_table(data, drop):
    # 12.4--14.4 in 0.2 dex steps, then the requested partial bin to 14.5.
    edges = np.r_[np.round(np.arange(12.4, 14.41, 0.2), 8), 14.5]
    mass = np.log10(np.asarray(data['Mtot'], dtype=np.float64))
    total = np.histogram(mass, edges)[0]
    removed = np.histogram(mass[drop], edges)[0]
    return [dict(logmass_lo=float(a), logmass_hi=float(b), rows=int(n),
                 removed=int(k), fraction=float(k / n) if n else None)
            for a, b, n, k in zip(edges[:-1], edges[1:], total, removed)]


def velocity_statistics(data, drop, box=1024.0):
    """Unweighted distinct halo pairs, peculiar radial relative velocities."""
    rows = np.flatnonzero(np.log10(np.asarray(data['Mtot'], dtype=np.float64)) >= 12.4)
    pairs, distance, dr = pairs_within(data, 2.0, box, rows)
    use = distance >= 0.5
    pairs, distance, dr = pairs[use], distance[use], dr[use]
    dv = np.column_stack([np.asarray(data[k][pairs[:, 1]], dtype=np.float64) -
                          np.asarray(data[k][pairs[:, 0]], dtype=np.float64)
                          for k in ['vx', 'vy', 'vz']])
    vr = np.einsum('ij,ij->i', dv, dr) / distance
    kept = ~drop[pairs].any(axis=1)
    result = []
    for lo, hi in [(0.5, 1.0), (1.0, 1.5), (1.5, 2.0), (0.5, 2.0)]:
        item = dict(r_lo=lo, r_hi=hi)
        for name, mask in [('raw', np.ones(len(pairs), bool)), ('clean', kept)]:
            v = vr[(distance >= lo) & (distance < hi) & mask]
            item[name] = dict(pairs=len(v), mean_km_s=float(v.mean()) if len(v) else None,
                              sigma_km_s=float(v.std()) if len(v) else None,
                              rms_km_s=float(np.sqrt(np.mean(v*v))) if len(v) else None)
        result.append(item)
    return result


def audit(data, box=1024.0, velocities=False):
    n = len(data['x'])
    pairs, distance, _ = pairs_within(data, 0.5, box)
    flags = edge_flags(data, pairs, distance)
    parent, drop = components(n, pairs[flags['strict']])
    _, exact_drop = components(n, pairs[flags['exact']])
    selected = np.log10(np.asarray(data['Mtot'], dtype=np.float64)) >= 12.4
    pair_selected = selected[pairs].all(axis=1)
    _, selected_drop = components(n, pairs[flags['strict'] & pair_selected])
    _, selected_exact = components(n, pairs[flags['exact'] & pair_selected])
    counts = []
    for radius in [0.001, 0.05, 0.2, 0.5]:
        for scope, mask in [('all', np.ones(len(pairs), bool)), ('logMtot>=12.4', pair_selected)]:
            use = (distance < radius) & mask
            npair = int(use.sum())
            nsame = int((use & flags['same']).sum())
            counts.append(dict(radius_mpc_h=radius, scope=scope, pairs=npair,
                               same_bound_mass_count=nsame,
                               fraction_same=nsame / npair if npair else None,
                               strict_edges=int((use & flags['strict']).sum())))
    remaining = (distance < 0.2) & flags['same'] & ~drop[pairs].any(axis=1)
    result = dict(rows=n, selected_rows=int(selected.sum()), removed=int(drop.sum()),
                  exact_removed=int(exact_drop.sum()),
                  selected_strict_removed=int(selected_drop.sum()),
                  selected_exact_removed=int(selected_exact.sum()),
                  selected_removed_global_mask=int((drop & selected).sum()),
                  pair_counts=counts,
                  remaining_same_bound_pairs_below_0p2=int(remaining.sum()),
                  remaining_strict_pairs=int((remaining & flags['strict']).sum()),
                  remaining_ambiguous_pairs=[dict(row_i=int(i), row_j=int(j),
                        Nhalo_i=int(data['Nhalo'][i]), Nhalo_j=int(data['Nhalo'][j]),
                        separation=float(d), speed=float(s), dlogmass=float(m))
                      for (i, j), d, s, m in zip(pairs[remaining], distance[remaining],
                          flags['speed'][remaining], flags['dlogmass'][remaining])],
                  mass_bins=mass_table(data, drop))
    if velocities:
        result['pairwise_radial_velocity'] = velocity_statistics(data, drop, box)
    return result, parent, drop, pairs[flags['strict']], exact_drop
