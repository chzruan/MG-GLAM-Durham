"""Same-field legacy/v3 halo statistics; never modifies campaign inputs.

Run with micromamba run -n cosemu python3 -B measure.py --help.
ASCII columns are deliberately used for BOTH finders (the membership index
has a different column layout). Velocities are peculiar, in km/s.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import resource
import struct
import subprocess
import time

import numpy as np
from numba import njit, prange, set_num_threads
from scipy.fft import rfftn
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
CAMPAIGN = REPO / 'BDM-refine/validation/20260908-convergence'
BOX = 256.
VARIANTS = ('legacy', 'v3')
MASS_EDGES = np.arange(12.5, 15.50001, .25)
THRESHOLDS = (12.5, 13.)
SELECTIONS = ('mass12p5', 'mass13', 'rank12p5', 'rank13')
R_EDGES = np.geomspace(.01, 50., 29)
K_EDGES = np.array([.024, .05, .075, .1, .125, .15, .2])
BIAS_BANDS = np.array([[.05, .15], [.025, .1], [.1, .2]])
MOMENT_NAMES = ('mean_vr', 'sigma_r', 'sigma_t_1d', 'skew_r', 'excess_kurt_r')


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        while b := f.read(16 << 20):
            h.update(b)
    return h.hexdigest()


def write_json(path, value):
    tmp = Path(str(path) + '.partial')
    tmp.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    tmp.replace(path)


def source_hashes():
    return {str(p.relative_to(REPO)): sha(p) for p in
            [HERE/'measure.py', HERE/'checks.py', REPO/'PMP2linker.f90',
             REPO/'PMP2mod_density.f90']}


def load_pair(case, z):
    p = CAMPAIGN / f'replay-{case}-z{z}-ng2048.json'
    receipt = json.loads(p.read_text())
    path = Path(receipt['catalogue_arrays'])
    expected = receipt['catalogue_arrays_sha256']
    if sha(path) != expected:
        raise ValueError(f'Catalogue hash mismatch: {path}')
    with np.load(path, allow_pickle=False) as a:
        cats = {v: a[v].copy() for v in VARIANTS}
    for v, a in cats.items():
        if a.ndim != 2 or a.shape[1] != 24 or not np.isfinite(a).all():
            raise ValueError(f'Invalid catalogue {case} z{z} {v}')
        if (a[:, 6] <= 0).any() or (a[:, :3] < 0).any() or (a[:, :3] > BOX).any():
            raise ValueError('Mass/coordinate range')
        if len(a) != receipt['variants'][v]['catalogue']['rows']:
            raise ValueError('Row count mismatch')
        a[:, :3] %= BOX  # ASCII can round a periodic edge to L.
    prov = dict(receipt=str(p.relative_to(REPO)), receipt_sha256=sha(p),
                catalogue_arrays=str(path.relative_to(REPO)),
                catalogue_arrays_sha256=expected,
                rows={v: len(a) for v, a in cats.items()})
    return cats, receipt, prov


def selections(cats):
    result = {v: [] for v in VARIANTS}
    ranks = {v: np.argsort(-a[:, 6], kind='stable') for v, a in cats.items()}
    for threshold in THRESHOLDS:
        masks = {v: np.flatnonzero(a[:, 6] >= 10.**threshold) for v, a in cats.items()}
        for v in VARIANTS:
            result[v].append(masks[v])
    for it, threshold in enumerate(THRESHOLDS):
        n = min(len(result[v][it]) for v in VARIANTS)
        for v in VARIANTS:
            result[v].append(ranks[v][:n])
    return result


def catalogue_statistics():
    out = dict(mass_edges=MASS_EDGES, variants=np.array(VARIANTS),
               cases=np.array(list('ABCDEFT')), redshifts=np.array([0, 1, 2]),
               profile_names=np.array(['bulk_sigma_1d', 'median_internal_vrms',
                                       'median_vmax', 'median_speed']))
    provenance = []
    for case in 'ABCDEFT':
        for z in (0, 1, 2):
            cats, _, prov = load_pair(case, z)
            provenance.append(prov)
            for v, a in cats.items():
                key = f'{case}_z{z}_{v}'
                logm = np.log10(a[:, 6])
                bins = np.searchsorted(MASS_EDGES, logm, side='right') - 1
                count = np.histogram(logm, MASS_EDGES)[0]
                octant = (a[:, :3] >= BOX/2).astype(int) @ np.array([1, 2, 4])
                octcounts = np.array([np.histogram(logm[octant == o], MASS_EDGES)[0]
                                      for o in range(8)])
                profiles = np.full((len(count), 4), np.nan)
                for b in range(len(count)):
                    t = a[bins == b]
                    if len(t):
                        vel = t[:, 3:6]
                        profiles[b] = [np.sqrt(np.mean((vel - vel.mean(axis=0))**2)),
                                       np.median(t[:, 9]), np.median(t[:, 10]),
                                       np.median(np.linalg.norm(vel, axis=1))]
                out[key+'_counts'] = count
                out[key+'_hmf'] = count / BOX**3 / np.diff(MASS_EDGES)
                out[key+'_octant_counts'] = octcounts
                out[key+'_profiles'] = profiles
                out[key+'_total_rows'] = np.array(len(a))
    path = HERE/'catalogue-statistics.npz'
    np.savez_compressed(path, **out)
    write_json(HERE/'catalogue-statistics.json', dict(
        completed=True, source_hashes=source_hashes(), inputs=provenance,
        output_sha256=sha(path), box_mpc_h=BOX,
        selection='Published catalogues, bound mass, no additional duplicate mask',
        velocity='Stored peculiar km/s; same staggered output convention in each pair'))


@njit(cache=False)
def _pair_accumulate(pos, vel, pairs, edges, box):
    nb = len(edges) - 1
    count = np.zeros(nb, np.int64)
    sums = np.zeros((nb, 7))  # r, vr, central2/3/4, transverse norm/2, raw vr2
    under = 0
    zero = 0
    for q in range(len(pairs)):
        i, j = pairs[q]
        dr = pos[j] - pos[i]
        dr -= box * np.rint(dr/box)
        r = np.sqrt(np.sum(dr*dr))
        if r == 0:
            zero += 1
        if r < edges[0]:
            under += 1
        b = np.searchsorted(edges, r, side='right') - 1
        if b < 0 or b >= nb or r == 0:
            continue
        dv = vel[j] - vel[i]
        vr = np.sum(dv*dr)/r
        count[b] += 1
        sums[b, 0] += r
        sums[b, 1] += vr
        sums[b, 5] += max(0., np.sum(dv*dv) - vr*vr)/2
        sums[b, 6] += vr*vr
    means = np.zeros(nb)
    for b in range(nb):
        if count[b]:
            means[b] = sums[b, 1]/count[b]
    # Second pass centres moments before powers; avoids raw-moment cancellation.
    for q in range(len(pairs)):
        i, j = pairs[q]
        dr = pos[j] - pos[i]
        dr -= box * np.rint(dr/box)
        r = np.sqrt(np.sum(dr*dr))
        b = np.searchsorted(edges, r, side='right') - 1
        if b < 0 or b >= nb or r == 0:
            continue
        d = np.sum((vel[j]-vel[i])*dr)/r - means[b]
        sums[b, 2] += d*d
        sums[b, 3] += d*d*d
        sums[b, 4] += d*d*d*d
    return count, sums, under, zero


def pair_statistics(pos, vel, edges=R_EDGES, box=BOX, threads=1):
    from pycorr import TwoPointCorrelationFunction
    pos = np.ascontiguousarray(pos % box, dtype=np.float64)
    vel = np.ascontiguousarray(vel, dtype=np.float64)
    corr = TwoPointCorrelationFunction(
        's', edges, data_positions1=pos, position_type='pos',
        boxsize=box, engine='corrfunc', nthreads=threads, estimator='natural')
    pairs = cKDTree(pos, boxsize=box).query_pairs(edges[-1], output_type='ndarray')
    count, sums, under, zero = _pair_accumulate(pos, vel, pairs, edges, box)
    ordered = np.asarray(corr.D1D2.ncounts, dtype=np.int64)
    if not np.array_equal(2*count, ordered):
        raise ValueError(f'Independent pair-count mismatch: {2*count} vs {ordered}')
    n = len(pos)
    rr = n*(n-1)*4*np.pi/3*np.diff(edges**3)/box**3
    xi = np.divide(ordered, rr, out=np.full(len(rr), np.nan), where=rr>0) - 1
    if not np.allclose(xi, corr.corr, rtol=1e-11, atol=1e-10, equal_nan=True):
        raise ValueError('Analytical periodic RR normalization mismatch')
    avg = np.divide(sums, count[:, None], out=np.full_like(sums, np.nan),
                    where=count[:, None]>0)
    moments = np.full((len(count), 5), np.nan)
    moments[:, 0] = avg[:, 1]
    moments[:, 1] = np.sqrt(avg[:, 2])
    moments[:, 2] = np.sqrt(avg[:, 5])
    nonzero = avg[:, 2] > 0
    moments[nonzero, 3] = avg[nonzero, 3]/avg[nonzero, 2]**1.5
    moments[nonzero, 4] = avg[nonzero, 4]/avg[nonzero, 2]**2 - 3
    return dict(count=count, xi=xi, r_mean=avg[:, 0], moments=moments,
                rr_ordered=rr, underflow=np.array(under), zero=np.array(zero))


def wavevectors(box=BOX, kmax=K_EDGES[-1]):
    nmax = int(np.floor(kmax*box/(2*np.pi)))
    a = np.arange(-nmax, nmax+1)
    v = np.array(np.meshgrid(a, a, a, indexing='ij')).reshape(3, -1).T
    # Exactly one member of each +/- pair, with nonnegative rFFT x index.
    keep = (v[:, 0] > 0) | ((v[:, 0] == 0) & (v[:, 1] > 0)) | (
        (v[:, 0] == 0) & (v[:, 1] == 0) & (v[:, 2] > 0))
    k = np.linalg.norm(v, axis=1)*2*np.pi/box
    return v[keep & (k < kmax)]


def block_mean(field, factor):
    n = field.shape[0]//factor
    return field.reshape(n, factor, n, factor, n, factor).mean(
        axis=(1, 3, 5), dtype=np.float64)


def matter_modes(field, vectors, fine_ng, box=BOX, threads=1, fine_cic=True):
    """Low-k physical Fourier coefficients from [z,y,x] block averages.

    Correct exact block-centre phase and discrete box-average window. The
    source is a CIC contrast field; its fine-grid CIC window is also removed.
    Halo modes below use a direct sum, with no halo assignment window.
    """
    ng = field.shape[0]
    field -= np.mean(field, dtype=np.float64)
    ft = rfftn(field, workers=threads, norm='forward', overwrite_x=True)
    modes = ft[vectors[:, 2] % ng, vectors[:, 1] % ng, vectors[:, 0]].copy()
    factor = fine_ng//ng
    offset = (factor-1)/2 * box/fine_ng
    phase = np.exp(-2j*np.pi*np.sum(vectors, axis=1)*offset/box)
    average_window = np.prod(np.sinc(vectors/ng)/np.sinc(vectors/fine_ng), axis=1)
    cic_window = np.prod(np.sinc(vectors/fine_ng)**2, axis=1) if fine_cic else 1.
    return modes * phase / average_window / cic_window


def stream_density(receipt, coarse=512, threads=1):
    d = receipt['density']
    path = Path(d['path'])
    if not d['shared_bit_identical_density']:
        raise ValueError('Different old/new density fields')
    hasher = hashlib.sha256()
    started = time.monotonic()
    with path.open('rb') as f:
        header = f.read(28)
        hasher.update(header)
        ng, npart, step, a, box = struct.unpack('>qqiff', header)
        if ng != 2048 or ng % coarse or box != BOX or step != receipt['step']:
            raise ValueError('Density metadata mismatch')
        factor = ng//coarse
        grid = np.empty((coarse, coarse, coarse), dtype=np.float64)
        slab_bytes = factor*ng*ng*4
        for iz in range(coarse):
            raw = f.read(slab_bytes)
            if len(raw) != slab_bytes:
                raise ValueError('Truncated density tape')
            hasher.update(raw)
            slab = np.frombuffer(raw, dtype='>f4').reshape(factor, coarse, factor, coarse, factor)
            grid[iz] = slab.mean(axis=(0, 2, 4), dtype=np.float64)
        if f.read(1):
            raise ValueError('Trailing density bytes')
    actual = hasher.hexdigest()
    if actual != d['sha256']:
        raise ValueError('Density SHA256 mismatch')
    if not np.isfinite(grid).all():
        raise ValueError('Nonfinite matter grid')
    vectors = wavevectors()
    lowgrid = block_mean(grid, 2)
    mean = float(grid.mean())
    low = matter_modes(lowgrid, vectors, ng, threads=threads)
    del lowgrid
    modes = matter_modes(grid, vectors, ng, threads=threads)
    return vectors, modes, low, dict(
        density_path=str(path.relative_to(REPO)), density_sha256=actual,
        fine_ng=ng, coarse_ng=coarse, comparison_ng=coarse//2,
        nparticles=npart, step=step, a=a, box=box, mean_delta=mean,
        elapsed_seconds=time.monotonic()-started)


@njit(parallel=True, cache=False)
def halo_modes(pos, vectors, box=BOX):
    """Direct point-sample Fourier transform, one independent +/- mode."""
    out = np.empty(len(vectors), np.complex128)
    for m in prange(len(vectors)):
        re, im = 0., 0.
        for i in range(len(pos)):
            phase = (2*np.pi/box)*(pos[i, 0]*vectors[m, 0] +
                                   pos[i, 1]*vectors[m, 1] +
                                   pos[i, 2]*vectors[m, 2])
            re += np.cos(phase)
            im -= np.sin(phase)
        out[m] = (re + 1j*im)/len(pos)
    return out


def bias_statistics(h, m, vectors, box=BOX):
    k = np.linalg.norm(vectors, axis=1)*2*np.pi/box
    cross = np.real(h*np.conj(m))
    auto = np.abs(m)**2
    vals = []
    for lo, hi in zip(K_EDGES[:-1], K_EDGES[1:]):
        use = (k >= lo) & (k < hi)
        vals.append([k[use].mean(), use.sum(), box**3*auto[use].mean(),
                     box**3*cross[use].mean(), cross[use].sum()/auto[use].sum()])
    bands = []
    for lo, hi in BIAS_BANDS:
        use = (k >= lo) & (k < hi)
        bands.append(cross[use].sum()/auto[use].sum())
    return np.array(vals), np.array(bands)


def measure_epoch(z, threads):
    start = time.monotonic()
    before = source_hashes()
    cats, receipt, provenance = load_pair('F', z)
    masks = selections(cats)
    vectors, matter, matter_low, density = stream_density(receipt, threads=threads)
    out = dict(r_edges=R_EDGES, k_edges=K_EDGES, bias_bands=BIAS_BANDS,
               wavevectors=vectors, matter_modes=matter, matter_modes_ng256=matter_low,
               selections=np.array(SELECTIONS), moment_names=np.array(MOMENT_NAMES))
    for v in VARIANTS:
        for s, ids in zip(SELECTIONS, masks[v]):
            a = cats[v][ids]
            key = v+'_'+s
            print(f'z{z} {key}: {len(a)} haloes', flush=True)
            out[key+'_rows'] = ids
            out[key+'_minimum_mass'] = np.min(a[:, 6])
            stat = pair_statistics(a[:, :3], a[:, 3:6], threads=threads)
            for name, value in stat.items():
                out[key+'_'+name] = value
            h = halo_modes(a[:, :3], vectors)
            out[key+'_halo_modes'] = h
            out[key+'_bias_k'], out[key+'_bias'] = bias_statistics(h, matter, vectors)
            out[key+'_bias_k_ng256'], out[key+'_bias_ng256'] = bias_statistics(h, matter_low, vectors)
    if source_hashes() != before:
        raise RuntimeError('Analysis sources changed while job was running')
    output = HERE/f'F-z{z}-statistics.npz'
    np.savez_compressed(output, **out)
    packages = {p: importlib.metadata.version(p) for p in
                ('numpy', 'scipy', 'numba', 'pycorr', 'Corrfunc')}
    write_json(HERE/f'F-z{z}-statistics.json', dict(
        completed=True, input=provenance, density=density, source_hashes=before,
        output_sha256=sha(output), packages=packages, threads=threads,
        job_id=os.environ.get('SLURM_JOB_ID'), hostname=os.uname().nodename,
        git_head=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=REPO, text=True).strip(),
        elapsed_seconds=time.monotonic()-start,
        maxrss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        cpu_seconds=resource.getrusage(resource.RUSAGE_SELF).ru_utime +
                    resource.getrusage(resource.RUSAGE_SELF).ru_stime,
        pair_count_crosscheck='All 8 selections agree exactly: 2*KDTree unordered = Corrfunc ordered',
        bias='sum Re(delta_h delta_m*) / sum |delta_m|^2; independent +/- modes; no shot subtraction',
        uncertainties='Single periodic realization; no independent-volume uncertainty estimate'))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--catalogues', action='store_true')
    p.add_argument('--redshift', type=int, choices=[0, 1, 2])
    p.add_argument('--threads', type=int, default=1)
    args = p.parse_args()
    set_num_threads(args.threads)
    if args.catalogues:
        catalogue_statistics()
    if args.redshift is not None:
        measure_epoch(args.redshift, args.threads)
