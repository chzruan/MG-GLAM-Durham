#!/usr/bin/env python3
"""Verify a DEGRACE companion run directory (pilot or production box).

usage: check_run.py <run_dir> [--node Node_002] [--aout docs/degrace_ecosmog/aout.json]
                    [--planes all|sample] [--job JOBID]

Checks: the 42-output schedule (exact flags, classes 17 + 25, epochs equal to
float32(1/(1+z)) of Init.dat and to the pinned ECOSMOG list), PowerDM and
LensPlanes files for every output (sizes, meta epochs, plane means = 1),
BDM catalogues at the class-1 outputs, the ic2pm velocity epoch, the
timing.log breakdown (per-step cost, output overhead) and, with sacct,
elapsed / MaxRSS / CPU efficiency.  Exit 0 = every check passed.
"""
import argparse, glob, json, os, re, struct, subprocess, sys
import numpy as np

def f32(x): return struct.unpack('f', struct.pack('f', x))[0]

def read_meta(f):
    m = {}
    for line in open(f):
        if '=' in line:
            k, v = line.split('=', 1); m[k.strip()] = v.strip()
    return m

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run'); ap.add_argument('--node', default=None); ap.add_argument('--aout', default=None)
    ap.add_argument('--planes', default='sample', choices=['all', 'sample']); ap.add_argument('--job', default=None)
    ap.add_argument('--box', type=int, default=1)
    a = ap.parse_args()
    run = a.run.rstrip('/'); R1 = os.path.join(run, 'Run%d' % a.box)
    fails, warns = [], []
    def check(ok, msg):
        print(('PASS ' if ok else 'FAIL ') + msg)
        if not ok: fails.append(msg)
    def warn(msg): print('WARN ' + msg); warns.append(msg)

    # ---- Init.dat redshift list
    init = open(os.path.join(run, 'Init.dat')).read().splitlines()
    i0 = [i for i, l in enumerate(init) if l.startswith('#outputs')][0]
    nout = int(init[i0].split('=')[1].split()[0])
    zs = []
    for l in init[i0 + 1:]:
        try: zs += [float(v) for v in l.split()]
        except ValueError: break
    zs = zs[:abs(nout)]
    check(nout < 0 and len(zs) == abs(nout), 'Init.dat: #outputs = %d with %d redshifts (exact mode)' % (nout, len(zs)))
    a_req = sorted(f32(1.0 / (1.0 + z)) for z in zs)

    # ---- OutputSchedule.dat
    fs = os.path.join(R1, 'OutputSchedule.dat')
    check(os.path.exists(fs), 'OutputSchedule.dat present')
    sched = []
    if os.path.exists(fs):
        for l in open(fs):
            if l.startswith('#'): continue
            p = l.split(); sched.append(dict(step=int(p[0]), a=float(p[1]), z=float(p[2]), chi=float(p[3]), exact=int(p[4]), cls=int(p[5])))
    steps = [s['step'] for s in sched]
    check(len(sched) == abs(nout), 'schedule has %d outputs (expected %d)' % (len(sched), abs(nout)))
    check(steps == sorted(steps) and len(set(steps)) == len(steps), 'schedule steps strictly increasing')
    check(all(s['exact'] == 1 for s in sched), 'all outputs flagged exact')
    a_sch = [s['a'] for s in sched]
    check(len(a_sch) == len(a_req) and all(abs(x - y) <= 3e-7 * y for x, y in zip(a_sch, a_req)),
          'schedule epochs equal float32(1/(1+z)) of the Init.dat list within 3 ulp (max rel diff %.1e)' %
          (max([abs(x - y) / y for x, y in zip(a_sch, a_req)] or [0])))
    ncls = {c: sum(1 for s in sched if s['cls'] == c) for c in sorted(set(s['cls'] for s in sched))}
    check(ncls.get(1, 0) == 17 and ncls.get(2, 0) == 25, 'classes: %s (expected 17 primary + 25 stations)' % ncls)
    if a.aout and a.node:
        ao = sorted(json.load(open(a.aout))[a.node])
        d = max(abs(x - y) / y for x, y in zip(a_sch, ao)) if len(ao) == len(a_sch) else 1
        check(len(ao) == len(a_sch) and d < 3e-7, 'schedule epochs match the pinned ECOSMOG %s list (%d values, max rel diff %.1e)' % (a.node, len(ao), d))
    # primary redshifts sit on the 0.125 grid
    prim = [s for s in sched if s['cls'] == 1]
    check(all(abs(round(s['z'] * 8) / 8 - s['z']) < 2e-6 for s in prim), 'class-1 outputs are the 17 primary redshifts (multiples of 0.125, 0..2)')

    # ---- per-output products
    npow = nlens = nmeta = ncat = ncatbin = 0; bad = []
    sizes = set(); plane_means = []
    sample = set([steps[0], steps[len(steps) // 2], steps[-1]]) if steps else set()
    for s in sched:
        st = s['step']; tag = '%04d.%04d' % (st, a.box)
        fp = os.path.join(R1, 'PowerDM.R.%s.dat' % tag)
        if os.path.exists(fp):
            npow += 1
            hdr = open(fp).read(400).splitlines()[1]
            ah = float(re.search(r'Aexpn *= *([0-9.Ee+-]+)', hdr).group(1))
            if abs(ah - s['a']) > 1e-5: bad.append('PowerDM %s epoch %g vs schedule %g' % (tag, ah, s['a']))
        fb = os.path.join(R1, 'LensPlanes.%s.bin' % tag); fm = fb[:-4] + '.meta'
        if os.path.exists(fb) and os.path.exists(fm):
            nlens += 1; m = read_meta(fm)
            ns, nl = int(m['N_slab']), int(m['N_lens'])
            exp = 3 * ns * nl * nl * 4; sizes.add(os.path.getsize(fb))
            if os.path.getsize(fb) != exp: bad.append('LensPlanes %s size %d != %d' % (tag, os.path.getsize(fb), exp))
            if abs(float(m['a']) - s['a']) > 1e-6 * s['a']: bad.append('LensPlanes %s meta a %s vs schedule %g' % (tag, m['a'], s['a']))
            if a.planes == 'all' or st in sample:
                pl = np.fromfile(fb, dtype='<f4').reshape(3 * ns, nl, nl)
                means = pl.reshape(3, ns, -1).mean(axis=(1, 2))        # per axis over the slabs
                plane_means.append((st, means, pl.min(), pl.max()))
                if np.abs(means - 1).max() > 1e-5 or pl.min() < 0: bad.append('LensPlanes %s axis means %s min %g' % (tag, means, pl.min()))
        fc = os.path.join(R1, 'CATALOGS', 'CatshortV.%s.DAT' % tag)      # text and/or binary (BDM.config CatalogueFormat)
        fb = os.path.join(R1, 'CATALOGS', 'CatshortV.%s.bin' % tag)
        fcm = os.path.join(R1, 'CATALOGS', 'CatshortV.%s.meta' % tag)
        has_text, has_bin = os.path.exists(fc), os.path.exists(fb)
        if has_bin != os.path.exists(fcm): bad.append('catalogue %s: .bin and .meta must come together' % tag)
        if has_text or has_bin:
            ncat += 1
            if s['cls'] != 1: bad.append('catalogue at a station output %s' % tag)
            if has_bin and os.path.exists(fcm):
                ncatbin += 1
                m = read_meta(fcm)
                rows = int(m.get('rows', '-1').split()[0])
                if os.path.getsize(fb) != 96 * rows: bad.append('catalogue %s: .bin size %d != 96 x rows %d' % (tag, os.path.getsize(fb), rows))
                if abs(float(m.get('a', '0').split()[0]) - s['a']) > 1e-6 * s['a']: bad.append('catalogue %s meta a %s vs schedule %g' % (tag, m.get('a'), s['a']))
                tc = (m.get('text_catalogue', '') or '').split()
                if (bool(tc) and tc[0] != 'none') != has_text: bad.append('catalogue %s: text_catalogue = %s but .DAT %s' % (tag, tc[:1], 'present' if has_text else 'absent'))
                if has_text:
                    with open(fc, 'rb') as fh: nl = sum(buf.count(b'\n') for buf in iter(lambda: fh.read(1 << 22), b''))
                    if nl - 8 != rows: bad.append('catalogue %s: text rows %d != binary rows %d' % (tag, nl - 8, rows))
        elif s['cls'] == 1:
            bad.append('missing catalogue at primary output %s (z = %.4f)' % (tag, s['z']))
    n = len(sched)
    check(npow == n, 'PowerDM.R files: %d of %d, epochs match' % (npow, n))
    check(nlens == n and len(sizes) == 1, 'LensPlanes bin+meta: %d of %d, one size (%s bytes)' % (nlens, n, sizes))
    check(ncat == ncls.get(1, 0), 'catalogues at primary outputs: %d of %d (%d with a binary .bin + .meta, checked: size = 96 x rows, meta a, text_catalogue, row counts)' % (ncat, ncls.get(1, 0), ncatbin))
    stale = glob.glob(os.path.join(R1, 'CATALOGS', '.Catshort*'))
    check(not stale, 'no staged catalogue temporaries (.Catshort*) left in CATALOGS/' if not stale else '%d staged catalogue temporaries (.Catshort*) in CATALOGS/: a finder died between staging and publication; delete them when no finder is running' % len(stale))
    for st, means, mn, mx in plane_means:
        print('     planes step %4d: axis means %s  min %.3g  max %.3g' % (st, np.array2string(means, precision=7), mn, mx))
    check(not bad, 'per-output consistency (%d problems)' % len(bad))
    for b in bad[:12]: print('     ' + b)

    # ---- ic2pm
    fl = os.path.join(R1, 'ic2pm.log')
    if os.path.exists(fl):
        t = open(fl).read()
        ep = re.search(r'velocity epoch *= *(\w+)', t)
        check(ep is not None and ep.group(1) == 'half', 'ic2pm velocity epoch = %s' % (ep.group(1) if ep else '?'))
        for key in ('rms at a_init', 'rms written'):
            l = [x for x in t.splitlines() if key in x]
            if l: print('     ' + l[0].strip())
    # ---- timing.log
    ft = os.path.join(R1, 'timing.log')
    if os.path.exists(ft):
        rows = []
        for l in open(ft):
            p = l.split()
            if len(p) >= 10 and p[0].isdigit(): rows.append([float(v) for v in p])
        T = np.array(rows)
        if len(T):
            step, tot = T[:, 0].astype(int), T[:, 1]
            cols = ['Force', 'Move', 'Density', 'IO', 'Bias', 'BDM', 'Analysis', 'MG']
            isout = np.isin(step, steps)
            base = np.median(tot[~isout]) if (~isout).any() else np.nan
            over = (tot[isout] - base).sum()
            print('timing.log: %d steps, total %.1f min, median non-output step %.3f min; output steps %d: extra %.1f min (%.1f%% of total)'
                  % (len(T), tot.sum(), base, isout.sum(), over, 100 * over / tot.sum()))
            print('     column totals [min]: ' + ', '.join('%s %.1f' % (c, T[:, 2 + i].sum()) for i, c in enumerate(cols)))
            pr = np.isin(step, [s['step'] for s in prim]); stn = isout & ~pr
            if pr.any(): print('     primary outputs: mean step %.2f min (BDM %.2f, Analysis %.2f, IO %.2f)' % (tot[pr].mean(), T[pr, 7].mean(), T[pr, 8].mean(), T[pr, 5].mean()))
            if stn.any(): print('     station outputs: mean step %.2f min (Analysis %.2f, IO %.2f)' % (tot[stn].mean(), T[stn, 8].mean(), T[stn, 5].mean()))
            check(len(T) >= steps[-1] if steps else True, 'timing.log covers the last scheduled step (%d of %d)' % (len(T), steps[-1] if steps else 0))
    # ---- sacct
    job = a.job
    if job is None:
        outs = glob.glob(os.path.join(run, 'slurm_*.out'))
        if outs: job = re.search(r'slurm_(\d+)\.out', sorted(outs)[-1]).group(1)
    if job:
        try:
            o = subprocess.run(['sacct', '-j', job, '-X', '-P', '-n', '-o', 'JobID,State,Elapsed,AllocCPUS,ReqMem,MaxRSS,TotalCPU,NodeList'],
                               capture_output=True, text=True, timeout=60).stdout.strip()
            o2 = subprocess.run(['sacct', '-j', job + '.batch', '-P', '-n', '-o', 'MaxRSS,TotalCPU'],
                                capture_output=True, text=True, timeout=60).stdout.strip()
            print('sacct %s: %s | batch step: %s' % (job, o, o2))
            p = o.split('|')
            def secs(t):
                d, _, hms = t.rpartition('-'); h, m, s = hms.split(':'); return (int(d) * 86400 if d else 0) + int(h) * 3600 + int(m) * 60 + float(s)
            if p[1] == 'COMPLETED':
                el, ncpu = secs(p[2]), int(p[3]); tc = secs(o2.split('|')[1]) if o2 else 0
                print('     elapsed %.2f h, %d cores, CPU efficiency %.0f%%, MaxRSS %s' % (el / 3600, ncpu, 100 * tc / (el * ncpu), o2.split('|')[0] if o2 else '?'))
            else:
                warn('job %s state %s (not COMPLETED yet)' % (job, p[1]))
        except Exception as e:
            warn('sacct unavailable: %s' % e)
    print('\n%d FAIL, %d WARN' % (len(fails), len(warns)))
    sys.exit(1 if fails else 0)

if __name__ == '__main__':
    main()
