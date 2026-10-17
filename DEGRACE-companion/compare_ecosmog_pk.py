#!/usr/bin/env python3
"""Compare a companion run's PowerDM.R spectra with the ECOSMOG twin's POWMES spectra.

usage: compare_ecosmog_pk.py <run_dir> <node e.g. 002> [--pkdir DIR] [--aout docs/degrace_ecosmog/aout.json]
                             [--box 1024] [--ngrid 2048] [--kmax 1.0]

ECOSMOG files PK_FR_N<node>_S<nn>.dat (POWMES, 1024^3 grid, folding): S01 is the IC
(z = 49), S02..S43 the 42 outputs in ascending a (the pinned aout list); columns
col1 = integer k index (k = col1 2 pi / L), col4 = P / L^3 (docs/degrace_ecosmog/README.md).
MG-GLAM PowerDM.R: 'bin k N Preal', raw CIC on the Ngrid mesh, deconvolved here with the
shell-averaged CIC window 1 - (2/3) sin^2(pi k / 2 k_Ny).  Both runs start from the same
IC files, so at linear k the ratio is a mode-exact test of the IC ingestion, growth and
epochs (expected 1 to < 1 per cent); a LCDM run against the f(R) twin shows the f(R)
enhancement growing at late times (|fR0| = 5.6e-6 for node 2: < 0.5 per cent at z > 3).
"""
import argparse, glob, json, math, os, re
import numpy as np

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run'); ap.add_argument('node')
    ap.add_argument('--pkdir', default='/cosma7/data/dp004/bl267/Runs/DEGRACE/Runs/L1024/FR/power_spectra/Node_%s')
    ap.add_argument('--aout', default='docs/degrace_ecosmog/aout.json')
    ap.add_argument('--box', type=float, default=1024.0); ap.add_argument('--ngrid', type=int, default=2048)
    ap.add_argument('--kmax', type=float, default=1.0); ap.add_argument('--realization', type=int, default=1)
    a = ap.parse_args()
    node = a.node.zfill(3); pkdir = a.pkdir % node if '%s' in a.pkdir else a.pkdir
    aout = sorted(json.load(open(a.aout))['Node_' + node])
    kf = 2 * math.pi / a.box; kny = math.pi * a.ngrid / a.box
    edges = [0.01, 0.03, 0.1, 0.3, 1.0, 3.0]
    edges = [e for e in edges if e <= a.kmax * 1.0001] + ([a.kmax] if a.kmax not in edges else [])
    files = sorted(glob.glob(os.path.join(a.run, 'Run%d' % a.realization, 'PowerDM.R.*.%04d.dat' % a.realization)))
    print('%d PowerDM.R files in %s; ECOSMOG spectra from %s' % (len(files), a.run, pkdir))
    print('   step      a       z   S   ' + '  '.join('%8s' % ('%.2f-%.2f' % (lo, hi)) for lo, hi in zip(edges[:-1], edges[1:])) + '   (P_GLAM/P_ECOSMOG, mean over bins; N bins)')
    for f in files:
        hdr = open(f).read(300).splitlines()[1]
        ag = float(re.search(r'Aexpn *= *([0-9.Ee+-]+)', hdr).group(1)); step = int(re.search(r'Step= *(\d+)', hdr).group(1))
        i = int(np.argmin(np.abs(np.array(aout) - ag)))
        if abs(aout[i] - ag) > 2e-5 * ag:
            print('  %5d %8.5f  no ECOSMOG output at this epoch (closest %.5f)' % (step, ag, aout[i])); continue
        fe = os.path.join(pkdir, 'PK_FR_N%s_S%02d.dat' % (node, i + 2))
        if not os.path.exists(fe):
            print('  %5d %8.5f  missing %s' % (step, ag, fe)); continue
        rows = [l.split() for l in open(f) if len(l.split()) == 4 and l.split()[0].isdigit()]
        kg = np.array([float(r[1]) for r in rows]); pg = np.array([float(r[3]) for r in rows])
        pg = pg / (1 - 2.0 / 3.0 * np.sin(math.pi * kg / (2 * kny)) ** 2)        # CIC deconvolution (shell average)
        E = np.loadtxt(fe); ke = E[:, 0] * kf; pe = E[:, 3] * a.box ** 3
        sel = (ke > 0) & (pe > 0)
        pe_at = np.exp(np.interp(np.log(kg), np.log(ke[sel]), np.log(pe[sel])))
        out = []
        for lo, hi in zip(edges[:-1], edges[1:]):
            s = (kg >= lo) & (kg < hi)
            out.append('%8.4f' % (pg[s] / pe_at[s]).mean() if s.any() else '       -')
        nb = [(kg >= lo).sum() - (kg >= hi).sum() for lo, hi in zip(edges[:-1], edges[1:])]
        print('  %5d %8.5f %6.3f  %02d  ' % (step, ag, 1 / ag - 1, i + 2) + '  '.join(out) + '   ' + ' '.join(str(int(x)) for x in nb))

if __name__ == '__main__':
    main()
