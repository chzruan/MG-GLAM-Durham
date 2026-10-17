#!/usr/bin/env python3
"""Pin the metadata of Baojiu Li's DEGRACE ECOSMOG suite (96 nodes) for
MG-GLAM companion runs.

Reads, per node:
  - the design row (DEGRACE_L512.csv, cosma7) and the CAMB .ini actually run,
  - 2LPTic.param, log_2LPTic.txt and inputspec_ics.txt of the ICs (cosma8 copy),
  - the ECOSMOG namelist (output list, f(R) parameters) and the RAMSES info files
    of the first and last output (cosma8 LCDM runs),
and writes nodes.csv (one row per node) and aout.json (the 42 output expansion
factors per node). It also checks that every derived quantity reproduces the
files on disk and that the 25 secondary outputs sit at chi = 256 k Mpc/h in the
node's own matter+Lambda background.

Usage (cosemu env, ~1 min, login-node light):
  micromamba run -n cosemu python3 -I extract_metadata.py --out-dir .
"""
import argparse, csv, glob, json, os, re
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

BASE8 = '/cosma8/data/dp203/bl267/Projects/Ongoing/DEGRACE'
RUNS = f'{BASE8}/Runs/LCDM'
ICS = f'{BASE8}/ICs/IC_data/L1024'
PYC = '/cosma7/data/dp004/bl267/Runs/DEGRACE/ICs/CAMB/Parameter_files/Python_codes'
DESIGN = f'{PYC}/DEGRACE_L512.csv'
C_H0 = 2997.92458          # Mpc/h
OMBH2, OMNUH2 = 0.02237, 0.0006442
PRIMARY_A = 1.0 / (1.0 + np.arange(0.0, 2.0001, 0.125))
IC_KEYS = ('Nmesh', 'Nsample', 'Box', 'Omega', 'OmegaLambda', 'OmegaBaryon',
           'HubbleParam', 'Redshift', 'Sigma8', 'Seed', 'GlassTileFac',
           'PrimordialIndex', 'WhichSpectrum')


def chi(a, om):
    f = lambda x: 1.0 / (x * x * np.sqrt(om / x**3 + 1.0 - om))  # noqa: E731
    return C_H0 * quad(f, a, 1.0, epsabs=0, epsrel=1e-11)[0]


def ramses_info(outdir):
    f = glob.glob(f'{outdir}/info_*.txt')
    d = {}
    if f:
        for line in open(f[0]):
            if line.strip().startswith('ordering'):
                break
            if '=' in line:
                k, v = line.split('=', 1)
                d[k.strip()] = v.strip()
    return d


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out-dir', default='.')
    args = ap.parse_args()
    design = [[float(x) for x in r] for r in csv.reader(open(DESIGN))]
    rows, aouts, failures = [], {}, []
    for n in range(1, 97):
        node = f'Node_{n:03d}'
        om, h, s8_, ns, mg = design[n - 1][:5]
        ob = OMBH2 / h**2
        sig8_target = s8_ / np.sqrt(om / 0.3)
        fr0 = -10.0**mg
        c1oc22 = (-fr0) * (3.0 * (1.0 + 4.0 * (1.0 - om) / om))**2
        dens_v_on = 80.0 if -fr0 < 1e-6 else (120.0 if -fr0 < 1e-5 else 160.0)
        r = dict(node=node, Omega_m=om, h=h, S_8=s8_, n_s=ns, log10_fR0=mg,
                 Omega_b=round(ob, 7), sigma8_target=round(sig8_target, 7),
                 fR0=fr0, fr_c1oc22_expected=round(c1oc22, 12), dens_V_on_expected=dens_v_on)
        # --- ICs
        ic = {}
        for line in open(f'{ICS}/{node}/2LPTic.param'):
            line = line.split('%')[0].strip()
            if line:
                k, *v = line.split()
                if k in IC_KEYS:
                    ic[k] = v[0]
        m = re.search(r'Sigma8 = ([\d.eE+-]+)\s*\n\s*Normalization adjusted to\s+Sigma8=([\d.eE+-]+)\s+\(Normfac=([\d.eE+-]+)\)',
                      open(f'{ICS}/{node}/log_2LPTic.txt').read())
        first = open(f'{ICS}/{node}/inputspec_ics.txt').readline().split()
        r.update(sigma8_ic=float(ic['Sigma8']), Omega_b_ic=float(ic['OmegaBaryon']),
                 seed=int(ic['Seed']), z_init=float(ic['Redshift']), Np=int(ic['Nsample']),
                 Nmesh=int(ic['Nmesh']), Box=float(ic['Box']), Dplus_z49=float(first[1]),
                 normfac_2lptic=float(m.group(3)), ic_nfiles=len(glob.glob(f'{ICS}/{node}/ics.*')))
        # --- CAMB input actually used
        ini = open(f'{PYC}/params_n{n:03d}.ini').read()
        g = lambda k: float(re.search(r'^\s*' + re.escape(k) + r'\s*=\s*([-+0-9.eE]+)', ini, re.M).group(1))  # noqa: E731
        r.update(omch2_camb=g('omch2'), ombh2_camb=g('ombh2'), omnuh2_camb=g('omnuh2'),
                 H0_camb=g('hubble'), n_s_camb=g('scalar_spectral_index(1)'), A_s=g('scalar_amp(1)'))
        r['Omega_m_camb'] = round((r['ombh2_camb'] + r['omch2_camb'] + r['omnuh2_camb']) / h**2, 5)
        # --- ECOSMOG namelist and run status (LCDM twin on cosma8)
        t = open(f'{RUNS}/{node}/namelist/cosmo.nml').read()
        aout = [float(x) for x in re.search(r'aout\s*=\s*([^\n/]+)', t).group(1).split(',') if x.strip()]
        r['fr_c1oc22'] = float(re.search(r'fr_c1oc22\s*=\s*([^\n/,]+)', t).group(1))
        r['dens_V_on'] = float(re.search(r'dens_V_on\s*=\s*([^\n/,]+)', t).group(1))
        outs = sorted(glob.glob(f'{RUNS}/{node}/data/output_*'))
        i1, il = ramses_info(outs[0]), ramses_info(outs[-1])
        r.update(n_outputs=len(outs), n_gadget_outputs=len(glob.glob(f'{RUNS}/{node}/gadget_data/output_*')),
                 has_gadget_data_old=os.path.isdir(f'{RUNS}/{node}/gadget_data_old'),
                 first_aexp=float(i1['aexp']), last_aexp=float(il['aexp']),
                 ramses_H0=float(i1['H0']), ramses_omega_m=float(i1['omega_m']))
        # --- output list: 17 primary + 25 stations at chi = 256 k
        isprim = np.array([np.any(np.isclose(a, PRIMARY_A, rtol=0, atol=1e-9)) for a in aout])
        sec = sorted([a for a, p in zip(aout, isprim) if not p], reverse=True)
        pred = [brentq(lambda x: chi(x, om) - 256.0 * k, 0.005, 0.9999999) for k in range(1, len(sec) + 1)]
        r.update(n_aout=len(aout), n_primary=int(isprim.sum()), n_secondary=len(sec),
                 station_max_reldiff=float(max(abs(p - s) / s for p, s in zip(pred, sec))))
        a_sorted = np.sort(aout)
        r['n_close_pairs_lt_0.5pct'] = int((((a_sorted[1:] - a_sorted[:-1]) / a_sorted[1:]) < 0.005).sum())
        aouts[node] = aout
        # --- consistency checks
        checks = {
            'Omega_m(2LPTic)': abs(float(ic['Omega']) - om) < 1e-6,
            'h(2LPTic)': abs(float(ic['HubbleParam']) - h) < 1e-6,
            'Omega_b(2LPTic)=0.02237/h^2': abs(r['Omega_b_ic'] - ob) < 1e-6,
            'sigma8(2LPTic) within 1e-4 of S8*sqrt(0.3/Om)': abs(r['sigma8_ic'] - sig8_target) < 1.01e-4,
            'sigma8(2LPTic)==2LPTic log target': abs(float(m.group(2)) - r['sigma8_ic']) < 1e-6,
            'fr_c1oc22 formula': abs(r['fr_c1oc22'] - c1oc22) < 1e-9,
            'dens_V_on rule': r['dens_V_on'] == dens_v_on,
            'CAMB omch2 == Omega_m*h^2 (sic)': abs(r['omch2_camb'] - om * h * h) < 1e-5,
            'CAMB H0, n_s': abs(r['H0_camb'] - 100 * h) < 0.01 and abs(r['n_s_camb'] - ns) < 1e-4,
            'CAMB ombh2, omnuh2 fixed': r['ombh2_camb'] == OMBH2 and abs(r['omnuh2_camb'] - OMNUH2) < 1e-9,
            'RAMSES header cosmology': abs(r['ramses_H0'] / 100 - h) < 1e-6 and abs(r['ramses_omega_m'] - om) < 1e-6,
            'seed 2026, z_init 49, 1024^3, Box 1024': (r['seed'], r['z_init'], r['Np'], r['Box']) == (2026, 49.0, 1024, 1024.0),
            '43 outputs, a from 0.02 to 1': r['n_outputs'] == 43 and r['first_aexp'] == 0.02 and r['last_aexp'] == 1.0,
            '42 aout = 17 primary + 25 stations': (r['n_aout'], r['n_primary'], r['n_secondary']) == (42, 17, 25),
            'stations at chi=256k to 1e-6': r['station_max_reldiff'] < 1e-6,
        }
        for k, ok in checks.items():
            if not ok:
                failures.append((node, k))
        rows.append(r)
    os.makedirs(args.out_dir, exist_ok=True)
    with open(os.path.join(args.out_dir, 'nodes.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    json.dump(aouts, open(os.path.join(args.out_dir, 'aout.json'), 'w'), indent=0)
    print(f'{len(rows)} nodes written to {args.out_dir}/nodes.csv and aout.json')
    print('consistency failures:', failures if failures else 'none')
    dm = np.array([r['Omega_m_camb'] - r['Omega_m'] for r in rows])
    print(f'CAMB total Omega_m exceeds the simulated Omega_m by {dm.min():.4f} .. {dm.max():.4f}')
    print('close pairs (<0.5% in a) per node:', sorted(set(r['n_close_pairs_lt_0.5pct'] for r in rows)))


if __name__ == '__main__':
    main()
