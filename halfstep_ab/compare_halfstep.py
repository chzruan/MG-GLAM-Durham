#!/usr/bin/env python3
"""Half-step velocity A/B analysis (see 2LPTIC_Gui/HANDOFF-ic2pm-halfstep.md section 5.2).
Run: micromamba run -n cosemu python3 compare_halfstep.py
T1  R(z)=P_sync/P_half at equal da, low k.   pred 1+0.8*eps*[1-(a0/a)^2.5]
T2  G(z)=[P(z)/P(step1)]/[D(a)/D(a1)]^2, low k.  pred half: 1.000+-0.5%; sync: excess
T3  z=0 growth-corrected P(da8)/P(da4), P(da16)/P(da4), half and sync; plus new da8_sync vs old conv_da8.
"""
import glob, re, os, sys, numpy as np
from scipy.integrate import quad
ROOT = "/cosma8/data/dp203/dc-ruan1/mgglam_claude/MG-GLAM"
AB = f"{ROOT}/halfstep_ab"
Om, OmL = 0.3, 0.7
A0 = 0.02
KLIN = 0.05                      # linear-test k limit [h/Mpc]
Hr = lambda a: np.sqrt(Om/a**3 + OmL)
Du = lambda a: 2.5*Om*Hr(a)*quad(lambda x: 1.0/(x*Hr(x))**3, 1e-8, a)[0]   # from plot_conv.py
D0 = Du(1.0)
DA = {4: 4e-4, 8: 8e-4, 16: 1.6e-3}

def eps(da):                      # first-kick momentum excess of synchronous ICs
    return (A0**1.5 - (A0-da/2)**1.5)/(A0+da/2)**1.5
def R_pred(a, da):                # T1 prediction
    return 1 + 0.8*eps(da)*(1 - (A0/a)**2.5)
def G_sync_pred(a, a1, da):       # T2 sync prediction: excess growth from step-1 output to a
    return R_pred(a, da)/R_pred(a1, da)

def read_pk(f):
    l = open(f).readlines()[1]
    a = float(re.search(r'Aexpn\s*=\s*([-0-9.]+)', l).group(1))
    z = float(re.search(r'Redshift=\s*([-0-9.]+)', l).group(1))
    st = int(re.search(r'Step=\s*(\d+)', l).group(1))
    d = np.loadtxt(f, skiprows=3, usecols=(1, 2, 3))
    return dict(a=a, z=z, step=st, k=d[:, 0], n=d[:, 1], P=d[:, 2], file=os.path.basename(f))

def load_run(d):
    fs = sorted(glob.glob(f"{d}/Run1/PowerDM.R.*.0001.dat"))
    out = [read_pk(f) for f in fs]
    return sorted(out, key=lambda r: r['step'])

def lowk_ratio(r1, r2, kmax=KLIN):
    m = (r1['k'] <= kmax)
    w = r1['n'][m]
    rat = r1['P'][m]/r2['P'][m]
    return np.sum(w*rat)/np.sum(w), rat.min(), rat.max()

runs = {}
for n in (4, 8, 16):
    for e in ('sync', 'half'):
        name = f"da{n}_{e}"
        if os.path.isdir(f"{AB}/{name}/Run1"):
            rr = load_run(f"{AB}/{name}")
            if rr: runs[name] = rr
print("runs found:", {k: (len(v), f"z {v[0]['z']:.2f}..{v[-1]['z']:.3f}") for k, v in runs.items()})

# ---------------- T1 ----------------
print("\n=== T1: R(z) = P_sync/P_half at same da, k <=", KLIN, "(mode-weighted mean; min..max over bins)")
for n in (4, 8, 16):
    s, h = runs.get(f"da{n}_sync"), runs.get(f"da{n}_half")
    if not (s and h): continue
    da = DA[n]
    print(f"--- da={da:.1e}  eps={eps(da):.4f}  late-time pred R={1+0.8*eps(da):.4f}")
    print("  step     z      R_meas   R_pred   (R_meas-1)/(R_pred-1)   [min..max over bins]")
    for rs, rh in zip(s, h):
        assert rs['step'] == rh['step'] and abs(rs['a']-rh['a']) < 1e-6, (rs['file'], rh['file'])
        if not (rs['z'] > 1.5 or rs['step'] == 1): continue
        R, lo, hi = lowk_ratio(rs, rh)
        Rp = R_pred(rs['a'], da)
        print(f"  {rs['step']:4d}  {rs['z']:6.2f}   {R:.5f}  {Rp:.5f}   {(R-1)/(Rp-1):7.3f}          [{lo:.5f}..{hi:.5f}]")

# ---------------- T2 ----------------
print("\n=== T2: G(z) = [P(z)/P(step1)] / [D(a)/D(a1)]^2 at k <=", KLIN)
print("   pred: half = 1.000 +- 0.5%;  sync = G_sync_pred (handoff table: +1.10/+2.06/+3.61% for da 4/8/16)")
summary = {}
for name, rr in runs.items():
    n = int(re.match(r'da(\d+)', name).group(1)); da = DA[n]
    r1 = rr[0]; assert r1['step'] == 1, r1['file']
    D1 = Du(r1['a'])
    print(f"--- {name}  (step-1 output a1={r1['a']:.5f} z1={r1['z']:.3f})")
    print("  step     z       G_meas    G_pred_sync   [min..max over bins]")
    vals = []
    for r in rr[1:]:
        if not (1.5 <= r['z'] <= 12): continue
        g = (Du(r['a'])/D1)**2
        G, lo, hi = lowk_ratio(r, r1)
        G /= g; lo /= g; hi /= g
        print(f"  {r['step']:4d}  {r['z']:6.2f}   {G:.5f}    {G_sync_pred(r['a'], r1['a'], da):.5f}      [{lo:.5f}..{hi:.5f}]")
        if 2 <= r['z'] <= 10: vals.append(G)
    if vals:
        summary[name] = (np.mean(vals), np.min(vals), np.max(vals))
print("\n  T2 summary over 2<=z<=10: name  mean(G)  min  max   |  pred excess (sync)")
for name, (m, lo, hi) in summary.items():
    n = int(re.match(r'da(\d+)', name).group(1)); da = DA[n]
    print(f"   {name:10s}  {m:.5f}  {lo:.5f}  {hi:.5f}   |  sync pred late-time G = {1+0.8*eps(da):.5f}/R(a1)")

# Contamination-free T2 (added 2026-09-17 after review 2): G_sync/G_half removes the nonlinear
# drift common to both epochs; compare with the discrete-leapfrog model (REVIEW-1 App. A):
# +1.118..1.141% / +2.126..2.171% / +3.840..3.929% over 2<=z<=10 for da 4/8/16.
print("\n  T2' G_sync/G_half - 1 [%] over 2<=z<=10 (min..max), k <=", KLIN)
for n in (4, 8, 16):
    s_, h_ = runs.get(f"da{n}_sync"), runs.get(f"da{n}_half")
    if not (s_ and h_): continue
    r1s, r1h = s_[0], h_[0]
    vals = []
    for rs, rh in zip(s_[1:], h_[1:]):
        if not (2 <= rs['z'] <= 10): continue
        gs = lowk_ratio(rs, r1s)[0]; gh = lowk_ratio(rh, r1h)[0]
        vals.append(100*(gs/gh - 1))
    if vals: print(f"   da{n}: {min(vals):+.3f} .. {max(vals):+.3f}")

# ---------------- T3 ----------------
print("\n=== T3: z=0 growth-corrected P ratios vs da4 (same epoch)")
print("   NOTE (2026-09-17, review 2): the linear-D^2 correction of the da16 run's last output (a=1.0102)")
print("   to a=1 is invalid at nonlinear k; see the matched-epoch band table (T3') below for the real error.")
def z0(rr):
    r = min(rr, key=lambda r: abs(r['z']))
    return r if abs(r['z']) < 0.05 else None
kk = [0.02, 0.05, 0.1, 0.2, 0.4, 0.8, 1.5, 2.5]
for e in ('sync', 'half'):
    ref = runs.get(f"da4_{e}")
    if not ref: continue
    ref0 = z0(ref)
    if ref0 is None: continue
    Pr = ref0['P']*(D0/Du(ref0['a']))**2
    print(f"--- epoch {e}: ratios to da4_{e} ({ref0['file']}, z={ref0['z']:.4f})")
    print("   k[h/Mpc]: " + " ".join(f"{k:7.2f}" for k in kk))
    for n in (8, 16):
        rr = runs.get(f"da{n}_{e}")
        if not rr: continue
        r0 = z0(rr)
        if r0 is None: print(f"   da{n}_{e}: no z=0 output"); continue
        P = r0['P']*(D0/Du(r0['a']))**2
        rat = P/np.interp(r0['k'], ref0['k'], Pr)
        print(f"   da{n}_{e}: " + " ".join(f"{rat[np.argmin(abs(r0['k']-k))]:7.4f}" for k in kk) + f"   (z={r0['z']:.4f})")

# ---------------- new da8_sync vs old conv_da8 ----------------
print("\n=== Binary check: new da8_sync vs old conv_da8/Run1 (July binaries), same step outputs")
old = load_run(f"{ROOT}/conv_da8")
new = runs.get("da8_sync")
if new and old:
    om = {r['step']: r for r in old}
    print("  step     z     P_new/P_old(k<=0.05)  [min..max]   P_new/P_old(k~1)")
    for r in new:
        if r['step'] in om:
            o = om[r['step']]
            m, lo, hi = lowk_ratio(r, o)
            i1 = np.argmin(abs(r['k']-1.0))
            print(f"  {r['step']:4d}  {r['z']:6.2f}   {m:.5f}   [{lo:.5f}..{hi:.5f}]   {r['P'][i1]/o['P'][i1]:.5f}")

# ---------------- T3' matched epochs (added 2026-09-17 after review 2) ----------------
# Outputs da4 step 221 / da8 step 111 / da16 step 56 lie within |dln a| <= 0.0016 of each other
# (a ~ 0.847), so P ratios there need only a tiny growth correction and are valid at nonlinear k.
MATCH = {4: 221, 8: 111, 16: 56}
BANDS = [(0, 0.05), (0.05, 0.1), (0.1, 0.2), (0.2, 0.5), (0.5, 1.0), (1.0, 2.5)]
def band_means(r, ref):
    g = (Du(ref['a'])/Du(r['a']))**2          # linear correction over the tiny |dln a|
    rat = r['P']*g/np.interp(r['k'], ref['k'], ref['P'])
    return [100*(np.sum(r['n'][m]*rat[m])/np.sum(r['n'][m]) - 1)
            for m in [(r['k'] > lo) & (r['k'] <= hi) for lo, hi in BANDS]]
print("\n=== T3': matched-epoch (a~0.847) band means of P(daX)/P(da4) - 1 [%], N-weighted")
print("   bands:", BANDS)
for e in ('sync', 'half'):
    ref = runs.get(f"da4_{e}")
    if not ref: continue
    ref0 = next((r for r in ref if r['step'] == MATCH[4]), None)
    if ref0 is None: continue
    for n in (8, 16):
        rr = runs.get(f"da{n}_{e}")
        if not rr: continue
        r0 = next((r for r in rr if r['step'] == MATCH[n]), None)
        if r0 is None: continue
        bm = band_means(r0, ref0)
        print(f"   {e:4s} da{n:<2d} (a={r0['a']:.5f} vs {ref0['a']:.5f}): " + " ".join(f"{x:+6.2f}" for x in bm))
