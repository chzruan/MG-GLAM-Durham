"""Test T1 of the in-line Rockstar: the in-line catalogue against standalone MPI-Rockstar on the same particles.
usage: python3 t1_compare.py INLINE_BASE REF_DIR BOX [--strict] [--shift S] [--npmin N]
  INLINE_BASE  .../Rockstar.<step>.<box> (without .bin/.meta), or a second reference directory (its
               out_0.parents.list), to measure the standalone-vs-standalone spread the thresholds come from
  REF_DIR      directory with the reference out_0.parents.list (find_parents output)
  --strict     the reference ran with STRICT_SO_MASSES = 1: its Mvir_all and R are the strict ones, which the
               in-line m_strict / r_strict (and pid) match; otherwise compare m and all halos only
  --shift S    the in-line positions were shifted by S Mpc/h (test T3)
Thresholds (calibrated on the standalone-vs-standalone spread, docs/rockstar-inline-audit-20261026.md 6.3):
  host mass-function ratios within Poisson (2 sigma) and < 0.5 % in bins with > 1000 hosts... reported as
  |ratio - 1| <= max(2 sigma, 0.005); matched fraction >= 0.93 (hosts with Np >= 100 within 0.1 R200m);
  |dln M| < 0.02 for >= 90 % of the matches.
Part of the MG-GLAM in-line Rockstar adapter; GNU GPL v3 or later."""
import argparse
import os
import sys

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rsbin import read_rsbin  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("inline_base")
ap.add_argument("ref_dir")
ap.add_argument("box", type=float)
ap.add_argument("--strict", action="store_true")
ap.add_argument("--shift", type=float, default=0.0)
ap.add_argument("--npmin", type=int, default=100)
a = ap.parse_args()
box = a.box
RHO_M_PER_OM = 2.77519737e11

def load_ref(d, shift=0.0):
    # out_*.list + find_parents' PID: 2 Mvir (= mgrav, BOUND_PROPS = 1), 7 Np, 8-10 x y z,
    # 19 Mvir_all (= m, the strict m with STRICT_SO_MASSES), 20 M200b, last PID
    r = np.loadtxt(os.path.join(d, "out_0.parents.list"), comments="#")
    return {"M": r[:, 2], "Mall": r[:, 19], "Np": r[:, 7], "pos": np.mod(r[:, 8:11] - shift, box), "M200b": r[:, 20],
            "host": r[:, -1] == -1, "n": len(r)}


R = load_ref(a.ref_dir)
if os.path.isdir(a.inline_base):
    I = load_ref(a.inline_base, a.shift)
    om = 0.3089
    print(f"second reference {a.inline_base}: {I['n']} haloes, {I['host'].sum()} hosts")
else:
    cat, meta = read_rsbin(a.inline_base)
    om = float(meta["Om"])
    I = {"M": cat["mgrav"].astype(float), "Mall": (cat["m_strict"] if a.strict else cat["m"]).astype(float),
         "Np": cat["num_p"].astype(float), "pos": np.mod(np.c_[cat["x"], cat["y"], cat["z"]] - a.shift, box),
         "M200b": cat["alt_m_1"].astype(float), "host": cat["pid"] == -1, "n": len(cat)}
    print(f"in-line {a.inline_base}: {I['n']} haloes, {I['host'].sum()} hosts (pid radius {meta.get('pid_radius')})")
ref = np.empty(R["n"])
print(f"reference {a.ref_dir}: {len(ref)} haloes, {R['host'].sum()} hosts")
ok_all = True


def mf(label, mi, mr, edges):
    global ok_all
    ni, _ = np.histogram(mi, edges)
    nr, _ = np.histogram(mr, edges)
    print(f"  {label}: bin, reference, in-line, ratio, allowed")
    for k in range(len(edges) - 1):
        if nr[k] == 0:
            continue
        r = ni[k] / nr[k]
        allow = max(2 * np.sqrt(1 / nr[k] + 1 / max(ni[k], 1)), 0.005)
        flag = "" if abs(r - 1) <= allow else "  <-- outside"
        ok_all &= abs(r - 1) <= allow
        print(f"    [{edges[k]:.0e},{edges[k+1]:.0e}) {nr[k]:8d} {ni[k]:8d}  {r:.4f}  +-{allow:.4f}{flag}")


edges = np.array([2e11, 5e11, 1e12, 3e12, 1e13, 3e13, 1e14, 1e16])
mf("all haloes, bound Mvir", I["M"], R["M"], edges)
if a.strict:
    mf("hosts, bound Mvir", I["M"][I["host"]], R["M"][R["host"]], edges)
    mf("hosts, strict Mvir_all", I["Mall"][I["host"]], R["Mall"][R["host"]], edges)
else:
    mf("all haloes, Mvir_all", I["Mall"], R["Mall"], edges)

# host-by-host (hosts with Np >= npmin of the reference matched to the nearest in-line host within 0.1 R200m)
sel = R["host"] & (R["Np"] >= a.npmin) if a.strict else (R["Np"] >= a.npmin)
pool = I["host"] if a.strict else np.ones(I["n"], bool)
r200 = (3 * R["M200b"][sel] / (4 * np.pi * 200 * om * RHO_M_PER_OM)) ** (1 / 3)
tree = cKDTree(I["pos"][pool], boxsize=box)
d, j = tree.query(R["pos"][sel], k=1)
okm = d < 0.1 * r200
dln = np.abs(np.log(I["M"][pool][j[okm]] / R["M"][sel][okm]))
frac, f02 = okm.mean(), np.mean(dln < 0.02)
print(f"  {'hosts' if a.strict else 'haloes'} with Np >= {a.npmin}: {sel.sum()}; matched within 0.1 R200m "
      f"{frac:.4f} (>= 0.93); |dln Mvir| < 0.02 for {f02:.4f} (>= 0.90), < 0.05 for {np.mean(dln < 0.05):.4f}; "
      f"median {np.median(dln):.4f}")
ok_all &= frac >= 0.93 and f02 >= 0.90
print("T1", "PASS" if ok_all else "FAIL")
sys.exit(0 if ok_all else 1)
