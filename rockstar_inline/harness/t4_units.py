"""Test T4 of the in-line Rockstar: positions and velocities in the right units and at the right epoch.
usage: python3 t4_units.py RUN_DIR STEP INLINE_BASE [--npmin N] [--box-index B]
  RUN_DIR      the run directory with PMcrd.<step>.DAT (+ PMcrs*) and CATALOGS/CatshortV.<step>.<box>
1. Particle level (the test): for hosts with Np >= npmin, the mean velocity of every snapshot particle within
   the halo radius r, converted independently by glamio (km/s with velocity_scale('a_v'), positions in Mpc/h),
   against the halo's bulk velocity (bulkvel). The median of the per-halo ratio q = v_bulk . v_mean / |v_mean|^2
   (halos with |v_mean| > 200 km/s) must lie within tol of 1 and farther than tol from a_v / AEXPN, the ratio a
   conversion with AEXPN would give; tol = max(3 errors, 20 % of |1 - a_v/AEXPN|). The floor allows for the
   small difference between Rockstar's particles (bulkvel) and a sphere of all particles.
2. BDM (information): BDM converts with AEXPN, so the median of the same ratio against BDM's halo velocity of the
   nearest distinct BDM halo (Mtot > 1e13, within 0.3 Rvir) is expected near AEXPN / a_v.
3. Descriptive: the slope of log vrms on log Mvir (bound) and the median T/|U| by mass.
Part of the MG-GLAM in-line Rockstar adapter; GNU GPL v3 or later."""
import argparse
import os
import sys

import numpy as np
from scipy.spatial import cKDTree

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
sys.path.insert(0, os.path.join(here, "..", "..", "python"))
import glamio  # noqa: E402
from glamio.snapshot import read_snapshot  # noqa: E402
from rsbin import read_rsbin  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("run")
ap.add_argument("step", type=int)
ap.add_argument("inline_base")
ap.add_argument("--npmin", type=int, default=1000)
ap.add_argument("--box-index", type=int, default=0)
a = ap.parse_args()

cat, meta = read_rsbin(a.inline_base)
aexpn, a_v, box = float(meta["a"]), float(meta["a_v"]), float(meta["Box"])
snap = read_snapshot(os.path.join(a.run, f"PMcrd.{a.step:04d}.DAT"))
h = snap.header
if abs(h.aexpn - aexpn) > 1e-6 or abs(h.a_v - a_v) > 1e-6:
    sys.exit(f"the snapshot (a = {h.aexpn}, a_v = {h.a_v}) is not the catalogue's output (a = {aexpn}, a_v = {a_v})")
print(f"step {a.step}: a = {aexpn:.6f}, a_v = {a_v:.6f}; a wrong epoch (AEXPN) would give q = a_v/a = {a_v / aexpn:.5f}")


def median_ratio(vh, vref, vmin):
    """median of v_h . v_ref / |v_ref|^2 over rows with |v_ref| > vmin, and its error (1.253 sigma / sqrt n)"""
    n2 = np.sum(vref ** 2, axis=1)
    ok = n2 > vmin ** 2
    q = np.sum(vh[ok] * vref[ok], axis=1) / n2[ok]
    lo, hi = np.percentile(q, [16, 84])
    return np.median(q), 1.253 * 0.5 * (hi - lo) / np.sqrt(ok.sum()), ok.sum()


# 1. particle level
x, y, z, vx, vy, vz = snap.read(units="physical", velocity_epoch="a_v", max_particles=None)
pos = np.c_[x, y, z]
del x, y, z
pos[pos >= box] -= box
tree = cKDTree(pos, boxsize=box)
host = (cat["pid"] == -1) & (cat["num_p"] >= a.npmin)
hc = np.mod(np.c_[cat["x"], cat["y"], cat["z"]][host].astype(float), box)
hc[hc >= box] -= box
vmean = np.zeros((host.sum(), 3))
for k, members in enumerate(tree.query_ball_point(hc, cat["r"][host].astype(float) * 1e-3)):
    m = np.asarray(members)
    vmean[k] = vx[m].mean(), vy[m].mean(), vz[m].mean()
vbulk = np.c_[cat["bulkvel_x"], cat["bulkvel_y"], cat["bulkvel_z"]][host].astype(float)
q, err, nq = median_ratio(vbulk, vmean, 200.0)
tol = max(3 * err, 0.2 * abs(1 - a_v / aexpn))
ok1 = abs(q - 1) <= tol and abs(q - a_v / aexpn) > tol
print(f"1. hosts with Np >= {a.npmin}: {host.sum()} ({nq} with |v_mean| > 200 km/s); median v_bulk/v_mean "
      f"q = {q:.5f} +- {err:.5f} ({(q - 1) / err:+.1f} sigma from 1, {(q - a_v / aexpn) / err:+.1f} sigma from the "
      f"wrong epoch; tolerance {tol:.5f}): {'ok' if ok1 else 'FAIL'}")
del tree, pos, vx, vy, vz

# 2. BDM, robust
bdm = glamio.read_catalogue(os.path.join(a.run, "CATALOGS", f"CatshortV.{a.step:04d}.{a.box_index:04d}"))
sel = (bdm["distinct"] == 0) & (bdm["Mtot"] > 1e13)
pb = np.mod(np.c_[bdm["x"], bdm["y"], bdm["z"]][sel], box)
pb[pb >= box] -= box
vb = np.c_[bdm["vx"], bdm["vy"], bdm["vz"]][sel]
hall = cat["pid"] == -1
pi = np.mod(np.c_[cat["x"], cat["y"], cat["z"]][hall].astype(float), box)
pi[pi >= box] -= box
d, j = cKDTree(pi, boxsize=box).query(pb)
near = d < 0.3 * bdm["Rvir"][sel] * 1e-3
vi = np.c_[cat["bulkvel_x"], cat["bulkvel_y"], cat["bulkvel_z"]][hall][j[near]].astype(float)
q2, err2, n2 = median_ratio(vi, vb[near], 300.0)
print(f"2. BDM distinct Mtot > 1e13: {sel.sum()}, matched {near.sum()} ({n2} with |v| > 300 km/s): median "
      f"v_bulk/v_BDM {q2:.5f} +- {err2:.5f}; expected AEXPN/a_v = {aexpn / a_v:.5f} (BDM converts with AEXPN)")

# 3. descriptive
hh = hall & (cat["num_p"] >= 100) & (cat["mgrav"] > 0) & (cat["vrms"] > 0)
p = np.polyfit(np.log10(cat["mgrav"][hh]), np.log10(cat["vrms"][hh]), 1)
print(f"3. log vrms = {p[0]:.4f} log Mvir + {p[1]:.3f} over {hh.sum()} hosts")
for lo, hi in [(1e12, 3e12), (3e12, 1e13), (1e13, 1e14), (1e14, 1e16)]:
    b = hh & (cat["mgrav"] >= lo) & (cat["mgrav"] < hi)
    if b.sum():
        print(f"   T/|U| median of hosts [{lo:.0e},{hi:.0e}): {np.median(cat['kin_to_pot'][b]):.4f} (n {b.sum()})")
print("T4", "PASS" if ok1 else "FAIL")
sys.exit(0 if ok1 else 1)
