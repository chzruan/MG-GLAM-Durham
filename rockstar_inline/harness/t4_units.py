"""Test T4 of the in-line Rockstar: units, against BDM's independent conversion (PMP2linker.f90).
usage: python3 t4_units.py RUN_DIR STEP INLINE_BASE BOX
1. Bulk velocities of BDM distinct haloes (Mtot > 1e13) and the nearest in-line host (within 0.3 Rvir): BDM converts
   with a = AEXPN, the in-line finder with a_v = AEXPN - ASTEP/2, so the least-squares slope should be AEXPN / a_v.
2. The slope of log vrms against log Mvir for hosts with Np >= 100 (virial scaling: about 1/3).
3. Median T/|U| of hosts by mass (the same quantity as standalone Rockstar's; a velocity unit error by a factor of
   a would move it by a^2).
Part of the MG-GLAM in-line Rockstar adapter; GNU GPL v3 or later."""
import os
import sys

import numpy as np
from scipy.spatial import cKDTree

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
sys.path.insert(0, os.path.join(here, "..", "..", "python"))
import glamio  # noqa: E402
from rsbin import read_rsbin  # noqa: E402

run, step, base, box = sys.argv[1], int(sys.argv[2]), sys.argv[3], float(sys.argv[4])
cat, meta = read_rsbin(base)
a, a_v = float(meta["a"]), float(meta["a_v"])
bdm = glamio.read_catalogue(os.path.join(run, "CATALOGS", f"CatshortV.{step:04d}.0000"))
sel = (bdm["distinct"] == 0) & (bdm["Mtot"] > 1e13)
pb = np.mod(np.c_[bdm["x"], bdm["y"], bdm["z"]][sel], box)
vb = np.c_[bdm["vx"], bdm["vy"], bdm["vz"]][sel]
rv = bdm["Rvir"][sel] * 1e-3
host = cat["pid"] == -1
tree = cKDTree(np.mod(np.c_[cat["x"], cat["y"], cat["z"]][host], box), boxsize=box)
d, j = tree.query(pb)
ok = d < 0.3 * rv
vi = np.c_[cat["vx"], cat["vy"], cat["vz"]][host][j[ok]]
slope = np.sum(vi * vb[ok]) / np.sum(vb[ok] ** 2)
resid = np.linalg.norm(vi - slope * vb[ok], axis=1)
err = np.median(resid) / np.sqrt(np.sum(np.linalg.norm(vb[ok], axis=1) ** 2) / ok.sum()) / np.sqrt(ok.sum())
print(f"step {step}: a = {a:.6f}, a_v = {a_v:.6f}; expected v_inline / v_BDM = a / a_v = {a / a_v:.5f}")
print(f"1. BDM Mtot > 1e13: {sel.sum()}, matched {ok.sum()}: slope {slope:.5f} +- {err:.5f} "
      f"({(slope - a / a_v) / err:+.1f} sigma); residual median {np.median(resid):.0f} km/s")
h = host & (cat["num_p"] >= 100) & (cat["mgrav"] > 0) & (cat["vrms"] > 0)
p = np.polyfit(np.log10(cat["mgrav"][h]), np.log10(cat["vrms"][h]), 1)
print(f"2. log vrms = {p[0]:.4f} log Mvir + {p[1]:.3f} over {h.sum()} hosts (virial scaling ~ 1/3)")
for lo, hi in [(1e12, 3e12), (3e12, 1e13), (1e13, 1e14), (1e14, 1e16)]:
    b = h & (cat["mgrav"] >= lo) & (cat["mgrav"] < hi)
    if b.sum():
        print(f"3. T/|U| median of hosts [{lo:.0e},{hi:.0e}): {np.median(cat['kin_to_pot'][b]):.4f} (n {b.sum()})")
