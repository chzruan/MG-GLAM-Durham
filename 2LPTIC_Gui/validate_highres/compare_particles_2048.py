#!/usr/bin/env python3
"""Per-particle spot check: our 2048^3 FML IC vs Gui's HEFT 2048^3 IC.

Gui's files have broken IDs (all 1), so match particles by position with a
KD-tree (positions should agree to ~float32 precision). His velocities are
a factor 100*box/a = 5.12e6 too small; rescale before comparing.
Uses our snap file .0 and HEFT file .0 (x in [0,4) Mpc/h slab).
"""
import numpy as np, struct
from scipy.spatial import cKDTree

BOX = 1024.0
VEL_FIX = 100.0 * BOX / 0.02   # 5.12e6: correct Gadget vel_norm / Gui's vel_norm

def read_gadget1(fn):
    raw = open(fn, 'rb').read()
    n = struct.unpack_from('6i', raw, 4)[1]
    o = 4 + 256 + 4
    bs = struct.unpack_from('i', raw, o)[0]; assert bs == 12*n
    pos = np.frombuffer(raw, np.float32, 3*n, o+4).reshape(-1, 3)
    o += 8 + bs
    bs = struct.unpack_from('i', raw, o)[0]; assert bs == 12*n
    vel = np.frombuffer(raw, np.float32, 3*n, o+4).reshape(-1, 3)
    return pos, vel

pos_o, vel_o = read_gadget1('snap/IC_ours_Np2048_L1024.0')
pos_r, vel_r = read_gadget1(
    '/cosma8/data/dp203/bl267/Projects/Ongoing/HEFT/ICs/IC_highres/IC_Np1d_2048_L_1024_2LPT.0')
print(f"ours file0: {len(pos_o):,}  HEFT file0: {len(pos_r):,} particles")
print(f"x range ours [{pos_o[:,0].min():.2f},{pos_o[:,0].max():.2f}] "
      f"HEFT [{pos_r[:,0].min():.2f},{pos_r[:,0].max():.2f}]")

# interior of the HEFT slab, away from Eulerian domain edges
so = (pos_o[:,0] > 0.5) & (pos_o[:,0] < 3.5)
sr = (pos_r[:,0] > 0.5) & (pos_r[:,0] < 3.5)
po, vo = pos_o[so], vel_o[so]
pr, vr = pos_r[sr], vel_r[sr] * VEL_FIX
print(f"interior subsets: ours {so.sum():,}, HEFT {sr.sum():,}")

tree = cKDTree(pr)
dist, idx = tree.query(po, k=1, distance_upper_bound=0.01)
m = np.isfinite(dist)
print(f"matched within 0.01 Mpc/h: {m.sum():,} of {len(po):,} ({m.mean()*100:.3f}%)")

dp = po[m] - pr[idx[m]]
dv = vo[m] - vr[idx[m]]
u_rms = np.sqrt((vr[idx[m]]**2).mean())
print(f"POSITION diff: rms = {np.sqrt((dp**2).mean()):.3e} Mpc/h, max = {np.abs(dp).max():.3e}")
print(f"VELOCITY (HEFT rescaled x{VEL_FIX:.3g}): u_rms HEFT = {u_rms:.2f}, ours = {np.sqrt((vo[m]**2).mean()):.2f} km/s")
print(f"VELOCITY diff: rms = {np.sqrt((dv**2).mean()):.3e} km/s ({np.sqrt((dv**2).mean())/u_rms:.2e} of u_rms)")
