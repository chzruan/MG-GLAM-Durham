#!/usr/bin/env python3
"""Per-particle comparison of our FML IC vs the 2LPTic reference (1024^3, L=1024).

Reads our snap file .0 (Eulerian slab of rank 0) and the reference ics.0
(Lagrangian slab of 2LPTic task 0), determines the reference ID ordering
empirically, matches particles on lattice site, and compares positions and
velocities particle by particle.
"""
import numpy as np, struct, sys

NP = 1024          # particles per dimension
BOX = 1024.0       # Mpc/h
CELL = BOX / NP

def read_gadget1(fn, id_dtype):
    raw = open(fn, 'rb').read()
    n = struct.unpack_from('6i', raw, 4)[1]
    o = 4 + 256 + 4
    bs = struct.unpack_from('i', raw, o)[0]; assert bs == 12*n, (bs, n)
    pos = np.frombuffer(raw, np.float32, 3*n, o+4).reshape(-1, 3)
    o += 8 + bs
    bs = struct.unpack_from('i', raw, o)[0]; assert bs == 12*n
    vel = np.frombuffer(raw, np.float32, 3*n, o+4).reshape(-1, 3)
    o += 8 + bs
    bs = struct.unpack_from('i', raw, o)[0]
    assert bs == n*np.dtype(id_dtype).itemsize, f"id block {bs} vs n={n}"
    ids = np.frombuffer(raw, id_dtype, n, o+4)
    return pos, vel, ids

ours_fn = sys.argv[1] if len(sys.argv) > 1 else 'snap/IC_ours_Np1024_L1024.0'
ref_fn = '/cosma7/data/dp004/bl267/Runs/DEGRACE/ICs/IC_data/L1024/Node_002/ics.0'

pos_o, vel_o, ids_o = read_gadget1(ours_fn, np.int64)
pos_r, vel_r, ids_r = read_gadget1(ref_fn, np.uint32)
print(f"ours file: {len(ids_o):,} particles; ref file: {len(ids_r):,} particles")
print(f"u_rms(1d) ours(file0)={np.sqrt((vel_o**2).mean()):.2f}  ref(file0)={np.sqrt((vel_r**2).mean()):.2f} km/s")

# --- lattice offsets (fractional position of the unperturbed grid) ---
for name, p in (("ours", pos_o), ("ref", pos_r)):
    fr = (p / CELL) % 1.0
    fr = np.where(fr > 0.5, fr - 1.0, fr)
    print(f"{name}: lattice frac offset median = {np.median(fr, axis=0)}, rms = {fr.std():.3f} cells")

# --- determine reference ID ordering from small-displacement particles ---
lat_r = np.round(pos_r / CELL).astype(np.int64) % NP
disp_r = pos_r - lat_r * CELL
disp_r -= BOX * np.round(disp_r / BOX)
small = (np.abs(disp_r) < 0.35 * CELL).all(axis=1)
ix, iy, iz = lat_r[small].T
idr = ids_r[small].astype(np.int64)
conventions = {
    'x-fastest id=1+ix+N(iy+N iz)': 1 + ix + NP*(iy + NP*iz),
    'z-fastest id=1+iz+N(iy+N ix)': 1 + iz + NP*(iy + NP*ix),
    'x-fastest id=0+...':           ix + NP*(iy + NP*iz),
    'z-fastest id=0+...':           iz + NP*(iy + NP*ix),
}
match_frac = {name: (pred == idr).mean() for name, pred in conventions.items()}
best = max(match_frac, key=match_frac.get)
print(f"ref ID convention: {best} (match {match_frac[best]:.4f} on {small.sum():,} small-disp particles)")
assert match_frac[best] > 0.999, match_frac

# --- decode lattice site from IDs on both sides ---
def decode(ids, conv):
    a = ids.astype(np.int64) - (1 if '1+' in conv else 0)
    if conv.startswith('x-fastest'):
        return a % NP, (a // NP) % NP, a // (NP*NP)
    return a // (NP*NP), (a // NP) % NP, a % NP

ox, oy, oz = decode(ids_o, 'x-fastest id=1+ix+N(iy+N iz)')
rx, ry, rz = decode(ids_r, best)

# ours file .0 = Eulerian x in [0,32); use ref-side Lagrangian ix in [2,29]
# so every matched particle is guaranteed present in both single files
key_o = ox + NP*(oy.astype(np.int64) + NP*oz.astype(np.int64))
key_r = rx + NP*(ry.astype(np.int64) + NP*rz.astype(np.int64))
sel_r = (rx >= 2) & (rx <= 29)
kr = key_r[sel_r]
order_o = np.argsort(key_o)
posn = np.searchsorted(key_o, kr, sorter=order_o)
ok = posn < len(key_o)
posn[~ok] = 0
matched = ok & (key_o[order_o[posn]] == kr)
io = order_o[posn[matched]]
ir = np.flatnonzero(sel_r)[matched]
print(f"matched particles: {len(io):,} of {sel_r.sum():,} candidate ref sites")

dp = pos_o[io] - pos_r[ir]
dp -= BOX * np.round(dp / BOX)
dv = vel_o[io] - vel_r[ir]
disp_ref = pos_r[ir] - np.stack([rx[ir], ry[ir], rz[ir]], axis=1) * CELL
disp_ref -= BOX * np.round(disp_ref / BOX)
rms_disp = np.sqrt((disp_ref**2).mean())
print(f"\nreference displacement rms(1d) = {rms_disp:.4f} Mpc/h")
print(f"POSITION diff: rms = {np.sqrt((dp**2).mean()):.3e} Mpc/h "
      f"({np.sqrt((dp**2).mean())/CELL:.2e} cells, {np.sqrt((dp**2).mean())/rms_disp:.2e} of disp rms); "
      f"max|dp| = {np.abs(dp).max():.3e} Mpc/h")
u_rms = np.sqrt((vel_r[ir]**2).mean())
print(f"VELOCITY diff: rms = {np.sqrt((dv**2).mean()):.3e} km/s "
      f"({np.sqrt((dv**2).mean())/u_rms:.2e} of ref u_rms={u_rms:.1f}); "
      f"max|dv| = {np.abs(dv).max():.3e} km/s")
