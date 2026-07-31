#!/usr/bin/env python3
"""Validate our FML-generated IC against a reference IC set.
Computes at a common Nmesh grid: auto P_ours, auto P_ref, cross ->
r(k)=P_x/sqrt(P_ours*P_ref) and ratio P_ours/P_ref.
Run under MPI:  mpirun -n N python measure_validate.py <ours_glob> <ref_glob> <Nmesh> <out.npz>
"""
import sys, os, numpy as np
from nbodykit.source.catalog import Gadget1Catalog
from nbodykit.lab import FFTPower

OURS, REF = sys.argv[1], sys.argv[2]
NM = int(sys.argv[3])
OUT = sys.argv[4]

def make_mesh(path, Nm):
    cat = Gadget1Catalog(path, ptype=1)
    box = float(cat.attrs['BoxSize'])
    m = cat.to_mesh(Nmesh=Nm, BoxSize=box, resampler='cic',
                    compensated=True, interlaced=True, position='Position')
    return m, cat, box

def unpack(r):
    p = r.power
    return (np.asarray(p['k']), np.asarray(p['power'].real),
            np.asarray(p['modes']), float(p.attrs['shotnoise']))

mo, cato, box = make_mesh(OURS, NM)
mr, catr, _   = make_mesh(REF, NM)
comm = cato.comm
if comm.rank == 0:
    print(f"[val] MPI size={comm.size} ours csize={cato.csize} ref csize={catr.csize} box={box}", flush=True)

r_o = FFTPower(mo, mode='1d')
r_r = FFTPower(mr, mode='1d')
r_x = FFTPower(mo, mode='1d', second=mr)

if comm.rank == 0:
    k, Po, nmo, sno = unpack(r_o)
    _, Pr, _, snr   = unpack(r_r)
    _, Px, _, _     = unpack(r_x)
    np.savez(OUT, k=k, Pours=Po, Pref=Pr, Pcross=Px, modes=nmo,
             sn_ours=sno, sn_ref=snr, box=box, Nmesh=NM,
             csize_ours=int(cato.csize), csize_ref=int(catr.csize))
    r = Px/np.sqrt(Po*Pr)
    ratio = Po/Pr
    sel = k < 1.0
    print(f"[val] wrote {OUT}", flush=True)
    print(f"[val] r(k<1): min={r[sel].min():.6f}  ratio P_ours/P_ref (k<1): "
          f"median={np.median(ratio[sel]):.4f} min={ratio[sel].min():.4f} max={ratio[sel].max():.4f}", flush=True)
