"""Halo properties of the in-line Rockstar against standalone MPI-Rockstar (part of test T1/T4).
usage: python3 t1_props.py REF_DIR OTHER [OTHER2 ...]   (OTHER: an in-line .bin base or another reference directory)
For distinct haloes (Np >= 100; the in-line catalogue holds only those, a reference its rows with PID = -1) in
bins of bound Mvir: the medians of Vmax, Rs, spin, T/|U| and b/a, with the error of the reference median
(1.25 sigma / sqrt n). Use a reference run with STRICT_SO_MASSES = 0: its PIDs use the kernel radius R, as the
in-line selection does. Descriptive (no pass/fail): read the ratios against the error column and against the
spread between two references.
Part of the MG-GLAM in-line Rockstar adapter; GNU GPL v3 or later."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rsbin import read_rsbin  # noqa: E402

# out_*.list columns: 2 Mvir, 3 Vmax, 6 Rs, 7 Np, 17 Spin, 27 b_to_a, 37 T/|U|, last PID
COLS = {"Vmax": 3, "Rs": 6, "spin": 17, "b_to_a": 27, "T/|U|": 37}
FIELDS = {"Vmax": "vmax", "Rs": "rs", "spin": "spin", "b_to_a": "b_to_a", "T/|U|": "kin_to_pot"}


def load(path):
    if os.path.isdir(path):
        r = np.loadtxt(os.path.join(path, "out_0.parents.list"), comments="#")
        d = {k: r[:, c] for k, c in COLS.items()}
        d.update(M=r[:, 2], Np=r[:, 7], host=r[:, -1] == -1)
        return d
    c, _ = read_rsbin(path)
    d = {k: c[f].astype(float) for k, f in FIELDS.items()}
    d.update(M=c["mgrav"].astype(float), Np=c["num_p"].astype(float), host=np.ones(len(c), bool))
    return d


paths = sys.argv[1:]
cats = [load(p) for p in paths]
bins = [(1e12, 3e12), (3e12, 1e13), (1e13, 1e14)]
print("sources: " + "; ".join(f"[{i}] {p}" for i, p in enumerate(paths)))
for lo, hi in bins:
    print(f"hosts with Np >= 100, Mvir in [{lo:.0e}, {hi:.0e}):")
    sel = [(c["host"] & (c["Np"] >= 100) & (c["M"] >= lo) & (c["M"] < hi)) for c in cats]
    print("   n: " + "  ".join(f"[{i}] {s.sum()}" for i, s in enumerate(sel)))
    for k in COLS:
        ref = cats[0][k][sel[0]]
        m0 = np.median(ref)
        err = 1.25 * np.std(ref) / np.sqrt(max(len(ref), 1))
        others = "  ".join(f"[{i}] {np.median(c[k][s]) / m0:.4f}" for i, (c, s) in enumerate(zip(cats, sel)) if i)
        print(f"   {k:7s} median [0] {m0:10.4g} (+-{err / abs(m0):.4f} rel);  ratios {others}")
