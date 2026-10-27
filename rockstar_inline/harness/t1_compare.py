"""Test T1 of the in-line Rockstar: the in-line catalogue against standalone MPI-Rockstar on the same particles.
The in-line catalogue holds distinct haloes only (no parent in find_parents with the kernel radius), so the
comparison is of distinct haloes ("hosts"): the reference's rows with PID = -1.
usage: python3 t1_compare.py INLINE_BASE REF_DIR BOX --baseline DIR [--strict --hosts-from DIR
                            --baseline-hosts-from DIR] [--shift F] [--npmin N] [--max-gap G]
  INLINE_BASE  .../Rockstar.<step>.<box> (without .bin/.meta), or a second reference directory (its
               out_0.parents.list): two standalone runs compared with each other, halo by halo not judged
  REF_DIR      directory with the reference out_0.parents.list (find_parents output)
  --baseline   a second standalone run of the same output (another rank count, or input order, or the same
               shift as the in-line run); its halo-by-halo agreement with REF_DIR is the baseline
  --strict     the reference ran with STRICT_SO_MASSES = 1: its Mvir_all is the strict mass, compared with the
               in-line m_strict. Its PIDs use the strict radius, so its hosts are taken from --hosts-from, the
               same run with STRICT_SO_MASSES = 0 (the same haloes, matched by position and Np), whose PIDs use
               the kernel radius as the in-line selection does; --baseline-hosts-from (and --inline-hosts-from
               for a strict INLINE_BASE directory) do the same for the other strict runs
  --shift F    the in-line positions were shifted by F * BOX (test T3; rs_harness inline --shift F)
Criteria:
  mass functions: the two catalogues hold the same haloes but for those the finders treat differently, so the
    ratio's error comes from a paired bootstrap (both catalogues resampled together over 4^3 subvolumes, 300
    draws); |ratio - 1| <= max(3 sigma, 0.005) in every bin;
  halo by halo: the reference's hosts with Np >= npmin, matched to the nearest in-line host within 0.1 R200m;
    the matched fraction and, of the matched, the fraction with |dln M| < 0.02 must each be at most --max-gap
    (0.03) below the same figures for the baseline run against the reference. Standalone runs disagree with
    each other because Rockstar's own FoF merges groups in dense regions in an order that changes between
    runs; the adapter's exact FoF adds a gap of about 2-2.5 points to that spread at every redshift.
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
ap.add_argument("--hosts-from")
ap.add_argument("--inline-hosts-from")
ap.add_argument("--baseline")
ap.add_argument("--baseline-hosts-from")
ap.add_argument("--max-gap", type=float, default=0.03)
ap.add_argument("--shift", type=float, default=0.0)
ap.add_argument("--npmin", type=int, default=100)
a = ap.parse_args()
box = a.box
RHO_M_PER_OM = 2.77519737e11

def wrap(pos):
    """into [0, box): np.mod of a float32 can round up to box, which cKDTree(boxsize=box) refuses"""
    pos = np.mod(pos, box)
    pos[pos >= box] -= box
    return pos


def parents_list(d):
    return np.loadtxt(os.path.join(d, "out_0.parents.list"), comments="#")


def hosts_from(r, d):
    """PID = -1 of the same haloes in directory d (rows matched by x, y, z and Np; the order may differ)"""
    q = parents_list(d)
    ka, kb = np.lexsort((r[:, 7], r[:, 10], r[:, 9], r[:, 8])), np.lexsort((q[:, 7], q[:, 10], q[:, 9], q[:, 8]))
    if len(q) != len(r) or not np.array_equal(r[ka][:, 7:11], q[kb][:, 7:11]):
        sys.exit(f"{d} does not hold the same haloes")
    host = np.empty(len(r), bool)
    host[ka] = q[kb][:, -1] == -1
    return host


def load_ref(d, shift=0.0, hosts_dir=None):
    # out_*.list + find_parents' PID: 2 Mvir (= mgrav, BOUND_PROPS = 1), 7 Np, 8-10 x y z,
    # 19 Mvir_all (= m, the strict m with STRICT_SO_MASSES), 20 M200b, last PID
    r = parents_list(d)
    host = hosts_from(r, hosts_dir) if hosts_dir else r[:, -1] == -1
    return {"M": r[:, 2], "Mall": r[:, 19], "Np": r[:, 7], "pos": wrap(r[:, 8:11] - shift), "M200b": r[:, 20],
            "host": host, "n": len(r)}


if a.strict and not a.hosts_from:
    sys.exit("--strict needs --hosts-from (the reference's STRICT_SO_MASSES = 0 twin)")
inline_is_ref = os.path.isdir(a.inline_base)
if not inline_is_ref and not a.baseline:
    sys.exit("an in-line catalogue needs --baseline (a second standalone run of the same output)")
if a.strict and a.baseline and not a.baseline_hosts_from:
    sys.exit("--strict needs --baseline-hosts-from (the baseline's STRICT_SO_MASSES = 0 twin)")
R = load_ref(a.ref_dir, hosts_dir=a.hosts_from)
shift = a.shift * box
if inline_is_ref:
    if a.strict and not a.inline_hosts_from:
        sys.exit("a strict second reference needs --inline-hosts-from")
    I = load_ref(a.inline_base, shift, a.inline_hosts_from if a.strict else None)
    om = 0.3089
    print(f"second reference {a.inline_base}: {I['n']} haloes, {I['host'].sum()} hosts")
else:
    cat, meta = read_rsbin(a.inline_base)
    om = float(meta["Omega_m"])
    I = {"M": cat["mgrav"].astype(float), "Mall": (cat["m_strict"] if a.strict else cat["m"]).astype(float),
         "Np": cat["num_p"].astype(float), "pos": wrap(np.c_[cat["x"], cat["y"], cat["z"]].astype(float) - shift),
         "M200b": cat["alt_m_1"].astype(float), "host": np.ones(len(cat), bool), "n": len(cat)}
    print(f"in-line {a.inline_base}: {I['n']} distinct haloes")
print(f"reference {a.ref_dir}: {R['n']} haloes, {R['host'].sum()} hosts")
ok_all = True


NSUB, NBOOT = 4, 300
rng = np.random.default_rng(12345)
draws = rng.integers(0, NSUB ** 3, size=(NBOOT, NSUB ** 3))
weights = np.stack([np.bincount(d, minlength=NSUB ** 3) for d in draws])          # NBOOT x NSUB^3


def sub_counts(m, pos, edges):
    """haloes per (subvolume, mass bin)"""
    c = np.minimum((pos / box * NSUB).astype(int), NSUB - 1)
    sv = (c[:, 0] * NSUB + c[:, 1]) * NSUB + c[:, 2]
    h, _, _ = np.histogram2d(sv, m, bins=[np.arange(NSUB ** 3 + 1), edges])
    return h


def mf(label, sel_i, sel_r, mi, mr, edges):
    global ok_all
    ci = sub_counts(mi[sel_i], I["pos"][sel_i], edges)
    cr = sub_counts(mr[sel_r], R["pos"][sel_r], edges)
    ni, nr = ci.sum(0), cr.sum(0)
    bi, br = weights @ ci, weights @ cr                                              # NBOOT x bins
    print(f"  {label}: bin, reference, in-line, ratio, allowed (max(3 sigma paired, 0.005))")
    for k in range(len(edges) - 1):
        if nr[k] == 0:
            continue
        r = ni[k] / nr[k]
        with np.errstate(divide="ignore", invalid="ignore"):
            sig = np.nanstd(bi[:, k] / br[:, k])
        allow = max(3 * sig, 0.005)
        flag = "" if abs(r - 1) <= allow else "  <-- outside"
        ok_all &= abs(r - 1) <= allow
        print(f"    [{edges[k]:.0e},{edges[k+1]:.0e}) {int(nr[k]):8d} {int(ni[k]):8d}  {r:.4f}  +-{allow:.4f}{flag}")


edges = np.array([2e11, 5e11, 1e12, 3e12, 1e13, 3e13, 1e14, 1e16])
mf("hosts, bound Mvir", I["host"], R["host"], I["M"], R["M"], edges)
mf("hosts, strict Mvir_all" if a.strict else "hosts, Mvir_all", I["host"], R["host"], I["Mall"], R["Mall"], edges)

# halo by halo: the reference's hosts with Np >= npmin matched to the nearest host of the other run within 0.1 R200m
def match(other):
    sel = R["host"] & (R["Np"] >= a.npmin)
    r200 = (3 * R["M200b"][sel] / (4 * np.pi * 200 * om * RHO_M_PER_OM)) ** (1 / 3)
    d, j = cKDTree(other["pos"][other["host"]], boxsize=box).query(R["pos"][sel], k=1)
    okm = d < 0.1 * r200
    dln = np.abs(np.log(other["M"][other["host"]][j[okm]] / R["M"][sel][okm]))
    print(f"    hosts with Np >= {a.npmin}: {sel.sum()}; matched within 0.1 R200m {okm.mean():.4f}; of the matched, "
          f"|dln Mvir| < 0.02 for {np.mean(dln < 0.02):.4f} [{okm.mean() * np.mean(dln < 0.02):.4f} of all], "
          f"< 0.05 for {np.mean(dln < 0.05):.4f}; median {np.median(dln):.4f}")
    return okm.mean(), np.mean(dln < 0.02)


print("  halo by halo, " + ("second reference" if inline_is_ref else "in-line") + " vs reference:")
frac, f02 = match(I)
if a.baseline:
    B = load_ref(a.baseline, shift, a.baseline_hosts_from if a.strict else None)
    print(f"  halo by halo, baseline {a.baseline} vs reference:")
    bfrac, bf02 = match(B)
    ok_hh = frac >= bfrac - a.max_gap and f02 >= bf02 - a.max_gap
    print(f"  gaps to the baseline: matched {bfrac - frac:+.4f}, within 2 % {bf02 - f02:+.4f} "
          f"(each <= {a.max_gap}): {'ok' if ok_hh else 'FAIL'}")
    ok_all &= ok_hh
print("T1", "PASS" if ok_all else "FAIL")
sys.exit(0 if ok_all else 1)
