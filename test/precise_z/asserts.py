#!/usr/bin/env python3
"""Assertions for the exact-redshift snapshot tests (see run_tests.sh).

Config convention under test: '#outputs = N' is the legacy nearest-step
schedule; '#outputs = -N' requests the same N redshifts EXACTLY.

Independently replicates, in Real*4 (float32) arithmetic, the timestep-ladder
construction of PMP2main.f90:Initialize, the legacy zout nearest-step marking,
the SetExactSteps insert-or-merge processing, and the run stop rule (legacy
half-step exit / exact-mode NlastX).  Every PMcrd.NNNN.DAT header produced by
the runs is then checked against this oracle:

  - exact-target steps must equal float32(1/(1+z)) BITWISE,
  - all other output epochs must match the ladder to <= 2 ulp,
  - the set of produced snapshot numbers must match exactly (no missing
    outputs, no duplicates, no spurious extras).

Pure stdlib; no numpy required.
"""
import glob
import os
import re
import struct
import sys

# ---- configuration mirrored from run_tests.sh / Init.base.dat -------------
ZINIT = 49.0
DA0 = 8.0e-4
NSTEPS = 100        # steps piped into PMP2main (default)
ZTOL = 0.005        # merge tolerance in SetExactSteps (unpinned steps only)
DUPTOL = 1.0e-6     # targets closer than this (relative in a) are duplicates
LOG_ATOL = 2.0e-5   # precision of the A= column in main.log

SIM_CASES = {       # case dir -> configuration (zout list + sign of #outputs)
    "case5a_legacy": dict(zout=[2.00, 0.50], exact=False),
    "case5b_trailing": dict(zout=[2.00, 0.50], exact=False),
    "case1_merge": dict(zout=[2.00, 0.49337], exact=True),
    "case2_insert": dict(zout=[2.00, 0.512], exact=True),
    "case3_multi": dict(zout=[2.5, 0.512, 1.1], exact=True),
    # z=0 must be INSERTED after the step where the legacy half-step exit
    # rule already fires (regression test for the NlastX stop-index fix)
    "case6_tailinsert": dict(zout=[0.0], exact=True, da0=6.9e-4, nsteps=300),
    # the four DEGRACE pairs closer than ZTOL (0.49%, 0.19%, 0.10%, 0.045% in
    # a): every member must be scheduled, the later one as a short step
    "case7_closepairs": dict(zout=[2.0, 1.98541, 1.75, 1.74475, 1.00202, 1.0, 0.37562, 0.375], exact=True, nsteps=140),
    # first member merges onto a ladder step (pinned); the second, 0.09%
    # later, must be INSERTED and must not re-move the pinned step
    "case8_pinned": dict(zout=[0.49337, 0.4920], exact=True),
    # per-output classes (I5) must not change the schedule; they must reach
    # OutputSchedule.dat (I4) attached to the right redshift after sorting
    "case9_classes": dict(zout=[2.5, 0.512, 1.1], exact=True, classes=[2, 1, 2]),
    # lens planes (I2): 2 slabs per axis, plane grid 72 (= Ngrid) and 36 (2x2 block mean)
    "case10_lens72": dict(zout=[2.5, 0.512, 1.1], exact=True, lens=(2, 72)),
    "case10_lens36": dict(zout=[2.5, 0.512, 1.1], exact=True, lens=(2, 36)),
    # low-k mode dumps (I8): density from POWERfft5's transform (iPower=1), momentum on N_dep=36
    "case11_modes": dict(zout=[2.5, 0.512, 1.1], exact=True, modes=dict(kc=0.5, ndep=36), subm=1),
    # f(R) smoke test: BDM + f_R at halo centres (I10); the box may host few or no haloes
    "case12_frhalo": dict(zout=[1.1, 0.512], exact=True, frhalo=True),
}
BOX = 150.0
NGRID = 72
C_H0 = 2997.92458   # Mpc/h, comoving-distance convention of OutputSchedule.dat
ERR_CASES = {       # case dir -> (log file, required message, forbidden product)
    "case4a_zbig": ("init.log", "outside the simulation range", "Setup.dat"),
    "case4b_zneg": ("init.log", "outside the simulation range", "Setup.dat"),
    "case4c_dup": ("Run1/main.log", "are duplicates", "Run1/PMcrd.0*.DAT"),
    "case4d_range": ("Run1/main.log", "at or before the initial redshift", "Run1/PMcrd.0*.DAT"),
    "case4e_classes": ("init.log", "one output class per listed redshift", "Setup.dat"),
}


def f32(x):
    """Round a Python float to Real*4."""
    return struct.unpack("!f", struct.pack("!f", x))[0]


# ---- oracle ---------------------------------------------------------------
def build_ladder(da0):
    """Replicate the Alist/dAlist builder in Initialize (all ops in Real*4)."""
    aexpn0 = f32(1.0 / f32(1.0 + f32(ZINIT)))
    astep0 = f32(da0)
    stepfactor = f32(astep0 / aexpn0)
    a, da = aexpn0, astep0
    A, dA = [None], [None]  # 1-based
    while True:
        if da < f32(f32(stepfactor / f32(1.25)) * a) and a < f32(0.300):
            da = f32(1.5 * da)
        a = f32(a + da)
        A.append(a)
        dA.append(da)
        if a >= 1.0:
            break
    return aexpn0, A, dA


def legacy_marks(A, zouts):
    """Replicate the zout nearest-step marking loop."""
    n = len(A) - 1
    marked = set()
    for z in zouts:
        at = f32(1.0 / f32(1.0 + f32(z)))
        j = 2
        while j <= n:
            if at < A[j]:
                break
            j += 1
        aj = A[j] if j <= n else 0.0        # Alist is zero-filled past Ntotal
        da1 = f32(aj - at)
        da0_ = f32(at - A[j - 1])
        marked.add(j if da1 < da0_ else j - 1)
    return marked


def apply_exact(aexpn0, A, dA, zex):
    """Replicate SetExactSteps: insert-or-merge each target, ascending in a.

    Returns (A, dA, marks, xmarks) where xmarks maps step -> f32 target epoch.
    """
    A, dA = list(A), list(dA)
    marks = set()
    xmarks = {}
    for at in sorted(f32(1.0 / f32(1.0 + f32(z))) for z in zex):
        n = len(A) - 1
        atol = f32(ZTOL * at)
        j = 1
        while A[j] < at and j < n:
            j += 1
        jn = j
        if j > 1 and f32(at - A[j - 1]) < f32(A[j] - at):
            jn = j - 1
        if abs(f32(A[jn] - at)) < atol and jn not in xmarks:   # merge: move an unpinned step
            A[jn] = at
            dA[jn] = f32(at - (A[jn - 1] if jn > 1 else aexpn0))
            if jn < n:
                dA[jn + 1] = f32(A[jn + 1] - at)
            marks.add(jn)
            xmarks[jn] = at
        else:                               # insert a new step before j
            aprev = A[j - 1] if j > 1 else aexpn0
            A.insert(j, at)
            dA.insert(j, f32(at - aprev))
            dA[j + 1] = f32(A[j + 1] - at)
            marks = {m + 1 if m >= j else m for m in marks}
            xmarks = {(s + 1 if s >= j else s): v for s, v in xmarks.items()}
            marks.add(j)
            xmarks[j] = at
    return A, dA, marks, xmarks


def stop_index(A, dA, marks, exact):
    """First step where the run exits: legacy half-step rule, or NlastX."""
    n = len(A) - 1
    s = n
    for i in range(1, n + 1):
        if A[i] >= f32(1.0 - f32(dA[i] / 2.0)):
            s = i
            break
    if exact and marks:
        s = max(s, max(marks))              # NlastX: never skip a marked step
    return s


def expected_outputs(A, marks, nsteps, stop):
    """Numbered snapshots: marked steps up to the last executed one, plus the
    end-of-run iSave dump if that step is unmarked."""
    last = min(nsteps, stop)
    out = {i: A[i] for i in sorted(marks) if i <= last}
    out.setdefault(last, A[last])
    return out


# ---- actual outputs -------------------------------------------------------
def read_aexpn(path):
    """AEXPN from a PMcrd header (big-endian: marker | 45-char HEADER | reals)."""
    with open(path, "rb") as fh:
        d = fh.read(4 + 45 + 4)
    return struct.unpack("!f", d[4 + 45:4 + 45 + 4])[0]


def actual_outputs(rundir):
    out = {}
    for p in glob.glob(os.path.join(rundir, "PMcrd.[0-9]*.DAT")):
        step = int(re.search(r"PMcrd\.(\d+)\.DAT$", p).group(1))
        out[step] = read_aexpn(p)
    return out


def a_sequence(mainlog):
    """The deterministic A= column of the per-step log lines."""
    seq = []
    with open(mainlog) as fh:
        for line in fh:
            m = re.match(r"(?:\*{4} STEP=|Step =)\s*\d+\s+A=\s*([0-9.]+)", line)
            if m:
                seq.append(m.group(1))
    return seq


def read_schedule(path):
    """OutputSchedule.dat -> {step: (a, z, chi, exact, class)}."""
    out = {}
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        s, a, z, chi, ex, cl = line.split()
        out[int(s)] = (float(a), float(z), float(chi), int(ex), int(cl))
    return out


def setup_cosmology(path):
    om = oml = None
    for line in open(path):
        if "Omega Matter" in line:
            om = float(line.split()[0])
        elif "Omega Lambda" in line:
            oml = float(line.split()[0])
    return om, oml


def chi_of(a, om, oml):
    """Comoving distance, flat matter+Lambda, by Simpson in ln a (no scipy)."""
    import math
    n = 20000
    x0 = math.log(a)
    h = -x0 / n
    s = 0.0
    for i in range(n + 1):
        x = x0 + h * i
        f = 1.0 / (math.exp(x) * math.sqrt(om * math.exp(-3.0 * x) + oml))
        s += f if i in (0, n) else (4.0 * f if i % 2 else 2.0 * f)
    return C_H0 * s * h / 3.0


def read_planes(path, nslab, nlens):
    # raw little-endian float32; Fortran column-major (i,j) per plane -> numpy [plane, j, i]
    import numpy as np
    arr = np.fromfile(path, dtype="<f4")
    assert arr.size == 3 * nslab * nlens * nlens, f"{path}: {arr.size} values"
    return arr.reshape(3 * nslab, nlens, nlens)


def read_projection(path, ngrid):
    # ProjectDENSIT output: big-endian Fortran unformatted; record 1 = (AEXPN, ISTEP, NGRID),
    # record 2 = Den(NGRID, NGRID) real*4 column-major -> numpy [j, i]
    import numpy as np
    with open(path, "rb") as fh:
        raw = fh.read()
    off = 4 + 12 + 4 + 4
    return np.frombuffer(raw, dtype=">f4", count=ngrid * ngrid, offset=off).reshape(ngrid, ngrid)


def read_modes(path, n):
    # complex64 little-endian, Fortran C(0:n,-n:n,-n:n) -> numpy [ix+n, iy+n, iz]
    import numpy as np
    arr = np.fromfile(path, dtype="<c8")
    assert arr.size == (2 * n + 1) ** 2 * (n + 1), f"{path}: {arr.size} values"
    return arr.reshape(2 * n + 1, 2 * n + 1, n + 1)


def read_powerdm(path):
    # PowerDM.R: 'bin k N Preal' rows after the header lines
    import numpy as np
    rows = []
    for line in open(path):
        p = line.split()
        if len(p) == 4 and p[0].isdigit():
            rows.append([float(v) for v in p])
    return np.array(rows)


def hermitian_ok(C, n, tol=1e-5):
    # at iz = 0: C(-ix,-iy,0) == conj(C(ix,iy,0))
    import numpy as np
    plane = C[:, :, 0]
    scale = max(np.max(np.abs(plane)), 1e-30)
    return np.max(np.abs(plane[::-1, ::-1] - np.conj(plane))) < tol * scale


def read_subsample(fb, fm):
    # int16 arrays qx,qy,qz,vx,vy,vz -> positions [Mpc/h] in [0,L), velocities [km/s]
    import numpy as np
    meta = read_meta(fm)
    num = lambda k: float(meta[k].split()[0])          # values may carry a unit suffix
    nsel, L, dv = int(num("n_sel")), num("Box"), num("dv_quant")
    arr = np.fromfile(fb, dtype="<i2").reshape(6, nsel)
    x = (arr[:3].astype(np.int64) + 32768 + 0.5) * L / 65536.0
    v = arr[3:].astype(np.float64) * dv
    return meta, x, v


def cube_from_particles(x, w, nd, L):
    # CIC deposit of weights w at positions x (3, N) onto nd^3, then delta_k = N^-3 sum e^{-ikx}
    import numpy as np
    pos = x * nd / L
    i0 = np.floor(pos).astype(np.int64)
    d = pos - i0
    g = np.zeros((nd, nd, nd))
    for dx in (0, 1):
        for dy in (0, 1):
            for dz in (0, 1):
                wt = w * (d[0] if dx else 1 - d[0]) * (d[1] if dy else 1 - d[1]) * (d[2] if dz else 1 - d[2])
                np.add.at(g, ((i0[0] + dx) % nd, (i0[1] + dy) % nd, (i0[2] + dz) % nd), wt)
    return np.fft.fftn(g) / nd**3


def binned_power(C, n, box, nbins):
    # mean of |C|^2 box^3 in bins nint(|m|) over the full sphere (weight 2 for iz > 0)
    import numpy as np
    ix, iy, iz = np.meshgrid(np.arange(-n, n + 1), np.arange(-n, n + 1), np.arange(0, n + 1), indexing="ij")
    mm = np.sqrt(ix**2 + iy**2 + iz**2)
    w = np.where(iz > 0, 2.0, 1.0)
    b = np.rint(mm).astype(int)
    p2 = np.abs(C) ** 2 * box**3
    return np.array([np.sum(w[b == k] * p2[b == k]) / np.sum(w[b == k]) for k in range(1, nbins + 1)])


def read_meta(path):
    out = {}
    for line in open(path):
        if "=" in line:
            k, v = line.split("=", 1)
            out[k.strip()] = v.strip()
    return out


# ---- test driver ----------------------------------------------------------
class Report:
    def __init__(self):
        self.failures = 0
        self.checks = 0

    def check(self, ok, case, msg):
        self.checks += 1
        if not ok:
            self.failures += 1
        print(f"{'PASS' if ok else 'FAIL'}  {case}: {msg}")


def oracle(cfg):
    """Returns (expected {step: aexpn}, xmarks, stop index, #inserted steps)."""
    return oracle_full(cfg)[:4]


def oracle_full(cfg):
    """As oracle(), plus the set of marked steps and the full ladder A."""
    a0, A, dA = build_ladder(cfg.get("da0", DA0))
    n0 = len(A)
    if cfg["exact"]:
        A, dA, marks, xmarks = apply_exact(a0, A, dA, cfg["zout"])
    else:
        marks, xmarks = legacy_marks(A, cfg["zout"]), {}
    stop = stop_index(A, dA, marks, exact=cfg["exact"])
    exp = expected_outputs(A, marks, cfg.get("nsteps", NSTEPS), stop)
    return exp, xmarks, stop, len(A) - n0, marks, A


def main(work):
    rep = Report()

    results = {}
    for case, cfg in SIM_CASES.items():
        rundir = os.path.join(work, case, "Run1")
        act = actual_outputs(rundir)
        results[case] = act
        exp, xmarks, _, _ = oracle(cfg)

        rep.check(set(act) == set(exp), case,
                  f"snapshot set {sorted(act)} == expected {sorted(exp)}")
        for step in sorted(set(act) & set(exp)):
            if step in xmarks:  # requested epoch: must be bit-exact
                rep.check(act[step] == xmarks[step], case,
                          f"step {step}: AEXPN={act[step]:.9g} bit-exact at requested 1/(1+z)={xmarks[step]:.9g}")
            else:               # scheduled epoch: <= 2 ulp of the ladder
                ok = abs(act[step] - exp[step]) <= 4e-7 * exp[step]
                rep.check(ok, case, f"step {step}: AEXPN={act[step]:.9g} matches schedule {exp[step]:.9g}")

    # OutputSchedule.dat (I4): exactly the marked steps, epochs as scheduled,
    # exact flags, comoving distances and classes (I5)
    for case, cfg in SIM_CASES.items():
        path = os.path.join(work, case, "Run1", "OutputSchedule.dat")
        if not os.path.exists(path):
            rep.check(False, case, "OutputSchedule.dat written")
            continue
        sched = read_schedule(path)
        _, xmarks, _, _, marks, A = oracle_full(cfg)
        rep.check(set(sched) == marks, case,
                  f"OutputSchedule.dat steps {sorted(sched)} == marked steps {sorted(marks)}")
        om, oml = setup_cosmology(os.path.join(work, case, "Setup.dat"))
        ok_a = ok_x = ok_chi = ok_cl = True
        zclass = {f32(1.0 / f32(1.0 + f32(z))): c for z, c in zip(cfg["zout"], cfg.get("classes", [1] * len(cfg["zout"])))}
        for step, (a, z, chi, ex, cl) in sched.items():
            if step not in marks:
                continue
            target = xmarks.get(step)
            ok_a &= (f32(a) == target) if target is not None else abs(a - A[step]) <= 4e-7 * A[step]
            ok_x &= ex == (1 if target is not None else 0)
            ok_chi &= abs(chi - chi_of(a, om, oml)) <= 1e-5 * max(chi, 1.0)
            expected_class = zclass.get(target, 1) if cfg["exact"] else 1
            ok_cl &= cl == expected_class
        rep.check(ok_a, case, "schedule epochs match the headers (bit-exact on pinned steps)")
        rep.check(ok_x, case, "exact flags mark exactly the pinned steps")
        rep.check(ok_chi, case, f"comoving distances match flat matter+Lambda (Om={om}, OmL={oml})")
        rep.check(ok_cl, case, f"output classes as requested ({cfg.get('classes', 'all 1')})")
    # lens planes (I2)
    import numpy as np
    planes = {}
    for case in ("case10_lens72", "case10_lens36"):
        cfg = SIM_CASES[case]
        nslab, nlens = cfg["lens"]
        _, xmarks, _, _, marks, _ = oracle_full(cfg)
        rundir = os.path.join(work, case, "Run1")
        ok_files = ok_meta = ok_mean = ok_proj = True
        for step in sorted(marks):
            fbin = os.path.join(rundir, f"LensPlanes.{step:04d}.0001.bin")
            fmeta = os.path.join(rundir, f"LensPlanes.{step:04d}.0001.meta")
            if not (os.path.exists(fbin) and os.path.exists(fmeta)):
                ok_files = False
                continue
            pl = read_planes(fbin, nslab, nlens)
            planes[(case, step)] = pl
            meta = read_meta(fmeta)
            ok_meta &= (int(meta["N_slab"]) == nslab and int(meta["N_lens"]) == nlens and int(meta["step"]) == step
                        and f32(float(meta["a"])) == xmarks[step] and abs(float(meta["depth"]) - 150.0 / nslab) < 1e-6)
            for ax in range(3):                      # mass conservation per axis
                ok_mean &= abs(pl[ax * nslab:(ax + 1) * nslab].mean() - 1.0) < 2e-6
            if nlens == NGRID:                       # L2: slab average == ProjectDENSIT full-box projection
                for ax, tag in enumerate(("YZ", "XZ", "XY")):
                    ref = read_projection(os.path.join(rundir, f"Projection{tag}.{step:04d}.0001.dat"), NGRID)
                    ok_proj &= np.max(np.abs(pl[ax * nslab:(ax + 1) * nslab].mean(axis=0) - ref)) < 1e-5
        rep.check(ok_files, case, "LensPlanes .bin and .meta written for every output")
        rep.check(ok_meta, case, "meta: N_slab, N_lens, step, a (bit-exact), depth")
        rep.check(ok_mean, case, "each axis: mean over its slabs of mean(1+delta) = 1 (mass conservation)")
        if nlens == NGRID:
            rep.check(ok_proj, case, "slab average equals the ProjectDENSIT full-box projection (< 1e-5)")
    # low-k modes (I8): files, meta, Hermitian layout, and P(k) from the density cube
    # against the code's own PowerDM.R, bin by bin (GetPower bins by nint(|m|), each
    # mode counted once over the full sphere: weight 2 for iz > 0, 1 for iz = 0)
    cfg = SIM_CASES["case11_modes"]
    kc, ndep = cfg["modes"]["kc"], cfg["modes"]["ndep"]
    n = int(kc * BOX / (2 * np.pi))
    _, xmarks, _, _, marks, _ = oracle_full(cfg)
    rundir = os.path.join(work, "case11_modes", "Run1")
    ok_files = ok_meta = ok_herm = ok_pk = ok_k0 = True
    for step in sorted(marks):
        for prod in ("deltam", "mom_x", "mom_y", "mom_z"):
            fb = os.path.join(rundir, f"{prod}.{step:04d}.0001.bin")
            fm = os.path.join(rundir, f"{prod}.{step:04d}.0001.meta")
            if not (os.path.exists(fb) and os.path.exists(fm)):
                ok_files = False
                continue
            C = read_modes(fb, n)
            meta = read_meta(fm)
            ok_meta &= (int(meta["n"]) == n and int(meta["N_dep"]) == (NGRID if prod == "deltam" else ndep)
                        and f32(float(meta["a"])) == xmarks[step] and float(meta["k_c"]) == kc)
            ok_herm &= hermitian_ok(C, n)
            ok_k0 &= C[n, n, 0] == 0
            if prod == "deltam":
                ix, iy, iz = np.meshgrid(np.arange(-n, n + 1), np.arange(-n, n + 1), np.arange(0, n + 1), indexing="ij")
                mm = np.sqrt(ix**2 + iy**2 + iz**2)
                w = np.where(iz > 0, 2.0, 1.0)
                b = np.maximum(np.rint(mm).astype(int), 1)
                p2 = np.abs(C) ** 2 * BOX**3
                pdm = read_powerdm(os.path.join(rundir, f"PowerDM.R.{step:04d}.0001.dat"))
                for bi in range(1, n):                      # complete bins only
                    sel = (b == bi) & (mm > 0)
                    pk_cube = np.sum(w[sel] * p2[sel]) / np.sum(w[sel])
                    pk_code = pdm[pdm[:, 0] == bi, 3][0]
                    ok_pk &= abs(pk_cube / pk_code - 1.0) < 1e-4
    rep.check(ok_files, "case11_modes", "deltam and mom_x/y/z .bin and .meta written for every output")
    rep.check(ok_meta, "case11_modes", f"meta: n={n}, N_dep, a (bit-exact), k_c")
    rep.check(ok_herm, "case11_modes", "iz = 0 plane is Hermitian (layout and sign convention)")
    rep.check(ok_k0, "case11_modes", "k = 0 mode is zero")
    rep.check(ok_pk, "case11_modes", "P(k) from the density cube equals PowerDM.R bin by bin (< 1e-4), V2")

    # subsample (I9) in full mode (M = 1: every particle, quantised): one file per class-1
    # output with n_sel = N, and V6: the density and momentum spectra of the decoded particles
    # equal the deltam/mom_* dumps at the five lowest bins to 1e-3 (quantisation only; the
    # hashed selection itself is covered by test/extras/run_subsample_test.sh)
    msub = cfg["subm"]
    npart = 36**3
    ok_sfiles = ok_nsel = ok_v6 = True
    for step in sorted(marks):
        fb = os.path.join(rundir, f"Subsample.{step:04d}.0001.bin")
        fm = os.path.join(rundir, f"Subsample.{step:04d}.0001.meta")
        if not (os.path.exists(fb) and os.path.exists(fm)):
            ok_sfiles = ok_nsel = ok_v6 = False
            continue
        meta, x, v = read_subsample(fb, fm)
        nsel = x.shape[1]
        ok_nsel &= (nsel == npart if msub == 1 else abs(nsel - npart / msub) < 3 * np.sqrt(npart / msub)) \
                   and int(meta["M"].split()[0]) == msub
        Cd = cube_from_particles(x, np.full(nsel, NGRID**3 / nsel), NGRID, BOX)   # deltam comes from the Ngrid mesh
        Csub = np.zeros((2 * n + 1, 2 * n + 1, n + 1), complex)
        for i in range(-n, n + 1):
            for j in range(-n, n + 1):
                Csub[i + n, j + n, :] = Cd[i, j, 0:n + 1]
        pd_sub = binned_power(Csub, n, BOX, 5)
        pd_full = binned_power(read_modes(os.path.join(rundir, f"deltam.{step:04d}.0001.bin"), n), n, BOX, 5)
        ok_v6 &= np.all(np.abs(pd_sub / pd_full - 1.0) < 1e-3)
        Cm = cube_from_particles(x, v[0] * ndep**3 / nsel, ndep, BOX)
        Cmsub = np.zeros((2 * n + 1, 2 * n + 1, n + 1), complex)
        for i in range(-n, n + 1):
            for j in range(-n, n + 1):
                Cmsub[i + n, j + n, :] = Cm[i, j, 0:n + 1]
        pm_sub = binned_power(Cmsub, n, BOX, 5)
        pm_full = binned_power(read_modes(os.path.join(rundir, f"mom_x.{step:04d}.0001.bin"), n), n, BOX, 5)
        ok_v6 &= np.all(np.abs(pm_sub / pm_full - 1.0) < 1e-3)
    rep.check(ok_sfiles, "case11_modes", "Subsample .bin and .meta written for every output")
    rep.check(ok_nsel, "case11_modes", f"n_sel = N (M={msub}, full quantised snapshot)")
    rep.check(ok_v6, "case11_modes", "decoded particles reproduce the deltam and mom_x spectra at k <= 5 k_f to 1e-3 (V6)")

    # f_R at halo centres (I10): one float32 pair per catalogue row, finite, |f_R| <= |fR_bg| * max(u)^2
    _, _, _, _, marks12, _ = oracle_full(SIM_CASES["case12_frhalo"])
    rundir = os.path.join(work, "case12_frhalo", "Run1")
    ok_fr = ok_rows = True
    nrows_total = 0
    for step in sorted(marks12):
        fb = os.path.join(rundir, f"fR_halo.{step:04d}.0001.bin")
        fm = os.path.join(rundir, f"fR_halo.{step:04d}.0001.meta")
        cat = os.path.join(rundir, "CATALOGS", f"CatshortV.{step:04d}.0001.DAT")
        if not (os.path.exists(fb) and os.path.exists(fm) and os.path.exists(cat)):
            ok_fr = False
            continue
        vals = np.fromfile(fb, dtype="<f4").reshape(-1, 2)
        meta = read_meta(fm)
        ncat = 0
        for line in open(cat):                       # data rows: >= 20 numeric fields
            p = line.split()
            try:
                float(p[0]); float(p[1])
                ncat += len(p) >= 20
            except (ValueError, IndexError):
                pass
        ok_rows &= vals.shape[0] == ncat == int(meta["rows"])
        nrows_total += vals.shape[0]
        if vals.size:
            fRbg = float(meta["fR_bg(a_old)"])
            ok_fr &= bool(np.all(np.isfinite(vals))) and bool(np.all(vals <= 0)) and bool(np.all(np.abs(vals) < 10 * abs(fRbg)))
    rep.check(ok_fr, "case12_frhalo", f"fR_halo files per output, finite, sign of fR0 (rows in total: {nrows_total})")
    rep.check(ok_rows, "case12_frhalo", "fR_halo rows == CatshortV rows == meta rows")

    # 36 = 2x2 block mean of 72 at the same steps (same IC, same schedule)
    ok_blk = True
    for (case, step), pl in planes.items():
        if case != "case10_lens72":
            continue
        pl36 = planes.get(("case10_lens36", step))
        if pl36 is None:
            ok_blk = False
            continue
        blk = pl.reshape(pl.shape[0], 36, 2, 36, 2).mean(axis=(2, 4))
        ok_blk &= np.max(np.abs(blk - pl36)) < 1e-5
    rep.check(ok_blk, "case10_lens36", "planes equal the 2x2 block means of the full-resolution planes")

    # L3: plane power spectra vs the 3-D P(k) of the same output through the exact projection
    # formula of the periodic box (Limber-free): P_2D(k_perp) = (1/L) sum_kz P_3D(|k|) |W(kz)|^2 with
    # the discrete slab top-hat W(kz) = sin(kz nc dz/2) / (nc sin(kz dz/2)), nc cells of size dz.
    # PowerDM.R and the planes are both raw CIC quantities and the CIC window is isotropic to
    # O(k^2 dx^2), so neither side is deconvolved. The comparison is statistical: a 2-D ring holds a
    # few modes while the 3-D bins average over shells, so each ring must agree within 3 sigma
    # (sigma = 1/sqrt(N_eff), N_eff = 3 axes x independent modes in the ring) and the mean ratio over
    # 0.08 < k_perp < 0.5 h/Mpc (k_N = 1.5) within 10 per cent.
    ok_l3 = ok_l3mean = True
    nslab, nlens = SIM_CASES["case10_lens72"]["lens"]
    rundir = os.path.join(work, "case10_lens72", "Run1")
    kf = 2 * np.pi / BOX
    nc = NGRID // nslab
    dz = BOX / NGRID
    mz = np.fft.fftfreq(NGRID) * NGRID
    kz = kf * mz
    with np.errstate(divide="ignore", invalid="ignore"):
        W2 = np.where(mz == 0, 1.0, (np.sin(kz * nc * dz / 2) / (nc * np.sin(kz * dz / 2))) ** 2)
    mx, my = np.meshgrid(np.fft.fftfreq(nlens) * nlens, np.fft.fftfreq(nlens) * nlens, indexing="ij")
    mperp = np.sqrt(mx ** 2 + my ** 2)
    bperp = np.rint(mperp).astype(int)
    for (case, step), pl in sorted(planes.items()):
        if case != "case10_lens72":
            continue
        fp = os.path.join(rundir, f"PowerDM.R.{step:04d}.0001.dat")
        if not os.path.exists(fp):
            ok_l3 = ok_l3mean = False
            continue
        pdm = read_powerdm(fp)
        kt, pt = pdm[:, 1], pdm[:, 3]
        p3d = lambda k: np.exp(np.interp(np.log(np.clip(k, kt[0], kt[-1])), np.log(kt), np.log(pt)))
        p2d = np.zeros((nlens, nlens))
        for p in pl:                                    # N^-2 convention, P_2D = L^2 <|delta_2D|^2>
            p2d += np.abs(np.fft.fft2(p - 1.0) / nlens ** 2) ** 2 * BOX ** 2
        p2d /= len(pl)
        ratios = []
        for b in range(2, 13):                          # 0.08 < k_perp < 0.5 h/Mpc
            sel = bperp == b
            pred = np.mean([np.sum(p3d(np.sqrt(kp ** 2 + kz ** 2)) * W2) / BOX for kp in kf * mperp[sel]])
            ratio = p2d[sel].mean() / pred
            neff = 3 * sel.sum() / 2
            sig = 1.0 / np.sqrt(neff)
            ratios.append(ratio)
            if abs(ratio - 1) >= 3 * sig:
                ok_l3 = False
                print(f"      L3 step {step} ring b={b} (k_perp={kf * b:.3f}): P_2D/pred = {ratio:.3f}, sigma = {sig:.3f}")
        print(f"      L3 step {step}: P_2D/prediction over the {len(ratios)} rings: mean {np.mean(ratios):.3f}, "
              f"min {min(ratios):.3f}, max {max(ratios):.3f}")
        ok_l3mean &= abs(np.mean(ratios) - 1) < 0.10
    rep.check(ok_l3, "case10_lens72", "L3: plane P_2D(k_perp) within 3 sigma of the projected 3-D P(k) in every ring, 0.08 < k_perp < 0.5")
    rep.check(ok_l3mean, "case10_lens72", "L3: mean P_2D/prediction over the rings within 10 per cent of 1")

    # case 9: classes must not change the schedule of case 3 (same redshifts)
    rep.check(set(results["case9_classes"]) == set(results["case3_multi"]), "case9_classes",
              "snapshot set identical to case3_multi (classes do not alter the schedule)")

    # case 1: the near-coincident target is merged (only z=2 inserts a step),
    # and exactly one snapshot exists at the merged moment
    _, _, _, nins1 = oracle(SIM_CASES["case1_merge"])
    rep.check(nins1 == 1, "case1_merge",
              f"z=0.49337 merged onto the existing step (only 1 step inserted, got {nins1})")
    act = results["case1_merge"]
    near = [s for s, a in act.items() if abs(a - f32(1.0 / f32(1.0 + f32(0.49337)))) < 0.005]
    rep.check(len(near) == 1, "case1_merge", f"exactly one snapshot at the merged moment (got {len(near)})")

    # case 2: both targets are genuine insertions
    _, _, _, nins2 = oracle(SIM_CASES["case2_insert"])
    rep.check(nins2 == 2, "case2_insert",
              f"z=2.0 and z=0.512 both inserted as new steps (got {nins2})")

    # case 3: all targets present, ordered (ascending step = descending z)
    _, xm3, _, _ = oracle(SIM_CASES["case3_multi"])
    steps3 = sorted(xm3)
    rep.check(len(steps3) == 3, "case3_multi", f"all 3 requested moments scheduled: steps {steps3}")
    rep.check([xm3[s] for s in steps3] == sorted(xm3.values()), "case3_multi",
              "requested snapshots ordered: step number increases as z decreases")

    # case 6: z=0 inserted AFTER the legacy stop step must still be executed
    _, xm6, stop6, _ = oracle(SIM_CASES["case6_tailinsert"])
    z0step = max(xm6)
    rep.check(stop6 == z0step, "case6_tailinsert",
              f"stop index extends to the inserted z=0 step ({z0step})")
    act6 = results["case6_tailinsert"]
    rep.check(z0step in act6 and act6.get(z0step) == 1.0, "case6_tailinsert",
              f"z=0 snapshot produced at step {z0step} with AEXPN exactly 1.0")

    # cases 7/8: close pairs -> each member pinned on its own step, the
    # schedule strictly increasing, and the run really executed every
    # (short) step in the oracle's order
    for case in ("case7_closepairs", "case8_pinned"):
        cfg = SIM_CASES[case]
        a0, A, dA = build_ladder(DA0)
        A, dA, marks, xm = apply_exact(a0, A, dA, cfg["zout"])
        rep.check(len(xm) == len(cfg["zout"]), case,
                  f"all {len(cfg['zout'])} requested moments pinned on distinct steps (got {len(xm)})")
        gaps = [f32(A[i] - A[i - 1]) for i in range(2, len(A))]
        rep.check(min(gaps) > 0 and min(dA[1:]) > 0, case,
                  f"schedule strictly increasing, shortest step da={min(dA[1:]):.3e}")
        stop = stop_index(A, dA, marks, exact=True)
        last = min(cfg.get("nsteps", NSTEPS), stop)
        seq = [float(v) for v in a_sequence(os.path.join(work, case, "Run1", "main.log"))]
        rep.check(len(seq) == last, case, f"run executed {len(seq)} steps, oracle expects {last}")
        ok = len(seq) == last and all(abs(seq[i] - A[i + 1]) <= LOG_ATOL for i in range(last))
        rep.check(ok, case, "logged expansion factors follow the oracle schedule step by step")
    _, xm8, _, nins8 = oracle(SIM_CASES["case8_pinned"])
    rep.check(nins8 == 1, "case8_pinned",
              f"z=0.49337 merged, z=0.4920 inserted (1 step inserted, got {nins8})")
    rep.check(len(xm8) == 2 and sorted(xm8.values())[0] == f32(1.0 / f32(1.0 + f32(0.49337))),
              "case8_pinned", "the pinned step stayed on its own target")

    # case 5: backward compatibility
    trailing = open(os.path.join(work, "case5b_trailing", "Setup.dat")).read()
    rep.check("exact output redshifts" in trailing.splitlines()[-2], "case5b_trailing",
              "Setup.dat really carries stale trailing lines")
    rep.check(set(results["case5b_trailing"]) == set(results["case5a_legacy"]),
              "case5b_trailing", "output moments identical despite trailing lines")
    seq_a = a_sequence(os.path.join(work, "case5a_legacy", "Run1", "main.log"))
    seq_b = a_sequence(os.path.join(work, "case5b_trailing", "Run1", "main.log"))
    rep.check(len(seq_a) == NSTEPS and seq_a == seq_b, "case5b_trailing",
              "expansion-factor sequence identical with and without trailing lines")

    # case 4: invalid requests -> clear message, no products
    for case, (log, msg, product) in ERR_CASES.items():
        text = ""
        logpath = os.path.join(work, case, log)
        if os.path.exists(logpath):
            text = open(logpath, errors="replace").read()
        rep.check(msg in text, case, f'error message "{msg}" reported')
        leftovers = glob.glob(os.path.join(work, case, product))
        rep.check(not leftovers, case, f"no {product} produced")

    print(f"\n{rep.checks - rep.failures}/{rep.checks} checks passed"
          + ("" if rep.failures == 0 else f"  ({rep.failures} FAILED)"))
    return 1 if rep.failures else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "work"))
