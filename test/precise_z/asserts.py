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
ZTOL = 0.005        # merge tolerance in SetExactSteps

SIM_CASES = {       # case dir -> configuration (zout list + sign of #outputs)
    "case5a_legacy": dict(zout=[2.00, 0.50], exact=False),
    "case5b_trailing": dict(zout=[2.00, 0.50], exact=False),
    "case1_merge": dict(zout=[2.00, 0.49337], exact=True),
    "case2_insert": dict(zout=[2.00, 0.512], exact=True),
    "case3_multi": dict(zout=[2.5, 0.512, 1.1], exact=True),
    # z=0 must be INSERTED after the step where the legacy half-step exit
    # rule already fires (regression test for the NlastX stop-index fix)
    "case6_tailinsert": dict(zout=[0.0], exact=True, da0=6.9e-4, nsteps=300),
}
ERR_CASES = {       # case dir -> (log file, required message, forbidden product)
    "case4a_zbig": ("init.log", "outside the simulation range", "Setup.dat"),
    "case4b_zneg": ("init.log", "outside the simulation range", "Setup.dat"),
    "case4c_pair": ("Run1/main.log", "closer than the 0.5% merge tolerance", "Run1/PMcrd.0*.DAT"),
    "case4d_range": ("Run1/main.log", "at or before the initial redshift", "Run1/PMcrd.0*.DAT"),
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
        if abs(f32(A[jn] - at)) < atol:     # merge: move the scheduled step
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
    a0, A, dA = build_ladder(cfg.get("da0", DA0))
    n0 = len(A)
    if cfg["exact"]:
        A, dA, marks, xmarks = apply_exact(a0, A, dA, cfg["zout"])
    else:
        marks, xmarks = legacy_marks(A, cfg["zout"]), {}
    stop = stop_index(A, dA, marks, exact=cfg["exact"])
    exp = expected_outputs(A, marks, cfg.get("nsteps", NSTEPS), stop)
    return exp, xmarks, stop, len(A) - n0


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
