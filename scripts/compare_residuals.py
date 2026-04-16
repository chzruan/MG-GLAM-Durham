#!/usr/bin/env python3
"""
Task A parser: compare per-V-cycle residuals across OMP thread counts.

Input: one or more log files of the form verify/race_T{N}.log — each is the
full stdout of a PMP2main-bitmatch.exe run at OMP_NUM_THREADS=N. One file
must be the T=1 reference; that run's residual sequence is the baseline
all others are compared against.

Signal: PMP2extradof.f90:105 emits
    'The residual is equal to <F20.12> after <I5> V-cycles.'
one line per V-cycle. We extract the (residual, iter_count) sequence from
each log, index by (step_idx, iter_count), and compute max|Δ| and
max|Δ|/|ref| against the T=1 reference.

Pass bar: max|Δ| ≤ 1e-14. F20.12 caps printed precision at ~1e-12, so
"exact bit-match" shows up here as max|Δ| == 0.0. Any non-zero delta at
the printed precision flags a potential data race.

Exit 0 on pass, 1 on fail, 2 on parse/shape error.
"""

import argparse
import re
import sys
from pathlib import Path

RESID_RE = re.compile(
    r"The residual is equal to\s+([-+\d.Ee]+)\s+after\s+(\d+)\s+V-cycles\."
)

THREADS_RE = re.compile(r"race_T(\d+)\.log$")

BIT_MATCH_TOL = 1e-14


def parse_log(path: Path):
    """Return list of (step_idx, iter_count, residual).

    step_idx increments every time iter_count resets to 1 (each scalar_field
    call begins a fresh MG solve with iter_count=1 on the first V-cycle).
    """
    seq = []
    step_idx = -1
    prev_iter = None
    with path.open("r", errors="replace") as fh:
        for line in fh:
            m = RESID_RE.search(line)
            if not m:
                continue
            res = float(m.group(1))
            it = int(m.group(2))
            if prev_iter is None or it <= prev_iter:
                step_idx += 1
            seq.append((step_idx, it, res))
            prev_iter = it
    return seq


def threads_from_name(path: Path) -> int:
    m = THREADS_RE.search(path.name)
    if not m:
        raise ValueError(f"cannot extract thread count from filename: {path.name}")
    return int(m.group(1))


def compare(ref_seq, other_seq):
    """Pair by (step_idx, iter_count); return (max_abs, max_rel, pairs, mismatches)."""
    ref_map = {(s, i): r for (s, i, r) in ref_seq}
    other_map = {(s, i): r for (s, i, r) in other_seq}
    common = sorted(set(ref_map) & set(other_map))
    max_abs = 0.0
    max_rel = 0.0
    max_abs_key = None
    for k in common:
        d = abs(other_map[k] - ref_map[k])
        if d > max_abs:
            max_abs = d
            max_abs_key = k
        ref_val = abs(ref_map[k])
        if ref_val > 0.0:
            r = d / ref_val
            if r > max_rel:
                max_rel = r
    only_ref = set(ref_map) - set(other_map)
    only_other = set(other_map) - set(ref_map)
    mismatches = (only_ref, only_other)
    return max_abs, max_rel, len(common), mismatches, max_abs_key


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("logs", nargs="+", type=Path,
                    help="race_T{N}.log files; one must be T=1.")
    ap.add_argument("--tol", type=float, default=BIT_MATCH_TOL,
                    help=f"pass bar on max|Δres| (default {BIT_MATCH_TOL:g})")
    args = ap.parse_args()

    by_threads = {}
    for p in args.logs:
        t = threads_from_name(p)
        if t in by_threads:
            print(f"ERROR: duplicate thread count {t} in inputs", file=sys.stderr)
            return 2
        by_threads[t] = p

    if 1 not in by_threads:
        print("ERROR: need race_T1.log as the reference", file=sys.stderr)
        return 2

    ref_path = by_threads[1]
    ref_seq = parse_log(ref_path)
    if not ref_seq:
        print(f"ERROR: no residual lines parsed from {ref_path}", file=sys.stderr)
        return 2

    print(f"# Task A: residual bit-match regression")
    print(f"# reference: T=1 ({ref_path.name}), {len(ref_seq)} V-cycles")
    print(f"# tolerance: max|Δres| ≤ {args.tol:g}")
    print()
    print("| T | V-cycles | max|Δ| | max|Δ|/|ref| | status |")
    print("|---|---------:|-------:|-------------:|:------:|")

    worst_abs = 0.0
    failed = False
    for t in sorted(by_threads):
        if t == 1:
            print(f"| 1 | {len(ref_seq)} | 0 | 0 | ref |")
            continue
        seq = parse_log(by_threads[t])
        if not seq:
            print(f"| {t} | 0 | — | — | FAIL (empty) |")
            failed = True
            continue
        max_abs, max_rel, n_common, (only_ref, only_other), worst_key = compare(ref_seq, seq)
        ok = max_abs <= args.tol and not only_ref and not only_other
        status = "PASS" if ok else "FAIL"
        if not ok:
            failed = True
        worst_abs = max(worst_abs, max_abs)
        print(f"| {t} | {n_common} | {max_abs:.3e} | {max_rel:.3e} | {status} |")
        if only_ref or only_other:
            print(f"  # shape mismatch at T={t}: "
                  f"ref-only={len(only_ref)} other-only={len(only_other)}",
                  file=sys.stderr)
        if not ok and worst_key is not None:
            print(f"  # worst |Δ| at T={t}: step={worst_key[0]} iter={worst_key[1]} "
                  f"ref={dict(((s,i),r) for s,i,r in ref_seq)[worst_key]:.12e} "
                  f"other={dict(((s,i),r) for s,i,r in seq)[worst_key]:.12e}",
                  file=sys.stderr)

    print()
    print(f"# overall max|Δ| = {worst_abs:.3e}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
