#!/usr/bin/env python3
"""
Task C parser: extract per-step MG relaxation time from timing.log copies.

Input: files named verify/bench_T{N}_r{R}.log — each is a preserved copy of
one run's timing.log (PMP2mod_tools.f90:513–527, fmt '1p,12G13.4').

timing.log shape:
    line 1: 'Ngrid = ... Npart = ... Nthreads = ...'
    line 2: 'Step   Tot/min   Force   Move   Density   IO   Bias   BDM   Analysis'
    line 3+: '  <iStep>   <tot>   <force>   <move>   <density>  <io> <bias> <bdm> <analysis>  <mg>'
where column index 9 (after whitespace-split, 0-indexed) of the data row is
the MG column — the total time per step charged by every MG solver via
TimingMain(8,1) in PMP2mod_tools.f90:510. Values are MINUTES per step
(CPU/60. in the write format at line 522), so the reported number must be
labelled accordingly.

Processing:
  1. For each run, drop step 1 (warmup), take median of remaining steps'
     Density values.
  2. For each thread count T, take median across repeats.
  3. Emit markdown table with speedup vs T=16 (plan's baseline).

Usage:
    python3 bench_report.py verify/bench_T16_r1.log verify/bench_T16_r2.log ... > verify/bench_report.md
"""

import argparse
import re
import statistics
import sys
from collections import defaultdict
from pathlib import Path

FILENAME_RE = re.compile(r"bench_T(\d+)_r(\d+)\.log$")

# timing.log column 9 (0-indexed) of the data row is the MG solver time.
MG_COL = 9


def parse_timing_log(path: Path):
    """Return list of MG solver values (minutes/step), one per timestep."""
    mg = []
    with path.open("r", errors="replace") as fh:
        for line in fh:
            parts = line.split()
            if not parts:
                continue
            try:
                int(parts[0])
            except ValueError:
                continue
            if len(parts) <= MG_COL:
                continue
            try:
                d = float(parts[MG_COL])
            except ValueError:
                continue
            mg.append(d)
    return mg


def threads_repeat_from_name(path: Path):
    m = FILENAME_RE.search(path.name)
    if not m:
        raise ValueError(f"cannot parse T/r from filename: {path.name}")
    return int(m.group(1)), int(m.group(2))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("logs", nargs="+", type=Path,
                    help="verify/bench_T{N}_r{R}.log files")
    ap.add_argument("--baseline", type=int, default=16,
                    help="thread count used as the speedup baseline (default 16)")
    ap.add_argument("--drop-warmup", action="store_true", default=True,
                    help="drop step 1 from each run (default on)")
    args = ap.parse_args()

    # by_T: dict[int T] -> list of per-run medians
    by_T = defaultdict(list)
    # sample_counts to spot runs that aborted early
    sample_counts = defaultdict(list)

    for p in args.logs:
        try:
            T, r = threads_repeat_from_name(p)
        except ValueError as e:
            print(f"# WARN: {e}", file=sys.stderr)
            continue
        d = parse_timing_log(p)
        if not d:
            print(f"# WARN: no data rows in {p}", file=sys.stderr)
            continue
        sample_counts[T].append(len(d))
        if args.drop_warmup and len(d) > 1:
            d = d[1:]
        med = statistics.median(d)
        by_T[T].append((r, med))

    if not by_T:
        print("ERROR: no valid runs parsed", file=sys.stderr)
        return 2

    # Per-T median across repeats.
    per_T_med = {T: statistics.median(m for _, m in rs) for T, rs in by_T.items()}

    baseline = per_T_med.get(args.baseline)
    if baseline is None:
        print(f"# WARN: baseline T={args.baseline} missing; speedups skipped",
              file=sys.stderr)

    print("# Task C: MG relaxation benchmark (fid_F5, production PMP2main.exe)")
    print("# Column: median(MG solver) [min/step], step 1 dropped; outer median over repeats.")
    print()
    print("| T | repeats | steps/run | median T(MG) [min] | min | max | speedup vs T={} |".format(args.baseline))
    print("|---|--------:|----------:|----------------------:|----:|----:|----------------:|")
    for T in sorted(by_T):
        meds = [m for _, m in by_T[T]]
        med_of_meds = statistics.median(meds)
        smin = min(meds)
        smax = max(meds)
        steps = ",".join(str(s) for s in sample_counts[T])
        if baseline is not None and med_of_meds > 0:
            sp = baseline / med_of_meds
            sp_str = f"{sp:.2f}×"
        else:
            sp_str = "—"
        print(f"| {T} | {len(meds)} | {steps} | {med_of_meds:.4g} | {smin:.4g} | {smax:.4g} | {sp_str} |")

    print()
    print("# Raw per-repeat medians:")
    for T in sorted(by_T):
        per_run = ", ".join(f"r{r}={m:.4g}" for r, m in sorted(by_T[T]))
        print(f"# T={T}: {per_run}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
