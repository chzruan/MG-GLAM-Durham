# z=0.25 campaign gate and cost estimate

Scope: nominal z=0.25 / snapshot 137, LCDM fiducial imodel 0 boxes 1--100,
and imodel 1--64, boxes 1--5, for both LCDM and fRn1. All 740 ASCII inputs
exist. Their total size is 708.165 GB. Other redshifts are inventoried but
are not part of this campaign; the explicit scope question received no reply
before the z=0.25 assumption was announced.

The full write-and-validate box 1 pilot, job **11942890**, passed on
cosma8-serial/dp004: one CPU, 2 GiB, ten-minute limit; elapsed 22 s,
TotalCPU 20.796 s (94.5% CPU efficiency), Slurm MaxRSS 1334.73 MiB.
ReqTRES and AllocTRES both show **one billed CPU**. Actual core-hours 0.00611,
time-limit core-hours 0.1667. Its 3,378,809 retained rows are byte-identical
to the corresponding input lines in all 24 columns. Both retained-row SHA256
streams equal `ce7147a2c413fa855cae92eccea2941a6a532943d16c397c7e046264385407f6`.
The pilot output is preserved and verified on campaign resume.

Before submission, input-byte scaling of the 22 s pilot gives **4.52 core-hours**
for 740 catalogues, or **34 minutes at eight workers**. Allow 4.5--9 core-hours
and roughly 34--68 minutes for I/O variation. Eight independent one-core
workers, each processing a fixed strided subset sequentially, use the shared
queue with array concurrency capped at eight. The two-hour limit per worker
caps the campaign at **16 time-limit core-hours**. No exclusive nodes are used.

The largest input is 1.741 GB, 1.817 times box 1. Scaling measured Slurm MaxRSS
by that ratio and 1.5 headroom gives 3.55 GiB, rounded up to an explicit
**4 GiB per worker**. Peak requested aggregate resources are eight CPUs and
32 GiB. The array organizes eight workers; it does not request full-node
allocations. Completed outputs are never overwritten during restarts.

Outputs mirror source paths below `products/strict-v1/`: fiducial files retain
`DESI_MGx100/data/GR/Run<ibox>/CATALOGS/`; training files retain their
`mg_glam/DurMun_hmfemu_<gravity>_...model<imodel>.../Run<ibox>/CATALOGS/`
structure. Each output has one adjacent `.cleaning.hdf5` sidecar, containing
the full-length row-drop masks, original-row/ID groups, validation and receipt.
This uses two files per catalogue rather than separate files for each metadata
component. The expected new ASCII volume is approximately 707 GB, plus small
compressed sidecars. No raw source or existing mirror is modified.

The actual job IDs, exact script and tool hashes, dependencies, sacct accounting
and achieved utilisation are recorded in `provenance/`. The pilot gate requires
successful exit, verified hashes, bitwise preservation, and zero remaining
**strict-rule** edges. Ambiguous same-bound-count pairs are reported, not treated
as a failed job or deleted merely to make their count zero.

## Completed campaign accounting

Array **11942910**, submitted after the successful cleaning pilot with an
`afterok:11942890` dependency, completed all eight tasks with exit 0:0.
The longest worker ran 4,089 s (68 min 9 s), versus the approximate pre-launch
34--68 minute wall-time range. Actual array allocation/billing was **7.8864
core-hours**, within the 4.5--9 core-hour allowance and below the 16 core-hour
limit. TotalCPU was 4.4277 hours, giving 56.14% aggregate CPU efficiency.
I/O stalls account for much of the difference from the short pilot; a process
sample and concurrent Slurm statistics are retained in provenance.

Slurm batch-step MaxRSS ranged from 4083.95 to 4084.00 MiB under the explicit
4 GiB requests; this value is distinct from sampled process RSS and the
per-worker process high-water marks in catalogue receipts. All tasks completed
without OOM or restart. ReqTRES/AllocTRES confirm one allocated and billed CPU
per worker, and no exclusive-node allocation.

All 740 products passed the final `--require-complete` summary, with zero
missing/invalid receipts and zero surviving strict-rule edges. The pre-existing
pilot product was verified and reused. Final output volume is 706.890 GB of
ASCII plus 216.204 MB of compressed sidecars. Including audit pilot 11942773,
cleaning pilot 11942890 and policy replay 11942925, all eleven allocations used
**7.9003 billed core-hours**, compared with 16.5 time-limit core-hours.
See `provenance/all-jobs.sacct.txt`, `provenance/accounting_summary.json` and
[RESULTS.md](RESULTS.md) for the final tables and scientific qualifications.
