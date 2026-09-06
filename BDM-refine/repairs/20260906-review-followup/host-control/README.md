# Host-exclusion checker: full-catalogue positive controls

Job 11949318 applied the original production membership checker, unchanged, to
all three original catalogues and then to a temporary copy with one invalid
host pair injected at each redshift. All three originals passed; all three
corrupt copies failed the intended host-pair assertion. The injection moved a
lower-mass selected halo to the highest-mass host centre, changing the raw
coordinate fields and matching catalogue fields together so the check reached
the host rule. Original raw tapes, indexes and catalogue arrays retained their
SHA-256 hashes. Temporary large copies were removed.

The original zero count of higher-priority neighbours was a valid null result,
but did not by itself exercise rejection. These production-size positive
controls demonstrate that the checker can reject the specified violation;
they do not establish completeness against every possible catalogue defect.

`results.json` records each pair and hashes. `submission.json`,
`resource-plan.json` and `accounting.json` retain launch and resource evidence.
The job used one shared COSMA8 core for 80 s (0.02223 allocated core-hours),
63.297 CPU s. The process reported 445 MiB peak RSS; Slurm charged 4084 MiB,
including file cache, close to its 4 GiB memory request. Use 6 GiB for a repeat.
