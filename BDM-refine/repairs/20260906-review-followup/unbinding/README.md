# Unbinding work diagnostics

The finder records the number of unbinding passes and the sum of the active
particle population presented to those passes for every candidate. The latter
is a work proxy, including a pass rejected for a singular potential; it is not
an instruction count. Candidates rejected before unbinding have both counters
zero. Each direct recall clears its counters, and normal finder teardown frees
both arrays and balances their 12 bytes per candidate in memory accounting.

The run log reports the number of candidates entering unbinding, maximum pass
count, total active rows across passes and a histogram with bins 1, 2, 3–4,
5–8, 9–16, 17–32 and >32. The five existing status counts are preserved.
These measurements do not truncate convergence or change halo selection.

GNU checked and optimized comparisons against the uninstrumented source cover
an empty sphere, one-pass binding, complete unbinding, a singular centre and
the independent review's adversarial energy ladders. All properties and exact
memberships agree. The 200/400-row ladders retain 20 rows and take 91/191
passes, with 10,010/40,110 active rows presented across those passes. Every
case also checks that a later empty recall clears the counters.

`diagnostics-results.json` records 12 comparisons (24 native executions),
source/compiler/test hashes and commands. Python must be invoked through
`micromamba run -n cosemu python3 -B`.

The quadratic worst case is real. A limit that silently accepts an unfinished
population or drops an otherwise converging halo would alter the scientific
sample. The N1024 production replay below establishes that even published
haloes can need more than 32 passes, so this repair retains exact convergence
and reports its cost. It does not introduce a pass cap or claim to remove the
quadratic worst case. An exact accelerated algorithm would need separate
membership/physics validation against these controls.

`integration-results.json` repeats the 12 comparisons on SO-v3 against the
normalization-only commit `e2a327a`, which has the same physics and no counters.
The test now defaults to that baseline; the initial v2 receipt still records
its own source and baseline hashes. Both runs preserve properties and IDs.
The extracted counter fixture does not exercise `ReleaseMaxima`; the separate
full native two-call preflight and snapshot wrapper explicitly assert release
of the diagnostic arrays and the original particle workspace.

## Measured production tail

The final v3 finder on the 1024^3-particle, 2048^3-mesh snapshots gives:

| z | Candidates entering unbinding | 99th-percentile passes | Maximum passes | Candidates >32 | Published candidates >32 | Maximum published passes | Total work in candidates >32 / all work |
|---|---:|---:|---:|---:|---:|---:|---:|
| 2 | 656,674 | 6 | 58 | 2 | 0 | 32 | 0.119% |
| 1 | 870,892 | 7 | 97 | 63 | 7 | 40 | 1.041% |
| 0 | 623,054 | 9 | 144 | 265 | 7 | 48 | 6.249% |

The work fraction includes **every pass of candidates whose final pass count
exceeds 32**; it is not the work performed after pass 32 and cannot be used as
the prospective saving of a 32-pass cap. Total active-row work is 104,387,091,
502,120,464 and 2,084,143,350 at z=2, 1, 0. Maximum retained bound-set sizes
are 13,081, 37,628 and 183,465 rows respectively; a host/duplicate-rejected
candidate can retain IDs while its post-selection mass is zero.

The same-density v2 reference reaches 209 passes at z=0, with 278 candidates
above 32 and one published halo needing 41. The optimized-v2 computational
control has byte-identical unbinding diagnostics, catalogue and member tape.
Thus the measured v2/v3 tail changes are consequences of the normalization
change rather than an altered unbinding stopping rule. Both versions' full
histograms, percentiles and highest-work candidates are retained in
[`native/results-t64.json`](../native/results-t64.json); the independent
32-thread z=0 controls are in [`native/results-t32.json`](../native/results-t32.json).
