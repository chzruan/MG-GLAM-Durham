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
sample. Production pass distributions will inform the follow-up limit policy.

`integration-results.json` repeats the 12 comparisons on SO-v3 against the
normalization-only commit `e2a327a`, which has the same physics and no counters.
The test now defaults to that baseline; the initial v2 receipt still records
its own source and baseline hashes. Both runs preserve properties and IDs.
The extracted counter fixture does not exercise `ReleaseMaxima`; the separate
full native two-call preflight and snapshot wrapper explicitly assert release
of the diagnostic arrays and the original particle workspace.
