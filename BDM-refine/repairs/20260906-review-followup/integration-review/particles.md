# Independent particle/replay integration review

Reviewed `c68c22dfa7e312334d504da6d8efcb8cad7ba7a7` (PMP2linker SHA256
`39f7e514d1db9fb3b9c7545fc7f3269b4e20b5ee2d1366945f8916b997e4cd5d`).
No blocker was found in the controlled replay scripts, cache admission, or
telemetry allocation/release integration. No production source changed.

The current four driver controls pass. The historical fixed-density reference
is the intended d64 catalogue (147,042 rows); timing extraction uses the outer
stage label, and final density verification performs an uncached hash.

A separate GNU checked fixture uses the actual `Memory`, `List`,
`BdmHaloMembershipInit` and `ReleaseMaxima` routines. Default and near-limit cache
admission pass, insufficient headroom fails explicitly, cache release preserves
the accounting value and link order, and three telemetry allocation/reinit/release
cycles clear counters and return to the original memory value.

At the saved z=0 counts, the cache is 2.18610 GiB. The accounted components during
list construction sum to about 78.822 GiB against the 500 GiB default; telemetry
adds 9,534,840 bytes later. These figures are bookkeeping, not a complete RSS cap.
The replay wrappers and allocator overhead remain covered by the measured Slurm
resource request. Large replay validation is separate.

[`particles.json`](particles.json) preserves exact source/script hashes, the
complete generated fixture, compiler/run commands, results, and original receipt
paths/hashes. Loose temporary builds were removed automatically; the original
small receipts were removed only after verifying their embedded readback.
