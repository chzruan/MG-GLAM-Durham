# Triage of the five reviews of `bdm-binary-catalogue` (2026-10-18)

Reviews run as independent subagents from `REVIEW-PROMPTS.md` against commit
167fb1c; findings in `REVIEW-<k>-*-findings.md`. Reported totals: **Blocker 2
(one issue, found twice), Major 9 (six distinct issues), Minor 28, Nit 28.**
Recommendations: R1 and R4 "not ready", R2 and R5 "ready after the listed
fixes", R3 "ready". Every Blocker and Major was re-checked by the author and
fixed in the follow-up commit on the branch (this file is part of it).

| review | Blocker | Major | Minor | Nit | scope |
|---|---|---|---|---|---|
| R1 Fortran writer, parser, publication | 1 | 2 | 3 | 5 | PMP2linker.f90, PMP2mod_extras.f90 |
| R2 format and sidecar semantics | 0 | 2 | 7 | 3 | definitions, units, epochs, prose |
| R3 glamio package | 0 | 0 | 8 | 10 | python/ |
| R4 tests, executed in a worktree | 1 | 1 | 5 | 5 | test/precise_z (33-mutation battery) |
| R5 downstream impact and hygiene | 0 | 4 | 5 | 5 | consumers, defaults, docs |

## Blocker, verified and fixed

| ID | Verified | What | Fix |
|---|---|---|---|
| R1-B1 = R4-F1 | yes, by R4's probe run and by the format arithmetic (`10x` + 15 characters reaches column 25; `T20` rewinds) | The default `BDM.config` the finder writes when none exists spelled the new key as `Catalogue = both`; the next `ReadParameters` (second catalogue output, or the first after a restart) rejected it as an unrecognized parameter. New in this commit: on `cz` the longest key written through these formats was 9 characters. | Format 30 is now `(10x,a,' = ',a,T45,a)`; new suite case14_catdefault runs with no `BDM.config` and two catalogue outputs and checks the written file spells the full key and is re-read. |

## Majors, verified and fixed

| ID | Verified | What | Fix |
|---|---|---|---|
| R1-M1 = R2-F1 | yes (`GetHalo`: `aperture = rso + grid*min(Rext/(rso/grid)^SlopeR, 0.75)`, `Rvir = aperture`, `Mvir = n*mass` over the bound survivors of `radii <= rso`; R2 measured a median enclosed overdensity of 99 inside the reported Rvir in the pilot) | The sidecar `centre =` line, the deck and the plan called Rvir the overdensity radius and Mbound the bound mass within it. | Sidecar line replaced by `definitions =` stating R_SO, the resolution correction, Mtot within Rvir, Mbound within R_SO, the bound-population quantities, the Vcirc/Cvir = 0 sentinel; `Rext`, `SlopeR`, `cell`, `a_x`, `MG_flag`, `MG_model` keys added; deck, plan, python README and `COLUMN_UNITS` corrected. |
| R1-M2 = R5-01 (+ R5-07) | yes (`case default` in the parser; `BDM.config` was parsed only inside `BDM`) | `gen_inputs.py` wrote the key unconditionally, and a pre-merge executable (the frozen pilot binaries) would reject it only at the first catalogue step, half-way through a run; any typo in `BDM.config` had the same late failure. | `CheckBdmConfig` (parse + validate, no side effects) runs at start-up when `Find BDM halos = 1` (`PMP2main.f90` after `ReadSetup`); the parser block was factored into `ReadBdmConfig`. `gen_inputs.py` writes the key only when `--catalogue-format` is given (help text names the minimum build). case4g now asserts the rejection happens before any time step. |
| R5-02 | yes | `check_run.py` knew only the `.DAT`: a `binary` run failed every catalogue check and a `both` run had its `.bin`/`.meta` unverified. | The checker accepts text and/or binary; for a binary it checks `.bin` size = 96 × rows, meta `a` against the schedule, `text_catalogue` against the presence of the `.DAT`, text row count = binary rows, and reports hidden `.Catshort*` staged temporaries. Verified on the ΛCDM pilot: 0 FAIL, 0 WARN. |
| R5-03 = R2-F3 | yes (`stat`: the z = 0 text catalogue is 220 MB apparent, 114 MB allocated; a float32 lens plane 201 → 185 MB) | The storage figures mixed apparent and on-disk bytes: this Lustre compresses the text 1.9× and float32 only 1.09×, so the binary saves ≈ 40 % of the catalogue space on disk, not 3×. | Plan sections 9, 11 (item 12) and 16, the deck and the companion README now give apparent and on-disk sizes (≈ 1.3 GB per 2048³ box, ≈ 0.25 TB over the suite). |
| R5-04 = R2-F2 | yes (`Vscale = 100 Xscale/AEXPN` vs `/a_v` in the extras: dividing the same momentum by the larger `a` gives the smaller number) | The README, plan section 3 and two deck lines said the catalogue velocities are "high by a/a_v, apply a_v/a"; the sidecar and glamio (× a/a_v) were right. The quoted 1.1 %/0.5 % were not the catalogue's numbers. | Prose corrected everywhere: low by a_v/a, 0.8 % at z = 2 (first catalogue output), 0.3 % at z = 0; Vrms, the virial ratio and λ carry the same offset (sidecar `velocity_scale` line extended). |
| R4-F2 | yes | No test exercised the default-config path, so 194/194 said nothing about the configuration every run without a `BDM.config` uses. | case14_catdefault (above). |

## Decision taken by the author (flag for the user)

**Compiled default changed from `both` to `text`** (R5-09, supported by
R5-03): every existing campaign directory lacks the key, the filesystem is
near 80 %, and `both` would have added ≈ 60 % to the on-disk catalogue space
of every BDM run silently. The binary catalogue is now opt-in through
`BDM.config`; `gen_inputs.py --catalogue-format both` is the recommended
setting for the companions until halocat reads the binary, then `binary`.
This reverses the "default both" of the 2026-10-17 proposal; it is a one-word
change in `PMP2linker.f90` (`CatalogueFormat='text'`) if the user prefers
`both`.

## Minors and Nits: handled

- R3-01 row self-check in the binary reader; R3-05 `columns`/`dtypes` lines
  required; R3-02/R2-F6 text and binary meta keys unified (`da`, `rows`,
  `realization`); R3-08 `mmap=True`; R3-07 version in one place (pyproject
  `dynamic`); R3-09/R2-F7 units and sentinels in `COLUMN_UNITS`, the sidecar
  and the README; R3-10/11/12/13/14/15/18 (repr, unreachable branch,
  assert → RuntimeError, case-insensitive `.dat`, text fallback when only
  `.DAT` exists, non-finite a/a_v rejected, exports and metadata); R3-06 two
  new unit tests (byte-swapped file, missing layout line).
- R4-F3 glamio exceptions caught as test failures; R4-F4 layout lines
  asserted; R4-F5 row totals must be positive; R4-F6 case15_cattext (text
  mode equals case13's binary rows column by column); R4-F7 hidden
  temporaries and snapshots in the forbidden globs; R4-F9/R5-14 trailing
  newline added to `test/BDM.config`; R4-F10 every column checked in case13;
  R4-F11 loadtxt warning suppressed; R4-F12 `catalogue_bin` of the fR_halo
  sidecar checked.
- R1-m1/R5-06 stale staged temporaries: documented in the plan and the
  companion README, reported by `check_run.py` (no automatic deletion: a
  concurrent finder could own them). R1-m2 `text_catalogue` taken from the
  staged path. R2-F4/F8/F10 deck numbers (effective mass floor, offsets,
  Δ_vir for Ω_m = 0.3). R2-F9 `a_x`, `MG_flag`, `MG_model` added (not `fR0`:
  not visible to the finder module). R5-11/12 plan figures and sources.

## Minors and Nits: left as they are

- R3-03/R3-17 the text reader's float32 cast (documented; `np.loadtxt` for
  the printed decimals). R3-16 byte-order-naive dtype comparison (COSMA is
  little-endian). R2-F5/R5-08 the fR_halo sidecar names the `.DAT`/`.bin`
  unconditionally and hard-codes the `V` label (pre-existing pattern,
  labelled "when present"). R2-F11 truncated text caption `RadRMS/k`
  (pre-existing). R2-F12 nine printed digits for real*4 epochs (same as the
  other sidecars). R5-10 `gen_inputs.py` overwrites `Run1/BDM.config` as it
  overwrites `Init.dat`. R5-13 the prompts file says one commit (the second
  commit was the prompts file itself).

## Verification after the fixes

- `test/precise_z/run_tests.sh`: 215/215 checks (was 194; cases 14 and 15
  and the hardened checks added), 46 s, 4 threads, 68 MB.
- `python/tests`: 7 passed (pytest and script mode).
- `test/extras/run_modes_test.sh`: V1 and fR_halo PASS (unchanged code).
- `DEGRACE-companion/check_run.py` on the ΛCDM pilot: 0 FAIL, 0 WARN, 17 of
  17 catalogues, no staged temporaries.
- `gen_inputs.py` without `--catalogue-format` writes no key; with `both`
  writes `CatalogueFormat = both`; provenance records the choice.
- Deck rebuilt (19 frames); the two changed frames inspected.
