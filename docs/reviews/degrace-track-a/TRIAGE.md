# Triage of the six reviews of `degrace-track-a` (2026-10-17)

Reviews run as independent subagents from `REVIEW-PROMPTS.md`; findings in
`REVIEW-<k>-*-findings.md`. Totals: **0 Blocker, 8 Major, 47 Minor, 38 Nit**.
Every reviewer recommends "ready after the listed fixes". All eight Majors were
re-checked against the code by the author (this file); none is a wrong result
in the pilot runs that exist, four are latent misuse traps, two are test gaps,
and three are documentation that does not match the implementation.

## Majors, verified

| ID | Verified | What | Fix |
|---|---|---|---|
| R3-M1 | yes: `ic2pm.f90:244` reads pos, vel and skips the ID block; `PMP2mod_subsample.f90:184` prints the PMP2start x-fastest lattice formula | The subsample hashes the array index; selection is reproducible across twins, threads and restarts, but Lagrangian positions cannot be recovered as the `.meta` implies: Baojiu's 2LPTic files are z-fastest and multi-rank FML output is in neither order | Decision needed: (a) make `ic2pm` read the Gadget IDs and scatter particles so index = ID (8 GB of IDs at 1024³, one permutation pass), then state the generator's lattice formula in the `.meta`; or (b) keep index order and write `lattice = IC file order (unknown)` plus drop the Lagrangian claim from the plan |
| R2-M1 | yes: `PMP2init.f90:762-810` returns the default when the fixed-width read fails (`ios /= 0`) | Optional `Init.dat` lines that are present but unparseable (comment inside the 10/12-character window, `=` in a comment) silently take the default; the legacy parser aborts. The pilots are unaffected (checked: the f(R) pilot's `Setup.dat` has every value) | Default only on EOF or no `=`; `Stop` with the line text on a read error |
| R2-M2 | yes: `PMP2init.f90:282`, `PMP2linker.f90:209` gate on `MG_flag` only | DGP/symmetron/k-mouflage/csf with `fR at halos = 1` would write an `fR_halo` product from the wrong scalar | Add `MG_model == 1` to both gates; write `MG_model` into the sidecars (R1-14) |
| R5-M1 | yes, observed: two epochs skipped in the ΛCDM comparison | `compare_ecosmog_pk.py` matches the `f8.5`-rounded `Aexpn` with a 2e-5 relative tolerance; below a = 0.25 the rounding exceeds it (259 of 4032 epochs suite-wide) | Take the epoch from `OutputSchedule.dat` by step |
| R5-M2 | yes, by construction | No suite test is sensitive to axis order or the sign of Im in the production `deltam`/`mom_*` dumps (V2, V6 compare power; V1 exercises `FFT3Dforward` on a synthetic field). `check_glam_deltaL.py` covers `deltam` end to end but is outside the main suite and does not cover `mom_*` | In case 11 (M = 1) compare the particle cubes already built for V6 with the dumps mode by mode |
| R6-M1 | yes | Plan section 7 "Runtime flags" table: wrong defaults (`FR_AT_HALOS` on vs off, `SUBSAMPLE_M` 100 vs 0, `N_DEP` N_g vs N_g/2) and three flags that do not exist | Replace the table with the implemented `Init.dat` keys and defaults |
| R6-M2 | yes (`PMP2linker.f90:537`, `Vscale = 100 Xscale/AEXPN`) | Plan section 3 says every velocity product uses the same (x, v) pair and records a_v; the BDM catalogue converts with a, not a_v (1.1 % at the first pilot output, 0.5 % at z = 0) | Decision needed: fix the legacy catalogue conversion (changes legacy output) or document it as a known offset in plan section 3, the README and the catalogue reader |
| R6-M3 | yes (`interlaced = false` always) | Plan I8 requires interlacing when N_dep < N_g and k_c > k_N,dep/4; the companions (N_dep = 1024, k_c,mom = 1.0 > 0.785) are in that regime and nothing interlaces | Decision needed: implement interlacing for the momentum deposit, or lower `kc_mom` to 0.75 at N_dep = N_g/2, or deposit momentum on N_g (memory); record the choice |

## Minors worth fixing with the Majors (grouped)

- Robustness of the parameter path: `ReadSetup` treats EOF and parse errors
  alike and resets six flags silently (R2-m3); no validation of the block it
  reads (R2-m2); `g12.5` input scales `kc_delta = 1` to 1e-5 (R2-m1);
  `Nout > 1000` overflows `zout` (R2-m5, pre-existing); `PMP2analyze.exe`
  inherits the hooks and overwrites a run's products (R2-m6).
- f_R product: hard-coded `CatshortV` label in the meta (R1-05); no
  `fR_halo` file for an empty catalogue (R1-07); written unstaged (R1-08).
- Memory and cost: separate N_dep³ buffer (34 GB at N_g = 4096) instead of
  FI reuse (R1-09, R6-m4); 8 atomics per particle in `MomentumDeposit`,
  untimed at 2048³ (R1-01; the f(R) pilot will time it); every FML rank holds
  the full δ_L cube (R4-14).
- δ_L: D_gen from the Lua cosmology only, no cross-check against the P(k)
  table or the GLAM run (R4-1); rank-0 open failure deadlocks (R4-5); the
  production Lua rerun would rewrite the 2048³ IC (R4-7); `run_cosma.sh`
  uses 192 ranks but FFTWGrid needs Nmesh divisible by the rank count (R4
  open question); V3 tolerances looser than the plan's criterion (R4-2).
- Tests: `test/PkTable.dat` and `test/TableSeeds.dat` are untracked but
  required (R4-3, R6-m3); `Main.o` has no header dependency (R4-4); case 12
  passes vacuously with no haloes (R5-m2); velocity epoch not independently
  asserted for the subsample (R3-m1); L3's N_eff ignores the two slabs per
  axis (R5-m1); `check_run.py` hard-codes 17 + 25 classes (R5-m3, R6-m6);
  `gen_inputs.py` field widths for Nrow ≥ 10000 / checkpoints ≥ 1000 (R5-m4);
  provenance JSON lacks generator, LPT order, seed and the ic2pm epoch
  (R5-m5, R6-m5); `__pycache__/` and `test/precise_z/work/` not ignored.
- Docs: stale status lines, L2 claim vs the specified test, open questions
  Q8/Q9 already decided, README product lists and build steps (R6-m1, m2,
  m7, m10, m11, n7); `rms_deltaL_cube` undefined (R4-9).

## Verified-OK highlights (do not re-review)

FFTPACK5 packing and the sign of Im derived independently (R1); `FFT3Dforward`
is a faithful copy of `POWERfft5`; lens planes orientation, normalisation,
mass conservation, thread ownership; momentum deposit mapping and a_v epoch
consistent with `MOVE`; f_R mapping matches `POTENT_FOR_MG`; growth and
comoving-distance quadratures to 1e-10; no unit-number collisions (R1).
Classes follow the exact-z permutation; `Nclass` is rebuilt deterministically
on restart; legacy path identical to `cz` apart from `OutputSchedule.dat`, one
console line and the appended `Setup.dat` block (R2). splitmix64 constants
against published vectors; 64-bit counters (2048³ safe); two-pass compaction
deterministic; exact quantisation/decoding; GLAM never permutes particle
arrays (R3). δ_L dump index- and convention-exact against the FML sources;
Simpson growth to 2e-11; `llround` provably correct (R4). Pilot numbers,
commit hashes and `.meta` key sets as documented (R6).
