# halfstep_ab — ic2pm half-step velocity shift: unit checks and L1024 A/B test

Task spec: `2LPTIC_Gui/HANDOFF-ic2pm-halfstep.md` (2026-09-16). Branch `ic2pm-halfstep` off `cz` (13b1afd).
Change: `ic2pm.f90` only — new 3rd CLI argument `[half|sync]` (default `half`) that rescales the
2LPTic (synchronous) velocities to a_v = a_init - ASTEP/2 with PMP2start's growing-mode factor
`(a_v/a)^1.5 * F(a_v)/F(a)`, F = sqrt(Om + OmL a^3), ASTEP from the PM header (= ASTEP0 in ../Setup.dat).

## Binaries (all built 2026-09-16 with intel_comp/2024.2.0 / ifx, `env.sh`)
- `bin/ic2pm_legacy.exe`  — `make ic2pm` at cz HEAD 13b1afd, before editing (legacy, unshifted).
- `../ic2pm.exe`, `../PMP2main.exe`, `../PMP2init.exe` — working tree = 13b1afd + uncommitted `ic2pm.f90` edit.
  md5 in `bin/md5_final.txt`.
- Note (pre-existing, out of scope): the PM header field `AEXP0` is an implicitly typed, never-initialised
  local of `Tools::WriteDataPM`, so it is garbage in every PM header (legacy and new, also PMP2start);
  `PMP2main` never reads it (only commented-out lines 440–441).

## 5.1 unit / regression checks (`unit128/`, 128^3 test IC, Box 256, da = 8e-4) — PASS
- Config: `conv_da8/{PkTable,lcdm,TableSeeds}.dat` + `Init.dat` with Box 256 / Nrow 128 / Ngrid 256; `PMP2init.exe`
  output identical to `conv_da8/pmp2init.log` apart from Box/Nrow/Ngrid.
- U1 legacy vs new `sync`: `PMcrs0.DAT` byte-identical (md5 95c99c05…); `PMcrd.DAT` differs only in the 4
  garbage bytes of `AEXP0` (offsets 53–56). All positions and velocities identical (numpy check).
- U2 new `half` vs new `sync`: positions bit-identical; every velocity component ratio = 0.97014997
  (printed factor) to max |ratio/factor − 1| = 1.3e-7 (float32 rounding). Printed diagnostics:
  a_v = 0.0196, rms v_pec 50.331 km/s at a_init → 49.825 km/s at a_v.
- U3: bad epoch string → usage message, exit 1; ASTEP ≥ a_init → FATAL, exit 1
  (`stop 1` also used in `die`). [Corrected 2026-09-17, review 1 F3: at ee44100 ONLY `die` and a bad epoch
  argument exited 1; no-argument usage, missing ../Setup.dat or IC file, NROW mismatch and the particle-count
  checks still exited 0 (`stop 'text'` returns 0 under ifx). Fixed in the post-review commit below.]
- PMcrs0.DAT files (50 MB each) deleted after the checks; logs and PMcrd.DAT kept.

## 5.2 L1024 A/B (IC `/cosma7/data/dp004/bl267/Runs/DEGRACE/ICs/IC_data/L1024/Node_002/ics.`, S_vel 1.0)
Configs `da{4,8,16}_{sync,half}/` copied from `conv_da{N}/`; `Setup.dat` regenerated with the new
`PMP2init.exe` (identical to `conv_da*/Setup.dat` except the exact-z comment on line 21). Account dp004.

### Job log
| job | what | shape | result |
|---|---|---|---|
| 12006223 | pilot `da8_half`, 300 steps | cosma8-serial, 64 cores, 160G, OMP_PROC_BIND=true/OMP_PLACES=cores | TIMEOUT at 40 min after 33 steps: ~66 s/step vs ~1.5 s/step in July jobs 11555300–02 (cosma8 exclusive, 128 threads, no binding). UserCPU 39 h, SystemCPU 12 min (so not cgroup reclaim). MaxRSS 98 GB. PM files removed by hand. |
| 12006379 | diag 5 steps `da8_half` | cosma8-serial, 64 cores, 160G, NO binding env | COMPLETED 4:59 wall; still ~35 s/step (density 13 s, POTENT 15 s per step); UserCPU 3h29, SystemCPU 6 min; MaxRSS 58 GB (PMP2main). Binding was only a 2x factor. |
| 12006409 | diag 5 steps `da8_half` | cosma8-serial, 128 cores (mem 1000G default), no binding | COMPLETED 3:13 wall (ic2pm ~55 s, PMP2main 2:10); ~19.5 s/step (density 7 s, POTENT 7.5 s, move/init 5 s); UserCPU 4h28, Sys 9 min; MaxRSS 98 GB (ic2pm) / 59 GB (PMP2main). |

(continued below as results arrive)

Compute-free checks while diagnosing the slowdown (2026-09-16):
- Source diff cz@2026-07-02 (fa0fca2) .. HEAD: only `PMP2main.f90` (+212/-21, exact-z mode) and
  `PMP2mod_tools.f90` (+20) changed; density / FFT / power modules and the FFLAGS are identical.
- The IC files `ics.*` date from 2026-03-11 (unchanged since the July runs).
- July jobs 11555300-02: cosma8 exclusive, 128 threads, no binding env, ~1.45 s/step; CPU efficiency
  TotalCPU/(Elapsed*128) = 3h57/(2:58*128) = 63%.
- Nodes report Sockets=8 CoresPerSocket=16 ThreadsPerCore=2, CPUTot=128 (Slurm cpu unit = 1 core).

### Resolution of the 'slowdown' (2026-09-16 21:00)
The handoff's reference jobs 11555300-02 (3.0-3.6 min) were NOT the convergence runs: their PowerDM
files stop at step 1. The `conv_da*` outputs were written by the resubmitted jobs 11555311-13
(cosma8 exclusive, 128 threads): conv_da8 1:53:47 for 123 steps (~55 s/step incl. 33 analyses),
conv_da16 1:04:06 for 62 steps, conv_da4 TIMEOUT at 3:00 (PowerDM up to step 122 at 07:47; later files
09:16-09:20 — origin not traced). So the current binary is ~2-3x FASTER than July, not slower, and the
handoff's '~1 node-hour' estimate was ~30x too low.

Measured per-step cost (current binary, cosma8-serial, no binding):
| threads | s/step | core-s/step |
|---|---|---|
| 64  (160G) | 35   | 2240 |
| 128 (1000G) | 19.5 | 2500 |

Cost of the planned six variants (244+123+62 = 429 steps per epoch, 858 total, + 43 analyses each):
~4.5 node-hours. Spent so far on the pilot + diagnostics: 42.7 + 5.3 + 6.9 = 55 core-h (0.43 node-h).
STOPPED here for user sign-off (handoff: check in before compute beyond ~1 node-hour).

### Production A/B submits (2026-09-16 21:10, user-authorised; cosma8-serial, 64 cores, 160G, no binding, account dp004)
| job | variants | steps | time limit | dependency |
|---|---|---|---|---|
| 12006891 | A: da8_sync, da8_half | 123 each (to z=0) | 04:00 | — |
| 12006892 | B: da4_sync, da4_half | 244 each | 07:00 | afterok:12006891 |
| 12006893 | B: da16_sync, da16_half | 62 each | 03:00 | afterok:12006891 |
Estimated cost at 2240 core-s/step + ~15 s per analysis output + ~90 s ic2pm per variant:
da8 pair ~180 core-h, da4 pair ~330 core-h, da16 pair ~105 core-h; total ~615 core-h (~4.8 node-h);
time-limit ceiling (4+7+3 h x 64) = 896 core-h.

### Results, da8 pair (job 12006891; da8_sync 1:27:13 wall, 33 outputs; k <= 0.05 h/Mpc, mode-weighted)
- Binary check: new `da8_sync` vs July `conv_da8` (jobs 11555312): P ratio 1.00000 at every output, all k
  (files differ only in trailing digits at high k).
- **T1** R = P_sync/P_half: 1.00223 at step 1 (pred 1.00217); 1.0235 at z=9.4 rising to 1.0239 at z=1.5
  (pred 1.0227-1.0232). Measured offset / predicted offset = 1.03 at all z. PASS (criterion: within 15%).
- **T2** G(z) = [P(z)/P(step 1)] / [D(a)/D(a1)]^2, 2 <= z <= 10:
  `half` 0.9995 (z=9.4) -> 0.9968 (z=2.05): within +-0.5% throughout. PASS.
  `sync` 1.0207 (z=9.4) -> 1.0183 (z=2.05); pred +2.05..2.10%. PASS (within 15% of predicted excess at z>4).
  Both epochs share the same slow -0.35% drift from z=9 to z=2 (and widening bin spread), so it is not an
  epoch effect. [Explained 2026-09-17, review 2 R2-2: it is a common NONLINEAR term, ∝ k^2 D^2 (fit
  NL = -α k^2 (D/D0)^2 with α ≈ 10 (Mpc/h)^2 for all three da; -0.25..-0.27% from z≈9.3 to z≈2 for every da),
  plus a small da^2 leapfrog term that is k-independent and matches the discrete linear model
  (G_half-1 at z=0: -0.013 / -0.051 / -0.196% for da 4/8/16). Absolute G at z≈2 is therefore not a clean
  linear test; the contamination-free quantity is G_sync/G_half (T2' in compare_halfstep.py).]
- Job 12006891 accounting: Elapsed 2:53:17, TotalCPU 155.9 core-h of 185 allocated (84% CPU efficiency),
  MaxRSS 59 GB (PMP2main; ic2pm peak ~98 GB); da8_half PMP2main 1:23:15. PM files removed by the script.
- B jobs started 01:48 (2026-09-17): 12006892 da4 pair on m8376, 12006893 da16 pair on m8411.

### Results, da16 pair (job 12006893; da16_sync 50:03 wall, 30 outputs)
- **T1** R = 1.00831 at step 1 (pred 1.00784); 1.0470 (z=9.3) -> 1.0478 (z=2.1) vs pred 1.0439 -> 1.0448;
  offset ratio 1.07. PASS.
- **T2** `half` 0.9988 (z=9.3) -> 0.9954 (z=2.06), mean 0.9973: within +-0.5%. PASS.
  `sync` mean 1.0362 over 2<=z<=10 vs pred +3.58..3.66%. PASS.
- The common-mode drift of G with decreasing z is mostly the nonlinear k^2 D^2 term above (same for all da);
  the da-dependent part (da16-da4 = -0.085% at z≈9.3, flat in k) is the da^2 leapfrog term, generated mostly by
  the EARLY steps, not "the late, larger steps" as first written. [Corrected 2026-09-17, review 2 R2-2.]
- Job 12006893 (da16 pair) accounting: COMPLETED, Elapsed 1:44:09, TotalCPU 90.3 core-h of 111 allocated (81% efficiency);
  (85% efficiency); da16_half PMP2main 50 min. PM files removed by the script.

### Results, da4 pair (job 12006892; da4_sync 2:35:48 wall, 39 outputs) — partial, da4_half pending
- **T2** `sync`: 1.0109 (z=9.3) -> 1.0103 (z=4.4) vs pred +1.10..1.12%. PASS.
- **T3** sync epoch, z=0 growth-corrected P/P(da4_sync) at k = 0.02/0.05/0.1/0.2/0.4/0.8/1.5/2.5 h/Mpc:
  da8_sync  1.0113 1.0112 1.0112 1.0123 1.0145 1.0160 1.0142 1.0111
  da16_sync 1.0330 1.0324 1.0327 1.0373 1.0464 1.0531 1.0469 1.0352
  = the July conv_da* result (+1.1-1.6% / +3-5%), as expected since the sync outputs are reproduced.

### Final results (all six variants complete, 2026-09-17 07:00; full tables in `results_final.txt`)
- Job 12006892 (da4 pair): COMPLETED, Elapsed 5:11:07, TotalCPU 280 core-h of 332 allocated (84%);
  da4_half PMP2main 2:32:45. PM files removed by the script.
- Total compute for the task, as BILLED (allocation, not TotalCPU; corrected 2026-09-17, review 2 R2-5):
  production 184.8 + 331.9 + 111.1 = 628 core-h at 81-84% CPU efficiency, plus pilot+diagnostics 55 =
  683 core-h (≈5.3 cosma8 node-h). (TotalCPU used: 156 + 280 + 90 = 526 core-h.)

| da | T1 R=P_sync/P_half (z 9.3 / z 2) | pred | offset ratio | T2 half mean G (min..max), 2<=z<=10 | T2 sync mean G | pred sync |
|---|---|---|---|---|---|---|
| 4e-4   | 1.0118 / 1.0120 | 1.0116 / 1.0118 | 1.017 | 0.9988 (0.9974..0.9997) | 1.0101 | 1.0110..1.0112 |
| 8e-4   | 1.0235 / 1.0239 | 1.0227 / 1.0232 | 1.03  | 0.9984 (0.9968..0.9995) | 1.0199 | 1.0205..1.0210 |
| 1.6e-3 | 1.0470 / 1.0478 | 1.0439 / 1.0448 | 1.07  | 0.9973 (0.9954..0.9988) | 1.0362 | 1.0358..1.0366 |

**T3** z=0, growth-corrected, P/P(da4) at k = 0.02/0.05/0.1/0.2/0.4/0.8/1.5/2.5 h/Mpc:
  half: da8  0.9996 0.9996 0.9996 0.9995 0.9993 0.9992 0.9991 0.9989
        da16 0.9980 0.9979 0.9979 0.9989 1.0009 1.0025 1.0016 0.9985
  sync: da8  1.0113 1.0112 1.0112 1.0123 1.0145 1.0160 1.0142 1.0111
        da16 1.0330 1.0324 1.0327 1.0373 1.0464 1.0531 1.0469 1.0352
=> [Corrected 2026-09-17, review 2 R2-1.] The linear-D^2 correction of the da16 run's last output (a=1.0102)
to a=1 is invalid at nonlinear k, so the rows above overstate da16's accuracy at k ≳ 0.3. Compared at matched
epochs (a≈0.847: steps da4 221 / da8 111 / da16 56; T3' in `compare_halfstep.py`, band means in %):
  k band:      ≤0.05  0.05-0.1  0.1-0.2  0.2-0.5  0.5-1  1-2.5
  half da8     -0.04   -0.04     -0.04    -0.07   -0.11  -0.19
  half da16    -0.18   -0.18     -0.22    -0.40   -0.67  -1.19
  sync da8     +1.13   +1.12     +1.17    +1.41   +1.58  +1.22
  sync da16    +3.32   +3.27     +3.40    +4.03   +4.35  +2.94
With the half-step shift the genuine time-stepping error is <=0.1% for da8 at k ≲ 1 (-0.2% at 1<k<2.5) and
<=0.3% for da16 only at k ≲ 0.2 (-0.5..-1.2% at 0.5<k<2.5); at low k the residuals (-0.04% / -0.18%) scale as
da^2 and equal the discrete linear-leapfrog model. The July +1.1-1.6% / +3-5% "step-count convergence"
differences were the unshifted-velocity offset at k ≲ 0.2, but at 0.5<k<2.5 about 0.5-1.2 points of the da16
difference is real stepping error. The sync da16 "+5.3%" upper end is +4.5% at matched epochs.

**Verdict:** T1 and T2 pass for every da (T2 sync for da4 on the 2<=z<=10 mean; per output it falls below the
±15% band at z ≲ 2.6 because of the nonlinear drift — review 2 R2-3; the contamination-free G_sync/G_half is
+1.12..1.14% / +2.13..2.16% / +3.84..3.92%, matching the discrete-leapfrog model to ≤0.02%). T3 passes at low
k as stated above. `half` is the correct default. Predictions: the handoff's 0.8·eps under-predicts the offset
by 1.7 / 3.4 / 7.0%; the discrete model (REVIEW-1 App. A) gives +1.20 / +2.40 / +4.80% and reproduces the data.

### Sign-off (2026-09-17)
- Code committed: ee44100 (`ic2pm.f90`) on branch `ic2pm-halfstep`. Docs committed: c01255e (VALIDATION.md errata,
  README.md usage, AGENTS.md line). `DEGRACE-highres/INGESTION_SCOPE.md` (untracked dir) edited in place.
- Follow-ups NOT started: (1) re-run `fid2LPTIC_L512Np2048Ng4096{,_da4,_da6}` with `half` (~3 x 10 node-h) and
  re-derive the ZA-vs-2LPT P(k) offsets and HMF ratios; (2) re-check the da0=6e-4 recommendation from the high-k
  side; (3) tell the collaborator to use ic2pm at ee44100 with the default `half` and NOT to edit the Gadget files;
  (4) merge `ic2pm-halfstep` into `cz`. MG (f(R)) runs use the same leapfrog — no separate change.


### Post-review fixes (Reviews 1-3, 2026-09-17)
See `2LPTIC_Gui/REVIEW-{1,2,3}-*.md` and `2LPTIC_Gui/FIX-PROMPTS-ic2pm-halfstep.md`. Code: F1 (half-mode PM
HEADER 'ic2pm: 2LPTic ingest, v at a-da/2'; sync header unchanged; NOTE the half-mode A/B files of jobs 12006891-93
were written by an ee44100 build and carry the sync string - their mode is in ic2pm.log only), F2 (S_vel written as
5.168609e6; strict parse), F3 (every ic2pm check exits 1; runtime I/O errors exit with ifx codes), F4/F5 (header wording). Docs: R2-1..R2-5, R3-1..R3-12 (VALIDATION.md,
README.md, AGENTS.md, INGESTION_SCOPE.md, this file, compare_halfstep.py T2'/T3').
Tests on the scratch build (halfstep_ab/unit128/RunFix{Sync,Half,Default}/, 128^3 IC, OMP 4):
- T1 sync: PMcrs0.DAT md5 95c99c05be40b46648ed5a54eff2d01a (unchanged); PMcrd.DAT vs legacy differs only at
  bytes 53-56 (AEXP0).
- T2 half and no-third-argument: PMcrs0.DAT md5 c84f28f54e1df92b8262b14aeeb6242a (unchanged); PMcrd.DAT vs the
  pre-fix half output differs only at bytes 24-36 (HEADER: 'ic2pm: 2LPTic ingest, v at a-da/2') and 53-56.
- T3 exit status: 1 for no arguments, bad/empty/mixed-case epoch, S_vel given as `sync`, "", "5.12e6/0.99059529",
  "5.12e6,1", -1, abc, missing ../Setup.dat, missing IC file, NROW mismatch, Om mismatch, ASTEP0 >= a_init (half);
  0 for the success case (S_vel 5.168609e6 prints "S_vel = 5.1686E+06"). ASTEP0 >= a_init with `sync` is not
  checked (the shift is not applied), as specified.
- T4 `ifx -warn all -c ic2pm.f90`: no warnings.
Commit: 6b629f1 (2026-09-17) on `ic2pm-halfstep`. Root `ic2pm.exe` rebuilt from it: md5 4a1f071415ff95c6607532c3ad6e2f04
(pre-fix, used for jobs 12006891-93) -> fc5135980e8097df735b54f4c63031d7. Re-run with the root binary
(`unit128/RunRoot{Sync,Half}/`): PMcrs0 md5 95c99c05… (sync) / c84f28f5… (half), PMcrd diffs as in T1/T2 above.
`bin/md5_final.txt` records the pre-fix md5s of the A/B binaries.
Verification (Prompt V, `2LPTIC_Gui/VERIFY-6b629f1-findings.md`): accept with fixes; nits N1-N8 (+N10) applied in
commit 54a8872. Root `ic2pm.exe` rebuilt: md5 fc5135980e8097df735b54f4c63031d7 -> b757c920cee6db8fc967f433acb53aba; re-run (`unit128/RunRoot2{Sync,Half}/`):
PMcrs0 md5 95c99c05… / c84f28f5… unchanged. S_vel validator now rejects "1-2", "1e", "1e1e1", "1.0.0" (accepts 1.0d0,
1e-2, +2.5E+01).
