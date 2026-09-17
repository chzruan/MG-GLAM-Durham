# Validation of the Intel-built FML LPT generator vs existing DEGRACE ICs

Date: 2026-07-28. Binary: `FML/FML/LPT/example/test` (mpiicpx,
intel_comp/2024.2.0, see README.md). All jobs account dp203.

## Reference IC sets

- **Low-res**: `/cosma7/data/dp004/bl267/Runs/DEGRACE/ICs/IC_data/L1024/Node_002/ics.{0..27}`
  1024^3, L=1024 Mpc/h, z=49, seed 2026 — generated with **original 2LPTic**
  (param: Nmesh=Nsample=1024, sigma8=0.8000650 renormalized at z=0, spectrum
  shape `PK_2LPTic_z00_002.dat`, growth D(0)/D(49)=38.9492 from
  `inputspec_ics.txt`; glass1_le 1-particle glass tiled 1024^3 = regular grid).
- **High-res**: the `gadget_z49.000.*` files in `.../L1024h/Node_002/` NO LONGER
  EXIST (only `white_noise_real.*` dumps remain). Used instead the complete
  2048^3 FML set `/cosma8/data/dp203/bl267/Projects/Ongoing/HEFT/ICs/IC_highres/IC_Np1d_2048_L_1024_2LPT.{0..255}`.

## Input spectra

- `validate_lowres/pofk_lowres_exact_z49.txt`: exact effective P(k,z=49) of the
  low-res run = `PK_2LPTic_z00_002.dat` shape (log10 k, log10 Delta^2)
  renormalized to sigma8=0.8000650, /38.9492^2. Cross-checked against 2LPTic's
  own `inputspec_ics.txt` dump: agreement 2e-5 median, 3e-6 scatter.
- `pofk_bli_z49.txt` (Gui's, used for the HEFT 2048^3 set) is a **flat x1.0181
  in P** above the low-res normalization — exactly the ~1.8% offset seen in the
  2026-07-01 low-vs-high phase test. Used for our 2048^3 run to match HEFT.

## Low-res validation (1024^3): PASSED

Generation: job 11658937 (cosma8-serial, 32 ranks x 2 threads, 500G, 99 s).
Internal check: P_ini/P_input = 0.9987 (0.1<k<1), max 1LPT displacement
0.816 cells (2LPTic log: 0.780).

- **P(k) cross-validation** (job 11658967, nbodykit Nmesh=1024, CIC+interlacing,
  `pk_validate_lowres.npz`): over all 511 bins to k_Ny = 3.14 h/Mpc:
  - r(k) = P_x/sqrt(P_ours P_ref) >= 0.99999816 (1-r < 1.9e-6) — **mode-by-mode
    phase identity** between our FML IC and the original-2LPTic IC.
  - P_ours/P_ref in [0.9997, 1.0003], median 0.99998 — **amplitude exact**.
- **Per-particle** (`compare_particles.py`, our file .0 vs ics.0, 29.4M matched
  particles via IDs):
  - positions: rms diff 2.2e-5 Mpc/h (= float32 storage precision), max 1.2e-4.
  - velocities: rms diff 0.088 km/s = 2.4e-4 of u_rms; u_rms ours 361.7 vs
    ref 362.0 km/s — **confirms our Gadget velocity-normalization fix**.
  - 2LPTic ID convention (for ic2pm etc.): **z-fastest, id = 1 + iz + N(iy + N ix)**;
    both codes place the unperturbed lattice at cell corners (offset 0).

## High-res validation (2048^3): PASSED

Generation: job 11658977 (cosma8-shm mad05, 64 ranks, 3000G, 8 min;
first attempt 11658973 failed — FML requires Nmesh % NTasks == 0).
Output: 64 files, exactly 2048^3 particles (275 GB), internal
P_ini/P_input = 0.9987, max 1LPT displacement 1.67 cells.

- **P(k) cross-validation** (job 11659092, nbodykit Nmesh=2048,
  `pk_validate_highres.npz`): all 1023 bins to k_Ny = 6.28 h/Mpc:
  - r(k) >= 0.999999847 (1-r < 1.6e-7)
  - P_ours/P_ref = 1.0000000 (within 5e-7), median exactly 1.
- **Per-particle spot check** (`compare_particles_2048.py`, our file .0 vs
  HEFT file .0, 25.2M particles matched by KD-tree, 100.000% match rate):
  - positions **bit-identical** (rms = max = 0 exactly), across different
    compilers (icpx vs g++), MPI libs, and rank counts (64 vs 256).
  - velocities: after fixing Gui's 5.12e6 scale bug, ours/HEFT = 1.00951
    = 1/0.99059529 exactly — the HEFT run used growth_rate1=0.99059529
    instead of f(z=49)=1.0000, i.e. **HEFT velocities are an additional
    ~0.95% low** (our f=1.0 velocities match original 2LPTic to 2.4e-4,
    so ours are the correct ones). Residual rms after this factor: consistent
    with pure scaling.

## Verdict

The Intel-compiled FML LPT generator is **fully validated**: it reproduces
the original-2LPTic low-res IC mode-by-mode and particle-by-particle (float32
precision), and Gui's high-res FML IC bit-for-bit in positions, while fixing
three defects of the existing HEFT 2048^3 set (IDs all 1; velocities x5.12e6
too small; velocities additionally x0.9906 low). Our regenerated
`validate_highres/snap/IC_ours_Np2048_L1024.{0..63}` (seed 2026, HEFT
normalization; positions identical to the HEFT set) is a drop-in replacement
with correct velocities and IDs.

## End-to-end GLAM run vs the fiducial suite (2026-07-29)

Full-chain test: fid-cosmology 2LPTIC IC (L=512, 2048^3, seed 2026, 2LPT
@ z=49; job 11659350) -> ic2pm (S_vel=1.0) -> PMP2main 123 steps to z=0
(job 11659352, cosma8/dp004, 11h40m; config `fid2LPTIC_L512Np2048Ng4096`,
identical to `fid_LCDM_L512Np2048Ng4096` except z_init=49, da=8e-4 = same
da/a=0.04). Compared with fid Run1-5 (GLAM ZA ICs @ z=100, seeds 1-5, 157
steps): `compare_fid.py`, `pk_fid_compare.png`.

Result: P_2LPTIC/<P_fid> = +4-6% at 0.3<k<1 (z=0), k-dependent, present
unchanged from the first common output (z=2.5) onward. Diagnosis (all
measured, not assumed):
- our chain is self-consistent: the Gadget IC realizes its input to 0.04%,
  and GLAM's z=2.5 spectrum matches the linearly-evolved measured IC
  spectrum + expected nonlinear growth;
- the fid boxes' own step-0 P(k) sits 1.1% BELOW their input table, and
  between IC and z=2.5 they grow by only 1.003x linear at k~0.07 (vs our
  1.023x = linear + nonlinear), losing progressively more toward high k
  (at z=2.5 the deficit reaches ~20% at k~10, washing down to 2-4% by z=0);
- i.e. the offset is the well-known ZA-transient under-growth of the
  z=100 Zel'dovich starts (GLAM's "Pk tune = 1.005" boosts the PMP2start IC
  amplitude by 1.005, i.e. +1% in P, still undercompensating), not an error
  in the 2LPTIC/ic2pm chain.
  [ERRATUM 2026-09-17: the 1.023x was NOT "linear + nonlinear" — it is the
  ic2pm half-step velocity offset (+2.4% at linear k for da=8e-4). The z=2.5
  "match" in the first bullet therefore hid that excess, and ≈ +3% of the
  +4-6% at 0.3<k<1 IS an error in the ic2pm chain (estimate). See the
  2026-09-17 erratum section at the end of this file. [corrected 2026-09-17 after review]]

Note on "Pk tune" (traced 2026-07-29): it is BiasPars(10), used in exactly
one place — PMP2init multiplies the box rms delta rho/rho written to
Setup.dat, which PMP2start uses as the overall IC normalization AMPLT
(PMP2init.f90:265, PMP2mod_tools.f90:211, PMP2start.f90:345). PMP2main
never consumes it (AMPLT is carried as header metadata only), so it is
INERT for runs whose ICs come from ic2pm. For hygiene, 2LPTIC-based
configs should set Pk tune = 1.0 in Init.dat; legacy ZA configs must keep
1.005 to reproduce the existing suites.

Consistency verdict: the two pipelines agree at low k within realization
scatter; the systematic difference at quasi-linear/nonlinear k is
attributable to the IC method, with the 2LPT @ z=49 start being the more
accurate. Any suite mixing legacy ZA-started GLAM boxes with 2LPTIC-started
ones must expect (and correct for) this few-% k-dependent offset.
[ERRATUM 2026-09-17: ≈ +3% of this offset at 0.3<k<1 is the ic2pm half-step
velocity offset, not the IC method; the ZA-transient part there is ≲1%
(estimate). The "few-% k-dependent offset" is withdrawn pending re-measurement
with `half` ICs; see the erratum section, item 1.]

Pk-tune accounting in the comparison (added 2026-07-29): the fid boxes'
ICs were deliberately boosted by Pk tune = 1.005 (x1.010 in P). Their
measured step-0 amplitude of 0.989 x intended table ALREADY includes that
boost (the raw GLAM draw+normalization alone would sit ~0.979). Decomposing
the as-run z=0 ratio at 0.3<k<1: ours/fid = 1.048 = [IC amplitude offset
1/0.989 = 1.011] x [pure ZA-transient growth deficit ~1.037] (linear
response; the amplitude part propagates slightly superlinearly at
nonlinear k, so ~1.037 is an upper bound on the pure-growth part there).
Without the tune, the as-run offset would have been ~1% larger (~+6%).
The dashed reference line in pk_fid_compare.png marks the 1.011 level.
[ERRATUM 2026-09-17: ≈1.033 of the "~1.037 pure ZA-transient" factor at
0.3<k<1 is our own half-step offset (1.4× its linear +2.4%; estimate); see
the erratum section. [corrected 2026-09-17 after review]]

## Time-stepping convergence at Ng=4096 (2026-07-30)

Motivated by GLAMdoc.pdf Sec. 5 (eq. 30-32: allowed step scales with force
resolution, beta = v dt / dx < 1; recommended late-time da/a = (0.75-1)e-2).
Schedules (recomputed exactly with PMP2init): fid 157 steps (z=100,
da=4e-4, late da/a=1.0e-2); our first run 123 steps (z=49, da=8e-4, late
da/a=1.4e-2 — ABOVE the recommendation at dx=0.125 Mpc/h); convergence run
244 steps (z=49, da=4e-4, late da/a=0.68e-2, early band also 2x finer;
job 11661548, 10h15m, config `fid2LPTIC_L512Np2048Ng4096_da4`, Pk tune=1.0).

P(123 steps)/P(244 steps), same IC, z~0: +1.1-1.3% for 0.05<k<0.8 (matches
the 2026-07-02 Ng=2048 test), crossing zero at k~2, then -2.4% (k~4) and
-4.3% (k~10): the coarse late-time step suppresses halo-scale power at this
force resolution — the resolution-dependent regime of eq. 32 that the
Ng=2048 test could not see. Lesson: at dx~0.125 Mpc/h use da=4e-4-class
schedules (late da/a<1e-2); re-evaluate per config via beta.
[ERRATUM 2026-09-17: the +1.1-1.3% at k<0.8 is the difference of the
half-step offsets (1.0240/1.0120 at linear k, ≈1.016 at 0.4<k<0.8; estimate),
not a stepping error; the high-k suppression is genuine and ≈1 point larger
net of the offset. See the erratum section. [corrected 2026-09-17 after review]]

CONVERGED comparison vs fid Run1-5 (244-step run): +1.1% (0.1<k<0.2),
+2.5% (0.2-0.4), +4.1% (0.4-0.8), +4.6% (0.8-1.6), +4.7% (1.6-3.2),
+5.6% (3.2-6.4); 0.3<k<1 median +3.9% = [1.011 IC amplitude offset incl
tune] x [~1.028 pure ZA transient]. At k>3 the fid boxes' own late-time
stepping (da/a=1.0e-2, beta>1 in clusters at this dx) also suppresses
their P, so the high-k gap mixes IC transients with fid stepping error;
the 244-step 2LPTIC run is the most accurate box of the set.
[ERRATUM 2026-09-17: every number in this paragraph contains the 244-step
run's half-step offset (≈ +1.2% at k<0.2, ≈ +1.7% at 0.3<k<1, ≈ +1.0-1.6% at
k>0.8; estimates), so "~1.028 pure ZA transient" is ≈1.011 and "the most
accurate box of the set" holds only after that correction. See the erratum
section, item 2.]

## Default-schedule run (163 steps) — stepping ladder complete (2026-07-30)

Third run with the same IC: da0=6e-4 from z=49 reproduces the production
default exactly where it matters (final constant da=1.0252e-2 identical to
fid's; 66 steps at z<2 identical; 163 steps total; early band slightly
finer, max da/a=0.035; config `fid2LPTIC_L512Np2048Ng4096_da6`,
job 11663514, 8h22m).

P(N steps)/P(244 steps), same IC, z~0 (pure stepping error):
  k range        123 steps   163 default
  0.05-0.8       +1.1-1.3%   +0.6-1.1%
  0.8-3.2        +0.6..-0.6% +0.6-0.9%
  3.2-6.4        -2.4%       -0.1%
  6.4-12.6       -4.3%       -0.7%
  12.6-25        -5.2%       -1.2%

The default schedule is accurate to ~1% everywhere at this force
resolution (dx=0.125 Mpc/h): mild quasi-linear over-growth +0.6-1.1%
(residual of 244 itself estimated +0.5-0.8% by da^2 scaling, so absolute
errors are slightly smaller), and high-k suppression only -0.7..-1.2%
(vs -4..-5% for the 123-step run). The 123-step (da0=8e-4) schedule is NOT
recommended at dx~0.125.
[ERRATUM 2026-09-17: the "+0.6-1.1% quasi-linear over-growth" of the 163-step
run relative to 244 steps, and its +0.6-0.9% at 0.8<k<3.2, are the difference
of the half-step offsets (1.0180/1.0120 = 1.006 at linear k, ≈1.008 at
0.4<k<0.8; estimate), not stepping errors, so the table's "(pure stepping
error)" label holds only at k>3. The "residual +0.5-0.8% by da^2 scaling" of
the 244-step run was its own +1.2% offset. Net of the offset the high-k
suppressions are ≈0.5 point (163) and ≈1 point (123) larger; the da0=6e-4
recommendation at dx~0.125 stands. See the erratum section, item 3. [corrected 2026-09-17 after review]]

Stepping-matched IC comparison — P(2LPTIC, 163-step default)/<P_fid,
157-step default>: +1.7% (0.1<k<0.2), +3.5% (0.2-0.4), +5.2% (0.4-0.8),
+5.6% (0.8-1.6), +5.3% (1.6-3.2), +5.6% (3.2-6.4). This is the practical
"switch the IC method in production, keep default stepping" offset.
[ERRATUM 2026-09-17: contains the 163-step run's half-step offset, ≈ +1.8% at
linear k, ≈ +2.2% at 0.2<k<0.4 and ≈ +1.5% at 3.2<k<6.4 (estimates); see the
erratum section, item 3.]

## Halo mass function comparison (2026-07-30)

BDM Mtot HMFs at z=0 (all runs interpolated in ln n vs ln a to a=1 using
their last two outputs; `compare_hmf.py`, `hmf_fid_compare.png`).
Total z~0 halo counts: 2LPTIC runs 2.032-2.041M (nearly step-independent);
fid boxes 1.68-1.75M — the 2LPTIC boxes have ~18% more haloes overall.

n_2LPTIC/<n_fid> at z=0 (244-step converged run):
  log10M   11.45  12.05  12.55  13.05  13.55  14.05  14.55  14.95
  ratio    1.39   1.08   1.04   1.01   1.02   1.02   1.07   1.4(noisy)
  (fid single-box scatter: 0.5-1% at 12-14, ~10-17% in the tail)

Key findings:
- The low-mass excess (+8% at 1e12 rising to +39% at 1.4e11 ~ 100
  particles) is IDENTICAL across 123/163/244 steps -> an IC-method
  effect, not stepping: ZA@z=100 transients suppress the small-scale
  structure that BDM needs to detect marginally-resolved haloes.
  Matters for any HOD/emulator use near the resolution limit.
- Well-resolved intermediate masses (~1e13, ~7500 particles): agreement
  to ~1% — the chain is fully consistent where the HMF is robust.
- Cluster masses (13.5-14.5): +2% (converged) with mild step dependence
  (123-step: +4-5%), consistent with the P(k) excess propagated through
  the HMF response; extreme tail consistent within realization scatter.
  [ERRATUM 2026-09-17: both the +2% and its step dependence contain the
  half-step P offset; see the erratum section, item 4.]

## ERRATUM (2026-09-17): `ic2pm` half-step velocity offset in every run above

(Section corrected in place 2026-09-17 after the independent reviews
`REVIEW-1-ee44100-findings.md`, `REVIEW-2-halfstep-ab-findings.md`,
`REVIEW-3-c01255e-errata-findings.md`: k-dependence of the offset, the T3
claim, and the quasi-linear estimates.)

All GLAM runs in or cited by this file that started from `ic2pm`-converted
ICs (`fid2LPTIC_L512Np2048Ng4096{,_da4,_da6}`; `conv_da{4,8,16}` = "the
2026-07-02 Ng=2048 test"; `ic2pm_val_L1024`), and any other ΛCDM or MG run
converted with `ic2pm` before commit ee44100, were made with the unshifted
converter, which passed
the synchronous 2LPTic velocities (at a_init) straight into GLAM's PM files.
GLAM's kick-then-drift leapfrog (PMP2main::MOVE) expects momenta at
a_init - da/2, which PMP2start supplies (AEXPV = AEXPN - ASTEP/2). The
missing shift over-boosts every momentum by ~0.75 da/a_init at the first
kick; in linear theory 40% of that feeds the growing mode, so at linear
scales the late-time P(k) is high by 0.8 × 0.75 da/a_init = 0.6 da/a_init
(an exact linear model of PMP2main's discrete leapfrog agrees to 0.1%; the
continuum 0.8·eps with eps normalised at a_init + da/2 is 2-7% lower). The
offset is NOT k-independent once modes are nonlinear: in the A/B at z=0 it
is ≈1.2× the linear value at 0.2<k<0.4 and ≈1.4× at 0.3<k<1 (at z=2.6:
≈1.25× at 0.3<k<1, ≈1.7× at 1.6<k<3.2).

| da0 (z_init=49) | steps | predicted P excess (discrete leapfrog model) | measured (L1024/Ng2048 A/B, k ≤ 0.05 h/Mpc, z≈9 .. z≈2) |
|---|---|---|---|
| 4e-4 | 244 | +1.20% | +1.18% (z 9) .. +1.20% (z 2) |
| 6e-4 | 163 | +1.80% | (not run; ≈ +1.8% interpolated) |
| 8e-4 | 123 | +2.40% | +2.35% .. +2.39% |
| 1.6e-3 | 62 | +4.80% | +4.70% .. +4.78% |

At 0.3<k<1 and z=0 the same A/B gives +1.67% / ≈+2.5% (interpolated) /
+3.35% / +6.70%; these, not the linear values, are the offsets to remove from
the quasi-linear ratios quoted in this file.

Fix: `ic2pm.f90` commit ee44100 (branch `ic2pm-halfstep`): new third
argument `[half|sync]`, default `half` = rescale the velocities to
a_v = a_init - ASTEP/2 with PMP2start's growing-mode factor
(a_v/a)^1.5 F(a_v)/F(a), F = sqrt(Om + OmL a^3); `sync` reproduces the old
particle files byte-for-byte (PMcrd.DAT differs only in the uninitialised
AEXP0 field). A/B test (`halfstep_ab/README.md`, an untracked directory at
/cosma8/data/dp203/dc-ruan1/mgglam_claude/MG-GLAM/halfstep_ab/, jobs
12006891-93, same Node_002 IC as conv_da*, all six da x epoch variants,
k <= 0.05 h/Mpc):
- T2 (decisive): from the step-1 output to z=2, the `half` runs grow as
  linear theory to 0.9973-0.9988 (min 0.9954), the `sync` runs 1.0101 /
  1.0199 / 1.0362 for da 4e-4 / 8e-4 / 1.6e-3 averaged over 2≤z≤10
  (predicted 1.011 / 1.021 / 1.036; per output the da4 value falls to 1.0088
  by z=2.25 because both epochs share a ≈ -0.25% nonlinear drift at
  k ≤ 0.05). The contamination-free ratio G_sync/G_half is 1.011-1.014 /
  1.021-1.022 / 1.038-1.039, matching an exact linear model of the leapfrog
  to ≤0.02%.
- T3: z=0 growth-corrected P(da8)/P(da4) and P(da16)/P(da4) at k <= 2.5:
  `sync` +1.1-1.6% / +3.3-5.3% with a linear-D^2 growth correction (= the
  2026-07-02 conv_da* result; +3.3-4.5% for da0=1.6e-3 when compared at
  matched epochs); `half` 0.9989-0.9996 / 0.998-1.0025 with the same
  correction, which is not valid at nonlinear k for the da0=1.6e-3 run's
  last output (a=1.0102). Compared at matched epochs (a≈0.847), the genuine
  time-stepping error at dx=0.5 Mpc/h is -0.04% (da0=8e-4) and -0.18%
  (da0=1.6e-3) at k ≲ 0.2, scaling as da^2 as a linear model of the leapfrog
  predicts, but -0.1..-0.2% (da0=8e-4, 1<k<2.5) and -0.5..-1.2%
  (da0=1.6e-3, 0.5<k<2.5) at nonlinear k. The 1-5% previously reported was
  the offset at k ≲ 0.2 but not entirely at k ≳ 0.5.

Consequences for the statements above (earlier text left as written). The
numbers below are ESTIMATES, not re-runs: they divide the Ng=4096
(dx=0.125 Mpc/h) ratios by the offsets measured at the same k and z in the
L1024/Ng2048 (dx=0.5 Mpc/h) A/B, which assumes the nonlinear response
transfers between resolutions (least safe at k>1, unmeasured at k>6). None
of the fid2LPTIC_L512 runs has been repeated with `half`.
1. "End-to-end GLAM run vs the fiducial suite": the 123-step run's 1.023x
   growth relative to linear from the IC to z=2.5 at k~0.07 is the offset
   (predicted 1.023), not "linear + nonlinear". The fid ZA boxes' 1.003x is
   approximately correct. In the decomposition ours/fid = 1.048 =
   1.011 (IC amplitude) x 1.037 at 0.3<k<1 (z=0), the 123-step offset at
   those k is ≈ +3.3% (1.4× its linear +2.4%), so the pure ZA-transient part
   is only ≈ 1.048/(1.011 × 1.033) ≈ +0.3% (estimate; the 244-step pair in
   item 2 gives ≈ +1.1%), not ~3.7%. At 0.3<k<1 a ZA-transient deficit is
   therefore ≲1% and not established by this comparison: the "few-%
   k-dependent offset" between ZA- and 2LPTIC-started boxes is withdrawn
   until re-measured with `half` ICs. "2LPT@z=49 is the more accurate start"
   remains supported at k ≳ 1 (≈ +2 to +3.5% at z=0 after removing the
   offset and the 1.1% amplitude; estimate), by the ~20% deficit at k~10 at
   z=2.5, and by the low-mass HMF excess.
2. "Time-stepping convergence at Ng=4096": P(123)/P(244) = +1.1-1.3% at
   0.05<k<0.8 is essentially the offset ratio: 1.0240/1.0120 = 1.012 at
   linear k, rising to ≈ +1.4% at 0.2<k<0.4 and ≈ +1.6% at 0.4<k<0.8 at z=0
   (estimate). At k ≲ 0.2 it is not stepping; at 0.4<k<0.8 the measured
   value sits ≈0.3-0.5 point below the offset, i.e. the 123-step run's
   genuine stepping error is already slightly negative there (estimate).
   The -2.4% (k~4) / -4.3% (k~10) suppression is genuine and, net of the
   offset difference (≈ +1.0% at 3.2<k<6.4 in the A/B; unmeasured at k>6),
   larger: ≈ -3.4% at k~4 and probably ≈ -5% at k~10 (estimates); the lesson
   about late da/a at dx~0.125 stands. The converged 2LPTIC-vs-fid figure
   +3.9% at 0.3<k<1 contains the 244-step offset at those k, ≈ +1.7% (not
   the linear 1.2%), so ≈ 1.011 x 1.011: a pure ZA-transient part of ≈ +1%
   (estimate).
3. "Default-schedule run (163 steps)": +0.6-1.1% at k<0.8 vs 244 steps is
   mostly the offset ratio: 1.0180/1.0120 = 1.006 at linear k, ≈ +0.7% at
   0.2<k<0.4 and ≈ +0.8% at 0.4<k<0.8 at z=0 (estimate), leaving ≲ +0.3%
   for stepping; the 0.8<k<3.2 row (+0.6-0.9%) is likewise ≈ the offset
   (+0.6-0.8%); the "residual +0.5-0.8% by da^2 scaling" was the 244-run's
   own 1.2% offset (the error scaled as da, not da^2). The high-k comparison
   and hence the da0=6e-4 recommendation at dx~0.125 stand (the offset is
   positive at every measured k, so net of it the suppressions relative to
   244 steps are ≈0.5 point larger for 163 steps and ≈1 point larger for
   123 steps; estimates, offset unmeasured at k>6). "Accurate to ~1%
   everywhere" should therefore read "within ~1-2% of the 244-step run at
   k>6". The stepping-matched IC-switch offset (+3.5% at k~0.3, +5.6% at
   k~5) contains ≈ +2.2% (k~0.3) and ≈ +1.5% (k~5) from this bug
   (estimates), leaving ≈ +1.3% and ≈ +4%.
4. "Halo mass function": the mild step dependence at cluster masses
   (+4-5% at 123 steps vs +2% converged) is partly the ≈1.2% P offset
   difference propagated through the HMF (Sheth-Tormen response: ≈ +0.4% at
   log10M=13.55, +1.0% at 14.05, +2.0% at 14.55; estimate), and the
   converged +2% itself contains the 244-step run's own offset (≈ +0.5%,
   +1.1%, +2.2% at the same masses), so the cluster-mass excess over the fid
   boxes is ≲1% after correction and not established. The step-independent
   low-mass excess (+8..+39%) is unaffected (response ≲0.2%, of opposite
   sign).

Follow-ups (not started): re-run the three fid2LPTIC_L512 comparisons with
`half` (~3 x 10 node-h) and re-derive the ZA-vs-2LPT offsets and the HMF
ratios; re-check the da0 recommendation from the high-k side only. To
reproduce any pre-fix run exactly, convert with `ic2pm.exe <IC> <S_vel> sync`
(S_vel = 1.0 for every run listed above; their submit scripts pass no third
argument and must have `sync` appended).
