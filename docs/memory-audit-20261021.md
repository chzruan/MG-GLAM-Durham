# Memory audit of the f(R) L1024 Np2048 Ng4096 run with all extra outputs

Audit and reduced-size measurement of 2026-10-20/21 for the brief
`docs/memory-budget-brief-20261020.md` (step 1 and 2 of its first prompt).
Code: branch `cz` at 9b88e4a (Fortran last changed in 56e0d07). Sections 1-4 were written before any code change; sections 5-6 record the implementation on branch `memory-cuts-20261021`.

**Units.** Everything below is in GiB (2^30 B). The brief mixes units: its array
sizes are decimal GB (4096^3 x 4 B = 275 GB = 256 GiB) while its measured numbers
are GiB (`sinfo` 1 030 000 MiB = 1005.9 GiB; job 11448712 MaxRSS 1 033 154 828 KiB
= 985.3 GiB). In one unit the "985 < 1031" paradox disappears: particles 192 +
3 x 256 = 960 GiB nominal, 985.3 GiB measured, and all three grids are resident.

## 1. What the code does (verified sites)

- Step (`PMP2main.f90` 41-97): DENSIT -> scalar_field -> POTENTfft5 -> POTENT_FOR_MG
  -> MOVE -> ADDTIME. Checkpoint and analysis outputs are bracketed by
  `MG_scratch_release` / `MG_scratch_ensure` (lines 73-86): FI3 is deallocated
  before WriteDataPM and Analysis and re-allocated afterwards. FI2 is kept across
  steps (warm start: `scalar_field` rescales it by u_bg at the start of each step
  and never resets it). Confirms the brief.
- The f(R) solver allocates nothing (`grep allocat` in PMP2extradof.f90 and
  PMP2MGsolver_fR.f90 is empty); POTENTfft5 uses per-thread real*8 work of 16384.
- FI3 at the fine level: `calculate_residual_fR(0)` (PMP2MGsolver_fR.f90 374-562)
  writes FI3(1:N,1:N,1:N) = residual and reduces RES. It is called twice per
  V-cycle (V_cycle(0) before the restriction, and in scalar_field after the cycle
  for the convergence print). `restrict_residual_fR(1)` (607-700) recomputes the
  restricted residual from FI and FI2 and never reads FI3; nothing else reads FI3
  at level 0 for f(R). The full-grid store is dead data (4 x 256 GiB of writes per
  step at 4096); only the scalar RES is used.
- FI3 coarse layout (level l >= 1, n = N/2^l, ioff = joff = n, koff = N - 2n,
  koff2 = N - n), six n^3 sections per level, named here by content:
  restricted residual FI3(i, j, koff+k); restricted density FI3(n+i, j, koff+k);
  working scalar field FI3(n+i, n+j, koff+k); restricted-copy of the finer field
  FI3(i, n+j, koff+k); physical RHS FI3(i, n+j, koff2+k); level residual
  FI3(n+i, j, koff2+k). Level 1 spans i, j, k up to N (k <= N/2 filled by four
  sections; k > N/2 holds the RHS and residual in two quadrants, levels >= 2 nest
  in the quadrant i, j <= N/2 and the quadrant i, j > N/2 is never used).
  Used: 6/7 N^3 = 219.4 GiB at 4096. Never used: N^3/7 = 36.6 GiB.
  => The brief's "coarse levels need ~N^3/7, FI3 could shrink to ~40 GB, saving
  ~235 GB" counts one field per level; there are six. Realistic FI3 savings:
  36.6 GiB (relayout), up to 146 GiB (relayout + three sections removed, below).
  The restricted residual is dead once the RHS is formed (RHS reads it
  cell-locally, `calculate_physical_right_hand_side_fR` 860); the restricted copy
  is R(u_finer) with u_finer unchanged while on coarser levels (recomputable at
  `correct_solution`); the level residual is consumed only by the next
  restriction.
- THP is `always` on the COSMA image, so leaving the unused octant untouched would
  not reliably keep it non-resident; shrinking FI3 needs a relayout.
- Primary output (`PMP2mod_analyze.f90` Analysis, OutClass = 1): LensPlanes
  (acc Nl^2 x Ns real*8 + plane), GetPower(0) in place (POWERfft5 dumps `deltam`
  with C(0:244,-244:244,-244:244) complex*8 = 0.44 GiB at k_c = 1.5),
  ExtraModes (one G(Nd^3) real*4 buffer, three momentum components in turn,
  in-place FFT, C = 0.06 GiB), WriteSubsample (q(nsel) int16), BDM.
  Station output (OutClass = 2): LensPlanes and GetPower only.
- BDM (`PMP2linker.f90` 136-237, rewrite of 2026-09-06 c68c22d, so the node-2
  pilots used the same structure): DENSIT; FindMaxima (22 words per maximum);
  deallocate FI; RescaleCoords moves the originals to BdmPM* and allocates a NEW
  rescaled set (+6 Np words); AddBuffer allocates imageCounts (1 B/particle), then
  the buffered set (6 real*4 + int64 id = 32 B x iPartMax) and frees the rescaled
  set; SizeList allocates Lst(iPartMax) + Label(cells) int64 after the check
  `8 (iPartMax + cells)/2^30 + Memory(0) > MaxMemory -> error stop`, MaxMemory
  defaults to 500 GiB (ReadBdmConfig 271) and gen_inputs.py does not write it;
  List (zCell16 transient); ParametersDistinct keeps BoundParticleIds%ids (8 B per
  bound particle, per maximum, overlapping spheres counted repeatedly) until
  ReleaseMaxima; WriteFiles (+ fR_halo reads FI2); RemoveBuffer; FI re-allocated.
  dBuffer = max(5, 15 x 2 Box/Ng + 0.75 Box/Ng) = 7.6875 Mpc/h at the 0.25 Mpc/h
  cell; image factor (1 + 2 dB/Box)^3 = 1.0457 at L1024, 1.1912 at L256; cells
  (ceil((Box+dB)/Cell) + 1 - floor(-dB/Cell) + 2)^3 = 2083^3 at L1024, Cell = 0.5.

## 2. Resident arrays at Ng = 4096, Np = 2048^3 (GiB)

Force step (peak inside the MG solve; nothing else is allocated):

| array | GiB | note |
|---|---|---|
| XPAR, YPAR, ZPAR, VX, VY, VZ | 192.0 | 6 x 2048^3 x 4 B |
| FI (density, then potential) | 256.0 | |
| FI2 (scalar field, warm start) | 256.0 | |
| FI3 (multigrid scratch) | 256.0 | dead fine-level store over the whole array; coarse sections 219.4; never used 36.6 |
| sum | **960.0** | measured 985.3 (sacct, cgroup incl. page cache, old BDM); the excess is the finder's resident residue (5-6 GiB at 128 threads, section 5) plus page cache charged to the cgroup |

Primary output (OutClass = 1), resident total after each change, in program order:

| phase | change | total |
|---|---|---|
| entry, FI3 released | -256 | 704.0 |
| LensPlanes acc(4096,4096,4) real*8 + plane | +0.56, freed | 704.6 |
| GetPower, POWERfft5 in place, deltam C | +0.44, freed | 704.4 |
| ExtraModes G(2048^3) real*4 (N_dep_mom = 0) | +32.0, freed | 736.0 |
| WriteSubsample q(nsel) int16 | +0.16, freed | 704.2 |
| BDM FindMaxima, 22 words x ~1.1e7 maxima | +0.9 | 704.9 |
| BDM deallocate FI | -256 | 448.9 |
| BDM RescaleCoords, second particle set | +192.0 | 640.9 |
| BDM AddBuffer imageCounts | +8.0 | 648.9 |
| BDM AddBuffer buffered set, 32 B x 1.0457 Np | +267.7 | **916.6** (peak of the output step) |
| free scratch and the rescaled set | -200.0 | 716.6 |
| SizeList check: 134.3 + 716.6 = 850.9 > MaxMemory 500 | | **error stop** (LCDM: 594.9 > 500, also stops) |
| Lst(iPartMax) + Label(2083^3) int64 | +134.3 | 850.9 |
| List zCell16 | +16.7, freed | 867.6 |
| ParametersDistinct bound ids, 8 B x Sum n over maxima | +15.5 (distinct haloes, 24 % of Np) to ~30 (overlap, from the +0.5 GiB measured at L256) | ~866 to ~881 |
| swap to the maxima lists | -134.3 | ~737 |
| WriteFiles + fR_halo (reads FI2) | | ~737 |
| RemoveBuffer, FI, FI3 | | 448.9 -> 704.9 -> 960.9 |

Station output (OutClass = 2): 704.0 + 0.56 (LensPlanes) = 704.6; GetPower in place.
Checkpoint: 704.0 + 0.02 (Xb buffers).

Verdict on the brief: FI3 is free at the output step (confirmed); the run's peak is the
force step (confirmed for the current code, 960 vs 916.6 nominal), but the output step
is within 44 GiB of it because of the finder's transient third particle set, not because
the extras are small; the extras (G 32 GiB) sit 180 GiB below the peak and N_dep_mom is
not a lever on the peak. The LCDM twin's output step (660.6) exceeds its force step (448).

## 3. Measurement: f(R) twin of the L256 pilot

`DEGRACE-MGGLAM/L256Np512Ng1024/Planck2015_FR`, job 12193314 (cosma8-serial, 32
cores, 40 GB, 1 h 22 min, CPU efficiency 93 %), realisation 0 (`Run0`, tag 0000),
same ICs, `gen_inputs.py ... --gravity fr --fR0 -1e-5 --fr-halos`, binaries of the
LCDM pilot plus the realisation-aware ic2pm of 9b88e4a; `check_run.py` 0 FAIL 0 WARN
(80 189 haloes at z = 0; MG solver 47 of 81 min). PMP2main's VmRSS/VmHWM sampled
every 2 s (`Run0/rss.log`), aligned with main.log's `time =` markers
(`phase_rss.py`, session scratchpad). Two earlier attempts (12193293, 12193309)
died in the submit script's ic2pm guard: the usage line wraps at 80 characters, so
`grep realisation` never matches; both submit scripts now test for `0..9999`.

| phase (L256, f(R)) | nominal | code counter (Run.log) | measured VmRSS |
|---|---|---|---|
| force step, before the first primary | 15.00 | 15.001 | 15.03 |
| force step, after the 17th primary | 15.00 | 15.001 | 16.30 (= VmHWM = time -v) |
| output entry (FI3 released), z = 0 | 11.00 | 11.001 | 12.31 |
| ExtraModes with G(512^3) | 11.50 | not counted | 12.80 |
| BDM AddBuffer transient (0.5 s) | 14.89 | 14.909 | not resolved; VmHWM shows it stayed below 16.30 |
| BDM list phase + zCell16 | 14.49 | 14.495 | 15.2 to 15.3 |
| BDM ParametersDistinct (23 s) | 14.2 + ids | not counted | 15.42 -> 15.70 (ids grow by ~0.5) |
| station output plateau | 11.00 | 11.001 | 11.0 + residue |
| sacct batch MaxRSS | | | 22.1 (cgroup, includes page cache) |

Findings: (i) the array sum predicts the counter to 0.02 GiB (AddBuffer 14.909 vs
14.89; lists 2.412 vs 2.411, also in the LCDM pilot); (ii) the solve itself has no
overhead (15.03 vs 15.00); (iii) every BDM call leaves a resident residue that is never
released: 15.03 -> 15.09, 15.42, ..., 16.30 over the 17 primaries (the last
plateau is the one after the 16th catalogue, steps 174-183; no force step
follows the 17th), +1.27 GiB in total, growing with the halo population, and
present during every later force step. The node-2 f(R) pilot (64 threads, 8x the particles, 10x the haloes) shows
+3.0 GiB (123.0 measured vs 120.0 nominal). Extrapolated to production: +7 (same
ratio) to +24 (linear in particles) GiB on top of 960, i.e. a force-step RSS of
967 to 984 GiB on a 1005.9 GiB node.

Residue mechanism, first attempt (superseded by section 5: the residue is the
cache of libiomp5's allocator and is size-driven; the per-thread scaling below
was coincidental). Diagnostic rerun 12193635, `Planck2015_FR_diagmalloc`, same
inputs and binary, submit script adding `MALLOC_MMAP_THRESHOLD_=65536
MALLOC_TRIM_THRESHOLD_=65536 MALLOC_TOP_PAD_=0`): the RSS trace reproduces the
base run to 0.01 GiB at every output (time -v 16.30 GiB; force-step plateau after
each primary 15.41, 15.47, 15.49, ..., 16.22, 16.31), so the residue is neither
top-of-heap padding nor large blocks held below the mmap threshold. It is not a
Fortran leak either: PMP2linker.f90 has 23 allocate and 24 deallocate sites,
ReleaseMaxima frees the 22 per-maximum arrays and BoundParticleIds (whose ids
components go with it), the output step frees acc, C, G, dens and the Xb..VZb
copies, and the finder has no automatic arrays. What is left is that the
per-maximum temporaries of GetHalo (rows, radii, potential, ids) are allocated
inside OpenMP parallel loops, i.e. in per-thread glibc arenas (glibc 2.34, 64 MiB
heaps, THP `always` on the image): chunks freed inside those heaps are returned
to the kernel only by an explicit `malloc_trim(0)`, and a partially used huge
page stays resident. Consistently, the residue scales with the thread count and
not with the particles or haloes: 1.30 GiB / 32 threads = 41 MB per thread at
L256, 3.0 GiB / 64 threads = 48 MiB per thread in the node-2 pilot (8x the
particles, 10x the haloes). At 128 threads that would be 5 to 6 GiB,
so the force step at production is ~967 GiB, and the "+24, linear in particles"
bound is unlikely. The rest of the baseline's 985.3 is page cache charged to the
cgroup, which `sacct` MaxRSS includes (22.1 vs 16.3 GiB here). The decisive test
is cut 1 itself (malloc_trim after ReleaseMaxima, rerun L256); the same rerun
should sample AnonHugePages and the arena mappings from `/proc/<pid>/smaps_rollup`
to locate whatever remains.

## 4. Proposed cuts, in order of return on the run's peak at 4096 (GiB)

0. **BDM `MaxMemory`** (0 GiB, mandatory): write `MaxMemory` (node size) into
   BDM.config from gen_inputs.py, or default it to the node; without it both
   production twins stop at their first primary output (850.9 / 594.9 > 500).
1. **Release the finder's heap residue** after each BDM call: about -7 on the
   force step at 128 threads (41 to 48 MiB per thread measured), no numerical
   effect. Implemented as `glam_kmp_alloc.c` (Fortran ALLOCATE through glibc
   instead of libiomp5's caching allocator, see section 5) plus TrimHeap
   (`malloc_trim(0)` after the finder); the submit-script tunables tested in
   12193635 and malloc_trim alone change nothing.
2. **FI3 relayout** into per-level coarse arrays and removal of the dead fine-level
   store: -36.6 (force step 923.4). Bit-identical by construction (same arithmetic,
   different addresses); touches every index expression of the f(R) routines and
   the shared restriction/prolongation.
3. **Coarse-section elimination**, -36.6 each, up to -109.7: (a) RHS formed in
   place of the restricted residual; (b) restricted copy recomputed at the
   prolongation; (c) level residual fused into the restriction (explicit real*4
   rounding to stay bit-identical). Force step 813.7 with all three. Beyond
   -44 GiB the output step's AddBuffer (916.6) caps the gain unless 4 is done.
4. **BDM rescale straight into the buffered arrays** (no intermediate second
   particle set): -192 at AddBuffer (916.6 -> 724.6); the output-step peak then
   becomes the list/ParametersDistinct phase at ~866-881 (+ residue). Catalogues
   bit-identical if the same expressions are used.
5. `N_dep_mom = 1024` (-28 at the ExtraModes point, 736 -> 708) or depositing the
   momentum into FI by sequence association (-32): 0 on the peak; harmless.
6. Lens planes, deltam, subsample: < 1 GiB; nothing to gain.

Resulting force-step / output-step peaks (nominal, plus the residue of 1): today
960 / 917; with 0+1 960 / 917 minus the residue; with 0+1+2 923 / 917; with
0+1+2+3+4 814 / ~880.

## 5. Implementation (branch `memory-cuts-20261021`, 2026-10-21)

User decision 2026-10-21: implement cuts 0, 1, 2 and 4 (not 3, not 5). Commits
on cz 9b88e4a: ed6cb3a (the cuts), 966f74c (review fixes, no numerical change),
d9bdf95 (the allocator fix of cut 1 and this note), then the round-2 review
fixes of section 6.

- Cut 0 (`PMP2linker.f90` ReadBdmConfig, DefaultMaxMemory): without a
  `MaxMemory` line in BDM.config the admission limit is the job's memory
  (`SLURM_MEM_PER_NODE`, or `SLURM_MEM_PER_CPU` x `SLURM_CPUS_ON_NODE`), else
  the node's MemTotal, else the old 500 GiB; the start-up check prints the
  resolved value. `gen_inputs.py --bdm-max-memory GIB` writes the key for
  older builds. A generated default BDM.config carries a commented template
  line instead of the resolved value.
- Cut 1 (`PMP2linker.f90` TrimHeap): `malloc_trim(0)` through `bind(C)` after
  RemoveBuffer, before FI is re-allocated.
- Cut 2 (`PMP2mod_tools.f90`, `PMP2MGsolver_fR.f90`, `PMP2extradof.f90`,
  `PMP2main.f90`): type `MGLevelArrays` with `res, rho, ucopy, u, rhs, lres`
  (the six old FI3 sections) and `MGL(1:levelmax-levelmin)`, allocated and
  zeroed by `MG_levels_allocate` (first touch in parallel), freed by
  `MG_levels_release`; `MG_scratch_ensure/release` dispatch on MG_model, so
  the checkpoint and output steps release the per-level arrays as they
  released FI3. Every coarse-level loop of the f(R) solver is an internal
  procedure with explicit-shape dummies; the PM-grid loops are untouched
  except that the dead residual store is gone. The shared restriction and
  prolongation routines got an `MG_model.EQ.1` branch; models 2-5 keep FI3
  and their old code. `levelmax/levelmin` are computed once in Tools
  (`MG_set_level_bounds`) for both the allocation and `initialise_extradof`.
  Side finding of the review: `tmpt` in `correct_solution` was a routine-scope
  variable left shared under `DEFAULT(SHARED)`, a formal data race in the old
  code for every model (benign in practice when the compiler keeps it in a
  register); it is private now in all three branches.
- Cut 4 (`PMP2linker.f90` RetainParticles, AddBuffer(rescaleRetained),
  BdmScaledCoordinate/Velocity): the inline finder moves the simulation arrays
  to BdmPM* and builds the buffered set straight from them; the rescaled copy
  (6 words per particle) no longer exists on the output step. RescaleCoords +
  AddBuffer without the flag remain for the standalone fixtures.

Verification:
- `scratch/memory-cuts-20261021/addbuffer_check`: a driver linked against the
  new objects runs both AddBuffer paths on 300,000 pseudo-random particles
  (with edge values) at 1 and 4 threads: buffered sets bit-identical, exact
  restoration, memory counter balanced. PASS.
- Tiny twins (L = 256, 128^3 particles, Ng = 256, FML 2LPT ICs seed 2026
  in `scratch/memory-cuts-20261021/ICs_L256Np128_seed2026`, one thread so
  that a run is deterministic; job 12197344, old = cz 9b88e4a build, new =
  ed6cb3a build): LCDM: all 146 binary products and 230 text files identical;
  the 17 finder logs (CATALOGS/outputB.*.dat) differ in their timing lines only.
  f(R): identical through step 108 (two catalogues, checkpoints at 50 and
  100, the scalar-field dump of step 1 with MG_test = 1), then last-bit
  differences that grow: P(k) agrees to 6 significant digits at all 42 outputs
  (the last printed digit of one column differs in 32 of 9576 rows; 60 of the
  84 PowerDM files are bitwise identical), lens planes at z = 0 within a
  maximum relative pixel difference of 5.7e-6, the z = 0 catalogue with the
  same 1505 haloes, the same particle counts and masses, positions within one
  float32 ulp (1.5e-5 Mpc/h). Diagnosis (jobs 12197470, 12197471): both
  binaries are deterministic (repeat runs identical), insensitive to
  MALLOC_PERTURB_ (no uninitialised reads) and to the mmap threshold (no
  placement sensitivity); the new binary without TrimHeap and with the old
  RescaleCoords+AddBuffer path reproduces the new results exactly, so the
  difference is cut 2 alone: the same arithmetic compiled on explicit-shape
  dummies instead of FI3 sections is not code-generation-identical under
  `-O3 -fp-model fast=1` (vector SVML vs scalar libm in remainders, FMA
  contraction). Strict test (job 12197553): cz 9b88e4a and the branch
  (966f74c) both built with `-O2 -fp-model strict -no-vec -no-simd` (no
  reassociation, no FMA contraction, scalar libm only) give bit-identical
  results for the whole tiny f(R) run: all 161 binary products and the
  particle checkpoint PMcrs0.DAT identical, all text products identical; the
  only differing data is the 4-byte AEXP0 word of the checkpoint header
  PMcrd.DAT, which GLAM never sets (documented in ic2pm.f90), besides the
  timing lines of the finder logs. So the source arithmetic of cut 2 is identical; the
  last-bit differences of the production build are code generation, of the
  same kind as the difference between the production and the strict build of
  cz itself (every product differs there).
- L256 Np512 Ng1024 twins with the new binary against the cz reference runs
  (f(R) job 12197341 vs 12193314, Run0; LCDM job 12197342 vs 12192772, Run1;
  `Planck2015_*_cuts20261021`, 32 threads): P(k) agrees to 6 significant
  digits at all 42 outputs in both models (the last printed digit of one
  column differs in 427 of 26,964 rows for f(R) and 393 for LCDM, maximum
  relative difference 1.5e-5 and 1.9e-5; 21 of 84 PowerDM files bitwise
  identical in each); lens planes at z = 0 within a maximum relative pixel
  difference of 3.3e-4 (f(R), rms of the relative differences 6.6e-6) and
  9.7e-5 (LCDM, rms 4.1e-6); z = 0 catalogues 80,189 vs 80,188 haloes (f(R))
  and 78,046 vs 78,045 (LCDM), 99.99-100 % matched within 0.05 Mpc/h,
  98.9-99.4 % with identical particle counts, 0.23 % (f(R)) and 0.13 % (LCDM)
  with bound masses differing by more than 1 %, N(> 1e13) and N(> 1e14)
  identical for both mass definitions. This is the thread-nondeterminism level
  of the code (atomic CIC), the same criterion as the confirmation pilot of
  2026-10-18. Wall time: f(R) 83 vs 82 min; LCDM 75 vs 28 min because node
  m8177 was contended (every phase, including the untouched LCDM force code,
  was slower from step 1; timing.log: force 46 vs 11 min, BDM 23 vs 11).

Memory at L256 (GiB, VmRSS from the 2 s sampler; nominal = counted arrays):

| phase | old (cz) | new (cuts) | nominal new |
|---|---|---|---|
| f(R) force step, early plateau (steps 20-50) | 15.03 | 14.45-14.50 | 14.43 |
| f(R) force steps 103-104, right before the first catalogue | 15.09 | 14.51 | 14.43 |
| f(R) force step after the 16th catalogue (steps 174-183) | 16.30 | 15.72 | 14.43 + residue |
| f(R) output step, maximum (list / ParametersDistinct phase, z = 0) | 15.70 | 15.71 | |
| f(R) AddBuffer transient at the first catalogue | 14.89 (nominal) | 11.9 (nominal) | |
| f(R) time -v maximum RSS | 16.30 | 15.71 | |
| LCDM force step before / after the catalogues | - / - | 7.09 / 8.22 | 7.00 |
| LCDM output step, maximum | 12.1 (time -v) | 11.59 | |

Cut 2 delivers the predicted 4/7 GiB (0.57) at L256, i.e. 36.6 GiB at
Ng = 4096. Cut 4 removes the rescaled copy (3 GiB at L256, 192 at
production) from the AddBuffer transient; at L256 the output-step maximum is
the list phase, so the peak itself is unchanged there, as predicted. Cut 1
as first implemented (malloc_trim after the finder) does NOT remove the
residue: the force-step plateau still rises by +0.34 after the first
catalogue (14.51 -> 14.85 at step 107) and by ~0.05 per later one (+1.22 over
the run, old +1.27). Cut 2 measured: 0.53 (early plateau) to 0.58 (steps
103-104) against 0.57 predicted.

Residue mechanism (LCDM L256 diagnostics, 2026-10-21):
- Mapping classes sampled from /proc/<pid>/smaps every 5 s over the whole run
  (job 12197761, `Planck2015_LCDM_smaps20261021`): the growth is almost all
  in ordinary anonymous mappings (11.18 GiB of 11.57 at the end); glibc arenas
  (64 MiB-aligned heaps) hold 0.24 GiB at the end, thread-stack-sized mappings
  0.11, the brk heap 0.03, the main stack nothing. Full smaps snapshots before and after the first
  catalogues (jobs 12198124 and 12198125, `Planck2015_LCDM_smapsdump20261021`
  and `_heaparrays20261021`): after the first catalogue a set of anonymous
  read-write regions of fixed sizes stays fully resident (128, 100, 24, 24,
  18, 16, 14 MiB), and later catalogues add or enlarge them (136, 100, 40
  (36 resident), 24 x 3, 18, 17, 16, 14, 12 MiB after four; the exact set
  differs from run to run); building the finder with `-heap-arrays` changes
  nothing, so it is not stack temporaries.
- `strace -f -k` on mmap (job 12198211, `Planck2015_LCDM_strace20261021`,
  `diagnostics/strace/strace.txt`): the anonymous regions that appear at the
  first catalogue are mmaps whose recorded stack is `for_allocate_handle`
  (the Fortran ALLOCATE, libifcoremt) -> `scalable_aligned_malloc` ->
  `rml::internal::MemoryPool::getFromLLOCache` -> `mallocLargeObject` ->
  ... -> `__mmap`, all inside libiomp5.so (the `kmp_aligned_malloc` stub in
  between is a tail jump and leaves no frame). The link from libifcoremt to
  libiomp5 is the symbol table: libifcoremt imports kmp_aligned_malloc,
  kmp_realloc and kmp_free as weak references, libiomp5 exports them, and
  libiomp5 does not export the allocator's clean-up call
  (scalable_allocation_command). In that run the regions still mapped after
  the first catalogue were mid-sized blocks of up to 144 MiB (128, 64, 60,
  32, 22, 14 MiB resident); the largest allocations of the finder (576 MiB to
  4.5 GiB: Lst, Label, the buffered particle set) were unmapped again. So the
  residue is the allocator's large-object cache holding freed mid-sized
  blocks (per-halo gathers and the finder's smaller work arrays), which it
  ages out only on later allocation activity that the force steps do not
  provide; the blocks are reused by the next catalogue, so the plateau only
  grows when a later catalogue needs bigger ones. glibc's malloc_trim never
  sees this memory.
- Fix (commit d9bdf95; since the round-2 fixes linked into every Fortran
  executable, before that only into the OBJ targets, PMP2MG and
  PMP2spectrum): `glam_kmp_alloc.c` defines kmp_aligned_malloc / kmp_malloc /
  kmp_calloc / kmp_realloc / kmp_free forwarding to glibc (posix_memalign,
  malloc, calloc, realloc, free). By ELF interposition the executable's
  definitions win over libiomp5's, so ALLOCATE goes to glibc: large arrays
  are mmap'd and unmapped on deallocation, the rest is released by TrimHeap.
  libiomp5 does not call these entry points internally (verified in review
  round 2 by objdump: no PLT call or relocation to the five names). Results
  are unchanged (allocation placement was shown above not to affect the
  numerics).
- Measured effect (job 12198268, `Planck2015_LCDM_glibcalloc20261021`, same
  snapshots): VmRSS 7.10 GiB at step 100, 7.09 GiB at steps 106, 109, 111
  and 114, i.e. no residue after four catalogues (the libiomp5 allocator in
  job 12198124: 7.09 -> 7.50 GiB); glibc arenas hold 17 MiB, the brk heap 40
  MiB; finder timings unchanged within node noise (ParametersDistinct 16.1 s
  vs 15.5 s, from the deleted diagnostic logs). The tiny single-thread f(R)
  twin with this build (md5 4b97717f, the repo-root build of 19:10 from the
  sources committed unchanged in d9bdf95 at 19:36; job 12198269) is bitwise
  identical to the ed6cb3a build (Run_new) in all 163 binary and 247 text
  products (the 17 finder logs differ in timing lines), so 966f74c and the
  allocator together change no result.

Expected peaks at L1024 Np2048 Ng4096 with cuts 0, 1, 2, 4 (GiB, node
1005.9): force step 923.4 with no finder residue (the per-thread solver and
OpenMP overheads measured at L256 are below 0.1 GiB); output step 724.6 at
AddBuffer and ~866-881 in the list / ParametersDistinct phase (unchanged by
the cuts, per-halo bound ids included); checkpoint 704. The run's peak is the
force step, 82 GiB below the node, against 960 + 7-24 before the cuts. The
production-size test (~2300 core-h) has not been run.

Not tested on this branch: the other gravity models (MG_model 2-5 keep the
FI3 branches untouched but no run was made); restart from a checkpoint with
the per-level arrays; the MPI programs (PMP2mainMPI, PMP2startMPI link the
allocator object but were not built); the standalone finder PMP2BDM; f(R)
with the glibc-allocator build at more than one thread (the 32-thread memory
measurement is LCDM, the f(R) checks are single-thread tiny runs); 128
threads; the `--mem-per-cpu` default; the production size. Open risks: the
interposition relies on libifcoremt's weak references and a dynamic libiomp5
of intel_comp/2024.2.0; a static Intel runtime or another toolchain restores
the caching allocator (static linking fails loudly, another toolchain does
not; `make check-alloc-imports` after every toolchain change).

How to re-run the verification (all under `scratch/memory-cuts-20261021`):
`tiny_job.sh` (old vs new single-thread twins), `strict_job.sh` with
`build_strict.sh` (strict builds), `tiny_glibc_job.sh` / `tiny_round2_job.sh`
(a new build against Run_new), `diagnostics/*/submit.sh` (the residue
diagnostics), `scripts/compare_runs.py A B` (bitwise product comparison),
`scripts/compare_l256.py A B` (P(k), lens planes, catalogue at z = 0),
`scripts/phase_rss.py RUNDIR` (per-phase memory from rss.log),
`scripts/smaps_diff.py BEFORE AFTER`, `scripts/strace_attr.py`; the suite
twins' submit.sh under `DEGRACE-MGGLAM/L256Np512Ng1024/Planck2015_*_cuts20261021`.

## 6. Pre-merge reviews (2026-10-21, four fresh sessions, prompts in `scratch/memory-cuts-20261021/review_prompts.md`)

Reviewer A (solver relayout): no Blocker or Major. Verified the F- and W-cycle
drivers, MG_test paths, restart and end-of-run allocation lifecycle, grid sizes
with odd factors (1536, 768), NGRID = 2 and 4, the OpenMP private lists, the
explicit-shape dummies (no copy-in/out, no aliasing) and that nothing reads
the dropped fine-level store. Nits: the per-ensure zeroing of the per-level
arrays was unnecessary work (219 GiB per output step at Ng = 4096) and doubled
the counter lines in Run.log; MGL(nlev)%lres is never used (negligible).
Fixed: the arrays are zeroed at the first allocation only, with one counter
line per allocation and release. Provenance caveat recorded: the production
build of the branch is not bitwise continuous with cz (section 5), so a run
resumed from a cz checkpoint with the new binary, or compared with a cz
reference, differs at the float32 last bit; record the binary per segment.

Reviewer B (finder): one Major. The BDM-refine harnesses extract finder
routines by name from PMP2linker.f90 and did not know RetainParticles and the
two scaling functions, so they no longer built. Fixed: the three extraction
lists (`BDM-refine/repairs/20260906/tests/particles_tests.py`,
`.../20260906-review-followup/particles/run_checks.py`,
`.../analysis/full-audit-20260906/run_audit.py`) now include them when
present; with gfortran 14.1 the two particle harnesses pass (104 and 256
experiments, the follow-up one comparing the pinned baseline and the branch
bit for bit at 1, 2 and 4 threads). The AddBuffer bit-identity driver is now
tracked as `test/addbuffer_check/` (PASS at 1 and 4 threads).
`BDM-refine/repairs/20260906-review-followup/build-config/regression.py`
stays as an archival harness: it asserts the former 500 GiB default and
cannot extract functions; it was already unbuildable before this branch.
Minors fixed: the start-up line and the two admission stops now print the
source of the default (`SLURM_MEM_PER_NODE`, `SLURM_MEM_PER_CPU x CPUs`,
`/proc/meminfo MemTotal` or the compiled 500) and the counted and required
GiB, and say that the finder counts nominal arrays only; gen_inputs writes
MaxMemory with `%g` and its help text says which builds accept the key
(2026-09-06 to 2026-10-20 default to 500 GiB, earlier ones reject it);
DEGRACE-companion/README.md notes that the pilots' `--mem` is now their
limit. The glibc dependency moved out of the Fortran: TrimHeap calls
`glam_trim_heap` in glam_kmp_alloc.c (malloc_trim under glibc, no-op
elsewhere).

Reviewer C (allocator interposition): one Major. Nothing checked that the
installed Fortran runtime imports only allocator entry points the interposer
covers; a newer runtime could allocate through another name and the first
DEALLOCATE would die in the wrong allocator. Fixed: `make
check-alloc-imports` lists the allocator imports of the libifcoremt the
executable links (readelf) and fails on any name outside the known set
(kmp_aligned_malloc, kmp_realloc, kmp_free, which are interposed;
__kmpc_aligned_alloc, __kmpc_alloc, __kmpc_free, kmpc_sharable_malloc,
kmpc_sharable_free, used only for `!$omp allocate` directives, which MG-GLAM
does not use; _mm_malloc, _mm_free of libintlc); it also checks that the
executable defines the five kmp_* symbols. Result today: all known. Verified
by the reviewer on the binaries: the executable's definitions are GLOBAL
DEFAULT and first in the lookup scope, libiomp5 has no DT_SYMBOLIC, its
exported kmp_* are stubs onto internal kmpc_* and it never calls the five
names through the PLT, the runtime's I/O buffers and temporaries use glibc
directly, and libifcoremt's `FOR_DISABLE_KMP_MALLOC` switch is read but does
not disable the path. Minors fixed: glam_kmp_alloc.o is now linked into every
Fortran executable (PMP2init, PMP2start, PMP2density, PMP2correlation,
gadget2pm, ic2pm were missing); the TrimHeap comment no longer blames glibc
arenas; the C file documents the five symbols, what does not go through them,
the static-linking failure mode and the check. Noted, not changed: every
ALLOCATE now takes glibc's posix_memalign path (alignment 32 requested by the
runtime) instead of the scalable allocator's thread caches, and blocks above
32 MiB are page-faulted fresh at every catalogue; the only timing so far is
ParametersDistinct 16.1 s against 15.5 s at L256 (4 %, within node noise).
Check the finder's timing at the first production-size catalogue.

Reviewer D (documentation and claims): one Blocker and four Majors, all
in the text, none in the code. "P(k) identical to the printed digits" was
false (the last printed digit of one column differs in a few hundred rows;
corrected above with the counts); kB had been written as GB for the
glibc-allocator measurement (corrected to GiB); the strace paragraph had
overstated what the trace shows (corrected: the trace gives the libiomp5
chain, the symbol tables give the libifcoremt link, and the cached blocks are
mid-sized, not the largest); the LCDM twin's 75 vs 28 min wall time was
unexplained (node contention, now stated); "linked into every executable"
was only true after round 2 (stated). Minors fixed: the force-step windows
are now named with their measured values, the lens-plane numbers are labelled
as maximum relative pixel differences, the mass definition of the 1 %
fraction is named, the finder logs are excluded from "bitwise identical",
section 3's per-thread model is marked as superseded, the 5-6 GiB
arithmetic, the 594.9 and 48 MiB figures, the gen_inputs help text (a build
with cut 0 rather than "cz after 2026-10-21"), and the "not tested / how to
re-run" paragraphs above. The analysis scripts are now kept under
`scratch/memory-cuts-20261021/scripts`.

Re-verification after the round-2 changes: build 024ca0e6; `make
check-alloc-imports` passes; `test/addbuffer_check/run.sh` PASS; the two
gfortran harnesses pass; the single-thread tiny f(R) twin of build 024ca0e6
against the previous build (job 12198565, Run_round2 vs Run_glibcalloc): all
161 binary products, the particle checkpoint and all text products identical;
the only difference is the unset AEXP0 word of the checkpoint header. Run.log
now carries one counter line per scratch allocation (474 lines instead of
1020 for the run).
