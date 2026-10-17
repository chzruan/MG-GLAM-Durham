# Review prompts for branch `degrace-track-a` (15 commits, `cz..degrace-track-a`)

Six independent, read-only reviews. Each prompt is self-contained: paste the
COMMON PREAMBLE followed by one REVIEW block into a fresh session or subagent.
Findings go to `docs/reviews/degrace-track-a/REVIEW-<k>-<topic>-findings.md`.

---------------------------------------------------------------------------
## COMMON PREAMBLE (prepend to every review)

You are reviewing code on the COSMA HPC system (Durham). Repository:
/cosma8/data/dp203/dc-ruan1/mgglam_claude/MG-GLAM (MG-GLAM, a Fortran
particle-mesh N-body code for LCDM and modified gravity, plus tools). Branch
`degrace-track-a` is checked out; the mainline is `cz` (NOT `main`). Review
range: `git diff cz..degrace-track-a` (merge base 69c9261, 15 commits, see
`git log --oneline cz..degrace-track-a`). The branch implements the "DEGRACE
companion" outputs specified in docs/degrace_companion_implementation_plan.md:
section 5 specifies items I1-I11, section 7 the Fourier and file conventions,
section 8 the tests, section 9 the storage budget; section 16 is the author's
session log. Treat the plan's statements about what was verified as claims to
check, not as evidence.

Hard rules (shared system; the user is personally accountable for what you do):
- Read-only review. Do not edit tracked files, do not commit, do not rebuild
  binaries, do not run simulations or the test suites, do not start background
  processes. Write exactly one file: your findings (path given below). Small
  single-threaded Python snippets to check arithmetic are fine:
  `micromamba run -n cosemu python3 -I -c '...'` (seconds, no numpy threads).
- Never run recursive grep/find/du/ls -R over /cosma*, /snap*, the home
  directory or the data root /cosma8/data/dp203/dc-ruan1. Search only inside
  the repository paths named below (`git grep` in the repo, or `grep -rn` on a
  named subdirectory). Do not read or touch ~/.ssh.
- The FML build directory 2LPTIC_Gui/FML/ is gitignored but present on disk;
  read it when the review needs the FML library sources.
- Use `git diff cz..degrace-track-a -- <paths>` and `git show` to see exactly
  what changed, and read the surrounding unchanged code when correctness
  depends on it.

Deliverable: write the findings file with these sections. (1) Summary, 3-6
lines, ending with a merge recommendation: ready / ready after the listed
fixes / not ready. (2) Findings table with columns ID, Severity, File:line,
Description, Evidence, Suggested fix. Severity: Blocker = wrong results, data
loss or a crash in production use; Major = wrong results in some
configurations, silent misuse, resource blow-up; Minor = robustness or
maintainability; Nit. (3) Verified OK: a bullet list of the specific
properties you checked and found correct, so nobody re-reviews them.
(4) Open questions for the author. Be concrete: cite file:line, quote the
offending expression, give a failing input where you can. Do not pad. In your
final message return only the findings file path, the counts per severity and
one line per Blocker or Major.

---------------------------------------------------------------------------
## REVIEW 1: Fortran extras module (lens planes, Fourier dumps, momentum, f_R at haloes)

Findings file: docs/reviews/degrace-track-a/REVIEW-1-extras-fortran-findings.md

Scope: PMP2mod_extras.f90 (new, module ExtraOutputs) and its call sites in
PMP2mod_power.f90 (DumpModes before the squaring in POWERfft5),
PMP2mod_analyze.f90 (LensPlanes, ExtraModes, WriteSubsample, BDM guard) and
PMP2linker.f90 (fRHaloOpen/Row/Close inside BDM's WriteFiles). Reference code
you will need: POWERfft5 and the 1-D FFT passes in PMP2mod_power.f90, the
FFTPACK5 module (`git grep -n -i 'module fft5' -- '*.f90'`), DENSIT in
PMP2mod_density.f90, the Tools module PMP2mod_tools.f90 (globals NGRID, NROW,
Box, AEXPN, ASTEP, Om, OmL, hubble, Nrealization, ISTEP, particle arrays
XPAR/YPAR/ZPAR/VX/VY/VZ and their coordinate convention [1, NGRID+1)), and
for f_R the MG background and FI2 definitions in PMP2MG_subroutines.f90 /
PMP2mod_MGbackground.f90.

Check at least:
- DumpModes / PackedBasis / FFT3Dforward. Derive independently the
  half-complex packing produced by POWERfft5's three 1-D real FFT passes plus
  the transposition (which index is p=1 / p=2m / p=2m+1 / p=N, the sqrt(2)
  factors on the first and Nyquist entries, the normalisation) and confirm
  that the weight tables in PackedBasis and the loop `val = val +
  wx(i)*wy(j)*wz(k)*G(pz(k), py(j), px(i))` yield delta_k = N^-3 sum_x
  delta(x) exp(-i k.x) with node i at x = (i-1) Box/N, the right sign of the
  imaginary part, the Hermitian relation at iz = 0, and the correct array index
  order of G after the transposition. Check the n <= N/2-1 limit and the k = 0
  handling. Is FFT3Dforward a faithful copy of POWERfft5's passes (compare
  line by line), and does it work on the transient N_dep grid when N_dep /=
  NGRID (work-array sizes, wsave initialisation per N, threadprivate arrays)?
- LensPlanes. Slab and pixel indexing per axis, the documented orientation
  (X-planes indexed (y,z) etc.) versus the actual loops, block-averaging
  factor f and wnorm, mass conservation, real*8 accumulation, OpenMP ownership
  of rows (no races), whether FI holds delta or 1+delta at the call site, the
  write order of the Fortran array versus the .meta description, and the cost
  at NGRID = 2048.
- MomentumDeposit. Mapping of coordinates in [1, NGRID+1) to the N_dep cells,
  CIC weights and periodic wrap, the velocity conversion (compare with how
  GLAM's own analysis converts VX to km/s and at which epoch, a or a - da/2),
  the (1+delta) v normalisation, OpenMP reduction correctness (private
  copies, atomics or critical sections), and memory of the transient grids at
  N_dep = 1024.
- ExtraModes: when it runs (class 1 only), how it avoids recomputing the
  density (mDENSIT, iPowerDone), epochs recorded as a_x and a_v.
- fRHaloOpen/Row/Close: f_R = fR_bg(a_old) * FI2**(n+1), fR_bg = fR0 (R0 /
  R(a))**(n+1), R = 3 (Om/a^3 + 4 OmL). Verify against the code's definition
  of FI2 (what it stores and in which units), the epoch at which FI2 is valid
  when BDM runs (a_old = a - da?), trilinear interpolation with periodic wrap,
  the coordinate units of the catalogue rows passed in (grid units or Mpc/h),
  the fRHaloActive flag lifetime, and row/meta counts.
- GrowthFactor / ComovingDistance / WriteOutputSchedule: formulas, numerical
  accuracy (Simpson intervals, lower limit), real*4/real*8 mixing, c/H0 =
  2997.92458, the class attached to the right step after the permutation.
- Fortran unit numbers 29/30/31 versus other units opened anywhere in GLAM
  (`git grep -n -i 'open *(' -- '*.f90'`), file-name collisions, status =
  'replace' on restart replays, implicit none / intent / deallocation, error
  messages.

---------------------------------------------------------------------------
## REVIEW 2: scheduler, parameter plumbing, build, legacy bit-identity

Findings file: docs/reviews/degrace-track-a/REVIEW-2-scheduler-params-findings.md

Scope: PMP2main.f90 (Nclass array, OutClass, SetExactSteps with the iord
permutation, WriteOutputSchedule call, legacy marks), PMP2init.f90 (optional
trailing Init.dat lines via iParseLineOpt/ParseLineOpt, validation, Setup.dat
block writer, CheckInit template, the print-format fix), PMP2mod_tools.f90
(new globals and the optional "extra outputs" block in ReadSetup),
PMP2mod_analyze.f90 (hook ordering), PMP2linker.f90 (hook), makefile.

Check at least:
- Legacy identity. With every new flag at its default (no "extra outputs"
  block or all zeros; all outputs class 1) and a positive `#outputs`, is the
  simulation bit-identical to cz? Trace OutClass and Nclass defaults, the new
  `If(iBDM==1 .and. OutClass==1)` guard, whether DENSIT is now called where
  it was not before (mDENSIT semantics), the Power hook condition, and
  PMP2init's console/Setup.dat output. Also with `#outputs < 0` and no
  classes.
- Exact-z scheduler. Are classes attached to the correct redshift after the
  sort/permutation? Edge cases: duplicate z, z >= z_init, class count /=
  output count, more than 1000 outputs, inserted versus moved ladder steps,
  the NlastX stop index. Is `Nclass(NstepM)` reconstructed deterministically
  on a restart from a checkpoint, or does a restart lose it? What happens to
  LensPlanes/extras files when a restart replays an output step?
- PMP2init parsing. Behaviour of iParseLineOpt/ParseLineOpt at end of file,
  on a present but malformed line, and on a value line whose trailing comment
  contains `=` (ParseLine keys on the LAST `=`). Does the Setup.dat block that
  PMP2init writes contain any `=` inside comments that ReadSetup would then
  mis-parse? Are validation messages and stops correct? Is the format
  `'(a,/,(10f10.5))'` right for 1..1000 values?
- ReadSetup. Detection of the optional block via `index(Line, 'extra
  outputs')`, iostat handling, defaults when the block is absent (Setup.dat
  files of campaigns running now must still be read), partial blocks, stale
  trailing lines, kinds and defaults of the new globals (KcDelta = 1.5, KcMom
  = 1.0, NdepMom, SubM, FRhalo, LensSlab, LensGrid, Nclasses, zclass(1000)).
- makefile. `CC = icx`, CFLAGS, the object order for PMP2main, PMP2MG,
  PMP2spectrum and PMP2start (Fortran module dependency order: Tools -> fft5
  -> Density -> ExtraOutputs -> Subsample -> Power -> ...), explicit dependency
  lines (PMP2linker.o on PMP2mod_extras.o; what about PMP2mod_analyze.o on the
  subsample and extras modules, PMP2mod_power.o on extras?), and whether a
  clean `make` with the documented toolchain (module purge; module load
  intel_comp/2024.2.0; export I_MPI_F90=ifx; module load compiler-rt tbb
  compiler mpi; -march=core-avx2) would succeed. Reason from the makefile;
  do not build.
- Did anything on the branch change the content or file names of the legacy
  BDM catalogues or snapshots for users who do not enable the new flags?

---------------------------------------------------------------------------
## REVIEW 3: hashed quantised particle subsample (I9)

Findings file: docs/reviews/degrace-track-a/REVIEW-3-subsample-findings.md

Scope: PMP2mod_subsample.f90, glam_hash.c, the glam_hash.o rules in the
makefile, the WriteSubsample call in PMP2mod_analyze.f90, and the tests
test/extras/test_subsample.f90, test/extras/hash_ref.c,
test/extras/check_subsample.py, test/extras/run_subsample_test.sh, plus
read_subsample / the V6 check in test/precise_z/asserts.py. Spec: plan
section 5, item I9 (splitmix64 hash of the particle ID with SALT
0x5EED00DE6C0A2026, keep 1 in M, int16 quantisation of positions and
velocities, deterministic across gravity models and thread counts).

Check at least:
- What "ID" is hashed? GLAM stores particles in plain arrays. If the array
  index is used, verify that ic2pm.f90 preserves the Gadget ID order and that
  PMP2main never reorders particles (MOVE, checkpoint read/write, iDumpPart,
  parallel regions), so the same physical particle is selected in f(R)/LCDM
  twins, across restarts, and across thread counts. Also check the legacy
  PMP2start path if it is reachable.
- C interoperability: bind(C) interfaces (c_int64_t / c_int32_t, value
  attributes, unsigned arithmetic emulated in Fortran?), the splitmix64
  implementation against the published constants (0x9E3779B97F4A7C15,
  0xBF58476D1CE4E5B9, 0x94D049BB133111EB, shifts 30/27/31), and whether
  glam_hash.c, hash_ref.c and the Python reference in check_subsample.py are
  independent implementations or copies of one idea. Selection rule: modulo
  bias versus threshold, signedness of the comparison.
- Integer kinds and overflow: 1024^3 = 2^30 particles fit int32, 2048^3 =
  2^33 do not. Check every counter, index, allocation size, loop variable,
  meta count and OpenMP reduction, including the existing GLAM variables the
  module relies on (Nparticles, loop bounds).
- Quantisation: the position step (Box/65536?), rounding, the periodic edge,
  the velocity scale and clipping bound, the logged n_clip, the exact
  decoding formulas in the .meta and whether check_subsample.py and
  asserts.read_subsample invert them exactly; endianness and the stream
  layout (array order of qx,qy,qz,vx,vy,vz).
- Determinism and order: is the output order fixed by particle index under
  OpenMP, and does a test actually check order? How is n_sel bounded and the
  output allocated (two passes, or an upper bound)?
- Velocity epoch: are subsample velocities at a_v = a - da/2 like the
  momentum products, converted with the same factor, and is a_v recorded?
- Test quality: does run_subsample_test.sh test 1e7 IDs as claimed, and
  could it pass with a broken writer (does it compare to an independent
  reference, or to the code's own output)? Cleanup of its work files.

---------------------------------------------------------------------------
## REVIEW 4: delta_L dump in the FML IC generator (I7, route 1)

Findings file: docs/reviews/degrace-track-a/REVIEW-4-deltaL-generator-findings.md

Scope: 2LPTIC_Gui/fml_intel_port/DeltaLDump.h (new), the Main.cpp hunks
(`git diff cz..degrace-track-a -- 2LPTIC_Gui/fml_intel_port/Main.cpp`: an
include, three Lua parameters, the dump block, and a separate pre-existing
ID-rounding hunk floor -> llround), APPLY.md, input_smoketest.lua,
input_production_L1024_Np2048.lua, Makefile.intel; the tests
test/extras/check_deltaL.py, test/extras/check_glam_deltaL.py,
test/extras/run_deltaL_test.sh; and the "delta_L dump" section of
2LPTIC_Gui/README.md. FML sources on disk (gitignored):
2LPTIC_Gui/FML/FML/FFTWGrid/FFTWGrid.h, FML/RandomFields/GaussianRandomField.h
(or RandomFields/*.h), FML/Global/Global.h, FML/MPIParticles/MPIParticles.h,
and the built copy FML/FML/LPT/example/Main.cpp (should equal the port copy).

Check at least:
- In FFTWGrid.h: the normalisation of the real-to-complex transform (is the
  Fourier field divided by N^3, so that delta_k = N^-3 sum_x delta(x)
  exp(-i k.x)?), the local complex index `(ix_local*Nmesh + iy)*(Nmesh/2+1)
  + iz`, get_local_x_start / get_local_nx, any extra slices in the Fourier
  grid, the half-complex axis, and the integer wavevector mapping `ix <= N/2
  ? ix : ix - N`. In the random-field generator: how the Gaussian field is
  drawn (power per mode, box-volume factors, cell-centre shift, any filter),
  whether the field is already at z_ini with the input P(k) and no growth
  factor, and whether fix_amplitude / phase reversal act before the dump.
- MPI: all ranks call MPI_Reduce with identical counts; chunking by 1<<28
  floats and the int cast; MPI_IN_PLACE on task 0; ranks holding no cube
  modes; the n > Nmesh/2-1 clamp; memory per rank at production (n = 244,
  192 ranks).
- Growth factor: growth_raw = sqrt(Om/a^3 + OL) * int e^{2.5x} / (Om + OL
  e^{3x})^{1.5} dx. Show that this is proportional to D(a) = H(a) int_0^a
  da' / (a' H(a'))^3 and compute independently whether 4000 Simpson intervals
  from ln(1e-8) give D(0.02) to better than 1e-6 relative. Note what is
  omitted (radiation, neutrinos).
- Meta correctness: the lattice statement and ID formula versus Main.cpp
  (the ID hunk uses llround(q*Np)); k_c versus n; the `N_dep = N_mesh` key;
  the definition of rms_deltaL_cube; the axis_map claim (does the FML x axis,
  the slab-decomposed one, become Gadget pos[0]?).
- The ID rounding hunk: is llround correct for lattice coordinates q in
  [0,1) at every Np (edge q*Np = Np-0.5 cannot occur? prove), and does it
  ever change an ID where floor was right? Is it documented truthfully in
  APPLY.md?
- The Python checks: derive independently the sign in delta(k) = -i k.Psi(k)
  for the stated convention and the half-cell deposit phase factor; judge
  whether the assertions would catch a transposed axis, a conjugation, a
  cell shift and a wrong D (and state why); check the Gadget-1 header offsets
  used (72, 96, 124-136, 136-160); check run_deltaL_test.sh for failure
  propagation (set -u without -e; `fail` set inside subshells or pipelines),
  cleanup and the paths it assumes.
- Production: time and memory of the gather, the 0.47 GB file size claim
  ((2n+1)^2 (n+1) * 8 B), and whether `write_deltaL = true` in the production
  Lua could break the existing IC jobs (file name collisions, disk).

---------------------------------------------------------------------------
## REVIEW 5: Python tools and the test suites

Findings file: docs/reviews/degrace-track-a/REVIEW-5-tools-tests-findings.md

Scope: DEGRACE-companion/gen_inputs.py, check_run.py, compare_ecosmog_pk.py,
README.md; test/precise_z/run_tests.sh and asserts.py (`git diff
cz..degrace-track-a -- test/precise_z`), test/precise_z/README.md;
test/extras/test_modes.f90, check_modes.py, run_modes_test.sh,
check_frhalo.py. Reference: docs/degrace_ecosmog/README.md,
docs/degrace_ecosmog/aout.json and nodes.csv (the pinned ECOSMOG metadata),
AGENTS.md (provenance rules), PMP2init.f90 (how PkTable.dat's header is
turned into Om, OmL, h, sigma8 and how `Pk tune`/`sig8` renormalise).

Check at least:
- gen_inputs.py: the PkTable header (ombh2, omch2 = Om h^2 - ombh2 so that
  PMP2init's h = sqrt((ombh2+omch2)/Om) is right), who normalises the table
  to sigma8 and to which sigma8 (`sigma8_ic` in nodes.csv), the input P(k)
  file and units (k/h, (Mpc/h)^3), the z list derived from aout.json
  (ordering, rounding to 7 decimals, the resulting float32 epoch error),
  class assignment (primary = multiples of 0.125), every generated Init.dat
  line against GLAM's parser (ParseLine keys on the LAST `=`: no `=` in
  comments; fixed widths), step da, z_init, Nrow/Ngrid, BDM.config, Run1 /
  CATALOGS creation, the f(R) parameters (fR0 sign and magnitude for
  log10|fR0| = -5.25 -> -5.6234e-6; fR_n), the provenance JSON (generator,
  IC method and epoch, seed) versus AGENTS.md, and `--check`.
- check_run.py and compare_ecosmog_pk.py: the ECOSMOG assumptions (S01 = IC,
  S02..S43 in ascending a; POWMES columns col1 = integer k index, col4 =
  P/L^3) against docs/degrace_ecosmog/README.md; the CIC deconvolution 1 -
  (2/3) sin^2(pi k / 2 k_Ny) (is this the shell-averaged CIC window? cite);
  PowerDM.R parsing; the 3-ulp tolerance; sacct parsing; failure modes on
  partial runs.
- asserts.py / run_tests.sh: do cases 9-12 and the L3 block test what
  test/precise_z/README.md claims? False-pass risks: missing files that
  silently pass, exceptions swallowed, tolerances that cannot fail (judge
  L3's 3 sigma with N_eff = 3 x ring modes / 2, V2's 1e-4, V6's 1e-3), the
  mirrored constants ZINIT / DA0 / NSTEPS / redshift lists between
  run_tests.sh and asserts.py, hard-coded 150.0 versus BOX, the `np.rint`
  mode binning versus GLAM's PowerDM binning, the big-endian ProjectDENSIT
  record offsets, cube_from_particles (CIC on NGRID) and binned_power.
- test/extras: what test_modes.f90 + check_modes.py assert (is the V1 layout
  test specific enough to pin the sign of the imaginary part and the ±m
  symmetry?), check_frhalo.py, the link lines and cleanup in the run scripts.
- Portability and hygiene: absolute paths, environment assumptions, exit
  codes, leftover work directories.

---------------------------------------------------------------------------
## REVIEW 6: documentation, provenance and repository hygiene

Findings file: docs/reviews/degrace-track-a/REVIEW-6-docs-provenance-findings.md

Scope: docs/degrace_companion_implementation_plan.md (`git diff
cz..degrace-track-a -- docs/`), DEGRACE-companion/README.md,
2LPTIC_Gui/README.md, 2LPTIC_Gui/fml_intel_port/APPLY.md,
test/precise_z/README.md, .gitignore and DEGRACE-companion/.gitignore, with
AGENTS.md as the rule book. Read the code only as far as needed to test a
documentation claim.

Check at least:
- Do the plan's sections 5-9 and 16 describe what the code does? Spot-check
  the file formats and .meta keys of LensPlanes, deltam / mom_*, Subsample,
  fR_halo and deltaL against the writers (PMP2mod_extras.f90,
  PMP2mod_subsample.f90, 2LPTIC_Gui/fml_intel_port/DeltaLDump.h); the
  defaults (k_c 1.5 / 1.0, N_dep = Ngrid/2, M); the storage numbers (201 MB
  planes, 0.47 GB deltaL, subsample sizes); the velocity-epoch statements
  (a_v = a - da/2); the claim that the code is bit-identical to cz with the
  flags off.
- Are the section-14 checklists consistent with section 16 and with git
  (`git show --stat <hash>` for every hash cited: does it exist and touch
  what is claimed)?
- AGENTS.md compliance: IC method, generator version, LPT order, z_init, seed
  and `ic2pm` velocity epoch recorded per campaign (gen_inputs.py provenance,
  the pilot submit scripts under DEGRACE-companion/runs/, which are
  gitignored but on disk); no silent substitution of PMP2start; cleanup and
  file-count guidance (are test work directories ignored, are generated files
  left in the tree?).
- README accuracy versus the actual CLIs (argument names and defaults in
  gen_inputs.py, check_run.py, compare_ecosmog_pk.py; the run scripts'
  options and runtimes).
- Open items listed versus actually open (L3, V9, I7 routes, pilot status),
  and anything a fresh session would need that is missing (rebuild steps,
  how to run each test, where the pilots and their logs live).
