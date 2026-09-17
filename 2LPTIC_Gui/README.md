# 2LPTIC (FML LPT) IC generator — Intel build for MG-GLAM

## Default IC workflow

**2LPTIC is the project default for new production, validation and convergence
simulations.** Use the corrected FML build described here with `lpt_order = 2`,
then convert its Gadget-format output with [`ic2pm`](../ic2pm.f90) before GLAM
evolution. Agent instructions are in the root [`AGENTS.md`](../AGENTS.md).
GLAM's native first-order `PMP2start` remains available for explicit legacy
reproduction or IC-method controls; existing campaigns retain their recorded
IC method. See [VALIDATION.md](VALIDATION.md) for the checked conventions and
end-to-end comparisons.

## `ic2pm` velocity epoch — erratum/update (2026-09-17)

Usage is now `ic2pm.exe <IC_basename> [S_vel] [half|sync]` (run from
`Run<box>/`, reads `../Setup.dat`). `<IC_basename>` includes the trailing
'.', and the PM files are written to the current directory. The epoch is
positional, so S_vel must come before it (`ic2pm.exe <IC> 1.0 sync`);
`ic2pm.exe <IC> sync` is rejected. Every check in `ic2pm` (usage error,
bad S_vel, missing `../Setup.dat` or IC file 0, header/Setup mismatch,
particle-count mismatch) prints a message and exits with status 1 (since
commit 6b629f1; at ee44100 only `die` aborts and a bad epoch did). Fortran
runtime I/O errors (a missing or truncated later IC file, an unreadable
Setup.dat value, an unwritable run directory) exit non-zero with the ifx
error code; a Setup.dat without the checkpoint block is rejected inside
Tools with a bare `stop` (exit 0), so job scripts should also check that
`PMcrd.DAT` exists. The third argument selects the epoch of the velocities
written to the PM files:

- `half` (**default**): the synchronous 2LPTic/Gadget velocities (at
  a_init) are rescaled to a_v = a_init − ASTEP/2 with PMP2start's
  growing-mode factor (a_v/a)^1.5 F(a_v)/F(a), F = sqrt(Om + OmL a^3),
  ASTEP being the PM-header step (= ASTEP0 from `Setup.dat`). This is what
  GLAM's kick-then-drift leapfrog expects. Do **not** edit the Gadget IC
  files themselves; keep them standard (synchronous). The shift depends on
  the GLAM run's ASTEP0, not on the IC (one IC serves runs with different
  da), and an edited file converted with the default `half` would be
  shifted twice.
- `sync`: no shift; particle files (PMcrs*.DAT) byte-identical to the
  pre-fix converter (commit 00ea3df), PMcrd.DAT identical apart from the
  uninitialised AEXP0 header field. Only for reproducing the runs made
  before the fix (`conv_da*`, `fid2LPTIC_*`, `ic2pm_val_L1024`). Their
  submit scripts (`conv_da{4,8,16}/Run1/submit_conv.sh`,
  `ic2pm_val_L1024/Run1/submit_val.sh`,
  `fid2LPTIC_L512Np2048Ng4096{,_da4,_da6}/Run1/submit.sh`) call
  `ic2pm.exe <IC> 1.0` without a third argument and would now get `half`;
  add `sync` explicitly to reproduce them. For files written by builds of
  commit 6b629f1 or later the PM header string records the mode:
  `ic2pm: 2LPTic ingest, v at a-da/2` (half) vs `ic2pm: 2LPTic ingest`
  (sync). Every file written by an older build carries the sync string
  whatever its mode (including the half-mode A/B runs of jobs 12006891-93,
  made with an ee44100 build); for those the mode is only in the
  `velocity epoch =` line of the ic2pm log.

Without the shift the late-time P(k) is high by ≈ 0.8 × 0.75 da/a_init at
linear scales (+1.2 / +1.8 / +2.4 / +4.8% for da0 = 4e-4 / 6e-4 / 8e-4 /
1.6e-3 from z = 49; measured at k ≤ 0.05 h/Mpc: +1.20 / – / +2.39 / +4.78%
at the z ≈ 2-3 peak, +1.18 / – / +2.36 / +4.72% at z = 0)
and by ≈1.4× that at 0.3<k<1 by z=0. This is what several earlier numbers
in [VALIDATION.md](VALIDATION.md) measured; see its 2026-09-17 erratum and
`../halfstep_ab/README.md` for the A/B test. `S_vel` values: 1.0 for ICs
from this build and for the original DEGRACE `ics.*`; `5.168609e6`
(= 5.12e6/0.99059529, but write the number: the parser rejects `/`) for
Gui's old HEFT files (5.12e6 = 100·box/a: see "Changes vs Gui's
Main_2LPT_lua.cpp" below; 0.99059529: VALIDATION.md, High-res validation).
Fix commit: ee44100 on branch `ic2pm-halfstep` (not yet merged into `cz`);
`ic2pm` binaries built before it have no epoch argument and always write
synchronous velocities.

## Generator identity

Gui Brando's "2LPTic" IC generator (email thread: `Gmail - Fw_ 2LPTic.html`)
is Hans Winther's **FML** library LPT example (C++ + Lua parameter file),
not Scoccimarro's classic 2LPTic. It uses the same GSL `ranlxd1` random
generator / N-GenIC-style seed table as NGenIC/2LPTic/monofonIC, so the same
seed (2026) gives the same large-scale phases as the existing DEGRACE
1024^3 2LPTic ICs.

## Layout

- `Gmail - Fw_ 2LPTic.html` — the email thread introducing the code (kept out of version control)
- `FML_LPT_Lua_Setup_Guide.txt`, `Main_2LPT_lua.cpp`, `Makefile`,
  `input_file.lua`, `setup_cosma.sh`, `run_cosma.sh` — Gui's original
  attachments (unmodified, for reference)
- `FML/` — our clone of https://github.com/HAWinther/FML (shallow, cloned
  2026-07-28) with the changes listed below
- `FML/FML/LPT/example/` — the working build directory:
  executable `test`, `Makefile` (Intel), `setup_env.sh`, `input.lua`
  (128^3 smoke test), `input_production_L1024_Np2048.lua`, `run_cosma.sh`,
  `pofk_bli_z49.txt` (the DEGRACE linear P(k) at z=49, copied from Gui's dir)

## Changes vs upstream FML (clone @ 2026-07-28)

1. `FML/ParticleTypes/SimpleParticle.h` — added Gui's `HEFTParticle`
   struct (copied verbatim from `/cosma8/data/dp203/dc-bran2/new_FML`).
2. `FML/MPIParticles/MPIParticles.h` — Gui's version (hardcodes ndim=3 in
   `create_particle_grid`; needed because `HEFTParticle`'s default ctor is
   not constexpr, which breaks `constexpr GetNDIM(T())`).

All other FML files are upstream. Gui's remaining tree changes are additive
(white-noise dump helpers etc.) and not needed. Crucially his
`RandomFields/GaussianRandomField.h` is purely additive, so upstream
generates the identical realization.

## Changes vs Gui's Main_2LPT_lua.cpp (two bug fixes)

`FML/FML/LPT/example/Main.cpp` = `Main_2LPT_lua.cpp` + two fixes:

1. **All particle IDs were written as 1.** IDs are computed from the
   Lagrangian position `q`, but `create_particle_grid()` only fills `pos`,
   never `q`, so `floor(q*Np)=0` for every particle. Fix: copy pos -> q
   right after `create_particle_grid` (before displacement).
2. **Velocities were a factor 100*box/a (~5e6 for L1024, z=49) too small.**
   Gui passed `vel_norm = 1/sqrt(a)`; the Gadget convention (cf. FML's own
   COLASolver `output_gadget`) is `vel_norm = 100*box/a^1.5` for internal
   velocities `a^2 (H/H0) f Psi` in box units.

**Both bugs are present in Gui's generated ICs** (verified directly on his
512^3 test output in `/cosma8/data/dp203/dc-bran2/new_FML/FML/LPT/example/snap/`:
all IDs = 1, max|u| = 6e-4 km/s), including the 2048^3 production IC.
Positions/phases in his files are fine; velocities and IDs are not.

## Build (Intel toolchain, same as MG-GLAM)

```sh
cd FML/FML/LPT/example
source setup_env.sh     # module purge; module load intel_comp/2024.2.0 compiler-rt tbb compiler mpi
make clean && make -j 8
```

Compiler `mpiicpx` (icpx + Intel MPI 2021.13), `-march=core-avx2`
(cosma8 is AMD Zen2 — never use `-x*` flags), `-DGADGET_LONG_INT_IDS`
(64-bit IDs, required for 2048^3 > 2^32 particles).
Libraries: FFTW `/cosma/local/fftw/intel_2024.2.0_intel_mpi_2024.2.0/3.3.10-epyc`
(MPI ABI matches Intel MPI), GSL 2.7.1, Gui's public static Lua 5.4
(`/cosma/apps/dp203/dc-bran2/lua-5.4.0/install`).

## Run

The code always reads `input.lua` from the cwd (the argv parameter is
ignored). Small interactive test (~2 s):

```sh
source setup_env.sh
export LD_LIBRARY_PATH=/cosma/local/fftw/intel_2024.2.0_intel_mpi_2024.2.0/3.3.10-epyc/lib:/cosma/local/gsl/2.7.1/lib:$LD_LIBRARY_PATH
mpirun -np 4 ./test input.lua
```

Production 2048^3: `cp input_production_L1024_Np2048.lua input.lua`, then
`sbatch run_cosma.sh` (3 cosma8 nodes; see memory estimate in the script).

Output: Gadget-1 format, one file per MPI rank (`<prefix>.<rank>`), float
pos [Mpc/h] / float vel [Gadget u = v_pec/sqrt(a), km/s] / int64 IDs, plus
`pofks.txt` (k, P_lin, P_particle, r_cross, P_ini/P_input) as a built-in
sanity check.

## Validation performed (2026-07-28)

- 128^3 smoke test (4 ranks): IDs unique 1..128^3; v_pec rms(1d) = 50 km/s
  at z=49; P_ini/P_input = 1 within bin noise; r_cross ~ 1.
- Realization reproducibility: ran our Intel build with Gui's exact 512^3
  config (seed 2026) — `pofks.txt` matches Gui's gnu/OpenMPI output to
  every printed digit in all 256 bins (both P_ini and P_2LPT), i.e. the
  realization is independent of compiler, MPI library and rank count.
