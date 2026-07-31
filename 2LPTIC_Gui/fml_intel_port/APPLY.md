# FML Intel port — how to rebuild the 2LPTIC generator from scratch

The working FML clone (`2LPTIC_Gui/FML/`, 103 MB, its own git repo) is kept
out of MG-GLAM version control. This directory holds every file we authored
or modified, so the build is reproducible from a fresh clone:

```sh
cd 2LPTIC_Gui
git clone https://github.com/HAWinther/FML.git      # cloned 2026-07-28
cp fml_intel_port/SimpleParticle.h  FML/FML/ParticleTypes/SimpleParticle.h
cp fml_intel_port/MPIParticles.h    FML/FML/MPIParticles/MPIParticles.h
cp fml_intel_port/Main.cpp          FML/FML/LPT/example/Main.cpp
cp fml_intel_port/Makefile.intel    FML/FML/LPT/example/Makefile
cp fml_intel_port/setup_env.sh      FML/FML/LPT/example/
cd FML/FML/LPT/example
source setup_env.sh && make clean && make -j 8      # -> ./test
```

File provenance:

- `SimpleParticle.h` — upstream + Gui's `HEFTParticle` struct (verbatim from
  `/cosma8/data/dp203/dc-bran2/new_FML`).
- `MPIParticles.h` — Gui's version (ndim=3 hardcode in `create_particle_grid`,
  needed because `HEFTParticle`'s default ctor is not constexpr).
- `Main.cpp` — Gui's emailed `Main_2LPT_lua.cpp` + our two bug fixes
  (q <- pos so IDs are correct; Gadget `vel_norm = 100*box/a^1.5`).
- `Makefile.intel` — Intel oneAPI toolchain (mpiicpx, -march=core-avx2,
  Intel-MPI FFTW, GSL 2.7.1, Gui's static Lua). Copy it to `Makefile`
  (that name is gitignored repo-wide, hence the suffix here). NOTE: it
  contains the absolute FML_INCLUDE path — adjust if you clone elsewhere.
- `input_smoketest.lua` (128^3, 2 s on 4 ranks), `input_production_L1024_Np2048.lua`,
  `run_cosma.sh` — run configs. The binary always reads `input.lua` from cwd.

Upstream FML pin: if upstream has moved, any commit from 2026-07 works; the
only load-bearing upstream behaviour we rely on (N-GenIC-compatible seed
table in GaussianRandomField.h) is stable. Number of MPI ranks must divide
Nmesh.
