# Unit and integration tests of the DEGRACE companion extras

All scripts load the Intel modules themselves (`SKIP_MODULES=1` if already
loaded), run on a login node within the light tier (≤ 4 threads or 4 MPI
ranks, seconds to ~1.5 min), write only under `test/extras/work/`
(gitignored; delete afterwards) and exit 0 when every check passes.

| script | what | needs | time |
|---|---|---|---|
| `run_modes_test.sh` | V1 layout of `DumpModes`/`FFT3Dforward` on synthetic cos/sin fields (`test_modes.f90`, `check_modes.py`); f_R trilinear interpolation (`check_frhalo.py`) | repo-root `.o` files from `make PMP2main` | seconds |
| `run_subsample_test.sh` | V5: splitmix64 hash on 1e7 IDs vs an independent numpy transcription, selection order, int16 quantisation and clipping (`test_subsample.f90`, `hash_ref.c`, `check_subsample.py`) | `glam_hash.o`, Intel compilers | ~1 min |
| `run_deltaL_test.sh` | I7: FML δ_L dump determinism (1/2/4 ranks), layout against the IC displacements (`check_deltaL.py`), end to end through `ic2pm` and `PMP2main` (`check_glam_deltaL.py`) | built FML generator `2LPTIC_Gui/FML/FML/LPT/example/test`, `PMP2init.exe`, `PMP2main.exe`, `ic2pm.exe`, `test/PkTable.dat`, `test/TableSeeds.dat` | ~75 s |
| `run_ic2pm_test.sh` | `ic2pm` places particles at index = Gadget ID, detects the ID lattice order (`IDz` 2LPTic-like int32 sequential, `IDx` FML-like int64 shuffled) and refuses a duplicate-ID set; verified through a one-step GLAM run with `Subsample M = 1` (`make_gadget_lattice.py`, `check_ic2pm_order.py`) | `PMP2init.exe`, `PMP2main.exe`, `ic2pm.exe` | seconds |

The scheduler, lens-plane, mode-dump, subsample (V2, V6, L1–L3) and f_R
smoke cases run in `../precise_z/run_tests.sh` (222 checks, ~1 min).
