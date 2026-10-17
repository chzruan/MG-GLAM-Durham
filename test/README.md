# MG-GLAM symmetron example (MG_model = 3)

A small, quick-to-run example configured for the **symmetron** modified-gravity
model. It exercises the full MG path (scalar-field multigrid solve + fifth
force) at low resolution so it finishes fast.

## Configuration
- **Cosmology:** box = 750 Mpc/h, 128³ particles, 256³ force mesh, z_init = 100 → z = 0
  (Ωm = 0.3089, ΩΛ = 0.6911, h = 0.6774, σ8 = 0.8159).
- **Gravity:** symmetron — `MG flag = 1`, `MG_model = 3`, with MG-GLAM's default
  symmetron parameters:
  - `a_*  = 0.5`   — symmetry-breaking scale factor (⇒ z_SSB = 1)
  - `xi   = 1e-3`  — range / Compton-wavelength parameter
  - `beta_* = 0.1` — matter coupling (enters the EOM as ∝ β² — set ~1 for a
    stronger, more visible fifth force)

  Edit these in `Setup.dat` (and `Init.dat`, see below) to change the model.

## Files
| file | role |
|------|------|
| `Setup.dat`      | **Parameter file the executables read at runtime — authoritative.** |
| `Init.dat`       | Human-readable input; `PMP2init.exe` reads it to (re)generate `Setup.dat`. |
| `PkTable.dat`, `lcdm.dat` | z = 100 linear P(k) used to seed the initial conditions. |
| `TableSeeds.dat` | RNG seed table; the realization index selects a seed (reproducible ICs). |
| `BDM.config`     | BDM halo-finder config (used only if `Find BDM halos = 1`). |
| `arun.sh`        | SLURM submit script (cosma8). |

`Setup.dat` and `Init.dat` are kept consistent — both already set
`MG flag = 1`, `MG_model = 3`. The executables only read `Setup.dat`; `Init.dat`
matters only if you regenerate `Setup.dat` with `PMP2init.exe`.

## Build (from the repo root)
```bash
module purge
module load intel_comp/2024.2.0
export I_MPI_F90=ifx
module load compiler-rt tbb compiler mpi

make PMP2start PMP2MG        # IC generator + modified-gravity evolution binary
# make PMP2init             # optional: only if regenerating Setup.dat from Init.dat
```
`PMP2MG.exe` is the modified-gravity binary. (`PMP2main.exe` is built from the
same objects and also honors `MG_flag` at runtime, but an MG run conventionally
uses `PMP2MG.exe`.)

## Run
```bash
cd test
bash arun.sh
```
For each `box` (realization index, default `2`) `arun.sh` creates `Run<box>/`
and submits a SLURM job that:
1. `PMP2start.exe <<< <box>` — generates initial conditions; the box index
   selects the RNG seed from `TableSeeds.dat`, so the ICs are reproducible.
2. `PMP2MG.exe <<< 157` — evolves **157** timesteps under symmetron gravity.

To run interactively instead of through SLURM, `cd Run<box>/` and execute the
two `../../PMP2*.exe <<< …` lines from `arun.sh` yourself (after the module
loads and `OMP_*` exports).

## Changing parameters
Edit `Setup.dat` directly, **or** edit `Init.dat` and run `PMP2init.exe` in this
directory (it reads `Init.dat`, writes `Setup.dat`). Keep the two in sync.
