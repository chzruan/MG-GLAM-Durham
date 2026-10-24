# DEGRACE companion runs

MG-GLAM companions of Baojiu Li's DEGRACE ECOSMOG suite (same ICs, same
cosmology, same 42 output epochs) and of new cosmologies. Plan:
`docs/degrace_companion_implementation_plan.md`; suite reference:
`docs/degrace_ecosmog/README.md`.

## Inputs

```bash
# S2 check: the epoch generator reproduces every ECOSMOG list and the fiducial list
micromamba run -n cosemu python3 gen_inputs.py --check
# companion of node 2 (design centre), LCDM, Ng 2048 on Baojiu's 1024^3 ICs,
# BDM at the 17 primary outputs only, lens planes 4 x 2048^2
micromamba run -n cosemu python3 gen_inputs.py --node 2 --gravity lcdm \
    --classes primary --lens-slabs 4 --lens-grid 2048 --out runs/Node_002_LCDM_L1024Np1024Ng2048
```

`gen_inputs.py` writes `PkTable.dat` (5-line cosmology header + Baojiu's
CAMB z=0 spectrum), `Init.dat` (`#outputs = -42` exact mode, `Pk tune = 1`,
`iSave = 0`, the optional trailing block with the output classes and the
lens-plane geometry) and `inputs_provenance.json`. Defaults: z_init 49,
da0 = 6e-4 (163 steps, the legacy late-time ladder), Nrow 1024, Ngrid 2048.
`--cosmology OMEGA_M H SIGMA8 --pk FILE [--fR0 ...]` serves new cosmologies
(stations computed from χ(a), matter + Λ only). `--redshifts Z ...` replaces
the 42 epochs by your own list of full outputs and `--stations Z ...` adds
lens-plane-and-P(k)-only outputs; `--run N` writes `Run<N>/` (default 1).
[`examples/twin_L128Np256Ng512/`](../examples/twin_L128Np256Ng512/README.md)
runs a small LCDM / f(R) twin with these options from ICs to packed outputs.

Extra on-the-fly products (plan items I8–I10) are switched on with
`--extra` (low-k density modes `deltam` from the P(k) transform and
momentum modes `mom_x/y/z` on an `N_dep` grid at class-1 outputs;
`--kc-delta 1.5 --kc-mom 1.0 --ndep-mom 0`), `--subsample-m 100` (1 %
hashed quantised particle subsample, once I9 is integrated) and
`--fr-halos` (f_R at the BDM halo centres, f(R) runs). Every product has a
`.meta` sidecar; see the plan's section 3 for the conventions.

## Running

A run directory holds `Init.dat`, `PkTable.dat`, `bin/` (frozen copies of
`PMP2init.exe`, `ic2pm.exe`, `PMP2main.exe` with md5 and git hash) and
`Run1/` (`BDM.config`, outputs). `submit.sh` runs `PMP2init` (once), then
in `Run<box>/`: `ic2pm.exe <ics.> 1.0 half <box>` (the fourth argument, default 1,
is the realisation number stamped on every product name) and
`PMP2main.exe <<< 300`, and
deletes the checkpoint files at the end. `PMP2init` exits with status 1
when it refuses a configuration (`error stop`; a plain `stop`, status 0,
before 2026-10-18), so `submit.sh` stops before `ic2pm`. Outputs per analysis step:
`PowerDM.*`, `LensPlanes.<step>.<box>.{bin,meta}`, `CATALOGS/CatshortV.*`
(class-1 outputs only; `.DAT` text and/or `.bin` + `.meta` binary, chosen
by `CatalogueFormat` in `Run1/BDM.config`, which `gen_inputs.py` writes
only when `--catalogue-format` is given; the compiled default is `binary`
since 2026-10-18, so a run without the key writes `.bin` + `.meta` only;
pass `both` or `text` while a consumer still needs the `.DAT`), and
`OutputSchedule.dat` (step, a, z, χ, exact, class) written at start-up.
`PMP2main` validates `BDM.config` at start-up. Read the catalogues with
`python/glamio` (`glamio.read_catalogue(stem)`); `check_run.py` verifies
text and/or binary catalogues and flags stale `CATALOGS/.Catshort*.tmp.*`
files left by a killed finder.

`runs/` is git-ignored: record each campaign's inputs, job IDs and results
in the plan's section 16.

## Verifying a run

```sh
micromamba run -n cosemu python3 DEGRACE-companion/check_run.py \
    DEGRACE-companion/runs/Node_002_LCDM_L1024Np1024Ng2048 --node Node_002 \
    --aout docs/degrace_ecosmog/aout.json [--planes all] [--job JOBID]
```

Checks the 42-output schedule (exact flags, 17 + 25 classes, epochs equal to
float32(1/(1+z)) of `Init.dat` and to the pinned ECOSMOG list), one
`PowerDM.R` and one `LensPlanes` pair per output (size, epoch, axis means
= 1 on a sample or all files), BDM catalogues at the class-1 outputs, the
`ic2pm` velocity epoch, the `timing.log` breakdown (median step, output
overhead) and, when the job is finished, sacct elapsed / MaxRSS / CPU
efficiency. Exit 0 = all checks passed. The output classes are compared with
`Init.dat`, so any redshift list works (the primary-redshift check applies to
the 42-output design only), and `--box` defaults to the run's only `Run<N>/`.
Packed runs are checked the same way
(the files are read through glamio), plus their `glamio.pack status`
(`packed` or `pruned` passes; an unfinished pack or prune warns).

## Packing a finished run

After the job has finished and `check_run.py` has passed, pack each box
(about 400-430 files) into one HDF5 file per output step: `Output.<step>.<box>.h5`,
plus `RunInfo.<box>.h5` and `packed.md5`, 44 files for 42 outputs. The bytes
are unchanged, and `glamio` reads the packed files with the same calls
(`python/README.md`, "Packed runs").

```sh
cd <run dir>     # the directory holding Init.dat and Run<box>/
sbatch -J pack_<name> <repo>/DEGRACE-companion/pack_run.sh "$PWD/Run0"                        # pack, byte check, delete the originals
sbatch -J pack_<name> --export=ALL,PRUNE=no <repo>/DEGRACE-companion/pack_run.sh "$PWD/Run0"  # pack and byte check only
```

`pack_run.sh` packs in place, rebuilds every original from the packed files
in scratch, and requires an identical file list and `cmp` equality of every
file (the rebuilt copy is removed afterwards, also on failure). It then
deletes the originals (unless `PRUNE=no`), through `glamio.pack prune`, which
refuses while any original is in no packed file and re-verifies every md5
first; this is the default since 2026-10-24, to cut the file count. Any
failure stops the job before
anything is deleted. A resubmission resumes from `python -m glamio.pack
status Run0`: a finished pack is re-checked against the originals, an
interrupted prune is finished (after a byte check of the originals still
present; not with `PRUNE=no`), and a box that is `busy` (a
`.partial` exists) or `inconsistent` is left alone. Run `check_run.py` (it
also reads packed runs) and any script that opens native files by name before
pruning; afterwards
`python -m glamio.pack unpack Run0 --out DIR` rebuilds the native files byte
for byte (75 GB for an L1024 box: run it under Slurm). For an L1024 Np2048
box (75 GB) on one core, pack takes 20 min, the check 15-30 min and the
prune 30 min. Copy a packed box with its run directory (`Init.dat`, `bin/`,
provenance) and check it with `md5sum -c packed.md5`.

## Comparing with the ECOSMOG twin

```sh
micromamba run -n cosemu python3 DEGRACE-companion/compare_ecosmog_pk.py \
    DEGRACE-companion/runs/Node_002_FR_L1024Np1024Ng2048 002 [--kmax 1.0]
```

Prints P_GLAM / P_ECOSMOG per output in the bands 0.01-0.03, 0.03-0.1, 0.1-0.3
and 0.3-1 h/Mpc, matching our `PowerDM.R` epochs to the ECOSMOG snapshots
S02..S43 (S01 = IC) through the pinned `aout.json`. ECOSMOG only has f(R)
spectra (POWMES, `/cosma7/data/dp004/bl267/Runs/DEGRACE/Runs/L1024/FR/power_spectra/`),
so a ΛCDM companion compared this way shows the f(R) enhancement at late
times; the f(R) companion is the real twin test. Our raw CIC spectra are
deconvolved with the shell-averaged window 1 − (2/3) sin²(πk/2k_Ny); at
k ≳ 0.5 and z ≳ 3 the comparison is limited by discreteness (P ~ a few
(Mpc/h)³), so judge the bands below 0.3 h/Mpc there. Both runs use the same
IC files, so the low-k ratio is mode-exact (no cosmic variance).

## Known conventions and offsets

- Velocities in the new products (`mom_*`, `Subsample`) are at a_v = a − da/2
  and the `.meta` records `a_v`. The legacy BDM catalogue (`CatshortV`)
  converts with 1/a instead: its velocities (and V_rms, the virial ratio and
  λ) are low by a_v/a (0.8 % at z = 2, the first catalogue output; 0.3 % at
  z = 0). Multiply by a/a_v when combining them (the binary catalogue's
  `.meta` records `a`, `a_v` and `da`;
  `glamio.read_catalogue(stem).velocities_at_av()` does it), or wait for
  the decision recorded in the plan (section 11). (Direction corrected
  2026-10-18; the earlier text said "high by a/a_v".)
- `kc_mom` must stay below k_Nyquist(N_dep_mom)/4 (0.785 h/Mpc at N_dep =
  1024, L = 1024; default 0.75): the momentum deposit is not interlaced and
  `PMP2init` refuses larger values. The node-2 f(R) pilot predates this rule
  (k_c,mom = 1.0): use its `mom_*` below k = 0.78 h/Mpc.
- Particle IDs: `ic2pm` places each particle at array index = Gadget ID and
  tags the lattice order in the PM header (`IDz` for Baojiu's 2LPTic files,
  `IDx` for our FML generator); see the `lattice` key of `Subsample.*.meta`.

Note (2026-10-21, MG-GLAM branch `memory-cuts-20261021`): the finder's memory admission limit
`MaxMemory` now defaults to the job's memory (`--mem`, via SLURM_MEM_PER_NODE) or the node's MemTotal
instead of 500 GiB, so the pilots' `--mem=200G/400G` become their limits (their counted memory is far
below). `gen_inputs.py --bdm-max-memory GIB` writes an explicit value for builds between 2026-09-06
and 2026-10-20, which default to 500 GiB.
