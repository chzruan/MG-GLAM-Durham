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
(stations computed from χ(a), matter + Λ only).

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
in `Run1/`: `ic2pm.exe <ics.> 1.0 half` and `PMP2main.exe <<< 300`, and
deletes the checkpoint files at the end. Outputs per analysis step:
`PowerDM.*`, `LensPlanes.<step>.<box>.{bin,meta}`, `CATALOGS/CatshortV.*`
(class-1 outputs only), and `OutputSchedule.dat` (step, a, z, χ, exact,
class) written at start-up.

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
efficiency. Exit 0 = all checks passed.
