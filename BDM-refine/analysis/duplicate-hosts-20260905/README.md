# BDM duplicate-host repair and catalogue products

This branch diagnoses and repairs the equal-bound-mass duplicate failure in
`PMP2linker.f90`, and creates new catalogue products for the nominal z=0.25
fiducial and training sample. **All 740 catalogues are complete and validated**;
see [RESULTS.md](RESULTS.md) for final incidence, accounting and cleanup.
See [DIAGNOSIS.md](DIAGNOSIS.md) for the source
chain, membership limitations and the box 1 reproduction, and
[CAMPAIGN.md](CAMPAIGN.md) for scope, pilot resources and cost estimates.

The scientific rule is `strict-v1`: periodic distance <0.2 Mpc/h, equal stored
Mbound and Nparticles, relative speed <5 km/s, and |dlog10 Mtot|<0.005.
Connected components retain their lowest original row/candidate index.
This deliberately matches Freyja's strict criterion. It does not introduce
an additional Rvir-fraction cut or claim particle membership from particle
counts. In the finder the same comparisons use unrounded in-memory quantities;
rounding near thresholds can differ from a historical ASCII decision.

The finder separates the immutable unequal-mass host decision from the
numerical-duplicate graph, uses minimum-image geometry, and applies deletions
after decisions. The stable index is the original candidate-array index among
host-test survivors. Candidate generation and its separate FI race are outside
this repair. Full-run bitwise invariance across thread counts cannot be claimed.
Making the unequal-mass host pass read-only and periodic can also change which
pre-write candidates survive that pass relative to an old racy execution;
the patch does not establish full historical host-set equivalence. No survivor
measurement is recomputed by either removal pass.
The reference replay establishes the **catalogue-stage** mask and survivor
values; deleted snapshots prevent a historical particle-to-catalogue replay.

## Products and validation

New catalogues mirror the upstream paths under `products/strict-v1/`:

```
DESI_MGx100/data/GR/Run<ibox>/CATALOGS/CatshortV.0137.<ibox:04d>.DAT
mg_glam/DurMun_hmfemu_<gravity>_wide_sample_first_64_model<imodel>_L1024Np2048Ng4096/
    Run<ibox>/CATALOGS/CatshortV.0137.<ibox:04d>.DAT
```

Each `.DAT` has an adjacent `.cleaning.hdf5` file. A sidecar is published only
after its catalogue passes validation, and serves as the completion receipt.
It stores one cleaning list for that catalogue, without changing the original
row frame. The original eight header lines and retained data lines are copied
verbatim; a ninth `# BDM-cleaning ...` provenance line is inserted. Nhalo values
are not renumbered. Halocat's `numpy.loadtxt(..., skiprows=8)` remains compatible
because its default comment handling skips the added line. Consumers that
assume a numeric ninth line without comment handling must be adjusted.

Sidecar datasets:

| Dataset | Meaning |
|---|---|
| `drop_mask` | Boolean array with one entry per original data row; True means remove. |
| `exact_drop_mask` | Independent exact-defining-fields-only sensitivity bracket. |
| `drop_row_indices`, `dropped_Nhalo` | Zero-based original row indices and original output IDs of strict-rule removals. |
| `dropped_representative_row_index` | Original row index retained for each dropped row, aligned with `drop_row_indices`. |
| `groups_json` | UTF-8 JSON as a compressed uint8 array: all member row indices/Nhalo and retained representative, one record per component. |
| `validation_json` | Per-catalogue pair counts, remaining ambiguous pairs, mass bins, and box 1 velocity statistics. |
| `receipt_json` | Input/output SHA256, original row count, grid identity, rule/parameters, exact tool hashes and git commit, resources and survivor byte-identity checks. |

The source SHA256, rule and tool commit also appear in HDF5 attributes and
apply to every group in that sidecar. Membership validation is explicitly null.

`catalogue_validation.csv` contains one row per successfully validated
catalogue, including gravity, imodel, ibox, nominal redshift, snapnum, incidence,
paths and hashes. `mass_bin_validation.csv` contains the per-catalogue mass
tables: 0.2 dex bins from 12.4 to 14.4 and the requested partial final bin
14.4--14.5. `campaign_summary.json` reports completeness, failures and model
totals. Consult its expected/completed counts before consuming the whole grid.
Receipt `maxrss_kib` values are process high-water marks across a worker's
catalogues. Slurm batch-step MaxRSS and CPU efficiency are recorded separately
in the accounting provenance.

For box 1, the rule removes 6,191 of 3,385,000 full-catalogue rows. Restricted
to the mass-selected sample before grouping, it removes 5,967 of 1,882,717
rows (0.31694%), exactly reproducing Freyja. The full-source mask removes 5,969
above that threshold because two groups cross the mass boundary. Always use
the full original-row mask first when using these delivered products.
The archived Freyja HOD test, which uses no minimum mass and `log10Mtot<14.5`,
also removes exactly **6,191** box 1 rows under the strict rule (509 under the
exact rule). Its selected catalogue has 3,383,222 rows. Thus the rule and box 1
removals reconcile in the actual HOD selection as well as the >=12.4 incidence
audit; see `provenance/freyja-reconciliation.json` for all five HOD boxes.
The delivered full masks remove eight additional rows at masses >=14.5 across
boxes 3--5. Restricted to Freyja's HOD mass range, removed-row counts agree in
all five boxes: 6,191, 5,961, 5,975, 5,972 and 6,038.

The box 1 cleaner preserves all columns of all 3,378,809 surviving rows
byte-for-byte. The compiled Fortran predicate/component policy gives exactly
the same 6,191-row mask at stored precision; see `provenance/policy-replay.json`.
The 0.5--2 Mpc/h radial peculiar-velocity standard deviation rises from
186.950110 to 187.013296 km/s. This statistic uses all mass>=12.4 distinct
halo pairs and subtracts the pair mean; it is not the radial RMS.

The strict graph has zero surviving edges. The broader same-bound-mass/count
test still has 1,041 box 1 pairs within 0.2 Mpc/h. They are explicitly listed
as ambiguous; removing them would broaden the scientific rule and would not
be supported by the missing particle-membership evidence.

## Downstream use

Run Python with the required environment:

```bash
micromamba run -n cosemu python3 your_script.py
```

To consume a cleaned ASCII catalogue, use the `output` path in the validation
table or sidecar receipt. To apply its mask to the original ASCII rows:

```python
import sys
import numpy as np
sys.path.insert(0, '/cosma8/data/dp203/dc-ruan1/mgglam_claude/MG-GLAM/'
                   'BDM-refine/analysis/duplicate-hosts-20260905/tools')
from load_cleaning import load_cleaning

drop, receipt = load_cleaning('/absolute/path/to/CatshortV.0137.0001.cleaning.hdf5')
raw = np.loadtxt(receipt['source'], skiprows=8)
clean = raw[~drop]
original_row_of_clean_row = np.flatnonzero(~drop)
```

`load_cleaning` verifies the source SHA256 by default. A pre-existing HDF5
mirror can use the same mask only after its original-row alignment with that
ASCII file has been established; matching counts alone are insufficient.
The cleaner also accepts a flat HDF5 input directly and preserves dataset
dtypes, values and attributes, using a mask tied to that HDF5 source hash.
Diagnostic arithmetic (positions, velocities, logarithmic masses, selections
and mass bins) uses float64 even when HDF5 storage is float32. Exact-field
equality uses the original values. This does not restore precision already lost
when a mirror was written; format equivalence requires identical numeric inputs.

No Freyja configuration, HOD realisation, raw mirror or trained emulator has
been changed. The consistent downstream calculation still requires remeasuring
HMF, halo correlations and velocity targets, then validating/retraining affected
emulators against these cleaned inputs. Genuine major-merger classification and
production particle-set validation remain open scientific questions.

## Reproduce and resume

```bash
micromamba run -n cosemu python3 -B -m unittest discover \
  -s BDM-refine/analysis/duplicate-hosts-20260905/tests -p 'test_*.py' -v
micromamba run -n cosemu python3 BDM-refine/analysis/duplicate-hosts-20260905/tests/source_reproduction.py \
  PMP2linker.f90 BDM-refine/analysis/duplicate-hosts-20260905/tests/fixed --fixed
micromamba run -n cosemu python3 BDM-refine/analysis/duplicate-hosts-20260905/tests/finder_cases.py
```

The actual full finder was compiled in an isolated build directory with
`intel_comp/2024.2.0`, `compiler-rt`, `tbb`, `compiler`, and the repository's
unchanged makefile flags. It requires the Intel runtime modules to run.
The retained binary is `bin/PMP2BDM.fixed.exe`, with SHA256 recorded in
`provenance/finder-binary.sha256`; disposable compilation files were removed.
The source-based tests use gfortran bounds checking and 1, 2, 4 OpenMP threads.

The campaign script assigns a fixed subset to each of eight shared-queue
workers. A restart verifies the source and output hashes for each existing
receipt and preserves completed files. It never overwrites an existing
catalogue or sidecar. If a catalogue exists without its sidecar, the cleaner
recomputes the expected result from the source and verifies the existing file:
complete bytes for ASCII; all dataset values, dtypes, shapes and attributes for
HDF5, including non-row metadata. Only a match permits creating the missing
receipt, with the original catalogue inode and bytes preserved. Mismatches stop
without modifying that output. Recovery receipts identify the current validating
tool revision and leave the original producer commit unknown.

For new outputs, the complete sidecar is staged before either file is published;
the sidecar is linked into place last. Unique temporary names permit retries in
one worker, and handled failures remove their staging files. See
[REVIEW_FIXES.md](REVIEW_FIXES.md) for the fault-injection and precision checks.
The historical 740 ASCII products and their original receipts remain valid and
unchanged. Float32 HDF5 products made with the older arithmetic need a fresh run
to a new destination; completed receipts are not silently replaced.

After successful completion, regenerate summary tables:

```bash
micromamba run -n cosemu python3 BDM-refine/analysis/duplicate-hosts-20260905/tools/summarize_campaign.py \
  BDM-refine/analysis/duplicate-hosts-20260905/provenance/z0p25_inventory.json \
  BDM-refine/analysis/duplicate-hosts-20260905/products/strict-v1 \
  BDM-refine/analysis/duplicate-hosts-20260905 --require-complete
```

Other redshifts are inventoried in `provenance/catalogue_inventory.json` but
require their own labelled campaign and redshift mapping. Across all snapshots
there are 8,460 input catalogues totalling 7.43 TB, versus 740 / 708 GB in this
z=0.25 campaign.
