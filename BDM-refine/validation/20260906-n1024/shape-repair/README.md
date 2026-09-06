# Empirically corrected transverse-axis ordering

Baseline: `cz` commit `b6e96669b6f7a8cba2f3dcaa9853278aac7edf11`.

The N128 catalogue plot review reported 7 of 1215 refined rows with `c/a > b/a`
beyond catalogue rounding precision. This repair reproduces the cause using
the production `GetHalo` code and controlled particle configurations.

`GetHalo` first orders the raw reduced-inertia tensor eigenvalues. It then
applies different empirical powers to the transverse ratios. Writing
`t = max(concentration_proxy - 0.4, 0)`, these powers are
`p_b = 1 + 2*t + (5.7*t)^3` and `p_c = 1 + 2*t + (5.5*t)^3`.
For `t > 0`, `p_b > p_c`; ordered inputs need not remain ordered after these
different transformations. For example, raw ratios `(0.75, 0.749)` and proxy
`0.65` give corrected values approximately `(0.282530, 0.305791)`.

The repair retains both corrected values and reports their maximum as the
intermediate-axis ratio `b/a` and their minimum as the minor-axis ratio `c/a`.
The empirical coefficients and powers remain unchanged. Since both corrected
values are in `[0, 1]`, this relabeling preserves the major axis and its stored
direction. The catalogue column layout is unchanged. An inverted pair is
swapped; an already ordered pair retains its values. There is no stored
transverse direction vector requiring a corresponding swap.

## Verification

`axis_order_regression.py` extracts the actual `GetHalo` and its dependencies
from the baseline and working source. It also extracts the actual correction
block for direct boundary and near-degeneracy tests. Generated sources,
binaries, and particle inputs are temporary and are removed after each run.

Both GNU Fortran 14.1.0 and Intel ifx 2024.2.0 passed **574 cases each**, with
checked and optimized builds. Each build covers:

- 281 correction-block cases, including an explicit near-degenerate reversal,
  zero and unit ratios, equal transverse ratios, and values on both sides of
  the empirical correction threshold. The baseline reverses 71 of these
  pairs. Every repaired pair satisfies `0 <= c/a <= b/a <= 1`; the unordered
  pair and major-axis direction remain bitwise identical to the baseline.
- Six complete `GetHalo` fixtures, each retaining all 512 input particles:
  near-degenerate transverse axes, its rotated counterpart, a nearly spherical
  halo, well-separated transverse axes, a halo below the empirical correction
  threshold, and a rank-one configuration. The near-degenerate fixture changes
  `(b/a, c/a)` from `(0.0964788, 0.1173751)` to `(0.1173751, 0.0964788)`.
  All other tested fields, including masses, energies, velocities, spin,
  position, shape major direction, status, and exact membership bytes, are
  unchanged. An independently accumulated reduced tensor confirms that the
  stored major direction remains a unit principal eigenvector.

These bounded checks validate the repair without a new production simulation.
The main validation workflow preserves the earlier pilot separately and must
rebuild the finder with this commit before subsequent production runs.

Recorded evidence: `gnu-results.json` and `ifx-results.json`, including compiler
commands, baseline/source/test SHA-256 hashes, per-fixture values, and member
hashes. Source SHA-256:
`1353532b32f6f3f80457208f5d023f8f359097a35a24edfc41856453a21b154c`.

From the repository root:

```bash
micromamba run -n cosemu python3 -B BDM-refine/validation/20260906-n1024/shape-repair/axis_order_regression.py --compiler gfortran --output /tmp/bdm-axis-gnu.json
```

For the native Intel check, load the same compiler/runtime modules and expose
their library path to the harness so the compiler and generated executable
use the matching Intel runtime:

```bash
module purge
module load intel_comp/2024.2.0
module load compiler-rt tbb compiler
export BDM_AUDIT_NATIVE_LIBS="$LD_LIBRARY_PATH"
micromamba run -n cosemu python3 -B BDM-refine/validation/20260906-n1024/shape-repair/axis_order_regression.py --compiler ifx --output /tmp/bdm-axis-ifx.json
```
