# Legacy versus refined BDM: halo properties

This analysis compares the **published legacy and final v3 catalogues on the
same saved density field**, without changing the finder, rerunning a simulation,
or applying an additional duplicate mask. The parent is `cz` at `3fcab81`;
the working branch is `analysis/bdm-legacy-properties-20260909`.

Open the local [nine-page comparison](figures/bdm_legacy_comparison.pdf) and
[numerical results](RESULTS.md). Plot PDFs stay uncommitted, as requested.
Scripts, small plot-ready measurements, results and provenance are committed.
The frozen convergence campaign remains untouched.

## Sample and scope

The main plots use run **F**, with 1024³ particles, a 4096³ evolution mesh,
256 Mpc/h periodic box, and outputs at z=0, 1, 2. Both finders use the same
2048³ density field at each output. Legacy is `a8c7715`; refined v3 is the
audited production finder with SHA256
`39f7e514d1db9fb3b9c7545fc7f3269b4e20b5ee2d1366945f8916b997e4cd5d`.
This is the retained **native first-order ZA, z_initial=100** campaign, not
a new default-2LPTIC simulation. The source campaign records the cosmology,
IC exception, resolution changes and simulations in
[its README](../../validation/20260908-convergence/README.md).

HMF and internal-velocity profiles are also measured for all seven cases
A/B/C/D/E/F/T. Spatial and pairwise velocity statistics and matter-referenced
bias are measured for F at all three epochs. These are same-field finder
differences, not an extension of the numerical convergence pass intervals.

Use reported **bound mass** throughout. ASCII-derived `catalogues-pair.npz`
arrays give identical storage precision to the two methods. Their zero-based
columns are position 0:3 (Mpc/h), bulk peculiar velocity 3:6 (km/s), bound mass
6 (Msun/h), internal Vrms 9 and Vmax 10 (km/s). Column 11 is halo ID, not Vmax;
the separate binary membership-index layout is different. No post hoc
cleaning is imposed on the legacy sample, since doing so would change the
requested comparison.

Selections for clustering, bias and velocity moments:

- Two common thresholds: bound mass >=10^12.5 and >=10^13 Msun/h. The first
  is above the 2.5e12 Msun/h publication limit and contains at least 2362
  particles at F's resolution. Mass limits are inclusive.
- Supporting equal-density selections: above each threshold let N be the
  smaller catalogue count, then retain the N largest bound masses from each
  full catalogue. Stable original row order resolves mass ties. This is rank
  selection, **not** individual-halo matching; the realised cuts are retained.

## Estimators

**HMF:** dn/dlog10(M) = count/(L³ Δlog10M), with 0.25-dex bins beginning at
log10M=12.5. Curves retain every nonempty bin. Displayed HMF and internal
velocity residuals require at least 30 haloes/bin in *both* catalogues.
All counts, eight octant counts and unmasked profiles are saved. No octant
error is advertised as an independent-volume uncertainty.

**Real-space xi_hh:** natural periodic estimator DD/RR−1, with ordered
RR=N(N−1) (4π/3)(r_hi³−r_lo³)/L³. Pair counting uses
[pycorr/Corrfunc's periodic estimator](https://py2pcf.readthedocs.io/en/latest/api/api.html),
with `data_positions2` omitted for an autocorrelation. Self-pairs are excluded;
distinct catalogue entries at the same position are not silently removed.
Bins span 0.01–50 Mpc/h logarithmically, below L/2. Every measured selection
is independently checked against unordered `scipy.spatial.cKDTree` pairs,
with ordered counts equal to twice the unordered counts. The plot uses a
symmetric-log xi axis to show exclusion (xi=−1) as well as positive clustering.
Its residual is **100 Δxi/(1+xi_legacy)**, the relative change in pair
probability, not Δxi/xi. Residuals need at least 50 unordered pairs/bin in
both catalogues; all unmasked results are retained.

**Bias:** effective halo–matter cross bias, b_hm(k)=P_hm/P_mm. Use the
retained CIC matter contrast FI, **not** FI−1 or FI/mean(FI). Stream its
big-endian 2048³ float32 tape into a 512³ float64 block-averaged grid while
checking its full SHA256; subtract the tiny measured mean. Fields are
indexed [z,y,x]. Correct the block-centre phase, the exact discrete block
average window, and the original fine CIC window. In terms of integer mode
n per axis, the latter two windows are
`sinc(n/N_coarse)/sinc(n/N_fine)` and `sinc(n/N_fine)^2`. The block-centre
shift is `(factor−1)L/(2 N_fine)` on each axis. Direct point-sample sums
`delta_h(k)=sum exp(−i k.x)/N` avoid any halo mesh-assignment window.

Use one member of each ±k pair, exclude DC, and retain k<0.2 h/Mpc. Power
normalisation is L³ |delta(k)|²; the bias ratio cancels this factor. Within
each bin/band take **sum Re(delta_h delta_m*) / sum |delta_m|²**, not the
mean of individual ratios. There is no auto-halo shot-noise subtraction or
assumed theory bias. The primary effective-bias band is **0.05–0.15 h/Mpc**;
0.025–0.1 and 0.1–0.2 are retained as scale-band checks. This finite-band
cross bias is not a measurement of the k→0 limit. Repeat the matter
measurement on 256³ using the same input; saved Fourier coefficients and
changes quantify the coarsening sensitivity. This 256³/512³ estimator check
does not test the 2048³ *finder* mesh.

**Velocity moments:** use equal-weight unordered halo pairs and the same
periodic real-space separation bins as xi. Define
`v_r = (v_j−v_i).(x_j−x_i)_minimum_image / r`; negative means infall.
No Hubble flow or redshift-space displacement is added. Report the mean
v12, centred radial dispersion sigma_r, transverse one-dimensional RMS
`sigma_t = sqrt(<|dv|²−v_r²>/2)`, radial skewness mu3/sigma_r³, and radial
excess kurtosis mu4/sigma_r⁴−3. Central moments use a second accumulation
pass after determining the mean. The transverse RMS describes the two
tangential components about their isotropic zero-mean convention; it is
not the standard deviation of the positive transverse speed. Coincident
and underflow pairs are counted separately, not assigned a radial direction.

The source bulk velocities use the simulation's staggered output convention;
legacy/v3 within each epoch share it. The mean, skewness and kurtosis plots
show **absolute differences**, avoiding ratios of quantities that can cross
zero; dispersions show percentage differences. Curves need at least 50
pairs/bin. Additional profiles show median published internal Vrms and Vmax
at fixed mass (not matched individual haloes).

## Interpretation limits

These are paired measurements in one volume. No independent-volume covariance
or formal significance claim is supplied; pairs and Fourier modes are not
assumed to be independent error samples. A 50-pair display floor does not
make skewness/kurtosis precise. The grey ±1% band is a visual reference,
not a physics tolerance or a convergence pass criterion. Empty or excluded
residual bins remain gaps, not zero changes.

Fixed-mass differences include mass reassignment, threshold crossings,
additions/removals, and changes of centres or velocities. Rank selection
controls count differences but does not isolate those causes. No legacy
membership tapes exist for exact member matching. Smaller or larger HMF,
bias, correlation or velocity moments do not by themselves establish physical
superiority. The prior algorithm/invariant validation supports correctness;
these plots show its effect on the published observables.

## Reproduce

From this directory (always use the specified environment):

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 micromamba run -n cosemu python3 -B checks.py
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 micromamba run -n cosemu python3 -B measure.py --catalogues
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 micromamba run -n cosemu python3 -B measure.py --redshift 0 --threads 2
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 micromamba run -n cosemu python3 -B measure.py --redshift 1 --threads 2
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 micromamba run -n cosemu python3 -B measure.py --redshift 2 --threads 2
micromamba run -n cosemu python3 -B report.py
micromamba run -n cosemu python3 -B plot.py
```

The three epoch commands each read a 32-GiB density tape; run sequentially
with appropriate resources. Per-epoch receipts bind inputs, measurement
sources, output hashes, versions, CPU time and peak RSS. `checks.json`
records seven independent controls. `jobs.json` and `accounting.json` retain
the cancelled queued pilot, while `execution.json` documents the actual
low-priority login-node runs under the user's idle-node authorization.
`cleanup.json` records removal/consolidation of task-owned scratch only.
The retained `launch-logs.tar.gz` preserves launch/log provenance; from this
directory, `tar -xzf launch-logs.tar.gz` restores its original `work/` paths.
