# Measured legacy versus refined BDM differences

Main sample: F (1024³ particles, 4096³ evolution mesh, L=256 Mpc/h), z=0, 1, 2. Both finders act on the same 2048³ density field. The [nine-page PDF](figures/bdm_legacy_comparison.pdf) contains all four requested statistics and supporting selection/internal-velocity checks.

**These plots quantify the effect of the repairs; they do not demonstrate a new convergence advantage or supply independent-volume uncertainties.** Full definitions and display cuts are in [README.md](README.md).

## Catalogue abundance

| z | Legacy rows | Refined rows | Difference | Change (%) |
|---|---:|---:|---:|---:|
| 0 | 28,262 | 28,226 | -36 | -0.127 |
| 1 | 20,794 | 20,777 | -17 | -0.082 |
| 2 | 8,674 | 8,754 | +80 | +0.922 |

Total rows use the native publication limit (2.5e12 Msun/h); the plots and clustering samples start at 10^12.5 Msun/h. Changes in the HMF are mass-dependent:

| z | Mass bin, log10(Mbound/[Msun/h]) | Legacy | Refined | HMF change (%) |
|---|---|---:|---:|---:|
| 0 | 12.50–12.75 | 9551 | 9547 | -0.042 |
| 0 | 12.75–13.00 | 5897 | 5876 | -0.356 |
| 0 | 13.00–13.25 | 3400 | 3390 | -0.294 |
| 0 | 13.25–13.50 | 1915 | 1921 | +0.313 |
| 0 | 13.50–13.75 | 1029 | 1020 | -0.875 |
| 0 | 13.75–14.00 | 558 | 554 | -0.717 |
| 0 | 14.00–14.25 | 264 | 265 | +0.379 |
| 0 | 14.25–14.50 | 124 | 126 | +1.613 |
| 0 | 14.50–14.75 | 54 | 54 | +0.000 |
| 1 | 12.50–12.75 | 8109 | 8106 | -0.037 |
| 1 | 12.75–13.00 | 4131 | 4128 | -0.073 |
| 1 | 13.00–13.25 | 1966 | 1956 | -0.509 |
| 1 | 13.25–13.50 | 948 | 949 | +0.105 |
| 1 | 13.50–13.75 | 382 | 381 | -0.262 |
| 1 | 13.75–14.00 | 149 | 153 | +2.685 |
| 1 | 14.00–14.25 | 43 | 46 | +6.977 |
| 2 | 12.50–12.75 | 3685 | 3711 | +0.706 |
| 2 | 12.75–13.00 | 1523 | 1536 | +0.854 |
| 2 | 13.00–13.25 | 542 | 554 | +2.214 |
| 2 | 13.25–13.50 | 184 | 189 | +2.717 |
| 2 | 13.50–13.75 | 41 | 41 | +0.000 |

Only bins with at least 30 haloes in both catalogues appear above. Sparse-bin fractional changes should be read alongside their counts.

The common mass cut does not represent the same particle count in the supporting simulations. At 10^12.5 Msun/h, A has about 37 particles per halo, B/C/D about 296 and E/F/T about 2362. The lowest-bin (12.50–12.75) HMF changes are:

| Case | Particle floor at common mass cut | z=0 change (%) | z=1 change (%) | z=2 change (%) |
|---|---:|---:|---:|---:|
| A | 37 | -18.777 | -18.672 | -18.724 |
| B | 296 | +2.109 | +2.058 | +1.428 |
| C | 296 | +0.198 | +0.264 | +0.057 |
| D | 296 | -0.530 | -0.619 | -0.138 |
| E | 2362 | +0.650 | +0.704 | +1.701 |
| F | 2362 | -0.042 | -0.037 | +0.706 |
| T | 2362 | -0.115 | -0.049 | +0.659 |

The much larger A response is in a poorly resolved mass range. The small F differences should not be generalised to this sample.

## Matter-referenced effective bias

The band estimator is sum P_hm / sum P_mm over **0.05 ≤ k < 0.15 h/Mpc**. This is a finite-band cross bias, not a fitted k→0 limit.

| z | Selection | N legacy / refined | b legacy | b refined | Change (%) |
|---|---|---:|---:|---:|---:|
| 0 | M ≥ 10^12.5 | 22,806 / 22,767 | 1.09412 | 1.08351 | -0.9702 |
| 0 | M ≥ 10^13 | 7,358 / 7,344 | 1.39869 | 1.38827 | -0.7449 |
| 0 | equal n, 12.5 cut | 22,767 / 22,767 | 1.09386 | 1.08351 | -0.9460 |
| 0 | equal n, 13 cut | 7,344 / 7,344 | 1.39815 | 1.38827 | -0.7068 |
| 1 | M ≥ 10^12.5 | 15,738 / 15,729 | 2.04536 | 2.02736 | -0.8797 |
| 1 | M ≥ 10^13 | 3,498 / 3,495 | 2.77809 | 2.76743 | -0.3836 |
| 1 | equal n, 12.5 cut | 15,729 / 15,729 | 2.04593 | 2.02736 | -0.9073 |
| 1 | equal n, 13 cut | 3,495 / 3,495 | 2.78054 | 2.76743 | -0.4712 |
| 2 | M ≥ 10^12.5 | 5,979 / 6,035 | 3.53228 | 3.50954 | -0.6438 |
| 2 | M ≥ 10^13 | 771 / 788 | 4.93971 | 4.91848 | -0.4298 |
| 2 | equal n, 12.5 cut | 5,979 / 5,979 | 3.53228 | 3.51833 | -0.3949 |
| 2 | equal n, 13 cut | 771 / 771 | 4.93971 | 4.93857 | -0.0231 |

Changing the matter measurement grid from 512³ to 256³ changes the effective bias by at most 0.00703% and the paired bias change by at most 0.00050 percentage points. Direct halo Fourier sums use no halo assignment grid. This validates measurement coarsening, not finder-mesh independence.

## Clustering and pairwise velocities

The following are **maximum absolute per-bin changes at 5–30 Mpc/h**, using geometric bin centres and at least 50 unordered pairs/bin in both catalogues. They are descriptive extrema, not error bars. The correlation column is the change in **1+xi**, not the percentage change in xi.

| z | log10 mass cut | Δ(1+xi)/(1+xi) (%) | Δv12 (km/s) | Δsigma_r (%) | Δsigma_t (%) | Δskewness | Δexcess kurtosis |
|---|---|---:|---:|---:|---:|---:|---:|
| 0 | 12.5 | 1.046 | 2.962 | 0.785 | 1.795 | 0.0115 | 0.1598 |
| 0 | 13.0 | 0.386 | 2.469 | 3.195 | 2.677 | 0.0541 | 0.3694 |
| 1 | 12.5 | 0.908 | 1.719 | 1.739 | 3.111 | 0.0180 | 0.1142 |
| 1 | 13.0 | 1.019 | 1.711 | 2.111 | 3.167 | 0.0755 | 0.2432 |
| 2 | 12.5 | 0.617 | 1.867 | 1.581 | 2.328 | 0.0188 | 0.0650 |
| 2 | 13.0 | 4.498 | 8.870 | 2.953 | 3.461 | 0.1678 | 0.2764 |

The plots retain the radial dependence, smaller-scale changes and gaps from the pair-count floor. `results.json` also reports 0.5–1 and 1–5 Mpc/h extrema, both equal-density controls, and ordinary Δxi/xi where |xi_legacy| > 0.1. All unmasked bins, raw pair counts, moments and low-k Fourier coefficients are in the three `F-z*-statistics.npz` files.

**Below 1 Mpc/h the velocity changes are substantially larger.** For the 10^12.5 Msun/h cut, the 0.5–1 Mpc/h bins meeting the same 50-pair floor give:

| z | Bins | Minimum pairs/catalogue/bin | Max abs Δv12 (km/s) | Max abs Δsigma_r (%) | Max abs Δsigma_t (%) |
|---|---:|---:|---:|---:|---:|
| 0 | 2 | 182 | 32.351 | 21.470 | 31.507 |
| 1 | 2 | 123 | 73.986 | 14.199 | 40.270 |
| 2 | 1 | 102 | 33.977 | 13.340 | 13.715 |

### Very close halo pairs

The spikes at the smallest separations in the legacy xi curve are produced by just a few pairs. The table uses exact periodic centre distances below 0.1 Mpc/h. A plotted bin centred below 0.1 can include pairs above 0.1; these counts are independent of the plotting bins:

| z | Mass cut | Legacy pairs with r < 0.1 Mpc/h | Refined pairs |
|---|---|---:|---:|
| 0 | 10^12.5 | 1 | 0 |
| 0 | 10^13 | 0 | 0 |
| 1 | 10^12.5 | 5 | 0 |
| 1 | 10^13 | 0 | 0 |
| 2 | 10^12.5 | 4 | 0 |
| 2 | 10^13 | 1 | 0 |

The refinement removes these very close pairs and produces a clear exclusion region. Centre distances alone do not establish that a pair had exactly identical particle membership; the legacy memberships were not retained. The spikes are not a well-sampled clustering signal.

## Internal halo velocities

| z | Profile | Maximum absolute change of bin median (%) |
|---|---|---:|
| 0 | internal_vrms | 3.188 |
| 0 | vmax | 1.072 |
| 1 | internal_vrms | 4.155 |
| 1 | vmax | 1.664 |
| 2 | internal_vrms | 2.927 |
| 2 | vmax | 1.195 |

These use fixed mass bins with 30 haloes per catalogue; they do not match individual haloes.

## Validation and provenance

- All 21 catalogue-pair files were checked against the frozen replay hashes.
- All three 32-GiB matter tapes were hashed during streaming and matched their replay receipts.
- Seven independent small controls passed; all 24 measured selections have exact agreement between Corrfunc ordered and KDTree unordered pair counts.
- The summary independently reconstructs the periodic RR normalisation and band bias from retained counts/Fourier modes; input, output and source hashes are checked.
- The single shared Slurm pilot was cancelled while pending, with zero allocated/billed time. Sequential low-priority login-node execution used the user’s idle-node authorization; actual time/memory are in `execution.json`.
- No production source, simulation input or frozen convergence result was changed.
