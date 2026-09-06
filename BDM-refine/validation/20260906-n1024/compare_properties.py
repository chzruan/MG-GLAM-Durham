#!/usr/bin/env python3
"""Same-snapshot BDM catalogue comparison; run with micromamba's cosemu Python.

Inputs are old_z*/new_z* native 24-column Catshort arrays and provenance JSON.
No candidate ID matching, inferred Rmax, or resolution-convergence claim is made.
See PLOTTING_NOTES.md for definitions, selection effects, and input units.
"""
from __future__ import annotations

import argparse
import hashlib
from importlib.metadata import version
import json
import os
from pathlib import Path
import re
import sys
from datetime import datetime, timezone

# Keep BLAS and tree queries bounded even when called inside an OMP job.
for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_name] = "1"
import numpy as np
from scipy.spatial import cKDTree


COLUMNS = ["x_mpc_h", "y_mpc_h", "z_mpc_h", "vx_km_s", "vy_km_s", "vz_km_s",
           "mass_msun_h", "total_mass_msun_h", "radius_kpc_h", "vrms_km_s",
           "vmax_km_s", "catalogue_id", "concentration", "particle_count",
           "distinct_sub_flag", "offset_over_radius", "two_K_over_Ep_minus_one",
           "spin", "radial_rms_kpc_h", "b_over_a", "c_over_a",
           "major_axis_x", "major_axis_y", "major_axis_z"]
# (column, domain, difference). Zero Vmax/Cvir/Rmax is a sentinel, not nonfinite.
PROPERTIES = {
    "mass": (6, "positive", "percent"),
    "total_mass": (7, "positive", "percent"),
    "radius": (8, "positive", "percent"),
    "vmax": (10, "sentinel", "percent"),
    "rmax": (None, "sentinel", "percent"),
    "concentration": (12, "sentinel", "percent"),
    "spin": (17, "nonnegative", "difference"),
    "b_over_a": (19, "axis_ratio", "difference"),
    "c_over_a": (20, "axis_ratio", "difference"),
    "offset": (15, "nonnegative", "difference"),
    "radial_rms": (18, "nonnegative", "percent"),
    "vrms": (9, "nonnegative", "percent"),
    "bulk_velocity": (None, "vector", "vector_norm"),
}
STATUS_NAMES = ["matched", "invalid_position_or_radius", "no_eligible_counterpart",
                "non_mutual_nearest", "outside_radius_cut", "ambiguous_neighbour"]
LIMITATIONS = [
    "HMF uses all published rows with finite positive catalogued mass, including unmatched rows; it is not a matched-sample HMF.",
    "Matches are geometric associations, not particle-identity matches. Candidate/catalogue IDs are not persistent.",
    "Matched trends are conditioned on publication by both finders, conservative positional matching, and per-property validity; they exclude disrupted/unmatched objects and selection crossings.",
    "The mass threshold applies to each finder's own published mass definition. Common MassMin does not imply an identical physical sample.",
    "The repaired bound mass, bulk velocity, energy, spin, shape and RMS use converged bound membership. The old implementation mixed populations/normalizations; differences are not numerical errors about a common definition.",
    "Reported Rvir is the extended aperture (Rext correction); it is not the unextended SO crossing radius. Xoff uses the reported aperture radius.",
    "Concentration is a mass/radius/Vmax diagnostic with a radius-based fallback, not an independently fitted profile concentration.",
    "Separate empirical corrections can make the reported c/a exceed b/a. Axis-order anomalies are counted separately; in-range calibrated axis ratios remain in matched-property statistics.",
    "Catshort contains no Rmax. Direct row-aligned auxiliary values may be supplied, but Rmax is never inferred from concentration.",
    "Quality fractions describe published catalogue rows only. Candidate rejection/status fractions require finder logs; missing candidates are not counted as valid or invalid published rows.",
    "HMF error bars are marginal sqrt(N) counting scales, not errors on the correlated old/new difference. Trend bands are 16th–84th percentiles of halo differences, not confidence intervals.",
    "Native ASCII serialization limits precision: positions have four decimal places, bulk velocities two, masses four significant digits and radii five. Zero displayed changes do not establish bitwise equality of internal properties.",
    "Low-count bins are flagged/masked, not evidence for convergence. These same-snapshot tests do not establish force, mass, time-step, or cosmological-volume convergence.",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def quantiles(values: np.ndarray) -> dict:
    finite = np.asarray(values)[np.isfinite(values)]
    result = {"n": int(len(finite))}
    for q in (1, 16, 50, 84, 99):
        result[f"p{q}"] = float(np.percentile(finite, q)) if len(finite) else None
    return result


def validate_table(table: np.ndarray, key: str) -> np.ndarray:
    if table.ndim != 2 or table.shape[1] != 24 or table.dtype.kind not in "fiu":
        raise ValueError(f"{key}: expected numeric (N,24) native Catshort table")
    return np.asarray(table, dtype=np.float64)


def match_catalogues(old: np.ndarray, new: np.ndarray, box: float, fraction: float,
                     ambiguity_ratio: float, tie_atol: float) -> dict:
    """Mutual periodic NN, d <= f min(Rold,Rnew), isolated on both sides.

    Each side's second-neighbour ambiguity criterion is d1 >= ratio*d2 OR
    d2-d1 <= tie_atol. It includes exact duplicate centres and rounding ties.
    All returned row indices refer to the input rows, never catalogue IDs.
    """
    tables = [old, new]
    eligible = [np.flatnonzero(np.all(np.isfinite(t[:, :3]), axis=1) &
                               np.isfinite(t[:, 8]) & (t[:, 8] > 0)) for t in tables]
    statuses = [np.full(len(t), 1, dtype=np.int8) for t in tables]
    result = {"old_index": np.empty(0, dtype=np.int64),
              "new_index": np.empty(0, dtype=np.int64),
              "separation_mpc_h": np.empty(0), "separation_over_min_radius": np.empty(0)}
    if not len(eligible[0]) or not len(eligible[1]):
        for status, indices in zip(statuses, eligible):
            status[indices] = 2
    else:
        positions = [np.mod(t[ix, :3], box) for t, ix in zip(tables, eligible)]
        trees = [cKDTree(p, boxsize=box) for p in positions]
        queries = [trees[1-s].query(positions[s], k=2, workers=1) for s in (0, 1)]
        distances = [q[0] for q in queries]
        nearest = [q[1][:, 0] for q in queries]
        ambiguous = [(d[:, 0] >= ambiguity_ratio * d[:, 1]) |
                     (d[:, 1] - d[:, 0] <= tie_atol) for d in distances]
        accepted = []
        for s in (0, 1):
            target = nearest[s]
            mutual = nearest[1-s][target] == np.arange(len(eligible[s]))
            radius = np.minimum(tables[s][eligible[s], 8],
                                tables[1-s][eligible[1-s][target], 8]) / 1000.0
            within = distances[s][:, 0] <= fraction * radius
            isolated = ~(ambiguous[s] | ambiguous[1-s][target])
            status = np.full(len(target), 3, dtype=np.int8)
            status[mutual] = 4
            status[mutual & within] = 5
            keep = mutual & within & isolated
            status[keep] = 0
            statuses[s][eligible[s]] = status
            accepted.append(keep)
        oi = eligible[0][accepted[0]]
        ni = eligible[1][nearest[0][accepted[0]]]
        assert len(oi) == np.count_nonzero(accepted[1])
        separation = distances[0][accepted[0], 0]
        result.update(old_index=oi, new_index=ni, separation_mpc_h=separation,
                      separation_over_min_radius=separation / (np.minimum(old[oi, 8], new[ni, 8]) / 1000.0))
    result.update(old_status=statuses[0], new_status=statuses[1])
    return result


def property_values(table: np.ndarray, name: str, auxiliary: np.ndarray | None) -> np.ndarray:
    column = PROPERTIES[name][0]
    if name == "rmax":
        return np.full(len(table), np.nan) if auxiliary is None else auxiliary
    if name == "bulk_velocity":
        return table[:, 3:6]
    return table[:, column]


def valid_values(values: np.ndarray, domain: str) -> np.ndarray:
    if domain == "vector":
        return np.all(np.isfinite(values), axis=1)
    valid = np.isfinite(values)
    if domain == "positive":
        valid &= values > 0
    elif domain in ("nonnegative", "sentinel"):
        valid &= values >= 0
    elif domain == "axis_ratio":
        valid &= (values >= 0) & (values <= 1.0005)  # native ASCII rounding
    return valid


def binned_quantiles(mass: np.ndarray, change: np.ndarray, edges: np.ndarray) -> tuple:
    finite = np.isfinite(mass) & (mass > 0) & np.isfinite(change)
    logmass = np.full(len(mass), np.nan)
    logmass[finite] = np.log10(mass[finite])
    bins = np.searchsorted(edges, logmass, side="right") - 1
    # Histogram convention includes the final right edge.
    bins[logmass == edges[-1]] = len(edges) - 2
    count = np.zeros(len(edges) - 1, dtype=np.int64)
    bands = np.full((len(count), 3), np.nan)
    for b in range(len(count)):
        values = change[finite & (bins == b)]
        count[b] = len(values)
        if len(values):
            bands[b] = np.percentile(values, [16, 50, 84])
    return count, bands


def analyse(catalogue_path: Path, metadata_path: Path, args: argparse.Namespace) -> tuple:
    metadata = json.loads(metadata_path.read_text())
    box = float(metadata["box_mpc_h"])
    nrow = int(metadata["nrow"])
    if not np.isfinite(box) or box <= 0 or nrow <= 0:
        raise ValueError("Metadata requires positive box_mpc_h and nrow")
    prepared = {}
    report = {"schema_version": 1, "created_utc": datetime.now(timezone.utc).isoformat(),
              "input_npz": str(catalogue_path.resolve()), "input_npz_sha256": sha256(catalogue_path),
              "metadata_sha256": sha256(metadata_path), "metadata": metadata,
              "script_sha256": sha256(Path(__file__)), "columns": COLUMNS,
              "software": {"python": sys.version.split()[0], "numpy": np.__version__,
                           "scipy": version("scipy"), "matplotlib": version("matplotlib")},
              "limitations": LIMITATIONS, "status_codes": dict(enumerate(STATUS_NAMES)),
              "configuration": {"box_mpc_h": box, "nrow": nrow, "min_count": args.min_count,
                                "match_fraction": args.match_fraction,
                                "ambiguity_ratio": args.ambiguity_ratio,
                                "tie_atol_mpc_h": args.tie_atol,
                                "expect_ordered_new": bool(getattr(args, "expect_ordered_new", False))}, "epochs": {}}
    report["finder_revisions"] = {
        "old": metadata.get("old_finder_commit", metadata.get("old_source_commit", "not supplied")),
        "new": metadata.get("new_finder_commit", metadata.get("source_commit", metadata.get("new_finder_build_commit", "not supplied")))}
    with np.load(catalogue_path, allow_pickle=False) as archive:
        tags = sorted((k[4:] for k in archive.files if re.fullmatch(r"old_z[0-9]+(?:p[0-9]+)?", k)),
                      key=lambda s: float(s[1:].replace("p", ".")))
        if not tags:
            raise ValueError("No old_z*/new_z* catalogue pairs found")
        if {k[4:] for k in archive.files if re.fullmatch(r"new_z[0-9]+(?:p[0-9]+)?", k)} != set(tags):
            raise ValueError("Old/new epoch sets differ; do not fabricate a missing epoch")
        tables = {f"{side}_{tag}": validate_table(archive[f"{side}_{tag}"], f"{side}_{tag}")
                  for tag in tags for side in ("old", "new")}
        masses = np.concatenate([t[np.isfinite(t[:, 6]) & (t[:, 6] > 0), 6] for t in tables.values()])
        if len(masses):
            low = np.floor(np.log10(masses.min()) * 10) / 10
            high = max(low + .1, np.ceil(np.log10(masses.max()) * 10) / 10)
        else:
            # Empty catalogues have no measured mass extent; this display range
            # is explicitly labelled and makes empty results reviewable.
            low, high = 12., 15.
        edges = np.linspace(low, high, args.mass_bins + 1)
        prepared["mass_edges_log10"] = edges
        report["mass_bins_log10"] = edges.tolist()
        report["mass_range_is_empty_display_default"] = len(masses) == 0
        for tag in tags:
            old, new = tables[f"old_{tag}"], tables[f"new_{tag}"]
            matched = match_catalogues(old, new, box, args.match_fraction, args.ambiguity_ratio, args.tie_atol)
            epoch = {"redshift": float(tag[1:].replace("p", ".")), "properties": {}, "quality": {},
                     "matching": {}, "catalogues": {}}
            report["epochs"][tag] = epoch
            for key, value in matched.items():
                prepared[f"{tag}__{key}"] = value
            oi, ni = matched["old_index"], matched["new_index"]
            epoch["matching"]["pairs"] = len(oi)
            epoch["matching"]["separation_mpc_h"] = quantiles(matched["separation_mpc_h"])
            epoch["matching"]["separation_over_min_radius"] = quantiles(matched["separation_over_min_radius"])
            prepared[f"{tag}__pair_old_mass"] = old[oi, 6]
            for side, table in (("old", old), ("new", new)):
                status = matched[f"{side}_status"]
                counts = {name: int(np.count_nonzero(status == i)) for i, name in enumerate(STATUS_NAMES)}
                epoch["matching"][side] = counts
                epoch["matching"][f"{side}_unmatched"] = len(table) - len(oi)
                good_mass = np.isfinite(table[:, 6]) & (table[:, 6] > 0)
                hmf_counts = np.histogram(np.log10(table[good_mass, 6]), edges)[0]
                matched_counts = np.histogram(np.log10(table[good_mass & (status == 0), 6]), edges)[0]
                prepared[f"{tag}__{side}_hmf_counts"] = hmf_counts
                prepared[f"{tag}__{side}_matched_mass_counts"] = matched_counts
                flags, flag_counts = np.unique(table[np.isfinite(table[:, 14]), 14], return_counts=True)
                epoch["catalogues"][side] = {"rows": len(table), "finite_positive_mass_rows": int(good_mass.sum()),
                    "nonfinite_any_rows": int(np.count_nonzero(~np.all(np.isfinite(table), axis=1))),
                    "position_outside_primary_rows": int(np.count_nonzero(np.any((table[:, :3] < 0) | (table[:, :3] >= box), axis=1))),
                    "flag_counts": {str(f): int(c) for f, c in zip(flags, flag_counts)},
                    "hmf_counts": hmf_counts.tolist(), "hmf_low_count_bins": int(np.count_nonzero((hmf_counts > 0) & (hmf_counts < args.min_count)))}
                epoch["quality"][side] = {}
                ordered_domain = np.isfinite(table[:, 19]) & np.isfinite(table[:, 20])
                order_anomaly = ordered_domain & (table[:, 20] > table[:, 19] + .0005)
                epoch["quality"][side]["axis_order"] = {
                    "available": True, "denominator_published_rows": len(table),
                    "invalid_count": int(order_anomaly.sum()),
                    "invalid_fraction": float(order_anomaly.mean()) if len(table) else None,
                    "zero_sentinel_count": None, "zero_sentinel_fraction": None,
                    "interpretation": "c/a > b/a + 5e-4 after independent empirical corrections; ordering anomaly, not a nonfinite value"}
                if side == "new" and getattr(args, "expect_ordered_new", False) and np.any(order_anomaly):
                    raise ValueError(f"{tag}: {int(order_anomaly.sum())} repaired rows violate expected c/a <= b/a ordering")
            for name, (_, domain, difference) in PROPERTIES.items():
                values = []
                available = []
                for side, table in (("old", old), ("new", new)):
                    auxiliary = None
                    key = f"{side}_rmax_{tag}"
                    if name == "rmax" and key in archive.files:
                        if metadata.get("rmax_units") != "comoving_kpc_h":
                            raise ValueError("Auxiliary Rmax requires metadata rmax_units='comoving_kpc_h'; raw finder Rmax is Mpc/h")
                        auxiliary = np.asarray(archive[key], dtype=np.float64)
                        if auxiliary.shape != (len(table),):
                            raise ValueError(f"{key}: Rmax must align exactly with catalogue rows")
                    available.append(name != "rmax" or auxiliary is not None)
                    value = property_values(table, name, auxiliary)
                    values.append(value)
                    valid = valid_values(value, domain)
                    unresolved = (valid & (value == 0)) if domain == "sentinel" else np.zeros(len(table), dtype=bool)
                    epoch["quality"][side][name] = {
                        "available": available[-1], "denominator_published_rows": len(table),
                        "invalid_count": int(np.count_nonzero(~valid)) if available[-1] else None,
                        "invalid_fraction": float(np.mean(~valid)) if len(table) and available[-1] else None,
                        "zero_sentinel_count": int(unresolved.sum()) if available[-1] and domain == "sentinel" else None,
                        "zero_sentinel_fraction": float(unresolved.mean()) if len(table) and available[-1] and domain == "sentinel" else None}
                a, b = values[0][oi], values[1][ni]
                valid = valid_values(a, domain) & valid_values(b, domain)
                change = np.full(len(oi), np.nan)
                if difference == "percent":
                    valid &= a > 0
                    if domain == "sentinel":
                        valid &= b > 0
                    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
                        change[valid] = 100. * (b[valid] / a[valid] - 1.)
                elif difference == "vector_norm":
                    change[valid] = np.linalg.norm(b[valid] - a[valid], axis=1)
                    prepared[f"{tag}__bulk_velocity_components_delta"] = b - a
                else:
                    change[valid] = b[valid] - a[valid]
                valid &= np.isfinite(change)
                change[~valid] = np.nan
                count, bands = binned_quantiles(old[oi, 6], change, edges)
                for field, value in (("old", a), ("new", b), ("change", change), ("valid", valid),
                                     ("bin_counts", count), ("bin_percentiles_16_50_84", bands)):
                    prepared[f"{tag}__{name}__{field}"] = value
                epoch["properties"][name] = {"available_both": all(available), "difference": difference,
                    "total_position_pairs": len(oi), "valid_property_pairs": int(valid.sum()),
                    "excluded_property_pairs": int((~valid).sum()), "change_percentiles": quantiles(change),
                    "mass_binned_pairs": int(count.sum()), "bins_at_min_count": int(np.count_nonzero(count >= args.min_count))}
    return prepared, report


def tex_escape(text: str) -> str:
    replacements = {"\\": r"\textbackslash{}", "&": r"\&", "%": r"\%", "$": r"\$",
                    "#": r"\#", "_": r"\_", "{": r"\{", "}": r"\}"}
    return "".join(replacements.get(c, c) for c in text)


def plot_figures(data: dict, report: dict, output: Path, usetex: bool = True) -> list[Path]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.colors import Normalize
    from house_style import use_house_style, PALETTE, save_figure
    use_house_style(usetex=usetex)
    plt.rcParams.update({"font.size": 12, "axes.labelsize": 13, "legend.fontsize": 11})
    cfg = report["configuration"]
    tags = list(report["epochs"])
    colors = {tag: PALETTE[i % len(PALETTE)] for i, tag in enumerate(tags)}
    edges = data["mass_edges_log10"]
    centres = .5 * (edges[:-1] + edges[1:])
    minimum = cfg["min_count"]
    prefix = "Small integration data" if cfg["nrow"] < 1024 else "Same-snapshot finder comparison"
    custom = report["metadata"].get("sample_label", prefix)
    context = (tex_escape(custom) if usetex else custom) + rf"; $N_\mathrm{{p}}={cfg['nrow']}^3$, $L={cfg['box_mpc_h']:g}\,h^{{-1}}\mathrm{{Mpc}}$"
    files = []

    def finish(fig, name, note):
        fig.text(.5, .965, context, ha="center", va="top", fontsize=11)
        revisions = report["finder_revisions"]
        source_text = f"Finder revisions: {str(revisions['old'])[:12]} to {str(revisions['new'])[:12]}"
        fig.text(.5, .93, tex_escape(source_text) if usetex else source_text, ha="center", va="top", fontsize=10)
        fig.text(.5, .012, note, ha="center", va="bottom", fontsize=10)
        fig.tight_layout(rect=(0, .065, 1, .90), h_pad=1.2)
        path = output / name
        save_figure(fig, path)
        plt.close(fig)
        files.append(path)

    def epoch_handles():
        return [Line2D([], [], color=colors[t], label=rf"$z={report['epochs'][t]['redshift']:g}$") for t in tags]

    # Full published populations; no positional/matching selection in the HMF.
    fig, (ax, residual) = plt.subplots(2, 1, figsize=(10.8, 6.5), sharex=True,
                                     gridspec_kw={"height_ratios": [2.8, 1.3]})
    positive = []
    for tag in tags:
        old = data[f"{tag}__old_hmf_counts"]
        new = data[f"{tag}__new_hmf_counts"]
        for side, counts, style in (("old", old, "--"), ("new", new, "-")):
            density = counts / (cfg["box_mpc_h"]**3 * np.diff(edges))
            good, low = counts >= minimum, (counts > 0) & (counts < minimum)
            positive.extend(density[counts > 0])
            ax.plot(centres, np.where(good, density, np.nan), style, color=colors[tag], lw=1.8 if side == "new" else 1.1)
            ax.errorbar(centres[good], density[good], yerr=np.sqrt(counts[good]) / (cfg["box_mpc_h"]**3 * np.diff(edges)[good]),
                        fmt="o" if side == "new" else "s", ms=3.5, mfc=colors[tag] if side == "new" else "white",
                        mec="k", mew=.4, color=colors[tag], elinewidth=.7, capsize=0)
            ax.plot(centres[low], density[low], "o" if side == "new" else "s", ms=4, mfc="none", mec=colors[tag])
        good = (old >= minimum) & (new >= minimum)
        change = np.full(len(old), np.nan)
        change[good] = 100. * (new[good] / old[good] - 1.)
        residual.plot(centres, change, "o-", color=colors[tag], ms=3.5, mec="k", mew=.4)
    ax.set_yscale("log")
    if not positive:
        ax.set_ylim(1.e-10, 1.e-3)
        ax.text(.5, .5, "No finite positive masses", transform=ax.transAxes, ha="center")
    ax.set_ylabel(r"$\mathrm{d}n/\mathrm{d}\log_{10}M\;[h^3\mathrm{Mpc}^{-3}]$")
    handles = epoch_handles() + [Line2D([], [], color=".25", ls="--", label="Pre-audit"),
                                Line2D([], [], color=".25", ls="-", label="Repaired")]
    ax.legend(handles=handles, loc="best", ncol=min(3, len(handles)))
    residual.axhline(0, color="k", lw=.8)
    residual.set_ylabel(r"$100(n_\mathrm{new}/n_\mathrm{old}-1)$")
    residual.set_xlabel(r"$\log_{10}(M_\mathrm{cat}/[h^{-1}M_\odot])$")
    residual.grid(True, ls=":", alpha=.3)
    finish(fig, "hmf.pdf", rf"All published populations. Open low-count bins: $0<N<{minimum}$; ratios require both $N\geq {minimum}$. Bars: $\sqrt{{N}}$ only.")

    labels = {
        "mass": r"$\Delta M_\mathrm{cat}/M_\mathrm{cat,old}\;[\%]$",
        "total_mass": r"$\Delta M_\mathrm{total}/M_\mathrm{total,old}\;[\%]$",
        "radius": r"$\Delta R_\mathrm{ap}/R_\mathrm{ap,old}\;[\%]$",
        "vmax": r"$\Delta V_\mathrm{max}/V_\mathrm{max,old}\;[\%]$",
        "rmax": r"$\Delta R_\mathrm{max}/R_\mathrm{max,old}\;[\%]$",
        "concentration": r"$\Delta c/c_\mathrm{old}\;[\%]$", "spin": r"$\Delta\lambda$",
        "b_over_a": r"$\Delta(b/a)$", "c_over_a": r"$\Delta(c/a)$", "offset": r"$\Delta X_\mathrm{off}$",
        "radial_rms": r"$\Delta R_\mathrm{rms}/R_\mathrm{rms,old}\;[\%]$",
        "vrms": r"$\Delta V_\mathrm{rms}/V_\mathrm{rms,old}\;[\%]$",
        "bulk_velocity": r"$|\boldsymbol{V}_\mathrm{new}-\boldsymbol{V}_\mathrm{old}|\;[\mathrm{km\,s}^{-1}]$"}
    groups = [("matched_primary.pdf", ["mass", "radius", "vmax", "concentration", "spin", "bulk_velocity"]),
              ("matched_structure.pdf", ["b_over_a", "c_over_a", "offset", "radial_rms", "vrms", "rmax"])]
    for filename, names in groups:
        fig, axes = plt.subplots(2, 3, figsize=(11.8, 6.4), sharex=True)
        for ax, name in zip(axes.flat, names):
            ax.axhline(0, color="k", lw=.8)
            available = any(report["epochs"][tag]["properties"][name]["available_both"] for tag in tags)
            if not available:
                ax.text(.5, .52, "Unavailable in paired Catshort\nNo value inferred", transform=ax.transAxes,
                        ha="center", va="center", fontsize=11)
                ax.set_yticks([])
            else:
                any_bin = False
                for tag in tags:
                    count = data[f"{tag}__{name}__bin_counts"]
                    bands = data[f"{tag}__{name}__bin_percentiles_16_50_84"].copy()
                    bands[count < minimum] = np.nan
                    any_bin |= np.any(count >= minimum)
                    ax.fill_between(centres, bands[:, 0], bands[:, 2], color=colors[tag], alpha=.16)
                    ax.plot(centres, bands[:, 1], "o-", color=colors[tag], ms=3.5, mec="k", mew=.4)
                if not any_bin:
                    ax.text(.5, .5, f"No bins with {minimum} valid pairs", transform=ax.transAxes, ha="center", fontsize=10)
            ax.set_ylabel(labels[name], fontsize=12)
            ax.grid(True, ls=":", alpha=.2)
            ax.tick_params(labelsize=10)
        axes[0, 0].legend(handles=epoch_handles(), loc="best", fontsize=10)
        for ax in axes[-1]:
            ax.set_xlabel(r"$\log_{10}(M_\mathrm{cat,old}/[h^{-1}M_\odot])$", fontsize=12)
        finish(fig, filename, rf"Conservative positional matches only. New minus old; median and 16--84 percentiles; $N_\mathrm{{valid}}\geq {minimum}$ per bin.")

    # Retain Mtotal as a separate lightweight figure rather than an overcrowded panel.
    fig, ax = plt.subplots(figsize=(7.5, 5.2))
    for tag in tags:
        bands = data[f"{tag}__total_mass__bin_percentiles_16_50_84"].copy()
        bands[data[f"{tag}__total_mass__bin_counts"] < minimum] = np.nan
        ax.fill_between(centres, bands[:, 0], bands[:, 2], color=colors[tag], alpha=.16)
        ax.plot(centres, bands[:, 1], "o-", color=colors[tag], ms=3.5, mec="k", mew=.4)
    ax.axhline(0, color="k", lw=.8)
    ax.set_ylabel(labels["total_mass"])
    ax.set_xlabel(r"$\log_{10}(M_\mathrm{cat,old}/[h^{-1}M_\odot])$")
    ax.legend(handles=epoch_handles())
    finish(fig, "matched_total_mass.pdf", "Aperture mass; matched objects only. Bands are halo-to-halo percentile ranges.")

    fig = plt.figure(figsize=(11.8, 6.5))
    grid = fig.add_gridspec(2, 2, height_ratios=[2.1, 1.25])
    coverage, distance = fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])
    table_ax = fig.add_subplot(grid[1, :]); table_ax.axis("off")
    table_rows, row_labels = [], []
    for tag in tags:
        for side, style in (("old", "--"), ("new", "-")):
            total = data[f"{tag}__{side}_hmf_counts"]
            selected = data[f"{tag}__{side}_matched_mass_counts"]
            fraction = np.full(len(total), np.nan)
            good = total >= minimum
            fraction[good] = selected[good] / total[good]
            coverage.plot(centres, fraction, style, color=colors[tag], lw=1.6)
            counts = report["epochs"][tag]["matching"][side]
            table_rows.append([f"{report['epochs'][tag]['catalogues'][side]['rows']:,}"] +
                              [f"{counts[name]:,}" for name in STATUS_NAMES])
            row_labels.append(f"{'Old' if side == 'old' else 'New'}, z={report['epochs'][tag]['redshift']:g}")
        separation = data[f"{tag}__separation_over_min_radius"]
        if len(separation):
            probability = np.linspace(0, 1, min(401, max(2, len(separation))))
            distance.plot(np.quantile(separation, probability), probability, color=colors[tag])
    coverage.set_ylim(-.02, 1.02)
    coverage.set_ylabel("Matched fraction of published rows")
    coverage.set_xlabel(r"$\log_{10}(M_\mathrm{cat}/[h^{-1}M_\odot])$")
    coverage.legend(handles=handles, ncol=2, fontsize=9, loc="best")
    distance.set_xlabel(r"$d_\mathrm{periodic}/\min(R_\mathrm{ap,old},R_\mathrm{ap,new})$")
    distance.set_ylabel("Cumulative matched fraction")
    distance.set_xscale("symlog", linthresh=1.e-4)
    distance.set_xlim(0, cfg["match_fraction"] * 1.02)
    distance.set_ylim(0, 1.02)
    distance.text(.96, .07, r"Linear for $d/R<10^{-4}$", transform=distance.transAxes,
                  ha="right", fontsize=9)
    table = table_ax.table(cellText=table_rows, rowLabels=row_labels,
                          colLabels=["Total", "Matched", "Invalid", "No target", "Non-mutual", "Too far", "Ambiguous"],
                          loc="center", cellLoc="center")
    table.auto_set_font_size(False); table.set_fontsize(10); table.scale(1, 1.3)
    for cell in table.get_celld().values():
        cell.set_edgecolor(".8"); cell.set_linewidth(.4)
    finish(fig, "matching.pdf", "Table categories are exclusive and sum to the published row count. Unmatched = all categories except Matched.")

    names = list(PROPERTIES) + ["axis_order"]
    short = [r"$M_\mathrm{cat}$", r"$M_\mathrm{tot}$", r"$R_\mathrm{ap}$", r"$V_\mathrm{max}$", r"$R_\mathrm{max}$",
             r"$c$", r"$\lambda$", r"$b/a$", r"$c/a$", r"$X_\mathrm{off}$", r"$R_\mathrm{rms}$", r"$V_\mathrm{rms}$", r"$\boldsymbol{V}$", r"$c>b$"]
    fig, axes = plt.subplots(2, 1, figsize=(11.8, 6.4))
    for ax, metric, title in zip(axes, ("invalid_fraction", "zero_sentinel_fraction"),
                                 ("Invalid values / axis-order anomalies [percent of published rows]", "Zero unresolved sentinels [percent of published rows]")):
        matrix = np.array([[report["epochs"][tag]["quality"][side][name][metric]
                            for name in names] for tag in tags for side in ("old", "new")], dtype=float) * 100
        finite = matrix[np.isfinite(matrix)]
        maximum = max(1., np.max(finite)) if len(finite) else 1.
        cmap = plt.get_cmap("viridis").copy(); cmap.set_bad(".94")
        im = ax.pcolormesh(np.arange(matrix.shape[1]+1)-.5, np.arange(matrix.shape[0]+1)-.5,
                          np.ma.masked_invalid(matrix), cmap=cmap, norm=Normalize(0, maximum), shading="flat")
        ax.set_xlim(-.5, matrix.shape[1]-.5); ax.set_ylim(matrix.shape[0]-.5, -.5)
        ax.set_xticks(np.arange(len(names))); ax.set_xticklabels(short, fontsize=11)
        ax.set_yticks(np.arange(len(row_labels))); ax.set_yticklabels(row_labels, fontsize=10)
        ax.set_title(title, fontsize=12, loc="left")
        for i, j in np.ndindex(matrix.shape):
            value = matrix[i, j]
            ax.text(j, i, "--" if not np.isfinite(value) else f"{value:.2g}", ha="center", va="center",
                    color=".35" if not np.isfinite(value) else ("white" if value < maximum * .55 else "black"), fontsize=9)
        fig.colorbar(im, ax=ax, fraction=.025, pad=.012)
    finish(fig, "quality.pdf", r"Grey: unavailable/not sentinel-defined. $c>b$: separate axis-order anomaly. Candidate rejections require finder logs.")
    return files


def self_check() -> dict:
    import tempfile
    def table(positions, radius=1000.):
        t = np.zeros((len(positions), 24), dtype=float)
        t[:, :3] = np.asarray(positions).reshape(-1, 3); t[:, 8] = radius; t[:, 6] = 1.e13
        return t
    # Independent dense periodic oracle on a tiny random cloud, plus exact edges.
    rng = np.random.default_rng(57)
    old = table(rng.uniform(0, 10, (24, 3)))
    permutation = rng.permutation(len(old))
    new = old[permutation].copy(); new[:, :3] = np.mod(new[:, :3] + .001, 10)
    result = match_catalogues(old, new, 10., .25, .5, 1.e-4)
    delta = abs(old[:, None, :3] - new[None, :, :3]); delta = np.minimum(delta, 10 - delta)
    oracle = np.sqrt(np.sum(delta**2, axis=2))
    assert np.array_equal(result["new_index"], np.argmin(oracle, axis=1))
    assert len(result["old_index"]) == 24
    periodic = match_catalogues(table([[.01, 5, 5]]), table([[9.99, 5, 5]]), 10., .25, .5, 1.e-4)
    assert np.allclose(periodic["separation_mpc_h"], .02)
    duplicate = match_catalogues(table([[0, 0, 0], [0, 0, 0]]), table([[0, 0, 0]]), 10., .25, .5, 1.e-4)
    assert not len(duplicate["old_index"]) and np.count_nonzero(duplicate["new_status"] == 5) == 1
    empty = match_catalogues(table([]), old, 10., .25, .5, 1.e-4)
    assert np.all(empty["new_status"] == 2)
    invalid = table([[np.nan, 0, 0]])
    assert match_catalogues(invalid, new, 10., .25, .5, 1.e-4)["old_status"][0] == 1
    far = match_catalogues(table([[0, 0, 0]], 1000), table([[.05, 0, 0]], 100), 10., .25, .5, 1.e-4)
    assert far["old_status"][0] == 4  # min radius, with kpc -> Mpc conversion
    count, bands = binned_quantiles(np.array([10., 100., 1000.]), np.array([1., 2., 3.]), np.array([1., 2., 3.]))
    assert np.array_equal(count, [1, 2]) and bands[1, 1] == 2.5
    assert np.array_equal(valid_values(np.array([0., -1., np.nan]), "sentinel"), [True, False, False])
    with tempfile.TemporaryDirectory(prefix="bdm-comparison-self-check-") as temporary:
        directory = Path(temporary)
        meta = directory / "metadata.json"
        fixture = directory / "catalogues.npz"
        meta.write_text(json.dumps({"box_mpc_h": 10, "nrow": 8}))
        options = argparse.Namespace(min_count=2, match_fraction=.25, ambiguity_ratio=.5,
                                     tie_atol=1.e-4, mass_bins=4, expect_ordered_new=False)
        np.savez(fixture, old_z0=table([]), new_z0=table([]))
        _, empty_report = analyse(fixture, meta, options)
        assert empty_report["epochs"]["z0"]["matching"]["pairs"] == 0
        one = table([[1., 1., 1.]])
        one[:, 19] = .4; one[:, 20] = .6
        np.savez(fixture, old_z0=one, new_z0=one, old_rmax_z0=[100.], new_rmax_z0=[200.])
        try:
            analyse(fixture, meta, options)
        except ValueError as error:
            assert "rmax_units" in str(error)
        else:
            raise AssertionError("Rmax units must be declared")
        meta.write_text(json.dumps({"box_mpc_h": 10, "nrow": 8, "rmax_units": "comoving_kpc_h"}))
        _, one_report = analyse(fixture, meta, options)
        assert one_report["epochs"]["z0"]["properties"]["rmax"]["change_percentiles"]["p50"] == 100.
        assert one_report["epochs"]["z0"]["quality"]["new"]["vmax"]["zero_sentinel_count"] == 1
        options.expect_ordered_new = True
        try:
            analyse(fixture, meta, options)
        except ValueError as error:
            assert "violate expected" in str(error)
        else:
            raise AssertionError("Expected axis order must be checked independently")
    return {"passed": True, "checks": ["dense periodic nearest-neighbour oracle", "periodic face",
            "duplicate-centre ambiguity", "empty counterpart", "invalid position", "minimum-radius cutoff and units",
            "rightmost mass-bin inclusion", "zero sentinel distinct from invalid", "empty full analysis",
            "required auxiliary Rmax units", "direct Rmax residual", "published zero-sentinel fraction",
            "independent repaired axis-order assertion"]}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    inputs = parser.add_mutually_exclusive_group()
    inputs.add_argument("--catalogues", type=Path)
    inputs.add_argument("--plot-ready", type=Path, help="Replot saved arrays without repeating catalogue matching")
    parser.add_argument("--metadata", type=Path)
    parser.add_argument("--summary", type=Path, help="Existing comparison-summary.json for --plot-ready")
    parser.add_argument("--outdir", type=Path)
    parser.add_argument("--mass-bins", type=int, default=16)
    parser.add_argument("--min-count", type=int, default=20)
    parser.add_argument("--match-fraction", type=float, default=.25)
    parser.add_argument("--ambiguity-ratio", type=float, default=.5)
    parser.add_argument("--tie-atol", type=float, default=1.e-4, help="Mpc/h; native ASCII position rounding scale")
    parser.add_argument("--expect-ordered-new", action="store_true", help="Fail if a repaired catalogue still reports c/a>b/a beyond ASCII tolerance")
    parser.add_argument("--skip-plots", action="store_true")
    parser.add_argument("--no-tex", action="store_true", help="Explicit mathtext fallback if a LaTeX installation is unavailable")
    parser.add_argument("--self-check", action="store_true")
    args = parser.parse_args()
    if args.self_check:
        print(json.dumps(self_check(), indent=2)); return
    if args.outdir is None or (args.catalogues is None and args.plot_ready is None):
        parser.error("Provide --outdir and --catalogues/--metadata or --plot-ready/--summary")
    if args.mass_bins < 1 or args.min_count < 1 or not 0 < args.match_fraction < 1 or not 0 < args.ambiguity_ratio < 1 or not np.isfinite(args.tie_atol) or args.tie_atol < 0:
        parser.error("Require positive bins/count, cutoff and ambiguity ratios in (0,1), and finite nonnegative tie tolerance")
    output = args.outdir.resolve(); output.mkdir(parents=True, exist_ok=True)
    if args.catalogues is not None:
        if args.metadata is None:
            parser.error("--catalogues requires --metadata")
        data, report = analyse(args.catalogues, args.metadata, args)
        data_path = output / "comparison-plot-ready.npz"
        summary_path = output / "comparison-summary.json"
        np.savez_compressed(data_path, **data)
        report["plot_ready_sha256"] = sha256(data_path)
        summary_path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    else:
        if args.summary is None:
            parser.error("--plot-ready requires --summary")
        report = json.loads(args.summary.read_text())
        if sha256(args.plot_ready) != report["plot_ready_sha256"]:
            raise ValueError("Saved plot-ready data hash disagrees with summary")
        with np.load(args.plot_ready, allow_pickle=False) as archive:
            data = {key: archive[key] for key in archive.files}
        data_path, summary_path = args.plot_ready, args.summary
    paths = [] if args.skip_plots else plot_figures(data, report, output, usetex=not args.no_tex)
    receipt = {"summary": str(summary_path), "summary_sha256": sha256(summary_path),
               "plot_ready": str(data_path), "plot_ready_sha256": sha256(data_path),
               "script_sha256": sha256(Path(__file__)), "tex_enabled": not args.no_tex,
               "style_sha256": {name: sha256(Path(__file__).with_name(name))
                                for name in ("house_style.py", "chz-paper.mplstyle")},
               "figures_sha256": {path.name: sha256(path) for path in paths}}
    (output / "comparison-figure-receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({"output": str(output), "epochs": {tag: {"old": epoch["catalogues"]["old"]["rows"],
          "new": epoch["catalogues"]["new"]["rows"], "matched": epoch["matching"]["pairs"]}
          for tag, epoch in report["epochs"].items()}, "figures": list(receipt["figures_sha256"])}, indent=2))


if __name__ == "__main__":
    main()
