"""Reusable Matplotlib helpers for Cheng-Zong Ruan's publication style."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Sequence

import matplotlib as mpl
import matplotlib.pyplot as plt

STYLE_PATH = Path(__file__).with_name("chz-paper.mplstyle")
PALETTE = ("#1d4e89", "#20b2aa", "#edae49", "#d1495b")


def use_house_style(*, usetex: bool = True) -> None:
    """Apply the bundled rcParams, optionally falling back to mathtext."""
    plt.style.use(STYLE_PATH)
    mpl.rcParams["text.usetex"] = bool(usetex)
    if usetex:
        mpl.rcParams["text.latex.preamble"] = (
            r"\usepackage{amsmath}\usepackage{physics}"
        )


def model_data_axes(
    *,
    figsize: tuple[float, float] = (5.2, 5.3),
    height_ratios: Sequence[float] = (3.0, 1.35),
) -> tuple[mpl.figure.Figure, mpl.axes.Axes, mpl.axes.Axes]:
    """Create vertically joined main and residual axes."""
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(
        2,
        1,
        height_ratios=height_ratios,
        hspace=0.0,
    )
    ax = fig.add_subplot(gs[0])
    axr = fig.add_subplot(gs[1], sharex=ax)
    ax.tick_params(labelbottom=False)
    return fig, ax, axr


def data_errorbar_style(
    color: str,
    *,
    marker_size: float = 3.6,
    edge_width: float = 0.35,
) -> dict[str, object]:
    """Return the standard measured-data error-bar keyword arguments."""
    return {
        "fmt": "o",
        "ms": marker_size,
        "color": color,
        "mfc": color,
        "mec": "k",
        "mew": edge_width,
        "lw": 0,
        "elinewidth": 0.8,
        "capsize": 0,
        "zorder": 5,
    }


def style_residual_axis(
    ax: mpl.axes.Axes,
    *,
    band: tuple[float, float] = (-1.0, 1.0),
    ylim: tuple[float, float] | None = None,
    grid: bool = True,
) -> None:
    """Add the standard zero line, tolerance band, limits, and grid."""
    lo, hi = sorted(float(value) for value in band)
    ax.axhline(0.0, color="k", lw=0.8)
    ax.axhspan(lo, hi, color="0.88", alpha=0.7)
    if ylim is not None:
        ax.set_ylim(*ylim)
    if grid:
        ax.grid(True, ls=":", alpha=0.35)


def add_excluded_region(
    ax: mpl.axes.Axes,
    xmin: float,
    xmax: float,
    *,
    label: str | None = None,
    label_y: float | None = None,
) -> None:
    """Shade and optionally label an excluded or noise-dominated range."""
    ax.axvspan(xmin, xmax, color="lightsteelblue", alpha=0.22)
    if label is not None and label_y is not None:
        xmid = (xmin * xmax) ** 0.5 if xmin > 0 else 0.5 * (xmin + xmax)
        ax.text(
            xmid,
            label_y,
            label,
            fontsize=8,
            ha="center",
            va="center",
            color="0.4",
        )


def save_figure(
    fig: mpl.figure.Figure,
    outputs: str | Path | Iterable[str | Path],
    *,
    pad_inches: float = 0.028,
) -> None:
    """Save one or more tightly cropped outputs."""
    paths = [outputs] if isinstance(outputs, (str, Path)) else list(outputs)
    for output in paths:
        path = Path(output)
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, bbox_inches="tight", pad_inches=pad_inches)
