"""Matplotlib plotting for CosMx QC.

Style matches the original Python QC visualisations:

- Bar / box charts: figsize=(16, 8), title 14pt bold, axis labels 12pt,
  grid alpha=0.3 axis=y, bar alpha 0.7-0.8.
- FOV-position scatter: figsize=(12, 10), tab20 cmap by FOV index,
  CosMx-standard 0.51 mm circles, equal aspect, inverted y for stage
  coordinates, colorbar on the side.
- Per-target colours: Gene #E74C3C, Negative #5DADE2, SystemControl #52C785.

Each plot type ships in two flavours:
  * `*_ax(ax, df, ...)` draws into a caller-provided Axes. Used by the
    Combined view, which runs `plt.subplots(1, N, sharey=True)` so the
    Y axis is *physically* shared across sample panels.
  * `*(df, ...)` is a thin wrapper that opens its own Figure and delegates
    to the `_ax` variant. Used by the per-sample tabs.
"""

from __future__ import annotations

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.figure import Figure

# Public colour palette - keyed by metric kind
COLORS = {
    "gene": "#E74C3C",  # red - real genes
    "negative": "#5DADE2",  # blue - negative probes
    "falsecode": "#52C785",  # green - system-control / false codes
    "assigned": "green",
    "unassigned": "red",
    "detected": "green",
    "undetected": "red",
    "neutral": "#A0A0A0",  # medium grey — same visual weight as boxplot
}

FIGSIZE_BAR = (16, 8)
FIGSIZE_POSITIONS = (12, 10)
TITLE_FONTSIZE = 14
LABEL_FONTSIZE = 12
BOX_TITLE_FONTSIZE = 16
BOX_LABEL_FONTSIZE = 14


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

_BAR_WIDTH = 0.8


def _set_fov_xticks(
    ax: Axes,
    fovs: np.ndarray,
    rotation: int = 90,
    ha: str = "right",
    va: str = "center",
    fontsize: int = 8,
) -> np.ndarray:
    """Place evenly-spaced integer FOV labels on `ax`. Returns x positions.

    Style is unified across bar and box plots: 90° rotation, 8-pt,
    anchored on the right edge of the rotated label so the top of each
    vertical label sits directly under its tick. `rotation_mode="anchor"`
    is essential — without it matplotlib rotates around the bbox centre
    and labels appear to drift off the ticks for multi-digit FOVs.

    Also pins x-limits to a uniform symmetric padding so bar and box
    panels share identical left/right whitespace.
    """
    xs = np.arange(len(fovs))
    ax.set_xticks(xs)
    ax.set_xticklabels(
        [str(int(v)) for v in fovs],
        rotation=rotation,
        ha=ha,
        va=va,
        fontsize=fontsize,
        rotation_mode="anchor",
    )
    ax.set_xlim(-0.5 - _BAR_WIDTH / 2, len(xs) - 0.5 + _BAR_WIDTH / 2)
    return xs


def _bar_finish(ax: Axes, title: str, ylabel: str, show_ylabel: bool = True) -> None:
    ax.set_xlabel("FOV", fontsize=LABEL_FONTSIZE)
    if show_ylabel:
        ax.set_ylabel(ylabel, fontsize=LABEL_FONTSIZE)
    if title:
        ax.set_title(title, fontsize=TITLE_FONTSIZE, fontweight="bold")
    ax.grid(True, alpha=0.3, axis="y")


def _box_finish(ax: Axes, title: str, ylabel: str, show_ylabel: bool = True) -> None:
    ax.set_xlabel("FOV", fontsize=BOX_LABEL_FONTSIZE)
    if show_ylabel:
        ax.set_ylabel(ylabel, fontsize=BOX_LABEL_FONTSIZE)
    if title:
        ax.set_title(title, fontsize=BOX_TITLE_FONTSIZE, fontweight="bold")
    ax.grid(True, alpha=0.3, axis="y")


# ---------------------------------------------------------------------------
# FOV positions (Profiling Position)
# ---------------------------------------------------------------------------


def plot_fov_positions(fov_pos: pd.DataFrame, title: str = "") -> Figure:
    """Scatter FOV index on (X_mm, Y_mm) with CosMx-standard 0.51 mm circles.

    Colour by FOV number (tab20), labels at circle centre, equal aspect,
    inverted y-axis to match microscope-stage coordinate convention.
    """
    fig, ax = plt.subplots(figsize=FIGSIZE_POSITIONS)

    x_range = fov_pos["X_mm"].max() - fov_pos["X_mm"].min()
    y_range = fov_pos["Y_mm"].max() - fov_pos["Y_mm"].min()
    pad_x = max(x_range * 0.1, 1.0)
    pad_y = max(y_range * 0.1, 1.0)
    ax.set_xlim(fov_pos["X_mm"].min() - pad_x, fov_pos["X_mm"].max() + pad_x)
    ax.set_ylim(fov_pos["Y_mm"].min() - pad_y, fov_pos["Y_mm"].max() + pad_y)

    cmap = plt.cm.tab20
    norm = plt.Normalize(vmin=fov_pos["FOV"].min(), vmax=fov_pos["FOV"].max())
    radius = 0.51 / 2  # CosMx FOV physical diameter is 0.51 mm

    for _, row in fov_pos.iterrows():
        c = cmap(norm(row["FOV"]))
        ax.add_patch(
            patches.Circle(
                (row["X_mm"], row["Y_mm"]),
                radius=radius,
                facecolor=c,
                edgecolor=c,
                alpha=0.6,
                linewidth=1,
            )
        )
        ax.text(
            row["X_mm"],
            row["Y_mm"],
            str(int(row["FOV"])),
            fontsize=8,
            ha="center",
            va="center",
            fontweight="bold",
            color="black",
        )

    ax.set_xlabel("X Global (mm)", fontsize=LABEL_FONTSIZE)
    ax.set_ylabel("Y Global (mm)", fontsize=LABEL_FONTSIZE)
    if title:
        ax.set_title(title, fontsize=TITLE_FONTSIZE, fontweight="bold")
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.3)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, shrink=0.6)
    cbar.set_label("FOV Number", fontsize=LABEL_FONTSIZE)

    ax.invert_yaxis()
    return fig


# ---------------------------------------------------------------------------
# bar charts (single-colour and stacked) -- ax + wrapper variants
# ---------------------------------------------------------------------------


def plot_count_per_fov_ax(
    ax: Axes,
    df: pd.DataFrame,
    *,
    y_col: str,
    ylabel: str,
    title: str = "",
    color: str = "#888888",
    y_max: float | None = None,
    show_ylabel: bool = True,
) -> None:
    """Single-colour bar chart drawn into `ax`."""
    fovs = df["fov"].values
    xs = _set_fov_xticks(ax, fovs)
    ax.bar(xs, df[y_col].values, width=0.8, color=color, alpha=0.8)
    if y_max is not None:
        ax.set_ylim(0, y_max)
    _bar_finish(ax, title, ylabel, show_ylabel=show_ylabel)


def plot_count_per_fov(df: pd.DataFrame, **kwargs) -> Figure:
    fig, ax = plt.subplots(figsize=FIGSIZE_BAR)
    plot_count_per_fov_ax(ax, df, **kwargs)
    return fig


def plot_stacked_two_per_fov_ax(
    ax: Axes,
    df: pd.DataFrame,
    *,
    bottom_col: str,
    top_col: str,
    bottom_label: str,
    top_label: str,
    bottom_color: str = "green",
    top_color: str = "red",
    title: str = "",
    show_ylabel: bool = True,
    show_legend: bool = False,
) -> None:
    """Generic stacked two-class percent bar drawn into `ax`.

    `df[bottom_col]` and `df[top_col]` are values in [0, 1]; the bar
    height is rendered in 0–100 % space.
    """
    fovs = df["fov"].values
    xs = _set_fov_xticks(ax, fovs)
    bot = df[bottom_col].values * 100.0
    top = df[top_col].values * 100.0
    ax.bar(xs, bot, width=0.8, label=bottom_label, color=bottom_color, alpha=0.7)
    ax.bar(xs, top, width=0.8, bottom=bot, label=top_label, color=top_color, alpha=0.7)
    ax.set_ylim(0, 100)
    if show_legend:
        ax.legend(loc="lower left", fontsize=10)
    _bar_finish(ax, title, "Percent (%)", show_ylabel=show_ylabel)


def plot_assignment_per_fov_ax(ax: Axes, df: pd.DataFrame, **kwargs) -> None:
    """Assigned (green) / Unassigned (red) stacked-percent bar."""
    plot_stacked_two_per_fov_ax(
        ax,
        df,
        bottom_col="assigned",
        top_col="unassigned",
        bottom_label="Assigned",
        top_label="Unassigned",
        bottom_color=COLORS["assigned"],
        top_color=COLORS["unassigned"],
        **kwargs,
    )


def plot_assignment_per_fov(df: pd.DataFrame, **kwargs) -> Figure:
    fig, ax = plt.subplots(figsize=FIGSIZE_BAR)
    plot_assignment_per_fov_ax(ax, df, **kwargs)
    return fig


def plot_detection_per_fov_ax(ax: Axes, df: pd.DataFrame, **kwargs) -> None:
    """Detected (green) / Undetected (red) stacked-percent bar."""
    plot_stacked_two_per_fov_ax(
        ax,
        df,
        bottom_col="detected",
        top_col="undetected",
        bottom_label="Detected",
        top_label="Undetected",
        bottom_color=COLORS["detected"],
        top_color=COLORS["undetected"],
        **kwargs,
    )


def plot_detection_per_fov(df: pd.DataFrame, **kwargs) -> Figure:
    fig, ax = plt.subplots(figsize=FIGSIZE_BAR)
    plot_detection_per_fov_ax(ax, df, **kwargs)
    return fig


def plot_stacked_kinds_per_fov_ax(
    ax: Axes,
    df: pd.DataFrame,
    *,
    title: str = "",
    percent: bool = False,
    y_max: float | None = None,
    show_ylabel: bool = True,
    show_legend: bool = False,
) -> None:
    """Stacked Gene / Negative / SystemControl bar drawn into `ax`.

    Legend defaults to off because the report places it as markdown
    under the section title; pass `show_legend=True` for stand-alone use.
    """
    fovs = df["fov"].values
    xs = _set_fov_xticks(ax, fovs)

    gene = df["gene"].values.astype(float)
    neg = df["negative"].values.astype(float)
    fc = df["falsecode"].values.astype(float)
    if percent:
        total = gene + neg + fc
        total = np.where(total == 0, 1.0, total)
        gene = gene / total * 100.0
        neg = neg / total * 100.0
        fc = fc / total * 100.0
        ylabel = "Percent (%)"
    else:
        ylabel = "Count"

    ax.bar(xs, gene, width=0.8, label="Gene", color=COLORS["gene"], alpha=0.8)
    ax.bar(xs, neg, width=0.8, bottom=gene, label="Negative", color=COLORS["negative"], alpha=0.8)
    ax.bar(
        xs,
        fc,
        width=0.8,
        bottom=gene + neg,
        label="SystemControl",
        color=COLORS["falsecode"],
        alpha=0.8,
    )

    if percent:
        ax.set_ylim(0, 100)
    elif y_max is not None:
        ax.set_ylim(0, y_max)
    if show_legend:
        ax.legend(loc="upper right", fontsize=10)
    _bar_finish(ax, title, ylabel, show_ylabel=show_ylabel)


def plot_stacked_kinds_per_fov(df: pd.DataFrame, **kwargs) -> Figure:
    fig, ax = plt.subplots(figsize=FIGSIZE_BAR)
    plot_stacked_kinds_per_fov_ax(ax, df, **kwargs)
    return fig


# ---------------------------------------------------------------------------
# box plots -- ax + wrapper variants
# ---------------------------------------------------------------------------


def plot_box_per_fov_ax(
    ax: Axes,
    df: pd.DataFrame,
    *,
    y_col: str,
    ylabel: str,
    title: str = "",
    color: str = "lightgray",
    y_max: float | None = None,
    show_ylabel: bool = True,
) -> None:
    """Boxplot of `y_col` grouped by FOV drawn into `ax`.

    Positions and widths match the bar plots (0..N-1, width 0.8) so
    side-by-side bar and box panels share identical x-axis spacing.
    """
    fovs = np.array(sorted(df["fov"].unique()))
    groups = [df.loc[df["fov"] == f, y_col].dropna().values for f in fovs]
    positions = np.arange(len(fovs))

    ax.boxplot(
        groups,
        positions=positions,
        widths=_BAR_WIDTH,
        patch_artist=True,
        boxprops=dict(facecolor=color, alpha=0.7),
        medianprops=dict(color="black", linewidth=2),
        whiskerprops=dict(color="black"),
        capprops=dict(color="black"),
        flierprops=dict(marker="o", markerfacecolor="gray", markersize=2, alpha=0.3),
    )
    _set_fov_xticks(ax, fovs)
    if y_max is not None:
        ax.set_ylim(0, y_max)
    _box_finish(ax, title, ylabel, show_ylabel=show_ylabel)


def plot_box_per_fov(df: pd.DataFrame, **kwargs) -> Figure:
    fig, ax = plt.subplots(figsize=FIGSIZE_BAR)
    plot_box_per_fov_ax(ax, df, **kwargs)
    return fig


# ---------------------------------------------------------------------------
# Combined-figure helper (sharey)
# ---------------------------------------------------------------------------


def make_combined_sharey(
    by_sample: dict[str, pd.DataFrame], plot_ax_fn, *, figsize: tuple | None = None, **plot_kwargs
) -> Figure:
    """Build a single Figure with N panels (one per sample), `sharey=True`.

    The Combined Figure has the same overall size as a single-sample
    Figure (default `FIGSIZE_BAR`) so it occupies the same screen real
    estate; each sub-panel is therefore ~1/N as wide. Panels touch
    (`wspace=0`); only the leftmost shows Y axis labels and ticks.

    `plot_ax_fn(ax, df, *, title, show_ylabel, **plot_kwargs)` draws
    into each axes.
    """
    n = len(by_sample)
    if figsize is None:
        figsize = FIGSIZE_BAR
    fig, axes = plt.subplots(1, n, figsize=figsize, sharey=True)
    if n == 1:
        axes = [axes]
    fig.subplots_adjust(wspace=0.05)
    for i, (name, df) in enumerate(by_sample.items()):
        plot_ax_fn(axes[i], df, title=name, show_ylabel=(i == 0), **plot_kwargs)
    return fig


# ---------------------------------------------------------------------------
# expression heatmap (kept multi-panel; per-sample-only — no sharey across samples)
# ---------------------------------------------------------------------------


def plot_expression_heatmap(by_kind: dict[str, pd.DataFrame], title: str = "") -> Figure:
    """3-panel viridis heatmap of per-FOV mean expression.

    `by_kind` is a dict mapping {"gene", "negative", "falsecode"} → a
    DataFrame indexed by FOV with target genes as columns. Empty / missing
    kinds are skipped. Each panel shares y-axis (FOV) but has its own
    target-column extent. Color is fixed to viridis with vmin=0, vmax=5
    (matching the original Python QC output).
    """
    label_map = {"gene": "Gene", "negative": "Negative", "falsecode": "SystemControl"}
    order = ["negative", "falsecode", "gene"]
    available = [k for k in order if k in by_kind and not by_kind[k].empty]
    if not available:
        fig, ax = plt.subplots(figsize=FIGSIZE_BAR)
        ax.text(0.5, 0.5, "No expression data", ha="center", va="center")
        ax.set_axis_off()
        return fig

    fig, axes = plt.subplots(1, len(available), figsize=(5 * len(available), 8))
    if len(available) == 1:
        axes = [axes]

    cmap = plt.cm.viridis.copy()
    cmap.set_bad(color="white", alpha=1.0)

    all_fovs = set()
    for k in available:
        all_fovs.update(by_kind[k].index)
    fovs = sorted(all_fovs)

    im = None
    for idx, kind in enumerate(available):
        df = by_kind[kind].reindex(fovs)
        ax = axes[idx]
        im = ax.imshow(df.values, cmap=cmap, aspect="auto", interpolation="nearest", vmin=0, vmax=5)
        ax.set_title(label_map[kind], fontsize=14, fontweight="bold")
        if idx == 0:
            ax.set_yticks(range(len(fovs)))
            ax.set_yticklabels([str(int(f)) for f in fovs], fontsize=6)
            ax.set_ylabel("FOV", fontsize=LABEL_FONTSIZE)
        else:
            ax.set_yticks([])
        ax.set_xticks([])

    if title:
        fig.suptitle(title, fontsize=BOX_TITLE_FONTSIZE, fontweight="bold")

    fig.subplots_adjust(right=0.9, top=0.92 if title else 0.95, wspace=0.05)
    cbar_ax = fig.add_axes([0.92, 0.15, 0.015, 0.7])
    fig.colorbar(im, cax=cbar_ax)
    return fig
