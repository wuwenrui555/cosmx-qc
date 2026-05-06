"""Pure-data QC metric functions. Each takes a SampleData (or dict thereof)
and returns a DataFrame.  No plotting, no IO."""

from __future__ import annotations

import logging
import time
from contextlib import contextmanager
from typing import Literal

import numpy as np
import pandas as pd

from cosmx_qc.io import SampleData

logger = logging.getLogger("cosmx_qc")

TargetKind = Literal["all", "gene", "negative", "falsecode"]
META_COLS = {"fov", "cell_ID", "cell"}  # known non-expression columns in exprMat


@contextmanager
def _timed(label: str):
    """Log start + duration around a block. Used by every metric."""
    t0 = time.monotonic()
    logger.info("%s start", label)
    try:
        yield
    finally:
        logger.info("%s done %.1fs", label, time.monotonic() - t0)


def classify_expr_columns(exprmat: pd.DataFrame) -> dict[str, list[str]]:
    """Split exprMat columns into gene / negative / falsecode."""
    gene, negative, falsecode = [], [], []
    for c in exprmat.columns:
        if c in META_COLS:
            continue
        if c.startswith("Negative"):
            negative.append(c)
        elif c.startswith("SystemControl"):
            falsecode.append(c)
        else:
            gene.append(c)
    return {"gene": gene, "negative": negative, "falsecode": falsecode}


def classify_target(target: pd.Series) -> pd.Series:
    """Element-wise classification of a tx-file `target` series.

    Most call sites should use `SampleData.tx_kind` (cached) instead.
    """
    out = pd.Series("gene", index=target.index, dtype=object)
    out[target.str.startswith("Negative", na=False)] = "negative"
    out[target.str.startswith("SystemControl", na=False)] = "falsecode"
    return out


def transcripts_per_fov(sd: SampleData, kind: TargetKind = "all") -> pd.DataFrame:
    """Count transcripts per FOV from tx_file, optionally filtered by kind."""
    with _timed(f"[{sd.name}] transcripts_per_fov(kind={kind})"):
        df = sd.tx
        if kind != "all":
            df = df[sd.tx_kind == kind]
        counts = df.groupby("fov").size().reset_index(name="count")
        return counts.sort_values("fov").reset_index(drop=True)


def unassigned_transcripts_per_fov(sd: SampleData) -> pd.DataFrame:
    """Count transcripts with cell_ID == 0 per FOV."""
    with _timed(f"[{sd.name}] unassigned_transcripts_per_fov"):
        df = sd.tx[sd.tx["cell_ID"] == 0]
        counts = df.groupby("fov").size().reset_index(name="unassigned")
        # ensure all FOVs present (zero-fill missing)
        all_fovs = pd.DataFrame({"fov": sorted(sd.tx["fov"].unique())})
        out = all_fovs.merge(counts, on="fov", how="left").fillna({"unassigned": 0})
        out["unassigned"] = out["unassigned"].astype(int)
        return out


def assignment_ratio_per_fov(sd: SampleData) -> pd.DataFrame:
    """Per-FOV (assigned, unassigned) proportions, each in [0, 1]."""
    with _timed(f"[{sd.name}] assignment_ratio_per_fov"):
        total = sd.tx.groupby("fov").size()
        unassigned = sd.tx[sd.tx["cell_ID"] == 0].groupby("fov").size()
        out = pd.DataFrame({"fov": total.index, "_total": total.values})
        u = out["fov"].map(unassigned).fillna(0).astype(float).values
        out["unassigned"] = u / total.values
        out["assigned"] = 1.0 - out["unassigned"]
        return out[["fov", "assigned", "unassigned"]].sort_values("fov").reset_index(drop=True)


def transcripts_breakdown_per_fov(sd: SampleData) -> pd.DataFrame:
    """Per-FOV counts split by transcript kind (gene / negative / falsecode)."""
    with _timed(f"[{sd.name}] transcripts_breakdown_per_fov"):
        grp = pd.crosstab(sd.tx["fov"], sd.tx_kind)
        for k in ("gene", "negative", "falsecode"):
            if k not in grp.columns:
                grp[k] = 0
        return (
            pd.DataFrame(
                {
                    "fov": grp.index,
                    "gene": grp["gene"].values,
                    "negative": grp["negative"].values,
                    "falsecode": grp["falsecode"].values,
                }
            )
            .sort_values("fov")
            .reset_index(drop=True)
        )


def unique_transcripts_per_fov(sd: SampleData) -> pd.DataFrame:
    """Count distinct target species seen in each FOV."""
    with _timed(f"[{sd.name}] unique_transcripts_per_fov"):
        df = sd.tx.groupby("fov")["target"].nunique().reset_index(name="unique")
        return df.sort_values("fov").reset_index(drop=True)


def panel_detection_ratio_per_fov(sd: SampleData) -> pd.DataFrame:
    """Per-FOV detection ratio: distinct targets seen / panel size."""
    with _timed(f"[{sd.name}] panel_detection_ratio_per_fov"):
        cls = classify_expr_columns(sd.exprmat)
        panel_size = len(cls["gene"]) + len(cls["negative"]) + len(cls["falsecode"])
        df = sd.tx.groupby("fov")["target"].nunique().reset_index(name="n_detected")
        if panel_size > 0:
            df["detected"] = df["n_detected"] / panel_size
        else:
            df["detected"] = 0.0
        df["undetected"] = 1.0 - df["detected"]
        df["percent"] = df["detected"] * 100.0
        df["panel_size"] = panel_size
        return df.sort_values("fov").reset_index(drop=True)


def expression_heatmap_data(sd: SampleData) -> dict[str, pd.DataFrame]:
    """Per-FOV × per-target mean expression from exprMat, split by transcript
    kind (gene / negative / falsecode).

    Each value is the mean exprMat count across all cells of that FOV for that
    target. Despite the function name, this is a *mean expression* metric — it
    has nothing to do with target uniqueness; the "heatmap_data" suffix only
    indicates the consumer (the per-kind heatmap panels)."""
    with _timed(f"[{sd.name}] expression_heatmap_data"):
        cls = classify_expr_columns(sd.exprmat)
        out = {}
        for kind, cols in cls.items():
            if not cols:
                out[kind] = pd.DataFrame()
                continue
            out[kind] = sd.exprmat.groupby("fov")[cols].mean()
        return out


def cell_count_per_fov(sd: SampleData) -> pd.DataFrame:
    """Count cells per FOV from metadata."""
    with _timed(f"[{sd.name}] cell_count_per_fov"):
        df = sd.metadata.groupby("fov").size().reset_index(name="n_cells")
        return df.sort_values("fov").reset_index(drop=True)


def cell_area_per_fov(sd: SampleData) -> pd.DataFrame:
    """Per-cell area distribution by FOV (long-form for boxplot)."""
    with _timed(f"[{sd.name}] cell_area_per_fov"):
        return sd.metadata[["fov", "Area.um2"]].copy()


def error_rates_per_fov(sd: SampleData) -> pd.DataFrame:
    """Per-FOV negative_rate and false_code_rate from tx_file."""
    with _timed(f"[{sd.name}] error_rates_per_fov"):
        grp = pd.crosstab(sd.tx["fov"], sd.tx_kind)
        for k in ("gene", "negative", "falsecode"):
            if k not in grp.columns:
                grp[k] = 0
        out = pd.DataFrame(
            {
                "fov": grp.index,
                "negative_rate": grp["negative"] / grp["gene"].replace(0, np.nan),
                "false_code_rate": grp["falsecode"] / grp["gene"].replace(0, np.nan),
            }
        ).reset_index(drop=True)
        return out


def transcripts_per_cell_per_fov(sd: SampleData, kind: TargetKind = "all") -> pd.DataFrame:
    """Per-cell sum of selected expression columns, with fov/cell_ID."""
    with _timed(f"[{sd.name}] transcripts_per_cell_per_fov(kind={kind})"):
        cls = classify_expr_columns(sd.exprmat)
        if kind == "all":
            cols = cls["gene"] + cls["negative"] + cls["falsecode"]
        else:
            cols = cls[kind]
        out = sd.exprmat[["fov", "cell_ID"]].copy()
        out["count"] = sd.exprmat[cols].sum(axis=1) if cols else 0
        return out


def error_rates_per_cell_per_fov(sd: SampleData) -> pd.DataFrame:
    """Per-cell negative_rate and false_code_rate from exprmat."""
    with _timed(f"[{sd.name}] error_rates_per_cell_per_fov"):
        cls = classify_expr_columns(sd.exprmat)
        gene = sd.exprmat[cls["gene"]].sum(axis=1) if cls["gene"] else 0
        neg = sd.exprmat[cls["negative"]].sum(axis=1) if cls["negative"] else 0
        fc = sd.exprmat[cls["falsecode"]].sum(axis=1) if cls["falsecode"] else 0
        out = sd.exprmat[["fov", "cell_ID"]].copy()
        gene_safe = gene.replace(0, np.nan)
        out["negative_rate"] = neg / gene_safe
        out["false_code_rate"] = fc / gene_safe
        return out


def global_view_table(sd: SampleData) -> pd.DataFrame:
    """One-row-per-slide summary."""
    with _timed(f"[{sd.name}] global_view_table"):
        md = sd.metadata
        rows = []
        for slide, sub in md.groupby("slide_ID"):
            rows.append(
                {
                    "slide_ID": slide,
                    "Run_Tissue_name": sub["Run_Tissue_name"].iloc[0]
                    if "Run_Tissue_name" in sub.columns
                    else None,
                    "n_fovs": sub["fov"].nunique(),
                    "n_cells": len(sub),
                    "n_transcripts": int((sd.tx["fov"].isin(sub["fov"])).sum()),
                }
            )
        return pd.DataFrame(rows)
