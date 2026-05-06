"""IO for CosMx flatfile samples."""

from __future__ import annotations

import gzip
import io
import logging
import time
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path

import pandas as pd
import polars as pl

logger = logging.getLogger("cosmx_qc")

# Columns we actually use from each flatfile. Loading only these:
#   1. cuts memory drastically on the 200M-row tx file (3 of 8 cols);
#   2. avoids polars schema-inference traps on metadata's RNA_quantile_*
#      columns, which are integer-valued in the first ~30k rows of
#      large samples and only later contain floats.
_TX_COLS = ["fov", "cell_ID", "target"]
_METADATA_COLS = ["fov", "cell_ID", "Area.um2", "nFeature_RNA", "slide_ID", "Run_Tissue_name"]


def _read_csv_gz(
    path: Path, columns: list[str] | None = None, n_rows: int | None = None
) -> pd.DataFrame:
    """Fast multi-threaded CSV reader for gzipped files via polars.

    Polars handles gzip auto-detection from the `.gz` extension and
    streams parsing into Arrow buffers; it is dramatically faster and
    more memory-efficient than pandas' C engine on large CosMx flatfiles.
    `columns`, if given, projects to a subset before materialising —
    essential for the 200M-row tx_file.

    `n_rows`, if given, limits the read to the first N rows. Polars
    cannot break out of compressed streams early, so we stream-decompress
    via `gzip.open` and feed only the requested header + N data rows to
    polars — this turns Pembro7 head-slice loads from "decompress all 7 GB"
    into "decompress ~50 MB".

    `infer_schema_length=10000` (vs polars' default of 100) covers the
    common case where CosMx metadata columns have integer-looking values
    in the first dozens of rows but fractional values later (e.g.
    `RNA_quantile_0.99`).
    """
    if n_rows is not None:
        # Stream-decompress just enough to satisfy n_rows + 1 (header).
        with gzip.open(str(path), "rb") as f:
            buf = bytearray()
            buf += f.readline()  # header
            for _ in range(n_rows):
                line = f.readline()
                if not line:
                    break
                buf += line
        df = pl.read_csv(io.BytesIO(bytes(buf)), columns=columns, infer_schema_length=10000)
    else:
        df = pl.read_csv(str(path), columns=columns, infer_schema_length=10000)
    return df.to_pandas()


REQUIRED_GLOBS = [
    "*_exprMat_file.csv.gz",
    "*_metadata_file.csv.gz",
    "*_tx_file.csv.gz",
    "*_fov_positions_file.csv.gz",
]


@dataclass
class SampleData:
    name: str
    exprmat: pd.DataFrame
    metadata: pd.DataFrame
    tx: pd.DataFrame
    fov_pos: pd.DataFrame

    @cached_property
    def tx_kind(self) -> pd.Series:
        """Element-wise gene/negative/falsecode for tx['target'].

        Computed once on first access; reused across all metric functions
        that filter or group by transcript class.
        """
        target = self.tx["target"]
        out = pd.Series("gene", index=target.index, dtype=object)
        out[target.str.startswith("Negative", na=False)] = "negative"
        out[target.str.startswith("SystemControl", na=False)] = "falsecode"
        return out


def _find_one(dir: Path, pattern: str) -> Path:
    matches = sorted(dir.glob(pattern))
    if not matches:
        raise FileNotFoundError(f"No file matching {pattern!r} in {dir}")
    if len(matches) > 1:
        raise RuntimeError(f"Multiple files match {pattern!r} in {dir}: {matches}")
    return matches[0]


def load_sample(name: str, dir: Path, n_rows: int | None = None) -> SampleData:
    """Load the four required flatfiles from a sample dir.

    `dir` may be the flatfiles directory directly, an AtoMx sample
    root (`<sample>/AtoMx/flatFiles/<run>/`), or a parent of either —
    `resolve_flatfiles_dir` finds the right place.

    `n_rows`, if given, applies to the two large files (exprMat and
    tx); metadata and fov_pos are always read in full because they
    are small. Used for fast style iteration on a head slice of a
    Pembro-scale sample.
    """
    target = resolve_flatfiles_dir(Path(dir))

    if n_rows is not None:
        logger.info("[%s] n_rows=%d (head-slice mode)", name, n_rows)

    t_load = time.monotonic()
    t0 = time.monotonic()
    exprmat = _read_csv_gz(_find_one(target, "*_exprMat_file.csv.gz"), n_rows=n_rows)
    logger.info(
        "[%s] exprMat read: %d rows × %d cols in %.1fs", name, *exprmat.shape, time.monotonic() - t0
    )

    t0 = time.monotonic()
    metadata = _read_csv_gz(_find_one(target, "*_metadata_file.csv.gz"), columns=_METADATA_COLS)
    logger.info(
        "[%s] metadata read: %d rows × %d cols in %.1fs",
        name,
        *metadata.shape,
        time.monotonic() - t0,
    )

    t0 = time.monotonic()
    tx = _read_csv_gz(_find_one(target, "*_tx_file.csv.gz"), columns=_TX_COLS, n_rows=n_rows)
    logger.info("[%s] tx read: %d rows × %d cols in %.1fs", name, *tx.shape, time.monotonic() - t0)

    t0 = time.monotonic()
    fov_pos = _read_csv_gz(_find_one(target, "*_fov_positions_file.csv.gz"))
    fov_pos = fov_pos.rename(columns={"x_global_mm": "X_mm", "y_global_mm": "Y_mm"})
    logger.info(
        "[%s] fov_positions read: %d rows × %d cols in %.1fs",
        name,
        *fov_pos.shape,
        time.monotonic() - t0,
    )

    logger.info("[%s] loaded total %.1fs", name, time.monotonic() - t_load)

    return SampleData(name=name, exprmat=exprmat, metadata=metadata, tx=tx, fov_pos=fov_pos)


def _has_all_flatfiles(dir: Path) -> bool:
    """True iff `dir` contains all four required CosMx flatfile globs."""
    return all(list(dir.glob(p)) for p in REQUIRED_GLOBS)


def resolve_flatfiles_dir(dir: Path) -> Path:
    """Find the actual flatfiles directory under `dir`.

    Accepts any of:
    - the flatfiles directory itself (e.g. `Pembro7/RNA_flatfiles/`)
    - the AtoMx-style sample root that contains
      `AtoMx/flatFiles/<run_subdir>/` (e.g. `cosmx_backup/<sample>/`)
    - a parent that has exactly one immediate child satisfying either
      of the above (e.g. `Pembro7/`)

    Raises FileNotFoundError if nothing matches.
    """
    dir = Path(dir)
    if not dir.exists():
        raise FileNotFoundError(f"Sample directory does not exist: {dir}")
    if not dir.is_dir():
        raise NotADirectoryError(f"Not a directory: {dir}")

    logger.info("resolving flatfiles under %s", dir)

    # 1) flatfiles right here
    if _has_all_flatfiles(dir):
        logger.info("  step 1 hit: flatfiles in %s", dir)
        return dir
    logger.info("  step 1 miss: no flatfiles directly in %s", dir)

    # 2) AtoMx/flatFiles/<run_subdir>
    atomx_candidates = sorted(dir.glob("AtoMx/flatFiles/*"))
    for c in atomx_candidates:
        if c.is_dir() and _has_all_flatfiles(c):
            logger.info("  step 2 hit: AtoMx/flatFiles/%s", c.name)
            return c
    if atomx_candidates:
        logger.info(
            "  step 2 miss: %d AtoMx/flatFiles/* candidates, none had all 4 files",
            len(atomx_candidates),
        )
    else:
        logger.info("  step 2 miss: no AtoMx/flatFiles/ directory")

    # 3) one level down (e.g. user passed sample-root for Pembro-style)
    for c in sorted(p for p in dir.iterdir() if p.is_dir()):
        if _has_all_flatfiles(c):
            logger.info("  step 3 hit: child %s", c.name)
            return c
    logger.info("  step 3 miss: no immediate subdir had all 4 files")

    raise FileNotFoundError(
        f"{dir} (and AtoMx/flatFiles/*, immediate subdirs) does not "
        f"contain a directory with all four required flatfiles."
    )


def validate_sample_dir(dir: Path) -> list[str]:
    """Return list of missing required flatfile globs (empty list = OK).

    Resolves the dir first via `resolve_flatfiles_dir`, so callers can
    pass either the flatfiles dir or its sample-root parent.
    """
    dir = Path(dir)
    if not dir.exists():
        raise FileNotFoundError(f"Sample directory does not exist: {dir}")
    if not dir.is_dir():
        raise NotADirectoryError(f"Not a directory: {dir}")
    try:
        target = resolve_flatfiles_dir(dir)
    except FileNotFoundError:
        # Report which globs are missing relative to `dir` itself for
        # the original "user pointed at the flatfiles dir but it's
        # incomplete" use case.
        return [p for p in REQUIRED_GLOBS if not list(dir.glob(p))]
    return [p for p in REQUIRED_GLOBS if not list(target.glob(p))]
