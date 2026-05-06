# cosmx-qc Report Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or executing-plans-test-first to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the `cosmx-qc` Python CLI that consumes one or more CosMx flatfile sample directories and produces a single self-contained HTML QC report mirroring the 13-section R tutorial structure, rendered via Quarto driving Python.

**Architecture:** Flat module layout under `src/cosmx_qc/` with strict separation of `io` (load flatfiles → `SampleData`), `metrics` (pure DataFrame transformations), `plots` (matplotlib figures), `cli` (argparse + validation), `render` (tempfile JSON + subprocess to Quarto), and `report.qmd` (Quarto template). The CLI passes runtime config to Quarto via a JSON tempfile pointed at by the `COSMX_QC_CONFIG` env var; the `.qmd` reads it once and invokes `metrics`/`plots` per section.

**Tech Stack:** Python 3.10+, pandas, numpy, matplotlib, PyYAML, argparse (stdlib), Quarto CLI ≥ 1.4 (system dep), pytest.

**Spec:** `docs/superpowers/specs/2026-05-05-cosmx-qc-report-design.md`

---

## File structure

**Create:**
- `src/cosmx_qc/__init__.py` (rewrite — minimal exports)
- `src/cosmx_qc/io.py`
- `src/cosmx_qc/metrics.py`
- `src/cosmx_qc/plots.py`
- `src/cosmx_qc/cli.py`
- `src/cosmx_qc/render.py`
- `src/cosmx_qc/report.qmd`
- `tests/__init__.py`
- `tests/__init__.py`
- `tests/conftest.py` (synthetic SampleData fixtures)
- `tests/test_io.py`
- `tests/test_metrics.py`
- `tests/test_plots.py`
- `tests/test_cli.py`
- `tests/test_render.py`
- `tests/fixtures/mini_sample/RNA_mini_exprMat_file.csv.gz`
- `tests/fixtures/mini_sample/RNA_mini_metadata_file.csv.gz`
- `tests/fixtures/mini_sample/RNA_mini_tx_file.csv.gz`
- `tests/fixtures/mini_sample/RNA_mini_fov_positions_file.csv.gz`
- `tests/fixtures/build_mini_sample.py` (script that produces the four `.csv.gz` files)

**Modify:**
- `pyproject.toml` (add deps + entry point)
- `README.md` (replace "TODO: QC report" line)

**Delete (after new code is in place and tested):**
- `src/cosmx_qc/cosmx_unified_viz.py`
- `src/cosmx_qc/htmldesign.py`
- `src/cosmx_qc/getfilepath.py`
- `src/cosmx_qc/fov_metrics.py`
- `src/cosmx_qc/cell_metrics.py`
- `src/cosmx_qc/file_paths.csv`
- `src/cosmx_qc/output_plots_unified/`
- `src/cosmx_qc/unified_output_plots/`
- `src/cosmx_qc/unified_output_plots_2/`
- `README_unified.md`
- `examples/cosmx_unified_examples.py`

---

## Task 1: Project setup — deps, entry point, test scaffolding

**Files:**
- Modify: `pyproject.toml`
- Create: `tests/__init__.py`, `tests/__init__.py`

- [ ] **Step 1: Update pyproject.toml**

Replace the `[project]` `dependencies` and `[project.scripts]` blocks:

```toml
[project]
name = "cosmx-utils"
version = "0.1.3"
description = "utilities for CosMx"
authors = [
    {name = "Wenrui Wu", email = "wuwenruiwwr@outlook.com"},
]
readme = "README.md"
requires-python = ">=3.10"
dependencies = [
    "ipykernel>=7.0.1",
    "matplotlib>=3.7",
    "numpy>=1.24",
    "pandas>=2.3.3",
    "paramiko<4",
    "pysftp>=0.2.9",
    "pyyaml>=6.0",
    "tqdm>=4.67.1",
]

[project.optional-dependencies]
dev = ["pytest>=7", "pytest-mpl>=0.16"]

[project.scripts]
atomx-transfer = "cosmx_utils.atomx.transfer:main"
cosmx-qc       = "cosmx_qc.cli:main"
cosmx-utils    = "cosmx_utils.cli:main"
```

- [ ] **Step 2: Create empty test packages**

```bash
mkdir -p tests/fixtures/mini_sample
touch tests/__init__.py tests/__init__.py
```

- [ ] **Step 3: Install dev deps and verify pytest runs**

```bash
uv pip install -e '.[dev]'
pytest --collect-only
```

Expected: pytest exits 0, "no tests ran" or similar (no test files yet).

- [ ] **Step 4: Commit**

```bash
git add pyproject.toml tests/__init__.py tests/__init__.py
git commit -m "chore: add matplotlib/pyyaml deps and qc-cli entry point"
```

---

## Task 2: Build mini sample test fixture

**Files:**
- Create: `tests/fixtures/build_mini_sample.py`
- Create: `tests/fixtures/mini_sample/*.csv.gz` (generated)

The fixture is a tiny synthetic CosMx flatfile sample: 2 FOVs, 5 cells (3 in fov 1, 2 in fov 2), 4 expression columns (1 real gene `GeneA`, 1 negative `Negative1`, 1 false code `SystemControl1`, 1 metadata col already excluded).

- [ ] **Step 1: Write the fixture generator**

Create `tests/fixtures/build_mini_sample.py`:

```python
"""Generate a tiny synthetic CosMx flatfile sample for tests.

Layout: 2 FOVs, 5 cells, 1 real gene + 1 negative probe + 1 false code.
Run from repo root:  python tests/fixtures/build_mini_sample.py
"""
from pathlib import Path
import pandas as pd

OUT = Path(__file__).parent / "mini_sample"
OUT.mkdir(parents=True, exist_ok=True)

# 5 cells across 2 FOVs
exprmat = pd.DataFrame({
    "fov":          [1, 1, 1, 2, 2],
    "cell_ID":      [1, 2, 3, 1, 2],
    "GeneA":        [10, 20, 30, 5,  15],
    "Negative1":    [1,  0,  2,  0,  1],
    "SystemControl1":[0, 1,  0,  1,  0],
})

metadata = pd.DataFrame({
    "fov":          [1, 1, 1, 2, 2],
    "cell_ID":      [1, 2, 3, 1, 2],
    "Area.um2":     [100.0, 200.0, 150.0, 80.0, 120.0],
    "nFeature_RNA": [3, 2, 3, 2, 2],
    "slide_ID":     ["slide_A"] * 5,
    "Run_Tissue_name": ["mini"] * 5,
})

# 13 transcripts, some unassigned (cell_ID = 0)
tx = pd.DataFrame({
    "fov":     [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 1, 2, 2],
    "cell_ID": [1, 1, 2, 3, 0, 1, 1, 2, 2, 0, 0, 0, 0],
    "target":  ["GeneA","Negative1","GeneA","SystemControl1","GeneA",
                "GeneA","GeneA","Negative1","GeneA","GeneA",
                "GeneA","Negative1","SystemControl1"],
    "CellComp":["Cytoplasm"] * 13,
    "x_local_px":[100,110,120,130,140, 200,210,220,230,240, 150,250,260],
    "y_local_px":[100,110,120,130,140, 200,210,220,230,240, 150,250,260],
})

fov_pos = pd.DataFrame({
    "FOV":          [1, 2],
    "x_global_mm":  [0.0, 5.0],
    "y_global_mm":  [0.0, 0.0],
})

prefix = "RNA_mini"
exprmat.to_csv(OUT / f"{prefix}_exprMat_file.csv.gz", index=False, compression="gzip")
metadata.to_csv(OUT / f"{prefix}_metadata_file.csv.gz", index=False, compression="gzip")
tx.to_csv(OUT / f"{prefix}_tx_file.csv.gz", index=False, compression="gzip")
fov_pos.to_csv(OUT / f"{prefix}_fov_positions_file.csv.gz", index=False, compression="gzip")
print(f"Wrote 4 files to {OUT}")
```

- [ ] **Step 2: Run it to produce the fixture files**

```bash
python tests/fixtures/build_mini_sample.py
ls tests/fixtures/mini_sample/
```

Expected:
```
RNA_mini_exprMat_file.csv.gz
RNA_mini_fov_positions_file.csv.gz
RNA_mini_metadata_file.csv.gz
RNA_mini_tx_file.csv.gz
```

- [ ] **Step 3: Commit (include the generated `.csv.gz` files for reproducibility)**

```bash
git add tests/fixtures/
git commit -m "test: add mini synthetic CosMx flatfile fixture"
```

---

## Task 3: io.SampleData + validate_sample_dir

**Files:**
- Create: `src/cosmx_qc/io.py`
- Create: `tests/test_io.py`

- [ ] **Step 1: Write the failing test**

Create `tests/test_io.py`:

```python
from pathlib import Path
import pytest
from cosmx_qc import io as qc_io

FIXTURE = Path(__file__).parent / "fixtures" / "mini_sample"


def test_validate_sample_dir_returns_empty_for_complete_fixture():
    missing = qc_io.validate_sample_dir(FIXTURE)
    assert missing == []


def test_validate_sample_dir_lists_missing(tmp_path):
    # only drop one of the four; validate flags it
    (tmp_path / "RNA_x_exprMat_file.csv.gz").write_bytes(b"")
    (tmp_path / "RNA_x_metadata_file.csv.gz").write_bytes(b"")
    (tmp_path / "RNA_x_tx_file.csv.gz").write_bytes(b"")
    missing = qc_io.validate_sample_dir(tmp_path)
    assert "*_fov_positions_file.csv.gz" in missing


def test_validate_sample_dir_rejects_nonexistent_dir(tmp_path):
    bogus = tmp_path / "does_not_exist"
    with pytest.raises(FileNotFoundError):
        qc_io.validate_sample_dir(bogus)
```

- [ ] **Step 2: Run test to verify it fails**

```bash
pytest tests/test_io.py -v
```

Expected: ImportError or "module has no attribute" (qc_io.validate_sample_dir undefined).

- [ ] **Step 3: Implement validate_sample_dir + SampleData skeleton**

Create `src/cosmx_qc/io.py`:

```python
"""IO for CosMx flatfile samples."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import pandas as pd

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


def validate_sample_dir(dir: Path) -> list[str]:
    """Return list of missing required flatfile globs (empty list = OK)."""
    dir = Path(dir)
    if not dir.exists():
        raise FileNotFoundError(f"Sample directory does not exist: {dir}")
    if not dir.is_dir():
        raise NotADirectoryError(f"Not a directory: {dir}")
    missing = []
    for pattern in REQUIRED_GLOBS:
        if not list(dir.glob(pattern)):
            missing.append(pattern)
    return missing
```

- [ ] **Step 4: Run tests, expect pass**

```bash
pytest tests/test_io.py -v
```

Expected: 3 passed.

- [ ] **Step 5: Commit**

```bash
git add src/cosmx_qc/io.py tests/test_io.py
git commit -m "feat(qc): add SampleData dataclass and validate_sample_dir"
```

---

## Task 4: io.load_sample

**Files:**
- Modify: `src/cosmx_qc/io.py`
- Modify: `tests/test_io.py`

- [ ] **Step 1: Add failing test for load_sample**

Append to `tests/test_io.py`:

```python
def test_load_sample_returns_complete_data():
    sd = qc_io.load_sample("mini", FIXTURE)
    assert sd.name == "mini"
    assert {"fov", "cell_ID", "GeneA", "Negative1", "SystemControl1"}.issubset(sd.exprmat.columns)
    assert {"fov", "cell_ID", "Area.um2"}.issubset(sd.metadata.columns)
    assert {"fov", "cell_ID", "target", "CellComp"}.issubset(sd.tx.columns)
    assert {"FOV", "X_mm", "Y_mm"}.issubset(sd.fov_pos.columns)
    assert len(sd.exprmat) == 5
    assert len(sd.tx) == 13


def test_load_sample_normalizes_position_columns():
    sd = qc_io.load_sample("mini", FIXTURE)
    # x_global_mm/y_global_mm should be renamed to X_mm/Y_mm
    assert "x_global_mm" not in sd.fov_pos.columns
    assert "X_mm" in sd.fov_pos.columns
    assert "Y_mm" in sd.fov_pos.columns
```

- [ ] **Step 2: Run, confirm failures**

```bash
pytest tests/test_io.py -v
```

Expected: the two new tests fail (load_sample undefined).

- [ ] **Step 3: Implement load_sample**

Add to `src/cosmx_qc/io.py`:

```python
def _find_one(dir: Path, pattern: str) -> Path:
    matches = sorted(dir.glob(pattern))
    if not matches:
        raise FileNotFoundError(f"No file matching {pattern!r} in {dir}")
    if len(matches) > 1:
        raise RuntimeError(f"Multiple files match {pattern!r} in {dir}: {matches}")
    return matches[0]


def load_sample(name: str, dir: Path) -> SampleData:
    """Load the four required flatfiles from a sample dir."""
    dir = Path(dir)
    missing = validate_sample_dir(dir)
    if missing:
        raise FileNotFoundError(f"{dir} missing required files: {missing}")

    exprmat = pd.read_csv(_find_one(dir, "*_exprMat_file.csv.gz"))
    metadata = pd.read_csv(_find_one(dir, "*_metadata_file.csv.gz"))
    tx = pd.read_csv(_find_one(dir, "*_tx_file.csv.gz"))
    fov_pos = pd.read_csv(_find_one(dir, "*_fov_positions_file.csv.gz"))

    fov_pos = fov_pos.rename(columns={"x_global_mm": "X_mm", "y_global_mm": "Y_mm"})

    return SampleData(name=name, exprmat=exprmat, metadata=metadata, tx=tx, fov_pos=fov_pos)
```

- [ ] **Step 4: Run, expect pass**

```bash
pytest tests/test_io.py -v
```

Expected: 5 passed.

- [ ] **Step 5: Commit**

```bash
git add src/cosmx_qc/io.py tests/test_io.py
git commit -m "feat(qc): add load_sample to read CosMx flatfiles into SampleData"
```

---

## Task 5: metrics — target classification helpers + per-FOV transcript metrics

**Files:**
- Create: `src/cosmx_qc/metrics.py`
- Create: `tests/conftest.py`
- Create: `tests/test_metrics.py`

- [ ] **Step 1: Add a shared in-memory SampleData fixture**

Create `tests/conftest.py`:

```python
"""Shared pytest fixtures for qc tests."""
import pandas as pd
import pytest
from cosmx_qc.io import SampleData


@pytest.fixture
def mini_sample() -> SampleData:
    """Hand-built tiny SampleData mirroring the file fixture."""
    exprmat = pd.DataFrame({
        "fov":          [1, 1, 1, 2, 2],
        "cell_ID":      [1, 2, 3, 1, 2],
        "GeneA":        [10, 20, 30, 5,  15],
        "Negative1":    [1,  0,  2,  0,  1],
        "SystemControl1":[0, 1,  0,  1,  0],
    })
    metadata = pd.DataFrame({
        "fov":          [1, 1, 1, 2, 2],
        "cell_ID":      [1, 2, 3, 1, 2],
        "Area.um2":     [100.0, 200.0, 150.0, 80.0, 120.0],
        "nFeature_RNA": [3, 2, 3, 2, 2],
        "slide_ID":     ["slide_A"] * 5,
        "Run_Tissue_name": ["mini"] * 5,
    })
    tx = pd.DataFrame({
        "fov":     [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 1, 2, 2],
        "cell_ID": [1, 1, 2, 3, 0, 1, 1, 2, 2, 0, 0, 0, 0],
        "target":  ["GeneA","Negative1","GeneA","SystemControl1","GeneA",
                    "GeneA","GeneA","Negative1","GeneA","GeneA",
                    "GeneA","Negative1","SystemControl1"],
        "CellComp":["Cytoplasm"] * 13,
        "x_local_px":[100,110,120,130,140, 200,210,220,230,240, 150,250,260],
        "y_local_px":[100,110,120,130,140, 200,210,220,230,240, 150,250,260],
    })
    fov_pos = pd.DataFrame({"FOV":[1,2], "X_mm":[0.0,5.0], "Y_mm":[0.0,0.0]})
    return SampleData(name="mini", exprmat=exprmat, metadata=metadata, tx=tx, fov_pos=fov_pos)
```

- [ ] **Step 2: Write failing tests for classification helpers and per-FOV transcript metrics**

Create `tests/test_metrics.py`:

```python
import pandas as pd
from cosmx_qc import metrics as M


def test_classify_columns(mini_sample):
    cls = M.classify_expr_columns(mini_sample.exprmat)
    assert cls["gene"] == ["GeneA"]
    assert cls["negative"] == ["Negative1"]
    assert cls["falsecode"] == ["SystemControl1"]


def test_classify_targets(mini_sample):
    s = mini_sample.tx["target"]
    assert M.classify_target(s).tolist() == [
        "gene","negative","gene","falsecode","gene",
        "gene","gene","negative","gene","gene",
        "gene","negative","falsecode",
    ]


def test_transcripts_per_fov_all(mini_sample):
    df = M.transcripts_per_fov(mini_sample, kind="all")
    # 13 transcripts total: 6 in fov 1, 7 in fov 2
    expected = pd.DataFrame({"fov":[1,2], "count":[6,7]})
    pd.testing.assert_frame_equal(df.reset_index(drop=True), expected)


def test_transcripts_per_fov_gene(mini_sample):
    df = M.transcripts_per_fov(mini_sample, kind="gene")
    # gene tx in fov 1: 4 (rows with target=GeneA at indices 0,2,4,10)
    # gene tx in fov 2: 4 (indices 5,6,8,9)
    expected = pd.DataFrame({"fov":[1,2], "count":[4,4]})
    pd.testing.assert_frame_equal(df.reset_index(drop=True), expected)


def test_transcripts_per_fov_negative(mini_sample):
    df = M.transcripts_per_fov(mini_sample, kind="negative")
    # negative tx in fov 1: 1 (index 1); fov 2: 2 (indices 7, 11)
    expected = pd.DataFrame({"fov":[1,2], "count":[1,2]})
    pd.testing.assert_frame_equal(df.reset_index(drop=True), expected)


def test_transcripts_per_fov_falsecode(mini_sample):
    df = M.transcripts_per_fov(mini_sample, kind="falsecode")
    # falsecode (SystemControl1) in fov 1: 1 (index 3); fov 2: 1 (index 12)
    expected = pd.DataFrame({"fov":[1,2], "count":[1,1]})
    pd.testing.assert_frame_equal(df.reset_index(drop=True), expected)
```

- [ ] **Step 3: Run, expect failures (module missing)**

```bash
pytest tests/test_metrics.py -v
```

- [ ] **Step 4: Implement metrics.py — classifiers + transcripts_per_fov**

Create `src/cosmx_qc/metrics.py`:

```python
"""Pure-data QC metric functions. Each takes a SampleData (or dict thereof)
and returns a DataFrame.  No plotting, no IO."""
from __future__ import annotations

from typing import Literal

import numpy as np
import pandas as pd

from cosmx_qc.io import SampleData

TargetKind = Literal["all", "gene", "negative", "falsecode"]
META_COLS = {"fov", "cell_ID", "cell"}  # known non-expression columns in exprMat


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
    """Element-wise classification of a tx-file `target` series."""
    out = pd.Series("gene", index=target.index, dtype=object)
    out[target.str.startswith("Negative", na=False)] = "negative"
    out[target.str.startswith("SystemControl", na=False)] = "falsecode"
    return out


def transcripts_per_fov(sd: SampleData, kind: TargetKind = "all") -> pd.DataFrame:
    """Count transcripts per FOV from tx_file, optionally filtered by kind."""
    df = sd.tx
    if kind != "all":
        df = df[classify_target(df["target"]) == kind]
    counts = df.groupby("fov").size().reset_index(name="count")
    return counts.sort_values("fov").reset_index(drop=True)
```

- [ ] **Step 5: Run, expect pass**

```bash
pytest tests/test_metrics.py -v
```

Expected: 6 passed.

- [ ] **Step 6: Commit**

```bash
git add src/cosmx_qc/metrics.py tests/conftest.py tests/test_metrics.py
git commit -m "feat(qc): add target classifiers and per-FOV transcript counts"
```

---

## Task 6: metrics — unassigned, unique, cell counts, cell area

**Files:**
- Modify: `src/cosmx_qc/metrics.py`
- Modify: `tests/test_metrics.py`

- [ ] **Step 1: Add failing tests**

Append to `tests/test_metrics.py`:

```python
def test_unassigned_transcripts_per_fov(mini_sample):
    df = M.unassigned_transcripts_per_fov(mini_sample)
    # cell_ID == 0 transcripts: 2 in fov 1 (indices 4, 10); 3 in fov 2 (indices 9, 11, 12)
    expected = pd.DataFrame({"fov":[1,2], "unassigned":[2,3]})
    pd.testing.assert_frame_equal(df.reset_index(drop=True), expected)


def test_unique_transcripts_per_fov(mini_sample):
    df = M.unique_transcripts_per_fov(mini_sample)
    # fov 1 sees {GeneA, Negative1, SystemControl1} = 3
    # fov 2 sees {GeneA, Negative1, SystemControl1} = 3
    expected = pd.DataFrame({"fov":[1,2], "unique":[3,3]})
    pd.testing.assert_frame_equal(df.reset_index(drop=True), expected)


def test_cell_count_per_fov(mini_sample):
    df = M.cell_count_per_fov(mini_sample)
    expected = pd.DataFrame({"fov":[1,2], "n_cells":[3,2]})
    pd.testing.assert_frame_equal(df.reset_index(drop=True), expected)


def test_cell_area_per_fov(mini_sample):
    df = M.cell_area_per_fov(mini_sample)
    # all rows passed through, with fov + Area.um2
    assert set(df.columns) == {"fov", "Area.um2"}
    assert len(df) == 5
    assert df.loc[df["fov"]==1, "Area.um2"].mean() == pytest.approx(150.0)
```

Also at the top of the test file add:
```python
import pytest
```

- [ ] **Step 2: Run, expect failures**

- [ ] **Step 3: Implement the four metrics**

Append to `src/cosmx_qc/metrics.py`:

```python
def unassigned_transcripts_per_fov(sd: SampleData) -> pd.DataFrame:
    """Count transcripts with cell_ID == 0 per FOV."""
    df = sd.tx[sd.tx["cell_ID"] == 0]
    counts = df.groupby("fov").size().reset_index(name="unassigned")
    # ensure all FOVs present (zero-fill missing)
    all_fovs = pd.DataFrame({"fov": sorted(sd.tx["fov"].unique())})
    out = all_fovs.merge(counts, on="fov", how="left").fillna({"unassigned": 0})
    out["unassigned"] = out["unassigned"].astype(int)
    return out


def unique_transcripts_per_fov(sd: SampleData) -> pd.DataFrame:
    """Count distinct target species seen in each FOV."""
    df = sd.tx.groupby("fov")["target"].nunique().reset_index(name="unique")
    return df.sort_values("fov").reset_index(drop=True)


def cell_count_per_fov(sd: SampleData) -> pd.DataFrame:
    """Count cells per FOV from metadata."""
    df = sd.metadata.groupby("fov").size().reset_index(name="n_cells")
    return df.sort_values("fov").reset_index(drop=True)


def cell_area_per_fov(sd: SampleData) -> pd.DataFrame:
    """Per-cell area distribution by FOV (long-form for boxplot)."""
    return sd.metadata[["fov", "Area.um2"]].copy()
```

- [ ] **Step 4: Run, expect pass**

```bash
pytest tests/test_metrics.py -v
```

Expected: 10 passed.

- [ ] **Step 5: Commit**

```bash
git add src/cosmx_qc/metrics.py tests/test_metrics.py
git commit -m "feat(qc): add unassigned, unique, cell_count, cell_area metrics"
```

---

## Task 7: metrics — error rates per FOV + global summary table

**Files:**
- Modify: `src/cosmx_qc/metrics.py`
- Modify: `tests/test_metrics.py`

- [ ] **Step 1: Add failing tests**

```python
def test_error_rates_per_fov(mini_sample):
    df = M.error_rates_per_fov(mini_sample)
    # fov 1 tx kinds: gene=4, negative=1, falsecode=1
    #   neg_rate = 1/4 = 0.25 ; false_code_rate = 1/4 = 0.25
    # fov 2 tx kinds: gene=4, negative=2, falsecode=1
    #   neg_rate = 2/4 = 0.5  ; false_code_rate = 1/4 = 0.25
    assert df.loc[df["fov"]==1, "negative_rate"].iloc[0] == pytest.approx(0.25)
    assert df.loc[df["fov"]==1, "false_code_rate"].iloc[0] == pytest.approx(0.25)
    assert df.loc[df["fov"]==2, "negative_rate"].iloc[0] == pytest.approx(0.5)
    assert df.loc[df["fov"]==2, "false_code_rate"].iloc[0] == pytest.approx(0.25)


def test_global_view_table(mini_sample):
    df = M.global_view_table(mini_sample)
    # one row per slide_ID
    assert len(df) == 1
    row = df.iloc[0]
    assert row["slide_ID"] == "slide_A"
    assert row["n_fovs"] == 2
    assert row["n_cells"] == 5
    assert row["n_transcripts"] == 13
```

- [ ] **Step 2: Run, expect failures**

- [ ] **Step 3: Implement**

```python
def error_rates_per_fov(sd: SampleData) -> pd.DataFrame:
    """Per-FOV negative_rate and false_code_rate from tx_file."""
    df = sd.tx.copy()
    df["kind"] = classify_target(df["target"])
    grp = df.groupby(["fov", "kind"]).size().unstack(fill_value=0)
    for k in ("gene", "negative", "falsecode"):
        if k not in grp.columns:
            grp[k] = 0
    out = pd.DataFrame({
        "fov": grp.index,
        "negative_rate": grp["negative"] / grp["gene"].replace(0, np.nan),
        "false_code_rate": grp["falsecode"] / grp["gene"].replace(0, np.nan),
    }).reset_index(drop=True)
    return out


def global_view_table(sd: SampleData) -> pd.DataFrame:
    """One-row-per-slide summary."""
    md = sd.metadata
    rows = []
    for slide, sub in md.groupby("slide_ID"):
        rows.append({
            "slide_ID": slide,
            "Run_Tissue_name": sub["Run_Tissue_name"].iloc[0]
                if "Run_Tissue_name" in sub.columns else None,
            "n_fovs": sub["fov"].nunique(),
            "n_cells": len(sub),
            "n_transcripts": int((sd.tx["fov"].isin(sub["fov"])).sum()),
        })
    return pd.DataFrame(rows)
```

- [ ] **Step 4: Run, expect pass**

```bash
pytest tests/test_metrics.py -v
```

Expected: 12 passed.

- [ ] **Step 5: Commit**

```bash
git add src/cosmx_qc/metrics.py tests/test_metrics.py
git commit -m "feat(qc): add error rates per FOV and global summary table"
```

---

## Task 8: metrics — per-cell-per-FOV transcript counts and error rates

**Files:**
- Modify: `src/cosmx_qc/metrics.py`
- Modify: `tests/test_metrics.py`

These metrics are computed from `exprmat` (per cell), not `tx` (per transcript).

- [ ] **Step 1: Add failing tests**

```python
def test_transcripts_per_cell_per_fov_all(mini_sample):
    df = M.transcripts_per_cell_per_fov(mini_sample, kind="all")
    # cells: per-cell sum of all expression cols
    # fov 1, cell 1: 10+1+0=11
    # fov 1, cell 2: 20+0+1=21
    # fov 1, cell 3: 30+2+0=32
    # fov 2, cell 1: 5+0+1=6
    # fov 2, cell 2: 15+1+0=16
    assert set(df.columns) >= {"fov", "cell_ID", "count"}
    assert df.loc[(df["fov"]==1) & (df["cell_ID"]==3), "count"].iloc[0] == 32


def test_transcripts_per_cell_per_fov_gene(mini_sample):
    df = M.transcripts_per_cell_per_fov(mini_sample, kind="gene")
    # only GeneA values
    expected_counts = {(1,1):10,(1,2):20,(1,3):30,(2,1):5,(2,2):15}
    for (fov, cid), v in expected_counts.items():
        assert df.loc[(df["fov"]==fov) & (df["cell_ID"]==cid), "count"].iloc[0] == v


def test_error_rates_per_cell_per_fov(mini_sample):
    df = M.error_rates_per_cell_per_fov(mini_sample)
    # per cell: neg_rate = neg_count / gene_count, false_code_rate = sc / gene
    assert set(df.columns) >= {"fov","cell_ID","negative_rate","false_code_rate"}
    # fov 1 cell 1: GeneA=10, Neg=1, SC=0 -> 0.1, 0.0
    row = df[(df["fov"]==1) & (df["cell_ID"]==1)].iloc[0]
    assert row["negative_rate"] == pytest.approx(0.1)
    assert row["false_code_rate"] == pytest.approx(0.0)
```

- [ ] **Step 2: Run, expect failures**

- [ ] **Step 3: Implement**

```python
def transcripts_per_cell_per_fov(sd: SampleData, kind: TargetKind = "all") -> pd.DataFrame:
    """Per-cell sum of selected expression columns, with fov/cell_ID."""
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
    cls = classify_expr_columns(sd.exprmat)
    gene = sd.exprmat[cls["gene"]].sum(axis=1) if cls["gene"] else 0
    neg = sd.exprmat[cls["negative"]].sum(axis=1) if cls["negative"] else 0
    fc = sd.exprmat[cls["falsecode"]].sum(axis=1) if cls["falsecode"] else 0
    out = sd.exprmat[["fov", "cell_ID"]].copy()
    gene_safe = gene.replace(0, np.nan)
    out["negative_rate"] = neg / gene_safe
    out["false_code_rate"] = fc / gene_safe
    return out
```

- [ ] **Step 4: Run, expect pass**

Expected: 15 passed.

- [ ] **Step 5: Commit**

```bash
git add src/cosmx_qc/metrics.py tests/test_metrics.py
git commit -m "feat(qc): add per-cell-per-FOV transcript and error-rate metrics"
```

---

## Task 9: plots — style constants + combined-layout helper

**Files:**
- Create: `src/cosmx_qc/plots.py`
- Create: `tests/test_plots.py`

- [ ] **Step 1: Failing test for combine_axes**

Create `tests/test_plots.py`:

```python
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pytest
from cosmx_qc import plots as P
from cosmx_qc import metrics as M


def test_make_combined_axes_returns_n_axes():
    fig, axes = P.make_combined_axes(3)
    assert len(axes) == 3
    plt.close(fig)


def test_make_combined_axes_single():
    fig, axes = P.make_combined_axes(1)
    assert len(axes) == 1
    plt.close(fig)
```

- [ ] **Step 2: Run, expect failures**

- [ ] **Step 3: Implement plots.py skeleton**

Create `src/cosmx_qc/plots.py`:

```python
"""Matplotlib plotting for CosMx QC. Each metric has a `_for_sample`
helper that draws into a caller-provided ax, plus a `_combined` helper
that arranges N samples horizontally."""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.figure import Figure

_STYLE = {
    "palette": plt.get_cmap("viridis"),
    "fig_height_in": 4.0,
    "panel_width_in": 4.5,
    "title_fontsize": 11,
    "label_fontsize": 9,
    "boxplot_kwargs": dict(showfliers=False, patch_artist=True),
}


def make_combined_axes(n: int) -> tuple[Figure, list[Axes]]:
    """Return (fig, [ax,...]) with N horizontal panels and constrained_layout."""
    width = max(1, n) * _STYLE["panel_width_in"]
    fig, axes = plt.subplots(
        1, n, figsize=(width, _STYLE["fig_height_in"]),
        constrained_layout=True, squeeze=False,
    )
    return fig, list(axes[0])
```

- [ ] **Step 4: Run, expect pass (2 tests)**

- [ ] **Step 5: Commit**

```bash
git add src/cosmx_qc/plots.py tests/test_plots.py
git commit -m "feat(qc): add plots module with style constants and combined-axes helper"
```

---

## Task 10: plots — per-FOV chart functions (bar/box, FOV positions)

**Files:**
- Modify: `src/cosmx_qc/plots.py`
- Modify: `tests/test_plots.py`

Smoke tests only: confirm each plot returns or fills an ax without raising.

- [ ] **Step 1: Add failing smoke tests**

```python
@pytest.fixture
def mini_for_plot(mini_sample):
    return mini_sample


def test_plot_fov_positions_for_sample(mini_for_plot):
    fig, axes = P.make_combined_axes(1)
    P.plot_fov_positions_for_sample(mini_for_plot.fov_pos, axes[0], label="mini")
    plt.close(fig)


def test_plot_fov_positions_combined(mini_for_plot):
    fig = P.plot_fov_positions_combined({"mini": mini_for_plot.fov_pos})
    assert isinstance(fig, plt.Figure)
    plt.close(fig)


def test_plot_count_per_fov_for_sample(mini_for_plot):
    df = M.transcripts_per_fov(mini_for_plot, kind="all")
    fig, axes = P.make_combined_axes(1)
    P.plot_count_per_fov_for_sample(df, axes[0], y_col="count", ylabel="N tx")
    plt.close(fig)


def test_plot_count_per_fov_combined(mini_for_plot):
    df = M.transcripts_per_fov(mini_for_plot, kind="all")
    fig = P.plot_count_per_fov_combined({"mini": df}, y_col="count", ylabel="N tx")
    assert isinstance(fig, plt.Figure)
    plt.close(fig)


def test_plot_box_per_fov_combined(mini_for_plot):
    df = M.cell_area_per_fov(mini_for_plot)
    fig = P.plot_box_per_fov_combined({"mini": df}, y_col="Area.um2", ylabel="Area (um2)")
    assert isinstance(fig, plt.Figure)
    plt.close(fig)
```

- [ ] **Step 2: Run, expect failures**

- [ ] **Step 3: Implement**

Append to `src/cosmx_qc/plots.py`:

```python
def plot_fov_positions_for_sample(fov_pos: pd.DataFrame, ax: Axes,
                                   label: str = "") -> None:
    """Scatter FOV index on (X_mm, Y_mm). Annotate each point with its FOV id."""
    ax.scatter(fov_pos["X_mm"], fov_pos["Y_mm"], s=120, alpha=0.6, color="#2c7fb8")
    for _, r in fov_pos.iterrows():
        ax.annotate(int(r["FOV"]), (r["X_mm"], r["Y_mm"]),
                    ha="center", va="center", fontsize=_STYLE["label_fontsize"])
    ax.set_xlabel("X (mm)"); ax.set_ylabel("Y (mm)")
    ax.set_aspect("equal", adjustable="datalim")
    if label:
        ax.set_title(label, fontsize=_STYLE["title_fontsize"])


def plot_fov_positions_combined(by_sample: dict[str, pd.DataFrame]) -> Figure:
    fig, axes = make_combined_axes(len(by_sample))
    for ax, (name, df) in zip(axes, by_sample.items()):
        plot_fov_positions_for_sample(df, ax, label=name)
    return fig


def plot_count_per_fov_for_sample(df: pd.DataFrame, ax: Axes,
                                   y_col: str, ylabel: str,
                                   label: str = "") -> None:
    """Bar chart: x=fov, y=y_col."""
    ax.bar(df["fov"].astype(str), df[y_col], color="#2c7fb8")
    ax.set_xlabel("FOV"); ax.set_ylabel(ylabel)
    if label:
        ax.set_title(label, fontsize=_STYLE["title_fontsize"])
    for tick in ax.get_xticklabels():
        tick.set_rotation(45)


def plot_count_per_fov_combined(by_sample: dict[str, pd.DataFrame],
                                 y_col: str, ylabel: str) -> Figure:
    fig, axes = make_combined_axes(len(by_sample))
    for ax, (name, df) in zip(axes, by_sample.items()):
        plot_count_per_fov_for_sample(df, ax, y_col=y_col, ylabel=ylabel, label=name)
    return fig


def plot_box_per_fov_for_sample(df: pd.DataFrame, ax: Axes,
                                 y_col: str, ylabel: str,
                                 label: str = "") -> None:
    """Boxplot of `y_col` grouped by fov. `df` is long-form (one row per cell)."""
    fovs = sorted(df["fov"].unique())
    data = [df.loc[df["fov"] == f, y_col].dropna().values for f in fovs]
    ax.boxplot(data, tick_labels=[str(f) for f in fovs], **_STYLE["boxplot_kwargs"])
    ax.set_xlabel("FOV"); ax.set_ylabel(ylabel)
    if label:
        ax.set_title(label, fontsize=_STYLE["title_fontsize"])


def plot_box_per_fov_combined(by_sample: dict[str, pd.DataFrame],
                               y_col: str, ylabel: str) -> Figure:
    fig, axes = make_combined_axes(len(by_sample))
    for ax, (name, df) in zip(axes, by_sample.items()):
        plot_box_per_fov_for_sample(df, ax, y_col=y_col, ylabel=ylabel, label=name)
    return fig
```

- [ ] **Step 4: Run, expect pass**

```bash
pytest tests/test_plots.py -v
```

Expected: 7 passed.

- [ ] **Step 5: Commit**

```bash
git add src/cosmx_qc/plots.py tests/test_plots.py
git commit -m "feat(qc): add per-FOV bar/box and FOV-positions plotting helpers"
```

---

## Task 11: cli — parse_sample_arg

**Files:**
- Create: `src/cosmx_qc/cli.py`
- Create: `tests/test_cli.py`

- [ ] **Step 1: Failing test**

Create `tests/test_cli.py`:

```python
from pathlib import Path
import pytest
from cosmx_qc import cli as C


def test_parse_sample_arg_ok():
    name, p = C.parse_sample_arg("Pembro7=/some/path")
    assert name == "Pembro7"
    assert p == Path("/some/path")


def test_parse_sample_arg_missing_equals():
    with pytest.raises(ValueError, match="must be NAME=PATH"):
        C.parse_sample_arg("Pembro7")


def test_parse_sample_arg_empty_name():
    with pytest.raises(ValueError, match="empty name"):
        C.parse_sample_arg("=/path")


def test_parse_sample_arg_path_with_equals():
    # path may legitimately contain '='; we split on the first
    name, p = C.parse_sample_arg("S1=/a/b=c/d")
    assert name == "S1"
    assert p == Path("/a/b=c/d")
```

- [ ] **Step 2: Run, expect failures**

- [ ] **Step 3: Implement**

Create `src/cosmx_qc/cli.py`:

```python
"""cosmx-qc command-line interface."""
from __future__ import annotations

import argparse
from pathlib import Path


def parse_sample_arg(s: str) -> tuple[str, Path]:
    """'NAME=PATH' → (name, Path). Raises ValueError on malformed input."""
    if "=" not in s:
        raise ValueError(f"--sample must be NAME=PATH, got {s!r}")
    name, _, path = s.partition("=")
    if not name:
        raise ValueError(f"--sample has empty name in {s!r}")
    if not path:
        raise ValueError(f"--sample has empty path in {s!r}")
    return name, Path(path)


def main():  # pragma: no cover  (filled in in later tasks)
    raise NotImplementedError
```

- [ ] **Step 4: Run, expect pass**

```bash
pytest tests/test_cli.py -v
```

Expected: 4 passed.

- [ ] **Step 5: Commit**

```bash
git add src/cosmx_qc/cli.py tests/test_cli.py
git commit -m "feat(qc): add CLI sample-arg parser"
```

---

## Task 12: cli — argparse, validation, config-file loader

**Files:**
- Modify: `src/cosmx_qc/cli.py`
- Modify: `tests/test_cli.py`

- [ ] **Step 1: Failing tests**

Append to `tests/test_cli.py`:

```python
import yaml


def test_load_config_file(tmp_path):
    cfg = tmp_path / "samples.yaml"
    cfg.write_text(yaml.safe_dump({
        "samples": {"A": "/p1", "B": "/p2"},
        "title": "X",
    }))
    samples, title = C.load_config(cfg)
    assert samples == {"A": Path("/p1"), "B": Path("/p2")}
    assert title == "X"


def test_load_config_missing_samples_key(tmp_path):
    cfg = tmp_path / "bad.yaml"
    cfg.write_text("title: only\n")
    with pytest.raises(ValueError, match="samples"):
        C.load_config(cfg)


def test_validate_samples_ok(tmp_path):
    """Sample dirs that contain all 4 globs pass."""
    d = tmp_path / "s1"
    d.mkdir()
    for pat in ["x_exprMat_file.csv.gz", "x_metadata_file.csv.gz",
                "x_tx_file.csv.gz", "x_fov_positions_file.csv.gz"]:
        (d / pat).write_bytes(b"")
    # should not raise
    C.validate_samples({"s1": d})


def test_validate_samples_missing_files(tmp_path):
    d = tmp_path / "s1"; d.mkdir()
    (d / "x_exprMat_file.csv.gz").write_bytes(b"")  # only one of four
    with pytest.raises(SystemExit):
        C.validate_samples({"s1": d})


def test_validate_samples_dir_does_not_exist(tmp_path):
    with pytest.raises(SystemExit):
        C.validate_samples({"s1": tmp_path / "nope"})


def test_check_quarto_installed_uses_shutil(monkeypatch):
    monkeypatch.setattr("shutil.which", lambda _: None)
    with pytest.raises(SystemExit) as ei:
        C.check_quarto_installed()
    assert "quarto" in str(ei.value)


def test_argparse_mutually_exclusive_sample_and_config(capsys):
    parser = C.build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(["report", "--sample", "A=/p", "--config", "x.yaml"])


def test_argparse_requires_one_of(capsys):
    parser = C.build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(["report"])


def test_argparse_repeated_sample_collected():
    parser = C.build_parser()
    args = parser.parse_args(["report",
                              "--sample", "A=/a",
                              "--sample", "B=/b",
                              "-o", "out.html"])
    assert args.samples == [("A", Path("/a")), ("B", Path("/b"))]


def test_duplicate_sample_name_rejected():
    args_samples = [("A", Path("/a")), ("A", Path("/b"))]
    with pytest.raises(SystemExit):
        C.dedupe_samples(args_samples)
```

- [ ] **Step 2: Run, expect failures**

- [ ] **Step 3: Implement**

Replace `src/cosmx_qc/cli.py` with:

```python
"""cosmx-qc command-line interface."""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import yaml

from cosmx_qc.io import validate_sample_dir


def parse_sample_arg(s: str) -> tuple[str, Path]:
    if "=" not in s:
        raise ValueError(f"--sample must be NAME=PATH, got {s!r}")
    name, _, path = s.partition("=")
    if not name:
        raise ValueError(f"--sample has empty name in {s!r}")
    if not path:
        raise ValueError(f"--sample has empty path in {s!r}")
    return name, Path(path)


def load_config(path: Path) -> tuple[dict[str, Path], str]:
    """Load YAML config. Returns (samples_dict, title)."""
    with open(path) as f:
        data = yaml.safe_load(f) or {}
    if "samples" not in data or not isinstance(data["samples"], dict):
        raise ValueError(f"{path} missing 'samples' mapping")
    samples = {k: Path(v) for k, v in data["samples"].items()}
    title = data.get("title", "CosMx QC Report")
    return samples, title


def dedupe_samples(pairs: list[tuple[str, Path]]) -> dict[str, Path]:
    """Convert [(name, path), ...] to dict, exiting on duplicate names."""
    seen: dict[str, Path] = {}
    for name, p in pairs:
        if name in seen:
            sys.exit(f"error: duplicate sample name {name!r}")
        seen[name] = p
    return seen


def validate_samples(samples: dict[str, Path]) -> None:
    """Exit with a message listing any missing files across all samples."""
    errors: list[str] = []
    for name, d in samples.items():
        try:
            missing = validate_sample_dir(d)
        except (FileNotFoundError, NotADirectoryError) as e:
            errors.append(f"  - {name}: {e}")
            continue
        if missing:
            errors.append(f"  - {name} ({d}): missing {missing}")
    if errors:
        sys.exit("error: invalid sample input(s):\n" + "\n".join(errors))


def check_quarto_installed() -> str:
    """Return path to `quarto` binary, exit with install hint if missing."""
    p = shutil.which("quarto")
    if p is None:
        sys.exit(
            "error: 'quarto' not found in PATH.\n"
            "Install Quarto: https://quarto.org/docs/get-started/\n"
            "  macOS: brew install --cask quarto\n"
            "  Linux: download .deb/.rpm from quarto.org"
        )
    return p


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="cosmx-qc")
    sub = p.add_subparsers(dest="cmd", required=True)

    rep = sub.add_parser("report", help="Generate a CosMx QC HTML report")
    g = rep.add_mutually_exclusive_group(required=True)
    g.add_argument("--sample", action="append", dest="samples",
                   type=parse_sample_arg, metavar="NAME=PATH",
                   help="Sample name and flatfile directory (repeatable)")
    g.add_argument("--config", type=Path,
                   help="YAML config: {samples: {NAME: PATH, ...}, title: ...}")
    rep.add_argument("--output", "-o", type=Path, default=Path("report.html"))
    rep.add_argument("--title", default="CosMx QC Report")
    return p


def cmd_report(args: argparse.Namespace) -> None:  # pragma: no cover until Task 14
    raise NotImplementedError


def main() -> None:  # pragma: no cover until Task 14
    parser = build_parser()
    args = parser.parse_args()
    if args.cmd == "report":
        cmd_report(args)
    else:
        parser.error(f"unknown command {args.cmd!r}")
```

- [ ] **Step 4: Run, expect pass**

```bash
pytest tests/test_cli.py -v
```

Expected: 14 passed.

- [ ] **Step 5: Commit**

```bash
git add src/cosmx_qc/cli.py tests/test_cli.py
git commit -m "feat(qc): add CLI argparse, validation, and config loader"
```

---

## Task 13: render — render_report (tempfile JSON + subprocess to Quarto)

**Files:**
- Create: `src/cosmx_qc/render.py`
- Create: `tests/test_render.py`

- [ ] **Step 1: Failing tests**

Create `tests/test_render.py`:

```python
import json
import subprocess
from pathlib import Path
from unittest.mock import patch
from cosmx_qc import render as R


def test_render_report_writes_config_and_calls_quarto(tmp_path):
    output = tmp_path / "out.html"
    captured = {}

    def fake_run(cmd, env=None, check=False, **kw):
        # capture env var pointing to JSON, read it
        cfg_path = env["COSMX_QC_CONFIG"]
        captured["cmd"] = cmd
        captured["cfg"] = json.loads(Path(cfg_path).read_text())
        # simulate quarto producing the output
        Path(output).write_text("<html></html>")
        return subprocess.CompletedProcess(cmd, 0)

    with patch("cosmx_qc.render.subprocess.run", side_effect=fake_run):
        R.render_report(
            samples={"A": tmp_path / "a", "B": tmp_path / "b"},
            output=output,
            title="T",
        )

    assert captured["cmd"][0] == "quarto"
    assert "render" in captured["cmd"]
    assert captured["cfg"]["samples"] == {"A": str(tmp_path / "a"), "B": str(tmp_path / "b")}
    assert captured["cfg"]["title"] == "T"


def test_render_report_creates_output_parent(tmp_path):
    output = tmp_path / "deeper" / "report.html"
    with patch("cosmx_qc.render.subprocess.run") as run:
        run.return_value = subprocess.CompletedProcess(["quarto"], 0)
        R.render_report(samples={"A": tmp_path}, output=output, title="X")
    assert output.parent.exists()
```

- [ ] **Step 2: Run, expect failures**

- [ ] **Step 3: Implement**

Create `src/cosmx_qc/render.py`:

```python
"""Quarto render orchestration: write JSON config to tempfile, call quarto."""
from __future__ import annotations

import json
import os
import subprocess
import tempfile
from datetime import datetime
from importlib.resources import files
from pathlib import Path


def render_report(samples: dict[str, Path], output: Path,
                  title: str = "CosMx QC Report") -> None:
    """Render report.qmd to `output` HTML via Quarto."""
    output = Path(output).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)

    qmd = files("cosmx_qc").joinpath("report.qmd")
    payload = {
        "samples": {k: str(v) for k, v in samples.items()},
        "title": title,
        "generated_at": datetime.now().isoformat(timespec="seconds"),
    }

    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as f:
        json.dump(payload, f)
        cfg_path = f.name
    try:
        env = {**os.environ, "COSMX_QC_CONFIG": cfg_path}
        subprocess.run(
            ["quarto", "render", str(qmd),
             "--to", "html", "--output", str(output)],
            env=env, check=True,
        )
    finally:
        os.unlink(cfg_path)
```

- [ ] **Step 4: Run, expect pass**

```bash
pytest tests/test_render.py -v
```

Expected: 2 passed.

- [ ] **Step 5: Commit**

```bash
git add src/cosmx_qc/render.py tests/test_render.py
git commit -m "feat(qc): add render_report wrapping quarto render with tempfile config"
```

---

## Task 14: cli — wire cmd_report through to render

**Files:**
- Modify: `src/cosmx_qc/cli.py`
- Modify: `tests/test_cli.py`

- [ ] **Step 1: Add failing tests for cmd_report integration**

Append to `tests/test_cli.py`:

```python
from unittest.mock import patch


def test_cmd_report_inline_samples(tmp_path):
    # build a valid sample dir
    d = tmp_path / "S"
    d.mkdir()
    for pat in ["x_exprMat_file.csv.gz", "x_metadata_file.csv.gz",
                "x_tx_file.csv.gz", "x_fov_positions_file.csv.gz"]:
        (d / pat).write_bytes(b"")
    parser = C.build_parser()
    args = parser.parse_args(["report", "--sample", f"S={d}",
                              "-o", str(tmp_path / "r.html")])
    with patch("cosmx_qc.cli.check_quarto_installed", return_value="/usr/bin/quarto"), \
         patch("cosmx_qc.cli.render_report") as rr:
        C.cmd_report(args)
    rr.assert_called_once()
    kwargs = rr.call_args.kwargs or rr.call_args[1]
    assert kwargs["samples"] == {"S": d}
    assert kwargs["output"] == tmp_path / "r.html"
    assert kwargs["title"] == "CosMx QC Report"


def test_cmd_report_config_file(tmp_path):
    d = tmp_path / "S"
    d.mkdir()
    for pat in ["x_exprMat_file.csv.gz", "x_metadata_file.csv.gz",
                "x_tx_file.csv.gz", "x_fov_positions_file.csv.gz"]:
        (d / pat).write_bytes(b"")
    cfg = tmp_path / "c.yaml"
    cfg.write_text(yaml.safe_dump({"samples": {"S": str(d)}, "title": "Cfg"}))
    parser = C.build_parser()
    args = parser.parse_args(["report", "--config", str(cfg),
                              "-o", str(tmp_path / "r.html")])
    with patch("cosmx_qc.cli.check_quarto_installed", return_value="q"), \
         patch("cosmx_qc.cli.render_report") as rr:
        C.cmd_report(args)
    kwargs = rr.call_args.kwargs or rr.call_args[1]
    assert kwargs["samples"] == {"S": d}
    assert kwargs["title"] == "Cfg"
```

- [ ] **Step 2: Run, expect failures**

- [ ] **Step 3: Implement cmd_report**

Replace the `cmd_report` and `main` placeholders in `src/cosmx_qc/cli.py`:

```python
from cosmx_qc.render import render_report  # add at top with other imports


def cmd_report(args: argparse.Namespace) -> None:
    if args.samples is not None:
        samples = dedupe_samples(args.samples)
        title = args.title
    else:
        samples, title = load_config(args.config)
        # CLI --title overrides only if user explicitly set it (not default)
        # argparse default is the same string, so we accept config's value as primary
    if str(args.output).endswith(".html") is False:
        sys.exit(f"error: --output must end in '.html' (got {args.output})")
    validate_samples(samples)
    check_quarto_installed()
    render_report(samples=samples, output=args.output, title=title)


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    if args.cmd == "report":
        cmd_report(args)
    else:
        parser.error(f"unknown command {args.cmd!r}")
```

- [ ] **Step 4: Run all tests**

```bash
pytest tests/ -v
```

Expected: all green (16 cli tests + others).

- [ ] **Step 5: Commit**

```bash
git add src/cosmx_qc/cli.py tests/test_cli.py
git commit -m "feat(qc): wire cmd_report through validation to render_report"
```

---

## Task 15: report.qmd — header, setup chunk, and tabset helpers

**Files:**
- Create: `src/cosmx_qc/report.qmd`

This file is hand-tested via Quarto, so no pytest. It must read
`os.environ["COSMX_QC_CONFIG"]`, load each sample once, and expose helpers
that emit a `panel-tabset` with a "Combined" tab plus one tab per sample.

The tabset technique uses Quarto's `#| output: asis` to print raw markdown
(including base64-embedded PNG images), so Quarto re-parses the printed
content as part of the document. This sidesteps the limitation that
`panel-tabset` needs markdown headings — we generate them dynamically.

- [ ] **Step 1: Create report.qmd with YAML, setup chunk, and helpers (no sections yet)**

Create `src/cosmx_qc/report.qmd`:

````qmd
---
title: "CosMx QC Report"
format:
  html:
    embed-resources: true
    toc: true
    toc-depth: 2
    toc-location: left
    theme: cosmo
    code-fold: true
    code-tools: false
execute:
  echo: false
  warning: false
  message: false
---

```{python}
#| label: setup
import base64
import io as _io
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from cosmx_qc import io, metrics as M, plots as P

cfg_path = os.environ["COSMX_QC_CONFIG"]
with open(cfg_path) as f:
    config = json.load(f)

samples = {name: io.load_sample(name, path)
           for name, path in config["samples"].items()}
sample_names = list(samples.keys())


def _fig_to_md_img(fig) -> str:
    """Render a matplotlib Figure to an inline base64 markdown image."""
    buf = _io.BytesIO()
    fig.savefig(buf, format="png", dpi=120, bbox_inches="tight")
    plt.close(fig)
    b64 = base64.b64encode(buf.getvalue()).decode()
    return f"![](data:image/png;base64,{b64})"


def _single_axis(plot_for_sample_fn, *args, label: str, **kwargs):
    """Wrap a *_for_sample plotter with its own 1-axis figure."""
    fig, axes = P.make_combined_axes(1)
    plot_for_sample_fn(*args, ax=axes[0], label=label, **kwargs)
    return fig


def emit_tabset(combined_factory, single_factory) -> None:
    """Print markdown for a panel-tabset block.

    `combined_factory()` returns a Figure with all samples side-by-side.
    `single_factory(name)` returns a Figure for one sample.
    Must be called from a `#| output: asis` chunk.
    """
    print("\n::: {.panel-tabset}\n")
    print("### Combined\n")
    print(_fig_to_md_img(combined_factory()))
    print()
    for name in sample_names:
        print(f"### {name}\n")
        print(_fig_to_md_img(single_factory(name)))
        print()
    print(":::\n")
```
````

- [ ] **Step 2: Smoke-render the skeleton (if Quarto is installed)**

If Quarto is not on this host, skip — validated in Task 18.

```bash
cat > /tmp/cfg.json <<'EOF'
{"samples": {"mini": "tests/fixtures/mini_sample"}, "title": "Smoke", "generated_at": "2026-05-05T00:00:00"}
EOF
COSMX_QC_CONFIG=/tmp/cfg.json quarto render src/cosmx_qc/report.qmd --to html -o /tmp/smoke.html
ls /tmp/smoke.html
```

Expected: `/tmp/smoke.html` exists. The body will be near-empty (no
sections yet), but the YAML header and setup chunk must execute without
error.

- [ ] **Step 3: Commit**

```bash
git add src/cosmx_qc/report.qmd
git commit -m "feat(qc): add report.qmd setup chunk and tabset emission helpers"
```

---

## Task 16: report.qmd — fill in all 13 sections with combined + per-sample tabsets

**Files:**
- Modify: `src/cosmx_qc/report.qmd`

Each metric section uses the `emit_tabset(combined_factory, single_factory)`
helper from Task 15. The two factories close over the section's metric +
plot choice. The pattern is:

```python
#| output: asis
data = {n: <metric_fn>(sd, **opts) for n, sd in samples.items()}
emit_tabset(
    combined_factory=lambda: P.<plot_combined>(data, y_col=..., ylabel=...),
    single_factory =lambda name: _single_axis(
        P.<plot_for_sample>, data[name], y_col=..., ylabel=..., label=name),
)
```

`<plot_combined>` is `plot_count_per_fov_combined` for bar charts and
`plot_box_per_fov_combined` for boxplots; the matching `_for_sample`
follows the same suffix.

- [ ] **Step 1: Append all 13 sections to report.qmd**

Append (after the setup chunk):

````qmd
# CosMx introduction

CosMx Spatial Molecular Imager (SMI) is a high-plex spatial transcriptomics
platform that detects individual RNA molecules in fixed tissue with
sub-cellular resolution. This report summarises per-FOV and per-cell
quality metrics across the configured samples. Each section provides a
*Combined* cross-sample view followed by per-sample tabs.

# Overview (per FOV)

## Global view of entire run

```{python}
rows = []
for name, sd in samples.items():
    t = M.global_view_table(sd)
    t.insert(0, "sample", name)
    rows.append(t)
pd.concat(rows, ignore_index=True)
```

## Profiling Position (FOV index)

```{python}
#| output: asis
data = {n: sd.fov_pos for n, sd in samples.items()}
emit_tabset(
    combined_factory=lambda: P.plot_fov_positions_combined(data),
    single_factory =lambda name: _single_axis(
        P.plot_fov_positions_for_sample, data[name], label=name),
)
```

## Unassigned transcripts

```{python}
#| output: asis
data = {n: M.unassigned_transcripts_per_fov(sd) for n, sd in samples.items()}
emit_tabset(
    combined_factory=lambda: P.plot_count_per_fov_combined(
        data, y_col="unassigned", ylabel="Unassigned transcripts"),
    single_factory =lambda name: _single_axis(
        P.plot_count_per_fov_for_sample, data[name],
        y_col="unassigned", ylabel="Unassigned transcripts", label=name),
)
```

## Unique transcripts

```{python}
#| output: asis
data = {n: M.unique_transcripts_per_fov(sd) for n, sd in samples.items()}
emit_tabset(
    combined_factory=lambda: P.plot_count_per_fov_combined(
        data, y_col="unique", ylabel="Unique transcript species"),
    single_factory =lambda name: _single_axis(
        P.plot_count_per_fov_for_sample, data[name],
        y_col="unique", ylabel="Unique transcript species", label=name),
)
```

# Transcripts of All (per FOV)

```{python}
#| output: asis
data = {n: M.transcripts_per_fov(sd, kind="all") for n, sd in samples.items()}
emit_tabset(
    combined_factory=lambda: P.plot_count_per_fov_combined(
        data, y_col="count", ylabel="Transcripts (all)"),
    single_factory =lambda name: _single_axis(
        P.plot_count_per_fov_for_sample, data[name],
        y_col="count", ylabel="Transcripts (all)", label=name),
)
```

# Transcripts of Real Genes (per FOV)

```{python}
#| output: asis
data = {n: M.transcripts_per_fov(sd, kind="gene") for n, sd in samples.items()}
emit_tabset(
    combined_factory=lambda: P.plot_count_per_fov_combined(
        data, y_col="count", ylabel="Transcripts (real genes)"),
    single_factory =lambda name: _single_axis(
        P.plot_count_per_fov_for_sample, data[name],
        y_col="count", ylabel="Transcripts (real genes)", label=name),
)
```

# Transcripts of Negative Probes (per FOV)

```{python}
#| output: asis
data = {n: M.transcripts_per_fov(sd, kind="negative") for n, sd in samples.items()}
emit_tabset(
    combined_factory=lambda: P.plot_count_per_fov_combined(
        data, y_col="count", ylabel="Transcripts (negative probes)"),
    single_factory =lambda name: _single_axis(
        P.plot_count_per_fov_for_sample, data[name],
        y_col="count", ylabel="Transcripts (negative probes)", label=name),
)
```

# Transcripts of False Codes (per FOV)

```{python}
#| output: asis
data = {n: M.transcripts_per_fov(sd, kind="falsecode") for n, sd in samples.items()}
emit_tabset(
    combined_factory=lambda: P.plot_count_per_fov_combined(
        data, y_col="count", ylabel="Transcripts (false codes)"),
    single_factory =lambda name: _single_axis(
        P.plot_count_per_fov_for_sample, data[name],
        y_col="count", ylabel="Transcripts (false codes)", label=name),
)
```

# Overall Error Rates (per FOV)

## Negative rate

```{python}
#| output: asis
data = {n: M.error_rates_per_fov(sd) for n, sd in samples.items()}
emit_tabset(
    combined_factory=lambda: P.plot_count_per_fov_combined(
        data, y_col="negative_rate", ylabel="Negative rate"),
    single_factory =lambda name: _single_axis(
        P.plot_count_per_fov_for_sample, data[name],
        y_col="negative_rate", ylabel="Negative rate", label=name),
)
```

## False-code rate

```{python}
#| output: asis
emit_tabset(
    combined_factory=lambda: P.plot_count_per_fov_combined(
        data, y_col="false_code_rate", ylabel="False-code rate"),
    single_factory =lambda name: _single_axis(
        P.plot_count_per_fov_for_sample, data[name],
        y_col="false_code_rate", ylabel="False-code rate", label=name),
)
```

# Cell number and size

## Total number of cells (per FOV)

```{python}
#| output: asis
data = {n: M.cell_count_per_fov(sd) for n, sd in samples.items()}
emit_tabset(
    combined_factory=lambda: P.plot_count_per_fov_combined(
        data, y_col="n_cells", ylabel="Cells per FOV"),
    single_factory =lambda name: _single_axis(
        P.plot_count_per_fov_for_sample, data[name],
        y_col="n_cells", ylabel="Cells per FOV", label=name),
)
```

## Size of cells (per cell per FOV)

```{python}
#| output: asis
data = {n: M.cell_area_per_fov(sd) for n, sd in samples.items()}
emit_tabset(
    combined_factory=lambda: P.plot_box_per_fov_combined(
        data, y_col="Area.um2", ylabel="Cell area (um^2)"),
    single_factory =lambda name: _single_axis(
        P.plot_box_per_fov_for_sample, data[name],
        y_col="Area.um2", ylabel="Cell area (um^2)", label=name),
)
```

# Transcripts of All (per cell per FOV)

```{python}
#| output: asis
data = {n: M.transcripts_per_cell_per_fov(sd, kind="all") for n, sd in samples.items()}
emit_tabset(
    combined_factory=lambda: P.plot_box_per_fov_combined(
        data, y_col="count", ylabel="Transcripts/cell (all)"),
    single_factory =lambda name: _single_axis(
        P.plot_box_per_fov_for_sample, data[name],
        y_col="count", ylabel="Transcripts/cell (all)", label=name),
)
```

# Transcripts of Real Genes (per cell per FOV)

```{python}
#| output: asis
data = {n: M.transcripts_per_cell_per_fov(sd, kind="gene") for n, sd in samples.items()}
emit_tabset(
    combined_factory=lambda: P.plot_box_per_fov_combined(
        data, y_col="count", ylabel="Transcripts/cell (real genes)"),
    single_factory =lambda name: _single_axis(
        P.plot_box_per_fov_for_sample, data[name],
        y_col="count", ylabel="Transcripts/cell (real genes)", label=name),
)
```

# Transcripts of Negative Probes (per cell per FOV)

```{python}
#| output: asis
data = {n: M.transcripts_per_cell_per_fov(sd, kind="negative") for n, sd in samples.items()}
emit_tabset(
    combined_factory=lambda: P.plot_box_per_fov_combined(
        data, y_col="count", ylabel="Transcripts/cell (negative)"),
    single_factory =lambda name: _single_axis(
        P.plot_box_per_fov_for_sample, data[name],
        y_col="count", ylabel="Transcripts/cell (negative)", label=name),
)
```

# Transcripts of False Codes (per cell per FOV)

```{python}
#| output: asis
data = {n: M.transcripts_per_cell_per_fov(sd, kind="falsecode") for n, sd in samples.items()}
emit_tabset(
    combined_factory=lambda: P.plot_box_per_fov_combined(
        data, y_col="count", ylabel="Transcripts/cell (false codes)"),
    single_factory =lambda name: _single_axis(
        P.plot_box_per_fov_for_sample, data[name],
        y_col="count", ylabel="Transcripts/cell (false codes)", label=name),
)
```

# Overall Error Rates (per cell per FOV)

## Negative rate

```{python}
#| output: asis
data = {n: M.error_rates_per_cell_per_fov(sd) for n, sd in samples.items()}
emit_tabset(
    combined_factory=lambda: P.plot_box_per_fov_combined(
        data, y_col="negative_rate", ylabel="Negative rate (per cell)"),
    single_factory =lambda name: _single_axis(
        P.plot_box_per_fov_for_sample, data[name],
        y_col="negative_rate", ylabel="Negative rate (per cell)", label=name),
)
```

## False-code rate

```{python}
#| output: asis
emit_tabset(
    combined_factory=lambda: P.plot_box_per_fov_combined(
        data, y_col="false_code_rate", ylabel="False-code rate (per cell)"),
    single_factory =lambda name: _single_axis(
        P.plot_box_per_fov_for_sample, data[name],
        y_col="false_code_rate", ylabel="False-code rate (per cell)", label=name),
)
```
````

- [ ] **Step 2: Add `cosmx_qc.report.qmd` to package data**

Modify `pyproject.toml` to include the `.qmd` file in the installed package:

```toml
[tool.setuptools.package-data]
"cosmx_qc" = ["report.qmd"]
```

(Add this section after the existing `[tool.setuptools.package-dir]` block.)

- [ ] **Step 3: Reinstall and verify report.qmd is packaged**

```bash
uv pip install -e '.[dev]'
python -c "from importlib.resources import files; print(files('cosmx_qc').joinpath('report.qmd'))"
```

Expected: prints a path that exists.

- [ ] **Step 4: Commit**

```bash
git add src/cosmx_qc/report.qmd pyproject.toml
git commit -m "feat(qc): fill in 13-section report.qmd and package the template"
```

---

## Task 17: rewrite __init__.py + delete old qc files + update README

**Files:**
- Modify: `src/cosmx_qc/__init__.py`
- Delete: see "File structure" → "Delete" list above
- Modify: `README.md`

- [ ] **Step 1: Rewrite __init__.py with minimal exports**

Replace `src/cosmx_qc/__init__.py` with:

```python
"""CosMx quality-control report tooling.

Public CLI entry: `cosmx-qc report --sample NAME=PATH ... -o report.html`
Programmatic entry: `cosmx_qc.render.render_report`.
"""
from cosmx_qc.io import SampleData, load_sample, validate_sample_dir
from cosmx_qc.render import render_report

__all__ = ["SampleData", "load_sample", "validate_sample_dir", "render_report"]
```

- [ ] **Step 2: Delete obsolete files**

```bash
git rm \
  src/cosmx_qc/cosmx_unified_viz.py \
  src/cosmx_qc/htmldesign.py \
  src/cosmx_qc/getfilepath.py \
  src/cosmx_qc/fov_metrics.py \
  src/cosmx_qc/cell_metrics.py \
  src/cosmx_qc/file_paths.csv \
  README_unified.md \
  examples/cosmx_unified_examples.py 2>/dev/null || true

git rm -rf \
  src/cosmx_qc/output_plots_unified \
  src/cosmx_qc/unified_output_plots \
  src/cosmx_qc/unified_output_plots_2 2>/dev/null || true
```

If any path doesn't exist (e.g., the examples file), the `|| true` keeps the
script going; verify with `git status`.

- [ ] **Step 3: Update README.md**

Edit `README.md`. Replace the line that says `- TODO: QC report generation for downloaded datasets.` with a real usage block:

```markdown
- QC report generation for downloaded datasets via `cosmx-qc`.
```

After the `### atomx module` section, add a new section:

```markdown
### `qc` module

Generate an HTML QC report from one or more CosMx flatfile sample directories.
Requires Quarto CLI (≥ 1.4) installed system-wide:
<https://quarto.org/docs/get-started/>.

```bash
cosmx-qc report \
    --sample Pembro7=/path/to/Pembro7/RNA_flatfiles \
    --sample Pembro8=/path/to/Pembro8/RNA_flatfiles \
    --output qc_report.html
```

Or with a config file:

```yaml
# samples.yaml
samples:
  Pembro7: /path/to/Pembro7/RNA_flatfiles
  Pembro8: /path/to/Pembro8/RNA_flatfiles
title: "CosMx QC: Pembro cohort"
```

```bash
cosmx-qc report --config samples.yaml --output qc_report.html
```
```

- [ ] **Step 4: Run full test suite to confirm nothing is broken**

```bash
pytest tests/ -v
```

Expected: all green.

- [ ] **Step 5: Commit**

```bash
git add src/cosmx_qc/__init__.py README.md
git commit -m "refactor(qc): drop legacy qc modules, rewrite __init__, update README"
```

---

## Task 18: Manual end-to-end on Pembro7 + Pembro8

**Files:** none (manual verification)

This task verifies the entire pipeline against real data. Quarto must be
installed on the host (`quarto --version` succeeds). On a host without
Quarto, install per:
<https://quarto.org/docs/get-started/>.

- [ ] **Step 1: Verify Quarto present**

```bash
quarto --version
```

Expected: a version ≥ 1.4. If Quarto is missing, install before continuing.

- [ ] **Step 2: Reinstall package and run end-to-end**

```bash
uv pip install -e '.[dev]'

cosmx-qc report \
    --sample Pembro7=/mnt/nfs/storage/rthuang/data/EC_Pembro/CosMx/Pembro7/RNA_flatfiles \
    --sample Pembro8=/mnt/nfs/storage/rthuang/data/EC_Pembro/CosMx/Pembro8/RNA_flatfiles \
    --output /tmp/cosmx_qc_pembro.html
```

Expected: command exits 0, prints Quarto's progress, and produces
`/tmp/cosmx_qc_pembro.html`. The file should be non-trivial in size
(likely 20–200 MB given two large samples and embedded resources).

- [ ] **Step 3: Spot-check the output**

```bash
ls -lh /tmp/cosmx_qc_pembro.html
grep -c "panel-tabset\|Pembro7\|Pembro8" /tmp/cosmx_qc_pembro.html
```

Open in a browser; confirm:
- TOC shows all 13 H1 sections.
- Each metric section has a "Combined" tab plus one tab per sample (Pembro7, Pembro8).
- The "Combined" tab shows side-by-side panels for the two samples.
- Each per-sample tab shows a single-axis figure for that sample.
- Cell-area boxplots render without error.
- No tracebacks visible in any code chunk.

- [ ] **Step 4: Run with --config form to confirm the alternative path**

```bash
cat <<'EOF' > /tmp/pembro_samples.yaml
samples:
  Pembro7: /mnt/nfs/storage/rthuang/data/EC_Pembro/CosMx/Pembro7/RNA_flatfiles
  Pembro8: /mnt/nfs/storage/rthuang/data/EC_Pembro/CosMx/Pembro8/RNA_flatfiles
title: "CosMx QC: Pembro cohort"
EOF
cosmx-qc report --config /tmp/pembro_samples.yaml --output /tmp/cosmx_qc_pembro_cfg.html
diff <(grep -c panel /tmp/cosmx_qc_pembro.html) <(grep -c panel /tmp/cosmx_qc_pembro_cfg.html)
```

Expected: both files render successfully and have similar content.

- [ ] **Step 5: Document the result**

If the output is correct, no commit needed. If you discovered any
issue (missing column, wrong axis label, etc.), open an issue or a
follow-up task — do not bundle hot-fixes into this task; pick them up
in a new branch.

---

## Self-review notes (for plan author, not engineer)

Spec coverage check:
- ✅ CLI both `--sample` and `--config` (Tasks 11, 12)
- ✅ Pre-flight validation (Task 12)
- ✅ Quarto-not-found behavior (Task 12)
- ✅ JSON tempfile + env var (Task 13)
- ✅ Single self-contained HTML (Task 15 YAML header)
- ✅ All 13 sections (Task 16)
- ✅ Per-FOV transcripts × 4 kinds (Task 5)
- ✅ Per-cell transcripts × 4 kinds (Task 8)
- ✅ Error rates per FOV / per cell (Tasks 7, 8)
- ✅ Global view table (Task 7)
- ✅ Cell count / area (Task 6)
- ✅ Migration / cleanup (Task 17)
- ✅ E2E (Task 18)

Deferred to v2 (per spec): polygon viz, multi-slide, `--no-embed`,
`--keep-plots`. Per-sample tabsets are implemented in v1 (Tasks 15–16).
