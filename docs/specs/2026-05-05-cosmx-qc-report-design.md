# CosMx QC Report — Design

- Date: 2026-05-05
- Status: Approved (pending user spec review)
- Related code: `src/cosmx_qc/`

## Goal

Build a Python CLI tool (`cosmx-qc`) that consumes one or more CosMx flatfile
sample directories and produces a single self-contained HTML quality-control
report. The report mirrors the 13-section structure of the R tutorial at
[SizunJiangLab/Tutorials CosMx/02_QC](https://github.com/SizunJiangLab/Tutorials/tree/main/CosMx/02_QC),
but is rendered via Quarto driving Python (no R toolchain).

## Non-goals (v1)

- Cell-segmentation polygon visualisation (`*-polygons.csv.gz` is not consumed
  in v1; deferred).
- The AtoMx `fovs` log file (used in R but not present in flatfile-only
  exports such as Pembro7/8 — omitted).
- Multi-slide-per-directory inputs. One sample directory ↔ one slide.
- Persisting individual PNG/CSV artefacts alongside the HTML.
  (`--keep-plots` is a possible v2 extension; not in v1.)

## Inputs

Each sample directory must contain the four CosMx flatfiles below. The
filename prefix is project-specific and matched by glob:

| Glob pattern                       | Required | Used for |
|------------------------------------|----------|----------|
| `*_exprMat_file.csv.gz`            | yes      | per-cell expression matrix |
| `*_metadata_file.csv.gz`           | yes      | cell metadata (`Area.um2`, `nFeature_RNA`, `slide_ID`, …) |
| `*_tx_file.csv.gz`                 | yes      | individual transcript records (`fov`, `cell_ID`, `target`, `CellComp`) |
| `*_fov_positions_file.csv.gz`      | yes      | FOV physical positions (`X_mm`, `Y_mm`) |
| `*-polygons.csv.gz`                | no       | not used in v1 |

## CLI

Two mutually exclusive ways to specify samples:

```bash
# Inline flag form (repeatable)
cosmx-qc report \
    --sample Pembro7=/mnt/.../Pembro7/RNA_flatfiles \
    --sample Pembro8=/mnt/.../Pembro8/RNA_flatfiles \
    --output report.html

# Config-file form
cosmx-qc report --config samples.yaml --output report.html
```

Config-file format (`yaml`):

```yaml
samples:
  Pembro7: /mnt/.../Pembro7/RNA_flatfiles
  Pembro8: /mnt/.../Pembro8/RNA_flatfiles
title: "CosMx QC: Pembro7 vs Pembro8"
```

CLI options:

| Flag | Type | Default | Notes |
|------|------|---------|-------|
| `--sample NAME=PATH` | str | — | Repeatable. Mutually exclusive with `--config`. |
| `--config PATH` | path | — | YAML file. Mutually exclusive with `--sample`. |
| `--output`, `-o` | path | `report.html` | Output HTML path (must end in `.html`). |
| `--title` | str | `"CosMx QC Report"` | Report title. |

Pre-flight validation (before invoking Quarto):

1. Each sample directory exists and is a directory.
2. Each sample directory contains all four required flatfiles (glob match).
3. `quarto` is on `PATH`. If not, exit with installation instructions.
4. Sample names are unique. With `action="append"`, argparse collects all
   occurrences in order; the CLI rejects the invocation if any name appears
   more than once.

## Architecture and data flow

```
$ cosmx-qc report --sample Pembro7=/...  --sample Pembro8=/...  -o report.html
        │
        ▼
┌─────────────── cli.py ───────────────┐
│ 1. Parse args (--sample / --config)  │
│ 2. Validate sample dirs              │
│ 3. shutil.which("quarto")            │
│ 4. Write JSON to NamedTemporaryFile  │
│ 5. env["COSMX_QC_CONFIG"] = path     │
│ 6. subprocess.run(["quarto","render",│
│        report.qmd, "--output", ...]) │
│ 7. Tempfile auto-cleaned             │
└──────────────────────────────────────┘
        │
        ▼
┌─────────────── report.qmd ───────────┐
│ Setup chunk: read env var → config   │
│ for name, dir in config["samples"]:  │
│     data[name] = io.load_sample(...) │
│ # 13 sections call metrics + plots   │
└──────────────────────────────────────┘
        │
        ▼
   report.html  (self-contained, embed-resources)
```

Each `SampleData` is loaded exactly once per render and reused by every
section.

JSON config written to tempfile:

```json
{
  "samples": {"Pembro7": "/path1", "Pembro8": "/path2"},
  "title": "CosMx QC Report",
  "generated_at": "2026-05-05T12:34:56"
}
```

## Module layout

```
src/cosmx_qc/
├── __init__.py        # minimal public exports
├── cli.py             # argparse, validation, subprocess to quarto
├── io.py              # load_sample(name, dir) → SampleData
├── metrics.py         # pure DataFrame in / DataFrame out functions
├── plots.py           # matplotlib plotting (metric → Figure / ax)
├── render.py          # tempfile JSON + quarto render wrapper
└── report.qmd         # Quarto template (13 sections)
```

### `io.py`

```python
@dataclass
class SampleData:
    name: str
    exprmat: pd.DataFrame   # cell × gene; cols: fov, cell_ID, <gene...>
    metadata: pd.DataFrame  # fov, cell_ID, Area.um2, nFeature_RNA, slide_ID
    tx: pd.DataFrame        # fov, cell_ID, target, CellComp, x_*, y_*
    fov_pos: pd.DataFrame   # FOV, X_mm, Y_mm  (renamed from x_global_mm)

def load_sample(name: str, dir: Path) -> SampleData: ...
def validate_sample_dir(dir: Path) -> list[str]:
    """Return list of missing required flatfile globs (empty list = OK)."""
```

Column names are normalised on load:

- `*fov*` columns → `fov`
- `*cell.*ID*` columns → `cell_ID`
- `x_global_mm` / `y_global_mm` → `X_mm` / `Y_mm`

### `metrics.py`

One function per chart (table) in the report. All functions are pure:
`SampleData` (or `dict[str, SampleData]`) in, `DataFrame` out, no side
effects, no plotting.

```python
# Section 2 — Overview (per FOV)
def global_view_table(sd: SampleData) -> pd.DataFrame
def unassigned_transcripts_per_fov(sd: SampleData) -> pd.DataFrame
def unique_transcripts_per_fov(sd: SampleData) -> pd.DataFrame

# Sections 3–6 — Transcripts per FOV
TargetKind = Literal["all", "gene", "negative", "falsecode"]
def transcripts_per_fov(sd: SampleData, kind: TargetKind) -> pd.DataFrame

# Section 7 — Error rates per FOV
def error_rates_per_fov(sd: SampleData) -> pd.DataFrame
    # columns: fov, negative_rate, false_code_rate

# Section 8 — Cell number and size
def cell_count_per_fov(sd: SampleData) -> pd.DataFrame
def cell_area_per_fov(sd: SampleData) -> pd.DataFrame

# Sections 9–12 — Transcripts per cell per FOV
def transcripts_per_cell_per_fov(sd: SampleData, kind: TargetKind) -> pd.DataFrame

# Section 13 — Error rates per cell per FOV
def error_rates_per_cell_per_fov(sd: SampleData) -> pd.DataFrame
```

Target classification (mirroring R):

- Column names / `target` values starting with `Negative` (case-sensitive
  prefix match, as in CosMx flatfiles) → negative probes
- Column names / `target` values starting with `SystemControl` → false codes
- All other expression columns → real genes
- `all` = gene + negative + falsecode

Error-rate definitions (used in sections 7 and 13):

- `negative_rate = sum(negative_counts) / sum(gene_counts)`
  per FOV (section 7) or per cell per FOV aggregated (section 13)
- `false_code_rate = sum(systemcontrol_counts) / sum(gene_counts)`
  computed analogously

### `plots.py`

Each metric has two plotting helpers:

```python
def plot_<metric>_for_sample(metric_df: pd.DataFrame, ax: Axes, **style) -> None
def plot_<metric>_combined(metric_by_sample: dict[str, pd.DataFrame]) -> Figure
```

`_for_sample` accepts an external `ax` so the caller controls the figure
(used inside `_combined` and inside per-sample tabs in `.qmd`).

`_combined` arranges N samples horizontally:
`fig, axes = plt.subplots(1, N, constrained_layout=True, figsize=...)` then
calls `_for_sample` per axis.

Visual style is captured in a single `_STYLE` dict at the top of `plots.py`
(palette, font sizes, dpi, boxplot defaults). Initial values are extracted
from the existing `fov_metrics.py` / `cell_metrics.py` so the matplotlib
look is preserved.

### `cli.py`

```python
def main():
    p = argparse.ArgumentParser(prog="cosmx-qc")
    sub = p.add_subparsers(dest="cmd", required=True)

    rep = sub.add_parser("report", help="Generate a CosMx QC HTML report")
    g = rep.add_mutually_exclusive_group(required=True)
    g.add_argument("--sample", action="append",
                   type=parse_sample_arg, dest="samples",
                   metavar="NAME=PATH",
                   help="Sample name and flatfile directory (repeatable)")
    g.add_argument("--config", type=Path,
                   help="YAML config file with samples and metadata")
    rep.add_argument("--output", "-o", type=Path,
                     default=Path("report.html"))
    rep.add_argument("--title", default="CosMx QC Report")
    rep.set_defaults(func=cmd_report)

    args = p.parse_args()
    args.func(args)

def parse_sample_arg(s: str) -> tuple[str, Path]:
    """'Pembro7=/path' → ('Pembro7', Path('/path')). Raises on malformed."""
```

`cosmx-qc` is registered in `pyproject.toml` under `[project.scripts]`.

### `render.py`

```python
def render_report(samples: dict[str, Path],
                  output: Path,
                  title: str = "CosMx QC Report") -> None:
    qmd = files("cosmx_qc").joinpath("report.qmd")
    payload = {
        "samples": {k: str(v) for k, v in samples.items()},
        "title": title,
        "generated_at": datetime.now().isoformat(),
    }
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as f:
        json.dump(payload, f); cfg_path = f.name
    try:
        env = {**os.environ, "COSMX_QC_CONFIG": cfg_path}
        subprocess.run(
            ["quarto", "render", str(qmd),
             "--to", "html", "--output", str(output.resolve())],
            env=env, check=True)
    finally:
        os.unlink(cfg_path)
```

`importlib.resources.files()` is used so the `.qmd` is found whether the
package is installed or run from source.

### `report.qmd`

YAML header:

```yaml
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
---
```

Each of the 13 sections follows the same shape:

```qmd
# Cell number and size

## Total number of cells (per FOV)

::: {.panel-tabset}

### Combined
```{python}
fig = plots.plot_cell_count_combined(
    {n: metrics.cell_count_per_fov(d) for n, d in data.items()}
)
fig
```

### Pembro7
```{python}
fig, ax = plt.subplots(constrained_layout=True)
plots.plot_cell_count_for_sample(metrics.cell_count_per_fov(data["Pembro7"]), ax)
fig
```

### Pembro8
```{python}
...
```

:::
```

Per-sample tabs are emitted programmatically in the setup chunk via Quarto's
inline expression mechanism so the `.qmd` does not hard-code sample names.
The exact emission technique (inline Markdown via `Markdown()` from
IPython.display, or templated strings) is to be chosen during
implementation; the constraint is that adding a new sample requires no edit
to `report.qmd`.

## Section list (mirroring R)

1. CosMx introduction (text only)
2. Overview (per FOV) — global table, FOV index map, unassigned, unique
3. Transcripts of All (per FOV)
4. Transcripts of Real Genes (per FOV)
5. Transcripts of Negative Probes (per FOV)
6. Transcripts of False Codes (per FOV)
7. Overall Error Rates (per FOV)
8. Cell number and size
9. Transcripts of All (per cell per FOV)
10. Transcripts of Real Genes (per cell per FOV)
11. Transcripts of Negative Probes (per cell per FOV)
12. Transcripts of False Codes (per cell per FOV)
13. Overall Error Rates (per cell per FOV)

## Migration / cleanup

**Delete** (recoverable from git history):

- `src/cosmx_qc/cosmx_unified_viz.py` (~2213 lines)
- `src/cosmx_qc/htmldesign.py` (~1232 lines)
- `src/cosmx_qc/getfilepath.py` (~395 lines)
- `src/cosmx_qc/fov_metrics.py` (~1575 lines)
- `src/cosmx_qc/cell_metrics.py` (~450 lines)
- `src/cosmx_qc/file_paths.csv`
- `src/cosmx_qc/output_plots_unified/`
- `src/cosmx_qc/unified_output_plots/`
- `src/cosmx_qc/unified_output_plots_2/`
- `README_unified.md`
- Any examples under `examples/` that import the deleted modules
  (`cosmx_unified_examples.py`, `cosmx_viz_examples.py` if present)

**Extract before deleting** (visual-style reference for new `plots.py`):

- Colour palette / viridis usage
- Boxplot, bar-chart, scatter default kwargs
- Font sizes, dpi, figure size conventions
- FOV-position panel arrangement

These are copied as constants into `_STYLE` in `plots.py`; no logic from
the old files is preserved.

**Add** to `pyproject.toml`:

```toml
[project]
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

[project.scripts]
atomx-transfer = "cosmx_utils.atomx.transfer:main"
cosmx-qc       = "cosmx_qc.cli:main"
cosmx-utils    = "cosmx_utils.cli:main"
```

**System dependency**: Quarto CLI ≥ 1.4. The README is updated with an
install pointer.

## Testing

`tests/` (new):

1. `test_io.py` — round-trip a tiny synthetic sample (handcrafted
   `tests/fixtures/mini_sample/*.csv.gz`, < 1 KB each); assert column names
   and dtypes match the `SampleData` contract.
2. `test_metrics.py` — primary test surface. Synthetic `SampleData`
   built from in-memory DataFrames; assert metric values exactly.
3. `test_plots.py` — smoke tests only: each `plot_*` returns a `Figure`
   without raising. No pixel comparison.
4. `test_cli.py` — `parse_sample_arg`, argument parsing, validation
   error messages.

End-to-end (running `quarto render` against a fixture) is **deferred**:
Quarto is a system dependency that is awkward to provision in CI. v1 ships
with manual e2e on Pembro7/8.

## Boundary cases

| Situation | Behaviour |
|-----------|-----------|
| Sample dir missing a required flatfile | CLI exits with a list of missing files before invoking Quarto |
| `--sample foo=/no/such/dir` | CLI exits with "directory not found" |
| `quarto` not on PATH | CLI exits with install instructions |
| Duplicate `--sample NAME=...` (same NAME) | CLI exits with "duplicate sample name" |
| `--output` does not end in `.html` | CLI exits with a clear error |
| `--output` parent dir does not exist | CLI creates it (`mkdir -p`) |
| Sample with zero cells / empty file | metrics return empty DataFrames; plots render an empty axis with a "no data" note |

## Open questions deferred to v2

- Whether to add `--keep-plots` to also dump individual PNGs.
- Whether to add `--no-embed` if HTML files become unwieldy (>100 MB).
- Polygon-based cell visualisations.
- Multi-slide samples.
