# cosmx-qc

[![CI](https://github.com/wuwenrui555/cosmx-qc/actions/workflows/ci.yml/badge.svg)](https://github.com/wuwenrui555/cosmx-qc/actions/workflows/ci.yml)

Quarto-based QC report tool for CosMx spatial-transcriptomics flatfiles.

Mirrors the 13-section structure of the [SizunJiangLab CosMx QC tutorial](https://github.com/SizunJiangLab/Tutorials/tree/main/CosMx/02_QC), rendered in Python via Quarto. Single command, multiple samples, one self-contained HTML.

## Requirements

- Python ≥ 3.10
- [Quarto CLI](https://quarto.org/docs/get-started/) ≥ 1.4 on `PATH`

## Install

```bash
git clone git@github.com:wuwenrui555/cosmx-qc.git
cd cosmx-qc
pip install -e '.[dev]'
```

## Usage

### Expected sample layout

The canonical CosMx AtoMx export looks like this (server-side, as exported):

```text
<sample>/
└── AtoMx/
    └── flatFiles/
        └── <run>/
            ├── <run>_exprMat_file.csv.gz
            ├── <run>_metadata_file.csv.gz
            ├── <run>_tx_file.csv.gz
            └── <run>_fov_positions_file.csv.gz
```

You can hand `cosmx-qc` either of these two paths and it will resolve to the same flatfiles. It logs which step matched.

| Input | Resolves via |
| - | - |
| The sample root, `<sample>/` | descends `AtoMx/flatFiles/<run>/` |
| The full flatfiles directory, `<sample>/AtoMx/flatFiles/<run>/` | used directly |

(There is also a fallback that looks one level deep when the four `*.csv.gz` live in an immediate subdirectory of the input — e.g. a hand-organized `<sample>/<some_name>/`. The two layouts above are the recommended ones.)

### Inline samples

```bash
cosmx-qc report \
    --sample data_1=/path/to/data_1 \
    --sample data_2=/path/to/data_2 \
    --output qc.html
```

Each `/path/to/data_*` is a sample root containing `AtoMx/flatFiles/<run>/...`.

### Config file

```yaml
# samples.yaml
samples:
  data_1: /path/to/data_1
  data_2: /path/to/data_2
title: "CosMx QC"
```

```bash
cosmx-qc report --config samples.yaml --output qc.html
```

### Fast iteration

```bash
cosmx-qc report --sample data_1=/path/to/data_1 --n-rows 100000 -o head.html
```

`--n-rows N` reads only the first `N` rows of the two large flatfiles (`exprMat` and `tx`) via streaming gzip. Use it for quick debugging — adjusting plot styles, sanity-checking a new sample layout, etc. — without waiting for the full data to load.

### Save plot data

```bash
cosmx-qc report --sample data_1=/path/to/data_1 --save-data ./qc_data/ -o qc.html
```

`--save-data DIR` writes the underlying DataFrame for each report section as Parquet, organised as `DIR/<sample>/<metric>.parquet` (e.g. `qc_data/data_1/assignment_ratio_per_fov.parquet`). The mean-expression heatmap is split into three files per sample (`mean_expression_gene.parquet`, `_negative.parquet`, `_falsecode.parquet`). Useful for downstream analysis or custom plotting without re-running the metric computations.

### Limit thread usage

```bash
cosmx-qc report --sample data_1=/path/to/data_1 --threads 4 -o qc.html
```

`--threads N` (default `4`) caps the worker pools that the render uses internally — `polars` for CSV reading, OpenBLAS / MKL via `numpy` for groupby/sum aggregations, and OpenMP. Without this cap, those libraries default to one thread per core and can saturate a shared machine. The cap is applied to the Quarto subprocess environment only; your shell's existing thread settings are untouched.

## Development

```bash
pip install -e '.[dev]'
pre-commit install
pytest tests/
```

`pre-commit install` wires up the hooks in `.pre-commit-config.yaml`: ruff (lint + format) for Python and markdownlint-cli for Markdown. Run them ad-hoc with `pre-commit run --all-files`.

Specs and plans live in `docs/specs/` and `docs/plans/` — they document the original design rationale.

## License

Apache 2.0 — see [LICENSE](LICENSE).
