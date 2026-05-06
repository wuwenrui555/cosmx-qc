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

The `AtoMx/` subtree is the raw export from AtoMx — `cosmx-qc` does not modify it and assumes it follows that exact structure. The sample-root directory *above* `flatFiles/` is named per our lab's data standard:

```text
<sample> = <YYYYMMDD>_<user>_<project>_<custom_name>_<atomx_version>
```

Putting it together:

```text
<sample>/
└── AtoMx/            ← everything under here is the raw AtoMx export, untouched
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
| The sample root, `<sample>/` | descends `<sample>/AtoMx/flatFiles/<run>/` |
| The full flatfiles directory, `<sample>/AtoMx/flatFiles/<run>/` | used directly |

### Flags

A `report` invocation needs samples (via `--sample` or `--config`); everything else is optional. Synopsis:

```text
cosmx-qc report (--sample NAME=PATH... | --config FILE)
                [-o FILE] [--title STRING] [--n-rows N]
                [--save-data DIR] [--threads N]
```

#### `--sample NAME=PATH` (repeatable) / `--config FILE`

Mutually exclusive; one is required. Passing both causes argparse to reject the call.

`--sample` is convenient for one-off runs; repeat the flag once per sample. `NAME` is the label shown in the report tabs; `PATH` is either of the two inputs accepted in the [Expected sample layout](#expected-sample-layout) table above.

```bash
cosmx-qc report \
    --sample data_1=/path/to/data_1 \
    --sample data_2=/path/to/data_2 \
    -o qc.html
```

`--config FILE` reads the same information from a YAML file. Useful when you have many samples or want a reproducible invocation.

```yaml
# samples.yaml
samples:
  data_1: /path/to/data_1
  data_2: /path/to/data_2
title: "CosMx QC Report"
```

```bash
cosmx-qc report --config samples.yaml -o qc.html
```

#### `--output FILE` / `-o FILE`

Path to write the rendered HTML. Must end in `.html`. Default: `report.html` in the current directory. The Quarto progress log is written to the same path with `.log` extension.

#### `--title STRING`

Override the report title shown in the HTML header. Resolution order, highest priority first:

1. `--title` on the CLI
2. `title:` field in the YAML config (only when `--config` is used)
3. `"CosMx QC Report"` (built-in default)

#### `--n-rows N`

Read only the first `N` rows of the two large flatfiles (`exprMat` and `tx`) via streaming gzip. Use it for quick debugging — adjusting plot styles, sanity-checking a new sample layout, etc. — without waiting for the full data to load.

```bash
cosmx-qc report --sample data_1=/path/to/data_1 --n-rows 100000 -o head.html
```

#### `--save-data DIR`

Also write the underlying DataFrame for each report section as [Parquet](https://parquet.apache.org/), organized as `DIR/<sample>/<metric>.parquet` (e.g. `qc_data/data_1/assignment_ratio_per_fov.parquet`).  Useful for downstream analysis or custom plotting without re-running the metric computations.

```bash
cosmx-qc report --sample data_1=/path/to/data_1 --save-data ./qc_data/ -o qc.html
```

#### `--threads N`

Cap the worker pools that the render uses internally — `polars` for CSV reading, OpenBLAS / MKL via `numpy` for groupby/sum aggregations, and OpenMP. Without this cap, those libraries default to one thread per core and can saturate a shared machine. Default: `8`. Applied to the Quarto subprocess environment only; your shell's existing thread settings are untouched.

```bash
cosmx-qc report --sample data_1=/path/to/data_1 --threads 16 -o qc.html
```

## Development

```bash
pip install -e '.[dev]'
pre-commit install
pytest tests/
```

`pre-commit install` wires up the hooks in `.pre-commit-config.yaml`: ruff (lint + format) for Python and markdownlint-cli for Markdown. Run them ad-hoc with `pre-commit run --all-files`.

## License

Apache 2.0 — see [LICENSE](LICENSE).
