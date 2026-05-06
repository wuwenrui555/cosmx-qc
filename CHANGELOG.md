# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2026-05-06

### Added

- `cosmx-qc` CLI: `report` subcommand with mutually-exclusive `--sample NAME=PATH` (repeatable) and `--config <yaml>`, plus `--output`, `--title`, and `--n-rows` for fast iteration.
- Auto-resolver in `io.resolve_flatfiles_dir` that accepts the flatfiles directory itself, an AtoMx-style sample root (descends into `AtoMx/flatFiles/<run>/`), or a one-level parent.
- Streaming gzip reader for `--n-rows` mode (polars cannot break out of compressed inputs early).
- 13 pure DataFrame metrics aligned with the SizunJiangLab CosMx 02_QC R tutorial: `transcripts_per_fov`, `unassigned_transcripts_per_fov`, `assignment_ratio_per_fov`, `unique_transcripts_per_fov`, `panel_detection_ratio_per_fov`, `transcripts_breakdown_per_fov`, `expression_heatmap_data`, `cell_count_per_fov`, `cell_area_per_fov`, `error_rates_per_fov`, `transcripts_per_cell_per_fov`, `error_rates_per_cell_per_fov`, `global_view_table`.
- `SampleData.tx_kind` cached property so transcript classification runs once per sample.
- Plot module matching the original Python QC styles: 16x8 bars/box, 12x10 FOV positions, Gene #E74C3C / Negative #5DADE2 / SystemControl #52C785, neutral grey #A0A0A0.
- Combined view as a single matplotlib Figure with `sharey=True` plus per-sample tabs; PIL hstack fallback for FOV positions and the 3-panel expression heatmap.
- Stacked bar variants for Unassigned (Assigned/Unassigned) and Detection ratio (Detected/Undetected) with a shared `plot_stacked_two_per_fov` helper.
- Quarto orchestration via `render.render_report`: writes a JSON config tempfile, hands off via `COSMX_QC_CONFIG`, copies the packaged `report.qmd` into a temp dir, runs `quarto render` with `cwd=temp_dir` so auxiliary files land alongside the output, then `shutil.move`s the HTML back to `--output`.
- Per-metric `_timed` context manager and per-section plot-emit timing; final summary chunk logs total render time + load/sections split + peak RSS.
- Pre-commit hooks: ruff (lint --fix + format), markdownlint-cli; CI workflow on Python 3.10 / 3.11 / 3.12.
- 46 unit tests across io, metrics, plots, cli, render.
- Apache 2.0 license.
