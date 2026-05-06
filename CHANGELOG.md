# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - 2026-05-06

### Added

- `--save-data DIR` CLI flag: also write each report section's underlying DataFrame as Parquet, organised as `DIR/<sample>/<metric>.parquet`. The mean-expression heatmap is split into three files per sample (gene / negative / falsecode).
- `--threads N` CLI flag (default `4`): cap the worker pools the render uses internally (`polars`, `OpenBLAS` / `MKL` via `numpy`, OpenMP). Set on the Quarto subprocess env only; the parent shell is untouched.
- Per-section description paragraphs in the QC report (each section now explains what its data is).
- README: previously-undocumented `--output` and `--title` flags now have their own subsections, including title resolution order.

### Changed

- Report section title `Unique transcripts — panel detection ratio` → `Panel detection ratio`.
- Report section title `Unique transcripts — expression heatmap` → `Mean expression per target (per FOV)`. The previous title implied uniqueness; the metric has always been per-FOV mean expression per target. `expression_heatmap_data` docstring updated to match.
- README: replaced the numbered usage subsections (`Inline samples` / `Config file` / `Fast iteration` / `Save plot data` / `Limit thread usage`) with a single `Flags` reference section plus a POSIX-style synopsis. Sample-layout section now states explicitly that the `AtoMx/` subtree is the raw export and documents the lab data standard for the sample-root directory name.

### Removed

- **Breaking**: `io.resolve_flatfiles_dir` no longer does the one-level subdirectory fallback. It now accepts only the flatfiles directory itself or an AtoMx-style sample root containing `AtoMx/flatFiles/<run>/`. Hand-organised layouts that relied on the fallback need to be reorganised.

## [0.1.1] - 2026-05-06

### Fixed

- CI: markdownlint now skips `docs/specs/` and `docs/plans/`. Those directories hold historical design artefacts (the original spec and implementation plan) and are not actively maintained against current lint rules; making them gate the build was a CI footgun.

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
