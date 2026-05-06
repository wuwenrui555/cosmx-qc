"""cosmx-qc command-line interface."""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import yaml

from cosmx_qc.io import validate_sample_dir
from cosmx_qc.render import render_report


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
    g.add_argument(
        "--sample",
        action="append",
        dest="samples",
        type=parse_sample_arg,
        metavar="NAME=PATH",
        help="Sample name and flatfile directory (repeatable)",
    )
    g.add_argument(
        "--config", type=Path, help="YAML config: {samples: {NAME: PATH, ...}, title: ...}"
    )
    rep.add_argument("--output", "-o", type=Path, default=Path("report.html"))
    rep.add_argument("--title", default=None)
    rep.add_argument(
        "--n-rows",
        type=int,
        default=None,
        dest="n_rows",
        help="Limit exprMat and tx reads to first N rows (fast iteration on Pembro-scale samples)",
    )
    rep.add_argument(
        "--save-data",
        type=Path,
        default=None,
        dest="save_data",
        metavar="DIR",
        help="Also write per-section plot data to <DIR>/<sample>/<metric>.parquet",
    )
    return p


def cmd_report(args: argparse.Namespace) -> None:
    if args.samples is not None:
        samples = dedupe_samples(args.samples)
        title = args.title if args.title is not None else "CosMx QC Report"
    else:
        samples, cfg_title = load_config(args.config)
        title = args.title if args.title is not None else cfg_title
    if not str(args.output).endswith(".html"):
        sys.exit(f"error: --output must end in '.html' (got {args.output})")
    validate_samples(samples)
    check_quarto_installed()
    render_report(
        samples=samples,
        output=args.output,
        title=title,
        n_rows=args.n_rows,
        save_data=args.save_data,
    )


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    if args.cmd == "report":
        cmd_report(args)
    else:
        parser.error(f"unknown command {args.cmd!r}")
