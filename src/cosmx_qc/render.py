"""Quarto render orchestration: write JSON config to tempfile, call quarto."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime
from importlib.resources import files
from pathlib import Path


def render_report(
    samples: dict[str, Path],
    output: Path,
    title: str = "CosMx QC Report",
    n_rows: int | None = None,
    save_data: Path | None = None,
    threads: int = 4,
) -> None:
    """Render report.qmd to `output` HTML via Quarto.

    The .qmd is copied into a temp dir and Quarto runs there with
    cwd=temp_dir so all aux files (`*_files/libs/...`) land alongside
    the output. The final HTML is then moved to the user's path.

    When `save_data` is given, the report also writes each section's
    plot data to `<save_data>/<sample>/<metric>.parquet`.

    `threads` caps the worker pools that report execution spawns
    (polars CSV reader, OpenBLAS / MKL via numpy, OpenMP) so the render
    does not claim every core on a shared machine. Set in the Quarto
    subprocess env only — the parent shell's environment is untouched.
    """
    output = Path(output).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    if save_data is not None:
        save_data = Path(save_data).resolve()
        save_data.mkdir(parents=True, exist_ok=True)

    qmd_src = files("cosmx_qc").joinpath("report.qmd")
    payload = {
        "samples": {k: str(v) for k, v in samples.items()},
        "title": title,
        "n_rows": n_rows,
        "generated_at": datetime.now().isoformat(timespec="seconds"),
    }

    log_path = output.with_suffix(".log")
    print(f"[cosmx-qc] Progress log: {log_path}", file=sys.stderr, flush=True)
    print(f"[cosmx-qc]   tail -f {log_path}", file=sys.stderr, flush=True)

    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as f:
        json.dump(payload, f)
        cfg_path = f.name
    try:
        env = {
            **os.environ,
            "COSMX_QC_CONFIG": cfg_path,
            "COSMX_QC_LOG": str(log_path),
            # Standard env-var names that polars / OpenBLAS / MKL / OpenMP
            # themselves read on import. They cannot be renamed — these are
            # the keys the libraries look for. Set on the child env only,
            # so the parent shell's thread settings are not modified.
            "POLARS_MAX_THREADS": str(threads),
            "OPENBLAS_NUM_THREADS": str(threads),
            "MKL_NUM_THREADS": str(threads),
            "OMP_NUM_THREADS": str(threads),
        }
        if save_data is not None:
            env["COSMX_QC_SAVE_DATA"] = str(save_data)
        with tempfile.TemporaryDirectory(prefix="cosmx_qc_") as render_dir:
            render_dir = Path(render_dir)
            qmd_local = render_dir / "report.qmd"
            shutil.copy(str(qmd_src), str(qmd_local))
            subprocess.run(
                ["quarto", "render", "report.qmd", "--to", "html", "-M", f"title={title}"],
                env=env,
                cwd=str(render_dir),
                check=True,
            )
            shutil.move(str(render_dir / "report.html"), str(output))
    finally:
        os.unlink(cfg_path)
