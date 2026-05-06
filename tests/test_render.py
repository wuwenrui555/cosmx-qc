import json
import os
import subprocess
from pathlib import Path
from unittest.mock import patch

from cosmx_qc import render as R


def _fake_quarto_run(captured: dict):
    """Return a fake subprocess.run that writes report.html into the
    cwd Quarto was launched in (the temp render dir)."""

    def fake_run(cmd, env=None, check=False, cwd=None, **kw):
        cfg_path = env["COSMX_QC_CONFIG"]
        captured["cmd"] = cmd
        captured["cwd"] = cwd
        captured["cfg"] = json.loads(Path(cfg_path).read_text())
        # render.py copies report.qmd into cwd before invoking quarto;
        # quarto would normally produce report.html alongside it.
        (Path(cwd) / "report.html").write_text("<html></html>")
        return subprocess.CompletedProcess(cmd, 0)

    return fake_run


def test_render_report_writes_config_and_calls_quarto(tmp_path):
    output = tmp_path / "out.html"
    captured: dict = {}

    with patch("cosmx_qc.render.subprocess.run", side_effect=_fake_quarto_run(captured)):
        R.render_report(
            samples={"A": tmp_path / "a", "B": tmp_path / "b"},
            output=output,
            title="T",
        )

    assert captured["cmd"][0] == "quarto"
    assert "render" in captured["cmd"]
    assert captured["cfg"]["samples"] == {"A": str(tmp_path / "a"), "B": str(tmp_path / "b")}
    assert captured["cfg"]["title"] == "T"
    assert "-M" in captured["cmd"]
    title_idx = captured["cmd"].index("-M") + 1
    assert captured["cmd"][title_idx] == "title=T"
    # quarto runs with cwd=temp render dir; report.qmd is the bare basename
    assert captured["cwd"] is not None
    assert captured["cmd"][2] == "report.qmd"
    # final HTML moved to user's requested location
    assert output.exists()


def test_render_report_creates_output_parent(tmp_path):
    output = tmp_path / "deeper" / "report.html"
    with patch("cosmx_qc.render.subprocess.run", side_effect=_fake_quarto_run({})):
        R.render_report(samples={"A": tmp_path}, output=output, title="X")
    assert output.parent.exists()
    assert output.exists()


def _capture_env_quarto_run(captured: dict):
    """Like `_fake_quarto_run` but captures the env passed to subprocess."""

    def fake_run(cmd, env=None, check=False, cwd=None, **kw):
        captured["env"] = env
        (Path(cwd) / "report.html").write_text("<html></html>")
        return subprocess.CompletedProcess(cmd, 0)

    return fake_run


_THREAD_ENV_KEYS = (
    "POLARS_MAX_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "OMP_NUM_THREADS",
)


def test_render_report_passes_save_data_env_var(tmp_path):
    output = tmp_path / "out.html"
    save_dir = tmp_path / "qc_data"
    captured: dict = {}
    with patch(
        "cosmx_qc.render.subprocess.run", side_effect=_capture_env_quarto_run(captured)
    ):
        R.render_report(
            samples={"A": tmp_path / "a"},
            output=output,
            title="T",
            save_data=save_dir,
        )
    assert captured["env"]["COSMX_QC_SAVE_DATA"] == str(save_dir.resolve())
    assert save_dir.exists()


def test_render_report_omits_save_data_env_var_by_default(tmp_path):
    output = tmp_path / "out.html"
    captured: dict = {}
    with patch(
        "cosmx_qc.render.subprocess.run", side_effect=_capture_env_quarto_run(captured)
    ):
        R.render_report(samples={"A": tmp_path / "a"}, output=output, title="T")
    assert "COSMX_QC_SAVE_DATA" not in captured["env"]


def test_render_report_thread_env_vars_default_4(tmp_path):
    output = tmp_path / "out.html"
    captured: dict = {}
    with patch(
        "cosmx_qc.render.subprocess.run", side_effect=_capture_env_quarto_run(captured)
    ):
        R.render_report(samples={"A": tmp_path / "a"}, output=output, title="T")
    for k in _THREAD_ENV_KEYS:
        assert captured["env"][k] == "4"


def test_render_report_thread_env_vars_explicit(tmp_path):
    output = tmp_path / "out.html"
    captured: dict = {}
    with patch(
        "cosmx_qc.render.subprocess.run", side_effect=_capture_env_quarto_run(captured)
    ):
        R.render_report(
            samples={"A": tmp_path / "a"}, output=output, title="T", threads=4
        )
    for k in _THREAD_ENV_KEYS:
        assert captured["env"][k] == "4"


def test_render_report_thread_env_does_not_leak_to_parent(tmp_path):
    """Setting threads inside the subprocess env must not mutate os.environ."""
    before = {k: os.environ.get(k) for k in _THREAD_ENV_KEYS}
    with patch(
        "cosmx_qc.render.subprocess.run", side_effect=_capture_env_quarto_run({})
    ):
        R.render_report(
            samples={"A": tmp_path / "a"},
            output=tmp_path / "out.html",
            title="T",
            threads=2,
        )
    after = {k: os.environ.get(k) for k in _THREAD_ENV_KEYS}
    assert before == after
