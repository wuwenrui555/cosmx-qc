import json
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
