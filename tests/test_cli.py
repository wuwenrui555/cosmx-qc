from pathlib import Path
from unittest.mock import patch

import pytest
import yaml

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


def test_load_config_file(tmp_path):
    cfg = tmp_path / "samples.yaml"
    cfg.write_text(
        yaml.safe_dump(
            {
                "samples": {"A": "/p1", "B": "/p2"},
                "title": "X",
            }
        )
    )
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
    for pat in [
        "x_exprMat_file.csv.gz",
        "x_metadata_file.csv.gz",
        "x_tx_file.csv.gz",
        "x_fov_positions_file.csv.gz",
    ]:
        (d / pat).write_bytes(b"")
    # should not raise
    C.validate_samples({"s1": d})


def test_validate_samples_missing_files(tmp_path):
    d = tmp_path / "s1"
    d.mkdir()
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
    args = parser.parse_args(["report", "--sample", "A=/a", "--sample", "B=/b", "-o", "out.html"])
    assert args.samples == [("A", Path("/a")), ("B", Path("/b"))]


def test_duplicate_sample_name_rejected():
    args_samples = [("A", Path("/a")), ("A", Path("/b"))]
    with pytest.raises(SystemExit):
        C.dedupe_samples(args_samples)


def test_cmd_report_inline_samples(tmp_path):
    # build a valid sample dir
    d = tmp_path / "S"
    d.mkdir()
    for pat in [
        "x_exprMat_file.csv.gz",
        "x_metadata_file.csv.gz",
        "x_tx_file.csv.gz",
        "x_fov_positions_file.csv.gz",
    ]:
        (d / pat).write_bytes(b"")
    parser = C.build_parser()
    args = parser.parse_args(["report", "--sample", f"S={d}", "-o", str(tmp_path / "r.html")])
    with (
        patch("cosmx_qc.cli.check_quarto_installed", return_value="/usr/bin/quarto"),
        patch("cosmx_qc.cli.render_report") as rr,
    ):
        C.cmd_report(args)
    rr.assert_called_once()
    kwargs = rr.call_args.kwargs or rr.call_args[1]
    assert kwargs["samples"] == {"S": d}
    assert kwargs["output"] == tmp_path / "r.html"
    assert kwargs["title"] == "CosMx QC Report"


def test_cmd_report_config_file(tmp_path):
    d = tmp_path / "S"
    d.mkdir()
    for pat in [
        "x_exprMat_file.csv.gz",
        "x_metadata_file.csv.gz",
        "x_tx_file.csv.gz",
        "x_fov_positions_file.csv.gz",
    ]:
        (d / pat).write_bytes(b"")
    cfg = tmp_path / "c.yaml"
    cfg.write_text(yaml.safe_dump({"samples": {"S": str(d)}, "title": "Cfg"}))
    parser = C.build_parser()
    args = parser.parse_args(["report", "--config", str(cfg), "-o", str(tmp_path / "r.html")])
    with (
        patch("cosmx_qc.cli.check_quarto_installed", return_value="q"),
        patch("cosmx_qc.cli.render_report") as rr,
    ):
        C.cmd_report(args)
    kwargs = rr.call_args.kwargs or rr.call_args[1]
    assert kwargs["samples"] == {"S": d}
    assert kwargs["title"] == "Cfg"


def test_argparse_save_data_parses(tmp_path):
    parser = C.build_parser()
    args = parser.parse_args(
        ["report", "--sample", "A=/a", "--save-data", str(tmp_path), "-o", "out.html"]
    )
    assert args.save_data == tmp_path


def test_cmd_report_passes_save_data(tmp_path):
    d = tmp_path / "S"
    d.mkdir()
    for pat in [
        "x_exprMat_file.csv.gz",
        "x_metadata_file.csv.gz",
        "x_tx_file.csv.gz",
        "x_fov_positions_file.csv.gz",
    ]:
        (d / pat).write_bytes(b"")
    save_dir = tmp_path / "qc_data"
    parser = C.build_parser()
    args = parser.parse_args(
        [
            "report",
            "--sample",
            f"S={d}",
            "--save-data",
            str(save_dir),
            "-o",
            str(tmp_path / "r.html"),
        ]
    )
    with (
        patch("cosmx_qc.cli.check_quarto_installed", return_value="q"),
        patch("cosmx_qc.cli.render_report") as rr,
    ):
        C.cmd_report(args)
    kwargs = rr.call_args.kwargs or rr.call_args[1]
    assert kwargs["save_data"] == save_dir


def test_cmd_report_default_save_data_none(tmp_path):
    d = tmp_path / "S"
    d.mkdir()
    for pat in [
        "x_exprMat_file.csv.gz",
        "x_metadata_file.csv.gz",
        "x_tx_file.csv.gz",
        "x_fov_positions_file.csv.gz",
    ]:
        (d / pat).write_bytes(b"")
    parser = C.build_parser()
    args = parser.parse_args(
        ["report", "--sample", f"S={d}", "-o", str(tmp_path / "r.html")]
    )
    with (
        patch("cosmx_qc.cli.check_quarto_installed", return_value="q"),
        patch("cosmx_qc.cli.render_report") as rr,
    ):
        C.cmd_report(args)
    kwargs = rr.call_args.kwargs or rr.call_args[1]
    assert kwargs["save_data"] is None


def test_argparse_threads_default_4():
    parser = C.build_parser()
    args = parser.parse_args(["report", "--sample", "A=/a", "-o", "out.html"])
    assert args.threads == 4


def test_argparse_threads_explicit():
    parser = C.build_parser()
    args = parser.parse_args(
        ["report", "--sample", "A=/a", "--threads", "4", "-o", "out.html"]
    )
    assert args.threads == 4


def test_cmd_report_passes_threads(tmp_path):
    d = tmp_path / "S"
    d.mkdir()
    for pat in [
        "x_exprMat_file.csv.gz",
        "x_metadata_file.csv.gz",
        "x_tx_file.csv.gz",
        "x_fov_positions_file.csv.gz",
    ]:
        (d / pat).write_bytes(b"")
    parser = C.build_parser()
    args = parser.parse_args(
        [
            "report",
            "--sample",
            f"S={d}",
            "--threads",
            "2",
            "-o",
            str(tmp_path / "r.html"),
        ]
    )
    with (
        patch("cosmx_qc.cli.check_quarto_installed", return_value="q"),
        patch("cosmx_qc.cli.render_report") as rr,
    ):
        C.cmd_report(args)
    kwargs = rr.call_args.kwargs or rr.call_args[1]
    assert kwargs["threads"] == 2


def test_cmd_report_explicit_title_overrides_config(tmp_path):
    d = tmp_path / "S"
    d.mkdir()
    for pat in [
        "x_exprMat_file.csv.gz",
        "x_metadata_file.csv.gz",
        "x_tx_file.csv.gz",
        "x_fov_positions_file.csv.gz",
    ]:
        (d / pat).write_bytes(b"")
    cfg = tmp_path / "c.yaml"
    cfg.write_text(yaml.safe_dump({"samples": {"S": str(d)}, "title": "FromCfg"}))
    parser = C.build_parser()
    args = parser.parse_args(
        ["report", "--config", str(cfg), "--title", "FromCLI", "-o", str(tmp_path / "r.html")]
    )
    with (
        patch("cosmx_qc.cli.check_quarto_installed", return_value="q"),
        patch("cosmx_qc.cli.render_report") as rr,
    ):
        C.cmd_report(args)
    kwargs = rr.call_args.kwargs or rr.call_args[1]
    assert kwargs["title"] == "FromCLI"
