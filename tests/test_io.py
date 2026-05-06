from pathlib import Path

import pytest

from cosmx_qc import io as qc_io

FIXTURE = Path(__file__).parent / "fixtures" / "mini_sample"


def test_validate_sample_dir_returns_empty_for_complete_fixture():
    missing = qc_io.validate_sample_dir(FIXTURE)
    assert missing == []


def test_validate_sample_dir_lists_missing(tmp_path):
    # only drop one of the four; validate flags it
    (tmp_path / "RNA_x_exprMat_file.csv.gz").write_bytes(b"")
    (tmp_path / "RNA_x_metadata_file.csv.gz").write_bytes(b"")
    (tmp_path / "RNA_x_tx_file.csv.gz").write_bytes(b"")
    missing = qc_io.validate_sample_dir(tmp_path)
    assert "*_fov_positions_file.csv.gz" in missing


def test_validate_sample_dir_rejects_nonexistent_dir(tmp_path):
    bogus = tmp_path / "does_not_exist"
    with pytest.raises(FileNotFoundError):
        qc_io.validate_sample_dir(bogus)


def test_load_sample_returns_complete_data():
    sd = qc_io.load_sample("mini", FIXTURE)
    assert sd.name == "mini"
    assert {"fov", "cell_ID", "GeneA", "Negative1", "SystemControl1"}.issubset(sd.exprmat.columns)
    assert {"fov", "cell_ID", "Area.um2"}.issubset(sd.metadata.columns)
    # tx is projected to the 3 columns we actually use; CellComp/x_*/y_* dropped
    assert set(sd.tx.columns) == {"fov", "cell_ID", "target"}
    assert {"FOV", "X_mm", "Y_mm"}.issubset(sd.fov_pos.columns)
    assert len(sd.exprmat) == 5
    assert len(sd.tx) == 13


def test_load_sample_normalizes_position_columns():
    sd = qc_io.load_sample("mini", FIXTURE)
    # x_global_mm/y_global_mm should be renamed to X_mm/Y_mm
    assert "x_global_mm" not in sd.fov_pos.columns
    assert "X_mm" in sd.fov_pos.columns
    assert "Y_mm" in sd.fov_pos.columns
