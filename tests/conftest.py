"""Shared pytest fixtures for qc tests."""

import pandas as pd
import pytest

from cosmx_qc.io import SampleData


@pytest.fixture
def mini_sample() -> SampleData:
    """Hand-built tiny SampleData mirroring the file fixture."""
    exprmat = pd.DataFrame(
        {
            "fov": [1, 1, 1, 2, 2],
            "cell_ID": [1, 2, 3, 1, 2],
            "GeneA": [10, 20, 30, 5, 15],
            "Negative1": [1, 0, 2, 0, 1],
            "SystemControl1": [0, 1, 0, 1, 0],
        }
    )
    metadata = pd.DataFrame(
        {
            "fov": [1, 1, 1, 2, 2],
            "cell_ID": [1, 2, 3, 1, 2],
            "Area.um2": [100.0, 200.0, 150.0, 80.0, 120.0],
            "nFeature_RNA": [3, 2, 3, 2, 2],
            "slide_ID": ["slide_A"] * 5,
            "Run_Tissue_name": ["mini"] * 5,
        }
    )
    tx = pd.DataFrame(
        {
            "fov": [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 1, 2, 2],
            "cell_ID": [1, 1, 2, 3, 0, 1, 1, 2, 2, 0, 0, 0, 0],
            "target": [
                "GeneA",
                "Negative1",
                "GeneA",
                "SystemControl1",
                "GeneA",
                "GeneA",
                "GeneA",
                "Negative1",
                "GeneA",
                "GeneA",
                "GeneA",
                "Negative1",
                "SystemControl1",
            ],
            "CellComp": ["Cytoplasm"] * 13,
            "x_local_px": [100, 110, 120, 130, 140, 200, 210, 220, 230, 240, 150, 250, 260],
            "y_local_px": [100, 110, 120, 130, 140, 200, 210, 220, 230, 240, 150, 250, 260],
        }
    )
    fov_pos = pd.DataFrame({"FOV": [1, 2], "X_mm": [0.0, 5.0], "Y_mm": [0.0, 0.0]})
    return SampleData(name="mini", exprmat=exprmat, metadata=metadata, tx=tx, fov_pos=fov_pos)
