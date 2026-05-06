"""Generate a tiny synthetic CosMx flatfile sample for tests.

Layout: 2 FOVs, 5 cells, 1 real gene + 1 negative probe + 1 false code.
Run from repo root:  python tests/qc/fixtures/build_mini_sample.py
"""

from pathlib import Path

import pandas as pd

OUT = Path(__file__).parent / "mini_sample"
OUT.mkdir(parents=True, exist_ok=True)

# 5 cells across 2 FOVs
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

# 13 transcripts, some unassigned (cell_ID = 0)
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

fov_pos = pd.DataFrame(
    {
        "FOV": [1, 2],
        "x_global_mm": [0.0, 5.0],
        "y_global_mm": [0.0, 0.0],
    }
)

prefix = "RNA_mini"
exprmat.to_csv(OUT / f"{prefix}_exprMat_file.csv.gz", index=False, compression="gzip")
metadata.to_csv(OUT / f"{prefix}_metadata_file.csv.gz", index=False, compression="gzip")
tx.to_csv(OUT / f"{prefix}_tx_file.csv.gz", index=False, compression="gzip")
fov_pos.to_csv(OUT / f"{prefix}_fov_positions_file.csv.gz", index=False, compression="gzip")
print(f"Wrote 4 files to {OUT}")
