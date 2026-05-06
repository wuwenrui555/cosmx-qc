import pandas as pd
import pytest

from cosmx_qc import metrics as M


def test_classify_columns(mini_sample):
    cls = M.classify_expr_columns(mini_sample.exprmat)
    assert cls["gene"] == ["GeneA"]
    assert cls["negative"] == ["Negative1"]
    assert cls["falsecode"] == ["SystemControl1"]


def test_classify_targets(mini_sample):
    s = mini_sample.tx["target"]
    assert M.classify_target(s).tolist() == [
        "gene",
        "negative",
        "gene",
        "falsecode",
        "gene",
        "gene",
        "gene",
        "negative",
        "gene",
        "gene",
        "gene",
        "negative",
        "falsecode",
    ]


def test_transcripts_per_fov_all(mini_sample):
    df = M.transcripts_per_fov(mini_sample, kind="all")
    # 13 transcripts total: 6 in fov 1, 7 in fov 2
    expected = pd.DataFrame({"fov": [1, 2], "count": [6, 7]})
    pd.testing.assert_frame_equal(df.reset_index(drop=True), expected)


def test_transcripts_per_fov_gene(mini_sample):
    df = M.transcripts_per_fov(mini_sample, kind="gene")
    # gene tx in fov 1: 4 (rows with target=GeneA at indices 0,2,4,10)
    # gene tx in fov 2: 4 (indices 5,6,8,9)
    expected = pd.DataFrame({"fov": [1, 2], "count": [4, 4]})
    pd.testing.assert_frame_equal(df.reset_index(drop=True), expected)


def test_transcripts_per_fov_negative(mini_sample):
    df = M.transcripts_per_fov(mini_sample, kind="negative")
    # negative tx in fov 1: 1 (index 1); fov 2: 2 (indices 7, 11)
    expected = pd.DataFrame({"fov": [1, 2], "count": [1, 2]})
    pd.testing.assert_frame_equal(df.reset_index(drop=True), expected)


def test_transcripts_per_fov_falsecode(mini_sample):
    df = M.transcripts_per_fov(mini_sample, kind="falsecode")
    # falsecode (SystemControl1) in fov 1: 1 (index 3); fov 2: 1 (index 12)
    expected = pd.DataFrame({"fov": [1, 2], "count": [1, 1]})
    pd.testing.assert_frame_equal(df.reset_index(drop=True), expected)


def test_unassigned_transcripts_per_fov(mini_sample):
    df = M.unassigned_transcripts_per_fov(mini_sample)
    # cell_ID == 0 transcripts: 2 in fov 1 (indices 4, 10); 3 in fov 2 (indices 9, 11, 12)
    expected = pd.DataFrame({"fov": [1, 2], "unassigned": [2, 3]})
    pd.testing.assert_frame_equal(df.reset_index(drop=True), expected)


def test_unique_transcripts_per_fov(mini_sample):
    df = M.unique_transcripts_per_fov(mini_sample)
    # fov 1 sees {GeneA, Negative1, SystemControl1} = 3
    # fov 2 sees {GeneA, Negative1, SystemControl1} = 3
    expected = pd.DataFrame({"fov": [1, 2], "unique": [3, 3]})
    pd.testing.assert_frame_equal(df.reset_index(drop=True), expected)


def test_cell_count_per_fov(mini_sample):
    df = M.cell_count_per_fov(mini_sample)
    expected = pd.DataFrame({"fov": [1, 2], "n_cells": [3, 2]})
    pd.testing.assert_frame_equal(df.reset_index(drop=True), expected)


def test_cell_area_per_fov(mini_sample):
    df = M.cell_area_per_fov(mini_sample)
    # all rows passed through, with fov + Area.um2
    assert set(df.columns) == {"fov", "Area.um2"}
    assert len(df) == 5
    assert df.loc[df["fov"] == 1, "Area.um2"].mean() == pytest.approx(150.0)


def test_error_rates_per_fov(mini_sample):
    df = M.error_rates_per_fov(mini_sample)
    # fov 1 tx kinds: gene=4, negative=1, falsecode=1
    #   neg_rate = 1/4 = 0.25 ; false_code_rate = 1/4 = 0.25
    # fov 2 tx kinds: gene=4, negative=2, falsecode=1
    #   neg_rate = 2/4 = 0.5  ; false_code_rate = 1/4 = 0.25
    assert df.loc[df["fov"] == 1, "negative_rate"].iloc[0] == pytest.approx(0.25)
    assert df.loc[df["fov"] == 1, "false_code_rate"].iloc[0] == pytest.approx(0.25)
    assert df.loc[df["fov"] == 2, "negative_rate"].iloc[0] == pytest.approx(0.5)
    assert df.loc[df["fov"] == 2, "false_code_rate"].iloc[0] == pytest.approx(0.25)


def test_global_view_table(mini_sample):
    df = M.global_view_table(mini_sample)
    # one row per slide_ID
    assert len(df) == 1
    row = df.iloc[0]
    assert row["slide_ID"] == "slide_A"
    assert row["n_fovs"] == 2
    assert row["n_cells"] == 5
    assert row["n_transcripts"] == 13


def test_transcripts_per_cell_per_fov_all(mini_sample):
    df = M.transcripts_per_cell_per_fov(mini_sample, kind="all")
    # cells: per-cell sum of all expression cols
    # fov 1, cell 1: 10+1+0=11
    # fov 1, cell 2: 20+0+1=21
    # fov 1, cell 3: 30+2+0=32
    # fov 2, cell 1: 5+0+1=6
    # fov 2, cell 2: 15+1+0=16
    assert set(df.columns) >= {"fov", "cell_ID", "count"}
    assert df.loc[(df["fov"] == 1) & (df["cell_ID"] == 3), "count"].iloc[0] == 32


def test_transcripts_per_cell_per_fov_gene(mini_sample):
    df = M.transcripts_per_cell_per_fov(mini_sample, kind="gene")
    # only GeneA values
    expected_counts = {(1, 1): 10, (1, 2): 20, (1, 3): 30, (2, 1): 5, (2, 2): 15}
    for (fov, cid), v in expected_counts.items():
        assert df.loc[(df["fov"] == fov) & (df["cell_ID"] == cid), "count"].iloc[0] == v


def test_error_rates_per_cell_per_fov(mini_sample):
    df = M.error_rates_per_cell_per_fov(mini_sample)
    # per cell: neg_rate = neg_count / gene_count, false_code_rate = sc / gene
    assert set(df.columns) >= {"fov", "cell_ID", "negative_rate", "false_code_rate"}
    # fov 1 cell 1: GeneA=10, Neg=1, SC=0 -> 0.1, 0.0
    row = df[(df["fov"] == 1) & (df["cell_ID"] == 1)].iloc[0]
    assert row["negative_rate"] == pytest.approx(0.1)
    assert row["false_code_rate"] == pytest.approx(0.0)
