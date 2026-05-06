"""Smoke tests for plots module: each function returns a Figure without raising."""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from cosmx_qc import metrics as M
from cosmx_qc import plots as P


def test_plot_fov_positions(mini_sample):
    fig = P.plot_fov_positions(mini_sample.fov_pos, title="mini")
    assert isinstance(fig, plt.Figure)
    plt.close(fig)


def test_plot_count_per_fov(mini_sample):
    df = M.transcripts_per_fov(mini_sample, kind="gene")
    fig = P.plot_count_per_fov(
        df, y_col="count", ylabel="Count", color=P.COLORS["gene"], title="mini"
    )
    assert isinstance(fig, plt.Figure)
    plt.close(fig)


def test_plot_assignment_per_fov(mini_sample):
    df = M.assignment_ratio_per_fov(mini_sample)
    fig = P.plot_assignment_per_fov(df, title="mini")
    assert isinstance(fig, plt.Figure)
    plt.close(fig)


def test_plot_stacked_kinds_per_fov_count(mini_sample):
    df = M.transcripts_breakdown_per_fov(mini_sample)
    fig = P.plot_stacked_kinds_per_fov(df, title="mini", percent=False)
    assert isinstance(fig, plt.Figure)
    plt.close(fig)


def test_plot_stacked_kinds_per_fov_percent(mini_sample):
    df = M.transcripts_breakdown_per_fov(mini_sample)
    fig = P.plot_stacked_kinds_per_fov(df, title="mini", percent=True)
    assert isinstance(fig, plt.Figure)
    plt.close(fig)


def test_plot_box_per_fov(mini_sample):
    df = M.cell_area_per_fov(mini_sample)
    fig = P.plot_box_per_fov(
        df, y_col="Area.um2", ylabel="Area (um^2)", color=P.COLORS["neutral"], title="mini"
    )
    assert isinstance(fig, plt.Figure)
    plt.close(fig)


def test_plot_box_per_fov_per_cell(mini_sample):
    df = M.transcripts_per_cell_per_fov(mini_sample, kind="all")
    fig = P.plot_box_per_fov(df, y_col="count", ylabel="Count", color=P.COLORS["gene"])
    assert isinstance(fig, plt.Figure)
    plt.close(fig)
