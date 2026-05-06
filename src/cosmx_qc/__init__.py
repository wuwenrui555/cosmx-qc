"""CosMx quality-control report tooling.

Public CLI entry: `cosmx-qc report --sample NAME=PATH ... -o report.html`
Programmatic entry: `cosmx_qc.render.render_report`.
"""

from cosmx_qc.io import SampleData, load_sample, validate_sample_dir
from cosmx_qc.render import render_report

__all__ = ["SampleData", "load_sample", "validate_sample_dir", "render_report"]
