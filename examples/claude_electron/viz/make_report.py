# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""PDF report for examples/claude_electron: the Claude electron (H44)

Run from the repository root:

    python examples/claude_electron/viz/make_report.py examples/claude_electron/claude_electron_report.pdf

The body is examples/claude_electron/README.md rendered as written (internal/pdfreport/mdreport.py),
followed by any figures in viz/ it does not show and the head of run_output.txt.
Nothing is re-simulated. Layout from internal/pdfreport/kronreport.py.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
EX = os.path.normpath(os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report  # noqa: E402
from mdreport import add_figures, add_run_excerpt, render_markdown  # noqa: E402


def build(out):
    r = Report("examples/claude_electron", "The Claude electron: the toroidal photon held by a knotted ether", "H23's toroidal photon held by H29's nonlinear ether as an axisymmetric Faddeev knot", "register entry H44")
    used = render_markdown(r, os.path.join(EX, "README.md"), base=EX)
    add_figures(r, os.path.join(EX, "viz"), skip=used)
    add_run_excerpt(r, os.path.join(EX, "run_output.txt"))
    r.h1("Reproducing the results")
    r.code("go test ./examples/claude_electron\ngo run ./examples/claude_electron")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(EX, "claude_electron_report.pdf"))
