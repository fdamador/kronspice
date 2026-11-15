# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""PDF report for examples/void_shells: void shells (H50)

Run from the repository root:

    python examples/void_shells/viz/make_report.py examples/void_shells/void_shells_report.pdf

The body is examples/void_shells/README.md rendered as written (internal/pdfreport/mdreport.py),
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
    r = Report("examples/void_shells", "Void shells: light pushing out, ether pushing in, built up in shells", "the round light bubble, excited modes, prepared shells and onions, several photons, and the g bound", "proposed register entry H50 (examples/void_shells/CLAUDE.md)")
    used = render_markdown(r, os.path.join(EX, "README.md"), base=EX)
    add_figures(r, os.path.join(EX, "viz"), skip=used)
    add_run_excerpt(r, os.path.join(EX, "run_output.txt"))
    r.h1("Reproducing the results")
    r.code("go test ./examples/void_shells\ngo run ./examples/void_shells -viz examples/void_shells/viz\npython3 examples/void_shells/viz/plot_shells.py")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(EX, "void_shells_report.pdf"))
