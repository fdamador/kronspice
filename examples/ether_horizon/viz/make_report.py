# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""PDF report for examples/ether_horizon: the ether horizon (H49), Parts 1-11

Run from the repository root:

    python examples/ether_horizon/viz/make_report.py examples/ether_horizon/ether_horizon_report.pdf

The body is examples/ether_horizon/README.md rendered as written (internal/pdfreport/mdreport.py),
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
    r = Report("examples/ether_horizon", "Ether horizon: a void and a horizon in a nonlinear ether, on the photon electron", "void, horizon, sink, fountain, round void, Meyl comparison and g − 2 on the photon ring (Parts 1–11)", "register entry H49 (examples/ether_horizon/CLAUDE.md)")
    used = render_markdown(r, os.path.join(EX, "README.md"), base=EX)
    add_figures(r, os.path.join(EX, "viz"), skip=used)
    add_run_excerpt(r, os.path.join(EX, "run_output.txt"))
    r.h1("Reproducing the results")
    r.code("go test ./examples/ether_horizon\ngo run ./examples/ether_horizon\ngo run ./examples/ether_horizon -viz examples/ether_horizon/viz")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(EX, "ether_horizon_report.pdf"))
