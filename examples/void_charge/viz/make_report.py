# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""PDF report for examples/void_charge: the charged void (H51)

Run from the repository root:

    python examples/void_charge/viz/make_report.py examples/void_charge/void_charge_report.pdf

The body is examples/void_charge/README.md rendered as written (internal/pdfreport/mdreport.py),
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
    r = Report("examples/void_charge", "The charged void: electrical properties for the void, tested", "Tewari's charged void re-read; charge radius against hydrogen and LEP; exact charge; electron and positron voids meeting; the void's vibrations", "proposed register entry H51 (examples/void_charge/CLAUDE.md)")
    used = render_markdown(r, os.path.join(EX, "README.md"), base=EX)
    add_figures(r, os.path.join(EX, "viz"), skip=used)
    add_run_excerpt(r, os.path.join(EX, "run_output.txt"))
    r.h1("Reproducing the results")
    r.code("go test ./examples/void_charge\ngo run ./examples/void_charge -viz examples/void_charge/viz\npython3 examples/void_charge/viz/plot_charge.py")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(EX, "void_charge_report.pdf"))
