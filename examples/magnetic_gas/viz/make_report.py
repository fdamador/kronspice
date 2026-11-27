# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""PDF report for examples/magnetic_gas: argon, H₂, O₂ and HHO at 500 psi,
their breakdown and magnetism, a Marx generator with gas gaps and a plasma
varactor charge pump, with ECE2 (H48) beside the standard run.

Run from the repository root:

    python3 examples/magnetic_gas/viz/make_report.py examples/magnetic_gas/magnetic_gas_report.pdf

The body is examples/magnetic_gas/README.md rendered as written
(internal/pdfreport/mdreport.py), so the PDF and the README never disagree,
followed by any figures in viz/ it does not show and the head of
run_output.txt. Nothing is re-simulated. Layout from
internal/pdfreport/kronreport.py.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
EX = os.path.normpath(os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report  # noqa: E402
from mdreport import add_figures, add_run_excerpt, render_markdown  # noqa: E402


def build(out):
    r = Report("examples/magnetic_gas",
               "Magnetic gas: argon, hydrogen, oxygen and HHO at 500 psi",
               "Breakdown, magnetism after one electron is removed, compression, a four-stage Marx generator with gas "
               "spark gaps and a plasma-varactor charge pump; ECE2 beside the standard prediction",
               "CLAUDE.md §3g; register entry H48 for the ECE2 runs")
    used = render_markdown(r, os.path.join(EX, "README.md"), base=EX)
    add_figures(r, os.path.join(EX, "viz"), skip=used)
    add_run_excerpt(r, os.path.join(EX, "run_output.txt"))
    r.h1("Reproducing the results")
    r.code("go test ./examples/magnetic_gas\n"
           "go run ./examples/magnetic_gas -viz examples/magnetic_gas/viz/data > examples/magnetic_gas/run_output.txt\n"
           "python3 examples/magnetic_gas/viz/plot.py\n"
           "python3 examples/magnetic_gas/viz/schematic.py\n"
           "python3 examples/magnetic_gas/viz/make_report.py examples/magnetic_gas/magnetic_gas_report.pdf")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(EX, "magnetic_gas_report.pdf"))
