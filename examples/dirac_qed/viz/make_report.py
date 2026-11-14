# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""PDF report for examples/dirac_qed: the standard electron benchmark

Run from the repository root:

    python examples/dirac_qed/viz/make_report.py examples/dirac_qed/dirac_qed_report.pdf

The body is examples/dirac_qed/README.md rendered as written (internal/pdfreport/mdreport.py),
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
    r = Report("examples/dirac_qed", "Dirac + QED: the standard electron, built as the benchmark", "spin ½ and g = 2 from Dirac's equation, hydrogen, QED's g − 2, and the Penning-trap measurement as circuits", "standard benchmark (not a hypothesis)")
    used = render_markdown(r, os.path.join(EX, "README.md"), base=EX)
    add_figures(r, os.path.join(EX, "viz"), skip=used)
    add_run_excerpt(r, os.path.join(EX, "run_output.txt"))
    r.h1("Reproducing the results")
    r.code("go test ./examples/dirac_qed\ngo run ./examples/dirac_qed -viz examples/dirac_qed/viz\npython3 examples/dirac_qed/viz/plot_dirac.py")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(EX, "dirac_qed_report.pdf"))
