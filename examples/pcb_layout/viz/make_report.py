# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""pcb_layout report: a microstrip pair evolved from Maxwell's equations (M9i, CLAUDE.md §3c).

python examples/pcb_layout/viz/make_report.py examples/pcb_layout/pcb_layout_report.pdf

Numbers come from README.md (default GA run and the gap sweep), doc/M9i.md
(extraction validation, shorter GA run) and
export/pcb_layout_pair_w3_s2_results.md (torn against MNA solve). The example
has no figures. Layout from internal/pdfreport/kronreport.py.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401


def build(out):
    r = Report("examples/pcb_layout", "Evolving a microstrip pair from Maxwell's equations",
               "Can a GA find a 50 Ω, low-crosstalk trace pair when every candidate's L and C come from a field "
               "solve? — milestone M9i, CLAUDE.md §3c",
               "CLAUDE.md §3c, M9i", "population 20, 12 generations, seed 1")

    r.h1("Summary")
    r.p("Two parallel microstrip traces on FR-4 are evolved for a 50 Ω characteristic impedance and near-end "
        "crosstalk (NEXT) no worse than −30 dB at 100 MHz, on the least board width. Nothing about the traces "
        "is looked up: for every candidate the per-unit-length inductance and capacitance matrices are extracted "
        "from the cross-section by the quasi-static field solver of milestone M9i, the pair becomes a "
        "coupled-line ladder in the circuit, and an AC analysis measures the crosstalk onto the quiet trace.")
    r.bullets([
        "<b>Evolved design (218 evaluations):</b> width w = 2.766 mm, gap s = 3.020 mm, board width 8.55 mm.",
        "<b>It meets the specification:</b> Z<sub>0</sub> = 49.82 Ω and NEXT = −33.1 dB at 100 MHz.",
        "<b>Extracted line parameters:</b> ε<sub>eff</sub> = 3.257, L<sub>11</sub> = 301 nH/m, C<sub>11</sub> = 121.3 pF/m, "
        "coupling K = L<sub>12</sub>/L<sub>11</sub> = 0.0857.",
        "<b>Coupling falls with the gap</b> for a 3 mm trace: K = 0.287, 0.128, 0.030 and NEXT = −20.8, −29.1, −43.3 dB "
        "at s = 0.5, 2, 6 mm.",
        "<b>The field solver is validated</b> against coax (+0.07 %), Cohn's stripline (−0.006 % extrapolated), the "
        "two-wire line (0.1 %) and Hammerstad–Jensen microstrip (Z<sub>0</sub> −1.8 %, ε<sub>eff</sub> −0.4 %).",
    ])

    r.h1("Background")
    r.p("CLAUDE.md §3c asks that circuit elements be traceable to Maxwell's equations, and that <i>em/</i> supply "
        "lumped parameters from geometry so evolved circuits can include physical layout and parasitic coupling. "
        "Milestone M9i (<i>doc/M9i.md</i>) builds that extraction; this example is its demonstration: a genetic "
        "algorithm that evolves physical trace geometry, scored by a circuit simulation whose element values come "
        "from a field solve. No hypothesis is involved.")

    r.h1("The model")
    r.bullets([
        "<b>Extraction</b> (<i>em.Section</i>): strips of zero thickness on top of a 1.6 mm substrate, ε<sub>r</sub> = 4.4, "
        "air above, inside a box twelve substrate thicknesses wide. The equation div(ε grad φ) = 0 is solved with each "
        "conductor driven in turn, giving the 2×2 capacitance matrix C. The same solve with ε = 1 gives C<sub>0</sub>, and "
        "L = μ<sub>0</sub>ε<sub>0</sub>C<sub>0</sub><super>−1</super> (the TEM relation).",
        "<b>Circuit</b> (<i>em.CoupledLadder</i>): the 10 cm pair as 20 sections of the coupled telegrapher's equations, "
        "mutual inductance as a K element and mutual capacitance as a capacitor between the traces. The aggressor "
        "trace is driven through Z<sub>0</sub>, and both traces are terminated in Z<sub>0</sub> at both ends.",
        "<b>Measurement:</b> NEXT = |v(v_in)/v(a_in)| in dB from a one-point AC analysis at 100 MHz (v_in is the near "
        "end of the quiet trace, a_in the near end of the driven one).",
        "<b>Genome:</b> trace width w (0.3–6 mm) and gap s (0.2–10 mm).",
        "<b>Objectives (weighted sum):</b> |Z<sub>0</sub> − 50|/50; crosstalk above −30 dB; and 0.01 × board width "
        "(2w + s) in units of 10 mm. Each candidate's extraction is memoized across the objectives.",
    ])

    r.page_break()
    r.h1("Results")
    r.h2("Evolved pair")
    r.table([
        ["Run", "Evals", "w", "s", "Board", "Z<sub>0</sub>", "ε<sub>eff</sub>", "K = L<sub>12</sub>/L<sub>11</sub>", "NEXT at 100 MHz"],
        ["Default (pop 20, 12 gens, seed 1)", "218", "2.766 mm", "3.020 mm", "8.55 mm", "49.82 Ω", "3.257", "0.0857", "−33.1 dB"],
        ["Shorter (doc/M9i.md, -gens 8 -pop 16)", "114", "2.769 mm", "8.913 mm", "14.45 mm", "49.86 Ω", "3.281", "0.0130", "−50.8 dB"],
    ], [0.24, 0.09, 0.09, 0.09, 0.09, 0.09, 0.07, 0.11, 0.13])
    r.p("The default run also reports L<sub>11</sub> = 301 nH/m and C<sub>11</sub> = 121.3 pF/m. Both runs land on almost "
        "the same width, which sets Z<sub>0</sub>. The shorter run settled on a wider gap; the board-width weight is small, "
        "so either design satisfies the specification.")
    r.h2("Coupling against gap")
    r.p("For a 3 mm trace (checked by the example's test):")
    r.table([
        ["Gap s", "K = L<sub>12</sub>/L<sub>11</sub>", "NEXT at 100 MHz"],
        ["0.5 mm", "0.287", "−20.8 dB"],
        ["2 mm", "0.128", "−29.1 dB"],
        ["6 mm", "0.030", "−43.3 dB"],
    ], [0.3, 0.35, 0.35])

    r.h1("Validation")
    r.h2("Field extraction (doc/M9i.md, em/extract_test.go)")
    r.table([
        ["Case", "Reference", "Result"],
        ["Coax, ε<sub>r</sub> = 2.1", "2πε/ln(b/a)", "+0.07 % (C, L, Z<sub>0</sub>)"],
        ["Stripline, zero thickness", "Cohn 1954, exact K(k)/K(k′)",
         "−0.006 % after Richardson extrapolation (b/40, b/80)"],
        ["Two-wire line", "πε<sub>0</sub>/acosh(D/2a)", "0.1 % in a 120 mm box"],
        ["Microstrip, FR-4, w = 3 mm", "Hammerstad–Jensen 1980", "Z<sub>0</sub> −1.8 %, ε<sub>eff</sub> −0.4 %"],
    ], [0.3, 0.3, 0.4])
    r.p("The telegrapher ladder is checked in <i>em/tline_test.go</i>: a 2 m microstrip-like line gives Z<sub>0</sub> "
        "51.85 Ω and a transient delay of 9.67 ns (0.690 c), matching √(L/C) and √(LC).")
    r.h2("Kron's torn solution against MNA")
    r.p("The exported w = 3 mm, s = 2 mm pair (185 elements, 83 nodes; <i>export/pcb_layout_pair_w3_s2_results.md</i>) is "
        "solved at 100 MHz both by full MNA and by Kron's diakoptics, torn automatically into 2 pieces of 2 and 162 "
        "unknowns joined by 1 link (the source resistor RS, carrying 0.0104559 A at phase 0.01°). The largest node-voltage "
        "difference is 1.2×10<super>−15</super> of the largest node voltage, and the power taken by the R, L, C branches "
        "is 0.0104559 at phase −0.01° in both (power invariant under tearing, Roth 1959 §5). In that solve the quiet trace's "
        "near end is 0.0176376 V at phase 68.20° against 0.50019 V at phase −0.01° on the driven trace.")
    r.h2("Determinism")
    r.p("The GA test runs without a per-evaluation timeout and fails if any candidate timed out, so the evolved result "
        "is the same on every machine (CLAUDE.md §1 goal 5).")

    r.h1("Limits of the model")
    r.bullets([
        "Strips have zero thickness.",
        "Quasi-static (TEM) only: dispersion and radiation are left out (<i>em/limits.go</i> applies).",
        "The ladder's own λ/10 rule sets the section count (20 sections for 10 cm).",
        "Crosstalk is measured at a single frequency, 100 MHz.",
    ])

    r.h1("Reproducing the results")
    r.code("go run ./examples/pcb_layout [-gens 12] [-pop 20]\n"
           "go test ./examples/pcb_layout          # -short skips the GA run\n"
           "KRONSPICE_EXPORT=1 go test ./examples/pcb_layout -run TestExport")
    r.small("Sources: C. R. Paul, <i>Analysis of Multiconductor Transmission Lines</i>, 2nd ed. (2008); E. Hammerstad, "
            "Ø. Jensen, \"Accurate models for microstrip computer-aided design,\" <i>IEEE MTT-S</i> (1980); doc/M9i.md.")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "pcb_layout_report.pdf"))
