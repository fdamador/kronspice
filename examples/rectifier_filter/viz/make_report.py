# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""rectifier_filter report: evolve a half-wave rectifier's filter for mean output and ripple (M5, CLAUDE.md §4).

python examples/rectifier_filter/viz/make_report.py examples/rectifier_filter/rectifier_filter_report.pdf

Numbers come from README.md, config.json, circuit.cir, doc/M5.md and
export/rectifier_filter_circuit_results.md. The example has no figures and
no run_output.txt. Layout from internal/pdfreport/kronreport.py.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401


def build(out):
    r = Report("examples/rectifier_filter", "Half-wave rectifier: filter for mean output and ripple",
               "The first config-driven GA run that is nonlinear and transient — value genome, milestone M5 "
               "(with the M4 diode), CLAUDE.md §4",
               "CLAUDE.md §4, M5", "population 30, seed 11")

    r.h1("Summary")
    r.p("A 10 V peak, 1 kHz sine drives a diode into a filter capacitor and a load resistor. A genetic algorithm chooses "
        "the capacitor, the load resistor and the diode's series resistance so that, once the start-up transient has "
        "settled, the output averages 9 V with no more than 0.2 V of peak-to-peak ripple. It is the first config-driven "
        "example that needs the nonlinear Newton solve and a transient analysis for every candidate.")
    r.bullets([
        "<b>Result (seed 11):</b> C1 = 820 µF, R1 = 270 Ω, diode RS = 0.119 Ω, reached after generation 13 (394 evaluations).",
        "<b>Mean output 9.007 V</b> over 4–8 ms, 0.08 % from the 9 V target, re-simulated from the written <i>best.cir</i>.",
        "<b>Ripple 0.154 V</b> peak to peak, under the 0.2 V limit.",
        "<b>Why it works:</b> R1·C1 = 221 ms against the 1 ms period holds the ripple down; the mean is then set by the crest "
        "voltage, 10 V less the diode drop at the charging peak.",
        "<b>A .model parameter evolves like a part value:</b> the diode's RS sets the charging current's peak and the drop at the crest.",
    ])

    r.h1("Background")
    r.p("CLAUDE.md §4 describes the value genome (fixed topology, component values on a log scale, optional E-series "
        "snapping); M5 adds genes that target a <i>.model</i> parameter. This example exercises that together with the "
        "M4 diode and its Newton fallback chain (doc/M4.md). It is a standard engineering exercise; no hypothesis is "
        "involved. The other config-driven M5 examples are <i>rc_lowpass_values</i> and <i>negative_resistance_gain</i>.")

    r.h1("Circuit, genes and objectives")
    r.p("Template netlist (<i>circuit.cir</i>), with the starting values:")
    r.code("Half-wave rectifier: evolve the filter capacitor and load for ripple and mean output\n"
           "V1 in 0 SIN(0 10 1k)\n"
           "D1 in out dr\n"
           "R1 out 0 1k\n"
           "C1 out 0 1u\n"
           ".model dr D(IS=1e-14 N=1 RS=0.5)\n"
           ".tran 10u 8m\n"
           ".end")
    r.table([
        ["Gene", "Range", "Series"],
        ["C1", "100 nF … 1 mF", "E12"],
        ["R1", "100 Ω … 100 kΩ", "E24"],
        ["dr RS (a .model parameter)", "0.05 … 50 Ω", "continuous (log scale)"],
    ], [0.35, 0.35, 0.3])
    r.table([
        ["Objective", "Window", "Target", "Error"],
        ["mean of V(out)", "4 ms – 8 ms", "9.0 V", "|mean − 9| / 9"],
        ["peak-to-peak ripple of V(out)", "4 ms – 8 ms", "≤ 0.2 V", "max(0, ripple − 0.2) / 0.2"],
    ], [0.3, 0.18, 0.14, 0.38])
    r.p("Both objectives weigh 1; the ripple is a constraint, costing nothing once it is met. GA settings: population 30, "
        "up to 25 generations, seed 11, a 5 s deadline per evaluation and <i>target_fitness</i> 0.002.")

    r.h1("Result")
    r.table([
        ["Part / measure", "Value"],
        ["C1", "820 µF"],
        ["R1", "270 Ω"],
        ["dr RS", "0.119 Ω"],
        ["Mean V(out), 4–8 ms", "9.007 V (0.08 % from target)"],
        ["Ripple, 4–8 ms", "0.154 V (limit 0.2 V)"],
        ["Time constant R1·C1", "221 ms (period 1 ms)"],
        ["Generations / evaluations", "target reached after generation 13 / 394"],
    ], [0.5, 0.5])
    r.p("doc/M5.md summarizes the run as \"within 0.08 %, ripple met\".")

    r.h1("Validation")
    r.bullets([
        "Every transient is checked by the energy ledger as it runs: a candidate whose books do not balance is penalized, "
        "never reported. Per doc/M5.md, each point's residual is judged against the whole run's peak power after the run, "
        "and gross imbalances (&gt; 10<super>−3</super>) stop a run immediately.",
        "The figures above are from re-simulating the written <i>best.cir</i>, not from the GA's own record.",
        "Exported template (<i>export/rectifier_filter_circuit_results.md</i>, 4 elements, 2 nodes): operating-point sheet only. "
        "Kron's torn solution is shown for linear circuits; this one contains a diode, so the sheet reports the MNA solution "
        "(the torn Newton solve, <i>kron.JacobianTearing</i>, is checked in kron's tests).",
        "Run outputs go to <i>runs/&lt;timestamp&gt;-rectifier_filter/</i> (<i>stats.csv</i>, <i>best/genNNNN.cir</i>, <i>best.cir</i>, "
        "<i>checkpoint.json</i>); <i>best.cir</i> opens in ngspice, LTspice and QSpice as it is.",
    ])

    r.h1("Limits")
    r.bullets([
        "The diode is the SPICE3 level-1 model with IS = 10<super>−14</super> A and N = 1; only RS evolves.",
        "Mean and ripple are measured over the 4–8 ms window of an 8 ms run.",
        "The source is ideal; source impedance, diode recovery and capacitor ESR are not modelled.",
    ])

    r.h1("Reproducing the results")
    r.code("go run ./cmd/evolve -config examples/rectifier_filter/config.json   # about 20 s\n"
           "KRONSPICE_EXPORT=1 go test ./internal/export                         # writes export/")
    r.small("See also doc/M5.md (the value-genome engine) and doc/M4.md (the diode and the Newton fallback chain).")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "rectifier_filter_report.pdf"))
