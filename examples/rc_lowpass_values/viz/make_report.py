# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""rc_lowpass_values report: evolve two standard values for a 1 kHz RC cutoff (M5, CLAUDE.md §4).

python examples/rc_lowpass_values/viz/make_report.py examples/rc_lowpass_values/rc_lowpass_values_report.pdf

Numbers come from README.md, config.json, circuit.cir, doc/M5.md and
export/rc_lowpass_values_circuit_results.md. The example has no figures and
no run_output.txt. Layout from internal/pdfreport/kronreport.py.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401


def build(out):
    r = Report("examples/rc_lowpass_values", "RC low-pass: two standard values for a 1 kHz cutoff",
               "The first and simplest genetic-algorithm experiment — value genome, milestone M5, CLAUDE.md §4",
               "CLAUDE.md §4, M5", "population 40, seed 1")

    r.h1("Summary")
    r.p("A first-order RC low-pass filter is driven by a 1 V AC source. A genetic algorithm chooses the resistor from the "
        "E24 series and the capacitor from the E12 series so that the −3 dB frequency of V(out)/V(in) lands at 1 kHz. "
        "Everything is in a config file: the netlist is plain SPICE, and <i>cmd/evolve</i> reads the genes, the objective "
        "and the GA settings from <i>config.json</i>.")
    r.bullets([
        "<b>Result:</b> R1 = 160 Ω, C1 = 1 µF, giving f<sub>c</sub> = 1/(2πRC) = 994.7 Hz, an error of 0.53 %.",
        "<b>Found at once:</b> seed 1 stops after generation 0 (40 evaluations), because the random initial population "
        "already contains a design under the 0.6 % target.",
        "<b>The limit is the value grid, not the search:</b> 0.53 % is the closest any E24 × E12 pair gets to 1 kHz (the "
        "other such pair, 1.6 kΩ with 100 nF, has the same product).",
        "<b>The written design reproduces:</b> <i>evo.TestRCLowpassExample</i> re-parses <i>best.cir</i> and measures the same "
        "cutoff from it, within 10 Hz of the target.",
    ])

    r.h1("Background")
    r.p("CLAUDE.md §4 makes the value genome the starting encoding: fixed topology, genes that are component values on a "
        "log scale, optionally snapped to E-series. Milestone M5 is done when this RC low-pass example hits its target cutoff "
        "within 1 %. It is a standard engineering exercise; no hypothesis is involved. The other two config-driven M5 "
        "examples are <i>negative_resistance_gain</i> and <i>rectifier_filter</i>.")

    r.h1("Circuit, genes and objective")
    r.p("Template netlist (<i>circuit.cir</i>), with the starting values:")
    r.code("RC low-pass, values to be evolved\n"
           "V1 in 0 AC 1\n"
           "R1 in out 10k\n"
           "C1 out 0 10n\n"
           ".ac dec 100 10 100k\n"
           ".end")
    r.table([
        ["Gene", "Range", "Series"],
        ["R1", "100 Ω … 1 MΩ", "E24"],
        ["C1", "100 pF … 10 µF", "E12"],
    ], [0.3, 0.4, 0.3])
    r.table([
        ["Objective", "Target", "Error"],
        ["upper −3 dB edge of V(out)", "1000 Hz", "|f − 1000| / 1000"],
    ], [0.4, 0.2, 0.4])
    r.p("GA settings: population 40, up to 40 generations, seed 1, a 2 s deadline per evaluation, and <i>target_fitness</i> "
        "0.006, which stops the run as soon as a design is within 0.6 % of the target.")

    r.h1("Result")
    r.table([
        ["Part", "Value"],
        ["R1", "160 Ω"],
        ["C1", "1 µF"],
        ["Cutoff f<sub>c</sub> = 1/(2πRC)", "994.7 Hz (error 0.53 %)"],
        ["Generations / evaluations", "stopped after generation 0 / 40"],
    ], [0.5, 0.5])
    r.p("doc/M5.md records the same outcome: R1 = 1.6 kΩ / 160 Ω with C1 = 100 nF / 1 µF (E24 × E12), f<sub>c</sub> = 994.7 Hz "
        "(0.53 %, the closest any E24 × E12 pair gets).")

    r.h1("Validation")
    r.bullets([
        "<i>evo.TestRCLowpassExample</i> runs this config exactly as <i>cmd/evolve</i> does, checks the cutoff is within 10 Hz, "
        "re-parses the written <i>best.cir</i> and measures the same cutoff from it.",
        "Exported template (<i>export/rc_lowpass_values_circuit_results.md</i>, 3 elements, 2 nodes, AC at 1000 Hz): Kron's "
        "diakoptic solve tears it automatically into 2 pieces joined by 1 link (R1, carrying 5.32018×10<super>−5</super> A at "
        "phase 57.86°). The torn and full MNA solves differ by 5.6×10<super>−17</super> of the largest node voltage, and the power "
        "taken by the R, L, C branches is the same in both (Roth 1959 §5). V(out) of the template is 0.846733 V at phase −32.14°.",
        "Run outputs go to <i>runs/&lt;timestamp&gt;-rc_lowpass_values/</i>: a copy of the config, <i>stats.csv</i> (best, mean and "
        "worst fitness per generation), <i>best/genNNNN.cir</i>, <i>best.cir</i> and <i>checkpoint.json</i>. <i>best.cir</i> opens in "
        "ngspice, LTspice and QSpice as it is.",
    ])

    r.h1("Limits")
    r.bullets([
        "Only E24 resistors and E12 capacitors are allowed, so 0.53 % is the best achievable error; a finer series would be "
        "needed to get closer.",
        "The objective is the cutoff alone; component tolerance, source and load impedance are not modelled.",
    ])

    r.h1("Reproducing the results")
    r.code("go run ./cmd/evolve -config examples/rc_lowpass_values/config.json   # about 1 s\n"
           "go test ./evo -run TestRCLowpassExample                              # runs this config end to end\n"
           "KRONSPICE_EXPORT=1 go test ./internal/export                          # writes export/")
    r.small("See also doc/M5.md (the value-genome engine, run directories, penalties).")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "rc_lowpass_values_report.pdf"))
