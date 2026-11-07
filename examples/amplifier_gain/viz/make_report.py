# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""PDF report for examples/amplifier_gain (2N3904 stage evolved for gain and bandwidth, M7).

python examples/amplifier_gain/viz/make_report.py examples/amplifier_gain/amplifier_gain_report.pdf

Numbers come from examples/amplifier_gain/README.md, config.json and circuit.cir,
and the seed table of doc/M7.md. The example has no run_output.txt and no figures.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401


def build(out):
    r = Report("examples/amplifier_gain", "Evolving a 2N3904 amplifier",
               "A common-emitter stage evolved for 30 dB, 100 Hz–200 kHz and its bias point "
               "— CLAUDE.md §4, milestone M7", "CLAUDE.md §4, M7",
               "seed 1, population 100, 300 generations")

    r.h1("Summary")
    r.p("This is the milestone M7 example and the first run of the value-genome genetic algorithm (§4) on a "
        "nonlinear circuit. A common-emitter stage around a 2N3904 transistor (Fairchild datasheet Gummel–Poon model) "
        "is driven from a 600 Ω source into a 10 kΩ load on 12 V. Evolution chooses eight standard component values "
        "(resistors from E24, capacitors from E12) to meet a mid-band gain, both −3 dB band edges and a bias point.")
    r.bullets([
        "<b>Seed 1 meets the spec:</b> 29.92 dB at 10 kHz, band edges at 100.0 Hz and 199.8 kHz, V<sub>CE</sub> 4.9 V, "
        "re-simulated from the written best.cir.",
        "<b>10 of 12 seeds pass</b> the acceptance test (30 ± 0.5 dB, edges within 5 %, V<sub>CE</sub> ≥ 3 V). Seeds 3 "
        "and 10 miss the gain by less than 0.1 dB beyond the 0.5 dB band (29.44 and 29.41 dB).",
        "<b>The E-series grid sets a floor:</b> the weighted error settles near 0.02, and no single-part change of the "
        "best design improves it.",
        "<b>Energy ledgers balance:</b> the acceptance test requires the AC and OP ledgers to close to 10<super>−9</super>.",
        "<b>Run time:</b> about 25 s for the full run; results are bit-identical for any worker count.",
    ])

    r.h1("The circuit")
    r.p("circuit.cir is a plain SPICE deck: VCC 12 V; VIN with RS 600 Ω coupled through CIN into the base; the bias "
        "divider R1/R2; RC to the collector; RE1 in series with the bypassed RE2, so RE1 alone sets the mid-band gain "
        "(≈ RC ∥ RL / RE1); CF from collector to base, the Miller capacitor that sets the upper edge; COUT 10 µF into "
        "RL 10 kΩ; and .ac dec 20 10 100meg.".replace("∥", "||"))
    r.p("The transistor is the datasheet's SPICE model (IS 6.734 fA, BF 416.4, VAF 74.03 V, CJC 3.638 pF, CJE 4.493 pF, "
        "TF 301 ps, RB 10 Ω). The datasheet print \"Ise=6.734\" is a known misprint for 6.734f and is corrected in the "
        "netlist. The engine solves the SPICE3 Gummel–Poon equations of device/bjt.go (forward and reverse Early effect, "
        "high injection, non-ideal base currents, depletion and diffusion charges, pnjlim limiting). The operating point "
        "is found by Newton with the M4 fallback chain, and the AC response is taken at a polished operating point (two "
        "extra Newton iterations so the small-signal transconductances are exact to round-off).")

    r.h1("Genes, objectives and GA settings")
    r.table([
        ["gene", "range", "series", "starting value (circuit.cir)"],
        ["R1", "10 kΩ – 1 MΩ", "E24", "47 kΩ"],
        ["R2", "1 kΩ – 100 kΩ", "E24", "10 kΩ"],
        ["RC", "470 Ω – 47 kΩ", "E24", "4.7 kΩ"],
        ["RE1", "10 Ω – 2.2 kΩ", "E24", "100 Ω"],
        ["RE2", "100 Ω – 10 kΩ", "E24", "1 kΩ"],
        ["CF", "1 pF – 1 nF", "E12", "10 pF"],
        ["CIN", "100 nF – 100 µF", "E12", "1 µF"],
        ["CE", "1 µF – 1 mF", "E12", "100 µF"],
    ], [0.15, 0.3, 0.15, 0.4])
    r.table([
        ["objective", "target", "weight"],
        ["gain at 10 kHz, V(o)/V(in)", "30 dB", "1"],
        ["upper −3 dB edge", "200 kHz", "1"],
        ["lower −3 dB edge", "100 Hz", "1"],
        ["collector voltage", "6 V", "0.2"],
        ["emitter voltage", "1.5 V (so V<sub>CE</sub> ≈ 4.5 V)", "0.2"],
    ], [0.45, 0.4, 0.15])
    r.p("GA settings: population 100, 300 generations, tournament 5, mutation rate 0.25 with σ = 0.1 decade, seed 1, "
        "2 s timeout per evaluation, stop at a weighted error of 0.02. Operators (as implemented in M5): tournament "
        "selection, BLX-0.5 crossover on log magnitudes (Eshelman &amp; Schaffer 1993), Gaussian mutation in decades, "
        "snapping to the E-series, elitism.")

    r.h1("Results")
    r.h2("Seed 1 after 300 generations")
    r.table([
        ["part", "evolved value"],
        ["R1 / R2", "30 kΩ / 6.8 kΩ"],
        ["RC", "5.1 kΩ"],
        ["RE1 / RE2", "68 Ω / 1.3 kΩ"],
        ["CF", "33 pF"],
        ["CIN / CE", "4.7 µF / 18 µF"],
    ], [0.4, 0.6])
    r.p("Re-simulated from the written best.cir: 29.92 dB at 10 kHz, edges at 100.0 Hz and 199.8 kHz, V<sub>CE</sub> "
        "4.9 V. best.cir opens in ngspice, LTspice and QSpice as it is.")
    r.h2("Seeds 1–12 (doc/M7.md)")
    r.p("The acceptance test (evo/amplifier_test.go) holds 30 ± 0.5 dB, edges within 5 % and V<sub>CE</sub> ≥ 3 V, the "
        "tolerance of the 5 % parts the design is made of.")
    r.table([
        ["seed", "gain", "lower edge", "upper edge", "V<sub>CE</sub>", "spec"],
        ["1", "29.92 dB", "100.0 Hz", "199.8 kHz", "4.91 V", "met"],
        ["2", "30.37 dB", "100.1 Hz", "201.3 kHz", "4.62 V", "met"],
        ["3", "29.44 dB", "100.0 Hz", "197.9 kHz", "5.04 V", "missed (gain)"],
        ["4", "30.18 dB", "100.5 Hz", "200.4 kHz", "4.40 V", "met"],
        ["5", "29.92 dB", "100.0 Hz", "199.9 kHz", "4.35 V", "met"],
        ["6", "29.96 dB", "99.9 Hz", "199.3 kHz", "3.97 V", "met"],
        ["7", "30.24 dB", "100.6 Hz", "199.6 kHz", "4.97 V", "met"],
        ["8", "29.90 dB", "100.0 Hz", "199.3 kHz", "4.25 V", "met"],
        ["9", "30.24 dB", "99.9 Hz", "200.0 kHz", "4.27 V", "met"],
        ["10", "29.41 dB", "101.1 Hz", "199.7 kHz", "4.58 V", "missed (gain)"],
        ["11", "29.80 dB", "100.2 Hz", "197.7 kHz", "4.45 V", "met"],
        ["12", "29.90 dB", "99.6 Hz", "197.6 kHz", "4.63 V", "met"],
    ], [0.1, 0.16, 0.16, 0.16, 0.14, 0.28], bold_rows=(1,))

    r.h2("How the search behaves")
    r.bullets([
        "<b>The E-series floor.</b> With E24/E12 values in about 4 % steps and five coupled targets, the weighted error "
        "settles near 0.02. No single-part change of the best design improves it, so that is the floor of the grid "
        "rather than of the search. A local search that moves several parts at once (the bias resistors together) would "
        "likely close the remaining gap; it is noted under M8 as future work.",
        "<b>Penalties, never a crash.</b> A candidate that fails to converge, gives a singular matrix, fails validation, "
        "breaks the ledger or times out gets a penalty fitness from the M5 classes (no-measurement &lt; no-convergence "
        "&lt; singular &lt; invalid &lt; ledger &lt; timeout &lt; other &lt; panic), so the run keeps going.",
        "<b>Determinism.</b> One PCG stream per generation derived from the seed, and results stored by individual "
        "index, so a run is bit-identical for any worker count and resumes exactly from its checkpoint (§1 goal 5).",
    ])

    r.h1("Validation")
    r.bullets([
        "TestAmplifierExample (evo/amplifier_test.go) is the M7 acceptance test: the spec above plus AC and OP energy "
        "ledgers balanced to 10<super>−9</super>.",
        "The export sheet for the starting circuit is MNA only: Kron's torn solution is shown for linear circuits, and the "
        "torn Newton solve (kron.JacobianTearing) is checked in kron's own tests.",
    ])

    r.h1("Limits")
    r.bullets([
        "The design is a single stage evaluated at 27 °C; the temperature coefficients of the model are unused.",
        "Results are for 5 % E-series parts; the spec band (±0.5 dB, 5 % edges) is set to that tolerance.",
        "Two of twelve seeds miss the gain target by under 0.1 dB beyond the band; a multi-part local search is future work.",
    ])

    r.h1("Reproducing the results")
    r.code("go run ./cmd/evolve -config examples/amplifier_gain/config.json   # about 25 s\n"
           "go test ./evo -run TestAmplifierExample                            # the M7 acceptance test\n"
           "python examples/amplifier_gain/viz/make_report.py examples/amplifier_gain/amplifier_gain_report.pdf")
    r.p("Outputs go to runs/&lt;timestamp&gt;-amplifier_gain/: a copy of the config, stats.csv (best, mean and worst "
        "fitness per generation), best.cir per generation and a final checkpoint.json that -resume continues from.")
    r.small("Sources: Fairchild 2N3904 datasheet (1997); Gummel &amp; Poon, Bell Syst. Tech. J. 49, 827 (1970); Nagel, "
            "ERL Memo M520 (1975), ch. 4; Massobrio &amp; Antognetti, Semiconductor Device Modeling with SPICE, 2nd ed., "
            "ch. 2; SPICE3f3 User's Manual §4.3; Eshelman &amp; Schaffer, FOGA 2 (1993).")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "amplifier_gain_report.pdf"))
