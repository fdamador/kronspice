# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Report: a resistive divider asked for gain, where evolution finds negative resistance, as a PDF.

    python examples/negative_resistance_gain/viz/make_report.py \
        examples/negative_resistance_gain/negative_resistance_gain_report.pdf

Numbers come from examples/negative_resistance_gain/README.md (seed 5 result),
config.json, circuit.cir and export/negative_resistance_gain_circuit_results.md
(template circuit, MNA against Kron's torn solve). No run_output.txt and no
figures. Layout from internal/pdfreport/kronreport.py.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401


def build(out):
    r = Report("examples/negative_resistance_gain", "A divider with gain",
               "Evolution with signed values finds negative resistance; its power is booked as model-supplied — "
               "CLAUDE.md §3c, §4",
               "CLAUDE.md §3c, §4", "seed 5, population 30 × 30 generations")

    r.h1("Summary")
    r.p("A two-resistor divider is asked for a voltage gain of 2 (6.02 dB) at 1 kHz. No passive divider can do it: "
        "the only solution is a lower leg of −2 kΩ, so the search has to cross zero and settle on a negative value. "
        "The example, from milestone M5 (value genome), shows the <i>sign: either</i> gene and the energy ledger's "
        "model-supplied column.")
    r.bullets([
        "<b>Result:</b> after 30 generations (842 evaluations, seed 5) R2 = −2000.86 Ω, a gain of 6.017 dB "
        "(error 1.9×10<super>−4</super>). The exact answer is R2 = −2 kΩ.",
        "<b>The negative resistor delivers power:</b> v²/|R| = 2²/2000 = 1 mW; the run reports model-supplied power "
        "0.000999 W.",
        "<b>It is never booked as inflow from the exterior.</b> An ideal negative resistance stands for a tunnel diode, "
        "Gunn diode or op-amp negative-impedance converter whose bias supply the model leaves out.",
        "<b>The template circuit</b> (both legs 1 kΩ, gain ½) gives the same answer by MNA and by Kron's torn solve, "
        "with power invariant under tearing.",
    ])

    r.h1("Why gain needs a negative leg")
    r.p("The divider's transfer function is H = R2/(R1 + R2). With R1 = 1 kΩ and both resistors positive, H is "
        "always between 0 and 1. Setting H = 2 gives R2 = −2 R1 = −2000 Ω: H = −2000/(1000 − 2000) = 2. So the target "
        "is reachable only if the sign of R2 can change.")
    r.p("CLAUDE.md §3c separates two transformations of a resistor. <i>Duality</i> (exchanging voltage and current, as "
        "a gyrator does) maps R to 1/R, still positive. <i>Sign inversion</i> (time reversal, or a negative-impedance "
        "converter) maps R to −R: dissipation becomes generation. A negative resistance is accepted in netlists and "
        "genomes, as in ngspice, and delivers v²/|R|. Its power is <b>model-supplied</b> — listed per element, balancing "
        "the ledger, never counted as inflow — and its thermal port is not modeled.")

    r.h1("Circuit, gene and objective")
    r.code("V1 in  0   DC 1 AC 1\n"
           "R1 in  out 1k\n"
           "R2 out 0   1k          gene: 10 Ω … 100 kΩ in magnitude, sign either\n"
           ".ac lin 1 1k 1k")
    r.table([
        ["objective", "target", "error"],
        ["gain of V(out)/V(in) at 1 kHz", "6.0206 dB (×2)", "|gain − 6.0206| in dB / 20"],
    ], [0.4, 0.25, 0.35])
    r.p("With <i>sign: either</i> the sign is its own gene, flipped by mutation, so the log-scale magnitude and the sign "
        "evolve separately. GA settings (config.json): population 30, 30 generations, seed 5, mutation step σ = 0.3 "
        "decades, target fitness 10<super>−4</super> (not reached, so the run goes the full 30 generations).")

    r.h1("Result")
    r.table([
        ["quantity", "evolved (seed 5)", "exact"],
        ["R2", "−2000.86 Ω", "−2000 Ω"],
        ["gain at 1 kHz", "6.017 dB", "6.0206 dB"],
        ["fitness error", "1.9×10<super>−4</super>", "0"],
        ["evaluations", "842", ""],
        ["model-supplied power", "0.000999 W", "2²/2000 = 1 mW"],
    ], [0.35, 0.35, 0.3])
    r.p("The run prints:")
    r.code("model-supplied power 0.000999 W (idealized active elements: a real circuit needs a supply for it)")
    r.p("The evolved best.cir opens in ngspice, LTspice and QSpice as it is (ngspice accepts negative resistors too).")

    r.h1("Ledger and twin check")
    r.bullets([
        "The negative resistor's power is booked in the ledger's model-supplied column, per element. A fitness that "
        "measures harvested energy must count this column as an input; <i>evo.ModelSuppliedLimit</i> is the objective "
        "that bounds it.",
        "Template circuit (both legs 1 kΩ) at 1000 Hz, from the export results sheet: V(in) = 1 at phase 0.00°, "
        "V(out) = 0.5 at phase 0.00° by both MNA and Kron's torn solve.",
        "Kron's solve tears it automatically (ground as common bus, Happ Ch. XII §4) into 2 pieces of [2 1] unknowns "
        "joined by 1 link (R1, current 0.0005 at phase 0.00° A, voltage 0.5 at phase 0.00° V); the interconnection system has "
        "order 1 against 3 unknowns for the whole circuit.",
        "Torn against full solve: largest node-voltage difference 0; power V·I* taken by the R, L, C branches "
        "0.0005 at phase 0.00° (torn) against 0.0005 at phase 0.00° (untorn) — power is invariant under tearing (Roth 1959 §5).",
    ])

    r.h1("Limits")
    r.bullets([
        "The negative resistor is ideal: it has no bias supply, no voltage or current range and no frequency limit. A "
        "real tunnel diode, Gunn diode or negative-impedance converter takes the 1 mW from its supply.",
        "One frequency point (1 kHz); the circuit has no reactive elements, so the gain does not depend on frequency.",
    ])

    r.h1("Reproducing the results")
    r.code("go run ./cmd/evolve -config examples/negative_resistance_gain/config.json   # about 2 s\n"
           "go run ./cmd/kronspice runs/<run>/best.cir                                   # re-simulate the result\n"
           "KRONSPICE_EXPORT=1 go test ./internal/export                                 # export/ sheets")
    r.small("See also doc/M5.md, \"Negative resistance and nonlinear elements\", and the other config-driven M5 "
            "examples rc_lowpass_values and rectifier_filter.")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "negative_resistance_gain_report.pdf"))
