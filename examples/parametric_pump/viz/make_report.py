# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Report: an LC tank whose capacitor is modulated at twice resonance, as a PDF.

    python examples/parametric_pump/viz/make_report.py examples/parametric_pump/parametric_pump_report.pdf

Numbers come from examples/parametric_pump/README.md (envelope growth and the
energy ledger of the run) and pump.cir. The example has no run_output.txt and
no figures. Layout from internal/pdfreport/kronreport.py.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401


def build(out):
    r = Report("examples/parametric_pump", "Parametric pump",
               "An LC tank whose capacitor is modulated at twice resonance: envelope growth against Mathieu's rate, "
               "with the pump's work on the ledger — CLAUDE.md §3a, M9c.1",
               "CLAUDE.md §3a, M9c.1")

    r.h1("Summary")
    r.p("A one-netlist example of kronspice's time-varying capacitor (the SCALE= extension). A 1 mH, 1 µF tank "
        "(f<sub>0</sub> = 5.033 kHz) starts with 1 V on the capacitor, and the capacitance is driven as "
        "C(t) = 1 µF·(1 − 0.1·sin 2ω<sub>0</sub>t): it falls while |v| peaks, so the tank's energy grows. It is the "
        "capacitive twin of the pumped inductance in examples/kron_machine, and the simplest form of the parametric "
        "circuits of H16 (examples/parametric_lcr). It gives the standard prediction beside which H16's \"energy from "
        "the time domain\" is read: the agent that varies C pays.")
    r.bullets([
        "<b>Envelope growth:</b> v(a) grows from 1.02 V at the first peak to 4.44 V at 1.992 ms, an exponential rate of "
        "740.8 /s.",
        "<b>Mathieu's rate</b> for pump depth m at 2ω<sub>0</sub>, less the series damping, is m·ω<sub>0</sub>/4 − "
        "R/(2L) = 790.6 − 50 = 740.6 /s (Kouřil &amp; Vrba §10.2.2, eq. 10.101).",
        "<b>The pump pays:</b> net pump work 10.36 µJ (11.70 µJ in through the mechanical port, 1.33 µJ returned) equals "
        "the stored-energy gain (0.50 → 10.25 µJ) plus the heat in R1 (0.61 µJ).",
        "<b>Ledger:</b> integration mismatch 1.2 nJ over the 2 ms; worst point residual 0.",
    ])

    r.h1("The netlist")
    r.code("C1 a 0 1u SCALE=SIN(1 0.1 10.066k 0 0 180) IC=1\n"
           "L1 a b 1m\n"
           "R1 b 0 0.1\n"
           ".tran 1u 2m UIC")
    r.p("SCALE=wave(...) multiplies the element's charge curve by s(t): the charge is q = s(t)·C·v, so the current is "
        "i = C·s·dv/dt + C·v·ds/dt. Here s(t) = 1 − 0.1·sin(2ω<sub>0</sub>t), at 10.066 kHz = 2f<sub>0</sub>, with a "
        "180° phase. The work of changing s, ṡ·∫Q, is the work of whatever agent squeezes the capacitor, and it is "
        "booked on the element's <b>mechanical port</b> (an EXT= channel can redirect it).")
    r.p("The 180° phase makes s fall as |v| peaks, so the pump puts energy in. With 0° the tank decays and the port "
        "takes energy out instead (analysis.TestParametricPumping).")

    r.h1("Result")
    r.table([
        ["quantity", "value"],
        ["envelope of v(a), first peak", "1.02 V"],
        ["envelope of v(a) at 1.992 ms", "4.44 V"],
        ["measured exponential rate", "740.8 /s"],
        ["Mathieu first-order rate m·ω<sub>0</sub>/4", "790.6 /s"],
        ["series damping R/(2L)", "50 /s"],
        ["predicted net rate", "740.6 /s"],
    ], [0.55, 0.45])
    r.p("analysis.TestParametricSelfExcitation holds the engine to this rate.")

    r.h1("Energy ledger over the 2 ms")
    r.table([
        ["entry", "energy"],
        ["in through the mechanical port", "11.70 µJ"],
        ["returned to the port", "1.33 µJ"],
        ["dissipated in R1", "0.61 µJ"],
        ["stored, start → end", "0.50 µJ → 10.25 µJ"],
        ["integration mismatch", "1.2 nJ"],
    ], [0.55, 0.45])
    r.p("The net work of the pump, 10.36 µJ, equals the stored-energy gain plus the heat: nothing comes from nowhere, "
        "and the ledger's worst point residual is 0. This is the standard prediction beside which H16 is read "
        "(CLAUDE.md §1): the energy the tank gains enters through the capacitor's mechanical port, from the agent that "
        "varies C.")

    r.h1("Validation")
    r.bullets([
        "analysis.TestParametricSelfExcitation: envelope growth rate against Mathieu's first-order rate.",
        "analysis.TestParametricPumping: 180° phase pumps energy in; 0° phase decays, with the port taking energy out.",
        "export/ holds the netlist, a KiCad schematic and the operating-point sheet (MNA only: Kron's torn solve on the "
        "sheet is linear-only, and C1 is a nonlinear capacitor).",
    ])

    r.h1("Limits")
    r.bullets([
        "The modulation s(t) is imposed; the agent that varies C is ideal and appears only through its work on the port.",
        "The growth rate is compared with Mathieu's first-order rate (pump depth 0.1).",
    ])

    r.h1("Reproducing the results")
    r.code("go run ./cmd/kronspice examples/parametric_pump/pump.cir\n"
           "go run ./cmd/kronspice -o pump.csv examples/parametric_pump/pump.cir     # waveforms and ledger columns\n"
           "go test ./analysis -run 'TestParametricPumping|TestParametricSelfExcitation'")
    r.small("Sources: F. Kouřil, K. Vrba, <i>Non-linear and Parametric Circuits</i>, Ellis Horwood 1988, §10.2.2; "
            "L. Mandelstam, N. Papalexi, <i>Z. Phys.</i> 73, 223, 1932; doc/M9c1.md (the pumped inductance against "
            "Mathieu's rate).")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "parametric_pump_report.pdf"))
