# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""PDF report for examples/bearing_motor (Huber ball-bearing motor and the SPG, H13 and H14).

python examples/bearing_motor/viz/make_report.py examples/bearing_motor/bearing_motor_report.pdf

Numbers come from examples/bearing_motor/README.md, the register entries H13
and H14 in CLAUDE.md, and export/bearing_motor_spg_loop_results.md. The example
has no run_output.txt and no figures.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401


def build(out):
    r = Report("examples/bearing_motor", "The ball-bearing motor and the Space Power Generator",
               "A measurement registry for the Huber motor, and whether bearing torque can explain the SPG's "
               "figures — register entries H13 and H14", "register entries H13, H14")

    r.h1("Summary")
    r.p("A current through a rotating ball bearing drives the shaft with a torque that scales as I², aids the "
        "rotation and is zero at rest (the Huber effect). kronspice carries it as the P element (device/bearing.go). "
        "This example holds the public record of what has been measured about such motors, the fitting code that "
        "would set the element's two coefficients, and one question: could bearing torque, aiding the drive motor, "
        "account for the published figures of Tewari's Space Power Generator (SPG, register H13)? Marinov's claim of "
        "heat beyond the electrical input is register entry H14.")
    r.bullets([
        "<b>No coefficients are invented.</b> No measured (I, ω, T) or back-EMF table is open access, so the dataset list "
        "is empty and the P element refuses to run without fitted κ and ω<sub>0</sub>. The program lists the nine papers "
        "whose tables would fix them.",
        "<b>The SPG's torque gap is 24.8 N·m (7.01 kW):</b> the claimed output needs 41.5 N·m at the shaft, the drive "
        "motor supplied 16.7 N·m.",
        "<b>EM mode (the measured mode) leaves 1.37 V:</b> the bearings' back-EMF takes 2.04 V of the 3.41 V no-load EMF, "
        "so the rest of the loop gets 4.71 kW, not 11.72 kW. That is exactly the loaded-EMF bound of the SPG ledger.",
        "<b>Thermal mode is excluded:</b> at the measured efficiency η = 0.10 the bearings would need 70 kW of electrical "
        "input, 20.4 V across them, more than the 3.41 V the machine makes.",
        "<b>The engine run agrees:</b> whatever the contact resistance (1–100 µΩ), the load gets at most "
        "(E<sub>0</sub> − back-EMF)·I, with ledger residual 0.",
    ])

    r.h1("What is being tested")
    r.h2("H13 — Tewari's Space Power Generator")
    r.p("The SPG test (Tewari, <i>Universal Principles of Space and Matter</i>, pp. 159–160) reports 2700 r/min, 3438 A, "
        "a no-load EMF of 3.41 V and a drive motor rising from 2943 W to 8620 W at η = 0.83. <b>Standard prediction:</b> the generator's output is paid for by its shaft; "
        "the loaded EMF is the unmeasured quantity on which the claim rests (predicted ≤ 1.37 V, assumed 3.41 V). This "
        "example asks whether a known, measured mechanism, the bearing motor, could close the gap between the shaft "
        "work supplied and the electrical output claimed.")
    r.h2("H14 — Marinov's excess heat")
    r.p("Marinov (<i>Electronics &amp; Wireless World</i>, April 1989) reports, from a calorimeter, about 8 % more heat "
        "with rotation than without. This is a <b>claim</b> with no data shown, and it is not used in any calculation. "
        "<b>Standard prediction</b>, given by the P element's ledger: heat = V·I − ω·T.")

    r.h1("The registry")
    r.p("Every number is a record with its source and a status: <i>measured</i> (reported from an experiment), "
        "<i>estimate</i> (no instrument named) or <i>claim</i> (asserted, no data shown). Only measured values enter a calculation.")
    r.table([
        ["finding", "source", "status"],
        ["Steady currents 43.5–70.15 A; self-rotation at high current, small efficiency; a drive motor needed at low current",
         "Moyssides &amp; Hatzikonstantinou, IEEE Trans. Magn. 26, 1274 (1990), abstract", "measured"],
        ["Back-EMF that depends on shaft speed; no generator action", "same", "measured"],
        ["Self-rotation 70–140 A, aided 18–140 A (two discs)", "Moyssides 1992, via Lauterbach, Soong &amp; Abbott (2001)", "measured, secondary"],
        ["Torque proportional to I²", "Choo, Soong &amp; Abbott, Proc. SPIE 5649 (2005), per McDonald (2020)", "measured"],
        ["Mechanical output / electrical input = 0.10", "Marinov, Electronics &amp; Wireless World, April 1989 (calorimeter)", "measured, data unpublished"],
        ["Heat 8 % higher with rotation than without", "Marinov 1989", "claim (H14), not used"],
        ["1000 r/min (Milroy 1959); \"at least 1500 rev/min\" (Ogden 1989)", "Milroy, J. Appl. Mech. 34, 525 (1967); Ogden 1989", "measured / estimate"],
        ["Needs a hand spin to start; thermal effects start it, the Lorentz force cannot", "Harrison (electricstuff); McDonald (2020) fn. 6", "measured"],
        ["Bearing temperature 80 °C at 30 A", "Univ. Adelaide 2022", "measured"],
    ], [0.45, 0.37, 0.18])

    r.h1("The model")
    r.p("The P element's torque law is T = κ·I²·tanh(ω/ω<sub>0</sub>), where I is the current through the bearing, ω the "
        "shaft speed, κ a torque coefficient and ω<sub>0</sub> a speed scale. It runs in one of two modes: in EM mode the "
        "element has a back-EMF E with E·I = ω·T, so electrical power is converted to shaft work; in thermal mode there "
        "is no back-EMF: the element is a heat engine on the contact heating, with shaft work at most η times that heat. FitTorque fits κ and ω<sub>0</sub> to a measured (I, ω, T) table and FitBackEMF to a "
        "measured back-EMF table; in EM mode the two fits must agree.")
    r.p("Because no such table is open access, the dataset list is empty and the element refuses to run without fitted "
        "KAPPA and W0. The program lists the nine papers whose tables would fix the coefficients: Gruenberg 1978, "
        "Moyssides 1990 and 1997, Watson 1993 and 1999, Choo 2005, Shen 1999, Polivanov 1973 and Huber 1959. "
        "TestFitRecovers checks the fitting code on synthetic points (κ and ω<sub>0</sub> recovered to 10<super>−6</super>), "
        "not on a bearing.")

    r.h1("Results: could bearing torque explain the SPG?")
    r.p("All inputs are measured: the SPG test figures above and the bearing efficiency η = 0.10.")
    r.table([
        ["quantity", "value"],
        ["shaft torque the claimed E<sub>0</sub>·I needs", "41.5 N·m"],
        ["shaft torque the drive motor supplied", "16.7 N·m"],
        ["gap a bearing torque would have to fill", "24.8 N·m, 7.01 kW"],
        ["κ·ψ(ω) needed if all 3438 A cross the bearings", "≥ 2.1×10<super>−6</super> N·m/A²"],
    ], [0.6, 0.4])
    r.bullets([
        "<b>EM mode</b> (the measured mode: a back-EMF exists). The bearings then drop 2.04 V of the 3.41 V, and the rest "
        "of the loop gets 4.71 kW, not 11.72 kW. The 1.37 V left is exactly the loaded-EMF bound of the SPG ledger: a "
        "bearing motor moves shaft work around and adds none.",
        "<b>Thermal mode</b> at the measured η = 0.10: the bearings would need 70 kW of electrical input, 20.4 V across "
        "them, more than the 3.41 V the machine makes. Excluded.",
    ])
    r.h2("Engine run (SimulateSPGLoop)")
    r.p("One circuit holds the no-load EMF, the P element at the SPG's speed with the derived κψ, and the rest of the loop "
        "sized for 3438 A. The unknown bearing contact resistance R<sub>b</sub> is swept:")
    r.table([
        ["R<sub>b</sub>", "load", "bearing shaft work", "E<sub>0</sub>·I", "ledger residual"],
        ["1 µΩ", "4.70 kW", "7.01 kW", "11.72 kW", "0"],
        ["10 µΩ", "4.59 kW", "7.01 kW", "11.72 kW", "0"],
        ["100 µΩ", "3.53 kW", "7.01 kW", "11.72 kW", "0"],
    ], [0.18, 0.18, 0.22, 0.2, 0.22])
    r.p("The load never gets more than (E<sub>0</sub> − back-EMF)·I, whatever R<sub>b</sub>. The exported netlist and results sheet "
        "(export/bearing_motor_spg_loop_results.md) are solved by MNA only; Kron's torn solution is shown for linear circuits.")

    r.h1("Validation")
    r.bullets([
        "main_test.go: fit recovery on synthetic data, every registry entry sourced, the SPG bound (24.8 N·m, 1.37 V) and "
        "the swept engine run.",
        "Ledger residual 0 at every contact resistance in the sweep.",
    ])

    r.h1("Limits")
    r.bullets([
        "κ and ω<sub>0</sub> are not known for any real bearing; the analysis uses only the product κ·ψ(ω) the SPG would need, "
        "and asks whether energy could balance, not whether a bearing could produce that torque.",
        "The SPG figures are Tewari's published numbers; the loaded EMF was not measured.",
        "Marinov's 8 % excess heat is a claim with no data and is not used.",
    ])

    r.h1("Status and next tests")
    r.p("<b>H13</b> stays <i>untested</i>: the SPG test rests on an unmeasured loaded EMF (≤ 1.37 V predicted, 3.41 V "
        "assumed); bearing torque could close the 24.8 N·m gap only in EM mode, whose back-EMF leaves the same 1.37 V, "
        "and is excluded in thermal mode at η ≈ 0.10. <b>H14</b> stays <i>untested</i>: no data published; the "
        "standard prediction is the P element's ledger, heat = V·I − ω·T.")
    r.p("Decisive bench tests: insulate the SPG bearings (ceramic balls or insulated housings) and repeat the loaded run; "
        "measure the voltage across each bearing under load; put a torque transducer between the drive motor and the rotor. "
        "For H14, a calorimetric repeat of Marinov's run with the ledger above as the standard prediction.")

    r.h1("Reproducing the results")
    r.code("go run ./examples/bearing_motor\n"
           "go test ./examples/bearing_motor\n"
           "python examples/bearing_motor/viz/make_report.py examples/bearing_motor/bearing_motor_report.pdf")
    r.small("Sources: Huber, Elektrotechnik und Maschinenbau 76(8), 169 (1959); Milroy, J. Appl. Mech. 34, 525 (1967); "
            "Gruenberg, Am. J. Phys. 46, 1213 (1978); Moyssides &amp; Hatzikonstantinou, IEEE Trans. Magn. 26, 1274 (1990) "
            "and 33, 4566 (1997); Watson et al., IEE Proc. A 140 (1993), IEE Proc. Sci. Meas. Technol. 146, 83 (1999); "
            "Choo, Soong &amp; Abbott, Proc. SPIE 5649, 700 (2005); Marinov and Ogden, Electronics &amp; Wireless World, "
            "April 1989; McDonald, \"Ball-Bearing Motor\", Princeton (2020); Tewari, Universal Principles of Space and "
            "Matter (2017), §6.7.")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "bearing_motor_report.pdf"))
