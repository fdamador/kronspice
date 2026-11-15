# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""H15 + H48 report: Eckardt's homopolar generator with the ECE and ECE2
spin-connection terms switched on, as a PDF.

python examples/eckardt_homopolar/viz/make_report_ece.py examples/eckardt_homopolar/eckardt_homopolar_ECE_report.pdf

Numbers come from `go run ./examples/eckardt_homopolar` (run 4) and the
README. Layout from internal/pdfreport/kronreport.py.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401
import kronreport  # noqa: E402

kronreport.H1.keepWithNext = 1
kronreport.H2.keepWithNext = 1


def build(out):
    r = Report("examples/eckardt_homopolar", "Eckardt's homopolar generator with ECE and ECE2",
               "The self-excited twin-disk machine rerun with the hypothetical spin-connection terms on — "
               "register entries H15 and H48, CLAUDE.md §3c, §3f", "register entries H15, H48")

    r.h1("Summary")
    r.p("The self-excited twin-disk homopolar machine of H. Eckardt (A.I.A.S./UPITEC, 2012) was rerun in kronspice "
        "with the hypothetical ECE and ECE2 spin-connection terms switched on. The plain run is in the companion "
        "report <i>eckardt_homopolar_report.pdf</i>. Each term was run in both forms, during spin-up under a constant "
        "drive and during the paper's flywheel discharge, with the full energy ledger.")
    r.bullets([
        "<b>The terms change the loop resistance, and so the resonance speed.</b> The loop resistance becomes "
        "R/(1 + s·κ) + s·c·ω<sub>0</sub>·L. ω<sub>res</sub> moves from 29.84 rad/s (Maxwell) to 37.30 (ECE torsion), "
        "14.92 (ECE2 torsion), 28.42 (ECE scalar) and 33.16 rad/s (ECE2 scalar).",
        "<b>The steady current does not change.</b> The back torque still sets it: α = K1·I + K2·I² gives "
        "I = 9.949 A in every case. Only the speed at which the machine holds that current moves, and the drive "
        "pays α·ω<sub>res</sub>.",
        "<b>ECE2 torsion is a negative resistance.</b> R heats at 29.7 W while the drive supplies 14.9 W. The other "
        "14.8 W is booked to the spin connection as <i>model-supplied</i>, never as inflow. In the flywheel run R "
        "takes twice the flywheel's energy, and the shaft is driven backwards.",
        "<b>Every run's books close</b> to 2.6×10<super>−5</super> J over about 100 J (driven) and "
        "4×10<super>−9</super> J (flywheel). Any gain is exactly what the declared terms supply.",
    ])

    r.h1("What is being tested")
    r.p("<b>H15.</b> ECE theory says homopolar resonances draw energy from spacetime. The plain run found that the "
        "paper's surplus is the work of the back torque its eq. (12) leaves out. With that torque included, the "
        "machine is an ordinary series dynamo, and its ledger closes.")
    r.p("<b>H48.</b> In Evans' ECE theory the spin connection enters the field–potential relations: "
        "E = −∇φ − ∂A/∂t − cω<sub>0</sub>A + ωφ (ECE; Eckardt, <i>ECE Unified Field Theory</i> Vol. 1, eq. 4.211). "
        "ECE2 has the same terms multiplied by −2 (eq. 6.189). Integrating them along a component gives two circuit "
        "terms. This mapping is a proposal of kronspice (<i>circuit/ece.go</i>), not taken from the book. With s = +1 "
        "for ECE and s = −2 for ECE2:")
    r.code("torsion (ω0·A):  v = L·di/dt + s·c·ω0·L·i          series resistance s·c·ω0·L in each coil\n"
           "scalar  (ω·φ):   body sees (1 + s·κ)·φ at each terminal   acts as R/(1 + s·κ) on a resistor")
    r.p("<b>Standard prediction.</b> ω and ω<sub>0</sub> are frame quantities with no field of their own, so the "
        "terms are zero and the plain run holds. A real energy gain needs a declared source.")

    r.h1("How the run is set up")
    r.p("The circuit engine attaches the terms to R, L, C and K elements only, not to the N machine element. In the "
        "plain run the coil's self-inductance (L = 1 mH) sits inside winding 3 of the N element. For this run it is "
        "moved out into a series L element. The loop current is the same everywhere, so the circuit is unchanged. "
        "The plain case reproduces the original run: drive work 84.75 J, heat in R 80.1 J, the same final speed and "
        "current, and books closing to 2.6×10<super>−5</super> J (<i>TestECEPlainTwin</i>).")
    r.p("No measurement of ω<sub>0</sub> or κ exists, so both are declared. They are sized as in "
        "<i>examples/ece_rlc</i>:")
    r.table([
        ["term", "declared value", "in this circuit", "ECE (s = +1)", "ECE2 (s = −2)"],
        ["torsion", "ω<sub>0</sub> = 2.50×10<super>−7</super> m<super>−1</super> on the coil",
         "c·ω<sub>0</sub>·L = 0.075 Ω (¼ of R)", "+0.075 Ω", "−0.150 Ω"],
        ["scalar", "κ = 0.05 on both terminals of R", "R/(1 + s·κ)", "0.286 Ω", "0.333 Ω"],
    ], [0.12, 0.3, 0.24, 0.17, 0.17])
    r.p("Parameters are otherwise those of the paper (p. 4): R 0.3 Ω, L 1 mH, N 800, l 5 cm, r 10 cm, μ<sub>r</sub> 100, "
        "B<sub>m</sub> 0.1 T, J 0.01 kg·m². That gives K1 = ½r²B<sub>m</sub> = 0.0005 V·s and "
        "K2 = ½r²μ<sub>0</sub>μ<sub>r</sub>N/l = 0.01005 V·s/A, and ω<sub>res</sub> = R<sub>loop</sub>/K2.")

    r.h1("Results")
    r.h2("Driven: α = 1 N·m from rest, 3 s")
    r.p("Energies in joules. <i>Spin-up</i> is the energy stored in the flywheel. <i>Model</i> is the energy the ECE "
        "terms supply (+) or absorb (−). <i>Books</i> is the ledger mismatch.")
    r.table([
        ["case", "ω<sub>res</sub> rad/s", "ω at 3 s", "I at 3 s", "drive", "spin-up", "heat in R", "model", "books"],
        ["plain (Maxwell)", "29.84", "30.33", "9.68 A", "84.75", "4.60", "80.10", "0", "2.6e−5"],
        ["ECE torsion", "37.30", "36.99", "9.08 A", "104.5", "6.84", "78.12", "−19.53", "2.6e−5"],
        ["ECE2 torsion", "14.92", "16.41", "8.11 A", "43.51", "1.35", "84.27", "+42.14", "1.9e−5"],
        ["ECE scalar", "28.42", "28.99", "10.04 A", "80.92", "4.20", "76.67", "+3.83", "2.6e−5"],
        ["ECE2 scalar", "33.16", "33.33", "9.06 A", "93.61", "5.56", "88.01", "−8.80", "2.6e−5"],
    ], [0.17, 0.11, 0.09, 0.1, 0.1, 0.1, 0.11, 0.1, 0.12])
    r.p("Each case settles toward its own ω<sub>res</sub> with the same current, 9.949 A. ECE2 torsion is still "
        "settling at 3 s (16.4 against 14.9 rad/s). The steady-state powers follow from the loop resistances:")
    r.table([
        ["case", "drive α·ω<sub>res</sub>", "heat in R", "ECE terms", "reading"],
        ["plain", "29.84 W", "29.69 W", "0", "the shaft pays for the heat"],
        ["ECE torsion", "37.30 W", "29.69 W", "−7.42 W", "an added loss in the coil"],
        ["ECE2 torsion", "14.92 W", "29.69 W", "+14.85 W", "a negative resistance: half the heat from the model"],
        ["ECE scalar", "28.42 W", "28.28 W (as R/(1 + κ))", "booked in ledger", "R appears smaller"],
        ["ECE2 scalar", "33.16 W", "32.99 W (as R/(1 − 2κ))", "booked in ledger", "R appears larger"],
    ], [0.15, 0.17, 0.18, 0.15, 0.35])

    r.h2("Flywheel: α = 0, disk spun to 40 rad/s, then connected, 0.5 s")
    r.table([
        ["case", "ω at 0.5 s", "peak I", "flywheel gave", "heat in R", "model"],
        ["plain (Maxwell)", "19.54 rad/s", "32.34 A", "6.092 J", "6.092 J", "0"],
        ["ECE torsion", "33.94 rad/s", "9.44 A", "2.241 J", "1.793 J", "−0.448 J"],
        ["ECE2 torsion", "−10.19 rad/s", "79.31 A", "7.481 J", "14.96 J", "+7.481 J"],
        ["ECE scalar", "16.72 rad/s", "36.80 A", "6.603 J", "6.603 J", "+0.330 J"],
        ["ECE2 scalar", "26.07 rad/s", "22.00 A", "4.603 J", "4.603 J", "−0.460 J"],
    ], [0.2, 0.16, 0.14, 0.17, 0.16, 0.17])
    r.p("With ECE torsion the loop resistance (0.375 Ω) holds the machine closer to self-excitation, so the "
        "discharge is weak and slow. With ECE2 torsion (0.150 Ω) the discharge is strong: the current peaks at "
        "79 A, R receives twice what the flywheel gave, and the dynamo drives the shaft through zero to "
        "−10 rad/s. The extra 7.48 J is booked to the spin connection. In the scalar cases the heat in R equals what "
        "the flywheel gave, and the model column is the scalar term's share at R's terminals.")

    r.h1("Ledger and validation")
    r.bullets([
        "<i>TestECEPlainTwin</i>: with no terms, the circuit with the coil as a series L element equals the original "
        "machine (final ω and I to 10<super>−3</super>), and it books no model-supplied energy.",
        "<i>TestECEBooks</i>: all eight ECE / ECE2 runs close their ledgers to 10<super>−5</super> of the energy "
        "throughput, and the ECE terms are the only model-supplied power.",
        "<i>TestECE2TorsionLowersResonance</i>: ECE2 ω<sub>res</sub> = 14.92 rad/s, and R heats beyond what the drive "
        "and the flywheel supply, with the difference booked as model-supplied.",
        "The power the terms deliver is labelled hypothetical (H48), booked per element as model-supplied, and never "
        "counted as inflow (CLAUDE.md §3c).",
    ])

    r.h1("Limits")
    r.bullets([
        "ω<sub>0</sub> and κ are declared, not measured. The results scale with them: the shift in ω<sub>res</sub> is "
        "linear in c·ω<sub>0</sub>·L/K2.",
        "The terms act only on the R and L elements. The N element (the disks and their speed voltages) carries "
        "none, because the circuit-engine mapping does not yet define them for a moving winding.",
        "The coil field is linear. The 20 T steady state is far above iron saturation, as in the plain run.",
        "UFT 456's conductivities (k0, kl0) are not used, and the two ECE-only resonances of UFT 107 are still "
        "not modelled, because no equations are supplied for them.",
    ])

    r.h1("Status")
    r.p("H15 stays <b>untested</b> and H48 <b>untested</b>. The ECE and ECE2 terms move the self-excitation speed, but "
        "they do not change the current the back torque sets. Every joule the machine gains beyond the shaft's "
        "work is supplied by the declared spin-connection terms, and the ledger shows that amount exactly. A bench "
        "test would have to measure the drive's shaft power against the heat in R at the resonance speed.")

    r.h1("Reproducing the results")
    r.code("go run ./examples/eckardt_homopolar          # run 4\n"
           "go test ./examples/eckardt_homopolar -run ECE\n"
           "python examples/eckardt_homopolar/viz/make_report_ece.py examples/eckardt_homopolar/eckardt_homopolar_ECE_report.pdf")
    r.small("Sources: H. Eckardt, “Simulation of resonances in the homopolar generator,” A.I.A.S./UPITEC, 13 February "
            "2012 (hypothesis source, H15); H. Eckardt, <i>ECE Unified Field Theory</i> Vol. 1 (2022), eqs. 4.211, "
            "6.189 (hypothesis source, H48); E. C. Bullard, Proc. Camb. Phil. Soc. 51, 744 (1955); G. Kron, "
            "Equivalent Circuits of Electric Machinery, Wiley 1951, Epilogue.")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "eckardt_homopolar_ECE_report.pdf"))
