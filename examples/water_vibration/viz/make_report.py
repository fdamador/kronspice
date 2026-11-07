# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""PDF report for examples/water_vibration: water's bonds as a circuit.

Run from the repository root:

    python examples/water_vibration/viz/make_report.py examples/water_vibration/water_vibration_report.pdf

Numbers come from examples/water_vibration/run_output.txt (the program's
output) and README.md (construction, sources); the netlist excerpt is from
netlists/h2o.cir. The example has no figures. Layout from
internal/pdfreport/kronreport.py.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401


def build(out):
    r = Report("examples/water_vibration", "Water's bonds as a circuit",
               "Does kronspice's AC analysis of a mass–spring circuit give the infrared fundamentals of H<sub>2</sub>O, "
               "and predict its isotopologues? — doc/chemistry-plan.md",
               "doc/chemistry-plan.md")

    r.h1("Summary")
    r.p("The vibration spectrum of a water molecule is computed by kronspice's AC analysis: atoms are mass capacitors, "
        "bonds are spring inductors, and Wilson's B matrix is an E/F (controlled-source) transformer between the bond "
        "frame and the atom frame. Three force constants are fitted to H<sub>2</sub>O's three measured fundamentals; "
        "the same constants then predict D<sub>2</sub>O, HDO and H<sub>2</sub><super>18</super>O with no refit, since "
        "isotopes change masses, not bonds. This is the first molecule of the chemistry plan.")
    r.bullets([
        "<b>The circuit is the molecule's equations:</b> its AC peaks equal the GF eigenvalues (Wilson's method) to "
        "10<super>−6</super>, for every mode of every isotopologue.",
        "<b>Fitted force field</b> (to the H<sub>2</sub>O fundamentals): f<sub>r</sub> 7.6712, f<sub>rr</sub> −0.0942, "
        "f<sub>a</sub> 0.7030 mdyn/Å; reproduces 1595, 3657 and 3756 cm<super>−1</super> exactly.",
        "<b>D<sub>2</sub>O is predicted within 1.32 %</b> (−1.32 %, −0.89 %, −1.29 %), <b>HDO within 1.09 %</b> "
        "(−0.48 %, −1.09 %, +0.05 %).",
        "<b>The remaining error is anharmonicity:</b> a harmonic field fitted to H<sub>2</sub>O's fundamentals comes out "
        "soft and under-predicts the deuterium modes by about 1 %; the O–H stretch of HDO is within 0.1 %.",
        "<b>H<sub>2</sub><super>18</super>O</b> (1589.15, 3647.51, 3740.70 cm<super>−1</super>) is a prediction only: "
        "no verified measurement has been checked yet.",
    ])

    r.h1("The construction")
    r.p("The mobility (force–current) analogy of Firestone 1933, the one Kron and Branin use for network analogies "
        "(CLAUDE.md §3a.1):")
    r.table([
        ["Mechanics", "Circuit"],
        ["velocity of a coordinate", "node voltage"],
        ["force", "current"],
        ["mass m", "capacitor C = m to ground"],
        ["spring of stiffness k", "inductor L = 1/k"],
        ["damper", "resistor to ground"],
    ], [0.5, 0.5])
    r.bullets([
        "The six in-plane Cartesian coordinates of O, H1 and H2 are nodes with the atomic mass as a capacitor.",
        "The three bond coordinates (two O–H stretches and the H–O–H bend) are springs. Each spring's rate of change is a "
        "linear combination of the atomic velocities, s′ = Σ B<sub>kj</sub> v<sub>j</sub> (Wilson's B matrix), built from "
        "a chain of voltage-controlled voltage sources (E); the spring force flows in an inductor; and that force acts "
        "back on every atom through current-controlled current sources (F). The E/F pair is an ideal multiport "
        "transformer, Kron's connection tensor between the two frames.",
        "The force constants form the inductance matrix L = F<super>−1</super>, so the stretch–stretch interaction "
        "f<sub>rr</sub> is a mutual inductance (K element).",
        "Circuit units: capacitance in daltons, inductance in Å/mdyn, so the circuit's ω = 1/√(LC) maps to the "
        "molecule's through √(100 N/m / 1 u); the netlist header gives f<sub>c</sub> [Hz] × 8185.679725 = wavenumber "
        "[cm<super>−1</super>].",
        "Geometry r<sub>e</sub> = 0.958 Å, H–O–H angle 104.4776° (Hoy &amp; Bunker 1979). Three force constants "
        "(f<sub>r</sub>, f<sub>rr</sub>, f<sub>a</sub>; no stretch–bend term) are fitted by Newton's method to "
        "H<sub>2</sub>O's three measured fundamentals.",
    ])
    r.p("An excerpt of <i>netlists/h2o.cir</i>, the first O–H stretch (masses in u, a large damping resistor on each "
        "node, the B-matrix chain of E sources, the spring inductor and the F sources that return its force):")
    r.code("CM_ox ox 0 15.9949146196\n"
           "RD_ox ox 0 1000000\n"
           "E_r1_ox r1_s0 0 ox 0 -0.790569884316532\n"
           "E_r1_oy r1_s1 r1_s0 oy 0 -0.612371829864621\n"
           "E_r1_h1x r1_s2 r1_s1 h1x 0 0.790569884316532\n"
           "E_r1_h1y r1_s3 r1_s2 h1y 0 0.612371829864621\n"
           "VS_r1 r1_s3 r1_l 0\n"
           "L_r1 r1_l 0 0.130378191754178\n"
           "F_r1_ox ox 0 VS_r1 -0.790569884316532\n"
           "...\n"
           "K_rr L_r1 L_r2 0.0122798366009168\n"
           "I_fx 0 h1x AC 1\n"
           "I_fy 0 h1y AC 1\n"
           ".ac lin 4001 0.0610822823 0.549740541")

    r.h1("Results")
    r.h2("Force fields")
    r.table([
        ["Fitted to", "f<sub>r</sub> [mdyn/Å]", "f<sub>rr</sub> [mdyn/Å]", "f<sub>a</sub> [mdyn/Å]"],
        ["H<sub>2</sub>O fundamentals (used)", "7.6712", "−0.0942", "0.7030"],
        ["H<sub>2</sub>O harmonic frequencies", "8.4388", "−0.1191", "0.7514"],
    ], [0.4, 0.2, 0.2, 0.2])
    r.h2("Vibration frequencies")
    r.p("The circuit's AC peaks (coarse sweep, fine sweep, parabola through the top three points) against the GF "
        "eigenvalues and the measured gas-phase fundamentals (NIST CCCBDB: Shimanouchi; Huber &amp; Herzberg):")
    r.table([
        ["Molecule", "Mode", "Circuit AC peak [cm<super>−1</super>]", "GF eigenvalue", "Measured", "Error"],
        ["H<sub>2</sub>O", "bend", "1595.00", "1595.00", "1595.0", "−0.00 % (fitted)"],
        ["H<sub>2</sub>O", "stretch", "3657.00", "3657.00", "3657.0", "−0.00 % (fitted)"],
        ["H<sub>2</sub>O", "stretch", "3756.00", "3756.00", "3756.0", "−0.00 % (fitted)"],
        ["D<sub>2</sub>O", "bend", "1162.51", "1162.51", "1178.0", "−1.32 %"],
        ["D<sub>2</sub>O", "stretch", "2647.26", "2647.26", "2671.0", "−0.89 %"],
        ["D<sub>2</sub>O", "stretch", "2752.05", "2752.05", "2788.0", "−1.29 %"],
        ["HDO", "bend", "1395.34", "1395.34", "1402.0", "−0.48 %"],
        ["HDO", "stretch", "2697.36", "2697.36", "2727.0", "−1.09 %"],
        ["HDO", "stretch", "3708.98", "3708.98", "3707.0", "+0.05 %"],
        ["H<sub>2</sub><super>18</super>O", "bend", "1589.15", "1589.15", "—", "prediction only"],
        ["H<sub>2</sub><super>18</super>O", "stretch", "3647.51", "3647.51", "—", "prediction only"],
        ["H<sub>2</sub><super>18</super>O", "stretch", "3740.70", "3740.70", "—", "prediction only"],
    ], [0.14, 0.12, 0.24, 0.16, 0.14, 0.2])
    r.bullets([
        "H<sub>2</sub>O is fitted (0 % by construction); D<sub>2</sub>O, HDO and H<sub>2</sub><super>18</super>O are "
        "predictions with the same bonds.",
        "H<sub>2</sub><super>18</super>O is not in CCCBDB and its band origins have not been checked against a primary "
        "source, so it is shown as a prediction.",
        "The remaining error is anharmonicity. A measured fundamental sits below the harmonic frequency by a shift that "
        "shrinks with heavier atoms (ω·x scales as 1/μ against ω as 1/√μ, with μ the reduced mass and x the anharmonicity "
        "constant). A harmonic force field fitted to H<sub>2</sub>O's fundamentals absorbs hydrogen's larger shift, so "
        "it comes out soft and under-predicts the deuterium modes by about 1 %.",
    ])

    r.h1("Validation")
    r.bullets([
        "AC peaks equal the GF eigenvalues to 10<super>−6</super>: the circuit is the molecule's equations, not an "
        "approximation of them.",
        "Tests (<i>main_test.go</i>): the fit reproduces H<sub>2</sub>O; AC peaks equal eigenvalues; isotopes within "
        "1.5 %; the three zero modes (two translations, one in-plane rotation) are present; the force field is "
        "physically sensible (stretch stiffer than bend, a weak negative stretch–stretch interaction).",
        "Each isotopologue is a plain SPICE netlist (E, F, K elements) in <i>netlists/</i>, runnable in kronspice or any "
        "SPICE; <i>export/</i> holds the same netlists with KiCad schematics and results sheets.",
    ])

    r.h1("Limits")
    r.bullets([
        "Harmonic springs only: anharmonicity is the remaining ~1 % error. The next step in the plan is Morse "
        "(nonlinear) inductors.",
        "No stretch–bend coupling term; three constants for three H<sub>2</sub>O lines.",
        "In-plane motion only; gas-phase single molecule.",
    ])

    r.h1("Reproducing the results")
    r.code("go run ./examples/water_vibration [-netlists dir]\n"
           "go test ./examples/water_vibration")
    r.small("Sources: F. A. Firestone, J. Acoust. Soc. Am. 4, 249 (1933); E. B. Wilson, J. C. Decius, P. C. Cross, "
            "Molecular Vibrations (1955), Ch. 4; A. R. Hoy, P. R. Bunker, J. Mol. Spectrosc. 74, 1 (1979); NIST CCCBDB "
            "(T. Shimanouchi; K. P. Huber, G. Herzberg, 1979); NIST \"Atomic Weights and Isotopic Compositions\".")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "water_vibration_report.pdf"))
