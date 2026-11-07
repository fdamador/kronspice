# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Kaluza's fifth dimension and the Kron machine, as a PDF report.

Run from the repository root:
    python examples/kaluza_kron/viz/make_report.py examples/kaluza_kron/kaluza_kron_report.pdf

Numbers come from examples/kaluza_kron/run_output.txt (the program's own
output) and README.md; the figure is viz/kaluza_kron.png. Layout from
internal/pdfreport/kronreport.py.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401


def build(out):
    r = Report("examples/kaluza_kron", "Kaluza's fifth dimension in a Kron machine",
               "Is electric charge motion along a fifth dimension, and is Kron's machine the same construction? "
               "CLAUDE.md §3a and §3f",
               "CLAUDE.md §3a, §3f",
               "c = G = 1 (Gaussian) for the geodesic; declared machine values")

    r.h1("Summary")
    r.p("Th. Kaluza (1921) added a fifth dimension to Einstein's space-time and found that one five-dimensional "
        "geometry holds both gravity and electromagnetism: a charge moving freely in 5-D is seen in 4-D as gravity "
        "plus the Lorentz force, and electric charge is momentum along the fifth dimension. Kaluza assumes nothing "
        "depends on the fifth coordinate (his \"cylinder condition\"); Hoffmann (1949) uses the same words for Kron's "
        "machines, whose winding charges never appear in the inductances. This example checks Kaluza's equations, "
        "solves the exact 5-D geodesic of a charge in a magnetic field, reproduces his difficulty with the electron "
        "and Klein's answer, and builds the same construction as a two-coil Kron machine solved three ways.")
    r.bullets([
        "<b>Kaluza's equations hold as printed.</b> His closed forms (eq. 3) agree with the three-index symbols of his "
        "metric to 6.8×10<super>−27</super> with his α and 1.8×10<super>−13</super> with α = 0.7. His coupling "
        "α = √(κ/2) = 3.0548×10<super>−14</super> (he printed 3.06×10<super>−14</super>).",
        "<b>Charge is fifth-dimension motion.</b> The exact geodesic is a Lorentz-force orbit with "
        "q/m = k·p<sub>w</sub>/√(1 + p<sub>w</sub>²), k = 2√G; orbit radius and frequency match the Lorentz force "
        "to all printed figures for p<sub>w</sub> from 0.1 to 10.",
        "<b>A ceiling on charge.</b> q/(m√G) can never exceed 2; the electron needs 2.041×10<super>21</super> and the "
        "proton 1.112×10<super>18</super>, so no real particle fits on a classical 5-D geodesic (Kaluza's own "
        "\"serious difficulty\", p. 971).",
        "<b>Klein's circle</b> has radius R = 2√G·ħ/(c·e) = 2 l<sub>P</sub>/√α = 3.784×10<super>−34</super> m; its "
        "first charged mode weighs 5.21×10<super>17</super> GeV, near the Planck mass, not an electron.",
        "<b>The Kron machine is the same construction.</b> Kron's geodesic in (q1, q2, θ), the rotor alone with the "
        "flux linkages held, and the kronspice machine element agree: 1.1×10<super>−10</super> rad/s and "
        "3.7×10<super>−4</super> rad/s largest speed differences; kronspice ledger mismatch 3×10<super>−11</super>.",
        "<b>No energy from the extra dimension.</b> The 5-D geodesic keeps its norm and the machine its energy "
        "(energy drift 6×10<super>−15</super>).",
    ])

    r.h1("Background")
    r.p("Kaluza's metric (his eq. 2) frames Einstein's 4-D metric with the electromagnetic potential. The modern form "
        "used for the geodesic is ĝ = (η + k²AA, kA; kA, 1), with k = 2√G, η the flat 4-D metric, A the vector "
        "potential and w the fifth coordinate. Because nothing depends on w, the momentum p<sub>w</sub> along it is "
        "conserved; the 4-D projection of the free-fall path is then a charged particle's orbit. Kaluza's eq. 8a, "
        "ρ0/μ0 = 2α·u<super>0</super>, is the small-momentum part of this relation.")
    r.p("In a Kron machine the coordinates are the winding charges and the rotor angle θ, and the inductances (the "
        "metric) depend on θ only. The charges are therefore cyclic coordinates, and their momenta, the flux "
        "linkages, are conserved, exactly as p<sub>w</sub> is. Eliminating them (Routh reduction) leaves the rotor "
        "alone in a potential set by the conserved fluxes: the machine's version of the force the fifth dimension "
        "leaves behind. This is not a hypothesis-register entry; it is a check of a historical theory and of the "
        "structure of Kron's machine equations.")

    r.h1("Results")
    r.h2("1. Kaluza's equations")
    r.table([
        ["quantity", "computed", "Kaluza / today"],
        ["α = √(κ/2), κ = 8πG/c²", "3.0548×10<super>−14</super>", "3.06×10<super>−14</super> (printed)"],
        ["eq. 3 closed forms vs symbols of eq. 1 metric (his α)", "largest difference 6.8×10<super>−27</super>", "holds as printed"],
        ["eq. 3, with α = 0.7", "largest difference 1.8×10<super>−13</super>", "holds as printed"],
        ["electron e/m", "1.77×10<super>7</super> is e/m in electromagnetic units per gram", "today 1.7588×10<super>7</super>"],
    ], [0.42, 0.33, 0.25])

    r.h2("2. A charge in a magnetic field as free motion in five dimensions")
    r.p("Exact 5-D geodesic, c = G = 1, Gaussian charge, uniform B = 0.05 along z, speed 0.3 c.")
    r.table([
        ["p<sub>w</sub>", "q/m", "orbit radius", "Lorentz radius", "ω", "Lorentz ω", "p<sub>w</sub> drift"],
        ["0.10", "0.1990", "31.60540", "31.60540", "0.009492", "0.009492", "6.7×10<super>−12</super>"],
        ["0.30", "0.5747", "10.94442", "10.94442", "0.027411", "0.027411", "2.7×10<super>−13</super>"],
        ["1.00", "1.4142", "4.44750", "4.44750", "0.067454", "0.067454", "1.7×10<super>−14</super>"],
        ["3.00", "1.8974", "3.31497", "3.31497", "0.090499", "0.090499", "9.6×10<super>−15</super>"],
        ["10.00", "1.9901", "3.16054", "3.16054", "0.094920", "0.094920", "6.8×10<super>−15</super>"],
    ], [0.1, 0.12, 0.15, 0.15, 0.14, 0.14, 0.2])

    r.h2("3. Kaluza's difficulty with the electron, and Klein's answer")
    r.table([
        ["quantity", "value"],
        ["ceiling on q/(m√G) for any 5-D geodesic", "k = 2"],
        ["electron e/(m√G)", "2.041×10<super>21</super>"],
        ["proton e/(m√G)", "1.112×10<super>18</super>"],
        ["Klein radius R = 2√G·ħ/(c·e) = 2 l<sub>P</sub>/√α", "3.784×10<super>−34</super> m (about 23 Planck lengths)"],
        ["first charged mode ħ/(Rc)", "5.21×10<super>17</super> GeV"],
    ], [0.55, 0.45])
    r.p("Klein (1926) makes the fifth dimension a tiny circle, so quantum momentum around it, p<sub>w</sub> = nħ/R, "
        "makes charge come in whole steps n·e. The files record a correction: an earlier written review gave "
        "1.3×10<super>−33</super> m, which paired the Heaviside coupling √(16πG) with the Gaussian charge, a factor "
        "√(4π) too large.")

    r.h2("4. The same construction in a Kron machine")
    r.p("A stator coil and a rotor coil with mutual inductance M·cos θ and lossless, shorted (superconducting) "
        "windings. Declared values: L1 = L2 = 1.0 H, M = 0.5 H, J = 0.01 kg·m²; start i1 = 1.0 A, i2 = 0.5 A, "
        "θ = 0.3 rad, at rest. Conserved flux linkages ψ1 = 1.2388 Wb, ψ2 = 0.9777 Wb (Kron geodesic: drift "
        "4×10<super>−15</super>; energy drift 6×10<super>−15</super>).")
    r.table([
        ["t (s)", "Kron's geodesic in (q1, q2, θ)", "rotor alone, fluxes held", "kronspice machine element"],
        ["0.1", "−0.707296 rad/s", "−0.707296 rad/s", "−0.706973 rad/s"],
        ["0.5", "−0.923503 rad/s", "−0.923503 rad/s", "−0.923217 rad/s"],
        ["1.0", "1.439010 rad/s", "1.439010 rad/s", "1.438928 rad/s"],
        ["1.5", "−1.340890 rad/s", "−1.340890 rad/s", "−1.341042 rad/s"],
        ["2.0", "0.666266 rad/s", "0.666266 rad/s", "0.666265 rad/s"],
    ], [0.1, 0.3, 0.3, 0.3])
    r.figure(os.path.join(HERE, "kaluza_kron.png"),
             "Figure 1. Left: 5-D free fall seen in 4-D, circular orbits in the magnetic field for p<sub>w</sub> = "
             "0.3, 1, 3 and 10. Centre: charge-to-mass ratio against momentum along the fifth dimension; the exact "
             "geodesic 2p/√(1 + p²) (blue) saturates at the ceiling 2 (dotted), while Kaluza's small-momentum eq. 8a, "
             "2p (dashed red), does not. Right: the Kron machine's rotor speed solved three ways, overlapping.")
    r.bullets([
        "The rotor-alone solution matches Kron's geodesic to 1.1×10<super>−10</super> rad/s, and kronspice to "
        "3.7×10<super>−4</super> rad/s (its windings are shorted through 10<super>−7</super> Ω).",
        "Because the charges never appear in the inductances (the cylinder condition), their momenta, the flux "
        "linkages, are conserved, like p<sub>w</sub>. Kron's machine is a working Kaluza model.",
    ])

    r.h1("Ledger and validation")
    r.bullets([
        "Momentum p<sub>w</sub> conserved to between 6.8×10<super>−15</super> and 6.7×10<super>−12</super> over the "
        "geodesic runs.",
        "Machine: flux linkages conserved to 4×10<super>−15</super>, energy to 6×10<super>−15</super>; kronspice "
        "ledger mismatch 3×10<super>−11</super>.",
        "Tests (<i>main_test.go</i>): eq. 3, geodesic = Lorentz force, Klein numbers, machine solved three ways.",
    ])

    r.h1("Limits and where it stands")
    r.bullets([
        "The geodesic uses the modern form of Kaluza's metric with the corner fixed (no scalar field).",
        "The Kaluza–Klein mathematics lives on as a standard tool in string theory and supergravity.",
        "No extra dimension has been observed: gravity follows Newton's inverse-square law down to 52 µm "
        "(Lee et al., Phys. Rev. Lett. 124, 101101, 2020), and a circle of Klein's size is far beyond any experiment.",
        "For kronspice the result is structural: Kaluza's geometry is exactly the structure of Kron's machine "
        "equations, and the machine shows it working. No energy comes from the extra dimension.",
    ])

    r.h1("Reproducing the results")
    r.code("go run ./examples/kaluza_kron -viz examples/kaluza_kron/viz   # about 1 second\n"
           "python3 examples/kaluza_kron/viz/render.py\n"
           "go test ./examples/kaluza_kron")
    r.small("Sources: Th. Kaluza, Sitzungsber. Preuss. Akad. Wiss. Berlin (Math. Phys.) 1921, 966–972 (V. T. Toth's "
            "translation); O. Klein, Z. Phys. 37, 895 (1926); B. Hoffmann, Rev. Mod. Phys. 21, 535 (1949); J. G. Lee "
            "et al., Phys. Rev. Lett. 124, 101101 (2020).")
    for f in r.s:
        if getattr(getattr(f, "style", None), "name", "") in ("h1", "h2"):
            f.keepWithNext = 1
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "kaluza_kron_report.pdf"))
