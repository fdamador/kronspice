# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""rate_of_change report: where higher time derivatives enter electromagnetism (water's Debye response; H19).

python examples/rate_of_change/viz/make_report.py examples/rate_of_change/rate_of_change_report.pdf

Numbers come from run_output.txt (the program's report) and README.md (the
construction, Kramers–Kronig and torn-solve checks); export/water_cell_debye_results.md
for the exported cell. The example has no figures. Layout from
internal/pdfreport/kronreport.py.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401
from reportlab.platypus import KeepTogether  # noqa: E402


def build(out):
    r = Report("examples/rate_of_change", "Higher rates of change in electromagnetism",
               "Where do second and higher time derivatives enter: in the material (water's Debye response) or in "
               "the field equations (Bopp–Podolsky)? — register entry H19",
               "register entry H19", "deterministic")

    r.h1("Summary")
    r.p("The question comes from a premise recorded in <i>doc/vassilatos-references.md</i>: that a fast pulse acts on space "
        "the way a fast slap acts on water, so that above some rate the medium answers like a solid. In kronspice the "
        "answer lives in two places, kept apart. In the <b>material</b> (standard physics): a permittivity that depends on "
        "how fast the field changes, here water's Debye relaxation, solved in time, in frequency and as a kronspice "
        "circuit. In the <b>field equations</b> (hypothesis H19): Bopp–Podolsky electrodynamics, which adds a fourth-order "
        "term with one new length a, printed beside Maxwell (a = 0) on every line.")
    r.bullets([
        "<b>Water stiffens when hit fast.</b> A field ramped up in 1 µs sees ε = 78.6; in 10 ps, 36.1; in 1 ps, 9.6 (the "
        "high-frequency limit is 5.4). After settling it is always 78.6. The slap-on-water picture is real, and it is Debye relaxation.",
        "<b>Three forms of one model agree:</b> the kronspice AC circuit twin matches ε(ω) from 1 kHz to 100 GHz to "
        "7.9×10<super>−13</super>; the ramp integrator meets the analytic value to 10<super>−6</super>; Kramers–Kronig recovers 73.2 "
        "(= ε(0) − ε<sub>∞</sub>).",
        "<b>Bopp–Podolsky (hypothetical):</b> a point charge's potential is finite at its centre and falls below Coulomb's by "
        "the fraction e<super>−r/a</super>, so the effect is gone a few a away.",
        "<b>Electron self-energy</b> becomes finite, q²/(8πε<sub>0</sub>a), and equals m<sub>e</sub>c² at a = 1.409×10<super>−15</super> m.",
        "<b>A Coulomb-law test</b> good to 1 part in 10<super>6</super> at 1 cm only bounds a &lt; 0.000724 m; Bhabha scattering "
        "bounds it far more tightly (value to be entered from Bufalo et al. 2014). <b>Status of H19: untested.</b>",
    ])

    r.h1("What is being asked")
    r.h2("In the material (standard)")
    r.p("Maxwell's constitutive relation D = (1/4π)K E (<i>Treatise</i> Art. 619, CLAUDE.md §3c) with a permittivity K that "
        "depends on the rate of change of the field. A Debye term carries the first rate of change of the polarization P, "
        "τ·dP/dt + P = ε<sub>0</sub>Δε·E; a Lorentz term the second, d²P/dt² + γ·dP/dt + ω<sub>0</sub>²P = "
        "ε<sub>0</sub>Δε·ω<sub>0</sub>²E; and D = ε<sub>0</sub>ε<sub>∞</sub>E + ΣP<sub>k</sub> (Debye 1929; Lorentz 1909; Taflove "
        "&amp; Hagness ch. 9). In the frequency domain they are Δε/(1 + jωτ) and Δε·ω<sub>0</sub>²/(ω<sub>0</sub>² − ω² + jγω).")
    r.h2("In the field equations (hypothesis H19)")
    r.p("Bopp–Podolsky generalized electrodynamics (F. Bopp 1940; B. Podolsky, <i>Phys. Rev.</i> 62, 68, 1942) adds "
        "(a²/2)∂<sub>μ</sub>F<super>μσ</super>∂<super>ν</super>F<sub>νσ</sub> to Maxwell's Lagrangian, giving the fourth-order "
        "field equation (1 + a²□)□A = μ<sub>0</sub>J in Lorenz gauge, with one new length a. <b>Standard prediction beside "
        "every result:</b> a = 0 — Coulomb's potential, infinite point-charge self-energy, one wave mode.")

    r.h1("The model")
    r.bullets([
        "<b>Water at 25 °C</b> is one Debye term: ε(0) = 78.6, ε<sub>∞</sub> = 5.4 (Elton, arXiv:1704.01667 §II.B), τ = 8.3 ps "
        "(Barthel et al. 1990 / Buchner et al. 1999 via Elton Table I). Valid below about 100 GHz: the faster relaxation near "
        "1 ps and the infrared bands are left out.",
        "<b>Time domain:</b> an exponential integrator, exact for a field varying linearly between grid points "
        "(<i>em/dispersion.go</i>).",
        "<b>Circuit twin</b> (Kron/Branin network analogy, §3a.1): a plate cell of area A and gap d filled with the material is "
        "C<sub>∞</sub> = ε<sub>0</sub>ε<sub>∞</sub>A/d in parallel with one branch per term — a series R–C for Debye "
        "(C<sub>k</sub> = ε<sub>0</sub>ΔεA/d, R<sub>k</sub> = τ/C<sub>k</sub>) and a series R–L–C for Lorentz. For the 1 cm², 1 mm "
        "cell: C<sub>∞</sub> = 4.78 pF, C<sub>D</sub> = 64.8 pF, R<sub>D</sub> = 0.128 Ω. Its admittance is exactly "
        "jω(ε<sub>0</sub>A/d)·ε(ω), so the netlist solved by kronspice AC is a third form of the same model.",
        "<b>Podolsky</b> (<i>em/podolsky.go</i>): static point-charge potential and field, finite self-energy, the two free-wave "
        "modes and the Coulomb-test bound.",
    ])

    r.h1("Results")
    r.h2("1. Water at 25 °C (single Debye: ε(0) 78.6, ε∞ 5.4, τ 8.3 ps)")
    r.table([
        ["Frequency", "ε′", "ε″", "Loss tangent"],
        ["1 MHz", "78.60", "0.00", "0.000"],
        ["100 MHz", "78.60", "0.38", "0.005"],
        ["1 GHz", "78.40", "3.81", "0.049"],
        ["10 GHz", "62.95", "30.01", "0.477"],
        ["19.2 GHz (1/2πτ)", "41.95", "36.60", "0.872"],
        ["50 GHz", "14.79", "24.47", "1.655"],
        ["100 GHz", "8.00", "13.54", "1.693"],
    ], [0.3, 0.2, 0.2, 0.3])
    r.p("\"Slapping\" the water with an electric field ramped from 0 to 1 over a rise time T:")
    r.table([
        ["Ramp time T", "ε seen at the end of the ramp", "After it settles"],
        ["1 µs", "78.6", "78.6"],
        ["1 ns", "78.0", "78.6"],
        ["100 ps", "72.5", "78.6"],
        ["10 ps", "36.1", "78.6"],
        ["1 ps", "9.6", "78.6"],
    ], [0.3, 0.4, 0.3])
    r.p("Slow, the molecules follow (ε ≈ 78.6). Fast, they cannot turn in time and the water answers like a stiff solid "
        "(ε → 5.4). The analytic value at the end of the ramp is ε<sub>∞</sub> + Δε·(1 − (τ/T)(1 − e<super>−T/τ</super>)), which "
        "the integrator meets to 10<super>−6</super>.")

    n0 = len(r.s)
    r.h2("2. Bopp–Podolsky (hypothetical, H19)")
    r.p("Point charge: φ = q/(4πε<sub>0</sub>r)·(1 − e<super>−r/a</super>), finite at the centre (q/(4πε<sub>0</sub>a)). For an "
        "electron with a = 1 fm (illustration):")
    r.table([
        ["r / a", "Maxwell φ [V] (standard)", "Podolsky φ [V] (hypothetical)", "Fraction below Coulomb"],
        ["0.0", "infinite", "1.43996×10<super>6</super>", "1"],
        ["0.1", "1.43996×10<super>7</super>", "1.37031×10<super>6</super>", "0.905"],
        ["1.0", "1.43996×10<super>6</super>", "910231", "0.368"],
        ["3.0", "479988", "456091", "0.0498"],
        ["10.0", "143996", "143990", "4.54×10<super>−5</super>"],
        ["30.0", "47998.8", "47998.8", "9.36×10<super>−14</super>"],
    ], [0.15, 0.28, 0.3, 0.27])
    r.s[n0:] = [KeepTogether(r.s[n0:])]
    r.bullets([
        "<b>Self-energy.</b> Maxwell: infinite. Podolsky: q²/(8πε<sub>0</sub>a), equal to m<sub>e</sub>c² when a = 1.409×10<super>−15</super> m "
        "(half the classical electron radius) — Podolsky's motivation.",
        "<b>Free waves.</b> The ordinary Maxwell wave ω = ck, and a massive mode ω = c√(k² + 1/a²) starting at ω = c/a = "
        "3×10<super>23</super> rad/s (4.77×10<super>22</super> Hz) for a = 1 fm. The Podolsky mass is ħ/(ac).",
        "<b>Bounds.</b> A Coulomb-law test good to 1 part in 10<super>6</super> at 1 cm bounds only a &lt; r/ln(1/ε) = 0.000724 m. "
        "Measured Bhabha (e<super>+</super>e<super>−</super>) scattering bounds the Podolsky mass far more tightly (Bufalo, "
        "Pimentel &amp; Soto, <i>Phys. Rev. D</i> 90, 085012, 2014); that value is to be entered from the paper, not assumed.",
    ])

    r.h1("Validation")
    r.table([
        ["Comparison", "Result"],
        ["Circuit twin (kronspice AC) against ε(ω), 1 kHz to 100 GHz", "worst difference 7.9×10<super>−13</super>"],
        ["Time-domain ramp against the analytic end value", "10<super>−6</super>"],
        ["Time-domain response against ε(ω) (em/dispersion_test.go)", "twin test passes"],
        ["Causality, Kramers–Kronig: ε′(0) − ε<sub>∞</sub> = (2/π)∫ε″/ω dω", "73.2 recovered"],
        ["Exported water cell, torn automatically into two pieces joined by the Debye branch, against the full solve",
         "2.2×10<super>−16</super>; same power in both (Roth 1959)"],
    ], [0.62, 0.38])
    r.p("Tests: <i>TestCircuitTwin</i> (the water cell and a Lorentz test material solved by kronspice AC equal "
        "jω(ε<sub>0</sub>A/d)ε(ω) to 10<super>−9</super>, §3a rule 4); <i>TestSlap</i> (ramp responses equal the analytic end value to "
        "10<super>−6</super> and settle to ε(0)); <i>em/dispersion_test.go</i> (water's ε(ω), Kramers–Kronig, time–frequency twin; "
        "<i>TestPodolsky</i>: potential, field, self-energy and modes); <i>TestRuns</i> (the whole report runs).")

    r.h1("Limits of the model")
    r.bullets([
        "Water is a single Debye term, valid below about 100 GHz; the faster relaxation near 1 ps and the infrared bands are "
        "left out, so the 1 ps ramp value is at the edge of the model.",
        "Podolsky is computed for static point charges and free waves only; a = 1 fm is an illustration, not a fitted value.",
        "The measured Bhabha-scattering bound on a is not yet entered.",
    ])

    r.h1("Status")
    r.p("The material part is standard physics and is reproduced three ways. Register entry H19 (Bopp–Podolsky) is "
        "<b>untested</b>: the standard prediction a = 0 (Maxwell) is printed beside every result, and the effect falls below "
        "Coulomb by e<super>−r/a</super>, so only very short-range experiments can see it. The Vassilatos slap-on-water analogy "
        "is history only, not a physics source.")

    r.h1("Reproducing the results")
    r.code("go run ./examples/rate_of_change\n"
           "go test ./examples/rate_of_change ./em\n"
           "KRONSPICE_EXPORT=1 go test ./examples/rate_of_change -run TestExport   # the water cell's netlist and results sheet")
    r.small("Sources: J. C. Maxwell, <i>Treatise</i> (1873) Vol. II Art. 619; P. Debye, <i>Polar Molecules</i> (1929); H. A. Lorentz, "
            "<i>The Theory of Electrons</i> (1909); A. Taflove, S. C. Hagness, <i>Computational Electrodynamics</i>, 3rd ed. (2005), "
            "ch. 9; D. C. Elton, arXiv:1704.01667; J. Barthel et al., <i>Chem. Phys. Lett.</i> 165, 369 (1990); R. Buchner et al., "
            "<i>Chem. Phys. Lett.</i> 306, 57 (1999); F. Bopp, <i>Ann. Phys.</i> 430, 345 (1940); B. Podolsky, <i>Phys. Rev.</i> 62, 68 "
            "(1942); R. Bufalo, B. M. Pimentel, D. E. Soto, <i>Phys. Rev. D</i> 90, 085012 (2014); CODATA 2018; G. Vassilatos, "
            "<i>Secrets of Cold War Technology</i>, p. 27 (history only).")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "rate_of_change_report.pdf"))
