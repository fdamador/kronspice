# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""H13 report: Tewari's space-vortex electron on the test bench, as a PDF.

python examples/space_vortex/viz/make_report.py examples/space_vortex/space_vortex_report.pdf

Numbers come from examples/space_vortex/README.md (the example has no
run_output.txt and no figures). Layout from internal/pdfreport/kronreport.py.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401


def build(out):
    r = Report("examples/space_vortex", "Tewari's space-vortex electron",
               "The electron as a vortex of a space fluid around a void, compared with measurement — register entry H13",
               "register entry H13", "Gaussian units; measured values CODATA 2018 and NIST ASD")

    r.h1("Summary")
    r.p("P. Tewari's 2018 paper (<i>Physics Essays</i> 31, 108) describes the electron as an irrotational vortex of an "
        "incompressible, massless “space fluid” whose speed reaches the speed of light c at the surface of a central "
        "spherical void of radius r<sub>e</sub>. From that void it states the electron's charge, mass, angular "
        "momentum and magnetic moment, and an electron orbit for hydrogen. This example puts those statements on the "
        "electron test bench of <i>examples/kelvin_ring</i> and compares each with measurement, with the Dirac "
        "electron and QED as the standard counterpart beside every row.")
    r.bullets([
        "<b>The g-factor is a real, parameter-free prediction:</b> g = 15/8 = 1.875 from eqs. 15–16 for any void "
        "radius, against the measured 2.00232, 6.4 % off a value known to about 10<super>−13</super>. With the "
        "standard ring moment ½qvρ it is 15/16 = 0.9375.",
        "<b>The proton charge</b> comes out as 12 e (eq. 30); atoms are neutral to ~10<super>−21</super> e.",
        "<b>The fine-structure constant</b> comes out as 15π²/16 = 9.253, 1268 times the measured "
        "7.297×10<super>−3</super>.",
        "<b>A density limit</b> of 3420 g/cm³ that white dwarfs (Sirius B, ~2×10<super>6</super> g/cm³) exceed.",
        "<b>The hydrogen orbit gives the wrong light:</b> speed 0.69 c and orbit frequency 3.405×10<super>18</super> Hz "
        "(X-rays), 1381 times the Lyman-α frequency.",
        "<b>Status: untested.</b> The paper gives no Zeeman prediction, the bench's discriminating test.",
    ])

    r.h1("The hypothesis and the standard prediction")
    r.p("<b>Hypothesis (H13).</b> The electron is a vortex of an incompressible space fluid around a void "
        "(Tewari, Physics Essays 31, 108, 2018). His reactionless generator, the other half of H13, is analysed in "
        "<i>examples/space_power_generator</i>.")
    r.p("<b>Standard prediction.</b> The Dirac electron with QED corrections, which gives the measured values in every "
        "row below: g = 2.00232, spin angular momentum ħ/2, and hydrogen with ground-state speed αc and ultraviolet "
        "Lyman lines.")

    r.h1("The construction")
    r.p("The paper's equations, in its own units (charge in cm³/s, mass in cm⁴/s):")
    r.code("q_e = pi^2 r_e^2 c          (eq. 4)\n"
           "m_e = (4 pi/3) r_e^3 c      (eq. 6)\n"
           "L   = (4/5) m_e c r_e       (eq. 15, angular momentum)\n"
           "mu  = (3/4) q_e c r_e       (eq. 16, magnetic moment)")
    r.bullets([
        "r<sub>e</sub> is fixed from the measured charge by eq. 4, as in the paper: r<sub>e</sub> = "
        "4.029×10<super>−11</super> cm.",
        "The paper defines its gram from the measured electron mass (eqs. 10–12), so charge and mass are both inputs.",
        "Eq. 16 as printed has the dimensions of an SI moment; the Gaussian moment of the same loop is that divided by "
        "c, and the g-factor is the same in either reading.",
        "Measured values: CODATA 2018 and the NIST Atomic Spectra Database (Gaussian units).",
    ])

    r.page_break()
    r.h1("Results")
    r.table([
        ["quantity", "model", "measured", "ratio", "note"],
        ["charge (eq. 4)", "4.803×10<super>−10</super> statC", "4.803×10<super>−10</super>", "1",
         "input: r<sub>e</sub> fitted to e"],
        ["void radius vs classical radius", "4.029×10<super>−11</super> cm", "2.818×10<super>−13</super>", "143",
         "no measured size; scattering bounds below about 10<super>−18</super> cm"],
        ["angular momentum (eq. 15) vs ħ/2", "8.802×10<super>−28</super> erg·s", "5.273×10<super>−28</super>", "1.669", ""],
        ["magnetic moment (eq. 16)", "1.451×10<super>−20</super> erg/G", "9.285×10<super>−21</super>", "1.563",
         "|μ<sub>e</sub>| = (g/2) μ<sub>B</sub>"],
        ["g-factor (eqs. 15–16)", "1.875 = 15/8", "2.00232", "0.9364", "independent of r<sub>e</sub>"],
        ["g with the ring moment ½qvρ", "0.9375 = 15/16", "2.00232", "0.4682",
         "p. 115 integrates dq·v·ρ; a ring current's moment is half that"],
        ["proton charge (eq. 30)", "12 e", "e", "12", "atoms neutral to ~10<super>−21</super> e"],
        ["fine-structure constant (UPSM eq. 5.2)", "15π²/16 = 9.253", "7.297×10<super>−3</super>", "1268",
         "book p. 130 blames experiment; α is measured to 1.5×10<super>−10</super>"],
        ["density limit (UPSM §2.8) vs Sirius B", "3420 g/cm³", "~2×10<super>6</super>", "0.0017", "white dwarfs exceed it"],
        ["hydrogen orbit radius (eq. 35) vs a<sub>0</sub>", "9.67×10<super>−10</super> cm", "5.292×10<super>−9</super>",
         "0.183", ""],
        ["hydrogen orbit speed vs αc", "2.069×10<super>10</super> cm/s (0.69 c)", "2.188×10<super>8</super>", "94.6", ""],
        ["orbit frequency vs Lyman-α", "3.405×10<super>18</super> Hz", "2.466×10<super>15</super>", "1381",
         "a classical orbit radiates at v/2πr"],
    ], [0.22, 0.2, 0.14, 0.1, 0.34], bold_rows=(5,))
    r.small("Ratio = model / measured. Measured values in the same units as the model column. a<sub>0</sub>: Bohr radius; α: fine-structure "
            "constant; μ<sub>B</sub>: Bohr magneton.")
    r.bullets([
        "<b>The g-factor is a real prediction.</b> g = 2m<sub>e</sub>cμ/(eL) from eqs. 15–16, and r<sub>e</sub>, q and "
        "m all cancel, so 15/8 (or 15/16 with the ring moment) holds for any void radius. It is 6.4 % from the "
        "measured value, which is known to about 10<super>−13</super>.",
        "<b>The hydrogen orbit gives the wrong light.</b> A charge orbiting at 0.69 c on 24 r<sub>e</sub> radiates "
        "X-rays; hydrogen's lines are ultraviolet and its ground state has v = αc.",
        "<b>Zeeman effect: not treated in the paper.</b> The bench's discriminating test needs the model's prediction "
        "for the splitting and polarisation of a line in a field B (measured: Δω = eB/2m<sub>e</sub>c per unit of "
        "m<sub>l</sub>, σ± circularly polarised).",
    ])

    r.h1("Validation")
    r.p("Before recording the comparisons, the test reproduces the paper's own numbers where it states them: "
        "r<sub>e</sub> ≈ 4×10<super>−11</super> cm, eq. 4 giving 4.7326×10<super>−10</super> at that radius, "
        "1 g = 8.8×10<super>6</super> cm⁴/s by eq. 11, and a void volume of 2.67×10<super>−31</super> cm³. "
        "<i>main_test.go</i> also checks g = 15/8 and α = 15π²/16.".replace("cm⁴/s", "cm<super>4</super>/s"))

    r.h1("Limits of the analysis")
    r.bullets([
        "The bench takes the paper's equations as printed; r<sub>e</sub> is fitted to the measured charge and the "
        "mass enters through the paper's own definition of the gram, so charge and mass are inputs, not predictions.",
    ])

    r.h1("Status")
    r.p("Register status: <b>untested</b>, pending a Zeeman prediction from the model. The register records g = 15/8 "
        "(15/16 with the standard ring moment) against 2.0023, proton charge 12e, α = 15π²/16 against 1/137, and a "
        "density limit white dwarfs exceed, as the values the model must also meet.")

    r.h1("Reproducing the results")
    r.code("go run ./examples/space_vortex\n"
           "go test ./examples/space_vortex")
    r.small("Sources: P. Tewari, Physics Essays 31(1), 108–129 (2018), doi:10.4006/0836-1398-31.1.108; P. Tewari, "
            "<i>Universal Principles of Space and Matter</i> (2017), eq. 5.2 and §2.8; CODATA 2018 (E. Tiesinga et al., "
            "Rev. Mod. Phys. 93, 025010, 2021); NIST Atomic Spectra Database. Notes and page checks in "
            "doc/tewari-references.md.")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "space_vortex_report.pdf"))
