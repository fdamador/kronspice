# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""PDF report for examples/double_hull (a craft's double hull as a ring waveguide).

python examples/double_hull/viz/make_report.py examples/double_hull/double_hull_report.pdf

Numbers come from examples/double_hull/run_output.txt, with the method,
tests and sources from examples/double_hull/README.md. The example has no
figures.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401


def build(out):
    r = Report("examples/double_hull", "A craft's double hull as a ring waveguide",
               "Gap voltage per watt of RF, DC bias, vacuum multipactor and the Kilpatrick limit, and the least RF "
               "power for a 20 kV gap — CLAUDE.md §3i", "CLAUDE.md §3i",
               "hull radius 5 m, height 10 m, gap 0.10 m, copper, mode 4, 1 kW")

    r.h1("Summary")
    r.p("The gap between a craft's inner and outer conducting hulls is a closed parallel-plate waveguide around the craft. "
        "Radiators in the gap drive it at one of its ring resonances, and a DC bias between the hulls is added to the RF. "
        "The questions are what gap voltage a given RF power reaches, whether the vacuum gap multipactors (a resonant "
        "secondary-electron discharge, §3i) on the way up, and what the least RF power is for a target peak voltage. This "
        "is the ring-tower experiment of examples/ring_towers moved inside a double wall.")
    r.bullets([
        "<b>1 kW gives about 1.7 kV:</b> at 37.685 MHz (mode 4, Q 9410) the travelling-wave gap amplitude is 1681–1683 V "
        "against the formula's 1680 V; a standing wave reaches 2379 V at its antinodes.",
        "<b>Almost every watt goes to the walls:</b> radiators 1000.0 W = wall loss 999.9 W + bias feed 0.141 W (travelling); "
        "39.7 mJ of RF energy sits in the gap.",
        "<b>Unbiased ring-up multipactors:</b> passing 760 V the worst growth per impact is 2.02 on as-received copper "
        "(1.13 conditioned).",
        "<b>A DC bias suppresses it:</b> the least bias on the scanned grid that clears the whole ring-up is 500 V.",
        "<b>Far below sparking:</b> Kilpatrick's limit at 37.7 MHz is 8.1 MV/m, against 0.049 MV/m for the 2.5 kV-biased "
        "standing wave.",
        "<b>Evolved 20 kV design:</b> a 0.470 m gap at mode 1 (9.09 MHz, Q 21729) with 451 W of RF and 12496 V of bias "
        "reaches 19963 V peak, multipactor-free (worst growth 0.06).",
    ])

    r.h1("What is being tested")
    r.p("Standard physics: a transmission-line resonator, the two-surface multipactor model, and Kilpatrick's empirical RF "
        "sparking limit. A DC bias across a coaxial gap is the measured cure for multipacting in the LEP2 352 MHz couplers "
        "(2.5 kV) and the KEKB 509 MHz couplers (2 kV); the example asks whether the same works for a hull-sized gap. The "
        "same multipactor model drives the electron pumps of examples/photon_electron.")

    r.h1("The model")
    r.bullets([
        "<b>The gap as a line.</b> The parallel-plate TEM line of the telegrapher's equations (Paul §4.2): per unit length "
        "around the hull L′ = μ<sub>0</sub>·d/w, C′ = ε<sub>0</sub>·w/d and Z<sub>0</sub> = η<sub>0</sub>·d/w, with d the gap "
        "and w the hull height. Wall loss R′ = 2·R<sub>s</sub>/w with R<sub>s</sub> = √(π·f·μ<sub>0</sub>·ρ) the skin-effect "
        "surface resistance of each copper wall, so α = R<sub>s</sub>/(η<sub>0</sub>·d). N sections (24 per wavelength) close "
        "the ring; here 96 sections and Z<sub>0</sub> = 3.77 Ω.",
        "<b>Solve.</b> The AC analysis gives the steady state at the drive frequency, the operating point gives the bias, and "
        "the gap voltage is their sum. Expected travelling-wave amplitude at mode n: √(Z<sub>0</sub>·P<sub>in</sub>·Q/(π·n)); "
        "a standing wave doubles it at the antinodes.",
        "<b>Bias.</b> The DC bias charges the outer hull through a high resistance, stores ½·C<sub>gap</sub>·V<sub>dc</sub>² "
        "once and draws no steady DC current. The gap then swings between V<sub>dc</sub> ± V<sub>rf</sub>.",
        "<b>Multipactor check</b> (em/multipactor.go; Hatch &amp; Williams 1954, Vaughan 1988). Electrons launched from either "
        "wall at every RF phase with 2 eV are tracked through the combined RF + DC field impact after impact, multiplying "
        "Furman's yield δ = δ<sub>max</sub>·s·x/(s − 1 + x<super>s</super>) at each one, with measured copper values (Baglin "
        "et al. 2002: as received δ<sub>max</sub> 2.03 at 262 eV, s 1.39; conditioned 1.13 at 318 eV, s 1.35). A phase whose "
        "mean yield per impact stays above 1 is a growing discharge. Ring-up passes through every amplitude below the final "
        "one, so each is checked.",
        "<b>Kilpatrick's limit</b> f[MHz] = 1.64·E²·exp(−8.5/E), E in MV/m, bounds the peak gap field.",
    ])

    r.h1("Results")
    r.h2("Example hull: radius 5 m, height 10 m, gap 0.10 m, copper, mode 4, 1 kW")
    r.p("All four runs resonate at f = 37.685 MHz with Q = 9410 and Z<sub>0</sub> = 3.77 Ω.")
    r.table([
        ["run", "RF gap amplitude", "formula (travelling)", "peak gap (RF + bias)", "wall loss", "bias feed", "bias energy"],
        ["travelling, 0 V", "1681–1683 V", "1680 V", "1683 V", "999.9 W", "0.141 W", "0 mJ"],
        ["travelling, 2500 V", "1681–1683 V", "1680 V", "4183 V", "999.9 W", "0.141 W", "87.794 mJ"],
        ["standing, 0 V", "0–2379 V", "1680 V", "2379 V", "999.7 W", "0.283 W", "0 mJ"],
        ["standing, 2500 V", "0–2379 V", "1680 V", "4879 V", "999.7 W", "0.283 W", "87.794 mJ"],
    ], [0.18, 0.15, 0.14, 0.15, 0.12, 0.12, 0.14])
    r.p("In every run the radiators deliver 1000.0 W, which the wall loss and the bias feed take, and 39.7 mJ of RF energy "
        "is held in the gap (Q·P/ω). The 2.5 kV bias stores its 87.794 mJ once.")

    r.h2("Multipactor in the vacuum gap (f·d = 3.77 MHz·m)")
    r.p("Worst-phase mean secondary yield per impact, as-received copper, with the impact energy; values above 1 (marked "
        "susceptible) grow.")
    r.table([
        ["V<sub>rf</sub>", "bias 0 V", "bias 500 V", "bias 2500 V", "bias 5000 V"],
        ["100 V", "0.37 at 14 eV", "0.07 at 25 eV", "0.06 at 126 eV", "0.06 at 248 eV"],
        ["300 V", "<b>1.10 at 53 eV, susceptible</b>", "0.07 at 22 eV", "0.06 at 115 eV", "0.06 at 238 eV"],
        ["1000 V", "<b>1.78 at 630 eV, susceptible</b>", "0.11 at 177 eV", "0.07 at 88 eV", "0.06 at 206 eV"],
        ["1700 V", "0.29 at 732 eV", "0.16 at 398 eV", "0.07 at 60 eV", "0.06 at 171 eV"],
        ["3400 V", "0.08 at 341 eV", "0.08 at 334 eV", "0.08 at 14 eV", "0.06 at 113 eV"],
        ["10000 V", "0.06 at 831 eV", "0.06 at 835 eV", "0.07 at 376 eV", "0.07 at 673 eV"],
    ], [0.12, 0.28, 0.2, 0.2, 0.2])
    r.h2("Ring-up")
    r.table([
        ["copper", "ring-up to", "worst growth, no bias", "least bias that clears it"],
        ["as received", "1680 V (travelling)", "2.02 at 760 V", "500 V"],
        ["as received", "2375 V (standing antinode)", "2.02 at 759 V", "500 V"],
        ["conditioned", "1680 V (travelling)", "1.13 at 760 V", "500 V"],
        ["conditioned", "2375 V (standing antinode)", "1.13 at 759 V", "500 V"],
    ], [0.2, 0.3, 0.25, 0.25])
    r.p("Kilpatrick's RF limit at 37.7 MHz is 8.1 MV/m; the 2.5 kV-biased standing wave peaks at 0.049 MV/m.")

    r.h2("Evolved design")
    r.p("The GA chooses gap, mode, RF power and bias for a 20000 V peak, at least 5000 V of it RF, with the least RF power, "
        "a multipactor-free ring-up and a field under Kilpatrick's limit.")
    r.table([
        ["quantity", "best of 332 evaluations"],
        ["gap, mode", "0.470 m, mode 1 (9.09 MHz, Q 21729)"],
        ["RF power, bias", "451 W, 12496 V"],
        ["peak gap voltage", "RF 7467 V + bias 12496 V = 19963 V"],
        ["multipactor", "ring-up worst growth 0.06 (at 7247 V)"],
        ["field", "0.042 MV/m against Kilpatrick's 5.3 MV/m"],
        ["ledger", "radiators 451.0 W = wall loss 448.3 W + bias feed 2.720 W"],
        ["bias supply", "delivers C·V² = 0.97 J once: 0.48 J stored, 0.48 J lost in the feed resistor"],
    ], [0.3, 0.7])

    r.h1("Ledger and validation")
    r.bullets([
        "Power balance in every run: radiator input = wall loss + bias feed; no steady DC current (no gap leakage is modelled).",
        "Tests (main_test.go) hold the resonance to n·c/ℓ times the ladder's lattice factor to 10<super>−9</super>, the "
        "amplitude to the formula within 1 %, Furman's yield, the multipactor zones and bias suppression, and Kilpatrick's inverse.",
        "Kron twin: the exported 193-node circuit, torn into 6 pieces, equals MNA to 5.6×10<super>−9</super> with the same power.",
    ])

    r.h1("Limits of the model")
    r.bullets([
        "Space charge and the spread of emission energies are left out of the multipactor check, so it marks where a "
        "discharge can start, not how large it grows.",
        "Gap leakage is not modelled, so the bias draws no steady current.",
        "The bias scan is on a grid (0, 500, 2500, 5000 V): 500 V is the least grid value that clears ring-up, not a threshold.",
        "The hull is an ideal parallel-plate ring line with copper walls; real hull penetrations and supports are not modelled.",
    ])

    r.h1("Reproducing the results")
    r.code("go run ./examples/double_hull\n"
           "go test ./examples/double_hull\n"
           "python examples/double_hull/viz/make_report.py examples/double_hull/double_hull_report.pdf")
    r.small("Sources: Paul, Analysis of Multiconductor Transmission Lines, 2nd ed., §4.2; Hatch &amp; Williams (1954); "
            "Vaughan, IEEE Trans. Electron Devices 35, 1172 (1988); Baglin et al., CERN LHC Project Report 472 (2002); "
            "Kindermann et al., EPAC 1996; Kijima et al., EPAC 2000; Kilpatrick, Rev. Sci. Instrum. 28, 824 (1957); "
            "Boyd, LANL LA-UR-82-1166 (1982).")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "double_hull_report.pdf"))
