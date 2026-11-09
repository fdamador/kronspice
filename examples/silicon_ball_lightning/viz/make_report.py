# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Report: silicon ball lightning, the chemical model and the welding arc.

python examples/silicon_ball_lightning/viz/make_report.py

Numbers from run_output.txt; figure from viz/render.py. Layout from
internal/pdfreport/kronreport.py.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report  # noqa: E402
from reportlab.platypus import KeepTogether  # noqa: E402


def build(out):
    r = Report("examples/silicon_ball_lightning", "Silicon ball lightning: the chemical model and the welding arc",
               "Can burning silicon nanoparticles make a ball that lasts seconds, and can a welding arc make them?",
               "standard counterpart to register entry H26", "about 4 seconds; deterministic")

    r.h1("Summary")
    r.p("Abrahamson and Dinniss (Nature 403, 519, 2000) propose that lightning striking soil reduces its silica and "
        "ejects a network of silicon nanoparticles whose slow oxidation in air gives a ball of light lasting seconds. "
        "Cen, Yuan and Xue (Phys. Rev. Lett. 112, 035001, 2014) recorded the first spectrum of a natural ball: silicon, "
        "iron and calcium lines from the soil throughout its 1.64 s. This example tests the chemical model's energy "
        "and timing, simulates a welding arc at the operating point of the laboratory silicon balls, and simulates the "
        "35 kV power line beside the Qinghai ball, all with kronspice and its energy ledger.")
    r.bullets([
        "<b>Energy:</b> burning silicon releases 32.4 kJ/g, 1.84 times the 17.6 kJ/g it costs to vaporize; an average "
        "ball (100 W for 10 s) needs only 31 mg. Iron returns 0.89 of its cost, manganese 1.24, calcium 3.08. Oxides "
        "(sand, rutile flux) store nothing.",
        "<b>Timing:</b> Deal–Grove kinetics burn 5–100 nm silicon particles in minutes to days at 1100–1500 K; "
        "10 nm particles burn in 1.6–10 s only at about 2200–2900 K (extrapolated), consistent with the 3140 K "
        "laboratory estimate.",
        "<b>Welding arc (kronspice):</b> 70 V open circuit behind 1.447 mH gives 120 A rms at 23 V, 2.5 kW, ledger "
        "2×10<super>−5</super>; 1–10 % of one second of arc vaporizing silicon stores 46–460 J, enough for a 100 W "
        "ball of 0.5–4.6 s.",
        "<b>Power line (kronspice):</b> the line's coupling makes any power in the ball pulse at 100 Hz (observed light "
        "99.4 Hz) but delivers at most 60 µW to a 15 cm ball against the 100 W it shines: it can flicker the ball, "
        "not power it.",
    ])

    r.h1("Results")
    r.figure(os.path.join(HERE, "silicon_ball_lightning.png"),
             "Figure 1. Top left: heat released on burning against the cost to vaporize, per gram. Top right: time for "
             "silicon particles to burn through (Deal–Grove, air), against the observed 1–10 s; grey, the measured range "
             "of the kinetics. Bottom left: the kronspice welder and arc. Bottom right: power the 35 kV line can drive "
             "into a ball 20 m away, against the 100 W it shines.", max_height_mm=170)

    r.h2("1. The chemical store")
    r.table([
        ["Element → oxide", "Burns (kJ/g)", "Vaporize (kJ/g)", "Return", "Mass for 1 kJ"],
        ["Si → SiO<sub>2</sub>", "32.43", "17.60", "1.84", "31 mg"],
        ["Fe → Fe<sub>2</sub>O<sub>3</sub>", "7.38", "8.32", "0.89", "136 mg"],
        ["Mn → MnO", "7.01", "5.65", "1.24", "143 mg"],
        ["Ca → CaO", "15.84", "5.14", "3.08", "63 mg"],
        ["Al → Al<sub>2</sub>O<sub>3</sub>", "31.05", "14.14", "2.20", "32 mg"],
        ["Ti → TiO<sub>2</sub>", "19.72", "11.79", "1.67", "51 mg"],
    ], [0.26, 0.18, 0.2, 0.14, 0.22])
    r.p("Soil holds silicon as SiO<sub>2</sub>: the strike must first reduce it, paying the same 32.4 kJ/g, so the "
        "ball releases energy the lightning stored. In a welding rod the steel core returns less than it costs; the "
        "silicon and manganese deoxidizers return more.")

    r.h2("2. How fast the particles burn")
    r.p("Deal and Grove's dry-oxidation constants for silicon (B = 772 µm²/h·exp(−1.23 eV/kT), "
        "B/A = 6.23×10<super>6</super> µm/h·exp(−2.00 eV/kT), J. Appl. Phys. 36, 3770, 1965), scaled by the oxygen partial "
        "pressure of air, with a particle burned through when its shell has consumed its radius, give the temperature "
        "at which a particle burns through in the observed lifetime:")
    r.table([
        ["Particle radius", "1.64 s", "3.5 s", "10 s"],
        ["10 nm", "2940 K", "2568 K", "2198 K"],
        ["25 nm", "4582 K", "3694 K", "2917 K"],
        ["50 nm", "above 5000 K", "above 5000 K", "4030 K"],
    ], [0.28, 0.24, 0.24, 0.24])
    r.p("At the 1100–1500 K where the kinetics were measured the same particles take minutes to days. Seconds need a "
        "hot ball (2000–3000 K) of fine particles, which the laboratory estimate (3140 K) and the white colour of the "
        "Qinghai ball support; above 1473 K these values are extrapolations.")

    r.h2("3. The welding arc in kronspice")
    r.p("An AC stick welder is a transformer whose leakage inductance sets the current; the arc is a current-"
        "controlled element whose voltage stays near 22.5 V. With 70 V open circuit, 1.447 mH gives 120.0 A rms at "
        "23.0 V rms and 2494 W, inside Paiva et al.'s 20–25 V, 100–140 A, 60 Hz; the ledger closes to "
        "2.1×10<super>−5</super>.")
    r.table([
        ["Metal", "Arc fraction vaporizing", "Vapor", "Chemical store", "100 W ball lasts"],
        ["Si", "1 % / 10 %", "1.4 / 14 mg", "46 / 460 J", "0.5 / 4.6 s"],
        ["Fe", "1 % / 10 %", "3.0 / 30 mg", "22 / 221 J", "0.2 / 2.2 s"],
        ["Mn", "1 % / 10 %", "4.4 / 44 mg", "31 / 310 J", "0.3 / 3.1 s"],
    ], [0.12, 0.26, 0.18, 0.22, 0.22])

    r.h2("4. The power line beside the Qinghai ball")
    r.p("A 35 kV line (20.2 kV to ground, 50 Hz) 20 m away couples through a mutual capacitance 4πε<sub>0</sub>a²/d; "
        "the ball has 4πε<sub>0</sub>a to ground. The power into its resistance pulses at 100.00 Hz in kronspice (light "
        "observed at 99.4 Hz). The most it can be is the matched-load value ¼ωC<sub>m</sub>²V²/(C<sub>m</sub> + C<sub>b</sub>): "
        "6.0×10<super>−5</super> W for a 15 cm ball, 0.25 W even for the full 5 m luminous size, against 100 W of light.")

    r.h1("Ledger and validation")
    r.bullets([
        "Welder: kronspice ledger mismatch 2.1×10<super>−5</super> of the source energy.",
        "Power line: the sweep maximum equals the matched-load formula to 2 % (TestLine).",
        "Deal–Grove B at 1473 K = 0.0477 µm²/h, the textbook value (TestDealGrove).",
        "Thermochemistry: NIST-JANAF / CRC values; vaporization cost uncertain by about ±10 % (mean heat capacities).",
    ])

    r.h1("Limits of the model")
    r.bullets([
        "Deal–Grove is planar-wafer kinetics measured at 973–1473 K; nanoparticles and higher temperatures are extrapolations.",
        "The fraction of arc energy that vaporizes the electrode is declared (1–10 %), not computed.",
        "The ball's radiation, convection, buoyancy and particle network structure are not modelled.",
        "The line-coupling capacitances are point-sphere estimates.",
    ])

    r.h1("Status")
    r.p("The chemical model passes the energy test easily and the timing test if the ball is about 2000–3000 K with "
        "particles of order 10 nm; the Qinghai flicker fits the nearby power line as a modulator, not a source. It is the "
        "standard counterpart to the field-held plasmoid of H26 (Force-Free_Harmonic_Plasmoids), which lasts "
        "microseconds to milliseconds. A welding arc can make the vapor: the bench test is the electrode fraction "
        "vaporized, the particle size and temperature from the spectrum, and silicon-bearing against plain steel electrodes.")

    n0 = len(r.s)
    r.h1("Reproducing the results")
    r.code("go run ./examples/silicon_ball_lightning -viz examples/silicon_ball_lightning/viz\n"
           "python examples/silicon_ball_lightning/viz/render.py\n"
           "python examples/silicon_ball_lightning/viz/make_report.py\n"
           "go test ./examples/silicon_ball_lightning")
    r.small("Sources: J. Abrahamson, J. Dinniss, <i>Nature</i> 403, 519 (2000); J. Cen, P. Yuan, S. Xue, <i>Phys. Rev. Lett.</i> "
            "112, 035001 (2014); G. S. Paiva et al., <i>Phys. Rev. Lett.</i> 98, 048501 (2007); K. D. Stephan, N. Massey, "
            "<i>J. Atmos. Sol.-Terr. Phys.</i> 70, 1589 (2008); V. Dikhtyar, E. Jerby, <i>Phys. Rev. Lett.</i> 96, 045002 (2006); "
            "B. E. Deal, A. S. Grove, <i>J. Appl. Phys.</i> 36, 3770 (1965); NIST-JANAF tables; CRC Handbook.")
    r.s[n0:] = [KeepTogether(r.s[n0:])]
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "silicon_ball_lightning_report.pdf"))
