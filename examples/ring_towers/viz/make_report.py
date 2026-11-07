# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Report: ring towers, standing waves on a closed ring, as a PDF.

python examples/ring_towers/viz/make_report.py examples/ring_towers/ring_towers_report.pdf

Numbers come from examples/ring_towers/run_output.txt (the program's output)
and examples/ring_towers/README.md and viz/README.md. The example has no
static figures (its viz/ pages are interactive HTML), so the report has
tables only. Layout from internal/pdfreport/kronreport.py.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401


def build(out):
    r = Report("examples/ring_towers", "Ring towers",
               "Standing, traveling and cancelling waves driven by towers on a closed ring — CLAUDE.md §3e",
               "CLAUDE.md §3e", "ring scaled to MHz; each state run to steady state by transient analysis")

    r.h1("Summary")
    r.p("A thought experiment as a circuit: towers at points around the Earth–ionosphere ring drive the same "
        "frequency into it, each tower sending one wave east and one west, with phases chosen so that the waves "
        "add, cancel or circulate. The ring is scaled to MHz so each lap takes microseconds; the patterns depend "
        "only on the number of wavelengths around the ring and on the phases. Every joule in the ring is supplied "
        "by the towers, and the energy ledger checks this (CLAUDE.md §3c).")
    r.bullets([
        "<b>Two towers</b> make a traveling wave (envelope min/max 0.925, 0.155 W flowing east), a standing wave "
        "(min/max 0.018, almost no net flow) or near-total cancellation (peak 0.167 V against 3.3 V), by phase alone.",
        "<b>Synchronized pairs add.</b> Two standing-wave pairs phased together double the amplitude (peak 6.656 V "
        "against 3.328 V) and quadruple the stored energy (0.951 µJ against 0.238 µJ) and the tower power.",
        "<b>Eight synchronized towers</b> give 11.95–13.22 V around the ring from 0.5006 W and store 3.8 µJ, sixteen "
        "times one pair; the 45° set in antiphase cancels everything to 0.32 V.",
        "<b>DC bias:</b> on a continuous ring two pairs at 100 V and 50 V settle at 75.0 V with 50 µA between the "
        "supplies; with 10 µF blocking capacitors each sector holds its own level (100.0 / 50.0 V) with no DC "
        "current, and the RF rides unchanged on top.",
        "<b>Ledger:</b> tower power equals ring loss in every state; the mismatch is at most 1.1×10<super>−9</super> "
        "of the energy put in.",
    ])

    r.h1("The construction")
    r.p("The ring is a lumped transmission line (telegrapher's equations; Steinmetz Lecture VI, Paul §2): 96 sections "
        "of series inductance L = 1 µH and resistance R with shunt capacitance C = 1 nF, the last section closing back "
        "on the first.")
    r.table([
        ["quantity", "value"],
        ["sections", "96, closed into a ring"],
        ["per section", "L = 1 µH, C = 1 nF, R = 0.132 Ω"],
        ["characteristic impedance Z0", "31.6 Ω"],
        ["loss per lap", "0.20 Np"],
        ["drive frequency", "1.314 MHz, the fourth resonance (four wavelengths around)"],
        ["tower", "two 10 mA current sources 6 sections (λ/4) apart"],
        ["steady state", "8 amplitude-decay times of transient analysis"],
    ], [0.35, 0.65])
    r.p("Because the ring is four wavelengths around, waves returning after each lap add to the new ones, and the "
        "build-up is set by the loss per lap. Each tower is a two-element end-fire array: with the eastern element "
        "lagging 90° it launches east only, leading 90° west only. Tower B sits at the antipode of A.")
    r.bullets([
        "<b>Traveling:</b> both towers launch east, B phased to the arrival of A's wave, so one wave circulates with "
        "the same amplitude everywhere.",
        "<b>Standing:</b> A launches east and B west, so counter-running waves form fixed nodes and antinodes with no "
        "net flow.",
        "<b>Cancel:</b> both launch east with B 180° from A's arrival, so the waves subtract.",
        "<b>Four towers:</b> pair C (90°) and D (270°) is added. With four wavelengths around, every quarter point is an "
        "antinode of A–B's cos(kx) pattern, all swinging in step, so phasing C–D to the same reference puts its "
        "pattern on top of A–B's. Contrasts: C–D in antiphase, and C–D shifted by an eighth of a wavelength (45°).",
        "<b>Eight towers</b> add the 45°-rotated pairs E–F and G–H.",
        "<b>DC bias:</b> a supply per pair through 1 MΩ feeds. A sectored ring puts 10 µF blocking capacitors at 45°, "
        "135°, 225° and 315° (reactance 0.012 Ω against ωL = 8.3 Ω per section).",
    ])

    r.h1("Results: two towers")
    r.table([
        ["state", "launcher phases", "peak [V]", "min / max", "east flow", "towers deliver = ring loss",
         "stored", "ledger mismatch"],
        ["traveling", "A1 +0°, A2 −90°, B1 −0°, B2 −90°", "3.294", "0.925", "0.155 W",
         "0.06255 W = 0.06254 W", "0.475 µJ", "1.1×10<super>−9</super> of 6.69×10<super>−6</super> J"],
        ["standing", "A1 +0°, A2 −90°, B1 −0°, B2 +90°", "3.328", "0.018", "0.000163 W",
         "0.03136 W = 0.03135 W", "0.238 µJ", "4.3×10<super>−10</super> of 3.43×10<super>−6</super> J"],
        ["cancel", "A1 +0°, A2 −90°, B1 +180°, B2 +90°", "0.167", "0.904", "0.0003 W",
         "0.000161 W = 0.0001609 W", "0.00122 µJ", "−4.5×10<super>−12</super> of 2.06×10<super>−7</super> J"],
    ], [0.09, 0.18, 0.07, 0.09, 0.1, 0.17, 0.1, 0.2])
    r.small("min/max: ratio of the smallest to the largest envelope amplitude around the ring; east flow: mean "
            "eastward power flow.")
    r.p("The two towers share the load equally in the traveling and cancel states (A 0.0313 W, B 0.0313 W; and "
        "8.05×10<super>−5</super> W each) and nearly equally when standing (A 0.0155 W, B 0.0159 W).")

    r.h1("Results: four towers")
    r.p("A is at 0°, C at 90°, B at 180°, D at 270°; each pair is an east launcher and a west launcher. |V| is the "
        "steady-state amplitude at the quarter points.")
    r.table([
        ["state", "|V| 0° [V]", "90°", "180°", "270°", "peak", "towers", "loss", "stored", "in step (spread)"],
        ["A–B standing", "3.012", "3.324", "3.328", "3.008", "3.328", "0.03136 W", "0.03135 W", "0.238 µJ", "yes (10 %)"],
        ["C–D standing", "3.008", "3.012", "3.324", "3.328", "3.328", "0.03136 W", "0.03135 W", "0.238 µJ", "yes (10 %)"],
        ["both, synchronized", "6.023", "6.339", "6.656", "6.340", "6.656", "0.1252 W", "0.1252 W", "0.951 µJ", "yes (10 %)"],
        ["both, C–D antiphase", "0.004", "0.315", "0.004", "0.323", "0.323", "0.0002405 W", "0.0002415 W", "0.00183 µJ", "no (99 %)"],
        ["both, C–D shifted 45°", "5.127", "5.433", "5.669", "5.338", "6.088", "0.1069 W", "0.1069 W", "0.812 µJ", "yes (10 %)"],
    ], [0.19, 0.08, 0.07, 0.07, 0.07, 0.07, 0.11, 0.11, 0.11, 0.12], bold_rows=(3,))
    r.small("Voltages in volts. “In step” means the voltages at the four quarter points swing with the same sign; "
            "spread is their relative difference.")
    r.p("Synchronized pairs give twice the amplitude and four times the stored energy and tower power. In antiphase "
        "C–D cancels A–B at the tower points, leaving about 0.3 V at 90° and 270°. Shifting C–D by an eighth of a "
        "wavelength moves the antinodes off the tower points: the peak (6.088 V) is no longer at a quarter point.")

    r.h1("Results: eight towers")
    r.p("A 0°, E 45°, C 90°, G 135°, B 180°, F 225°, D 270°, H 315°.")
    r.table([
        ["state", "|V| at 0°, 45°, 90°, 135°, 180°, 225°, 270°, 315° [V]", "towers = loss", "stored"],
        ["E–F standing", "2.99, 2.99, 3.30, 3.30, 3.30, 3.31, 2.99, 2.98", "0.03136 W = 0.03135 W", "0.238 µJ"],
        ["eight towers, synchronized", "11.95, 12.26, 12.58, 12.90, 13.22, 12.90, 12.58, 12.26", "0.5006 W = 0.5005 W", "3.8 µJ"],
        ["eight towers, 45° set in antiphase", "0.00, 0.32, 0.00, 0.32, 0.00, 0.32, 0.00, 0.32", "0.0002414 W = 0.0002414 W", "0.00232 µJ"],
    ], [0.26, 0.42, 0.2, 0.12], bold_rows=(2,))

    r.h1("Results: DC bias per pair")
    r.p("Synchronized RF, A–B supply 100 V, C–D supply 50 V, 1 MΩ feeds. The ring is one conductor at DC, so on a "
        "continuous ring the pairs cannot hold different levels: it settles at the average and "
        "(V<sub>AB</sub> − V<sub>CD</sub>)/R<sub>feed</sub> flows between the supplies.")
    r.table([
        ["ring", "DC at 0° / 90° / 180° / 270°", "RF on top", "peak", "DC current between supplies", "dissipated in feeds"],
        ["continuous", "75.0 / 75.0 / 75.0 / 75.0 V", "5.97 / 6.29 / 6.60 / 6.29 V", "81.7 V", "50 µA", "2.5 mW"],
        ["sectored (10 µF caps at 45°, 135°, 225°, 315°)", "100.0 / 50.0 / 100.0 / 50.0 V", "5.97 / 6.29 / 6.60 / 6.29 V",
         "106.7 V", "2.8×10<super>−8</super> µA", "1.4×10<super>−9</super> mW"],
    ], [0.22, 0.2, 0.18, 0.09, 0.16, 0.15])
    r.p("In both rings the RF rides unchanged on top of the DC level; the blocking capacitors let each sector hold "
        "its own bias while the RF passes.")

    r.h1("Animations")
    r.p("<i>viz/</i> holds self-contained HTML pages: <i>four_towers.html</i> (each pair alone, synchronized, "
        "antiphase, shifted), <i>phase_control.html</i> (sliders for timing and position phase with the quarter-point "
        "voltages and tower power), <i>eight_towers.html</i> (all four pairs with on/off per pair) and "
        "<i>dc_bias.html</i> (continuous and sectored ring). They are interactive, so they are not reproduced here. "
        "Each embeds <i>launchers.json</i>, the steady-state phasors of each launcher simulated alone; the ring is "
        "linear, so the browser sums them with the chosen phases, matching full simulations to 0.25 % "
        "(TestFourTowers, TestBias).")

    r.h1("Ledger and validation")
    r.bullets([
        "In every state the towers' delivered power equals the ring's loss (e.g. 0.06255 W = 0.06254 W traveling); "
        "the ledger mismatch is at most 1.1×10<super>−9</super> of the energy put in.",
        "<i>main_test.go</i> checks two, four and eight towers and the bias, each against the expected pattern with "
        "the ledger closed.",
        "The DC levels and the current between supplies are checked against the kronspice operating point (TestBias).",
        "Netlists, KiCad schematics and results sheets for the four-tower RF and sectored-bias cases are in "
        "<i>export/</i>.",
    ])

    r.h1("Limits of the model")
    r.bullets([
        "This is a thought experiment scaled to MHz: the ring is a uniform lumped LC ladder with a declared loss "
        "of 0.20 Np per lap, standing in for the Earth–ionosphere ring.",
        "Each tower is two ideal 10 mA current sources (declared values).",
        "The patterns depend only on the number of wavelengths around the ring and the phases, so they carry over to "
        "other scales only in that sense. The same experiment inside a craft's double wall is "
        "<i>examples/double_hull</i>.",
    ])

    r.h1("Reproducing the results")
    r.code("go run ./examples/ring_towers [-json out.json] [-phasors viz/launchers.json]\n"
           "go test ./examples/ring_towers\n"
           "go run ./examples/ring_towers -phasors examples/ring_towers/viz/launchers.json\n"
           "python3 examples/ring_towers/viz/build_viz.py examples/ring_towers/viz/launchers.json examples/ring_towers/viz")
    r.small("Sources: C. P. Steinmetz, <i>Elementary Lectures on Electric Discharges, Waves and Impulses</i>, Lecture "
            "VI (1911); C. R. Paul, <i>Analysis of Multiconductor Transmission Lines</i>, 2nd ed., §2.")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "ring_towers_report.pdf"))
