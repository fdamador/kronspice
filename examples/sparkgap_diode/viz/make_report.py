# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Report: a water gap that acts as a diode (spark-gap diode and diode bridge), as a PDF.

python examples/sparkgap_diode/viz/make_report.py examples/sparkgap_diode/sparkgap_diode_report.pdf

Numbers come from examples/sparkgap_diode/run_output.txt (the program's
output) and examples/sparkgap_diode/README.md; figures are
viz/sparkgap_diode.png, viz/pulses.png and viz/frequency_trend.png.
Layout from internal/pdfreport/kronreport.py.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401


def build(out):
    r = Report("examples/sparkgap_diode", "A water gap as a diode",
               "Aluminium valve cells and polarity spark gaps as rectifiers — CLAUDE.md §3g, §3h",
               "CLAUDE.md §3g, §3h", "declared cell and gap parameters")

    r.h1("Summary")
    r.p("Can a water gap act like a diode because of its two electrode materials? Yes, and it was done over a "
        "century ago. A <b>valve metal</b> (aluminium) in water grows an insulating oxide when it is the positive "
        "electrode and conducts freely when it is the negative one. Paired with an ordinary metal (lead, iron) it is "
        "the aluminium electrolytic rectifier of Pollak (1896) and Nodon (1899): about 140 V per cell, about "
        "75–80 % efficient (<i>Encyclopaedia Britannica</i> 1911, “Transformers”). At kilovolts the water itself "
        "breaks down in a spark; there the one-way action comes from the <b>shape</b> of the electrodes (a sharp "
        "positive point breaks down at a lower voltage than a sharp negative one), the polarity spark gap. The "
        "example models both, builds rectifiers from them, and compares them in a pulse test and a frequency sweep.")
    r.bullets([
        "<b>Aluminium against any ordinary metal rectifies:</b> forward/reverse 1.37×10<super>4</super> at 100 V for "
        "lead, iron or silver-plated iron against aluminium; two ordinary metals are symmetric; aluminium against "
        "aluminium blocks both ways.",
        "<b>Best rectifier:</b> a bridge of four aluminium–lead cells at 100 V peak gives 88.0 V with 5.21 V ripple "
        "at 85.0 % efficiency, at the top of the 1911 account's 75–80 %.",
        "<b>The oxide film limits frequency.</b> It is a 4.38 µF capacitor; the bridge falls to 17.2 % efficiency at "
        "20 kHz (84.1 % at 1 kHz).",
        "<b>Spark-gap bridges need a polarity ratio above 2.</b> With the declared ratio 2 the bridge never strikes "
        "(4 kV) or breaks down backwards (6 kV); with ratio 3 it rectifies at 67.9 %. Both breakdown voltages are "
        "extrapolated and the ratio is declared.",
        "<b>Fastest amps:</b> the same 10 J gives 4247 A in 0.108 µs through a 3 mm water gap, against 470 A in "
        "9.22 µs through a lead–aluminium cell, but the cell hands over all the energy and the gap leaves 8.10 J in "
        "the capacitor.",
        "<b>Spark gaps work only at mains frequency:</b> the water-gap bridge gives 0.247 A/V at 50 Hz and 0 from "
        "300 Hz up; the cells hold about 0.459 A/V to 100 kHz, but their efficiency falls to 34 %.",
    ])

    r.h1("Model")
    r.bullets([
        "<b>Water cell</b> (<i>cell.go</i>): plates 100 cm², 1 cm apart, water 0.1 Ω (phosphate solution), "
        "decomposition onset 2.0 V; aluminium's oxide blocks the reverse direction and sparks over at 140 V; the "
        "film is a 4.38 µF capacitor. Each I–V curve is checked against a kronspice DC sweep of the same element.",
        "<b>Rectifiers</b> (<i>rectifier.go</i>): kronspice .TRAN with the ledger on; load 50 Ω with a 2,200 µF "
        "filter, source resistance 0.5 Ω.",
        "<b>Polarity spark gap</b> (<i>gap.go</i>, <i>gapsim.go</i>): 3 mm point–plane gap in pure water, breakdown "
        "+2415 V with the point positive (extrapolated from Stygar et al. 2006, sub-microsecond data, to 60 Hz) and "
        "−4830 V with it negative (ratio 2, declared from J. C. Martin's constants, to be checked); leak 57.9 MΩ; "
        "arc 50 V + 0.5 Ω. A spark gap is a switch with memory, so it is stepped by a small nodal solver of its own, "
        "which matches kronspice below breakdown (TestGapSolverTwin).",
        "<b>Recovery:</b> declared as 1 − e<super>−t/1 ms</super>; a hot gap never re-strikes below the Paschen "
        "minimum of water vapour, 480 V (Škoro et al. 2011).",
    ])

    r.h1("Part 1 — electrode pairs in water")
    r.p("Current from the first electrode to the second.")
    r.table([
        ["pair", "+5 V", "−5 V", "+100 V", "−100 V", "−200 V", "fwd/rev at 100 V", "twin"],
        ["lead | water | aluminium", "30 A", "−0.00357 A", "980 A", "−0.0714 A", "−28.7 A", "1.37×10<super>4</super>", "6.8e-17"],
        ["iron | water | aluminium", "30 A", "−0.00357 A", "980 A", "−0.0714 A", "−28.7 A", "1.37×10<super>4</super>", "6.8e-17"],
        ["copper with silver tips | water | aluminium", "3 A", "−0.00357 A", "98 A", "−0.0714 A", "−20.1 A", "1.37×10<super>3</super>", "3.1e-16"],
        ["silver-plated iron | water | aluminium", "30 A", "−0.00357 A", "980 A", "−0.0714 A", "−28.7 A", "1.37×10<super>4</super>", "6.8e-17"],
        ["iron | water | copper", "30 A", "−30 A", "980 A", "−980 A", "−1.98×10<super>3</super> A", "—", "0.0"],
        ["lead | water | iron", "30 A", "−30 A", "980 A", "−980 A", "−1.98×10<super>3</super> A", "—", "0.0"],
        ["aluminium | water | aluminium", "0.00357 A", "−0.00357 A", "0.0714 A", "−0.0714 A", "−28.7 A", "—", "4.7e-15"],
    ], [0.22, 0.1, 0.11, 0.1, 0.1, 0.13, 0.12, 0.1])
    r.p("A valve metal against an ordinary metal rectifies; any ordinary metal, silver included, works as the "
        "partner, because the one-way action comes from the aluminium's oxide. Two ordinary metals do not rectify; "
        "they make a weak battery (iron | water | copper 0.789 V open circuit; lead | water | iron −0.321 V). Two "
        "valve metals block both ways. Past the oxide's sparking voltage (140 V) the blocked direction sparks over, "
        "like a Zener diode.")

    r.h1("Part 2 — aluminium–lead water cells as rectifiers")
    r.table([
        ["circuit", "Vdc", "ripple", "load", "efficiency", "reverse peak", "ledger"],
        ["half-wave, 1 cell, 100 V peak, 60 Hz", "76.2 V", "27.28 V", "118.3 W", "20.5 %", "167 V (sparks over)", "4.6e-09"],
        ["half-wave, 2 cells in series, 100 V peak", "83.3 V", "10.91 V", "139.0 W", "84.0 %", "183 V", "−1.6e-09"],
        ["bridge, 1 cell per arm, 100 V peak, 60 Hz", "88.0 V", "5.21 V", "154.8 W", "85.0 %", "93 V", "1.1e-08"],
        ["bridge, 1 cell per arm, 170 V peak (120 V rms)", "145.1 V", "9.69 V", "421.1 W", "40.1 %", "153 V (sparks over)", "1.2e-07"],
        ["bridge, 2 cells per arm, 170 V peak", "146.5 V", "8.36 V", "429.4 W", "85.6 %", "156 V", "3.5e-09"],
        ["bridge, 1 cell per arm, 100 V peak, 1 kHz", "88.3 V", "0.32 V", "155.9 W", "84.1 %", "91 V", "2.5e-08"],
        ["bridge, 1 cell per arm, 100 V peak, 20 kHz", "85.9 V", "0.02 V", "147.5 W", "17.2 %", "89 V", "2.9e-09"],
        ["bridge, copper with silver tips | aluminium", "79.2 V", "4.04 V", "125.6 W", "77.5 %", "88 V", "7.9e-09"],
        ["bridge, silver-plated iron | aluminium", "88.0 V", "5.21 V", "154.8 W", "85.0 %", "93 V", "1.1e-08"],
    ], [0.34, 0.09, 0.09, 0.1, 0.1, 0.17, 0.11], bold_rows=(3,))
    r.small("Efficiency = load power ÷ source power over the last cycles. 60 Hz unless stated.")
    r.bullets([
        "A filtered half-wave rectifier sees twice the peak in reverse, so one cell sparks over; two in series hold, "
        "as in the historical practice of “two or more valves in series for higher pressures”.",
        "At 20 kHz the oxide film's capacitance passes current both ways and the efficiency collapses, the 1911 "
        "article's warning that these rectifiers are unsuitable for high frequencies.",
        "85 % sits at the top of the 75–80 % the 1911 account gives; the cell's leakage and resistance are declared.",
    ])
    r.h2("Silver against aluminium")
    r.p("Both silver versions rectify, but when a cell conducts the silver is the positive electrode, and a positive "
        "metal electrode in water is eaten away. If all the current dissolved silver (Faraday's law, the worst case):")
    r.table([
        ["silver electrode (declared size)", "silver", "charge it passes", "lasts at the bridge's current"],
        ["copper with 1 mm solid silver tips on 10 % of the plate", "10.5 g", "9.38×10<super>3</super> C", "3.3 h (0.79 A per cell)"],
        ["iron with 10 µm silver plating", "1.05 g", "938 C", "18 min (0.88 A per cell)"],
    ], [0.45, 0.12, 0.18, 0.25])
    r.p("After that the copper or iron underneath carries the current and corrodes. This is why Pollak used lead: as "
        "the positive electrode it forms a stable lead-dioxide skin instead of dissolving.")

    r.h1("Part 3 — polarity spark gaps in water as rectifiers")
    r.p("A one-way spark needs a forward voltage above +Vb and a reverse voltage below −Vb·ratio. With a filter "
        "capacitor a half-wave rectifier sees 2·Vpeak in reverse, and a bridge needs 2·Vb forward (two gaps in "
        "series) while each blocked gap sees the full output: both need a polarity ratio above 2.")
    r.table([
        ["circuit", "ratio", "Vdc", "ripple", "efficiency", "reverse peak", "strikes f/b", "result"],
        ["half-wave, resistive load, 4 kV peak", "2.0", "900 V", "4014 V", "96.1 %", "3932 V", "8617/0", "rectifies"],
        ["half-wave, filter capacitor, 4 kV peak", "2.0", "366 V", "6014 V", "7.0 %", "4830 V", "21/20", "breaks down backwards"],
        ["bridge, filter capacitor, 4 kV peak", "2.0", "0 V", "0 V", "0.0 %", "2029 V", "0/0", "never strikes"],
        ["bridge, filter capacitor, 6 kV peak", "2.0", "1837 V", "7173 V", "0.1 %", "9474 V", "540/520", "breaks down backwards"],
        ["bridge, filter capacitor, 6 kV peak", "3.0", "3107 V", "4519 V", "67.9 %", "5965 V", "6/0", "rectifies"],
    ], [0.27, 0.06, 0.08, 0.08, 0.1, 0.11, 0.1, 0.2], bold_rows=(5,))
    r.p("A spark gap is not a low-drop diode: it conducts only once the voltage passes its breakdown (2.4 kV here), "
        "so the output sags until the source exceeds the load voltage by the gaps' breakdown, then jumps back up. "
        "<b>A bench measurement</b> of both breakdown voltages and their ratio decides whether a water spark-gap "
        "bridge can work.")
    r.figure(os.path.join(HERE, "sparkgap_diode.png"),
             "Figure 1. Top left: I–V curves of the seven electrode pairs (symmetric-log current axis); aluminium "
             "against an ordinary metal conducts forward and blocks to 140 V in reverse, aluminium against aluminium "
             "blocks both ways. Top right: output of the bridge of four aluminium–lead cells at 100 V peak. Bottom "
             "left: polarity spark gap, half-wave, resistive load, 4 kV peak (pulsed output). Bottom right: spark-gap "
             "bridge with polarity ratio 3 at 6 kV peak: the output decays until the gaps strike again.",
             max_height_mm=135)

    r.h1("Part 4 — pulse test: which gives the fastest amps?")
    r.p("The same 10 J from a charged capacitor through each device, a 1 µH loop and a 0.1 Ω load (declared).")
    r.table([
        ["device", "charge", "C", "peak", "rise (10–90 %)", "to load", "in device", "left in C"],
        ["lead | water | aluminium", "100 V", "0.002 F", "470 A", "9.22 µs", "4.80 J", "5.19 J", "0.00 J"],
        ["silver-plated iron | water | aluminium", "100 V", "0.002 F", "470 A", "9.22 µs", "4.80 J", "5.19 J", "0.00 J"],
        ["copper with silver tips | water | aluminium", "100 V", "0.002 F", "195 A", "2.22 µs", "0.86 J", "8.82 J", "0.32 J"],
        ["water polarity gap, 3 mm", "41000 V", "1.19×10<super>−8</super> F", "4247 A", "0.108 µs", "0.31 J", "1.59 J", "8.10 J"],
        ["mist-filled air gap, 1 mm", "8000 V", "3.13×10<super>−7</super> F", "3510 A", "0.505 µs", "1.09 J", "5.55 J", "3.36 J"],
    ], [0.3, 0.09, 0.12, 0.08, 0.11, 0.09, 0.1, 0.11], bold_rows=(4,))
    r.p("Spark gaps give the highest and fastest current because they run at kilovolts. The cells give hundreds of "
        "amps in microseconds, limited by the water's resistance, and hand over all the energy. A gap's arc goes out "
        "at the first current zero, so energy still ringing in the capacitor stays there. The water gap's charge "
        "follows Stygar's relation at 1 µs, inside its measured range; the mist gap uses the fog factor from "
        "<i>water_lamp</i> (extrapolated).")
    r.figure(os.path.join(HERE, "pulses.png"),
             "Figure 2. Current through each device for the same 10 J pulse (log time axis, µs). The water gap peaks "
             "at about 4.2 kA and the mist gap at 3.5 kA within a microsecond; the cells rise to hundreds of amps over "
             "about 10 µs.", max_height_mm=95)

    r.h1("Part 5 — frequency sweep: amps per volt")
    r.p("Every device in the same bridge (source 0.1 Ω, load 1 Ω, no filter), 50 Hz to 100 kHz. Amps per volt = "
        "mean load current ÷ source peak voltage. A perfect diode bridge gives 0.579 A/V at every frequency.")
    r.table([
        ["device", "50 Hz", "100", "300", "1k", "3k", "10k", "30k", "100k"],
        ["lead | aluminium (100 V), A/V", "0.459", "0.459", "0.459", "0.459", "0.460", "0.463", "0.481", "0.508"],
        ["   efficiency", "73 %", "73 %", "73 %", "73 %", "73 %", "71 %", "60 %", "34 %"],
        ["silver-plated iron | aluminium (100 V)", "same as lead", "", "", "", "", "", "", ""],
        ["copper with silver tips | aluminium (100 V), A/V", "0.192", "0.192", "0.192", "0.192", "0.193", "0.193", "0.195", "0.191"],
        ["   efficiency", "31 %", "31 %", "31 %", "30 %", "30 %", "29 %", "22 %", "8 %"],
        ["water gap, ratio 3, 6 kV fixed, A/V", "0.247", "0.198", "0", "0", "0", "0", "0", "0"],
        ["   efficiency", "47 %", "47 %", "0 %", "0 %", "0 %", "0 %", "0 %", "0 %"],
        ["water gap, ratio 3, voltage raised with f, A/V", "0.239", "0.240", "0", "0", "0", "0", "0", "0"],
        ["   efficiency", "47 %", "47 %", "0 %", "0 %", "0 %", "0 %", "0 %", "0 %"],
    ], [0.36, 0.1, 0.08, 0.065, 0.065, 0.065, 0.065, 0.065, 0.075])
    r.bullets([
        "<b>Water cells</b> hold steady to about 3 kHz. Above that the oxide film passes current both ways: the mean "
        "current holds up but the efficiency collapses (34 % at 100 kHz).",
        "<b>Water spark gaps</b> work only at mains frequencies. At a fixed 6 kV the breakdown voltage rises with "
        "frequency (E·τ<super>1/3</super>, extrapolated) until, above 100 Hz, the gaps never strike. With the "
        "voltage raised to keep them striking, above about 100 Hz the gaps that just carried an arc have not "
        "recovered and break down backwards: the bridge shorts. Spark gaps are pulse switches, not high-frequency "
        "rectifiers.",
    ])
    r.figure(os.path.join(HERE, "frequency_trend.png"),
             "Figure 3. Amps per volt against frequency for every device in the same bridge. The cells stay flat "
             "(lead and silver-plated iron overlap); the water gaps drop to no output above 100 Hz (shaded band). "
             "Dotted line: perfect diode bridge, 0.58 A/V.", max_height_mm=110)

    r.h1("Ledger and validation")
    r.bullets([
        "Every rectifier run's energy ledger closes to about 10<super>−7</super> or better (largest listed "
        "1.2×10<super>−7</super>).",
        "Each cell I–V curve matches a kronspice DC sweep of the same element (twin ≤ 4.7×10<super>−15</super>).",
        "The spark-gap solver matches kronspice below breakdown (TestGapSolverTwin).",
        "The water and the metals add no energy: every watt at the load comes from the AC source; the rest is heat "
        "and, for the cells, some electrolysis.",
    ])

    r.h1("Limits of the model")
    r.bullets([
        "Cell leakage and resistance are declared values; oxide sparking at 140 V follows the 1911 figure per cell.",
        "The water gap's positive breakdown is extrapolated from sub-microsecond data to 60 Hz; the polarity ratio "
        "of 2 is declared (to be checked); the mist gap's fog factor is extrapolated.",
        "Gap recovery (1 ms) is declared.",
        "Silver lifetimes are a worst case (all current dissolving silver); part of it may evolve oxygen or form a "
        "film instead.",
    ])

    r.h1("Reproducing the results")
    r.code("go run ./examples/sparkgap_diode -viz examples/sparkgap_diode/viz   # a few seconds\n"
           "python3 examples/sparkgap_diode/viz/render.py\n"
           "go test ./examples/sparkgap_diode")
    r.small("Bridges that combine these cells and gaps for particular jobs (surge protection, a spark pulser, a voltage "
            "doubler, the oxide film as a varactor, a frequency doubler) are in <i>examples/sparkgap_diodebridge</i>. "
            "Sources: Encyclopaedia Britannica 11th ed. (1911) vol. 27 “Transformers”; W. A. Stygar et al., Phys. Rev. "
            "ST Accel. Beams 9, 070401 (2006); M. Škoro et al., Phys. Rev. E 84, 055401 (2011); L. Liu et al. (2018); "
            "CRC Handbook electrochemical series.")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "sparkgap_diode_report.pdf"))
