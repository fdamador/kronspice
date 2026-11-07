# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Report: Tesla's two-disk homopolar with its load current through the ball bearings, as a PDF.

python examples/homopolar_bearing/viz/make_report.py examples/homopolar_bearing/homopolar_bearing_report.pdf

Numbers come from examples/homopolar_bearing/run_output.txt and README.md (the
example has no figures). Layout from internal/pdfreport/kronreport.py.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401
import kronreport  # noqa: E402

kronreport.H1.keepWithNext = 1  # keep headings with the text that follows
kronreport.H2.keepWithNext = 1


def build(out):
    r = Report("examples/homopolar_bearing", "Homopolar generator through its ball bearings",
               "Least drive torque per load watt when Tesla's two-disk homopolar takes its current through the shaft "
               "bearings — CLAUDE.md §3f, register entries H13, H14", "CLAUDE.md §3f, H13, H14",
               "15 cm copper disks, 0.8 T; GA 702 evaluations per case")

    r.h1("Summary")
    r.p("Tesla's two-disk homopolar generator (U.S. Patent 406,968) is built with its output current taken through "
        "the shaft ball bearings instead of brushes, so that each bearing is also a Huber ball-bearing motor (the P "
        "element of <i>bearing_motor</i>). This is a circuit simulation only: every part is a lumped element and the "
        "shafts are mechanical nodes (across variable ω, through variable torque). The question is the least drive "
        "torque per watt delivered to an external resistor or capacitor. No coefficient is assumed; unmeasured ones "
        "are swept.")
    r.bullets([
        "<b>Reference unit</b> (40 mm bearings × 2 per shaft, 600 rad/s, 2 Ω): 41.80 W into the load from 77.29 W at "
        "the shaft, efficiency 0.541, 0.00308 N·m/W — but the bearing current is 7.62× the 0.3 A/mm² damage limit.",
        "<b>At the damage limit</b> a bearing bank passes only 0.238–1.248 A (bores 10–120 mm), so a through-bearing "
        "unit is a watts-to-tens-of-watts machine (2.22–9.50 W).",
        "<b>The Huber torque never helps the torque per watt:</b> in EM mode it lowers the drive torque only by lowering "
        "the output through its back-EMF; in thermal mode it leaves 41.797 W at 0.0031 N·m/W.",
        "<b>Insulated bearings do better:</b> 56.26 W at efficiency 0.661 with a shaft collector, against 41.80 W "
        "at 0.541 through the bearings.",
        "<b>Evolved:</b> through the bearings, 21.13 W at 0.0101 N·m/W (efficiency 0.134); bearings insulated, "
        "689.28 W at 0.00135 N·m/W (efficiency 0.998), 7.4× less torque per watt.",
        "<b>Ledger:</b> residual 0 at every operating point; no free torque (shaft power is never below load plus "
        "measured losses).",
    ])

    r.h1("Background")
    r.p("Register entries H13 and H14 concern rotating machines whose ledgers have been questioned: Tewari's "
        "space-power generator (H13), whose 24.8 N·m torque gap was compared with the bearing (Huber) torque, and "
        "Marinov's claim that a running ball-bearing motor yields heat beyond its electrical input (H14). The P "
        "element models the Huber effect from measured data only. This example asks a practical engineering "
        "question with the same element: if the homopolar current must pass through the bearings, what does it cost "
        "in drive torque? The standard prediction for every design is that the drive torque per load watt is at least "
        "1/ω, with equality only at efficiency 1.")

    r.h1("The trade-off and the measured inputs")
    r.p("The bearing bore sets the axle radius r<sub>i</sub>, so bigger bearings remove magnetic area from each disk, "
        "EMF = ½Bω(r<sub>o</sub>² − r<sub>i</sub>²) (B: axial field, ω: shaft speed, r<sub>o</sub>: disk radius). But "
        "bigger, more heavily loaded bearings carry more current before they are damaged, and have more friction.")
    r.table([
        ["input", "value", "source"],
        ["bearing voltage under DC", "almost constant, 0.73 V average at 4 A (0.5–1.5 V reported), from repeated "
         "breakdown of the oxide and lubricant film; carried by the P element as VF",
         "Safdarzadeh, Farahi, Binder et al., Lubricants 12, 148 (2024); greased 51208 bearings, 40 mm bore"],
        ["damage limit", "0.1–0.3 A/mm² of Hertzian contact area; A<sub>Hz</sub> = 1.0 mm² at 400 N and 2.5 mm² at "
         "1600 N; other sizes by Hertz theory, A proportional to (F·R)<super>2/3</super>", "same paper; Johnson, Contact Mechanics"],
        ["friction torque", "M = ½·f·F·d with f = 0.0015 (F: load, d: bore)", "SKF simplified model"],
        ["rim-speed limit", "hoop stress (3+ν)/8·ρv² at the yield of annealed copper, 69 MPa, safety factor 1.5: "
         "111 m/s", "Timoshenko; ASM Handbook"],
        ["swept (unmeasured)", "Huber coefficient κψ in EM and thermal modes; belt rolling-contact resistance; shaft "
         "collector resistance", "—"],
    ], [0.18, 0.52, 0.3])
    r.p("<i>TestHertzArea</i> recovers the paper's second contact-area point from the first within 2 %.")

    r.h1("Results (15 cm copper disks, 0.8 T)")
    r.h2("Reference unit")
    r.p("Bore 40 mm × 2 per shaft, 400 N per bearing, disk r<sub>o</sub> 0.150 m, r<sub>i</sub> 0.020 m, thickness "
        "5.0 mm, ω 600 rad/s (rim 90 m/s), load R<sub>L</sub> 2 Ω. EMF 10.608 V; load 4.571 A, 41.80 W; shaft "
        "77.29 W; efficiency 0.541; drive torque 0.0644 + 0.0644 N·m = 0.00308 N·m/W. Bearings: contact heat 6.70 W, "
        "friction 0.0240 N·m per shaft, Huber torque 0. The bearing current is 7.62× the 0.3 A/mm² limit (22.86× the "
        "0.1 A/mm² one). Ledger residual 0.")
    r.h2("Bearing bore against disk magnetic area (load matched to the bank's current limit)")
    r.table([
        ["bore", "area kept", "EMF", "I<sub>max</sub>", "P<sub>load</sub>", "N·m/W", "current / limit"],
        ["10 mm", "99.9 %", "10.788 V", "0.238 A", "2.22 W", "0.0073", "1.00"],
        ["20 mm", "99.6 %", "10.752 V", "0.378 A", "3.51 W", "0.0088", "1.00"],
        ["40 mm", "98.2 %", "10.608 V", "0.600 A", "5.49 W", "0.0107", "1.00"],
        ["60 mm", "96.0 %", "10.368 V", "0.786 A", "7.00 W", "0.0122", "1.00"],
        ["80 mm", "92.9 %", "10.032 V", "0.952 A", "8.16 W", "0.0137", "1.00"],
        ["100 mm", "88.9 %", "9.600 V", "1.105 A", "8.99 W", "0.0153", "1.00"],
        ["120 mm", "84.0 %", "9.072 V", "1.248 A", "9.50 W", "0.0171", "1.00"],
    ], [0.12, 0.13, 0.14, 0.13, 0.13, 0.13, 0.22])
    r.h2("Huber torque, κψ swept (reference unit)")
    r.table([
        ["κψ (N·m/A²)", "mode", "Huber torque", "drive torque", "P<sub>load</sub>", "N·m/W"],
        ["0", "EM / thermal", "0.0000 N·m", "0.1288 N·m", "41.797 W", "0.0031"],
        ["10<super>−3</super>", "EM", "0.0062 N·m", "0.0978 N·m", "24.738 W", "0.0040"],
        ["10<super>−3</super>", "thermal", "0.0005 N·m", "0.1278 N·m", "41.797 W", "0.0031"],
        ["10<super>−2</super>", "EM", "0.0065 N·m", "0.0551 N·m", "2.614 W", "0.0211"],
        ["10<super>−2</super>", "thermal", "0.0006 N·m", "0.1277 N·m", "41.797 W", "0.0031"],
        ["10<super>−1</super>", "EM", "0.0011 N·m", "0.0484 N·m", "0.044 W", "1.1124"],
        ["10<super>−1</super>", "thermal", "0.0006 N·m", "0.1277 N·m", "41.797 W", "0.0031"],
    ], [0.17, 0.17, 0.17, 0.17, 0.16, 0.16])
    r.p("In EM mode the Huber torque lowers the drive torque only by lowering the output through its back-EMF, so it "
        "never improves the torque per watt. In thermal mode it recovers at most η·H, about 10 % of the bearing heat. "
        "<i>TestNoFreeTorque</i> checks that the shaft power is never below the load plus the measured losses.")
    r.h2("Bearing-voltage range and the insulated alternative")
    r.table([
        ["configuration", "P<sub>load</sub>", "efficiency", "N·m/W"],
        ["through bearings, 0.50 V each", "46.11 W", "0.578", "0.0029"],
        ["through bearings, 0.73 V each", "41.80 W", "0.541", "0.0031"],
        ["through bearings, 1.50 V each", "28.91 W", "0.418", "0.0040"],
        ["insulated bearings, collector 10<super>−5</super> Ω", "56.26 W", "0.661", "0.0025"],
        ["insulated bearings, collector 10<super>−4</super> Ω", "56.25 W", "0.661", "0.0025"],
        ["insulated bearings, collector 10<super>−3</super> Ω", "56.15 W", "0.661", "0.0025"],
    ], [0.46, 0.18, 0.18, 0.18])
    r.h2("Capacitor load: 1 F through R<sub>L</sub>")
    r.p("The field is stepped (1 ms) or ramped, then the circuit settles for 5·R<sub>L</sub>·C.")
    r.table([
        ["field ramp", "stored", "shaft", "friction", "conversion", "mismatch"],
        ["0 s (stepped)", "41.28 J", "384.36 J", "288.00 J", "42.8 %", "−1.6×10<super>−3</super> J"],
        ["20 s", "41.80 J", "927.73 J", "864.00 J", "65.6 %", "7.3×10<super>−5</super> J"],
        ["100 s", "41.96 J", "3225.22 J", "3168.00 J", "73.3 %", "7.0×10<super>−5</super> J"],
    ], [0.18, 0.15, 0.17, 0.17, 0.15, 0.18])
    r.p("The bearing drops cost 2·V<sub>b</sub>·Q whatever the ramp (V<sub>b</sub>: voltage per bearing, Q: charge "
        "delivered), so electrical conversion is at best ½CV²/(½CV² + 2V<sub>b</sub>CV) = 75.8 % (V = 9.15 V), while "
        "friction runs throughout.")
    r.h2("Evolved for the least drive torque per watt at ≥ 20 W")
    r.p("Genes: bore, preload, disk thickness, speed up to the rim limit, load and bearings per shaft; limits: the "
        "bearing current limit and the rim speed; 702 evaluations each.")
    r.table([
        ["", "through the bearings", "bearings insulated (collector 5×10<super>−5</super> Ω declared)"],
        ["design", "108 mm × 6 per shaft, 185 N each, r<sub>i</sub> 0.054 m, t 4.4 mm",
         "10 mm × 1 per shaft, 50 N each, r<sub>i</sub> 0.005 m, t 6.4 mm"],
        ["speed; load", "739 rad/s (rim 111 m/s); 4.84 Ω", "739 rad/s (rim 111 m/s); 0.256 Ω"],
        ["EMF", "11.575 V", "13.292 V"],
        ["load", "21.13 W at 2.089 A (1.00× the current limit)", "689.28 W at 51.897 A"],
        ["shaft power; efficiency", "157.62 W; 0.134", "690.38 W; 0.998"],
        ["drive torque per watt", "<b>0.0101 N·m/W</b>", "<b>0.00135 N·m/W</b> (7.4× less)"],
        ["ledger residual", "0", "0"],
    ], [0.24, 0.36, 0.4], bold_rows=(6,))
    r.p("The insulated unit's torque per watt is 1/(η·ω); in every case the drive torque per load watt is at least "
        "1/ω.")

    r.h1("Ledger and validation")
    r.bullets([
        "Ledger residual 0 at every operating point; transient mismatch below 10<super>−4</super> of the shaft energy.",
        "<i>main_test.go</i>: Hertz scaling (<i>TestHertzArea</i>), no free torque (<i>TestNoFreeTorque</i>), bore "
        "monotonicity, the capacitor conversion ceiling.",
        "<i>export_test.go</i>: reference and capacitor netlists, KiCad schematics and results sheets in <i>export/</i>.",
    ])

    r.h1("Limits and bench data still needed")
    r.bullets([
        "Circuit simulation only: lumped elements, shafts as mechanical nodes.",
        "The Huber coefficient κψ, the belt's rolling-contact resistance and the real collector's contact drop are "
        "unmeasured and swept. A carbon brush drops about 1 V, which would narrow the insulated design's lead; a "
        "liquid-metal collector would not.",
        "The damage limit is given as a range (0.1–0.3 A/mm²); the designs use 0.3 A/mm².",
    ])

    r.h1("Reproducing the results")
    r.code("go run ./examples/homopolar_bearing [-p 20] [-gens 25] [-pop 30]\n"
           "go test ./examples/homopolar_bearing\n"
           "python examples/homopolar_bearing/viz/make_report.py "
           "examples/homopolar_bearing/homopolar_bearing_report.pdf")
    r.small("Sources: O. Safdarzadeh, A. Farahi, A. Binder et al., Lubricants 12(5), 148 (2024), DOI "
            "10.3390/lubricants12050148; K. L. Johnson, Contact Mechanics, Cambridge 1985, §4.2; T. A. Harris, Rolling "
            "Bearing Analysis; SKF rolling-bearing catalogue; S. Timoshenko, Strength of Materials II, §44; ASM "
            "Handbook vol. 2 (C10100 annealed); N. Tesla, U.S. Patent 406,968 (1889); doc/bearing-motor.md.")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "homopolar_bearing_report.pdf"))
