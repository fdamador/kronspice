# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""H13 report: Tewari's Space Power Generator test through the energy ledger, as a PDF.

python examples/space_power_generator/viz/make_report.py examples/space_power_generator/space_power_generator_report.pdf

Numbers come from examples/space_power_generator/README.md (the example has
no run_output.txt and no figures) and the H13 register entry in CLAUDE.md.
Layout from internal/pdfreport/kronreport.py.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401


def build(out):
    r = Report("examples/space_power_generator", "Tewari's Space Power Generator test",
               "Does the published loaded test show an efficiency above one? The energy ledger of UPSM pp. 159–160 "
               "— register entry H13", "register entry H13", "published test values (UPSM §6.7)")

    r.h1("Summary")
    r.p("P. Tewari's Space Power Generator (SPG) is a homopolar machine: an electromagnet on a rotating iron core, a "
        "static yoke, and output brushes on slip rings, with the Faraday-disk EMF V = ½ωBR² (UPSM eq. 6.16; ω the "
        "angular speed, B the field, R the radius). His book reports one loaded test and reads it as an efficiency "
        "above one. This example runs that test, number by number, through the energy ledger.")
    r.bullets([
        "<b>The arithmetic holds</b> apart from a typo: 3.41 V × 3438 A = 11 724 W (eq. 6.22 prints 11 727 W); the "
        "ratios 2.065 and 2.49 follow.",
        "<b>Counting all inputs</b> (drive motor 8620 W + excitation 650 W = 9270 W) the claimed ratio falls to 1.26; "
        "the book leaves the excitation out.",
        "<b>The whole surplus rests on one assumption:</b> that the EMF under 3438 A load equals the 3.41 V measured "
        "at no load. Energy balance with the book's own shaft-power figures allows at most 1.37 V (1.65 V even with a "
        "lossless drive motor).",
        "<b>Torque:</b> the drive motor delivered 16.7 N·m extra under load; the claimed 3.41 V at 3438 A needs "
        "41.5 N·m of back torque.",
        "<b>Status: untested.</b> The book reports no loaded-EMF measurement. The decisive measurement is the sum of "
        "voltage drops around the loaded loop, or a torque transducer on the shaft.",
    ])

    r.h1("The claim and the standard prediction")
    r.p("<b>Hypothesis (H13, SPG part).</b> Tewari's reactionless generator delivers more electrical power than the "
        "mechanical power that drives it; the book calls the result a violation of energy conservation (p. 161).")
    r.p("<b>Standard prediction.</b> The mechanical power the rotor takes from the shaft equals E<sub>load</sub>·I "
        "(Poynting's theorem; the back torque is k·I when E = k·ω, wherever the reaction lands). Under load the EMF "
        "falls through armature reaction, saturation and speed sag at the rotor: the load loop carries 3438 "
        "ampere-turns against the excitation.")

    r.h1("The published test (UPSM pp. 159–160)")
    r.p("A dc drive motor (DM) runs the SPG at 2700 r/min.")
    r.table([
        ["condition", "reported"],
        ["no load", "DM input W<sub>0</sub> = 2943 W; excitation 650 W; EMF E = 3.41 V"],
        ["loaded, same speed", "I = 3438 A; DM input W<sub>L</sub> = 8620 W"],
    ], [0.3, 0.7])
    r.p("From these the book computes:")
    r.table([
        ["eq.", "quantity", "value"],
        ["6.20", "W<sub>r</sub> = W<sub>L</sub> − W<sub>0</sub>", "5677 W"],
        ["6.21", "W<sub>E</sub> = E·I", "11 724 W"],
        ["6.22", "incremental efficiency", "2.066"],
        ["6.23", "SPG alone, after dividing by the DM's efficiency 0.83", "2.49"],
    ], [0.12, 0.6, 0.28])

    r.h1("The ledger")
    r.h2("1. Arithmetic")
    r.p("3.41 V × 3438 A = 11 724 W; eq. 6.22 prints 11 727 W, a typo. The ratios 2.065 and 2.49 follow.")
    r.h2("2. All inputs")
    r.p("DM 8620 W + excitation 650 W = 9270 W in, against the claimed 11 724 W out: ratio 1.26. The excitation is "
        "left out of the book's ratio.")
    r.h2("3. The assumption")
    r.p("W<sub>E</sub> multiplies the load current by the EMF measured at <i>no load</i>. The book's own numbers give "
        "the extra shaft power under load: 0.83 × 5677 W = 4712 W, a torque increase of 16.7 N·m at 282.7 rad/s. "
        "Since the rotor's shaft power must equal E<sub>load</sub>·I, this bounds the loaded EMF:")
    r.table([
        ["", "loaded EMF", "back torque at 3438 A"],
        ["standard prediction", "≤ 1.37 V (≤ 1.65 V even with a lossless DM)", "16.7 N·m, as the drive motor delivered"],
        ["the claim", "3.41 V, unchanged from no load", "needs k·I = 41.5 N·m"],
    ], [0.22, 0.4, 0.38], bold_rows=(1,))
    r.p("So the standard prediction is that the EMF falls by at least 60 % at 3438 A. The implied loop resistance "
        "E<sub>0</sub>/I is 0.992 mΩ. The book reports no loaded-EMF measurement; it says the drops “can also be” "
        "measured (p. 160).")

    r.h1("Validation")
    r.bullets([
        "<i>main_test.go</i> reproduces the book's numbers and the bound (1.3706 V, and 1.6513 V with a lossless "
        "drive motor), with torque consistency to 10<super>−12</super>.",
        "<i>main.go</i> holds the published values and the analysis of eqs. 6.20–6.23, all inputs and the loaded-EMF "
        "bound.",
    ])

    r.h1("Limits of the analysis")
    r.bullets([
        "It uses only the published figures; the shaft power under load is inferred from the drive motor's input and "
        "the book's assumed efficiency 0.83, not measured.",
        "It is a ledger of one reported test, not a field model of the machine.",
        "Whether the shaft bearings could supply the missing torque is asked in <i>examples/bearing_motor</i>; the "
        "register records that bearing (Huber) torque could close the 24.8 N·m torque gap only in EM mode, whose "
        "back-EMF leaves the same 1.37 V, and is excluded in thermal mode at the measured η ≈ 0.10. The electron "
        "model of the same hypothesis is on the bench in <i>examples/space_vortex</i>.",
    ])

    r.h1("Status")
    r.p("Register status: <b>untested</b>. The claim rests on an unmeasured loaded EMF: ≤ 1.37 V predicted, 3.41 V "
        "assumed. <b>Decisive measurement, not in the book:</b> the sum of the voltage drops around the loaded loop "
        "(brushes, rings, conductor, shunt, load) at 3438 A, or the shaft torque from a transducer rather than "
        "inferred from the drive motor's input and an assumed efficiency. A closed-loop, self-running test would "
        "settle it outright.")

    r.h1("Reproducing the results")
    r.code("go run ./examples/space_power_generator\n"
           "go test ./examples/space_power_generator")
    r.small("Sources: P. Tewari, <i>Universal Principles of Space and Matter</i>, Grain and Chaff, Toronto, 2017 (1st "
            "ed. 2009), §6.7, pp. 149–161 (hypothesis source); P. Tewari, Physics Essays 31, 108 (2018); reading notes "
            "in doc/tewari-references.md.")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "space_power_generator_report.pdf"))
