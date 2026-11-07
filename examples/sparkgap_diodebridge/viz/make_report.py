# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""PDF report for examples/sparkgap_diodebridge: water cells and spark gaps combined for six jobs.

Run from the repository root:

    python examples/sparkgap_diodebridge/viz/make_report.py examples/sparkgap_diodebridge/sparkgap_diodebridge_report.pdf

Numbers come from examples/sparkgap_diodebridge/run_output.txt (the program's
output) and README.md (circuits, declared values, sources); the figure is
viz/diodebridge.png. Layout from internal/pdfreport/kronreport.py.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401


def build(out):
    r = Report("examples/sparkgap_diodebridge", "Water diodes and spark gaps in bridges",
               "What useful circuits can be built from aluminium water cells, polarity spark gaps in water "
               "and air spark gaps? — CLAUDE.md §3g, §3h",
               "CLAUDE.md §3g, §3h")

    r.h1("Summary")
    r.p("The companion example <i>sparkgap_diode</i> built two one-way devices from water: the <b>water cell</b> "
        "(aluminium against lead or iron in a phosphate solution, the 1890s rectifier of Pollak and Nodon, one-way up "
        "to about 140 V) and the <b>polarity spark gap</b> (a point–plane gap in water that breaks down at a lower "
        "voltage with the point positive). This example combines them, with plain air spark gaps, into circuits for "
        "six jobs, each with an energy ledger.")
    r.bullets([
        "<b>A. DC supply:</b> four lead–aluminium cells give 87.9 V DC from 100 V peak, 85 % efficient; kronspice "
        "gives 87.96 V for the same circuit (−0.01 %).",
        "<b>B. Surge protection:</b> under a 4 kV, 1.2/50 µs surge the largest bridge-cell voltage falls from 221 V "
        "(unprotected, above the 140 V spark-over) to 117 V with an air-gap arrester and an aluminium–aluminium clamp.",
        "<b>C. Spark pulser:</b> only an output gap fired at each mains zero crossing works: 120 pulses/s, 505 A peak, "
        "1.13 J per pulse, 50.3 % of mains power to the load. Self-breaking gaps reach 0.6 % and 1.8 %.",
        "<b>D. Voltage doubler:</b> a Greinacher doubler needs two cells per position: 174 V DC at 85 % (one cell: "
        "151 V, 17 %, sparks over).",
        "<b>E. Varactor with memory:</b> the oxide film, set by its forming voltage (10 V to 140 V), tunes a 1 mH tank "
        "from 643 Hz to 2,405 Hz, matching 1/(2π√LC) to 1 Hz.",
        "<b>F. Frequency doubler:</b> 60 Hz in, 120 Hz out; 95 % of the output's AC power is at 120 Hz with cells, 89 % "
        "with water spark gaps.",
        "<b>No energy is added</b> by cells or gaps: every joule at a load comes from the mains or the surge source.",
    ])

    r.h1("How it is solved")
    r.p("Circuits with only cells run in kronspice (D, E and the twin of A). A spark gap is a switch with memory: it "
        "strikes at its breakdown voltage and stays on while its arc carries current, which kronspice's static elements "
        "cannot hold. Circuits with gaps therefore run in the example's own small nodal solver: backward Euler for "
        "capacitors and inductors, Newton on the cells' curves, gap states updated after each step. It is checked "
        "against an exact RC discharge and against kronspice on the supply. Its energy ledger splits what capacitors "
        "and inductors take in into stored energy and backward Euler's numerical damping, and closes to a millionth "
        "of the energy supplied or better.")
    r.p("Spark gaps (declared unless cited):")
    r.bullets([
        "Water polarity gap: 3 mm point–plane; breakdown with the point positive from Stygar's water relation "
        "(extrapolated to 60 Hz), negative ×3; arc 50 V + 0.5 Ω, out below 0.1 A after 10 µs.",
        "After an arc the breakdown recovers as 1 − e<super>−t/1 ms</super>, and never below the Paschen minimum of "
        "water vapour, 480 V (Škoro et al. 2011).",
        "Air gaps: Paschen's law for air (Raizer), arc 30 V + 0.2 Ω, never below air's Paschen minimum, about 330 V.",
    ])

    r.h1("Results")
    r.h2("A. Mains DC supply")
    r.p("Four lead–aluminium cells, 100 V peak at 60 Hz, 0.5 Ω source, 2,200 µF, 50 Ω load: <b>87.9 V DC</b>, ripple "
        "5.2 V, efficiency 85 %. kronspice gives 87.96 V (difference −0.01 %), which checks the example's own solver. "
        "Ledger residual 3.7×10<super>−9</super> J of 74.8 J supplied.")

    r.figure(os.path.join(HERE, "diodebridge.png"),
             "Figure 1. B (top left): largest voltage across any bridge cell during the 4 kV surge, unprotected (red, "
             "peak 221 V) and protected (blue, peak 117 V), against the 140 V oxide spark-over (dotted). C (top right): "
             "capacitor voltage of the spark pulser for the three output switches; only the zero-crossing-fired gap "
             "(green) charges and dumps cleanly every half cycle. E (bottom left): tank resonance against forming "
             "voltage, kronspice AC points on the 1/(2π√LC) curve, labelled with the film capacitance. F (bottom right): "
             "DC, 120 Hz and 240 Hz amplitudes divided by the source peak for the three bridges.",
             max_height_mm=125)

    r.h2("B. Surge-protected supply")
    r.p("A 4 kV, 1.2/50 µs surge from a 2 Ω source (the IEC 61000-4-5 test-wave shape and impedance; its 8/20 µs current "
        "shape is not modelled) arrives on the mains line at the mains peak. Protection has two stages: a 0.2 mm air gap "
        "across the line (breaks down at 1,423 V, Paschen's law) takes the bulk of the surge; after 10 µH and 0.5 Ω of "
        "line, an aluminium–aluminium water cell across the bridge input blocks both ways up to 140 V and then sparks "
        "over, a two-way clamp like a pair of Zener diodes.")
    r.table([
        ["", "Largest bridge-cell voltage", "Extra energy into bridge cells", "Arrester", "Clamp", "DC output peak", "Ledger residual"],
        ["Unprotected", "221 V (sparks over at 140 V)", "15.42 J", "0.0 J (0 strike)", "0.000 J", "131 V", "3×10<super>−10</super> J"],
        ["Protected", "117 V", "0.64 J", "22.9 J (1 strike)", "0.006 J", "95 V", "2×10<super>−10</super> J"],
    ], [0.13, 0.19, 0.17, 0.14, 0.09, 0.13, 0.15], bold_rows=(2,))

    r.h2("C. Spark pulser")
    r.p("A 6 kV, 60 Hz transformer (1 kΩ) feeds a bridge of four water polarity gaps (polarity ratio 3, the ratio that "
        "rectified in <i>sparkgap_diode</i>), which charges 0.1 µF. An air gap discharges the capacitor through 1 µH into "
        "a load.")
    r.table([
        ["Output switch and load", "Pulses per s", "Peak", "To load per pulse", "Load / mains", "Bridge backward sparks", "What happens"],
        ["Self-breaking 1 mm air gap, 1 Ω (rings)", "120", "1,200 A", "0.73 J", "88 W of 13,960 W (0.6 %)",
         "2,564 (15.26 J)", "capacitor swings negative, bridge gaps re-strike backwards: the transformer is shorted"],
        ["Self-breaking 1 mm air gap, 8 Ω (overdamped)", "120", "482 A", "2.11 J", "253 W of 14,081 W (1.8 %)",
         "0 (0.00 J)", "the mains keeps the output arc lit (power follow): gap conducting 0.065 s of 0.1 s"],
        ["2 mm air gap fired at each mains zero crossing, 8 Ω", "120", "505 A", "1.13 J", "135 W of 269 W (50.3 %)",
         "2,804 (0.02 J)", "the transformer cannot follow at its zero; the pulser works"],
    ], [0.2, 0.09, 0.08, 0.1, 0.15, 0.13, 0.25], bold_rows=(3,))
    r.p("Firing at a fixed point of the mains cycle is what Tesla-coil builders do with a synchronous rotary spark gap; "
        "here the trigger is ideal (declared). In the working design the bridge gaps still make tiny backward sparks, "
        "0.02 J in total, because the transformer and the full capacitor nearly balance. Ledger residuals: "
        "8×10<super>−7</super> J, 8×10<super>−7</super> J and 6×10<super>−6</super> J.")

    r.h2("D. Voltage doubler")
    r.p("Greinacher's circuit (a pump capacitor and two diodes; 2,200 µF, 200 Ω) with lead–aluminium cells, 100 V peak, "
        "solved in kronspice:")
    r.table([
        ["Cells per position", "Output", "Efficiency", "Reverse peak", "Ledger mismatch"],
        ["1", "151 V", "17 % (sparks over)", "160 V", "−0.0002"],
        ["2", "174 V", "85 %", "182 V", "−1×10<super>−5</super>"],
    ], [0.2, 0.15, 0.25, 0.2, 0.2], bold_rows=(2,))
    r.p("Each position sees twice the peak in reverse, so it needs two cells, as in the historical \"two or more valves in "
        "series for higher pressures\".")

    r.h2("E. A varactor with memory")
    r.p("The aluminium plate's oxide film is a capacitor whose thickness is set when it is formed: about 1.3 nm per volt "
        "(declared, typical barrier anodic alumina). Formed at a higher voltage, the film is thicker and its capacitance "
        "smaller. A 1 mH tank across a cell biased in its blocking direction, solved by kronspice AC:")
    r.table([
        ["Formed at", "Film", "Resonance (kronspice)", "1/(2π√LC)", "Q"],
        ["10 V", "61.30 µF", "643 Hz", "643 Hz", "24"],
        ["20 V", "30.65 µF", "909 Hz", "909 Hz", "33"],
        ["40 V", "15.32 µF", "1,285 Hz", "1,286 Hz", "48"],
        ["70 V", "8.76 µF", "1,701 Hz", "1,701 Hz", "62"],
        ["140 V", "4.38 µF", "2,405 Hz", "2,405 Hz", "87"],
    ], [0.18, 0.18, 0.24, 0.22, 0.18])
    r.p("Biased forward (−5 V) the conducting cell shorts the tank and there is no resonance. <b>This is not a "
        "semiconductor varactor:</b> its capacitance does not follow the bias and stays fixed below the forming voltage. "
        "It is set by the highest voltage the film has seen, a capacitor with memory. Re-forming higher thickens the film "
        "at once; thinning it again takes the slow dissolution of the oxide, which is not modelled. Q is low at low forming "
        "voltage because the film then leaks more (declared leakage).")

    r.h2("F. Frequency doubler")
    r.p("A full bridge's output repeats twice per cycle. Amplitudes divided by the source peak, resistive load:")
    r.table([
        ["Bridge", "DC", "2f (120 Hz)", "4f (240 Hz)", "Share of AC power at 2f"],
        ["Perfect diode bridge", "0.637", "0.424", "0.085", "95 %"],
        ["Lead–aluminium cell bridge, 100 V", "0.580", "0.411", "0.081", "95 %"],
        ["Water spark-gap bridge, 6 kV", "0.495", "0.534", "0.108", "89 %"],
    ], [0.34, 0.12, 0.16, 0.16, 0.22])
    r.p("A filter tuned to 120 Hz picks out the doubled frequency. The spark-gap bridge's sudden jumps put more into "
        "120 Hz and the higher harmonics than a smooth bridge.")

    r.h1("Ledger and validation")
    r.bullets([
        "Every run has an energy ledger; the gap solver's ledger closes to a millionth of the energy supplied or better "
        "(A: 3.7×10<super>−9</super> J of 74.8 J; B: 2–3×10<super>−10</super> J; C: ≤ 6×10<super>−6</super> J).",
        "Twin: the example's solver and kronspice agree on the DC supply to 0.01 %; the solver is also checked against "
        "an exact RC discharge.",
        "Varactor resonances from kronspice AC agree with 1/(2π√LC) to within 1 Hz.",
        "Tests in <i>main_test.go</i>: model match, solver checks, the twin, and each job's result.",
    ])

    r.h1("Limits of the model")
    r.bullets([
        "The water gap's breakdown extrapolates Stygar's pulsed-power relation to 60 Hz; the negative/positive polarity "
        "ratio, arc constants and 1 ms recovery are declared.",
        "The surge's 8/20 µs current shape is not modelled; the zero-crossing trigger in C is ideal.",
        "The oxide's growth rate (1.3 nm/V) and leakage are declared; slow dissolution of the film is not modelled.",
    ])

    r.h1("Reproducing the results")
    r.code("go run ./examples/sparkgap_diodebridge -viz examples/sparkgap_diodebridge/viz   # about 20 seconds\n"
           "python3 examples/sparkgap_diodebridge/viz/render.py\n"
           "go test ./examples/sparkgap_diodebridge")
    r.small("Sources: Encyclopaedia Britannica, 11th ed. (1911), vol. 27, \"Transformers\"; H. Greinacher, Z. Phys. 4, "
            "195 (1921); IEC 61000-4-5; Yu. P. Raizer, Gas Discharge Physics (1991); W. A. Stygar et al., Phys. Rev. ST "
            "Accel. Beams 9, 070401 (2006); Škoro et al., Phys. Rev. E 84, 055401 (2011).")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "sparkgap_diodebridge_report.pdf"))
