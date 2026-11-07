# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""pulsed_bifilar report: pulse-driven resonance in ordinary and Tesla bifilar coils (CLAUDE.md §3e).

python examples/pulsed_bifilar/viz/make_report.py examples/pulsed_bifilar/pulsed_bifilar_report.pdf

Numbers come from run_output.txt (coil sweep and ledgers), pickup_output.txt
(the tuned pickup) and README.md (construction); figures are viz/*.png made
by viz/render.py. Layout from internal/pdfreport/kronreport.py.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401
from reportlab.platypus import KeepTogether  # noqa: E402


def build(out):
    r = Report("examples/pulsed_bifilar", "Pulse-driven ordinary and bifilar coils",
               "When do one-sign pulses build a coil into a large standing wave? — CLAUDE.md §3e step 3",
               "CLAUDE.md §3e", "100 V, 50 Ω one-sign pulses; 20 sections per coil")

    r.h1("Summary")
    r.p("Four coils, each wound from about 126 m of the same 0.5 mm copper wire, are driven by one-sign "
        "100 V pulses through a 5-turn primary: an ordinary solenoid, a Tesla bifilar solenoid (US 512,340), "
        "a bifilar cone and a bifilar flat spiral. Each coil is a network of 20 magnetically coupled sections, "
        "not a lumped L and C. The example asks at which pulse rates the energy builds up into a large "
        "standing wave, and then takes power off the ordinary coil with a low-voltage, high-current pickup.")
    r.bullets([
        "<b>Build-up only in step with resonance.</b> A one-sign train at rate r contains r, 2r, 3r, …, so a coil "
        "builds up only when one of these lands on its first resonance f<sub>1</sub>: at r = f<sub>1</sub>, "
        "f<sub>1</sub>/2, f<sub>1</sub>/3. At 0.9, 1.1, 1.5 or 2 × f<sub>1</sub> none does.",
        "<b>Ordinary solenoid:</b> f<sub>1</sub> = 1.167 MHz (wire = 0.489 wavelengths, a free half-wave line), Q ≈ 210; "
        "pulses at f<sub>1</sub> give 7,426.4 V between the coil's ends from 100 V pulses (×74).",
        "<b>Bifilar windings resonate 28 to 33 times lower</b> (34.93–41.61 kHz) because of the 12 nF between their "
        "strands, with Q ≈ 83–96, and reach 154.6–281.6 V at f<sub>1</sub>.",
        "<b>V/2 between neighbouring turns.</b> At r = f<sub>1</sub> neighbouring turns differ by 141.3 V in the bifilar "
        "solenoid, 77.9 V in the cone and 127.0 V in the spiral (half the coil voltage), against 37.9 V in the ordinary solenoid.",
        "<b>Every joule comes from the pulse source.</b> At f<sub>1</sub> the ledger closes to 1.6×10<super>−5</super> J "
        "(ordinary) and 10<super>−6</super>–10<super>−7</super> J (bifilar).",
        "<b>The pickup converts voltage to current; it adds no energy.</b> A tuned 4-turn pickup at the coil's middle "
        "delivers 3.828 A at 0.766 V (1.4010 W), 66.5 % of the power that reaches the primary past the source's 50 Ω; "
        "the best variant reaches 91.3 %.",
    ])

    r.h1("What is being tested")
    r.p("CLAUDE.md §3e step 3 states that, for pulse drive, a coil's reflected wave adds to the drive when the "
        "pulse repeats at the round-trip time of the winding (or a whole multiple of it). Tesla's US 512,340 "
        "bifilar winding joins the end of one strand to the start of the other, so neighbouring turns sit half "
        "a wire length apart in potential (the patent's V/2). This example checks both statements on full "
        "distributed models and books the energy of every run. It is standard circuit physics; no hypothesis "
        "is involved. The bifilar solenoid's first mode is also the load used by <i>water_lamp</i>, and the same "
        "coil formulas serve <i>coiled_tether</i> and <i>bial_bifilar</i>.")

    r.h1("The model")
    r.h2("Wire and geometry")
    r.p("0.5 mm copper at 0.7 mm between wire centres (0.2 mm of enamel, ε = 3 declared). Each winding holds about "
        "the same length of wire as the 400-turn solenoid:")
    r.table([
        ["Coil", "Turns", "Radius", "Height", "Wire", "Surroundings"],
        ["Ordinary cylindrical solenoid", "400", "5.0 cm", "27.9 cm", "125.7 m", "grounded can of twice its radius"],
        ["Tesla bifilar solenoid (US 512,340)", "400", "5.0 cm", "27.9 cm", "125.7 m", "same can"],
        ["Tesla bifilar cone, 45° opening upward", "248", "2.0 to 14.2 cm", "12.2 cm", "126.4 m", "grounded plane 5 cm below the base"],
        ["Tesla bifilar flat (Archimedean) spiral", "214", "2.0 to 16.9 cm", "0", "127.1 m", "grounded plane 5 cm below"],
    ], [0.30, 0.08, 0.14, 0.10, 0.10, 0.28])
    r.p("In a bifilar coil the wire runs out along the even positions (strand 1) and back along the odd ones "
        "(strand 2), the end of the first joined to the start of the second, so neighbouring turns belong to the two strands.")
    r.h2("Distributed network")
    r.bullets([
        "<b>Inductance:</b> 20 sections per coil. The section inductance matrix is summed turn by turn from Maxwell's "
        "formula for coaxial circles (<i>Treatise</i> art. 701, complete elliptic integrals), with a single turn's "
        "high-frequency self-inductance μ<sub>0</sub>R(ln(8R/r<sub>w</sub>) − 2) (Grover 1946). Sections are coupled in "
        "kronspice with K elements.",
        "<b>Capacitance:</b> between adjacent turns, the two-wire formula π·ε<sub>0</sub>ε<sub>r</sub>·πD/acosh(p/d) with each "
        "pair's mean diameter D (12.1 nF summed over all neighbouring turns); to the surroundings, the winding as a "
        "sheet facing the nearest grounded surface, split half to each end node of a section.",
        "<b>Losses:</b> Terman's skin-effect resistance at f<sub>1</sub>; a 1 GΩ bleeder gives the floating coil a DC reference.",
        "<b>Drive:</b> 100 V, 50 Ω source into a 5-turn primary of 1 mm wire at 2 mm pitch just below the base, coupled to "
        "every section by the same coaxial-circle mutuals. Pulse edges are 2 % and width 10 % of 1/f<sub>1</sub>; rates "
        "r/f<sub>1</sub> = 1/3, 1/2, 0.9, 0.97, 1, 1.03, 1.1, 1.5 and 2; each run lasts min(600, max(60, 5Q/π + 20)) "
        "cycles of f<sub>1</sub> (about five ring-up time constants) at 50 steps per period.",
        "<b>f<sub>1</sub> and Q</b> come from an AC sweep of the terminal voltage per volt of drive, 1 kHz to 100 MHz.",
        "<b>Declared approximations:</b> sections are equipotential inside; the bifilar junction lead and radiation are "
        "left out; proximity effect is not in the wire resistance; absolute tolerances are loosened to the circuit's "
        "scale (1 µA, 1 mV, 1 pC).",
    ])

    r.h1("Results")
    r.h2("The coils")
    r.table([
        ["Coil", "f<sub>1</sub>", "Wire / λ", "AC gain at top", "Q", "Inductance", "Pulse width", "Cycles"],
        ["Ordinary solenoid", "1.167 MHz", "0.489", "×338.6", "210", "4.853 mH (Nagaoka 4.873 mH)", "85.7 ns", "354"],
        ["Bifilar solenoid", "41.61 kHz", "0.017", "×12.1", "83", "4.853 mH", "2.40 µs", "151"],
        ["Bifilar cone", "35.81 kHz", "0.015", "×6.6", "96", "6.542 mH", "2.79 µs", "172"],
        ["Bifilar flat spiral", "34.93 kHz", "0.015", "×10.8", "94", "6.842 mH", "2.86 µs", "169"],
    ], [0.18, 0.11, 0.08, 0.11, 0.06, 0.24, 0.12, 0.08])
    r.p("The two solenoids have the same inductance to 10<super>−12</super> (the strands of a series bifilar coil carry "
        "the same current the same way), within 0.4 % of Nagaoka's current-sheet formula; the flat spiral is within "
        "1.7 % of Wheeler's flat-spiral formula. The ordinary solenoid's wire is 0.489 wavelengths at f<sub>1</sub>: it rings "
        "as a free half-wave line, the first resonance Knight reports for a free single-layer coil, with voltage maxima "
        "at its ends and the current maximum at its middle.")

    r.h2("Pulse-rate sweep")
    r.p("Voltage between the coil's ends at the end of each run (100 V pulses from 50 Ω):")
    r.table([
        ["Coil", "f<sub>1</sub>/3", "f<sub>1</sub>/2", "0.9 f<sub>1</sub>", "0.97 f<sub>1</sub>", "f<sub>1</sub>",
         "1.03 f<sub>1</sub>", "1.1 f<sub>1</sub>", "1.5 f<sub>1</sub>", "2 f<sub>1</sub>"],
        ["Ordinary solenoid", "2631.5 V", "3868.1 V", "454.8 V", "892.8 V", "<b>7426.4 V</b>", "657.4 V", "213.9 V", "479.5 V", "356.4 V"],
        ["Bifilar solenoid", "97.4 V", "143.4 V", "15.9 V", "54.9 V", "<b>281.6 V</b>", "58.5 V", "19.5 V", "6.6 V", "5.1 V"],
        ["Bifilar cone", "53.4 V", "78.6 V", "7.9 V", "26.5 V", "<b>154.6 V</b>", "28.4 V", "9.9 V", "3.9 V", "3.1 V"],
        ["Bifilar flat spiral", "86.8 V", "127.9 V", "12.8 V", "43.6 V", "<b>252.1 V</b>", "46.3 V", "15.5 V", "5.7 V", "4.7 V"],
    ], [0.19] + [0.09] * 9)
    r.p("Voltage between neighbouring turns at r = f<sub>1</sub>: 37.9 V (ordinary), 141.3 V (bifilar solenoid), "
        "77.9 V (cone), 127.0 V (spiral). Every bifilar winding has half its coil voltage between neighbouring turns, "
        "as the patent says.")
    r.figure(os.path.join(HERE, "pulse_sweep.png"),
             "Figure 1. Voltage between the coil's ends against pulse rate ÷ f<sub>1</sub> (log scale), for the four coils. "
             "Peaks appear only at f<sub>1</sub>/3, f<sub>1</sub>/2 and f<sub>1</sub>, where one harmonic of the one-sign "
             "train lands on the resonance.", 0.9)
    r.figure(os.path.join(HERE, "ring_up.png"),
             "Figure 2. Coil-end voltage envelope building up pulse by pulse with pulses at f<sub>1</sub>. The ordinary "
             "solenoid settles near 7.4 kV within about 300 µs; the bifilar coils settle at 150–280 V over a few "
             "milliseconds.", 0.9)
    r.figure(os.path.join(HERE, "coil_voltage.png"),
             "Figure 3. The two solenoids at f<sub>1</sub>, at the moment of largest voltage. Top: every turn coloured by its "
             "voltage. Bottom: voltage along the wire (blue) and between neighbouring turns (red). The ordinary coil "
             "shows a half-wave pattern with a small turn-to-turn voltage (32 V at this instant); in the bifilar coil every "
             "neighbouring pair differs by about half the coil voltage (141 V).", 0.95)
    r.figure(os.path.join(HERE, "coil_voltage_shapes.png"),
             "Figure 4. The bifilar cone and flat spiral at f<sub>1</sub>, drawn as in Figure 3. Neighbouring turns differ "
             "by 77 V (cone) and 127 V (spiral), about half the end-to-end voltage.", 0.95)

    r.h2("Low-voltage, high-current pickup")
    r.p("A 4-turn pickup of 2 mm copper (radius 6.5 cm, inside the can; L = 4.082 µH, 0.072 Ω at f<sub>1</sub>) is wound "
        "round the ordinary coil's middle, where the current is largest, and tuned with a series capacitor "
        "C = 1/((2πf<sub>1</sub>)²L) = 4.56 nF. The load (and, where stated, the capacitor) is chosen from a grid for the "
        "largest load power per V² of drive; the pulse rate is the coupled system's best frequency, searched from "
        "0.7 to 1.4 f<sub>1</sub>. \"Past 50 Ω\" is the power that reaches the primary after the source's own resistance; "
        "the share is load power divided by it.")
    r.table([
        ["Pickup", "C (nF)", "Load (Ω)", "Rate / f<sub>1</sub>", "I load (A)", "V load (V)", "P load (W)", "Past 50 Ω (W)", "Share", "Coil (V)"],
        ["Middle, tuned 1/(ω<sub>1</sub>²L)", "4.56", "0.2", "1.2975", "3.828", "0.766", "1.4010", "2.1055", "66.5 %", "3651"],
        ["Middle, capacitor optimized", "1.14", "10", "0.9717", "0.667", "6.668", "1.1996", "1.3138", "91.3 %", "3442"],
        ["Middle, untuned", "none", "1", "1.1073", "1.702", "1.702", "1.3051", "1.6369", "79.7 %", "3502"],
        ["Top (high-voltage end), tuned", "4.56", "0.5", "1.0487", "2.038", "1.019", "1.0064", "1.1590", "86.8 %", "2230"],
    ], [0.20, 0.08, 0.07, 0.08, 0.08, 0.08, 0.09, 0.10, 0.08, 0.08])
    r.bullets([
        "The pickup's load never gets more than the power that reaches the primary past the source's resistance; the "
        "pickup converts voltage into current and does not add energy.",
        "A tightly coupled tuned pickup splits the coil's resonance in two (about 0.83 and 1.30 f<sub>1</sub>), so the best "
        "pulse rate moves with it (1.2975 f<sub>1</sub> for the tuned middle pickup).",
    ])

    r.h1("Energy ledger")
    r.p("Energy from the pulse source over each run at r = f<sub>1</sub>, split into heat and energy stored at the end:")
    r.table([
        ["Coil, at f<sub>1</sub>", "In", "Heat", "Stored", "Time-stepping error"],
        ["Ordinary solenoid", "0.0023518 J", "0.0022037 J", "0.00013172 J", "1.6×10<super>−5</super> J"],
        ["Bifilar solenoid", "0.079822 J", "0.079705 J", "0.0001151 J", "1.2×10<super>−6</super> J"],
        ["Bifilar cone", "0.1083 J", "0.10827 J", "3.4799×10<super>−5</super> J", "1.4×10<super>−7</super> J"],
        ["Bifilar flat spiral", "0.10853 J", "0.10844 J", "9.2976×10<super>−5</super> J", "9.9×10<super>−8</super> J"],
    ], [0.26, 0.17, 0.17, 0.2, 0.2])
    r.p("The bifilar coils take far more energy from the source (their wide pulses at low f<sub>1</sub> spend most of it "
        "as heat in the primary circuit) and store less: their low-frequency, low-Q mode builds up less than the "
        "ordinary coil's half-wave resonance. Pickup ledgers:")
    r.table([
        ["Pickup", "In", "Heat (of it, load)", "Stored", "Time-stepping error"],
        ["Middle, tuned", "0.003514 J", "0.003469 J (0.0003756 J)", "4.23×10<super>−5</super> J", "2.7×10<super>−6</super> J"],
        ["Middle, capacitor optimized", "0.002382 J", "0.002349 J (0.0003185 J)", "2.72×10<super>−5</super> J", "3.4×10<super>−6</super> J"],
        ["Middle, untuned", "0.00323 J", "0.003187 J (0.0003571 J)", "2.85×10<super>−5</super> J", "1.5×10<super>−5</super> J"],
        ["Top, tuned", "0.002516 J", "0.002493 J (0.0002707 J)", "2.21×10<super>−5</super> J", "6.5×10<super>−7</super> J"],
    ], [0.26, 0.15, 0.27, 0.15, 0.17])
    r.h2("Tests")
    r.table([
        ["Test", "Checks"],
        ["TestCoilNetwork", "every winding uses each position once; L within 1 % of Nagaoka and equal for both solenoids; "
         "symmetric capacitance matrix; the bifilar coil holds (N − 1) strand-to-strand gaps"],
        ["TestFlatSpiralInductance", "the spiral against Wheeler's a²N²/(8a + 11c), held to 10 %"],
        ["TestPulseBuildUp", "pulses at f<sub>1</sub> build more than three times higher than at 1.5 f<sub>1</sub>; ledger closes to 1 %"],
        ["TestPickup", "the tuning capacitor resonates the pickup at f<sub>1</sub> to 10<super>−12</super>; ledger to 1 %; "
         "load power below the power past the source"],
        ["TestExport", "the four coils swept around resonance, and the tuned pickup, as netlists and results sheets"],
    ], [0.28, 0.72])

    n0 = len(r.s)
    r.h1("Limits of the model")
    r.bullets([
        "Twenty sections per coil, each equipotential inside; finer structure along the wire is not resolved.",
        "The bifilar junction lead and radiation are not modelled; proximity effect is not in the wire resistance.",
        "Turn-to-turn capacitance uses the two-wire formula; capacitance to the surroundings treats the winding as a sheet.",
        "Absolute tolerances are loosened to the circuit's scale (1 µA, 1 mV, 1 pC).",
    ])
    r.s[n0:] = [KeepTogether(r.s[n0:])]

    r.h1("Reproducing the results")
    r.code("go run ./examples/pulsed_bifilar                    # about 10 minutes, writes viz/*.csv\n"
           "go run ./examples/pulsed_bifilar -only-pickup       # the pickup study alone\n"
           "python3 examples/pulsed_bifilar/viz/render.py       # pictures and animations\n"
           "go test ./examples/pulsed_bifilar\n"
           "KRONSPICE_EXPORT=1 go test ./examples/pulsed_bifilar -run TestExport")
    r.small("Sources: N. Tesla, US 512,340 (1894) and US 593,138 (1897); J. C. Maxwell, <i>Treatise</i> Vol. II art. 701; "
            "F. W. Grover, <i>Inductance Calculations</i> (1946); H. Nagaoka, <i>J. Coll. Sci. Tokyo</i> 27, 1 (1909); "
            "H. A. Wheeler, <i>Proc. IRE</i> 16, 1398 (1928); Ramo, Whinnery &amp; Van Duzer, ch. 4; F. E. Terman, "
            "<i>Radio Engineers' Handbook</i> (1943); D. W. Knight, self-resonance of solenoid coils (2016), §13.")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "pulsed_bifilar_report.pdf"))
