# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""PDF report for examples/water_lamp: the flashing neon lamp with water instead of neon.

Run from the repository root:

    python examples/water_lamp/viz/make_report.py examples/water_lamp/water_lamp_report.pdf

Numbers come from examples/water_lamp/run_output.txt (the program's output) and
examples/water_lamp/README.md (construction, declared values, sources); figures
are viz/lamps.png and viz/coil.png. Layout from internal/pdfreport/kronreport.py.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401


def build(out):
    r = Report("examples/water_lamp", "The flashing lamp with water",
               "Can water replace the neon and the radioactive starter of a relaxation-oscillator lamp, "
               "and how much of each flash can it pass into the Tesla bifilar coil? — CLAUDE.md §3g, §3j, H5",
               "CLAUDE.md §3g, §3j, register entry H5", "GA seed 1, 24 designs × 12 generations")

    r.h1("Summary")
    r.p("The classic flashing-lamp (relaxation) oscillator charges a capacitor C through a resistor R from a "
        "supply V<sub>s</sub> until the lamp strikes and empties it. A neon NE-2 lamp strikes reliably because "
        "a trace of radioactive krypton-85 supplies the first free electron; in the dark, without it, the lamp "
        "waits for a cosmic ray. This example replaces the neon and the starter with water in three ways "
        "(a low-pressure water-vapour glow tube, a spark gap in pure liquid water, an air spark gap with water "
        "mist), lets each lamp drive the first mode of the Tesla bifilar coil of <i>pulsed_bifilar</i>, tunes "
        "the flash circuit to the coil as a spark-gap Tesla coil is tuned, and tries gaps and capacitor plates "
        "in parallel. The misted gap is the §3j front end; its faster recovery is register hypothesis H5, "
        "carried here as an untested parameter.")
    r.bullets([
        "<b>The dark effect is the same for water vapour as for neon.</b> Without a starter, the water-vapour tube "
        "flashed 54 times in 36,000 s, every 662 s ± 108 % (expected wait 792 s); with an ultraviolet LED starter it "
        "flashes every 0.174 s ± 5 % at 492 ± 4 V.",
        "<b>The liquid-water gap strikes at 1,076 ± 91 V</b> from a 3,000 V supply, every 0.00442 s ± 11 %; its "
        "strike voltage is an extrapolation of a relation measured for sub-microsecond pulses.",
        "<b>Mist raises the air gap's strike voltage</b> from 6,798 ± 894 V (dry) to 8,224 ± 643 V and lengthens "
        "the flash interval from 0.0118 s to 0.0181 s. Faster recovery with mist (H5) is not modelled; it is untested.",
        "<b>Untuned, almost nothing reaches the coil:</b> 0.2 % of each liquid-water spark and 0.6 % of each mist "
        "spark; the glow lamps deliver about 10<super>−11</super> J.",
        "<b>Tuned and quenched, an evolved mist-gap circuit puts 74.7 % of each flash into the coil</b> "
        "(0.249 µF, 20 primary turns, quench at 0.65 of the transfer time; 4.35 J per flash, coil peak "
        "5.157×10<super>4</super> V, 25.5 % of the supply energy). The coil voltage is far above what enamelled "
        "wire holds.",
        "<b>Parallel gaps do not add spark current</b> (438 A for 1, 2 and 4 gaps, 436 A for 8); parallel capacitor "
        "plates do (438 A to 2,818 A for 1 to 8 plates), but the share reaching the coil falls from 74.7 % to 18.0 %.",
        "<b>Every ledger closes</b> (residual ≤ 10<super>−5</super> relative; the GA design to 2.9×10<super>−6</super>). "
        "Every joule comes from the supply; the water adds no energy.",
    ])

    r.h1("Background")
    r.p("A gas lamp strikes only when a free electron is present to start the avalanche. In the model free "
        "electrons appear at a declared rate S, and once V is past breakdown V<sub>b</sub> each one starts a "
        "discharge with probability 1 − V<sub>b</sub>/V (declared form), so the mean wait at the supply voltage is "
        "1/(S·(1 − V<sub>b</sub>/V<sub>s</sub>)). Cosmic rays give about 2 ion pairs per cm³ per second in air at "
        "1 atm (declared), scaled with gas density. With no statistical wait the period is the textbook "
        "R·C·ln((V<sub>s</sub> − V<sub>m</sub>)/(V<sub>s</sub> − V<sub>b</sub>)) plus the short on-time, where "
        "V<sub>m</sub> is the maintaining voltage. All randomness comes from one seeded generator, so runs repeat exactly.")
    r.p("<b>H5 (hypothesis):</b> water mist speeds spark-gap recovery and sharpens the pulse train. "
        "<b>Standard prediction used here:</b> the misted gap recovers as the dry gap (1 ms declared); only the "
        "fog's measured rise in breakdown voltage (Liu 2018) is applied. No measurement of faster recovery "
        "with mist has been found, so H5 stays untested.")

    r.h1("The model")
    r.p("The lamps, with cited or declared values (from the README):")
    r.table([
        ["Lamp", "Supply, R, C", "Strike", "Maintain", "Starting electrons"],
        ["Neon NE-2, radioactive starter", "120 V, 10 MΩ, 10 nF", "76 V (TekWiki)", "58 V", "10<super>4</super>/s (declared)"],
        ["Neon NE-2, no starter, dark", "same", "76 V", "58 V", "cosmic rays in ~0.03 cm³ at ~60 Torr: 0.0047/s"],
        ["Water-vapour tube, dark", "600 V, 22 MΩ, 10 nF",
         "480 V, Paschen minimum of water vapour at 0.6 Torr·cm (Škoro 2011)", "0.76 × strike (declared)",
         "cosmic rays in 2 cm³ at 1.2 Torr: 0.0063/s"],
        ["Water-vapour tube, UV LED starter", "same", "480 V", "365 V", "10³/s photoelectrons (declared)"],
        ["Spark gap in pure liquid water", "3,000 V, 1 MΩ, 10 nF",
         "E·τ<super>1/3</super> = 0.135 MV/cm·µs<super>1/3</super> ± 12 % over 1 mm (Stygar 2006)", "50 V arc", "none needed"],
        ["Air spark gap, dry", "10,000 V, 1 MΩ, 10 nF", "Paschen, 1 mm at 1 atm: 5,035 V (Raizer)", "30 V arc",
         "10³/s room-air ions; 1 ms recovery (declared)"],
        ["Air spark gap with water mist", "same", "dry × 29.4/21.7 = 6,822 V (Liu 2018)", "30 V", "as dry; recovery as dry"],
    ], [0.22, 0.17, 0.27, 0.13, 0.21])
    r.bullets([
        "Paschen's law for air: V = B·pd/ln(A·pd/ln(1 + 1/γ)), A = 15 /(cm·Torr), B = 365 V/(cm·Torr) "
        "(Raizer, Table 4.1 and eq. 7.3), γ = 0.01 declared. For water vapour the curve V = B·pd/ln(pd/x<sub>0</sub>) is "
        "pinned to the measured minimum, 480 V at 0.6 Torr·cm.",
        "The liquid gap's breakdown relation was measured for pulses of 0.011 to 0.6 µs; the millisecond charging here "
        "is far outside that range, so its strike voltage is an <b>extrapolation</b>. Pure water conducts "
        "(0.055 µS/cm), so a leak current splits water: the ledger books 1.229 V × charge on the chemical channel "
        "and the rest as heat.",
        "The fog factor is Liu's 5 cm rod–plane result applied to 1 mm, also an <b>extrapolation</b>.",
        "One RK4 integrator carries the state [V, I<sub>primary</sub>, I<sub>coil</sub>, V<sub>coil</sub>, supply "
        "energy, R heat, leak, lamp energy, coil heat]; the energies are states, so every ledger closes to RK4 accuracy.",
        "The coil's first mode: L<sub>s</sub> 4.853 mH, f<sub>1</sub> 41.61 kHz, Q 83 (from <i>pulsed_bifilar</i>), "
        "5-turn primary L<sub>p</sub> 5.97 µH, mutual 28.1 µH (k = 0.165). A one-mode model, checked against a "
        "kronspice transient of the same network.",
    ])

    r.h1("Results")
    r.h2("1. Seven lamps in the same circuit")
    r.table([
        ["Lamp", "Supply", "Simulated", "Flashes", "Strike [V]", "Flash interval", "Energy per flash"],
        ["Neon NE-2, radioactive starter", "120 V", "2 s", "49", "78 ± 1", "0.039 s ± 6 %", "1.35×10<super>−5</super> J"],
        ["Neon NE-2, no starter, dark", "120 V", "36,000 s", "67", "120 ± 0", "540 s ± 100 % (expected 576 s)", "5.51×10<super>−5</super> J"],
        ["Water-vapour tube, dark", "600 V", "36,000 s", "54", "600 ± 0", "662 s ± 108 % (expected 792 s)", "0.00113 J"],
        ["Water-vapour tube, UV LED starter", "600 V", "2 s", "10", "492 ± 4", "0.174 s ± 5 %", "0.000545 J"],
        ["Spark gap in liquid water (pure)", "3,000 V", "2 s", "452", "1,076 ± 91", "0.00442 s ± 11 %", "0.00582 J"],
        ["Air spark gap, dry", "10,000 V", "2 s", "168", "6,798 ± 894", "0.0118 s ± 27 %", "0.235 J"],
        ["Air spark gap with water mist", "10,000 V", "2 s", "110", "8,224 ± 643", "0.0181 s ± 24 %", "0.34 J"],
    ], [0.23, 0.09, 0.11, 0.09, 0.11, 0.21, 0.16])
    r.figure(os.path.join(HERE, "lamps.png"),
             "Figure 1. Capacitor voltage of the flashing-lamp circuit; each drop is one flash. Left column: neon with its "
             "radioactive starter (top), the water-vapour tube with a UV LED starter (middle) and the spark gap in liquid "
             "water (bottom, scatter from the ±12 % breakdown relation). Right column: neon and the water-vapour tube "
             "without a starter in the dark, each vertical line one flash over 10 hours (random cosmic-ray waits), and the "
             "air gap with water mist (bottom).", max_height_mm=150)
    r.bullets([
        "The leak through the liquid-water gap splits 0.000166 J worth of water into hydrogen and oxygen in 2 s "
        "(7×10<super>−10</super> mol H<sub>2</sub>); the rest of the 0.1025 J leakage is heat.",
        "Example ledger (liquid-water gap): supply 14.33 J = heat in R 11.6 + lamp 2.631 + leakage 0.1025 + stored "
        "0.0007132 J, residual −3.5×10<super>−11</super> J. The cosmic-ray runs close to ≤ 1.0×10<super>−5</super> J.",
    ])

    r.h2("2. Driving the bifilar coil (untuned)")
    r.p("The reduced one-mode coil model agrees with the kronspice transient for the same pulse to 0.0008 of the peak.")
    r.table([
        ["Lamp", "Flashes", "Coil peak [V]", "Energy to coil", "Of flash energy"],
        ["Neon NE-2, radioactive starter", "11", "0.03036", "1.35×10<super>−11</super> J", "0.0 %"],
        ["Water-vapour tube, UV LED starter", "7", "0.04821", "2.28×10<super>−11</super> J", "0.0 %"],
        ["Spark gap in liquid water (pure)", "41", "43.79", "9.56×10<super>−5</super> J", "0.2 %"],
        ["Air spark gap with water mist", "14", "312.9", "0.00189 J", "0.6 %"],
    ], [0.36, 0.1, 0.16, 0.2, 0.18])
    r.p("Glow discharges are too slow and gentle. Sparks ring the coil, but the 10 nF flash capacitor rings with the "
        "6 µH primary at about 650 kHz, 16 times the coil's 41.6 kHz, so almost nothing reaches it.")

    r.figure(os.path.join(HERE, "coil.png"),
             "Figure 2. Coil voltage. Top: the liquid-water spark gap rings the coil to about 40 V per flash (untuned). "
             "Middle: the misted air gap, about 300 V per flash (untuned). Bottom: one flash of the evolved, tuned and "
             "quenched mist-gap circuit, ringing the coil to about 50 kV and decaying with the coil's Q.",
             max_height_mm=120)

    r.h2("3. Tuning the flash circuit to the coil")
    r.p("Tuning sets L<sub>p</sub>·C<sub>p</sub> = L<sub>s</sub>·C<sub>s</sub>, as in a spark-gap Tesla coil: with the "
        "5-turn primary C<sub>p</sub> = 2.45 µF, and energy swings from primary to coil after half a beat, "
        "1/(2k·f<sub>1</sub>) = 72.8 µs. An ordinary arc keeps conducting through current zeros until the current has "
        "stayed small for its de-ionization time (20 µs declared), so energy swings back and burns in the arc; a "
        "quenched gap opens at the first current zero after the transfer and leaves the energy in the coil.")
    r.table([
        ["Gap and circuit", "Flashes", "Spark [A]", "To coil per flash", "Of flash", "Of supply", "Coil peak"],
        ["Liquid water, untuned (10 nF)", "41", "47.62", "2.33×10<super>−6</super> J", "0.2 %", "0.00 %", "43.79 V"],
        ["Liquid water, tuned, ordinary arc", "47", "173.4", "0.00622 J", "4.4 %", "0.23 %", "2,401 V"],
        ["Liquid water, tuned, quenched gap", "47", "174.6", "0.00622 J", "4.4 %", "0.23 %", "2,421 V"],
        ["Mist, untuned (10 nF)", "14", "339.8", "0.000135 J", "0.6 %", "0.01 %", "312.9 V"],
        ["Mist, tuned, ordinary arc", "5", "3,490", "1.9 J", "3.3 %", "1.08 %", "6.853×10<super>4</super> V"],
        ["Mist, tuned, quenched gap", "5", "3,486", "5.41 J", "9.5 %", "3.16 %", "6.844×10<super>4</super> V"],
        ["Mist, evolved: 0.249 µF, 20 turns, quench at 0.65 of transfer", "", "438", "4.35 J", "74.7 %", "25.5 %",
         "5.157×10<super>4</super> V"],
    ], [0.3, 0.09, 0.09, 0.14, 0.09, 0.1, 0.19], bold_rows=(7,))
    r.p("The genetic algorithm (<i>evo.Run</i>, 24 designs × 12 generations, seed 1) chose the flash capacitor "
        "(1 nF to 20 µF), the primary turns (1 to 20) and the quench time (0.2 to 3 beat halves) for the mist gap, "
        "minimizing 1 − share with at least three flashes in six charging time constants. Its ledger closes to "
        "2.9×10<super>−6</super>. 20 turns is the top of the allowed range. Charging a capacitor through a resistor from "
        "a fixed supply loses half the supply energy in the resistor, so \"of supply\" is at most 50 %.")
    r.p("<b>Caution:</b> in a bifilar coil neighbouring turns carry about half the coil voltage, here about 26 kV, far "
        "above what enamelled wire holds (a few kV, declared). The design would need much thicker insulation or a "
        "smaller drive (§3e rule).")
    r.h2("4. Plates in parallel")
    r.p("(a) N mist gaps in parallel across the evolved capacitor; only the first to break down carries the flash.")
    r.table([
        ["Gaps", "Strike [V]", "Flash interval", "Spark [A]", "To coil per flash", "Of flash"],
        ["1", "6,914 ± 29", "0.532 s ± 1 %", "438", "4.35 J", "74.7 %"],
        ["2", "5,922 ± 2,394", "0.468 s ± 60 %", "438", "3.23 J", "71.5 %"],
        ["4", "4,254 ± 2,726", "0.24 s ± 140 %", "438", "1.34 J", "59.7 %"],
        ["8", "3,602 ± 2,452", "0.15 s ± 192 %", "436", "0.677 J", "45.3 %"],
    ], [0.1, 0.2, 0.2, 0.14, 0.2, 0.16])
    r.p("More gaps give more starting electrons, so a gap re-strikes before the others have recovered: smaller, more "
        "frequent flashes, and the spark current does not add.")
    r.p("(b) n capacitor sections in parallel, primary turns re-chosen to stay tuned.")
    r.table([
        ["Plates", "C [µF]", "Turns", "Detuning", "Spark [A]", "Flash energy", "To coil per flash", "Of flash"],
        ["1", "0.249", "20", "−0.5 %", "438", "5.82 J", "4.35 J", "74.7 %"],
        ["2", "0.498", "13", "−1.0 %", "846", "11.5 J", "6.66 J", "57.8 %"],
        ["4", "0.997", "8", "4.7 %", "1,695", "23.6 J", "7.85 J", "33.3 %"],
        ["8", "1.99", "6", "−5.4 %", "2,818", "46.2 J", "8.3 J", "18.0 %"],
    ], [0.1, 0.1, 0.09, 0.12, 0.12, 0.15, 0.18, 0.14])
    r.p("More plates raise the energy per flash (×n) and the spark current nearly as much; with fewer primary turns the "
        "coupling falls, so the share reaching the coil drops. A bigger bank needs the coil and primary redesigned together.")

    r.h1("Ledger and validation")
    r.bullets([
        "Every ledger: supply = heat in R + lamp (light and heat) + leakage + coil + stored; residual ≤ 10<super>−5</super> "
        "relative (≤ 10<super>−10</super> J where no cosmic-ray run is involved).",
        "Reduced coil model vs kronspice transient: largest difference 0.0008 of the peak.",
        "Tests: <i>TestPaschen</i> (water vapour's 480 V minimum at 0.6 Torr·cm; air at 1 atm, 1 mm near 4.5 to 5 kV), "
        "<i>TestRelaxationPeriod</i> (primed lamp follows R·C·ln(...); ledgers close), <i>TestDarkWait</i> (mean dark wait "
        "equals 1/(S·(1 − V<sub>b</sub>/V<sub>s</sub>))), <i>TestCoil</i>, <i>TestTuning</i>, <i>TestParallel</i>.",
    ])

    r.h1("Limits of the model")
    r.bullets([
        "The liquid-water strike voltage extrapolates Stygar's relation from ≤ 0.6 µs pulses to millisecond charging.",
        "The fog factor extrapolates Liu's 5 cm rod–plane result to a 1 mm gap; recovery with mist is set equal to dry.",
        "Starting-electron rates for the radioactive starter, UV LED and room-air ions, the maintaining-voltage ratio of "
        "the water-vapour tube, γ, the arc de-ionization time and the insulation limit are declared, not measured.",
        "The coil is reduced to its first mode; the quenched gap is named in the code without a citation.",
    ])

    r.h1("Status")
    r.p("H5 (water mist speeds spark-gap recovery) remains <b>untested</b>: here the mist only raises breakdown, as "
        "measured, and recovery is a parameter set equal to the dry gap. No energy enters from the water.")

    r.h1("Reproducing the results")
    r.code("go run ./examples/water_lamp -viz examples/water_lamp/viz   # about 3 minutes\n"
           "python3 examples/water_lamp/viz/render.py\n"
           "go test ./examples/water_lamp\n"
           "KRONSPICE_EXPORT=1 go test ./examples/water_lamp -run TestExport")
    r.small("Sources: Škoro et al., Phys. Rev. E 84, 055401(R) (2011); TekWiki, \"Neon Lamps\"; Stygar et al., "
            "Phys. Rev. ST Accel. Beams 9, 070401 (2006); Raizer, Gas Discharge Physics (1991); Liu et al., Energies 11, "
            "2251 (2018); Bhoomireddy, PhD thesis, HBNI/BARC (2016); Terman, Radio Engineers' Handbook (1943); CODATA 2018.")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "water_lamp_report.pdf"))
