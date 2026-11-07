# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""PDF report for examples/coiled_tether (a five-mile tether, coiled and bifilar, in orbit).

python examples/coiled_tether/viz/make_report.py examples/coiled_tether/coiled_tether_report.pdf

Numbers come from examples/coiled_tether/run_output.txt, with the declared
inputs, method and tests from examples/coiled_tether/README.md. The figure is
examples/coiled_tether/ring_up.png.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401


def build(out):
    r = Report("examples/coiled_tether", "A five-mile tether, coiled and bifilar, in orbit",
               "Does coiling a tether, or winding it as Tesla's bifilar coil, make it an energy source? "
               "— CLAUDE.md §3e, §3f", "CLAUDE.md §3e, §3f",
               "orbit 7.7 km/s, B = 30 µT, 816 W solar array (declared)")

    r.h1("Summary")
    r.p("A 5-mile (8,046.72 m) copper tether of 2 mm wire (43.0 Ω) in low Earth orbit is compared in three windings: "
        "straight (the ordinary electrodynamic tether), an ordinary single-layer coil 2 m across with 1,281 turns, and the "
        "same coil wound as Tesla's bifilar coil (US 512,340), in which each turn lies next to the turn half a coil away. "
        "The question is whether \"multiplication between windings\" makes the coil an energy source. Energy comes only "
        "from declared sources: the Sun, and Earth's magnetic and electric fields.")
    r.bullets([
        "<b>Coiling throws the field's voltage away:</b> the straight tether gives 1858.792 V, the ordinary coil 0.888 V and "
        "the bifilar coil 0.889 V, which is E′ = v×B = 0.231 V/m times the coil's own 3.84 m length.",
        "<b>The Sun can drive the coil at resonance:</b> 816 W at 565.7 Hz (Q 112.5) ring 265.1 V up to 29,784 V peak, ×112.",
        "<b>Bifilar moves the energy between the turns, not into the books:</b> 25.78 J of the 25.783 J stored sits between "
        "the bifilar turns against 6.283×10<super>−5</super> J for the ordinary coil, the same total.",
        "<b>Every joule is from the Sun:</b> in 309.92 J = heat 284.13 J + stored 25.781 J (integration error 0.0064 J); a "
        "matched load takes half the input (406.9 W) and never more.",
        "<b>The bifilar turn gap is overstressed:</b> 14.9 MV/m between turns against the declared 10 MV/m limit.",
        "<b>No multipactor at 566 Hz</b> (growth per impact 0.063–0.071); avalanche only at 1 GHz and 3 GHz, where the 8 km "
        "wire is many wavelengths long.",
    ])

    r.h1("What is being tested")
    r.p("The idea under test is that winding a tether into a coil, or as a bifilar coil whose adjacent turns differ by half "
        "the coil voltage, multiplies what it can harvest. <b>Standard prediction:</b> in a uniform field the EMF between a "
        "wire's ends is (v×B)·(r<sub>end</sub> − r<sub>start</sub>) whatever the path, so coiling cancels it; a resonant coil "
        "multiplies voltage by Q but stores only what its drive delivers, and the bifilar winding changes where the energy "
        "is stored, not how much. This example computes all of it with the energy ledger on. The multipactor check between "
        "turns is the one built for double_hull (§3i), and the burst-speed limit is the rotor rule of §3l.")

    r.h1("Declared inputs")
    r.table([
        ["input", "value", "note"],
        ["wire", "5 international miles of 2 mm copper", "8,049 m wound (1,281 whole turns), 43.0 Ω DC"],
        ["orbit speed", "7,700 m/s", "relative to the co-rotating field"],
        ["Earth's field B", "30 µT", "order of magnitude in LEO; a real orbit takes B from IGRF"],
        ["solar array", "2 m² at 30 %, TSI 1,360.8 W/m²", "Kopp &amp; Lean 2011; 816 W to the drive"],
        ["coil", "2 m across, pitch 3 mm", "1 mm bare gap between turns, vacuum (ε = 1); 3.84 m long"],
        ["gap limit", "10 MV/m", "to be replaced by the chosen insulation's data-sheet value"],
    ], [0.18, 0.36, 0.46])

    r.h1("The model")
    r.bullets([
        "<b>Windings</b> (em.Coil, the lumped first approximation of §3e): Wheeler inductance, 1.362 H for either winding "
        "(the two strands of a series bifilar coil carry the same current the same way); the two-wire capacitance between "
        "adjacent turns; Terman's skin-effect resistance; and the patent's turn voltage, V/N ordinary and V/2 bifilar. "
        "Referred to the terminals, the turn-to-turn self-capacitance is (N − 1)·C<sub>tt</sub>/N² for the ordinary coil "
        "(1.42×10<super>−13</super> F) and (N − 1)·C<sub>tt</sub>/4 for the bifilar coil (5.81×10<super>−8</super> F). Their "
        "ratio is Tesla's N²/4, 410,240 for 1,281 turns (the patent's 250,000 at 1,000 turns).",
        "<b>Earth's field.</b> The motional EMF ∫(v×B)·dl is integrated along the actual wire path with em.MotionalEMF "
        "(the Lorentz force by grade selection, §3f): a straight line, a helix of 36 points per turn, and for the bifilar "
        "coil two interleaved strands with the return lead. A fourth path folds the straight tether in half, out and back "
        "1 cm apart. In the craft's frame the same coupling is the electric field E′ = v×B, not a second source.",
        "<b>The Sun.</b> A sine source on the light channel in series with the wire resistance, the coil inductance and the "
        "capacitance across the coil (the bifilar coil's own, or a tank capacitor of the same 5.81×10<super>−8</super> F for "
        "the ordinary coil), with an optional series load. The amplitude √(2(R<sub>wire</sub> + R<sub>load</sub>)·P<sub>array</sub>) "
        "makes the loop draw the array's whole 816 W at resonance. Each run is a kronspice transient at 200 points per period "
        "for seven ring-up time constants (τ = 2L/R); peaks and mean powers are taken over the last 20 cycles.",
        "<b>Turn gap.</b> Voltage, field and stored energy between neighbouring turns; em.Multipactor over every amplitude "
        "the gap sees during ring-up, with measured as-received copper yield (Baglin et al. 2002), keeping the worst growth.",
    ])

    r.h1("Results")
    r.h2("1. Earth's field: the EMF along the real wire path")
    r.table([
        ["winding", "EMF", "reason"],
        ["straight 5-mile tether", "1858.792 V", "E′ × 8 km"],
        ["coiled (ordinary)", "0.888 V", "E′ × the coil's 3.84 m length: the turns cancel"],
        ["coiled (Tesla bifilar)", "0.889 V", "both strands run the same way through the field"],
        ["straight, folded out and back", "1.14×10<super>−13</super> V", "out and back cancel"],
    ], [0.3, 0.2, 0.5])
    r.p("The dipole field's gradient across the orbit adds at most N·A·v·3B/r ≈ 0.00041 V, and all of it is quasi-DC: the "
        "field seen by the craft changes over an orbit (~90 min), far below any LC resonance. A coil spun at 1 rev/s in the "
        "field gives 0.759 V peak at 1 Hz, paid for by the spin, which slows.")
    r.h2("2. The Sun drives the coil at resonance")
    r.p("L 1.362 H, C 5.81×10<super>−8</super> F, f 565.7 Hz, wire R 43.0 Ω; the wire is 0.015 λ, so the lumped model holds "
        "(valid below 0.1).")
    r.table([
        ["", "ordinary coil + tank capacitor", "Tesla bifilar", "bifilar, matched load"],
        ["Q", "112.5", "112.5", "56.2"],
        ["drive amplitude", "265.1 V", "265.1 V", "374.9 V"],
        ["coil voltage (peak)", "29,784 V (×112)", "29,784 V (×112)", "21,064 V (×56)"],
        ["coil current", "6.15 A", "6.15 A", "4.35 A"],
        ["between neighbouring turns", "23.3 V, 2.33×10<super>4</super> V/m", "14,892.2 V, 1.49×10<super>7</super> V/m",
         "10,532.1 V, 1.05×10<super>7</super> V/m"],
        ["energy between turns / total stored", "6.283×10<super>−5</super> J / 25.783 J", "25.78 J / 25.783 J", "12.89 J / 12.895 J"],
        ["multipactor growth per impact (worst)", "0.063", "0.070", "0.070"],
        ["steady state", "Sun in 815.2 W = heat 814.2 W", "Sun in 815.2 W = heat 814.2 W",
         "Sun in 815.2 W = heat 406.9 W + load 406.9 W"],
        ["whole run: in = heat + stored (+ error)", "309.92 J = 284.13 + 25.781 J (+0.0064 J)", "same",
         "155.13 J = 142.23 + 12.894 J (+0.0032 J)"],
    ], [0.25, 0.25, 0.25, 0.25])
    r.figure(os.path.join(HERE, "..", "ring_up.png"),
             "Figure 1. The bifilar coil driven by 816 W of sunlight at 566 Hz. Top: coil voltage ringing up from a 265 V "
             "drive to about 30 kV (×112). Bottom: energy in from the Sun, energy turned to heat in the wire, and energy "
             "stored in the coil; the stored energy levels off at about 26 J.", 0.78)
    r.p("Multiplication is real for voltage (by Q), and the bifilar winding puts the whole stored energy between its turns, "
        "N²/4 times the ordinary coil's turn-to-turn energy, with the same total. All of it comes from the Sun's array; "
        "the matched load takes half the input and never more. The bifilar turn field exceeds the declared 10 MV/m gap "
        "limit, so it needs insulation or wider spacing.")
    r.h2("3. Driving the resonance from Earth's field instead")
    r.p("The coil would have to spin at 566 rev/s: rim speed 3554 m/s against 1881 m/s, the burst speed of T1000G carbon "
        "fibre (the strongest rotor material in spacetime/). It would fly apart.")
    r.h2("4. Where electrons would multiply in the 1 mm vacuum gap")
    r.table([
        ["frequency", "worst growth per impact", "at amplitude", "verdict"],
        ["566 Hz", "0.071", "480 V", "no avalanche"],
        ["1 MHz", "0.067", "2000 V", "no avalanche"],
        ["100 MHz", "0.079", "2000 V", "no avalanche"],
        ["1 GHz", "1.657", "120 V", "avalanche"],
        ["3 GHz", "2.029", "480 V", "avalanche"],
    ], [0.2, 0.3, 0.2, 0.3])
    r.p("Worst amplitude scanned up to 2 kV. At gigahertz the electrons take their energy from the drive, and the 8 km wire "
        "is many wavelengths long, so it is not a coil resonance.")

    r.h1("Ledger and validation")
    r.bullets([
        "TestMotionalEMF: straight = vBL to 10<super>−12</super>; each coil = vB × its axial length; folded = 0.",
        "TestResonance: step-up = Q within 1 %; Sun in = array power = heat within 1 %; ledger closes to 10<super>−3</super>; "
        "turn-energy ratio = N²/4 to 10<super>−9</super> with equal totals; matched load takes half the input.",
        "TestRuns: the whole report runs.",
        "Whole-run integration error 0.0064 J of 309.92 J (ordinary and bifilar) and 0.0032 J of 155.13 J (matched load).",
        "The export sheets carry only the MNA column: automatic tearing finds no valid cut of this single loop (every torn "
        "piece is singular, Aitchison 1983 Lemma 2), which the sheet reports instead of solving.",
    ])

    r.h1("Limits of the model")
    r.bullets([
        "B = 30 µT is an order of magnitude; a real orbit takes B from IGRF.",
        "The coil is the lumped first approximation of §3e (turn-to-turn capacitance model, Wheeler inductance); the "
        "distributed and helical-line models are not used here.",
        "The gap limit (10 MV/m) is a placeholder for the chosen insulation's data-sheet value.",
    ])

    r.h1("Reproducing the results")
    r.code("go run ./examples/coiled_tether\n"
           "go test ./examples/coiled_tether\n"
           "KRONSPICE_EXPORT=1 go test ./examples/coiled_tether -run TestExport   # netlists, schematics, results sheets\n"
           "python examples/coiled_tether/viz/make_report.py examples/coiled_tether/coiled_tether_report.pdf")
    r.small("Sources: Kopp &amp; Lean, Geophys. Res. Lett. 38, L01706 (2011); Tesla, US 512,340 (1894), p. 2, ll. 3–7; "
            "Wheeler, Proc. IRE 16, 1398 (1928); Terman, Radio Engineers' Handbook (1943); Ramo, Whinnery &amp; Van Duzer, "
            "Fields and Waves, ch. 4; Joot, Geometric Algebra for Electrical Engineers, Thm 3.16; Baglin et al., LHC "
            "Project Report 472 (2002); Hatch &amp; Williams (1954); Vaughan, IEEE Trans. Electron Devices 35, 1172 (1988); "
            "Toray T1000G data sheet (2024); Aitchison, Quart. Appl. Math. 41, 265 (1983).")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "coiled_tether_report.pdf"))
