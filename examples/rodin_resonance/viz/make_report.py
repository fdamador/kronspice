# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""H43 report: the Rodin coil in quadrature, its field, its resonance and its
output coils, as a PDF.

python examples/rodin_resonance/viz/make_report.py examples/rodin_resonance/rodin_resonance_report.pdf

Numbers come from examples/rodin_resonance/run_output.txt (the program's
output) and README.md; figures are viz/geometry.png, fields_xy.png,
fields_xz.png, energy.png, modes.png, map.png, bank.png and toroidal.png
(plot_all.py).
Layout from internal/pdfreport/kronreport.py.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401
import kronreport  # noqa: E402

kronreport.H1.keepWithNext = 1
kronreport.H2.keepWithNext = 1


def fig(name):
    return os.path.join(HERE, name)


def build(out):
    r = Report("examples/rodin_resonance", "The Rodin coil in quadrature",
               "A counter-wound Rodin pair driven by a sine and a cosine: its field, its resonance, the winding as a "
               "line on a knot, and output coils charging a capacitor bank — register entry H43",
               "register entry H43", "GA seeds 43 and 44, 64 × 50")

    r.h1("Summary")
    r.p("The Rodin coil winds its wire over a torus and through the hole, 150° round the axis per pass between 12 "
        "poles: 1→6→11→4→9→2→7→12→5→10→3→8→1. Seen from above it is a {12/5} star; in space one circuit goes five "
        "times round the axis and twelve times through the hole, a torus-knot path of the (5, 12) class. This "
        "example builds the coil characterized by Ramírez Galindo et al. (<i>Rev. Mex. Fís.</i> 70, 030901, 2024) as "
        "a counter-wound pair, drives the windings with a sine and a cosine, maps the field, builds it up by "
        "resonance, and evolves output coils, together with the drive frequency, that charge a capacitor bank. "
        "Every joule is booked against the two drive sources.")
    r.bullets([
        "<b>Validated against the paper:</b> one fitted turn count (12 circuits per winding of 0.6 mm wire) gives "
        "3.948 and 3.084 mT at 1.60 and 1.25 A against 4.13 and 3.05 mT measured (−4.4 %, +1.1 %). At 0.60 A the "
        "model gives 1.480 mT against 1.16: the paper's own data are not linear in the current.",
        "<b>Sine and cosine turn the field vector but do not move the pattern:</b> ellipticity 0.09–0.11 in the "
        "planes, while the 12-fold pattern round the axis advances by at most 3×10<super>−13</super>° per cycle.",
        "<b>The sources pass power through the coupling:</b> 0.4401 W from A to B (toroidal sense) and 0.4546 W from "
        "B to A (dipole sense), equal to ω·M·I<sub>A</sub>·I<sub>B</sub>·sin φ to four figures.",
        "<b>Energy accumulates only by resonance:</b> a tuned winding rises to its Q-limited level, Q = 76.4 from the "
        "rise against ω<sub>0</sub>L/R = 77.3, paid by the source. With k = 0.39 the pair cannot be tuned as a "
        "pair; tuned to one mode its currents lock in phase (−0.2°) under a sine/cosine drive.",
        "<b>A line on a knot near 2 MHz:</b> first self-resonance 1.966 MHz; the series counter-wound loop closed on "
        "itself is photon_electron's crossed ring, with 720° modes for the (5, 12) and the (2, 1) path alike.",
        "<b>Output coils:</b> the free search builds a transformer (one coaxial coil, 4.1 kHz; bank 0.176 of the "
        "source energy); two coils under the aperture at 29 kHz charge the bank with 0.321. The books close to the "
        "transient's integration error. Status: <b>bounded</b>.",
    ])

    r.h1("What is being tested")
    r.p("<b>Hypothesis (H43).</b> Hobbyist and promotional sources (M. Rodin; forum posts) say the Rodin winding makes a "
        "“vortex” field that accumulates energy, delivers more than its input, and follows a "
        "“vortex-based mathematics” (the doubling sequence 1-2-4-8-7-5). No independent measurement is "
        "given, and the mathematics gives no equation that a simulation could check.")
    r.p("<b>Standard prediction.</b> Two coupled inductors: a transformer, and near self-resonance a transmission "
        "line. The field is the Biot–Savart field of the currents; energy builds up only in a resonance, to the "
        "level where loss per cycle equals input; output is at most input, and the ledger closes.")

    r.h1("The relation to the photon electron (H23)")
    r.p("examples/photon_electron tests Williamson and van der Mark's electron, one photon of the Compton wavelength "
        "on a (2, 1) torus-knot path whose field returns to itself only after 720°. It builds that property as a "
        "circuit — a two-wire ring line with its conductors crossed once, whose modes are half-integer and whose "
        "pulse returns inverted after one trip — and adds the circuit mirror of the fdamador neutrino: the neutrino's "
        "field (E round the ring, B through the hole) is a series RLC loop, and the electron's mirror is a parallel "
        "RLC tank.")
    r.p("The Rodin winding is the same family of paths with other winding numbers, (5, 12) instead of (2, 1). It "
        "carries a current in a wire, not light, so the wire does the confining. At the drive frequencies used here "
        "(1 kHz, 100 kHz) the 19.2 m of wire per winding is far shorter than a wavelength, so the winding is a lumped "
        "inductor. Only near its self-resonance (1.966 MHz) does it become a transmission line laid on a knot. The "
        "counter-wound pair connected in series and closed on itself <i>is</i> the crossed ring: its lowest modes "
        "come back inverted after one trip through one winding. They do so for the (5, 12) Rodin winding and for a "
        "(2, 1) winding on the same former alike, so the 720° property belongs to the connection, not to the knot's "
        "winding numbers. In the mirror language, a winding tuned by a series capacitor is the series RLC loop.")
    r.figure(fig("geometry.png"),
             "Figure 1. Left: the Rodin pair from above — winding A (1→6→11…) and its mirror B, the 12 poles, and the "
             "Part 3 output coils beside the toroid (green: free search; violet: distributed groups). Middle: one circuit of each "
             "winding in 3-D. Third: a (2, 1) winding, the class of photon_electron's path. Right: the Part 3t "
             "toroidal output winding, 12 sections of 39 turns, every turn round the tube through the hole.", max_height_mm=70)

    r.h1("Part 1. Geometry and validation")
    r.p("Taut-wire path: from above, each pass is the chord from pole P<sub>k</sub> to P<sub>k+5</sub>, the first "
        "half on the top surface of the torus and the second underneath, changing side through the hole. Former "
        "R<sub>0</sub> = 3.2125 cm, a = 1.7625 cm. The paper gives no turn count or gauge: 0.6 mm copper is "
        "declared, each winding lays its path Nc times in a round bundle, and Nc is fitted to the measured centre "
        "field (A and B in series, the paper's connection): Nc = 12, 2.4672 mT/A.")
    r.table([
        ["I [A]", "paper theory", "paper FEM", "measured", "kronspice", "vs measured", "vs theory"],
        ["1.60", "3.90", "4.08", "4.13", "3.948", "−4.4 %", "+1.2 %"],
        ["1.25", "3.01", "3.00", "3.05", "3.084", "+1.1 %", "+2.5 %"],
        ["0.60", "1.30", "1.10", "1.16", "1.480", "+27.6 %", "+13.9 %"],
    ], [0.1, 0.14, 0.13, 0.13, 0.14, 0.18, 0.18])
    r.table([
        ["check", "result"],
        ["5 cm ring vs Kirchhoff μ<sub>0</sub>R(ln(8R/r) − 7/4)", "3.42106 vs 3.42147 ×10<super>−7</super> H (−1.2×10<super>−4</super>)"],
        ["L<sub>A</sub>, L<sub>B</sub>, M, k", "263.44 µH, 310.53 µH, +111.97 µH, 0.3915"],
        ["series inductance, toroidal / dipole sense", "797.91 / 350.03 µH"],
        ["reciprocity M<sub>AB</sub> vs M<sub>BA</sub>", "111.9694754 µH both (2.4×10<super>−16</super>)"],
        ["flux of A's field through a 1 cm circle vs Neumann", "516.443 vs 516.388 nH (1.1×10<super>−4</super>)"],
        ["the same, circle round the tube (3 mm from the bundles)", "−1123.68 vs −1090.02 nH (3.1 %)"],
        ["field energy outside a 2.5 mm tube vs ½L<sub>ext</sub>I², (2, 1) circuit", "125.02 vs 125.7 nJ (−0.54 %)"],
        ["the same, one Rodin circuit (tubes overlap at crossings)", "695.14 vs 743.79 nJ (−6.54 %)"],
        ["wire resistance A / B", "1.1414 Ω (19.209 m) / 1.3212 Ω (22.235 m)"],
    ], [0.52, 0.48])

    r.h1("Part 2. Sine on A, cosine on B")
    r.p("Lumped two-winding circuit in kronspice, sources 3.383 V peak at 1 kHz through 1 Ω (A alone carries 1.25 A). "
        "<b>Toroidal sense:</b> the through-the-hole parts of A and B add (field in the tube). <b>Dipole sense:</b> the "
        "round-the-axis parts add (field through the hole, the paper's connection).")
    r.table([
        ["drive", "|I<sub>A</sub>| A", "|I<sub>B</sub>| A", "P<sub>A</sub> W", "P<sub>B</sub> W", "A→B via M, W",
         "ωMI<sub>A</sub>I<sub>B</sub>sin φ", "W<sub>m</sub> µJ"],
        ["A only", "1.2499", "0", "1.6728", "0", "0", "0", "0.022–205.8"],
        ["B only", "0", "1.1156", "0", "1.4446", "0", "0", "0.026–193.2"],
        ["in phase, toroidal", "1.0706", "0.9521", "1.2544", "1.0253", "−0.0270", "−0.0268", "0.26–405.3"],
        ["in phase, dipole", "1.4327", "1.2944", "2.1860", "1.9569", "+0.0126", "+0.0122", "0.07–322.9"],
        ["sin/cos, toroidal", "1.4624", "0.8908", "1.8500", "1.3613", "+0.4401", "+0.4401", "104.3–260"],
        ["sin/cos, dipole", "1.0297", "1.3373", "1.5904", "1.6209", "−0.4546", "−0.4548", "117.1–247.2"],
    ], [0.2, 0.09, 0.09, 0.09, 0.09, 0.12, 0.13, 0.19])
    r.p("Every case closes its books over the last cycle (sources = copper + source resistors to ≤ "
        "1.6×10<super>−5</super> W) and over the run (ledger mismatch ≤ 9×10<super>−7</super> of the energy in).")
    r.table([
        ["drive", "centre peak mT", "tube, pole: peak mT (ellipticity)", "largest peak z=0 / y=0 mT",
         "mean B²/2μ<sub>0</sub> z=0 / y=0 J/m³", "ellipticity z=0 / y=0"],
        ["A only", "1.9339", "1.1292 (0.000)", "2.2098 / 2.1852", "0.1137 / 0.06669", "0 / 0"],
        ["in phase, toroidal", "0.7860", "1.8229 (0.002)", "3.5089 / 3.2838", "0.2606 / 0.08905", "0.003 / 0.007"],
        ["in phase, dipole", "3.4075", "0.1892 (0.090)", "3.4925 / 3.6579", "0.05043 / 0.1572", "0.008 / 0.002"],
        ["sin/cos, toroidal", "2.6128", "1.3403 (0.065)", "2.8371 / 2.8258", "0.1676 / 0.1249", "0.090 / 0.094"],
        ["sin/cos, dipole", "2.3243", "1.2433 (0.078)", "2.4479 / 2.7285", "0.1434 / 0.1214", "0.108 / 0.098"],
    ], [0.18, 0.12, 0.2, 0.17, 0.18, 0.15])
    r.p("The quadrature pair: the toroidal and dipole parts of the field peak a quarter cycle apart, so the field at "
        "a point traces an ellipse; the pattern itself keeps the poles' 12-fold symmetry at every instant and does "
        "not travel round the axis. A linear coil's field is the superposition of the drives at each instant; "
        "nothing accumulates from cycle to cycle without a resonance.")
    r.figure(fig("fields_xz.png"),
             "Figure 2. The field in the vertical plane y = 0 (±7 cm) for four drives: |B| peak over a cycle with "
             "the direction of the ellipse's major axis, rms, time-averaged energy density, and ellipticity (0: the "
             "field only swings; 1: it turns at constant size). viz/fields_xy.png shows the mid-plane; "
             "viz/cycle.gif animates one cycle of the sine/cosine drive.", max_height_mm=165)
    r.table([
        ["case at 100 kHz", "|I<sub>A</sub>| A", "|I<sub>B</sub>| A", "φ<sub>B</sub> − φ<sub>A</sub>",
         "mean W<sub>m</sub> µJ", "P<sub>in</sub> W", "Q = πf<sub>0</sub>τ", "ω<sub>0</sub>L/R"],
        ["A only, tuned to L<sub>A</sub>", "1.2480", "0", "—", "102.6", "1.6686", "76.4", "77.3"],
        ["A only, not tuned", "0.0162", "0", "—", "0.0172", "0.0003", "—", "—"],
        ["sin/cos, toroidal, each tuned to own L", "0.0392", "0.0368", "−90.2°", "0.2065", "0.0032", "—", "—"],
        ["sin/cos, toroidal, not tuned", "0.0204", "0.0175", "132.8°", "0.03753", "0.0008", "—", "—"],
        ["in phase, toroidal, tuned to L + M", "1.1957", "1.1958", "−0.1°", "285.2", "3.1930", "111", "—"],
        ["sin/cos, toroidal, tuned to L + M", "0.8594", "0.8327", "−0.2°", "142.5", "1.5969", "112", "—"],
    ], [0.3, 0.09, 0.09, 0.1, 0.1, 0.09, 0.12, 0.11])
    r.p("One tuned winding builds up over τ = 2L/R to the level where loss per cycle equals input: 5,965 times the "
        "untuned winding's stored energy at the same voltage, Q = 76.4 against 77.3 predicted. With k = 0.39 the "
        "pair has two modes, L ± M; tuning each winding to its own L leaves the drive between them, and tuned to "
        "the L + M mode the currents lock in phase under the sine/cosine drive: a resonant pair carries its mode, "
        "and the turning vector is gone.")
    r.figure(fig("energy.png"),
             "Figure 3. Left: stored magnetic energy over two cycles at 1 kHz. Middle: B(t) in the tube at pole 1 "
             "under the sine/cosine drive (toroidal sense). Right: the build-ups at 100 kHz (cycle maxima).",
             max_height_mm=60)

    r.h1("Part 2b. The winding as a line on a knot")
    r.p("Sectional model: every pass of every circuit is a section (144 per winding), the inductance matrix from the "
        "geometry, declared capacitances (188.1 pF/m between neighbouring wires of a bundle, A–B contact ×0.1, "
        "5.54 pF to the surroundings). Lossless eigenmodes; in brackets the trip parity (−1: inverted after one trip "
        "through one winding, the 720° modes).")
    r.table([
        ["winding", "open (far end open)", "series, toroidal", "series, dipole", "plain rings", "crossed ring"],
        ["Rodin (5, 12), 19.209 m", "1.966, 13.66 MHz", "0.5314, 1.266", "0.7052, 1.578",
         "9.534 (−1.00), 9.814 (−1.00)", "1.27 (−1.00), 6.258 (−1.00)"],
        ["(2, 1) knot, 5.091 m", "8.039, 50.51 MHz", "3.578, 4.091", "2.279, 10.34",
         "38.85 (−1.00), 40 (−1.00)", "4.144 (−1.00), 25.2 (−1.00)"],
    ], [0.19, 0.15, 0.14, 0.14, 0.19, 0.19])
    r.figure(fig("modes.png"), "Figure 4. The first six modes of each network, coloured by trip parity.", max_height_mm=60)

    r.h1("Part 3. Output coils round the toroid, the frequency evolved with them")
    r.p("Sine on A, cosine on B, 3.383 V peak each. Output coils of 0.5 mm copper, each tuned to the drive frequency by "
        "a series capacitor, coupled to every section of the pair through mutuals computed from the geometry; they "
        "load the drive in every solution. The sectional networks are solved by a complex modified-nodal solver "
        "written for the example (mutuals of either sign), twinned against a kronspice transient of the same "
        "network to 3×10<super>−3</super>.")
    r.h2("Where round the toroid, at which frequency")
    r.table([
        ["coil, sense", "EMF max/min over angle: 1 kHz", "100 kHz", "1.86 MHz", "5.28 MHz", "18.9 MHz",
         "largest load power"],
        ["tube, toroidal", "1.04 (m 12)", "1.05 (m 12)", "1.09 (m 12)", "4.13 (m 1)", "26.4 (m 5)", "0.03586 W, 5.68 kHz"],
        ["tube, dipole", "1.05 (m 12)", "1.05 (m 12)", "1.43 (m 2)", "112 (m 5)", "5.01 (m 1)", "0.03851 W, 4.51 kHz"],
        ["top, toroidal", "1.02 (m 12)", "1.05 (m 12)", "2.51 (m 5)", "106 (m 2)", "64.9 (m 7)", "0.03318 W, 18.9 MHz"],
        ["top, dipole", "1.07 (m 12)", "1.04 (m 12)", "1.24 (m 3)", "12.5 (m 7)", "57 (m 7)", "0.0331 W, 18.9 MHz"],
    ], [0.15, 0.17, 0.11, 0.11, 0.11, 0.11, 0.24])
    r.p("Below the first self-resonance every pass carries the same current, so the field keeps the 30° symmetry of "
        "the poles and a coil sees the same EMF within 2–7 % at every angle (harmonic m = 12): where a group sits "
        "round the toroid does not matter, and the sine/cosine pair gives no angle-dependent phase, since the "
        "pattern does not turn. Above it the current varies along the wire as a standing wave; the wire goes round "
        "the axis 5 times per circuit, 150° per pass, so the field round the toroid takes the low harmonics of that "
        "wave (m = 1, 2, 5, 7) and the EMF changes by up to 112× with angle, the best angles moving with "
        "frequency. These modes depend on the declared capacitances.")
    r.figure(fig("map.png"),
             "Figure 5. log<sub>10</sub> of the open-circuit EMF (top, full sectional net) and of the load power "
             "(bottom, coil in the network, reduced net) of a reference coil against the drive frequency and its "
             "angle round the toroid; tube coil and top coil, both connection senses.", max_height_mm=110)
    r.h2("The genetic searches")
    r.p("Genes: drive frequency (1 kHz–30 MHz), number of groups (1–8), first angle and spacing, coil tilt, position, "
        "radius (3–60 mm), turns (1–200), combination (each into its own bridge, all in series into one, all in "
        "parallel into one), load, connection sense. Constraints: coils clear of the bundles by 2 mm and of each "
        "other, shorter than λ/20, tuning capacitors ≥ 30 pF. Fitness: load power in the AC solution times a smooth "
        "bridge factor (V/(V + 0.6 V))². Searches run on 2 sections per circuit; the winners give the same power on "
        "4, 6 and 12.")
    r.table([
        ["", "free search", "distributed groups (≥ 2, radius ≤ 20 mm)"],
        ["fitness, generation 0 → evolved", "1.664×10<super>−7</super> → 0.181 W", "1.265×10<super>−7</super> → 0.1172 W"],
        ["drive frequency, sense", "4,134 Hz, dipole", "29.38 kHz, dipole"],
        ["groups", "1 coaxial coil round the whole coil, r 53.3 mm, 200 turns, 14.4 mm above the mid-plane",
         "2 at 79.7° and 251.1°, r 20 mm, 200 turns, tilted 50.1°, under the aperture (z −34.8 mm), in series"],
        ["coil L, R, tuning C", "7677 µH, 5.733 Ω, 193.1 nF", "1894 µH, 2.15 Ω, 15.5 nF"],
        ["AC: sources = load + copper + source R",
         "0.4708 + 0.2793 = 0.2522 + 0.307 + 0.1909 W", "0.2061 + 0.1462 = 0.2018 + 0.126 + 0.02447 W"],
        ["reaction: sources without → with the coils", "0.5877 + 0.03219 → 0.4708 + 0.2793 W",
         "0.058 − 0.04325 → 0.2061 + 0.1462 W"],
        ["separate / series / parallel bridges", "one coil: 0.2522 W", "0.1598 / 0.2018 / 0.1053 W"],
    ], [0.3, 0.35, 0.35])
    r.p("Left free, the search builds a transformer: one coil coaxial with the toroid, linking the dipole flux through "
        "the hole. Made to use several small coils, it puts two under the aperture where that flux leaves, in phase "
        "(138.3° both) and in series; in parallel they give half as much, on separate bridges 21 % less. Neither "
        "search goes above the self-resonance: the angle matters there, but the winding's input impedance leaves "
        "less power to take.")
    r.h2("The winners charging the bank")
    r.table([
        ["kronspice transient", "free search", "distributed groups"],
        ["run, bank", "400 cycles (96.77 ms), 2200 µF", "400 cycles (13.62 ms), 680 µF"],
        ["bank voltage at the end", "3.367 V", "1.86 V"],
        ["bank energy / source energy A + B", "12.47 mJ / 47.04 + 23.65 mJ", "1.176 mJ / 2.393 + 1.266 mJ"],
        ["charging efficiency", "0.1764", "0.3214"],
        ["copper / source R / diodes", "45.08 / 28.96 / 6.47 %", "34.76 / 8.10 / 21.65 %"],
        ["books: sources − (bank + losses + stored)", "35.6 µJ", "4.95 µJ"],
        ["ledger integration mismatch", "35.8 µJ (4.0×10<super>−4</super>)", "4.97 µJ (1.1×10<super>−3</super>)"],
        ["radiation of the output coil", "3.59×10<super>−15</super> Ω vs 5.73 Ω", "1.82×10<super>−13</super> Ω vs 2.15 Ω"],
    ], [0.36, 0.32, 0.32])
    r.figure(fig("bank.png"),
             "Figure 6. Bank voltage and energy of the two winners in the kronspice transient (Schottky bridges), and "
             "the best fitness of each genetic search by generation.", max_height_mm=60)

    r.h1("Part 3t — a toroidal output winding in 12 sections")
    r.p("The output winding wound round the torus as on a toroidal inductor core: every turn over the outside, "
        "through the hole and back under, laid on the two Rodin windings, in 12 sections 30° apart (arc ≤ 29°, "
        "0.5 mm copper, up to 6 layers; a 29° section holds 39 turns). Each turn's mutual to every Rodin segment is "
        "tabulated over 30° and carried round by the windings' 30° symmetry (pass q of A → q + 5, of B → q + 7). "
        "Checks: interpolation 1.6×10<super>−4</super>; a uniform 114-turn toroid 158.59 vs the current sheet "
        "158.76 µH; rings vs an explicit helix with its return lead, mutuals 0.25 %, self 3 %.")
    r.p("A toroidal turn collects only the toroidal flux of the through-the-hole parts. Under sine/cosine the two "
        "windings' parts are a quarter cycle apart, so a turn sees √(A² + B²) in either sense (11.10 mV per turn at "
        "1 kHz); driven in phase in the toroidal sense they add, A + B (15.65 mV). Round the torus the EMF varies by "
        "3–4 % (m = 12) below the self-resonance, so where the sections sit hardly matters (0.16 % in load power).")
    r.table([
        ["", "sine on A, cosine on B", "drive phase evolved too"],
        ["frequency, phase of B, sense", "1,575 Hz, +90°, dipole", "1,557 Hz, 356.4° (in phase), toroidal"],
        ["winding", "12 × 39 turns, arc 28.2°, series, one bridge", "12 × 39 turns, arc 29°, series, one bridge"],
        ["L, R of the whole winding", "2,974 µH, 6.15 Ω", "2,973 µH, 6.15 Ω"],
        ["AC: sources = load + copper + source R", "2.334 = 0.297 + 1.156 + 0.881 W", "1.462 = 0.570 + 0.562 + 0.331 W"],
        ["AC efficiency", "0.127", "0.390"],
        ["bank (kronspice, 4,700 µF, ~255 ms)", "4.300 V, 43.45 mJ", "6.097 V, 87.37 mJ"],
        ["charging efficiency", "0.0723", "0.2262"],
        ["ledger integration mismatch", "142 µJ (2.2×10<super>−4</super>)", "154 µJ (3.7×10<super>−4</super>)"],
    ], [0.34, 0.33, 0.33])
    r.p("The search fills every section with as many turns as fit and puts them in series into one bridge; separate "
        "or parallel bridges lose 4–50×, alternating sections cancel. Under sine/cosine the toroidal winding takes "
        "more power than Part 3's best outside coil (0.297 vs 0.252 W) at a third of its efficiency; left free to "
        "choose, the search drives the two windings in phase in the toroidal sense, where the same winding takes "
        "0.570 W at 0.390 and banks 7× the coaxial coil's energy. A toroidal winding cannot use the rotation of "
        "the field vector that sine and cosine produce; it uses only the toroidal flux.")
    r.figure(fig("toroidal.png"),
             "Figure 7. Part 3t: (a) the 12 sections on the pair, (b) a cross-section, (c) the EMF of one turn against "
             "angle at 1 kHz for every drive, (d) its (frequency, angle) map, (e) bank charging, (f) the searches.",
             max_height_mm=120)

    r.h1("Status of H43")
    r.table([
        ["claim", "what the runs show", "standard prediction"],
        ["a “vortex” field", "the vector turns at each point (ellipticity 0.09–0.11); the pattern stands still",
         "two Biot–Savart fields a quarter cycle apart"],
        ["the field accumulates", "only by resonance, Q 76.4 vs 77.3, paid by the source; a tuned pair locks in phase",
         "a series RLC tank"],
        ["over-unity", "every ledger closes; the bank holds 0.07–0.32 of the source energy (toroidal winding 0.07 under sine/cosine, 0.23 in phase)", "transformer: output ≤ input"],
        ["the (5, 12) path is special", "it sets L, k and frequencies; 720° modes for (5, 12) and (2, 1) alike",
         "a property of the connection"],
        ["vortex-based mathematics", "no equation given", "—"],
    ], [0.22, 0.48, 0.3])
    r.p("<b>Status: bounded.</b> Computed under the standard model with closed books, the Rodin pair is a pair of "
        "coupled inductors (k = 0.39) and, near 2 MHz, a line laid on a (5, 12) knot. What would move the entry: a "
        "bench measurement in which the energy delivered to a load or bank exceeds the energy measured at the drive "
        "terminals, with every input counted. Bench test: the free-search winner (or the paper's coil with a coaxial "
        "secondary), a power analyser on both drive channels, and the bank's ½CV².")

    r.h1("Approximations")
    r.bullets([
        "Turn count and gauge are not in the paper: 0.6 mm copper, 12 circuits per winding fitted to the measured "
        "centre field; each winding a round bundle on one path.",
        "Taut-wire path; crossings and curvature of a real winding idealized.",
        "Skin effect of an isolated wire, no proximity effect; air core.",
        "Sectional model with declared capacitances; results above the first self-resonance depend on them.",
        "Searches on 2 sections per circuit with a resistive load and a smooth bridge factor; the kronspice "
        "transient replaces them with real diodes.",
        "Radiation estimated (small loop), not modelled.",
        "Toroidal output winding: turns as circles at their azimuths, layers at their mean radius, leads twisted; "
        "no winding capacitance of the output winding.",
    ])

    r.h1("Run it")
    r.code("go run ./examples/rodin_resonance\n"
           "python3 examples/rodin_resonance/viz/plot_all.py\n"
           "python3 examples/rodin_resonance/viz/plot_cycle.py\n"
           "go test -short ./examples/rodin_resonance")

    r.h1("Sources")
    r.bullets([
        "A. D. Ramírez Galindo, H. A. Pérez Olivas, G. Basurto Islas, A. Hernández Rayas, F. González López, "
        "T. Córdova Fraga, “Vortex magnetic induction: Mathematical, geometric and experimental "
        "characterization”, <i>Rev. Mex. Fís.</i> 70, 030901 (2024), doi:10.31349/RevMexFis.70.030901.",
        "M. Rodin, promotional material on the Rodin coil and vortex-based mathematics: the hypothesis source.",
        "J. G. Williamson, M. B. van der Mark, <i>Ann. Fond. Louis de Broglie</i> 22, 133 (1997); "
        "examples/photon_electron (H23).",
        "F. W. Grover, <i>Inductance Calculations</i>, 1946; A. E. Ruehli, <i>IBM J. Res. Dev.</i> 16, 470 (1972); "
        "J. C. Maxwell, <i>Treatise</i> Vol. II art. 701; J. D. Jackson, <i>Classical Electrodynamics</i>, §5.2–5.3.",
        "C.-W. Ho, A. Ruehli, P. Brennan, modified nodal analysis, <i>IEEE Trans. Circuits Syst.</i>, 1975; "
        "G. H. Golub, C. F. Van Loan, <i>Matrix Computations</i>, §8.4.",
        "C. A. Balanis, <i>Antenna Theory</i>; F. E. Terman, <i>Radio Engineers' Handbook</i>, 1943.",
    ])
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "rodin_resonance_report.pdf"))
