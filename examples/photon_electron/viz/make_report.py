# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""H23 report: the toroidal-photon electron and three electron pumps, as a PDF.

python examples/photon_electron/viz/make_report.py examples/photon_electron/photon_electron_report.pdf

Numbers come from examples/photon_electron/run_output.txt (the program's
output) and README.md; figures are viz/round_trip.png, viz/knot_light.png,
viz/pumps.png and viz/magnetron.png. Layout from internal/pdfreport/kronreport.py.
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
    r = Report("examples/photon_electron", "The electron as a photon on a double loop",
               "Williamson &amp; van der Mark's toroidal-photon electron on the test bench, and the electron "
               "pumps circuits can build — register entry H23, CLAUDE.md §3i, §3j",
               "register entry H23", "magnetron GA seed 2450, 40 × 30")

    r.h1("Summary")
    r.p("Williamson and van der Mark (1997) proposed that the electron is one photon of the Compton wavelength "
        "going twice round a torus. This example puts the paper's numbers beside measurement, builds its "
        "720° double loop as a circuit (a two-wire ring line with its conductors crossed once), and tests the "
        "paper's central assumption — that the light stays in the loop — with the field solver "
        "<i>em/fdtd.go</i>. It then builds three electron pumps that real circuits can make, each with an "
        "energy ledger, and evolves a multipactor pump driven by a microwave-oven magnetron.")
    r.bullets([
        "<b>The model's charge is 9.0 % low:</b> q = 0.910318 e, so its magnetic anomaly (g − 2)/2 = 0.000962436 is "
        "17.0 % below the measured 0.00115965218.",
        "<b>The 720° property is reproduced as a circuit:</b> the crossed ring line has half-integer modes "
        "(54,548 Hz = ½ wavelength round the ring) and a pulse returns inverted after one trip (−0.91 of launch) "
        "and upright after two (+0.70).",
        "<b>Light does not stay in the loop:</b> under linear Maxwell equations the energy within 2R falls from "
        "98.93 % to 5.773 % after a quarter period, and 99.99 % has left the region after 2 periods "
        "(ledger residual −1.9×10<super>−13</super>). The paper assumes a confinement it does not supply.",
        "<b>Single-electron pump:</b> 1.000000000 electrons per cycle at 1 MHz (I = e·f = 0.1602 pA), the gates "
        "paying the heat.",
        "<b>Resonant ground pump:</b> 10 V at 100 kHz → 199.9 V peak → bank −1194.7 V against the ideal −1199.5 V; "
        "the earth lends 7.46×10<super>14</super> electrons and gets the same number back.",
        "<b>Magnetron-driven multipactor pump, evolved:</b> 7.475 µA → 659.4 µA collected (176.6 → 2.002 mains W per "
        "µA). Arrays of gaps reach 1.95 A (LG 2M214, 1.32 kW mains) and 216 A (Stellant CWM-100kW, 119 kW mains).",
    ])

    r.h1("What is being tested")
    r.p("<b>Hypothesis (H23).</b> J. G. Williamson and M. B. van der Mark, “Is the electron a photon with toroidal "
        "topology?”, <i>Ann. Fond. Louis de Broglie</i> 22, 133 (1997). One photon of wavelength λ = h/(m<sub>e</sub>c) "
        "is confined to a closed path that goes twice round a torus. The paper assumes the confinement; it adds no "
        "field equation. From the field in a sphere of diameter λ/2 it states q = (1/2π)·√(3ε<sub>0</sub>ħc) (eq. 5), "
        "the coupling α′ = q²/(4πε<sub>0</sub>ħc) = (q/e)²α (eq. 8) and g = 2(1 + α′/2π) (eq. 12). Spin ½ comes from "
        "the double loop: the field returns to itself only after 720°.")
    r.p("<b>Standard prediction.</b> The Dirac electron with QED corrections, which gives the measured charge, "
        "magnetic moment and spin. Scattering shows the electron point-like below about 10<super>−18</super> m. "
        "Under linear Maxwell equations nothing holds light in a loop, so it leaves at the speed of light.")

    r.h1("Part 1. The paper's numbers against measurement")
    r.table([
        ["quantity", "model", "measured (CODATA 2018; Fan et al. 2023)"],
        ["charge, q = (1/2π)·√(3ε<sub>0</sub>ħc) (eq. 5)", "0.910318 e", "1 e (exact)"],
        ["coupling α′ = (q/e)²·α (eq. 8)", "0.00604716", "α = 0.00729735"],
        ["magnetic anomaly (g − 2)/2 = α′/2π (eq. 12)", "0.000962436", "0.00115965218"],
        ["… with q set to e (= α/2π, QED first order)", "0.001161410", "off by 0.0015"],
        ["spin", "1/2 (double loop)", "1/2"],
        ["size: Compton wavelength λ", "2.426×10<super>−12</super> m", "point-like to ~10<super>−18</super> m"],
        ["loop radius λ/4π", "1.931×10<super>−13</super> m", ""],
    ], [0.45, 0.22, 0.33])
    r.p("The model anomaly is 17.0 % below the measured one because the charge is 9.0 % below e.")

    r.h1("Part 2. The double loop as a circuit: plain and crossed ring lines")
    r.p("A two-wire ring transmission line of 200 LC sections is an LC ladder closed on itself (a closed path in "
        "Kron's sense). It is built once plain and once with its two conductors swapped at one point. On the plain "
        "ring a wave returns to itself after one trip, so its modes hold a whole number of wavelengths. On the "
        "crossed ring the differential voltage changes sign at the crossing, so a wave matches itself only after two "
        "trips and the modes hold half-integer numbers of wavelengths. Both rings are solved in kronspice: an AC "
        "sweep finds the modes (the twin is the ladder's analytic dispersion) and a transient sends a 2.93 V pulse "
        "round. One trip takes 9.17 µs.")
    r.table([
        ["#", "plain: AC peak", "ladder theory", "rel. diff", "λ round ring",
         "crossed: AC peak", "ladder theory", "rel. diff", "λ round ring"],
        ["1", "109102 Hz", "109104 Hz", "−2.66×10<super>−5</super>", "1.0",
         "54548.1 Hz", "54553.9 Hz", "−1.06×10<super>−4</super>", "0.5"],
        ["2", "218181 Hz", "218182 Hz", "−6.63×10<super>−6</super>", "2.0",
         "163646 Hz", "163648 Hz", "−1.18×10<super>−5</super>", "1.5"],
        ["3", "327205 Hz", "327206 Hz", "−2.95×10<super>−6</super>", "3.0",
         "272701 Hz", "272702 Hz", "−4.24×10<super>−6</super>", "2.5"],
        ["4", "436148 Hz", "436149 Hz", "−1.65×10<super>−6</super>", "4.0",
         "381688 Hz", "381689 Hz", "−2.17×10<super>−6</super>", "3.5"],
    ], [0.06, 0.12, 0.12, 0.13, 0.09, 0.12, 0.12, 0.13, 0.09])
    r.table([
        ["pulse back at the start", "trip 1", "trip 2", "trip 3", "trip 4"],
        ["plain ring", "+2.65 V (+0.90)", "+2.02 V (+0.69)", "+1.88 V (+0.64)", "+1.66 V (+0.57)"],
        ["crossed ring", "−2.68 V (−0.91)", "+2.06 V (+0.70)", "−1.72 V (−0.59)", "+1.52 V (+0.52)"],
    ], [0.28, 0.18, 0.18, 0.18, 0.18])
    r.p("Fractions in brackets are of the 2.93 V launch. Ledgers: plain ring, in 3.438×10<super>−8</super> J = heat "
        "1.096×10<super>−8</super> J + stored 2.341×10<super>−8</super> J; crossed ring, in 3.438×10<super>−8</super> J "
        "= heat 1.117×10<super>−8</super> J + stored 2.32×10<super>−8</super> J; time-stepping error "
        "1.3×10<super>−11</super> J in each.")
    r.figure(os.path.join(HERE, "round_trip.png"),
             "Figure 1. Differential voltage at the launch point. Top: the plain ring returns the pulse the same way up "
             "every trip (360°). Bottom: the crossed ring returns it inverted after one trip and upright after two "
             "(720°), the paper's spin-½ property in circuit form.", 0.9)

    r.h1("Part 3. Does light stay in the loop?")
    r.p("One wavelength of light is set on the paper's (2, 1) torus-knot path — twice round a ring of radius "
        "R = λ/4π while once round a tube of radius R/2 — as a travelling wave cos(ks − ωt) in a Gaussian tube of "
        "width R/3, with E and B across the path so the energy flows along it. E and B are built as discrete curls of "
        "sampled potentials, so the start is exactly divergence-free. Grid: λ = 1 m, R = 0.07958 m, cell 0.00663 m "
        "(R/12), region ±6R with an absorbing layer, path length 1.035 λ. The equations are linear, so lengths scale "
        "out. The solver is validated in <i>em/fdtd_test.go</i> (TM<sub>110</sub> cavity to 10<super>−7</super>, "
        "Hertzian dipole within 0.5 % of Larmor, flying doughnut within 3.8 %).")
    r.table([
        ["time (periods T)", "0", "0.125", "0.25", "0.5", "1", "1.5", "2"],
        ["energy within 2R", "98.93 %", "65.11 %", "5.773 %", "1.594×10<super>−6</super> %",
         "0.0001033 %", "0.00173 %", "5.549×10<super>−7</super> %"],
        ["in the region", "100 %", "100 %", "100 %", "90.42 %", "0.3495 %", "0.07582 %", "0.01209 %"],
        ["left the region", "0 %", "0 %", "4.038×10<super>−11</super> %", "9.584 %", "99.65 %", "99.92 %", "99.99 %"],
    ], [0.19, 0.09, 0.09, 0.13, 0.14, 0.12, 0.11, 0.13])
    r.p("Ledger residual −1.9×10<super>−13</super> of the starting energy; div B 7×10<super>−13</super> (round-off).")
    r.figure(os.path.join(HERE, "knot_light.png"),
             "Figure 2. Left: the paper's path, twice round the ring (dotted: the ring's centre line). Right: share of "
             "the starting energy still within 2R of the loop (red) and that has left through the absorber at ±6R "
             "(blue); the dashed line is what the model needs, 100 % staying for ever.", 0.9)

    r.h1("Part 4. Three electron pumps that work")
    r.p("The model's own frequency, m<sub>e</sub>c²/h = 1.24×10<super>20</super> Hz, is about 10<super>8</super> "
        "times beyond any circuit, so these pumps run at frequencies circuits reach. Electrons come from metal or "
        "ground; energy comes from a declared source.")
    r.h2("4.1 Single-electron pump (Pothier et al. 1992)")
    r.p("Two aluminium islands between three tunnel junctions of 200 aF and 200 kΩ, two 50 aF gates driven 90° "
        "apart, at 20 mK; charging energy e²/2C = 0.178 meV (2.1 K). Solved with the orthodox theory of "
        "single-electron tunnelling: a master equation over the islands' electron numbers with rates "
        "Γ = −ΔF/(e²R(1 − e<super>ΔF/kT</super>)), ΔF being the free-energy change of a tunnelling event.")
    r.table([
        ["drive", "electrons per cycle", "current", "heat per electron (paid by the gates)", "ledger"],
        ["1 MHz", "1.000000000", "0.1602 pA", "12.7 µeV", "−4.1×10<super>−10</super>"],
        ["10 MHz", "0.999999977", "1.602 pA", "53.1 µeV", "−4.1×10<super>−8</super>"],
        ["30 MHz", "0.999458493", "4.804 pA", "93.2 µeV", "−3.7×10<super>−7</super>"],
        ["100 MHz", "0.950531648", "15.23 pA", "160 µeV", "−4.1×10<super>−6</super>"],
        ["300 MHz", "0.667354712", "32.08 pA", "204 µeV", "−1.0×10<super>−5</super>"],
    ], [0.13, 0.2, 0.15, 0.32, 0.2])
    r.p("Errors grow when the drive outruns the tunnelling (RC ≈ 40 ps). Uphill against 20 µV at 10 MHz it still "
        "pumps 0.999998 electrons per cycle; the gates supply 72.8 µeV per cycle = 20.0 µeV lifted into the bias "
        "source + 52.8 µeV heat (ledger −5.6×10<super>−8</super>).")
    r.h2("4.2 Multipactor pump in a 1 GHz resonator (Farnsworth US 2,091,439; Kishek &amp; Lau 1998)")
    r.p("A 1 mm, 1 cm² copper gap (as-received yield, Baglin 2002) is the capacitor of a resonator (C = 0.885 pF, "
        "L = 28.6 nH, Q = 1000) driven by 0.779 mA for 140 V unloaded; a 10 % grid passes electrons to a collector. "
        "The resonator alone agrees with a kronspice transient to 0.0031 of V<sub>0</sub>. In step (140 V, 1000 seed "
        "electrons) the population grows to 1.98×10<super>8</super>, pulls the gap voltage down to 39 V and stops "
        "multiplying; the resonator recovers and a stray electron (one per 20 ns while the gap is empty, declared; 7 "
        "used) restarts it, giving bursts as Kishek and Lau found. Out of step (300 V) no electrons are left after 40 "
        "periods. Over 2000 periods 2.16×10<super>8</super> electrons are collected, 15.2 µA on average, at 1426 eV "
        "of RF drive per collected electron.")
    r.p("Ledger [nJ]: drive 49.32 = resonator heat 22.98 + impacts on the plates 31.43 + collector 1.78 + stored "
        "change −6.876; residual −3.5×10<super>−7</super>. The stored change is negative because the bursts also drain "
        "the energy the resonator held at the start.")
    r.h2("4.3 Ground charge pump with a resonant drive (§3j; Cockcroft &amp; Walton 1932)")
    r.p("10 V at 100 kHz into a series resonator (Q = 20) → 3-stage negative ladder → 100 nF bank with a 10 MΩ "
        "load, 600 cycles, solved by kronspice. Step-up 10 V → 199.9 V peak (×20.0); bank −1194.7 V against the "
        "ideal −2·n·peak = −1199.5 V. The earth lends 7.46×10<super>14</super> electrons to the bank's negative plate "
        "(plus 4.02×10<super>12</super> through the load) and the bank's earthed plate returns the same number, so the "
        "earth's net charge is unchanged. Ledger: source 3.559 J = heat 3.425 J + stored 0.133 J (bank 0.0714 J) + "
        "time-stepping error 2.7×10<super>−5</super> J.")
    r.figure(os.path.join(HERE, "pumps.png"),
             "Figure 3. (1) Single-electron pump: electrons moved per cycle against drive frequency; exactly one until "
             "the drive is too fast. (2) Multipactor pump: electrons in the gap (red, log scale) and gap voltage "
             "(blue) — bursts that drain the resonator and restart. (3) Ground pump: bank voltage approaching the "
             "ladder's limit −2·n·peak.")

    r.h1("Part 5. The multipactor pump on a microwave-oven magnetron, evolved")
    r.p("An LG 2M214 oven magnetron (2450 MHz, 900 W RF; efficiency 70 % declared, as the sheet gives it only as a "
        "chart; 34.6 W heater; 1320 W from the mains) feeds a circulator, a coupling port and a copper re-entrant "
        "cavity (Q<sub>0</sub> = 2000) whose capacitive gap multipacts, then a grid and a collector. Power the "
        "cavity does not take goes to the circulator's water load. The cavity's circuit twin is C = ε<sub>0</sub>A/d, "
        "L resonating it, wall loss R<sub>0</sub> = Q<sub>0</sub>/(ωC), port R<sub>ext</sub> = R<sub>0</sub>/β and a "
        "Norton source of the available power. The genetic algorithm (<i>evo</i>, value genome, seed 2450, 40 designs "
        "× 30 generations) chooses gap, gap area (side ≤ λ/10, so the gap stays lumped), coupling β, grid "
        "transparency and surface, for the fewest mains watts per collected microamp, below Kilpatrick's spark "
        "limit and 10 kV (non-relativistic). The search ran at ≥ 200 steps per period over 1500 periods; results were "
        "re-run at 800 steps over 2000 periods.")
    r.table([
        ["", "best of generation 0", "evolved"],
        ["gap; area", "2.761 mm; 0.376 cm²", "1.519 mm; 1.500 cm²"],
        ["coupling β; grid transparency", "0.000139; 0.03", "0.00153; 0.43"],
        ["surface", "copper, as received", "copper, as received"],
        ["gap amplitude, unloaded; peak field / Kilpatrick", "1038 V; 0.009", "1278 V; 0.020"],
        ["collector current", "7.475 µA", "659.4 µA"],
        ["mains W per µA; electrons per mains joule", "176.6; 3.53×10<super>10</super>",
         "2.002; 3.12×10<super>12</super>"],
        ["into the cavity: wall heat + electrons", "0.332 W: 0.297 + 0.0633", "2.05 W: 1.79 + 0.816"],
        ["to the water load", "899.7 W", "897.9 W"],
        ["ledger residual", "4.0×10<super>−6</super>", "3.9×10<super>−6</super>"],
    ], [0.44, 0.28, 0.28], bold_rows=(5,))
    r.p("The collector current varies by about ±10 % with step size and run length (burst timing); the order of "
        "magnitude is stable. An earlier search found a false 0.3 A from time steps too coarse for fast electrons; "
        "the engine now takes at least 40 steps per gap crossing. The gap needs only 2.1 W, so one magnetron could "
        "feed about 438 such gaps (321 µA per watt into the cavity).")
    r.figure(os.path.join(HERE, "magnetron.png"),
             "Figure 4. Left: best fitness (mains watts per collected µA) per generation, as measured during the "
             "search. Right: the evolved design's electrons in the gap (red) and gap voltage (blue) — bursts. The "
             "plotted search fitness (about 0.80 W/µA at the end) and the panel title (about 0.85 mA) are not the "
             "re-run values; the program's re-run at 800 steps per period gives 2.002 W/µA and 659.4 µA (table above).")

    r.h1("Part 6. Scaling to amps: many gaps on one magnetron")
    r.p("The genetic algorithm also chooses how many identical gaps (1 to 1,000,000) share the magnetron through an "
        "ideal, lossless, isolated power divider (declared), each taking P<sub>out</sub>/N.")
    r.table([
        ["", "LG 2M214 (oven)", "Stellant CWM-100kW (industrial)"],
        ["frequency, RF power", "2450 MHz, 1 kW", "915 MHz, 100 kW"],
        ["gaps in parallel", "2686", "786127"],
        ["gap; area per gap", "0.56 mm; 1.50 cm²", "1.52 mm; 3.96 cm²"],
        ["total gap area", "0.403 m²", "311 m²"],
        ["gap voltage (unloaded)", "191 V", "196 V"],
        ["collected current (per gap)", "<b>1.95 A</b> (726 µA)", "<b>216 A</b> (275 µA)"],
        ["into the gaps / water load", "0.695 / 0.205 kW", "76.4 / 23.6 kW"],
        ["mains", "1.32 kW", "119 kW"],
        ["amps per kW of mains; per m² of gap", "1.48; 4.84", "1.83; 0.696"],
        ["ledger residual (one gap)", "6.7×10<super>−7</super>", "6.1×10<super>−7</super>"],
    ], [0.4, 0.3, 0.3])
    r.p("Scaling the better design: 1 A needs about 0.548 kW from the mains and 1.44 m² of gap plates; 10 A about "
        "5.48 kW and 14.4 m²; 100 A about 54.8 kW and 144 m²; 300 A about 164 kW and 431 m². The searches were still "
        "improving when they stopped, so these are lower bounds. Feeding thousands of gaps needs a real "
        "power-distribution structure (a large multi-gap cavity rather than a divider), which is not modelled. The "
        "electrons come from the copper; every watt comes from the magnetron.")

    r.page_break()
    r.h1("The circuit mirror of the fdamador neutrino")
    r.p("The fdamador neutrino (H28, examples/fdamador_neutrino) is this example's field turned inside out: E is the "
        "ring and B circles it. In circuit form it is a series RLC loop: E round the ring is the capacitor, B through "
        "the hole the inductor, radiation a series resistor. Exchanging E and B is circuit duality (voltage and "
        "current, L and C, series and parallel exchange, and R becomes R0<super>2</super>/R), so the photon "
        "electron's mirror is a <b>parallel RLC tank</b>: B round the ring (toroidal flux) is the inductor, E through "
        "the hole the capacitor, radiation a shunt resistor.")
    r.table([
        ["element", "neutrino: series loop", "electron mirror: parallel tank"],
        ["inductor", "L = 21.56 nH (B through the hole)", "L = 21.56 nH (B round the ring)"],
        ["capacitor", "C = 0.2767 pF (E round the ring)", "C = 0.2767 pF (E through the hole)"],
        ["radiation", "R<sub>rad</sub> = 6.862 Ω in series", "R<sub>p</sub> = R0<super>2</super>/R<sub>rad</sub> "
                      "= 11.36 kΩ in shunt"],
        ["dielectric loss", "R<sub>d</sub> = 27.9 mΩ in series", "R0<super>2</super>/R<sub>d</sub> = 2.79 MΩ in shunt"],
        ["resonance, Q", "f0 = 2.0606 GHz, Q = 40.51", "f0 = 2.0606 GHz, Q = 40.51 (equal to 10<super>−12</super>)"],
    ], [0.2, 0.38, 0.42])
    r.p("The scale is R0 = √(L/C) = 279.1 Ω of the neutrino's bench loop (a = 10 mm, εr = 100; values from "
        "fdamador_neutrino's run_output.txt), which keeps L, C, f0 and Q. In the series loop every element carries "
        "the same current; in the tank every element sees the same voltage.")
    r.figure(os.path.join(HERE, "circuit.png"),
             "Figure 5. The fdamador neutrino's series RLC loop, the duality, and the photon electron's parallel RLC "
             "mirror, with values (viz/plot_circuit.py). An animation of one cycle, the light travelling along the "
             "(2, 1) path beside the mirror tank's energy moving between C and L, is viz/cycle.gif.")
    r.p("What the mirror is and is not. It is exact circuit algebra between the two <b>models</b>, not a statement "
        "about the measured electron or neutrino. The real dielectric torus is not dual (εr = 100, μr = 1): its "
        "B-ring mode rings at 4.82 GHz with Q 112 in the field solver (fdamador_neutrino Part 7b, grid r/2), not at "
        "the E-ring's 2.06 GHz. The H23 electron itself is a travelling wave on a ring line, with E and B rising and "
        "falling together (Part 2), whereas in the tank they alternate; a lumped tank stands for one of its modes at "
        "a time, and under linear Maxwell its light leaves within a quarter to half a period (Part 3). Neither "
        "circuit has mass or the weak interaction, and neither derives spin ½.")

    r.h1("Ledger and validation")
    r.table([
        ["test", "checks"],
        ["TestPaperNumbers", "eqs. 5, 8 and 12 of the paper"],
        ["TestRingModes, TestRoundTrips", "whole- and half-wave modes against the ladder; upright and inverted returns; ledgers close"],
        ["TestKnotLightLeaves", "on a coarse grid the light leaves, the field stays divergence-free, the ledger closes"],
        ["TestSingleElectronPump", "one electron per cycle at low drive, uphill too, the gates paying"],
        ["TestMultipactorPump", "resonator matches kronspice; growth in step, death out of step; ledger"],
        ["TestGroundPump", "the bank reaches −2·n·peak to 1 %; ledger"],
        ["TestMagnetronEvolution", "cavity matches kronspice; a short evolution improves on its start inside the limits"],
        ["TestGapArray", "N gaps are N copies of one gap fed P<sub>out</sub>/N; the 100 kW tube's declared numbers"],
    ], [0.32, 0.68])

    r.h1("Limits of the model")
    r.bullets([
        "Part 3 follows one initial field configuration under linear Maxwell equations; it shows the configuration is "
        "not self-confined, not that no confining mechanism could exist.",
        "The ring lines reproduce the 720° topology of the double loop, not the photon or the electron.",
        "The magnetron efficiency (70 %) and the lossless power divider are declared; the divider stands in for a "
        "multi-gap cavity that is not modelled. Collector currents vary by about ±10 % with step size.",
        "The stray-electron rate of the multipactor pump (one per 20 ns while the gap is empty) is declared.",
    ])

    r.h1("Status")
    r.p("H23 is <b>bounded</b>. The model charge is 0.910 e and its anomaly 17 % low against Fan et al. (2023); spin ½ "
        "from the double loop is reproduced as the crossed ring's half-integer modes and its sign flip per trip; the "
        "self-confinement is assumed, not derived, and under linear Maxwell the configuration leaves in less than "
        "half a period. The model needs a confining mechanism it does not give. Standard counterpart: the Dirac "
        "electron with QED.")

    r.h1("Reproducing the results")
    r.code("go run ./examples/photon_electron -viz examples/photon_electron/viz   # about 10 minutes\n"
           "python3 examples/photon_electron/viz/render.py                         # pictures\n"
           "go test ./examples/photon_electron                                    # about 1 minute\n"
           "KRONSPICE_EXPORT=1 go test ./examples/photon_electron -run TestExport  # netlists, results sheets\n"
           "python examples/photon_electron/viz/make_report.py examples/photon_electron/photon_electron_report.pdf")
    r.small("Sources: Williamson &amp; van der Mark, Ann. Fond. Louis de Broglie 22, 133 (1997); Fan et al., PRL 130, "
            "071801 (2023); CODATA 2018; Yee (1966); Hellwarth &amp; Nouchi, PRE 54, 889 (1996); Pothier et al., "
            "Europhys. Lett. 17, 249 (1992); Grabert &amp; Devoret (1992) ch. 2; Farnsworth US 2,091,439; Kishek &amp; "
            "Lau, PRL 80, 193 (1998); Baglin et al., LHC Project Report 472 (2002); Kilpatrick, Rev. Sci. Instrum. 28, "
            "824 (1957); Cockcroft &amp; Walton, Proc. R. Soc. A 136, 619 (1932); LG 2M214 and Stellant CWM-100kW "
            "sheets.")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "photon_electron_report.pdf"))
