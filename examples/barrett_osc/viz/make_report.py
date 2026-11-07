# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""H36 report: Barrett's oscillator-shuttle-circuits (Tesla OSC) as a PDF.

python examples/barrett_osc/viz/make_report.py examples/barrett_osc/barrett_osc_report.pdf

Numbers come from examples/barrett_osc/run_output.txt (the program's output)
and examples/barrett_osc/README.md; figures are viz/bench.png and
viz/tesla_coil.png. Layout from internal/pdfreport/kronreport.py.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401


def build(out):
    r = Report("examples/barrett_osc", "Barrett's oscillator-shuttle-circuits",
               "Are Tesla's OSCs nonlinear beyond Kirchhoff and Ohm? — register entry H36",
               "register entry H36", "component values declared in the code (the sources give none)")

    r.h1("Summary")
    r.p("T. W. Barrett (<i>Ann. Fond. Louis de Broglie</i> 16, 23, 1991; US 5,493,691, 1996) says that Tesla's "
        "circuits, which shuttle energy among several floating “grounds”, are nonlinear in a way Kirchhoff's and "
        "Ohm's laws cannot describe. He says they mix a pump and a signal into an idler that obeys the "
        "Manley–Rowe relations, phase-conjugate the signal, and carry the algebra of quaternions (SU(2)) and higher "
        "Cayley–Dickson numbers. This example builds his Figs. 2–5 from the drawings, with every component value "
        "declared, then adds the nonlinear parts a real build has: a Metglas core and a spark gap.")
    r.bullets([
        "<b>The circuits as drawn are linear.</b> In all seven OSCs of Fig. 2 (A–G) the transient load current "
        "equals the AC superposition of pump and signal to ≤ 1.9×10<super>−4</super>, the largest mixing product is "
        "≤ 1.1×10<super>−14</super> (round-off), and the ledger closes to ≤ 1×10<super>−7</super>.",
        "<b>A nonlinear store makes the idler.</b> With a varactor top-load, an idler at fp − fs = 370 kHz appears at "
        "4.063×10<super>−5</super> of the pump line; Manley–Rowe sums are 4.0×10<super>−7</super> and "
        "4.8×10<super>−7</super>; a +60° signal phase moves the idler to −62.71° (phase conjugation).",
        "<b>Standard scaling, not the paper's.</b> The idler goes as Vs<super>1.010</super>·Vp<super>1.027</super> "
        "(standard χ<super>(2)</super> mixing, Es*·Ep); the paper's Eq. 24 has |Es|²·Ep.",
        "<b>Fig. 3 is a working pulse former:</b> load lobes −6.71, +3.42, −1.37 kV; ledger residual "
        "1.0×10<super>−10</super> J on 150.00 mJ stored.",
        "<b>The coil networks of Fig. 4 are ordinary two-port algebra:</b> cascades do not commute (XY vs YX 0.406) "
        "but are associative (1.1×10<super>−10</super>) with det ABCD = 1, so they cannot carry Cayley "
        "non-associativity or zero divisors.",
        "<b>Real nonlinear parts do mix.</b> Air-core pancakes give no idler (≤ 6.91×10<super>−15</super> V); a "
        "biased Metglas core gives a 1.09 mV idler; a timed spark gap multiplies a 1 V top-load signal into "
        "sidebands every 1 kHz (0.43× the signal rms, ordinary secondary; 0.66×, bifilar).",
    ])

    r.h1("The claim and the standard prediction")
    r.p("<b>Hypothesis (H36).</b> OSCs, with coils shuttling energy among floating grounds, are nonlinear beyond "
        "Kirchhoff and Ohm; they mix pump and signal into a Manley–Rowe idler, phase-conjugate, and carry SU(2) "
        "and higher Cayley–Dickson algebra.")
    r.p("<b>Standard prediction.</b> Coils and capacitors, however their grounds float, form a linear network, and "
        "a linear network makes no new frequency. Parametric mixing, the Manley–Rowe relations and phase conjugation "
        "come from a nonlinear element (a varactor, a saturable core, a switch). Cascaded two-ports multiply as "
        "2×2 matrices, which do not commute but do associate.")
    r.p("Neither the paper nor the patent gives a component value, so every value below is declared in the code. "
        "Spark gaps (Fig. 3) use the gap model and solver of <i>examples/sparkgap_diodebridge</i>; everything else "
        "runs in kronspice.")

    r.h1("Fig. 2 — the OSCs A–G as drawn")
    r.p("Pump 10 V at 1 MHz from 10 Ω; pump coil α 10 µH, shuttle coil β 30 µH, load 50 Ω with its coil γ 5 µH, "
        "2 Ω per coil; A's store 1 nF; each floating terminal a 10 pF top-load with a 1 GΩ leak; coupling k = 0.3 "
        "(0.6 for G's pancakes); F's capacitor 1 nF. The signal (1 V, 630 kHz) enters at the shuttle's terminal.")
    r.table([
        ["Fig.", "load |I(fp)|", "load |I(fs)|", "transient vs AC (fp, fs)", "largest mixing product", "ledger"],
        ["A", "5.870×10<super>−2</super> A", "2.384×10<super>−6</super> A", "9.8e-05, 1.2e-05", "9.7e-15 at (+3,+0)", "7e-08"],
        ["B", "6.209×10<super>−2</super> A", "1.041×10<super>−3</super> A", "5.5e-05, 1.7e-05", "5.9e-15 at (−1,+3)", "1e-08"],
        ["C", "5.865×10<super>−2</super> A", "4.951×10<super>−6</super> A", "1.3e-04, 2.2e-05", "6.9e-15 at (+3,+3)", "1e-08"],
        ["D", "5.832×10<super>−2</super> A", "1.450×10<super>−5</super> A", "7.9e-05, 4.3e-05", "1.1e-14 at (+2,+3)", "1e-08"],
        ["E", "5.461×10<super>−2</super> A", "1.457×10<super>−6</super> A", "1.4e-04, 2.1e-05", "3.0e-15 at (+3,−1)", "7e-09"],
        ["F", "3.508×10<super>−2</super> A", "1.541×10<super>−6</super> A", "1.9e-04, 4.3e-05", "1.0e-14 at (+2,+0)", "1e-07"],
        ["G", "5.853×10<super>−2</super> A", "9.894×10<super>−6</super> A", "6.8e-05, 4.8e-05", "9.8e-15 at (+3,+1)", "1e-08"],
    ], [0.07, 0.17, 0.17, 0.2, 0.25, 0.14])
    r.p("Every circuit gives the pump and signal lines that AC analysis gives for each source alone, and no mixing "
        "product above round-off. <b>Bench test:</b> drive the pump, feed a signal to G2, and look for fp − fs in "
        "the load; the standard prediction is only what the parts' own nonlinearity gives (ferrite cores, ceramic "
        "dielectrics; air-core coils and film or vacuum capacitors give least).")
    r.h2("With a nonlinear store (varactor)")
    r.p("Fig. 2A with the terminal top-load a varactor, q = Ct·(v + v²/(2·20 V)), pump 2 V, signal 0.2 V:")
    r.bullets([
        "idler at fp − fs = 370 kHz at the varactor terminal, 4.063×10<super>−5</super> of the pump line;",
        "Manley–Rowe sums for the power into the varactor (relative, harmonics to 4th): 4.0×10<super>−7</super> and "
        "4.8×10<super>−7</super> (zero for a lossless nonlinear reactance);",
        "idler proportional to Vs<super>1.010</super> · Vp<super>1.027</super>: standard χ<super>(2)</super> mixing, Es*·Ep. "
        "The paper's Eq. 24 (patent Eq. 15) has |Es|²·Ep, which would give exponent 2 in the signal;",
        "signal phase +60° moves the idler phase to −62.71° (the exact conjugate is −60°): phase conjugation.",
    ])
    r.p("The mixing, Manley–Rowe and phase conjugation Barrett describes are those of any parametric mixer; they "
        "come from the nonlinear element, not from the floating grounds.")

    r.h1("Fig. 3 — switched −, +, − capacitors (spark gaps)")
    r.p("10 kV on 1 nF each, 1 µH into a 50 Ω load, 10 mm air gaps (breakdown 35.5 kV, so they hold the 20 kV "
        "between neighbours), triggered at 0, 100 and 200 ns.")
    r.table([
        ["quantity", "value"],
        ["load-voltage lobes, in order", "−6.71, +3.42, −1.37 kV"],
        ["strikes per gap", "49, 3, 1"],
        ["stored", "150.00 mJ"],
        ["to the load", "66.27 mJ"],
        ["in the gaps (arcs)", "76.49 mJ"],
        ["left in the capacitors", "3.802 mJ"],
        ["solver damping", "3.442 mJ"],
        ["ledger residual", "1.0×10<super>−10</super> J"],
    ], [0.5, 0.5])
    r.p("This is a lumped frozen-wave generator, the impulse source of later ultrawideband radar. When a gap fires, "
        "the next capacitor shares charge with the one before it with no inductance between them; that loss goes "
        "to the arcs, which is why the + and second − lobes are smaller. Transmission-line segments in place of "
        "capacitors avoid it.")

    r.h1("Fig. 4 — the coil networks (A quaternion, B Cayley, C sedenion-like)")
    r.p("α…ψ = 10…70 µH (1 mΩ each), 1 pF per node, 10 pF per open terminal, 5 MHz (declared). One-terminal "
        "networks are driven at their terminal against G1:")
    r.table([
        ["one-terminal network", "Y / (jωC<sub>term</sub>)", "largest coil current / terminal current"],
        ["A", "1.209029 − 1.60×10<super>−8</super>i", "1.73×10<super>−1</super>"],
        ["B(i)", "1.360694 − 5.44×10<super>−8</super>i", "2.65×10<super>−1</super>"],
        ["C(i)", "1.619236 − 1.80×10<super>−7</super>i", "3.82×10<super>−1</super>"],
    ], [0.3, 0.35, 0.35])
    r.p("With no capacitance at the inner nodes every coil current is zero: a network with one terminal is only its "
        "terminal's capacitance, so the coils (and the α, β, γ “products”) act only through stray capacitance, or as "
        "antennas in an outside field.")
    r.table([
        ["two-ports X = B(ii), Y = C(ii), Z = B(ii) with values reversed", "result"],
        ["cascade XY in kronspice vs ABCD(X)·ABCD(Y)", "1.9×10<super>−16</super>"],
        ["XY vs YX (order matters: “αβ ≠ βα”)", "0.406"],
        ["XYZ vs (XY)Z and X(YZ) (associative)", "1.1×10<super>−10</super>"],
        ["|det ABCD − 1| (reciprocal: invertible, no zero divisors)", "1.1×10<super>−16</super>"],
        ["S†S − I (lossless: S unitary, an element of U(2))", "1.8×10<super>−6</super>"],
    ], [0.7, 0.3])
    r.p("Order dependence is real and is ordinary two-port theory. The SU(2) Barrett names is the group of a "
        "lossless two-port's scattering matrix, up to phase. Box H's “(αβ)γ ≠ γ(αβ)” is a commutator, not an "
        "associator.")

    r.h1("Fig. 5 — the quaternion “idler”")
    r.p("x = x0 + x1·i + x2·j + x3·k with x0 = Es*·e<super>−iωs t</super>, x1 = Es·e<super>−iβs t</super>, "
        "x2 = √Ep·e<super>iωp t</super>, x3 = √Ep·e<super>−iβp t</super>. Each component's line at fp − fs is "
        "4.1×10<super>−17</super>, 4.1×10<super>−17</super>, 7.2×10<super>−17</super>, 7.2×10<super>−17</super> "
        "(relative): a sum of pump and signal holds no idler; one needs a product, i.e. a nonlinear element. "
        "x1, x3 and Eq. 24 also write β·t, but β is a propagation constant (rad/m), so β·t is not a phase; the "
        "travelling-wave phase is β·z.")

    r.h1("Bench version of Fig. 2G: bifilar pancakes and a Metglas core")
    r.p("Pump coil α and shuttle coil β are the flat spirals of <i>pulsed_bifilar</i> (214 turns of 0.5 mm copper, "
        "127 m): L 6.84 mH and 10.9 Ω each; self-capacitance 6.5 pF ordinary (first resonance 757 kHz) against "
        "3.03 nF bifilar (34.9 kHz). The bench therefore runs below 35 kHz: pump 10 kHz 30 V, signal 6.3 kHz 1 V "
        "(50 Ω each) into a 10 nF store, idler 3.7 kHz, load 50 Ω, coupling k 0.6 (declared). Metglas 2605SA1 "
        "(declared): magnetizing inductance ×4, saturating from 20 mA, core Q 50; ribbon k<sub>me</sub> 0.3, Q 300, "
        "0.54 m long for 3.7 kHz.")
    r.table([
        ["case", "load at fp", "fs [V]", "idler fp−fs", "2fp−fs [V]", "Manley–Rowe", "ledger"],
        ["ordinary pancakes, air", "2.989 V", "1.18e-02", "8.41e-16 V", "5.78e-16", "—", "3e-09"],
        ["bifilar pancakes, air", "2.571 V", "1.22e-02", "6.91e-15 V", "1.01e-14", "—", "5e-09"],
        ["bifilar, Metglas core", "0.461 V", "1.57e-02", "1.41e-07 V", "2.49e-03", "2e-03, −2e-05", "5e-04"],
        ["bifilar, Metglas, DC bias", "0.595 V", "1.59e-02", "1.09e-03 V", "1.81e-03", "1e-04, −1e-04", "2e-04"],
        ["bifilar, Metglas, bias, ribbon 3.7 kHz", "0.656 V", "1.54e-02", "2.88e-05 V", "1.84e-03", "1e-04, −2e-04", "3e-04"],
    ], [0.28, 0.1, 0.1, 0.13, 0.1, 0.17, 0.09], bold_rows=(4,))
    r.small("Values as printed in run_output.txt. The DC bias is 13 mA (README).")
    r.bullets([
        "<b>Air core:</b> ordinary or bifilar, the pancakes are linear and make no idler; the bifilar winding changes "
        "the self-capacitance 470×, so the lines differ, but adds no new frequency.",
        "<b>Unbiased Metglas:</b> the core is symmetric and makes only odd orders; the strongest product is the "
        "four-wave line 2fp − fs.",
        "<b>Biased Metglas:</b> a DC bias makes the even orders too, including the three-wave idler fp − fs. "
        "Manley–Rowe holds for the core's flux curve.",
        "<b>Ribbon tuned to the idler</b> cuts it 38× at the load: a magnetostrictive resonance is an impedance peak "
        "of the core branch, a trap at that frequency. To help the idler, tune it elsewhere.",
    ])

    r.figure(os.path.join(HERE, "bench.png"),
             "Figure 1. Fig. 2G bench. Top: load voltage with air-core bifilar pancakes and with a biased Metglas "
             "core. Middle: what is left after subtracting the pump (with harmonics) and the signal — round-off for "
             "air, a few mV of mixing products for the biased core. Bottom: load spectrum; only the Metglas core "
             "makes new lines (fp−fs, 2fp−fs, fp+fs, …). The ~10<super>−6</super> V floor of the Metglas spectrum "
             "is numerical, from the core curve's piecewise-linear table.", max_height_mm=150)
    r.h1("A spark-gap Tesla coil with Barrett's floating terminal")
    r.p("Secondary: the 400-turn solenoid of <i>pulsed_bifilar</i>, L 4.85 mH, 15 pF top-load. Primary (declared): "
        "10 kV DC through 10 kΩ charges 20 nF; a rotary gap fires 1000 times a second; k 0.15. A 1 V signal from "
        "50 Ω enters in series with the top-load at 0.9 of the secondary's frequency. The gap is a timed switch "
        "(0.2 Ω for 100 µs per firing), so the coil is a linear, periodically switched circuit. The solver "
        "<i>tcsim.go</i> steps it by the trapezoidal rule; TestTeslaTwin checks it against kronspice with the gap "
        "held closed (top-load voltage within 0.1 %, ledger to 10<super>−9</super>). Runs with the signal at phase "
        "0 and π separate what is linear in the signal; subtracting the signal-alone run leaves what the switching "
        "does to it.")
    r.table([
        ["", "ordinary solenoid", "Tesla bifilar solenoid"],
        ["self-capacitance; ring frequency", "3.84 pF; 526.4 kHz", "3.01×10<super>3</super> pF; 41.5 kHz"],
        ["primary Lp; signal frequency", "4.6 µH; 473.5 kHz", "735.1 µH; 37.5 kHz"],
        ["top-load peak with sparks", "289 kV", "22 kV"],
        ["signal alone at the top-load", "3.38 V", "0.022 V"],
        ["sparks × signal, rms", "1.02 V (0.43 × signal)", "0.0103 V (0.66 × signal)"],
        ["sidebands fs + k·1 kHz, k = −3…3 [V]", "2.16e-01, 2.45e-01, 2.68e-01, 2.85e-01, 2.92e-01, 2.91e-01, 2.80e-01",
         "3.22e-04, 4.13e-04, 5.29e-04, 6.86e-04, 9.20e-04, 1.34e-03, 2.50e-03"],
        ["ledger residual", "−4.7×10<super>−10</super> J of 43.1 J", "6.2×10<super>−11</super> J of 49.5 J"],
        ["gap heat; firings", "7.77 J; 13 in 13 ms", "0.276 J; 20 in 20 ms"],
    ], [0.3, 0.35, 0.35])
    r.figure(os.path.join(HERE, "tesla_coil.png"),
             "Figure 2. Spark-gap Tesla coil, ordinary secondary (left) and Tesla bifilar secondary (right). Top: "
             "top-load voltage with sparks and signal over three firings. Middle: the signal alone and the part the "
             "switching adds (sparks × signal). Bottom: spectra around the signal; the switching puts sidebands every "
             "1 kHz (the break rate), flat for the ordinary coil and rising toward the bifilar coil's 41.5 kHz "
             "resonance.", max_height_mm=150)
    r.p("This is real mixing, made by a time-varying element, of the kind Barrett's own Tesla coils would have shown. "
        "With the full arc model (strikes at thresholds, re-strikes while the channel is hot) the ordinary coil is "
        "hypersensitive: the phase-flip difference is 4.67×10<super>3</super> V rms (1957× the signal at the top), "
        "with 1020 strikes in 13 ms. The 1 V signal flips individual re-strikes in a ring-down of hundreds of kV; a "
        "real gap's random jitter (ns–µs) swamps that, so it is not usable mixing.")

    r.h1("Ledger and validation")
    r.bullets([
        "Fig. 2 ledgers 7×10<super>−9</super> to 1×10<super>−7</super>; transient equals AC superposition to "
        "≤ 1.9×10<super>−4</super>.",
        "Fig. 3 ledger residual 1.0×10<super>−10</super> J of 150 mJ stored.",
        "Tesla-coil solver ledgers −4.7×10<super>−10</super> J of 43.1 J and 6.2×10<super>−11</super> J of 49.5 J; "
        "twinned with kronspice to 0.1 % (TestTeslaTwin).",
        "Metglas runs close to 2×10<super>−4</super>–5×10<super>−4</super>; air-core bench runs to "
        "3×10<super>−9</super>–5×10<super>−9</super>.",
        "Spectra use a Hann window: a floating terminal charges slowly when a nonlinear store rectifies (leak time "
        "constant about a second), and the β–top-load ring has not quite died out after warm-up.",
    ])

    r.h1("Limits of the model")
    r.bullets([
        "All component values are declared; the sources give none.",
        "Coils are lumped equivalents (L, wire R, self-capacitance); the Metglas core is a piecewise-linear "
        "flux curve with a core-loss resistor, the ribbon a parallel tank in series with it.",
        "The rotary gap is a timed switch; the full arc model is used only to show the coil's sensitivity.",
        "Claims beyond these (noise rejection, lower transmission loss from “higher symmetry”) are untested here: "
        "the paper gives no measurement and no equation to simulate.",
    ])

    r.h1("Status")
    r.p("Register status: <b>untested</b>. Figs. 2–5 as drawn are linear and obey Kirchhoff's laws; Barrett's mixing "
        "and phase conjugation are real effects of a nonlinear element, with standard scaling and the Manley–Rowe "
        "relations; the SU(2) and order effects are ordinary two-port algebra; Fig. 3 is a working pulse former. "
        "The bench test is fp − fs in the load of an air-core Fig. 2 build.")

    r.h1("Reproducing the results")
    r.code("go run ./examples/barrett_osc                  # about 1.5 minutes\n"
           "go run ./examples/barrett_osc -data DIR        # also the figure data (CSV)\n"
           "python3 examples/barrett_osc/viz/plot_bench.py DIR viz/bench.png\n"
           "python3 examples/barrett_osc/viz/plot_tesla.py DIR viz/tesla_coil.png\n"
           "go test ./examples/barrett_osc")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "barrett_osc_report.pdf"))
