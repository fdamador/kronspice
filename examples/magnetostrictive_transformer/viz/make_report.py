# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""H37 report: the magnetostrictive transformer with two ultrasonic transducers, as a PDF.

    python examples/magnetostrictive_transformer/viz/make_report.py \
        examples/magnetostrictive_transformer/magnetostrictive_transformer_report.pdf

Numbers come from examples/magnetostrictive_transformer/run_output.txt (the
program's output) and README.md (model, materials, bench notes); figures are
viz/setup.png and viz/lines.png. Layout from internal/pdfreport/kronreport.py.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401


def build(out):
    r = Report("examples/magnetostrictive_transformer",
               "Magnetostrictive transformer with two ultrasonic transducers",
               "Does a magnetostrictive core carrying ultrasound mix it with the electrical drive, and does the "
               "swap test show it? — register entry H37",
               "register entry H37")

    r.h1("Summary")
    r.p("A transformer is wound on a magnetostrictive core — Terfenol-D or Galfenol, shaped as a rod, a ring or a "
        "helical spring — on a DC bias, and two ultrasonic transducers set up standing waves in the core. The primary "
        "is driven at 5 kHz (f<sub>e</sub>) and the transducers at the core's modes (f<sub>a</sub>, f<sub>b</sub>). "
        "The question (H37, proposed after Barrett's H36) is whether the core mixes the sound with the electrical drive, "
        "and whether driving the secondary instead (the swap test) shows the nonlinearity. The standard prediction is "
        "computed with an energy-consistent magnetoelastic model.")
    r.bullets([
        "<b>Drive × sound products appear on the secondary:</b> f<sub>a</sub> ± f<sub>e</sub> and f<sub>b</sub> ± "
        "f<sub>e</sub> up to 2.9 mA (Terfenol-D ring, f<sub>b</sub> + f<sub>e</sub>), and only with both drive and "
        "sound on.",
        "<b>Sound × sound products are about 100× weaker</b> (f<sub>a</sub> + f<sub>b</sub>, f<sub>b</sub> − "
        "f<sub>a</sub>: 0.017 mA at most).",
        "<b>Reciprocal at small drive:</b> the transformer is reciprocal to 2.1×10<super>−10</super> or better, the "
        "check that the magnetoelastic coupling is energy-consistent.",
        "<b>Not reciprocal at full drive for the products:</b> swapped, the Galfenol ring's summed products fall from "
        "3.04 to 1.71 mA (×0.56) and the rod's from 0.388 to 0.290 mA (×0.75), while the transformer line stays "
        "nearly unchanged.",
        "<b>The Galfenol helix reads the ultrasound best</b> (20.0 mA at f<sub>a</sub>); the Terfenol-composite helix "
        "is weakest (H/a = 0.02). The spring shape does not amplify.",
        "<b>Ledger residuals</b> 8×10<super>−5</super> to 1×10<super>−2</super> (first-order time stepping).",
    ])

    r.h1("What is being tested")
    r.p("<b>Hypothesis H37</b> (proposed 2026, after Barrett H36): a Terfenol-D or Galfenol core (rod, ring, helical "
        "spring) carrying standing waves from transducers at its ends mixes them with the electrical drive, and the "
        "swap test (drive the secondary, watch the primary) shows the nonlinearity. It is a mechanical relative of "
        "Barrett's oscillator-shuttle mixer (examples/barrett_osc).")
    r.p("<b>Standard prediction</b> (computed here): an anhysteretic magnetoelastic material mixes drive and sound "
        "through the curvature of its magnetization and strain curves. Products fall to zero at zero bias, where "
        "λ′ = 0, and vanish when either the drive or the sound is off. The linear transfer must be reciprocal; the "
        "products need not be.")

    r.figure(os.path.join(HERE, "setup.png"),
             "Figure 1. The set-up. Top: the three core shapes — (1) Terfenol-D rod, 50 mm × Ø10 mm, with primary (blue) "
             "and secondary (red) interleaved, a return yoke and prestress clamp, and transducers A and B at the ends "
             "(modes 1 and 3, 18 and 54 kHz); (2) ring, Ø40 mm, wound all the way round, radial transducers driving the "
             "breathing mode; (3) helix spring, 20 turns, Ø20 mm coil of Ø2 mm wire (1.26 m), windings along the wire. "
             "The DC bias winding is dashed. Bottom: the circuit — 20 V, 5 kHz source with 2 Ω on the 200-turn primary, "
             "400-turn secondary into the 20 Ω load, DC bias coil, and two generators driving the transducers; the "
             "spectrum analyser watches the secondary.")
    r.h1("The model")
    r.h2("Material")
    r.p("An anhysteretic magnetization with quadratic magnetostriction, derived from one potential:")
    r.code("M = Ms·tanh(Ht/a),   λ = λm·tanh²(Ht/a),   Ht = H + Hbias\n"
           "σ = E·(ε − λ),       B = μ0·(Ht + M) + σ·λ′(Ht)\n"
           "Φ(Ht, ε) = ½E(ε−λ)² − μ0·(Ht²/2 + Ms·a·ln cosh(Ht/a)),   σ = ∂Φ/∂ε,  B = −∂Φ/∂Ht")
    r.p("Here M is the magnetization, M<sub>s</sub> its saturation value, a the field scale of the curve, λ the "
        "magnetostrictive strain, σ the stress, ε the strain and E Young's modulus. Because both laws come from one "
        "potential, the magnetic and mechanical sides exchange energy exactly; the stored energy is Φ + H<sub>t</sub>·B. "
        "Hysteresis and eddy currents are not modelled. The bias puts every material at H<sub>b</sub> = 0.66·a, where "
        "λ′ = dλ/dH is largest.")
    r.table([
        ["material (declared)", "d = λ′ at bias", "μ<sub>r</sub>", "k", "E", "c", "Q", "bias", "sources"],
        ["Terfenol-D", "13.3 nm/A", "7.1", "0.77", "30 GPa", "1801 m/s", "10", "57 kA/m",
         "Etrema: μ<sub>r</sub> 3–10, E 25–35 GPa, 800–1200 ppm"],
        ["Terfenol-D/epoxy composite (51 %)", "5.4 nm/A", "3.7", "0.31", "15 GPa", "1698 m/s", "15", "66 kA/m",
         "measured d<sub>33</sub> ≈ 5.5 nm/A, k<sub>33</sub> ≈ 0.34"],
        ["Galfenol (Fe<sub>81</sub>Ga<sub>19</sub>)", "20.9 nm/A", "87.1", "0.49", "60 GPa", "2774 m/s", "100",
         "6 kA/m", "Etrema: d<sub>33</sub> 15–30 nm/A, μ<sub>r</sub> 75–100, B<sub>s</sub> 1.5 T, 7800 kg/m³"],
    ], [0.17, 0.1, 0.06, 0.06, 0.08, 0.09, 0.06, 0.09, 0.29])
    r.small("k is the magnetomechanical coupling factor, c the sound speed, Q the mechanical quality factor at the first "
            "tone (Kelvin–Voigt damping).")
    r.h2("Cores, windings and numerics")
    r.bullets([
        "<b>Rod:</b> Terfenol-D, 50 mm × Ø10 mm; flux closes through a return yoke; a prestress clamp holds it. "
        "Modelled as a free–free line of segments, transducers at the two ends.",
        "<b>Ring:</b> Ø40 mm, 5 × 5 mm; Terfenol as bonded segments or machined Galfenol. Modelled by its breathing "
        "mode, driven by radial transducers.",
        "<b>Helix:</b> 20 turns, Ø20 mm, of Ø2 mm wire (1.26 m), Galfenol or moulded Terfenol composite; a free–free "
        "line along the wire.",
        "<b>Windings:</b> uniform, so H is uniform along the core and the coils link its mean flux. Primary 200 turns, "
        "secondary 400 turns, each with 5 % leakage; ports 2 Ω (primary) and 20 Ω (load).",
        "<b>Drive:</b> 20 V at 5 kHz; transducers 5 N each. Tones: the rod's modes 1 and 3; the ring's breathing mode "
        "and 1.37× it; the helix wire's odd modes near 20 and 40 kHz. Even free–free modes and a ring's "
        "circumferential modes have no net strain, so uniform coils cannot feel them.",
        "<b>Numerics:</b> coil currents are stiff (leakage time constants can be nanoseconds), so they step by "
        "backward Euler and the mechanics by symplectic Euler. With a 4× finer step the lines move by 1–2 % and the "
        "ledger residual falls 4×, as a first-order scheme should.",
    ])

    r.h1("Results")
    r.h2("Forward run: drive the primary, watch the secondary")
    r.p("All lines are load current in mA (from run_output.txt). H/a is the peak field swing relative to the curve's "
        "field scale.")
    r.table([
        ["core", "tones f<sub>a</sub>, f<sub>b</sub> [kHz]", "H peak/a", "f<sub>e</sub>", "f<sub>a</sub>",
         "f<sub>b</sub>", "ledger"],
        ["Terfenol-D rod", "18.0, 54.0", "0.07", "1397.810", "3.1619", "1.7808", "1e-03"],
        ["Terfenol-D ring", "14.3, 19.6", "0.11", "1199.246", "6.9186", "13.6145", "1e-02"],
        ["Galfenol ring", "22.1, 30.2", "0.18", "1334.504", "3.4827", "1.9444", "5e-03"],
        ["Terfenol-D composite helix", "20.9, 39.9", "0.02", "13.139", "0.7920", "0.7493", "8e-05"],
        ["Galfenol helix", "21.0, 40.8", "0.17", "271.884", "20.0115", "4.3118", "2e-03"],
    ], [0.28, 0.16, 0.1, 0.12, 0.12, 0.12, 0.1])
    r.p("The mixing products on the same runs (mA):")
    r.table([
        ["core", "f<sub>a</sub>+f<sub>e</sub>", "f<sub>a</sub>−f<sub>e</sub>", "f<sub>b</sub>+f<sub>e</sub>",
         "f<sub>b</sub>−f<sub>e</sub>", "f<sub>a</sub>+f<sub>b</sub>", "f<sub>b</sub>−f<sub>a</sub>"],
        ["Terfenol-D rod", "2.2e-01", "3.4e-02", "7.8e-02", "5.7e-02", "2.1e-03", "1.4e-03"],
        ["Terfenol-D ring", "7.7e-01", "2.1e-01", "2.9e+00", "1.2e-01", "1.7e-02", "5.1e-03"],
        ["Galfenol ring", "5.6e-01", "1.2e-01", "1.5e-01", "2.2e+00", "4.1e-03", "2.2e-03"],
        ["Terfenol-D composite helix", "2.8e-04", "1.7e-04", "4.6e-04", "3.6e-04", "3.4e-07", "1.1e-07"],
        ["Galfenol helix", "1.7e+00", "1.3e+00", "4.4e-01", "4.1e-01", "5.3e-03", "3.0e-03"],
    ], [0.28, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12])
    r.h2("What the secondary sees")
    r.bullets([
        "<b>Transformer line (f<sub>e</sub>):</b> 1.2–1.4 A for the rod and rings. The helix's long, thin magnetic "
        "path gives a far smaller inductance and field (13.139 mA composite, 271.884 mA Galfenol).",
        "<b>Ultrasound itself (f<sub>a</sub>, f<sub>b</sub>):</b> the coil reads it through the Villari effect, "
        "1–20 mA.",
        "<b>Products of drive and sound (f<sub>a</sub> ± f<sub>e</sub>, f<sub>b</sub> ± f<sub>e</sub>):</b> the "
        "nonlinear part, up to 2.9 mA. They appear only with both drive and sound on (TestRodAndRing), and are largest "
        "where the field swings furthest on the curve (H/a) and the coupling is strong.",
        "<b>Sound × sound products:</b> about 100× weaker. In this model the magnetization is linear in strain, so these "
        "products need the magnetic field as a go-between.",
        "<b>Rings</b> (closed flux path, no yoke) give the strongest products per amp. The Terfenol ring's "
        "f<sub>b</sub> line (13.6 mA, off resonance) beats its f<sub>a</sub> line (6.9 mA, at the breathing "
        "frequency): with the coils loaded the ring's resonance moves up by the ΔE effect, 1/√(1 − k²) = 1.57 for "
        "k = 0.77.",
        "<b>The spring shape does not amplify.</b> Uniform strain along the wire changes the coil's diameter, not its "
        "length. A spring's amplified motion comes from twisting the wire (the Wiedemann effect), which needs a current "
        "through the helix itself; that is not modelled here and is the next step.",
        "<b>Ring standing waves:</b> on a ring wound all the way round, circumferential standing waves cannot reach the "
        "coils; only the breathing mode couples. Sector windings would change that.",
    ])
    r.figure(os.path.join(HERE, "lines.png"),
             "Figure 2. Top: forward run, secondary (load) current per line on a log scale for the five cores — the "
             "transformer line f<sub>e</sub>, the ultrasound read by the coil (f<sub>a</sub>, f<sub>b</sub>) and, in "
             "the shaded band, the products of drive and sound and of sound with sound. Bottom: swap test, drive × "
             "ultrasound products summed, forward (blue) and swapped (red), with the swapped/forward ratio above each pair.")
    r.h2("Swap test: drive the secondary, watch the primary")
    r.p("The same source is moved to the secondary port, in series with its 20 Ω. Reciprocity is checked at small "
        "drive; the line currents are at full drive.")
    r.table([
        ["core", "reciprocity (small drive)", "f<sub>e</sub> fwd → rev [mA]",
         "products f<sub>a</sub>±f<sub>e</sub>, f<sub>b</sub>±f<sub>e</sub> summed: fwd → rev [mA]", "ratio"],
        ["Terfenol-D rod", "3.8e-12", "1397.810 → 1397.814", "3.88e-01 → 2.90e-01", "×0.75"],
        ["Terfenol-D ring", "4.8e-12", "1199.246 → 1199.771", "4.01e+00 → 4.06e+00", "×1.01"],
        ["Galfenol ring", "1.2e-13", "1334.504 → 1334.521", "3.04e+00 → 1.71e+00", "×0.56"],
        ["Terfenol-D composite helix", "2.0e-11", "13.139 → 13.138", "1.27e-03 → 1.27e-03", "×1.00"],
        ["Galfenol helix", "2.1e-10", "271.884 → 268.662", "3.79e+00 → 3.95e+00", "×1.04"],
    ], [0.22, 0.17, 0.2, 0.31, 0.1])
    r.small("Ratios as printed in the README and on Figure 2.")
    r.p("At small drive the transformer is reciprocal to 1e-10 or better. At full drive the transformer line stays "
        "almost reciprocal, but the products need not be: the Galfenol ring gives 44 % less going backwards, the rod "
        "25 % less. The two directions put the core at different points of its curve (the windings' current, and so "
        "the field swing, differ) — a clean sign of the nonlinearity, which is what the swap test looks for.")

    r.h1("Ledger and validation")
    r.bullets([
        "Energy ledger residuals per core: rod 1e-3, Terfenol ring 1e-2, Galfenol ring 5e-3, composite helix 8e-5, "
        "Galfenol helix 2e-3 (first-order scheme; the residual falls 4× with a 4× finer step).",
        "Lines converged to 1–2 % under a 4× finer step.",
        "Small-signal reciprocity 1.2e-13 to 2.1e-10: the magnetoelastic coupling is energy-consistent.",
        "TestRodAndRing: the drive × sound products appear only with both the drive and the sound on.",
    ])

    r.h1("Limits of the model")
    r.bullets([
        "Anhysteretic: no hysteresis loss and no eddy currents. Eddy currents are why the helix uses Terfenol composite; "
        "a solid Terfenol core would need laminating.",
        "Uniform windings only, so only modes with net strain couple; sector windings are not modelled.",
        "No Wiedemann (twist) effect: a current through the helix wire itself is not modelled.",
        "Material values are declared from data sheets and the cited measurement, not fitted to a bench core.",
        "First-order time stepping (ledger 1e-4 to 1e-2).",
    ])

    r.h1("Bench notes")
    r.bullets([
        "<b>Mounting:</b> Terfenol-D is brittle (28 MPa in tension). Pre-compress it (a few to 10 MPa) and bond the "
        "transducers under that clamp.",
        "<b>Galfenol</b> can be machined into a ring or bent into a spring.",
        "<b>Bias:</b> set it with a DC coil or permanent magnets. The products vanish at zero bias, where λ′ = 0.",
        "<b>Spectrum analyser:</b> watch f<sub>a</sub> ± f<sub>e</sub> and f<sub>b</sub> ± f<sub>e</sub> on the "
        "secondary. Repeat with the transducers off (products gone) and with the drive off (only f<sub>a</sub>, "
        "f<sub>b</sub>).",
    ])

    r.h1("Status")
    r.p("<b>H37: bounded.</b> The standard prediction is computed: drive × sound products of up to 2.9 mA on the "
        "secondary, absent without sound; sound × sound 100× weaker; reciprocal at small drive, products not reciprocal "
        "at full drive (Galfenol ring ×0.56, rod ×0.75 when swapped). The spring shape does not amplify; amplification "
        "would need the Wiedemann twist with current through the helix, not modelled. A bench build compares its "
        "spectrum against these lines.")

    r.h1("Reproducing the results")
    r.code("go run ./examples/magnetostrictive_transformer            # about 2 s\n"
           "go run ./examples/magnetostrictive_transformer -data DIR  # also the figure data\n"
           "python3 examples/magnetostrictive_transformer/viz/plot_lines.py DIR viz/lines.png\n"
           "python3 examples/magnetostrictive_transformer/viz/setup.py viz/setup.png\n"
           "go test ./examples/magnetostrictive_transformer")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else
          os.path.join(HERE, "..", "magnetostrictive_transformer_report.pdf"))
