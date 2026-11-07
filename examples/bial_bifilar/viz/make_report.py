# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""PDF report for examples/bial_bifilar (Bi-Al capacitor element with a bifilar coil, M5a).

python examples/bial_bifilar/viz/make_report.py examples/bial_bifilar/bial_bifilar_report.pdf

Numbers come from examples/bial_bifilar/README.md and from the run report it
links, doc/experiments/bial_bifilar_run2.md (seed 1, 5 weighted-sum runs per
variant, population 40, 30 generations). The example has no run_output.txt and
no PNG figures (the equivalent-circuit drawings in doc/experiments are SVG).
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401


def build(out):
    r = Report("examples/bial_bifilar", "Bi–Al capacitor element with a bifilar coil",
               "The project's first experiment: Q at resonance against stored energy per volume, with every "
               "joule traced to its port — CLAUDE.md §3e, milestone M5a", "CLAUDE.md §3e, M5a",
               "run 2, seed 1, 5 weighted-sum runs × population 40 × 30 generations per variant")

    r.h1("Summary")
    r.p("A bismuth–aluminium capacitor element and a bifilar coil (Tesla, US 512,340) are connected in parallel and "
        "left self-driven: the element is the only energy source, and the energy ledger must trace every joule stored "
        "or dissipated to the element's declared port. Evolution trades the tank's quality factor Q at resonance against "
        "the energy stored per volume, for three element variants: dry (Seebeck thermocouple, thermal port), wet "
        "(Al|electrolyte|Bi galvanic cell, chemical port) and thermopile (n Bi–Al pairs in series, thermal port).")
    r.bullets([
        "<b>Dry: microscopic storage.</b> Best energy per volume 1.254 µJ/m³ (analytic ½·ε·(EMF/gap)² bound 1.255 µJ/m³); "
        "highest Q 1,050 at 311.7 kHz with 5.4 pJ/m³. The 7.55 mV Seebeck EMF only charges the plates.",
        "<b>Wet: up to 89.3 J/m³,</b> mostly ½·L·I² from a 4.079 A steady cell current that a real aluminium electrode "
        "would not sustain (passivation): an upper bound until a cell is measured.",
        "<b>Thermopile: 31.25 J/m³ with 5,000 pairs</b> (EMF 37.75 V, EMF² scaling), but holding ΔT = 100 K across the legs "
        "conducts 47.83 W of heat to store 1.4 mJ.",
        "<b>Bifilar against ordinary:</b> the bifilar winding raises self-capacitance by about N²/4 and lowers the resonance "
        "(on 13 of 59 dry designs by more than 0.1 %), but leaves the stored energy unchanged and raises Q on no dry design.",
        "<b>The ledger closes:</b> worst per-point residual 1.6×10<super>−14</super> of peak power, worst integration mismatch "
        "7.9×10<super>−4</super> of the energy taken in, model-supplied energy 0 J.",
        "<b>The coil model is superseded for prediction</b> (Knight 2016): the Q and resonance values must be re-run with "
        "the helical-line model before they are relied on.",
    ])

    r.h1("Background")
    r.p("Tesla's US 512,340 winds two conductors side by side so that adjacent turns differ by half the coil voltage "
        "rather than V/N, which raises the coil's built-in capacitance; at 1,000 turns the idealized turn-to-turn energy "
        "ratio is (N/2)² = 250,000 (reproduced by TestTeslaComparison). This experiment asks what such a coil does when the "
        "tank is charged only by a Bi–Al element, and keeps the Tesla comparison separate from the realistic results, as "
        "§3e requires.")
    r.p("The element is also the bench for the temporal-field hypothesis <b>H4</b> (Jack 2003, em/temporal.go), which "
        "predicts, under rapid heating with an imposed anisotropy, an EMF from the time variation of temperature set by the "
        "Thomson coefficients and conductivities of Bi and Al. <b>Standard prediction:</b> the Seebeck EMF from the "
        "instantaneous temperature difference only, which is what this example computes. The difference is the measurement "
        "target; H4 itself is not modelled here and remains <i>untested</i>.")

    r.h1("The model")
    r.table([
        ["variant", "element model", "energy enters through"],
        ["dry", "Bi and Al plates with a solid dielectric; the Bi–Al leads form a thermocouple, EMF "
         "(S<sub>Al</sub> − S<sub>Bi</sub>)·ΔT = 75.5 µV/K·ΔT in series with the plate capacitance", "thermal port"],
        ["thermopile", "the dry element charged by n Bi–Al pairs in series (n·75.5 µV/K·ΔT); genes add n (1–5000) and each "
         "leg's diameter and length; legs add series resistance, volume and a heat leak n·(k<sub>Bi</sub> + k<sub>Al</sub>)·A/ℓ·ΔT; "
         "legs shorter than 5 diameters are rejected", "thermal port"],
        ["wet", "Al|electrolyte|Bi galvanic cell (≤ 1.970 V), electrolyte resistance, Randles electrodes (double layer || "
         "charge transfer), geometric plate capacitance; plate capacity by Faraday's law", "chemical port"],
    ], [0.13, 0.67, 0.2])
    r.bullets([
        "<b>Coil</b> (em/coil.go): Wheeler inductance; turn-to-turn capacitance from the two-wire formula with V/N per gap "
        "(ordinary) or V/2 (bifilar); skin effect at resonance; wire longer than λ/10 at resonance is rejected as not quasi-static.",
        "<b>Circuit:</b> the coil as series R and L with its self-capacitance across its terminals, in parallel with the "
        "element, between node a and ground. A 1 A AC probe current into a makes v(a) the tank impedance. For the switch-on "
        "transient the EMF steps on at t = 0 (idealized).",
        "<b>Resonance:</b> where the impedance phase falls through zero from inductive to capacitive; Q = (f<sub>0</sub>/2)·|dφ/df| "
        "from the phase slope, f<sub>0</sub> fixed by bisection to 10<super>−12</super>. A tank with no crossing is booked as no-measurement.",
        "<b>Genes:</b> turns 10–2000, coil diameter 5–150 mm, wire diameter 0.1–2 mm, pitch 1.05–3 wire diameters, insulation "
        "ε 2–6, plate side 5–150 mm, gap 10 µm–2 mm (wet 50 µm–5 mm), dielectric ε 2–100 (wet: set by the electrolyte); "
        "thermopile adds pairs 1–5000, leg diameter 50 µm–2 mm, leg length 1–50 mm.",
        "<b>Objectives:</b> Q at resonance, and energy per volume as the error log<sub>10</sub>(1 MJ/m³ / E<sub>v</sub>). Five "
        "weighted-sum runs per variant trace the trade-off; non-dominated designs are reported. Every front design is "
        "re-simulated with an ordinary coil of the same geometry.",
    ])
    r.h2("Declared assumptions")
    r.table([
        ["assumption", "value", "source"],
        ["Bi–Al Seebeck coefficient", "75.5 µV/K (S<sub>Al</sub> +3.5, S<sub>Bi</sub> −72 µV/K vs Pt)", "Electronics Cooling (Nov 2006)"],
        ["dry: plates above coil connections", "ΔT = 100 K, EMF 7.55 mV", "assumption"],
        ["wet: cell EMF", "1.97 V (Al<super>3+</super>/Al −1.662 V, Bi<super>3+</super>/Bi +0.308 V), the thermodynamic upper bound", "CRC electrochemical series"],
        ["wet: electrolyte", "ρ = 0.21 Ω·m (sea water, S = 35, 20 °C), ε = 80", "UNESCO 1981"],
        ["wet: double layer, charge transfer", "20 µF/cm², 10 Ω·cm² per electrode", "Bard &amp; Faulkner ch. 13 (typical)"],
        ["plates", "1 mm Bi + 1 mm Al", "assumption"],
        ["coil", "copper; Wheeler L; two-wire turn-to-turn C; skin effect at f<sub>0</sub>", "Wheeler 1928; Tesla; Terman 1943"],
        ["spectral line", "Cu I 324.754 nm (air), 922.9 THz", "NIST"],
    ], [0.3, 0.45, 0.25])

    r.h1("Results (run 2)")
    r.p("Each variant: 5,710 evaluations. Dry and wet are bit-identical to run 1 (same seed); run 2 adds the thermopile.")
    r.table([
        ["variant", "penalized", "front", "highest Q", "highest energy per volume"],
        ["dry", "271", "59 designs", "39 turns, 150 mm, 311.7 kHz, Q 1,050, 5.4 pJ/m³",
         "15 turns, 5 mm, 10 µm gap, 78.9 kHz, Q 0.75, 1.254 µJ/m³"],
        ["wet", "599", "38 designs", "23 turns, 141 mm, 11.16 kHz, Q 6.4, 4.7 mJ/m³", "10 turns, 5.0 mm, 15.63 kHz, Q 2.4, 89.3 J/m³"],
        ["thermopile", "423", "19 designs", "3 pairs, 11 turns, 21.15 MHz, Q 1,250, 86 pJ/m³",
         "5,000 pairs, 19 turns, 27.75 MHz, Q 43.3, 31.25 J/m³"],
    ], [0.12, 0.1, 0.11, 0.33, 0.34])
    r.small("Penalty reasons. Dry: coil too short for Wheeler 171, wire too thick 66, gap not small 34, no resonance 8, "
            "wire &gt; λ/10 3. Wet: wire too thick 262, no resonance 216, coil too short 85, gap 46. Thermopile: legs shorter "
            "than 5 diameters 168, coil too short 147, wire &gt; λ/10 59, gap 44, wire too thick 11.")
    r.h2("Extreme designs in detail (doc/experiments/bial_bifilar_run2.md)")
    r.table([
        ["design", "L", "self-C bifilar (ordinary)", "steady state", "switch-on transient ledger"],
        ["dry, highest Q", "225.5 µH", "867.8 pF (2.282 pF; ratio 380.2)", "no current; plates charge in 121.7 ps",
         "in 613.7 pJ, dissipated 613.7 pJ, stored 8.216 fJ; residual 7.4×10<super>−16</super>"],
        ["dry, highest energy", "1.264 µH", "3.449 pF (61.32 fF; 56.25)", "no current; τ 987.3 ns",
         "in 1.645 nJ, dissipated 1.591 nJ, stored 55.69 pJ; residual 1.6×10<super>−14</super>"],
        ["wet, highest Q", "82.66 µH", "169.1 pF (1.279 pF; 132.2)", "24.48 mA, 48.22 mW; plate lasts 8h13m",
         "in 118.2 µJ, dissipated 91.19 µJ, stored 4.819 µJ; residual 8.2×10<super>−16</super>"],
        ["wet, highest energy", "238.1 nH", "14.56 pF (582.4 fF; 25)", "4.079 A, 8.035 W; plate lasts 8h18m",
         "in 13.1 mJ, dissipated 11.07 mJ, stored 792.1 µJ; residual 3×10<super>−16</super>"],
        ["thermopile, highest Q", "1.469 µH", "34.66 pF (1.146 pF; 30.25)", "no current; heat leak 5.463 W",
         "in 160.8 fJ, dissipated 159.8 fJ, stored 996.5 aJ; residual 2.5×10<super>−16</super>"],
        ["thermopile, highest energy", "1.921 µH", "17.12 pF (189.7 fF; 90.25)", "no current; τ 33.02 ms; heat leak 47.83 W",
         "in 395.7 nJ, dissipated 395.7 nJ, stored 32.86 pJ; residual 6.8×10<super>−16</super>"],
    ], [0.17, 0.1, 0.21, 0.22, 0.30])
    r.small("Each ratio of bifilar to ordinary self-capacitance equals N²/4 at that turn count. Transients are shorter than "
            "the plate charging time where stated, so they end with the plates partly charged. For the two thermopile "
            "extremes the ordinary-coil re-simulation is outside the lumped model (wire longer than λ/10 at the ordinary "
            "coil's higher resonance), which the report prints in place of the comparison.")

    r.h2("What the fronts say")
    r.bullets([
        "<b>Dry.</b> The only store is the field in the plate gap, ½·ε·(EMF/gap)² diluted by the plates: 1.255 µJ/m³ by that "
        "formula against 1.254 µJ/m³ simulated for the best design. With a millivolt Seebeck EMF the energy per volume is "
        "necessarily microscopic. At the steady state no current flows (the plates block DC).",
        "<b>Wet.</b> The store is mostly ½·L·I² from the steady cell current, set by the assumed charge-transfer resistance: "
        "4.079 A (96.9 mA/cm²) in the best design. The plate's charge lasts the same few hours (8 h 18 min) for every design "
        "because capacity and current both scale with plate area.",
        "<b>Thermopile.</b> 5,000 pairs raise the EMF to 37.75 V and the stored energy per volume from 1.25 µJ/m³ to 31 J/m³, "
        "the EMF² scaling (5,000² = 2.5×10<super>7</super>). Holding ΔT = 100 K conducts 47.83 W of heat (outside the electrical "
        "ledger) to store 1.4 mJ, and the plates charge through 16.86 kΩ with a 33 ms time constant.",
        "<b>Bifilar against ordinary.</b> At DC no voltage stands across the turn-to-turn capacitance, so the stored energy "
        "is the same for both windings. Lower resonance by more than 0.1 %: dry 13 of 59, wet 0 of 38, thermopile 4 of 19. "
        "Q raised by more than 0.1 %: dry 0, wet 0, thermopile 2.",
        "<b>Edges of the gene ranges.</b> The highest-energy optimum of every variant sits at the edge of some allowed ranges "
        "(e.g. dry: coil and wire diameter, insulation ε, gap, dielectric ε), so it is set by the chosen limits.",
    ])

    r.h1("Ledger and validation")
    r.bullets([
        "No evaluated design failed the operating-point ledger check. Over the six switch-on transients the worst per-point "
        "residual is 1.6×10<super>−14</super> of peak power, the worst integration mismatch 7.9×10<super>−4</super> of the energy "
        "taken in, and model-supplied energy is 0 J: every joule comes from the element's declared port.",
        "Tests: TestTeslaComparison (250,000 at 1,000 turns), TestResonanceOfKnownTank (f<sub>0</sub>, Q of a parallel RLC), "
        "TestDrySteadyState, TestWetSteadyState, TestThermopileScaling (EMF × n, stored × n²), TestOutsideValidityIsRejected "
        "(λ/10), TestQuickExperiment, TestExport.",
        "Resonances are also reported as the nearest whole-number subharmonic of the Cu I line; adjacent subharmonics are "
        "about 1/n apart with n ~ 10<super>7</super>–10<super>11</super>, so this is reported, not a constraint.",
        "The neutrino channel is listed as not modelled in this run.",
    ])

    r.h1("Limits of the model")
    r.bullets([
        "<b>Superseded coil model.</b> Knight (2016, §6, §13 point 8) shows turn-to-turn capacitance theories have no "
        "predictive power for single-layer solenoids (Medhurst 1947 already found the model underestimates an ordinary coil's "
        "self-capacitance). The V/N model is kept to reproduce Tesla's idealized comparison; ordinary-coil self-capacitance "
        "and self-resonance belong to em/helix.go, and the strand coupling to the coupled-line model. The Q and resonance "
        "values above must be re-run before they are relied on (§3e).",
        "The wet results are upper bounds: the charge-transfer resistance is assumed, and aluminium passivates well below "
        "the current densities found.",
        "ΔT = 100 K, the 1.97 V cell EMF and the plate thicknesses are assumptions; the EMF steps on instantly, while a real "
        "thermal or chemical source turns on over milliseconds or longer.",
    ])

    r.page_break()
    r.h1("Reproducing the results")
    r.code("go run ./examples/bial_bifilar                 # about 1 min: report in runs/<time>-bial_bifilar/\n"
           "go run ./examples/bial_bifilar -quick          # a few seconds\n"
           "go run ./examples/bial_bifilar -deltaT 300 -emf-factor 0.6 -line Al\n"
           "go run ./examples/bial_bifilar -netlists doc/experiments   # netlists + schematics of the best designs\n"
           "go test ./examples/bial_bifilar\n"
           "python examples/bial_bifilar/viz/make_report.py examples/bial_bifilar/bial_bifilar_report.pdf")
    r.p("The kept runs, their fronts (front_dry_run1.csv, front_wet_run1.csv, front_thermopile_run2.csv) and the "
        "netlists, LTspice files and equivalent-circuit drawings of the best designs are in doc/experiments/.")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "bial_bifilar_report.pdf"))
