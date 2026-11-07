# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Bubble collapse (CLAUDE.md §3h, register entry H30) as a PDF report.

    python examples/bubble_collapse/viz/make_report.py [out.pdf]

Numbers come from examples/bubble_collapse/run_output.txt (the recorded
`go run ./examples/bubble_collapse`), README.md and the text of the interactive
report viz/report.html (left unchanged; its charts are drawn in the browser from
embedded data, so they have no image files and are not reproduced here).
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401


def main(out):
    r = Report("examples/bubble_collapse", "Bubble collapse",
               "Is a sonoluminescing bubble a large-scale picture of the knot's implosion, and does its light "
               "carry the knot's topology? — CLAUDE.md §3h, register entry H30",
               "CLAUDE.md §3h, register entry H30",
               "argon R0 = 4.5 µm, 1.3 atm at 25.02 kHz; recorded run: run_output.txt")

    r.h1("Summary")
    r.p("A single sonoluminescing argon bubble is put on the matched water flask of <i>acoustic_drive</i> as an "
        "equivalent circuit (the bubble returns to scope only as a cited equivalent-circuit model, CLAUDE.md §1 "
        "non-goals), and compared with the implosion of knotted light in <i>knotted_light</i> (H29, H30). The "
        "user's question (2026-11-05): is the energy at the centre of a collapsing bubble a large-scale picture of "
        "the knot's implosion, and does its light carry the knot's topology? The run follows every joule of one "
        "acoustic cycle, loads the flask's drive chain with the bubble, tests reversibility, estimates the flash "
        "and its spectrum against measurement, and measures the flash's circular polarization in the field solver.")
    r.bullets([
        "<b>One cycle, ledger closed.</b> 14.2923 nJ in from the sound field = 1.04039 nJ viscous heat + "
        "13.2519 nJ radiated sound (Δ stored −4.99×10<super>−18</super> J); residual 6.2×10<super>−10</super> "
        "of the input. 91 % of the radiated sound leaves within ±10 ns of the minimum.",
        "<b>Circuit and Keller–Miksis agree on the slow motion</b> (R<sub>max</sub> 32.60 against 32.73 µm, minimum "
        "at phase 0.5170 against 0.5181) and differ at the collapse (wall Mach 0.781 against 0.516).",
        "<b>The bubble is a small load on the flask:</b> 357.66 µW against 6.95727 mW into the flask's mode; the "
        "generator EMF for 1.3 atm rises from 1.8981 V to 1.9240 V.",
        "<b>Reversible only without losses:</b> the lossless bubble returns to |ΔR|/R = 2.5×10<super>−12</super>; "
        "with sound radiation it misses by 0.90.",
        "<b>The flash estimate:</b> 12.7746 fJ, about 6.2×10<super>4</super> photons, FWHM 180 ps, against "
        "measured 10<super>5</super>–10<super>7</super> photons and 60–250 ps; the width matches, the colour and "
        "intensity do not.",
        "<b>No topology in the flash:</b> the hopfion gives P<sub>c</sub> = +0.9996 (exact +0.9999); the thermal "
        "core gives mean +0.0058, rms 0.0184 over six seeds, consistent with zero.",
    ])

    r.h1("What is being tested")
    r.p("<b>The idea (H30 context, user question).</b> The knot's implosion and explosion (H30, "
        "<i>knotted_light</i>) might have a large-scale counterpart in the collapse of a sonoluminescing bubble, "
        "and the bubble's light might carry the knot's topology. <b>Standard prediction.</b> The collapse is a "
        "nonlinear, irreversible, spherically symmetric flow of matter stopped by gas pressure; its light is "
        "thermal (Hilgenfeldt, Grossmann &amp; Lohse, <i>Nature</i> 398, 402, 1999; Brenner, Hilgenfeldt &amp; "
        "Lohse, <i>Rev. Mod. Phys.</i> 74, 425, 2002), with degree of circular polarization P<sub>c</sub> = 0 "
        "within counting noise. A knot of field energy is one helicity and gives P<sub>c</sub> of order 1.")

    r.h1("The model")
    r.h2("The bubble as a circuit")
    r.p("In the acoustic domain the across variable is pressure and the through variable is volume velocity "
        "U = dV/dt. From the far liquid to the gas:")
    r.code("p∞(t) ── [ M(V) || R_rad(V) ] ── R_visc(V) ── wall ── gas p_g(V) − 2σ/R")
    r.table([
        ["Element", "Value", "Source"],
        ["liquid inertance", "M = ρ/(4πR), so ½MU² = 2πρR³Ṙ²", "Rayleigh 1917"],
        ["radiation resistance", "R<sub>rad</sub> = ρc/(4πR²)",
         "pulsating sphere: M || R<sub>rad</sub> = (ρc/4πR²)·jkR/(1+jkR), Kinsler et al. §7.2"],
        ["wall viscosity", "R<sub>visc</sub> = μ/(πR³) (= 4μṘ/R per U)", "Newtonian liquid"],
        ["gas", "argon, adiabatic γ = 5/3, hard core R0/h = 8.86", "Brenner, Hilgenfeldt &amp; Lohse 2002"],
        ["surface", "2σ/R", "IAPWS R1-76(2014) σ = 71.69 mN/m at 300 K"],
    ], [0.2, 0.37, 0.43])
    r.p("M depends on the volume, the circuit's own \"charge\", so the inertance is a Kron inductance on a "
        "coordinate. With momentum Π = M·U<sub>M</sub>, dΠ/dt = Δp + ½(dM/dV)·U·U<sub>M</sub>, where the ½M′ term "
        "is the Christoffel (speed-voltage) term of Kron 1951 (Epilogue eqs. 3–4). Without the two losses this is "
        "the Rayleigh–Plesset equation exactly. R<sub>rad</sub> adds the liquid's compressibility to first order in "
        "Ṙ/c, the job Keller–Miksis does; Keller–Miksis runs beside it as the twin. Integration is Dormand–Prince "
        "5(4) at tolerance 10<super>−10</super>.")
    r.p("Parameters: argon bubble R0 = 4.5 µm in water at 300 K, drive 1.30 atm at 25.02 kHz (the matched flask's "
        "peak); ρ = 996.556 kg/m³, c = 1501.52 m/s (IAPWS-95), μ = 853.8 µPa·s (IAPWS R12-08), "
        "σ = 71.69 mN/m.")
    r.h2("The ledger")
    r.p("It closes by construction (d(½Π²/M)/dt = U<sub>M</sub>·Δp, and U<sub>M</sub>·Δp = U·Δp − "
        "Δp²/R<sub>rad</sub>):")
    r.code("∫ −p_ac·U dt = Δ[½MU_M² + p_g(V−Vh)/(γ−1) + 4πσR² + P0·V] + ∫R_visc·U² dt + ∫Δp²/R_rad dt\n"
           "in from the sound field = Δ stored + viscous heat + sound radiated")

    r.h1("Results")
    r.h2("Part 1: one acoustic cycle (the sixth from rest)")
    r.table([
        ["", "Circuit", "Keller–Miksis"],
        ["R<sub>max</sub>", "32.60 µm", "32.73 µm"],
        ["R<sub>min</sub>", "0.730 µm", "0.769 µm"],
        ["minimum at phase [cycles]", "0.5170", "0.5181"],
        ["max wall Mach |Ṙ|/c", "0.781", "0.516"],
    ], [0.4, 0.3, 0.3])
    r.bullets([
        "Peak gas: 2.35×10<super>9</super> Pa and 14985 K (adiabatic); R<sub>max</sub>/R<sub>min</sub> = 44.7.",
        "Ledger over the cycle: in from the sound field 14.2923 nJ = Δ stored −4.99483×10<super>−18</super> J + "
        "viscous 1.04039 nJ + radiated sound 13.2519 nJ; residual 6.2×10<super>−10</super> of the input. The "
        "collapse energy P0·V<sub>max</sub> is 14.7069 nJ.",
        "91 % of the radiated sound leaves within ±10 ns of the minimum: the collapse turns a 25 kHz input into a "
        "nanosecond outgoing pulse (peak 2.7 kW for a few nanoseconds, per the interactive report).",
        "Energy density goes from 3.86 J/m³ in the sound field to 9.03×10<super>9</super> J/m³ in V<sub>min</sub>, "
        "a concentration of 2.3×10<super>9</super>.",
        "The two models agree on the slow quantities (R<sub>max</sub> to 0.4 %, the minimum's phase to 0.1 %) and "
        "differ at the collapse, where Ṙ/c is no longer small: the circuit is first order in Ṙ/c, while "
        "Keller–Miksis keeps (1 − Ṙ/c) factors that slow the wall more. Neither is exact there (Brenner et al. "
        "discuss Gilmore and full gas dynamics).",
    ])

    r.h2("Part 2: the bubble as a load on the drive chain")
    r.p("The bubble's fundamental admittance Y<sub>1</sub> = −U<sub>1</sub>/P<sub>1</sub> = "
        "4.123×10<super>−14</super> + 2.825×10<super>−14</super>j m³/(Pa·s) is put at node p of the matched water "
        "chain as a resistor and a capacitor. This harmonic balance is exact for power, because only the "
        "fundamental of U does work against the sinusoidal centre pressure: ½Re(Y<sub>1</sub>)P<sub>a</sub>² = "
        "357.66 µW equals the cycle work × f = 357.66 µW.")
    r.table([
        ["", "No bubble", "With bubble"],
        ["peak", "25.025 kHz", "25.024 kHz"],
        ["pressure per volt", "69.40 kPa/V", "68.46 kPa/V"],
        ["generator EMF for 1.3 atm", "1.8981 V", "1.9240 V"],
        ["source power", "26.1816 mW", "26.5416 mW"],
        ["into the liquid (mode)", "6.95727 mW", "6.95727 mW"],
        ["into the bubble", "—", "357.66 µW"],
        ["kronspice ledger residual", "0.0", "0.0"],
    ], [0.4, 0.3, 0.3])
    r.p("The bubble takes about 5 % of what the flask's own losses take. export/bubble_on_flask.cir is the chain "
        "with the bubble element (RBUBBLE, CBUBBLE), with its KiCad sheet and the MNA against torn results.")

    r.h2("Part 3: reversibility")
    r.p("The knot's implosion is its explosion run backwards (Maxwell's equations are time-symmetric). Here the "
        "bubble runs from R<sub>max</sub> through the collapse and back out (2 × the fall time); then its momentum "
        "is flipped and the drive run backwards.")
    r.table([
        ["Viscosity", "Radiation", "|ΔR|/R on return", "|ΔΠ|/Π<sub>max</sub>"],
        ["off", "off", "2.5×10<super>−12</super>", "4.6×10<super>−13</super>"],
        ["on", "off", "6.5×10<super>−2</super>", "1.3×10<super>−2</super>"],
        ["off", "on", "9.0×10<super>−1</super>", "1.8×10<super>−2</super>"],
        ["on", "on", "8.0×10<super>−1</super>", "8.1×10<super>−3</super>"],
    ], [0.2, 0.2, 0.3, 0.3])
    r.p("Without losses the bubble is a conservative (Hamiltonian) system and retraces its collapse, like the knot. "
        "The radiated sound makes it irreversible, because the outgoing pulse does not come back; running a bubble "
        "backwards would need that pulse to converge on it again, the job H30's Huygens surface did for the knot.")

    r.h2("Part 4: the light and its spectrum")
    r.p("Free–free (bremsstrahlung) emission from argon ionized per the Saha equation, assuming a uniform, "
        "optically thin core (Rybicki &amp; Lightman eqs. 5.14b, 5.15b; Gaunt factor 1.2; NIST argon levels):")
    r.bullets([
        "energy per flash 12.7746 fJ (8.7×10<super>−7</super> of the collapse energy), peak 63.831 µW, "
        "ionization 5.8×10<super>−3</super> at the peak;",
        "flash FWHM 180 ps, about 6.2×10<super>4</super> photons (energy / kT at the peak);",
        "measured: 10<super>5</super>–10<super>7</super> photons per flash (about 10<super>−12</super> J), widths "
        "60–250 ps (Gompf et al., <i>PRL</i> 79, 1405, 1997). The estimate is low by one to two orders in photons "
        "and right in width.",
    ])
    r.table([
        ["Quantity", "Model", "Measured", "Agreement"],
        ["Width, UV band (250–300 nm) / red band (600–700 nm)", "164 / 176 ps", "equal (Gompf et al. 1997)", "yes"],
        ["Spectral shape", "peaks near 475 nm, falls into the UV",
         "rises into the UV to water's edge near 200 nm", "no"],
        ["Planck-fit temperature of the shape (250–700 nm)", "6489 K", "of order 10<super>4</super> K and above",
         "too cool"],
    ], [0.34, 0.22, 0.28, 0.16])
    r.p("Time-integrated spectrum (energy per nm, relative): 250 nm 0.55, 325 nm 0.81, 400 nm 0.95, 475 nm 1.00, "
        "550 nm 0.99, 625 nm 0.95, 700 nm 0.90. An optically thin gas at the model's 15,000 K peak emits most in "
        "the blue-green. The processes left out (electron–neutral bremsstrahlung, recombination, heat loss, the "
        "lowered ionization potential) are where the standard account finds the ultraviolet. None of the "
        "shortfall points to a non-thermal source. The light is far below the other ledger terms and is not taken "
        "out of the dynamics.")

    r.h2("Part 5: does the flash carry the knot's topology?")
    r.p("A knot of field energy carries net optical helicity; thermal light from random currents carries none. "
        "The bench measurement is the whole flash's degree of circular polarization, summed over every frequency "
        "and over 600 probes on a sphere (V &gt; 0 is spin along the direction of travel):")
    r.code("P_c = ΣV/ΣI,  V = −2 Im(Ê_θ*·Ê_φ),  I = |Ê_θ|² + |Ê_φ|²")
    r.p("The sources run in em/fdtd.go, scaled to a central wavelength λ<sub>0</sub> = 450 nm at 8 cells: the "
        "collapsed core (R<sub>min</sub> ≈ 0.79 µm) is 1.75 λ<sub>0</sub>, and the probes sit at 6 λ<sub>0</sub>.")
    r.table([
        ["Source", "P<sub>c</sub>", "Field ledger residual", "Energy absorbed"],
        ["hopfion of knotted_light (exact field at the probes +0.9999)", "+0.9996", "1.6×10<super>−13</super>",
         "0.937"],
        ["chiral dipole pair (p || m, quadrature)", "+0.9797", "5.7×10<super>−15</super>", "0.992"],
        ["rotating electric dipole", "+0.0007", "4.3×10<super>−15</super>", "0.992"],
        ["thermal core, seed 1", "+0.0087", "1.2×10<super>−15</super>", "0.999"],
        ["thermal core, seed 2", "−0.0015", "2.9×10<super>−15</super>", "0.999"],
        ["thermal core, seed 3", "−0.0268", "3.5×10<super>−14</super>", "0.999"],
        ["thermal core, seed 4", "+0.0296", "1.3×10<super>−14</super>", "0.999"],
        ["thermal core, seed 5", "+0.0069", "1.3×10<super>−14</super>", "0.999"],
        ["thermal core, seed 6", "+0.0176", "4.1×10<super>−15</super>", "0.999"],
        ["<b>thermal core, 6 seeds</b>", "mean +0.0058, rms 0.0184", "≤ 4×10<super>−14</super>", ""],
    ], [0.46, 0.2, 0.18, 0.16], bold_rows=(10,))
    r.bullets([
        "<b>The knot is one helicity.</b> Rañada's hopfion in the Bateman form satisfies i∂F/∂t = curl F, the "
        "condition for one helicity; the measurement sees +0.9996 and the exact field +0.9999.",
        "<b>An electric dipole alone cannot carry helicity</b>, even when circular on its axis: the two hemispheres "
        "cancel. Helicity needs electric and magnetic parts together, being the charge of electric–magnetic "
        "duality (Calkin, <i>Am. J. Phys.</i> 33, 958, 1965; Cameron, Barnett &amp; Yao, <i>New J. Phys.</i> 14, "
        "053050, 2012). Random hot-gas currents are electric currents.",
        "<b>The thermal core gives P<sub>c</sub> consistent with zero</b>: rms 0.0184 against the independent-"
        "radiator estimate 1/√(2N) = 0.0408. A real flash of 10<super>5</super>–10<super>7</super> photons has a "
        "counting floor 1/√N of 3×10<super>−3</super> to 3×10<super>−4</super>.",
    ])
    r.small("The interactive report (viz/report.html) shows these results as browser-drawn charts and 3-D scenes: "
            "the bubble beside the hopfion over one cycle, the radius and ledger curves, the polarization globe "
            "and field slices. Those charts are generated in the browser from embedded data and have no image "
            "files, so this PDF gives their numbers in the tables above.")

    r.h1("Validation")
    r.table([
        ["Test", "Checks", "Held to"],
        ["TestRayleighPlessetLimit", "lossless circuit = Rayleigh–Plesset", "10<super>−9</super>"],
        ["TestLinearAdmittance", "small-amplitude Y<sub>1</sub> = linear acoustic admittance", "10<super>−4</super>"],
        ["TestCycleLedgerAndTwin", "ledger; ½Re(Y<sub>1</sub>)P<sub>a</sub>² = cycle work × f; Keller–Miksis "
         "R<sub>max</sub>, R<sub>min</sub>, phase; flash plausible",
         "10<super>−8</super>; 10<super>−6</super>; 1 %, 10 %, 0.5 %"],
        ["TestChainWithBubble", "bubble element takes ½Re(Y<sub>1</sub>)P<sub>a</sub>², mode heat unchanged, "
         "kronspice ledger", "10<super>−3</super>; 10<super>−9</super>"],
        ["TestReversal", "lossless return, lossy failure", "10<super>−9</super>; &gt; 10<super>−2</super>"],
        ["TestStokesSign, TestIonization", "sign of V; Saha limits", "exact"],
        ["TestHelicity", "hopfion = exact, chiral &gt; 0.9, rotating and thermal ≈ 0, ledgers",
         "0.01; 10<super>−10</super>"],
        ["em.TestFDTDSharedEdgeSources", "two sources on one edge act as their sum, ledger closed (the solver fix "
         "this example needed)", "10<super>−10</super>"],
    ], [0.3, 0.5, 0.2])

    r.h1("Limits of the model")
    r.bullets([
        "The circuit is first order in Ṙ/c; near the collapse neither it nor Keller–Miksis is exact (wall Mach "
        "0.781 and 0.516).",
        "The light is the simplest standard estimate: uniform, optically thin core with Saha ionization and "
        "electron–ion free–free emission only. It is low by one to two orders in photons and too red.",
        "The polarization test uses scaled field-solver sources (λ<sub>0</sub> = 450 nm at 8 cells), 300 random "
        "current elements for the thermal core and six seeds.",
        "Part 2 is a harmonic balance at the fundamental: exact for power, not for the waveform.",
    ])

    r.h1("Verdict and status")
    r.bullets([
        "<b>Where the analogy holds:</b> the bubble is a driven implosion and explosion with its inward supply "
        "booked, one cycle per acoustic period: 14.3 nJ in per cycle from the flask's standing wave and 13.3 nJ "
        "out as a sound pulse. This is the structure H30 found for the pulsing knot, with the flask's standing "
        "wave as the declared medium.",
        "<b>Where it fails:</b> the collapse is a nonlinear, irreversible flow of matter, stopped by gas "
        "pressure; its flow is spherically symmetric and carries no topology. The light is thermal in every test "
        "made: its width does not depend on colour, and its polarization is zero within the speckle of random "
        "currents, against +0.9996 for the knot.",
        "<b>Bench test:</b> a polarimeter on single flashes. A knotted field would give P<sub>c</sub> of order 1; "
        "thermal light gives at most the counting floor. The two are separated by about three orders of "
        "magnitude.",
    ])
    r.p("H30 is recorded as <i>bounded</i> in the register; this example is its bench case, the sonoluminescing "
        "bubble as the same driven cycle with its inward supply booked.")

    r.h1("Reproducing the results")
    r.code("go run ./examples/bubble_collapse            # about 1 min\n"
           "go run ./examples/bubble_collapse -nofield   # parts 1–4 only\n"
           "go test ./examples/bubble_collapse           # -short skips the field solver\n"
           "go run ./examples/bubble_collapse -data report_data.json\n"
           "go run ./examples/bubble_collapse -report examples/bubble_collapse/viz/report.html\n"
           "python examples/bubble_collapse/viz/make_report.py examples/bubble_collapse/bubble_collapse_report.pdf")
    r.small("Sources: Brenner, Hilgenfeldt &amp; Lohse, <i>Rev. Mod. Phys.</i> 74, 425 (2002); Rayleigh, <i>Phil. "
            "Mag.</i> 34, 94 (1917); Plesset, <i>J. Appl. Mech.</i> 16, 277 (1949); Keller &amp; Miksis, <i>J. "
            "Acoust. Soc. Am.</i> 68, 628 (1980); Kinsler et al., <i>Fundamentals of Acoustics</i>, 4th ed., §7.2; "
            "Kron (1951), Epilogue eqs. 3–4; Hilgenfeldt, Grossmann &amp; Lohse, <i>Nature</i> 398, 402 (1999); "
            "Rybicki &amp; Lightman (1979); Gompf et al., <i>PRL</i> 79, 1405 (1997); IAPWS R1-76(2014), R12-08; "
            "Calkin (1965); Cameron, Barnett &amp; Yao (2012); Rañada (1989); Kedia et al., <i>PRL</i> 111, "
            "150404 (2013).")
    r.build(out)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "bubble_collapse_report.pdf"))
