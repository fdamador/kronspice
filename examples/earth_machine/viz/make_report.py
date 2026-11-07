# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""PDF report for examples/earth_machine (the Earth as three Kron networks).

python examples/earth_machine/viz/make_report.py examples/earth_machine/earth_machine_report.pdf

Numbers come from examples/earth_machine/run_output.txt (layers 1-5 and the
studies) and, for layer 7 (nonlinear ionosphere, not in run_output.txt), from
examples/earth_machine/README.md. Figures are the PNGs in viz/.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401


def fig(name):
    return os.path.join(HERE, name)


def build(out):
    r = Report("examples/earth_machine", "The Earth as an electromagnetic machine",
               "Geodynamo, global atmospheric circuit and Schumann cavity as Kron networks "
               "— CLAUDE.md §3a, §3f", "CLAUDE.md §3a, §3f",
               "dimensionless Rikitake μ = 1, A = 1")

    r.h1("Summary")
    r.p("The planet is modelled as three Kron networks built from kronspice elements. Layer 1 is the core "
        "dynamo as Rikitake's two-disk dynamo (two cross-coupled homopolar machines, machine elements N). "
        "Layer 2 is the global atmospheric electric circuit: storms charge the ionosphere through leaky "
        "fair-weather air, modelled as 60 shells of R || C. Layer 3 is the Earth–ionosphere cavity as a "
        "transmission line in colatitude, which carries the Schumann resonances. The layers are then linked "
        "in one network, the dynamo is put through convergence, parameter and drive studies, a stacked "
        "two-disk design is tested, and a nonlinear (electron-heated) ionosphere is added.")
    r.bullets([
        "<b>The Rikitake dynamo reverses chaotically:</b> 71 reversals in t = 50…1000, intervals 3.92–70.36, mean "
        "12.17, CV 0.92 (Poisson: 1). Equilibrium I<sub>1</sub>² = 1.618034 (the golden ratio), natural period 4.2018.",
        "<b>Energy closes:</b> drive work 2581 = heat 2578 + Δstored 3.078, residual 1.8×10<super>−5</super> %; "
        "engine ledger mismatch 1.1×10<super>−5</super> J (4×10<super>−9</super> of the drive work).",
        "<b>Global circuit:</b> calibrated on 250 kV, 2 pA/m² and 120 V/m, it gives 1020 A, 244.9 Ω, 255 MW, "
        "an Earth charge of −5.42×10<super>5</super> C, and a relaxation time ε<sub>0</sub>/σ<sub>0</sub> = 531 s. "
        "Half the balancing positive charge sits below 1.42 km.",
        "<b>Schumann cavity:</b> the perfect-cavity ladder equals Legendre's c/(2πa)·√(n(n+1)) to 4×10<super>−4</super>. "
        "A two-height ionosphere fitted to the five measured modes reproduces them within 2 % (rms error 0.0095), "
        "with predicted Q<sub>2</sub>–Q<sub>5</sub> of 4.7–4.2.",
        "<b>One network for layers 2 and 3:</b> the same elements give 249.59–250.27 kV at DC and ELF peaks "
        "7.2, 13.0, 18.9, 24.9, 30.9 Hz with no fitted electric height (measured 7.8, 14.3, 20.8, 27.3, 33.8 Hz).",
        "<b>Stacked disks do not reverse:</b> any self-coupling drains the speed difference A = ω<sub>1</sub> − ω<sub>2</sub> "
        "and the field locks to one polarity; no gap from 0.53 to 6 r<sub>o</sub> reverses.",
        "<b>No layer powers another:</b> every drive study (paper107 variations a–e) closes its ledger; none adds energy.",
    ])

    r.h1("What is being modelled")
    r.p("This example is standard physics, not a register hypothesis. It applies Kron's network method (§3a) "
        "and the homopolar machine of §3f to the Earth. The paper107 drive studies (Evans, Amador &amp; Eckardt) "
        "give the standard Maxwell–Lorentz counterpart of the ECE reading, which stays in the register as H15 and H16. "
        "The dynamo results are dimensionless: there is no claim about reversal times in years.")
    r.figure(fig("fig0_rikitake_circuit.png"),
             "Figure 1. Equivalent circuit of the Rikitake dynamo: two loops, each with a disk (speed EMF "
             "Mω<sub>i</sub>·I<sub>j</sub>), a resistance R and a field coil L whose field threads the other disk. "
             "Each shaft is driven by a torque G.", 0.85)

    r.h1("Layer 1 — the Rikitake two-disk dynamo")
    r.h2("Model")
    r.table([
        ["part", "built as"],
        ["disk <i>i</i>", "machine element with two windings: its radial path and the other loop's field coil; "
         "speed matrix G[0][1] = ±M gives the speed EMF M·ω<sub>i</sub>·I<sub>coil</sub> and the back torque M·I<sub>1</sub>·I<sub>2</sub>"],
        ["shaft <i>i</i>", "mechanical node, inertia C (capacitor), drive torque G (current source on the mechanical channel)"],
        ["loop <i>i</i>", "disk path → R → field coil L (threads the other disk)"],
    ], [0.18, 0.82])
    r.p("The equations are L dI<sub>1</sub>/dt = −R I<sub>1</sub> + M ω<sub>1</sub> I<sub>2</sub>, "
        "L dI<sub>2</sub>/dt = −R I<sub>2</sub> + M ω<sub>2</sub> I<sub>1</sub>, C dω<sub>i</sub>/dt = G<sub>i</sub> − M I<sub>1</sub> I<sub>2</sub>. "
        "With all constants 1 this is Rikitake's dimensionless system with μ = 1, A = 1. The linearization about the "
        "equilibrium is (λ + 2)(λ² + k² + k<super>−2</super>): a neutral oscillation of period 4.20, the floor of the reversal intervals.")
    r.h2("Four cases")
    r.table([
        ["case", "reversals", "intervals (min / mean / max)", "CV", "energy residual (trace)", "engine ledger"],
        ["1 co-rotating (Rikitake)", "71 (RK4: 49)", "3.92 / 12.17 / 70.36", "0.92", "1.8×10<super>−5</super> %", "1.1×10<super>−5</super> J"],
        ["2 disk 2 reversed, same wiring", "0", "—", "—", "6.7×10<super>−5</super> %", "−4.2×10<super>−5</super> J"],
        ["3 disk 2 reversed, coil 2 rewound", "71 (RK4: 49)", "3.92 / 12.17 / 70.36", "0.92", "1.8×10<super>−5</super> %", "1.1×10<super>−5</super> J"],
        ["4 co-rotating, coil 2 rewound (Tesla's opposed fields)", "0", "—", "—", "6.7×10<super>−5</super> %", "−4.1×10<super>−5</super> J"],
    ], [0.30, 0.12, 0.18, 0.07, 0.16, 0.17], bold_rows=(1,))
    r.bullets([
        "Case 1: the circuit tracks an independent RK4 integration to 0.0041 for t ≤ 10; after that chaos separates them, "
        "so only the statistics are comparable. Final shaft speeds ω<sub>1</sub> = 2.47, ω<sub>2</sub> = 1.47.",
        "Case 2: the feedback is negative; |I<sub>1</sub>| after t = 50 is 6.67×10<super>−14</super> and both shafts spin up "
        "without load (ω = ±10³).",
        "Case 3 equals case 1 with I<sub>2</sub> and ω<sub>2</sub> negated (test TestCounterRotatingRewoundMirrors, to 10<super>−9</super>): "
        "the sense of the feedback loop matters, not the direction of spin.",
        "Case 4: about I = 0, λ = −1 ± i; the field dies. Tesla's magnets excite his disks separately.",
    ])
    r.figure(fig("fig1_timeseries.png"),
             "Figure 2. Case 1 transient from kronspice: I<sub>1</sub>, I<sub>2</sub> and ω<sub>1</sub> against "
             "dimensionless time (t = 50…1000). Polarity reversals appear as sign changes of I<sub>1</sub>.", 0.8)
    r.figure(fig("fig2_attractor.png"),
             "Figure 3. The Rikitake attractor traced by the kronspice transient in (I<sub>1</sub>, I<sub>2</sub>, ω<sub>1</sub>).", 0.55)
    r.figure(fig("fig3_reversals.png"),
             "Figure 4. Histogram of the 71 reversal intervals of case 1 against an exponential (Poisson) density "
             "with the same mean, 12.2. No interval is shorter than about 4, the oscillation period.", 0.65)

    r.h2("Layer 1 studies: convergence and the (μ, A) plane")
    r.p("Counts at a fixed time depend on the integrator; the interval distribution does not. RK4 over t = 20 000:")
    r.table([
        ["step h", "reversals", "mean", "CV", "KS vs h = 0.0005: D, p"],
        ["0.004", "1385", "14.35", "0.99", "0.062, 0.01"],
        ["0.002", "1442", "13.78", "0.99", "0.027, 0.65"],
        ["0.001", "1470", "13.52", "1.00", "0.022, 0.88"],
        ["0.0005", "1473", "13.49", "1.01", "—"],
    ], [0.15, 0.15, 0.15, 0.15, 0.4])
    r.p("The kronspice circuit at maximum step 0.02, 0.01 and 0.005 (t = 4000) gives 269, 307 and 291 reversals and KS p = "
        "0.13, 0.38, 0.39 against the finest RK4. A parametric bootstrap KS test (400 resamples) rejects exponential "
        "intervals at every (μ, A) point tried (p = 0.0025, the floor). μ = 0.5 with A ≥ 2 is periodic (CV 0.00). "
        "Longest/mean is at most 6.0, against of order 100 or more in the polarity record (the Cretaceous Normal Superchron, ≈ 38 Myr, "
        "against a Cenozoic mean of 0.2–0.3 Myr). Where the machine is chaotic (μ = 1, 2), the last swing before a reversal "
        "is 3–23 % above the median swing (ratio 1.03–1.23), whereas the record shows the dipole weakening before it reverses "
        "(Valet &amp; Meynadier 1993).")
    r.figure(fig("fig6_reversal_statistics.png"),
             "Figure 5. (a) Fraction of intervals longer than a given multiple of the mean, Rikitake μ = 1, A = 1 "
             "(RK4, t = 20 000), against an exponential: a floor and no superchrons. (b) Longest/mean interval over "
             "the (μ, A) plane; the polarity record has of order 100 or more.", 0.9)

    r.h2("paper107 drive variations (task 4)")
    r.table([
        ["variation", "built as", "result"],
        ["a. modulated torque G = 1 + 0.3 sin(ω<sub>d</sub>t)", "SIN current sources on the shafts",
         "reversal rate 0.0326–0.0865 per unit time over the run_output sweep (0.0747 unmodulated); RMS I<sub>1</sub> 1.313–1.343; "
         "circuit vs RK4 |ΔI<sub>1</sub>| ≤ 0.0076"],
        ["b. each disk on its own speed drive", "voltage sources on the mechanical nodes",
         "growth (−R + M√(ω<sub>1</sub>ω<sub>2</sub>))/L; circuit = exact to 4 decimals; opposed speed modulation decays at the old threshold"],
        ["c. feedback coil around each disk", "G[0][0] = −K·M",
         "K ≥ +0.25: no reversals; K = −0.5: 21 reversals against 1438; circuit vs RK4 |ΔI<sub>1</sub>| ≤ 0.0014"],
        ["d. Tesla's opposed fields", "case 4", "λ = −1 ± i about I = 0; |I<sub>1</sub>| after t = 60 is 1.7×10<super>−27</super>"],
        ["e. rotating-frame terms", "one disk at ω(t) = 5 sin 10t",
         "machine = exact to 7×10<super>−6</super>; for a 1 T, 10 cm, 50 rev/s disk the Stewart–Tolman, centrifugal and "
         "Coriolis terms are 6×10<super>−11</super>, 2×10<super>−9</super>, 4×10<super>−9</super> of the 1.571 V motional EMF"],
    ], [0.25, 0.22, 0.53])
    r.p("None adds energy and every ledger closes. The feedback coil changes how the drive's energy splits between "
        "field and heat: its EMF K·M·ω·I is paid by the back torque K·M·I² on the shaft.")
    r.figure(fig("fig7_paper107_studies.png"),
             "Figure 6. (a) Reversal rate against torque modulation frequency ω<sub>d</sub>/ω<sub>n</sub> (a finer sweep "
             "than the run_output table; peak near 0.8 ω<sub>n</sub>). (b) Growth rate of the field with both disks on "
             "speed drives: it grows only where ω<sub>1</sub>ω<sub>2</sub> &gt; 1. (c) Reversals in t = 20 000 against "
             "the feedback-coil gain K.", 1.0)

    r.h1("Layer 2 — the global atmospheric circuit")
    r.p("Wilson's circuit (Rycroft, Israelsson &amp; Price 2000): storms lift positive charge to the ionosphere, and the "
        "return current flows down through fair-weather air. Each spherical shell of air is R || C (conduction plus "
        "displacement current). 60 shells, thin near the ground, form a ladder; the storms are one current source "
        "booked on the thermal channel. The conductivity σ(z) = σ<sub>0</sub>·exp(z/H) is calibrated on three measured "
        "anchors, V<sub>I</sub> = 250 kV, J = 2 pA/m² and E<sub>0</sub> = 120 V/m, which give σ<sub>0</sub> = "
        "1.67×10<super>−14</super> S/m and H = 2.08 km.")
    r.table([
        ["quantity", "result"],
        ["global current J·4πR²", "1020 A"],
        ["global resistance", "244.9 Ω"],
        ["V<sub>I</sub>, circuit vs exact ladder", "249.84 kV, rel. diff 9×10<super>−6</super> (sphericity −0.07 % from the flat column)"],
        ["Earth's charge, circuit vs −ε<sub>0</sub>E<sub>0</sub>·4πR²", "−5.42×10<super>5</super> C, both"],
        ["balancing positive charge", "5.42×10<super>5</super> C in the air, half below 1.42 km; 1.7×10<super>−6</super> C at the ionosphere"],
        ["capacitance", "geometric shell to 60 km 0.076 F; effective Q/V<sub>I</sub> 2.17 F"],
        ["relaxation time", "R·C<sub>eff</sub> = ε<sub>0</sub>/σ<sub>0</sub> = 531 s (8.9 min)"],
        ["power", "255 MW from the storms, all to heat; ledger residual −4×10<super>−6</super>"],
        ["switch-on", "63 % of V<sub>I</sub> in 229 s; circuit vs exact ladder within 8×10<super>−5</super> of V<sub>I</sub>"],
        ["Carnegie curve (±15 % daily storm current)", "V<sub>I</sub> follows with amplitude ratio 0.999691, lag 4.4 min"],
    ], [0.38, 0.62])
    r.p("The charge that balances the Earth's sits in the lowest kilometres of air, not on the ionosphere, so the "
        "spherical-capacitor picture understates the capacitance about 30 times.")
    r.p("<b>Refinements</b> (declared, illustrative inputs; the three global anchors still calibrate):")
    r.bullets([
        "<b>Two-scale profile</b> (boundary layer 1.5 km, free-air scale height 5 km; solved H<sub>BL</sub> = 1.09 km): "
        "V<sub>I</sub> 249.67 kV, E<sub>0</sub> 119.8 V/m, half the space charge below 0.72 km; τ = 531 s unchanged "
        "(R·C<sub>eff</sub> = ε<sub>0</sub>/σ<sub>0</sub> for any profile).",
        "<b>Land (29 %, σ<sub>0</sub> 1×10<super>−14</super> S/m) and ocean (2.2×10<super>−14</super> S/m) in parallel:</b> "
        "V<sub>I</sub> 224.8 kV; land 160 A, ocean 860 A; E<sub>0</sub> 108 V/m over both.",
        "<b>Storm generators with their own paths</b> (up, down, local return), with the generator current set so that "
        "1020 A still reaches the ionosphere; the rest returns locally:",
    ])
    r.table([
        ["R<sub>up</sub> [Ω]", "R<sub>down</sub> [Ω]", "R<sub>cloud</sub> [Ω]", "I<sub>s</sub> [A]", "reaching ionosphere",
         "P<sub>gen</sub> [MW]", "in storms [MW]", "fair weather [MW]", "ledger"],
        ["50", "100", "1000", "1423", "71.7 %", "573", "318", "255", "3×10<super>−7</super>"],
        ["100", "100", "500", "1928", "52.9 %", "875", "620", "255", "5×10<super>−7</super>"],
        ["100", "500", "250", "4468", "22.8 %", "3851", "3596", "255", "1×10<super>−7</super>"],
    ], [0.1, 0.1, 0.1, 0.1, 0.14, 0.1, 0.12, 0.12, 0.12])
    r.figure(fig("fig4_global_circuit.png"),
             "Figure 7. (a) The calibrated conductivity profiles (one scale, H = 2.08 km; two-scale). (b) Fraction of "
             "the positive space charge below height z: half lies below 1.42 km. (c) V<sub>I</sub> after the storms "
             "switch on: kronspice transient against the exact ladder sum.", 1.0)

    r.h1("Layer 3 — the Earth–ionosphere cavity")
    r.p("The shell between ground and ionosphere is a line in colatitude (Madden &amp; Thompson 1965): V = E<sub>r</sub>·h, "
        "I = H<sub>φ</sub>·2πa·sinθ, L′ = μ<sub>0</sub>h/(2πa sinθ), C′ = ε<sub>0</sub>·2πa sinθ/h, so the telegrapher's "
        "equations become Legendre's with f<sub>n</sub> = c/(2πa)·√(n(n+1)). The ladder has 181 bands of 1°; a 1 A AC "
        "current source into the pole band stands for a vertical lightning stroke.")
    r.table([
        ["n", "perfect ladder [Hz]", "Legendre [Hz]", "rel.", "constant height: f [Hz], Q", "two-height: f [Hz], Q", "measured [Hz]", "two-height / measured"],
        ["1", "10.5910", "10.5913", "−3×10<super>−5</super>", "7.80, 5.0", "7.96, 5.0", "7.8", "1.020"],
        ["2", "18.3432", "18.3446", "−8×10<super>−5</super>", "13.55, 8.7", "14.27, 4.7", "14.3", "0.998"],
        ["3", "25.9393", "25.9432", "−0.0002", "19.18, 12.2", "20.64, 4.5", "20.8", "0.992"],
        ["4", "33.4841", "33.4926", "−0.0003", "24.77, 15.8", "27.11, 4.3", "27.3", "0.993"],
        ["5", "41.0044", "41.0198", "−0.0004", "30.34, 19.4", "33.68, 4.2", "33.8", "0.996"],
    ], [0.05, 0.13, 0.12, 0.11, 0.16, 0.15, 0.12, 0.16])
    r.bullets([
        "<b>Constant heights, fitted to mode 1 only</b> (7.8 Hz, Q = 5; Sentman 1995): v/c = √(h<sub>e</sub>/h<sub>m</sub>) = 0.7402 "
        "(h<sub>e</sub> = 49.3 km at h<sub>m</sub> = 90 km), α = 4.93 s<super>−1</super>. Modes 2–5 come out up to 10 % low "
        "(0.948, 0.922, 0.907, 0.898 of measured) and Q grows with f.",
        "<b>Two-scale-height ionosphere</b> (Greifinger &amp; Greifinger 1978; Sentman 1990), with the ladder rebuilt at every "
        "frequency (781 ladders, AC ledger residual ≤ 0): h<sub>e</sub> = 52.3 km, h<sub>m</sub> = 90.2 km at 8 Hz, "
        "ζ<sub>e</sub> = 0.55 km, ζ<sub>m</sub> = 9.94 km. Q<sub>2</sub>–Q<sub>5</sub> are predictions and fall inside the "
        "reported 4–6. The frequencies are well constrained; the split between the two scale heights is not.",
        "An equator station sees only the even modes (P<sub>n</sub>(0) = 0 for odd n): peaks at 13.58, 24.78 and 35.91 Hz.",
    ])
    r.figure(fig("fig5_schumann.png"),
             "Figure 8. (a) Perfect cavity: ladder/Legendre error against band width, second order. (b) Spectra at "
             "the pole and the equator with constant heights, and at the pole with the two-height ionosphere; dotted "
             "lines are the measured modes. (c) Model/measured frequency per mode (solid) and Q (dotted, right axis; "
             "band: reported 4–6).", 1.0)

    r.h1("Linking the layers")
    r.bullets([
        "<b>One network for layers 2 and 3:</b> 61 bands × 40 shells to 100 km (two-scale profile), series path at "
        "h<sub>m</sub> = 90 km. DC: ionosphere 249.59–250.27 kV across all bands, 255 MW, ledger 2×10<super>−9</super>. "
        "ELF, same elements: peaks 7.2, 13.0, 18.9, 24.9, 30.9, 37.0 Hz at the pole and 13.5, 25.3, 37.2 Hz at the equator. "
        "The electric height comes from σ(z) = ωε<sub>0</sub> in the profile: 45.6 km at 8 Hz, 52.7 km at 33 Hz (layer-3 fit 52.3 km).",
        "<b>Dynamo → atmosphere:</b> the dipole shields the air from galactic cosmic rays, and σ is proportional to √q. With the unshielded "
        "ion production declared as 2× today's, V<sub>I</sub> falls to 176.8 kV (circuit 176.7 kV) and τ to 376 s when the "
        "dipole vanishes. Over a Rikitake run V<sub>I</sub> has a mean of 219 kV and dips to 177 kV at each reversal. The link "
        "is quasi-static: minutes for the air against millennia for a reversal.",
        "<b>One ledger:</b> the dynamo on the mechanical channel (geodynamo ohmic losses estimated at 0.2–0.5 TW, Christensen "
        "&amp; Tilgner 2004), the storms on the thermal channel (0.26–3.9 GW with the storm paths); the cavity draws on "
        "lightning. No layer supplies another with energy; the dynamo only modulates the air's conductivity.",
    ])
    r.figure(fig("fig8_layers_linked.png"),
             "Figure 9. (a) One network of air columns plus cavity, with no fitted h<sub>e</sub>: ELF spectrum at the stroke "
             "point (pole) and at the equator; dotted lines are measured modes. (b) Rikitake I<sub>1</sub> (blue) and the "
             "resulting ionospheric potential V<sub>I</sub> (red; unshielded ion production ×2, declared).", 0.95)

    r.h1("The stacked design")
    r.p("The two disks (eddies) sit on one axis, each inside the coil fed by the other. The couplings come from the copper "
        "machine of examples/homopolar scaled to the disk radius (coil 1.2 r<sub>o</sub>, 0.4 r<sub>o</sub> long). KS is each coil's "
        "stray field on its own disk, as a fraction of M; KM is the coil–coil coupling. In the circuit the coils are "
        "inductors with a Mutual element; it equals RK4 to 0.0012 at the 0.53 r<sub>o</sub> gap and the ledger closes "
        "(−3×10<super>−6</super> J).")
    r.table([
        ["gap [r<sub>o</sub>]", "KS", "KM", "reversals (t = 20 000)", "RMS I<sub>1</sub>", "mean ω<sub>1</sub>, ω<sub>2</sub>"],
        ["side by side", "0", "0", "1438 (mean 13.76, CV 0.99)", "1.343", "—"],
        ["0.53", "0.5190", "0.3870", "0", "0.811", "0.66, 0.66"],
        ["1.00", "0.2630", "0.1932", "0", "0.890", "0.79, 0.79"],
        ["2.00", "0.0813", "0.0621", "0", "0.962", "0.92, 0.92"],
        ["6.00", "0.0052", "0.0042", "0", "0.998", "1.00, 0.99"],
    ], [0.16, 0.12, 0.12, 0.28, 0.14, 0.18])
    r.p("Rikitake's chaos rests on a conserved speed difference A = ω<sub>1</sub> − ω<sub>2</sub>: side by side, both disks "
        "feel the same back torque M·I<sub>1</sub>·I<sub>2</sub>. Stacked, each disk also feels KS·M·I<sub>i</sub>², "
        "which differs between the two disks, so A drains away (gap 6 r<sub>o</sub>: A = 1.00 at t = 0, 0.77 at 50, 0.35 at 200, "
        "0.005 at 1000) and the field locks to one polarity. Holding A with unequal drives G = 1 ± δ brings reversals "
        "back only at gap 0.53 r<sub>o</sub> with δ = 0.5 (5119 reversals), and then like a clock: interval 3.87, CV 0.13.")
    r.figure(fig("fig10_stacked.png"),
             "Figure 10. (a) Two eddies stacked on one axis, each in the other's coil (KS = 0.52, KM = 0.39). (b) I<sub>1</sub>: "
             "side by side it reverses, stacked it holds one polarity. (c) The speed difference A drains when stacked.", 1.0)
    r.figure(fig("fig9_fluid_machine.png"),
             "Figure 11. The side-by-side machine drawn as two eddies in the liquid core, each one's current "
             "building the field that threads the other (schematic; streamline shading shows flow speed).", 0.75)

    r.h1("Layer 7 — a nonlinear ionosphere")
    r.p("The numbers in this section are from README.md; layer 7 is not part of run_output.txt. The layer-3 cavity (fitted to "
        "mode 1, Q 5) has its ionospheric loss made field-dependent by electron heating, G(E) = G<sub>0</sub>/(1 + (E/E<sub>p</sub>)²) "
        "(Gurevich 1978), the mechanism behind the Luxembourg effect. The heating field E<sub>p</sub> ≈ 0.02 V/m is declared for "
        "the night-time D region (T 200 K, δ 10<super>−3</super>, ν 10<super>6</super> s<super>−1</super>). The field at a band "
        "is taken as V/h<sub>e</sub>, an upper bound, so V<sub>h</sub> = E<sub>p</sub>·h<sub>e</sub> ≈ 1 kV. Two transmitters "
        "at 7.8 and 14.3 Hz drive the source band; the station is 30° away.")
    r.table([
        ["source peak / V<sub>h</sub>", "drive", "station f<sub>1</sub>, f<sub>2</sub>", "2f<sub>1</sub>−f<sub>2</sub> (1.3 Hz)",
         "2f<sub>2</sub>−f<sub>1</sub> (20.8 Hz)", "linear cavity at those lines"],
        ["0.01", "1 A", "2.9 V", "7×10<super>−6</super> V", "2×10<super>−5</super> V", "1×10<super>−6</super>, 8×10<super>−6</super> V"],
        ["0.1", "10 A", "29 V", "3.6 mV", "3.2 mV", "4×10<super>−4</super>, 3×10<super>−4</super> V"],
        ["0.5", "49 A", "150 V", "0.48 V", "0.40 V", "1×10<super>−3</super>, 9×10<super>−4</super> V"],
    ], [0.14, 0.09, 0.15, 0.18, 0.18, 0.26])
    r.p("With heating the third-order (four-wave) products appear, and 2f<sub>2</sub> − f<sub>1</sub> = 20.8 Hz falls on the "
        "measured third Schumann frequency. The 7.8 Hz line changes by +3 % at 0.5 V<sub>h</sub>; its sign follows from the "
        "assumption that the loss is proportional to σ(E). The ledger closes to 3×10<super>−5</super> in both runs. V<sub>h</sub> corresponds to "
        "about 20 mV/m over the whole cavity height; natural Schumann fields are far smaller and real ELF transmitters stay "
        "in the linear regime, so at ELF the cavity is linear in practice.")
    r.figure(fig("ionosphere_nl.png"),
             "Figure 12. Two transmitters (7.8 + 14.3 Hz) at 0.5 V<sub>h</sub>, station 30° away. Top: station voltage with "
             "a linear and a nonlinear ionosphere. Middle: their difference (the products of the two waves). Bottom: "
             "station spectrum; new lines (2f<sub>1</sub>−f<sub>2</sub>, 2f<sub>2</sub>−f<sub>1</sub>, 3f<sub>1</sub>, …) appear "
             "only with the nonlinear ionosphere.", 0.65)

    r.h1("Ledger and validation")
    r.bullets([
        "Rikitake case 1: drive in 2581 = heat 2578 + Δstored 3.078 (residual 0.00047, 1.8×10<super>−5</super> %); engine "
        "ledger 1.1×10<super>−5</super> J. An earlier start-up mismatch of −3.75×10<super>4</super> J from the UIC point is "
        "fixed in analysis/tran.go (TestUICImposedPointNotIntegrated).",
        "Circuit vs RK4: case 1 0.0041 (t ≤ 10); task 4a 0.0076; task 4c 0.0014; stacked 0.0012; task 4e 7×10<super>−6</super> of peak.",
        "Global circuit: V<sub>I</sub> circuit vs exact ladder 9×10<super>−6</super>; switch-on 8×10<super>−5</super>; ledger −4×10<super>−6</super>; "
        "storm-path ledgers 1–5×10<super>−7</super>.",
        "Cavity: perfect ladder vs Legendre ≤ 4×10<super>−4</super>, second order in band width; AC ledger residual ≤ 0 "
        "in every sweep; linked network DC ledger 2×10<super>−9</super>.",
    ])

    r.h1("Limits of the model")
    r.bullets([
        "The dynamo is dimensionless; there is no claim about reversal times in years.",
        "Rikitake gives no superchrons (longest/mean ≤ 6 against of order 100 or more) and no field weakening before reversals; "
        "self-consistent 3-D dynamos (Glatzmaier &amp; Roberts 1995) are the benchmark it does not reach.",
        "Reversal counts at a fixed time are sensitive to the integrator; only interval statistics are comparable.",
        "The one-scale conductivity profile is an approximation; storm generator paths, land/ocean values and the "
        "unshielded ion production are declared, illustrative inputs.",
        "The two-height fit trades ζ<sub>e</sub> against h<sub>e0</sub>; ζ<sub>e</sub> = 0.55 km is small for a lower-ionosphere scale height.",
        "Layer 7 takes the cavity's vertical field as an upper bound on what reaches the plasma, and the sign of the "
        "line shift rests on the loss being proportional to σ(E).",
    ])

    r.h1("Reproducing the results")
    r.code("go run ./examples/earth_machine                 # everything (several minutes)\n"
           "go run ./examples/earth_machine -layer N        # 1 geodynamo, 2 global circuit, 3 cavity,\n"
           "                                                # 4 layer-1 studies, 5 the links, 6 stacked,\n"
           "                                                # 7 nonlinear ionosphere\n"
           "go run ./examples/earth_machine -data viz/data  # figure data (CSV), then: python viz/plot_layers.py\n"
           "go test ./examples/earth_machine\n"
           "python examples/earth_machine/viz/make_report.py examples/earth_machine/earth_machine_report.pdf")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "earth_machine_report.pdf"))
