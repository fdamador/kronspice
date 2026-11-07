# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""PDF report for examples/dna_vibration (DNA bases and Watson-Crick pairs as circuits).

python examples/dna_vibration/viz/make_report.py examples/dna_vibration/dna_vibration_report.pdf

Numbers come from examples/dna_vibration/run_output.txt (bases and pairs),
summary/measured.csv (model against measurement), stacks/summary*.csv (stacks)
and README.md (method). Figures are the PNGs in summary/.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401

NH2 = "NH<sub>2</sub>"


def fig(name):
    return os.path.join(HERE, "..", "summary", name)


def build(out):
    r = Report("examples/dna_vibration", "DNA bases and base pairs as circuits",
               "Infrared absorption and release of the four bases and the Watson–Crick pairs from kronspice's AC "
               "analysis — doc/chemistry-plan.md", "doc/chemistry-plan.md",
               "AMBER OL15, harmonic, gas phase")

    r.h1("Summary")
    r.p("The four DNA bases (methyl-capped) and the two Watson–Crick pairs A·T and G·C are built as electric circuits "
        "with the mobility analogy: atoms are capacitors (mass), every bond and angle of the AMBER force field is a spring "
        "inductor, and Kron's connection tensor joins the bond frame to the atom frame. Light enters as a current source at "
        "every atom. kronspice's AC analysis then gives the infrared light each molecule absorbs, and two independent routes "
        "give the light it releases. Stacks of 2–4 pairs show how stacking splits each line into a band.")
    r.bullets([
        "<b>All six molecules are true minima</b> (RMS gradient ≤ 4.3×10<super>−10</super> kcal/mol/Å, all 3N−6 vibrations "
        "real) and electrically small (d/λ = 2.1×10<super>−4</super> to 4.4×10<super>−4</super>, lumped by the λ/10 rule).",
        "<b>Thermal infrared release at 300 K</b> is 1.28×10<super>−20</super> W (T) to 5.20×10<super>−20</super> W (G·C) per molecule.",
        "<b>N–H stretches sit 3.5–6.0 % below the measurements,</b> the harmonic AMBER value against anharmonic data; C=O "
        "stretches range from +1.5 % (T) to −11.9 % (G).",
        "<b>Pairing red shifts are not reproduced:</b> with fixed charges the H-bonded stretches move up on pairing (T N3–H "
        "+27.6, G N2–H21 +126.3 cm<super>−1</super>), where measurements show shifts of −39 to −186 cm<super>−1</super>.",
        "<b>Modes living in the pairing</b> (≥ 25 % of their energy in the inter-base forces): A·T at 24, 45, 54, 72 and "
        "123 cm<super>−1</super>; G·C at 62, 63, 116 and 154 cm<super>−1</super>. The A·T H-bond stretch at 123.5 cm<super>−1</super> "
        "(one pair, vacuum) compares with 94 cm<super>−1</super> measured in a duplex in water (+31 %).",
        "<b>Stacking splits lines into bands</b> and adds rigid-pair modes below 200 cm<super>−1</super>; the lowest mode falls "
        "from 12.7 cm<super>−1</super> (one A·T) to 3.4 cm<super>−1</super> (four).",
    ])

    r.h1("What is being modelled")
    r.p("This is the chemistry branch's step after water_vibration (doc/chemistry-plan.md): standard molecular mechanics, "
        "not a register hypothesis. The question is whether a circuit built term by term from a published force field gives "
        "the same vibrations, and how close those are to measured gas-phase, matrix and helium-droplet spectra.")

    r.h1("The model")
    r.bullets([
        "<b>Molecules.</b> 9-methyladenine, 1-methylthymine, 9-methylguanine and 1-methylcytosine (the methyl stands in for "
        "the sugar, as in gas-phase spectroscopy) and the pairs A·T and G·C, cut from the all-atom B-DNA dodecamer "
        "d(CGCGAATTCGCG) and minimized in the AMBER DNA force field OL15 (parm99/bsc0 base terms). prep/prep.py writes "
        "molecules.json from OpenMM's DNA.OL15.xml and bdna.xyz. The methyl caps are this example's choice.",
        "<b>Circuit</b> (circuit.go). Atoms are capacitors (mass); each AMBER bond and angle is a spring branch: a chain of "
        "VCVS forms the rate of change of its internal coordinate (Wilson's B row), an inductor carries the spring force and "
        "CCCS hand it back to each atom. Everything else (torsions, out-of-plane, electrostatics, van der Waals, pre-stress) "
        "is the residual stiffness R = K − Σ k·BᵀB, applied through displacement integrators and one VCCS per non-zero entry. "
        "Forces between the two bases of a pair form their own group G<sub>HB</sub>; in a stack, forces between bases of "
        "different pairs form G<sub>ST</sub>. The circuit's stiffness equals K up to entries below 10<super>−9</super> of the largest.",
        "<b>Outputs</b> (analysis.go, spectra.go). Absorption is the power the field puts in during the AC sweep "
        "(0–3700 cm<super>−1</super>, 2 cm<super>−1</super> steps, averaged over x, y, z). Release per mode comes by two "
        "routes that must agree: the Einstein A coefficient from the IR intensity, and Maxwell's radiated power of the "
        "vibrating dipole (em.LarmorPower, itself checked against the Poynting flux). Release spectra are given for one "
        "quantum per mode and for thermal equilibrium at 300 K (Kirchhoff).",
    ])

    r.h1("Results")
    r.h2("Bases and pairs (run_output.txt)")
    r.p("A = 9-methyladenine, T = 1-methylthymine, G = 9-methylguanine, C = 1-methylcytosine.")
    r.table([
        ["molecule", "atoms", "vibrations", "range [cm<super>−1</super>]", "RMS gradient [kcal/mol/Å]", "springs",
         "pairing sources", "d/λ", "thermal release at 300 K"],
        ["A", "18", "48", "49–3367", "7.5×10<super>−13</super>", "49", "0", "2.5×10<super>−4</super>", "2.58×10<super>−20</super> W"],
        ["T", "18", "48", "32–3314", "4.3×10<super>−10</super>", "48", "0", "2.1×10<super>−4</super>", "1.28×10<super>−20</super> W"],
        ["G", "19", "51", "30–3363", "2.1×10<super>−13</super>", "52", "0", "2.5×10<super>−4</super>", "2.90×10<super>−20</super> W"],
        ["C", "16", "42", "13–3366", "2.6×10<super>−13</super>", "41", "0", "2.2×10<super>−4</super>", "2.65×10<super>−20</super> W"],
        ["A·T pair", "36", "102", "13–3365", "4.5×10<super>−13</super>", "97", "6154", "4.4×10<super>−4</super>", "3.68×10<super>−20</super> W"],
        ["G·C pair", "35", "99", "18–3370", "3.6×10<super>−13</super>", "93", "5785", "4.2×10<super>−4</super>", "5.20×10<super>−20</super> W"],
    ], [0.11, 0.08, 0.11, 0.11, 0.14, 0.09, 0.1, 0.1, 0.16])
    r.p("Examples of the strongest absorption lines, which are also the strongest release lines: A at 3367 cm<super>−1</super> "
        "(N6–H stretches; 324 km/mol, A = 460 s<super>−1</super>); G at 717 cm<super>−1</super> (torsion/out-of-plane/nonbonded; "
        "539 km/mol, A = 35 s<super>−1</super>); G·C at 878 cm<super>−1</super> (581 km/mol, A = 56 s<super>−1</super>) and "
        "3368 cm<super>−1</super> (C N4–H41 stretch; 461 km/mol, A = 654 s<super>−1</super>). Here A is the Einstein coefficient "
        "of spontaneous release.")
    r.figure(fig("base_vibrations.png"),
             "Figure 1. Model lines of the four bases (bar height: strength of each vibration) with the laboratory "
             "measurements (triangles). Bands from left: ring bends out of plane; ring stretches, C=O and N–H bends; C–H "
             "and N–H stretches.", 0.8)

    r.page_break()
    r.h2("Model against measurement (summary/measured.csv)")
    r.table([
        ["band", "measured [cm<super>−1</super>]", "model [cm<super>−1</super>]", "difference", "measurement conditions"],
        ["A " + NH2 + " antisymmetric stretch", "3569", "3367.1", "−5.7 %", "9H-adenine, gas, jet-cooled IR–UV"],
        ["A " + NH2 + " symmetric stretch", "3452", "3246.8", "−5.9 %", "9H-adenine, gas, jet-cooled IR–UV"],
        ["T N3–H stretch", "3435", "3314.1", "−3.5 %", "1-methylthymine, gas"],
        ["T C4=O stretch", "1720", "1745.7", "+1.5 %", "thymine, gas, approximate"],
        ["G " + NH2 + " antisymmetric stretch", "3543", "3363.2", "−5.1 %", "9-methylguanine, Ar matrix 14 K"],
        ["G N1–H stretch", "3440", "3308.7", "−3.8 %", "9-methylguanine, Ar matrix"],
        ["G " + NH2 + " symmetric stretch", "3430", "3244.1", "−5.4 %", "9-methylguanine, Ar matrix"],
        ["G C6=O stretch", "1752", "1544.0", "−11.9 %", "9-methylguanine, Ar matrix"],
        ["C " + NH2 + " antisymmetric stretch", "3572.7", "3365.9", "−5.8 %", "cytosine keto, He nanodroplets 0.37 K"],
        ["C " + NH2 + " symmetric stretch", "3451.7", "3244.4", "−6.0 %", "cytosine keto, He nanodroplets"],
        ["C C2=O stretch", "1740", "1685.0", "−3.2 %", "cytosine keto, gas GC-FTIR 300 °C"],
        ["A·T inter-base H-bond stretch", "94", "123.5", "+31 %", "measured: A·T duplex in water; model: one pair in vacuum"],
    ], [0.27, 0.12, 0.11, 0.12, 0.38])
    r.figure(fig("model_vs_measured.png"),
             "Figure 2. Model (dots) against laboratory values (triangles) for the eleven base bands, with the relative "
             "difference.", 0.75)
    r.bullets([
        "The N–H stretches sit 3–6 % below the data: the harmonic AMBER value against anharmonic measurements.",
        "Measured red shifts of the N–H stretches on pairing (−39 to −186 cm<super>−1</super>, gas-phase A·T, a non-Watson–Crick "
        "isomer) are not reproduced. With fixed charges the model's stretches move up: A N6–H61 −1.8, T N3–H3 +27.6, "
        "G N1–H1 +21.5, G N2–H21 +126.3, C N4–H41 +2.2 cm<super>−1</super>.",
        "UV absorption and fluorescence maxima are listed as measured only; they are electronic and not in the vibration circuit.",
    ])

    r.h2("Pairing")
    r.table([
        ["pair", "hydrogen bond", "H···A [Å]", "D···A [Å]", "angle"],
        ["A·T", "A:N6–H61···T:O4", "1.87", "2.89", "178°"],
        ["A·T", "T:N3–H3···A:N1", "1.90", "2.91", "177°"],
        ["G·C", "G:N1–H1···C:N3", "1.92", "2.92", "173°"],
        ["G·C", "G:N2–H21···C:O2", "1.81", "2.82", "176°"],
        ["G·C", "C:N4–H41···G:O6", "1.82", "2.84", "180°"],
    ], [0.12, 0.38, 0.17, 0.17, 0.16])
    r.p("Donor–acceptor distances of 2.82–2.92 Å are near the crystal values (Saenger 1984). Modes with at least 25 % of their "
        "energy in the inter-base forces, with the share in the H-bonds: A·T 24 (45 %), 45 (47 %), 54 (74 %), 72 (82 %), "
        "123 (82 %) cm<super>−1</super>; G·C 62 (29 %), 63 (79 %), 116 (63 %), 154 (80 %) cm<super>−1</super>.")
    r.figure(fig("base_pairs.png"),
             "Figure 3. Top: minimized geometry of the Watson–Crick pairs (T–A and C–G are A–T and G–C seen from the other "
             "strand: same circuit and spectrum). Bottom: frequency response of each pair's circuit — absorbed power, release "
             "for one quantum per mode and thermal release at 300 K; insets show 0–250 cm<super>−1</super> with the H-bond "
             "modes dashed. Lines drawn 6 cm<super>−1</super> wide.", 1.0)

    r.h2("Stacks (stacks/summary*.csv)")
    r.table([
        ["stack", "pairs", "atoms", "vibrations", "lowest [cm<super>−1</super>]", "modes below 200 cm<super>−1</super>",
         "of which stacking", "twist start → minimized"],
        ["AT", "1", "36", "102", "12.7", "12", "0", "—"],
        ["AT2", "2", "72", "210", "11.6", "30", "5", "45.73° → 12.67°"],
        ["AT3", "3", "108", "318", "5.9", "48", "8", "43.64°, 42.85° → 20.97°, 34.19°"],
        ["AT4", "4", "144", "426", "3.4", "66", "14", "37.80°, 43.48°, 37.81° → 17.34°, 74.48°, 17.34°"],
        ["GC", "1", "35", "99", "17.7", "11", "0", "—"],
        ["GC2", "2", "70", "204", "10.9", "28", "5", "45.90° → 8.76°"],
        ["GC3", "3", "105", "309", "9.2", "45", "10", "43.08°, 41.56° → 65.98°, 18.21°"],
        ["GC4", "4", "140", "414", "5.4", "60", "15", "38.55°, 42.39°, 38.27° → 28.41°, 20.92°, 28.41°"],
        ["ApT", "2", "72", "210", "9.1", "30", "6", "47.54° → 34.54°"],
        ["polyA2", "2", "72", "210", "11.6", "30", "5", "45.85° → 12.67°"],
        ["polyA3", "3", "108", "318", "11.1", "48", "10", "43.51°, 42.06° → 15.76°, 11.11°"],
        ["polyA4", "4", "144", "426", "8.3", "66", "14", "38.30°, 43.73°, 36.62° → 19.22°, 13.23°, 10.50°"],
    ], [0.09, 0.07, 0.08, 0.11, 0.1, 0.13, 0.1, 0.32])
    r.p("Stacking turns each single-pair line into a band of nearby lines and adds rigid-pair (stacking) modes below "
        "200 cm<super>−1</super>. With the backbone absent, the rise and twist of the free stacks move away from the duplex's "
        "on minimization (for example AT2 twist 45.7° → 12.7°), which the summaries record.")
    r.figure(fig("stack_band_splitting.png"),
             "Figure 4. Ring/C=O stretch region (top) and N–H stretch region (bottom) for 1–4 stacked A·T pairs: each "
             "single-pair line splits into a cluster as pairs are added.", 0.8)
    r.figure(fig("stack_low_modes.png"),
             "Figure 5. Low-frequency vibrations of stacked A·T (top) and G·C (bottom) pairs. Thick lines: the pairs move as "
             "solid blocks against each other (stacking motions); thin lines: motion inside the pairs. Dotted: far-infrared "
             "bands of poly(dA)·poly(dT) (Powell et al. 1987); dashed: inter-strand H-bond mode of A·T DNA in water "
             "(González-Jiménez et al. 2016).", 0.8)

    r.h1("Validation")
    r.bullets([
        "Every molecule is at a true minimum and lumped by the λ/10 rule.",
        "The circuit's AC solution and a direct solve of the same admittance give the same absorbed power on and off "
        "resonance; the Einstein and Larmor release routes agree for every mode (dna_test.go).",
        "dna_test.go and ff_test.go also check the stiffness split, passivity, hydrogen bonds, the modal twin, the stack "
        "parts, and the analytic gradient against finite differences.",
        "In the pairs, 3610 (A·T) and 4088 (G·C) stiffness entries below 10<super>−9</super> of the largest were dropped "
        "(largest dropped 8.4×10<super>−9</super> and 8.8×10<super>−9</super> mdyn/Å).",
    ])

    r.h1("Limits of the model")
    r.bullets([
        "Harmonic force field: anharmonic shifts are absent, which is why N–H stretches come out 3–6 % low.",
        "Fixed charges: the model cannot give the measured red shifts of H-bonded stretches on pairing.",
        "Gas phase, methyl caps instead of sugars, no backbone and no water: free stacks unwind on minimization, and the "
        "single-pair H-bond mode is compared with a duplex measured in water.",
        "Electronic (UV) absorption and fluorescence are outside the vibration circuit.",
    ])

    r.h1("Reproducing the results")
    r.code("go run ./examples/dna_vibration [-out examples/dna_vibration] [-stacks=false] [-only-stacks] [-stack-keys AT2,GC2]\n"
           "go test ./examples/dna_vibration\n"
           "python examples/dna_vibration/viz/make_report.py examples/dna_vibration/dna_vibration_report.pdf")
    r.p("Outputs: netlists/ (each base and pair as a SPICE netlist), spectra/, modes/, summary/, stacks/ and sound/ (the "
        "spectra as audio, 1 cm<super>−1</super> → 3 Hz).")
    r.small("Sources: AMBER OL15 (Cornell et al., JACS 117, 5179, 1995; Pérez et al., Biophys. J. 92, 3817, 2007; Zgarbová "
            "et al., JCTC 11, 5723, 2015); Plützer et al., PCCP 3, 5466 (2001) and ChemPhysChem 4, 838 (2003); Nowak et al., "
            "PCCP 21, 22857 (2019); Min et al., Bull. Korean Chem. Soc. 30, 3039 (2009); Nir et al., PCCP 4, 732 (2002); "
            "González-Jiménez et al., Nat. Commun. 7, 11799 (2016); Cavaluzzi &amp; Borer, Nucleic Acids Res. 32, e13 "
            "(2004); Onidas et al., J. Phys. Chem. B 106, 11367 (2002); Saenger (1984), Ch. 6; Wilson, Decius &amp; Cross, "
            "Molecular Vibrations (1955).")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "dna_vibration_report.pdf"))
