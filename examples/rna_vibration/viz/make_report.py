# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""PDF report for examples/rna_vibration (RNA bases, Watson-Crick and wobble pairs as circuits).

python examples/rna_vibration/viz/make_report.py examples/rna_vibration/rna_vibration_report.pdf

Numbers come from examples/rna_vibration/run_output.txt (bases and pairs),
summary/measured.csv (model against measurement), stacks/summary.csv and
stacks/run_output.txt (stacks) and README.md (method). Figures are the PNGs
in summary/.
"""
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401

NH2 = "NH<sub>2</sub>"
CM = "cm<super>−1</super>"


def fig(name):
    return os.path.join(HERE, "..", "summary", name)


def sci(s):
    """'2.52e-20' -> 2.52×10<super>−20</super>."""
    m, e = s.lower().split("e")
    return "%s×10<super>%s</super>" % (m, str(int(e)).replace("-", "−"))


def build(out):
    r = Report("examples/rna_vibration", "RNA bases and base pairs as circuits",
               "Infrared absorption and release of the four RNA bases, the Watson–Crick pairs and the G·U wobble pair "
               "from kronspice's AC analysis — doc/chemistry-plan.md", "doc/chemistry-plan.md",
               "AMBER RNA OL3, harmonic, gas phase")

    r.h1("Summary")
    r.p("The RNA twin of dna_vibration. The four RNA bases (methyl-capped), the Watson–Crick pairs A·U and G·C and RNA's "
        "wobble pair G·U are built as electric circuits with the mobility analogy: atoms are capacitors (mass), every bond "
        "and angle of the AMBER force field is a spring inductor, and Kron's connection tensor joins the bond frame to the "
        "atom frame. Light enters as a current source at every atom. kronspice's AC analysis gives the infrared light each "
        "molecule absorbs, and two independent routes give the light it releases. Stacks of 2–4 pairs are placed on the "
        "A-form step of the NMR duplex 1EKA.")
    r.bullets([
        "<b>Starting structure:</b> PDB 1EKA, the solution NMR structure of r(GAGUGCUC)<sub>2</sub> with all hydrogens "
        "(Chen et al., Biochemistry 39, 8970, 2000): A·U = A2·U7, G·C = G1·C8, G·U = G5·U4; everything minimized in Go.",
        "<b>All seven molecules are true minima</b> (RMS gradient ≤ 7.0×10<super>−10</super> kcal/mol/Å, all 3N−6 "
        "vibrations real) and electrically small (d/λ = 2.0×10<super>−4</super> to 4.3×10<super>−4</super>).",
        "<b>Thermal infrared release at 300 K</b> is 1.24×10<super>−20</super> W (U) to 5.16×10<super>−20</super> W (G·C) "
        "per molecule; the wobble pair G·U gives 4.13×10<super>−20</super> W.",
        "<b>N–H stretches sit 3.4–6.0 % below the measurements</b> (harmonic AMBER against anharmonic data); uracil's "
        "C2=O comes out +3.0 %, its C4=O has no clean line (best mode 1458.0 " + CM + ", 19 % of the bond's energy, −14.2 %).",
        "<b>The wobble pair keeps its two hydrogen bonds</b> (G N1–H1···O2 U, U N3–H3···O6 G) and its stiffest pairing "
        "mode, 131 " + CM + ", sits beside A·U's 134 " + CM + " and below G·C's 152 " + CM + ".",
        "<b>Pairing red shifts are not reproduced:</b> with fixed charges the H-bonded stretches move up on pairing.",
    ])

    r.h1("The model")
    r.bullets([
        "<b>Molecules.</b> 9-methyladenine, 1-methyluracil, 9-methylguanine and 1-methylcytosine (the methyl stands in for "
        "the ribose) and the pairs A·U, G·C and G·U, cut from PDB 1EKA and minimized in the AMBER RNA force field OL3 "
        "(Cornell et al. base terms, bsc0, OL3 χ; OpenMM's amber14/RNA.OL3.xml). prep/prep.py writes molecules.json; the "
        "methyl caps are made exactly as in dna_vibration.",
        "<b>Circuit</b> (circuit.go, as in dna_vibration). Each AMBER bond and angle is a spring branch (VCVS chain for "
        "Wilson's B row, an inductor for the spring force, CCCS back to the atoms); everything else is the residual stiffness "
        "R = K − Σ k·BᵀB through displacement integrators and one VCCS per non-zero entry. Forces between the two bases of a "
        "pair form the group G<sub>HB</sub>; in a stack, forces between bases of different pairs form G<sub>ST</sub>.",
        "<b>Outputs.</b> Absorption is the power the field puts in during the AC sweep (0–3700 " + CM + ", 2 " + CM +
        " steps, averaged over x, y, z). Release per mode comes from the Einstein A coefficient and from Maxwell's radiated "
        "power of the vibrating dipole (em.LarmorPower); they must agree. Release spectra are given for one quantum per mode "
        "and for thermal equilibrium at 300 K.",
        "<b>Stacks.</b> 1EKA has no AA/UU, GG/CC or CG step, so stacks repeat its real GA/UC step (G1·C8 → A2·U7) in a "
        "simplified base-pair frame (Olson et al. 2001); U·A and C·G pairs are A·U and G·C turned 180° about the dyad. The "
        "tandem wobble stack GU2 is cut from 1EKA as it is.",
    ])

    r.h1("Results")
    r.h2("Bases and pairs (run_output.txt)")
    r.table([
        ["molecule", "atoms", "vibrations", "range [" + CM + "]", "RMS gradient [kcal/mol/Å]", "springs",
         "pairing sources", "d/λ", "thermal release at 300 K"],
        ["A", "18", "48", "48–3367", sci("3.8e-13"), "49", "0", sci("2.5e-04"), sci("2.52e-20") + " W"],
        ["U", "15", "39", "52–3313", sci("1.6e-13"), "39", "0", sci("2.0e-04"), sci("1.24e-20") + " W"],
        ["G", "19", "51", "31–3363", sci("5.0e-12"), "52", "0", sci("2.5e-04"), sci("2.97e-20") + " W"],
        ["C", "16", "42", "45–3366", sci("2.1e-13"), "41", "0", sci("2.2e-04"), sci("2.52e-20") + " W"],
        ["A·U pair", "33", "93", "12–3365", sci("7.5e-13"), "88", "5157", sci("4.3e-04"), sci("3.60e-20") + " W"],
        ["G·C pair", "35", "99", "18–3370", sci("4.7e-11"), "93", "5787", sci("4.2e-04"), sci("5.16e-20") + " W"],
        ["G·U pair", "34", "96", "12–3363", sci("7.0e-10"), "91", "5436", sci("4.3e-04"), sci("4.13e-20") + " W"],
    ], [0.11, 0.08, 0.11, 0.11, 0.14, 0.09, 0.1, 0.1, 0.16])
    r.p("A = 9-methyladenine, U = 1-methyluracil, G = 9-methylguanine, C = 1-methylcytosine. Against dna_vibration: the "
        "A, G and C N–H lines agree with the DNA run to within half a wavenumber (A " + NH2 + " antisymmetric 3367.0 against "
        "3367.1 " + CM + "); uracil releases 1.24×10<super>−20</super> W against thymine's 1.28×10<super>−20</super> W, and "
        "A·U 3.60×10<super>−20</super> W against A·T's 3.68×10<super>−20</super> W.")
    r.figure(fig("base_vibrations.png"),
             "Figure 1. Model lines of the four RNA bases (bar height: IR intensity) with the laboratory band positions "
             "(triangles).", 0.8)

    r.page_break()
    r.h2("Model against measurement (summary/measured.csv)")
    rows = [["band", "measured [" + CM + "]", "model [" + CM + "]", "difference", "conditions"]]
    short = {"Adenine": "A", "Uracil": "U", "Guanine": "G", "Cytosine": "C"}
    for m in csv.DictReader(open(os.path.join(HERE, "..", "summary", "measured.csv"), encoding="utf-8")):
        if not m["model_cm-1"]:
            continue
        cond = m["conditions"].split(";")[0]
        d = float(m["difference_pct"])
        rows.append([short[m["molecule"]] + " " + esc(m["band"]).replace("NH₂", NH2), m["measured_cm-1"], m["model_cm-1"],
                     ("+" if d > 0 else "−") + "%.1f %%" % abs(d), esc(cond)])
    r.table(rows, [0.3, 0.11, 0.1, 0.11, 0.38])
    r.p("Adenine, guanine and cytosine values and sources are reused from dna_vibration (the same molecules); the uracil "
        "values are new: 1-methyluracil in an argon matrix (Ivanov et al. 2002) and uracil in the gas phase (Colarusso et al. "
        "1997). No gas-phase or droplet IR spectrum of an A·U or G·U pair, and no measured low-frequency mode of an A-form "
        "RNA duplex, was found, so neither is compared.")
    r.figure(fig("model_vs_measured.png"),
             "Figure 2. Model (dots) against laboratory values (triangles), with the relative difference; uracil rows in red.", 0.8)

    r.h2("Pairing and the wobble pair")
    hb = list(csv.DictReader(open(os.path.join(HERE, "..", "summary", "hbonds.csv"), encoding="utf-8")))
    pname = {"AU": "A·U", "GC": "G·C", "GU": "G·U"}
    t = [["pair", "hydrogen bond", "H···A [Å]", "D···A [Å]", "angle"]]
    for h in hb:
        t.append([pname[h["pair"]], "%s–%s···%s" % (h["donor"], h["hydrogen"].split(":")[1], h["acceptor"]),
                  "%.2f" % float(h["h_acceptor_A"]), "%.2f" % float(h["donor_acceptor_A"]), "%.0f°" % float(h["angle_deg"])])
    r.table(t, [0.12, 0.38, 0.17, 0.17, 0.16])
    r.p("Modes with at least 25 % of their energy in the inter-base forces: A·U 28, 45, 53, 56, 72, 134 " + CM + "; G·C 61, "
        "62, 115, 152 " + CM + "; G·U 21, 43, 52, 55, 68, 131 " + CM + ". The wobble pair is held by the two wobble bonds only "
        "(G N1–H1···O2 U and U N3–H3···O6 G, a test); its strongest H-bond mode at 131 " + CM + " (79 % in the H-bonds) "
        "matches A·U's 134 " + CM + " (84 %) rather than G·C's 152 " + CM + " (80 %). On pairing, G N1–H1 moves +30.4 " + CM +
        " in G·U (+21.8 in G·C) and U N3–H3 +3.1 " + CM + " in G·U (+26.3 in A·U). Measured H-bonded N–H stretches move "
        "down; the fixed-charge model cannot give that red shift.")
    r.figure(fig("base_pairs.png"),
             "Figure 3. Top: minimized geometry of A·U, G·C and the G·U wobble pair, hydrogen bonds dashed. Bottom: frequency "
             "response of each pair's circuit — absorbed power, release for one quantum per mode and thermal release at 300 K; "
             "insets show 0–250 " + CM + " with the H-bond modes dashed. Lines drawn 6 " + CM + " wide.", 1.0)

    r.h2("Stacks (stacks/summary.csv)")
    t = [["stack", "pairs", "atoms", "vibrations", "lowest [" + CM + "]", "modes below 200 " + CM, "of which stacking",
          "twist start → minimized"]]
    for s in csv.DictReader(open(os.path.join(HERE, "..", "stacks", "summary.csv"), encoding="utf-8")):
        tw = "—"
        if s["twist_start_deg"]:
            tw = ", ".join(x + "°" for x in s["twist_start_deg"].split()) + " → " + ", ".join(x + "°" for x in s["twist_min_deg"].split())
        t.append([s["stack"], s["pairs"], s["atoms"], s["vibrations"], s["lowest_cm-1"], s["modes_below_200"],
                  s["stacking_modes_below_200"], tw])
    r.table(t, [0.09, 0.07, 0.08, 0.11, 0.1, 0.13, 0.1, 0.32])
    r.p(STACKS_TEXT)
    r.figure(fig("stack_band_splitting.png"),
             "Figure 4. Ring/C=O stretch region (top) and N–H stretch region (bottom) for 1–4 stacked A·U pairs: each "
             "single-pair line splits into a cluster as pairs are added.", 0.8)
    r.figure(fig("stack_low_modes.png"),
             "Figure 5. Low-frequency vibrations of stacked A·U (top) and G·C (bottom) pairs. Thick lines: the pairs move as "
             "solid blocks against each other; thin lines: motion inside the pairs. No measured A-form RNA band was found to "
             "mark.", 0.8)

    r.h1("Validation")
    r.bullets([
        "Every molecule is at a true minimum and lumped by the λ/10 rule.",
        "The circuit's AC solution and a direct solve of the same admittance give the same absorbed power on and off "
        "resonance; the Einstein and Larmor release routes agree for every mode (rna_test.go).",
        "rna_test.go and ff_test.go also check the stiffness split, passivity, hydrogen bonds (the wobble pair's two by "
        "name), the modal twin, the stack parts and the analytic gradient against finite differences.",
        "For the 2-pair stacks the modal spectrum agrees with the direct solve and with kronspice's AC analysis of the "
        "circuit (stacks/run_output.txt).",
    ])

    r.h1("Limits of the model")
    r.bullets([
        "Harmonic force field: anharmonic shifts are absent, which is why N–H stretches come out 3–6 % low.",
        "Fixed charges: the model cannot give the measured red shifts of H-bonded stretches on pairing.",
        "Gas phase, methyl caps instead of riboses, no backbone, no 2′-OH and no water: free stacks change their rise and "
        "twist on minimization.",
        "Stacks repeat one real step of 1EKA (GA/UC) for every sequence, a stated approximation; only GU2 has its own "
        "measured geometry.",
        "Electronic (UV) absorption and fluorescence are outside the vibration circuit.",
    ])

    r.h1("Reproducing the results")
    r.code("python examples/rna_vibration/prep/prep.py RNA.OL3.xml 1eka.pdb examples/rna_vibration/molecules.json\n"
           "go run ./examples/rna_vibration -out examples/rna_vibration -stacks=false\n"
           "go run ./examples/rna_vibration -out examples/rna_vibration -only-stacks\n"
           "go test ./examples/rna_vibration\n"
           "python examples/rna_vibration/viz/make_report.py examples/rna_vibration/rna_vibration_report.pdf")
    r.small("Sources: AMBER RNA OL3 (Cornell et al., JACS 117, 5179, 1995; Pérez et al., Biophys. J. 92, 3817, 2007; "
            "Zgarbová et al., JCTC 7, 2886, 2011); Chen et al., Biochemistry 39, 8970 (2000), PDB 1EKA; Ivanov et al., "
            "Biophysical Bulletin 2(11), 14 (2002); Colarusso et al., Chem. Phys. Lett. 269, 39 (1997); Plützer et al., PCCP "
            "3, 5466 (2001); Nowak et al., PCCP 21, 22857 (2019); Min et al., Bull. Korean Chem. Soc. 30, 3039 (2009); Nir et "
            "al., PCCP 4, 732 (2002); Cavaluzzi &amp; Borer, Nucleic Acids Res. 32, e13 (2004); Varani &amp; McClain, EMBO "
            "Rep. 1, 18 (2000); Olson et al., J. Mol. Biol. 313, 229 (2001); Wilson, Decius &amp; Cross, Molecular "
            "Vibrations (1955).")
    r.build(out)


STACKS_TEXT = ("Stacks: A·U (strand I AA, AAU, AAUU), G·C (CG, CGC, CGCG), poly(A)·poly(U) (AA, AAA, AAAA), ApU and the tandem "
               "wobble pairs GU2. Built stacks repeat 1EKA's GA/UC step, so AU2 and polyA2 are the same molecule. Stacking turns "
               "each single-pair line into a band and adds rigid-pair modes below 200 cm<super>−1</super>; the lowest mode falls "
               "from 12.5 cm<super>−1</super> (one A·U) to 6.4 (AAUU) and 4.8 cm<super>−1</super> (AAAA), against 12.7 → 3.4 "
               "cm<super>−1</super> for A·T in dna_vibration. With the backbone absent the free stacks leave the A-form step on "
               "minimization (AU2 rise 4.10 → 5.83 Å, GC2 twist 35.00° → 77.18°; simplified centroid and C1′–C1′ measures).")

if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "rna_vibration_report.pdf"))
