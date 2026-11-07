# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""H25 report: Forward's vacuum-fluctuation battery as a Kron circuit, as a PDF.

python examples/forward_casimir/viz/make_report.py examples/forward_casimir/forward_casimir_report.pdf

Numbers come from examples/forward_casimir/run_output.txt and README.md; the
figure is viz/forward_casimir.png. Layout from internal/pdfreport/kronreport.py.
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
    r = Report("examples/forward_casimir", "Forward's vacuum-fluctuation battery",
               "Can charged leaves pulled together by the Casimir force deliver electrical energy, and is it a "
               "continuous source? — register entry H25", "register entry H25",
               "100 nm Al leaves, 1 µm → 20 nm at 1 mm/s")

    r.h1("Summary")
    r.p("R. L. Forward (<i>Phys. Rev. B</i> 30, 1700, 1984) proposed a battery driven by the Casimir force, the "
        "attraction between close metal plates that comes from the electromagnetic zero-point fluctuations of empty "
        "space. A stack of thin aluminium leaves carries equal charges, so their electric repulsion pushes them "
        "apart; a fast feedback (servo) supply holds that repulsion just below the Casimir pull, so the leaves creep "
        "together, and the work the Casimir force does goes into the supply as electrical energy. Pushing the leaves "
        "apart again “recharges” it. This example computes the pull for real aluminium, models one leaf pair as a "
        "Kron (q, a) system with a full ledger, twins its electrical side in kronspice, and runs a full cycle.")
    r.bullets([
        "<b>Real aluminium pulls less than perfect mirrors:</b> at 20 nm the Lifshitz energy is 0.230 of Forward's "
        "ideal value (1.245×10<super>−5</super> against 5.417×10<super>−5</super> J/m²).",
        "<b>One closing stroke gains energy, once:</b> 12.3 µJ/m² for aluminium (54.0 µJ/m² for perfect mirrors), "
        "equal to the Casimir work less servo losses.",
        "<b>The supply pays first:</b> during the stroke it puts in 4.22 mJ/m² (hundreds of times the Casimir work), "
        "stored in the field; the gain appears only when the leaves are stopped and discharged.",
        "<b>A full cycle is a net loss:</b> reopening costs 12.6 µJ/m², so close + reopen gives −0.28 µJ/m² "
        "(aluminium). The Casimir force is conservative.",
        "<b>kronspice twin agrees:</b> supply energy 4.21981 mJ against 4.22084 mJ from the model's own solution; the "
        "leaves' work on the capacitor equals the Casimir work (12.4543 against 12.4542 µJ).",
        "<b>Scale:</b> 0.023 J per kg of leaf metal per stroke, against about 900,000 J/kg for a lithium-ion cell.",
    ])

    r.h1("What is being tested")
    r.p("<b>Hypothesis (H25).</b> Electrical energy can be extracted from the vacuum by the Casimir cohesion of "
        "charged foliated conductors (Forward 1984). Forward himself compared it to a hydroelectric dam: the force is "
        "a catalyst that releases stored energy, here the surface energy of the separated leaves.")
    r.p("<b>Standard prediction.</b> The Casimir force is conservative, so one stroke releases stored surface energy "
        "once, consistent with energy conservation (Cole &amp; Puthoff, <i>Phys. Rev. E</i> 48, 1562, 1993); "
        "reopening the leaves costs at least what closing gave. The energy cannot exceed the leaves' surface energy, "
        "about 1.1 J/m² per aluminium surface (Tyson &amp; Miller 1977).")

    r.h1("1. How strong is the pull?")
    r.p("Forward used the formulas for perfect mirrors (his eqs. 1–2). Real aluminium stops reflecting above its "
        "plasma frequency, so the Lifshitz formula is used with aluminium's plasma wavelength λ<sub>p</sub> = 107 nm "
        "(Lambrecht &amp; Reynaud 2000). With λ<sub>p</sub> set to zero the calculation reproduces Forward's eq. 1 to "
        "10<super>−5</super>. Energies and pressures are per m² of plate; negative values mean attraction.")
    r.table([
        ["gap", "energy, ideal", "energy, Al", "Al ÷ ideal", "pressure, ideal", "pressure, Al"],
        ["10 nm", "−4.334×10<super>−4</super> J", "−5.861×10<super>−5</super> J", "0.135",
         "−1.300×10<super>5</super> Pa", "−1.284×10<super>4</super> Pa"],
        ["20 nm", "−5.417×10<super>−5</super> J", "−1.245×10<super>−5</super> J", "0.230",
         "−8.126×10<super>3</super> Pa", "−1.421×10<super>3</super> Pa"],
        ["100 nm", "−4.334×10<super>−7</super> J", "−2.547×10<super>−7</super> J", "0.588",
         "−1.300×10<super>1</super> Pa", "−6.586 Pa"],
        ["1000 nm", "−4.334×10<super>−10</super> J", "−4.056×10<super>−10</super> J", "0.936",
         "−1.300×10<super>−3</super> Pa", "−1.191×10<super>−3</super> Pa"],
    ], [0.12, 0.18, 0.18, 0.12, 0.2, 0.2])
    r.p("Below λ/2π the force turns to the van der Waals form (Forward eq. 4); the plasma model includes it.")

    r.h1("2. One leaf pair as a Kron circuit")
    r.p("Kron's coordinates are the charge q on the leaves and the gap a, as in his machines (Hoffmann 1949; see "
        "<i>hoffmann_motor</i>). Both leaves sit at the supply voltage V, so there is no field in the gap; each faces a "
        "grounded neighbour at distance L, giving C(a) = 4ε<sub>0</sub>/(L − a) per m². The field outside pushes the "
        "leaves apart: Forward's repulsion. The equations are Lagrange–Maxwell with Casimir's energy as a potential "
        "(r: supply resistance, m: leaf mass per m²):")
    r.code("r·dq/dt + q/C(a) = V\n½m·d²a/dt² = (electric push) − (Casimir pull)")
    r.p("The servo sets V at every instant so that the push equals the pull, plus just enough to move the leaves at a "
        "steady 1 mm/s. Declared values: 100 nm aluminium leaves, 2 µm to the grounded neighbours, closing from 1 µm "
        "to 20 nm.")
    r.table([
        ["per m² of gap", "perfect mirrors", "aluminium"],
        ["work done by the Casimir force", "5.4173×10<super>−5</super> J", "1.2454×10<super>−5</super> J"],
        ["servo voltage at 20 nm", "63.8 V", "26.7 V"],
        ["energy the supply gives during the stroke", "2.4162×10<super>−2</super> J", "4.2208×10<super>−3</super> J"],
        ["energy stored in the field at 20 nm", "2.4216×10<super>−2</super> J", "4.2332×10<super>−3</super> J"],
        ["<b>net gain after discharging at the stop</b>", "<b>5.4023×10<super>−5</super> J</b>",
         "<b>1.2315×10<super>−5</super> J</b>"],
        ["losses: heat; motion; ledger residual", "2.9×10<super>−8</super>; 1.5×10<super>−10</super>; 1.2×10<super>−7</super> J",
         "2.4×10<super>−8</super>; 7.8×10<super>−11</super>; 1.2×10<super>−7</super> J"],
        ["charge sent during the stroke", "7.587×10<super>−4</super> C", "3.172×10<super>−4</super> C"],
        ["electrons returned at discharge", "4.74×10<super>15</super>", "1.98×10<super>15</super>"],
        ["Casimir energy carried out per electron", "0.071 eV", "0.039 eV"],
    ], [0.44, 0.28, 0.28], bold_rows=(5,))
    r.bullets([
        "<b>During the stroke the supply pays.</b> The Casimir pull grows as 1/a<super>4</super>, so the servo must raise "
        "the voltage steeply; the supply puts in hundreds of times the Casimir work, and it is stored in the field.",
        "<b>The gain appears only at the end.</b> When the leaves are stopped and discharged, all the charge comes back "
        "(net charge moved over the whole operation: zero), each electron carrying a few hundredths of an eV more "
        "than it went in with. Net, the supply gains the Casimir work less tiny servo losses.",
        "<b>The ledger closes</b> to about 1 % of the Casimir work and 3×10<super>−5</super> of the energy flows "
        "(residual 1.2×10<super>−7</super> J/m²).",
    ])

    r.h1("3. The same electrical side in kronspice")
    r.p("kronspice's varying capacitor C(t) = 4ε<sub>0</sub>/(L − a(t)) follows the gap of the stroke and is driven by "
        "the servo voltage through the supply resistance (1 µΩ·m²).")
    r.table([
        ["", "kronspice", "own solution", "ledger mismatch"],
        ["perfect mirrors: supply energy", "2.41528×10<super>−2</super> J", "2.41615×10<super>−2</super> J",
         "4×10<super>−4</super>"],
        ["perfect mirrors: charge sent", "7.58463×10<super>−4</super> C", "7.58679×10<super>−4</super> C", ""],
        ["perfect mirrors: leaves' work on the capacitor", "5.41731×10<super>−5</super> J",
         "Casimir 5.41731×10<super>−5</super> J", ""],
        ["aluminium: supply energy", "4.21981×10<super>−3</super> J", "4.22084×10<super>−3</super> J",
         "2×10<super>−4</super>"],
        ["aluminium: charge sent", "3.17151×10<super>−4</super> C", "3.17207×10<super>−4</super> C", ""],
        ["aluminium: leaves' work on the capacitor", "1.24543×10<super>−5</super> J",
         "Casimir 1.24542×10<super>−5</super> J", ""],
    ], [0.4, 0.2, 0.24, 0.16])

    r.h1("4. A full cycle: close, then reopen")
    r.table([
        ["per m²", "closing gain", "reopening cost", "cycle"],
        ["perfect mirrors", "+5.4023×10<super>−5</super> J", "−5.4311×10<super>−5</super> J",
         "<b>−2.884×10<super>−7</super> J</b>"],
        ["aluminium", "+1.2315×10<super>−5</super> J", "−1.2592×10<super>−5</super> J",
         "<b>−2.771×10<super>−7</super> J</b>"],
    ], [0.25, 0.25, 0.25, 0.25])
    r.p("Reopening costs what closing gave, and the servo's losses make every cycle a small net loss. Forward's "
        "battery is a one-shot release of stored surface energy, like a spring — not a continuous source.")
    r.figure(os.path.join(HERE, "forward_casimir.png"),
             "Figure 1. Left: attractive Casimir pressure against gap for perfect mirrors (Forward eq. 2, grey) and "
             "aluminium (Lifshitz, λ<sub>p</sub> = 107 nm, blue); real metal pulls less at small gaps. Right: the "
             "servo voltage needed to balance the pull, closing (solid) and reopening (dashed), for perfect mirrors "
             "and aluminium.")

    r.h1("5. Scale, and electron pumps")
    r.bullets([
        "<b>Energy per kilogram.</b> One stroke of real aluminium yields 1.23×10<super>−5</super> J per m² of gap, "
        "0.023 J per kg of leaf metal (100 nm leaves, two per gap). A lithium-ion cell stores about 900,000 J/kg.",
        "<b>The upper limit.</b> Pressing the leaves into full contact cannot release more than their surface energy, "
        "about 1.1 J/m² per aluminium surface (Tyson &amp; Miller 1977) — the same energy that had to be spent to make "
        "the separate leaves.",
        "<b>As an electron pump.</b> One stroke is a variable-capacitor charge pump driven by the Casimir force: it "
        "lifts 2.0×10<super>15</super> electrons per m² by 0.039 eV each, once; reopening costs that back. A "
        "continuous pump needs a continuous source of work, as the pumps of <i>photon_electron</i> have.",
        "<b>What it is good for.</b> The servo-controlled variable capacitor itself, and the Casimir force as a "
        "nanoscale actuator. It is not a source of energy.",
        "<b>Measurements.</b> The Casimir force between parallel metal plates has been measured (Bressi et al., "
        "<i>PRL</i> 88, 041804, 2002, 0.5–3 µm, 15 %). Forward noted that no measurements existed for metal plates "
        "at the 1–20 nm spacings his battery needs.",
    ])

    r.h1("Validation and limits")
    r.p("<i>main_test.go</i> checks the perfect-mirror limit (Forward eq. 1 to 10<super>−5</super>), the stroke "
        "ledger and its kronspice twin, and the cycle. Scope: one leaf pair; aluminium in the plasma model "
        "(λ<sub>p</sub> = 107 nm); a servo moving the leaves at a steady 1 mm/s; declared leaf thickness and spacing.")

    r.h1("Status")
    r.p("H25: a one-shot release of surface (Casimir) energy, consistent with energy conservation (Cole &amp; Puthoff "
        "1993); 12 µJ/m² per aluminium stroke, 0.02 J/kg; a close–reopen cycle is a net loss (computed) — not a "
        "continuous source. The ledger closes: what the supply gains is the leaves' surface energy.")

    r.h1("Reproducing the results")
    r.code("go run ./examples/forward_casimir -viz examples/forward_casimir/viz   # about 3 seconds\n"
           "python3 examples/forward_casimir/viz/render.py\n"
           "go test ./examples/forward_casimir\n"
           "python examples/forward_casimir/viz/make_report.py examples/forward_casimir/forward_casimir_report.pdf")
    r.small("Sources: R. L. Forward, Phys. Rev. B 30, 1700 (1984); H. B. G. Casimir, Proc. K. Ned. Akad. Wet. 51, 793 "
            "(1948); E. M. Lifshitz, Sov. Phys. JETP 2, 73 (1956); A. Lambrecht, S. Reynaud, Eur. Phys. J. D 8, 309 "
            "(2000); D. C. Cole, H. E. Puthoff, Phys. Rev. E 48, 1562 (1993); G. Bressi et al., PRL 88, 041804 (2002); "
            "W. R. Tyson, W. A. Miller, Surface Science 62, 267 (1977).")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "forward_casimir_report.pdf"))
