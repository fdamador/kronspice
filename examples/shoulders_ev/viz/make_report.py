# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""H27 report: Ken Shoulders' EVs tested with standard physics, as a PDF.

python examples/shoulders_ev/viz/make_report.py examples/shoulders_ev/shoulders_ev_report.pdf

Numbers come from examples/shoulders_ev/README.md (the example has no
run_output.txt); the figure is viz/shoulders_ev.png. Layout from
internal/pdfreport/kronreport.py.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401


def build(out):
    r = Report("examples/shoulders_ev", "Ken Shoulders' EVs",
               "Charge clusters of 10<super>11</super> electrons checked against Coulomb, Gauss, the Casimir force and "
               "charge and energy budgets — register entry H27",
               "register entry H27", "bench values declared")

    r.h1("Summary")
    r.p("K. R. Shoulders (<i>EV: A Tale of Discovery</i>, 1987) reports micrometre beads of electric charge thrown off "
        "by fast discharges from sharp cathodes, which he calls EVs (<i>electrum validum</i>). He reports about "
        "10<super>11</super> electrons per 1 µm bead, near-solid electron density with almost no ions, motion at about "
        "c/10, and he proposes that the zero-point fluctuations of empty space hold them together. This example takes "
        "his numbers at face value and checks them against Coulomb's and Gauss's laws, the known containers, and the "
        "charge and energy budgets of his sources; then it runs his source and an EV strike as kronspice circuits "
        "with energy ledgers. Page numbers below are the book's chapter-page numbers.")
    r.bullets([
        "<b>Energy:</b> a 1 µm bead of 10<super>11</super> electrons needs <b>2.77 J</b> to assemble, 3 million times "
        "his smallest source's 1 µJ and 1,400 times his larger source's 2.0 mJ.",
        "<b>Left bare it explodes:</b> outward pressure 1.5×10<super>18</super> Pa; the radius doubles in "
        "1.7×10<super>−15</super> s, 20,000 times shorter than his reported lifetime of 3.3×10<super>−11</super> s.",
        "<b>No known container:</b> the Casimir pull (1.3×10<super>−3</super> Pa) is 10<super>21</super> times too weak, "
        "and on a sphere the zero-point stress pushes outward; a magnetic trap would need 2×10<super>5</super> T.",
        "<b>Gauss's law:</b> his kilovolt readings need a “container” of at least 14.4 cm radius (1 kV), or a net "
        "charge of at most about 3×10<super>5</super> electrons.",
        "<b>Charge budget:</b> the 1 pF source holds 8.8×10<super>9</super> electrons; one EV bead is 11 times that. "
        "In kronspice the source puts exactly C·V (8.83×10<super>9</super> electrons) through the scope, ledger "
        "closed to 3×10<super>−5</super>.",
        "<b>Test:</b> a one-shot Faraday-cup charge measurement against C·V of the source. Standard physics predicts "
        "it never exceeds C·V.",
    ])

    r.h1("The claim and the standard prediction")
    r.p("<b>Hypothesis (H27), from Shoulders 1987:</b>")
    r.bullets([
        "<b>Size and charge:</b> about 1 µm across, often in chains and rings; one bead carries “around "
        "10<super>11</super> charges” (8-49).",
        "<b>Content:</b> electrons at a density near that of a solid (6.6×10<super>23</super> per cm³), with fewer than "
        "1 ion per 100,000 electrons (4-63).",
        "<b>Motion and life:</b> about c/10, living about 10<super>−11</super> s; they bore pits, leave ring-shaped "
        "marks and give X-rays where they hit.",
        "<b>Binding:</b> the zero-point fluctuations of empty space (Preface; 9-2).",
        "<b>Source:</b> a cathode of about 1 pF storing less than 10<super>−6</super> J (2-7).",
    ])
    r.p("<b>Standard prediction.</b> A bead of 10<super>11</super> electrons in 1 µm is held by nothing known and "
        "explodes in femtoseconds; the charge collected from one shot can never exceed C·V of the source, and the "
        "energy delivered never exceeds what the source stored. The measured counterpart is G. A. Mesyats' ectons: "
        "explosive-emission micro-explosions at cathode spots.")

    r.h1("1. A bead of 10<super>11</super> electrons, 1 µm across")
    r.table([
        ["quantity", "value"],
        ["charge", "1.6×10<super>−8</super> C"],
        ["density", "1.9×10<super>23</super> electrons/cm³ (his 3.5×10<super>14</super> electrons in 10 µm gives his "
         "6.6×10<super>23</super>)"],
        ["energy needed to assemble it", "<b>2.77 J</b>"],
        ["field at its surface", "5.8×10<super>14</super> V/m"],
        ["outward pressure", "1.5×10<super>18</super> Pa (10<super>13</super> atmospheres)"],
        ["kinetic energy of all its electrons at c/10", "4.1×10<super>−5</super> J (the field holds 67,000 times more)"],
    ], [0.4, 0.6])
    r.p("<b>Left bare, it flies apart.</b> The outer electrons are pushed off at nearly the speed of light: the radius "
        "doubles in 1.7×10<super>−15</super> s and grows tenfold in 1.5×10<super>−14</super> s, and each outer electron "
        "leaves with about 260 MeV. Shoulders' lifetime (1 mm of track at c/10) is 3.3×10<super>−11</super> s, 20,000 "
        "times longer than the doubling time. Something would have to hold it.")

    r.h1("2. What could hold it?")
    r.table([
        ["candidate", "pressure", "against the 1.5×10<super>18</super> Pa push"],
        ["magnetic attraction of electrons moving together at c/10", "1 % of the repulsion", "the push still wins"],
        ["a magnetic trap at this density (Brillouin limit)", "needs 2×10<super>5</super> T",
         "the strongest field made indoors is about 1,200 T"],
        ["Casimir pull between ideal plates 1 µm apart", "1.3×10<super>−3</super> Pa", "10<super>21</super> times too weak"],
        ["zero-point stress on a 1 µm conducting sphere", "1.9×10<super>−3</super> Pa, <b>outward</b>",
         "pushes the wrong way (Boyer 1968)"],
        ["positive ions", "would cancel the charge", "his own measurement: under 1 ion per 100,000 electrons"],
    ], [0.4, 0.25, 0.35])
    r.p("The zero-point field is the container Shoulders proposes. Its known, measured effect, the Casimir force, is 21 "
        "orders of magnitude too weak, and for a sphere it pushes outward.")

    r.h1("3. Gauss's law and his kilovolt measurements")
    r.p("Shoulders states the puzzle himself (9-2): an electron leaving a bare 1 µm bead of 10<super>11</super> "
        "electrons gains 288 MeV, but his probes saw “only kilovolts at best”. His answer is a container that hides "
        "the core's field. Gauss's law rules this out: outside any closed surface, the potential of a net charge q is "
        "k·q/r whatever is inside. For escaping electrons to gain only V volts the “container” must be at least:")
    r.table([
        ["escape energy seen", "minimum container radius"],
        ["1 kV", "14.4 cm"],
        ["10 kV", "1.4 cm"],
    ], [0.5, 0.5])
    r.p("Those sizes are larger than the millimetre gaps of his apparatus. The alternative is a net charge of at most "
        "about 3×10<super>5</super> electrons for a 1 µm bead seen at 1 kV, 300,000 times less than "
        "10<super>11</super>, unless positive charge cancels it, which his ion measurement excludes.")

    r.h1("4. Charge and energy budgets")
    r.table([
        ["source (2-7)", "electrons it holds (C·V)", "energy"],
        ["smallest: 1 pF, &lt; 1 µJ (at most 1.4 kV)", "8.8×10<super>9</super>", "1.0 µJ"],
        ["earlier: 40 pF at 10 kV", "2.5×10<super>12</super>", "2.0 mJ"],
    ], [0.45, 0.3, 0.25])
    r.table([
        ["what he reports", "electrons", "÷ small source", "÷ larger source"],
        ["residual charge of a 3 µm EV on an electrode (4-63)", "2×10<super>10</super>", "2.3", "0.008"],
        ["one EV bead (8-49)", "10<super>11</super>", "11", "0.04"],
        ["emitted by a 10 µm EV over 1 mm of track (4-9)", "3.5×10<super>14</super>", "40,000", "<b>140</b>"],
    ], [0.5, 0.17, 0.16, 0.17])
    r.bullets([
        "<b>Charge</b> is conserved exactly. A 10<super>11</super>-electron bead cannot come from the 1 pF source, and "
        "the 3.5×10<super>14</super> electrons he infers from camera brightness cannot come from either source. That "
        "count is probably the inference to check first: it rests on a phosphor calibration and on isotropic emission.",
        "<b>Energy:</b> the 2.77 J assembly energy is 3 million times the small source's energy and 1,400 times the "
        "larger one's.",
    ])

    r.h1("5. kronspice: the source and a strike, with ledgers")
    r.p("The bench (declared values): a 1 pF anode, the witness plate or a collector, feeding a 50 Ω scope line.")
    r.h2("a. His smallest source")
    r.p("A 1 pF cathode at 1.414 kV (1 µJ) discharges across the gap (a 20 Ω arc with 5 nH).")
    r.table([
        ["quantity", "value"],
        ["charge through the scope", "8.83×10<super>9</super> electrons (all of C·V)"],
        ["energy", "1.000 µJ stored → 0.638 µJ to the scope, 0.362 µJ heat in the arc"],
        ["ledger", "closes to 3×10<super>−5</super>"],
    ], [0.35, 0.65])
    r.h2("b. An EV strike")
    r.p("An ideal current source brings the EV's charge to the anode in 20 ps.")
    r.table([
        ["electrons", "peak at the scope", "energy the arriving charge must bring", "÷ source energy"],
        ["2×10<super>10</super>", "2.7 kV", "4.7 µJ", "5"],
        ["10<super>11</super>", "13.4 kV", "117 µJ", "117"],
    ], [0.2, 0.2, 0.4, 0.2])
    r.p("The ledger closes to 10<super>−5</super> in both.")
    r.figure(os.path.join(HERE, "shoulders_ev.png"),
             "Figure 1. Left: outward pressure of a 1 µm EV of 10<super>11</super> electrons against the magnetic pinch "
             "at c/10, the pressure of the record 1200 T magnet and the Casimir pull between plates 1 µm apart (log "
             "scale). Middle: electrons his sources hold (green) against what he reports (orange). Right: kronspice "
             "scope voltage on the 1 pF anode into 50 Ω for the 1 pF, 1.4 kV source discharge and for EV strikes of "
             "2×10<super>10</super> and 10<super>11</super> electrons.")

    r.h1("6. The rings, and the mainstream counterpart")
    r.bullets([
        "<b>Ring marks.</b> Ring-shaped witness marks (2-19) fit a toroidal current ring. A real electron ring needs an "
        "outside field to hold it: the collapsing ring in <i>examples/plasmoid_rings</i> needs a 2 T coil, and its "
        "electrons are 10<super>6</super> times sparser than an EV's.",
        "<b>Ectons.</b> Mesyats' ectons are explosive-emission micro-explosions at cathode spots, each lasting "
        "nanoseconds and releasing about 10<super>11</super>–10<super>12</super> electrons. They leave craters and fast "
        "current pulses and give X-rays.",
        "<b>The difference.</b> An ecton is a neutral plasma of cathode metal, with ions. Shoulders reports almost no "
        "ions, so his ion measurement is where the two pictures disagree. It is worth repeating.",
    ])

    r.h1("Validation")
    r.bullets([
        "<i>main_test.go</i> checks the electrostatics, the Coulomb explosion against the exact solution, the "
        "containment numbers, and the kronspice charge and ledger.",
        "kronspice ledgers close to 3×10<super>−5</super> (source discharge) and 10<super>−5</super> (strikes); the "
        "source delivers exactly its C·V.",
    ])

    r.h1("Limits of the model")
    r.bullets([
        "Shoulders' reported numbers are taken at face value; the bench (1 pF anode, 50 Ω scope, 20 Ω arc with 5 nH) "
        "is declared.",
        "The EV strike is an ideal current source delivering the charge in 20 ps; it shows the energy the charge "
        "must bring, not a model of how it is made.",
        "The observations themselves (pits, rings, chains, X-rays, very fast pulses from tiny discharges) are real and "
        "carefully documented; what is tested here is the EV interpretation.",
    ])

    r.h1("Status")
    r.p("At the reported size and charge an EV holds 1,400 to 3 million times its source's energy, flies apart in "
        "femtoseconds, shows a field at odds with Gauss's law and his own kilovolt measurements, and has no known "
        "container. What would change this: a one-shot Faraday-cup charge measurement exceeding C·V of the source, "
        "with this ledger on the shot. That is the H27 test. Shoulders' later patents (5,018,180; 5,054,046; "
        "5,054,047; 5,123,039; 5,148,461) and the Puthoff &amp; Piestrup review carry the same claims.")

    r.h1("Reproducing the results")
    r.code("go run ./examples/shoulders_ev -viz examples/shoulders_ev/viz   # under a second\n"
           "python3 examples/shoulders_ev/viz/render.py\n"
           "go test ./examples/shoulders_ev")
    r.small("Sources: K. R. Shoulders, <i>EV: A Tale of Discovery</i> (1987); H. E. Puthoff, M. A. Piestrup, "
            "arXiv:physics/0408114 (2004); T. H. Boyer, Phys. Rev. 174, 1764 (1968); H. B. G. Casimir, Proc. K. Ned. "
            "Akad. Wet. 51, 793 (1948); L. Brillouin, Phys. Rev. 67, 260 (1945); G. A. Mesyats, Plasma Phys. Control. "
            "Fusion 47, A109 (2005); D. Nakamura et al., Rev. Sci. Instrum. 89, 095106 (2018).")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "shoulders_ev_report.pdf"))
