# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""H26 report: Nachamkin's force-free time-harmonic plasmoids, harmonics and size.

python examples/Force-Free_Harmonic_Plasmoids/viz/make_report.py

Numbers come from run_output.txt (the program's report); the figures are
made by viz/render.py from the CSVs the program writes with -viz. Layout
from internal/pdfreport/kronreport.py.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report  # noqa: E402
from reportlab.platypus import KeepTogether  # noqa: E402


def build(out):
    r = Report("examples/Force-Free_Harmonic_Plasmoids", "Force-free time-harmonic plasmoids: harmonics and size",
               "Nachamkin's ball-lightning plasmoid studied as written, over the reported range of size, energy and lifetime",
               "register entry H26", "about 5 seconds (plus 25 s for the pictures); deterministic")

    r.h1("Summary")
    r.p("J. Nachamkin's 1992 Phillips Laboratory report (PL-TR-92-3044) models ball lightning as a sphere of plasma "
        "holding a standing electromagnetic wave whose electric and magnetic fields are parallel everywhere (a "
        "force-free field), with a vortical flow balancing its pressure and a \"resonant size\" at which it should not "
        "radiate. This example takes the paper as written and asks three questions at ball-lightning scale: which "
        "frequencies (harmonics) a plasmoid of each size may have, what the paper's own recipe gives over the reported "
        "range of size and energy, and how long such a plasmoid lasts. kronspice analyses its spectrum, its ring-down "
        "and its charge condition as circuits.")
    r.bullets([
        "<b>The field is sound.</b> Eq. 8 is force-free to 1.8×10<super>−10</super>; its energy integral is 35.868/α³, as the paper's B3 says.",
        "<b>Appendix B has an arithmetic slip.</b> Its field E<sub>0</sub> = 3.566×10<super>9</super> V/m follows from 35.87/α², "
        "not from its own B3 (α³): the field should be √α = 6.7 times larger. Every other printed number is reproduced.",
        "<b>Each size allows a comb of harmonics</b> f<sub>n</sub> = x<sub>n</sub>·c/2πR (x<sub>n</sub> the zeros of j<sub>1</sub>); at "
        "the mean ball-lightning radius, 12 cm, they start at 3.07 GHz.",
        "<b>The paper's recipe drives its electrons faster than light</b>, v/c = ξ·kR = 1.4 to 12 c, over the whole reported "
        "range; keeping them below 0.1 c moves the plasmoid into the terahertz.",
        "<b>Free, it lasts microseconds to milliseconds, not seconds.</b> Lasting the 1–10 s that witnesses report needs, at a "
        "12 cm radius, 7–74 free electrons per sea-level air molecule; fully stripped air has 14.5.",
        "<b>kronspice</b> returns the plasmoid's eight lowest resonances to 10<super>−15</super> in frequency and 2×10<super>−6</super> in Q, "
        "closes the ring-down ledger to 5×10<super>−6</super>, and shows the paper's charge condition (eq. 49) is the resonance "
        "of a plasma-loaded line.",
    ])

    r.h1("What is being tested")
    r.h2("The claim (hypothesis H26)")
    r.p("Long-lived plasmoids (ball lightning) are force-free, time-harmonic electromagnetic fields in a plasma, "
        "stabilized by vortical flow and non-radiating at a resonant size: j<sub>1</sub>(αR) = 0 and j<sub>1</sub>(kR) = 0 "
        "with k > α (eq. 35). The free-electron density is fixed by eq. 49, ρ<sub>e</sub> = (m/e)ε<sub>0</sub>c²(k² − α²).")
    r.h2("The standard prediction beside it")
    r.p("Eq. 49 is the cold-plasma dispersion relation ω² = ω<sub>p</sub>² + α²c². A plasma sphere fed by nothing outside "
        "has only leaky (radiating) modes; its TE and TM halves are separate modes needing different sizes, so the "
        "force-free mix of the two is not one free mode. Ball lightning as reported: diameter 10–50 cm (mean 24 cm), "
        "lifetime peaking at 2–5 s (mean about 10 s), energy 10 kJ to over 1 MJ (eyewitness statistics).")

    r.page_break()
    r.h1("The plasmoid")
    r.figure(os.path.join(HERE, "plasmoid_structure.png"),
             "Figure 1. Left: the l = 1 field of eq. 8 lies on nested tori (three shown, one field line on each); E and H "
             "follow the same lines. Middle: a cut through the axis, energy density |M|² and the tori, nesting around the "
             "current ring at 0.611 R. Right: over a cycle the energy moves between the electric field and the electrons' "
             "motion (plus a magnetic share (α/k)², negligible at Appendix B's k/α = 59.8); the total stays constant "
             "because ω<sub>p</sub>² + α²c² = ω².")
    r.p("Eq. 8 for l = 1: M = 2 cos θ · j<sub>1</sub>(αr)/(αr) e<sub>r</sub> − [j<sub>0</sub>(αr) − j<sub>1</sub>(αr)/(αr)] sin θ e<sub>θ</sub> ± "
        "j<sub>1</sub>(αr) sin θ e<sub>φ</sub> (unit vectors e), with E = E<sub>0</sub>M sin(ωt + φ) and H proportional to M cos(ωt + φ).")

    r.h1("Results")
    r.h2("1. Appendix B recomputed")
    r.table([
        ["Case", "E<sub>0</sub> (V/m)", "kR used", "f (GHz)", "ρ<sub>e</sub> (C/m³)", "Electrons", "v/c"],
        ["As printed, his k = 2686 /m", "3.566×10<super>9</super>", "268.60", "128.16", "−32.63", "8.532×10<super>17</super>", "2.69"],
        ["As printed, nearest zero", "3.566×10<super>9</super>", "265.46", "126.66", "−31.87", "8.333×10<super>17</super>", "2.65"],
        ["B3 as printed (α³)", "2.39×10<super>10</super>", "683.30", "326.0", "−211.2", "5.52×10<super>18</super>", "6.83"],
    ], [0.27, 0.13, 0.1, 0.1, 0.12, 0.16, 0.08])
    r.bullets([
        "Every printed value of the first row is reproduced; eq. 49 equals the cold-plasma dispersion to 4×10<super>−17</super>.",
        "His step from k = 2642 to 2686 /m skips the nearest zero of j<sub>1</sub> (265.46) for the next one.",
        "Scaled to 50 µm (his eq. 86, energy proportional to R): 500 J and λ = 1.17 µm as printed, but 4.27×10<super>14</super> electrons "
        "against the printed 4.3×10<super>13</super>, and 30 rather than \"about forty\" times sea-level density.",
    ])

    r.page_break()
    r.h2("2. Harmonics and size")
    r.figure(os.path.join(HERE, "harmonics_and_size.png"),
             "Figure 2. Left: the harmonics eq. 35 allows a 12 cm plasmoid, above the critical frequency αc/2π, with the free "
             "TE lifetime of each. Middle: free radiative lifetime against radius, by the paper's recipe (ξ = 0.01; electrons "
             "above c) and with the electrons kept at or below 0.1 c, for 10 kJ and 1 MJ, against the reported 1–10 s (gold). "
             "Right: the free-electron density a free plasmoid needs to last τ, against fully stripped sea-level air.",
             max_height_mm=90)
    r.table([
        ["Radius", "Lifetime wanted", "Frequency", "Free electrons ÷ sea-level air molecules"],
        ["5 cm", "1 – 10 s", "470 – 1485 THz", "102 – 1020"],
        ["12 cm (mean)", "1 – 10 s", "126 – 399 THz", "7.4 – 74"],
        ["25 cm", "1 – 10 s", "42 – 133 THz", "0.8 – 8.1"],
        ["50 cm", "1 – 10 s", "15 – 47 THz", "0.10 – 1.0"],
    ], [0.18, 0.2, 0.22, 0.4])
    r.p("For kR much larger than 1 the free energy lifetimes are τ<sub>TE</sub> = R(kR)²/(2c·x<sub>1</sub>²) and "
        "τ<sub>TM</sub> = R(kR)²|g′(x<sub>t</sub>)|/(2c·x<sub>t</sub>³), with g = 1 + x·j<sub>1</sub>′/j<sub>1</sub>, "
        "x<sub>1</sub> = 4.4934, x<sub>t</sub> = 2.7437; both agree with the complex root finder to four figures at kR = 268.6, "
        "1000 and 3000. So lasting τ fixes kR, and eq. 49 then fixes the electron density, whatever the energy.")

    r.h2("3. kronspice analysis")
    r.figure(os.path.join(HERE, "kronspice_spectrum.png"),
             "Figure 3. Left: the Appendix B plasmoid's four TE and four TM l = 1 resonances as a Foster chain of parallel "
             "RLC tanks in a kronspice AC sweep (TE1 is Nachamkin's 128.158 GHz mode). Right: a kronspice ring-down of one "
             "tank, its envelope drawn from the ledger's stored energy.", max_height_mm=80)
    r.bullets([
        "Spectrum: each peak's frequency matches the mode finder to ≤ 2×10<super>−15</super> and its half-power Q to 2×10<super>−6</super>.",
        "Ring-down of TE1 over 400 periods: Q = 4.7992×10<super>5</super> against 4.7985×10<super>5</super>; heat = energy lost to "
        "5×10<super>−6</super> of the stored energy; energy lifetime 0.60 µs.",
        "Eq. 49 as a circuit: a one-dimensional plasma line (series μ<sub>0</sub>Δ, shunt ε<sub>0</sub>Δ, electrons' kinetic "
        "inductance m/(ne²Δ) at every node), shorted half a wavelength π/α apart, resonates at ω² = ω<sub>p</sub>² + α²c²: "
        "kronspice equals the exact lattice to 7×10<super>−10</super>, and the lattice approaches eq. 49 as 1/N² "
        "(7.3×10<super>−5</super> at 50 cells, 4.6×10<super>−6</super> at 200).",
    ])

    r.h1("The reported seconds")
    r.p("Witnesses report ball lightning lasting seconds, and this example does not dispute the reports; it measures "
        "the paper's model against them. As written, a free plasmoid of the paper's kind radiates its field away in "
        "microseconds (the paper's recipe) to tens of milliseconds (electrons kept non-relativistic), at every reported "
        "size and energy. Lasting seconds would need more free electrons than the air contains at the mean size. So "
        "either the plasmoid is not free, or its energy is not mainly in this field:")
    r.bullets([
        "<b>Held from outside.</b> The paper itself (§IX) says the boundary must be matched by an external standing field "
        "or by surface currents; a plasmoid fed or held that way can last as long as its supply. The microwave-sustained "
        "cavity plasmoids the paper mentions are the laboratory case. The test: a measured external field or surface "
        "current at the resonant size.",
        "<b>Energy stored elsewhere.</b> Standard models that reach seconds store the energy chemically, not in the field: "
        "luminous balls from electrical discharges into silicon last several seconds in the laboratory (Abrahamson &amp; "
        "Dinniss, <i>Nature</i> 403, 519, 2000; Paiva et al., <i>Phys. Rev. Lett.</i> 98, 048501, 2007). These are context, "
        "not part of the paper's model.",
    ])

    r.h1("Ledger and validation")
    r.bullets([
        "kronspice ring-down: heat equals the energy lost to 5×10<super>−6</super> of the stored energy.",
        "Foster-chain AC: frequency to 10<super>−15</super>, Q to 2×10<super>−6</super> for all eight modes.",
        "Large-kR lifetimes equal the complex root finder (TE to 10<super>−3</super>, TM to 5×10<super>−3</super> in the tests).",
        "<i>main_test.go</i>: field and B3, Appendix B, harmonics, lifetimes, kronspice spectrum, ring-down and plasma line.",
    ])

    r.h1("Limits of the model")
    r.bullets([
        "Lifetimes are for a uniform cold-plasma sphere radiating into free space; the vortical flow, collisions, "
        "recombination and any external field are not modelled.",
        "The energy follows the paper's split (half field, half fluid) and its time-averaged field energy (B2).",
        "Only l = 1 is analysed in detail, as in the paper; the minimum radii of l = 2–4 are listed.",
        "Ball-lightning statistics are eyewitness-based and broad.",
    ])

    r.h1("Status")
    r.p("Register entry H26 stays <i>untested</i> as a ball-lightning model. Its field mathematics and charge condition "
        "are sound; its worked example has a √α slip and drives electrons faster than light; and a free plasmoid of its "
        "kind cannot last the seconds reported at any reported size without more free electrons than air has. What would "
        "move it: a measured external field or surface current holding a plasmoid at its resonant size.")

    n0 = len(r.s)
    r.h1("Reproducing the results")
    r.code("go run ./examples/Force-Free_Harmonic_Plasmoids -viz examples/Force-Free_Harmonic_Plasmoids/viz\n"
           "python examples/Force-Free_Harmonic_Plasmoids/viz/render.py      # pictures and plasmoid_cycle.gif\n"
           "python examples/Force-Free_Harmonic_Plasmoids/viz/make_report.py # this report\n"
           "go test ./examples/Force-Free_Harmonic_Plasmoids")
    r.small("Sources: J. Nachamkin, <i>Force-Free Time-Harmonic Plasmoids</i>, PL-TR-92-3044 (1992), DTIC AD-A257 765; "
            "S. Chandrasekhar, P. C. Kendall, <i>Astrophys. J.</i> 126, 457 (1957); J. A. Stratton, <i>Electromagnetic "
            "Theory</i> (1941), §9.22; F. F. Chen, <i>Introduction to Plasma Physics</i>; ball-lightning statistics: "
            "arXiv:physics/0312122, arXiv:1102.2829; J. Abrahamson, J. Dinniss, <i>Nature</i> 403, 519 (2000); "
            "G. S. Paiva et al., <i>Phys. Rev. Lett.</i> 98, 048501 (2007).")
    r.s[n0:] = [KeepTogether(r.s[n0:])]
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "Force-Free_Harmonic_Plasmoids_report.pdf"))
