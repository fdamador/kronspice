# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Knotted light (H29 step 1, H30) as a PDF report.

Run from the repository root:
    python examples/knotted_light/viz/make_report.py examples/knotted_light/knotted_light_report.pdf

Numbers come from examples/knotted_light/run_output.txt (the program's own
output) and README.md; figures are viz/knot_3d.png, viz/knot_time.png and
viz/h30_cycle.png. Layout from internal/pdfreport/kronreport.py.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401


def build(out):
    r = Report("examples/knotted_light", "Knotted light",
               "Does topology alone hold a knot of light together, and what does a pulsing knot cost? "
               "Register entries H29 (step 1) and H30",
               "register entries H29, H30",
               "10 cells per L0, region ±6 L0, absorber 12 cells")

    r.h1("Summary")
    r.p("Rañada's knotted light (the hopfion) is an exact solution of Maxwell's equations in vacuum whose "
        "magnetic field lines are closed loops, every pair linked once. It was set in the kronspice FDTD field "
        "solver (<i>em/fdtd.go</i>) and run against the exact solution, to ask whether topology alone holds a "
        "knot of field energy together (H29, step 1). The incoming, contracting half of the same knot was then fed "
        "in through a Huygens surface, so the ledger books every joule that enters and leaves (H30), and the "
        "measured QED (Euler–Heisenberg) vacuum was tested as a possible source of that inward wave.")
    r.bullets([
        "<b>The topology is kept.</b> The exact field's traced B lines link once (+1.0000) at every time to "
        "t = 4 L0/c; the solver's magnetic helicity follows the exact value (0.9873 against 0.9824 at t = 4).",
        "<b>The energy is not held.</b> Energy within 2 L0 falls from 0.9659 (solver; exact 0.9662) at t = 0 to "
        "0.0186 (exact 0.0200) by t = 3; the rms radius grows as √(1 + t²) (0.996 → 4.117 L0).",
        "<b>Every joule of the implosion and explosion is booked.</b> The surface puts in 0.9862 of a knot and takes "
        "out 0.9799; the difference, 6.6×10<super>−3</super>, leaks numerically into the absorber. Ledger residual "
        "−1.8×10<super>−13</super> of one knot's energy.",
        "<b>The knot survives the implosion.</b> At t = 0 two lines traced on the solver's own field link +1.0000 "
        "(gap 0.034 L0).",
        "<b>A pulsing knot costs one knot's energy per period from outside.</b> With period 8 L0/c it is rebuilt to "
        "0.9571, 0.9579, 0.9573 (exact 0.9668, 0.9674, 0.9668), linking +1.0000 each time; the net work over the "
        "middle cycle is +0.0046 knot energies, equal to the leak (in = out).",
        "<b>The measured QED vacuum returns nothing.</b> Of an explosion, 6.1×10<super>−28</super> comes back between "
        "the separated halves (round-off); at most about 6a ≈ 1.25×10<super>−6</super> of the field energy could "
        "ever be redirected, against the 1 knot per period needed.",
    ])

    r.h1("What is being tested")
    r.h2("H29: particles as knots of field energy")
    r.p("The user hypothesis (2026-11-05) is that particles are knots of field energy, held together by their "
        "topology, since a knot cannot untie without being cut. Kelvin proposed the same for atoms in 1867 (vortex "
        "knots in the ether). Step 1 asks the question every ring model in the repository has met: under ordinary, "
        "linear Maxwell electrodynamics, does a knotted field stay together?")
    r.p("<b>Standard prediction.</b> Linking is conserved, because Maxwell's equations carry the topology, but the "
        "energy is not held: the knot expands at c while staying linked. Derrick's theorem (1964) says no static, "
        "finite lump of a linear field exists in three dimensions; holding a knot needs a nonlinear term (Skyrme 1961; "
        "Faddeev and Niemi 1997).")
    r.h2("H30: implosion and explosion")
    r.p("The user hypothesis (2026-11-05, after Walter Russell's \"rhythmic balanced interchange\") is that a "
        "particle is a cycle of implosion and explosion, and that its energy books balance only when the inward "
        "half is counted with the outward one.")
    r.p("<b>Standard prediction.</b> Maxwell's equations are time-symmetric, so the contraction is the expansion run "
        "backwards. A knot rebuilt every period needs an incoming wave supplied from outside at one knot's energy per "
        "period, here a power U·c/(8 L0), where U is the knot's energy.")

    r.h1("The model")
    r.p("The field is the hopfion in the Bateman form of Kedia et al. (2013), with lengths in L0 and times in L0/c:")
    r.p("&nbsp;&nbsp;&nbsp;&nbsp;α = (r² − t² − 1 + 2iz) / (r² − (t − i)²)<br/>"
        "&nbsp;&nbsp;&nbsp;&nbsp;β = 2(x − iy) / (r² − (t − i)²)<br/>"
        "&nbsp;&nbsp;&nbsp;&nbsp;F = E + icB = (grad α) × (grad β) = curl(α grad β)")
    r.bullets([
        "It is set in <i>em/fdtd.go</i> as discrete curls of Re(α grad β) and Im(α grad β), so div E = div B = 0 exactly on the "
        "grid. The helicity sign of the formula is checked numerically against Maxwell's equations (+1).",
        "The solver also carries the vector potential A in the temporal gauge (A ← A − Δt·E after each step), so "
        "curl A = B at every step and the magnetic helicity ∫A·B dV, the integral measure of field-line linking, can be "
        "followed. For a null field such as the hopfion (E·B = 0, |E| = c|B|) it is conserved exactly.",
        "Linking of two field lines is the Gauss linking integral. Field-line seeds at (±0.6, 0, 0) L0.",
        "H30: a closed Huygens surface (total-field/scattered-field; Taflove and Hagness ch. 5), a cube of half-width "
        "5.75 L0, made of 160,080 edge currents and 158,700 magnetic face currents (<i>em.FaceSource</i>, with its own "
        "exact ledger term). The declared incoming wave stands in for whatever medium would have to supply it.",
        "One knot = 1.748×10<super>−10</super> J at E0 = 1 V/m, L0 = 1 m.",
    ])

    r.h1("Results")
    r.h2("H29 step 1: the knot in the solver")
    r.table([
        ["t (L0/c)", "within 2 L0, solver / exact", "absorbed", "rms radius (L0)", "helicity, solver / exact",
         "curl A − B", "linking exact / solver (closed)"],
        ["0.00", "0.9659 / 0.9662", "0.0000", "0.996", "1.0000 / 1.0000", "0", "+1.0000 / +1.0000 (yes)"],
        ["0.50", "0.9554 / 0.9564", "0.0001", "1.114", "0.9999 / 0.9999", "3.5×10<super>−15</super>", "+1.0000 / +1.0043 (yes)"],
        ["1.00", "0.9050 / 0.9102", "0.0001", "1.411", "0.9992 / 0.9995", "6.0×10<super>−15</super>", "+1.0000 / +1.0000 (yes)"],
        ["1.50", "0.7297 / 0.7513", "0.0001", "1.800", "0.9982 / 0.9988", "9.7×10<super>−15</super>", "+1.0000 / not closed"],
        ["2.00", "0.3990 / 0.3918", "0.0001", "2.233", "0.9971 / 0.9976", "1.3×10<super>−14</super>", "+1.0000 / not closed"],
        ["3.00", "0.0186 / 0.0200", "0.0003", "3.159", "0.9935 / 0.9933", "2.0×10<super>−14</super>", "+1.0000 / not closed"],
        ["4.00", "0.0012 / 0.0012", "0.0007", "4.117", "0.9873 / 0.9824", "2.9×10<super>−14</super>", "+1.0000 / not closed"],
    ], [0.08, 0.17, 0.1, 0.1, 0.17, 0.13, 0.25])
    r.small("Energies are fractions of the starting energy. Solver helicity at t = 0 is 0.9921 of the exact value "
            "(grid discretization at 10 cells per L0); the helicity columns are normalized to their own t = 0 values. "
            "Lines traced on the solver's field stop closing after t = 1 (grid errors turn near-closed loops into open "
            "ones), which is why helicity, an integral, is the solver's topology measure.")
    r.figure(os.path.join(HERE, "knot_3d.png"),
             "Figure 1. 24 magnetic field lines of the exact hopfion at t = 0, seeded on three nested rings "
             "(light: near the axis; dark: near the core). All close, and the first pair links once (1.0000).",
             width_frac=0.42)
    r.figure(os.path.join(HERE, "knot_time.png"),
             "Figure 2. The same family of exact field lines at t = 0, 1, 2 and 3 L0/c on one scale, seeds scaled by "
             "the knot's size √(1 + t²). The grey sphere marks the energy's rms radius (1.00, 1.41, 2.23, 3.16 L0). "
             "All 24 lines close at every time and the linking stays 1.000, while the knot spreads at the speed of light.",
             width_frac=0.78)
    r.bullets([
        "The slow fall of the helicity, in solver and exact solution alike, is the field's tail leaving the finite "
        "box, not the knot coming undone.",
        "Conclusion of step 1: linear Maxwell carries the knot's topology perfectly and does nothing to hold its "
        "energy. This is the standard prediction (Rañada 1989; Irvine and Bouwmeester 2008).",
    ])

    r.h2("H30 Part 1: one implosion, then one explosion")
    r.table([
        ["t (L0/c)", "within 2 L0, solver / exact", "inside surface", "put in / taken out", "leak (numerical)",
         "field error", "linking (gap, L0)"],
        ["−10", "0.0000 / 0.0000", "0.0041", "0.0000 / 0.0000", "0", "0.000", "—"],
        ["−6", "0.0000 / 0.0000", "0.8198", "0.8157 / 0.0000", "7.2×10<super>−8</super>", "0.001", "—"],
        ["−4", "0.0012 / 0.0012", "0.9883", "0.9842 / 0.0000", "8.5×10<super>−8</super>", "0.003", "—"],
        ["−2", "0.3773 / 0.3934", "0.9902", "0.9861 / 0.0000", "8.5×10<super>−8</super>", "0.030", "—"],
        ["0", "0.9567 / 0.9661", "0.9903", "0.9862 / 0.0000", "8.5×10<super>−8</super>", "0.043", "+1.0000 (0.034)"],
        ["2", "0.4129 / 0.3893", "0.9902", "0.9862 / 0.0001", "8.8×10<super>−8</super>", "0.060", "—"],
        ["4", "0.0012 / 0.0012", "0.9883", "0.9862 / 0.0020", "4.5×10<super>−7</super>", "0.004", "—"],
        ["6", "0.0000 / 0.0000", "0.8390", "0.9862 / 0.1498", "1.6×10<super>−3</super>", "0.001", "—"],
        ["10", "0.0000 / 0.0000", "0.0038", "0.9862 / 0.9799", "6.6×10<super>−3</super>", "0.080", "—"],
    ], [0.09, 0.2, 0.12, 0.18, 0.14, 0.11, 0.16])
    r.small("Energies as fractions of one knot's energy; selected rows of run_output.txt. Ledger residual "
            "−1.8×10<super>−13</super> of one knot's energy.")

    r.h2("H30 Part 2: a pulsing knot")
    r.p("Knots tightest at t = 0, 8 and 16 L0/c are fed through the same surface; the first starts inside the "
        "surface at t = −4, the next two come in through it.")
    r.table([
        ["peak at t (L0/c)", "within 2 L0, solver", "exact", "linking of two solver lines (gap)"],
        ["0", "0.9571", "0.9668", "+1.0000 (0.020 L0)"],
        ["8", "0.9579", "0.9674", "+1.0000 (0.035 L0)"],
        ["16", "0.9573", "0.9668", "+1.0000 (0.037 L0)"],
    ], [0.22, 0.22, 0.2, 0.36])
    r.figure(os.path.join(HERE, "h30_cycle.png"),
             "Figure 3. Top: Part 1, one knot fed in from t = −10 and out by t = 10. Bottom: Part 2, the pulsing "
             "knot. Blue: energy within 2 L0 (thick light line exact, dots solver). Dashed green and orange: the "
             "surface's cumulative work in and out. When an arriving and a departing knot overlap the columns show "
             "only the net flow, so the gross flow each way cannot be read from them.",
             width_frac=0.9)
    r.bullets([
        "Over the middle cycle (t = 4 … 12) the surface's net work is +0.0046 knot energies, equal to the energy that "
        "leaked past it: in equals out. Ledger residual −3.7×10<super>−13</super> of one knot's energy.",
        "Under linear Maxwell an imploding and exploding knot keeps its topology through both halves and can be kept "
        "pulsing indefinitely, at the cost of an inward wave of one knot's energy per period supplied from outside. "
        "The run books that inflow but supplies no source for it.",
    ])

    r.h2("H30 Part 3: the measured QED vacuum as the source of the inward wave")
    r.p("Euler–Heisenberg's vacuum is the one nonlinear vacuum with experimental support (light-by-light scattering, "
        "ATLAS 2017, CMS 2019). Waves act on one another only where they cross; a lone plane wave never acts on "
        "itself. Both measurements run on Kron's plane-wave network with the Euler–Heisenberg vacuum as its nonlinear "
        "primitive (<i>kron.EHLine</i>, validated on light-by-light phase shifts to 0.4 %), beside Maxwell's linear "
        "vacuum, at 0.08–0.09 of the Schwinger field E<sub>s</sub> (1.3×10<super>18</super> V/m). Here "
        "a = (α/90π)(E/E<sub>s</sub>)² = 2.09×10<super>−7</super> at 0.09 E<sub>s</sub>.")
    r.table([
        ["measurement", "Euler–Heisenberg vacuum", "needed for a knot rebuilt every period"],
        ["A. Explosion into two halves: nonlinear energy back between the separated halves",
         "6.1×10<super>−28</super> of the explosion (linear lattice alone: 3.4×10<super>−17</super>)", "1"],
        ["Nonlinear part of the field after the halves crossed (phase shifts, all outgoing)",
         "3.1×10<super>−13</super>", "—"],
        ["B. Probe crossing a counter-propagating field of 0.08 E<sub>s</sub> (a = 1.65×10<super>−7</super>): energy turned round",
         "2.4×10<super>−17</super> (Fresnel bound (8a)² = 1.7×10<super>−12</super>)", "1"],
        ["Nonlinear share of the field energy at the peak (the most the vacuum could redirect)",
         "1.25×10<super>−6</super> (≈ 6a)", "1"],
    ], [0.44, 0.36, 0.2])
    r.bullets([
        "Nothing comes back: what the vacuum gives the crossing halves is a phase shift that travels outward with "
        "them; once apart, each is a plane wave the vacuum no longer touches.",
        "In three dimensions the expanding shell is locally a plane wave within a few L0, so the bound on what the "
        "vacuum could ever turn round is about 10<super>−6</super> of the knot per period, against the whole knot needed.",
        "So this route is bounded. The Wheeler–Feynman absorber, the zero-point field and Born–Infeld remain untested.",
    ])

    r.h1("Ledger and validation")
    r.table([
        ["check", "result (run_output.txt)"],
        ["H29 energy ledger residual", "4.6×10<super>−14</super> of the starting energy"],
        ["max |div B| (relative)", "5.2×10<super>−15</super>"],
        ["curl A − B", "≤ 2.9×10<super>−14</super> of max B (t = 4)"],
        ["H30 Part 1 ledger residual", "−1.8×10<super>−13</super> of one knot's energy"],
        ["H30 Part 2 ledger residual", "−3.7×10<super>−13</super> of one knot's energy"],
    ], [0.45, 0.55])
    r.table([
        ["test", "checks", "held to"],
        ["TestLinkingKnownCurves", "Gauss linking integral: Hopf link = ±1, separated circles = 0", "1e-3"],
        ["TestHopfionExact", "field is null (F·F = 0); traced B lines close and link once at t = 0 and t = 3", "1e-12; 1e-3"],
        ["TestKnotInSolver", "coarse grid: energy spreads as exact; helicity kept; curl A = B; ledger", "0.02; 0.02; 1e-12; 1e-10"],
        ["TestKnotEnergy", "energy density at t = 0 spherically symmetric; radial total = 3-D sum", "1e-12; 2e-3"],
        ["TestImplosionThroughSurface", "≥ 0.9 of a knot enters; centre reaches exact energy; linking 1 at t = 0; in = out; ledger",
         "0.05; 0.02; 0.05; 1e-10"],
        ["TestVacuumInward", "EH halves return nothing (below linear lattice ripple); reflection below Fresnel bound; nonlinear share 6a",
         "1e-24; Fresnel; 1 %"],
        ["em.TestFDTDMagneticSourceLedger", "ledger closes with magnetic face sources; sources in the absorber refused", "1e-10"],
    ], [0.3, 0.5, 0.2])

    r.h1("Limits of the model")
    r.bullets([
        "The grid has 10 cells per L0; the solver's starting helicity is 0.9921 of the exact value, and after the wave "
        "has crossed 10 L0 of grid the solver's field is 4.3 % from the exact one within 2 L0.",
        "Field lines traced on the solver's field are fragile: they close only up to t = 1 in step 1, and in H30 each "
        "is closed at its nearest return (gaps 0.020–0.037 L0).",
        "The declared incoming wave in H30 stands in for a medium; the run gives no physical source for it.",
        "Part 3 is computed on a one-dimensional plane-wave network at the top of the Euler–Heisenberg weak-field "
        "range; the three-dimensional statement is an argument from the local plane-wave form of the expanding shell.",
    ])

    r.h1("Status")
    r.p("<b>H29 step 1: bounded.</b> Linear Maxwell keeps the topology and not the energy, which is the standard "
        "prediction; a holding knot needs a nonlinear term (step 2 is <i>examples/ether_knot</i>). "
        "<b>H30: bounded.</b> Linear Maxwell gives both halves by time symmetry; the pulsing knot needs an inward wave "
        "of one knot per period from outside; the measured QED vacuum is bounded as its source at about 10<super>−6</super> or less "
        "of a knot per period; the absorber, zero-point field and Born–Infeld routes are untested.")

    r.h1("Reproducing the results")
    r.code("go run ./examples/knotted_light               # step 1 and H30, about 3 min\n"
           "go run ./examples/knotted_light -h29only      # step 1 only\n"
           "go test ./examples/knotted_light              # -short skips the 3-D solve\n"
           "go run ./examples/knotted_light -viz examples/knotted_light/viz   # traced 3-D lines\n"
           "python3 examples/knotted_light/viz/render.py       # knot_3d.png, knot_time.png\n"
           "python3 examples/knotted_light/viz/render_h30.py   # h30_cycle.png from run_output.txt")
    r.small("Sources: Rañada, Lett. Math. Phys. 18, 97 (1989); Irvine and Bouwmeester, Nature Physics 4, 716 (2008); "
            "Kedia et al., PRL 111, 150404 (2013); Taflove and Hagness, Computational Electrodynamics, 3rd ed., ch. 5; "
            "Wheeler and Feynman, Rev. Mod. Phys. 17, 157 (1945); Born and Infeld (1934); Heisenberg and Euler (1936); "
            "Moffatt, J. Fluid Mech. 35, 117 (1969); Derrick (1964); Skyrme (1961); Faddeev and Niemi, Nature 387, 58 "
            "(1997); Kelvin, Phil. Mag. 34, 15 (1867); Russell (1926, 1947).")
    for f in r.s:
        if getattr(getattr(f, 'style', None), 'name', '') in ('h1', 'h2'):
            f.keepWithNext = 1
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "knotted_light_report.pdf"))
