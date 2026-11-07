# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""H29 step 2 report: a Faddeev knot held by a nonlinear ether, as a PDF.

python examples/ether_knot/viz/make_report.py examples/ether_knot/ether_knot_report.pdf

Numbers come from examples/ether_knot/run_output.txt (relaxation, controls,
3-D check), examples/ether_knot/vib_output.txt (vibrations) and README.md
(eigenmode box study, tests). The example has no figures. Layout from
internal/pdfreport/kronreport.py.
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
    r = Report("examples/ether_knot", "Ether knot",
               "Can a nonlinear ether hold a knot of field energy at a fixed size? — register entry H29, step 2",
               "register entry H29",
               "κ2 = κ4 = 1; arrested Newton flow, dt = 0.05 h")

    r.h1("Summary")
    r.p("Hypothesis H29 holds that particles are knots of field energy, held together by their topology. Step 1 "
        "(<i>knotted_light</i>) found that linear Maxwell keeps a knot's topology but not its energy: the knot links "
        "once for ever while flying apart at the speed of light. Step 2, here, asks whether a nonlinear, fluid-like "
        "ether can hold a knot by itself, at a size of its own. The ether chosen is a unit-vector field n(x) whose "
        "energy is its stiffness plus the magnetic energy of knotted light (Faddeev's model).")
    r.bullets([
        "<b>The knot settles at a fixed size:</b> extrapolated to zero grid spacing, energy E = 547.5·√(κ2κ4) and "
        "stiffness/magnetic ratio E2/E4 = 0.9972 (the virial theorem requires 1), core radius about 1.2 ℓ.",
        "<b>Topology intact in 3-D:</b> two pairs of preimage curves link at 1.0001 and 1.0002 — Hopf charge 1.",
        "<b>Controls behave as Derrick's theorem predicts:</b> the magnetic term alone spreads the knot (rms radius "
        "1.003 → 5.125 ℓ, energy 311.713 → 0.087); the stiffness alone shrinks it until it unties through the grid.",
        "<b>The knot has a tone:</b> struck, it rings by squashing at ω ≈ 1.36–1.42 c/ℓ with Q ≈ 6–7; breathing is "
        "broad (ω ≈ 0.43–0.50 c/ℓ, Q ≈ 2) and twisting hardly rings (Q ≈ 1–2).",
        "<b>Gap to the literature recorded:</b> Battye &amp; Sutcliffe's 100³ lattice gave E = 504 (spacing not "
        "stated), 8 % below this run.",
        "<b>Not fixed by the model:</b> κ2 and κ4; no measured ether stiffness sets them, so size, energy and "
        "frequency are ratios, not yet numbers.",
    ])

    r.h1("What is being tested")
    r.p("<b>Hypothesis (H29).</b> Particles are knots of field energy, held together by their topology (a knot cannot "
        "untie without being cut); Kelvin's vortex atoms (<i>Phil. Mag.</i> 34, 15, 1867) in field form. Step 2 asks "
        "whether a nonlinear ether holds a knot at rest without any inflow of energy.")
    r.p("<b>Standard prediction.</b> Derrick (1964): no static lump of a linear field exists in three dimensions. "
        "Under a change of size x → λx the stiffness term grows as λ and the magnetic term falls as 1/λ, so alone the "
        "magnetic term spreads the knot and the stiffness shrinks it; together they fix a size where the two are "
        "equal. A holding knot needs a nonlinear term (Skyrme 1961; Faddeev &amp; Niemi, <i>Nature</i> 387, 58, 1997).")
    r.h2("Which ether? (literature survey from the README)")
    r.table([
        ["candidate", "what the literature says", "holds a knot at rest?"],
        ["Light made nonlinear: Born–Infeld (1934), Euler–Heisenberg (1936, seen as light-by-light scattering by "
         "ATLAS, Nature Phys. 13, 852, 2017)",
         "the only finite-energy, source-free, stationary Maxwell–Born–Infeld solution is zero (Kiessling, "
         "arXiv:1107.2333, Theorem 1)", "no"],
        ["an ordinary superfluid (Gross–Pitaevskii), Kelvin's vortex atoms",
         "all 1,458 simulated vortex knots untie by reconnection at the vortex cores, where the field vanishes "
         "(Kleckner, Kauffman &amp; Irvine, Nature Phys. 12, 650, 2016)", "no"],
        ["a fluid whose order parameter never vanishes: a unit vector n(x); a charged two-component superfluid "
         "reduces to it (Babaev, Faddeev &amp; Niemi, PRB 65, 100512, 2002)",
         "Faddeev's knots; stable, computed by Battye &amp; Sutcliffe (PRL 81, 4798, 1998): Q = 1 energy 504 on "
         "their 100³ lattice", "yes: tested here"],
    ], [0.38, 0.47, 0.15])

    r.h1("The model")
    r.p("The ether's energy is")
    r.code("E = κ2 ∫ (∂_i n)² d³x  +  κ4 ∫ |B_n|² d³x,      B_n,k = ½ ε_kij n·(∂_i n × ∂_j n)")
    r.p("For Rañada's Hopf map, the starting field, B<sub>n</sub> is exactly Rañada's knotted magnetic field times −16 "
        "at every point (tested). So the second term is the knotted light's magnetic energy, the term step 1 already "
        "had, and the first term (E2) is the new ingredient: the ether's stiffness. Units: κ2 = κ4 = 1, so lengths "
        "are in ℓ = √(κ4/κ2) and energies in √(κ2κ4).")
    r.p("The Q = 1 knot is axially symmetric, with the field turning once round the axis: n = R<sub>z</sub>(φ) m(ρ, z). "
        "The energy is therefore minimized on the (ρ, z) half-plane, about a thousand times cheaper than in 3-D, with m "
        "kept a unit vector. Relaxation is arrested Newton flow (as Battye &amp; Sutcliffe): the field accelerates "
        "downhill and is stopped whenever the energy rises. The start is Rañada's map pulled in smoothly to reach the "
        "vacuum inside the box. The relaxed knot is then rebuilt on a 3-D grid to check its energy and topology.")

    r.h1("Results")
    r.h2("1. Stiffness + magnetic: grid and box study")
    r.table([
        ["box (ℓ)", "h (ℓ)", "steps", "E", "E2 (stiffness)", "E4 (magnetic)", "E2/E4", "core radius (ℓ)",
         "rms radius (ℓ)"],
        ["10", "0.1", "4000", "541.084", "261.147", "279.937", "0.9329", "1.150", "1.816"],
        ["16", "0.1", "6000", "540.244", "262.109", "278.135", "0.9424", "1.150", "1.870"],
        ["16", "0.07", "8000", "544.004", "268.062", "275.942", "0.9714", "1.225", "1.911"],
        ["16", "0.05", "12000", "545.701", "270.656", "275.046", "0.9840", "1.225", "1.929"],
        ["h → 0", "(h = 0.07, 0.05)", "", "<b>547.5</b>", "", "", "<b>0.9972</b>", "", ""],
    ], [0.08, 0.13, 0.08, 0.1, 0.13, 0.13, 0.09, 0.13, 0.13], bold_rows=(5,))
    r.bullets([
        "The box edge hardly matters (about 1 part in 500 from 10 ℓ to 16 ℓ); the grid does. As h → 0 the stiffness "
        "and magnetic terms become equal to 0.3 % — Derrick's balance, reached by the relaxation without being "
        "imposed.",
        "Battye &amp; Sutcliffe's 100³ lattice gave 504 (their spacing is not stated). Coarse grids give lower "
        "energies here too (540 at h = 0.1), so the 8 % gap is recorded, not explained away; a continuum value from "
        "the literature in the same normalization is still to be checked.",
    ])
    r.h2("2. Controls (box 10 ℓ, h 0.1 ℓ, from the same start)")
    r.table([
        ["magnetic only (κ2 = 0)", "step 0", "1500", "3000", "4500", "6000"],
        ["E", "311.713", "6.835", "0.600", "0.148", "0.087"],
        ["rms radius (ℓ)", "1.003", "3.229", "4.182", "5.124", "5.125"],
    ], [0.25, 0.15, 0.15, 0.15, 0.15, 0.15])
    r.table([
        ["stiffness only (κ4 = 0)", "step 0", "100", "200", "300", "500", "700", "1000"],
        ["E", "270.795", "236.362", "129.480", "66.072", "37.233", "14.569", "4.943"],
        ["core radius (ℓ)", "0.950", "0.750", "0.350", "2.650", "2.550", "5.650", "1.750"],
        ["m<sub>z</sub> at core", "−0.997", "−0.994", "−0.016", "+0.754", "+0.849", "+0.979", "+0.986"],
    ], [0.25, 0.11, 0.11, 0.11, 0.11, 0.11, 0.1, 0.1])
    r.p("Without the ether (κ2 = 0, knotted light alone) the knot spreads: step 1's result, now as an energy "
        "minimization. With stiffness alone it shrinks; by step 200 the core reaches the grid spacing, m<sub>z</sub> "
        "at the core leaves −1 and the knot unties through the lattice. In the continuum it would shrink to a point; "
        "a grid cannot hold topology below its own spacing.")
    r.h2("3. The relaxed knot in 3-D")
    r.p("The finest knot (box 16 ℓ, h = 0.05) was rebuilt on 101³ points over ±8 ℓ. Two pairs of preimage curves of n "
        "(the field lines of B<sub>n</sub>), each closed at its nearest return, link at <b>1.0001</b> (gap 0.019 ℓ) and "
        "<b>1.0002</b> (gap 0.011 ℓ): Hopf charge 1, the knot intact. The 3-D energy reads E2 297.28 + E4 255.40 = "
        "552.68 against 270.656 + 275.046 in the axisymmetric solve: the ±8 ℓ box cuts the tail of a knot relaxed in a "
        "16 ℓ box (stiffness high) and the 3-D grid (h = 0.16) is about three times coarser (magnetic term low by "
        "O(h²)). A knot relaxed in an 8 ℓ box and rebuilt the same way matched to 0.5 % in stiffness (258.7 against "
        "260.0), with the magnetic term 6 % low (README).")

    r.h2("4. Vibrations: the knot's own tone")
    r.p("Small motions x about the relaxed knot obey 2K ẍ = −H x, from Faddeev's Lagrangian (kinetic energy "
        "∫|ṅ|² + Σ|ṅ × ∂<sub>i</sub>n|²): H is the energy's Hessian on the sphere and K the kinetic form, both built "
        "for motions that keep the axial symmetry. Two motions cost nothing and come out at ω² ≈ 0 (sliding along z, "
        "turning the field about ẑ): 2.2×10<super>−4</super> and 9.9×10<super>−8</super> (box 16 ℓ, h 0.15), "
        "1.6×10<super>−4</super> and −3.5×10<super>−7</super> (box 24 ℓ), 1.5×10<super>−4</super> and "
        "1.6×10<super>−7</super> (box 16 ℓ, h 0.1), against the real modes' ~1.")
    r.p("Small waves of this ether travel at c with no lowest frequency, so a vibrating knot leaks energy and has no "
        "sealed tone. In a box the low eigenmodes are standing waves of the box: the lowest falls 0.428 → 0.324 → "
        "0.252 c/ℓ as the box grows 8 → 12 → 16 ℓ, and its share within 3 ℓ of the knot falls 42 % → 9 % → 2 % "
        "(README; vib_output.txt gives 0.2519 c/ℓ and 0.02 for box 16 ℓ). So the knot was struck instead, with three "
        "localized pushes (breathing, squashing, twisting, under exp(−r²/8ℓ²), zero modes removed), left to ring, and "
        "its motion fitted to damped cosines (matrix pencil) up to the return of the first echo from the wall. "
        "Strongest component of each push:")
    r.table([
        ["push", "box 16 ℓ, h 0.15 ℓ", "box 24 ℓ, h 0.15 ℓ", "box 16 ℓ, h 0.10 ℓ"],
        ["<b>squashing</b>", "<b>ω 1.416, Q 7.0</b>", "<b>ω 1.415, Q 7.0</b>", "<b>ω 1.385, Q 6.0</b>"],
        ["breathing", "ω 0.429, Q 1.9", "ω 0.447, Q 1.9", "ω 0.498, Q 2.0"],
        ["twisting", "ω 0.931, Q 1.7", "ω 0.859, Q 1.1", "ω 0.829, Q 1.3"],
    ], [0.19, 0.27, 0.27, 0.27])
    r.small("ω in c/ℓ. Records to t = 20 ℓ/c (box 16 ℓ) and 36 ℓ/c (box 24 ℓ). Weaker components are listed in "
            "vib_output.txt.")
    r.bullets([
        "<b>Squashing is the knot's tone.</b> A bigger box leaves it unchanged (1.416 → 1.415); the finer grid lowers "
        "it 2 %; extrapolated in h², about 1.36 c/ℓ. With Q ≈ 6–7 the energy falls by e in about one cycle.",
        "<b>Breathing</b>, the size change Derrick's balance fixes, is a broad resonance near 0.45–0.5 c/ℓ, gone in "
        "about two cycles (Q ≈ 2). <b>Twisting</b> hardly rings: no tone of its own.",
        "The frequency is set by the knot's size: f ≈ 1.4 c/(2πℓ), ω·(core radius)/c ≈ 1.7. For ℓ = 1 m it would be "
        "67 MHz, for ℓ = 1 fm, 6.7×10<super>22</super> Hz. The squashing-to-breathing ratio, about 3, needs no "
        "constant at all.",
    ])

    r.h1("Validation")
    r.table([
        ["test", "checks", "held to"],
        ["TestEtherCurvatureIsKnottedLight", "B<sub>n</sub> of Rañada's map = −16 × his knotted B, at four points",
         "10<super>−6</super>"],
        ["TestGradient", "analytic gradients of the 3-D and axisymmetric energies against finite differences",
         "10<super>−6</super>; 10<super>−5</super>"],
        ["TestAxiEqualsThreeD", "axisymmetric and 3-D energies of the same field (E4 converges as h²)", "1 %; 6 %"],
        ["TestVibrations", "coarse grid: two zero modes at ω² ≈ 0; a squashing push rings at ω ≈ 1.45 c/ℓ "
         "(zero crossings)", "0.02; 0.15"],
        ["TestEtherHoldsKnot", "coarse grid: both terms near the virial balance with a core ~ℓ; magnetic only spreads; "
         "stiffness only shrinks, then unties", "0.2; ×1.1; ×0.9, 5 %"],
    ], [0.3, 0.52, 0.18])

    r.h1("Limits of the model")
    r.bullets([
        "<b>κ2 and κ4 are not fixed.</b> The knot's size is ℓ = √(κ4/κ2) and its energy 547.5·√(κ2κ4); no measured "
        "“ether stiffness” exists to set them. Matching one particle's mass and size uses both numbers, so the "
        "first prediction would be a ratio between two particles (e.g. the Q = 2 knot against Q = 1: Battye &amp; "
        "Sutcliffe 835 against 504).",
        "<b>It is not electromagnetism.</b> n is an ordinary field with its own energy. B<sub>n</sub> has the form of a "
        "magnetic field, but how n would couple to electric charges, and whether such a knot carries charge or "
        "spin ½, is open.",
        "Only axially symmetric relaxation and vibrations are computed; the 8 % gap to Battye &amp; Sutcliffe is not "
        "resolved.",
    ])

    r.h1("Status")
    r.p("H29 step 2 is <b>reproduced</b>: a nonlinear ether whose order parameter never vanishes holds the knot at a "
        "fixed size and energy, with its linking number intact and no inflow needed. That is the standard result of "
        "Faddeev &amp; Niemi (1997), now reproduced in kronspice and tied to the knotted light of step 1. Next steps: "
        "the Q = 2 knot (a prediction with no new constant), vibrations that break the axial symmetry, and a Faddeev "
        "ether coupled to the field solver <i>em/fdtd.go</i>.")

    r.h1("Reproducing the results")
    r.code("go test -short ./examples/ether_knot     # field link, gradients, 2-D vs 3-D (< 1 s)\n"
           "go test ./examples/ether_knot            # + the relaxation and controls (~15 s)\n"
           "go run ./examples/ether_knot             # grid study, controls, 3-D check (~30 min on 2 cores)\n"
           "go run ./examples/ether_knot -ring examples/ether_knot/vib     # vibrations, ~8 min on 2 cores\n"
           "python3 examples/ether_knot/ring.py examples/ether_knot/vib/ring_*.json\n"
           "python3 examples/ether_knot/modes.py examples/ether_knot/vib/ops_b16_h0.15.txt\n"
           "python examples/ether_knot/viz/make_report.py examples/ether_knot/ether_knot_report.pdf")
    r.small("Sources: Faddeev &amp; Niemi, Nature 387, 58 (1997); Battye &amp; Sutcliffe, PRL 81, 4798 (1998); Babaev, "
            "Faddeev &amp; Niemi, PRB 65, 100512 (2002); Kleckner, Kauffman &amp; Irvine, Nature Phys. 12, 650 (2016); "
            "Kiessling, arXiv:1107.2333; ATLAS, Nature Phys. 13, 852 (2017); Derrick, J. Math. Phys. 5, 1252 (1964); "
            "Rañada, Lett. Math. Phys. 18, 97 (1989); Hua &amp; Sarkar, IEEE Trans. ASSP 38, 814 (1990); Volovik, The "
            "Universe in a Helium Droplet (2003).")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "ether_knot_report.pdf"))
