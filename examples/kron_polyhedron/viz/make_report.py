# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Kron's field network and his 1962 polyhedron (H20 step 1) as a PDF report.

Run from the repository root:
    python examples/kron_polyhedron/viz/make_report.py examples/kron_polyhedron/kron_polyhedron_report.pdf

Numbers come from examples/kron_polyhedron/README.md (the example has no
run_output.txt and no figures). Layout from internal/pdfreport/kronreport.py.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401


def build(out):
    r = Report("examples/kron_polyhedron", "Kron's field network and his 1962 polyhedron",
               "Kron's equivalent circuit of Maxwell's equations and his curve-fitting polyhedron, checked against "
               "every number he printed. The standard-physics half of register entry H20 (step 1)",
               "register entry H20",
               "")

    r.h1("Summary")
    r.p("Register entry H20 reads Kron's polyhedral waves (1962) as a description of a physical field region. Its "
        "first step is to build the structures and check them against standard physics and against Kron's own "
        "numbers. This example solves Kron's 1945 field network (a cavity and a waveguide lattice) and the "
        "curve-fitting polyhedron of his 1962 paper with the circuit engine, and compares the results with every "
        "number he printed (figures and tables transcribed in <i>doc/kron-figure-catalogue.md</i>).")
    r.bullets([
        "<b>The 1945 cavity matches his network analyzer</b> to a worst 5.9 % of the largest potential (mean 2 %); his "
        "own unbalance was 8 %.",
        "<b>The TE<sub>11</sub> box resonance</b> of the lattice equals the spatial part of Yee's dispersion to "
        "2×10<super>−14</super> and the continuum to 0.36 % (second order); Kron's solve equals MNA to "
        "4×10<super>−10</super>.",
        "<b>Kron 1962 Tables I–VI are reproduced:</b> given columns to 1.9×10<super>−6</super>, the first planar frame "
        "to 8.9×10<super>−6</super>, the oscillatory column to 2.8×10<super>−5</super>.",
        "<b>A planar reference frame is a vector of whole numbers.</b> The genetic algorithm recovers its nine integers "
        "from the printed column alone (17,702 circuit evaluations).",
        "<b>The oscillatory column returns its input.</b> Its frame fits any data, 1000 random sets to "
        "7×10<super>−15</super>, so its agreement shows no estimating power and does not by itself support the "
        "physical reading of H20.",
        "<b>Power input is the error-square</b> (Kron p. 67): a one-mesh circuit's dissipation equals the summed "
        "error-square, ledger residual 4×10<super>−13</super> W.",
    ])

    r.h1("What is being tested")
    r.p("<b>The hypothesis (H20).</b> In his 1962 curve-fitting paper Kron built a k-dimensional polyhedron whose "
        "q-networks are copies of his 1944 Maxwell network, even-dimensional ones magnetic and odd-dimensional ones "
        "dielectric, joined to a dual polyhedron. He wrote that a plasma-filled field region \"crystallizes\" into such "
        "a structure, that it can oscillate self-adaptively, and that crystals realize it. Kron himself called the "
        "physical interpretation \"only tentative\" (p. 66); the hypothesis is the physical reading the user placed in "
        "the register.")
    r.p("<b>Standard prediction.</b> In three dimensions the polyhedron is a discretization of Maxwell's equations "
        "(Kron 1944; Tonti's cell method) and agrees with FDTD and TLM to discretization error; a passive linear "
        "structure has fixed modes; X-ray-excited crystal waves follow Ewald–Laue dynamical diffraction. This "
        "example covers the standard half of step 1: does the network reproduce Kron's printed numbers and standard "
        "field results?")

    r.h1("What was checked against what")
    r.table([
        ["source", "kronspice", "result"],
        ["Kron 1945 Fig. 14, cavity potentials read on the network analyzer", "the lattice solved exactly "
         "(kron.Kron1945Cavity)", "worst 5.9 % of the largest potential, mean 2 %; his own unbalance was 8 %"],
        ["Kron 1945 Fig. 15, unbalanced currents", "recomputed from his Fig. 14 by his eq. (1)",
         "23 of 25 junctions within 0.33 A; the other two look like his slips (below)"],
        ["Balasubramanian 1968 Fig. 15, waveguide matching error", "the TE lattice (kron.Waveguide), with a "
         "closed-form match leaving ripple 1e-16", "+0.69, +1.50, +3.47, +6.43 % against about 1.0, 1.6, 3.4, 6.2 % "
         "read from his curve"],
        ["TE<sub>11</sub> box resonance", "lattice dispersion, the spatial part of Yee's",
         "equal to 2e-14; continuum to 0.36 %, second order; Kron's solve equals MNA to 4e-10"],
        ["Kron 1962 Tables I–VI, given columns", "CurveFit.Characteristics on his printed z", "worst 1.9e-6"],
        ["Tables I–VI, first planar frame", "the integer frame below", "worst 8.9e-6"],
        ["Tables I–VI, oscillatory column", "the closed-form frame below", "worst 2.8e-5"],
        ["Kron 1962 p. 67, power input is the error-square", "one-mesh circuit, MNA ledger and Kron's mesh form",
         "dissipation equals the summed error-square; ledger residual 4e-13 W"],
    ], [0.33, 0.32, 0.35])

    r.h1("The two findings")
    r.h2("A planar reference frame is a vector of whole numbers")
    r.p("Kron never states his frames. With 13 rows and 12 coefficients, a fit has one residual direction w, and "
        "ŷ = y − w(w<super>T</super>y)/(w<super>T</super>w) whichever 12 columns span the model (y is the data "
        "column, ŷ the estimate). The errors in the first planar column have rank one across the six tables, and that "
        "direction is")
    r.code("w = ( -3, 1, 1, 1 |  4, -1, -3, -4, -4, 2 over the branch lengths |  2, -10, -4 over 9 )")
    r.p("The genetic algorithm then recovers these nine integers from the printed column alone, with every candidate "
        "solved as a circuit (17,702 evaluations).")
    r.h2("The oscillatory column returns its input")
    r.p("The oscillatory column equals the given column for all six functions, to rounding. All 13 characteristics "
        "are linear in the four values at the points, so they span only a 4-dimensional subspace. A frame with w "
        "orthogonal to that subspace fits every function exactly, with zero power, and "
        "<i>kron.CurveFit.OscillatoryFrame</i> writes one down:")
    r.p("&nbsp;&nbsp;&nbsp;&nbsp;N = (9, M<sub>2</sub>/2, M<sub>1</sub>/2, −13, −10, M<sub>3</sub>/2), any even M")
    r.bullets([
        "That is why the oscillatory polyhedron \"can estimate satisfactorily even a set of random numbers\", of which "
        "Kron wrote \"the significance of this is not yet understood\" (p. 69): 1000 random sets here fit to "
        "7×10<super>−15</super>. It is also why six \"widely different\" functions fit without adjustment.",
        "A model with 12 coefficients is fitting 4 independent numbers, so this agreement shows no estimating power. "
        "It is a property of the frame, and it does not by itself support the physical reading of H20.",
        "This is an inference from the tables, with three limits: any member of the family matches the printed "
        "column, so which frame Kron used cannot be determined; he wrote that \"no exact fit is possible\", and a "
        "residual below about 10<super>−5</super> cannot be excluded from the printed figures; Kron 1959 and the 1960 "
        "General Electric report have not been read and might state otherwise.",
    ])

    r.h1("Where the GA did not succeed")
    r.p("On Kron's objective, minimum power input, the search over integer frames did not reach zero. In a sweep of "
        "200 × 300 runs, none of 17 completed runs (three mutation steps, up to six seeds each) got there. The best "
        "was 3.7×10<super>−4</super> W, against 90 W for frame 1. The sweep was cut off by its time limit before the "
        "0.2-decade step ran. Zero power is a linear condition, so it is solved by the closed form, and the GA adds "
        "nothing there. <i>TestGAMinimumPower</i> records the shortfall.")

    r.h1("What is not reproduced")
    r.bullets([
        "<b>The second planar frame.</b> Its function rows follow frame 1's residual times a constant 0.3032, but its "
        "divided differences are not one linear map of z. That means several misprints or a different frame per "
        "table; six tables cannot separate the two.",
        "<b>Least-square + linear.</b> With two corrections it is linear to 2×10<super>−5</super>. Two outer branches "
        "come out with the sign reversed and the third exact, but the inner rows have a rank-3 error map, not a "
        "single residual direction, so no frame of this form gives it.",
        "<b>The full 3-D cell with its ideal transformers</b> (Kron 1962 Fig. 1) is built in <i>kron/maxwell3d.go</i> "
        "and its tests, not in this example. Still open from stage B of the catalogue's plan: the TLM twin.",
        "<b>q-networks in general.</b> <i>kron.CurveFit</i> handles this one polyhedron. The second-difference formula "
        "S = −(Δz/L)(A/L) is read off the tables, not derived from Kron's construction.",
    ])

    r.h1("Apparent errors in the printed sources")
    r.p("These are inferences from internal consistency, not published errata.")
    r.table([
        ["where", "printed", "needed / computed", "note"],
        ["Table II, least-square + linear, second function row", "−12.339281", "−12.330283", ""],
        ["Table I, same column, fourth divided difference", "−6.4341227", "−0.64341227", "a slipped decimal point"],
        ["1945 Fig. 15, row 2 column 3", "2.15", "−2.16", "the figure's own e<sub>u</sub> is −0.51 (a dropped minus sign)"],
        ["1945 Fig. 15, row 1 column 0", "2.12", "−5.94 from the printed potentials", "cause unresolved"],
    ], [0.36, 0.14, 0.22, 0.28])
    r.p("Kron's p. 62 claim that a polyhedron \"improves always\" the least-square estimate fails in his own Table I: "
        "z<sub>1</sub> is linear, the least-square fit is exact, and frame 1 moves it from 9.000 to 8.879.")

    r.h1("Status")
    r.p("H20 remains <b>untested</b> in the register: the physical claims (plasma, crystals, self-adaptive "
        "oscillation) have no equations in any source and no test yet. The standard half of step 1 is partly done: the "
        "2-D lattice and the first planar frame reproduce Kron's numbers, the 3-D cell (in <i>kron/maxwell3d.go</i>) "
        "equals Yee's equations, and the TLM twin, the second planar frame and the least-square + linear column are "
        "still open. The oscillatory column's fit, which matches any values including random numbers, is not evidence "
        "for the physical reading.")

    r.h1("Reproducing the results")
    r.code("go run ./examples/kron_polyhedron\n"
           "go test ./examples/kron_polyhedron ./kron")
    r.small("Sources: G. Kron, J. Math. Anal. Appl. 5, 46–69 (1962); G. Kron, J. Appl. Phys. 16, 172–186 (1945); "
            "N. V. Balasubramanian, PhD thesis, Liverpool, 1968; H. Nicholson, ACSE Research Report 26, Sheffield, 1974; "
            "E. Tonti, PIER 32, 1–44 (2001). The example has no run_output.txt or figures; numbers are those recorded "
            "in its README.")
    for f in r.s:
        if getattr(getattr(f, "style", None), "name", "") in ("h1", "h2"):
            f.keepWithNext = 1
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "kron_polyhedron_report.pdf"))
