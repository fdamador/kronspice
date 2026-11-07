# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""H22 report: "quantumtricity" — circuit theory as U(1) gauge theory on Kron's network, then SU(3).

python examples/quantumtricity/viz/make_report.py examples/quantumtricity/quantumtricity_report.pdf

Numbers come from README.md (the example has no run_output.txt and no
figures) and the H22 register entry in CLAUDE.md. Layout from
internal/pdfreport/kronreport.py.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401
from reportlab.platypus import KeepTogether  # noqa: E402


def build(out):
    r = Report("examples/quantumtricity", "Circuit electricity and the quark model as gauge theories",
               "Is the quark model circuit electricity at the quantum level? — register entry H22 (the user's "
               "\"quantumtricity\" hypothesis)",
               "register entry H22", "work in progress")

    r.h1("Summary")
    r.p("Register entry H22 holds the user's hypothesis that the Ohm's-law wheel (P, V, I, R) and the Standard Model "
        "wheel describe one structure, the quark model being the principles of electricity at the quantum level. As "
        "stated it makes no prediction, so it is tested in a definite form. The one established link is that both are "
        "gauge theories on a graph: circuit theory is the U(1) (abelian) theory on Kron's network, and the lattice quark "
        "model is the SU(3) (non-abelian) theory on the same kind of graph. Step 1 checks the abelian dictionary; step 2 "
        "searches for any circuit quantity that composes the non-abelian way.")
    r.bullets([
        "<b>Step 1 passes (standard):</b> on a two-mesh RLC circuit in Kron's mesh frame, KVL is zero holonomy of the "
        "branch voltage (residual 2×10<super>−16</super>) and KCL is the rate of the Gauss generator (10<super>−17</super>).",
        "<b>The U(1) lattice field is Yee's FDTD:</b> on a 21×17×6 grid it equals <i>em.FDTD</i> edge by edge to "
        "6.6×10<super>−13</super> of max E over 4000 steps of the TM110 mode; energies agree to 8.6×10<super>−15</super> of W<sub>0</sub>.",
        "<b>Step 2 is null so far:</b> SU(3) mesh loops compose with defect 0.75 on average (min 0.37 over 200 random fields), "
        "0 for the abelian (Cartan) part; every circuit quantity composes with defect ≤ 3×10<super>−16</super>, the round-off of "
        "addition, in the linear circuit and in a diode pair.",
        "<b>The quantum level of V, I and R is electrons and photons:</b> K<sub>J</sub> = 2e/h = 483 597.848 4 GHz/V, "
        "R<sub>K</sub> = h/e² = 25 812.807 459 Ω, I = e·f — no quark content.",
        "<b>Status: bounded, work in progress.</b> The matter sector (quark ↔ electron at the nodes) is not yet modelled.",
    ])

    r.h1("What is being tested")
    r.h2("The hypothesis (H22)")
    r.p("The quark model and the principles of circuit electricity (V, I, R, P) are one structure, the latter being "
        "electricity at the quantum level (user hypothesis, 2026-11-01). The name \"quantumtricity\" is the user's and is "
        "used for this example only, never as an established term.")
    r.h2("The testable form and the standard prediction")
    r.bullets([
        "<b>Prediction of the hypothesis:</b> a circuit quantity exists that composes as SU(3) link variables do, i.e. one "
        "whose loop products around two meshes meeting at a node do not commute.",
        "<b>Standard prediction:</b> none exists. Circuit theory is the classical, abelian (U(1)) limit of a lattice gauge "
        "theory and the quark model is the non-abelian (SU(3)) one on the same kind of graph; every Kirchhoff sum commutes, "
        "so the resemblance ends at \"both are gauge theories on a network\". The quantum level of V, I and R is set by "
        "electrons and photons through e and h; conduction current is electrons only, and quarks fix only the lattice charge Z.",
    ])

    r.h1("The dictionary: both are gauge theories on a graph")
    r.p("A lattice gauge theory (Wilson 1974) has a group element on every link of a graph, a gauge transformation at every "
        "node, and gauge-invariant quantities on closed loops. Kron's network has the same three objects: the branch, the "
        "junction and the mesh. <i>kron/gauge.go</i> writes the correspondence out (A is Kron's incidence tensor, C his "
        "connection tensor):")
    r.table([
        ["Lattice gauge theory", "Kron's network (kron/)"],
        ["node, link, plaquette", "junction, branch, mesh"],
        ["Grad χ (node → link)", "E = A·v, the incidence tensor"],
        ["Curl a (link → plaquette), holonomy", "e′ = C<super>T</super>·e, the loop sum of the connection tensor"],
        ["Curl(Grad) = 0", "A<super>T</super>·C = 0 (connection.go)"],
        ["Div = Grad<super>T</super>", "A<super>T</super>·i, Kirchhoff's current law"],
        ["Curl<super>T</super>", "i = C·i′, branch currents from mesh currents"],
        ["Gauss generator Div Q − ρ", "node charge balance; its time derivative is KCL"],
        ["pure gauge a = Grad χ, zero holonomy", "the branch voltage, E = Grad(v<sub>node</sub>): KVL"],
        ["Hamiltonian Σ Q²/2C + Σ Φ<sub>p</sub>²/2L<sub>p</sub>", "the ledger's ½Cv² and ½Li² (mesh current i<sub>p</sub> = Φ<sub>p</sub>/L<sub>p</sub>)"],
        ["Wilson loop exp(iqΦ<sub>p</sub>/ħ)", "mesh flux linkage in flux quanta"],
    ], [0.45, 0.55])
    r.p("For the group U(1) this is Maxwell's electrodynamics on Kron's 1944 field network (Tonti 2001 is the modern "
        "statement), and in the lumped limit it is circuit theory. For SU(3) it is the lattice form of the quark model "
        "(Kogut &amp; Susskind 1975).")

    r.h1("Step 1 — abelian (standard): circuit theory is U(1) on Kron's network")
    r.p("A two-mesh RLC circuit (V1 R1 L1 L2 R2 C1, three fundamental meshes) is solved in Kron's mesh frame at 1 kHz, "
        "and the identities are read off the solution:")
    r.table([
        ["Identity", "Residual"],
        ["Curl(Grad) = 0 (A<super>T</super>·C = 0)", "0"],
        ["KVL as zero holonomy of the branch voltage", "2×10<super>−16</super>"],
        ["the voltage is pure gauge, E = Grad(v<sub>node</sub>)", "2×10<super>−16</super>"],
        ["mesh EMF = −jω × flux holonomy (Faraday with the flux in the inductors)", "1×10<super>−16</super>"],
        ["KCL as the rate of the Gauss generator, Div I = 0", "1×10<super>−17</super>"],
        ["continuity, Div I<sub>C</sub> = jω·Div Q = jω·ρ<sub>node</sub>", "1×10<super>−18</super>"],
    ], [0.75, 0.25])
    r.p("Mesh 0 links 6.1×10<super>−5</super> Wb at the phasor peak, 2.9×10<super>10</super> flux quanta: a Cooper pair carried "
        "round it would have its Wilson loop turned through 2π that many times. The Hamiltonian's two terms are the "
        "ledger's energies element by element (C1: ½C|v|², L1, L2: ½L|i|²).")
    r.h2("The same field on Kron's 1944 field network against FDTD")
    r.p("The same U(1) field on a 21×17×6 Yee grid closed by conducting walls is stepped beside <i>em.FDTD</i> on the TM110 "
        "cavity mode for 4000 steps:")
    r.table([
        ["Lattice vs FDTD", "Agreement"],
        ["E, edge by edge", "6.6×10<super>−13</super> of max E"],
        ["energy: lattice Hamiltonian vs FDTD ledger", "8.6×10<super>−15</super> of W<sub>0</sub>"],
        ["exact discrete energy drift", "2.0×10<super>−14</super>"],
        ["Gauss generator at interior nodes", "0"],
        ["ω against Yee's own dispersion / the continuum", "1.1×10<super>−7</super> / 5.2×10<super>−4</super>"],
    ], [0.6, 0.4])
    r.p("They agree to round-off because Hamilton's equations on the lattice, dQ/dt = Curl<super>T</super>(Φ/L) and "
        "da/dt = −Q/C, are Yee's leapfrog. The field network and the circuit are one construction: a branch capacitor is an "
        "edge's ε·(dual face)/length, a mesh inductance is a face's μ·(face)/(dual length).")

    n0 = len(r.s)
    r.h1("Step 2 — non-abelian (the hypothesis): SU(3) on the same graph")
    r.p("<i>kron.SU3Field</i> puts an SU(3) element on every branch of the same two-mesh graph. Gauge transformations "
        "U<sub>k</sub> → g<sub>from</sub> U<sub>k</sub> g<sub>to</sub><super>†</super> leave the Wilson-loop traces unchanged "
        "(6×10<super>−15</super>) and turn the Gauss generator covariantly (G<sub>n</sub> → g<sub>n</sub> G<sub>n</sub> "
        "g<sub>n</sub><super>†</super>, 10<super>−14</super>) but not invariantly: it changes by order 1, because it is a matrix, "
        "not a number. Two meshes meeting at a node compose by matrix multiplication; the composition defect is "
        "||W<sub>p</sub>W<sub>q</sub> − W<sub>q</sub>W<sub>p</sub>||/3 for loop elements W<sub>p</sub>, W<sub>q</sub>:")
    r.table([
        ["Field on the links", "Composition defect"],
        ["random SU(3), 200 fields", "mean 0.75, min 0.37"],
        ["diagonal (Cartan) SU(3), the abelian U(1)×U(1) part", "0"],
    ], [0.6, 0.4])
    r.s[n0:] = [KeepTogether(r.s[n0:])]
    r.p("A circuit quantity that composed as SU(3) currents do would show such a defect. The search covers every quantity "
        "the simulator carries on a branch:")
    r.table([
        ["Circuit", "Quantities", "Largest composition defect"],
        ["two-mesh RLC, Kron mesh frame", "voltage, current, charge, flux, complex power", "3×10<super>−16</super> (the round-off of addition)"],
        ["diode pair, operating point (D1 0.669 V 1.73 mA, D2 0.652 V 0.87 mA)", "voltage, current; KCL 9×10<super>−15</super>, KVL 0", "0"],
    ], [0.38, 0.34, 0.28])
    r.p("None found. Every Kirchhoff sum commutes; a nonlinear constitutive law (the diodes) changes the numbers, not the "
        "composition law. The difference between the two wheels is exactly the difference between U(1) and SU(3): the "
        "circuit is the abelian, classical limit.")

    r.h1("The quantum level of V, I and R")
    r.p("Since the 2019 SI the three units of the Ohm's-law wheel are realized from the exact e and h (computed in the run "
        "from CODATA 2018):")
    r.table([
        ["Unit", "Realization", "Source"],
        ["volt", "K<sub>J</sub> = 2e/h = 483 597.848 4 GHz/V", "Josephson 1962"],
        ["ohm", "R<sub>K</sub> = h/e² = 25 812.807 459 Ω", "von Klitzing, Dorda &amp; Pepper 1980"],
        ["ampere", "I = e·f; 1 nA = 6.24×10<super>9</super> electrons/s", "Pekola et al. 2013"],
    ], [0.15, 0.5, 0.35])
    r.bullets([
        "One lepton and one gauge boson of the Standard Model wheel; no quark.",
        "The only quark content a circuit sees is the lattice charge Z of its conductors (copper: 29 × uud = 29 e). "
        "Conduction current is electrons.",
        "The fractional charges e/3 and e/5 measured in the fractional quantum Hall effect (de-Picciotto et al. 1997; "
        "Saminadayar et al. 1997) are collective excitations of electrons, not quarks.",
    ])

    r.h1("Limits and open work")
    r.bullets([
        "<b>The matter sector.</b> The dictionary matches the gauge fields only (photon ↔ gluon). A quark is a charged matter "
        "field at the nodes; its circuit counterpart is the electron, entering the Gauss law as the node charge ρ. Step 2 "
        "carried the quark side only as abstract colour charges ρ<sub>n</sub>. A charge field on the nodes with its hopping "
        "current along the links is the next construction.",
        "<b>Where the wheel's terms land.</b> V is the electric field on a link, I the matter current along it, P = V·I the "
        "work term J·E between matter and field, and R has no counterpart: a gauge theory is lossless, and resistance is "
        "coupling to a heat bath. Still to be written into the dictionary with a test.",
        "<b>The TLM twin.</b> The register asks for agreement with <i>em/tlm.go</i> as well as <i>em/fdtd.go</i>; TLM is not "
        "yet written.",
    ])

    r.h1("Status")
    r.p("Register entry H22: <b>bounded</b>, work in progress. Step 1 passes; step 2 is a null result so far, which records a "
        "bound. The entry stays in the register and the example keeps growing. What would move it: a circuit quantity with "
        "a non-zero composition defect, i.e. one that carries an index transported along branches and whose node sums do "
        "not commute. None of the simulator's quantities does.")

    r.h1("Reproducing the results")
    r.code("go run ./examples/quantumtricity\n"
           "go test ./examples/quantumtricity")
    r.small("Sources: K. G. Wilson, <i>Phys. Rev. D</i> 10, 2445 (1974); J. Kogut, L. Susskind, <i>Phys. Rev. D</i> 11, 395 (1975); "
            "G. Kron, <i>Proc. IRE</i> 32, 289 (1944); E. Tonti, <i>PIER</i> 32, 1 (2001); K. S. Yee, <i>IEEE Trans. Antennas "
            "Propag.</i> 14, 302 (1966); U. Vool, M. Devoret, <i>Int. J. Circuit Theory Appl.</i> 45, 897 (2017); M. Gell-Mann, "
            "<i>Phys. Rev.</i> 125, 1067 (1962); D. Marcos et al., <i>Phys. Rev. Lett.</i> 111, 110504 (2013); B. D. Josephson, "
            "<i>Phys. Lett.</i> 1, 251 (1962); K. v. Klitzing et al., <i>Phys. Rev. Lett.</i> 45, 494 (1980); J. P. Pekola et al., "
            "<i>Rev. Mod. Phys.</i> 85, 1421 (2013); BIPM, <i>SI Brochure</i>, 9th ed. (2019); R. de-Picciotto et al., "
            "<i>Nature</i> 389, 162 (1997); L. Saminadayar et al., <i>Phys. Rev. Lett.</i> 79, 2526 (1997).")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "quantumtricity_report.pdf"))
