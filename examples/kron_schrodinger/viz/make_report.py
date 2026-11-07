# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Kron's Schrödinger network (H atom, H2+ bond) as a PDF report.

Run from the repository root:
    python examples/kron_schrodinger/viz/make_report.py examples/kron_schrodinger/kron_schrodinger_report.pdf

Numbers come from examples/kron_schrodinger/run_output.txt (the program's own
output), README.md and export/kron_schrodinger_h2plus_coarse_results.md. The
example has no figures. Layout from internal/pdfreport/kronreport.py.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401


def build(out):
    r = Report("examples/kron_schrodinger", "Kron's Schrödinger network",
               "The electron of a chemical bond (H<sub>2</sub><super>+</super>) solved as an electric network after "
               "Kron 1945, with the H atom as its check. CLAUDE.md §3a.1",
               "CLAUDE.md §3a.1",
               "atomic units; grids h = 0.2, 0.1, 0.05 a0")

    r.h1("Summary")
    r.p("Kron (<i>Phys. Rev.</i> 67, 39, 1945) wrote the Schrödinger equation as an electric network: series "
        "inductors carry the kinetic energy, shunt inductors to ground the potential energy, and shunt capacitors the "
        "energy eigenvalue. The wave function ψ is the junction voltage to ground, and an eigenstate is a resonance of "
        "the network. This example builds that network for H<sub>2</sub><super>+</super>, the simplest chemical bond "
        "(two protons, one electron), with the hydrogen atom as its check; both have exact solutions. It is the "
        "electronic side of the chemistry plan (the vibrational side is <i>water_vibration</i>).")
    r.bullets([
        "<b>Hydrogen atom:</b> E → −0.499728 Ha = −13.5983 eV as h → 0 (order 1.44), against the exact −0.5 Ha = "
        "−13.6057 eV.",
        "<b>H<sub>2</sub><super>+</super> at R = 2 a<sub>0</sub>:</b> electronic energy E → −1.102291 Ha (order 1.46), "
        "against the exact −1.1026342 Ha.",
        "<b>The bond:</b> length R<sub>e</sub> = 1.998 a<sub>0</sub> = 1.057 Å (exact 1.997 a<sub>0</sub> = 1.057 Å) "
        "and energy D<sub>e</sub> = 0.1024 Ha = 2.786 eV (exact 0.1026 Ha = 2.79 eV).",
        "<b>Circuit twin:</b> the same network as an LC circuit in kronspice resonates at E = −1.13848455 Ha on the "
        "coarse grid (h = 0.4 a<sub>0</sub>, 1950 junctions), equal to the driven-network eigenvalue −1.13848455 Ha "
        "(test tolerance 10<super>−7</super>).",
        "<b>Torn solve equals MNA:</b> the 1950-node circuit torn into 4 pieces joined by 179 links matches the full "
        "solve to 3.1×10<super>−15</super> of the largest node voltage, with the same power 0.302832 at phase −90.00°.",
    ])

    r.h1("The construction")
    r.p("The ground state is axially symmetric (m = 0), so the network is laid out on the (ρ, z) half-plane with the "
        "nuclei on the axis. In atomic units, the equation multiplied by ρ so the network is reciprocal is")
    r.p("&nbsp;&nbsp;&nbsp;&nbsp;−½ [∂<sub>ρ</sub>(ρ ∂<sub>ρ</sub>ψ) + ρ ∂<sub>z</sub>²ψ] + ρ V ψ = E ρ ψ,"
        "&nbsp;&nbsp;&nbsp;&nbsp;V = −Σ Z<sub>a</sub> / |r − R<sub>a</sub>|")
    r.p("where Z<sub>a</sub> and R<sub>a</sub> are the nuclear charges and positions. On a grid of spacing h with "
        "junctions at ρ<sub>i</sub> = (i + ½)h (no junction on the axis, so it needs no boundary condition) this is "
        "A ψ = E B ψ with B = diag(ρ<sub>i</sub>). In Kron's network: series branches of admittance "
        "ρ<sub>face</sub>/(2h²) between neighbours, a shunt branch ρ<sub>i</sub>·V<sub>i</sub> at each junction and a "
        "shunt capacitor ρ<sub>i</sub>. The outer boundary is a short circuit (ψ = 0) far from the nuclei, Kron's "
        "termination \"by a short circuit at a point of zero voltage\".")
    r.p("The eigenvalue is found two ways:")
    r.bullets([
        "<b>Kron's driven network.</b> The network is driven at a trial energy σ by sources B·x and its response is the "
        "next x (inverse iteration with the factors computed once); the response converges to the resonance nearest σ "
        "and the Rayleigh quotient gives E. Energies at h, h/2 and h/4 are extrapolated to h → 0 (Richardson, with the "
        "order estimated from the three runs).",
        "<b>A kronspice LC circuit.</b> The same network written as a SPICE netlist in its dual form (capacitors and "
        "inductors exchanged, because Kron's all-inductor network has loops of shorts at DC): series capacitors between "
        "junctions, shunt capacitors ρ(V + V0) and shunt inductors 1/ρ, with V0 shifting the potential so every "
        "element is positive. Every eigenstate is then a resonance at ω² = 1/(E + V0), the ground state the highest, "
        "and kronspice's AC analysis finds it.",
    ])

    r.h1("Results")
    r.p("Finite-difference energies lie below the exact ones (a difference scheme is not variational; the grid "
        "samples the Coulomb wells too deeply) and converge with an order set by the cusps, about 1.45 on these grids.")
    r.h2("Hydrogen atom (check)")
    r.table([
        ["h (a<sub>0</sub>)", "junctions", "E (Ha)"],
        ["0.200", "9,800", "−0.51087970"],
        ["0.100", "39,200", "−0.50382944"],
        ["0.050", "156,800", "−0.50123643"],
        ["h → 0 (order 1.44)", "", "−0.499728 = −13.5983 eV (exact −0.5 Ha = −13.6057 eV)"],
    ], [0.25, 0.2, 0.55], bold_rows=(4,))
    r.h2("H<sub>2</sub><super>+</super> at R = 2 a<sub>0</sub>, electronic energy")
    r.table([
        ["h (a<sub>0</sub>)", "junctions", "E (Ha)"],
        ["0.200", "7,800", "−1.11745397"],
        ["0.100", "31,200", "−1.10780517"],
        ["0.050", "124,800", "−1.10429632"],
        ["h → 0 (order 1.46)", "", "−1.102291 (exact −1.1026342)"],
    ], [0.25, 0.2, 0.55], bold_rows=(4,))
    r.h2("The bond curve")
    r.p("Total energy E(R) + 1/R (the 1/R term is the protons' repulsion), each R extrapolated to h → 0 from h = 0.2 "
        "and 0.1, with a least-squares cubic through the minimum:")
    r.table([
        ["R (a<sub>0</sub>)", "R (Å)", "total energy (Ha)"],
        ["1.6", "0.847", "−0.590490"],
        ["1.8", "0.953", "−0.599866"],
        ["2.0", "1.058", "−0.602291"],
        ["2.2", "1.164", "−0.600531"],
        ["2.4", "1.270", "−0.596271"],
        ["2.8", "1.482", "−0.584110"],
    ], [0.3, 0.3, 0.4])
    r.table([
        ["", "model", "exact"],
        ["bond length R<sub>e</sub>", "1.998 a<sub>0</sub> = 1.057 Å", "1.997 a<sub>0</sub> = 1.057 Å"],
        ["bond energy D<sub>e</sub>", "0.1024 Ha = 2.786 eV", "0.1026 Ha = 2.79 eV"],
    ], [0.3, 0.35, 0.35])

    r.h2("Circuit twin and torn solve")
    r.bullets([
        "Coarse grid, h = 0.4 a<sub>0</sub>, 1950 junctions: kronspice AC resonance E = −1.13848455 Ha; driven-network "
        "eigenvalue on the same grid −1.13848455 Ha. The test requires agreement to 10<super>−7</super>.",
        "The exported netlist (<i>kron_schrodinger_h2plus_coarse.cir</i>, 7831 elements, 1950 nodes) is solved at AC "
        "(0.078178 Hz) by automatic tearing into 4 pieces of 976, 974, 976 and 974 unknowns, joined by 179 links; the "
        "only coupled solve is of order 179, against 3900 unknowns for the whole circuit.",
        "Torn against full (MNA) solve: largest node-voltage difference 3.1×10<super>−15</super> of the largest node "
        "voltage; power taken by the R, L, C branches 0.302832 at phase −90.00° both torn and untorn (power is invariant "
        "under tearing, Roth 1959 §5).",
    ])

    r.h1("Validation")
    r.bullets([
        "<i>main_test.go</i>: convergence to the exact energies; circuit resonance equals the eigenvalue; cubic minimum.",
        "Exact references: Bates, Ledsham and Stewart (1953) and Peek (1965) for H<sub>2</sub><super>+</super>; CODATA "
        "2018 for the hartree and the bohr.",
    ])

    r.h1("Limits")
    r.bullets([
        "Only the axially symmetric (m = 0) ground state is treated, on a 2-D (ρ, z) grid.",
        "Finite-difference energies are not variational and lie below the exact ones; the h → 0 values rely on "
        "Richardson extrapolation with an order (about 1.45) set by the Coulomb cusps.",
        "One electron only. Water has ten interacting electrons: the next step is a self-consistent field "
        "(Hartree–Fock) on a 3-D network, with a Poisson network for the repulsion, which is not done here.",
    ])

    r.h1("Reproducing the results")
    r.code("go run ./examples/kron_schrodinger\n"
           "go test ./examples/kron_schrodinger")
    r.small("Sources: G. Kron, \"Electric circuit models of the Schrödinger equation,\" Phys. Rev. 67, 39–43 (1945); "
            "D. R. Bates, K. Ledsham, A. L. Stewart, Phil. Trans. R. Soc. A 246, 215 (1953); J. M. Peek, J. Chem. Phys. "
            "43, 3004 (1965); CODATA 2018.")
    for f in r.s:
        if getattr(getattr(f, "style", None), "name", "") in ("h1", "h2"):
            f.keepWithNext = 1
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "kron_schrodinger_report.pdf"))
