# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""PDF report for examples/vortex_ring: Kelvin's vortex atoms through his ring-electromagnet analogy.

Run from the repository root:

    python examples/vortex_ring/viz/make_report.py examples/vortex_ring/vortex_ring_report.pdf

Numbers come from examples/vortex_ring/README.md (Results, Γ = 1, R = 1). The
example has no run_output.txt and no figures. Layout from
internal/pdfreport/kronreport.py.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401


def build(out):
    r = Report("examples/vortex_ring", "Kelvin's vortex rings as ring electromagnets",
               "Can the engine's inductance and Biot–Savart code compute Kelvin's 1867 vortex atoms through his "
               "ring-electromagnet analogy? — CLAUDE.md §6 (Kelvin)",
               "CLAUDE.md §6 (Kelvin)", "Γ = 1, R = 1")

    r.h1("Summary")
    r.p("Kelvin's 1867 \"On Vortex Atoms\" proposes that atoms are closed vortex rings in a perfect liquid, and notes "
        "that the liquid's velocity round such a ring is the magnetic field of the \"corresponding ring electromagnet\", "
        "with the ring's impulse equal to that electromagnet's magnetic moment (pp. 9–10). This example takes the analogy "
        "literally: the circulation Γ plays the current I, the velocity plays B/μ<sub>0</sub>, and the ring's kinetic "
        "energy is ½ρΓ²·L/μ<sub>0</sub>, with ρ the liquid density and L a self-inductance. The engine's inductance and "
        "Biot–Savart code (milestone M9i) then computes the vortex ring directly. It sits beside the charged ring of "
        "<i>kelvin_ring</i> on the electron test bench.")
    r.bullets([
        "<b>Ring speed from the inductance matches Kelvin's formula</b> to −0.97 %, −0.26 % and −0.07 % for cores "
        "a/R = 0.20, 0.10 and 0.05; the difference falls as (a/R)², the order of the terms the thin-ring formula drops.",
        "<b>The analogy's identities hold:</b> the flow at the ring's centre is 0.50000 = Γ/2R (as the field at a loop's "
        "centre is μ<sub>0</sub>I/2R), and the impulse ρΓπR² is the counterpart of the moment IπR² = 3.14159.",
        "<b>Two coaxial rings leapfrog:</b> 21 passes in 60 time units, period 2.92, radii swinging 0.730 … 1.211, the "
        "pair advancing at 0.4990 against 0.3288 for a single ring.",
        "<b>Invariants hold:</b> impulse drift 2.3×10<super>−5</super>, energy drift 2.1×10<super>−5</super>, which also "
        "checks that the Biot–Savart velocities and the elliptic mutual inductance agree.",
        "<b>Same shape, different dynamics:</b> a vortex ring moves itself; a current loop with the identical field does "
        "not. The analogy is static.",
    ])

    r.h1("The construction")
    r.h2("1. Ring speed")
    r.p("The self-inductance of a ring conductor of radius R and round core radius a with uniform current is Maxwell's "
        "mutual inductance of coaxial circles (Treatise Art. 701) averaged over pairs of filaments in the core section "
        "(Gauss–Legendre in ρ², midpoint in φ). The ring's speed is U = dE/dP at fixed circulation and core volume "
        "(a²R constant), with energy E = ½ρΓ²L/μ<sub>0</sub> and impulse P = ρΓπR². It is compared with Kelvin's "
        "U = Γ/4πR·(ln 8R/a − ¼) (Vol. IV Art. 3, as quoted by Saffman §10.2).")
    r.h2("2. Leapfrogging")
    r.p("Two coaxial rings of equal Γ (Helmholtz 1858; Tait's smoke rings, Kelvin pp. 2, 11–12) are integrated with "
        "RK4; each moves at its self speed plus the other ring's Biot–Savart flow at its centre line. The energy is the "
        "self terms plus the mutual-inductance interaction, so impulse and energy must stay constant.")
    r.h2("3. Where the analogy stops")
    r.p("The flow and the field share their shape but not their dynamics; see below.")

    r.h1("Results (Γ = 1, R = 1)")
    r.h2("Self-inductance and ring speed")
    r.table([
        ["Core a/R", "L/μ<sub>0</sub>R (numeric)", "Thin ring ln 8R/a − 7/4", "U (numeric)", "Kelvin's U", "Difference"],
        ["0.20", "1.95837", "1.93888", "0.27100", "0.27366", "−0.97 %"],
        ["0.10", "2.63730", "2.63203", "0.32798", "0.32882", "−0.26 %"],
        ["0.05", "3.32624", "3.32517", "0.38369", "0.38398", "−0.07 %"],
    ], [0.12, 0.18, 0.22, 0.16, 0.16, 0.16])
    r.p("The numeric inductance approaches the thin-ring value ln(8R/a) − 7/4 as the core shrinks, and the speed "
        "difference falls as (a/R)².")
    r.h2("Leapfrogging rings")
    r.p("Two rings of core 0.1, 0.6 apart, integrated for 60 time units:")
    r.table([
        ["Quantity", "Value"],
        ["Passes", "21, period 2.92"],
        ["Radii swing", "0.730 … 1.211"],
        ["Pair advance", "0.4990 (single ring 0.3288)"],
        ["Impulse drift", "2.3×10<super>−5</super>"],
        ["Energy drift", "2.1×10<super>−5</super>"],
    ], [0.4, 0.6])
    r.h2("Same shape, different dynamics")
    r.p("A vortex ring travels because its vorticity is carried along by the flow it makes. A current loop with the "
        "identical magnetic field is not moved by that field: in Maxwell's equations the ring electromagnet stays put. "
        "A self-propelled toroidal <i>electromagnetic</i> structure has to come from the full field equations, not from "
        "this static analogy; the standard case is Hellwarth and Nouchi's flying doughnut (<i>em/doughnut.go</i>), which "
        "travels at c.")
    r.p("A vortex ring in a neutral liquid has no charge and no coupling to an external field, so the Zeeman bench of "
        "<i>kelvin_ring</i> cannot be applied to it until a charge–vortex coupling is specified.")

    r.h1("Validation")
    r.bullets([
        "Speed error against Kelvin's thin-ring formula shrinks as a² (test in <i>main_test.go</i>).",
        "Centre-flow identity: 0.50000 = Γ/2R.",
        "Leapfrog invariants: impulse and energy constant to 2.3×10<super>−5</super> and 2.1×10<super>−5</super> over "
        "60 time units.",
    ])

    r.h1("Limits")
    r.bullets([
        "Inviscid, incompressible liquid; round cores with uniform vorticity; coaxial rings only.",
        "The analogy maps the shape of the flow onto a magnetic field; it carries no electromagnetic dynamics and no "
        "charge, so it says nothing about vortex atoms' spectra.",
    ])

    r.h1("Reproducing the results")
    r.code("go run ./examples/vortex_ring\n"
           "go test ./examples/vortex_ring")
    r.small("Sources: W. Thomson (Lord Kelvin), \"On Vortex Atoms,\" Proc. R. Soc. Edinb. 6 (1867); Mathematical and "
            "Physical Papers Vol. IV (1910), Art. 1, pp. 1–12, Art. 3, p. 67; H. von Helmholtz, J. reine angew. Math. 55, "
            "25 (1858); P. G. Saffman, Vortex Dynamics (1992), §10.2; J. C. Maxwell, Treatise, Vol. II, Art. 701.")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "vortex_ring_report.pdf"))
