# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Kelvin's 1899 charged-ring molecule as a PDF report.

Run from the repository root:
    python examples/kelvin_ring/viz/make_report.py examples/kelvin_ring/kelvin_ring_report.pdf

Numbers come from examples/kelvin_ring/README.md (the example has no
run_output.txt and no figures). Layout from internal/pdfreport/kronreport.py.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401


def build(out):
    r = Report("examples/kelvin_ring", "Kelvin's charged-ring molecule",
               "Induction and Kelvin's own Zeeman test of his 1899 ring, with every charge under the Lorentz force. "
               "CLAUDE.md §6 (Kelvin); the electron test bench for H13, H21 and H23",
               "CLAUDE.md §6 (Kelvin)",
               "24 charges per half, RK4")

    r.h1("Summary")
    r.p("Lord Kelvin's \"Magnetism and Molecular Rotation\" (1899) pictures a molecule as a rigid charged ring filled "
        "with an oppositely charged fluid that can slide round it. This example builds that ring with standard "
        "electrodynamics only, every charge under the Lorentz force, and runs the two things Kelvin did with it: "
        "induction when a field is switched on, and his own Zeeman test, which made him give the model up. It is the "
        "electron test bench on which <i>space_vortex</i> (H13), <i>kron_electron</i> (H21) and "
        "<i>photon_electron</i> (H23) are measured.")
    r.bullets([
        "<b>Induction matches Kelvin and Larmor.</b> Tube and fluid gain moments of momentum "
        "−/+1.11207×10<super>−10</super> kg·m²/s = −/+eΦ/2π, and the fluid turns at Larmor's rate "
        "eB̄/2w = 4.44829×10<super>−6</super> rad/s.",
        "<b>The flux is checked two ways:</b> 6.98736×10<super>−4</super> Wb by Biot–Savart against "
        "6.98732×10<super>−4</super> Wb from Maxwell's elliptic-integral mutual inductance.",
        "<b>The ledger closes:</b> source work −4.94683×10<super>−16</super> J = kinetic 4.94683×10<super>−16</super> J "
        "+ mutual-field −9.89365×10<super>−16</super> J, residual −8×10<super>−16</super>.",
        "<b>Kelvin's gyrostat as written broadens the Zeeman line</b> (splitting (ω<sub>L</sub>/2)·cos θ; ensemble "
        "FWHM 0.0188), which is his stated reason for discarding the model.",
        "<b>The full ring gives sharp lines</b> (splitting 0.0243 … 0.0260 at every orientation, FWHM 0.0037), "
        "removing that objection, <b>but still fails Zeeman's test:</b> its light along B is nearly unpolarised "
        "(−0.02 / −0.03 against the measured ±1), and the splitting is set by the total mass w + w′, not by e/m of "
        "the electron.",
    ])

    r.h1("Background")
    r.p("The ring is a non-conducting tube of radius a carrying +e (mass w′), filled with an incompressible fluid "
        "carrying −e (mass w) that slides without friction. Kelvin used it in §§1–3 of his paper to describe "
        "induction, and in §§5–7 to ask whether such molecules give the Zeeman effect. He discarded it after "
        "Zeeman's splitting (p. 357). The standard classical model beside it is Lorentz's bound charge, which gives "
        "the normal Zeeman effect because all of its charge has one charge-to-mass ratio (the condition of Larmor's "
        "theorem).")

    r.h1("The model")
    r.bullets([
        "Twenty-four charges per half (tube and fluid), integrated with fourth-order Runge–Kutta (RK4).",
        "<b>Induction (Kelvin §§1–3).</b> A coaxial coil's current is ramped on, held and ramped off. Its field comes "
        "from Biot–Savart (<i>em.Loop</i>), the induced E from Faraday's law, and the force on each charge from "
        "<i>em.LorentzVector</i>. The flux is checked against Maxwell's elliptic-integral mutual inductance "
        "(<i>em.MutualCircles</i>).",
        "<b>Zeeman test (Kelvin §§5–7).</b> The ring is held by isotropic angular springs (ω<sub>0</sub> = 1) in a "
        "steady field, kicked, and the spectrum of its motion and of its radiating magnetic moment is taken over "
        "twelve orientations. Three models are compared: Lorentz's bound charge, Kelvin's gyrostat as he wrote it "
        "(fixed internal spin along the axis), and the full ring (fluid free to slip, rigid-body rotation with "
        "<i>hyper/quat</i>, nothing linearised).",
    ])

    r.h1("Results")
    r.h2("Part 1: induction")
    r.p("Ring a = 5 cm, ±1 µC, tube and fluid 10 g each; coil radius 10 cm, 5 cm below, 20 kA·turns, 1 s ramps.")
    r.table([
        ["quantity", "simulated", "expected"],
        ["flux through the ring", "6.98736×10<super>−4</super> Wb (Biot–Savart)", "6.98732×10<super>−4</super> Wb (Maxwell, elliptic)"],
        ["moment of momentum, tube / fluid", "−1.11207×10<super>−10</super> / +1.11207×10<super>−10</super> kg·m²/s",
         "−/+eΦ/2π = 1.11207×10<super>−10</super> (Kelvin's corrected eAM/2π)"],
        ["net angular impulse on the ring", "−5.5×10<super>−30</super>", "0: the halves cancel; the coil gets none back"],
        ["fluid angular velocity", "4.44829×10<super>−6</super> rad/s", "Larmor's eB̄/2w = 4.44829×10<super>−6</super> rad/s"],
        ["ring moment", "−1.11×10<super>−14</super> A·m² (diamagnetic), 1.44×10<super>−14</super> N repulsion",
         "Langevin diamagnetism"],
        ["ledger", "source work −4.94683×10<super>−16</super> J = kinetic 4.94683×10<super>−16</super> J + mutual-field "
         "−9.89365×10<super>−16</super> J", "residual −8×10<super>−16</super>"],
        ["after switch-off", "tube 1.5×10<super>−26</super>, fluid −1.5×10<super>−26</super> kg·m²/s; net source work "
         "−5×10<super>−32</super> J", "everything at rest"],
    ], [0.25, 0.42, 0.33])
    r.p("Kelvin's ring is a mechanical picture of Larmor's theorem. At atomic scale (electron charge and mass, 1 T) "
        "the same fluid turns at 8.794×10<super>10</super> rad/s.")

    r.h2("Part 2: Kelvin's Zeeman test")
    r.p("ω<sub>0</sub> = 1, ω<sub>L</sub> = eB/2w = 0.05, 12 orientations, spectral resolution 2π/T = 0.0026. "
        "θ is the angle between the ring's axis and B.")
    r.table([
        ["model", "splitting vs cos θ", "ensemble upper line (FWHM)", "circular polarisation along B, radiated (upper / lower)",
         "across B: centre parallel, sides perpendicular to B"],
        ["Lorentz bound charge", "0.0513 at every θ", "sharp, 0.0038", "−1.00 / +1.00", "100 % / 100 %"],
        ["Kelvin gyrostat as written", "(ω<sub>L</sub>/2)·cos θ: 0 … 0.0243", "broad, 0.0188", "+0.72 / −0.71", "67 % / 59 %"],
        ["Kelvin ring, full Lorentz force", "0.0243 … 0.0260 at every θ", "sharp, 0.0037",
         "−0.02 / −0.03 (motion −0.47 / +0.46; analytic 5/11)", "64 % / 71 %"],
    ], [0.2, 0.2, 0.16, 0.26, 0.18])
    r.small("In the last column, \"centre parallel, sides perpendicular\" is the share of the centre line polarised parallel to B and "
            "of the side lines polarised perpendicular to it, viewed across the field.")
    r.bullets([
        "<b>The gyrostat</b> reproduces Kelvin's reason for discarding the model (p. 357): with a fixed spin the "
        "splitting goes as cos θ, so randomly oriented molecules broaden the line instead of splitting it.",
        "<b>The full ring</b> is new. When the ring tilts, the flux through it changes and the fluid's spin follows, "
        "since its canonical angular momentum w·a²·u − eΦ/2π is conserved (u is the fluid's angular velocity). The "
        "gyroscopic vector becomes G = (ea²/2)(B − 2(n·B)n), the reflection of B in the ring's plane (n is the ring "
        "axis), whose length does not depend on orientation (<i>GyroVector</i>). Every ring splits by eB/2(w + w′) and "
        "the lines are sharp, which removes Kelvin's and Larmor's objection.",
        "<b>It still fails Zeeman's test.</b> The modes circle G, not B, so the light from the ring's own moment is "
        "nearly unpolarised along the field where Zeeman measured ±1, and the splitting is set by the total mass "
        "w + w′, not the electron's e/m. A ring of two opposite charges with different charge-to-mass ratios cannot "
        "give the normal Zeeman effect; Lorentz's single charge can.",
    ])

    r.h1("Consequence for the electron test bench")
    r.p("Any toroidal model of the electron must give a splitting set by e/m<sub>e</sub> and the observed circular "
        "polarisations before the anomalous Zeeman effect, g and spin ½ are even asked for. That is the bar applied "
        "to the ring models of H13, H21 and H23.")

    r.h1("Validation")
    r.bullets([
        "<i>main_test.go</i>: induction identities and ledgers; the Zeeman test for all three models; |G| constant.",
        "Twin check: flux by Biot–Savart against Maxwell's elliptic mutual inductance (6.98736 against "
        "6.98732 ×10<super>−4</super> Wb).",
        "Energy ledger residual −8×10<super>−16</super> J on 4.94683×10<super>−16</super> J of kinetic energy.",
    ])

    r.h1("Limits of the model")
    r.bullets([
        "Each half is represented by 24 point charges; the fluid is incompressible and frictionless.",
        "The Zeeman test uses isotropic angular springs (ω<sub>0</sub> = 1) to hold the ring, as a stand-in for "
        "whatever binds a molecule, and dimensionless units.",
        "The example has no run_output.txt or figures; the numbers above are those recorded in its README.",
    ])

    r.h1("Reproducing the results")
    r.code("go run ./examples/kelvin_ring\n"
           "go test ./examples/kelvin_ring          # -short skips the Zeeman ensemble")
    r.small("Sources: W. Thomson (Lord Kelvin), \"Magnetism and Molecular Rotation,\" Proc. R. Soc. Edinb. 22, 631–635 "
            "and Phil. Mag. 48, 236–239 (1899); Mathematical and Physical Papers Vol. V, Art. 175, pp. 354–358; notes in "
            "doc/kelvin-references.md. H. A. Lorentz, bound-charge model of the normal Zeeman effect; P. Zeeman, Phil. "
            "Mag. 43, 226 (1897); J. Larmor, Larmor's theorem (1897); J. C. Maxwell, Treatise Vol. II, mutual inductance "
            "of coaxial circles.")
    for f in r.s:
        if getattr(getattr(f, "style", None), "name", "") in ("h1", "h2"):
            f.keepWithNext = 1
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "kelvin_ring_report.pdf"))
