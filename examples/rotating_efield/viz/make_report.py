# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Report: can a rotating electric field make a magnetic field? as a PDF.

python examples/rotating_efield/viz/make_report.py examples/rotating_efield/rotating_efield_report.pdf

Numbers come from examples/rotating_efield/README.md (the example has no
run_output.txt and no figures). Layout from internal/pdfreport/kronreport.py.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401


def build(out):
    r = Report("examples/rotating_efield", "A rotating electric field",
               "Can a rotating E field make a B field? Real fields, and the magnetic-charge case labelled "
               "hypothetical — CLAUDE.md §3f, register entry H2",
               "CLAUDE.md §3f, H2", "milestone M9n")

    r.h1("Summary")
    r.p("This is the inverse of the homopolar dynamo (<i>examples/homopolar</i>): there a conductor spinning in a "
        "magnetic field B gives a voltage; what does a <i>rotating electric field</i> E give? The example keeps two "
        "questions apart. First the real physics: a magnetic field comes from moving charge (Rowland 1876) and from "
        "a changing electric field (the displacement current of the Ampère–Maxwell law), checked with the "
        "quaternion field operators of milestone M9h. Then, opt-in and labelled, the hypothetical physics: the force "
        "a magnetic charge would feel in these fields (the dual pair of Joot's Theorem 3.16), kept out of the real "
        "momentum ledger.")
    r.bullets([
        "<b>A rotating E pattern by itself makes no B.</b> B comes from the charge that moves, or from ∂E/∂t where "
        "the field really changes in time, exactly as the Ampère–Maxwell law says.",
        "<b>Rowland disk</b> (20,000 rpm): B<sub>z</sub> = 1.316×10<super>−9</super> T at the centre and "
        "4.496×10<super>−10</super> T at z = 5 cm, equal to the analytic values.",
        "<b>Charged ring</b>: the same static E pattern at rest and spinning, but B at the centre is 0 T at rest and "
        "1.26×10<super>−8</super> T spinning (the mirror of the Faraday-disk paradox).",
        "<b>Röntgen–Eichenwald disk</b> (dielectric spinning in a fixed E): B<sub>r</sub> = −3.46×10<super>−10</super> T "
        "against the thin-disk −3.5×10<super>−10</super> T. The 1/c² is why the electric mirror of the homopolar "
        "effect is nanotesla where the homopolar one is volts.",
        "<b>Hypothetical magnetic charge</b> (q<sub>m</sub> = 1 A·m, none ever observed): forces of order "
        "10<super>−9</super> N, booked only in the hypothetical column; the real net force stays 0.",
    ])

    r.h1("What is being asked")
    r.p("<b>Standard prediction.</b> By the Ampère–Maxwell law (Maxwell's Treatise, Vol. II, Art. 619: "
        "4π C = curl H in his quaternion notation, where the total current C includes the displacement current ∂D/∂t), any time-varying electric field produces a "
        "magnetic field, and moving charge produces one (Rowland). A field pattern that only appears to rotate while "
        "the charges stay still produces none.")
    r.p("<b>Hypothesis (register entry H2: magnetic charge, dual-symmetric Maxwell).</b> The dual pair of the Lorentz "
        "force would push a magnetic charge q<sub>m</sub> with q<sub>m</sub>(B − v×E/c²). No magnetic monopole has "
        "been observed, so in normal mode q<sub>m</sub> = 0 and the dual pair gives zero force. With "
        "<i>em.Hypothetical{MagneticCharge: true}</i> the example computes what such a charge would feel, labelled "
        "hypothetical and kept in a separate ledger column (CLAUDE.md §3f).")

    r.h1("The five real cases")
    r.table([
        ["case", "set-up", "result", "expected / check"],
        ["1a. Rowland disk", "σ = 10<super>−5</super> C/m², R = 10 cm, 20,000 rpm; B<sub>z</sub> on the axis by "
         "summing 400 current rings",
         "1.316×10<super>−9</super> T at z = 0; 4.496×10<super>−10</super> T at z = 5 cm",
         "analytic (μ<sub>0</sub>σω/2)[(R² + 2z²)/√(R² + z²) − 2|z|]: 1.316×10<super>−9</super>, "
         "4.496×10<super>−10</super> T"],
        ["1b. Rotating uniform E", "1 MV/m rotating at 1 MHz, as between quadrature-driven plate pairs",
         "|B| = 1.75×10<super>−6</super> T at r = 5 cm",
         "(μ<sub>0</sub>ε<sub>0</sub>/2)·E·ω·r; Ampère–Maxwell residual 3.3×10<super>−10</super>, Faraday residual "
         "4.2×10<super>−7</super> (the neglected (ωr/c)² term)"],
        ["1c. Charged ring", "1 µC on 5 cm, at rest and spinning at 1000 rev/s",
         "B at centre 0 T at rest, 1.26×10<super>−8</super> T spinning",
         "the E pattern is identical and static in both; B follows the moving charge"],
        ["1d. Rotating charge pair", "±1 nC on a 10 cm rod at 1000 rev/s",
         "|B| = 8.99×10<super>−12</super> T on the axis 2.5 cm up",
         "the E pattern really rotates; Ampère–Maxwell residual 1.3×10<super>−4</super> of the displacement term "
         "(stencil truncation)"],
        ["1e. Röntgen– Eichenwald", "ε<sub>r</sub> = 3 disk, R = 10 cm, 1 mm thick, 50 rev/s in a fixed 1 MV/m",
         "B<sub>r</sub> = −3.46×10<super>−10</super> T at r = 5 cm (ring sum)",
         "thin-disk −(ε<sub>r</sub> − 1)·E·v/c² = −3.5×10<super>−10</super> T"],
    ], [0.18, 0.26, 0.21, 0.35])
    r.h2("Case 1e: the electric mirror of the homopolar disk")
    r.p("A dielectric moving in an electric field is magnetized, M = P × v (Minkowski; Jackson §11.10). The ring sum "
        "uses only the bound face charge ε<sub>0</sub>(ε<sub>r</sub> − 1)E carried round as two opposite surface "
        "currents and Biot–Savart's μ<sub>0</sub>, so the 1/c² = μ<sub>0</sub>ε<sub>0</sub> comes out by itself, not "
        "by hand. Its deficit against the thin-disk value falls linearly in h/R (h the thickness): the test checks "
        "&lt; 1.5 % at h/R = 0.01 and that halving h halves it. That 1/c² is why the mirror effect is nanotesla where "
        "the homopolar one is volts: 1 T at 100 m/s gives 100 V/m, while 1 MV/m at 100 m/s gives about 1 nT.")
    r.h2("Answer")
    r.p("A rotating E pattern by itself makes no B. B comes from the charge that moves (1a, 1c, 1d, 1e) or from ∂E/∂t "
        "where the field really changes in time (1b, 1d), exactly as the Ampère–Maxwell law says; the quaternion "
        "residuals confirm the fields are solutions of Maxwell's equations.")

    r.h1("Hypothetical: a magnetic charge in these fields")
    r.p("Opt-in with <i>em.Hypothetical{MagneticCharge: true}</i>, magnetic charge q<sub>m</sub> = 1 A·m. No magnetic "
        "charge has ever been observed (CLAUDE.md §3f). These numbers go to the hypothetical column of "
        "<i>em.MomentumLedger</i> only.")
    r.table([
        ["situation", "hypothetical force", "real net force"],
        ["at rest at the Rowland disk's centre", "1.32×10<super>−9</super> N (q<sub>m</sub>·B)",
         "0 (the momentum ledger's real column stays empty)"],
        ["moving at 100 m/s through 1 MV/m", "1.11×10<super>−9</super> N (q<sub>m</sub>·v×E/c²)", "—"],
    ], [0.38, 0.3, 0.32])

    r.h1("Validation")
    r.bullets([
        "Each case is tested against its analytic value (<i>main_test.go</i>).",
        "Field solutions are checked with the M9h quaternion operators: Ampère–Maxwell and Faraday residuals as "
        "listed in the table (3.3×10<super>−10</super> and 4.2×10<super>−7</super> for 1b; 1.3×10<super>−4</super> "
        "of the displacement term for 1d, from stencil truncation).",
        "Röntgen–Eichenwald: deficit against the thin-disk value &lt; 1.5 % at h/R = 0.01, halving with h.",
    ])

    r.h1("Limits of the model")
    r.bullets([
        "Case 1b keeps the quasi-static field; the Faraday residual is the neglected (ωr/c)² term.",
        "The magnetic-charge results are hypothetical (H2, untested): they show what a monopole would feel, not "
        "that one exists.",
    ])

    r.h1("Status and reproducing the results")
    r.p("The real cases are standard electromagnetism and agree with their analytic values. Register entry H2 "
        "(magnetic charge, dual-symmetric Maxwell) stays <b>untested</b>: this example only books what a magnetic "
        "charge would feel, in the hypothetical column.")
    r.code("go run ./examples/rotating_efield\n"
           "go test ./examples/rotating_efield")
    r.small("Sources: H. A. Rowland, Am. J. Sci. 15, 30 (1878); W. C. Röntgen, Ann. Phys. 35, 264 (1888); "
            "A. Eichenwald, Ann. Phys. 11, 1 (1903); J. D. Jackson, <i>Classical Electrodynamics</i>, 3rd ed., §11.10; "
            "Panofsky &amp; Phillips, <i>Classical Electricity and Magnetism</i>, §18-3; J. C. Maxwell, <i>Treatise</i>, "
            "Vol. II, Art. 619; P. Joot, <i>Geometric Algebra for Electrical Engineers</i>, Theorem 3.16 (reference only).")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "rotating_efield_report.pdf"))
