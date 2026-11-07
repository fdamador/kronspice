# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Report: Einstein's 1928 distant parallelism beside Kron's machines, as a PDF.

python examples/einstein_teleparallel/viz/make_report.py examples/einstein_teleparallel/einstein_teleparallel_report.pdf

Numbers come from examples/einstein_teleparallel/run_output.txt and README.md;
the figure is viz/einstein_teleparallel.png. Layout from internal/pdfreport/kronreport.py.
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
    r = Report("examples/einstein_teleparallel", "Einstein's distant parallelism and Kron's machines",
               "Does Einstein's 1928 teleparallel geometry answer Kron's non-Riemannian machine theory? — "
               "CLAUDE.md §3a, M9c.1", "CLAUDE.md §3a, M9c.1",
               "120 V, 60 Hz motor, pushed to 40 rad/s, 4 s")

    r.h1("Summary")
    r.p("In 1928 Einstein built space-time from a <b>4-bein</b> — four frame legs at every point — with the rule that "
        "legs at different points are parallel (“distant parallelism”). That geometry has no curvature but has "
        "<b>torsion</b> Λ. He proposed that the trace of the torsion, φ<sub>μ</sub>, is the electromagnetic "
        "potential, and that to first order the theory gives both vacuum gravity and Maxwell's equations. This "
        "example asks whether that answers Kron's multi-dimensional, non-Riemannian machine work, using Hoffmann's "
        "single-phase motor (see <i>hoffmann_motor</i>) as the test object.")
    r.bullets([
        "<b>One torsion, three notations:</b> Einstein's Λ equals minus the frame's anholonomy to 6×10<super>−11</super>; "
        "Hoffmann's Ω equals Λ with one index lowered by the inductances to 3×10<super>−17</super>. Kron's "
        "non-Riemannian frame is a teleparallel 4-bein.",
        "<b>Einstein's potential belongs to the frame, not the machine:</b> φ = 0 in Kron's frame (largest "
        "component 7×10<super>−11</super>), but φ<sub>θ</sub> = 0.300 for the same motor in a frame whose rotor legs "
        "also grow as e<super>0.3θ</super>.",
        "<b>Curvature is the real difference:</b> Einstein's connection is flat (largest curvature component 0); the "
        "machine's own connection is curved (largest component 541).",
        "<b>Einstein's first approximation holds:</b> for plane waves his eqs. 7–9 hold to ≤ 2×10<super>−16</super>; "
        "on the light cone φ obeys Maxwell's vacuum equations.",
        "<b>A motor runs from the torsion alone:</b> built only from the 4-bein and its torsion, the kronspice motor "
        "reaches 376.8235 rad/s, equal to Kron's frame II(b) to 2×10<super>−12</super> rad/s. With the torsion "
        "removed it only coasts to 26.8128 rad/s.",
    ])

    r.h1("Background")
    r.p("A. Einstein, “New possibility for a unified field theory of gravitation and electricity,” "
        "<i>Sitzungsber. Preuss. Akad. Wiss.</i> 1928, 224–227 (A. Unzicker's translation). Kron's non-Riemannian "
        "dynamics of rotating machines (1934), presented for physicists by Hoffmann (<i>Rev. Mod. Phys.</i> 21, 535, "
        "1949), writes a machine's equations in non-holonomic frames — frames of current paths that are not the "
        "derivatives of any coordinates — and carries the motion in a torsion Ω (Hoffmann eq. 19). Hoffmann's "
        "motor has exactly four coordinates (three winding charges q1, q2, q3 and the rotor angle θ), so its change "
        "to stationary rotor coils (Hoffmann eq. 35) is a 4-bein, Einstein's case n = 4.")
    r.p("<b>Standard position.</b> Einstein's distant-parallelism unified theory was dropped by 1931; its gravity "
        "part survives as the teleparallel equivalent of general relativity, while electromagnetism from the torsion "
        "trace did not survive (Sauer 2006; Maluf 2013). This example is a comparison of formalisms, not a "
        "register hypothesis.")

    r.h1("Results")
    r.h2("1. Einstein's 4-bein is Kron's frame")
    r.table([
        ["the same torsion, computed as", "agreement"],
        ["Einstein's torsion Λ (his definition), in the frame's components",
         "equals minus the frame's anholonomy (sign convention), to 6×10<super>−11</super>"],
        ["Hoffmann's Ω (<i>kron.Torsion</i>, his eq. 19)",
         "equals Λ with one index lowered by the inductances, Ω<sub>λμν</sub> = L′<sub>νρ</sub>·Λ<super>ρ</super><sub>λμ</sub>, "
         "to 3×10<super>−17</super>"],
        ["torsion of Einstein's own connection", "2Λ, to 10<super>−10</super>"],
    ], [0.45, 0.55])
    r.h2("2. Einstein's electromagnetic potential on the motor")
    r.table([
        ["frame of the same motor", "φ<sub>μ</sub> = Λ<super>α</super><sub>μα</sub> (trace of the torsion, his eq. 2)"],
        ["Kron's stationary rotor coils (a rotation)",
         "(0, 0, 0, 7×10<super>−11</super>), i.e. 0 — a “pure gravitational field” to Einstein"],
        ["rotor legs also growing as e<super>0.3θ</super>", "φ<sub>θ</sub> = 0.300"],
    ], [0.45, 0.55])
    r.p("In a Kron machine the frame is a choice of current paths and both frames describe the same motor, so φ here "
        "is a property of the frame chosen, not a field the machine has. This echoes Einstein's own remark that "
        "“the separation of the gravitational and electromagnetic field seems arbitrary”.")
    r.h2("3. Curvature: the real difference")
    r.table([
        ["connection", "curvature"],
        ["Einstein's distant parallelism, either frame", "largest component 0 — flat; all its geometry is torsion"],
        ["the machine's own (Christoffel symbols of its inductances and inertia)",
         "largest component 541 (mixed units of henry and kg·m²) — not zero; the machine's space is curved"],
    ], [0.45, 0.55])
    r.p("Kron's torsion is that curved space seen from a moving frame (Hoffmann eq. 17); a change of frame cannot "
        "remove curvature. Einstein's torsion <i>replaces</i> curvature. The mathematics of frames is the same; the "
        "geometry is different.")
    r.h2("4. Einstein's first approximation, checked")
    r.p("His linearized field equation (eq. 5) was solved for plane waves e<super>p·x</super>, with x<super>4</super> "
        "imaginary as in his paper. R<sub>αβ</sub> is the Ricci tensor.")
    r.table([
        ["wave", "solutions of eq. 5", "largest φ", "eq. 8", "eq. 9", "eq. 7", "R<sub>αβ</sub> when φ = 0"],
        ["p² ≠ 0", "4", "6.3×10<super>−17</super>", "8×10<super>−17</super>", "9×10<super>−17</super>",
         "1×10<super>−16</super>", "2×10<super>−16</super>"],
        ["p² ≠ 0, complex", "4", "4.2×10<super>−16</super>", "1×10<super>−16</super>", "5×10<super>−17</super>",
         "2×10<super>−16</super>", "3×10<super>−16</super>"],
        ["p² = 0 (light-like)", "10", "0.5", "0", "0", "0", "0"],
    ], [0.19, 0.14, 0.13, 0.12, 0.12, 0.12, 0.18])
    r.p("Away from the light cone the only solutions are coordinate changes (φ = 0, R = 0). On it, φ ≠ 0 obeys ∂φ = 0 "
        "and □φ = 0 — Maxwell's vacuum equations — and with φ = 0 the metric obeys R<sub>αβ</sub> = 0. Einstein's "
        "first-approximation claims hold.")
    r.h2("5. kronspice: the motor built from Einstein's torsion")
    r.p("The motor's machine model was built using only the 4-bein and its torsion: constant inductances L′ (Hoffmann's L-bar, the inductances in the stationary-coil frame) and a "
        "speed matrix G = 2(Ω<sub>θλν</sub> + Ω<sub>νλθ</sub>), with no inductance derivatives. From the torsion "
        "G = [[0 0 0] [0 0 0.27] [−0.23 −0.23 0]], the same as <i>kron.RotatingFrame</i>. Run: 120 V, 60 Hz, pushed to "
        "40 rad/s, 4 s.")
    r.table([
        ["model", "final speed", "ledger mismatch"],
        ["frame I (coils turn with the rotor)", "376.7928 rad/s", "1×10<super>−6</super>"],
        ["Kron's frame II(b) (<i>kron.RotatingFrame</i>)", "376.8235 rad/s", "5×10<super>−8</super>"],
        ["<b>built from Einstein's torsion</b>", "<b>376.8235 rad/s</b>", "5×10<super>−8</super>"],
        ["same inductances, torsion removed", "26.8128 rad/s — coasts on friction only", "4×10<super>−7</super>"],
    ], [0.45, 0.35, 0.2], bold_rows=(3,))
    r.p("The torsion-built motor and Kron's frame II(b) agree to 2×10<super>−12</super> rad/s; frame I agrees to "
        "0.03 rad/s (time step). Without the torsion there is no torque, and the rotor coasts as "
        "40·e<super>−4 s·friction/J</super> = 26.81 rad/s.")
    r.figure(os.path.join(HERE, "einstein_teleparallel.png"),
             "Figure 1. Left: the nonzero torsion components of Kron's motor at θ = 0.7 (axes: 1 = stator, 2s, 3s = "
             "stationary rotor coils, θ = rotor) in three notations — Einstein's −Λ in frame components (blue), the "
             "frame's anholonomy (yellow) and Hoffmann's Ω = L′·Λ in henry (red). Right: rotor speed of the four "
             "kronspice models; frame I, frame II(b) and the torsion-built motor overlap, while the motor without "
             "torsion only coasts down.")

    r.h1("Validation")
    r.p("<i>main_test.go</i> checks the torsion three ways, φ and the curvature, the linear plane-wave check of "
        "eqs. 5–9 and R<sub>αβ</sub>, and the motor runs. Every motor run closes its energy ledger to "
        "≤ 1×10<super>−6</super>.")

    r.h1("Where it stands")
    r.bullets([
        "<b>Einstein's theory was dropped.</b> He worked on distant parallelism from summer 1928 to spring 1931 and "
        "gave it up because no single set of field equations could be justified (Sauer 2006).",
        "<b>Its gravity part survives</b> as the teleparallel equivalent of general relativity. Electromagnetism from "
        "the torsion trace did not survive.",
        "<b>For Kron.</b> Kron's non-Riemannian machine dynamics (1934; Hoffmann 1949) uses the same objects — "
        "n-beins, anholonomy and torsion — as a working engineering method. Its torsion is trace-free (φ = 0) and "
        "sits on a curved space, so it is not electromagnetism in Einstein's sense.",
        "<b>Not decided here:</b> whether Kron drew on Einstein's 1928–31 papers is a history question; it would take "
        "Kron's 1934 references to settle.",
    ])

    r.h1("Reproducing the results")
    r.code("go run ./examples/einstein_teleparallel -viz examples/einstein_teleparallel/viz   # about 2 seconds\n"
           "python3 examples/einstein_teleparallel/viz/render.py\n"
           "go test ./examples/einstein_teleparallel\n"
           "python examples/einstein_teleparallel/viz/make_report.py "
           "examples/einstein_teleparallel/einstein_teleparallel_report.pdf")
    r.small("Sources: A. Einstein, Sitzungsber. Preuss. Akad. Wiss. 1928, 217–221 and 224–227 (A. Unzicker's "
            "translation); B. Hoffmann, Rev. Mod. Phys. 21, 535 (1949); G. Kron, J. Math. Phys. 13, 103 (1934); "
            "T. Sauer, Historia Mathematica 33, 399 (2006); J. W. Maluf, Ann. Phys. 525, 339 (2013).")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "einstein_teleparallel_report.pdf"))
