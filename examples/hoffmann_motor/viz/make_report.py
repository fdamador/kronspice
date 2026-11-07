# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Report: Hoffmann's single-phase motor in Kron's non-Riemannian frames, as a PDF.

python examples/hoffmann_motor/viz/make_report.py examples/hoffmann_motor/hoffmann_motor_report.pdf

Numbers come from examples/hoffmann_motor/run_output.txt and README.md; the
figure is viz/hoffmann_motor.png. Layout from internal/pdfreport/kronreport.py.
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
    r = Report("examples/hoffmann_motor", "Hoffmann's single-phase motor in Kron's frames",
               "Kron's non-Riemannian electrodynamics (Hoffmann 1949): does the torsion of a non-holonomic frame "
               "carry the motor's torque? — CLAUDE.md §3a, M9c.1", "CLAUDE.md §3a, M9c.1",
               "120 V, 60 Hz, 5 s runs")

    r.h1("Summary")
    r.p("B. Hoffmann, “Kron's Non-Riemannian Electrodynamics,” <i>Rev. Mod. Phys.</i> 21, 535–540 (1949), explains "
        "for physicists what Kron did with rotating machines and works one example, a single-phase squirrel-cage "
        "induction motor. This example builds that motor and solves it three ways: with rotor coils that turn with the "
        "rotor (frame I), with Kron's stationary rotor coils and speed matrix (frame II(b)), and by integrating "
        "Hoffmann's eq. 5′ with the torsion term directly. New in kronspice for this: <i>kron.Torsion</i> (eq. 19), "
        "<i>kron.NewMetric</i> (eq. 15) and <i>kron.NonRiemannian</i> (eq. 5′), in <i>kron/hoffmann.go</i>.")
    r.bullets([
        "<b>Constant inductances, but torsion:</b> in the stationary-coil frame every inductance is constant, and the "
        "frame has torsion Ω(2s θ 3s) = −0.125, Ω(3s θ 1) = +0.115, Ω(3s θ 2s) = +0.125, independent of θ (spread "
        "10<super>−16</super>).",
        "<b>Eq. 5′ is exact:</b> at 200 random states it gives the frame-I forces to 5×10<super>−16</super>; without "
        "the torsion term it misses by up to 293 % of the largest force — Hoffmann's warning.",
        "<b>Three solutions agree:</b> pushed to +40 rad/s, the motor reaches 376.587 rad/s in frame I and by eq. 5′ "
        "(difference 4.1×10<super>−4</super> rad/s), 376.618 rad/s in frame II(b).",
        "<b>No start from rest; runs either way when pushed:</b> from rest the final speed is 0; pushed to −40 rad/s it "
        "runs to −376.587 rad/s.",
        "<b>Without the torsion it only coasts:</b> 40·e<super>−5 s·friction/J</super> = 24.261 rad/s.",
        "<b>Every ledger closes</b> to 2×10<super>−6</super> or better; every joule comes from the 120 V supply (or "
        "the push).",
    ])

    r.h1("The idea")
    r.p("A machine's currents and its rotor angle together are the coordinates of one moving point. The inductances "
        "and the inertia act as the geometry (the metric) of the space it moves in, and Lagrange's equations are the "
        "path of the point (Hoffmann eqs. 1–8).")
    r.bullets([
        "<b>Frame I: coils that turn with the rotor</b> (Fig. 2). The squirrel cage is replaced by two perpendicular "
        "coils fixed to the rotor. Their inductances to the stator change with the angle θ (eqs. 24, 42). The space is "
        "ordinary (Riemannian) and the equations are Lagrange–Maxwell with Christoffel symbols (eq. 5) "
        "(<i>kron.Christoffel</i>, <i>kron.Geodesic</i>).",
        "<b>Stationary coils</b> (Hoffmann's barred 2 and 3, written 2s and 3s here; Fig. 3, eq. 35). Current paths "
        "that stay still while the cage turns make every inductance a constant (eq. 29). The price is that this "
        "change of frame is <i>non-holonomic</i>: no set of coordinates goes with it (eq. 38). The space gains a "
        "torsion Ω (eq. 19), and the equations gain a term 2Ω·ẋ·ẋ (eq. 5′) that carries both the speed voltages "
        "and the torque.",
        "<b>Hoffmann's warning</b> (p. 540): putting the constant inductances into the ordinary Lagrange–Maxwell "
        "equations “would obtain wrong results”.",
    ])
    r.h2("The motor (all values declared)")
    r.table([
        ["part", "value", "part", "value"],
        ["stator coil L11", "0.25 H, 2 Ω", "inertia J", "0.001 kg·m²"],
        ["rotor coils L<sub>S</sub>", "0.25 H, 3 Ω each", "friction", "10<super>−4</super> N·m·s"],
        ["mutual M", "0.23 H", "supply", "120 V, 60 Hz"],
        ["salient pole L<sub>D</sub>", "0 (and 0.01, 0.02 H in some runs)", "poles; synchronous speed",
         "two; 377.0 rad/s"],
    ], [0.2, 0.3, 0.22, 0.28])

    r.h1("1. The two frames")
    r.p("In the stationary-coil frame the inductance matrix (henry; inertia in kg·m²) no longer depends on θ:")
    r.table([
        ["", "1", "2s", "3s", "θ"],
        ["1", "0.2500", "0.2300", "0.0000", "0.0000"],
        ["2s", "0.2300", "0.2500", "0.0000", "0.0000"],
        ["3s", "0.0000", "0.0000", "0.2500", "0.0000"],
        ["θ", "0.0000", "0.0000", "0.0000", "0.0010"],
    ], [0.2, 0.2, 0.2, 0.2, 0.2])
    r.p("The frame's torsion Ω<sub>λμν</sub> has non-zero components Ω(2s θ 3s) = −0.125, Ω(3s θ 1) = +0.115 and "
        "Ω(3s θ 2s) = +0.125, with Ω(μλν) = −Ω(λμν); it does not change with θ (spread 10<super>−16</super>). At 200 "
        "random states (angles, currents, their rates, speed, acceleration), eq. 5′ with the torsion gives the frame-I "
        "forces to 5×10<super>−16</super>. Ordinary Lagrange–Maxwell equations on the constant inductances miss by up "
        "to 293 % of the largest force.")

    r.h1("2. Starting (5 s runs)")
    r.table([
        ["run", "final speed", "energy from the mains", "ledger mismatch"],
        ["from rest, frame I", "<b>0.000 rad/s</b> — does not start", "1409 J = heat 1407 J + stored 2.2 J",
         "3×10<super>−7</super>"],
        ["pushed to +40 rad/s, frame I", "<b>376.587 rad/s</b>", "867 J = heat 796 J + stored 70.5 J",
         "2×10<super>−6</super>"],
        ["pushed to +40 rad/s, frame II(b)", "376.618 rad/s", "867 J = heat 796 J + stored 70.6 J",
         "1×10<super>−7</super>"],
        ["pushed to +40 rad/s, eq. 5′ (Runge–Kutta)", "<b>376.587 rad/s</b>", "—", "—"],
        ["pushed to −40 rad/s, frame I", "−376.587 rad/s", "867 J = heat 796 J + stored 70.5 J",
         "2×10<super>−6</super>"],
        ["pushed to +40 rad/s, no torsion", "24.261 rad/s — coasts on friction", "1409 J = heat 1407 J + stored 1.7 J",
         "3×10<super>−7</super>"],
    ], [0.3, 0.24, 0.32, 0.14])
    r.bullets([
        "<b>No start from rest.</b> A single-phase motor has no starting torque: its pulsating field is two equal "
        "fields turning in opposite directions. Real motors add a start winding or a capacitor.",
        "<b>Pushed, it runs either way</b>, up to just under synchronous speed (slip 0.107 % with no load).",
        "<b>The three solutions agree.</b> Frame I and eq. 5′ agree to 4.1×10<super>−4</super> rad/s. Eq. 5′ is "
        "integrated by this example's own Runge–Kutta code, which shares nothing with kronspice's machine element.",
        "<b>Hoffmann's warning.</b> Without the torsion there is no speed voltage and no torque; the pushed rotor "
        "follows friction alone.",
    ])

    r.h1("3. Torque–speed curve")
    r.p("With the shaft held at each speed, mean torque over the last 10 cycles (N·m). The curve is odd in speed "
        "(the negative half is the same with the sign changed).")
    r.table([
        ["rad/s", "0.0", "62.5", "125.0", "187.6", "250.1", "312.6", "375.1"],
        ["L<sub>D</sub> = 0", "0.0000", "0.0707", "0.1523", "0.2619", "0.4327", "0.7121", "0.0714"],
        ["L<sub>D</sub> = 0.02 H (salient pole)", "0.0000", "0.2309", "0.4845", "0.7854", "1.1438", "1.3829", "0.0845"],
    ], [0.24, 0.1, 0.1, 0.11, 0.11, 0.11, 0.11, 0.12])
    r.bullets([
        "<b>Shape.</b> Zero at standstill, odd in speed, largest near 91 % of synchronous speed. Frames I and II(b) "
        "agree to 1.7×10<super>−3</super> of the peak.",
        "<b>Salient pole.</b> It makes stationary coil 2s couple more tightly to the stator (stationary-frame L<sub>2s2s</sub> = "
        "L<sub>S</sub> − L<sub>D</sub> = M), which roughly doubles the torque at most speeds.",
        "<b>No-load slip</b> falls with L<sub>D</sub>: 376.587 rad/s (0.107 %) at 0 H, 376.694 rad/s (0.079 %) at "
        "0.01 H, 376.793 rad/s (0.053 %) at 0.02 H. The motor does not pull into synchronism at these values: a "
        "reluctance motor needs a much larger L<sub>D</sub>.",
    ])
    r.figure(os.path.join(HERE, "hoffmann_motor.png"),
             "Figure 1. Left: rotor speed in the start-up runs. From rest the motor stays at 0; pushed to +40 rad/s, "
             "frame I, frame II(b) and eq. 5′ overlap and run up to just under synchronous speed (dotted lines); pushed "
             "to −40 rad/s it runs backwards; without the torsion term it only coasts down. Right: mean torque against "
             "held shaft speed, frames I and II(b) overlapping, and with a salient pole L<sub>D</sub> = 0.02 H.")

    r.h1("Validation")
    r.p("<i>main_test.go</i> checks the random-state twin, the constancy of the torsion, start-up agreement, the "
        "coast-down and the torque curve; <i>kron/hoffmann_test.go</i> tests <i>Torsion</i>, <i>NewMetric</i> and "
        "<i>NonRiemannian</i>. Every kronspice run's energy ledger closes to 2×10<super>−6</super> or better. The "
        "motor adds no energy: every joule of shaft motion and heat comes from the 120 V supply (or the push).")
    r.p("<b>Limits.</b> All machine values are declared, not taken from a real motor; the rotor cage is represented "
        "by two coils, as in Hoffmann's paper; the motor runs with no load.")

    r.h1("Reproducing the results")
    r.code("go run ./examples/hoffmann_motor -viz examples/hoffmann_motor/viz   # about 10 seconds\n"
           "python3 examples/hoffmann_motor/viz/render.py\n"
           "go test ./examples/hoffmann_motor ./kron\n"
           "python examples/hoffmann_motor/viz/make_report.py examples/hoffmann_motor/hoffmann_motor_report.pdf")
    r.small("Sources: B. Hoffmann, Rev. Mod. Phys. 21, 535–540 (1949); G. Kron, “Non-Riemannian Dynamics of Rotating "
            "Electrical Machinery,” J. Math. Phys. 13, 103–194 (1934); G. Kron, Equivalent Circuits of Electric "
            "Machinery (Wiley, 1951), Epilogue — frames I and II(b) as used by kron.RotatingFrame.")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "hoffmann_motor_report.pdf"))
