# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Kron's moving inductors (milestone M9c.1) as a PDF report.

Run from the repository root:
    python examples/kron_machine/viz/make_report.py examples/kron_machine/kron_machine_report.pdf

Numbers come from examples/kron_machine/README.md and the results sheets in
examples/kron_machine/export/ (the example has no run_output.txt and no
figures). Layout from internal/pdfreport/kronreport.py.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401


def build(out):
    r = Report("examples/kron_machine", "Kron's moving inductors",
               "An induction motor started in Kron's frames I and II, and a parametric inductance pump against "
               "Mathieu's equation. CLAUDE.md §3a, §3l item 3; milestone M9c.1",
               "CLAUDE.md §3a, M9c.1",
               "transient 1.5 s, 0.5 ms step, RelTol 1e-5")

    r.h1("Summary")
    r.p("This is the milestone M9c.1 example for the N machine element (<i>device/machine.go</i>). One two-phase "
        "induction motor is started from rest twice: in Kron's holonomic frame I, where the rotor windings sit on "
        "the rotor and the inductances are trigonometric in the rotor angle θ, and in his non-holonomic rotor frame "
        "II(b), where the inductances are constant and speed voltages G·i·ω appear. The two runs must agree. Then the "
        "same element, with its coordinate driven at constant speed, is a time-varying inductance pumping an LC tank, "
        "and its growth is checked against Mathieu's equation.")
    r.bullets([
        "<b>The two frames agree.</b> Final speeds 313.568 rad/s (frame I) and 313.696 rad/s (frame II(b)), against "
        "synchronous 314.16 rad/s; over the whole start-up, speed within 1.1×10<super>−4</super> and stator current "
        "within 5.8×10<super>−4</super> of peak.",
        "<b>Both ledgers close:</b> supply 326.10 J = heat 227.29 J + stored 98.75 J (frame I, mismatch "
        "9.8×10<super>−5</super>); 326.21 J = 227.38 J + 98.83 J (frame II(b), mismatch −1.4×10<super>−6</super>).",
        "<b>Speed terms appear only in the rotor rows</b> of the transformed tensor G′, as Kron's frame II requires; "
        "a frame that leaves θ in L is refused.",
        "<b>The pumped tank grows at Mathieu's rate:</b> 78.550 /s against m·ω<sub>0</sub>/4 = 78.540 /s undamped, "
        "39.273 /s against 39.270 /s at half damping, and decays at −78.571 /s against −78.540 /s when damping is "
        "twice the pump rate.",
        "<b>The pump pays for everything the tank gains:</b> pump work = tank gain + heat to 6×10<super>−4</super>. "
        "A coefficient of performance computed from the capacitor's source alone, as in H16, omits that work.",
    ])

    r.h1("The construction")
    r.h2("The primitive machine")
    r.p("Kron's two-phase primitive machine has stator windings α, β and rotor windings a, b at angle θ, "
        "self-inductances Ls = Lr = 0.12 H and stator–rotor mutuals M·R(θ) with M = 0.11 H, so the inductance tensor "
        "is [[Ls·I, M·R(θ)], [M·R(θ)<super>T</super>, Lr·I]], with R the 2×2 rotation. In kronspice this is a "
        "<i>circuit.MachineSpec</i> with a first harmonic in θ.")
    r.h2("Frame I, holonomic")
    r.p("The element is solved with L(θ) as written: Kron's eq. 3, e<sub>α</sub> = R<sub>αβ</sub>·i<super>β</super> "
        "+ L<sub>αβ</sub>·di<super>β</super>/dt + [βγ, α]·i<super>β</super>·i<super>γ</super>, with the Christoffel "
        "symbols [βγ, α] of eq. 4 built from the metric (L(θ), J), and the torque ½·i<super>T</super>·(∂L/∂θ)·i "
        "(<i>Equivalent Circuits of Electric Machinery</i>, 1951, Epilogue pp. 254–255). <i>kron.Christoffel</i> and "
        "<i>kron.Geodesic</i> reproduce the element's voltages and torque from those symbols to 10<super>−12</super>.")
    r.h2("Frame II(b), the rotor frame")
    r.p("<i>kron.RotatingFrame</i> applies the connection tensor C(θ) of <i>kron.RotorFrame</i> to the rotor axes "
        "(power invariant, e′ = C<super>T</super>·e). L′ = C<super>T</super>·L·C becomes constant and the rotor rows "
        "gain the speed matrix G = C<super>T</super>·L·∂C/∂θ, Kron's non-holonomic object, so that "
        "e = R·i + p(L·i) + G·i·pθ (p = d/dt) with torque i<super>T</super>·G·i (Epilogue eq. 9, p. 257). The program "
        "prints the transformed tensors:")
    r.code("L' = [[0.12 0 0.11 0], [0 0.12 0 0.11], [0.11 0 0.12 0], [0 0.11 0 0.12]]\n"
           "G' = [[0 0 0 0], [0 0 0 0], [0 0.11 0 0.12], [-0.11 0 -0.12 0]]")
    r.h2("The motor circuit")
    r.table([
        ["part", "value"],
        ["supply", "balanced two-phase 100 V, 50 Hz (sin and cos sources), 0.5 Ω per phase"],
        ["rotor windings", "shorted through 0.8 Ω"],
        ["rotor inertia", "0.002 kg·m², as a capacitor on the mechanical node w"],
        ["friction load", "0.2 mN·m per rad/s, as a 5 kΩ resistor on node w (angular velocity across, torque through)"],
        ["synchronous speed", "314.16 rad/s"],
        ["transient", "1.5 s, 0.5 ms step, RelTol 1e-5"],
    ], [0.25, 0.75])
    r.h2("The parametric pump")
    r.p("One winding with L(t) = L0·(1 + m·cos 2ω<sub>0</sub>t), L0 = 1 mH, m = 0.05, made from a second harmonic in θ "
        "with the shaft turned at ω<sub>0</sub> by a voltage source on the mechanical node, so θ = ω<sub>0</sub>t and "
        "the inductance varies at 2ω<sub>0</sub>. It is in series with C = 1/(ω<sub>0</sub>²L0) (1 V initial charge, "
        "f0 = 1 kHz) and a resistor R = damp·(m·ω<sub>0</sub>/4)·2L0, where damp is the damping R/(2L0) as a multiple "
        "of the pump rate. The tank runs 60 periods at 100 steps per period; the growth rate is the slope of "
        "log(peak) over periods 30 to 59.")

    r.h1("Results")
    r.h2("Start-up in the two frames")
    r.table([
        ["frame", "final speed", "points", "supply", "heat", "stored (magnetic + kinetic)", "ledger mismatch"],
        ["I, holonomic, L(θ)", "313.568 rad/s", "7953", "326.10 J", "227.29 J", "98.75 J", "9.8×10<super>−5</super>"],
        ["II(b), rotor frame, G·i·ω", "313.696 rad/s", "8439", "326.21 J", "227.38 J", "98.83 J", "−1.4×10<super>−6</super>"],
    ], [0.22, 0.13, 0.08, 0.11, 0.11, 0.2, 0.15])
    r.bullets([
        "The test holds the two final speeds within 10<super>−3</super> of each other and the motor between "
        "300 rad/s and synchronous (<i>analysis.TestMachineFrameTwin</i>).",
        "The ledger mismatch is the trapezoidal integration error over the sources' breakpoints.",
    ])
    r.h2("The pumped tank against Mathieu")
    r.p("ω<sub>0</sub> = 2π·1 kHz, m = 0.05, pump rate m·ω<sub>0</sub>/4 = 78.54 /s.")
    r.table([
        ["damping R/2L0 (× pump rate)", "growth, run", "Mathieu m·ω<sub>0</sub>/4 − R/2L0", "pump work", "heat"],
        ["0", "78.550 /s", "78.540 /s", "0.07876 J", "5.1×10<super>−10</super> J"],
        ["0.5", "39.273 /s", "39.270 /s", "1.396 mJ", "0.707 mJ"],
        ["2", "−78.571 /s", "−78.540 /s", "4.38 µJ", "17.1 µJ"],
    ], [0.26, 0.17, 0.25, 0.16, 0.16])
    r.bullets([
        "The test holds the growth within 1 % of Mathieu's rate for the first two rows "
        "(<i>analysis.TestMachineParametric</i>).",
        "The energy the tank gains enters through the pump's mechanical port: pump work = tank gain + heat, mismatch "
        "≤ 6×10<super>−4</super>.",
        "This is the standard account of parametric excitation (Mandelstam and Papalexi 1932; Kouřil and Vrba "
        "§10.2.2, whose averaged Mathieu exponent σ = −α + √((mω<sub>0</sub>/4)² − Δω²) reduces to the column above "
        "at zero detuning). It is the standard counterpart of the parametric circuits of H16 "
        "(<i>parametric_lcr</i>): the agent that varies L pays for everything the tank gains.",
    ])
    r.h2("Exports")
    r.p("<i>export/</i> holds the frame-I and frame-II motors (.tran 0.5m 1.5 uic) and the damping-0.5 parametric "
        "tank (.tran 10u 60m uic) as netlists with the N element and .model … machine, KiCad schematics and results "
        "sheets. The sheets show the DC operating point by MNA only: the machine is nonlinear, so Kron's linear torn "
        "solution is not drawn, and the torn Newton solve is checked in <i>kron</i>'s own tests.")
    r.table([
        ["netlist", "elements, nodes", "DC operating point (MNA)"],
        ["kron_machine_motor_frame1.cir", "9, 7", "sa = 100 V; ia, ib, ra, rb, sb, w = 0"],
        ["kron_machine_motor_frame2.cir", "9, 7", "sa = 100 V; ia, ib, ra, rb, sb, w = 0"],
        ["kron_machine_parametric.cir", "4, 3", "w = 6283.19 (shaft speed ω<sub>0</sub>); a, b = 0"],
    ], [0.36, 0.18, 0.46])

    r.h1("Validation")
    r.bullets([
        "Frame twin: speeds agree to 10<super>−3</super>, near synchronous; speed within 1.1×10<super>−4</super> and "
        "stator current within 5.8×10<super>−4</super> of peak over the start-up (doc/M9c1.md).",
        "Christoffel form (Kron 1951 eqs. 3–4) equals the element to 10<super>−12</super> (<i>kron</i> tests).",
        "Mathieu growth within 1 %; pump-work ledger mismatch ≤ 6×10<super>−4</super>.",
    ])

    r.h1("Limits")
    r.bullets([
        "The machine has a first-harmonic winding distribution only (the motor) or a second harmonic only (the pump).",
        "The growth rate is fitted from peaks over periods 30 to 59 at 100 steps per period; the fit and the "
        "trapezoidal integrator set the small differences from Mathieu's rate.",
        "The example has no run_output.txt or figures; numbers are those recorded in its README and export sheets.",
    ])

    r.h1("Reproducing the results")
    r.code("go run ./examples/kron_machine\n"
           "go test ./examples/kron_machine\n"
           "KRONSPICE_EXPORT=1 go test ./examples/kron_machine -run TestExport   # netlists, KiCad, results sheets")
    r.small("Sources: G. Kron, Equivalent Circuits of Electric Machinery, Wiley, 1951, Epilogue pp. 253–265; G. Kron, "
            "J. Math. Phys. 13, 103–194 (1934); F. Kouřil and K. Vrba, Non-linear and Parametric Circuits, 1988, "
            "§10.2.2; L. Mandelstam and N. Papalexi, Z. Phys. 73, 223 (1932); doc/M9c1.md.")
    for f in r.s:
        if getattr(getattr(f, "style", None), "name", "") in ("h1", "h2"):
            f.keepWithNext = 1
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "kron_machine_report.pdf"))
