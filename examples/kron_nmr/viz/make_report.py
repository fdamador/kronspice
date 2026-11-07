# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Magnetic resonance as Kron's rotating frame, as a PDF report.

Run from the repository root:
    python examples/kron_nmr/viz/make_report.py examples/kron_nmr/kron_nmr_report.pdf

Numbers come from examples/kron_nmr/README.md (the example has no
run_output.txt). The only picture, viz/nmr.gif, is an animation whose first
frame shows an empty signal panel, so it is described but not drawn.
Layout from internal/pdfreport/kronreport.py.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401


def build(out):
    r = Report("examples/kron_nmr", "Magnetic resonance as Kron's rotating frame",
               "NMR and ESR solved as a three-winding Kron machine, and the bench measurement that discriminates the "
               "ring-electron hypothesis. CLAUDE.md §3b; register entry H21",
               "CLAUDE.md §3b, H21",
               "proton at 1 T, T1 = 5 µs, T2 = 2 µs, 400 points per Larmor period")

    r.h1("Summary")
    r.p("Bloch's equations for nuclear and electron magnetic resonance are a rotation of the magnetization about the "
        "field, and a rotation of currents at a fixed rate is Kron's frame II(b). This example solves NMR for the "
        "proton (and gives the ESR numbers for the electron) with the circuit engine's N machine element: three "
        "windings carry the magnetization components as currents, a skew-symmetric speed matrix carries the field, "
        "resistors carry the relaxation, and a fourth winding is the pickup coil. It is a companion to "
        "<i>kron_electron</i> (H21): the resonance frequency is the bench measurement that tells a ring electron from "
        "the measured one.")
    r.bullets([
        "<b>Free precession at 42.5766 MHz</b> against γB<sub>0</sub>/2π = 42.5775 MHz (−2×10<super>−5</super>, "
        "the trapezoidal phase error at 400 points per period).",
        "<b>Relaxation exact to the printed digits:</b> transverse decay 0.01832 at 4T2 (e<super>−4</super> = 0.01832); "
        "Mz recovery 0.79810 (1 − e<super>−4T2/T1</super> = 0.79810).",
        "<b>Rabi nutation and the π/2 pulse:</b> period 23.487 µs at B<sub>1</sub> = 1 mT, Mz − cos(γB<sub>1</sub>t) "
        "within 2.4×10<super>−4</super>, zero shaft power; a 5.87 µs pulse leaves M = (0, 1.0000, 0).",
        "<b>The continuous-wave line equals Bloch's steady state</b> to 5 digits at seven detunings.",
        "<b>Ledger closes:</b> source − heat − Δstored = −3×10<super>−14</super> J during free precession.",
        "<b>The discriminating measurement for H21:</b> a ring with the particle's charge and mass (g = 1) would "
        "resonate at 13.996 GHz (electron) and 7.623 MHz (proton) at 1 T, against the measured 28.025 GHz and "
        "42.578 MHz.",
    ])

    r.h1("The construction")
    r.p("Bloch's equations:")
    r.p("&nbsp;&nbsp;&nbsp;&nbsp;dM/dt = γ·M × B − (Mx, My, 0)/T2 − (0, 0, Mz − M0)/T1")
    r.p("where M is the magnetization, γ the gyromagnetic ratio, B the field, T1 and T2 the longitudinal and "
        "transverse relaxation times and M0 the equilibrium magnetization. Kron's frame II(b) (<i>Equivalent Circuits "
        "of Electric Machinery</i>, 1951, Epilogue eq. 9, p. 257) has constant inductances and a skew-symmetric speed "
        "matrix G, so that L·di/dt = −ω·G·i. The mapping:")
    r.bullets([
        "Three windings x, y, z carry Mx, My, Mz as their currents (L = 1 nH, so the winding voltages ωL·M are volts).",
        "G = L·γ·[B]<sub>×</sub>, the matrix of M × B, and the shaft turns at ω = 1. With B<sub>eff</sub> = "
        "(B<sub>1</sub>, 0, B<sub>0</sub> − ω<sub>rf</sub>/γ) the non-zero entries are G[x][y] = −LγB<sub>z</sub>, "
        "G[y][x] = LγB<sub>z</sub>, G[y][z] = −LγB<sub>1</sub>, G[z][y] = LγB<sub>1</sub>. A skew G has zero torque "
        "i<super>T</super>·G·i, so the precession takes no shaft power.",
        "The NMR rotating frame is Kron's rotating reference frame (Park's transform, CLAUDE.md §3b): at the RF "
        "frequency ω<sub>rf</sub> the RF field B<sub>1</sub> is a constant x-field and the static field becomes "
        "B<sub>0</sub> − ω<sub>rf</sub>/γ, so G is constant and resonance is where the effective z-field vanishes; "
        "ω<sub>rf</sub> = 0 is the laboratory frame.",
        "Relaxation is resistance on the windings: T2 = L/R<sub>2</sub> on x and y, T1 = L/R<sub>1</sub> on z with a "
        "series source R<sub>1</sub>·M<sub>0</sub> that restores equilibrium. The heat is on the ledger.",
        "A fourth, open winding (1 MΩ) coupled to x by M<sub>p</sub> = 0.1 nH is the pickup coil; its voltage "
        "M<sub>p</sub>·dMx/dt is the spectrometer's signal.",
    ])
    r.p("Proton run: γ<sub>p</sub> = 2.675 221 87×10<super>8</super> rad/s/T, B<sub>0</sub> = 1 T, M<sub>0</sub> = 1 "
        "in its own units, T1 = 5 µs and T2 = 2 µs (chosen short so the runs stay small), 400 time points per Larmor "
        "period. <i>BlochSteady</i> is Bloch's closed-form steady state at a detuning, the standard the "
        "continuous-wave run is read against.")

    r.h1("Results")
    r.h2("Proton at 1 T")
    r.table([
        ["check", "machine", "expected"],
        ["free precession", "42.5766 MHz", "γB<sub>0</sub>/2π = 42.5775 MHz (−2×10<super>−5</super>, the trapezoidal "
         "phase error (ωh)²/12 at 400 points per period; 136,252 points over 4T2)"],
        ["transverse decay at 4T2", "0.01832", "e<super>−4</super> = 0.01832"],
        ["Mz recovery at 4T2", "0.79810", "1 − e<super>−4T2/T1</super> = 0.79810"],
        ["Rabi nutation, B<sub>1</sub> = 1 mT, relaxation off",
         "period 23.487 µs; Mz − cos(γB<sub>1</sub>t) within 2.4×10<super>−4</super> over two periods",
         "2π/γB<sub>1</sub>; shaft power 0"],
        ["π/2 pulse (5.87 µs)", "M = (0, 1.0000, 0)", "onto the transverse plane"],
        ["free-induction decay", "42.5766 MHz at the pickup, peak 0.02676 V",
         "M<sub>p</sub>·γB<sub>0</sub>·M<sub>0</sub> = 0.02675 V, decaying with T2"],
        ["continuous-wave line, B<sub>1</sub> = 10 µT", "absorption, dispersion and Mz at seven detunings "
         "Δω·T2 = −3 … 3", "Bloch's steady state, to 5 digits"],
        ["ledger (free precession)", "source 8.02×10<super>−10</super> J, heat 9.83×10<super>−10</super> J, stored "
         "5.0×10<super>−10</super> → 3.19×10<super>−10</super> J",
         "source − heat − Δstored = −3×10<super>−14</super> J"],
    ], [0.25, 0.35, 0.4])

    r.h2("The continuous-wave line at B<sub>1</sub> = 10 µT")
    r.table([
        ["Δω·T2", "My (absorption)", "Bloch", "Mx (dispersion)", "Bloch", "Mz"],
        ["−3", "0.00054", "0.00054", "0.00161", "0.00161", "0.99999"],
        ["−1", "0.00268", "0.00268", "0.00268", "0.00268", "0.99996"],
        ["−0.5", "0.00428", "0.00428", "0.00214", "0.00214", "0.99994"],
        ["0", "0.00535", "0.00535", "0.00000", "0.00000", "0.99993"],
        ["0.5", "0.00428", "0.00428", "−0.00214", "−0.00214", "0.99994"],
        ["1", "0.00268", "0.00268", "−0.00268", "−0.00268", "0.99996"],
        ["3", "0.00054", "0.00054", "−0.00161", "−0.00161", "0.99999"],
    ], [0.12, 0.18, 0.16, 0.18, 0.16, 0.2])
    r.p("The absorption is the Lorentzian and the dispersion its odd partner, as Bloch's solution gives; at this "
        "B<sub>1</sub> the line is unsaturated (Mz stays at 0.9999).")

    r.h2("Resonance frequencies at 1 T: measured against the H21 ring")
    r.p("The H21 ring is a ring with the particle's charge and mass circulating together (g = 1, γ = e/2m).")
    r.table([
        ["", "measured", "ring (g = 1)"],
        ["proton NMR", "42.578 MHz (g<sub>p</sub> = 5.5857)", "7.623 MHz"],
        ["electron ESR", "28.025 GHz (g<sub>e</sub> = 2.0023)", "13.996 GHz"],
    ], [0.3, 0.4, 0.3])
    r.bullets([
        "A ring electron whose mass circulates with its charge would resonate at half the ESR frequency (its g ≤ 1 is "
        "derived in <i>kron_electron</i>), and a ring proton at a fifth of the NMR frequency.",
        "Standard prediction beside it: the Dirac electron with QED, g = 2.002 319 304 36, and the measured proton g.",
    ])

    r.h2("Animation")
    r.p("<i>viz/nmr.gif</i> (100 frames) shows the magnetization vector in the rotating frame (left) through the "
        "π/2 pulse and the decay that follows, with the pickup coil's signal drawn as it arrives (right). The frames "
        "are the solver's time points, the decay run converted to the rotating frame by Park's rotation. It is not "
        "reproduced here because a single frame does not show the motion; open the GIF to see it.")

    r.h1("Against existing measurements")
    r.p("The resonance frequencies are measured inputs, so the match is exact by construction: the proton's "
        "42.577 MHz/T is the CODATA value of γ<sub>p</sub>/2π, and the electron's 28.025 GHz/T follows from the "
        "measured g-factor. What the machine adds is the dynamics, and those are Bloch's 1946 equations, whose "
        "predictions (precession, T1 and T2 relaxation, the Lorentzian absorption and dispersion lines, saturation) "
        "were confirmed in the first NMR experiments of Purcell, Torrey and Pound and of Bloch, Hansen and Packard "
        "(both 1946) and in every spectrometer since; ESR was first observed by Zavoisky in 1945. The scaled "
        "relaxation times and pickup units are not measurements.")

    r.h1("Validation")
    r.table([
        ["test", "check", "held to"],
        ["TestFreePrecession", "frequency γB<sub>0</sub>; |Mxy| = e<super>−2</super> and Mz = 1 − e<super>−2T2/T1</super> "
         "at 2T2; ledger closed", "1e-4; mismatch 1e-3 of heat"],
        ["TestRabi", "Mz = cos(γB<sub>1</sub>t) at resonance; π/2 pulse gives (0, 1, 0)", "1e-3"],
        ["TestPickup", "pickup peak M<sub>p</sub>·γB<sub>0</sub>·M<sub>0</sub> at frequency γB<sub>0</sub>", "1e-3, 1e-4"],
        ["TestCWLine", "absorption, dispersion and Mz at Δω·T2 = −2, −0.5, 0, 0.5, 2 against BlochSteady", "2e-5"],
        ["TestExport", "free precession (.tran 50p 8u uic) and the CW run at resonance (.tran 40n 60u uic) exported "
         "with the N element", "—"],
    ], [0.22, 0.56, 0.22])

    r.h1("Limits")
    r.bullets([
        "T1 and T2 are chosen short (microseconds) so the runs stay small; they are not measured values.",
        "The resonance frequencies use measured g-factors as inputs; the model checks the dynamics, not the "
        "magnetic moments.",
        "Free precession carries a trapezoidal phase error of −2×10<super>−5</super> at 400 points per period.",
    ])

    r.h1("Status of H21 (context)")
    r.p("H21, the ring electron as a Kron machine, is <i>bounded</i> in the register: no ring with the electron's "
        "charge-to-mass ratio reaches g = 2. This example supplies the bench measurement that discriminates it: ESR at "
        "28.0 GHz/T measured against 14.0 GHz/T for a ring electron, and NMR at 42.6 MHz/T against 7.6 MHz/T for a "
        "ring proton.")

    r.h1("Reproducing the results")
    r.code("go run ./examples/kron_nmr -gif viz/nmr.gif\n"
           "go test ./examples/kron_nmr\n"
           "KRONSPICE_EXPORT=1 go test ./examples/kron_nmr -run TestExport   # netlists, KiCad, results sheets")
    r.small("Sources: F. Bloch, Phys. Rev. 70, 460 (1946); F. Bloch, W. W. Hansen, M. Packard, Phys. Rev. 69, 127 "
            "(1946); E. M. Purcell, H. C. Torrey, R. V. Pound, Phys. Rev. 69, 37 (1946); E. Zavoisky, J. Phys. USSR 9, "
            "211 (1945); G. Kron, Equivalent Circuits of Electric Machinery, 1951, Epilogue eq. 9; CODATA 2018 "
            "(Tiesinga et al., Rev. Mod. Phys. 93, 025010, 2021); D. Hanneke et al., PRL 100, 120801 (2008); X. Fan "
            "et al., PRL 130, 071801 (2023).")
    for f in r.s:
        if getattr(getattr(f, "style", None), "name", "") in ("h1", "h2"):
            f.keepWithNext = 1
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "kron_nmr_report.pdf"))
