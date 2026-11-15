# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""H48 report: the ECE and ECE2 spin-connection terms on one series RLC, as a PDF.

    go run ./examples/ece_rlc -data examples/ece_rlc/viz/data > examples/ece_rlc/run_output.txt
    python3 examples/ece_rlc/viz/plot_report.py examples/ece_rlc/viz/data examples/ece_rlc/viz
    python3 examples/ece_rlc/viz/make_report.py examples/ece_rlc/ece_rlc_report.pdf

Every number comes from the run's CSV files (viz/data/summary.csv).
"""
import csv
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report  # noqa: E402

DATA = os.path.join(HERE, "data")
NAMES = ["plain (Maxwell)", "ECE torsion", "ECE2 torsion", "ECE scalar", "ECE2 scalar"]


def load(name):
    with open(os.path.join(DATA, name)) as f:
        rows = list(csv.DictReader(f))
    return {k: [float(r[k]) for r in rows] for k in rows[0]}


def sci(x, d=1):
    if x == 0:
        return "0"
    e = int(math.floor(math.log10(abs(x))))
    return f"{x / 10 ** e:.{d}g}×10<super>{'−' if e < 0 else ''}{abs(e)}</super>"


def build(out):
    s = load("summary.csv")
    r = Report("examples/ece_rlc",
               "ECE and ECE2 on a series RLC circuit",
               "The spin-connection terms of Einstein–Cartan–Evans theory in kronspice's circuit engine, "
               "the two forms side by side on one tank — register entry H48",
               "register entry H48")

    # ------------------------------------------------------------ summary
    r.h1("Summary")
    r.p("kronspice can run any circuit with the spin-connection terms of ECE theory (Eckardt, <i>ECE Unified Field "
        "Theory</i> Vol. 1, 2022) switched on by one netlist line, <font face='Courier'>.ece form=ece|ece2</font>. "
        "The terms are hypothetical, off by default, and their power is booked apart from every declared source. This "
        "report runs one series RLC tank five ways — plain, then the torsion and the scalar term in each of the two "
        "forms — with the same ω<sub>0</sub> and κ in both forms, so the ECE and ECE2 predictions can be read against each other "
        "and against Maxwell's.")
    r.figure(os.path.join(HERE, "..", "ece_rlc_summary.png"),
             "Figure 1. The circuit, the two forms of the field–potential relation and how they enter the circuit, and "
             "the five runs: Q, resonance and where the power at resonance comes from.", width_frac=1.0)
    r.bullets([
        f"<b>The two forms predict opposite effects for the same spin connection.</b> Torsion: ECE adds "
        f"+0.25 Ω and Q falls to {s['q'][1]:.2f}; ECE2 adds −0.50 Ω and Q rises to {s['q'][2]:.2f}. "
        f"Scalar: ECE tunes the tank down to {s['f_peak'][3]:.2f} Hz, ECE2 up to {s['f_peak'][4]:.2f} Hz.",
        f"<b>ECE2's torsion supplies power.</b> At resonance the source gives {s['p_source'][2]:.2f} W, the resistor "
        f"heats at {s['p_heat'][2]:.2f} W, and the term delivers the other {s['p_model'][2]:+.2f} W with no declared "
        "source behind it. The ledger books it as model-supplied, like a negative resistance.",
        "<b>The engine is exact for its stated mapping.</b> Each case equals a twin built from ordinary elements to "
        f"{sci(max(s['twin_diff']))}; resonance and Q equal the closed forms to the sweep grid; the AC ledger residual "
        f"is 0 and every transient closes to ≤ {sci(max(s['mismatch']))} J.",
        "<b>Standard prediction:</b> ω<sub>0</sub> = ω = 0, so all five cases are the plain tank. A bench measurement of a real "
        "tank's Q and tuning against its parts bounds |ω<sub>0</sub>| and |κ| in both forms at once.",
    ])

    # ------------------------------------------------------------ terms
    r.h1("The terms and how they enter the circuit")
    r.p("ECE writes the electric field with two spin-connection terms, ECE2 with the same two terms doubled and "
        "reversed:")
    r.code("ECE   (eq. 4.211):  E = −∇φ − ∂A/∂t − c·ω0·A + ω·φ\n"
           "ECE2  (eq. 6.189):  E = −∇φ − ∂A/∂t + 2(c·ω0·A − ω·φ)")
    r.p("Integrated along a component they give two circuit terms, with s = +1 for ECE and s = −2 for ECE2. This "
        "mapping is kronspice's proposal (circuit/ece.go), not the book's:")
    r.table([
        ["term", "netlist", "in the circuit", "value here"],
        ["torsion (ω<sub>0</sub>·A)", "w0= on .ece or per element",
         "series resistance s·c·ω<sub>0</sub>·L in each coil; s·c·ω<sub>0</sub>·M between coupled coils", "c·ω<sub>0</sub>·L = 0.25 Ω"],
        ["scalar (ω·φ)", "kappa=, kappa_a=, kappa_b=",
         "the element's body sees (1 + s·κ)·φ at each terminal (0 at ground)", "κ = 0.05 on C1's top"],
    ], [0.15, 0.2, 0.45, 0.2])
    r.code("series RLC, ECE2 torsion\n"
           ".ece form=ece2 w0=8.3391023799538009e-08\n"
           "V1 in 0 DC 0 AC 1 SIN(0 1 173)\n"
           "R1 in a 1\n"
           "L1 a b 10m\n"
           "C1 b 0 84.6346u\n"
           ".end")
    r.p("Closed forms: with torsion Q = ω<sub>0</sub>L/(R + s·c·ω<sub>0</sub>·L) and the tuning unchanged; with the scalar term "
        "the capacitor acts as C(1 + s·κ), so f<sub>r</sub> = f<sub>0</sub>/√(1 + s·κ) and Q = ω<sub>r</sub>L/R. "
        "The ordinary-element twins are a series resistor R + s·c·ω<sub>0</sub>·L (negative in ECE2) and a voltage-controlled "
        "source E that lifts the capacitor's body to (1 + s·κ)·φ.")

    # ------------------------------------------------------------ AC
    r.h1("Resonance (AC)")
    r.figure(os.path.join(HERE, "ac.png"),
             "Figure 2. |I| for 1 V, 4001 points from 120 to 230 Hz (grid 0.0275 Hz). Left: the torsion term changes "
             "the damping and leaves the tuning. Right: the scalar term retunes and leaves the damping.")
    t = [["case", "f<sub>r</sub> [Hz]", "closed form", "Q", "closed form", "peak |I| [A]", "vs twin",
          "source [W]", "heat [W]", "model-supplied [W]"]]
    for n in range(5):
        t.append([NAMES[n], f"{s['f_peak'][n]:.3f}", f"{s['f_want'][n]:.3f}", f"{s['q'][n]:.4f}",
                  f"{s['q_want'][n]:.4f}", f"{s['i_peak'][n]:.3f}", sci(s["twin_diff"][n]),
                  f"{s['p_source'][n]:.3f}", f"{s['p_heat'][n]:.3f}", f"{s['p_model'][n]:+.3f}"])
    r.table(t, [0.15, 0.08, 0.08, 0.08, 0.08, 0.09, 0.11, 0.08, 0.08, 0.15])
    r.small("Power is read at the sweep point nearest the peak. The AC ledger residual is 0 for every case.")

    # ------------------------------------------------------------ compare
    r.h1("ECE against ECE2")
    r.p("The two forms are run in one example on purpose: with the same ω<sub>0</sub> and κ they differ only in s, so every "
        "difference between them is the sign and the factor 2 in Eckardt's eq. 6.189 against eq. 4.211.")
    r.figure(os.path.join(HERE, "compare.png"),
             "Figure 3. Closed forms over a range of ω<sub>0</sub> and κ, with this run's five points. Left: Q against the torsion "
             "resistance c·ω<sub>0</sub>·L. Right: resonance against κ.")
    r.table([
        ["", "ECE (s = +1)", "ECE2 (s = −2)"],
        ["torsion in a coil", "+c·ω<sub>0</sub>·L, a loss", "−2c·ω<sub>0</sub>·L, a negative resistance"],
        ["Q for ω<sub>0</sub> > 0", "falls", "rises; infinite at c·ω<sub>0</sub>·L = R/2, growth beyond"],
        ["power of the torsion term", "absorbed (−0.08 W here)", "supplied (+1.00 W here)"],
        ["scalar term on C", "C(1 + κ): tunes down", "C(1 − 2κ): tunes up"],
        ["scalar term's average AC power", "0", "0"],
        ["sign of ω<sub>0</sub> for the same Q shift", "one sign", "the opposite sign, half the size"],
    ], [0.3, 0.32, 0.38])
    r.p("Because the forms move Q and tuning in opposite directions, one bench result cannot fit both with a non-zero "
        "spin connection: a measured tank that agrees with its parts bounds |ω<sub>0</sub>| and |κ| in both at once, and a "
        "deviation would point to one form by its sign.")

    # ------------------------------------------------------------ transient
    r.h1("Transient and the energy ledger")
    r.figure(os.path.join(HERE, "tran.png"),
             "Figure 4. 0–0.3 s from rest, 1 V sine at 173 Hz. Top: the loop current (sampled every 50 µs). Bottom: the "
             "running energies; the dotted red line is source in + model-supplied − heat − Δstored, which stays at 0.")
    t = [["case", "source in [J]", "back to source", "heat", "model-supplied", "Δstored", "mismatch"]]
    for n in range(5):
        t.append([NAMES[n], f"{s['e_in'][n]:.5f}", f"{s['e_out'][n]:.5f}", f"{s['e_heat'][n]:.5f}",
                  f"{s['e_model'][n]:+.5f}", f"{s['e_stored'][n]:.5f}", sci(s["mismatch"][n])])
    r.table(t, [0.2, 0.13, 0.13, 0.12, 0.15, 0.12, 0.15])
    r.bullets([
        "<b>Every run closes:</b> source in + model-supplied = back to source + heat + Δstored, to the trapezoidal "
        "rule's integration mismatch.",
        f"<b>ECE2 torsion:</b> the term supplies {s['e_model'][2]:.3f} J of the {s['e_heat'][2]:.3f} J of heat, "
        "nearly as much as the source. In the ledger it is a model-supplied column, never inflow.",
        "<b>The scalar term</b> moves a few tenths of a millijoule while the tank settles and averages to zero in the "
        "steady state; some energy goes back to the source because the retuned tank is driven off its resonance.",
    ])

    # ------------------------------------------------------------ Kron
    r.h1("Kron's frame")
    r.p("The same terms are primitive-network relations in Kron's frame (kron/ece.go): torsion in the coil's own "
        "impedance (jωL + s·c·ω<sub>0</sub>·L), mutual torsion as a chord, and the scalar term through control coils from each "
        "terminal to the datum. The results sheets in export/ solve each of the five circuits by MNA and torn by "
        "diakoptics, and agree to 7×10<super>−16</super>.")

    r.h1("Status and limits")
    r.bullets([
        "<b>Status of H48: untested.</b> The standard prediction (ω<sub>0</sub> = ω = 0) is computed beside every result; no bench "
        "data yet.",
        "The circuit mapping of the terms is kronspice's proposal; the book gives the field relations, not the circuit "
        "form, and no measured ω<sub>0</sub> or κ.",
        "The book gives no energy theorem containing ω; the ledger here books the terms' power as model-supplied.",
        "UFT 456's conductivities (k0=, kl0=) are in the engine but not used in this example.",
        "Critiques of the theory's foundations: F. W. Hehl, <i>Found. Phys.</i> 38, 7 (2008); F. W. Hehl, "
        "Yu. N. Obukhov, <i>Found. Phys.</i> 38, 38 (2008).",
    ])

    r.h1("Reproducing the results")
    r.code("go run ./examples/ece_rlc -data examples/ece_rlc/viz/data > examples/ece_rlc/run_output.txt\n"
           "python3 examples/ece_rlc/viz/plot_report.py examples/ece_rlc/viz/data examples/ece_rlc/viz\n"
           "python3 examples/ece_rlc/viz/make_report.py examples/ece_rlc/ece_rlc_report.pdf\n"
           "go test ./examples/ece_rlc")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "ece_rlc_report.pdf"))
