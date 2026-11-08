# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""The fdamador neutrino as a circuit: a series RLC loop (H28).

    python3 examples/fdamador_neutrino/viz/plot_circuit.py [OUT.png]

Writes viz/circuit.png: the bench ring of the README (§1; a = 10 mm,
r = 2.5 mm, εr = 100) as the lumped loop kronspice solves (main.go, ERing),
with the values printed in run_output.txt:
  L     = 2.156e-08 H   ring self-inductance μ0·a(ln(8a/r) − 7/4)   (B through the hole)
  C     = 2.767e-13 F   the dielectric ring as series capacitance  (E round the ring)
  R_rad = 6.862 Ω       magnetic-dipole radiation resistance        (energy out, magnetic channel)
  R_d   = 0.02791 Ω     dielectric loss                              (heat, thermal channel)
f0 = 1/(2π√(LC)) and Q = ω0L/(R_rad + R_d) are recomputed here from them.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
L, C, RRAD, RD = 2.156e-08, 2.767e-13, 6.862, 0.02791
E_COL, B_COL, R_COL, K = "#d9480f", "#1c64b4", "#2b8a3e", "#222222"


def coil(ax, x0, x1, y, n=5, col=B_COL):
    ax.plot([x0, x0 + 0.15], [y, y], color=K, lw=2)
    ax.plot([x1 - 0.15, x1], [y, y], color=K, lw=2)
    t = np.linspace(0, np.pi * n, 400)
    xs = np.linspace(x0 + 0.15, x1 - 0.15, 400)
    ax.plot(xs, y + 0.22 * np.abs(np.sin(t)), color=col, lw=2.4)


def zigzag(ax, x, y0, y1, col, n=6, horiz=False):
    m = 0.15
    pts = np.linspace(y0 + m, y1 - m, 2 * n + 1)
    off = np.array([0] + [0.16 * (-1) ** k for k in range(2 * n - 1)] + [0])
    if horiz:
        ax.plot([y0, y0 + m], [x, x], color=K, lw=2); ax.plot([y1 - m, y1], [x, x], color=K, lw=2)
        ax.plot(pts, x + off, color=col, lw=2.2)
    else:
        ax.plot([x, x], [y0, y0 + m], color=K, lw=2); ax.plot([x, x], [y1 - m, y1], color=K, lw=2)
        ax.plot(x + off, pts, color=col, lw=2.2)


def cap(ax, x, y0, y1, col=E_COL):
    ym = (y0 + y1) / 2
    ax.plot([x, x], [y0, ym - 0.09], color=K, lw=2)
    ax.plot([x, x], [ym + 0.09, y1], color=K, lw=2)
    for yy in (ym - 0.09, ym + 0.09):
        ax.plot([x - 0.32, x + 0.32], [yy, yy], color=col, lw=3.2)


def main(out):
    w0 = 1 / np.sqrt(L * C)
    f0 = w0 / (2 * np.pi)
    Q = w0 * L / (RRAD + RD)
    Qrad = w0 * L / RRAD

    fig = plt.figure(figsize=(13, 5.4))
    fig.suptitle("The fdamador neutrino as a circuit: a series RLC loop (bench ring, a = 10 mm, εr = 100)",
                 fontsize=13)
    ax = fig.add_axes([0.0, 0.06, 0.5, 0.82])
    ax.set_xlim(-0.6, 6.4); ax.set_ylim(-0.9, 3.6); ax.set_aspect("equal"); ax.axis("off")
    x0, x1, y0, y1 = 0.5, 4.5, 0.4, 2.8
    coil(ax, x0, x1, y1)                           # top: L
    cap(ax, x1, y0, y1)                            # right: C
    zigzag(ax, y0, x0, x1, R_COL, horiz=True)      # bottom: R_rad
    zigzag(ax, x0, y0, y1, "#888888", n=4)         # left: R_d
    for xx, yy in ((x0, y0), (x1, y0), (x0, y1), (x1, y1)):
        ax.plot(xx, yy, "o", color=K, ms=4)
    ax.text(2.5, y1 + 0.42, f"L = {L * 1e9:.2f} nH", ha="center", fontsize=12, color=B_COL, weight="bold")
    ax.text(2.5, y1 + 0.78, "B through the hole: ½LI²", ha="center", fontsize=9.5, color=B_COL)
    ax.text(x1 + 0.45, (y0 + y1) / 2 + 0.15, f"C = {C * 1e12:.4f} pF", fontsize=12, color=E_COL, weight="bold")
    ax.text(x1 + 0.45, (y0 + y1) / 2 - 0.2, "E round the ring:\n½CV²", fontsize=9.5, color=E_COL, va="top")
    ax.text(2.5, y0 - 0.45, f"R_rad = {RRAD:.3f} Ω", ha="center", fontsize=12, color=R_COL, weight="bold")
    ax.text(2.5, y0 - 0.78, "radiation out (magnetic dipole), magnetic channel", ha="center", fontsize=9, color=R_COL)
    ax.text(x0 - 0.3, (y0 + y1) / 2 + 0.15, f"R_d =\n{RD * 1e3:.1f} mΩ", ha="right", fontsize=11, color="#666")
    ax.text(x0 - 0.3, (y0 + y1) / 2 - 0.45, "dielectric\nheat", ha="right", fontsize=9, color="#666", va="top")
    th = np.linspace(0.3, 2 * np.pi - 0.3, 100)
    ax.plot(2.5 + 0.55 * np.cos(th), 1.6 + 0.55 * np.sin(th), color=K, lw=1.2)
    ax.annotate("", xy=(2.5 + 0.55 * np.cos(th[-1]), 1.6 + 0.55 * np.sin(th[-1])),
                xytext=(2.5 + 0.55 * np.cos(th[-3]), 1.6 + 0.55 * np.sin(th[-3])),
                arrowprops=dict(arrowstyle="-|>", color=K, lw=1.2))
    ax.text(2.5, 1.6, "I", ha="center", va="center", fontsize=12, style="italic")

    tx = fig.add_axes([0.53, 0.06, 0.46, 0.82]); tx.axis("off")
    rows = [
        ("Resonance", f"f0 = 1/(2π√(LC)) = {f0 / 1e9:.4f} GHz"),
        ("Quality", f"Q = ω0L/(R_rad + R_d) = {Q:.2f}   (radiation alone {Qrad:.2f})"),
        ("Period, decay", f"T = {1e9 / f0:.3f} ns; energy time Q/ω0 = {Q / w0 * 1e9:.2f} ns"),
        ("Energy per cycle", f"E and B trade it each quarter period; {100 * (1 - np.exp(-2 * np.pi / Q)):.0f} % "
                             "leaves per cycle"),
        ("Where from", "L: Maxwell, Treatise Art. 706, μ0a(ln(8a/r) − 7/4)\n"
                       "R_rad: (η0π/6)(ka)⁴ (Balanis 5-24), ka = 0.432\n"
                       "C: set by the ring's εr = 100 so that f0 is the TE01δ mode"),
        ("Checked", "kronspice AC peak and ring-down agree with these values;\n"
                    "energy ledger closed (README §1)"),
        ("At neutrino size", "to last as long as the SN1987A neutrinos (proper\n"
                             "time 5.3×10⁴ s, Q = 8×10¹⁸) the loop needs an\n"
                             "effective vacuum permittivity of 3.4×10¹³\n"
                             "at a = 1.5×10⁻¹² m (Part 2)"),
        ("Not in the circuit", "mass and the weak interaction, by which every\nneutrino has been detected"),
    ]
    y = 0.98
    for k, v in rows:
        tx.text(0.0, y, k, fontsize=10, weight="bold", va="top", color="#1f3a5f")
        tx.text(0.27, y, v, fontsize=9, va="top")
        y -= 0.065 + 0.055 * v.count("\n") + 0.02
    fig.savefig(out, dpi=140)
    print(f"wrote {out}  (f0 {f0:.6e} Hz, Q {Q:.3f}, Qrad {Qrad:.3f})")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "circuit.png"))
