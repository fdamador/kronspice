# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""The photon electron as the circuit mirror of the fdamador neutrino.

    python3 examples/photon_electron/viz/plot_circuit.py [OUT.png]

Writes viz/circuit.png. The fdamador neutrino (H28) is this example's field
turned inside out: E round the ring and B through the hole, a series RLC
loop (examples/fdamador_neutrino/viz/circuit.png, values from its
run_output.txt). Exchanging E and B is circuit duality: voltage ↔ current,
L ↔ C, series ↔ parallel, R ↔ R0²/R. With the scale R0 = √(L/C) of the
neutrino loop the dual keeps L, C, f0 and Q, and gives the electron's
mirror: a parallel RLC tank, B round the ring (toroidal flux, the
inductor), E through the hole (the capacitor), radiation as a shunt
resistor R0²/R_rad.

What the mirror is and is not:
  - it is exact as circuit algebra (checked here: f0 and Q identical);
  - it is a mirror of the models, not of the measured electron or neutrino;
  - the real dielectric torus is not dual (εr = 100, μr = 1): its B-ring
    mode rings at 4.8245 GHz, Q 112 (field solver, fdamador_neutrino Part
    7b, grid r/2; 4.7021 GHz, Q 103 on r/3), not at the E-ring's 2.06 GHz;
  - the H23 electron itself is a travelling wave on a ring line (E and B in
    phase, README §2), which a lumped tank stands for one mode at a time,
    and in the field solver its light leaves within a quarter to half a
    period (README §3), i.e. Q of order 1.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
L, C, RRAD, RD = 2.156e-08, 2.767e-13, 6.862, 0.02791   # fdamador_neutrino run_output.txt
E_COL, B_COL, R_COL, K, G = "#d9480f", "#1c64b4", "#2b8a3e", "#222222", "#888888"


def coil(ax, p0, p1, col, n=4):
    p0, p1 = np.array(p0, float), np.array(p1, float)
    d = p1 - p0
    ln = np.hypot(*d)
    u, nrm = d / ln, np.array([-d[1], d[0]]) / ln
    s = np.linspace(0.12, ln - 0.12, 300)
    bump = 0.2 * np.abs(np.sin(np.pi * n * (s - 0.12) / (ln - 0.24)))
    pts = p0 + np.outer(s, u) + np.outer(bump, nrm)
    ax.plot(*zip(p0, p0 + 0.12 * u), color=K, lw=2); ax.plot(*zip(p1 - 0.12 * u, p1), color=K, lw=2)
    ax.plot(pts[:, 0], pts[:, 1], color=col, lw=2.4)


def zig(ax, p0, p1, col, n=5):
    p0, p1 = np.array(p0, float), np.array(p1, float)
    d = p1 - p0
    ln = np.hypot(*d)
    u, nrm = d / ln, np.array([-d[1], d[0]]) / ln
    s = np.linspace(0.12, ln - 0.12, 2 * n + 1)
    off = np.array([0] + [0.14 * (-1) ** k for k in range(2 * n - 1)] + [0])
    pts = p0 + np.outer(s, u) + np.outer(off, nrm)
    ax.plot(*zip(p0, p0 + 0.12 * u), color=K, lw=2); ax.plot(*zip(p1 - 0.12 * u, p1), color=K, lw=2)
    ax.plot(pts[:, 0], pts[:, 1], color=col, lw=2.2)


def cap(ax, p0, p1, col):
    p0, p1 = np.array(p0, float), np.array(p1, float)
    d = p1 - p0
    ln = np.hypot(*d)
    u, nrm = d / ln, np.array([-d[1], d[0]]) / ln
    m = (p0 + p1) / 2
    ax.plot(*zip(p0, m - 0.08 * u), color=K, lw=2); ax.plot(*zip(m + 0.08 * u, p1), color=K, lw=2)
    for sgn in (-1, 1):
        c = m + sgn * 0.08 * u
        ax.plot(*zip(c - 0.3 * nrm, c + 0.3 * nrm), color=col, lw=3.2)


def series(ax):
    x0, x1, y0, y1 = 0.3, 3.3, 0.5, 2.6
    coil(ax, (x0, y1), (x1, y1), B_COL)
    cap(ax, (x1, y1), (x1, y0), E_COL)
    zig(ax, (x1, y0), (x0, y0), R_COL)
    zig(ax, (x0, y0), (x0, y1), G, n=3)
    ax.text(1.8, y1 + 0.38, f"L = {L * 1e9:.2f} nH", ha="center", color=B_COL, weight="bold", fontsize=10.5)
    ax.text(1.8, y1 + 0.68, "B through the hole", ha="center", color=B_COL, fontsize=8.5)
    ax.text(x1 + 0.4, 1.7, f"C = {C * 1e12:.4f} pF", color=E_COL, weight="bold", fontsize=10.5)
    ax.text(x1 + 0.4, 1.42, "E round the ring", color=E_COL, fontsize=8.5)
    ax.text(1.8, y0 - 0.45, f"R_rad = {RRAD:.3f} Ω", ha="center", color=R_COL, weight="bold", fontsize=10.5)
    ax.text(x0 - 0.2, 1.55, f"R_d\n{RD * 1e3:.1f} mΩ", ha="right", color="#666", fontsize=9)
    th = np.linspace(0.3, 2 * np.pi - 0.3, 80)
    ax.plot(1.8 + 0.4 * np.cos(th), 1.55 + 0.4 * np.sin(th), color=K, lw=1)
    ax.annotate("", xy=(1.8 + 0.4 * np.cos(th[-1]), 1.55 + 0.4 * np.sin(th[-1])),
                xytext=(1.8 + 0.4 * np.cos(th[-3]), 1.55 + 0.4 * np.sin(th[-3])),
                arrowprops=dict(arrowstyle="-|>", color=K, lw=1))
    ax.text(1.8, 1.55, "I", ha="center", va="center", style="italic", fontsize=11)


def parallel(ax, rp, rpd):
    yt, yb = 2.9, 0.3
    xs = [0.5, 1.9, 3.4, 4.8]
    ax.plot([xs[0], xs[-1]], [yt, yt], color=K, lw=2)
    ax.plot([xs[0], xs[-1]], [yb, yb], color=K, lw=2)
    coil(ax, (xs[0], yt), (xs[0], yb), B_COL)
    cap(ax, (xs[1], yt), (xs[1], yb), E_COL)
    zig(ax, (xs[2], yt), (xs[2], yb), R_COL)
    zig(ax, (xs[3], yt), (xs[3], yb), G, n=3)
    for x in xs:
        ax.plot([x, x], [yt, yt], "o", color=K, ms=4); ax.plot(x, yb, "o", color=K, ms=4)
    ax.text(xs[0] - 0.3, 1.75, f"L = {L * 1e9:.2f} nH", ha="right", color=B_COL, weight="bold", fontsize=10.5)
    ax.text(xs[0] - 0.3, 1.47, "B round the ring\n(toroidal flux)", ha="right", va="top", color=B_COL, fontsize=8.5)
    ax.text(xs[1] + 0.38, 2.15, f"C =\n{C * 1e12:.4f} pF", color=E_COL, weight="bold", fontsize=10)
    ax.text(xs[1] + 0.38, 1.2, "E through\nthe hole", color=E_COL, fontsize=8.5)
    ax.text(xs[2] + 0.28, 2.15, f"R_p =\n{rp / 1e3:.2f} kΩ", color=R_COL, weight="bold", fontsize=10)
    ax.text(xs[2] + 0.28, 1.2, "radiation\n(R0²/R_rad)", color=R_COL, fontsize=8.5)
    ax.text(xs[3] + 0.28, 1.75, f"{rpd / 1e6:.2f} MΩ", color="#666", fontsize=9)
    ax.text(xs[3] + 0.28, 1.45, "(R0²/R_d)", color="#666", fontsize=8)
    ax.text(2.65, yt + 0.3, "V across every element", ha="center", fontsize=9)


def main(out):
    w0 = 1 / np.sqrt(L * C)
    f0 = w0 / (2 * np.pi)
    R0 = np.sqrt(L / C)
    rp, rpd = R0 ** 2 / RRAD, R0 ** 2 / RD
    Qs = w0 * L / (RRAD + RD)
    Gp = 1 / rp + 1 / rpd
    Qp = 1 / (w0 * L * Gp)                    # parallel: Q = R_total/(ω0 L)
    assert abs(Qp / Qs - 1) < 1e-12, (Qs, Qp)

    fig = plt.figure(figsize=(15, 5.6))
    fig.suptitle("The photon electron as the circuit mirror of the fdamador neutrino: "
                 "series RLC loop ↔ parallel RLC tank", fontsize=13)
    a1 = fig.add_axes([0.0, 0.22, 0.33, 0.68]); a1.set_xlim(-0.8, 5.0); a1.set_ylim(-0.3, 3.6)
    a2 = fig.add_axes([0.36, 0.22, 0.2, 0.68]); a2.set_xlim(0, 1); a2.set_ylim(0, 1)
    a3 = fig.add_axes([0.58, 0.22, 0.42, 0.68]); a3.set_xlim(-1.4, 6.2); a3.set_ylim(-0.3, 3.6)
    for a in (a1, a2, a3):
        a.axis("off")
    for a in (a1, a3):
        a.set_aspect("equal")
    series(a1)
    parallel(a3, rp, rpd)
    a1.set_title("fdamador neutrino (H28): series loop\nE is the ring, B circles it", fontsize=10.5)
    a3.set_title("photon electron (H23), its mirror: parallel tank\nB is the ring, E threads it", fontsize=10.5)
    a2.annotate("", xy=(0.95, 0.55), xytext=(0.05, 0.55), arrowprops=dict(arrowstyle="<|-|>", lw=2, color=K))
    a2.text(0.5, 0.95, "duality (E ↔ B)", ha="center", fontsize=11, weight="bold")
    a2.text(0.5, 0.86, "voltage ↔ current\nL ↔ C, series ↔ parallel\nR ↔ R0²/R", ha="center", va="top", fontsize=9.5)
    a2.text(0.5, 0.45, f"R0 = √(L/C) = {R0:.1f} Ω\nkeeps L, C, f0, Q", ha="center", va="top", fontsize=9.5)

    t = (f"Both: f0 = {f0 / 1e9:.4f} GHz, Q = {Qs:.2f} (series ω0L/(R_rad+R_d) = parallel R_p∥R_pd/(ω0L), equal to 10⁻¹²).  "
         "In the series loop the current is shared, in the tank the voltage.\n"
         "The mirror is exact as circuit algebra, for the models only. The real dielectric torus is not dual (εr = 100, μr = 1): "
         "its B-ring rings at 4.82 GHz, Q 112\n(field solver, fdamador_neutrino Part 7b, r/2), not at the E-ring's 2.06 GHz. "
         "The H23 electron itself is a travelling wave on a ring line (E and B in phase, README §2) at 1.24×10²⁰ Hz;\n"
         "a tank stands for one of its modes at a time. Under linear Maxwell its light leaves within ¼–½ period (README §3), "
         "Q of order 1. Neither circuit has mass, spin ½ from first principles or the weak interaction.")
    fig.text(0.02, 0.04, t, fontsize=9, va="bottom", color="#333")
    fig.savefig(out, dpi=130)
    print(f"wrote {out}  (f0 {f0:.6e} Hz, R0 {R0:.2f} Ω, R_p {rp:.1f} Ω, R_pd {rpd:.4e} Ω, Q {Qs:.4f} = {Qp:.4f})")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "circuit.png"))
