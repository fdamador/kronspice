# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""H49 Part 9: the round void in the two-pole fountain.

    go run ./examples/ether_horizon -viz examples/ether_horizon/viz
    python3 examples/ether_horizon/viz/plot_roundvoid.py [OUT.png]

Writes viz/round_void.png from viz/round_void_shapes.csv:
  (a) the void (radius b) at the stagnation point of the two-pole flow, with
      the straining flow round it, φ = E(r² + (2/3)b⁵/r³)P2(cos θ): in round
      the equator, out along the axis; the light inside pushing out (orange),
      the ether pushing in (grey arrows); the axis vortex's hollow tails;
  (b) the wall's shape at Weber numbers We = ρ(Eb)²b/σ = 0, 0.1, 0.3, 1
      (linearized Young–Laplace): the flow stretches the void along the axis.
"""
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))


def main(out):
    with open(os.path.join(HERE, "round_void_shapes.csv")) as f:
        rows = list(csv.DictReader(f))
    th = np.array([float(r["theta"]) for r in rows])
    fig, axs = plt.subplots(1, 2, figsize=(12.5, 6.2))

    ax = axs[0]
    b = 1.0
    x = np.linspace(-3, 3, 400); z = np.linspace(-3, 3, 400)
    X, Z = np.meshgrid(x, z)
    R = np.hypot(X, Z)
    S = np.abs(X)
    # Stokes stream function of the straining flow round a sphere:
    # ψ = −(3/2)·E·s²·z·(1 − b⁵/r⁵)   (u_z = (1/s)∂ψ/∂s outward along the axis)
    psi = 1.5 * S ** 2 * Z * (1 - b ** 5 / np.maximum(R, 1e-9) ** 5) * np.sign(X + 1e-12)
    psi[R < b] = np.nan
    ax.contour(X, Z, psi, levels=np.linspace(-12, 12, 33), colors="#378add", linewidths=0.8, linestyles="solid")
    for sx, sz, dx, dz in [(2.6, 0.15, -0.5, 0), (-2.6, 0.15, 0.5, 0), (2.6, -0.15, -0.5, 0), (-2.6, -0.15, 0.5, 0),
                           (0.25, 1.6, 0, 0.6), (-0.25, 1.6, 0, 0.6), (0.25, -1.6, 0, -0.6), (-0.25, -1.6, 0, -0.6)]:
        ax.annotate("", xy=(sx + dx, sz + dz), xytext=(sx, sz), arrowprops=dict(arrowstyle="->", color="#1c64b4", lw=1.6))
    t = np.linspace(0, 2 * np.pi, 200)
    ax.fill(np.cos(t), np.sin(t), color="#fff4e0", ec="#5f5e5a", lw=1.5)
    for a in np.linspace(0, 2 * np.pi, 12, endpoint=False):
        ax.annotate("", xy=(0.92 * np.cos(a), 0.92 * np.sin(a)), xytext=(0.55 * np.cos(a), 0.55 * np.sin(a)),
                    arrowprops=dict(arrowstyle="->", color="#d9480f", lw=1.2))
        ax.annotate("", xy=(1.08 * np.cos(a), 1.08 * np.sin(a)), xytext=(1.38 * np.cos(a), 1.38 * np.sin(a)),
                    arrowprops=dict(arrowstyle="->", color="#888", lw=1.0))
    ax.plot([0, 0], [1, 3], color="#5f5e5a", lw=3); ax.plot([0, 0], [-1, -3], color="#5f5e5a", lw=3)
    ax.text(0.12, 2.3, "hollow tail\n(axis vortex)", fontsize=8.5, color="#444")
    ax.text(0, -0.1, "void\nlight pushes out", ha="center", fontsize=9, color="#9a3412")
    ax.text(1.55, 1.1, "ether pushes in\n(p∞, σ)", fontsize=8.5, color="#555")
    ax.set_aspect("equal"); ax.set_xlim(-3, 3); ax.set_ylim(-3, 3)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title("(a) a round void at the centre of the two-pole fountain:\nin round the equator, out along the axis", fontsize=10)

    ax = axs[1]
    cols = {"we0": ("#222", "We = 0 (round)"), "we0.1": ("#378add", "We = 0.1"), "we0.3": ("#1d9e75", "We = 0.3"), "we1": ("#d1361f", "We = 1 (beyond linear)")}
    for k, (c, lab) in cols.items():
        r = np.array([float(row[k]) for row in rows])
        for sgn in (1, -1):
            ax.plot(sgn * r * np.sin(th), r * np.cos(th), color=c, lw=1.8 if k != "we1" else 1.2,
                    ls="-" if k != "we1" else "--", label=lab if sgn == 1 else None)
    ax.axhline(0, color="#ccc", lw=0.6); ax.axvline(0, color="#ccc", lw=0.6)
    ax.set_aspect("equal"); ax.set_xlim(-1.5, 1.5); ax.set_ylim(-1.5, 1.5)
    ax.set_xlabel("s / b"); ax.set_ylabel("z / b (axis)")
    ax.legend(fontsize=8.5, loc="lower right")
    ax.set_title("(b) the wall's shape: the flow stretches the void along the axis\n"
                 "(10 % longer than wide at We = ρ(Eb)²b/σ = 0.27)", fontsize=10)
    fig.suptitle("H49 Part 9: a round void held by its light against the nonlinear ether, in the fountain flow", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(out, dpi=110)
    print("wrote", out)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "round_void.png"))
