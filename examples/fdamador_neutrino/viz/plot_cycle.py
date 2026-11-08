# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""One cycle of the fdamador neutrino (H28), animated.

    python3 examples/fdamador_neutrino/viz/plot_cycle.py [OUT.gif]

Writes viz/cycle.gif (default): one period of the README's bench ring
(a = 10 mm, r = 2.5 mm, εr = 100, f0 = 2.0606 GHz, T = 0.485 ns), slowed so
the loop lasts about 3 s. Same geometry and B lines as viz/structure.png:

  (a) 3-D: E round the ring (orange arrows, length ∝ E) and B through the
      hole and round the tube (blue lines, strength ∝ |B|, arrow on the axis
      gives its sign);
  (b) cross-section: the B field as arrows, E into or out of the page at the
      tube;
  (c) E and B a quarter period apart, and the energy passing from electric
      (½CV², the ring's displacement) to magnetic (½LI², the moment) and back.

The decay is left out (Q = 40.5 loses 14 % of the energy per cycle); this
is the shape of one cycle, the ring-down is viz/structure.png panel (c).
Amplitudes are normalised.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.animation import FuncAnimation, PillowWriter  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_structure import A_MM, B_COL, E_COL, F0, R_MM, b_lines, psi  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
NF = 48


def main(out):
    rho, z, P, levels = b_lines()
    cs = plt.figure().add_subplot().contour(rho, z, P, levels=levels)
    segs = [s for sl in cs.allsegs for s in sl]
    plt.close("all")

    # cross-section B on a coarse grid (from ψ, as in plot_structure)
    xs = np.linspace(-26, 26, 27)
    zs = np.linspace(-18, 18, 19)
    XX, ZZ = np.meshgrid(xs, zs)
    h = 1e-3
    rr = np.maximum(np.abs(XX), 2e-3)
    Bz = (psi(rr + h, ZZ) - psi(rr - h, ZZ)) / (2 * h) / rr
    Bx = -np.sign(XX) * (psi(rr, ZZ + h) - psi(rr, ZZ - h)) / (2 * h) / rr
    mag = np.hypot(Bx, Bz)
    near = np.hypot(np.abs(XX) - A_MM, ZZ) < R_MM * 1.1
    Bx, Bz = np.where(near, np.nan, Bx), np.where(near, np.nan, Bz)
    # arrow length ∝ |B|^0.3, so the weak outer field stays visible
    g = (mag / np.nanmax(np.where(near, np.nan, mag))) ** 0.3 / (mag + 1e-30)
    Bx, Bz = Bx * g, Bz * g

    fig = plt.figure(figsize=(14, 5.4))
    fig.suptitle("The fdamador neutrino (H28): one cycle — E round the ring and B through it, "
                 "a quarter period apart", fontsize=12.5)
    ax3 = fig.add_subplot(1, 3, 1, projection="3d")
    axb = fig.add_subplot(1, 3, 2)
    axc = fig.add_subplot(1, 3, 3)
    T = 1 / F0
    tt = np.linspace(0, 1, 400)

    u, v = np.meshgrid(np.linspace(0, 2 * np.pi, 48), np.linspace(0, 2 * np.pi, 20))
    torus = ((A_MM + R_MM * np.cos(v)) * np.cos(u), (A_MM + R_MM * np.cos(v)) * np.sin(u), R_MM * np.sin(v))
    th = np.linspace(0, 2 * np.pi, 10, endpoint=False)

    def draw(k):
        ph = 2 * np.pi * k / NF
        e, b = np.cos(ph), np.sin(ph)

        ax3.cla()
        ax3.plot_surface(*torus, color="#f2c6a8", alpha=0.3, linewidth=0)
        for phi in np.linspace(0, 2 * np.pi, 8, endpoint=False):
            for s in segs:
                ax3.plot(s[:, 0] * np.cos(phi), s[:, 0] * np.sin(phi), s[:, 1], color=B_COL,
                         lw=0.4 + 1.2 * abs(b), alpha=0.08 + 0.75 * abs(b))
        ax3.quiver(A_MM * np.cos(th), A_MM * np.sin(th), np.zeros_like(th),
                   -np.sin(th) * 7 * e, np.cos(th) * 7 * e, np.zeros_like(th),
                   color=E_COL, lw=2.2, arrow_length_ratio=0.45)
        if abs(b) > 0.05:
            ax3.quiver(0, 0, -13 * b, 0, 0, 26 * b, color=B_COL, lw=2.4, arrow_length_ratio=0.15)
        ax3.set_xlim(-22, 22); ax3.set_ylim(-22, 22); ax3.set_zlim(-16, 16)
        ax3.set_box_aspect((1, 1, 0.75))
        ax3.set_xticks([]); ax3.set_yticks([]); ax3.set_zticks([])
        ax3.view_init(elev=24, azim=-60 + 360 * k / NF * 0.25)
        ax3.set_title("(a) E round the ring (orange), B through the hole (blue)\n"
                      "a = 10 mm, εr = 100", fontsize=9.5)

        axb.cla()
        axb.quiver(XX, ZZ, Bx * b, Bz * b, color=B_COL, scale=14, width=0.004, alpha=0.25 + 0.75 * abs(b))
        for sx in (1, -1):
            axb.add_patch(plt.Circle((sx * A_MM, 0), R_MM, color="#f2c6a8", ec=E_COL, lw=1.5))
            into = (sx * e) > 0
            ms = 4 + 14 * abs(e)
            if into:
                axb.plot(sx * A_MM, 0, "x", color=E_COL, ms=ms, mew=2.5)
            else:
                axb.plot(sx * A_MM, 0, "o", color=E_COL, ms=ms, mew=2.5, mfc="none")
                axb.plot(sx * A_MM, 0, ".", color=E_COL, ms=ms / 2.5)
        axb.set_aspect("equal"); axb.set_xlim(-28, 28); axb.set_ylim(-20, 20)
        axb.set_xlabel("ρ [mm] (cut through the axis)"); axb.set_ylabel("z [mm]")
        axb.set_title("(b) cross-section: B arrows; E into (×) or out of (○)\nthe page at the tube",
                      fontsize=9.5)

        axc.cla()
        axc.plot(tt, np.cos(2 * np.pi * tt), color=E_COL, lw=1.6, label="E round the ring")
        axc.plot(tt, np.sin(2 * np.pi * tt), color=B_COL, lw=1.6, label="B through the hole (moment m)")
        axc.axvline(k / NF, color="#333", lw=1)
        axc.plot(k / NF, e, "o", color=E_COL); axc.plot(k / NF, b, "o", color=B_COL)
        axc.bar([1.12, 1.26], [e * e, b * b], width=0.11, color=[E_COL, B_COL], alpha=0.8)
        axc.text(1.12, -0.13, "U_E", ha="center", fontsize=8.5)
        axc.text(1.26, -0.13, "U_B", ha="center", fontsize=8.5)
        axc.axvline(1.04, color="#bbb", lw=0.6)
        axc.set_xlim(0, 1.34); axc.set_ylim(-1.2, 1.25)
        axc.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
        axc.set_xlabel(f"time [periods; T = {T * 1e9:.3f} ns]")
        axc.legend(fontsize=7.5, loc="lower left")
        axc.grid(alpha=0.3)
        axc.set_title(f"(c) f0 = {F0 / 1e9:.4f} GHz: energy passes electric → magnetic\n"
                      f"→ electric; t = {k / NF:.2f} T", fontsize=9.5)
        return []

    fig.text(0.5, 0.015, "Bench ring of the README (§1). Decay left out (Q = 40.5 loses 14 % per cycle); "
             "slowed about 6×10⁹× (0.485 ns → 3 s). The average moment over a cycle is zero.", ha="center", fontsize=8.5, color="#444")
    fig.tight_layout(rect=(0, 0.04, 1, 0.93))
    anim = FuncAnimation(fig, draw, frames=NF)
    anim.save(out, writer=PillowWriter(fps=16), dpi=80)
    print("wrote", out)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "cycle.gif"))
