# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""One beat of two fdamador rings (Part 9, H42), animated: viz/pair_cycle.gif.

    python3 examples/fdamador_neutrino/viz/plot_pair_cycle.py [OUT.gif]

Two bench rings (a = 10 mm, r = 2.5 mm) on a common axis d = a apart.
Ring 0 starts with all the energy (its displacement, E round the ring, at a
maximum), ring 1 at rest. The state is the sum of the in-phase mode
(1.7957 GHz) and the antiphase mode (2.3056 GHz) of Part 9a, so the energy
passes to ring 1 and back once per beat, 1/Δf = 4.02 periods; one beat is
shown, slowed so it lasts about 5 s.

  (a) cross-section through the common axis: each ring's E (into or out of
      the page, symbol size ∝ |E|, colour intensity ∝ its energy) and the B
      lines of both rings' displacement currents together (contours of
      I0·ψ(z − d/2) + I1·ψ(z + d/2), ψ = ρ·A_φ as in plot_structure.py);
  (b) the energy in each ring and the cursor.

The decay is left out (the in-phase half radiates with Q 41.8, the
antiphase half with Q 917; Part 9a shows the ring-down: after a few beats
only the antiphase half is left and the exchange stops). Amplitudes are
normalised.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.animation import FuncAnimation, PillowWriter  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_structure import A_MM, B_COL, E_COL, R_MM, psi  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
D_MM = 10.0
F_IN, F_ANTI = 1.7957e9, 2.3056e9
NF = 120


def state(t):
    """Displacement q (∝ E) and current I (∝ B) of each ring at time t [s], unit start."""
    ws, wa = 2 * np.pi * F_IN, 2 * np.pi * F_ANTI
    q0 = 0.5 * (np.cos(ws * t) + np.cos(wa * t))
    q1 = 0.5 * (np.cos(ws * t) - np.cos(wa * t))
    i0 = -0.5 * (ws * np.sin(ws * t) + wa * np.sin(wa * t))
    i1 = -0.5 * (ws * np.sin(ws * t) - wa * np.sin(wa * t))
    wm = 0.5 * (ws + wa)
    e0 = 0.5 * (q0 ** 2 + (i0 / wm) ** 2)
    e1 = 0.5 * (q1 ** 2 + (i1 / wm) ** 2)
    return q0, q1, i0 / wm, i1 / wm, e0, e1


def main(out):
    tb = 1 / (F_ANTI - F_IN)
    ts = np.linspace(0, tb, NF, endpoint=False)
    tt = np.linspace(0, tb, 2000)
    _, _, _, _, E0, E1 = state(tt)
    tot = E0 + E1

    xs = np.linspace(-30, 30, 241)
    zs = np.linspace(-24, 24, 193)
    XX, ZZ = np.meshgrid(xs, zs)
    RR = np.maximum(np.abs(XX), 1e-3)
    P0 = psi(RR, ZZ - D_MM / 2)
    P1 = psi(RR, ZZ + D_MM / 2)
    pmax = np.max(np.abs(P0[np.hypot(np.abs(XX) - A_MM, ZZ - D_MM / 2) > R_MM]))
    levels = np.concatenate([-pmax * np.array([0.9, 0.6, 0.4, 0.25, 0.15, 0.08]),
                             pmax * np.array([0.08, 0.15, 0.25, 0.4, 0.6, 0.9])])

    fig = plt.figure(figsize=(12.5, 6.0))
    fig.suptitle("Two fdamador rings (H42), d = a: the energy passes to the other ring and back once per beat "
                 "(in-phase + antiphase mode)", fontsize=11.5)
    ax = fig.add_axes([0.07, 0.1, 0.42, 0.72])
    bx = fig.add_axes([0.56, 0.14, 0.41, 0.7])

    bx.plot(tt / tb, E0 / tot, color=E_COL, label="ring 0 (top, started)")
    bx.plot(tt / tb, E1 / tot, color=B_COL, label="ring 1 (bottom, at rest)")
    bx.set_xlabel("time [beats]  (1 beat = 1/Δf = 4.02 periods = 1.96 ns)")
    bx.set_ylabel("energy in the ring / total")
    bx.set_xlim(0, 1)
    bx.set_ylim(0, 1.08)
    bx.legend(fontsize=9, loc="center right", bbox_to_anchor=(1.0, 0.3))
    cur = bx.axvline(0, color="#444", lw=1)
    bx.set_title("(b) energy per ring: complete exchange (no decay shown)", fontsize=10)

    def draw(k):
        t = ts[k]
        q0, q1, i0, i1, e0, e1 = state(t)
        ax.clear()
        P = i0 * P0 + i1 * P1
        ax.contour(XX, ZZ, P, levels=levels, colors=B_COL, linewidths=0.9, linestyles="solid")
        for zc, q, e in ((D_MM / 2, q0, e0), (-D_MM / 2, q1, e1)):
            a = 0.15 + 0.85 * min(1, e / 0.5)
            for sx in (1, -1):
                ax.add_patch(plt.Circle((sx * A_MM, zc), R_MM, facecolor=E_COL, alpha=0.12 + 0.5 * a, ec=E_COL, lw=1.4, zorder=3))
                into = sx * q > 0
                ms = 3 + 13 * abs(q)
                ax.plot(sx * A_MM, zc, marker="x" if into else "o", color="#7a2400", ms=ms, mew=2,
                        mfc="#7a2400" if into else "none", zorder=4)
        ax.axvline(0, color="#999", lw=0.6, ls=":")
        ax.set_aspect("equal")
        ax.set_xlim(-30, 30)
        ax.set_ylim(-24, 24)
        ax.set_xlabel("ρ [mm] (cut through the common axis)")
        ax.set_ylabel("z [mm]")
        ax.set_title(f"(a) E (× into, ○ out of the page; size ∝ |E|) and B lines\nt = {t / tb:.2f} beat", fontsize=10)
        cur.set_xdata([t / tb, t / tb])
        return []

    anim = FuncAnimation(fig, draw, frames=NF, blit=False)
    anim.save(out, writer=PillowWriter(fps=24), dpi=70)
    print("wrote", out)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "pair_cycle.gif"))
