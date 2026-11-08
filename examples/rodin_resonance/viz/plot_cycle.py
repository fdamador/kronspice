# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""One cycle of the sine/cosine drive of the Rodin pair, animated.

    python3 examples/rodin_resonance/viz/plot_cycle.py [OUT.gif]

Writes viz/cycle.gif from viz/field_planes.csv (the per-ampere fields of
windings A and B on the planes z = 0 and y = 0) and viz/drives.csv (the
steady-state current phasors of each drive, from the kronspice transient
at 1 kHz). Columns: the quadrature drive in the toroidal sense and in the
dipole sense; rows: the mid-plane z = 0 and the vertical plane y = 0.
Colour is |B(t)|, arrows the field's direction in the plane at that
instant. One loop is one period (1 ms, slowed 3000× to 3 s).

What to look for: in the toroidal sense the field inside the tube (round
the axis) and the field through the hole peak a quarter cycle apart, so at
each point the arrow turns; the pattern as a whole does not travel round
the axis (it keeps the 12-fold symmetry of the poles at every instant).
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.animation import FuncAnimation, PillowWriter  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_all import INK, SEQ, drives, load_planes  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
NF = 36
CASES = [("A + B quadrature (sin/cos), toroidal sense", "sin/cos, toroidal sense"),
         ("A + B quadrature (sin/cos), dipole sense", "sin/cos, dipole sense")]


def main(out):
    P = load_planes()
    D = drives()
    fig, axs = plt.subplots(2, 2, figsize=(7.2, 7))
    panels = []
    for j, (case, label) in enumerate(CASES):
        ia, ib = D[case]
        for i, (pl, (u, v)) in enumerate((("xy", (0, 1)), ("xz", (0, 2)))):
            n, pts, bA, bB, near = P[pl]
            c = ia * bA + ib * bB
            peak = np.sqrt(np.sum(np.abs(c) ** 2, axis=1))
            vmax = np.percentile(peak[~near], 99) * 1e3
            ext = [pts[:, u].min() * 100, pts[:, u].max() * 100, pts[:, v].min() * 100, pts[:, v].max() * 100]
            ax = axs[i, j]
            im = ax.imshow(np.zeros((n, n)), origin="lower", extent=ext, cmap=SEQ, vmin=0, vmax=vmax, interpolation="nearest")
            s = slice(None, None, max(1, n // 16))
            X = pts[:, u].reshape(n, n).T * 100
            Y = pts[:, v].reshape(n, n).T * 100
            q = ax.quiver(X[s, s], Y[s, s], np.zeros_like(X[s, s]), np.zeros_like(X[s, s]), color=INK, scale=24, width=0.004)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title(f"{label}, {'z = 0' if pl == 'xy' else 'y = 0'}", fontsize=8)
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02, label="|B| [mT]")
            panels.append((im, q, c, near, n, u, v, s))
    txt = fig.suptitle("", fontsize=9)

    def frame(k):
        wt = 2 * np.pi * k / NF
        for im, q, c, near, n, u, v, s in panels:
            b = np.real(c * np.exp(1j * wt))
            mag = np.sqrt(np.sum(b ** 2, axis=1)) * 1e3
            mag[near] = np.nan
            im.set_data(mag.reshape(n, n).T)
            U = b[:, u].reshape(n, n).T
            V = b[:, v].reshape(n, n).T
            L = np.hypot(U, V) + 1e-30
            q.set_UVC((U / L)[s, s], (V / L)[s, s])
        txt.set_text(f"Rodin pair, sine on A, cosine on B, 1 kHz — t = {k}/{NF} of a period")
        return []

    anim = FuncAnimation(fig, frame, frames=NF, blit=False)
    fig.tight_layout()
    anim.save(out, writer=PillowWriter(fps=12), dpi=72)
    print("wrote", out)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "cycle.gif"))
