# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""The fountain of H49 Part 8a in the meridional plane.

    go run ./examples/ether_horizon -viz examples/ether_horizon/viz
    python3 examples/ether_horizon/viz/plot_fountain.py [OUT.png]

Writes viz/fountain.png from fountain_l{1,2}_psi.csv and fountain_l{1,2}_q2.csv
(part8.go): stream-function contours of the two steady swirling flows in a
sphere of radius a (Bragg–Hawthorne with C = √(C0² + k²ψ²)), arrows for the
poloidal flow, the axis vortex (the void) as a thick line, and the streamline
torus that winds twice round the axis per turn round its tube (q = 2, red),
fitted to tube/ring = ½; the photon_electron torus cross-section (ring R,
tube R/2) dashed for comparison.
"""
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))


def grid(path):
    with open(path) as f:
        rows = [[float(x) for x in r] for r in list(csv.reader(f))[1:]]
    a = np.array(rows)
    s = np.unique(a[:, 0]); z = np.unique(a[:, 1])
    psi = a[:, 2].reshape(len(s), len(z)).T
    return s, z, psi


def curve(path):
    with open(path) as f:
        rows = [[float(x) for x in r] for r in list(csv.reader(f))[1:]]
    return np.array(rows)


def panel(ax, l, title):
    s, z, psi = grid(os.path.join(HERE, f"fountain_l{l}_psi.csv"))
    q2 = curve(os.path.join(HERE, f"fountain_l{l}_q2.csv"))
    S, Z = np.meshgrid(s, z)
    m = np.nanmax(np.abs(psi))
    for sg in (1, -1):
        ax.contour(sg * S, Z, psi, levels=np.linspace(-m, m, 25), colors="#7b6fb0", linewidths=0.7, linestyles="solid")
    # Poloidal arrows from ψ: u_s = −(1/s)∂ψ/∂z, u_z = (1/s)∂ψ/∂s.
    dpz, dps = np.gradient(psi, z, s)
    with np.errstate(divide="ignore", invalid="ignore"):
        us, uz = -dpz / S, dps / S
    k = (slice(10, None, 14), slice(8, None, 14))
    for sign in (1, -1):
        ax.quiver(sign * S[k], Z[k], sign * us[k], uz[k], color="#378add", scale=60, width=0.004)
    t = np.linspace(0, 2 * np.pi, 300)
    ax.plot(np.cos(t), np.sin(t), color="#444", lw=1)
    ax.plot([0, 0], [-1, 1], color="#5f5e5a", lw=4, solid_capstyle="butt")
    for sign in (1, -1):
        ax.plot(sign * q2[:, 0], q2[:, 1], color="#d1361f", lw=2)
        if l == 2:
            ax.plot(sign * q2[:, 0], -q2[:, 1], color="#d1361f", lw=2)
    so = np.mean(q2[:, 0]); zo = np.mean(q2[:, 1])
    for sign in (1, -1):
        for zz in ([zo] if l == 1 else [zo, -zo]):
            ax.plot(sign * so + so / 2 * np.cos(t), zz + so / 2 * np.sin(t), color="#d1361f", lw=0.9, ls="--")
    ax.set_aspect("equal"); ax.set_xlim(-1.05, 1.05); ax.set_ylim(-1.05, 1.05)
    ax.set_xlabel("s / a (distance from the axis)"); ax.set_ylabel("z / a (along the axis)")
    ax.set_title(title, fontsize=10)


def main(out):
    fig, axs = plt.subplots(1, 2, figsize=(12.5, 6.4))
    panel(axs[0], 1, "(a) one cell (l = 1): in at the south pole, out at the north;\n"
                     "q = 2 torus (red) in the equatorial plane, a = 316 fm")
    panel(axs[1], 2, "(b) two cells (l = 2): in round the equator, out at both poles;\n"
                     "q = 2 tori above and below the equator, a = 352 fm")
    fig.suptitle("H49 Part 8a: the fountain — steady swirling ether flow with a void on the axis (grey); "
                 "red: the streamline torus that winds 2:1 like the photon ring", fontsize=11)
    fig.text(0.5, 0.015, "Purple: stream surfaces (cross-sections of nested tori); blue: poloidal flow; dashed red: the photon_electron "
             "torus (tube = ring/2) fitted to the same ring radius. Scaled so the ring radius is λ/4π = 193 fm.",
             ha="center", fontsize=8.5, color="#444")
    fig.tight_layout(rect=(0, 0.06, 1, 0.93))
    fig.savefig(out, dpi=110)
    print("wrote", out)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "fountain.png"))
