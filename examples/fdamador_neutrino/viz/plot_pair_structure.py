# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Two fdamador rings (Part 9, H42): viz/structure_pair.png.

    python3 examples/fdamador_neutrino/viz/plot_pair_structure.py [OUT.png]

Same geometry and drawing as viz/structure.png (bench ring a = 10 mm,
r = 2.5 mm), two rings on a common axis d = a apart:

  (a) in phase: E round both rings the same way at the same instant; the
      B lines (contours of ψ1 + ψ2, ψ = ρ·A_φ of each ring, as in
      plot_structure.py) close round both rings together: one larger
      magnetic dipole, which radiates more (superradiant);
  (b) antiphase: E round the rings in opposite senses; the B lines of
      each ring close on themselves and between the rings they cancel on
      the axis: the net dipole is zero and what is left is a quadrupole
      (subradiant);
  (c) radiation pattern (far field, power per solid angle, same current
      per ring) of one ring, the in-phase pair and the antiphase pair:
      sin²θ·|1 ± e^{ikd·cosθ}|² with kd at each mode's own frequency (Part 9a:
      1.7957 and 2.3056 GHz at d = a).

Amplitudes are normalised; geometry and frequencies are the README's.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_structure import A_MM, B_COL, E_COL, R_MM, psi  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
D_MM = 10.0                       # centre distance (d = a)
F_IN, F_ANTI = 1.7957e9, 2.3056e9  # Part 9a, coaxial d = a
C0 = 299792458.0


def panel(ax, sign, title):
    xs = np.linspace(-30, 30, 301)
    zs = np.linspace(-26, 26, 261)
    XX, ZZ = np.meshgrid(xs, zs)
    rr = np.maximum(np.abs(XX), 2e-3)
    h = 1e-3

    def P(r, z):
        return psi(r, z - D_MM / 2) + sign * psi(r, z + D_MM / 2)

    Bz = (P(rr + h, ZZ) - P(rr - h, ZZ)) / (2 * h) / rr
    Bx = -np.sign(XX) * (P(rr, ZZ + h) - P(rr, ZZ - h)) / (2 * h) / rr
    for zc in (D_MM / 2, -D_MM / 2):
        inside = np.hypot(np.abs(XX) - A_MM, ZZ - zc) < R_MM * 0.35
        Bx[inside] = np.nan
        Bz[inside] = np.nan
    # field lines as contours of ψ (exact, closed), on both halves of the cut
    RR = np.abs(XX)
    PP = P(np.maximum(RR, 1e-3), ZZ)
    lv = []
    for zc, s in ((D_MM / 2, 1), (-D_MM / 2, sign)):
        for dd in (1.0, 2.2, 4.0, 6.5, 10.0):
            lv.append(P(np.array([A_MM + dd]), np.array([zc]))[0])
        for dz in (6.0, 10.0, 15.0):
            lv.append(P(np.array([A_MM]), np.array([zc + np.sign(zc) * dz]))[0])
    for dd in (8, 13, 19):
        lv.append(P(np.array([A_MM + dd]), np.array([0.0]))[0])
    lv = np.array(lv)
    if sign < 0:
        lv = np.concatenate([lv, -lv])
    lv = np.unique(np.round(lv[np.abs(lv) > 1e-9], 9))
    mag = np.hypot(Bx, Bz)
    ax.contour(XX, ZZ, PP, levels=lv, colors=B_COL, linewidths=1.0, linestyles="solid")
    ax.quiver(XX[::20, ::20], ZZ[::20, ::20], (Bx / (mag + 1e-30) * (mag / np.nanmax(mag)) ** 0.25)[::20, ::20],
              (Bz / (mag + 1e-30) * (mag / np.nanmax(mag)) ** 0.25)[::20, ::20], color=B_COL, alpha=0.6, scale=28, width=0.004)
    for zc, s in ((D_MM / 2, 1), (-D_MM / 2, sign)):
        for sx in (1, -1):
            ax.add_patch(plt.Circle((sx * A_MM, zc), R_MM, color="#f2c6a8", ec=E_COL, lw=1.5, zorder=3))
            into = (sx * s) > 0
            ax.plot(sx * A_MM, zc, marker="x" if into else "o", color=E_COL, ms=10, mew=2.3,
                    mfc=E_COL if into else "none", zorder=4)
            if not into:
                ax.plot(sx * A_MM, zc, ".", color=E_COL, ms=5, zorder=4)
    ax.axvline(0, color="#999", lw=0.6, ls=":")
    ax.set_aspect("equal")
    ax.set_xlim(-30, 30)
    ax.set_ylim(-26, 26)
    ax.set_xlabel("ρ [mm] (cut through the common axis)")
    ax.set_ylabel("z [mm]")
    ax.set_title(title, fontsize=10)


def main(out):
    fig = plt.figure(figsize=(16, 5.8))
    fig.suptitle("Two fdamador rings (H42) on a common axis, d = a: B lines of the in-phase and antiphase modes",
                 fontsize=13)
    panel(fig.add_subplot(1, 3, 1), +1,
          f"(a) in phase ({F_IN / 1e9:.3f} GHz): the lines close round both rings;\n"
          "one larger dipole, superradiant (Q 41.8 against 48.9)")
    panel(fig.add_subplot(1, 3, 2), -1,
          f"(b) antiphase ({F_ANTI / 1e9:.3f} GHz): each ring's lines close on itself,\n"
          "net dipole 0, a quadrupole: subradiant (Q 917)")
    ax = fig.add_subplot(1, 3, 3, projection="polar")
    th = np.linspace(0, 2 * np.pi, 721)
    one = np.sin(th) ** 2
    kin = 2 * np.pi * F_IN / C0 * D_MM * 1e-3
    kan = 2 * np.pi * F_ANTI / C0 * D_MM * 1e-3
    pin = np.sin(th) ** 2 * np.abs(1 + np.exp(1j * kin * np.cos(th))) ** 2
    pan = np.sin(th) ** 2 * np.abs(1 - np.exp(1j * kan * np.cos(th))) ** 2
    ax.plot(th, pin / 4, color=B_COL, lw=3, alpha=0.7, label="pair in phase (÷4)")
    ax.plot(th, one, color="#222222", lw=1.2, ls="--", label="one ring (magnetic dipole)")
    ax.plot(th, pan * 10, color=E_COL, lw=1.6, label="pair antiphase (×10)")
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_rticks([])
    ax.legend(fontsize=8, loc="upper right", bbox_to_anchor=(1.18, 0.08))
    tot = lambda p: np.trapezoid(p * np.abs(np.sin(th)), th) / 2  # noqa: E731
    ax.set_title(f"(c) far-field power vs angle from the axis (same current\nper ring): antiphase total {tot(pan) / tot(one):.3f} of one ring,"
                 f" in phase {tot(pin) / tot(one):.2f}", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(out, dpi=110)
    print("wrote", out)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "structure_pair.png"))
