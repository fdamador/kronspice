# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""The dual-motor homopolar machine (dual.go), side by side and stacked.

Cross-sections through the axes at the real dimensions (disks r_o 15 cm,
coils r_c 18 cm, 6 cm long, 4 turns; stacked gap 8 cm). The field lines are
the coils' own field, from the exact field of a circular loop (complete
elliptic integrals), equal coil currents; colour = |B| (log). The numbers in
the titles are from `go run ./examples/homopolar -dual`.

Usage: python viz/plot_dual.py [out_dir]   (matplotlib, numpy)
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.patches import FancyArrowPatch, Rectangle

here = os.path.dirname(os.path.abspath(__file__))
out = sys.argv[1] if len(sys.argv) > 1 else here
MU0 = 4e-7 * np.pi
RO, RI, T = 0.15, 0.01, 0.005
RC, LC, N = 0.18, 0.06, 4
GAP = 0.08
TEAL, ORANGE, COPPER, INK, QUIET = "#127b8a", "#d9772b", "#b87333", "#1b1b1f", "#5b5b66"
BOX = dict(boxstyle="round,pad=0.25", facecolor="#fbfaf7", edgecolor="none", alpha=0.92)


def ellipke(m):
    """Complete elliptic integrals K(m), E(m) by the AGM."""
    a, b, c = np.ones_like(m), np.sqrt(1 - m), np.sqrt(m)
    s = 0.5 * c**2
    p = 0.5
    for _ in range(30):
        a, b, c = 0.5 * (a + b), np.sqrt(a * b), 0.5 * (a - b)
        p *= 2
        s += p * c**2 * 0.5
    K = np.pi / (2 * a)
    return K, K * (1 - s)


def loop_field(r, z, a, z0, I=1.0):
    """B_r, B_z of a circular loop of radius a at height z0 (axis r = 0)."""
    zz = z - z0
    rr = np.maximum(np.abs(r), 1e-6)
    m = 4 * a * rr / ((a + rr) ** 2 + zz**2)
    m = np.clip(m, 0, 1 - 1e-12)
    K, E = ellipke(m)
    den = np.sqrt((a + rr) ** 2 + zz**2)
    q = (a - rr) ** 2 + zz**2
    bz = MU0 * I / (2 * np.pi) / den * (K + (a * a - rr * rr - zz * zz) / q * E)
    br = MU0 * I / (2 * np.pi) * zz / (rr * den) * (-K + (a * a + rr * rr + zz * zz) / q * E)
    return np.sign(r) * br, bz


def coil_field(x, z, xc, zc, I=1.0):
    """Field in the x–z plane of an N-turn coil on the vertical axis x = xc."""
    bx = np.zeros_like(x)
    bz = np.zeros_like(x)
    for k in range(N):
        z0 = zc + LC * ((k + 0.5) / N - 0.5)
        brk, bzk = loop_field(x - xc, z, RC, z0, I)
        bx += brk
        bz += bzk
    return bx, bz


def draw_unit(ax, xc, zd, coil_col, coil_label, disk_label, loop_col, motor_up):
    # disk
    ax.add_patch(Rectangle((xc - RO, zd - T / 2 - 0.003), 2 * RO, T + 0.006, facecolor=COPPER, edgecolor=INK, lw=0.8, zorder=5))
    ax.text(xc - RO / 2, zd + 0.016, disk_label, ha="center", va="bottom", fontsize=8.5, fontweight="bold", color=INK, zorder=8, bbox=BOX)
    # coil cross-sections (both sides of the axis)
    for s in (-1, 1):
        ax.add_patch(Rectangle((xc + s * RC - 0.01, zd - LC / 2), 0.02, LC, facecolor=coil_col, edgecolor=INK, lw=0.6, alpha=0.95, zorder=5))
        for k in range(N):
            zz = zd - LC / 2 + LC * (k + 0.5) / N
            ax.plot(xc + s * RC, zz, marker="o" if s > 0 else "x", ms=3.2, mfc="white", mec=INK, mew=0.7, color=INK, zorder=6)
    ax.text(xc + RC + 0.022, zd, coil_label, ha="left", va="center", fontsize=7.5, color=coil_col, fontweight="bold", zorder=8, bbox=BOX)
    # shaft and motor
    sgn = 1 if motor_up else -1
    ax.plot([xc, xc], [zd, zd + sgn * 0.16], color=INK, lw=3, solid_capstyle="butt", zorder=4)
    ax.add_patch(Rectangle((xc - 0.045, zd + sgn * 0.16 - (0 if motor_up else 0.06)), 0.09, 0.06, facecolor="#e7e5df", edgecolor=INK, lw=0.8, zorder=5))
    ax.text(xc, zd + sgn * 0.16 + (0.03 if motor_up else -0.03), "motor", ha="center", va="center", fontsize=7.5, color=INK, zorder=6)


def field_panel(ax, xs, zs, coils):
    x, z = np.meshgrid(xs, zs)
    bx = np.zeros_like(x)
    bz = np.zeros_like(x)
    for xc, zc, I in coils:
        a, b = coil_field(x, z, xc, zc, I)
        bx += a
        bz += b
    mag = np.hypot(bx, bz)
    ax.streamplot(x, z, bx, bz, color=mag, cmap="Blues", norm=LogNorm(mag.min() * 3, mag.max() / 3), density=1.5,
                  linewidth=0.7, arrowsize=0.6, zorder=1)


fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 6.3), dpi=300, gridspec_kw={"width_ratios": [1.45, 1]})
for ax in (a1, a2):
    ax.set_facecolor("#fbfaf7")
fig.patch.set_facecolor("#fbfaf7")

# (a) side by side: unit 1 at x = −0.33, unit 2 at x = +0.33. Coil 2 (I2) on unit 1, coil 1 (I1) on unit 2.
X = 0.33
field_panel(a1, np.linspace(-0.72, 0.72, 300), np.linspace(-0.3, 0.38, 160), [(-X, 0, 1.0), (X, 0, 1.0)])
draw_unit(a1, -X, 0, ORANGE, "coil 2: I₂", "disk 1", TEAL, False)
draw_unit(a1, X, 0, TEAL, "coil 1: I₁", "disk 2", ORANGE, False)
# cross wiring
a1.add_patch(FancyArrowPatch((-X + RO, 0.0), (X - RC - 0.01, 0.02), connectionstyle="arc3,rad=-0.35", arrowstyle="-|>", mutation_scale=12, color=TEAL, lw=2.2, zorder=7))
a1.add_patch(FancyArrowPatch((X - RO, 0.0), (-X + RC + 0.01, 0.02), connectionstyle="arc3,rad=0.35", arrowstyle="-|>", mutation_scale=12, color=ORANGE, lw=2.2, zorder=7))
a1.text(0, 0.2, "disk 1's current (I₁) feeds coil 1 on unit 2,\ndisk 2's current (I₂) feeds coil 2 on unit 1", ha="center", fontsize=8, color=QUIET, zorder=8, bbox=BOX)
a1.set_title("(a) side by side: each coil reaches only its partner's disk\nself-excites at 1206 rad/s, rim 181 m/s — above the 157 m/s copper limit",
             fontsize=9.5, color=INK, loc="left")

# (b) stacked: disk 2 at z = 0 inside coil 1 (I1), disk 1 at z = GAP inside coil 2 (I2).
field_panel(a2, np.linspace(-0.36, 0.36, 180), np.linspace(-0.3, 0.38, 170), [(0, 0, 1.0), (0, GAP, 1.0)])
draw_unit(a2, 0, GAP, ORANGE, "coil 2: I₂", "disk 1", TEAL, True)
draw_unit(a2, 0, 0, TEAL, "coil 1: I₁", "disk 2", ORANGE, False)
a2.annotate("coil 1's stray field\nthreads disk 1 too:\nfeedback K = 0.52,\nno extra copper", xy=(-0.09, GAP - 0.004), xytext=(-0.355, 0.25),
            fontsize=7.5, color=TEAL, arrowprops=dict(arrowstyle="-|>", color=TEAL, lw=1.0), zorder=9, bbox=BOX)
# cross wiring, outside the coils: disk 1 (top) → coil 1 (bottom), disk 2 (bottom) → coil 2 (top)
a2.add_patch(FancyArrowPatch((RO, GAP), (RC + 0.012, -0.012), connectionstyle="arc3,rad=-0.9", arrowstyle="-|>", mutation_scale=12, color=TEAL, lw=2.2, zorder=7))
a2.add_patch(FancyArrowPatch((-RO, 0.0), (-RC - 0.012, GAP + 0.012), connectionstyle="arc3,rad=-0.9", arrowstyle="-|>", mutation_scale=12, color=ORANGE, lw=2.2, zorder=7))
a2.text(0.0, -0.285, "coils couple to each other: k = 0.39", ha="center", fontsize=7.5, color=QUIET, zorder=8, bbox=BOX)
a2.set_title("(b) stacked on one axis, 8 cm apart\nself-excites at 794 rad/s, rim 119 m/s — below the limit",
             fontsize=9.5, color=INK, loc="left")

for ax, xl in ((a1, 0.72), (a2, 0.36)):
    ax.set_xlim(-xl, xl)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
a1.set_ylim(-0.3, 0.38)
a2.set_ylim(-0.3, 0.38)
fig.suptitle("Dual-motor homopolar machine, cross-excited: field lines of the coils (blue, darker = stronger), copper disks, ⊙/× coil turns",
             fontsize=10, color=INK, x=0.01, ha="left")
fig.tight_layout()
fig.savefig(os.path.join(out, "dual_motor_layouts.png"), facecolor=fig.get_facecolor())
print("wrote dual_motor_layouts.png")
