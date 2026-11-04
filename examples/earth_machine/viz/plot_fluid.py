# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Figure 9: the Rikitake two-disk machine drawn as the core sees it.

Each disk is a convective eddy of liquid iron (a Lamb-Oseen vortex, top view,
streamlines coloured by speed); each eddy's current runs out along a ribbon to
the coil wrapped around the OTHER eddy, whose field threads its partner. The
two wirings cross: Rikitake's cross-excitation. An illustration, not a field
solve; the parts map one to one onto the circuit (README, Fig. 0).

Usage: python viz/plot_fluid.py [out_dir]   (matplotlib, numpy)
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BOX = dict(boxstyle="round,pad=0.25", facecolor="#fbfaf7", edgecolor="none", alpha=0.92)
from matplotlib.patches import Circle, FancyArrowPatch, Wedge

here = os.path.dirname(os.path.abspath(__file__))
out = sys.argv[1] if len(sys.argv) > 1 else here

# Two co-rotating eddies (disks), Lamb-Oseen profile.
c1, c2 = np.array([-2.2, 0.0]), np.array([2.2, 0.0])
rc, gamma = 0.9, 6.0
x, y = np.meshgrid(np.linspace(-5.4, 5.4, 540), np.linspace(-3.2, 3.2, 320))
u = np.zeros_like(x)
v = np.zeros_like(x)
for c, g in ((c1, gamma), (c2, 0.8 * gamma)):  # disk 2 slower: A = ω1 − ω2 > 0
    dx, dy = x - c[0], y - c[1]
    r2 = dx * dx + dy * dy + 1e-9
    vt = g / (2 * np.pi * r2) * (1 - np.exp(-r2 / rc**2))
    u += -vt * dy
    v += vt * dx
speed = np.hypot(u, v)

TEAL, ORANGE, COPPER, INK, QUIET = "#127b8a", "#d9772b", "#b0643a", "#1b1b1f", "#5b5b66"
fig, ax = plt.subplots(figsize=(10, 6.2), dpi=300)
fig.patch.set_facecolor("#fbfaf7")
ax.set_facecolor("#fbfaf7")

# The liquid: streamlines, width and colour by speed.
lw = 0.25 + 1.4 * speed / speed.max()
ax.streamplot(x, y, u, v, density=1.25, color=speed, cmap="magma_r", linewidth=lw,
              arrowsize=0.7, broken_streamlines=False, zorder=1)

# Disks: the eddies' cores (where the radial current is collected).
for c, lab in ((c1, "disk 1  (ω₁)"), (c2, "disk 2  (ω₂)")):
    ax.add_patch(Circle(c, 0.95, facecolor="#f1efe9", edgecolor=INK, lw=1.0, alpha=1.0, zorder=3))
    ax.text(c[0], c[1] - 0.15, lab, ha="center", va="center", fontsize=10, fontweight="bold", color=INK, zorder=6)

# Field coils: coil 1 (carries I1) wraps disk 2; coil 2 (carries I2) wraps disk 1.
def coil(center, colour, label, side):
    for k in range(4):  # a few turns
        ax.add_patch(Circle(center, 1.55 + 0.09 * k, fill=False, edgecolor=colour, lw=2.0, alpha=0.95, zorder=4))
    # current direction arrows on the coil
    for ang in (60, 200, 320):
        a = np.radians(ang)
        p = center + 1.69 * np.array([np.cos(a), np.sin(a)])
        t = np.array([-np.sin(a), np.cos(a)]) * 0.35
        ax.add_patch(FancyArrowPatch(p - t / 2, p + t / 2, arrowstyle="-|>", mutation_scale=13, color=colour, lw=0, zorder=5))
    ax.text(center[0] + side * 1.95, center[1] + 1.85, label, ha="center", fontsize=9, color=colour, fontweight="bold", zorder=7, bbox=BOX)


coil(c2, TEAL, "coil 1: carries I₁", +1)
coil(c1, ORANGE, "coil 2: carries I₂", -1)

# Field threading each disk (out of the page): ⊙ marks.
for c, colour, lab in ((c1, ORANGE, "B₂ threads disk 1"), (c2, TEAL, "B₁ threads disk 2")):
    for off in ((-0.45, 0.35), (0.45, 0.35), (0.0, -0.55)):
        p = c + np.array(off)
        ax.add_patch(Circle(p, 0.11, facecolor="none", edgecolor=colour, lw=1.3, zorder=6))
        ax.add_patch(Circle(p, 0.03, facecolor=colour, edgecolor="none", zorder=6))
    ax.text(c[0], c[1] + 0.95 + 0.12, lab, ha="center", va="bottom", fontsize=8.5, color=colour, zorder=7, bbox=BOX)

# Cross wiring: the two ribbons cross between the eddies (figure eight).
t = np.linspace(0, 1, 200)
def ribbon(p0, p1, bulge, colour):
    mid = (p0 + p1) / 2 + np.array([0, bulge])
    pts = np.outer((1 - t) ** 2, p0) + np.outer(2 * (1 - t) * t, mid) + np.outer(t**2, p1)
    ax.plot(pts[:, 0], pts[:, 1], color=colour, lw=4.2, alpha=0.9, solid_capstyle="round", zorder=5)
    k = 120
    ax.add_patch(FancyArrowPatch(pts[k - 6], pts[k], arrowstyle="-|>", mutation_scale=16, color=colour, lw=0, zorder=6))

ribbon(c1 + np.array([0.7, 0.62]), c2 + np.array([-1.45, -0.62]), -1.1, TEAL)   # I1: disk 1 → coil around disk 2
ribbon(c2 + np.array([-0.7, 0.62]), c1 + np.array([1.45, -0.62]), -1.1, ORANGE)  # I2: disk 2 → coil around disk 1
ax.text(0, -1.45, "the wirings cross:\neach eddy excites the other", ha="center", va="top", fontsize=8.5, color=QUIET, zorder=7, bbox=BOX)

# Drive torque from convection.
for c in (c1, c2):
    ax.add_patch(FancyArrowPatch(c + np.array([0.0, -2.35]), c + np.array([1.2, -2.0]),
                                 connectionstyle="arc3,rad=0.35", arrowstyle="-|>", mutation_scale=14,
                                 color=QUIET, lw=1.2, zorder=6))
ax.text(c1[0] + 0.6, -2.75, "drive torque G\n(core convection)", ha="center", va="top", fontsize=8, color=QUIET, zorder=7, bbox=BOX)
ax.text(c2[0] + 0.6, -2.75, "drive torque G\n(core convection)", ha="center", va="top", fontsize=8, color=QUIET, zorder=7, bbox=BOX)

ax.set_title("The Rikitake machine as two eddies in the liquid core: each one's current builds the other's field",
             fontsize=11, color=INK, loc="left", pad=10)
ax.text(-5.35, 3.05, "streamlines: flow speed of the liquid iron (darker = faster); disk 2 turns slower, A = ω₁ − ω₂ > 0",
        fontsize=7.5, color=QUIET, va="top", zorder=7, bbox=BOX)
ax.set_xlim(-5.4, 5.4)
ax.set_ylim(-3.2, 3.2)
ax.set_aspect("equal")
ax.axis("off")
fig.tight_layout()
fig.savefig(os.path.join(out, "fig9_fluid_machine.png"), facecolor=fig.get_facecolor())
print("wrote fig9_fluid_machine.png")
