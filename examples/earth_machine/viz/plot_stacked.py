# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Figure 10: the stacked design (stacked.go).

(a) The two eddies on one axis, side view, each inside the coil fed by the
    other; blue lines = the coils' magnetic field (exact circular-loop field,
    complete elliptic integrals, equal currents), geometry in disk radii as in
    examples/homopolar: coils 1.2 r_o, 0.4 r_o long, gap 0.53 r_o.
(b) I1(t): side by side it reverses; stacked it locks to one polarity.
(c) the speed difference A = ω1 − ω2 that Rikitake's chaos needs: conserved
    side by side, drained by the stray self-coupling when stacked.

Usage: python viz/plot_stacked.py [data_dir] [out_dir]   (matplotlib, numpy)
"""
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.patches import Ellipse, FancyArrowPatch, Rectangle

here = os.path.dirname(os.path.abspath(__file__))
data = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "data")
out = sys.argv[2] if len(sys.argv) > 2 else here
RC, LC, N, GAP = 1.2, 0.4, 4, 0.08 / 0.15
TEAL, ORANGE, INK, QUIET = "#127b8a", "#d9772b", "#1b1b1f", "#5b5b66"
BOX = dict(boxstyle="round,pad=0.25", facecolor="#fbfaf7", edgecolor="none", alpha=0.92)


def ellipke(m):
    a, b, c = np.ones_like(m), np.sqrt(1 - m), np.sqrt(m)
    s, p = 0.5 * c**2, 0.5
    for _ in range(30):
        a, b, c = 0.5 * (a + b), np.sqrt(a * b), 0.5 * (a - b)
        p *= 2
        s += p * c**2 * 0.5
    K = np.pi / (2 * a)
    return K, K * (1 - s)


def loop_field(r, z, a, z0):
    zz = z - z0
    rr = np.maximum(np.abs(r), 1e-6)
    m = np.clip(4 * a * rr / ((a + rr) ** 2 + zz**2), 0, 1 - 1e-12)
    K, E = ellipke(m)
    den = np.sqrt((a + rr) ** 2 + zz**2)
    q = (a - rr) ** 2 + zz**2
    bz = (K + (a * a - rr * rr - zz * zz) / q * E) / den
    br = zz / (rr * den) * (-K + (a * a + rr * rr + zz * zz) / q * E)
    return np.sign(r) * br, bz


def coil(x, z, zc):
    bx, bz = np.zeros_like(x), np.zeros_like(x)
    for k in range(N):
        a, b = loop_field(x, z, RC, zc + LC * ((k + 0.5) / N - 0.5))
        bx += a
        bz += b
    return bx, bz


def load(name):
    with open(os.path.join(data, name), newline="") as f:
        r = csv.reader(f)
        head = next(r)
        rows = np.array([[float(v) for v in row] for row in r])
    return {h: rows[:, i] for i, h in enumerate(head)}


fig = plt.figure(figsize=(12, 6.4), dpi=300)
fig.patch.set_facecolor("#fbfaf7")
gs = fig.add_gridspec(2, 2, width_ratios=[1, 1.25], height_ratios=[1.3, 1])
ax = fig.add_subplot(gs[:, 0])
bt = fig.add_subplot(gs[0, 1])
ba = fig.add_subplot(gs[1, 1], sharex=bt)
for a in (ax, bt, ba):
    a.set_facecolor("#fbfaf7")

# (a) the stacked eddies.
x, z = np.meshgrid(np.linspace(-2.3, 2.3, 260), np.linspace(-1.5, 2.05, 220))
b1x, b1z = coil(x, z, 0.0)   # coil 1 (I1) around eddy 2 at z = 0
b2x, b2z = coil(x, z, GAP)   # coil 2 (I2) around eddy 1 at z = GAP
bx, bz = b1x + b2x, b1z + b2z
mag = np.hypot(bx, bz)
ax.streamplot(x, z, bx, bz, color=mag, cmap="Blues", norm=LogNorm(mag.min() * 3, mag.max() / 4), density=1.35,
              linewidth=0.7, arrowsize=0.6, zorder=1)
for zc, col, eddy, lab in ((GAP, TEAL, "eddy 1", "coil 2: I₂"), (0.0, ORANGE, "eddy 2", "coil 1: I₁")):
    coil_col = ORANGE if lab.startswith("coil 2") else TEAL
    # the eddy: a rotating slab of fluid, swirl rings in perspective
    for k, w in enumerate((2.0, 1.55, 1.1, 0.65)):
        ax.add_patch(Ellipse((0, zc), w, w * 0.16, facecolor=col if k == 0 else "none", alpha=0.18 if k == 0 else 1,
                             edgecolor=col, lw=1.3 if k else 0.8, zorder=4))
    ax.add_patch(FancyArrowPatch((0.55, zc - 0.075), (-0.25, zc - 0.085), connectionstyle="arc3,rad=-0.08",
                                 arrowstyle="-|>", mutation_scale=11, color=col, lw=1.3, zorder=5))
    ax.text(-0.95, zc + 0.12, eddy, ha="center", fontsize=8.5, fontweight="bold", color=INK, zorder=8, bbox=BOX)
    for s in (-1, 1):
        ax.add_patch(Rectangle((s * RC - 0.06, zc - LC / 2), 0.12, LC, facecolor=coil_col, edgecolor=INK, lw=0.6, zorder=5))
    ax.text(RC + 0.1, zc, lab, ha="left", va="center", fontsize=7.5, color=coil_col, fontweight="bold", zorder=8, bbox=BOX)
# cross wiring: eddy 1's current (I1, teal) to coil 1 below; eddy 2's (I2, orange) to coil 2 above
ax.add_patch(FancyArrowPatch((1.0, GAP), (RC + 0.07, -0.12), connectionstyle="arc3,rad=-1.1", arrowstyle="-|>",
                             mutation_scale=12, color=TEAL, lw=2.2, zorder=7))
ax.add_patch(FancyArrowPatch((-1.0, 0.0), (-RC - 0.07, GAP + 0.12), connectionstyle="arc3,rad=-1.1", arrowstyle="-|>",
                             mutation_scale=12, color=ORANGE, lw=2.2, zorder=7))
ax.annotate("coil 1's stray field also\nthreads eddy 1: KS = 0.52\n(and coil 2's threads eddy 2)", xy=(0.35, GAP + 0.03),
            xytext=(-2.2, 1.55), fontsize=7.5, color=TEAL, arrowprops=dict(arrowstyle="-|>", color=TEAL, lw=1), zorder=9, bbox=BOX)
ax.text(0, -1.38, "rotation axis · coils couple to each other, KM = 0.39", ha="center", fontsize=7.5, color=QUIET, zorder=8, bbox=BOX)
ax.axvline(0, color=QUIET, lw=0.6, ls=(0, (4, 4)), zorder=2)
ax.set_xlim(-2.3, 2.3)
ax.set_ylim(-1.5, 2.05)
ax.set_aspect("equal")
ax.axis("off")
ax.set_title("(a) two eddies stacked on one axis, each in the other's coil", fontsize=9.5, color=INK, loc="left")

# (b), (c) the dynamics.
d = load("stacked_traces.csv")
bt.plot(d["t"], d["I1_side"], color=INK, lw=0.7, label="side by side: reverses")
bt.plot(d["t"], d["I1_near"], color=TEAL, lw=1.4, label="stacked, gap 0.53 r_o: locked")
bt.plot(d["t"], d["I1_far"], color=ORANGE, lw=1.2, ls="--", label="stacked, gap 3 r_o: locks once A drains")
bt.axhline(0, color=QUIET, lw=0.5)
bt.set_ylabel("I₁ (field of eddy 1's loop)")
bt.legend(fontsize=7, loc="lower right", frameon=False)
bt.set_title("(b) the field: side by side it reverses, stacked it holds one polarity", fontsize=9.5, color=INK, loc="left")
ba.plot(d["t"], d["A_side"], color=INK, lw=1.0, label="side by side: conserved")
ba.plot(d["t"], d["A_near"], color=TEAL, lw=1.4, label="stacked, 0.53 r_o")
ba.plot(d["t"], d["A_far"], color=ORANGE, lw=1.2, ls="--", label="stacked, 3 r_o")
ba.set_ylabel("A = ω₁ − ω₂")
ba.set_xlabel("Rikitake time (μ = 1, equal drives)")
ba.legend(fontsize=7, loc="upper right", frameon=False)
ba.set_title("(c) the speed difference the chaos needs drains when stacked", fontsize=9.5, color=INK, loc="left")
for a in (bt, ba):
    a.grid(alpha=0.3)
    a.spines["top"].set_visible(False)
    a.spines["right"].set_visible(False)
bt.set_xlim(0, 400)
fig.tight_layout()
fig.savefig(os.path.join(out, "fig10_stacked.png"), facecolor=fig.get_facecolor())
print("wrote fig10_stacked.png")
