"""H37 drawing: the magnetostrictive transformer set-ups (rod, ring, helix).

python3 examples/magnetostrictive_transformer/viz/setup.py OUT.png
"""
import math
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Rectangle, Wedge

CORE, PRI, SEC, BIAS, TX = "#b08d57", "#3b6ea5", "#c8553d", "0.45", "#2a9d8f"
out = sys.argv[1]
fig = plt.figure(figsize=(14, 9.2))
fig.suptitle("Magnetostrictive transformer with two ultrasonic transducers — the three core shapes", fontsize=14, y=0.985)


def label(ax, x, y, s, **k):
    k.setdefault("fontsize", 9)
    ax.text(x, y, s, **k)


def transducer(ax, x, y, w, h, name, rot=0):
    ax.add_patch(Rectangle((x - w / 2, y - h / 2), w, h, color=TX, zorder=5))
    for i in range(1, 4):
        if rot == 0:
            ax.plot([x - w / 2 + i * w / 4] * 2, [y - h / 2, y + h / 2], color="white", lw=0.8, zorder=6)
        else:
            ax.plot([x - w / 2, x + w / 2], [y - h / 2 + i * h / 4] * 2, color="white", lw=0.8, zorder=6)
    label(ax, x, y - h / 2 - 0.08, name, ha="center", va="top", color=TX, fontsize=8.5)


# ---- Rod ----
ax = fig.add_axes([0.02, 0.38, 0.32, 0.55])
ax.set_xlim(-1.6, 1.6)
ax.set_ylim(-1.2, 1.2)
ax.axis("off")
fig.text(0.18, 0.93, "1. Rod (Terfenol-D, 50 mm × Ø10 mm)", ha="center", fontsize=11)
ax.add_patch(Rectangle((-1.0, -0.12), 2.0, 0.24, color=CORE, zorder=3))
for i in range(22):
    x = -0.95 + i * 0.09
    ax.add_patch(Rectangle((x, -0.2), 0.035, 0.4, color=PRI if i % 2 == 0 else SEC, zorder=4))
ax.add_patch(Rectangle((-1.05, -0.38), 2.1, 0.76, fill=False, ls="--", color=BIAS, lw=1.3, zorder=2))
ax.add_patch(Rectangle((-1.45, -0.95), 2.9, 0.18, color="0.82", zorder=1))
ax.add_patch(Rectangle((-1.45, -0.95), 0.18, 0.95, color="0.82", zorder=1))
ax.add_patch(Rectangle((1.27, -0.95), 0.18, 0.95, color="0.82", zorder=1))
transducer(ax, -1.17, 0, 0.18, 0.4, "A (fa)", rot=1)
transducer(ax, 1.17, 0, 0.18, 0.4, "B (fb)", rot=1)
label(ax, 0, 0.55, "primary (blue) and secondary (red) interleaved along the rod", ha="center")
label(ax, 0, -0.62, "bias winding (DC), dashed", ha="center", color=BIAS)
label(ax, 0, -1.1, "return yoke + prestress clamp (grey)", ha="center", color="0.4")
ax.annotate("", (0.85, 0.85), (-0.85, 0.85), arrowprops=dict(arrowstyle="<->", color="k"))
label(ax, 0, 0.92, "standing wave: modes 1 and 3 (18, 54 kHz)", ha="center")
xs = [-1 + 2 * i / 100 for i in range(101)]
ax.plot(xs, [0.75 + 0.06 * math.cos(math.pi * (x + 1) / 2) for x in xs], color=TX, lw=1)
ax.plot(xs, [0.75 + 0.06 * math.cos(3 * math.pi * (x + 1) / 2) for x in xs], color=TX, lw=1, ls=":")

# ---- Ring ----
ax = fig.add_axes([0.35, 0.38, 0.30, 0.55])
ax.set_xlim(-1.5, 1.5)
ax.set_ylim(-1.75, 1.2)
ax.set_aspect("equal")
ax.axis("off")
fig.text(0.50, 0.93, "2. Ring (Terfenol-D segments, or Galfenol), Ø40 mm", ha="center", fontsize=11)
ax.add_patch(Wedge((0, 0), 0.82, 0, 360, width=0.18, color=CORE, zorder=3))
n = 16
for k in range(n):  # polygon joints of bonded Terfenol segments
    a = 2 * math.pi * k / n
    ax.plot([0.64 * math.cos(a), 0.82 * math.cos(a)], [0.64 * math.sin(a), 0.82 * math.sin(a)], color="#7a6038", lw=0.6, zorder=4)
for k in range(48):
    a = 2 * math.pi * (k + 0.5) / 48
    ax.plot([0.58 * math.cos(a), 0.88 * math.cos(a)], [0.58 * math.sin(a), 0.88 * math.sin(a)],
            color=PRI if k % 2 == 0 else SEC, lw=2.2, zorder=5)
ax.add_patch(Circle((0, 0), 1.0, fill=False, ls="--", color=BIAS, lw=1.3))
transducer(ax, -1.08, 0, 0.16, 0.36, "A (fa)", rot=1)
transducer(ax, 1.08, 0, 0.16, 0.36, "B (fb)", rot=1)
label(ax, 0, 0.05, "closed flux path\n(no yoke)", ha="center")
label(ax, 0, -0.3, "radial transducers drive\nthe breathing mode", ha="center", color=TX, fontsize=8.5)
label(ax, 0, -1.62, "windings all the way round (blue / red);\ncircumferential waves have no net strain:\nthe coils feel only the breathing mode", ha="center", fontsize=8.5)

# ---- Helix ----
ax = fig.add_axes([0.67, 0.38, 0.32, 0.55])
ax.set_xlim(-1.6, 1.6)
ax.set_ylim(-1.2, 1.2)
ax.axis("off")
fig.text(0.83, 0.93, "3. Helix spring (Galfenol wire, or Terfenol composite)", ha="center", fontsize=11)
N = 8
ts = [i / 400 for i in range(401)]
xs = [-1.1 + 2.2 * t for t in ts]
ys = [0.32 * math.sin(2 * math.pi * N * t) for t in ts]
zs = [math.cos(2 * math.pi * N * t) for t in ts]
for i in range(len(ts) - 1):
    ax.plot(xs[i:i + 2], ys[i:i + 2], color=CORE if zs[i] > 0 else "#d9c3a0", lw=4.5, zorder=3 if zs[i] > 0 else 1, solid_capstyle="round")
for i in range(0, len(ts), 6):
    if zs[i] > 0:
        ax.plot(xs[i], ys[i], "o", ms=3.2, color=PRI if (i // 6) % 2 == 0 else SEC, zorder=4)
transducer(ax, -1.28, 0, 0.16, 0.36, "A (fa)", rot=1)
transducer(ax, 1.28, 0, 0.16, 0.36, "B (fb)", rot=1)
ax.add_patch(Rectangle((-1.15, -0.5), 2.3, 1.0, fill=False, ls="--", color=BIAS, lw=1.3))
label(ax, 0, 0.62, "20 turns, Ø20 mm coil, Ø2 mm wire (1.26 m of wire)", ha="center")
label(ax, 0, -0.68, "primary and secondary wound along the wire (dots);\ntransducers launch waves along the wire from its ends", ha="center", fontsize=8.5)
label(ax, 0, -1.02, "bias field along the wire (DC), dashed", ha="center", color=BIAS, fontsize=8.5)

# ---- Electrical / acoustic connections ----
ax = fig.add_axes([0.04, 0.02, 0.92, 0.32])
ax.set_xlim(0, 10)
ax.set_ylim(0, 3)
ax.axis("off")


def box(x, y, w, h, s, c):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.03", fc="white", ec=c, lw=1.6))
    label(ax, x + w / 2, y + h / 2, s, ha="center", va="center", color=c, fontsize=9.5)


def wire(x1, y1, x2, y2, c="k"):
    ax.plot([x1, x2], [y1, y2], color=c, lw=1.2)


box(4.0, 1.1, 2.0, 0.9, "core\n(rod / ring / helix)", "#7a6038")
box(0.3, 2.0, 1.6, 0.7, "AC source 20 V, 5 kHz\n+ 2 Ω (primary port)", PRI)
box(8.1, 2.0, 1.6, 0.7, "load 20 Ω\n(secondary port)", SEC)
box(0.3, 0.2, 1.6, 0.7, "generator fa\n→ transducer A", TX)
box(8.1, 0.2, 1.6, 0.7, "generator fb\n→ transducer B", TX)
box(4.2, 2.35, 1.6, 0.55, "DC supply → bias coil", BIAS)
box(2.3, 2.0, 1.4, 0.7, "primary\n200 turns", PRI)
box(6.3, 2.0, 1.4, 0.7, "secondary\n400 turns", SEC)
wire(1.9, 2.35, 2.3, 2.35, PRI)
wire(3.7, 2.35, 4.0, 1.75, PRI)
wire(6.0, 1.75, 6.3, 2.35, SEC)
wire(7.7, 2.35, 8.1, 2.35, SEC)
wire(5.0, 2.35, 5.0, 2.0, BIAS)
wire(1.9, 0.55, 4.0, 1.25, TX)
wire(8.1, 0.55, 6.0, 1.25, TX)
label(ax, 5.0, 0.55, "spectrum analyser on the secondary: fe (transformer), fa, fb (ultrasound read by the coil),\n"
      "fa ± fe, fb ± fe, fa + fb, fb − fa (the products: the nonlinear part).  Swap test: move the source to the\n"
      "secondary port and the analyser to the primary.", ha="center", va="center", fontsize=9)
fig.savefig(out, dpi=120)
