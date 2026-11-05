"""3-D pictures of the knotted light (run `go run ./examples/knotted_light -viz examples/knotted_light/viz` first).

    python3 examples/knotted_light/viz/render.py

Writes
  viz/knot_3d.png      the 24 magnetic field lines of the hopfion at t = 0: three
                       nests of closed loops, every pair linked once
  viz/knot_time.png    the same family of lines at t = 0, 1, 2 and 3 (L0/c), seeds
                       scaled with the knot's size √(1 + t²), drawn to one scale with
                       a grey sphere at the energy's rms radius: the loops stay
                       closed and linked while the knot grows
Lines are traced on the exact solution (see README); colour marks the nest a
line's seed belongs to, from nearest the axis (light) to nearest the core (dark).
"""
import json
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
INK, INK2 = "#0b0b0b", "#52514e"
NEST = ["#86b6ef", "#2a78d6", "#0d366b"]  # one blue ramp, ordered by nest
plt.rcParams.update({"font.size": 14, "figure.facecolor": "#fcfcfb", "savefig.facecolor": "#fcfcfb",
                     "axes.facecolor": "#fcfcfb", "text.color": INK, "axes.labelcolor": INK2})

with open(os.path.join(HERE, "knot_lines.json")) as f:
    frames = json.load(f)


def sphere(ax, r):
    u, v = np.mgrid[0:2 * np.pi:28j, 0:np.pi:14j]
    ax.plot_wireframe(r * np.cos(u) * np.sin(v), r * np.sin(u) * np.sin(v), r * np.cos(v),
                      color="#b4b2a9", lw=0.4, alpha=0.5)


def draw(ax, frame, lim, lw=1.6, rms=False):
    if rms:
        sphere(ax, frame["rms"])
    for line in frame["lines"]:
        p = np.array(line["points"])
        if len(p) < 2:
            continue
        p = np.vstack([p, p[:1]])  # close the loop
        ax.plot(p[:, 0], p[:, 1], p[:, 2], color=NEST[line["torus"]], lw=lw, alpha=0.9)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_zlim(-lim, lim)
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=24, azim=38)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.set_facecolor("#f1f0ec")
        axis.pane.set_edgecolor("#e4e3df")
        axis.set_tick_params(colors=INK2, labelsize=10)
    ax.set_xlabel("x (L0)", fontsize=11)
    ax.set_ylabel("y (L0)", fontsize=11)
    ax.set_zlabel("z (L0)", fontsize=11)


# 1. The knot at t = 0, large.
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection="3d")
f0 = frames[0]
draw(ax, f0, 3.2, lw=2)
fig.suptitle("Knotted light at t = 0: 24 closed field lines, every pair linked once", fontsize=17, y=0.95)
ax.set_title(f"linking number of the first pair {f0['linking']:.4f}   ·   colour = nest of lines (light: near the axis, dark: near the core)",
             fontsize=11, color=INK2, pad=4)
out = os.path.join(HERE, "knot_3d.png")
fig.savefig(out, dpi=110, bbox_inches="tight")
print("wrote", out)

# 2. Four times, one scale.
fig = plt.figure(figsize=(12, 12.5))
lim = 4.5
for k, fr in enumerate(frames[:4]):
    ax = fig.add_subplot(2, 2, k + 1, projection="3d")
    draw(ax, fr, lim, lw=1.3, rms=True)
    ax.set_title(f"t = {fr['t']:g} L0/c   ·   rms radius {fr['rms']:.2f} L0   ·   linking {fr['linking']:.3f}",
                 fontsize=12, color=INK)
fig.suptitle("The same knot over time: still linked, spreading at the speed of light", fontsize=17, y=0.985)
fig.text(0.5, 0.945, "grey sphere = rms radius of the energy; seeds scaled with the knot's size √(1 + t²)", ha="center", fontsize=12, color=INK2)
fig.subplots_adjust(left=0.02, right=0.98, top=0.90, bottom=0.03, wspace=0.05, hspace=0.12)
out = os.path.join(HERE, "knot_time.png")
fig.savefig(out, dpi=100)
print("wrote", out)
