"""B(3) transformer: drawing of the N-turn (spring) core and its two windings.

python3 examples/b3_transformer/viz/core.py OUT.png

The core is a closed magnetic loop: a helix of N = 10 turns whose two ends are
joined by an outside return leg. Flux Phi runs along the core wire, so it goes
round the spring axis N times. A winding wrapped round the spring wall (each
turn passing through the bore) encloses n_c passes of the core wire and links
n_c * Phi per winding turn: primary round all 10 passes, secondary round 2.
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, FancyBboxPatch

CORE, PRI, SEC, RET = "#b08d57", "#3b6ea5", "#c8553d", "#8a6d3b"
N, R, P = 10, 1.0, 1.2           # core turns, helix radius, pitch (design: Ø20 mm, 12 mm pitch)
out = sys.argv[1]

fig = plt.figure(figsize=(15, 11))
fig.suptitle("B(3) transformer — 10-turn spring core (Ø20 mm, Ø4 mm core wire, 12 mm pitch, return leg 16× section): primary round 10 passes, secondary round 2",
             fontsize=14, y=0.975)

# ---------------- 3-D view ----------------
ax = fig.add_axes([0.0, 0.03, 0.52, 0.9], projection="3d")
th = np.linspace(np.pi / 2, np.pi / 2 + 2 * np.pi * N, 2000)
x, y, z = R * np.cos(th), R * np.sin(th), P * (th - np.pi / 2) / (2 * np.pi)
ax.plot(x, y, z, color=CORE, lw=10, solid_capstyle="round", zorder=2)
# return leg: from the top end out, down, and back in to the start (at phi = 90 deg)
L = 0.75
zt = z[-1]
ret = np.array([[0, R, zt], [0, R + L, zt], [0, R + L, 0], [0, R, 0]])
ax.plot(ret[:, 0], ret[:, 1], ret[:, 2], color=RET, lw=18, solid_capstyle="round")
# flux arrows along the core
for t in np.pi / 2 + 2 * np.pi * np.array([0.6, 2.6, 4.6, 6.6, 8.6]):
    p0 = np.array([R * np.cos(t), R * np.sin(t), P * (t - np.pi / 2) / (2 * np.pi)])
    d = np.array([-np.sin(t), np.cos(t), P / (2 * np.pi * R)])
    ax.quiver(*p0, *(0.45 * d / np.linalg.norm(d)), color="k", lw=1.6, arrow_length_ratio=0.5)
ax.quiver(0, R + L, zt * 0.6, 0, 0, -0.5, color="k", lw=1.6, arrow_length_ratio=0.5)

def winding(xc, z0, z1, color, nturns, half=0.22, yspan=0.5):
    """Rectangular turns in the x-z plane round the wall at x = xc, stacked along y."""
    for yy in np.linspace(-yspan / 2, yspan / 2, nturns):
        xs = [xc - half, xc + half, xc + half, xc - half, xc - half]
        zs = [z0, z0, z1, z1, z0]
        ax.plot(xs, [yy] * 5, zs, color=color, lw=1.4)

# passes of the core through phi = 0 (x = +R) sit at z = P*(k - 1/4), k = 1..10
zp0 = P * (np.arange(1, N + 1) - 0.25)
winding(R, zp0[0] - 0.35, zp0[-1] + 0.35, PRI, 9, half=0.32)
# passes through phi = 180 deg (x = -R) sit at z = P*(k + 1/4), k = 0..9; secondary round k = 4, 5
zpp = P * (np.arange(0, N) + 0.25)
winding(-R, zpp[4] - 0.35, zpp[5] + 0.35, SEC, 5, half=0.32)
ax.text(R + 0.35, -0.2, zt + 0.15, "primary\n(round 10 core passes)", color=PRI, fontsize=10)
ax.text(-R - 0.3, -0.2, zpp[5] + 0.55, "secondary\n(round 2 core passes)", color=SEC, fontsize=10, ha="right")
ax.text(0.15, R + L, zt * 0.8, "return leg (behind)", color=RET, fontsize=10)
ax.text(0.1, 0.1, -0.55, "flux Φ runs along the core wire (arrows): 10 times round the axis, back down the return leg",
        fontsize=9, ha="center")
ax.set_box_aspect((1, 1.15, 2.6))
lim = 1.7
ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim + 0.3); ax.set_zlim(-0.4, zt + 0.4)
ax.view_init(elev=16, azim=-78)
ax.set_axis_off()

# ---------------- section through the axis (plane y = 0) ----------------
ax2 = fig.add_axes([0.55, 0.08, 0.43, 0.82])
ax2.set_aspect("equal"); ax2.axis("off")
ax2.set_xlim(-2.7, 2.7); ax2.set_ylim(-0.75, zt + 0.75)
ax2.set_title("Section through the spring axis (looking along −y)", fontsize=11)
rw = 0.2
for zz in zp0:   # right wall: flux toward the viewer? helix tangent at phi=0 is +y -> into page
    ax2.add_patch(Circle((R, zz), rw, color=CORE, zorder=3))
    ax2.plot([R - 0.07, R + 0.07], [zz - 0.07, zz + 0.07], color="k", lw=1.2, zorder=4)
    ax2.plot([R - 0.07, R + 0.07], [zz + 0.07, zz - 0.07], color="k", lw=1.2, zorder=4)
for zz in zpp:   # left wall: tangent at phi=180 is -y -> out of page
    ax2.add_patch(Circle((-R, zz), rw, color=CORE, zorder=3))
    ax2.add_patch(Circle((-R, zz), 0.035, color="k", zorder=4))
ax2.plot([0, 0], [-0.5, zt + 0.5], color="0.6", ls="-.", lw=1)
ax2.text(0.05, zt + 0.55, "spring axis", color="0.45", fontsize=9)
ax2.add_patch(FancyBboxPatch((R - 0.32, zp0[0] - 0.35), 0.64, zp0[-1] - zp0[0] + 0.7,
                             boxstyle="round,pad=0.02", fill=False, ec=PRI, lw=2.2, zorder=5))
ax2.add_patch(FancyBboxPatch((-R - 0.32, zpp[4] - 0.35), 0.64, zpp[5] - zpp[4] + 0.7,
                             boxstyle="round,pad=0.02", fill=False, ec=SEC, lw=2.2, zorder=5))
ax2.text(R + 0.38, zt / 2 + 0.25, "primary turn\nencloses n_p = 10\ncore passes\n\nλ_p = N_p · 10 · Φ",
         color=PRI, fontsize=10, va="center")
ax2.text(-R - 0.38, zpp[5] + 0.1, "secondary turn\nencloses n_s = 2\ncore passes\n\nλ_s = N_s · 2 · Φ",
         color=SEC, fontsize=10, va="center", ha="right")
ax2.text(0, -0.62, "⊗ flux into page (front wall)    ⊙ flux out of page (back wall)\n"
         "one flux Φ in every pass, so  v_s / v_p = (N_s · n_s) / (N_p · n_p) = (N_s · 2) / (N_p · 10)",
         ha="center", fontsize=9.5)

fig.savefig(out, dpi=130)
print("wrote", out)
