# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Charts for Parts 5 and 6 of examples/fdamador_neutrino.

    go run ./examples/fdamador_neutrino -field -viz examples/fdamador_neutrino/viz
    python3 examples/fdamador_neutrino/viz/render_part56.py

Writes
  viz/implosion.png   Part 5: the E-ring imploded by the time-reversed ring-down and
                      exploding again, as a circuit (5a: drives compared, the pulsing
                      ring's books) and in the field solver (5b: energy inside the
                      Huygens surface, the field pattern at the peak)
  viz/flavours.png    Part 6 (H39): the three wave types on the ring line, the
                      two-state beat a defect makes, the three-state band-edge mixing
                      against the 3x3 formula, and implosion / explosion of each type

Every curve is read from the CSV files the Go program writes (-viz); nothing is
computed here except plotting.
"""
import csv
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
BLUE, ORANGE, GREEN, RED, PURPLE = "#2a78d6", "#eb6834", "#2f9e58", "#c93c3c", "#7d4fc4"
SEQ = ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]
plt.rcParams.update({"font.size": 13, "axes.titlesize": 14.5, "axes.titleweight": "normal",
                     "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb",
                     "savefig.facecolor": "#fcfcfb", "axes.spines.top": False, "axes.spines.right": False,
                     "legend.frameon": False, "legend.fontsize": 11})


def load(name):
    with open(os.path.join(HERE, name)) as f:
        rows = list(csv.reader(f))[1:]
    return np.array([[float(x) for x in r] for r in rows])


def exists(name):
    return os.path.exists(os.path.join(HERE, name))


def grid(ax):
    ax.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)


# ================================================================ implosion.png
runs = {}
with open(os.path.join(HERE, "p5a_runs.csv"), encoding="utf-8") as f:
    for line in f:
        i, rest = line.strip().split(",", 1)
        name, ab, pk = rest.rsplit(",", 2)
        runs[int(i)] = (name, float(ab), float(pk))

have5b = exists("p5b_plane.csv")
fig = plt.figure(figsize=(14, 21 if have5b else 13))
gs = fig.add_gridspec(4 if have5b else 2, 2, hspace=0.42, wspace=0.26, left=0.07, right=0.97, top=0.95, bottom=0.04)
fig.suptitle("fdamador neutrino, Part 5: implosion and explosion of the E-ring", fontsize=18, color=INK, y=0.985)

ax = fig.add_subplot(gs[0, 0])
sel = [(0, BLUE, "time-reversed ring-down (CPA rate)"), (5, ORANGE, "constant burst, 2 τ_E"),
       (6, GREEN, "constant burst, 8 τ_E"), (4, RED, "ring-down not reversed (decaying)")]
for i, col, lab in sel:
    d = load(f"p5a_run{i}.csv")
    ax.plot(d[:, 0], d[:, 1], color=col, lw=2, label=f"{lab}: peak {runs[i][2]:.3f}")
ax.set_xlabel("time (ring periods)")
ax.set_ylabel("stored energy / incident energy")
ax.set_title("5a. Circuit: same incident energy, different drives")
ax.legend(loc="upper left", fontsize=10)
grid(ax)

ax = fig.add_subplot(gs[0, 1])
d = load("p5a_run0.csv")
ax.plot(d[:, 0], d[:, 2], color=BLUE, lw=1.6, label="incident power (time-reversed drive)")
ax.plot(d[:, 0], d[:, 3], color=ORANGE, lw=1.6, label="outgoing power: ~0 while driven, the explosion after")
ax.set_xlabel("time (ring periods)")
ax.set_ylabel("power / peak incident power")
ax.set_title(f"5a. Coherent perfect absorption: absorbed {runs[0][1]:.5f}")
ax.legend(loc="upper left", fontsize=10)
grid(ax)

ax = fig.add_subplot(gs[1, 0])
for i, col, lab in [(5, ORANGE, "burst 2 τ_E"), (4, RED, "decaying drive")]:
    d = load(f"p5a_run{i}.csv")
    ax.plot(d[:, 0], d[:, 3], color=col, lw=1.4, label=f"{lab}: outgoing (absorbed {runs[i][1]:.3f})")
ax.set_xlabel("time (ring periods)")
ax.set_ylabel("outgoing power / peak incident")
ax.set_title("5a. Other drives are reflected while they arrive")
ax.legend(loc="upper right", fontsize=10)
grid(ax)

ax = fig.add_subplot(gs[1, 1])
tl = load("p5a_train_locked.csv")
x = tl[:, 0]
wbar = 0.26
ax.bar(x - wbar, tl[:, 1] * 1e12, wbar, color=BLUE, label="in (incident wave)")
ax.bar(x, tl[:, 2] * 1e12, wbar, color=ORANGE, label="out (reflection + explosion)")
ax.bar(x + wbar, (tl[:, 1] - tl[:, 2]) * 1e12 * 20, wbar, color=GREEN, label="in − out, ×20")
ax.plot(x + wbar, tl[:, 3] * 1e12 * 20, "k_", ms=18, mew=2, label="dielectric heat, ×20")
ax.set_xlabel("cycle (one implosion every 40 periods)")
ax.set_ylabel("energy per cycle (pJ)")
ax.set_title("5a. Pulsing ring: every cycle paid by the incoming wave")
ax.set_ylim(0, 300)
ax.legend(loc="upper center", fontsize=10, ncol=2)
grid(ax)

if have5b:
    ax = fig.add_subplot(gs[2, 0])
    summ = load("p5b_summary.csv")
    for k, c in enumerate(summ[:, 0]):
        col = SEQ[1 + k % 3]
        d = load(f"p5b_implosion_c{c:g}.csv")
        ax.plot(d[:, 0], d[:, 1], color=col, lw=2, label=f"inside the surface, grid r/{c:g}")
        ax.plot(d[:, 0], d[:, 2], color=col, lw=1.2, ls="--", label=f"outside + absorbed, r/{c:g}")
        if k == 0:
            ax.plot(d[:, 0], d[:, 4], color=INK2, lw=1, ls=":", label="declared incident energy so far")
    ax.set_xlabel("time (periods)")
    ax.set_ylabel("energy / incident energy")
    ax.set_title("5b. Field solver: the torus imploded, then exploding")
    ax.legend(loc="upper left", fontsize=9.5)
    grid(ax)

    ax = fig.add_subplot(gs[2, 1])
    for k, c in enumerate(summ[:, 0]):
        col = SEQ[1 + k % 3]
        d = load(f"p5b_burst_c{c:g}.csv")
        ax.plot(d[:, 0], d[:, 1], color=col, lw=2, label=f"burst: inside, r/{c:g}")
        ax.plot(d[:, 0], d[:, 2], color=col, lw=1.2, ls="--", label=f"burst: outside + absorbed, r/{c:g}")
    ax.set_xlabel("time (periods)")
    ax.set_ylabel("energy / incident energy")
    ax.set_title("5b. A constant burst of the same pattern")
    ax.legend(loc="upper left", fontsize=9.5)
    grid(ax)

    ax = fig.add_subplot(gs[3, 0])
    pl = load("p5b_plane.csv")
    keep = pl[:, 0] <= 0.020  # inside the Huygens surface (outside it the solver holds the scattered field only)
    pl = pl[keep]
    xs = np.unique(pl[:, 0])
    zs = np.unique(pl[:, 1])
    ey = pl[:, 2].reshape(len(xs), len(zs)).T
    vmax = np.abs(ey).max()
    im = ax.pcolormesh(xs * 1e3, zs * 1e3, ey, cmap="RdBu_r", vmin=-vmax, vmax=vmax, shading="auto")
    th = np.linspace(0, 2 * np.pi, 200)
    ax.plot(10 + 2.5 * np.cos(th), 2.5 * np.sin(th), color=INK, lw=1.2)
    ax.set_aspect("equal")
    ax.set_xlabel("x (mm), y = 0 plane")
    ax.set_ylabel("z (mm)")
    ax.set_title("5b. At the peak: E_φ fills the tube (grid r/3)")
    fig.colorbar(im, ax=ax, shrink=0.8, label="E_y = E_φ (arb.)")

    ax = fig.add_subplot(gs[3, 1])
    bx = pl[:, 3].reshape(len(xs), len(zs)).T
    bz = pl[:, 4].reshape(len(xs), len(zs)).T
    mag = np.hypot(bx, bz)
    ax.pcolormesh(xs * 1e3, zs * 1e3, mag, cmap="Blues", shading="auto")
    ax.streamplot(xs * 1e3, zs * 1e3, bx, bz, color=INK, density=1.2, linewidth=0.8, arrowsize=0.8)
    ax.plot(10 + 2.5 * np.cos(th), 2.5 * np.sin(th), color=ORANGE, lw=1.6)
    ax.set_aspect("equal")
    ax.set_xlabel("x (mm), y = 0 plane")
    ax.set_ylabel("z (mm)")
    ax.set_title("5b. T/4 later: B through the hole, round the tube")

out = os.path.join(HERE, "implosion.png")
fig.savefig(out, dpi=100)
print("wrote", out)

# ================================================================ flavours.png
fig = plt.figure(figsize=(14, 22))
gs = fig.add_gridspec(4, 3, hspace=0.48, wspace=0.3, left=0.07, right=0.97, top=0.95, bottom=0.04)
fig.suptitle("fdamador neutrino, Part 6 (H39): flavours as wave types on the ring line", fontsize=18, color=INK, y=0.985)

titles = {"C": "constructive (travelling)", "S": "standing (antinode at c)", "D": "destructive (node at c)"}
for col_i, ty in enumerate(["C", "S", "D"]):
    ax = fig.add_subplot(gs[0, col_i])
    d = load(f"p6_type_{ty}.csv")
    for s in range(1, d.shape[1]):
        ax.plot(d[:, 0], d[:, s], color=SEQ[min(s - 1, 3)] if s < 5 else INK, lw=1.8, marker="o", ms=3,
                label=f"t = {s - 1}T/8")
    ax.axvline(0.5, color=ORANGE, lw=1, ls=":")
    ax.set_xlabel("node j (c = ½ dotted)")
    ax.set_ylabel("u_j (V)")
    ax.set_title(titles[ty])
    if col_i == 0:
        ax.legend(fontsize=8.5, loc="upper center", ncol=3, frameon=True, framealpha=0.95, edgecolor=GRID)
        ax.set_ylim(-1.15, 1.6)
    grid(ax)

ax = fig.add_subplot(gs[1, :2])
deltas = load("p6_df.csv")[:, 0]
for i, dl in enumerate(deltas):
    d = load(f"p6_beat{i}.csv")
    ax.plot(d[:, 0], d[:, 1], color=SEQ[i], lw=1.8, label=f"CW start, δ = {dl:g}: P(CCW)")
for j, (s, col) in enumerate([(1, GREEN), (2, PURPLE), (3, RED)]):
    d = load(f"p6_beat{len(deltas) + j}.csv")
    ax.plot(d[:, 0], d[:, 1], color=col, lw=1.4, ls="--", label=f"standing start s = {s} (θ = {15 * s}°)")
ax.set_xlabel("time (beat periods 1/Δf, each run's exact Δf)")
ax.set_ylabel("conversion probability")
ax.set_title("Test 2: a defect splits the pair; P = sin²2θ · sin²(πΔf t)")
ax.legend(fontsize=9.5, loc="upper right", ncol=2)
ax.set_ylim(0, 1.25)
grid(ax)

ax = fig.add_subplot(gs[1, 2])
df = load("p6_df.csv")
ax.plot(df[:, 0], df[:, 3] / 1e3, color=INK2, lw=3, ls="--", label="first order (on the exact line)")
ax.plot(df[:, 0], df[:, 1] / 1e3, color=BLUE, lw=2, label="exact (eigenmodes)")
ax.plot(df[:, 0], df[:, 2] / 1e3, "o", color=ORANGE, ms=8, label="kronspice beat")
ax.set_xlabel("defect δ")
ax.set_ylabel("Δf (kHz)")
ax.set_title("Beat frequency vs defect")
ax.legend(fontsize=10)
grid(ax)

ax = fig.add_subplot(gs[2, :])
th = load("p6_three.csv")
with open(os.path.join(HERE, "p6_three_meta.csv")) as f:
    frm, e_, m_, t_ = [int(x) for x in f.read().strip().split(",")]
lab = {e_: "e", m_: "μ", t_: "τ"}
names = ["f1 (cos, κ = 11π/12)", "f2 (sin, κ = 11π/12)", "f3 (π-mode)"]
for b, col in zip(range(3), [BLUE, ORANGE, GREEN]):
    ax.plot(th[:, 0], th[:, 1 + b], color=col, lw=2.2, alpha=0.75, label=f"kronspice: {lab[b]} = {names[b]}")
    ax.plot(th[:, 0], th[:, 7 + b], color=INK, lw=0.8, ls=":")
    ax.plot(th[:, 0], th[:, 4 + b], color=col, lw=0.9, ls="--")
ax.plot([], [], color=INK, lw=0.8, ls=":", label="all ring modes (exact)")
ax.plot([], [], color=INK2, lw=0.9, ls="--", label="3×3 formula |Σ U_αi U_βi e^{iω_i t}|²")
ax.set_xlabel("time (periods)")
ax.set_ylabel("energy fraction")
ax.set_title(f"Test 2: three ring states mixed by two defects, started in {lab[frm]} (tuned toward NuFIT 6.0)")
ax.legend(fontsize=10, loc="upper right", ncol=2)
ax.set_ylim(0, 1.15)
grid(ax)

ax = fig.add_subplot(gs[3, 0])
for i, (ty, col) in enumerate([("S", BLUE), ("D", ORANGE), ("C", GREEN)]):
    d = load(f"p6_explode{i}.csv")
    ax.plot(d[:, 0], d[:, 1], color=col, lw=2, label=f"{ty}: energy")
d = load("p6_explode2.csv")
ax.plot(d[:, 0], d[:, 2], color=GREEN, lw=1.2, ls="--", label="C: fraction in D")
ax.set_xlabel("time (periods)")
ax.set_ylabel("energy / start, or fraction")
ax.set_title("Test 3: free explosion through one port")
ax.legend(fontsize=9.5, loc="center right")
grid(ax)

ax = fig.add_subplot(gs[3, 1])
for i, (lab_, col) in enumerate([("S, one port", BLUE), ("C, two ports", GREEN), ("C, one port", RED)]):
    d = load(f"p6_rebuild{i}.csv")
    ax.plot(d[:, 0], d[:, 1], color=col, lw=2, label=lab_)
ax.set_xlabel("time (periods)")
ax.set_ylabel("ring energy / incident energy")
ax.set_title("Test 3: rebuilt by the reversed explosion")
ax.legend(fontsize=9.5, loc="upper left")
grid(ax)

ax = fig.add_subplot(gs[3, 2])
d = load("p6_fexplode.csv")
ax.plot(d[:, 0], d[:, 2], color=INK2, lw=1.2, label="free explosion: P(F⊥)")
d = load("p6_rebuild3.csv")
ax.plot(d[:, 0], d[:, 2], color=PURPLE, lw=1.4, label="reversed explosion: P(F⊥)")
ax.plot(d[:, 0], d[:, 1], color=PURPLE, lw=2, ls="--", label="reversed explosion: energy")
d = load("p6_rebuild4.csv")
ax.plot(d[:, 0], d[:, 2], color=ORANGE, lw=1.2, label="single-tone drive: P(F⊥)")
ax.set_xlabel("time (periods)")
ax.set_ylabel("fraction")
ax.set_title("Test 3: a flavour superposition")
ax.legend(fontsize=9, loc="upper right")
ax.set_ylim(0, 1.3)
grid(ax)

out = os.path.join(HERE, "flavours.png")
fig.savefig(out, dpi=100)
print("wrote", out)
