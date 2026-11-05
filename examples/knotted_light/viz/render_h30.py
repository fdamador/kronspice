"""H30 chart from the run's own output (run `go run ./examples/knotted_light > examples/knotted_light/run_output.txt` first).

    python3 examples/knotted_light/viz/render_h30.py

Writes viz/h30_cycle.png: top, part 1 (one implosion then explosion) — energy within
2 L0 from the solver and the exact solution, and the energy the Huygens surface has put
in and taken out; bottom, part 2 (the pulsing knot) — the same, three cycles.
"""
import os
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
BLUE, BLUE_L, ORANGE, GREEN = "#2a78d6", "#86b6ef", "#c8601a", "#2f8f5b"
plt.rcParams.update({"font.size": 13, "figure.facecolor": "#fcfcfb", "savefig.facecolor": "#fcfcfb",
                     "axes.facecolor": "#fcfcfb", "text.color": INK, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.edgecolor": GRID})

text = open(os.path.join(HERE, "..", "run_output.txt")).read()
row = re.compile(r"^\s+(-?\d+)\s+\|\s+([\d.]+) / ([\d.]+)\s+\|\s+([\d.]+)\s+\|\s+([\d.]+) / ([\d.]+)\s+\|", re.M)


def part(name, nxt):
    body = text.split(name, 1)[1]
    if nxt:
        body = body.split(nxt, 1)[0]
    return [tuple(float(g) for g in m.groups()) for m in row.finditer(body)]


parts = [("Part 1: one knot, tightest at t = 0 — fed in from t = −10, out by t = 10", part("Part 1:", "Part 2:")),
         ("Part 2: a pulsing knot — tightest at t = 0, 8, 16 L0/c, the first starts inside, the next two come in through the boundary", part("Part 2:", None))]

fig, axes = plt.subplots(2, 1, figsize=(12, 10.5))
for ax, (title, rows) in zip(axes, parts):
    t = [r[0] for r in rows]
    ax.plot(t, [r[2] for r in rows], color=BLUE_L, lw=5, alpha=0.8, label="within 2 L0, exact")
    ax.plot(t, [r[1] for r in rows], color=BLUE, lw=2, marker="o", ms=4, label="within 2 L0, solver")
    ax.plot(t, [r[4] for r in rows], color=GREEN, lw=2, ls="--", label="put in by the surface (cumulative)")
    ax.plot(t, [r[5] for r in rows], color=ORANGE, lw=2, ls="--", label="taken out by the surface (cumulative)")
    ax.set_title(title, fontsize=13, loc="left", color=INK)
    ax.set_xlabel("time (L0/c)")
    ax.set_ylabel("energy / one knot's energy")
    ax.grid(color=GRID, lw=0.8)
    ax.axvline(0, color=INK2, lw=0.6)
    ax.legend(frameon=False, fontsize=11, loc="upper left")
fig.suptitle("Knotted light, contracting then expanding: every joule enters and leaves through the boundary",
             fontsize=15, y=0.995)
fig.tight_layout()
out = os.path.join(HERE, "h30_cycle.png")
fig.savefig(out, dpi=110)
print("wrote", out)
