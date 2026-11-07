"""H38 one-tank figure: concentric pipes.

go run ./examples/al_cell_pump -data DIR
python3 examples/al_cell_pump/viz/plot_pipes.py DIR OUT.png

(a) efficiency against bank voltage: plate stacks and the four kinds of one-tank pipe rectifier;
(b) seals on the best series tank: efficiency and the worst oxide against the leak between annuli;
(c) parallel annuli for current: the best single cell rewired with N parallel annuli, at three loads;
(d) the two tanks: parallel (alternating aluminium and iron on two bus bars, top view) and
    graded series (clad pipes, masked aluminium, sealed bottoms, cutaway).
"""
import csv
import math
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle

d, out = sys.argv[1], sys.argv[2]
rd = lambda n: list(csv.DictReader(open(f"{d}/{n}")))
f = lambda r, k: float(r[k])
plates = rd("pareto.csv")
pipes = rd("pipes_pareto.csv")
try:
    seal = rd("seal_sweep.csv")
except FileNotFoundError:
    seal = []
amps = rd("pipes_amps.csv")

fig, ax = plt.subplots(2, 2, figsize=(13, 10))
fig.suptitle("H38 — one tank per rectifier: concentric aluminium / iron pipes", fontsize=12)

# (a)
a = ax[0][0]
a.scatter([-f(r, "out_V") for r in plates], [f(r, "eff") for r in plates], s=28, color="#bbbbbb", label="plate stacks (separate cells)")
kinds = {"single cell": "#2e86c1", "parallel": "#16a085", "series": "#e67e22", "graded series": "#8e44ad"}
for k, col in kinds.items():
    sel = []
    for r in pipes:
        n, s = int(r["annuli"]), int(r["series"])
        kk = "single cell" if n == 1 else ("parallel" if s == 0 else ("graded series" if r.get("graded", "0") == "1" else "series"))
        if kk == k:
            sel.append(r)
    if sel:
        a.scatter([-f(r, "out_V") for r in sel], [f(r, "eff") for r in sel], s=44, color=col, edgecolor="k", linewidth=0.4, label=f"pipes: {k} ({len(sel)})")
a.set_xscale("log")
a.set_ylim(0, 1)
a.set_xlabel("bank voltage |V| [V]")
a.set_ylabel("catalytic efficiency η")
a.set_title("(a) Pareto fronts")
a.grid(alpha=0.3)
a.legend(fontsize=8, loc="lower left")

# (b)
a = ax[0][1]
if seal:
    R = [f(r, "rseal") for r in seal]
    a.plot(R, [f(r, "eff") for r in seal], "o-", color="#8e44ad", label="η")
    a.plot(R, [f(r, "rev") for r in seal], "s--", color="#c0392b", label="worst oxide / its forming voltage")
    a.axhline(1, color="#c0392b", lw=0.6, ls=":")
    a.axhline(0.9, color="#c0392b", lw=0.6, ls=":")
    a.set_xscale("log")
    a.set_xlabel("leak resistance between neighbouring annuli, Rseal [Ω]")
    a.set_ylim(0, 1.1)
    a.legend(fontsize=8, loc="lower right")
    a.set_title("(b) seals on the best series tank that needs its annuli")
    a.grid(alpha=0.3)
else:
    a.text(0.5, 0.5, "no series tank on the front", ha="center")
    a.axis("off")

# (c)
a = ax[1][0]
loads = sorted({f(r, "rload") for r in amps}, reverse=True)
for L, col in zip(loads, ["#2e86c1", "#16a085", "#e67e22"]):
    sel = [r for r in amps if f(r, "rload") == L]
    n = [int(r["n"]) for r in sel]
    a.plot(n, [f(r, "eff") for r in sel], "o-", color=col, label=f"load {L:.3g} Ω")
    for r in sel:
        a.annotate(f"{1e3*f(r, 'iout_A'):.3g} mA", (int(r["n"]), f(r, "eff")), fontsize=7, xytext=(3, 4), textcoords="offset points")
a.set_xlabel("parallel annuli N")
a.set_ylabel("η")
a.set_ylim(0, 1)
a.set_xticks([1, 2, 4, 6])
a.set_title("(c) parallel annuli for current (best single cell rewired)")
a.text(0.98, 0.04, "the ladder and drive set the current, and the 2 V\nelectrolysis onset (not the water) is the forward loss:\nextra annuli only add film capacitance",
       transform=a.transAxes, ha="right", va="bottom", fontsize=7.5, bbox=dict(fc="white", ec="#cccccc"))
a.legend(fontsize=8)
a.grid(alpha=0.3)

# (d) drawings
a = ax[1][1]
a.set_title("(d) parallel tank (top view)          graded series tank (cutaway)")
AL, FE, W = "#b0b7bd", "#4a4f55", "#d6eaf8"
# parallel top view
cx, cy = 2.2, 2.6
radii = [0.25, 0.75, 1.3]
a.add_patch(Circle((cx, cy), 2.0, color=W))
for i, r in enumerate(reversed(radii)):
    j = len(radii) - 1 - i
    a.add_patch(Circle((cx, cy), r + 0.08, color=AL if j % 2 == 0 else FE))
    a.add_patch(Circle((cx, cy), r - 0.08 if j else 0, color=W))
a.add_patch(Circle((cx, cy), 0.25, color=AL))
a.add_patch(Circle((cx, cy), 2.0, fill=False, ec=FE, lw=3))
a.plot([cx - 2.1, cx + 1.3], [cy + 2.3, cy + 2.3], color=AL, lw=4)
a.plot([cx - 2.1, cx + 2.0], [cy - 2.3, cy - 2.3], color=FE, lw=4)
for r in [0.0, 1.3]:
    a.plot([cx + r, cx + r], [cy, cy + 2.3], color=AL, lw=1.2)
for r in [0.75, 2.0]:
    a.plot([cx + r, cx + r], [cy, cy - 2.3], color=FE, lw=1.2)
a.text(cx - 2.1, cy + 2.45, "aluminium bus", fontsize=8)
a.text(cx - 2.1, cy - 2.65, "iron bus (tank wall is iron)", fontsize=8)
a.text(cx, cy - 3.2, "3 annuli, each Al facing Fe, all across the\nsame two buses: no end seals needed", ha="center", fontsize=7.5)

# graded series cutaway (half section, rod on the left)
x0, y0 = 5.6, 0.6
H = 4.2
a.add_patch(Rectangle((x0, y0), 4.4, H, color=W))
a.add_patch(Rectangle((x0 - 0.1, y0 - 0.3), 4.6, 0.35, color="#7d6608"))
a.text(x0 + 2.2, y0 - 0.55, "bonded seal (epoxy potting) under every pipe", ha="center", fontsize=7.5)
a.plot([x0, x0 + 4.4], [y0 + H, y0 + H], color="#2e86c1", lw=1, ls="--")
a.text(x0 + 4.45, y0 + H, "water line", fontsize=7, va="center")
xs = [x0, x0 + 1.3, x0 + 2.6, x0 + 3.9]
wet = [H * 0.9, H * 0.9 * 0.35 / 1.45, H * 0.9 * 0.35 / 2.75]
for i, x in enumerate(xs):
    top = y0 + H + 0.5
    if i == 0:
        a.add_patch(Rectangle((x, y0), 0.35, top - y0, color=AL))
        continue
    if i < 3:
        a.add_patch(Rectangle((x, y0), 0.12, top - y0, color=FE))
        a.add_patch(Rectangle((x + 0.12, y0), 0.2, top - y0, color="#e6b0aa"))
        a.add_patch(Rectangle((x + 0.12, y0), 0.2, wet[i], color=AL))
    else:
        a.add_patch(Rectangle((x, y0), 0.25, top - y0, color=FE))
a.text(x0 + 0.05, y0 + H + 0.65, "Al rod", fontsize=7)
a.text(x0 + 1.3, y0 + H + 0.65, "clad pipes: Fe inside, Al outside", fontsize=7)
a.text(x0 + 3.6, y0 + H + 0.95, "Fe tank", fontsize=7)
a.text(x0 + 2.2, y0 - 1.05, "pink: aluminium masked (epoxy) above the height\nthat gives every annulus the same area", ha="center", fontsize=7.5)
a.set_xlim(0, 10.8)
a.set_ylim(-1.4, 5.9)
a.set_aspect("equal")
a.axis("off")

fig.tight_layout()
fig.savefig(out, dpi=130)
print("wrote", out)
