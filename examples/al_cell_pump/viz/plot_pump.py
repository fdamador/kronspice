"""H38 figure: the evolved aluminium-cell ground pump.

go run ./examples/al_cell_pump -data DIR
python3 examples/al_cell_pump/viz/plot_pump.py DIR OUT.png

(a) the Pareto front: catalytic efficiency against bank voltage, coloured by load current;
(b) where the input goes for the two baselines and the most efficient evolved design;
(c) the most efficient design's bank charging from rest;
(d) a multi-plate stack: cells in series for voltage, plates in parallel for current.
"""
import csv
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle

d, out = sys.argv[1], sys.argv[2]
par = list(csv.DictReader(open(f"{d}/pareto.csv")))
base = list(csv.DictReader(open(f"{d}/baseline.csv")))
bank = list(csv.DictReader(open(f"{d}/best_bank.csv")))
f = lambda r, k: float(r[k])

fig, ax = plt.subplots(2, 2, figsize=(13, 9.5))
fig.suptitle("H38 — aluminium-cell negative ground pump, evolved for catalytic efficiency (η = P_load / P_in)", fontsize=12)

# (a) Pareto front
a = ax[0][0]
v = [-f(r, "out_V") for r in par]
e = [f(r, "eff") for r in par]
i = [1e3 * f(r, "iout_A") for r in par]
sc = a.scatter(v, e, c=i, cmap="viridis", norm=LogNorm(), s=46, edgecolor="k", linewidth=0.4)
fig.colorbar(sc, ax=a, label="load current [mA]")
for r in base:
    a.scatter(-f(r, "out_V"), max(f(r, "eff"), 1e-5), marker="x", s=70, color="#c0392b")
    a.annotate(r["name"] + " (baseline)", (-f(r, "out_V"), 0.02), fontsize=8, color="#c0392b", ha="center")
a.set_xscale("log")
a.set_xlabel("bank voltage |V| [V]  (bank is negative)")
a.set_ylabel("catalytic efficiency η")
a.set_ylim(0, 1)
a.set_title("(a) Pareto front: 38 designs, 40 × 25 NSGA-II")
a.grid(alpha=0.3)

# (b) loss breakdown
a = ax[0][1]
best = par[0]
rows = [(base[0]["name"], base[0], "pout"), (base[1]["name"], base[1], "pout"), ("evolved, most efficient", best, None)]
cats = [("load (useful)", "#2e86c1"), ("coil", "#e67e22"), ("forward (electrolysis + electrolyte)", "#8e44ad"), ("leakage", "#16a085"), ("rest", "#95a5a6")]
for k, (name, r, pk) in enumerate(rows):
    pin = f(r, "pin_W")
    pout = f(r, pk) if pk else f(r, "eff") * pin
    vals = [pout, f(r, "coil"), f(r, "forward"), f(r, "leak"), f(r, "rest")]
    left = 0
    for (lab, col), x in zip(cats, vals):
        a.barh(k, 100 * x / pin, left=left, color=col, label=lab if k == 0 else None)
        left += 100 * x / pin
    a.text(101, k, f"{pin:.3g} W in", va="center", fontsize=8)
a.set_yticks(range(3), [r[0] for r in rows])
a.set_xlim(0, 118)
a.set_xlabel("share of the input power [%]", labelpad=1)
a.set_title("(b) where the input goes")
a.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=3)

# (c) bank charging
a = ax[1][0]
t = [1e3 * f(r, "t") for r in bank]
b = [f(r, "bank") for r in bank]
a.plot(t, b, color="#2e86c1")
a.axhline(f(best, "out_V"), ls="--", color="k", lw=0.8)
a.text(t[-1], f(best, "out_V"), f"  {f(best, 'out_V'):.0f} V", va="bottom", ha="right", fontsize=8)
a.set_xlabel("time [ms]")
a.set_ylabel("bank voltage [V]")
a.set_title(f"(c) most efficient design charging: η {f(best, 'eff'):.3f}, {1e3*f(best, 'iout_A'):.3g} mA, {f(best, 'pin_W'):.3g} W in")
a.grid(alpha=0.3)

# (d) stack drawing: Ns cells in series, Np plates in parallel per electrode
a = ax[1][1]
a.set_title("(d) a stack: cells in series (voltage) × plates in parallel (current)")
ns, npl = 3, 3
W, H = 2.6, 3.0
for c in range(ns):
    x0 = 0.4 + c * (W + 0.6)
    a.add_patch(Rectangle((x0, 0.3), W, H, fill=True, color="#d6eaf8", ec="k"))
    a.text(x0 + W / 2, 0.05, f"cell {c+1}: Al and plate in water", ha="center", fontsize=7)
    for p in range(npl):
        y = 0.6 + p * 0.9
        a.add_patch(Rectangle((x0 + 0.25, y), 0.12, 0.6, color="#7f8c8d"))        # Al anode plate (oxide)
        a.add_patch(Rectangle((x0 + 0.22, y), 0.03, 0.6, color="#c0392b"))        # oxide film
        a.add_patch(Rectangle((x0 + W - 0.37, y), 0.12, 0.6, color="#34495e"))   # cathode plate
    a.plot([x0 + 0.31] * 2, [0.6, 0.6 + (npl - 1) * 0.9 + 0.6], color="k", lw=1)
    a.plot([x0 + W - 0.31] * 2, [0.6, 0.6 + (npl - 1) * 0.9 + 0.6], color="k", lw=1)
    if c < ns - 1:
        a.plot([x0 + W - 0.31, x0 + W + 0.6 + 0.31], [H + 0.5] * 2, color="k", lw=1)
        a.plot([x0 + W - 0.31] * 2, [0.6 + (npl - 1) * 0.9 + 0.6, H + 0.5], color="k", lw=1)
        a.plot([x0 + W + 0.6 + 0.31] * 2, [0.6 + (npl - 1) * 0.9 + 0.6, H + 0.5], color="k", lw=1)
a.text(0.2, H + 0.75, "anode (Al, formed oxide, red)", fontsize=8)
a.text(0.4 + ns * (W + 0.6) - 0.6, H + 0.75, "cathode", fontsize=8, ha="right")
a.text(0.2, -0.45,
       "genes: cells in series, plates in parallel, plate area, gap, forming voltage\n"
       "holds Ns·Vform; forward Ns·2 V + Ns·ρ·gap/(A·Np); film ε·A·Np/(k·Vform·Ns)",
       fontsize=8)
a.set_xlim(0, 0.4 + ns * (W + 0.6))
a.set_ylim(-0.8, H + 1.1)
a.axis("off")

fig.tight_layout()
fig.savefig(out, dpi=130)
print("wrote", out)
