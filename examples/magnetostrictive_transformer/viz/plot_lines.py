"""H37 figure: what reaches the secondary for each core, and the swap test.

go run ./examples/magnetostrictive_transformer -data DIR
python3 examples/magnetostrictive_transformer/viz/plot_lines.py DIR OUT.png
"""
import csv
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

d, out = sys.argv[1], sys.argv[2]
rows = list(csv.DictReader(open(f"{d}/lines.csv")))
cores = []
for r in rows:
    if r["core"] not in cores:
        cores.append(r["core"])
groups = [("fe", "transformer\n(fe)"), ("fa", "ultrasound read\nby the coil (fa)"), ("fb", "(fb)"),
          ("fa+fe", "fa+fe"), ("fa-fe", "fa−fe"), ("fb+fe", "fb+fe"), ("fb-fe", "fb−fe"), ("fa+fb", "fa+fb"), ("fb-fa", "fb−fa")]
colors = ["#3b6ea5", "#c8553d", "#6a4c93", "#2a9d8f", "#e9c46a"]
fig, ax = plt.subplots(2, 1, figsize=(12, 9), constrained_layout=True)
x = np.arange(len(groups))
w = 0.16
for i, c in enumerate(cores):
    fw = next(r for r in rows if r["core"] == c and r["direction"] == "forward")
    vals = [max(float(fw[g]) * 1e3, 1e-7) for g, _ in groups]
    ax[0].bar(x + (i - 2) * w, vals, w, label=c, color=colors[i])
ax[0].set_yscale("log")
ax[0].set_xticks(x)
ax[0].set_xticklabels([l for _, l in groups], fontsize=8)
ax[0].axvspan(2.5, 8.5, color="0.93", zorder=-1)
ax[0].text(5.5, 3e3, "products of the electrical drive and the ultrasound (the nonlinear part)", ha="center", fontsize=9)
ax[0].set_ylabel("secondary (load) current line [mA]")
ax[0].set_ylim(1e-7, 1e4)
ax[0].set_title("Forward: drive the primary at 5 kHz, transducers at the core's modes — what the secondary sees")
ax[0].legend(fontsize=8, ncol=1, loc="lower left")
prod = ["fa+fe", "fa-fe", "fb+fe", "fb-fe"]
fwd, rev = [], []
for c in cores:
    f = next(r for r in rows if r["core"] == c and r["direction"] == "forward")
    b = next(r for r in rows if r["core"] == c and r["direction"] == "reverse")
    fwd.append(sum(float(f[p]) for p in prod) * 1e3)
    rev.append(sum(float(b[p]) for p in prod) * 1e3)
xx = np.arange(len(cores))
ax[1].bar(xx - 0.18, fwd, 0.36, label="forward: drive primary, watch secondary", color="#3b6ea5")
ax[1].bar(xx + 0.18, rev, 0.36, label="swapped: drive secondary, watch primary", color="#c8553d")
ax[1].set_yscale("log")
ax[1].set_xticks(xx)
ax[1].set_xticklabels(cores, fontsize=9)
ax[1].set_ylabel("drive × ultrasound products, summed [mA]")
ax[1].set_title("Swap test: the linear transfer is reciprocal (to 1e-10); the nonlinear products need not be")
ax[1].legend(fontsize=9)
for i in range(len(cores)):
    ax[1].text(i, max(fwd[i], rev[i]) * 1.6, f"×{rev[i]/fwd[i]:.2f}", ha="center", fontsize=9)
ax[1].set_ylim(min(fwd + rev) / 3, max(fwd + rev) * 8)
fig.savefig(out, dpi=125)
