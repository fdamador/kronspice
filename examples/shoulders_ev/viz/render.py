"""Pictures for examples/shoulders_ev (run `go run . -viz viz` first).

    python3 viz/render.py

Writes viz/shoulders_ev.png: the pressure an EV must hold against the
candidate containers, the electron counts Shoulders reports against what
his sources hold, and the kronspice scope traces. Light background and
large type.
"""
import csv
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
plt.rcParams.update({"font.size": 15, "axes.titlesize": 16, "figure.facecolor": "white",
                     "axes.facecolor": "white", "savefig.facecolor": "white"})


def rows(name):
    return list(csv.DictReader(open(os.path.join(HERE, name))))


fig, axs = plt.subplots(1, 3, figsize=(24, 8))

ax = axs[0]
p = rows("pressure.csv")
names = [r["name"] for r in p]
vals = [float(r["pa"]) for r in p]
cols = ["#c0392b"] + ["#1f5fbf"] * (len(p) - 1)
ax.barh(names[::-1], vals[::-1], color=cols[::-1])
ax.set_xscale("log")
ax.set_xlabel("pressure (Pa)")
ax.set_title("A 1 µm EV of 10¹¹ electrons pushes outward\nfar harder than anything known can hold", loc="left")
for i, v in enumerate(vals[::-1]):
    ax.text(v * 1.5, i, "%.0e" % v, va="center", fontsize=14)
ax.set_xlim(1e-4, 1e21)

ax = axs[1]
b = rows("budget.csv")
names = [r["name"].replace(" of a 3 µm EV on an electrode", "\n(3 µm EV)").replace(" by a 10 µm EV over 1 mm of track", "\n(10 µm EV, 1 mm)") for r in b]
vals = [float(r["electrons"]) for r in b]
cols = ["#16a085", "#16a085", "#e67e22", "#e67e22", "#e67e22"]
ax.barh(names[::-1], vals[::-1], color=cols[::-1])
ax.set_xscale("log")
ax.set_xlabel("electrons")
ax.set_title("Charge: what his sources hold (green)\nand what he reports (orange)", loc="left")
for i, v in enumerate(vals[::-1]):
    ax.text(v * 1.5, i, "%.1e" % v, va="center", fontsize=14)
ax.set_xlim(1e9, 1e16)

ax = axs[2]
tr = {}
for r in rows("scope.csv"):
    tr.setdefault(r["case"], []).append((float(r["t_ns"]), float(r["v"])))
for (name, pts), col in zip(tr.items(), ["#16a085", "#e67e22", "#c0392b"]):
    ax.plot([q[0] for q in pts], [q[1] / 1e3 for q in pts], color=col, lw=4, label=name)
ax.set_xlim(0, 1)
ax.set_xlabel("time (ns)")
ax.set_ylabel("scope voltage (kV)")
ax.set_title("kronspice: 1 pF anode into 50 Ω\n(source discharge vs. EV strikes)", loc="left")
ax.legend()
ax.grid(alpha=0.3)

fig.tight_layout()
fig.savefig(os.path.join(HERE, "shoulders_ev.png"), dpi=110)
