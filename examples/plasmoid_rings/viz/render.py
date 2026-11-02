"""Pictures for examples/plasmoid_rings (run `go run . -viz viz` first).

    python3 viz/render.py

Writes viz/plasmoid_rings.png: field lines of the force-free plasmoid on
nested tori, the spheromak ring's energy decay at three temperatures, and
the collapsing electron ring. Light background and large type.
"""
import csv
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
plt.rcParams.update({"font.size": 15, "axes.titlesize": 16, "figure.facecolor": "white",
                     "axes.facecolor": "white", "savefig.facecolor": "white"})

fig = plt.figure(figsize=(24, 8))
ax = fig.add_subplot(1, 3, 1, projection="3d")
lines = {}
for r in csv.DictReader(open(os.path.join(HERE, "fieldlines.csv"))):
    lines.setdefault(r["line"], []).append((float(r["x"]) * 100, float(r["y"]) * 100, float(r["z"]) * 100))
for (name, pts), col in zip(lines.items(), ["#c0392b", "#1f5fbf", "#16a085"]):
    ax.plot([p[0] for p in pts], [p[1] for p in pts], [p[2] for p in pts], color=col, lw=1.2)
ax.set_xlabel("x (cm)")
ax.set_ylabel("y (cm)")
ax.set_zlabel("z (cm)")
ax.set_title("Force-free plasmoid: field lines on nested tori\naround the current ring", loc="left")

ax = fig.add_subplot(1, 3, 2)
dec = {}
for r in csv.DictReader(open(os.path.join(HERE, "decay.csv"))):
    dec.setdefault(r["T_eV"], []).append((float(r["t"]), float(r["energy"])))
for (T, pts), col in zip(dec.items(), ["#c0392b", "#e67e22", "#1f5fbf"]):
    ax.plot([p[0] * 1e6 for p in pts], [p[1] / 1e6 for p in pts], color=col, lw=4, label="%s eV" % T)
ax.set_xscale("log")
ax.set_xlabel("time (µs)")
ax.set_ylabel("energy left (MJ)")
ax.set_title("Spheromak ring (1 MJ, 10 cm) in kronspice:\nresistance drains it in µs–ms", loc="left")
ax.grid(alpha=0.3, which="both")
ax.legend(fontsize=13, title="electron temperature")

ax = fig.add_subplot(1, 3, 3)
rows = list(csv.DictReader(open(os.path.join(HERE, "ring.csv"))))
t = [float(r["t"]) * 1e6 for r in rows]
ax.plot(t, [float(r["radius"]) * 100 for r in rows], color="#1f5fbf", lw=5, label="ring radius (cm)")
ax.plot(t, [float(r["B"]) * 10 for r in rows], color="#7f8c8d", lw=3, ls="--", label="field × 10 (T)")
ax.plot(t, [float(r["gamma"]) for r in rows], color="#c0392b", lw=3, label="electron energy (γ)")
ax.set_xlabel("time (µs)")
ax.set_title("Collapsing electron ring:\nthe bank's rising field squeezes it", loc="left")
ax.grid(alpha=0.3)
ax.legend(fontsize=13)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "plasmoid_rings.png"), dpi=90)
print("wrote", os.path.join(HERE, "plasmoid_rings.png"))
