"""Pictures for examples/einstein_teleparallel (run `go run . -viz viz` first).

    python3 viz/render.py

Writes viz/einstein_teleparallel.png: the torsion components computed
three ways, and the motor built from Einstein's torsion against Kron's
frames. Light background and large type.
"""
import csv
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
plt.rcParams.update({"font.size": 15, "axes.titlesize": 16, "figure.facecolor": "white",
                     "axes.facecolor": "white", "savefig.facecolor": "white"})
names = ["1", "2s", "3s", "θ"]  # stator coil, stationary rotor coils 2̄ and 3̄, rotor angle

fig, axs = plt.subplots(1, 2, figsize=(18, 8))
rows = list(csv.DictReader(open(os.path.join(HERE, "torsion.csv"))))
ax = axs[0]
x = list(range(len(rows)))
labels = ["%s%s%s" % (names[int(r["lambda"])], names[int(r["mu"])], names[int(r["nu"])]) for r in rows]
ax.bar([i - 0.27 for i in x], [-float(r["einstein_frame"]) for r in rows], 0.27, color="#1f5fbf", label="Einstein 1928: −Λ (frame components)")
ax.bar(x, [float(r["kron_anholonomy"]) for r in rows], 0.27, color="#f1c40f", label="Kron's frame: anholonomy")
ax.bar([i + 0.27 for i in x], [float(r["hoffmann_omega"]) for r in rows], 0.27, color="#c0392b", label="Hoffmann's Ω (= L̄·Λ, henry)")
ax.set_xticks(x)
ax.set_xticklabels(labels, rotation=60, fontsize=12)
ax.axhline(0, color="black", lw=1)
ax.set_ylabel("torsion component")
ax.set_title("One torsion, three notations (Kron's motor, θ = 0.7)\n(axes: 1 = stator, 2s, 3s = stationary rotor coils, θ = rotor)", loc="left")
ax.grid(alpha=0.3, axis="y")
ax.legend(fontsize=12)

run = {}
for r in csv.DictReader(open(os.path.join(HERE, "speed.csv"))):
    run.setdefault(r["run"], []).append((float(r["t"]), float(r["speed"])))
ax = axs[1]
for (name, pts), col, ls, lw in zip(run.items(), ["#7f8c8d", "#1f5fbf", "#c0392b", "#8e44ad"], ["-", "-", ":", "-"], [9, 5, 3, 4]):
    ax.plot([p[0] for p in pts], [p[1] for p in pts], color=col, ls=ls, lw=lw, label=name)
ax.set_xlabel("time (s)")
ax.set_ylabel("rotor speed (rad/s)")
ax.set_title("kronspice: the motor built from Einstein's torsion", loc="left")
ax.grid(alpha=0.3)
ax.legend(fontsize=12, loc="center right")
fig.tight_layout()
fig.savefig(os.path.join(HERE, "einstein_teleparallel.png"), dpi=100)
print("wrote", os.path.join(HERE, "einstein_teleparallel.png"))
