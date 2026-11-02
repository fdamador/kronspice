"""Pictures for examples/hoffmann_motor (run `go run . -viz viz` first).

    python3 viz/render.py

Writes viz/hoffmann_motor.png: start-up speed in each frame and the
torque–speed curve. Light background and large type.
"""
import csv
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
plt.rcParams.update({"font.size": 15, "axes.titlesize": 16, "figure.facecolor": "white",
                     "axes.facecolor": "white", "savefig.facecolor": "white"})

runs = {}
for row in csv.DictReader(open(os.path.join(HERE, "speed.csv"))):
    runs.setdefault(row["run"], []).append((float(row["t"]), float(row["speed"])))

fig, axs = plt.subplots(1, 2, figsize=(18, 8))
ax = axs[0]
style = {
    "from rest, frame I": ("#7f8c8d", "-", 4),
    "pushed to +40 rad/s, frame I": ("#1f5fbf", "-", 9),
    "pushed to +40 rad/s, frame II(b)": ("#f1c40f", "--", 5),
    "pushed to +40 rad/s, eq. 5′ (Runge–Kutta)": ("#c0392b", ":", 3.5),
    "pushed to −40 rad/s, frame I": ("#16a085", "-", 4),
    "pushed to +40 rad/s, no torsion (Hoffmann's warning)": ("#8e44ad", "-", 4),
}
for name, pts in runs.items():
    col, ls, lw = style.get(name, ("black", "-", 2))
    ax.plot([p[0] for p in pts], [p[1] for p in pts], color=col, ls=ls, lw=lw, label=name)
ax.axhline(377, color="black", lw=1, ls=":")
ax.axhline(-377, color="black", lw=1, ls=":")
ax.set_xlabel("time (s)")
ax.set_ylabel("rotor speed (rad/s)")
ax.set_title("Start-up: three solutions agree; without torsion it only coasts", loc="left")
ax.grid(alpha=0.3)
ax.legend(fontsize=11, loc="lower right", bbox_to_anchor=(1, 0.08))

ax = axs[1]
rows = list(csv.DictReader(open(os.path.join(HERE, "torque.csv"))))
w = [float(r["speed"]) for r in rows]
ax.plot(w, [float(r["torque"]) for r in rows], "o-", color="#1f5fbf", lw=6, ms=10, label="frame I")
ax.plot(w, [float(r["torque_frame2"]) for r in rows], "s--", color="#f1c40f", lw=3, ms=6, label="frame II(b)")
ax.plot(w, [float(r["torque_salient"]) for r in rows], "^-", color="#c0392b", lw=3, ms=9, label="salient pole, L_D = 0.02 H")
ax.axhline(0, color="black", lw=1)
ax.axvline(0, color="black", lw=1)
ax.set_xlabel("rotor speed (rad/s)")
ax.set_ylabel("mean torque (N·m)")
ax.set_title("Torque–speed: zero at standstill, either direction runs", loc="left")
ax.grid(alpha=0.3)
ax.legend(fontsize=13, loc="upper left")
fig.tight_layout()
fig.savefig(os.path.join(HERE, "hoffmann_motor.png"), dpi=100)
print("wrote", os.path.join(HERE, "hoffmann_motor.png"))
