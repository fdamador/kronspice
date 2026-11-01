"""Pictures for examples/water_lamp (run `go run . -viz viz` first).

    python3 viz/render.py

Writes
  viz/lamps.png   capacitor voltage of each flashing-lamp circuit (Part 1)
  viz/coil.png    the bifilar coil's voltage when the spark gaps drive it (Part 2)
Light background and large type.
"""
import csv
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
plt.rcParams.update({"font.size": 15, "axes.titlesize": 16, "figure.facecolor": "white",
                     "axes.facecolor": "white", "savefig.facecolor": "white"})


def load(name):
    rows = list(csv.reader(open(os.path.join(HERE, name))))[1:]
    return [float(r[0]) for r in rows], [float(r[1]) for r in rows]


panels = [
    ("neon_primed", "Neon NE-2 with radioactive starter", "#8e44ad", 0.3, 1, "s"),
    ("neon_dark", "Neon NE-2, no starter, in the dark", "#8e44ad", 7200, 1 / 60, "min"),
    ("vapour_uv", "Water-vapour tube, UV LED starter", "#1f5fbf", 1.0, 1, "s"),
    ("vapour_dark", "Water-vapour tube, no starter, in the dark", "#1f5fbf", 7200, 1 / 60, "min"),
    ("liquid", "Spark gap in liquid water", "#16a085", 0.03, 1e3, "ms"),
    ("mist", "Air gap with water mist", "#c0392b", 0.1, 1e3, "ms"),
]
fig, axs = plt.subplots(3, 2, figsize=(16, 13))
for ax, (tag, title, colr, span, scale, unit) in zip(axs.flat, panels):
    ax.set_title(title, loc="left")
    ax.grid(alpha=0.3)
    if "dark" in tag:
        # hours-long waits: mark each flash
        fl = [float(r[0]) for r in list(csv.reader(open(os.path.join(HERE, tag + "_flashes.csv"))))[1:]]
        ax.vlines([x / 3600 for x in fl], 0, 1, color=colr, lw=2.5)
        ax.set_xlim(0, 10)
        ax.set_ylim(0, 1.2)
        ax.set_yticks([])
        ax.set_xlabel("time (hours) — each line is one flash")
        continue
    t, v = load(tag + "_trace.csv")
    keep = [i for i, x in enumerate(t) if x <= span]
    ax.plot([t[i] * scale for i in keep], [v[i] for i in keep], color=colr, lw=2)
    ax.set_xlabel(f"time ({unit})")
    ax.set_ylabel("capacitor volts")
fig.suptitle("Flashing-lamp circuit: each drop is one flash", fontsize=19)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "lamps.png"), dpi=100)

fig, axs = plt.subplots(3, 1, figsize=(15, 13))
for ax, (tag, title, colr) in zip(axs, [("liquid", "Spark gap in liquid water drives the bifilar coil (untuned)", "#16a085"),
                                        ("mist", "Air gap with water mist drives the bifilar coil (untuned)", "#c0392b"),
                                        ("evolved", "Mist gap, flash circuit tuned to the coil and quenched (evolved)", "#8e44ad")]):
    t, v = load(tag + "_coil.csv")
    if tag == "evolved":  # zoom on the first flash
        i0 = next(i for i, x in enumerate(v) if abs(x) > 1)
        t0 = t[i0]
        keep = [i for i, x in enumerate(t) if t0 - 0.2 <= x <= t0 + 3]
    else:
        keep = [i for i, x in enumerate(t) if x <= 60]
    ax.plot([t[i] for i in keep], [v[i] for i in keep], color=colr, lw=1.2)
    ax.set_title(title, loc="left")
    ax.set_ylabel("coil volts")
    ax.set_xlabel("time (ms)")
    ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "coil.png"), dpi=100)
print("wrote lamps.png and coil.png")
