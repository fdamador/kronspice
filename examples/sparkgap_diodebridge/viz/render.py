"""Pictures for examples/sparkgap_diodebridge (run `go run . -viz viz` first).

    python3 viz/render.py

Writes viz/diodebridge.png: the surge test, the spark pulser, the oxide-film
varactor and the frequency doubler. Light background and large type.
"""
import csv
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
plt.rcParams.update({"font.size": 15, "axes.titlesize": 16, "figure.facecolor": "white",
                     "axes.facecolor": "white", "savefig.facecolor": "white"})


def groups(name, key):
    out = {}
    for row in csv.DictReader(open(os.path.join(HERE, name))):
        out.setdefault(row[key], []).append(row)
    return out


fig, axs = plt.subplots(2, 2, figsize=(18, 12))

ax = axs[0][0]
for (name, rows), col, lw in zip(groups("surge.csv", "case").items(), ["#c0392b", "#1f5fbf"], [3, 5]):
    ax.plot([float(r["t_us"]) for r in rows], [float(r["v_arm"]) for r in rows], color=col, lw=lw, label=name)
ax.axhline(140, color="black", ls=":", lw=2)
ax.lines[-1].set_label("a cell's oxide sparks over (140 V)")
ax.set_xlabel("time (µs)")
ax.set_ylabel("largest volts across any bridge cell")
ax.set_title("B. 4 kV surge: arrester gap + aluminium clamp protect the bridge", loc="left")
ax.grid(alpha=0.3)
ax.legend(fontsize=12, loc="upper right")

ax = axs[0][1]
for (name, rows), col, lw in zip(groups("pulser.csv", "design").items(), ["#e67e22", "#8e44ad", "#16a085"], [2.5, 2.5, 5]):
    ax.plot([float(r["t_ms"]) for r in rows], [float(r["v_store"]) / 1e3 for r in rows], color=col, lw=lw, label=name)
ax.set_xlabel("time (ms)")
ax.set_ylabel("capacitor kilovolts")
ax.set_title("C. Spark pulser: only the zero-crossing-fired gap works", loc="left")
ax.grid(alpha=0.3)
ax.legend(fontsize=11, loc="lower left")

ax = axs[1][0]
rows = list(csv.DictReader(open(os.path.join(HERE, "varactor.csv"))))
vf = [float(r["vform"]) for r in rows]
ax.plot(vf, [float(r["f0_formula"]) for r in rows], color="black", lw=2, ls=":", label="1/(2π√LC)")
ax.plot(vf, [float(r["f0"]) for r in rows], "o", color="#1f5fbf", ms=14, label="kronspice AC resonance")
for r in rows:
    ax.annotate("%.1f µF" % float(r["film_uF"]), (float(r["vform"]), float(r["f0"])), textcoords="offset points",
                xytext=(10, -22), fontsize=12)
ax.set_xlabel("volts the oxide film was formed at")
ax.set_ylabel("tank resonance (Hz)")
ax.set_title("E. Oxide film as a varactor with memory (1 mH tank)", loc="left")
ax.grid(alpha=0.3)
ax.legend(fontsize=13, loc="lower right")

ax = axs[1][1]
rows = list(csv.DictReader(open(os.path.join(HERE, "doubling.csv"))))
labels = ["DC", "2f (120 Hz)", "4f (240 Hz)"]
w = 0.26
for k, (r, col) in enumerate(zip(rows, ["#7f8c8d", "#1f5fbf", "#c0392b"])):
    ax.bar([x + (k - 1) * w for x in range(3)], [float(r["dc"]), float(r["h2"]), float(r["h4"])], w, color=col, label=r["bridge"])
ax.set_xticks(range(3))
ax.set_xticklabels(labels)
ax.set_ylabel("amplitude ÷ source peak")
ax.set_title("F. Frequency doubler: 60 Hz in, 120 Hz out", loc="left")
ax.grid(alpha=0.3, axis="y")
ax.legend(fontsize=12)

fig.tight_layout()
fig.savefig(os.path.join(HERE, "diodebridge.png"), dpi=100)
print("wrote", os.path.join(HERE, "diodebridge.png"))
