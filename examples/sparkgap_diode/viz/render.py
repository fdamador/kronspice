"""Pictures for examples/sparkgap_diode (run `go run . -viz viz` first).

    python3 viz/render.py

Writes viz/sparkgap_diode.png: the I–V curves of electrode pairs in water,
the water-cell bridge's output, and the spark-gap rectifiers' output.
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


def rows(name):
    r = list(csv.reader(open(os.path.join(HERE, name))))
    return r[0], [[float(x) for x in row] for row in r[1:]]


fig, axs = plt.subplots(2, 2, figsize=(17, 12))
head, iv = rows("iv.csv")
ax = axs[0][0]
colours = ["#c0392b", "#e67e22", "#7f8c8d", "#f1c40f", "#1f5fbf", "#16a085", "#8e44ad"]
for k, name in enumerate(head[1:]):
    ax.plot([r[0] for r in iv], [r[k + 1] for r in iv], lw=[9, 6.5, 3, 4.5, 7, 4, 2.5][k], color=colours[k], label=name, alpha=0.85)
ax.set_yscale("symlog", linthresh=0.01)
ax.set_xlabel("volts (first electrode positive →)")
ax.set_ylabel("amps")
ax.set_title("Electrode pairs in water: aluminium makes it one-way", loc="left")
ax.grid(alpha=0.3)
ax.legend(fontsize=12, loc="upper left")

for ax, (fname, title, colr) in zip([axs[0][1], axs[1][0], axs[1][1]], [
        ("bridge_cells.csv", "Bridge of four aluminium–lead water cells (100 V peak)", "#1f5fbf"),
        ("gap_resistive.csv", "Polarity spark gap, half-wave, resistive load (4 kV peak)", "#c0392b"),
        ("gap_bridge_ratio3.csv", "Spark-gap bridge, polarity ratio 3 (6 kV peak)", "#8e44ad")]):
    _, d = rows(fname)
    keep = [r for r in d if r[0] <= 200]
    ax.plot([r[0] for r in keep], [r[1] for r in keep], color=colr, lw=1.8)
    ax.set_xlabel("time (ms)")
    ax.set_ylabel("output volts")
    ax.set_title(title, loc="left")
    ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "sparkgap_diode.png"), dpi=100)
print("wrote", os.path.join(HERE, "sparkgap_diode.png"))

# Frequency trend: amps per volt against frequency, every device on one chart.
groups = {}
for row in list(csv.DictReader(open(os.path.join(HERE, "frequency_trend.csv")))):
    groups.setdefault(row["device"], []).append((float(row["f_hz"]), float(row["a_per_v"])))
FLOOR = 2e-3
fig, ax = plt.subplots(figsize=(15, 9))
styles = [("#1f5fbf", "o", 9, "-"), ("#16a085", "s", 5, "--"), ("#e67e22", "D", 7, "-"),
          ("#c0392b", "^", 9, "-"), ("#8e44ad", "v", 9, "--")]
for (name, pts), (col, mk, lw, ls) in zip(groups.items(), styles):
    f = [p[0] for p in pts]
    y = [max(p[1], FLOOR) for p in pts]
    ax.plot(f, y, color=col, marker=mk, ms=12, lw=lw if lw < 8 else 4, ls=ls, label=name, alpha=0.9)
ideal = 2 / 3.141592653589793 / 1.1
ax.axhline(ideal, color="black", lw=2, ls=":", label="perfect diode bridge (%.2f A/V)" % ideal)
ax.axhspan(FLOOR * 0.7, FLOOR * 1.4, color="#f5b7b1", alpha=0.5)
ax.text(60, FLOOR * 1.6, "no output: never strikes, or breaks down backwards (shorts)", fontsize=14, color="#922b21")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_ylim(FLOOR * 0.7, 1.2)
ax.set_xlabel("frequency (Hz)")
ax.set_ylabel("amps per volt (mean load amps ÷ peak volts)")
ax.set_title("Water diodes and spark gaps in the same bridge: amps per volt vs frequency", loc="left")
ax.grid(alpha=0.3, which="both")
ax.legend(fontsize=13, loc="center left")
fig.tight_layout()
fig.savefig(os.path.join(HERE, "frequency_trend.png"), dpi=100)
print("wrote", os.path.join(HERE, "frequency_trend.png"))

# Pulse test: current in the first 50 µs of a 10 J pulse.
pul = {}
for row in csv.DictReader(open(os.path.join(HERE, "pulses.csv"))):
    pul.setdefault(row["device"], []).append((float(row["t_us"]), float(row["amps"])))
fig, ax = plt.subplots(figsize=(15, 8))
for (name, pts), (col, mk, lw, ls) in zip(pul.items(), styles):
    ax.plot([p[0] for p in pts], [p[1] for p in pts], color=col, lw=3, ls=ls, label=name)
ax.set_xscale("log")
ax.set_xlim(0.01, 50)
ax.set_xlabel("time (µs)")
ax.set_ylabel("amps")
ax.set_title("The same 10 J pulse through each device", loc="left")
ax.grid(alpha=0.3, which="both")
ax.legend(fontsize=13, loc="upper left")
fig.tight_layout()
fig.savefig(os.path.join(HERE, "pulses.png"), dpi=100)
print("wrote", os.path.join(HERE, "pulses.png"))
