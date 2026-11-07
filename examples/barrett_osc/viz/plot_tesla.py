"""Spark-gap Tesla coil figure: top-load voltage and the sparks × signal mixing.

go run ./examples/barrett_osc -data DIR
python3 examples/barrett_osc/viz/plot_tesla.py DIR OUT.png
"""
import csv
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

d, out = sys.argv[1], sys.argv[2]
fig, ax = plt.subplots(3, 2, figsize=(12, 9.5), constrained_layout=True)
for col, tag, name in [(0, "ordinary", "ordinary solenoid (526 kHz)"), (1, "bifilar", "Tesla bifilar solenoid (41.5 kHz)")]:
    w = list(csv.DictReader(open(f"{d}/tesla_{tag}_wave.csv")))
    t = [float(r["t"]) * 1e3 for r in w]
    n = int(len(t) * 0.25)
    ax[0][col].plot(t[:n], [float(r["both"]) / 1e3 for r in w[:n]], color="#c8553d", lw=0.5)
    ax[0][col].set_title(f"{name}: top-load voltage, sparks + signal")
    ax[0][col].set_ylabel("kV")
    m = int(len(t) * 0.12)
    ax[1][col].plot(t[:m], [float(r["signal"]) for r in w[:m]], color="#3b6ea5", lw=0.7, label="signal alone (no sparks)")
    ax[1][col].plot(t[:m], [float(r["mix"]) for r in w[:m]], color="#6a4c93", lw=0.7, label="sparks × signal (the multiplication)")
    ax[1][col].set_ylabel("V at the top-load")
    ax[1][col].set_xlabel("time [ms]")
    ax[1][col].legend(fontsize=8, loc="upper right")
    ax[1][col].set_title("the signal's part: alone, and what the switching adds")
    s = list(csv.DictReader(open(f"{d}/tesla_{tag}_spectrum.csv")))
    f = [float(r["f"]) / 1e3 for r in s]
    ax[2][col].semilogy(f, [max(float(r["signal"]), 1e-9) for r in s], color="#3b6ea5", lw=1.2, label="signal alone")
    ax[2][col].semilogy(f, [max(float(r["mix"]), 1e-9) for r in s], color="#6a4c93", lw=1.0, label="sparks × signal")
    ax[2][col].set_xlabel("frequency [kHz]")
    ax[2][col].set_ylabel("line amplitude [V]")
    ax[2][col].set_title("sidebands every 1 kHz (the break rate) around the signal")
    ax[2][col].legend(fontsize=8, loc="upper right")
fig.savefig(out, dpi=120)
