"""Fig. 2G bench figure: load voltage with air-core and Metglas-core bifilar pancakes.

go run ./examples/barrett_osc -data DIR
python3 examples/barrett_osc/viz/plot_bench.py DIR OUT.png
"""
import csv
import math
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

d, out = sys.argv[1], sys.argv[2]
w = list(csv.DictReader(open(f"{d}/bench_wave.csv")))
t = [float(r["t"]) for r in w]
names = ["ordinary pancakes, air", "bifilar pancakes, air", "bifilar, Metglas", "bifilar, Metglas, DC bias", "bifilar, Metglas, bias, ribbon at idler"]
FP, FS = 10e3, 6.3e3


def lines_removed(y):
    """Remove the DC, the pump and its harmonics, and the signal (least squares):
    what is left are the products that hold both waves, pump × signal."""
    import numpy as np
    tt = np.array(t)
    fr = [k * FP for k in range(1, 9)] + [FS, 2 * FS, 3 * FS]
    A = np.column_stack([np.ones_like(tt)] + [g(2 * math.pi * f * tt) for f in fr for g in (np.cos, np.sin)])
    c, *_ = np.linalg.lstsq(A, np.array(y), rcond=None)
    return np.array(y) - A @ c, A @ c


air = [float(r["case1"]) for r in w]
mg = [float(r["case3"]) for r in w]
res_air, fit_air = lines_removed(air)
res_mg, fit_mg = lines_removed(mg)
n = int(2e-3 / (t[1] - t[0]))
ms = [x * 1e3 for x in t[:n]]
fig, ax = plt.subplots(3, 1, figsize=(9, 9.5), constrained_layout=True)
ax[0].plot(ms, [a - sum(air) / len(air) for a in air[:n]], color="#3b6ea5", lw=1.4, label="bifilar pancakes, air core")
ax[0].plot(ms, [a - sum(mg) / len(mg) for a in mg[:n]], color="#c8553d", lw=1.0, label="bifilar pancakes, Metglas core, DC bias")
ax[0].set_ylabel("load voltage (AC) [V]")
ax[0].set_xlabel("time [ms]")
ax[0].set_title("Fig. 2G bench: pump 10 kHz + signal 6.3 kHz at the shuttle (load voltage)")
ax[0].legend(fontsize=9, loc="upper right")
ax[1].plot(ms, res_air[:n] * 1e3, color="#3b6ea5", lw=1.4, label="air core: residual (round-off)")
ax[1].plot(ms, res_mg[:n] * 1e3, color="#c8553d", lw=1.0, label="Metglas core: the mixing products")
ax[1].set_ylabel("pump × signal products [mV]")
ax[1].set_xlabel("time [ms]")
ax[1].set_title("After subtracting the pump (and its harmonics) and the signal: the 'multiplication'")
ax[1].legend(fontsize=9, loc="upper right")
s = list(csv.DictReader(open(f"{d}/bench_spectrum.csv")))
f = [float(r["f"]) / 1e3 for r in s]
cols = {"case1": ("#3b6ea5", 1.4), "case2": ("#999999", 1.0), "case3": ("#c8553d", 1.0)}
for k, (col, lw) in cols.items():
    ax[2].semilogy(f, [max(float(r[k]), 1e-12) for r in s], color=col, lw=lw, label=names[int(k[-1])])
for x, lab in [(3.7, "fp−fs"), (6.3, "fs"), (10, "fp"), (13.7, "2fp−fs"), (16.3, "fp+fs"), (20, "2fp"), (7.4, "2(fp−fs)")]:
    ax[2].axvline(x, color="0.8", lw=0.6, ls="--")
    ax[2].text(x, 3e-12, lab, fontsize=7, ha="center", rotation=90, va="bottom")
ax[2].set_xlim(0, 40)
ax[2].set_ylim(1e-12, 10)
ax[2].set_xlabel("frequency [kHz]")
ax[2].set_ylabel("load line amplitude [V]")
ax[2].set_title("Load spectrum: only the Metglas core makes new lines")
ax[2].legend(fontsize=8, loc="upper right")
fig.savefig(out, dpi=130)
