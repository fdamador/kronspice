"""Layer 7 figure: the cavity voltage with a linear and a nonlinear ionosphere.

go run ./examples/earth_machine -layer 7 -nldata DIR
python3 examples/earth_machine/viz/plot_ionosphere_nl.py DIR OUT.png
"""
import csv
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

d, out = sys.argv[1], sys.argv[2]
rows = list(csv.DictReader(open(f"{d}/ionosphere_nl.csv")))
t = [float(r["t"]) - float(rows[0]["t"]) for r in rows]
lin = [float(r["station_linear"]) for r in rows]
nl = [float(r["station_nonlinear"]) for r in rows]
sp = list(csv.DictReader(open(f"{d}/ionosphere_nl_spectrum.csv")))
f = [float(r["f"]) for r in sp]
slin = [max(float(r["station_linear"]), 1e-7) for r in sp]
snl = [max(float(r["station_nonlinear"]), 1e-7) for r in sp]

fig, ax = plt.subplots(3, 1, figsize=(9, 9.5), constrained_layout=True)
n = int(2.0 / (t[1] - t[0]))
ax[0].plot(t[:n], lin[:n], color="#3b6ea5", lw=1.6, label="linear ionosphere: f1 + f2 (sum)")
ax[0].plot(t[:n], nl[:n], color="#c8553d", lw=1.0, label="nonlinear ionosphere (heating)")
ax[0].set_ylabel("station voltage V = E·h [V]")
ax[0].set_title("Two transmitters, 7.8 Hz + 14.3 Hz, station 30° away (drive: source band at 0.5 Vh)")
ax[0].legend(loc="upper right", fontsize=9)
diff = [a - b for a, b in zip(nl, lin)]
ax[1].plot(t[:n], diff[:n], color="#6a4c93", lw=1.0)
ax[1].set_ylabel("nonlinear − linear [V]")
ax[1].set_title("What the nonlinearity adds: the products of the two waves (the 'multiplication')")
ax[2].semilogy(f, slin, color="#3b6ea5", lw=1.2, label="linear")
ax[2].semilogy(f, snl, color="#c8553d", lw=1.0, label="nonlinear")
for x, lab in [(1.3, "2f1−f2"), (20.8, "2f2−f1\n(= mode 3)"), (29.9, "2f1+f2"), (23.4, "3f1"), (36.4, "2f2+f1"), (42.9, "3f2"), (7.8, "f1"), (14.3, "f2")]:
    ax[2].axvline(x, color="0.75", lw=0.6, ls="--")
    ax[2].text(x, 2e-7 if lab.startswith("f") else 3e-6, lab, fontsize=8, ha="center")
ax[2].set_xlim(0, 45)
ax[2].set_ylim(1e-7, 1e3)
ax[2].set_xlabel("frequency [Hz]")
ax[2].set_ylabel("line amplitude [V]")
ax[2].set_title("Station spectrum: new lines appear only with the nonlinear ionosphere")
ax[2].legend(loc="upper right", fontsize=9)
fig.savefig(out, dpi=130)
