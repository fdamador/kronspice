"""Pictures for examples/zpf_hydrogen (run `go run . -viz viz` first).

    python3 viz/render.py

Writes viz/zpf_hydrogen.png: power absorbed from the zero-point field and
power radiated on Coulomb orbits, and the oscillator's energy in the
simulated field. Light background and large type.
"""
import csv
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
plt.rcParams.update({"font.size": 15, "axes.titlesize": 16, "figure.facecolor": "white",
                     "axes.facecolor": "white", "savefig.facecolor": "white"})

fig, axs = plt.subplots(1, 2, figsize=(18, 8))
rows = list(csv.DictReader(open(os.path.join(HERE, "balance.csv"))))
r = [float(x["r_over_a0"]) for x in rows]
ax = axs[0]
ax.plot(r, [float(x["absorbed"]) for x in rows], color="#1f5fbf", lw=6, label="absorbed from the zero-point field (eq. 18)")
ax.plot(r, [float(x["radiated"]) for x in rows], color="#c0392b", lw=4, ls="--", label="radiated, Larmor (eq. 19)")
ax.axvline(1, color="black", lw=1.5, ls=":")
ax.text(1.04, 3e-6, "balance at the Bohr radius a0", fontsize=14)
ax.set_yscale("log")
ax.set_xlabel("orbit radius ÷ Bohr radius a0")
ax.set_ylabel("power (W)")
ax.set_title("Puthoff 1987: the orbit gains inside a0, loses outside", loc="left")
ax.grid(alpha=0.3, which="both")
ax.legend(fontsize=13, loc="upper right")

ax = axs[1]
rows = list(csv.DictReader(open(os.path.join(HERE, "oscillator.csv"))))
t = [float(x["t"]) for x in rows]
e = [float(x["energy"]) for x in rows]
ax.plot(t, e, color="#16a085", lw=1.2, label="oscillator energy, one run")
mean = sum(e) / len(e)
ax.axhline(mean, color="#c0392b", lw=9, alpha=0.8, label="time average %.3f ħω0" % mean, zorder=4)
ax.axhline(0.5, color="black", lw=3, ls="--", label="½ ħω0", zorder=5)
ax.set_xlabel("time (scaled: 1/ω0)")
ax.set_ylabel("energy ÷ ħω0")
ax.set_title("The oscillator in the simulated zero-point field", loc="left")
ax.grid(alpha=0.3)
ax.legend(fontsize=13, loc="upper right")
fig.tight_layout()
fig.savefig(os.path.join(HERE, "zpf_hydrogen.png"), dpi=100)
print("wrote", os.path.join(HERE, "zpf_hydrogen.png"))
