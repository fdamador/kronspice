# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Resonance and harmonics of the stacked dual motor (resonance.go).

Data: go run ./examples/homopolar -data examples/homopolar/viz/data
Usage: python viz/plot_resonance.py [data_dir] [out_dir]   (matplotlib, numpy)
"""
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

here = os.path.dirname(os.path.abspath(__file__))
data = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "data")
out = sys.argv[2] if len(sys.argv) > 2 else here
TEAL, ORANGE, INK, QUIET = "#127b8a", "#d9772b", "#1b1b1f", "#5b5b66"


def load(name):
    with open(os.path.join(data, name), newline="") as f:
        r = csv.reader(f)
        head = next(r)
        rows = np.array([[float(v) for v in row] for row in r])
    return {h: rows[:, i] for i, h in enumerate(head)}


c, w, g = load("res_curve.csv"), load("res_wave.csv"), load("res_tank.csv")
fig, ax = plt.subplots(1, 3, figsize=(12, 3.8), dpi=300)
fig.patch.set_facecolor("#fbfaf7")
for a in ax:
    a.set_facecolor("#fbfaf7")
    a.grid(alpha=0.3)
    a.spines["top"].set_visible(False)
    a.spines["right"].set_visible(False)

ax[0].plot(c["f"], c["dw"], color=TEAL, lw=1.6)
ax[0].set_xlabel("torque oscillation [Hz] (±0.5 N·m on 20 N·m)")
ax[0].set_ylabel("speed oscillation δω [rad/s]", color=TEAL)
a2 = ax[0].twinx()
a2.plot(c["f"], c["di"] / 1e3, color=ORANGE, lw=1.1, ls="--")
a2.set_ylabel("current oscillation δI [kA]", color=ORANGE)
a2.spines["top"].set_visible(False)
ax[0].set_title("(a) no capacitor: inertia and coils resonate,\npeak 0.575 Hz (linear f_em 0.587 Hz), Q = 41", fontsize=9.5, color=INK, loc="left")

t = w["t"] - w["t"][0]
ax[1].plot(t, w["i1"] / 1e3, color=ORANGE, lw=1.2)
ax[1].set_xlabel("time [s]")
ax[1].set_ylabel("loop current I₁ [kA]", color=ORANGE)
a3 = ax[1].twinx()
a3.plot(t, w["w1"] * 60 / (2 * np.pi), color=TEAL, lw=1.0)
a3.set_ylabel("shaft speed [rpm]", color=TEAL)
a3.spines["top"].set_visible(False)
ax[1].set_title("(b) driven ±5 N·m at f_em: a non-sinusoidal\ncurrent, 2nd harmonic 32 %, 3rd 7.5 %", fontsize=9.5, color=INK, loc="left")

ax[2].plot(g["t"], g["i1_300F"] / 1e3, color=QUIET, lw=1.6, label="300 F: DC holds, 2.3 kW in the load")
ax[2].plot(g["t"], g["i1_1000F"] / 1e3, color=ORANGE, lw=1.0, label="1000 F: AC, 2.0 kW")
ax[2].plot(g["t"], g["i1_12000F"] / 1e3, color=TEAL, lw=1.0, label="12 000 F: AC, 1.0 kW")
ax[2].axhline(0, color=INK, lw=0.6)
ax[2].set_xlim(30, 60)
ax[2].set_xlabel("time [s]  (torque-driven, 20 N·m)")
ax[2].set_ylabel("loop current I₁ [kA]")
ax[2].legend(fontsize=7, frameon=False, loc="lower left")
ax[2].set_title("(c) C across the load: the tank locks in\nas AC past C ≈ 450 F (free shafts)", fontsize=9.5, color=INK, loc="left")
fig.tight_layout()
fig.savefig(os.path.join(out, "dual_motor_resonance.png"), facecolor=fig.get_facecolor())
print("wrote dual_motor_resonance.png")
