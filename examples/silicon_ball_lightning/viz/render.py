# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Figure for the silicon ball-lightning example.

    go run ./examples/silicon_ball_lightning -viz examples/silicon_ball_lightning/viz
    python examples/silicon_ball_lightning/viz/render.py
"""
import csv
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def read(name):
    with open(os.path.join(HERE, name), newline="") as f:
        rows = list(csv.DictReader(f))
    out = {}
    for k in rows[0]:
        try:
            out[k] = np.array([float(r[k]) for r in rows])
        except ValueError:
            out[k] = [r[k] for r in rows]
    return out


def main():
    el, bt, wd, ln = read("elements.csv"), read("burntime.csv"), read("weld.csv"), read("line.csv")
    fig, axs = plt.subplots(2, 2, figsize=(14, 10))

    ax = axs[0, 0]
    x = np.arange(len(el["element"]))
    ax.bar(x - 0.2, el["burn_kJ_per_g"], 0.4, label="released on burning in air (kJ/g)", color="C3")
    ax.bar(x + 0.2, el["vaporize_kJ_per_g"], 0.4, label="cost to make the vapor (kJ/g)", color="C0")
    for i, r in enumerate(el["return"]):
        ax.text(i, max(el["burn_kJ_per_g"][i], el["vaporize_kJ_per_g"][i]) + 0.8, f"×{r:.2f}", ha="center")
    ax.set_xticks(x, el["element"])
    ax.set_title("The chemical store per gram (return = burn ÷ vaporize)\nsilicon in soil; iron, manganese, silicon in welding rods")
    ax.legend(fontsize=8)

    ax = axs[0, 1]
    for r in (5, 10, 25, 50, 100):
        m = bt["r_nm"] == r
        ax.semilogy(bt["T_K"][m], bt["t_s"][m], label=f"particle radius {r:g} nm")
    ax.axhspan(1, 10, color="gold", alpha=0.35, label="observed lifetimes 1–10 s")
    ax.axvspan(973, 1473, color="0.85", alpha=0.6, label="Deal–Grove fitted range")
    ax.axvline(1687, color="k", ls=":", lw=1)
    ax.text(1700, 1e5, "Si melts", fontsize=8)
    ax.axvline(3140, color="C3", ls="--", lw=1)
    ax.text(3150, 1e5, "Stephan & Massey\n≈ 3140 K", fontsize=8, color="C3")
    ax.set_ylim(1e-2, 1e7)
    ax.set_xlabel("temperature (K)")
    ax.set_ylabel("time to burn through (s), air at sea level")
    ax.set_title("How long silicon nanoparticles burn (Deal–Grove kinetics)")
    ax.legend(fontsize=7, loc="lower left")

    ax = axs[1, 0]
    ax.plot(wd["t_ms"], wd["v_source"], color="0.6", label="welder open-circuit voltage (V)")
    ax.plot(wd["t_ms"], wd["v_arc"], "C3", label="arc voltage (V)")
    ax2 = ax.twinx()
    ax2.plot(wd["t_ms"], wd["i_arc"], "C0", label="arc current (A)")
    ax.set_xlabel("time (ms)")
    ax.set_ylabel("voltage (V)")
    ax2.set_ylabel("current (A)", color="C0")
    ax.set_title("kronspice: AC stick welder at Paiva et al.'s point\n(120 A rms, 23 V, 2.5 kW, 60 Hz)")
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=8, loc="lower left")

    ax = axs[1, 1]
    for a, c in ((0.15, "C0"), (2.5, "C3")):
        m = ln["a_m"] == a
        ax.loglog(ln["R_ohm"][m], ln["P_W"][m], c, label=f"ball radius {a:g} m, 20 m from a 35 kV line")
    ax.axhline(100, color="k", ls="--")
    ax.text(2, 140, "≈ 100 W the ball shines", fontsize=8)
    ax.set_ylim(1e-12, 1e3)
    ax.set_xlabel("ball resistance R (Ω)")
    ax.set_ylabel("mean power from the line (W)")
    ax.set_title("kronspice: what the power line can drive into the ball\n(it pulses at 100 Hz; observed light 99.4 Hz)")
    ax.legend(fontsize=8, loc="lower center")

    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "silicon_ball_lightning.png"), dpi=120)
    print("wrote silicon_ball_lightning.png")


if __name__ == "__main__":
    main()
