# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Charts for examples/magnetic_gas from the CSVs that
`go run ./examples/magnetic_gas -viz examples/magnetic_gas/viz/data` writes.

    python3 examples/magnetic_gas/viz/plot.py

Writes paschen.png, marx.png and pump.png beside this script.
"""
import csv
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
# Categorical slots in fixed order, one per gas (validated: CVD ΔE ≥ 9;
# aqua and yellow are below 3:1 on the surface, so every line is labelled).
GAS = [("argon", "Ar", "#2a78d6"), ("hydrogen", "H₂", "#eb6834"),
       ("oxygen", "O₂", "#1baf7a"), ("hho", "HHO (est.)", "#eda100")]
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
PSI500_TORR = 500 * 6894.757293168 / (101325 / 760)


def read(name):
    with open(os.path.join(DATA, name), newline="") as f:
        rows = list(csv.reader(f))
    head = rows[0]
    cols = {h: [] for h in head}
    for r in rows[1:]:
        for h, v in zip(head, r):
            cols[h].append(float(v) if v else float("nan"))
    return cols


def style(ax, title, sub):
    ax.set_facecolor(SURF)
    ax.figure.set_facecolor(SURF)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=9)
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.set_title(title, loc="left", fontsize=13, color=INK, pad=22)
    ax.text(0, 1.02, sub, transform=ax.transAxes, fontsize=9, color=INK2)


def paschen():
    d = read("paschen.csv")
    fig, ax = plt.subplots(figsize=(8, 5))
    pd = d["pd_torr_cm"]
    for key, label, col in GAS:
        ax.plot(pd, d[key], color=col, linewidth=2)
        ax.text(pd[-1] * 1.05, d[key][-1], label, color=INK, fontsize=9, va="center")
    x = PSI500_TORR * 0.02  # 0.2 mm gap
    ax.axvline(x, color=INK2, linewidth=1, linestyle=(0, (4, 3)))
    ax.text(x * 1.08, 150, "500 psi, 0.2 mm gap", color=INK2, fontsize=9)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.2, pd[-1] * 1.6)
    ax.set_xlabel("pressure × gap (Torr·cm)", color=INK2)
    ax.set_ylabel("breakdown voltage (V)", color=INK2)
    style(ax, "Argon sparks at the lowest voltage, oxygen at the highest",
          "Uniform-field breakdown from each gas's measured Paschen minimum (Naidu); HHO estimated. Upper bound at high pd.")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "paschen.png"), dpi=150)


def marx():
    fig, ax = plt.subplots(figsize=(8, 5))
    for key, label, col in GAS:
        d = read(f"marx_{key}.csv")
        t = [v * 1e9 for v in d["t_s"]]
        v = [x / 1e3 for x in d["v_out"]]
        ax.plot(t, v, color=col, linewidth=2)
        k = max(range(len(v)), key=v.__getitem__)
        ax.text(t[k] + 12, v[k], f"{label} {v[k]:.0f} kV", color=INK, fontsize=9, va="center")
    ax.set_xlim(0, 600)
    ax.set_ylim(0, None)
    ax.set_xlabel("time after the gaps fire (ns)", color=INK2)
    ax.set_ylabel("output into 50 Ω (kV)", color=INK2)
    style(ax, "Four stages give 3.06× the charge voltage; the gas sets the charge voltage",
          "Marx generator, 4 × 10 nF, 0.2 mm gaps at 500 psi charged to 90 % of breakdown, 1 µH loop (kronspice)")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "marx.png"), dpi=150)


def pump():
    d = read("pump_argon.csv")
    t = [x * 1e3 for x in d["t_s"]]
    fig, (a, b) = plt.subplots(1, 2, figsize=(10, 4.6), gridspec_kw={"width_ratios": [3, 2]})
    a.plot(t, d["v_out"], color="#2a78d6", linewidth=2)
    a.axhline(300, color=INK2, linewidth=1, linestyle=(0, (4, 3)))
    a.text(0.3, 306, "ideal: 3 × 100 V", color=INK2, fontsize=9)
    a.set_xlabel("time (ms)", color=INK2)
    a.set_ylabel("output (V)", color=INK2)
    a.set_ylim(0, 330)
    style(a, "Plasma varactor pumps 100 V up to 295 V",
          "Argon gap, C ×3 when the plasma is lit, 10 kHz, 1 MΩ load (kronspice)")
    t0 = t[-1] - 0.3
    k0 = next(i for i, x in enumerate(t) if x >= t0)
    tt = [x - t0 for x in t[k0:]]
    b.plot(tt, d["v_varactor"][k0:], color="#eb6834", linewidth=2)
    b.plot(tt, d["v_out"][k0:], color="#2a78d6", linewidth=2)
    b.text(0.01, 75, "varactor", color=INK, fontsize=9)
    b.text(0.01, 280, "output", color=INK, fontsize=9)
    b.set_xlabel("last 0.3 ms (ms)", color=INK2)
    b.set_ylim(0, 330)
    style(b, "Each cycle: fill at 100 V, lift to ~300 V", "C high (plasma on) → C low (off)")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "pump.png"), dpi=150)


if __name__ == "__main__":
    paschen()
    marx()
    pump()
