# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Chart for Part 7 (H40) of examples/fdamador_neutrino: a spinning ring.

    go run ./examples/fdamador_neutrino -only 7 -field -viz examples/fdamador_neutrino/viz
    python3 examples/fdamador_neutrino/viz/render_part7.py

Writes
  viz/spin.png   7a: the CP phase and the mixing angles a rotating ring line reaches
                 against its non-reciprocity (with the rim speed it needs and the
                 burst limit), and P(e->mu) against P(mu->e) in kronspice with and
                 without rotation; 7b: angular momentum, helicity and energy of the
                 spinning fields in the field solver

Every curve is read from the CSV files the Go program writes (-viz); nothing is
computed here except plotting (and the rim speed nu*c/2 for the top axis).
"""
import csv
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
BLUE, ORANGE, GREEN, RED, PURPLE, GREY = "#2a78d6", "#eb6834", "#2f9e58", "#c93c3c", "#7d4fc4", "#8a8985"
plt.rcParams.update({"font.size": 13, "axes.titlesize": 14.5, "axes.titleweight": "normal",
                     "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb",
                     "savefig.facecolor": "#fcfcfb", "axes.spines.top": False, "axes.spines.right": False,
                     "legend.frameon": False, "legend.fontsize": 11})
C0 = 299792458.0


def load(name):
    with open(os.path.join(HERE, name)) as f:
        rows = list(csv.reader(f))[1:]
    return np.array([[float(x) for x in r] for r in rows])


def exists(name):
    return os.path.exists(os.path.join(HERE, name))


def grid(ax):
    ax.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)


meta = {}
with open(os.path.join(HERE, "p7_meta.csv")) as f:
    for k, v in list(csv.reader(f))[1:]:
        meta[k] = float(v)

have7b = exists("p7b_runs.csv")
fig = plt.figure(figsize=(14, 19 if have7b else 12.5))
rows = 3 if have7b else 2
gs = fig.add_gridspec(rows, 2, hspace=0.42, wspace=0.26, left=0.07, right=0.97, top=0.93, bottom=0.05)
fig.suptitle("fdamador neutrino, Part 7 (H40): a spinning ring", fontsize=18, color=INK, y=0.985)

# ---- (a) delta_CP and J/Jmax against |nu| (rotation sense that gives delta > 180)
s = load("p7_scan.csv")
neg = s[s[:, 0] < 0]
neg = neg[np.argsort(-neg[:, 0])]
x = -neg[:, 0]
ax = fig.add_subplot(gs[0, 0])
grid(ax)
lo, hi = meta["delta"] - 41, meta["delta"] + 26
ax.axhspan(lo, hi, color=BLUE, alpha=0.10, lw=0)
ax.axhline(meta["delta"], color=BLUE, lw=1, ls="--")
ax.text(4e-4, 176, "measured 212° (+26/−41)", color=BLUE, fontsize=11)
keep = x <= 5e-3  # beyond, the cluster's modes re-order and the labels jump
ax.plot(x[keep], neg[keep, 1], color=INK, lw=2, label="δ_CP (Part 6 defects)")
ax.axvline(meta["burst_nu"], color=RED, lw=1.3)
ax.text(meta["burst_nu"] * 1.15, 232, "burst speed\n(T1000G, v_g = c)", color=RED, fontsize=11)
ax.plot([-meta["fit_nu"]], [meta["fit_delta"]], "o", color=GREEN, ms=8, label="fit, δ free in 1σ band")
ax.plot([-meta["pinned_nu"]], [meta["pinned_delta"]], "s", color=ORANGE, ms=8, label="fit, δ pinned at 212°")
ax.set_xscale("log")
ax.set_xlim(1e-6, 5e-3)
ax.set_ylim(172, 252)
ax.set_xlabel("non-reciprocity ν = Δω/ω of the pair")
ax.set_ylabel("δ_CP [°]")
ax.set_title("7a  CP phase from rotation")
ax.legend(loc="upper left", fontsize=10.5)
top = ax.secondary_xaxis("top", functions=(lambda v: v * C0 / 2, lambda r: 2 * r / C0))
top.set_xlabel("rim speed aΩ needed with v_g = c [m/s]", fontsize=11.5)

# ---- (b) angles and ratio, deviation in measured sigma
ax = fig.add_subplot(gs[0, 1])
grid(ax)
for col, key, sk, name, c in [(3, "th12", "s12", "θ12", BLUE), (4, "th23", "s23", "θ23", GREEN), (5, "th13", "s13", "θ13", ORANGE)]:
    ax.plot(x, (neg[:, col] - meta[key]) / meta[sk], color=c, lw=2, label=name)
ax.plot(x, np.log(neg[:, 6] / meta["ratio"]) / meta["sratio"], color=PURPLE, lw=2, label="splitting ratio")
ax.axhspan(-1, 1, color=GREY, alpha=0.12, lw=0)
ax.axvline(meta["burst_nu"], color=RED, lw=1.3)
ax.set_xscale("log")
ax.set_xlim(1e-6, 3e-2)
ax.set_yscale("symlog", linthresh=1)
ax.set_xlabel("non-reciprocity ν")
ax.set_ylabel("deviation from measured [σ]")
ax.set_title("7a  what rotation does to the angles (Part 6 defects)")
ax.legend(loc="upper left", fontsize=10.5)

# ---- (c) T violation in kronspice
ax = fig.add_subplot(gs[1, 0])
grid(ax)
t = load("p7_tviol.csv")
z = load("p7_tviol_zero.csv")
ax.plot(t[:, 0], t[:, 1], color=BLUE, lw=1.4, label="P(e→μ), rotating")
ax.plot(t[:, 0], t[:, 2], color=ORANGE, lw=1.4, label="P(μ→e), rotating")
ax.set_xlabel("time [periods]")
ax.set_ylabel("probability (kronspice)")
ax.set_title("7a  δ_CP ≈ 211°: P(e→μ) ≠ P(μ→e)")
ax.legend(loc="upper right", fontsize=10.5)
ax = fig.add_subplot(gs[1, 1])
grid(ax)
ax.plot(t[:, 0], t[:, 1] - t[:, 2], color=INK, lw=1.4, label="rotating (kronspice)")
ax.plot(t[:, 0], t[:, 3] - t[:, 4], color=RED, lw=1, ls="--", label="rotating (exact modes)")
ax.plot(z[:, 0], z[:, 1] - z[:, 2], color=GREY, lw=1.2, label="no rotation (readout only)")
ax.set_xlabel("time [periods]")
ax.set_ylabel("P(e→μ) − P(μ→e)")
ax.set_title("7a  the T-odd difference")
ax.legend(loc="lower left", fontsize=10.5)

if have7b:
    names = {}
    with open(os.path.join(HERE, "p7b_runs.csv"), encoding="utf-8") as f:
        for line in list(f)[1:]:
            i, rest = line.strip().split(",", 1)
            parts = rest.rsplit(",", 3)
            names[int(i)] = parts[0]
    runs = {i: load(f"p7b_run{i}.csv") for i in names}
    cells = load("p7b_grid.csv").ravel()[0] if exists("p7b_grid.csv") else 2
    gl = f"r/{cells:g}"

    def upto(r, frac=1e-3):
        k = np.where(r[:, 1] >= frac)[0]
        return r[: k[-1] + 1] if len(k) else r

    ax = fig.add_subplot(gs[2, 0])
    grid(ax)
    for i, c, ls in [(6, BLUE, "-"), (7, ORANGE, "-"), (5, GREY, "-"), (2, GREEN, "--"), (8, PURPLE, "-")]:
        r = upto(runs[i])
        ax.plot(r[:, 0], r[:, 5], color=c, lw=1.8, ls=ls, label=names[i])
    r = upto(runs[6])
    ax.plot(r[:, 0], r[:, 2], color=BLUE, lw=1.4, ls=":", label="travelling (+), Abraham ε0E×B")
    for y in (1, -1):
        ax.axhline(y, color=INK2, lw=0.8, ls=":")
    ax.set_xlabel("time [periods of the E-ring]")
    ax.set_ylabel("J_z·ω / U  (Minkowski D×B)")
    ax.set_ylim(-2.4, 2.4)
    ax.set_title(f"7b  angular momentum about the axis ({gl})")
    ax.legend(loc="upper left", fontsize=9.5, ncol=2)

    ax = fig.add_subplot(gs[2, 1])
    grid(ax)
    for i, c in [(2, BLUE), (3, ORANGE), (6, GREEN)]:
        r = runs[i]
        r = r[r[:, 0] <= 3]
        ax.plot(r[:, 0], r[:, 3], color=c, lw=1.8, label=names[i])
    ax.set_xlabel("time [periods of the E-ring], first 3")
    ax.set_ylabel("h·ω / U")
    ax.set_ylim(-6, 8)
    ax.set_title(f"7b  helicity: sign swings as E and B rings beat ({gl})")
    ax.legend(loc="upper right", fontsize=9.5)

out = os.path.join(HERE, "spin.png")
fig.savefig(out, dpi=100)
print("wrote", out)
