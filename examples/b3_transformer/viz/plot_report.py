# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""B(3) transformer report figures from the run's data.

    go run ./examples/b3_transformer -data examples/b3_transformer/viz/data
    python3 examples/b3_transformer/viz/plot_report.py examples/b3_transformer/viz/data examples/b3_transformer/viz

Writes wave.png (input and output, volts per turn, dΦ/dt per pass, against a
normal 100:100 transformer), freq.png (frequency response), ratio.png (the
turns ratio: spring against normal) and passes.png (dΦ/dt along the spring).
"""
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PRI, SEC, AX, NORM, CORE = "#3b6ea5", "#c8553d", "#2a9d8f", "0.45", "#b08d57"
data, out = sys.argv[1], sys.argv[2]


def load(name):
    with open(os.path.join(data, name)) as f:
        rows = list(csv.DictReader(f))
    return {k: [float(r[k]) for r in rows] for k in rows[0]}


# ---------- wave.png ----------
w = load("wave.csv")
t = w["t_ms"]
fig, ax = plt.subplots(3, 2, figsize=(13, 10), sharex=True)
fig.suptitle("Input and output at 1 kHz — spring core (left) against a normal 100:100 transformer (right)", fontsize=13)
a = ax[0, 0]
a.plot(t, w["v_in"], color=PRI, label="primary v_p (100 turns × 10 passes)")
a.plot(t, w["v_out"], color=SEC, label="secondary v_s (100 turns × 2 passes)")
a.plot(t, w["v_axis"], color=AX, label="bore wire (1 turn × 10 passes)")
a.set_ylabel("voltage [V]"); a.legend(fontsize=8, loc="upper right"); a.grid(alpha=0.3)
a.set_title("spring core: v_s / v_p = (100·2)/(100·10) = 0.2", fontsize=10)
a = ax[0, 1]
a.plot(t, w["normal_v_in"], color=PRI, label="primary v_p (100 turns)")
a.plot(t, w["normal_v_out"], color=SEC, ls="--", label="secondary v_s (100 turns)")
a.legend(fontsize=8, loc="upper right"); a.grid(alpha=0.3)
a.set_title("normal transformer: v_s / v_p = 100/100 = 1", fontsize=10)
a = ax[1, 0]
a.plot(t, w["vpt_in"], color=PRI, label="primary: volts per winding turn = 10·dΦ/dt")
a.plot(t, w["vpt_out"], color=SEC, label="secondary: volts per winding turn = 2·dΦ/dt")
a.set_ylabel("volts per winding turn [V]"); a.legend(fontsize=8, loc="upper right"); a.grid(alpha=0.3)
a.set_title("volts per turn differ by n_p/n_s = 5", fontsize=10)
a = ax[1, 1]
a.plot(t, w["normal_vpt_in"], color=PRI, label="primary: volts per turn = dΦ/dt")
a.plot(t, w["normal_vpt_out"], color=SEC, ls="--", label="secondary: volts per turn = dΦ/dt")
a.legend(fontsize=8, loc="upper right"); a.grid(alpha=0.3)
a.set_title("volts per turn equal in both windings", fontsize=10)
a = ax[2, 0]
a.plot(t, w["dphi_pass1"], color=CORE, lw=3, label="pass 1")
a.plot(t, w["dphi_pass6"], color="k", lw=1.2, ls="--", label="pass 6 (inside the secondary)")
a.plot(t, w["dphi_pass10"], color=PRI, lw=1, ls=":", label="pass 10")
a.set_ylabel("dΦ/dt in the core [Wb/s]"); a.set_xlabel("time [ms]"); a.legend(fontsize=8, loc="upper right"); a.grid(alpha=0.3)
a.set_title("the rate of change is the same in every core pass", fontsize=10)
a = ax[2, 1]
a.plot(t, w["i_in"], color=PRI, label="primary current")
a.plot(t, w["i_out"], color=SEC, label="secondary current")
a.set_ylabel("current [A]"); a.set_xlabel("time [ms]"); a.legend(fontsize=8, loc="upper right"); a.grid(alpha=0.3)
a.set_title("spring core: winding currents (secondary into its load)", fontsize=10)
fig.tight_layout(rect=(0, 0, 1, 0.96))
fig.savefig(os.path.join(out, "wave.png"), dpi=120)

# ---------- freq.png ----------
f = load("freq.csv")
fig, ax = plt.subplots(2, 2, figsize=(13, 8.5))
fig.suptitle("Frequency response (loaded, AC analysis in kronspice)", fontsize=13)
cores = [("nano", "nanocrystalline tape (design)", PRI), ("steel", "silicon-steel wire", CORE), ("galfenol", "Galfenol wire", SEC)]
for tag, name, col in cores:
    ax[0, 0].semilogx(f["f"], f["ratio_" + tag], color=col, label=name)
    ax[0, 1].semilogx(f["f"], f["phase_" + tag], color=col, label=name)
    ax[1, 0].loglog(f["f"], f["zin_" + tag], color=col, label=name)
    ax[1, 1].semilogx(f["f"], f["eta_" + tag], color=col, label=name)
ax[0, 0].axhline(0.2, color="k", ls=":", lw=1); ax[0, 0].text(1.2, 0.205, "ideal (N_s·n_s)/(N_p·n_p) = 0.2", fontsize=8)
ax[0, 0].set_ylabel("|v_s / v_p|"); ax[0, 0].set_ylim(0, 0.24)
ax[0, 1].set_ylabel("phase of v_s / v_p [°]")
ax[1, 0].set_ylabel("|Z_in| at the primary [Ω]")
ax[1, 1].set_ylabel("efficiency P_out / P_in"); ax[1, 1].set_ylim(0, 1.02)
for a in ax.flat:
    a.axvline(1e3, color="0.6", ls="--", lw=0.8)
    a.set_xlabel("frequency [Hz]"); a.grid(alpha=0.3, which="both"); a.legend(fontsize=8)
fig.tight_layout(rect=(0, 0, 1, 0.95))
fig.savefig(os.path.join(out, "freq.png"), dpi=120)

# ---------- ratio.png ----------
r = load("ratio.csv")
fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
fig.suptitle("The turns ratio: spring core against a normal transformer (100-turn windings)", fontsize=13)
i0 = [i for i, s in enumerate(r["sweep"]) if s == 0]
i1 = [i for i, s in enumerate(r["sweep"]) if s == 1]
a = ax[0]
a.plot([r["x"][i] for i in i0], [r["normal"][i] for i in i0], "s--", color=NORM, label="normal: N_s/N_p = 1")
a.plot([r["x"][i] for i in i0], [r["ideal"][i] for i in i0], "-", color="k", lw=1, label="(N_s·n_s)/(N_p·n_p)")
a.plot([r["x"][i] for i in i0], [r["ratio"][i] for i in i0], "o", color=PRI, label="spring core (network)")
a.set_xlabel("core passes enclosed by the secondary, n_s (primary round all 10)"); a.set_ylabel("v_s / v_p")
a.grid(alpha=0.3); a.legend(fontsize=9)
a = ax[1]
a.loglog([r["x"][i] for i in i1], [r["normal"][i] for i in i1], "s--", color=NORM, label="normal: N_s/N_p = 1")
a.loglog([r["x"][i] for i in i1], [r["ideal"][i] for i in i1], "-", color="k", lw=1, label="2/N")
a.loglog([r["x"][i] for i in i1], [r["ratio"][i] for i in i1], "o", color=PRI, label="spring core (network)")
a.set_xlabel("core turns N (primary round all N, secondary round 2)"); a.set_ylabel("v_s / v_p")
a.grid(alpha=0.3, which="both"); a.legend(fontsize=9)
fig.tight_layout(rect=(0, 0, 1, 0.93))
fig.savefig(os.path.join(out, "ratio.png"), dpi=120)

# ---------- passes.png ----------
p = load("passes.csv")
fig, a = plt.subplots(figsize=(9, 4.2))
for k, name, col, mk in [("nano_design", "nanocrystalline, design (a 2 mm, gap 8 mm, return 16×)", PRI, "o"),
                         ("steel_design", "silicon steel, design geometry", CORE, "s"),
                         ("nano_plain", "nanocrystalline, plain (a 1 mm, gap 1 mm, return 1×)", AX, "^"),
                         ("galfenol_design", "Galfenol, design geometry", SEC, "v")]:
    a.plot(p["pass"], p[k], mk + "-", color=col, label=name)
a.set_xlabel("core pass (turn of the spring)"); a.set_ylabel("dΦ/dt in the pass / largest")
a.set_title("Rate of change along the spring, primary round all 10 passes", fontsize=11)
a.grid(alpha=0.3); a.legend(fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(out, "passes.png"), dpi=120)
print("wrote wave.png, freq.png, ratio.png, passes.png to", out)

# ---------- ledger.png ----------
L = load("ledger.csv")
n = L["N"]
fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
fig.suptitle("Energy ledger against core turns N (20 V peak, 1 kHz, 20 cycles, fixed 632 Ω load)", fontsize=13)
a = ax[0]
mj = lambda k: [v * 1e3 for v in L[k]]
a.loglog(n, mj("e_in"), "o-", color=PRI, label="E in (delivered by the source)")
a.loglog(n, mj("e_back"), "s--", color=NORM, label="E back (returned by the magnetizing current)")
a.loglog(n, mj("e_net"), "^-", color="k", label="net in = E in − E back")
a.loglog(n, mj("e_load"), "v:", color=SEC, label="E load")
a.loglog(n, mj("e_copper"), "x-", color=CORE, label="E copper (windings + source R)")
a.set_xlabel("core turns N (primary round all N passes, secondary round 2)"); a.set_ylabel("energy over the run [mJ]")
a.grid(alpha=0.3, which="both"); a.legend(fontsize=8)
a = ax[1]
a.semilogx(n, [o / i for o, i in zip(L["p_out"], L["p_in"])], "o-", color=PRI, label="efficiency P_out / P_in")
a.semilogx(n, [1e6 * abs(r) / e for r, e in zip(L["residual"], L["e_in"])], "s--", color=SEC, label="ledger |residual| / E in (×1e6)")
a.set_ylim(0, 1.05); a.set_xlabel("core turns N"); a.grid(alpha=0.3, which="both")
tw = a.twinx()
tw.loglog(n, L["dphi_pk"], "^:", color=AX, label="peak dΦ/dt per pass [Wb/s]")
tw.set_ylabel("peak dΦ/dt per pass [Wb/s]", color=AX)
h1, l1 = a.get_legend_handles_labels(); h2, l2 = tw.get_legend_handles_labels()
a.legend(h1 + h2, l1 + l2, fontsize=8, loc="center left")
fig.tight_layout(rect=(0, 0, 1, 0.93))
fig.savefig(os.path.join(out, "ledger.png"), dpi=120)
print("wrote ledger.png")


def plain_ticks(axes, ticks):
    for a in axes:
        a.set_xticks(ticks)
        a.set_xticklabels([str(t) for t in ticks])
        a.minorticks_off()


plain_ticks([ax[0], ax[1]], [2, 4, 6, 10, 16, 20, 30, 40])
fig.savefig(os.path.join(out, "ledger.png"), dpi=120)
