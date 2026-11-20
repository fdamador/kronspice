# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Report figures for the three-coil example from the run's CSV data.

    python3 examples/three_coil_sync/viz/plot_report.py examples/three_coil_sync/viz/data examples/three_coil_sync/viz

Writes convergence.png, compare.png, wave.png, trends.png, losses.png and startup.png.
"""
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HVC, ROY, OUT, GOLD, INK, MUTED, GRID = "#2a78d6", "#eb6834", "#1baf7a", "#6250d6", "#0b0b0b", "#52514e", "#e4e3df"
data, out = sys.argv[1], sys.argv[2]
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
                     "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 1.8,
                     "legend.frameon": False})


def load(name):
    with open(os.path.join(data, name)) as f:
        rows = list(csv.DictReader(f))
    outd = {}
    for k in rows[0]:
        try:
            outd[k] = [float(r[k]) for r in rows]
        except ValueError:
            outd[k] = [r[k] for r in rows]
    return outd


def save(fig, name):
    fig.savefig(os.path.join(out, name), dpi=160, bbox_inches="tight")
    plt.close(fig)


# convergence
c = load("convergence.csv")
fig, ax = plt.subplots(figsize=(7, 3.4))
ax.plot(c["gen"], [100 / b - 100 for b in c["best"]], color=OUT, marker="o", ms=4, label="best (rating-penalized) power")
ax.set_xlabel("generation")
ax.set_ylabel("score as output power [W]")
ax.set_title("GA convergence: 100/fitness − 100 (equals P_out when no rating is exceeded)", fontsize=9, loc="left")
ax.set_ylim(bottom=0)
save(fig, "convergence.png")

# compare
cm = load("compare.csv")
names = cm["case"]
fig, axs = plt.subplots(1, 2, figsize=(11, 3.8), sharey=False)
for ax, tag in zip(axs, ("evolved", "baseline")):
    idx = [i for i, n in enumerate(names) if n.startswith(tag)]
    labels = [names[i].split(" ", 1)[1] for i in idx]
    x = range(len(idx))
    w = 0.26
    ax.bar([i - w - 0.02 for i in x], [cm["p_hv"][i] for i in idx], w, color=HVC, label="HV generator in")
    ax.bar([i for i in x], [cm["p_roy"][i] for i in idx], w, color=ROY, label="Royer supply in")
    ax.bar([i + w + 0.02 for i in x], [cm["p_out"][i] for i in idx], w, color=OUT, label="load out")
    for i, j in zip(x, idx):
        ax.text(i + w + 0.02, cm["p_out"][j], f"{cm['p_out'][j]:.1f} W", ha="center", va="bottom", fontsize=8, color=INK)
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels)
    ax.set_ylabel("mean power [W]")
    ax.axhline(0, color=MUTED, lw=0.8)
    ax.set_title(f"{tag} design", fontsize=10, loc="left")
fig.legend(*axs[0].get_legend_handles_labels(), fontsize=8, loc="lower center", ncol=3, bbox_to_anchor=(0.5, 1.0))
save(fig, "compare.png")

# wave
w = load("wave.csv")
t = w["t_us"]
fig, axs = plt.subplots(3, 1, figsize=(10, 7.5), sharex=True)
ax = axs[0]
ax.plot(t, w["v_hv"], color=HVC, label="HV terminal (coil 1 → coil 3)")
ax.plot(t, w["v_tank"], color=ROY, label="Royer tank c1 − c2")
ax.plot(t, w["v_out"], color=OUT, label="load voltage")
ax.set_ylabel("voltage [V]")
ax.legend(fontsize=8, loc="lower left", bbox_to_anchor=(0, 1.0), ncol=3)
ax = axs[1]
ax.plot(t, w["i_coil3"], color=OUT, label="coil 3 current")
ax.plot(t, w["i_tank"], color=ROY, label="coil 2 tank current (½(i2a + i2b))")
ax.plot(t, w["i_hvgen"], color=HVC, label="HV generator / amplifier current")
ax.set_ylabel("current [A]")
ax.legend(fontsize=8, loc="lower left", bbox_to_anchor=(0, 1.0), ncol=3)
ax = axs[2]
ax.plot(t, w["v_c3"], color=GOLD, label="across C3")
ax.set_ylabel("voltage [V]")
ax.set_xlabel("time in the last five HV periods [µs]")
ax.legend(fontsize=8, loc="upper right")
save(fig, "wave.png")

# trends
tr = load("trends.csv")
labels = {"detune": "HV frequency / evolved value", "gain": "HV amplifier gain / evolved gain",
          "load": "load / evolved load", "c3": "C3 / evolved C3", "gap": "coil 3 → coil 2 spacer [mm]"}
sync = "gain" in tr["sweep"]
labels["phase"] = ("HV drive phase relative to the Royer tank [°]" if sync else
                   "HV generator phase at switch-on [°]")
order = ["gain" if sync else "detune", "phase", "load", "c3", "gap"]
sweeps = [(k, labels[k]) for k in order]
fig, axs = plt.subplots(2, 3, figsize=(13, 7))
for ax, (s, xl) in zip(axs.flat, sweeps):
    idx = [i for i, n in enumerate(tr["sweep"]) if n == s]
    x = [tr["x"][i] for i in idx]
    ax.plot(x, [tr["p_out"][i] for i in idx], color=OUT, marker="o", ms=3.5, label="load out")
    ax.plot(x, [tr["p_hv"][i] for i in idx], color=HVC, lw=1.2, label="HV in")
    ax.plot(x, [tr["p_roy"][i] for i in idx], color=ROY, lw=1.2, label="Royer in")
    lk = [x[j] for j, i in enumerate(idx) if tr["locked"][i] == "true"]
    if lk and not sync:
        ax.axvspan(min(lk), max(lk), color=OUT, alpha=0.08, lw=0)
    if s in ("load", "c3"):
        ax.set_xscale("log")
    ax.set_xlabel(xl)
    ax.set_ylabel("mean power [W]")
axs[0, 0].legend(fontsize=8, loc="best")
ax = axs[1, 2]
idx = [i for i, n in enumerate(tr["sweep"]) if n == "gap"]
ax.plot([tr["x"][i] for i in idx], [tr["k23"][i] for i in idx], color=MUTED, marker="o", ms=3.5)
ax.set_xlabel("coil 3 → coil 2 spacer [mm]")
ax.set_ylabel("coupling k(2–3)")
fig.suptitle("Trends around the evolved design" + ("" if sync else " (shaded: Royer locked to the HV generator)"),
             fontsize=10, x=0.02, ha="left")
fig.tight_layout()
save(fig, "trends.png")

# losses
with open(os.path.join(data, "design.csv")) as f:
    dz = {r["name"]: float(r["value"]) for r in csv.DictReader(f)}
groups = [("load (useful)", dz["POut"], OUT),
          ("coil 2 copper", dz["loss_R2a"] + dz["loss_R2b"], ROY),
          ("tank capacitor ESR", dz["loss_RC2"], ROY),
          ("MOSFETs", dz["loss_M1"] + dz["loss_M2"], ROY),
          ("gate drive + diodes", dz["loss_RG1"] + dz["loss_RG2"] + dz["loss_D1"] + dz["loss_D2"], ROY),
          ("supply + choke", dz["loss_Rs"], ROY),
          ("HV generator source R", dz["loss_Rg"], HVC),
          ("coil 1 copper", dz["loss_R1p"] + dz["loss_R1s"], HVC),
          ("coil 3 copper", dz["loss_R3"], OUT)]
fig, ax = plt.subplots(figsize=(8, 3.8))
ys = range(len(groups))
ax.barh(list(ys), [g[1] for g in groups], color=[g[2] for g in groups], height=0.6)
for y, g in zip(ys, groups):
    ax.text(g[1], y, f"  {g[1]:.2f} W", va="center", fontsize=8, color=INK)
ax.set_yticks(list(ys))
ax.set_yticklabels([g[0] for g in groups])
ax.invert_yaxis()
ax.set_xlabel("mean power [W]")
ax.set_title(f"Where the input goes: HV {dz['PInHV']:.1f} W + Royer {dz['PInRoy']:.1f} W in, "
             f"η = {dz['Eta']:.3f}", fontsize=9, loc="left")
save(fig, "losses.png")

# startup
s = load("startup.csv")
fig, ax = plt.subplots(figsize=(8, 3.2))
ax.plot(s["t_ms"], s["v_tank_peak"], color=ROY, label="Royer tank peak |c1 − c2|")
ax.plot(s["t_ms"], s["v_out_peak"], color=OUT, label="load peak voltage")
ax.axvspan(2.0, s["t_ms"][-1], color=MUTED, alpha=0.08, lw=0)
ax.text(2.02, ax.get_ylim()[1] * 0.92, "measuring window", fontsize=8, color=MUTED)
ax.set_xlabel("time from switch-on [ms]")
ax.set_ylabel("peak voltage per 10 µs [V]")
ax.legend(fontsize=8, loc="lower right")
save(fig, "startup.png")
