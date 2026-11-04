# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Figures 4-9 of the earth_machine study, from the CSVs written by

    go run ./examples/earth_machine -data examples/earth_machine/viz/data

Usage: python viz/plot_layers.py [data_dir] [out_dir]   (matplotlib, numpy)
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
MEAS = [7.8, 14.3, 20.8, 27.3, 33.8]

plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.3, "savefig.dpi": 300,
                     "axes.spines.top": False, "axes.spines.right": False})


def load(name):
    with open(os.path.join(data, name), newline="") as f:
        r = csv.reader(f)
        head = next(r)
        rows = np.array([[float(x) for x in row] for row in r])
    return {h: rows[:, i] for i, h in enumerate(head)}


def save(fig, name):
    fig.tight_layout()
    fig.savefig(os.path.join(out, name))
    plt.close(fig)
    print("wrote", name)


# Fig. 4 — layer 2: the leaky atmosphere.
p, q, s = load("l2_profile.csv"), load("l2_charge.csv"), load("l2_switchon.csv")
fig, ax = plt.subplots(1, 3, figsize=(10, 3.4))
ax[0].semilogx(p["sigma_exp"], p["z"] / 1e3, label="one scale, H = 2.08 km")
ax[0].semilogx(p["sigma_two"], p["z"] / 1e3, "--", label="two-scale (BL 1.5 km, 5 km above)")
ax[0].set_xlabel("conductivity σ [S/m]"); ax[0].set_ylabel("height [km]"); ax[0].legend(fontsize=7)
ax[0].set_title("(a) calibrated profiles")
ax[1].plot(q["cum"] / q["cum"][-1], q["z"] / 1e3)
ax[1].axvline(0.5, color="gray", lw=0.6)
ax[1].set_ylim(0, 10); ax[1].set_xlabel("fraction of the positive space charge below z")
ax[1].set_ylabel("height [km]"); ax[1].set_title("(b) where the Earth's balancing charge sits")
ax[2].plot(s["t"], s["vi_circuit"] / 1e3, label="kronspice transient")
ax[2].plot(s["t"], s["vi_exact"] / 1e3, ":", lw=2, label="Σ I·R_k(1 − e^(−t/R_kC_k))")
ax[2].set_xlabel("time after storms switch on [s]"); ax[2].set_ylabel("V_I [kV]")
ax[2].set_title("(c) charging the ionosphere"); ax[2].legend(fontsize=7)
save(fig, "fig4_global_circuit.png")

# Fig. 5 — layer 3: Schumann resonances.
sp, md, cv, tw = load("l3_spectra.csv"), load("l3_modes.csv"), load("l3_convergence.csv"), load("l3_two.csv")
fig, ax = plt.subplots(1, 3, figsize=(10, 3.4))
for k in range(1, 6):
    ax[0].loglog(180 / cv["bands"], np.abs(cv["e%d" % k]), "o-", ms=3, label="n = %d" % k)
ax[0].loglog(180 / cv["bands"], 3e-4 * (180 / cv["bands"]) ** 2 / 1.0, "k:", lw=0.8, label="∝ Δθ²")
ax[0].set_xlabel("band width Δθ [°]"); ax[0].set_ylabel("|ladder / Legendre − 1|")
ax[0].set_title("(a) perfect cavity = Legendre"); ax[0].legend(fontsize=6)
ax[0].set_xticks([0.5, 1, 2, 6]); ax[0].set_xticklabels(["0.5", "1", "2", "6"]); ax[0].minorticks_off()
norm = lambda y: y / np.nanmax(y)
ax[1].plot(sp["f"], norm(sp["const_pole"]), label="constant heights, pole")
ax[1].plot(sp["f"], norm(sp["const_equator"]), lw=0.8, label="constant heights, equator")
ax[1].plot(tw["f"], norm(tw["pole"]), "--", color="C2", label="two-height, pole (rebuilt per f)")
for m in MEAS:
    ax[1].axvline(m, color="gray", lw=0.6, ls=":")
ax[1].set_xlabel("frequency [Hz]"); ax[1].set_ylabel("|V| (normalized)"); ax[1].set_xlim(3, 42)
ax[1].set_title("(b) spectra; dotted: measured modes"); ax[1].legend(fontsize=6)
n = md["n"]
ax[2].plot(n, md["const_f"] / md["measured"], "s-", label="constant heights (fit to n = 1)")
ax[2].plot(n, md["two_f"] / md["measured"], "o-", label="two-height (fit)")
ax[2].axhline(1, color="gray", lw=0.6)
ax2 = ax[2].twinx()
ax2.plot(n, md["const_q"], "s:", color="C0", alpha=0.6)
ax2.plot(n, md["two_q"], "o:", color="C1", alpha=0.6)
ax2.axhspan(4, 6, color="gray", alpha=0.12)
ax2.set_ylabel("Q (dotted; band: reported 4–6)")
ax[2].set_xlabel("mode n"); ax[2].set_ylabel("model / measured frequency")
ax[2].set_title("(c) modes against measurement"); ax[2].legend(fontsize=6, loc="center left")
save(fig, "fig5_schumann.png")

# Fig. 6 — layer 1 statistics.
sv, mua = load("l1_survival.csv"), load("l1_mua.csv")
fig, ax = plt.subplots(1, 2, figsize=(8, 3.4))
m = np.mean(sv["interval"])
ax[0].semilogy(sv["interval"] / m, sv["survival"], label="Rikitake μ = 1, A = 1 (RK4, t = 20 000)")
x = np.linspace(0, 6, 100)
ax[0].semilogy(x, np.exp(-x), "k:", label="exponential (Poisson reversals)")
ax[0].set_xlabel("polarity interval / mean"); ax[0].set_ylabel("fraction longer")
ax[0].set_title("(a) intervals: a floor, and no superchrons"); ax[0].legend(fontsize=7)
for mu in (0.5, 1, 2):
    sel = mua["mu"] == mu
    ax[1].plot(mua["A"][sel], mua["maxovermean"][sel], "o-", label="μ = %g" % mu)
ax[1].axhline(100, color="gray", ls="--", lw=0.8)
ax[1].text(0.6, 60, "polarity record: longest/mean ≳ 100", fontsize=7)
ax[1].set_yscale("log"); ax[1].set_xlabel("A = ω1 − ω2"); ax[1].set_ylabel("longest / mean interval")
ax[1].set_title("(b) the (μ, A) plane"); ax[1].legend(fontsize=7)
save(fig, "fig6_reversal_statistics.png")

# Fig. 7 — paper107 expansions.
dr, fb = load("l1_drive.csv"), load("l1_feedback.csv")
fig, ax = plt.subplots(1, 3, figsize=(10, 3.4))
ax[0].plot(dr["wd_over_wn"], dr["rate"], "o-", ms=3)
ax[0].axvline(1, color="gray", lw=0.6)
ax[0].set_xlabel("drive modulation ωd / ωn"); ax[0].set_ylabel("reversals per unit time")
ax[0].set_title("(a) torque 1 + 0.3 sin ωd·t")
w = np.linspace(-3, 3, 241)
W1, W2 = np.meshgrid(w, w)
lam = np.where(W1 * W2 > 0, -1 + np.sqrt(np.abs(W1 * W2)), -1.0)
im = ax[1].contourf(W1, W2, lam, levels=np.linspace(-1, 2, 13), cmap="RdBu_r")
ax[1].contour(W1, W2, lam, levels=[0], colors="k", linewidths=1)
fig.colorbar(im, ax=ax[1], label="growth rate of the field")
ax[1].set_xlabel("ω1 (speed drive)"); ax[1].set_ylabel("ω2 (speed drive)")
ax[1].set_title("(b) speed drives: grows if ω1ω2 > 1")
ax[2].plot(fb["K"], fb["reversals"], "o-")
ax[2].set_yscale("symlog", linthresh=1)
ax[2].set_xlabel("feedback coil K (aiding > 0)"); ax[2].set_ylabel("reversals in t = 20 000")
ax[2].set_title("(c) feedback coil")
save(fig, "fig7_paper107_studies.png")

# Fig. 8 — the layers linked.
un, dp = load("link_unified.csv"), load("link_dipole.csv")
w = dp["t"] <= 220
dp = {k: v[w] for k, v in dp.items()}
fig, ax = plt.subplots(1, 2, figsize=(8, 3.4))
ax[0].plot(un["f"], un["pole"] / un["pole"].max(), label="stroke point (pole)")
ax[0].plot(un["f"], un["equator"] / un["equator"].max(), label="equator")
for m in MEAS:
    ax[0].axvline(m, color="gray", lw=0.6, ls=":")
ax[0].set_xlabel("frequency [Hz]"); ax[0].set_ylabel("|V| (normalized)")
ax[0].set_title("(a) one network: columns + cavity (no fitted h_e)"); ax[0].legend(fontsize=7)
ax[1].plot(dp["t"], dp["I1"], lw=0.7, label="dipole ∝ I1 (layer 1)")
ax3 = ax[1].twinx()
ax3.plot(dp["t"], dp["VI"] / 1e3, color="C3", lw=0.7)
ax3.set_ylabel("V_I [kV] (red; unshielded q ×2, declared)", color="C3")
ax[1].set_xlabel("Rikitake time"); ax[1].set_ylabel("I1")
ax[1].set_title("(b) dynamo → cosmic-ray shielding → V_I")
save(fig, "fig8_layers_linked.png")
