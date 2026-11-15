# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""ECE / ECE2 series RLC report figures from the run's data.

    go run ./examples/ece_rlc -data examples/ece_rlc/viz/data > examples/ece_rlc/run_output.txt
    python3 examples/ece_rlc/viz/plot_report.py examples/ece_rlc/viz/data examples/ece_rlc/viz

Writes ac.png (resonance curves), compare.png (ECE against ECE2 over ω₀ and
κ: closed forms with the run's points), tran.png (transient current and the
running energy ledger), and ../ece_rlc_summary.png (one-page picture of the
circuit, the two forms and their results).
"""
import csv
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

data, out = sys.argv[1], sys.argv[2]
PLAIN, ECE, ECE2 = "0.35", "#3b6ea5", "#c8553d"
TAGS = ["plain", "ece_t", "ece2_t", "ece_s", "ece2_s"]
NAMES = ["plain (Maxwell)", "ECE torsion", "ECE2 torsion", "ECE scalar", "ECE2 scalar"]
COL = [PLAIN, ECE, ECE2, ECE, ECE2]
R, L, C, c0 = 1.0, 10e-3, 84.6346e-6, 299792458.0
F0 = 1 / (2 * math.pi * math.sqrt(L * C))
W0 = 0.25 / (c0 * L)
KAPPA = 0.05


def load(name):
    with open(os.path.join(data, name)) as f:
        rows = list(csv.DictReader(f))
    return {k: [float(r[k]) for r in rows] for k in rows[0]}


ac, tr, sm = load("ac.csv"), load("tran.csv"), load("summary.csv")

# ---------- ac.png ----------
fig, ax = plt.subplots(1, 2, figsize=(13, 4.8), sharey=True)
for a, idx, title in ((ax[0], (0, 1, 2), "torsion term: c·ω₀·L = 0.25 Ω changes the damping"),
                      (ax[1], (0, 3, 4), "scalar term: κ = 0.05 on the capacitor retunes")):
    for n in idx:
        a.plot(ac["f"], ac["mag_" + TAGS[n]], color=COL[n], lw=2 if n else 1.5,
               ls="--" if n == 0 else "-", label=f"{NAMES[n]}: f_r {sm['f_peak'][n]:.2f} Hz, Q {sm['q'][n]:.2f}")
    a.set_title(title, fontsize=10)
    a.set_xlabel("frequency [Hz]")
    a.grid(alpha=0.3)
    a.legend(fontsize=8, loc="upper right")
ax[0].set_ylabel("|I| [A] for 1 V")
ax[1].set_xlim(150, 205)
ax[0].set_xlim(150, 196)
fig.suptitle("Resonance of the series RLC, plain and with the ECE / ECE2 terms (kronspice AC, 4001 points)", fontsize=12)
fig.tight_layout()
fig.savefig(os.path.join(out, "ac.png"), dpi=150)
plt.close(fig)

# ---------- compare.png ----------
fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
rt = [x / 100 for x in range(-40, 101)]  # c·ω₀·L [Ω]
a = ax[0]
q = lambda s, r: 2 * math.pi * F0 * L / (R + s * r)
a.plot(rt, [q(1, r) for r in rt], color=ECE, lw=2, label="ECE (s = +1): Q = ω0L/(R + cω₀L)")
r2 = [r for r in rt if R - 2 * r > 0.02]
a.plot(r2, [q(-2, r) for r in r2], color=ECE2, lw=2, label="ECE2 (s = −2): Q = ω0L/(R − 2cω₀L)")
a.axhline(q(0, 0), color=PLAIN, ls="--", lw=1, label="plain Q0 = 10.87")
a.axvline(0.5, color=ECE2, ls=":", lw=1)
a.text(0.52, 80, "ECE2 undamped\nat cω₀L = R/2", color=ECE2, fontsize=8)
for n in (1, 2):
    a.plot(0.25, sm["q"][n], "o", color=COL[n], ms=8, mec="k")
a.set_ylim(0, 100)
a.set_xlabel("c·ω₀·L [Ω]   (ω₀ = 8.34×10⁻⁸ m⁻¹ per 0.25 Ω here)")
a.set_ylabel("Q")
a.set_title("torsion: the same ω₀ damps in ECE and undamps in ECE2", fontsize=10)
a.legend(fontsize=8, loc="upper left")
a.grid(alpha=0.3)
a = ax[1]
ks = [x / 1000 for x in range(-150, 151)]
fr = lambda s, k: F0 / math.sqrt(1 + s * k)
a.plot(ks, [fr(1, k) for k in ks], color=ECE, lw=2, label="ECE: f_r = f0/√(1 + κ)")
a.plot(ks, [fr(-2, k) for k in ks], color=ECE2, lw=2, label="ECE2: f_r = f0/√(1 − 2κ)")
a.axhline(F0, color=PLAIN, ls="--", lw=1, label="plain f0 = 173.00 Hz")
for n in (3, 4):
    a.plot(KAPPA, sm["f_peak"][n], "o", color=COL[n], ms=8, mec="k")
a.set_xlabel("κ on the capacitor's top terminal")
a.set_ylabel("resonance f_r [Hz]")
a.set_title("scalar: the same κ tunes down in ECE and up in ECE2", fontsize=10)
a.legend(fontsize=8, loc="upper left")
a.grid(alpha=0.3)
fig.suptitle("ECE against ECE2: closed forms (lines) and the kronspice runs (dots)", fontsize=12)
fig.tight_layout()
fig.savefig(os.path.join(out, "compare.png"), dpi=150)
plt.close(fig)

# ---------- tran.png ----------
fig, ax = plt.subplots(2, 3, figsize=(14, 7.5), sharex=True)
t = [x * 1e3 for x in tr["t"]]
for j, n in enumerate((0, 1, 2)):
    a = ax[0, j]
    a.plot(t, tr["i_" + TAGS[n]], color=COL[n], lw=0.6)
    a.set_title(NAMES[n], fontsize=10)
    a.grid(alpha=0.3)
    a.set_ylim(-2.3, 2.3)
    a = ax[1, j]
    tag = TAGS[n]
    ein, eh, em, st = (tr[k + tag] for k in ("ein_", "eheat_", "emodel_", "stored_"))
    a.plot(t, [x * 1e3 for x in ein], color="k", lw=1.5, label="source in")
    a.plot(t, [x * 1e3 for x in eh], color="#e09f3e", lw=1.5, label="heat in R")
    a.plot(t, [x * 1e3 for x in em], color="#7b2cbf", lw=1.5, label="model-supplied (term)")
    a.plot(t, [x * 1e3 for x in st], color="#2a9d8f", lw=1.2, label="stored ½Li² + ½Cv²")
    bal = [(a_ + b_ - c_ - d_) * 1e3 for a_, b_, c_, d_ in zip(ein, em, eh, [s - st[0] for s in st])]
    a.plot(t, bal, color="r", lw=0.8, ls=":", label="in + model − heat − Δstored")
    a.set_xlabel("time [ms]")
    a.grid(alpha=0.3)
    a.legend(fontsize=7, loc="upper left")
ax[0, 0].set_ylabel("loop current [A]")
ax[1, 0].set_ylabel("running energy [mJ]")
fig.suptitle("Transient, 1 V sine at 173 Hz from rest (kronspice, trapezoidal, 20 µs): the ledger closes in every case",
             fontsize=12)
fig.tight_layout()
fig.savefig(os.path.join(out, "tran.png"), dpi=150)
plt.close(fig)

# ---------- ece_rlc_summary.png ----------
fig = plt.figure(figsize=(14, 8.2))
fig.suptitle("ECE and ECE2 spin-connection terms on one series RLC  —  kronspice H48 (hypothetical, off by default)",
             fontsize=14, weight="bold")

# schematic
a = fig.add_axes([0.02, 0.50, 0.40, 0.40])
a.axis("off")
a.set_xlim(0, 10)
a.set_ylim(0, 6)
lw = 2


def wire(xs, ys):
    a.plot(xs, ys, color="k", lw=lw)


wire([1, 1], [1, 2.2]); wire([1, 1], [3.8, 5]); wire([1, 2.4], [5, 5])
a.add_patch(plt.Circle((1, 3), 0.8, fill=False, lw=lw))
a.text(1, 3, "~", ha="center", va="center", fontsize=22)
a.text(-0.1, 3, "V1\n1 V", ha="right", va="center", fontsize=9)
a.add_patch(plt.Rectangle((2.4, 4.7), 1.6, 0.6, fill=False, lw=lw))
a.text(3.2, 5.6, "R1  1 Ω", ha="center", fontsize=9)
wire([4.0, 5.0], [5, 5])
for k in range(4):
    a.add_patch(matplotlib.patches.Arc((5.3 + 0.6 * k, 5), 0.6, 0.8, theta1=0, theta2=180, lw=lw))
a.text(6.2, 5.6, "L1  10 mH", ha="center", fontsize=9)
a.text(6.2, 4.2, "torsion:  + s·c·ω₀·L in series", ha="center", fontsize=8, color="#7b2cbf")
wire([7.4, 8.6], [5, 5]); wire([8.6, 8.6], [5, 3.3]); wire([8.6, 8.6], [2.7, 1])
wire([8.0, 9.2], [3.3, 3.3]); wire([8.0, 9.2], [2.7, 2.7])
a.text(9.4, 3.0, "C1  84.63 µF", va="center", fontsize=9)
a.text(9.4, 2.4, "scalar: body sees\n(1 + s·κ)·φ on top", va="top", fontsize=8, color="#7b2cbf")
wire([1, 8.6], [1, 1])
a.text(4.8, 0.3, "f0 = 173.00 Hz,  Q0 = 10.87", ha="center", fontsize=10)

# equation box
a = fig.add_axes([0.02, 0.06, 0.40, 0.40])
a.axis("off")
a.add_patch(FancyBboxPatch((0.02, 0.02), 0.96, 0.96, boxstyle="round,pad=0.01", fc="#f6f4ef", ec="0.6"))
lines = [
    ("ECE  (Eckardt eq. 4.211):", "E = −∇φ − ∂A/∂t − c·ω₀·A + ω·φ", ECE),
    ("ECE2 (Eckardt eq. 6.189):", "E = −∇φ − ∂A/∂t + 2(c·ω₀·A − ω·φ)", ECE2),
]
y = 0.86
for head, eq, col in lines:
    a.text(0.06, y, head, fontsize=10, weight="bold", color=col)
    a.text(0.10, y - 0.08, eq, fontsize=10, family="monospace")
    y -= 0.22
a.text(0.06, 0.40, "In the circuit (proposed mapping, circuit/ece.go), s = +1 ECE, −2 ECE2:", fontsize=9)
a.text(0.10, 0.31, "torsion   v = L·di/dt + s·c·ω₀·L·i", fontsize=9, family="monospace")
a.text(0.10, 0.23, "scalar    body sees (1 + s·κ_a)φ_a − (1 + s·κ_b)φ_b", fontsize=9, family="monospace")
a.text(0.06, 0.10, "Power of the terms booked as model-supplied, never as inflow.\n"
       "Standard prediction: ω₀ = ω = 0, all cases = plain.", fontsize=9, style="italic")

# Q bars
a = fig.add_axes([0.50, 0.56, 0.22, 0.32])
a.bar(range(5), sm["q"], color=COL, edgecolor="k")
a.axhline(sm["q"][0], color=PLAIN, ls="--", lw=1)
for n in range(5):
    a.text(n, sm["q"][n] + 0.4, f"{sm['q'][n]:.2f}", ha="center", fontsize=8)
a.set_xticks(range(5))
a.set_xticklabels(["plain", "ECE\ntorsion", "ECE2\ntorsion", "ECE\nscalar", "ECE2\nscalar"], fontsize=8)
a.set_title("Q", fontsize=11)
a.set_ylim(0, 25)

a = fig.add_axes([0.77, 0.56, 0.21, 0.32])
a.bar(range(5), sm["f_peak"], color=COL, edgecolor="k")
for n in range(5):
    a.text(n, sm["f_peak"][n] + 0.6, f"{sm['f_peak'][n]:.1f}", ha="center", fontsize=8)
a.set_ylim(160, 188)
a.set_xticks(range(5))
a.set_xticklabels(["plain", "ECE\ntorsion", "ECE2\ntorsion", "ECE\nscalar", "ECE2\nscalar"], fontsize=8)
a.set_title("resonance f_r [Hz]", fontsize=11)

# power at resonance
a = fig.add_axes([0.50, 0.10, 0.48, 0.34])
x = range(5)
wdt = 0.27
a.bar([i - wdt for i in x], sm["p_source"], wdt, color="k", label="from source")
a.bar(list(x), sm["p_heat"], wdt, color="#e09f3e", label="heat in R")
a.bar([i + wdt for i in x], sm["p_model"], wdt, color="#7b2cbf", label="model-supplied (the term)")
a.axhline(0, color="k", lw=0.8)
a.set_xticks(list(x))
a.set_xticklabels(NAMES, fontsize=9)
a.set_ylabel("power [W]")
a.legend(fontsize=8, loc="upper left")
a.set_title("Power at resonance: ECE2 torsion's extra watt comes from the term itself", fontsize=10)
a.grid(alpha=0.3, axis="y")
for i in x:
    if sm["p_model"][i] != 0:
        a.text(i + wdt, sm["p_model"][i] + (0.05 if sm["p_model"][i] > 0 else 0.05), f"{sm['p_model'][i]:+.2f}",
               ha="center", fontsize=8, color="#7b2cbf")
fig.text(0.99, 0.01, "engine = ordinary-element twin to 4×10⁻¹⁵ · AC ledger residual 0 · transients close to ≤ 2.4×10⁻⁶ J",
         ha="right", fontsize=8, color="0.3")
fig.savefig(os.path.join(out, "..", "ece_rlc_summary.png"), dpi=150)
plt.close(fig)
