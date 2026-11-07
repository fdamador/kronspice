#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Picture of the RNA base pairs (A·U, G·C, G·U wobble): minimized geometry
and the frequency response of each pair's circuit.

    python prep/render_pairs.py . summary/base_pairs.png   (from examples/rna_vibration)
"""
import json, csv, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

base = sys.argv[1]; out = sys.argv[2]
atoms = json.load(open(base + "/summary/atoms.json"))
bonds = json.load(open(base + "/summary/bonds.json"))
def spec(m):
    r = list(csv.DictReader(open(base + "/spectra/%s.csv" % m.lower())))
    return (np.array([float(x["wavenumber_cm-1"]) for x in r]), np.array([float(x["absorption"]) for x in r]),
            np.array([float(x["release_one_quantum"]) for x in r]), np.array([float(x["release_thermal_300K"]) for x in r]))
modes = {m: list(csv.DictReader(open(base + "/modes/%s.csv" % m.lower()))) for m in ("AU", "GC", "GU")}

col = {"C": "#5b5b5b", "N": "#2f6fdb", "O": "#d93b3b", "H": "#d9d9d9"}
rad = {"C": 70, "N": 70, "O": 70, "H": 28}
panels = [("AU", "A–U", False), ("GC", "G–C", False), ("GU", "G·U wobble", False)]
fig = plt.figure(figsize=(13.5, 10), dpi=130)
fig.suptitle("RNA base pairs as circuits: geometry (AMBER RNA OL3, minimized from PDB 1EKA) and infrared light in and out", fontsize=13)
for k, (mol, label, flip) in enumerate(panels):
    A = [a for a in atoms if a["molecule"] == mol]
    X = np.array([[a["x"], a["y"], a["z"]] for a in A])
    X -= X.mean(0)
    # Plane of the pair: principal axes, longest axis along x.
    u, s, vt = np.linalg.svd(X)
    P = X @ vt.T
    # Put the first base (A or G) on the left for A–T/G–C, on the right for the flipped view.
    first = [i for i, a in enumerate(A) if a["name"].startswith(mol[0] + ":")]
    if (P[first, 0].mean() > 0) != flip:
        P[:, 0] *= -1
        P[:, 2] *= -1
    ax = fig.add_subplot(2, 3, k + 1, projection="3d")
    ax.view_init(elev=28, azim=-62)
    for b in bonds:
        if b["molecule"] != mol:
            continue
        p, q = P[b["i"]], P[b["j"]]
        if b["kind"] == "hbond":
            ax.plot(*zip(p, q), color="#e0a000", lw=2.2, ls=(0, (2, 2)), zorder=1)
        else:
            ax.plot(*zip(p, q), color="#444444", lw=1.6, zorder=1)
    for i, a in enumerate(A):
        ax.scatter(*P[i], s=rad[a["element"]], color=col[a["element"]], edgecolor="#222222", linewidth=0.4, depthshade=True, zorder=2)
    for i, a in enumerate(A):
        nm = a["name"].split(":")[1]
        if a["element"] in ("N", "O") or nm in ("C1p",):
            ax.text(*(P[i] + [0, 0, 0.35]), nm, fontsize=6.5, color="#222222")
    ax.set_title(label + " pair", fontsize=12, pad=0)
    lim = np.abs(P[:, :2]).max() * 0.6
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-lim, lim)
    ax.set_box_aspect((1, 1, 1))
    ax.set_axis_off()
    nb = len([b for b in bonds if b["molecule"] == mol and b["kind"] == "hbond"])
    ax.text2D(0.5, 0.02, "%d hydrogen bonds (dashed)" % nb, transform=ax.transAxes, ha="center", fontsize=8.5, color="#7a5a00")

    sx = fig.add_subplot(2, 3, k + 4)
    nu, a, q, t = spec(mol)
    sx.plot(nu, a, color="#2f6fdb", lw=1.1, label="absorbed (in)")
    sx.plot(nu, q, color="#d93b3b", lw=1.0, label="released, 1 quantum (out)")
    sx.plot(nu, t, color="#2a9d55", lw=1.0, label="released at 300 K (out)")
    for m in modes[mol]:
        if float(m["pairing_energy_share"]) >= 0.25:
            sx.axvline(float(m["wavenumber_cm-1"]), color="#e0a000", lw=0.8, ls="--")
    sx.set_xlim(0, 3700); sx.set_ylim(0, 1.08)
    sx.set_xlabel("wavenumber (cm⁻¹)", fontsize=9)
    if k == 0:
        sx.set_ylabel("relative power", fontsize=9)
    sx.tick_params(labelsize=8)
    sx.set_title(label + ": frequency response of the circuit", fontsize=10)
    if k == 0:
        sx.legend(fontsize=7.5, loc="upper center", frameon=False)
    ins = sx.inset_axes([0.12, 0.42, 0.3, 0.38])
    lo = nu <= 250
    ins.plot(nu[lo], a[lo] / max(a[lo].max(), 1e-30), color="#2f6fdb", lw=1.0)
    ins.plot(nu[lo], t[lo] / max(t[lo].max(), 1e-30), color="#2a9d55", lw=1.0)
    for m in modes[mol]:
        if float(m["pairing_energy_share"]) >= 0.25:
            ins.axvline(float(m["wavenumber_cm-1"]), color="#e0a000", lw=0.9, ls="--")
    ins.set_xlim(0, 250); ins.set_yticks([]); ins.tick_params(labelsize=6.5)
    ins.set_title("0–250 cm⁻¹: H-bond modes (dashed)", fontsize=6.8, color="#7a5a00")
handles = [Line2D([], [], marker="o", ls="", color=col[e], markeredgecolor="#222", label=n) for e, n in (("C", "carbon"), ("N", "nitrogen"), ("O", "oxygen"), ("H", "hydrogen"))]
handles.append(Line2D([], [], color="#e0a000", ls="--", lw=2, label="hydrogen bond"))
fig.legend(handles=handles, loc="upper right", fontsize=9, frameon=False, ncol=5, bbox_to_anchor=(0.98, 0.955))
fig.text(0.5, 0.005, "A·U and G·C are Watson–Crick pairs; G·U is the wobble pair (G N1–H···O2 U and U N3–H···O6 G). "
         "Methyl groups stand in for the riboses. Harmonic model, gas phase; lines drawn 6 cm⁻¹ wide.", ha="center", fontsize=8.5, color="#444")
fig.tight_layout(rect=(0, 0.02, 1, 0.95))
fig.savefig(out)
print("ok")
