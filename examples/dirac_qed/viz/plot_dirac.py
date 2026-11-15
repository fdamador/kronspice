# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Pictures for examples/dirac_qed, from the Go run's CSVs.

    go run ./examples/dirac_qed -viz examples/dirac_qed/viz
    python3 examples/dirac_qed/viz/plot_dirac.py

writes dirac_qed_graph.png (QED series, Landau levels and g = 2, the trapped
electron's dip in kronspice) and dirac_qed_hydrogen.png (the n = 1, 2 levels:
Bohr, Dirac, and the Lamb shift Dirac misses).
"""
import csv
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))


def rows(name):
    return list(csv.DictReader(open(os.path.join(HERE, name))))


def graph():
    fig, ax = plt.subplots(1, 3, figsize=(17, 5.2))
    # (a) QED series
    q = rows("qed_terms.csv")
    names = [r["term"].replace(" (", "\n(").replace(", computed", ",\ncomputed") for r in q]
    vals = [abs(float(r["value"])) for r in q]
    cols = ["#1f77b4" if "computed" in r["term"] else "#7f7f7f" for r in q]
    y = np.arange(len(vals))
    ax[0].barh(y, vals, color=cols)
    ax[0].set_xscale("log")
    ax[0].set_yticks(y)
    ax[0].set_yticklabels(names, fontsize=7)
    ax[0].invert_yaxis()
    ax[0].axvline(0.13e-12, color="m", ls=":", lw=1.4, label="measurement uncertainty 1.3×10⁻¹³")
    ax[0].axvline(3.4e-13, color="k", ls="--", lw=1, label="measured − theory (α from Rb) 3.4×10⁻¹³")
    ax[0].set_xlabel("|contribution to a_e = (g − 2)/2|")
    ax[0].set_title("(a) QED's series for the anomaly\n(blue: computed here; gray: literature)")
    ax[0].legend(fontsize=7, loc="lower right")
    # (b) Landau levels at b = 0.05
    L = rows("landau_levels.csv")
    for r in L:
        e = float(r["E"]) - 1
        n, up = int(r["n"]), r["up"] == "true"
        x0 = 0.1 if up else 1.1
        ax[1].plot([x0, x0 + 0.8], [e, e], color="#d62728" if up else "#1f77b4", lw=2.5)
        ax[1].text(x0 + 0.82, e, f"n={n}", fontsize=7, va="center")
    ax[1].set_xticks([0.5, 1.5])
    ax[1].set_xticklabels(["spin up", "spin down"])
    ax[1].set_ylabel("E − mc²  [mc²], field b = B/B_crit = 0.05")
    ax[1].set_title("(b) Dirac in a field: each spin-up level sits exactly\non the next spin-down level → g = 2")
    ax[1].set_xlim(0, 2.3)
    # (c) axial dip
    d = rows("axial_dip.csv")
    df = np.array([float(r["df_hz"]) for r in d])
    w = np.array([float(r["v_with"]) for r in d])
    wo = np.array([float(r["v_without"]) for r in d])
    ax[2].plot(df, wo, color="#7f7f7f", lw=1.5, label="detector tank alone")
    ax[2].plot(df, w, color="#1f77b4", lw=1.8, label="with one trapped electron")
    ax[2].set_xlabel("frequency − 114 MHz  [Hz]")
    ax[2].set_ylabel("|V| at the detector  [V per A]")
    ax[2].set_title("(c) kronspice AC: one electron as a series LC\nnotches the 114 MHz detector (Wineland–Dehmelt)")
    ax[2].legend(fontsize=8)
    fig.suptitle("dirac_qed: the standard electron, computed — spin ½ and g = 2 from Dirac, the anomaly from QED, the measurement as circuits",
                 fontsize=11)
    fig.tight_layout()
    out = os.path.join(HERE, "dirac_qed_graph.png")
    fig.savefig(out, dpi=125)
    plt.close(fig)
    print("wrote", out)


def hydrogen():
    H = {r["level"]: float(r["binding_eV"]) for r in rows("hydrogen.csv")}
    fig, ax = plt.subplots(1, 2, figsize=(13, 5.5), gridspec_kw={"width_ratios": [1, 1.4]})
    # left: n = 1, 2 overview (binding energies)
    for name, xs in (("1S½", 0), ("2S½", 0), ("2P½", 1), ("2P3/2", 2)):
        e = -H[name]
        ax[0].plot([xs, xs + 0.8], [e, e], lw=3, color="#1f77b4")
        ax[0].text(xs + 0.4, e + 0.3, name, ha="center", fontsize=9)
    ax[0].set_ylabel("energy [eV]")
    ax[0].set_xticks([])
    ax[0].set_title("Hydrogen from Dirac's equation\n(shooting = exact formula to double precision)")
    ax[0].set_ylim(-14.5, -2)
    # right: n = 2 zoom, Dirac vs with QED
    base = -H["2P½"]
    fs = (H["2P½"] - H["2P3/2"]) * 2.417989242e14 / 1e9  # eV → GHz
    ax[1].plot([0.0, 0.8], [0, 0], lw=3, color="#1f77b4")
    ax[1].text(0.4, 0.35, "2S½ = 2P½ (Dirac: exactly equal)", ha="center", fontsize=9)
    ax[1].plot([0.0, 0.8], [fs, fs], lw=3, color="#1f77b4")
    ax[1].text(0.4, fs + 0.35, f"2P3/2: +{fs:.3f} GHz", ha="center", fontsize=9)
    ax[1].plot([1.2, 2.0], [0, 0], lw=3, color="#d62728")
    ax[1].text(1.6, -0.6, "2P½", ha="center", fontsize=9)
    ax[1].plot([1.2, 2.0], [1.0578, 1.0578], lw=3, color="#d62728")
    ax[1].text(1.6, 1.4, "2S½: +1.0578 GHz (Lamb shift)", ha="center", fontsize=9)
    ax[1].plot([1.2, 2.0], [10.969, 10.969], lw=3, color="#d62728")
    ax[1].text(1.6, 11.3, "2P3/2: +10.969 GHz", ha="center", fontsize=9)
    ax[1].set_xticks([0.4, 1.6])
    ax[1].set_xticklabels(["Dirac (computed here)", "with QED (measured / QED theory)"])
    ax[1].set_ylabel("energy above 2P½  [GHz]")
    ax[1].set_ylim(-1.5, 13)
    ax[1].set_title("The n = 2 levels: Dirac gets the fine structure to 0.23 %,\nQED adds the Lamb shift and the rest")
    fig.tight_layout()
    out = os.path.join(HERE, "dirac_qed_hydrogen.png")
    fig.savefig(out, dpi=125)
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    graph()
    hydrogen()
