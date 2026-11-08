# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Part 9 chart (H42: two neutrinos, i.e. two rings): viz/pair.png.

    go run ./examples/fdamador_neutrino -only 9 -field -viz examples/fdamador_neutrino/viz
    python3 examples/fdamador_neutrino/viz/render_part9.py [DATA_DIR]

Panels: (1) normal-mode frequencies of a pair against distance, coaxial and
side by side; (2) beats: ring 0 started alone at d = a, the energy in each
ring and in both; (3) time-averaged force of the in-phase and antiphase
modes against distance (curve: analytic, dots: kronspice; crosses: random
relative phases after 30 periods); (4) Q of both modes against distance,
circuit and field solver; (5) chains: Q of the darkest and brightest mode
against N, with the N³ guide and the Q Part 2 needs; (6) the N = 16 chain's
darkest (band-edge) and brightest mode patterns.
"""
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = sys.argv[1] if len(sys.argv) > 1 else HERE
E_COL, B_COL, K_COL, G_COL = "#d9480f", "#1c64b4", "#222222", "#2b8a3e"


def read(name, strcols=("geom", "label")):
    p = os.path.join(DATA, name)
    if not os.path.exists(p):
        return None
    with open(p, encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    out = {}
    for k in rows[0]:
        if k in strcols:
            out[k] = np.array([r[k] for r in rows])
        else:
            out[k] = np.array([float(r[k]) for r in rows])
    return out


def meta():
    with open(os.path.join(DATA, "p9_meta.csv"), encoding="utf-8") as f:
        return {r["key"]: float(r["value"]) for r in csv.DictReader(f)}


def main():
    m = meta()
    cur = read("p9_curve.csv")
    pairs = read("p9_pairs.csv")
    fc = read("p9_forcecurve.csv")
    fk = read("p9_force.csv")
    ch = read("p9_chain.csv")
    cv = read("p9_chainvec.csv")
    fld = read("p9_field.csv")
    fs = m["f_single"] / 1e9

    fig, axs = plt.subplots(2, 3, figsize=(17, 10))
    fig.suptitle("Part 9 (H42): two rings (two model neutrinos) and chains of rings — coupled resonators", fontsize=14)

    # 1 splitting
    ax = axs[0, 0]
    for g, ls, name in (("coaxial", "-", "coaxial"), ("side", "--", "side by side")):
        s = cur["geom"] == g
        ax.plot(cur["d_over_a"][s], cur["f_in"][s] / 1e9, ls, color=B_COL, label=f"in phase, {name}")
        ax.plot(cur["d_over_a"][s], cur["f_anti"][s] / 1e9, ls, color=E_COL, label=f"antiphase, {name}")
    ax.axhline(fs, color=K_COL, lw=0.8, ls=":", label=f"one ring {fs:.4f} GHz")
    ax.set_xscale("log")
    ax.set_xlabel("centre distance d / a")
    ax.set_ylabel("mode frequency [GHz]")
    ax.legend(fontsize=8)
    ax.set_title("(1) normal-mode splitting (M > 0 coaxial, M < 0 side by side)", fontsize=10)

    # 2 beats
    ax = axs[0, 1]
    b = read("p9_beat0.csv")
    u0 = b["u"][0]
    ax.plot(b["t_periods"], b["e0"] / u0, color=E_COL, lw=1, label="ring 0 (started)")
    ax.plot(b["t_periods"], b["e1"] / u0, color=B_COL, lw=1, label="ring 1 (at rest)")
    ax.plot(b["t_periods"], b["u"] / u0, color=K_COL, lw=1.2, label="both (with the mutual M·I0·I1)")
    ax.set_xlim(0, 40)
    ax.set_xlabel("time [periods]")
    ax.set_ylabel("energy / starting energy")
    T = m["beat0_fmean"] / m["beat0_df"]
    ax.set_title(f"(2) beats, coaxial d = {m['beat0_d_over_a']:g}a: 1/Δf = {T:.2f} periods;\n"
                 f"the in-phase half leaves (Q {m['beat0_qin']:.0f}), the antiphase half stays (Q {m['beat0_qanti']:.0f})", fontsize=10)
    ax.legend(fontsize=8)

    # 3 force
    ax = axs[0, 2]
    for g, ls, name in (("coaxial", "-", "coaxial"), ("side", "--", "side by side")):
        s = fc["geom"] == g
        ax.plot(fc["d_over_a"][s], fc["f_in"][s], ls, color=B_COL, label=f"in phase, {name}")
        ax.plot(fc["d_over_a"][s], fc["f_anti"][s], ls, color=E_COL, label=f"antiphase, {name}")
    ax.plot(fk["d_over_a"], fk["f_in_k"], "o", color=B_COL, ms=4)
    ax.plot(fk["d_over_a"], fk["f_anti_k"], "o", color=E_COL, ms=4, label="kronspice (mode runs)")
    ax.errorbar(fk["d_over_a"], fk["f_rand_late"], yerr=fk["s_late"], fmt="x", color=G_COL, ms=6,
                label="random phases, after 30 periods")
    ax.errorbar(fk["d_over_a"], fk["f_rand_early"], yerr=fk["s_early"], fmt="+", color="#888888", ms=6,
                label="random phases, first period")
    ax.axhline(0, color=K_COL, lw=0.6)
    ax.set_xscale("log")
    ax.set_xlabel("d / a")
    ax.set_ylabel("time-averaged force [U/a]   (+ repels)")
    ax.legend(fontsize=7)
    ax.set_title("(3) force: only for correlated phases; random pairs\nrelax into the repelling dark mode", fontsize=10)

    # 4 Q vs d
    ax = axs[1, 0]
    for g, ls, name in (("coaxial", "-", "coaxial"), ("side", "--", "side by side")):
        s = cur["geom"] == g
        ax.plot(cur["d_over_a"][s], cur["q_in"][s], ls, color=B_COL, label=f"in phase, {name}")
        ax.plot(cur["d_over_a"][s], cur["q_anti"][s], ls, color=E_COL, label=f"antiphase, {name} (tan δ 1e-4)")
        ax.plot(cur["d_over_a"][s], cur["qrad_anti"][s], ls, color=E_COL, lw=0.7, alpha=0.6)
    s = pairs["geom"] == "coaxial"
    ax.plot(pairs["d_over_a"][s], pairs["q_anti_kron"][s], "o", color=E_COL, ms=4, label="kronspice")
    ax.plot(pairs["d_over_a"][s], pairs["q_in_kron"][s], "o", color=B_COL, ms=4)
    ax.axhline(m["q_single"], color=K_COL, lw=0.8, ls=":", label=f"one ring Q {m['q_single']:.1f}")
    if fld is not None:
        for lab, mk in (("in phase", "s"), ("antiphase", "D")):
            sel = np.array([l.startswith(lab) and "absorber" not in l for l in fld["label"]]) & (fld["cells"] == 2)
            if sel.any():
                ax.plot(fld["d_over_a"][sel], fld["q"][sel], mk, mfc="none", color=G_COL, ms=7,
                        label=f"field solver r/2, {lab} (lossless)")
        one = np.array([l == "one ring" for l in fld["label"]]) & (fld["cells"] == 2)
        if one.any():
            ax.axhline(fld["q"][one][0], color=G_COL, lw=0.8, ls=":", label=f"field: one ring Q {fld['q'][one][0]:.1f}")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("d / a")
    ax.set_ylabel("Q")
    ax.legend(fontsize=7)
    ax.set_title("(4) subradiance: the antiphase pair keeps its energy ~20× longer\n(thin line: radiation alone)", fontsize=10)

    # 5 chains
    ax = axs[1, 1]
    ax.loglog(ch["n"], ch["qrad_dark"], "-o", color=E_COL, ms=4, label="darkest mode, radiation only")
    ax.loglog(ch["n"], ch["q_dark"], "-", color=E_COL, lw=0.8, alpha=0.6, label="darkest, with tan δ = 1e-4")
    ax.loglog(ch["n"], ch["q_bright"], "-o", color=B_COL, ms=4, label="brightest mode")
    k = ch["q_dark_kron"] > 0
    ax.loglog(ch["n"][k], ch["q_dark_kron"][k], "s", mfc="none", color=K_COL, ms=8, label="kronspice")
    nn = np.logspace(0, np.log10(m["n_needed"]) + 0.2, 50)
    ax.loglog(nn, m["chain_coef"] * nn ** m["chain_exp"], ":", color=K_COL, lw=0.8,
              label=f"{m['chain_coef']:.1f}·N^{m['chain_exp']:.2f}")
    ax.axhline(m["q_needed"], color=G_COL, lw=1, ls="--", label=f"Q needed (Part 2) {m['q_needed']:.1e}")
    ax.plot([m["n_needed"]], [m["q_needed"]], "*", color=G_COL, ms=12)
    ax.text(m["n_needed"], m["q_needed"] / 30, f"N ≈ {m['n_needed']:.1e}", ha="center", fontsize=8, color=G_COL)
    ax.set_xlabel("rings in the chain N")
    ax.set_ylabel("Q")
    ax.legend(fontsize=7, loc="center left")
    ax.set_title(f"(5) chains (coaxial, d = {m['chain_d_over_a']:g}a): the dark mode compounds\nas a power of N, not exponentially", fontsize=10)

    # 6 chain mode patterns
    ax = axs[1, 2]
    s = cv["n"] == 16
    j = cv["j"][s]
    ax.bar(j - 0.2, cv["dark"][s], 0.4, color=E_COL, label="darkest (band edge)")
    ax.bar(j + 0.2, cv["bright"][s], 0.4, color=B_COL, label="brightest (all in phase)")
    ax.axhline(0, color=K_COL, lw=0.6)
    ax.set_xlabel("ring j (N = 16)")
    ax.set_ylabel("current amplitude (normalized)")
    ax.legend(fontsize=8)
    ax.set_title("(6) mode patterns: neighbours opposite under a smooth\nenvelope (net moment 0) against all in phase", fontsize=10)

    from matplotlib.ticker import FuncFormatter
    for ax in (axs[0, 0], axs[0, 2], axs[1, 0]):
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
        ax.xaxis.set_minor_formatter(FuncFormatter(lambda v, _: f"{v:g}" if v in (0.5, 2, 3, 5) else ""))
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = os.path.join(HERE, "pair.png")
    fig.savefig(out, dpi=110)
    print("wrote", out)


if __name__ == "__main__":
    main()
