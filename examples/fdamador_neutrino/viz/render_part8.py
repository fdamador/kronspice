# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Part 8 chart (H41: the flavours as three spin states): viz/triplet.png.

    go run ./examples/fdamador_neutrino -only 8 -field -viz examples/fdamador_neutrino/viz
    python3 examples/fdamador_neutrino/viz/render_part8.py

Panels: (1) Gram eigenvalues of the triplet m = +1, 0, −1 against Part 6's
types C, S, D; (2) the triplet's levels under an outside rotation Ω, with
the Ω that gives the measured Δm²31/Δm²21; (3) the reciprocal fit in
kronspice, started in e = m 0, against the exact mode expansion; (4) the
angular-momentum ledger of +1 → −1 under one defect: J of the field and
J(0) + ∫torque dt; (5) energy added per quantum by the outside rotation
against rim speed, Coriolis and centrifugal, bench and Part 2 rings;
(6) the dipoles each state carries in the field solver (per quantum).
"""
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
E_COL, B_COL, K_COL, G_COL = "#d9480f", "#1c64b4", "#222222", "#2b8a3e"
HBAR, EV, C0 = 1.054571817e-34, 1.602176634e-19, 299792458.0


def read(name):
    with open(os.path.join(HERE, name)) as f:
        rows = list(csv.DictReader(f))
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}


def meta():
    with open(os.path.join(HERE, "p8_meta.csv")) as f:
        return {r["key"]: float(r["value"]) for r in csv.DictReader(f)}


def main():
    m = meta()
    fig, axs = plt.subplots(2, 3, figsize=(16, 9.5))
    fig.suptitle("Part 8 (H41): the flavours as three spin states m = +1, 0, −1 of the ring", fontsize=14)

    # 1 Gram
    ax = axs[0, 0]
    gt = [m[f"gram_triplet_{a}"] for a in range(3)]
    gy = [m[f"gram_types_{a}"] for a in range(3)]
    x = np.arange(3)
    ax.bar(x - 0.18, gt, 0.34, color=B_COL, label="triplet m = +1, 0, −1")
    ax.bar(x + 0.18, gy, 0.34, color=E_COL, label="Part 6 types C, S, D")
    for xi, v in zip(x, gy):
        ax.text(xi + 0.18, v + 0.04, f"{v:.2f}", ha="center", fontsize=8)
    ax.set_xticks(x, ["λ₁", "λ₂", "λ₃"])
    ax.set_ylabel("Gram eigenvalue")
    ax.set_ylim(0, 2.4)
    ax.legend(fontsize=8)
    ax.set_title("(1) independence: the triplet spans 3 dimensions,\nthe types 2 (C = (S + iD)/√2)", fontsize=10)

    # 2 levels under rotation
    ax = axs[0, 1]
    w0, w1, om = m["w0"], m["w1"], m["om_ratio"]
    r = np.linspace(0, 0.03, 300)
    ax.plot(r, (w1 + r * w1) / w1, color=B_COL, lw=1.6, label="m = +1: ω₁ + Ω")
    ax.plot(r, (w1 - r * w1) / w1, color=B_COL, lw=1.6, ls="--", label="m = −1: ω₁ − Ω")
    ax.plot(r, np.full_like(r, w0 / w1), color=E_COL, lw=1.6, label="m = 0: ω₀ (no first-order shift)")
    ax.axvline(om / w1, color=K_COL, lw=1, ls=":")
    ax.text(om / w1 + 0.0005, 0.62, f"Δm²31/Δm²21 = {m['nufit_ratio']:.1f}\nat Ω/ω₁ = {om / w1:.2e}\n"
            f"rim {om / w1 * w1 * 0.01:.2g} m/s\n= {om * 0.01 / m['burst']:.0f}× burst", fontsize=8)
    ax.set_xlabel("Ω/ω₁ (outside rotation)")
    ax.set_ylabel("ω / ω₁ (lab frame)")
    ax.legend(fontsize=7.5, loc="center right")
    ax.set_title("(2) levels under outside rotation (field-solver ω₀, ω₁):\nthe ±1 pair splits, m = 0 stays", fontsize=10)

    # 3 kronspice run
    ax = axs[0, 2]
    tr = read("p8_triplet_run.csv")
    for k, lab, col in (("pe", "e (m = 0)", E_COL), ("pmu", "μ (m = +1)", B_COL), ("ptau", "τ (m = −1)", G_COL)):
        ax.plot(tr["t_periods"], tr[k], color=col, lw=1.2, label=lab + ", kronspice")
        ax.plot(tr["t_periods"], tr[k + "_full"], color=K_COL, lw=0.6, ls="--")
    ax.plot([], [], color=K_COL, lw=0.6, ls="--", label="exact mode expansion")
    ax.set_xlabel("time [periods]")
    ax.set_ylabel("share in the flavour state")
    ax.legend(fontsize=7.5, loc="upper right")
    ax.set_title(f"(3) reciprocal fit in kronspice (leak {100 * m['recip_leak']:.0f} %):\n"
                 f"θ23 = {m['recip_th23']:.1f}°, δ_CP = {m['recip_dcp']:.0f}° by μ–τ symmetry", fontsize=10)

    # 4 torque ledger
    ax = axs[1, 0]
    tq = read("p8_torque.csv")
    ax.plot(tq["t_beats"], tq["j"], color=B_COL, lw=1.8, label="J·ω/U₀ of the field")
    ax.plot(tq["t_beats"], tq["j"][0] + tq["torque_int"], color=E_COL, lw=1.0, ls="--", label="J(0) + ∫torque dt (from −∂H/∂α)")
    ax.plot(tq["t_beats"], tq["u"], color="#888", lw=0.8, label="field energy / U₀")
    ax.set_xlabel("time [beats]")
    ax.set_ylabel("J·ω₁/U₀")
    ax.legend(fontsize=7.5, loc="lower left")
    ax.set_title(f"(4) +1 → −1 under one defect: ΔJ = {m['torque_dj']:.3f}, ∫τ dt = {m['torque_int']:.3f}\n"
                 "(2ħ per quantum) go to the defect; a free ring keeps J", fontsize=10)

    # 5 rotation energy vs rim speed
    ax = axs[1, 1]
    v = np.logspace(-1, np.log10(C0), 400)
    rings = {}
    for i in range(12):
        k = f"rm{i}_a"
        if k not in m:
            break
        rings.setdefault(m[k], []).append(i)
    styles = ["-", "--"]
    for j, (a, ids) in enumerate(sorted(rings.items(), reverse=True)):
        i1 = [i for i in ids if m[f"rm{i}_m"] == 1][0]
        w, rho2 = m[f"rm{i1}_omega"], m[f"rm{i1}_rho2"]
        name = "bench ring a = 10 mm" if a > 1e-6 else f"Part 2 ring a = {a:.2g} m"
        ax.loglog(v, HBAR * v / a / EV, color=B_COL, ls=styles[j], lw=1.4, label=f"Coriolis ħΩ, {name}")
        ax.loglog(v, 0.5 * HBAR * w * rho2 * (v / C0) ** 2 / EV, color=E_COL, ls=styles[j], lw=1.4, label=f"centrifugal ½IΩ², {name}")
    ax.axhline(0.1, color=K_COL, lw=1, ls=":")
    ax.text(0.15, 0.13, "0.1 eV (declared)", fontsize=8)
    ax.axvline(m["burst"], color="#888", lw=1)
    ax.text(m["burst"] * 0.8, 1e-26, "burst\n1,881 m/s", fontsize=8, color="#666", ha="right")
    ax.axvline(C0, color="#888", lw=1, ls="--")
    ax.text(C0 * 0.6, 1e-26, "c", fontsize=9, color="#666", ha="right")
    ax.set_ylim(1e-34, 1e2)
    ax.set_xlabel("rim speed aΩ [m/s]")
    ax.set_ylabel("energy added per quantum [eV]")
    ax.legend(fontsize=6.8, loc="lower right")
    ax.set_title("(5) outside rotation: mħΩ (none for m = 0) and ½IΩ²;\nall of it work by the outside agent", fontsize=10)

    # 6 dipoles
    ax = axs[1, 2]
    try:
        fm = list(csv.DictReader(open(os.path.join(HERE, "p8_field_meta.csv"))))
    except FileNotFoundError:
        fm = []
    if fm:
        cells = max(float(r["cells"]) for r in fm)
        rows = [r for r in fm if float(r["cells"]) == cells]
        order = {"p1": 0, "0": 1, "m1": 2}
        rows.sort(key=lambda r: order[r["state"]])
        lab = ["m = +1", "m = 0", "m = −1"]
        x = np.arange(3)
        floor = 1e-6
        mz = [max(float(r["mz_uB"]), floor) for r in rows]
        mp = [max(float(r["mperp_uB"]), floor) for r in rows]
        pc = [max(float(r["p_uBc"]), floor) for r in rows]
        ax.bar(x - 0.25, mz, 0.24, color=B_COL, label="|m_z| (axial magnetic)")
        ax.bar(x, mp, 0.24, color="#7fa8d8", label="|m_⊥| (in-plane magnetic)")
        ax.bar(x + 0.25, pc, 0.24, color=E_COL, label="|p|·c (electric)")
        ax.set_yscale("log")
        ax.set_ylim(1e-6, 1e14)
        ax.set_xticks(x, lab)
        ax.set_ylabel("amplitude per quantum ħω [μB]")
        ax.legend(fontsize=7.5, loc="upper right")
        ax.set_title(f"(6) dipoles per quantum (field solver r/{cells:g}, bench ring);\n"
                     "time averages all 0 (Borexino static bound 2.8e-11 μB)", fontsize=10)
    else:
        ax.text(0.5, 0.5, "run with -field", ha="center")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = os.path.join(HERE, "triplet.png")
    fig.savefig(out, dpi=130)
    print("wrote", out)


if __name__ == "__main__":
    main()
