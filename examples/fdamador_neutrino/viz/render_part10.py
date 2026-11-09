# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Part 10 chart (H45: the weak interaction as the ring's own field): viz/weak.png.

    go run ./examples/fdamador_neutrino -only 10 -viz examples/fdamador_neutrino/viz
    python3 examples/fdamador_neutrino/viz/render_part10.py [DATA_DIR]

Panels: (1, 2) electron recoil spectra at E = 862 keV (7Be) and 10 MeV:
weak (nu_e, nu_mu, anti-nu_e) against the ring's dipole term (m = 0 point
and with the ring form factor, m = +-1 upper estimate) and a point dipole at
the XENONnT bound; (3) total cross-section against E (T > 1 keV); (4) the
ring's moments at the Part 2 scale against the bounds and the Standard
Model; (5) SNO: the deuteron E1 photodisintegration cross-section used and
the ring form factor; (6) kronspice twin: ring current, power into the
electron's held moment and the accumulated work.
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


def read(name):
    with open(os.path.join(DATA, name), encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}


def meta():
    with open(os.path.join(DATA, "p10_meta.csv"), encoding="utf-8") as f:
        return {r["key"]: float(r["value"]) for r in csv.DictReader(f)}


def main():
    m = meta()
    rec = read("p10_recoil.csv")
    sig = read("p10_sigma.csv")
    ffd = read("p10_formfactor.csv")
    e1 = read("p10_e1.csv")
    tw = read("p10_twin.csv")
    fig, ax = plt.subplots(2, 3, figsize=(16, 9.5))

    for k, ee in enumerate([0.862, 10.0]):
        a = ax[0, k]
        s = np.isclose(rec["e_nu"], ee)
        t = rec["t"][s] * 1e3
        a.loglog(t, rec["weak_e"][s], color=B_COL, lw=2, label="weak ν_e (SM)")
        a.loglog(t, rec["weak_mu"][s], color=B_COL, lw=1.2, ls="--", label="weak ν_μ,τ")
        a.loglog(t, rec["weak_ebar"][s], color=B_COL, lw=1, ls=":", label="weak ν̄_e")
        a.loglog(t, rec["ring0_point"][s], color=E_COL, lw=1, ls="--", label="ring m = 0, point dipole")
        a.loglog(t, rec["ring0"][s], color=E_COL, lw=2, label="ring m = 0, with ring form factor")
        a.loglog(t, rec["ring1"][s], color=G_COL, lw=1.5, label="ring m = ±1 (upper estimate)")
        a.loglog(t, rec["dip_xenon"][s], color=K_COL, lw=1, ls="-.", label="point dipole at XENONnT bound")
        a.set_xlabel("electron recoil T [keV]")
        a.set_ylabel("dσ/dT [cm²/MeV]")
        a.set_title(f"({k + 1}) recoil spectra, E_ν = {ee:g} MeV")
        a.grid(alpha=0.3, which="both")
        if k == 0:
            a.legend(fontsize=7.5, loc="center left")

    a = ax[0, 2]
    a.loglog(sig["e_nu"], sig["weak_e"], color=B_COL, lw=2, label="weak ν_e")
    a.loglog(sig["e_nu"], sig["weak_mu"], color=B_COL, lw=1.2, ls="--", label="weak ν_μ,τ")
    a.loglog(sig["e_nu"], sig["weak_ebar"], color=B_COL, lw=1, ls=":", label="weak ν̄_e")
    a.loglog(sig["e_nu"], sig["ring0"], color=E_COL, lw=2, label="ring m = 0 (= ν̄)")
    a.loglog(sig["e_nu"], sig["ring1"], color=G_COL, lw=1.5, label="ring m = ±1 (upper)")
    a.set_xlabel("E_ν [MeV]")
    a.set_ylabel("σ (T > 1 keV) [cm²]")
    a.set_title("(3) total cross-section: weak ∝ E, ring flat")
    a.grid(alpha=0.3, which="both")
    a.legend(fontsize=8)

    a = ax[1, 0]
    labels = ["μ m=0\n(rms)", "|pc| m=±1\n(upper)", "μ_eff m=0\nat 10 keV", "XENONnT\nbound", "Borexino\nbound", "SM Dirac\n(0.1 eV)"]
    vals = [m["mu0_rms"], m["pc1_amp"], m["eff_xe0"], m["bound_xe"], m["bound_bx"], 3.2e-20]
    cols = [E_COL, G_COL, E_COL, K_COL, K_COL, B_COL]
    a.bar(range(len(vals)), vals, color=cols)
    a.set_yscale("log")
    a.set_xticks(range(len(vals)))
    a.set_xticklabels(labels, fontsize=8)
    a.set_ylabel("dipole moment [μB]")
    a.set_title("(4) ring moments (Part 2 scale) vs bounds")
    a.text(0.98, 0.70, f"ring charge 0, ⟨r²⟩ 0,\nanapole 0 (exact)\nSM ⟨r²⟩ ≈ 4e-33 cm²", transform=a.transAxes, fontsize=8, ha="right", va="top")
    a.grid(alpha=0.3, axis="y", which="both")

    a = ax[1, 1]
    a.plot(e1["w"], e1["sigma_e1"] * 1e27, color=K_COL, lw=2, label="γ d → n p, E1 (Bethe–Peierls + eff. range)")
    a.set_xlabel("photon energy ω [MeV]")
    a.set_ylabel("σ_E1 [mb]")
    a.set_title("(5) SNO: deuteron break-up through the photon")
    a.grid(alpha=0.3)
    b = a.twinx()
    x = ffd["x"]
    b.loglog(x * 197.327 / (m["a_nu"] * 1e15), ffd["ff2"], color=E_COL, lw=1.5, label="ring ⟨F²(qa)⟩ vs q [MeV]")
    b.set_ylabel("ring form factor ⟨F²⟩", color=E_COL)
    b.set_xlim(e1["w"][0], e1["w"][-1])
    b.set_xscale("linear")
    a.legend(fontsize=7.5, loc="upper right")
    b.legend(fontsize=7.5, loc="center right")
    a.text(0.03, 0.05, f"10 MeV: ring m=0 σ ≈ {m['sno_ring0']:.1e} cm²\nSM NC ~ {m['sig_nc_ref']:.0e} cm²; CC: ring 0",
           transform=a.transAxes, fontsize=8)

    a = ax[1, 2]
    a.plot(tw["t_periods"], tw["i_ring"] / np.max(np.abs(tw["i_ring"])), color=E_COL, lw=3, label="ring current / peak")
    pe = tw["p_electron"]
    a.plot(tw["t_periods"], pe / np.max(np.abs(pe)), color=B_COL, lw=0.8, alpha=0.7, label="power into held moment / peak")
    w = tw["work"]
    a.plot(tw["t_periods"], w / np.max(np.abs(w)), color=K_COL, lw=1.5, label="work on electron / its peak (= M·I_e·ΔI_ring, so it\nlies on the current: no net gain)")
    a.set_xlabel("time [periods]")
    a.set_title("(6) kronspice twin: no static ring–electron energy")
    a.set_xlim(0, 12)
    a.grid(alpha=0.3)
    a.legend(fontsize=7.5)

    fig.suptitle("Part 10 (H45): the ring's own field on the electron, against the weak interaction", fontsize=13)
    fig.tight_layout()
    out = os.path.join(HERE, "weak.png")
    fig.savefig(out, dpi=110)
    print("wrote", out)


if __name__ == "__main__":
    main()
