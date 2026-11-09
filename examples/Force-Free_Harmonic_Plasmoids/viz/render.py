# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Pictures and an animation of Nachamkin's force-free plasmoid.

    go run ./examples/Force-Free_Harmonic_Plasmoids -viz examples/Force-Free_Harmonic_Plasmoids/viz
    python examples/Force-Free_Harmonic_Plasmoids/viz/render.py

Writes, next to this script:
  plasmoid_structure.png   field lines (eq. 8, l = 1), energy density, the energy cycle
  harmonics_and_size.png   the harmonic comb, lifetimes against size, the electrons seconds need
  kronspice_spectrum.png   the kronspice AC spectrum of the Appendix B plasmoid, ring-down
  plasmoid_cycle.gif       the plasmoid oscillating and fading, beside a clock from 1 ns to 30 s
The field is computed here from eq. 8; every other number comes from the
Go example's CSVs.
"""
import csv
import math
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter

HERE = os.path.dirname(os.path.abspath(__file__))
X1 = 4.493409457909064  # first zero of j1: αR


def read(name):
    with open(os.path.join(HERE, name), newline="") as f:
        rows = list(csv.DictReader(f))
    return {k: np.array([float(r[k]) if k != "kind" else 0 for r in rows]) if k != "kind" else [r[k] for r in rows] for k in rows[0]}


# ---------- eq. 8, l = 1, α = 1 (radius X1) ----------
def j0(x):
    return np.sin(x) / x


def j1(x):
    return np.sin(x) / x**2 - np.cos(x) / x


def field(p, sign=1.0):
    """Cartesian M at points p (..., 3), α = 1."""
    x, y, z = p[..., 0], p[..., 1], p[..., 2]
    r = np.sqrt(x * x + y * y + z * z) + 1e-12
    rho = np.sqrt(x * x + y * y) + 1e-12
    ct, st = z / r, rho / r
    cp, sp = x / rho, y / rho
    jx = np.where(r < 1e-3, 1 / 3, j1(r) / r)
    fr = 2 * jx * ct
    ft = -(j0(r) - jx) * st
    fp = sign * j1(r) * st
    return np.stack([fr * st * cp + ft * ct * cp - fp * sp, fr * st * sp + ft * ct * sp + fp * cp, fr * ct - ft * st], -1)


def trace(p0, n=2600, h=0.03):
    pts = [np.array(p0, float)]
    p = pts[0]
    for _ in range(n):
        def d(q):
            f = field(q)
            return f / (np.linalg.norm(f) + 1e-15)
        k1 = d(p); k2 = d(p + h / 2 * k1); k3 = d(p + h / 2 * k2); k4 = d(p + h * k3)
        p = p + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        pts.append(p)
    return np.array(pts)


AXIS = 0.6106 * X1  # the ring the tori nest around
STARTS = [0.40, 0.80, 0.95]  # three tori: inner (around the ring), middle, outer
CMAP = plt.get_cmap("plasma")
LINES = [trace([s * X1, 0, 0.0], n=4000) for s in STARTS]


def torus(L, nphi=48):
    """The torus a traced line lies on: its (ρ, z) points sorted by angle
    round the current ring, revolved about the axis (units of R)."""
    rho = np.hypot(L[:, 0], L[:, 1]) / X1
    z = L[:, 2] / X1
    ang = np.arctan2(z, rho - AXIS / X1)
    o = np.argsort(ang)
    rho, z = rho[o][::8], z[o][::8]
    rho, z = np.append(rho, rho[0]), np.append(z, z[0])
    ph = np.linspace(0, 2 * np.pi, nphi)
    return np.outer(rho, np.cos(ph)), np.outer(rho, np.sin(ph)), np.outer(z, np.ones_like(ph))


TORI = [torus(L) for L in LINES]
SHOW = [L[:900] for L in LINES]  # one stretch of field line per torus


def draw_plasmoid(ax, alpha=1.0, linecolor=None):
    for i, (X, Y, Z) in enumerate(TORI):
        c = CMAP(0.15 + 0.7 * i / (len(TORI) - 1))
        ax.plot_surface(X, Y, Z, color=c, alpha=0.12 * alpha, linewidth=0, shade=False)
        L = SHOW[i]
        ax.plot(L[:, 0] / X1, L[:, 1] / X1, L[:, 2] / X1, lw=1.0, color=linecolor or c, alpha=alpha)


def structure():
    fig = plt.figure(figsize=(16, 5.4))
    ax = fig.add_subplot(1, 3, 1, projection="3d")
    draw_plasmoid(ax)
    u, v = np.mgrid[0 : 2 * np.pi : 40j, 0 : np.pi : 20j]
    ax.plot_wireframe(np.cos(u) * np.sin(v), np.sin(u) * np.sin(v), np.cos(v), color="0.8", lw=0.3)
    ax.set_box_aspect((1, 1, 1))
    ax.set_xlim(-1, 1); ax.set_ylim(-1, 1); ax.set_zlim(-1, 1)
    ax.set_title("Eq. 8 (l = 1): three nested tori, one field line on each\nE and H both follow these lines (E ∥ H)")
    ax.set_xlabel("x/R"); ax.set_ylabel("y/R"); ax.set_zlabel("z/R")

    ax = fig.add_subplot(1, 3, 2)
    rr = np.linspace(1e-3, 1, 300)
    zz = np.linspace(-1, 1, 600)
    RHO, Z = np.meshgrid(rr, zz)
    P = np.stack([RHO * X1, np.zeros_like(RHO), Z * X1], -1)
    F = field(P)
    E2 = (F**2).sum(-1)
    r = np.sqrt(RHO**2 + Z**2)
    E2 = np.where(r <= 1, E2, np.nan)
    psi = np.where(r <= 1, r * X1 * j1(r * X1) * (RHO / np.maximum(r, 1e-9)) ** 2, np.nan)
    for sgn in (1, -1):
        im = ax.pcolormesh(sgn * RHO, Z, E2 / np.nanmax(E2), shading="auto", cmap="inferno")
        ax.contour(sgn * RHO, Z, psi, levels=10, colors="w", linewidths=0.6)
    th = np.linspace(0, 2 * np.pi, 400)
    ax.plot(np.cos(th), np.sin(th), "c", lw=1)
    ax.plot([AXIS / X1, -AXIS / X1], [0, 0], "c+", ms=10)
    ax.set_aspect("equal")
    ax.set_title("Cut through the axis: energy density |M|² (colour)\nand the tori (white); + = the current ring at 0.611 R")
    ax.set_xlabel("x/R"); ax.set_ylabel("z/R")
    fig.colorbar(im, ax=ax, shrink=0.8, label="|M|² ÷ max")

    ax = fig.add_subplot(1, 3, 3)
    t = np.linspace(0, 2, 400)
    for kOverA, ls in ((59.776, "-"), (1.5, "--")):
        we = np.sin(2 * np.pi * t) ** 2
        wk = (1 - 1 / kOverA**2) * np.cos(2 * np.pi * t) ** 2
        wm = (1 / kOverA**2) * np.cos(2 * np.pi * t) ** 2
        lab = f"k/α = {kOverA:g}" + (" (Appendix B)" if kOverA > 2 else "")
        ax.plot(t, we, "C3" + ls, label="electric field, " + lab)
        ax.plot(t, wk, "C0" + ls, label="electron motion, " + lab)
        ax.plot(t, wm, "C2" + ls, label="magnetic field, " + lab)
    ax.plot(t, np.ones_like(t), "k:", lw=1, label="total")
    ax.set_xlabel("time (field periods)")
    ax.set_ylabel("share of the plasmoid's energy")
    ax.set_title("Where the energy is over a cycle\n(electric ↔ electron motion + magnetic; ωp² + α²c² = ω²)")
    ax.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "plasmoid_structure.png"), dpi=130, bbox_inches="tight")
    plt.close(fig)


def harmonics_and_size():
    comb = read("comb.csv")
    life = read("lifetimes.csv")
    need = read("needed.csv")
    fig, axs = plt.subplots(1, 3, figsize=(17, 5.4))

    ax = axs[0]
    n = comb["n"][:12]
    ax.vlines(comb["f_Hz"][:12] / 1e9, 0, 1, color="C0")
    for i in range(len(n)):
        ax.text(comb["f_Hz"][i] / 1e9, 1.03, f"n={int(n[i])}\n{comb['tauTE'][i]*1e9:.1f} ns", ha="center", fontsize=7)
    f0 = X1 * 299792458 / 0.12 / (2 * np.pi) / 1e9
    ax.axvline(f0, color="C3", ls="--")
    ax.text(f0, 0.5, " critical\n frequency\n αc/2π", color="C3", fontsize=8)
    ax.set_ylim(0, 1.3); ax.set_yticks([])
    ax.set_xlabel("frequency (GHz)")
    ax.set_title("Harmonics allowed for R = 12 cm (mean ball-lightning size)\nj1(kR) = 0; free TE lifetime under each line")

    ax = axs[1]
    for wj, c in ((1e4, "C0"), (1e6, "C3")):
        m = life["W_J"] == wj
        R = life["R_m"][m] * 100
        ax.plot(R, life["tauTE_recipe"][m], c + "-", label=f"paper's recipe ξ = 0.01, {wj:.0e} J (electrons > c)")
        ax.plot(R, life["tauTE_nonrel"][m], c + "--", label=f"electrons ≤ 0.1 c, {wj:.0e} J")
    ax.axhspan(1, 10, color="gold", alpha=0.35, label="reported ball lightning: 1–10 s")
    ax.axvspan(5, 25, color="0.85", alpha=0.6, label="reported radius 5–25 cm")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("plasmoid radius R (cm)"); ax.set_ylabel("free radiative lifetime (s)")
    ax.set_title("How long a free plasmoid lasts (TE half; TM ≈ 2×)")
    ax.legend(fontsize=7, loc="lower right")
    ax.set_ylim(1e-8, 1e2)

    ax = axs[2]
    styles = {1e-6: ("C7", ":"), 1e-3: ("C2", ":"), 1: ("C1", "-"), 3: ("C3", "-"), 10: ("C4", "-")}
    for tau, (c, ls) in styles.items():
        m = need["tau_s"] == tau
        lab = f"to last {tau:g} s" if tau >= 1 else f"to last {tau*1e6:g} µs" if tau < 1e-3 else f"to last {tau*1e3:g} ms"
        ax.plot(need["R_m"][m] * 100, need["ne_over_loschmidt"][m], color=c, ls=ls, label=lab)
    ax.axhline(14.5, color="k", lw=1.2)
    ax.text(1.1, 17, "fully stripped sea-level air (14.5 electrons per molecule)", fontsize=7)
    ax.axhline(1, color="k", lw=0.8, ls="--")
    ax.text(1.1, 1.2, "one free electron per air molecule", fontsize=7)
    ax.axvspan(5, 25, color="0.85", alpha=0.6)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("plasmoid radius R (cm)"); ax.set_ylabel("free electrons needed ÷ air molecules (sea level)")
    ax.set_title("Electron density a free plasmoid needs to last τ\n(his eq. 49; independent of the energy)")
    ax.legend(fontsize=7, loc="upper right")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "harmonics_and_size.png"), dpi=130)
    plt.close(fig)


def spectrum():
    sp = read("spectrum.csv")
    md = read("modes.csv")
    rd = read("ringdown.csv")
    fig, (ax, ins) = plt.subplots(1, 2, figsize=(15, 5), gridspec_kw={"width_ratios": [3, 1]})
    f0 = 128.1e9
    ax.semilogy((sp["f_Hz"] - f0) / 1e6, sp["absZ_ohm_normalized"], lw=0.7)
    counts = {"TE": 0, "TM": 0}
    for k, f, q in zip(md["kind"], md["f_Hz"], md["Q"]):
        counts[k] += 1
        ax.annotate(f"{k}{counts[k]}\nQ {q:.2g}", ((f - f0) / 1e6, q * 1.4), ha="center", fontsize=8,
                    color="C3" if k == "TE" else "C2")
    ax.set_ylim(top=sp["absZ_ohm_normalized"].max() * 6)
    ax.set_xlabel("frequency − 128.1 GHz (MHz)")
    ax.set_ylabel("|Z| of the Foster chain (normalized Ω)")
    ax.set_title("kronspice AC: the Appendix B plasmoid's l = 1 resonances\n(TE1 is Nachamkin's 128.158 GHz mode)")
    ins.plot(rd["t_periods"], rd["v"], lw=0.6)
    ins.plot(rd["t_periods"], np.sqrt(2 * rd["stored"]), "C3", lw=1, label="envelope from the ledger's stored energy")
    ins.plot(rd["t_periods"], -np.sqrt(2 * rd["stored"]), "C3", lw=1)
    ins.set_xlabel("time (field periods)")
    ins.set_title("kronspice ring-down of one mode\n(Q drawn as 40 for visibility; the real\nQ 4.8×10⁵ empties it in 0.6 µs)", fontsize=9)
    ins.legend(fontsize=7, loc="upper right")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "kronspice_spectrum.png"), dpi=130)
    plt.close(fig)


def animation():
    # the case drawn: R = 12 cm, 1 MJ, electrons ≤ 0.1 c (the longest-lived row of map b)
    life = read("lifetimes.csv")
    m = life["W_J"] == 1e6
    i12 = np.argmin(abs(life["R_m"][m] - 0.12))
    tau_best = life["tauTE_nonrel"][m][i12]
    tau_rec = life["tauTE_recipe"][m][i12]
    cases = [("free plasmoid, paper's recipe (12 cm, 1 MJ)", tau_rec, "C0"),
             ("free plasmoid, electrons ≤ 0.1 c (12 cm, 1 MJ)", tau_best, "C3")]
    nf = 90
    tl = np.logspace(-9, math.log10(30), nf)
    fig = plt.figure(figsize=(12, 5.4))
    ax3 = fig.add_subplot(1, 2, 1, projection="3d")
    ax2 = fig.add_subplot(1, 2, 2)

    def draw(k):
        ax3.cla(); ax2.cla()
        t = tl[k]
        phase = 2 * np.pi * k / 9  # slowed: one field cycle per 9 frames
        we = math.sin(phase) ** 2  # electric share; the rest is electron motion
        left = math.exp(-t / tau_best)
        c = np.array(matplotlib.colors.to_rgb("orangered")) * we + np.array(matplotlib.colors.to_rgb("deepskyblue")) * (1 - we)
        draw_plasmoid(ax3, alpha=max(left, 0.03), linecolor=tuple(c))
        ax3.view_init(elev=22, azim=3 * k)
        ax3.set_xlim(-1, 1); ax3.set_ylim(-1, 1); ax3.set_zlim(-1, 1)
        ax3.set_box_aspect((1, 1, 1)); ax3.set_axis_off()
        ax3.set_title(f"orange: energy in the electric field ({we:.0%})\nblue: in the electrons' motion; brightness: energy left {left:.2g}", fontsize=9)
        for lab, tau, c in cases:
            ax2.plot(tl, np.exp(-tl / tau), color=c, label=f"{lab}: τ = {tau:.2g} s")
            ax2.plot([t], [math.exp(-t / tau)], "o", color=c)
        ax2.axvspan(1, 10, color="gold", alpha=0.35, label="reported ball lightning lasts 1–10 s")
        ax2.axvline(t, color="k", lw=0.8)
        ax2.set_xscale("log"); ax2.set_ylim(-0.03, 1.05); ax2.set_xlim(1e-9, 30)
        ax2.set_xlabel("time after formation (s), log scale"); ax2.set_ylabel("field energy left")
        ax2.set_title(f"t = {t:.2g} s  (the field oscillates at THz; shown slowed)")
        ax2.legend(fontsize=7, loc="lower left")
        ax2.text(1.5, 0.55, "reported\nball lightning", fontsize=8, ha="left")
        return []

    anim = FuncAnimation(fig, draw, frames=nf, blit=False)
    anim.save(os.path.join(HERE, "plasmoid_cycle.gif"), writer=PillowWriter(fps=10), dpi=70)
    plt.close(fig)


if __name__ == "__main__":
    structure()
    harmonics_and_size()
    spectrum()
    animation()
    print("wrote plasmoid_structure.png, harmonics_and_size.png, kronspice_spectrum.png, plasmoid_cycle.gif")
