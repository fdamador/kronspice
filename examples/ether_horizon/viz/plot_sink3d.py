# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""The sink electron (H49) in 3-D: the photon electron's light with a sink at its centre.

    go run ./examples/ether_horizon -viz examples/ether_horizon/viz   # writes sink_profile.csv
    python3 examples/ether_horizon/viz/plot_sink3d.py [OUT.gif]

Writes viz/sink3d.gif: one period of the photon electron's light, slowed to
3 s, with the H49 sink drawn at its centre:

  (a) 3-D: the (2, 1) path twice round the torus and the wave cos(ks − ωt)
      travelling along it, E (orange) and B (blue) across the path, taken
      from examples/photon_electron/viz/plot_structure.py (field.go's
      frame); at the centre the sink core (grey) that removes ether, the
      sphere where the inflow would reach the wave speed (red wireframe:
      the would-be horizon, ruled out by Earth's heat flow, README §6), and
      ether tracers falling in with speed ∝ 1/r² (r³ linear in time,
      constant volume flux);
  (b) Part 5's steady inflow into a localized loss (n = 0, r_s = 1 ξ) from
      viz/sink_profile.csv: density and Mach number against radius, at the
      critical loss γ_c (Mach 1 at the sink's edge) and at γ = 10 (a true
      horizon outside the sink, in loss-free ether); the dot rides the
      γ = 10 inflow.

Schematic: the core and the horizon are drawn at 40 and 58 fm so they can be
seen; a gravity-fitted sink's horizon would be below 10⁻¹⁴ m (README §7),
and the surviving sink has none. The ring geometry, frequency and the
plotted profile are the examples' own.
"""
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.animation import FuncAnimation, PillowWriter  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "photon_electron", "viz"))
from plot_structure import A, B_COL, E_COL, P_COL, R, arc, frame, knot  # noqa: E402

NF = 48
CORE, HORIZON, OUTER = 40.0, 58.0, 330.0  # fm, schematic
SINK_COL, HOR_COL, FLOW_COL = "#5f5e5a", "#e24b4a", "#378add"


def sphere(r, nu=24, nv=12):
    u, v = np.meshgrid(np.linspace(0, 2 * np.pi, nu), np.linspace(0, np.pi, nv))
    return r * np.cos(u) * np.sin(v), r * np.sin(u) * np.sin(v), r * np.cos(v)


def load_profile(path):
    with open(path) as f:
        rows = list(csv.DictReader(f))
    col = {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}
    return col


def main(out):
    prof = load_profile(os.path.join(HERE, "sink_profile.csv"))
    us = np.linspace(0, 4 * np.pi, 32, endpoint=False)
    ps, ts, n1, n2 = frame(us)
    s_us, Lp = arc(us)
    ppath, _ = knot(np.linspace(0, 4 * np.pi, 600))
    uu, vv = np.meshgrid(np.linspace(0, 2 * np.pi, 48), np.linspace(0, 2 * np.pi, 20))
    torus = ((R + A * np.cos(vv)) * np.cos(uu), (R + A * np.cos(vv)) * np.sin(uu), A * np.sin(vv))
    rng = np.random.default_rng(49)
    ntr = 70
    d = rng.normal(size=(ntr, 3))
    d /= np.linalg.norm(d, axis=1, keepdims=True)
    phase = rng.random(ntr)

    fig = plt.figure(figsize=(13.5, 5.6))
    fig.suptitle("H49: a sink at the centre of the photon electron — the light twice round the torus, "
                 "the ether falling in", fontsize=12.5)
    ax3 = fig.add_subplot(1, 2, 1, projection="3d")
    axp = fig.add_subplot(1, 2, 2)

    def draw(k):
        ph = 2 * np.pi * k / NF
        amp = np.cos(2 * np.pi * s_us / Lp - ph)
        ax3.cla()
        ax3.plot_surface(*torus, color="#d9d2e9", alpha=0.18, linewidth=0)
        ax3.plot(ppath[:, 0], ppath[:, 1], ppath[:, 2], color=P_COL, lw=1.6, alpha=0.7)
        e = (-amp)[:, None] * n2 * 75
        b = amp[:, None] * n1 * 75
        ax3.quiver(ps[:, 0], ps[:, 1], ps[:, 2], e[:, 0], e[:, 1], e[:, 2], color=E_COL, lw=1.4, arrow_length_ratio=0.3)
        ax3.quiver(ps[:, 0], ps[:, 1], ps[:, 2], b[:, 0], b[:, 1], b[:, 2], color=B_COL, lw=1.4, arrow_length_ratio=0.3)
        ax3.plot_surface(*sphere(CORE), color=SINK_COL, alpha=0.85, linewidth=0)
        ax3.plot_wireframe(*sphere(HORIZON, 18, 9), color=HOR_COL, lw=0.6, alpha=0.55)
        # Tracers: constant volume flux, so r³ falls linearly in time.
        s = (phase + k / NF) % 1.0
        r = (OUTER ** 3 - s * (OUTER ** 3 - CORE ** 3)) ** (1 / 3)
        p = d * r[:, None]
        ax3.scatter(p[:, 0], p[:, 1], p[:, 2], color=FLOW_COL, s=9 + 30 * (CORE / r) ** 1.5, alpha=0.8, depthshade=False)
        ax3.set_xlim(-330, 330); ax3.set_ylim(-330, 330); ax3.set_zlim(-200, 200)
        ax3.set_box_aspect((1, 1, 0.6))
        ax3.set_xticks([]); ax3.set_yticks([]); ax3.set_zticks([])
        ax3.view_init(elev=24, azim=-58 + 90 * k / NF)
        ax3.set_title(f"(a) light on the (2, 1) path, R = λ/4π = {R:.0f} fm; grey: sink core;\n"
                      "red: would-be horizon (ruled out); blue: ether falling in, v ∝ 1/r² (core, horizon not to scale)",
                      fontsize=9)

        axp.cla()
        x = prof["r_xi"]
        axp.axvspan(0, 1.0, color="#efefef")
        axp.text(0.5, 0.06, "sink\n(loss −iγψ)", ha="center", fontsize=8.5, color="#555", transform=axp.get_xaxis_transform())
        axp.plot(x, prof["rho_c"], color="#888", lw=1.4, ls="--", label="density, γ = γ_c")
        axp.plot(x, prof["rho_10"], color="#444", lw=1.6, label="density, γ = 10")
        axp.plot(x, prof["mach_c"], color=FLOW_COL, lw=1.4, ls="--", label="Mach, γ = γ_c")
        axp.plot(x, prof["mach_10"], color=FLOW_COL, lw=1.8, label="Mach, γ = 10")
        axp.axhline(1, color=HOR_COL, lw=0.9)
        hor = x[np.where((prof["mach_10"][:-1] - 1) * (prof["mach_10"][1:] - 1) < 0)[0][-1]]
        axp.axvline(hor, color=HOR_COL, lw=0.9, ls=":")
        axp.text(hor + 0.05, 2.75, f"horizon {hor:.2f} ξ", color=HOR_COL, fontsize=8.5)
        # Dot riding the γ = 10 inflow from 4 ξ to the sink edge, r³ linear in time.
        sr = (4.0 ** 2 - (k / NF) * (4.0 ** 2 - 1.0)) ** 0.5  # cylinder: r² linear in time
        axp.plot(sr, np.interp(sr, x, prof["mach_10"]), "o", color=FLOW_COL, ms=8)
        axp.set_xlim(0, 4); axp.set_ylim(0, 3)
        axp.set_xlabel("radius r [healing lengths ξ]"); axp.set_ylabel("density ρ/ρ∞,  Mach |v|/c_local")
        axp.legend(fontsize=8, loc="upper right")
        axp.grid(alpha=0.3)
        axp.set_title("(b) Part 5: steady inflow into a localized loss (n = 0, r_s = 1 ξ);\n"
                      "above γ_c the inflow passes Mach 1 outside the sink, a true horizon", fontsize=9)
        return []

    fig.text(0.5, 0.015, "Light: f = m_e c²/h = 1.24×10²⁰ Hz, slowed about 4×10²⁰× to 3 s (examples/photon_electron). "
             "Earth's heat flow (47 TW) forbids the horizon's Hawking radiation, so the sink that survives "
             "stays subsonic (README §6–7).", ha="center", fontsize=8.5, color="#444")
    fig.tight_layout(rect=(0, 0.06, 1, 0.93))
    FuncAnimation(fig, draw, frames=NF).save(out, writer=PillowWriter(fps=16), dpi=80)
    draw(NF // 4)
    fig.savefig(os.path.splitext(out)[0] + ".png", dpi=110)
    print("wrote", out)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "sink3d.gif"))
