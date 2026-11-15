# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Part 12 pictures: the donut void (the photon in the hollow core of an ether vortex ring).

    go run ./examples/ether_horizon -viz examples/ether_horizon/viz   # writes viz/donut_scan.csv
    python3 examples/ether_horizon/viz/plot_donut.py                  # donut_g.png, donut_void.png,
                                                                      # donut_void.mp4, donut_void.gif

donut_g.png is drawn from donut_scan.csv (the Go run). The image and the
animation are schematic: shape a/R = 0.5, the TE11 mode across the core, the
wave turning over once per circuit (so it repeats after two), the ether
swirling round the core.
"""
import csv
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import animation
from scipy.special import j1, jvp

HERE = os.path.dirname(os.path.abspath(__file__))
AE = 0.00115965218059
SCHW = 7.2973525643e-3 / (2 * np.pi)
X11 = 1.8411837813
COL = {0.0: "#1f77b4", 1e-4: "#2ca02c", 1e-2: "#ff7f0e", 1.0: "#d62728"}


def load():
    rows = list(csv.DictReader(open(os.path.join(HERE, "donut_scan.csv"))))
    data = {}
    for r in rows:
        data.setdefault(float(r["pi"]), []).append({k: float(v) for k, v in r.items()})
    return data


def graph(data):
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
    for pi, rs in data.items():
        t = [r["t"] for r in rs]
        lab = "Π = 0 (swirl only)" if pi == 0 else f"Π = {pi:g}"
        ax[0].plot(t, [r["g"] for r in rs], color=COL.get(pi), lw=2, label=lab)
        ax[2].plot(t, [r["M"] for r in rs], color=COL.get(pi), lw=2, label=lab)
    ax[0].axhline(2, color="k", lw=0.8, ls="--")
    ax[0].axhline(2 * (1 + AE), color="m", lw=1, ls=":", label="measured g = 2.00232")
    ax[0].set(xscale="log", xlabel="core / ring  a/R", ylabel="g = Mc² / E_light",
              title="(a) g: the swirl alone gives 2; pressure pulls it down")
    ax[0].set_ylim(1.2, 2.05)
    ax[0].legend(fontsize=8, loc="lower left")

    rs = data[0.0]
    t = np.array([r["t"] for r in rs])
    a = np.array([r["g"] for r in rs]) / 2 - 1
    ax[1].plot(t, a * 1e3, color=COL[0.0], lw=2.2, label="charge on the core wall (Π = 0)")
    ax[1].axhline(AE * 1e3, color="m", ls=":", lw=1.2, label="measured 1.15965×10⁻³")
    ax[1].axhline(SCHW * 1e3, color="gray", ls="--", lw=0.8, label="Schwinger α/2π")
    ax[1].axvspan(0.6, 1.0, color="#eeeeee", zorder=0)
    ax[1].text(0.62, 0.15, "thin-ring\nformulas\nstrained", fontsize=8, color="#666")
    ax[1].plot([1.0], [a[-1] * 1e3], "o", color="k", ms=6)
    ax[1].annotate(f"the edge, a = R (horn torus):\n{a[-1]*1e3:.4f}×10⁻³ = {a[-1]/AE*100:.1f}% of measured",
                   (1.0, a[-1] * 1e3), (0.025, 0.72), arrowprops=dict(arrowstyle="->", color="k"), fontsize=9)
    ax[1].set_ylim(0, 1.3)
    ax[1].set(xscale="log", xlabel="core / ring  a/R", ylabel="(g − 2)/2  [×10⁻³]",
              title="(b) the anomaly from the charge's static field")
    ax[1].legend(fontsize=8, loc="upper left", framealpha=0.95)

    for pi, rs in data.items():
        ax[2].plot([rs[-1]["t"]], [rs[-1]["M"]], "o", color=COL.get(pi), ms=5)
    ax[2].annotate("the edge: the hole closes\n(the double loop needs it)", (1.0, data[0.0][-1]["M"]), (0.05, 1.0),
                   arrowprops=dict(arrowstyle="->", color="k"), fontsize=9)
    ax[2].set(xscale="log", yscale="log", xlabel="core / ring  a/R", ylabel="mass M  [E0]",
              title="(c) the mass falls as the core swells:\nno thin donut is held")
    ax[2].legend(fontsize=8)
    fig.suptitle("H49 Part 12: the donut void — the double-loop photon in the hollow core of an ether vortex ring",
                 fontsize=12)
    fig.tight_layout()
    out = os.path.join(HERE, "donut_g.png")
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("wrote", out)


# ---------------------------------------------------------------- geometry
R, A = 1.0, 0.5


def torus(u, v, a=A):
    return ((R + a * np.cos(v)) * np.cos(u), (R + a * np.cos(v)) * np.sin(u), a * np.sin(v))


def te11(rho, phi, turn):
    """|E| of TE11 across the core, its pattern turned by angle `turn`."""
    k = X11 / A
    x = k * rho
    er = np.where(x > 1e-9, j1(x) / np.maximum(x, 1e-9), 0.5) * np.sin(phi - turn)
    ep = jvp(1, x) * np.cos(phi - turn)
    return np.hypot(er, ep)


def image():
    fig = plt.figure(figsize=(18, 6.6))
    ax = fig.add_subplot(1, 3, 1, projection="3d")
    draw3d(ax, phase=0.6)
    ax.set_title("The donut void: hollow vortex ring, photon in its core\n"
                 "(orange: the light's wave; arrows: its field turning over once per circuit;\n"
                 "blue: ether swirling round the core)", fontsize=9)

    bx = fig.add_subplot(1, 3, 2)
    n = 500
    xs = np.linspace(-2.1, 2.1, n)
    X, Z = np.meshgrid(xs, xs)
    # cross-section through the axis: two cores at x = ±R
    dens = np.ones_like(X)
    field = np.zeros_like(X)
    for c in (-R, R):
        d = np.hypot(X - c, Z)
        dens *= np.tanh(np.maximum(d - A, 0) / 0.12) ** 2
        ph = np.arctan2(Z, X - c)
        f = np.where(d < A, te11(d, ph, 0.0), 0)
        field = np.maximum(field, f / 0.58)
    bx.imshow(dens, extent=[-2.1, 2.1, -2.1, 2.1], origin="lower", cmap="Blues", vmin=0, vmax=1.4, alpha=0.9)
    bx.imshow(np.ma.masked_where(field < 0.02, field), extent=[-2.1, 2.1, -2.1, 2.1], origin="lower",
              cmap="hot", vmin=0, vmax=1.1)
    # ether swirl round each core (opposite senses in the cut: one vortex ring)
    for c, sgn in ((-R, -1), (R, 1)):
        for rr in (0.7, 0.95, 1.2):
            th = np.linspace(0, 2 * np.pi, 200)
            bx.plot(c + rr * np.cos(th), rr * np.sin(th), color="#1f4f8f", lw=0.7, alpha=0.6)
            for k in range(6):
                q = 2 * np.pi * k / 6
                bx.annotate("", (c + rr * np.cos(q + sgn * 0.25), rr * np.sin(q + sgn * 0.25)),
                            (c + rr * np.cos(q), rr * np.sin(q)),
                            arrowprops=dict(arrowstyle="->", color="#1f4f8f", lw=0.8))
        circ = plt.Circle((c, 0), A, fill=False, color="k", lw=1.2, ls="--")
        bx.add_patch(circ)
    bx.axvline(0, color="gray", lw=0.6, ls=":")
    bx.text(0.04, 1.95, "axis", fontsize=8, color="gray")
    bx.text(-R, -A - 0.2, "empty core\n(TE₁₁ light)", ha="center", va="top", fontsize=8)
    bx.text(1.55, 1.75, "ether\n(density, blue)", fontsize=8, color="#1f4f8f")
    bx.set(xlim=(-2.1, 2.1), ylim=(-2.1, 2.1), aspect="equal", xlabel="x / R", ylabel="z / R",
           title="Cut through the axis (a/R = 0.5): light fills the empty core,\n"
                 "the ether swirls round it; pressure and swirl must hold the wall")
    cx = fig.add_subplot(1, 3, 3)
    th = np.linspace(0, 2 * np.pi, 300)
    for t, col, lab in ((0.2, "#9ecae1", "a/R = 0.2"), (0.5, "#4292c6", "a/R = 0.5"), (1.0, "#d62728", "a/R = 1: the edge")):
        for c in (-1, 1):
            cx.fill(c + t * np.cos(th), t * np.sin(th), color=col, alpha=0.25 if t < 1 else 0.18)
            cx.plot(c + t * np.cos(th), t * np.sin(th), color=col, lw=2 if t == 1 else 1.3,
                    label=lab if c == 1 else None)
    cx.annotate("", (0.05, -1.35), (0.6, -1.35), arrowprops=dict(arrowstyle="->", lw=1.5))
    cx.text(0.32, -1.5, "the core swells", ha="center", fontsize=9)
    cx.plot([0], [0], "o", color="#d62728", ms=7)
    cx.annotate("the hole closes here:\nno room left for the\ndouble loop to wind", (0, 0), (-0.95, 1.35),
                arrowprops=dict(arrowstyle="->", color="#d62728"), fontsize=9, color="#d62728")
    cx.axvline(0, color="gray", lw=0.6, ls=":")
    cx.set(xlim=(-2.2, 2.2), ylim=(-2.0, 2.0), aspect="equal", xlabel="x / R", ylabel="z / R",
           title="Why it sits on the edge: the energy falls as the core swells,\n"
                 "until a = R; at the horn torus the anomaly is 95.9 % of measured")
    cx.legend(fontsize=8, loc="lower right")
    fig.text(0.5, 0.01, "Schematic. Part 12 results: the swirl's line energy is the constant-force string (g = 2), "
             "the core swells to the edge a = R, where the anomaly is 95.9 % of measured with nothing fitted; R ≈ 3.8 ƛ_C.",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    out = os.path.join(HERE, "donut_void.png")
    fig.savefig(out, dpi=120)
    plt.close(fig)
    print("wrote", out)


def draw3d(ax, phase, view=25):
    """phase: fraction of the double loop (0…1 is two circuits)."""
    ax.cla()
    u = np.linspace(0, 2 * np.pi, 90)
    v = np.linspace(0, 2 * np.pi, 40)
    U, V = np.meshgrid(u, v)
    x, y, z = torus(U, V)
    # the wave round the ring: wavelength 2 circumferences, packet at the phase position
    s = 2 * 2 * np.pi * phase                      # angle travelled (two turns per phase 0→1)
    wave = np.cos(0.5 * U - 0.5 * s)               # one wavelength per two turns
    env = np.exp(-((np.angle(np.exp(1j * (U - s)))) ** 2) / 0.9)
    c = plt.cm.YlOrRd(0.15 + 0.8 * np.clip(env * (0.5 + 0.5 * wave ** 2), 0, 1))
    c[..., 3] = 0.12 + 0.45 * env
    ax.plot_surface(x, y, z, facecolors=c, linewidth=0, antialiased=False, shade=False)
    # field arrows at the packet: TE11 field across the core, turned by half the angle travelled
    for du in np.linspace(-0.9, 0.9, 7):
        uu = s + du
        turn = 0.5 * uu                            # turns over (π) once per circuit
        cx, cy = R * np.cos(uu), R * np.sin(uu)
        nr = np.array([np.cos(uu), np.sin(uu), 0.0])
        ez = np.array([0.0, 0.0, 1.0])
        d = np.cos(turn) * nr + np.sin(turn) * ez
        w = np.exp(-du ** 2 / 0.5) * 0.55
        ax.quiver(cx, cy, 0, *(w * d), color="#b30000", lw=2.2, arrow_length_ratio=0.3)
        ax.quiver(cx, cy, 0, *(-w * d), color="#b30000", lw=2.2, arrow_length_ratio=0.3)
    # ether tracers swirling round the core (poloidal), outside the wall
    k = np.arange(48)
    uu = 2 * np.pi * (k % 12) / 12 + 0.2
    rr = A + 0.18 + 0.12 * (k // 12)
    vv = 2 * np.pi * (k * 0.618 % 1) + 2 * np.pi * phase * 2 / (rr / A) ** 2
    tx, ty, tz = torus(uu, vv, rr)
    ax.scatter(tx, ty, tz, s=7, color="#1f4f8f", alpha=0.8, depthshade=False)
    ax.plot([0, 0], [0, 0], [-1.3, 1.3], color="gray", lw=0.8, ls=":")
    lim = 1.35
    ax.set(xlim=(-lim, lim), ylim=(-lim, lim), zlim=(-lim * 0.7, lim * 0.7))
    ax.set_box_aspect((1, 1, 0.7))
    ax.view_init(elev=view, azim=-60 + 40 * phase)
    ax.set_axis_off()


def movie(frames=160):
    fig = plt.figure(figsize=(7.2, 6.2))
    ax = fig.add_subplot(projection="3d")
    title = fig.suptitle("", fontsize=10)

    def frame(i):
        ph = i / frames
        draw3d(ax, ph)
        turns = 2 * ph
        title.set_text("H49 Part 12: the donut void — photon in the hollow core of an ether vortex ring\n"
                       f"circuit {turns:4.2f} of 2: the field (red) has turned {180 * turns:5.1f}° "
                       "— it repeats only after two circuits (spin ½)")
        return []

    anim = animation.FuncAnimation(fig, frame, frames=frames, interval=60)
    out_mp4 = os.path.join(HERE, "donut_void.mp4")
    anim.save(out_mp4, writer=animation.FFMpegWriter(fps=24, bitrate=2400))
    print("wrote", out_mp4)
    out_gif = os.path.join(HERE, "donut_void.gif")
    anim.save(out_gif, writer=animation.PillowWriter(fps=16), dpi=60)
    print("wrote", out_gif)
    plt.close(fig)


if __name__ == "__main__":
    what = sys.argv[1:] or ["graph", "image", "movie"]
    if "graph" in what:
        graph(load())
    if "image" in what:
        image()
    if "movie" in what:
        movie()
